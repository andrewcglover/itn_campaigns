# mdc.R
# Functions for estimating mass distribution campaign timings

append_access_meanlife <- function(dataset) {
  area_ids <- unique(dataset$area_id)
  N_area_ids <- length(area_ids)
  dataset$prior_mean_access_meanlife <- rep(NA, dim(dataset)[1])
  for (i in 1:N_area_ids) {
    ids <- which(dataset$area_id == area_ids[i])
    dataset$prior_mean_access_meanlife[ids] <- prior_mean_access_meanlife[i]
  }
  return(dataset)
}

calculate_net_receipt_weights <- function(dataset) {
  # hv005 = dhs weighting
  # hml4 = months since net obtained
  dhs_weight <- dataset$hv005
  net_age <- dataset$hml4
  
  # Set error/not recorded values to NA for net ages
  net_age[which(net_age > 36)] <- NA
  
  growth_rate <- 1 / dataset$prior_mean_access_meanlife
  dataset$receipt_w <- dhs_weight * exp(net_age * growth_rate)
  return(dataset)
}

append_total_weights_by_interview_date <- function(dataset) {
  # net_data expected as dataset
  # date of interview hv008
  # Household sample weight hv005
  
  all_nets_in_area <- function(dataset, idx) {
    area_idx <<- dataset$area_id[idx]
    area_all_nets <<- all_net_data[which(all_net_data$area_id == area_idx),]
  }
  
  dataset$tot_dhs_weight <- rep(0, dim(dataset)[1])
  all_nets_in_area(dataset, 1)
  for (i in 1:dim(dataset)[1]) {
    if (i > 1) {
      if (dataset$area_id[i] != area_idx) {
        all_nets_in_area(dataset, i)
      }
    }
    all_surveyed_net_ids <- which(area_all_nets$hv008 == dataset$CMC[i])
    dataset$tot_dhs_weight[i] <- sum(area_all_nets$hv005[all_surveyed_net_ids])
  }
  return(dataset)
}

append_weight_window <- function(dataset) {
  
  if (dim(dataset)[1] != (N_areas * N_CMC)) {
    print("warning: dataset size mismatch")
    }

  max_window <- min(36, N_CMC)
  j_threshold <- N_CMC - max_window + 1
  k <- 1
  for (i in 1:N_areas) {
    if (dataset$area_id[k] != i) {print("warning: weight assignment mismatch")}
    area_tot_nets <- dataset[which(dataset$area_id == i),]
    for (j in 1:N_CMC) {
      if (j == 1) {
        run_weight <- sum(area_tot_nets$tot_dhs_weight[1:max_window])
      } else {
        run_weight <- run_weight - area_tot_nets$tot_dhs_weight[j-1]
        if (j <= j_threshold) {
          run_weight <- run_weight + area_tot_nets$tot_dhs_weight[j+max_window-1]
        }
      }
      dataset$weight_window[k] <- run_weight
      k <- k + 1
    }
  }
  return(dataset)
}

append_total_receipt_weights <- function(dataset) {
  dataset$tot_receipt_w <- rep(0, dim(dataset)[1])
  k <- 1
  for (i in 1:N_areas) {
    for (j in 1:N_CMC) {
      t <- CMC_series[j]
      ids <- which((all_net_data$area_id == i) &
                     (all_net_data$CMC_net_obtained == t))
      dataset$tot_receipt_w[k] <- sum(all_net_data$receipt_w[ids])
      k <- k + 1
    }
    print(paste0(i,"/",N_areas))
  }
  return(dataset)
}

append_adj_receipt_weight <- function(dataset) {
  dataset$adj_receipt_w <- dataset$tot_receipt_w / dataset$weight_window
  dataset$adj_receipt_w[is.na(dataset$adj_receipt_w)] <- 0
  dataset$adj_receipt_w[is.infinite(dataset$adj_receipt_w)] <- 0
  return(dataset)
}

combine_weights <- function(dataset, density_name) {
  
  dataset$urb_comb_w <- rep(NA, dim(dataset)[1])
  
  for (i in 1:N_ADM1) {
    
    # Identify rural and urban areas for admin i
    ISO2_ADM1 <- uni_ADM1_ISO2[i]
    rural_area <- paste(ISO2_ADM1, "rural", sep = " ")
    urban_area <- paste(ISO2_ADM1, "urban", sep = " ")
    rural_ids <- which(dataset$area == rural_area)
    urban_ids <- which(dataset$area == urban_area)
    
    # Calculate grand totals of dhs weights by region for weighted mean
    grand_rural_dhs_weight <- sum(dataset$tot_dhs_weight[rural_ids])
    grand_urban_dhs_weight <- sum(dataset$tot_dhs_weight[urban_ids])
    
    # Identify densities to be combined
    rural_density <- dataset[dataset$area == rural_area, density_name]
    urban_density <- dataset[dataset$area == urban_area, density_name]
    
    # Weighted mean
    comb_density <- ((rural_density * grand_rural_dhs_weight) +
                       (urban_density * grand_urban_dhs_weight)) /
                          (grand_rural_dhs_weight + grand_urban_dhs_weight)
    
    # Append combined density
    if (!identical(rural_ids, integer(0))) {
      dataset$urb_comb_w[rural_ids] <- comb_density
    }
    if (!identical(urban_ids, integer(0))) {
      dataset$urb_comb_w[urban_ids] <- comb_density
    }
    
  }
  
  return(dataset)
  
}

#-------------------------------------------------------------------------------
# Estimate MDC dates

estimate_MDC_timings <- function(dataset, net_density_name = NULL) {
  
  # mdc_nodes <- data.frame(area = character(),
  #                         area_id = integer(),
  #                         selected_mdc_nodes = integer())
  mdc_nodes <<- NULL
  
  dataset$smth_nets <- rep(NA, dim(dataset)[1])
  
  net_density <- dataset[, net_density_name]

  for(i in 1:N_areas) {
    area_ids <- which(dataset$area_id == i)
    area_net_density <- net_density[area_ids]
    kde_series <- ksmooth(CMC_series, area_net_density, 'normal',
                         bandwidth = ksmooth_bandwidth)
    dataset$smth_nets[area_ids] <- kde_series
    
    # Remaining code to the end of if N_modes > 0 statement copied from previous
    # ksmth_fun on 02/11/23
    kde_modes <- (kde_series > c(NA,kde_series[1:(length(kde_series)-1)])
                  & kde_series > c(kde_series[2:length(kde_series)],NA))
    max_kde_series <- max(kde_series)
    camp_net_modes <- (kde_modes & kde_series > max_kde_series * prop_max_kde_mdc)
    
    camp_net_modes_id <- which(camp_net_modes)
    camp_net_modes_kde <- kde_series[camp_net_modes_id]
    camp_net_modes_kde_order <- order(camp_net_modes_kde, decreasing = TRUE)
    N_modes <- length(camp_net_modes_kde)#min(c(length(camp_net_modes_kde), max_modes))
    
    selected_nodes <- rep(FALSE, N_CMC)
    
    if (N_modes > 0) {
      
      peak_checked_ordered_id <- c()
      N_peaked <- 0
      for (j in 1:N_modes) {
        i <- camp_net_modes_kde_order[j]
        lower_window <- max(1, camp_net_modes_id[i] - local_mode_window)
        lower_avg_kde <- mean(kde_series[lower_window:(camp_net_modes_id[i]-1)])
        upper_window <- min(N_CMC, camp_net_modes_id[i] + local_mode_window)
        upper_avg_kde <- mean(kde_series[(camp_net_modes_id[i]+1):upper_window])
        if ((kde_series[camp_net_modes_id[i]] > (lower_avg_kde*peak_window_ratio)
             & kde_series[camp_net_modes_id[i]] > (upper_avg_kde*peak_window_ratio))) {
          peak_checked_ordered_id <- c(peak_checked_ordered_id, camp_net_modes_id[i])
          N_peaked <- N_peaked + 1
        }
      }
      
      
      if (N_peaked > 0) {
        selected_nodes_id <- peak_checked_ordered_id[1]
        if (N_peaked > 1) {
          for (j in 2:N_peaked) {
            if (sum(abs(selected_nodes_id - peak_checked_ordered_id[j])
                    < min_kde_int_mdc) == 0) {
              selected_nodes_id <- c(selected_nodes_id, peak_checked_ordered_id[j])
            } else {
              node_included_after_adjustment <- FALSE
              for (k in 1:length(selected_nodes_id)) {
                if ( !node_included_after_adjustment &
                     (abs(peak_checked_ordered_id[j] - selected_nodes_id[k]) < min_kde_int_mdc) ) {
                  if (peak_checked_ordered_id[j] < selected_nodes_id[k]) {
                    #check adjacent points to mode here - CHECK
                    if ((CMC_series[peak_checked_ordered_id[j]] - min_kde_int_mdc) >= CMC_first) {
                      peak_checked_ordered_id[j] <- selected_nodes_id[k] - min_kde_int_mdc
                    }
                  } else {
                    if ((CMC_series[peak_checked_ordered_id[j]] + min_kde_int_mdc) <= CMC_last) {
                      peak_checked_ordered_id[j] <- selected_nodes_id[k] + min_kde_int_mdc
                    }
                  }
                  if (sum(abs(selected_nodes_id - peak_checked_ordered_id[j])
                          < min_kde_int_mdc) == 0) {
                    node_included_after_adjustment <- TRUE
                  }
                }
              }
              if (node_included_after_adjustment) {
                selected_nodes_id <- c(selected_nodes_id, peak_checked_ordered_id[j])
              }
            }
          }
        }
      }
      
      N_selected <- length(selected_nodes_id)
      k <- N_peaked
      while (N_selected > max_modes) {
        selected_nodes_id <- selected_nodes_id[which(selected_nodes_id != peak_checked_ordered_id[k])]
        N_selected <- N_selected - 1
        k <- k - 1
      }
      
      #selected_nodes <- rep(FALSE, N_CMC)
      selected_nodes[selected_nodes_id] <- TRUE
    }
    
    # Update global data frame summary of selected nodes
    N_sn <- length(selected_nodes)
    areas_mdc_nodes <- data.frame("ISO2" = rep(id_link$ISO2[i], N_sn),
                                  "ADM1" = rep(id_link$ADM1[i], N_sn),
                                  "area" = rep(id_link$area[i], N_sn),
                                  "area_id" = rep(id_link$new_area_id[i], N_sn),
                                  "selected_nodes_ids" = selected_nodes_id)
    mdc_nodes <<- rbind.data.frame(mdc_nodes, areas_mdc_nodes)
    
    # Add logical indicator of MDC to net data frame
    dataset$selected_nodes[area_ids] <- selected_nodes
    
  }
  
  return(dataset)
  
}