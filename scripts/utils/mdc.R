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
  # hv005 = dhs wing
  # hml4 = months since net obtained
  dhs_w <- dataset$hv005
  net_age <- dataset$hml4
  
  # Set error/not recorded values to NA for net ages
  net_age[which(net_age > 36)] <- NA
  
  growth_rate <- 1 / dataset$prior_mean_access_meanlife
  dataset$grw_w <- dhs_w * exp(net_age * growth_rate)
  return(dataset)
}

append_total_weights_by_interview_date <- function(dataset) {
  # net_data expected as dataset
  # date of interview hv008
  # Household sample w hv005
  
  all_nets_in_area <- function(dataset, idx) {
    area_idx <<- dataset$area_id[idx]
    area_all_nets <<- all_net_data[which(all_net_data$area_id == area_idx),]
  }
  
  dataset$tot_dhs_w <- rep(0, dim(dataset)[1])
  all_nets_in_area(dataset, 1)
  for (i in 1:dim(dataset)[1]) {
    if (i > 1) {
      if (dataset$area_id[i] != area_idx) {
        all_nets_in_area(dataset, i)
      }
    }
    all_surveyed_net_ids <- which(area_all_nets$hv008 == dataset$CMC[i])
    dataset$tot_dhs_w[i] <- sum(area_all_nets$hv005[all_surveyed_net_ids])
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
    if (dataset$area_id[k] != i) {print("warning: w assignment mismatch")}
    area_tot_nets <- dataset[which(dataset$area_id == i),]
    for (j in 1:N_CMC) {
      if (j == 1) {
        run_w <- sum(area_tot_nets$tot_dhs_w[1:max_window])
      } else {
        run_w <- run_w - area_tot_nets$tot_dhs_w[j - 1]
        if (j <= j_threshold) {
          run_w <- run_w + area_tot_nets$tot_dhs_w[j + max_window - 1]
        }
      }
      dataset$w_window[k] <- run_w
      k <- k + 1
    }
  }
  return(dataset)
}

append_total_receipt_weights <- function(dataset) {
  dataset$rcpt_grw_w <- rep(0, dim(dataset)[1])
  dataset$rcpt_dhs_w <- rep(0, dim(dataset)[1])
  k <- 1
  for (i in 1:N_areas) {
    for (j in 1:N_CMC) {
      t <- CMC_series[j]
      ids <- which((all_net_data$area_id == i) &
                     (all_net_data$CMC_net_obtained == t))
      dataset$rcpt_grw_w[k] <- sum(all_net_data$grw_w[ids])
      dataset$rcpt_dhs_w[k] <- sum(all_net_data$hv005[ids])
      k <- k + 1
    }
    print(paste0(i,"/",N_areas))
  }
  return(dataset)
}

append_adj_receipt_weights <- function(dataset) {
  scale_by_weight_window <- function(den) {
    scale_den <- den / dataset$w_window
    scale_den[is.na(scale_den) | is.infinite(scale_den)] <- 0
  }
  dataset$adj_rcpt_grw_w <- scale_by_weight_window(dataset$rcpt_grw_w)
  dataset$adj_rcpt_dhs_w <- scale_by_weight_window(dataset$rcpt_dhs_w)
  # dataset$adj_rcpt_grw_w <- dataset$rcpt_grw_w / dataset$w_window
  # dataset$adj_rcpt_grw_w[is.na(dataset$adj_rcpt_grw_w)] <- 0
  # dataset$adj_rcpt_grw_w[is.infinite(dataset$adj_rcpt_grw_w)] <- 0
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
    
    # Warning if no rural or urban densities identified
    if (identical(rural_ids, integer(0)) & identical(urban_ids, integer(0))) {
      print(paste0("warning: no rural or urban areas found for admin ", i,
                   ": ", ISO2_ADM1))
    }
    
    # Identify densities to be combined
    if (identical(rural_ids, integer(0))) {
      grand_rural_dhs_w <- 0
      rural_density <- rep(0, N_CMC)
    } else {
      grand_rural_dhs_w <- sum(dataset$tot_dhs_w[rural_ids])
      rural_density <- dataset[rural_ids, density_name]
    }
    if (identical(urban_ids, integer(0))) {
      grand_urban_dhs_w <- 0
      urban_density <- rep(0, N_CMC)
    } else {
      grand_urban_dhs_w <- sum(dataset$tot_dhs_w[urban_ids])
      urban_density <- dataset[urban_ids, density_name]
    }
    
    # wed mean
    comb_density <- ((rural_density * grand_rural_dhs_w) +
                       (urban_density * grand_urban_dhs_w)) /
                          (grand_rural_dhs_w + grand_urban_dhs_w)
    
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

mode_smoothing <- function(dataset, net_density_name = NULL) {
  
  # mdc_modes <- data.frame(area = character(),
  #                         area_id = integer(),
  #                         selected_mdc_modes = integer())
  mdc_modes <<- NULL
  
  dataset[, paste0("smth_", net_density_name)] <- rep(NA, dim(dataset)[1])
  
  net_density <- dataset[, net_density_name]
  
  #dataset$fit_nets <- net_density

  for(i in 1:N_areas) {
    area_ids <- which(dataset$area_id == i)
    area_net_density <- net_density[area_ids]
    k_out <- ksmooth(CMC_series, area_net_density, 'normal', 
                     bandwidth = ksmooth_bandwidth)
    kde_series <- k_out$y
    dataset[area_ids, paste0("smth_", net_density_name)] <- kde_series
    #dataset$smth_nets[area_ids] <- kde_series
    
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
    
    selected_modes <- rep(FALSE, N_CMC)
    
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
        selected_modes_id <- peak_checked_ordered_id[1]
        if (N_peaked > 1) {
          for (j in 2:N_peaked) {
            if (sum(abs(selected_modes_id - peak_checked_ordered_id[j])
                    < min_kde_int_mdc) == 0) {
              selected_modes_id <- c(selected_modes_id, peak_checked_ordered_id[j])
            } else {
              node_included_after_adjustment <- FALSE
              for (k in 1:length(selected_modes_id)) {
                if ( !node_included_after_adjustment &
                     (abs(peak_checked_ordered_id[j] - selected_modes_id[k]) < min_kde_int_mdc) ) {
                  if (peak_checked_ordered_id[j] < selected_modes_id[k]) {
                    #check adjacent points to mode here - CHECK
                    if ((CMC_series[peak_checked_ordered_id[j]] - min_kde_int_mdc) >= CMC_first) {
                      peak_checked_ordered_id[j] <- selected_modes_id[k] - min_kde_int_mdc
                    }
                  } else {
                    if ((CMC_series[peak_checked_ordered_id[j]] + min_kde_int_mdc) <= CMC_last) {
                      peak_checked_ordered_id[j] <- selected_modes_id[k] + min_kde_int_mdc
                    }
                  }
                  if (sum(abs(selected_modes_id - peak_checked_ordered_id[j])
                          < min_kde_int_mdc) == 0) {
                    node_included_after_adjustment <- TRUE
                  }
                }
              }
              if (node_included_after_adjustment) {
                selected_modes_id <- c(selected_modes_id, peak_checked_ordered_id[j])
              }
            }
          }
        }
      }
      
      N_selected <- length(selected_modes_id)
      k <- N_peaked
      while (N_selected > max_modes) {
        selected_modes_id <- selected_modes_id[which(selected_modes_id != peak_checked_ordered_id[k])]
        N_selected <- N_selected - 1
        k <- k - 1
      }
      
      #selected_modes <- rep(FALSE, N_CMC)
      selected_modes[selected_modes_id] <- TRUE
    }
    
    # Update global data frame summary of selected modes
    N_sn <- length(selected_modes_id)
    areas_mdc_modes <- data.frame("ISO2" = rep(id_link$ISO2[i], N_sn),
                                  "ADM1" = rep(id_link$ADM1[i], N_sn),
                                  "area" = rep(id_link$area[i], N_sn),
                                  "area_id" = rep(id_link$new_area_id[i], N_sn))
    areas_mdc_modes[,paste0("mode_ids_", net_density_name)] <- selected_modes_id
    mdc_modes <<- rbind.data.frame(mdc_modes, areas_mdc_modes)
    
    # Add logical indicator of MDC to net data frame
    dataset[area_ids, paste0("modes_", net_density_name)] <- selected_modes
    
  }
  
  return(dataset)
  
}

#-------------------------------------------------------------------------------
# Estimating MDC timing from reference data

identify_antimodes <- function(dataset, density_name) {
  
  # Density names
  smth_name <- paste0("smth_", density_name)
  mode_name <- paste0("modes_", density_name)
  antimode_name <- paste0("antimodes_", density_name)
  
  # Select all area densities and modes
  all_densities <- dataset[, smth_name]
  all_modes <- dataset[, mode_name]
  
  for (i in 1:N_areas){
    
    # Select area
    aid <- uni_area_ids[i]
    ids <- which(dataset$area_id == aid)
    
    if (length(ids) != N_CMC) {print("warning: unexpected country data size")}
    
    # Select density and modes for area i
    area_density <- all_densities[ids]
    area_modes <- all_modes[ids]
    mode_ids <- which(area_modes == TRUE)
    N_modes <- sum(area_modes)
    
    # Select antimodes - first and last time points also included
    area_antimodes <- rep(FALSE, N_CMC)
    area_antimodes[1] <- TRUE
    area_antimodes[N_CMC] <- TRUE
    if (N_modes > 1) {
      for (j in 2:N_modes) {
        ma <- mode_ids[j-1]
        mb <- mode_ids[j]
        period <- area_density[ma:mb]
        antimode_here <- min(period)
        prop_amid <- which(area_density == antimode_here)
        antimode_id_here <- prop_amid[which(prop_amid > ma & prop_amid < mb)]
        area_antimodes[antimode_id_here] <- TRUE
      }
    }
    
    # Combine antimode list into dataset
    if (length(area_antimodes) != N_CMC) {"warning: unexpected antimode length"}
    dataset[ids, antimode_name] <- area_antimodes
    
  }
  
  return(dataset)
  
}

# Function to select additional antimode around start of time series
additional_early_antimode <- function(dataset, density_name) {
  
  # Define variable names
  smth_name <- paste0("smth_", density_name)
  mode_name <- paste0("modes_", density_name)
  antimode_name <- paste0("antimodes_", density_name)
  
  # Define new dataset
  new_antimodes <- NULL
  
  for (i in 1:N_areas) {
    # Select area ids and smooth density
    aids <- which(dataset$area_id == i)
    area_data <- dataset[aids,]
    area_smth <- area_data[, smth_name]
    
    # Select first mode
    area_mode_ids <- which(area_data[, mode_name])
    first_mode_id <- min(area_mode_ids)
    first_mode_val <- area_smth[first_mode_id]
    
    # Select first antimode
    area_antimodes <- area_data[, antimode_name]
    area_antimode_ids <- which(area_data[, antimode_name])
    first_antimode_id <- min(area_antimode_ids)
    first_antimode_val <- area_smth[first_antimode_id]
    
    # Return warning if first antimode is not at first CMC value
    if(area_data[first_antimode_id, "CMC"] != CMC_first) {
      print(paste0("warning: first antimode mismatch with first CMC in area ",i))
    }
    
    # Candidate antimode
    early_smth <- area_smth[first_antimode_id:first_mode_id]
    candidate_val <- min(early_smth)
    candidate_id <- which(early_smth == candidate_val)
    
    # Decide whether to include
    avg_smth <- mean(area_smth)
    if ((first_antimode_val >= avg_smth * min_first_antimode_overall_prop) &
        (first_antimode_val >= candidate_val * min_first_antimode_min_ratio)) {
      area_antimodes[candidate_id] <- TRUE
    }
    
    # Update antimode vector
    if (length(area_antimodes) != N_CMC) {
      print("warning: unexpected antimode length")
    }
    new_antimodes <- c(new_antimodes, area_antimodes)
  }
  
  # Update antimodes in dataset and return
  dataset[, antimode_name] <- new_antimodes
  return(dataset)
}

# Function to deselect adjacent antimodes
deselect_adjacent_antimodes <- function(dataset, density_name) {
  
  antimode_name <- paste0("antimodes_", density_name)
  
  for (i in 1:N_areas) {
    ids <- which(dataset$area_id == i)
    for (j in 2:N_CMC) {
      a0 <- dataset[ids[j-1], antimode_name]
      a1 <- dataset[ids[j], antimode_name]
      if (a0 & a1) {dataset[ids[j-1], antimode_name] <- FALSE}
    }
  }
  
  return(dataset)
  
}

# Function to generate composite density
regression_based_compostie_density <- function(dataset,
                                       rec_name,
                                       ref_name,
                                       force_positive_gradient = FALSE,
                                       force_zero_intercept = FALSE,
                                       use_predefined_extreme_nets = FALSE) {
  
  composite_densities <- NULL
  
  for (i in 1:N_areas) {
    
    # Select reference and recorded densities for area
    area_data <- dataset[which(dataset$area_id == i),]
    rec_den <- area_data[, rec_name]
    ref_den <- area_data[, ref_name]
    
    # Identify extreme points
    min_id <- which(cumsum(rec_den) == 0) %>% max + 1
    max_id <- which(revcumsum(rec_den) == 0) %>% min - 1
    if (use_predefined_extreme_nets | is.infinite(min_id)) {
      min_rec <- c(CMC_first, extreme_nets$min_rec[i]) %>% max
      min_id <- which(CMC_series == min_rec)
    }
    if (use_predefined_extreme_nets | is.infinite(max_id)) {
      max_rec <- c(CMC_last, extreme_nets$max_rec[i]) %>% min
      max_id <- which(CMC_series == max_rec)
    }
    
    mid_rec <- rec_den[min_id:max_id]
    mid_ref <- ref_den[min_id:max_id]
    
    # Linear regression
    reg_df <- data.frame("x" = mid_ref, "y" = mid_rec)
    reg_mod <- lm(y ~ x, data = reg_df)
    intr <- reg_mod$coefficients[[1]]
    grad <- reg_mod$coefficients[[2]]
    if (grad < 0) {grad = 0} # force gradient to zero if negative
    
    # Scale lower and upper reference density
    a <- min_id - 1
    b <- max_id + 1
    if (a < 1) {lwr_ref <- NULL} else {lwr_ref <- grad*ref_den[1:a]+intr}
    if (b > N_CMC) {upr_ref <- NULL} else {upr_ref <- grad*ref_den[b:N_CMC]+intr}
    
    # Area composite density
    area_comp <- c(lwr_ref, mid_rec, upr_ref)

    # # Scale lower and upper reference density
    # a <- max(c(1, min_id - 1))
    # b <- min(c(max_id + 1, N_CMC))
    # lwr_ref <- intr + grad * ref_den[1:a]
    # upr_ref <- intr + grad * ref_den[b:N_CMC]
    # 
    # # Area composite density
    # area_comp <- NULL
    # if (a > 1) area_comp <- c(area_comp, lwr_ref)
    # area_comp <- c(area_comp, mid_rec)
    # if (b < N_CMC) area_comp <- c(area_comp, upr_ref)
    if (length(area_comp) != N_CMC) {
      print(paste0("warning: unexpected composite density length for area ", i))
    }
    
    # Append to composite density vector
    composite_densities <- c(composite_densities, area_comp)
    
  }
  
  # Combine with dataset
  dataset$comp_nets <- composite_densities
  return(dataset)
  
}

# Function to generate composite density
generate_compostie_density <- function(dataset,
                                       rec_name,
                                       ref_name,
                                       scale_from_means = TRUE,
                                       use_predefined_extreme_nets = FALSE) {
  
  composite_densities <- NULL
  
  for (i in 1:N_areas) {
    
    # Select reference and recorded densities for area
    area_data <- dataset[which(dataset$area_id == i),]
    rec_den <- area_data[, rec_name]
    ref_den <- area_data[, ref_name]
    
    # Identify extreme points
    min_id <- which(cumsum(rec_den) == 0) %>% max + 1
    max_id <- which(revcumsum(rec_den) == 0) %>% min - 1
    if (use_predefined_extreme_nets | is.infinite(min_id)) {
      min_rec <- c(CMC_first, extreme_nets$min_rec[i]) %>% max
      min_id <- which(CMC_series == min_rec)
    }
    if (use_predefined_extreme_nets | is.infinite(max_id)) {
      max_rec <- c(CMC_last, extreme_nets$max_rec[i]) %>% min
      max_id <- which(CMC_series == max_rec)
    }
    
    mid_rec <- rec_den[min_id:max_id]
    mid_ref <- ref_den[min_id:max_id]
    
    a <- min_id - 1
    b <- max_id + 1
    
    if (scale_from_means) {
      # Calculate scale factor
      scale_fac <- mean(mid_rec) / mean(mid_ref)
      # Scale lower and upper reference density
      if (a < 1) {lwr_ref <- NULL} else {lwr_ref <- scale_fac*ref_den[1:a]}
      if (b > N_CMC) {upr_ref <- NULL} else {upr_ref <- scale_fac*ref_den[b:N_CMC]}
    } else {
      # Linear regression
      reg_df <- data.frame("x" = mid_ref, "y" = mid_rec)
      reg_mod <- lm(y ~ x, data = reg_df)
      intr <- reg_mod$coefficients[[1]]
      grad <- reg_mod$coefficients[[2]]
      if (grad < 0) {grad = 0} # force gradient to zero if negative
      # Scale lower and upper reference density
      if (a < 1) {lwr_ref <- NULL} else {lwr_ref <- grad*ref_den[1:a]+intr}
      if (b > N_CMC) {upr_ref <- NULL} else {upr_ref <- grad*ref_den[b:N_CMC]+intr}
    }
    
    # Area composite density
    area_comp <- c(lwr_ref, mid_rec, upr_ref)
    if (length(area_comp) != N_CMC) {
      print(paste0("warning: unexpected composite density length for area ", i))
    }
    
    # Append to composite density vector
    composite_densities <- c(composite_densities, area_comp)
    
  }
  
  # Combine with dataset
  dataset$comp_nets <- composite_densities
  return(dataset)
  
}


# Function to estimate MDC timings
estimate_mdc_timings <- function(dataset, mdc_bounds_name, density_name) {
  
  all_MDCs <- NULL
  
  for (i in 1:N_areas) {
    
    # Declare MDC vector
    MDCs <- rep(FALSE, N_CMC)
    
    # Select area data, MDC bounds and density of interest
    area_data <- dataset[which(dataset$area_id == i),]
    density <- area_data[, density_name]
    bounds <- area_data[, mdc_bounds_name]
    bound_ids <- which(bounds)
    
    # MDCs
    N_MDCs <- sum(bounds) - 1
    
    for (j in 1:N_MDCs) {
      j0 <- bound_ids[j]
      j1 <- bound_ids[j + 1] - 1
      sub_density <- density[j0:j1]
      norm_sub <- sub_density / sum(sub_density)
      CMC_sub <- CMC_series[j0:j1]
      mean_mdc <- sum(norm_sub * CMC_sub) %>% round
      mdc_id <- which(CMC_series == mean_mdc)
      MDCs[mdc_id] <- TRUE
    }
    
    # t <- 1
    # 
    # for (j in 1:N_MDCs) {
    #   
    #   j0 <- bound_ids[j]
    #   j1 <- bound_ids[j + 1] - 1
    #   
    #   sub_density <- density[j0:j1]
    #   N_sub <- length(sub_density)
    #   norm_sub <- sub_density / sum(sub_density)
    #   
    #   mean_mdc <- 0
    #   
    #   for (k in 1:(N_sub-1)) {
    #     mean_mdc <- mean_mdc + norm_sub[k] * CMC_series[t]
    #     print(mean_mdc)
    #     t <- t + 1
    #   }
    #   
    #   mean_mdc %<>% round
    #   mdc_id <- which(CMC_series == mean_mdc)
    #   MDCs[mdc_id] <- TRUE
    #   
    # }
    
    all_MDCs <- c(all_MDCs, MDCs)
    
  }
  
  dataset$mdc <- all_MDCs
  return(dataset)
  
}


#-------------------------------------------------------------------------------

# Function to normalise densities
normalise_area_densities <- function(dataset,
                                     density_names,
                                     norm_over_net_rec_range = FALSE,
                                     time_unit = NULL) {
  
  calc_norm_fac <- function(t0 = CMC_first, tm = CMC_last) {
    
    time_range <- tm - t0
    
    # Normalisation factor depending on time units
    if (is.null(time_unit)) {
      nf = 1
    } else if (time_unit %in% c("month","months","Month","Months","m","M")) {
      nf = 1
    } else if (time_unit %in% c("year","years","Year","Years","y","Y")) {
      nf = 12
    } else {
      print("warning: unrecognised time units")
    }
    return(nf)
  }
  
  if (!norm_over_net_rec_range) {norm_fac <- calc_norm_fac()}
  
  # Loop over density names to be normalised
  for (j in 1:length(density_names)) {
    # Normalise densities for each area
    if (density_names[j] %in% colnames(dataset)) {
      raw_densities <- dataset[, density_names[j]]
      norm_densities <- rep(NA, dim(dataset)[1])
      for (i in 1:N_areas) {
        if (norm_over_net_rec_range) {
          norm_fac <- calc_norm_fac(extreme_nets$min_rec[i],
                                    extreme_nets$max_rec[i])
        }
        ids <- which(dataset$area_id == i)
        area_raw_density <- raw_densities[ids]
        area_sum_density <- sum(area_raw_density)
        area_norm_density <- norm_fac * area_raw_density / area_sum_density
        dataset[ids, paste0(density_names[j] , "_norm")] <- area_norm_density
      }
    } else {
      print(paste0("warning: column name ", density_name[j], " not found"))
    }
  }
  
  return(dataset)
  
}

# Function to estimate local MDCs between density minima
estimate_local_mdcs <- function