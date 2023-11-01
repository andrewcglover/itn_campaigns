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