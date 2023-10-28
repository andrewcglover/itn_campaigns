# nets_received.R
# Functions generating the distribution of when nets were received
# These functions are used to smooth the empirical d

# kernal smoothing function

ksmth_fun <- function(DHS_for_MDC, AMP_for_MDC,
                      scaled_adm_camp_nets, scaled_adm_dist_nets,
                      CMC_series, CMC_first, CMC_last,
                      min_kde_mode, prop_max_kde_mdc, max_modes,
                      local_mode_window, peak_window_ratio, min_kde_int_mdc,
                      dhs_bw, dst_bw) {
  
  smth_dhs <- ksmooth(CMC_series, scaled_adm_camp_nets, 'normal', bandwidth = dhs_bw)
  smth_dist <- ksmooth(CMC_series, scaled_adm_dist_nets, 'normal', bandwidth = dst_bw)
  
  if (DHS_for_MDC) {
    if (AMP_for_MDC) {
      smth_comb <- sqrt(smth_dhs$y * smth_dist$y)
    } else {
      smth_comb <- smth_dhs$y
    }
  } else {
    if (AMP_for_MDC) {
      smth_comb <- smth_dist$y
    } else {
      stop("Error: no source data specified for estimating MDC timings")
    }
  }
  
  #smth_comb <- smth_comb / sum(smth_comb, na.rm = TRUE)
  
  kde_series <- smth_comb
  
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
      # if (kde_series[camp_net_modes_id[i]] > (lower_avg_kde)
      #     & kde_series[camp_net_modes_id[i]] > (upper_avg_kde)
      #     &(kde_series[camp_net_modes_id[i]] > (lower_avg_kde*peak_window_ratio)
      #       | kde_series[camp_net_modes_id[i]] > (upper_avg_kde*peak_window_ratio))) {
      #   peak_checked_ordered_id <- c(peak_checked_ordered_id, camp_net_modes_id[i])
      #   N_peaked <- N_peaked + 1
      # }
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
  
  kde_df <- data.frame("comb_net_series" = kde_series,
                       "smth_dhs" = smth_dhs,
                       "smth_dist" = smth_dist,
                       "selected_nodes" = selected_nodes)
  
  kde_list <- list(kde_df, selected_nodes_id)
  
  return(kde_list)
}




ksmth_fun_0 <- function(scaled_adm_camp_nets, scaled_adm_dist_nets,
                        CMC_series, CMC_first, CMC_last,
                        min_kde_mode, prop_max_kde_mdc, max_modes,
                        local_mode_window, peak_window_ratio, min_kde_int_mdc,
                        dhs_bw, dst_bw) {
  
  smth_dhs <- ksmooth(CMC_series, scaled_adm_camp_nets, 'normal', bandwidth = dhs_bw)
  smth_dist <- ksmooth(CMC_series, scaled_adm_dist_nets, 'normal', bandwidth = dst_bw)
  
  smth_comb <- sqrt(smth_dhs$y * smth_dist$y)
  #smth_comb <- smth_comb / sum(smth_comb, na.rm = TRUE)
  
  kde_series <- smth_comb
  
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
      # if (kde_series[camp_net_modes_id[i]] > (lower_avg_kde)
      #     & kde_series[camp_net_modes_id[i]] > (upper_avg_kde)
      #     &(kde_series[camp_net_modes_id[i]] > (lower_avg_kde*peak_window_ratio)
      #       | kde_series[camp_net_modes_id[i]] > (upper_avg_kde*peak_window_ratio))) {
      #   peak_checked_ordered_id <- c(peak_checked_ordered_id, camp_net_modes_id[i])
      #   N_peaked <- N_peaked + 1
      # }
    }
    
    if (N_peaked > 0) {
      selected_nodes_id <- peak_checked_ordered_id[1]
      if (N_peaked > 1) {
        for (j in 2:N_peaked) {
          if (sum(abs(selected_nodes_id - peak_checked_ordered_id[j])
                  < min_kde_int_mdc) == 0) {
            selected_nodes_id <- c(selected_nodes_id, peak_checked_ordered_id[j])
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
  
  kde_df <- data.frame("comb_net_series" = kde_series,
                       "smth_dhs" = smth_dhs,
                       "smth_dist" = smth_dist,
                       "selected_nodes" = selected_nodes)
  
  kde_list <- list(kde_df, selected_nodes_id)
  
  return(kde_list)
}