#usage_access_fitting.R





#prepare data list for time-varying usage/access model in Stan
create_usage_access_list <- function(usage = TRUE) {
  if (usage) {
    usage_list <<- list(N_c = N_ISO2,
                        N_a = N_areas,
                        N_t = N_CMC,
                        cc = area_link$CTRY,
                        mu_n = est_mean_used_netlife,
                        sigma_n = est_sd_used_netlife,
                        N = dim(net_data)[1],
                        c = net_data$CTRY,
                        a = net_data$area_id,
                        t = net_data$CMC,
                        m = net_data$post_MDC,
                        u = net_data$used,
                        n = net_data$total,
                        s = net_data$source_rec,
                        z = net_data$camp_rec,
                        max_rho = max_rounds,
                        r_meas = MDC_matrix,
                        r_tau = MDC_tau_matrix,
                        #r_tau = matrix(3, N_a, max_rounds),
                        rho = binomial_df$MDC_round,
                        disc_rnge = 4,
                        N_disc = 20,
                        max_m = 72,
                        no_round = -9999)
  }
}