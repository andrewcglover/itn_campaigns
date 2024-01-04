#usage_access_fitting.R

#prepare data list for time-varying usage/access model in Stan
create_usage_access_list <- function(usage = TRUE) {
  base_list <- list(N_c = N_ISO2,
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
                    #u = net_data$used,
                    n = net_data$total,
                    s = net_data$source_rec,
                    z = net_data$camp_rec,
                    max_rho = max_rounds,
                    r_meas = MDC_matrix,
                    r_tau = MDC_tau_matrix,
                    #r_tau = matrix(3, N_a, max_rounds),
                    rho = net_data$MDC_round,
                    disc_rnge = 4,
                    N_disc = 20,
                    max_m = 72,
                    no_round = -9999)
  if (usage) {
    base_list$u <- net_data$used
    usage_list <<- base_list
  } else {
    base_list$u <- net_data$access
    access_list <<- base_list
  }
}

# Run Stan model

usage_access_stan_fit <- function(usage = TRUE) {
  if (usage) {
    usage_fit <<- stan('./scripts/stan/use_acc_reg_all_ccc.stan',
                       data = usage_list,
                       iter = 300,
                       warmup = 200,
                       chains = 1,
                       init_r = 1e-2
                       #control = list(adapt_delta = 0.99,
                       #stepsize = 0.5,
                       #max_treedepth = 15
                       #               )
                       )
  } else {
    access_fit <<- stan('./scripts/stan/use_acc_reg_all_ccc.stan',
                        data = access_list,
                        iter = 300,
                        warmup = 200,
                        chains = 4,
                        init_r = 1e-2
                        )
  }
}

