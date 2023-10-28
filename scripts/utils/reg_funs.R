#reg_funs.R

### Regression Model Functions ###

















stan_decay_fit <- function(nets_weighted, adm_net_link) {
  
  N_a <- length(unique(nets_weighted$area))
  N_c <- N_ISO2
  
  net_decay_dat <- list(N = dim(nets_weighted)[1],
                        N_a = N_a,
                        N_c = N_c,
                        a = nets_weighted$area_id,
                        c = nets_weighted$CTRY,
                        cc = adm_net_link$CTRY,
                        m = nets_weighted$months_since_obtained)
  
  net_decay_fit <- stan('exp_model_d2.stan',
                        data = net_decay_dat,
                        iter = 1000,
                        warmup = 500,
                        chains = 4,
                        #init_r = 1e-2,
                        control = list(adapt_delta = 0.95))
  
  net_decay_samples <- extract(net_decay_fit)
  
  return(net_decay_samples)
  
}



