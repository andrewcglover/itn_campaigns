#usage_access_fitting.R

#prepare data list for time-varying usage/access model in Stan
create_usage_access_list <- function(usage = TRUE) {
  base_list <- list(N_c = N_ISO2,
                    N_a = N_areas,
                    N_t = N_CMC,
                    cc = area_link$CTRY,
                    #mu_n = est_mean_used_netlife,
                    #sigma_n = est_sd_used_netlife,
                    N = dim(net_data)[1],
                    c = net_data$CTRY,
                    a = net_data$area_id,
                    t = net_data$CMC,
                    #m = net_data$post_MDC,
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
                    N_disc = 9,
                    max_m = 72,
                    no_round = -9999,
                    N_bb = N_bb)
  if (usage) {
    base_list$mu_n <- est_mean_used_netlife
    base_list$sigma_n <- est_sd_used_netlife
    base_list$u <- net_data$used
    usage_list <<- base_list
  } else {
    base_list$mu_n <- est_mean_access_netlife
    base_list$sigma_n <- est_sd_access_netlife
    base_list$u <- net_data$access
    access_list <<- base_list
  }
}

# Run Stan model

usage_access_stan_fit <- function(usage = TRUE) {
  if (usage) {
    usage_fit <<- stan('./scripts/stan/copy1_use_acc_reg_all_ccc.stan',
                       data = usage_list,
                       iter = 250,
                       warmup = 200,
                       chains = 8,
                       #algorithm = 'HMC',
                       init_r = 0.01
                       # control = list(adapt_delta = 0.99,
                       #                stepsize = 0.5,
                       #                max_treedepth = 15
                       #                )
                       )
  } else {
    access_fit <<- stan('./scripts/stan/copy1_use_acc_reg_all_ccc.stan',
                        data = access_list,
                        iter = 250,
                        warmup = 200,
                        chains = 8,
                        #algorithm = 'HMC',
                        init_r = 0.01
                        # control = list(adapt_delta = 0.99,
                        #                stepsize = 0.5,
                        #                max_treedepth = 15
                        #                )
                        )
  }
}

usage_access_cmdstanr_fit <- function(usage = TRUE) {
  
  #ua_stan_file <- './scripts/stan/ua_reg_cmdstanr.stan'
  ua_stan_file <- './scripts/stan/ua_reg_cmdstanr_malsim_branch_recov.stan'
  ua_mod <- cmdstan_model(ua_stan_file)
  
  if (usage) {
    usage_fit_raw <<- ua_mod$sample(data = usage_list,
                                    seed = Ucmd_seed,
                                    init = Ucmd_init,
                                    chains = Ucmd_chains,
                                    parallel_chains = Ucmd_parallel_chains,
                                    iter_warmup = Ucmd_warmup,
                                    iter_sampling = Ucmd_sampling,
                                    refresh = Ucmd_refresh
    )
  } else {
    access_fit_raw <<- ua_mod$sample(data = access_list,
                                     seed = Acmd_seed,
                                     init = Acmd_init,
                                     chains = Acmd_chains,
                                     parallel_chains = Acmd_parallel_chains,
                                     iter_warmup = Acmd_warmup,
                                     iter_sampling = Acmd_sampling,
                                     refresh = Acmd_refresh
    )
  }
}

#-------------------------------------------------------------------------------
# Extract samples

extract_time_series_draws <- function(cmdstanr = TRUE,
                                      usage = TRUE,
                                      access = TRUE) {
  
  # Extract usage samples
  if (usage) {
    
    # Extract usage parameters
    if (cmdstanr) {
      usage_draws <- usage_fit_raw$draws(format = "draws_df")
      Pbb_u <<- usage_draws %>% select(starts_with("u_tilde[")) %>% divide_by(N_bb)
      P_u <<- usage_draws %>% select(starts_with("P["))
      P0_u <<- usage_draws %>% select(starts_with("P0["))
      D_u <<- usage_draws %>% select(starts_with("D["))
      C_u <<- usage_draws %>% select(starts_with("C["))
      C0_u <<- usage_draws %>% select(starts_with("C0["))
      invlam_u <<- usage_draws %>% select(starts_with("inv_lambda["))
    } else {
      extracted_usage <- extract(usage_fit)
      Pbb_u <<- extracted_usage$u_tilde / N_bb
      P_u <<- extracted_usage$P
      P0_u <<- extracted_usage$P0
      D_u <<- extracted_usage$D
      C_u <<- extracted_usage$C
      C0_u <<- extracted_usage$C0
      invlam_u <<- extracted_usage$inv_lambda
    }
    
    # Time repetitions of inverse lambda
    invlamrep_u <<- invlam_u[, rep(seq_len(ncol(invlam_u)), each = N_CMC)]
    
    # Conditional usage
    PC_u <<- C_u / P_u
    
    # Time repetitions of inverse lambda
    invlamrep_u <<- invlam_u[, rep(seq_len(ncol(invlam_u)), each = N_CMC)]
    
    # Mean retention
    lamrep_u <<- 1 / invlamrep_u
    ret_u <<- 1 / (lamrep_u * (1-D_u))
  }
  
  if (access) {
    
    # Extract access parameters
    if (cmdstanr) {
      access_draws <- access_fit_raw$draws(format = "draws_df")
      Pbb_a <<- access_draws %>% select(starts_with("u_tilde[")) %>% divide_by(N_bb)
      P_a <<- access_draws %>% select(starts_with("P["))
      P0_a <<- access_draws %>% select(starts_with("P0["))
      D_a <<- access_draws %>% select(starts_with("D["))
      C_a <<- access_draws %>% select(starts_with("C["))
      C0_a <<- access_draws %>% select(starts_with("C0["))
      invlam_a <<- access_draws %>% select(starts_with("inv_lambda["))
    } else {
      extracted_access <- extract(access_fit)
      Pbb_a <<- extracted_access$u_tilde / N_bb
      P_a <<- extracted_access$P
      P0_a <<- extracted_access$P0
      D_a <<- extracted_access$D
      C_a <<- extracted_access$C
      C0_a <<- extracted_access$C0
      invlam_a <<- extracted_access$inv_lambda
    }
    
    # Proportion of accessible nets from campaigns
    PC_a <<- C_a / P_a
    
    # Time repetitions of inverse lambda
    invlamrep_a <<- invlam_a[, rep(seq_len(ncol(invlam_a)), each = N_CMC)]
    
    # Mean retention
    lamrep_a <<- 1 / invlamrep_a
    ret_a <<- 1 / (lamrep_a * (1-D_a))
  }
  
}

append_time_series_stats <- function(dataset,
                                     cmdstanr = TRUE,
                                     usage = TRUE,
                                     access = TRUE,
                                     lower_CrI1 = 0.025,
                                     upper_CrI1 = 0.975,
                                     lower_CrI2 = 0.1,
                                     upper_CrI2 = 0.9,
                                     lower_CrI3 = 0.25,
                                     upper_CrI3 = 0.75) {
  
  if (usage) {
    
    # Calculate mean values
    Pbb_u_mean <- Pbb_u %>% apply(2, mean)
    P_u_mean <- P_u %>% apply(2, mean)
    P0_u_mean <- P0_u %>% apply(2, mean)
    D_u_mean <- D_u %>% apply(2, mean)
    C_u_mean <- C_u %>% apply(2, mean)
    PC_u_mean <- PC_u %>% apply(2, mean)
    invlam_u_mean <- invlamrep_u %>% apply(2, mean)
    ret_u_mean <- ret_u %>% apply(2, mean)
    
    # Calculate lower bounds of credible intervals
    Pbb_u_LB1 <- Pbb_u %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    Pbb_u_LB2 <- Pbb_u %>% apply(2, quantile, probs = lower_CrI2, na.rm=TRUE)
    Pbb_u_LB3 <- Pbb_u %>% apply(2, quantile, probs = lower_CrI3, na.rm=TRUE)
    P_u_LB1 <- P_u %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    P0_u_LB1 <- P0_u %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    D_u_LB1 <- D_u %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    C_u_LB1 <- C_u %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    PC_u_LB1 <- PC_u %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    ret_u_LB1 <- ret_u %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    invlam_u_LB1 <- invlamrep_u %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    
    # Calculate upper bounds of credible intervals
    Pbb_u_UB1 <- Pbb_u %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    Pbb_u_UB2 <- Pbb_u %>% apply(2, quantile, probs = upper_CrI2, na.rm=TRUE)
    Pbb_u_UB3 <- Pbb_u %>% apply(2, quantile, probs = upper_CrI3, na.rm=TRUE)
    P_u_UB1 <- P_u %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    P0_u_UB1 <- P0_u %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    D_u_UB1 <- D_u %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    C_u_UB1 <- C_u %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    PC_u_UB1 <- PC_u %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    ret_u_UB1 <- ret_u %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    invlam_u_UB1 <- invlamrep_u %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    
    # Append to dataset
    dataset <- data.frame(dataset,
                          "Pbb_u_mean" = Pbb_u_mean,
                          "Pbb_u_LB1" = Pbb_u_LB1,
                          "Pbb_u_UB1" = Pbb_u_UB1,
                          "Pbb_u_LB2" = Pbb_u_LB2,
                          "Pbb_u_UB2" = Pbb_u_UB2,
                          "Pbb_u_LB3" = Pbb_u_LB3,
                          "Pbb_u_UB3" = Pbb_u_UB3,
                          "P_u_mean" = P_u_mean,
                          "P_u_LB1" = P_u_LB1,
                          "P_u_UB1" = P_u_UB1,
                          "P0_u_mean" = P0_u_mean,
                          "P0_u_LB1" = P0_u_LB1,
                          "P0_u_UB1" = P0_u_UB1,
                          "D_u_mean" = D_u_mean,
                          "D_u_LB1" = D_u_LB1,
                          "D_u_UB1" = D_u_UB1,
                          "C_u_mean" = C_u_mean,
                          "C_u_LB1" = C_u_LB1,
                          "C_u_UB1" = C_u_UB1,
                          "PC_u_mean" = PC_u_mean,
                          "PC_u_LB1" = PC_u_LB1,
                          "PC_u_UB1" = PC_u_UB1,
                          "invlam_u_mean" = invlam_u_mean,
                          "invlam_u_LB1" = invlam_u_LB1,
                          "invlam_u_UB1" = invlam_u_UB1,
                          "ret_u_mean" = ret_u_mean,
                          "ret_u_LB1" = ret_u_LB1,
                          "ret_u_UB1" = ret_u_UB1)
  }
  
  # Extract access samples
  if (access) {
    
    # Calculate mean values
    Pbb_a_mean <- Pbb_a %>% apply(2, mean)
    P_a_mean <- P_a %>% apply(2, mean)
    P0_a_mean <- P0_a %>% apply(2, mean)
    D_a_mean <- D_a %>% apply(2, mean)
    C_a_mean <- C_a %>% apply(2, mean)
    PC_a_mean <- PC_a %>% apply(2, mean)
    invlam_a_mean <- invlamrep_a %>% apply(2, mean)
    ret_a_mean <- ret_a %>% apply(2, mean)
    
    # Calculate lower bounds of credible intervals
    Pbb_a_LB1 <- Pbb_a %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    Pbb_a_LB2 <- Pbb_a %>% apply(2, quantile, probs = lower_CrI2, na.rm=TRUE)
    Pbb_a_LB3 <- Pbb_a %>% apply(2, quantile, probs = lower_CrI3, na.rm=TRUE)
    P_a_LB1 <- P_a %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    P0_a_LB1 <- P0_a %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    D_a_LB1 <- D_a %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    C_a_LB1 <- C_a %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    PC_a_LB1 <- PC_a %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    ret_a_LB1 <- ret_a %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    invlam_a_LB1 <- invlamrep_a %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
    
    # Calculate upper bounds of credible intervals
    Pbb_a_UB1 <- Pbb_a %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    Pbb_a_UB2 <- Pbb_a %>% apply(2, quantile, probs = upper_CrI2, na.rm=TRUE)
    Pbb_a_UB3 <- Pbb_a %>% apply(2, quantile, probs = upper_CrI3, na.rm=TRUE)
    P_a_UB1 <- P_a %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    P0_a_UB1 <- P0_a %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    D_a_UB1 <- D_a %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    C_a_UB1 <- C_a %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    PC_a_UB1 <- PC_a %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    ret_a_UB1 <- ret_a %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    invlam_a_UB1 <- invlamrep_a %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
    
    # Append to dataset
    dataset <- data.frame(dataset,
                          "Pbb_a_mean" = Pbb_a_mean,
                          "Pbb_a_LB1" = Pbb_a_LB1,
                          "Pbb_a_UB1" = Pbb_a_UB1,
                          "Pbb_a_LB2" = Pbb_a_LB2,
                          "Pbb_a_UB2" = Pbb_a_UB2,
                          "Pbb_a_LB3" = Pbb_a_LB3,
                          "Pbb_a_UB3" = Pbb_a_UB3,
                          "P_a_mean" = P_a_mean,
                          "P_a_LB1" = P_a_LB1,
                          "P_a_UB1" = P_a_UB1,
                          "P0_a_mean" = P0_a_mean,
                          "P0_a_LB1" = P0_a_LB1,
                          "P0_a_UB1" = P0_a_UB1,
                          "D_a_mean" = D_a_mean,
                          "D_a_LB1" = D_a_LB1,
                          "D_a_UB1" = D_a_UB1,
                          "C_a_mean" = C_a_mean,
                          "C_a_LB1" = C_a_LB1,
                          "C_a_UB1" = C_a_UB1,
                          "PC_a_mean" = PC_a_mean,
                          "PC_a_LB1" = PC_a_LB1,
                          "PC_a_UB1" = PC_a_UB1,
                          "invlam_a_mean" = invlam_a_mean,
                          "invlam_a_LB1" = invlam_a_LB1,
                          "invlam_a_UB1" = invlam_a_UB1,
                          "ret_a_mean" = ret_a_mean,
                          "ret_a_LB1" = ret_a_LB1,
                          "ret_a_UB1" = ret_a_UB1)
  }
  
  # Extract conditional usage
  if (usage & access) {
    N_u_samples <- dim(Pbb_u)[1]
    N_a_samples <- dim(Pbb_a)[1]
    if (N_u_samples == N_a_samples) {
      
      P_condu <- P_u / P_a
      P_condu_mean <- P_condu %>% apply(2, mean)
      P_condu_LB1 <- P_condu %>% apply(2, quantile, probs = lower_CrI1, na.rm=TRUE)
      P_condu_UB1 <- P_condu %>% apply(2, quantile, probs = upper_CrI1, na.rm=TRUE)
      
      # Append to dataset
      dataset <- data.frame(dataset,
                            "P_condu_mean" = P_condu_mean,
                            "P_condu_LB1" = P_condu_LB1,
                            "P_condu_UB1" = P_condu_UB1)
    } else {
      print(paste0("Warning: Conditional usage not calculated due to different",
                   "sample sizes for usage and access. There are ", N_u_samples,
                   " usage samples but ", N_a_samples, " access samples."))
    }
  }
  
  # Return dataset
  return(dataset)
  
}


# append_time_series_fits deprecated in favour of extract_time_series_draws and
# append_time_series_stats
append_time_series_fits <- function(dataset,
                                    cmdstanr = TRUE,
                                    usage = TRUE,
                                    access = TRUE,
                                    lower_CrI1 = 0.025,
                                    upper_CrI1 = 0.975,
                                    lower_CrI2 = 0.1,
                                    upper_CrI2 = 0.9,
                                    lower_CrI3 = 0.25,
                                    upper_CrI3 = 0.75) {
  
  # Extract usage samples
  if (usage) {
    
    # Extract usage parameters
    if (cmdstanr) {
      usage_draws <- usage_fit_raw$draws(format = "draws_df")
      Pbb_u <- usage_draws %>% select(starts_with("u_tilde[")) %>% divide_by(N_bb)
      P_u <- usage_draws %>% select(starts_with("P["))
      P0_u <- usage_draws %>% select(starts_with("P0["))
      D_u <- usage_draws %>% select(starts_with("D["))
      C_u <- usage_draws %>% select(starts_with("C["))
      C0_u <- usage_draws %>% select(starts_with("C0["))
      invlam_u <- usage_draws %>% select(starts_with("inv_lambda["))
    } else {
      extracted_usage <- extract(usage_fit)
      Pbb_u <- extracted_usage$u_tilde / N_bb
      P_u <- extracted_usage$P
      P0_u <- extracted_usage$P0
      D_u <- extracted_usage$D
      C_u <- extracted_usage$C
      C0_u <- extracted_usage$C0
      invlam_u <- extracted_usage$inv_lambda
    }
    
    # Conditional usage
    PC_u <- C_u / P_u
    
    # Time repetitions of inverse lambda
    invlam_u <- invlam_u[, rep(seq_len(ncol(invlam_u)), each = N_CMC)]
    
    # Mean retention
    lam_u <- 1 / invlam_u
    ret_u <- 1 / (lam_u * (1-D_u))

    # Calculate mean values
    Pbb_u_mean <- Pbb_u %>% apply(2, mean)
    P_u_mean <- P_u %>% apply(2, mean)
    P0_u_mean <- P0_u %>% apply(2, mean)
    D_u_mean <- D_u %>% apply(2, mean)
    C_u_mean <- C_u %>% apply(2, mean)
    PC_u_mean <- PC_u %>% apply(2, mean)
    invlam_u_mean <- invlam_u %>% apply(2, mean)
    ret_u_mean <- ret_u %>% apply(2, mean)

    # Sort usage parameters
    Pbb_u %<>% apply(2, sort)
    P_u %<>% apply(2, sort)
    P0_u %<>% apply(2, sort)
    D_u %<>% apply(2, sort)
    C_u %<>% apply(2, sort)
    PC_u %<>% apply(2, sort)
    invlam_u %<>% apply(2, sort)
    ret_u %<>% apply(2, sort)

    # Credible interval bounds
    N_u_samples <- dim(Pbb_u)[1]
    LB1_ID <- round(N_u_samples * lower_CrI1)
    UB1_ID <- round(N_u_samples * upper_CrI1)
    LB2_ID <- round(N_u_samples * lower_CrI2)
    UB2_ID <- round(N_u_samples * upper_CrI2)
    LB3_ID <- round(N_u_samples * lower_CrI3)
    UB3_ID <- round(N_u_samples * upper_CrI3)
    
    # Calculate lower bound of credible intervals
    Pbb_u_LB1 <- Pbb_u[LB1_ID,]
    P_u_LB1 <- P_u[LB1_ID,]
    P0_u_LB1 <- P0_u[LB1_ID,]
    D_u_LB1 <- D_u[LB1_ID,]
    C_u_LB1 <- C_u[LB1_ID,]
    PC_u_LB1 <- PC_u[LB1_ID,]
    invlam_u_LB1 <- invlam_u[LB1_ID,]
    ret_u_LB1 <- ret_u[LB1_ID,]
    Pbb_u_LB2 <- Pbb_u[LB2_ID,]
    Pbb_u_LB3 <- Pbb_u[LB3_ID,]
    
    # Calculate upper bound of credible intervals
    Pbb_u_UB1 <- Pbb_u[UB1_ID,]
    P_u_UB1 <- P_u[UB1_ID,]
    P0_u_UB1 <- P0_u[UB1_ID,]
    D_u_UB1 <- D_u[UB1_ID,]
    C_u_UB1 <- C_u[UB1_ID,]
    PC_u_UB1 <- PC_u[UB1_ID,]
    invlam_u_UB1 <- invlam_u[UB1_ID,]
    ret_u_UB1 <- ret_u[UB1_ID,]
    Pbb_u_UB2 <- Pbb_u[UB2_ID,]
    Pbb_u_UB3 <- Pbb_u[UB3_ID,]
    
    # Append to dataset
    dataset <- data.frame(dataset,
                          "Pbb_u_mean" = Pbb_u_mean,
                          "Pbb_u_LB1" = Pbb_u_LB1,
                          "Pbb_u_UB1" = Pbb_u_UB1,
                          "Pbb_u_LB2" = Pbb_u_LB2,
                          "Pbb_u_UB2" = Pbb_u_UB2,
                          "Pbb_u_LB3" = Pbb_u_LB3,
                          "Pbb_u_UB3" = Pbb_u_UB3,
                          "P_u_mean" = P_u_mean,
                          "P_u_LB1" = P_u_LB1,
                          "P_u_UB1" = P_u_UB1,
                          "P0_u_mean" = P0_u_mean,
                          "P0_u_LB1" = P0_u_LB1,
                          "P0_u_UB1" = P0_u_UB1,
                          "D_u_mean" = D_u_mean,
                          "D_u_LB1" = D_u_LB1,
                          "D_u_UB1" = D_u_UB1,
                          "C_u_mean" = C_u_mean,
                          "C_u_LB1" = C_u_LB1,
                          "C_u_UB1" = C_u_UB1,
                          "PC_u_mean" = PC_u_mean,
                          "PC_u_LB1" = PC_u_LB1,
                          "PC_u_UB1" = PC_u_UB1,
                          "invlam_u_mean" = invlam_u_mean,
                          "invlam_u_LB1" = invlam_u_LB1,
                          "invlam_u_UB1" = invlam_u_UB1,
                          "ret_u_mean" = ret_u_mean,
                          "ret_u_LB1" = ret_u_LB1,
                          "ret_u_UB1" = ret_u_UB1)
  }
    
  # Extract access samples
  if (access) {

    # Extract access parameters
    if (cmdstanr) {
      access_draws <- access_fit_raw$draws(format = "draws_df")
      Pbb_a <- access_draws %>% select(starts_with("u_tilde[")) %>% divide_by(N_bb)
      P_a <- access_draws %>% select(starts_with("P["))
      P0_a <- access_draws %>% select(starts_with("P0["))
      D_a <- access_draws %>% select(starts_with("D["))
      C_a <- access_draws %>% select(starts_with("C["))
      C0_a <- access_draws %>% select(starts_with("C0["))
      invlam_a <- access_draws %>% select(starts_with("inv_lambda["))
    } else {
      extracted_access <- extract(access_fit)
      Pbb_a <- extracted_access$u_tilde / N_bb
      P_a <- extracted_access$P
      P0_a <- extracted_access$P0
      D_a <- extracted_access$D
      C_a <- extracted_access$C
      C0_a <- extracted_access$C0
      invlam_a <- extracted_access$inv_lambda
    }
    
    # Proportion of accessible nets from campaigns
    PC_a <- C_a / P_a
    
    # Time repetitions of inverse lambda
    invlam_a <- invlam_a[, rep(seq_len(ncol(invlam_a)), each = N_CMC)]
    
    # Mean retention
    lam_a <- 1 / invlam_a
    ret_a <- 1 / (lam_a * (1-D_a))
    
    # Calculate mean values
    Pbb_a_mean <- Pbb_a %>% apply(2, mean)
    P_a_mean <- P_a %>% apply(2, mean)
    P0_a_mean <- P0_a %>% apply(2, mean)
    D_a_mean <- D_a %>% apply(2, mean)
    C_a_mean <- C_a %>% apply(2, mean)
    PC_a_mean <- PC_a %>% apply(2, mean)
    invlam_a_mean <- invlam_a %>% apply(2, mean)
    ret_a_mean <- ret_a %>% apply(2, mean)
    
    # Sort usage parameters
    Pbb_a %<>% apply(2, sort)
    P_a %<>% apply(2, sort)
    P0_a %<>% apply(2, sort)
    D_a %<>% apply(2, sort)
    C_a %<>% apply(2, sort)
    PC_a %<>% apply(2, sort)
    invlam_a %<>% apply(2, sort)
    ret_a %<>% apply(2, sort)
    
    # Credible interval bounds
    N_a_samples <- dim(Pbb_a)[1]
    LB1_ID <- round(N_a_samples * lower_CrI1)
    UB1_ID <- round(N_a_samples * upper_CrI1)
    LB2_ID <- round(N_a_samples * lower_CrI2)
    UB2_ID <- round(N_a_samples * upper_CrI2)
    LB3_ID <- round(N_a_samples * lower_CrI3)
    UB3_ID <- round(N_a_samples * upper_CrI3)
    
    # Calculate lower bound of credible intervals
    Pbb_a_LB1 <- Pbb_a[LB1_ID,]
    P_a_LB1 <- P_a[LB1_ID,]
    P0_a_LB1 <- P0_a[LB1_ID,]
    D_a_LB1 <- D_a[LB1_ID,]
    C_a_LB1 <- C_a[LB1_ID,]
    PC_a_LB1 <- PC_a[LB1_ID,]
    invlam_a_LB1 <- invlam_a[LB1_ID,]
    ret_a_LB1 <- ret_a[LB1_ID,]
    Pbb_a_LB2 <- Pbb_a[LB2_ID,]
    Pbb_a_LB3 <- Pbb_a[LB3_ID,]
    
    # Calculate upper bound of credible intervals
    Pbb_a_UB1 <- Pbb_a[UB1_ID,]
    P_a_UB1 <- P_a[UB1_ID,]
    P0_a_UB1 <- P0_a[UB1_ID,]
    D_a_UB1 <- D_a[UB1_ID,]
    C_a_UB1 <- C_a[UB1_ID,]
    PC_a_UB1 <- PC_a[UB1_ID,]
    invlam_a_UB1 <- invlam_a[UB1_ID,]
    ret_a_UB1 <- ret_a[UB1_ID,]
    Pbb_a_UB2 <- Pbb_a[UB2_ID,]
    Pbb_a_UB3 <- Pbb_a[UB3_ID,]
    
    # Append to dataset
    dataset <- data.frame(dataset,
                          "Pbb_a_mean" = Pbb_a_mean,
                          "Pbb_a_LB1" = Pbb_a_LB1,
                          "Pbb_a_UB1" = Pbb_a_UB1,
                          "Pbb_a_LB2" = Pbb_a_LB2,
                          "Pbb_a_UB2" = Pbb_a_UB2,
                          "Pbb_a_LB3" = Pbb_a_LB3,
                          "Pbb_a_UB3" = Pbb_a_UB3,
                          "P_a_mean" = P_a_mean,
                          "P_a_LB1" = P_a_LB1,
                          "P_a_UB1" = P_a_UB1,
                          "P0_a_mean" = P0_a_mean,
                          "P0_a_LB1" = P0_a_LB1,
                          "P0_a_UB1" = P0_a_UB1,
                          "D_a_mean" = D_a_mean,
                          "D_a_LB1" = D_a_LB1,
                          "D_a_UB1" = D_a_UB1,
                          "C_a_mean" = C_a_mean,
                          "C_a_LB1" = C_a_LB1,
                          "C_a_UB1" = C_a_UB1,
                          "PC_a_mean" = PC_a_mean,
                          "PC_a_LB1" = PC_a_LB1,
                          "PC_a_UB1" = PC_a_UB1,
                          "invlam_a_mean" = invlam_a_mean,
                          "invlam_a_LB1" = invlam_a_LB1,
                          "invlam_a_UB1" = invlam_a_UB1,
                          "ret_a_mean" = ret_a_mean,
                          "ret_a_LB1" = ret_a_LB1,
                          "ret_a_UB1" = ret_a_UB1)
  }
  
  # Extract conditional usage
  if (usage & access) {
    if (N_u_samples == N_a_samples) {
      
      P_condu <- P_u / P_a
      P_condu_mean <- P_condu %>% apply(2, mean)
      P_condu_LB1 <- P_condu[LB1_ID,]
      P_condu_UB1 <- P_condu[UB1_ID,]
      
      # Append to dataset
      dataset <- data.frame(dataset,
                            "P_condu_mean" = P_condu_mean,
                            "P_condu_LB1" = P_condu_LB1,
                            "P_condu_UB1" = P_condu_UB1)
    } else {
      print(paste0("Warning: Conditional usage not calculated due to different",
                   "sample sizes for usage and access. There are ", N_u_samples,
                   " usage samples but ", N_a_samples, " access samples."))
    }
  }
  
  # Return dataset
  return(dataset)
  
}
