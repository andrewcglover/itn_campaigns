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
                    N_disc = 20,
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
    usage_fit <<- stan('./scripts/stan/use_acc_reg_all_ccc.stan',
                       data = usage_list,
                       iter = 300,
                       warmup = 200,
                       chains = 4,
                       init_r = 1e-2
                       # control = list(adapt_delta = 0.99,
                       #                stepsize = 0.5,
                       #                max_treedepth = 15
                       #                )
                       )
  } else {
    access_fit <<- stan('./scripts/stan/use_acc_reg_all_ccc.stan',
                        data = access_list,
                        iter = 300,
                        warmup = 200,
                        chains = 4,
                        init_r = 1e-2
                        # control = list(adapt_delta = 0.99,
                        #                stepsize = 0.5,
                        #                max_treedepth = 15
                        #                )
                        )
  }
}

#-------------------------------------------------------------------------------
# Extract samples

append_time_series_fits <- function(dataset,
                                    usage = TRUE,
                                    access = TRUE,
                                    lower_CrI = 0.025,
                                    upper_CrI = 0.975) {
  
  # Extract usage samples
  if (usage) {
    extracted_usage <- extract(usage_fit)
    
    # Extract usage parameters
    Pbb_u <- extracted_usage$u_tilde / N_bb
    P_u <- extracted_usage$P
    P0_u <- extracted_usage$P0
    D_u <- extracted_usage$D
    C_u <- extracted_usage$C
    PC_u <- C_u / P_u
    invlam_u <- extracted_usage$inv_lambda
    
    # Calculate mean values
    Pbb_u_mean <- Pbb_u %>% apply(2, mean)
    P_u_mean <- P_u %>% apply(2, mean)
    P0_u_mean <- P0_u %>% apply(2, mean)
    D_u_mean <- D_u %>% apply(2, mean)
    C_u_mean <- C_u %>% apply(2, mean)
    PC_u_mean <- PC_u %>% apply(2, mean)
    invlam_u_mean <- invlam_u %>% apply(2, mean)
    
    # Sort usage parameters
    Pbb_u %<>% apply(2, sort)
    P_u %<>% apply(2, sort)
    P0_u %<>% apply(2, sort)
    D_u %<>% apply(2, sort)
    C_u %<>% apply(2, sort)
    PC_u %<>% apply(2, sort)
    invlam_u %<>% apply(2, sort)
    
    # Credible interval bounds
    N_u_samples <- dim(Pbb_u)[1]
    LB_ID <- round(N_u_samples * lower_CrI)
    UB_ID <- round(N_u_samples * upper_CrI)
    
    # Calculate lower bound of credible intervals
    Pbb_u_LB <- Pbb_u[LB_ID,]
    P_u_LB <- P_u[LB_ID,]
    P0_u_LB <- P0_u[LB_ID,]
    D_u_LB <- D_u[LB_ID,]
    C_u_LB <- C_u[LB_ID,]
    PC_u_LB <- PC_u[LB_ID,]
    invlam_u_LB <- invlam_u[LB_ID,]
    
    # Calculate upper bound of credible intervals
    Pbb_u_UB <- Pbb_u[UB_ID,]
    P_u_UB <- P_u[UB_ID,]
    P0_u_UB <- P0_u[UB_ID,]
    D_u_UB <- D_u[UB_ID,]
    C_u_UB <- C_u[UB_ID,]
    PC_u_UB <- PC_u[UB_ID,]
    invlam_u_UB <- invlam_u[UB_ID,]
    
    # Append to dataset
    dataset <- data.frame(dataset,
                          "Pbb_u_mean" = Pbb_u_mean,
                          "Pbb_u_LB" = Pbb_u_LB,
                          "Pbb_u_UB" = Pbb_u_UB,
                          "P_u_mean" = P_u_mean,
                          "P_u_LB" = P_u_LB,
                          "P_u_UB" = P_u_UB,
                          "P0_u_mean" = P0_u_mean,
                          "P0_u_LB" = P0_u_LB,
                          "P0_u_UB" = P0_u_UB,
                          "D_u_mean" = D_u_mean,
                          "D_u_LB" = D_u_LB,
                          "D_u_UB" = D_u_UB,
                          "C_u_mean" = C_u_mean,
                          "C_u_LB" = C_u_LB,
                          "C_u_UB" = C_u_UB,
                          "PC_u_mean" = PC_u_mean,
                          "PC_u_LB" = PC_u_LB,
                          "PC_u_UB" = PC_u_UB,
                          "invlam_u_mean" = invlam_u_mean,
                          "invlam_u_LB" = invlam_u_LB,
                          "invlam_u_UB" = invlam_u_UB)
  }
    
  # Extract access samples
  if (access) {
    extracted_access <- extract(access_fit)
    
    # Extract access parameters
    Pbb_a <- extracted_access$u_tilde / N_bb
    P_a <- extracted_access$P
    P0_a <- extracted_access$P0
    D_a <- extracted_access$D
    C_a <- extracted_access$C
    PC_a <- C_a / P_a
    invlam_a <- extracted_access$inv_lambda
    
    # Calculate mean values
    Pbb_a_mean <- Pbb_a %>% apply(2, mean)
    P_a_mean <- P_a %>% apply(2, mean)
    P0_a_mean <- P0_a %>% apply(2, mean)
    D_a_mean <- D_a %>% apply(2, mean)
    C_a_mean <- C_a %>% apply(2, mean)
    PC_a_mean <- PC_a %>% apply(2, mean)
    invlam_a_mean <- invlam_a %>% apply(2, mean)
    
    # Sort usage parameters
    Pbb_a %<>% apply(2, sort)
    P_a %<>% apply(2, sort)
    P0_a %<>% apply(2, sort)
    D_a %<>% apply(2, sort)
    C_a %<>% apply(2, sort)
    PC_a %<>% apply(2, sort)
    invlam_a %<>% apply(2, sort)
    
    # Credible interval bounds
    N_a_samples <- dim(Pbb_a)[1]
    LB_ID <- round(N_a_samples * lower_CrI)
    UB_ID <- round(N_a_samples * upper_CrI)
    
    # Calculate lower bound of credible intervals
    Pbb_a_LB <- Pbb_a[LB_ID,]
    P_a_LB <- P_a[LB_ID,]
    P0_a_LB <- P0_a[LB_ID,]
    D_a_LB <- D_a[LB_ID,]
    C_a_LB <- C_a[LB_ID,]
    PC_a_LB <- PC_a[LB_ID,]
    invlam_a_LB <- invlam_a[LB_ID,]
    
    # Calculate upper bound of credible intervals
    Pbb_a_UB <- Pbb_a[UB_ID,]
    P_a_UB <- P_a[UB_ID,]
    P0_a_UB <- P0_a[UB_ID,]
    D_a_UB <- D_a[UB_ID,]
    C_a_UB <- C_a[UB_ID,]
    PC_a_UB <- PC_a[UB_ID,]
    invlam_a_UB <- invlam_a[UB_ID,]
    
    # Append to dataset
    dataset <- data.frame(dataset,
                          "Pbb_a_mean" = Pbb_a_mean,
                          "Pbb_a_LB" = Pbb_a_LB,
                          "Pbb_a_UB" = Pbb_a_UB,
                          "P_a_mean" = P_a_mean,
                          "P_a_LB" = P_a_LB,
                          "P_a_UB" = P_a_UB,
                          "P0_a_mean" = P0_a_mean,
                          "P0_a_LB" = P0_a_LB,
                          "P0_a_UB" = P0_a_UB,
                          "D_a_mean" = D_a_mean,
                          "D_a_LB" = D_a_LB,
                          "D_a_UB" = D_a_UB,
                          "C_a_mean" = C_a_mean,
                          "C_a_LB" = C_a_LB,
                          "C_a_UB" = C_a_UB,
                          "PC_a_mean" = PC_a_mean,
                          "PC_a_LB" = PC_a_LB,
                          "PC_a_UB" = PC_a_UB,
                          "invlam_a_mean" = invlam_a_mean,
                          "invlam_a_LB" = invlam_a_LB,
                          "invlam_a_UB" = invlam_a_UB)
  }
  
  # Extract conditional usage
  if (N_u_samples == N_a_samples) {
    if (usage & access) {
      
      P_condu <- P_u / P_a
      P_condu_mean <- P_condu %>% apply(2, mean)
      P_condu_LB <- P_condu[LB_ID,]
      P_condu_UB <- P_condu[UB_ID,]
      
      # Append to dataset
      dataset <- data.frame(dataset,
                            "P_condu_mean" = P_condu_mean,
                            "P_condu_LB" = P_condu_LB,
                            "P_condu_UB" = P_condu_UB)
    }
  } else {
    print("Warning: different sample sizes for usage and access")
  }
  
  # Return dataset
  return(dataset)
  
}
