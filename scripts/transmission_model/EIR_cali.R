# EIR_cali.R

par_net_region_sequential_v4 <- function(param_list) {
  
  #library(magrittr)
  
  #N_reps <- length(param_list)
  
  
  # Extract parameters from parameter list
  site_pars <- param_list#[[1]] #Uncomment for testing
  mean_retu <- site_pars$mean_retu
  mean_reta <- site_pars$mean_reta
  npc_beta <- site_pars$npc_beta
  npc_gamma <- site_pars$npc_gamma
  net_type <- site_pars$net_type
  net_name <- site_pars$net_name
  net_strategy <- site_pars$net_strategy
  last_camp <- site_pars$last_camp
  mass_int_mn <- site_pars$mass_int_mn
  ISO2 <- site_pars$ISO2
  fs_area <- site_pars$fs_area
  ISO2 <- site_pars$ISO2
  fs_name_1 <- site_pars$fs_name_1
  urbanicity <- site_pars$urbanicity
  fs_area_id <- site_pars$fs_area_id
  N_species <- site_pars$N_species
  CMC_first <- site_pars$CMC_first
  CMC_Jan2000 <- site_pars$CMC_Jan2000
  N_CMC <- site_pars$N_CMC
  N_CMC_sim <- site_pars$N_CMC_sim
  sim_population <- site_pars$sim_population
  P_a_mid <- site_pars$P_a_mid
  P0_a_mid <- site_pars$P0_a_mid
  D_a_mid <- site_pars$D_a_mid
  lam_a_mid <- site_pars$lam_a_mid
  P_u_mid <- site_pars$P_u_mid
  P0_u_mid <- site_pars$P0_u_mid
  D_u_mid <- site_pars$D_u_mid
  lam_u_mid <- site_pars$lam_u_mid
  dn0_mat <- site_pars$dn0_mat
  rn_mat <- site_pars$rn_mat
  rnm_mat <- site_pars$rnm_mat
  gam_vec <- site_pars$gam_vec
  DOY_1st <- site_pars$DOY_1st
  DOY_mid <- site_pars$DOY_mid
  prevalence_rendering_min_ages <- site_pars$prevalence_rendering_min_ages
  prevalence_rendering_max_ages <- site_pars$prevalence_rendering_max_ages

  # if (biennial_reduction & (mass_int_mn < 25)) {
  #   net_strategy <- paste0(net_strategy, "_bien_costed")
  # } else if ((cost_factor < 0.9999) | (cost_factor > 1.0001)) {
  #   net_strategy <- paste0(net_strategy, "_costed")
  # }

  
  
  fs_area_undrscr <- gsub(" ", "_", fs_area)
  
  year = 365

  # Monthly rates
  decay_a <- 1 / mean_reta
  decay_u <- 1 / mean_retu
  decay_uga <- decay_u - decay_a
  decay_npc <- decay_a / npc_gamma
  
  # convert retention to days
  mean_retu_dy <- 365 * mean_retu / 12
  mean_reta_dy <- 365 * mean_reta / 12
  
  # Central time point for first regular mass campaign

  # Back-fill preliminary period
  N_early <- CMC_first - CMC_Jan2000
  
  # Times for fitting
  times_mn <- seq(1, N_CMC)
  times_yr <- rep(seq(0, ceiling(N_CMC_sim / 12)), each=12)
  times_1st_dy <- DOY_1st + (times_yr * year)
  times_mid_dy <- DOY_mid + (times_yr * year)
  input_net_times <- times_mid_dy[1:N_CMC_sim]    # usage for fitting
  output_net_times <- times_1st_dy[1:N_CMC_sim]   # distribution times for netz
  
  #-----------------------------------------------------------------------------
  # usage profiles
  
  # Extract values for selected sample
  P_ui <- P_u_mid
  P0_ui <- P0_u_mid
  D_ui <- D_u_mid
  lambda_ui <- lam_u_mid
  
  P_u_full <- c(rep(D_ui[1], N_early), P_ui)
  
  # Fit nets with no future mass campaigns
  cali_dist_list <- fit_usage_sequential_distributions(
    target_usage = P_u_full,
    target_usage_timesteps = input_net_times,
    distribution_timesteps = output_net_times,
    mean_retention = mean_retu_dy
  )
  
  cali_dist <- cali_dist_list[[1]]
  
  dn0_mat <-  site_pars$bednet_dn0[
    rep(seq_len(nrow(site_pars$bednet_dn0)), each = 12),
  ]
  rn_mat <-  site_pars$bednet_rn[
    rep(seq_len(nrow(site_pars$bednet_rn)), each = 12),
  ]
  rnm_mat <-  site_pars$bednet_rnm[
    rep(seq_len(nrow(site_pars$bednet_rnm)), each = 12),
  ]
  gamman_vec <- rep(site_pars$bednet_gamman, each = 12)
  
  # set bednets
  bednet_pars <- malariasimulation::set_bednets(site_pars,
                                                timesteps = output_net_times,
                                                coverages = cali_dist,
                                                retention = mean_retu_dy,
                                                dn0 = dn0_mat,
                                                rn = rn_mat,
                                                rnm = rnm_mat,
                                                gamman = gamman_vec)
  
  target <- adm_site$prevalence$pfpr
  weights = rep(1, length(target))
  weighted_target <- target * weights
  
  weighted_annual_pfpr_summary <- function(x, w = weights){
    year <- ceiling(x$timestep / 365)
    pfpr <- x$n_detect_lm_182_1825  / x$n_age_182_1825
    tapply(pfpr, year, mean) * w
  }
  
  
  EIR_out <- calibrate(
    parameters = bednet_pars,
    target = weighted_target,
    summary_function = weighted_annual_pfpr_summary,
    eq_prevalence = target[1]
  )
  
  
  
  bednet_pars$human_population <- 5000
  parameters <- set_equilibrium(bednet_pars, init_EIR = EIR_out)
  raw <- run_simulation(bednet_pars$timesteps + 100, parameters = bednet_pars)
  raw$pfpr <- raw$n_detect_lm_182_1825  / raw$n_age_182_1825
  
  
  ggplot() +
    geom_point(aes(x = 365 * (0:24 + 0.5), y = target), col = "dodgerblue", size = 4) +
    geom_line(data = raw, aes(x = timestep, y = pfpr), col = "deeppink", linewidth = 1) +
    ylim(0, 1) +
    ylab(expression(italic(Pf)*Pr[2-10])) +
    xlab("Time") +
    theme_bw()
  
  
  
  
  
  
  target <- dataset$rdt_prev[area_time_ids]
  target_alpha <- dataset$rdt_pos[area_time_ids] + 0.5
  target_beta <- dataset$rdt_neg[area_time_ids] + 0.5
  target_var <- (target_alpha * target_beta) /
    ((target_alpha + target_beta)^2 * (target_alpha + target_beta + 1))
  
  inc_months <- which(!is.nan(target))
  target_lo <- dataset$lo_prev[area_time_ids]
  target_hi <- dataset$hi_prev[area_time_ids]
  target_lo <- target_lo[inc_months]
  target_hi <- target_hi[inc_months]
  target <- target[inc_months]
  target_var <- target_var[inc_months]
  weights <- 1 / target_var
  weights <- weights / sum(weights)
  weighted_target <- target * weights
  #weights = rep(1, length(target))
  #weighted_target <- target * weights
  
  # annual_filtered_pfpr_summary <- function(x){
  #   month <- ceiling(x$timestep / (365/12))
  #   pfpr <- x$n_detect_lm_182_1825  / x$n_age_182_1825
  #   filtered_months <- which(month %in% inc_months)
  #   pfpr <- pfpr[filtered_months]
  #   month <- month[filtered_months]
  #   tapply(pfpr, year, mean)
  # }
  # 
  
  weighted_filtered_annual_pfpr_summary <- function(x, w = weights){
    month <- ceiling(x$timestep / (365/12))
    pfpr_6_59_mo <- (x$n_detect_lm_182_730 + x$n_detect_lm_730_1825) /
      (x$n_age_182_730 + x$n_age_730_1825)
    filtered_months <- which(month %in% inc_months)
    pfpr_6_59_mo <- pfpr_6_59_mo[filtered_months]
    month <- month[filtered_months]
    tapply(pfpr_6_59_mo, month, mean) * w
  }
  
  filtered_annual_pfpr_summary <- function(x){
    month <- ceiling(x$timestep / (365/12))
    pfpr_6_59_mo <- (x$n_detect_lm_182_730 + x$n_detect_lm_730_1825) /
      (x$n_age_182_730 + x$n_age_730_1825)
    filtered_months <- which(month %in% inc_months)
    pfpr_6_59_mo <- pfpr_6_59_mo[filtered_months]
    month <- month[filtered_months]
    tapply(pfpr_6_59_mo, month, mean)
  }
  
  set.seed(123)
  EIR_out <- cali::calibrate(
    parameters = bednet_pars,
    target = weighted_target,
    summary_function = weighted_filtered_annual_pfpr_summary,
    eq_prevalence = target[1],
    use_pfpr_6_59_mo = TRUE
  )
  pfpr1 <-cali:::pfpr_6_59_mo(eir = EIR1, ft = 0)
  EIR2 <- cali::calibrate(
    parameters = bednet_pars,
    target = weighted_target,
    summary_function = weighted_filtered_annual_pfpr_summary,
    eq_prevalence = pfpr1,
    use_pfpr_6_59_mo = TRUE
  )
  
  
  
  bednet_pars$human_population <- 5000
  bednet_pars <- malariasimulation::set_equilibrium(bednet_pars, init_EIR = EIR_out)
  raw <- malariasimulation::run_simulation(bednet_pars$timesteps, parameters = bednet_pars)
  raw$pfpr <- (raw$n_detect_lm_182_730 + raw$n_detect_lm_730_1825) /
    (raw$n_age_182_730 + raw$n_age_730_1825)

  
  
  ggplot() +
    geom_line(
      data = raw,
      aes(
        x = timestep,
        y = pfpr
      ),
      col = "deeppink",
      linewidth = 0.5
    ) +
    geom_point(
      aes(
        x = inc_months * 365 / 12,
        y = target
      ),
      col = "dodgerblue",
      size = 2
    ) +
    geom_errorbar(
      aes(
        x =  inc_months * 365 / 12,
        y = target,
        ymin = target_lo,
        ymax = target_hi
      ),
      col = "dodgerblue"
    ) +
    ylim(0, 1) +
    ylab(expression(P*italic(f)*Pr["6-59mo"])) +
    xlab("Time") +
    theme_bw()
  
  
  
  pfpr_alp <- dataset$rdt_pos + 0.5
  pfpr_bet <- dataset$rdt_neg + 0.5
  pfpr_alp <- pfpr_alp[inc_months]
  pfpr_bet <- pfpr_bet[inc_months]
  
  
  
  
  ggplot() +
    geom_point(aes(x = 365 * (0:24 + 0.5), y = target), col = "dodgerblue", size = 4) +
    geom_line(data = raw, aes(x = timestep, y = pfpr), col = "deeppink", linewidth = 1) +
    ylim(0, 1) +
    ylab(expression(italic(Pf)*Pr[2-10])) +
    xlab("Time") +
    theme_bw()
  
  
  
  
  
  
  # run simulation
  output <- malariasimulation::run_simulation(timesteps = bednet_pars$timesteps,
                                              parameters = bednet_pars)
  
  
  
  
  # Extend values
  P0_ui_end <- tail(P0_ui, n = 1)
  D_ui_end <- tail(D_ui, n = 1)
  P0_ui_long <- c(P0_ui, rep(P0_ui_end, projection_window_mn))
  D_ui_long <- c(D_ui, rep(D_ui_end, projection_window_mn))
  P0_ui_tail <- P0_ui_long[last_camp:proj_end]
  D_ui_tail <- D_ui_long[last_camp:proj_end]
  C0_ui_tail <- P0_ui_tail - D_ui_tail
  decay_ui_tail <- exp(-lambda_ui * m_tail)
  C_ui_tail <- C0_ui_tail * decay_ui_tail
  P_ui_tail <- C_ui_tail + D_ui_tail
  P_ui_long <- c(P_ui[1:(last_camp-1)], P_ui_tail)
  C_ui_long <- P_ui_long - D_ui_long
  
  # P0 and D over projection window
  P0_ui_proj <- tail(P0_ui_long, n = N_proj)
  D_ui_proj <- tail(D_ui_long, n = N_proj)
  C0_ui_proj <- P0_ui_proj - D_ui_proj
  
  # Back-fill first value for early usage
  D_ui_early <- rep(D_ui[1], N_early)
  P_ui_early <- rep(P_ui[1], N_early)
  C_ui_early <- P_ui_early - D_ui_early
  
  # Calculate routine usage
  D_ui_full <- c(D_ui_early, D_ui_long[1:(proj_camp_1-1)], D_ui_proj)
  
  # Calculate campaign usage
  decay_ui_proj <- exp(-lambda_ui * m_proj)
  C_ui_proj <- C0_ui_proj * decay_ui_proj
  C_ui_pre_proj <- c(C_ui_early, C_ui_long[1:(proj_camp_1-1)])
  C_ui_full <- c(C_ui_early, C_ui_long[1:(proj_camp_1-1)], C_ui_proj)
  
  # Calculate overall usage
  P_ui_full <- C_ui_full + D_ui_full
  
  # Usage with no future mass campaigns
  D_ui_proj_only <- c(P_ui_early,P_ui_long)
  
  # Fit nets with no future mass campaigns
  dist_used_no_future_mdc <- fit_usage_sequential_distributions(
    target_usage = D_ui_proj_only,
    target_usage_timesteps = input_net_times,
    distribution_timesteps = output_net_times,
    mean_retention = mean_retu_dy
  )
  
  # Fit nets with future mass campaigns
  dist_used_future_mdc <- fit_usage_sequential_distributions(
    target_usage = P_ui_full,
    target_usage_timesteps = input_net_times,
    distribution_timesteps = output_net_times,
    mean_retention = mean_retu_dy
  )
  
  # Future nets distributed through mass campaigns
  dist_used_future_mdc_only <- list(dist_used_future_mdc[[1]] - dist_used_no_future_mdc[[1]],
                                    dist_used_future_mdc[[2]] - dist_used_no_future_mdc[[2]],
                                    dist_used_future_mdc[[3]] - dist_used_no_future_mdc[[3]])
  
  #-----------------------------------------------------------------------------
  # Access profiles
  
  # Extract values for selected sample
  P_ai <- P_a_mid[sid,]
  P0_ai <- P0_a_mid[sid,]
  D_ai <- D_a_mid[sid,]
  lambda_ai <- lam_a_mid[sid]
  
  # Extend values
  P0_ai_end <- tail(P0_ai, n = 1)
  D_ai_end <- tail(D_ai, n = 1)
  P0_ai_long <- c(P0_ai, rep(P0_ai_end, projection_window_mn))
  D_ai_long <- c(D_ai, rep(D_ai_end, projection_window_mn))
  P0_ai_tail <- P0_ai_long[last_camp:proj_end]
  D_ai_tail <- D_ai_long[last_camp:proj_end]
  C0_ai_tail <- P0_ai_tail - D_ai_tail
  decay_ai_tail <- exp(-lambda_ai * m_tail)
  C_ai_tail <- C0_ai_tail * decay_ai_tail
  P_ai_tail <- C_ai_tail + D_ai_tail
  P_ai_long <- c(P_ai[1:(last_camp-1)], P_ai_tail)
  C_ai_long <- P_ai_long - D_ai_long
  
  # P0 and D over projection window
  P0_ai_proj <- tail(P0_ai_long, n = N_proj)
  D_ai_proj <- tail(D_ai_long, n = N_proj)
  C0_ai_proj <- P0_ai_proj - D_ai_proj
  
  # Back-fill first value for early usage
  D_ai_early <- rep(D_ai[1], N_early)
  P_ai_early <- rep(P_ai[1], N_early)
  C_ai_early <- P_ai_early - D_ai_early
  
  # Calculate routine usage
  D_ai_full <- c(D_ai_early, D_ai_long[1:(proj_camp_1-1)], D_ai_proj)
  
  # Calculate campaign usage
  decay_ai_proj <- exp(-lambda_ai * m_proj)
  C_ai_proj <- C0_ai_proj * decay_ai_proj
  C_ai_pre_proj <- c(C_ai_early, C_ai_long[1:(proj_camp_1-1)])
  C_ai_full <- c(C_ai_early, C_ai_long[1:(proj_camp_1-1)], C_ai_proj)
  
  # Calculate overall usage
  P_ai_full <- C_ai_full + D_ai_full
  
  # Usage with no future mass campaigns
  D_ai_proj_only <- c(P_ai_early,P_ai_long)
  
  # Fit nets with no future mass campaigns
  dist_access_no_future_mdc <- fit_usage_sequential_distributions(
    target_usage = D_ai_proj_only,
    target_usage_timesteps = input_net_times,
    distribution_timesteps = output_net_times,
    mean_retention = mean_reta_dy
  )
  
  # Fit nets with future mass campaigns
  dist_access_future_mdc <- fit_usage_sequential_distributions(
    target_usage = P_ai_full,
    target_usage_timesteps = input_net_times,
    distribution_timesteps = output_net_times,
    mean_retention = mean_reta_dy
  )
  
  # Future nets distributed through mass campaigns
  dist_access_future_mdc_only <- list(dist_access_future_mdc[[1]] - dist_access_no_future_mdc[[1]],
                                      dist_access_future_mdc[[2]] - dist_access_no_future_mdc[[2]],
                                      dist_access_future_mdc[[3]] - dist_access_no_future_mdc[[3]])
  
  #-----------------------------------------------------------------------------
  # nets per capita calculations
  
  # Default usage distribution
  u_dist <- dist_used_future_mdc[[1]]
  
  # Routine nets per capita
  D_a0 <- tail(dist_access_no_future_mdc[[2]], n = 1)
  D_a1 <- tail(dist_access_no_future_mdc[[3]], n = 1)
  D_r0 <- (D_a0 / npc_beta) ^ (1 / npc_gamma)
  D_r1 <- (D_a1 / npc_beta) ^ (1 / npc_gamma)
  D_npc_dist <- (D_r0 - D_r1) / (1 - D_r1)
  
  # Uncosted nets per capita
  P_a0 <- dist_access_future_mdc[[2]]
  P_a1 <- dist_access_future_mdc[[3]]
  P_r0 <- (P_a0 / npc_beta) ^ (1 / npc_gamma)
  P_r1 <- (P_a1 / npc_beta) ^ (1 / npc_gamma)
  P_npc_dist <- (P_r0 - P_r1) / (1 - P_r1)
  
  # Total uncosted new nets distributed
  total_new_npc <- sum(P_npc_dist[new_start:new_end])
  total_new_npc_D <- D_npc_dist * projection_window_mn
  total_new_npc_C <- total_new_npc - total_new_npc_D
  ann_npc_D <- 12 * total_new_npc_D / projection_window_mn
  
  # Costed nets per capita
  if (net_costings | biennial_reduction) {
    
    # Usage and use given access
    P_u0 <- dist_used_future_mdc[[2]]
    P_u1 <- dist_used_future_mdc[[3]]
    
    # Approximate use given access relationship
    uu <- tail(P_u0, n = projection_window_mn)
    aa <- tail(P_a0, n = projection_window_mn)
    uga_approx <- lm(uu ~ aa)
    
    # Total costed new nets distributed
    total_costed_new_npc <- total_new_npc * cost_factor
    total_costed_new_npc_C <- total_costed_new_npc - total_new_npc_D
    costed_camp_reduction <- total_costed_new_npc_C / total_new_npc_C
    if (biennial_reduction) {
      costed_camp_reduction <- costed_camp_reduction * 2 / 3
    }
    
    # Costed new nets per capita
    costed_u_dist <- dist_used_future_mdc[[1]]
    costed_P_u0 <- P_u0
    costed_P_u1 <- P_u1
    costed_P_a0 <- P_a0
    costed_P_a1 <- P_a1
    costed_P_r0 <- P_r0
    costed_P_r1 <- P_r1
    costed_P_uga <- P_u0 / P_a0
    costed_P_npc_dist <- P_npc_dist
    last_npc <- P_r1[new_start]
    for (t in new_start:new_end) {
      rho_D <- D_npc_dist
      rho_C <- P_npc_dist[t] - D_npc_dist
      costed_P_npc_dist[t] <- rho_D + rho_C * costed_camp_reduction
      costed_P_r0[t] <- costed_P_npc_dist[t] * (1 - costed_P_r1[t]) + costed_P_r1[t]
      costed_P_a0[t] <- npc_beta * costed_P_r0[t] ^ npc_gamma
      costed_P_u0[t] <- predict(uga_approx, data.frame(aa = costed_P_a0[t]))
      
      # Costed usage distributed
      costed_u_dist[t] <- (costed_P_u0[t] - costed_P_u1[t]) / (1 - costed_P_u1[t]) 
      
      # Update last values
      if (t < new_end) {
        costed_P_r1[t+1] <- costed_P_r0[t] * exp(-decay_npc)
        costed_P_a1[t+1] <- costed_P_a0[t] * exp(-decay_a)
        costed_P_u1[t+1] <- costed_P_u0[t] * exp(-decay_u)
      }
    }
    
    costed_P_npc_dist <- (costed_P_r0 - costed_P_r1) / (1 - costed_P_r1)
    
    # Update costed usage distributed
    u_dist <- costed_u_dist
    
  }
  
  #-----------------------------------------------------------------------------
  
  # Select net strategy
  if (no_future_nets) {
    all_output_nets <- u_dist
    all_output_nets[new_start:new_end] <- 0
    net_strategy <- "no future nets"
  } else if (routine_baseline) {
    #all_output_nets <- dist_used_no_future_mdc[[1]]
    D_topup <- mean(tail(dist_used_no_future_mdc[[1]], 12))
    all_output_nets <- u_dist
    all_output_nets[new_start:new_end] <- D_topup
    
    net_strategy <- paste(net_name, "routine only")
  } else {
    all_output_nets <- u_dist
    if (net_costings) {
      if (biennial_reduction) {
        net_strategy <- paste0(net_name, " ", mass_int_yr,
                               "-year campaigns net and freq costed")
      } else {
        net_strategy <- paste0(net_name, " ", mass_int_yr,
                               "-year campaigns net costed")
      }
    } else {
      if (biennial_reduction) {
        net_strategy <- paste0(net_name, " ", mass_int_yr,
                               "-year campaigns freq costed")
      } else {
        net_strategy <- paste0(net_name, " ", mass_int_yr,
                               "-year campaigns uncosted")
      }
    }
  }
  all_output_nets[all_output_nets < 0] <- 0
  all_output_nets[all_output_nets > 1] <- 1
  
  # Net costings
  if (routine_baseline) {
    avg_tail_nets <- D_npc_dist * 12
  } else if (net_costings | biennial_reduction) {
    avg_tail_nets <- sum(
      tail(
        costed_P_npc_dist * 12 / projection_window_mn,
        n = projection_window_mn
      )
    )
  } else {
    avg_tail_nets <- sum(
      tail(
        P_npc_dist * 12 / projection_window_mn,
        n = projection_window_mn
      )
    )
  }
  ann_cost_percap <- avg_tail_nets * new_net_cost
  
  # cost override
  if (override_cost) {
    if (override_mdc_only) {
      avg_tail_nets_routine <- sum(tail(dist_used_no_future_mdc * tail_pop, n = 6 * 12)) / 6
      avg_pop_adj_tail_nets_routine <- avg_tail_nets_routine / 1.8
      tail_net_cost_routine <- avg_pop_adj_tail_nets_routine * new_net_cost
      tail_mdc_cost <- tail_net_cost - tail_net_cost_routine
      mdc_budget <- max(0, override_cost_value - tail_net_cost_routine)
      override_cost_factor <- mdc_budget / tail_mdc_cost
      dist_used_future_mdc_only[new_net_range] <- dist_used_future_mdc_only[new_net_range] * override_cost_factor
      all_output_nets <- dist_used_no_future_mdc + dist_used_future_mdc_only
    } else {
      override_cost_factor <- override_cost_value / tail_net_cost
      all_output_nets[new_net_range] <- all_output_nets[new_net_range] * override_cost_factor
    }
  }
  
  # set bednets
  bednet_pars <- malariasimulation::set_bednets(site_pars,
                                                timesteps = output_net_times,
                                                coverages = all_output_nets,
                                                retention = mean_retu_dy,
                                                dn0 = dn0_mat,
                                                rn = rn_mat,
                                                rnm = rnm_mat,
                                                gamman = gam_vec)
  
  # run simulation
  output <- malariasimulation::run_simulation(timesteps = bednet_pars$timesteps,
                                              parameters = bednet_pars)
  
  N_timesteps <- bednet_pars$timesteps
  
  #TEST WITH NEW AGE RENDERING OUTPUTS FIRST!!!!
  #save outputs by future ITN years
  
  pfpr_min_bounds <- prevalence_rendering_min_ages[1:3]
  pfpr_max_bounds <- c(tail(prevalence_rendering_max_ages, 1),
                       prevalence_rendering_max_ages[3:4])
  
  n_age_cols <- paste("n_",
                      prevalence_rendering_min_ages, "_",
                      prevalence_rendering_max_ages, "_tot",
                      sep = "")
  detect_age_cols <- paste("n_detect_",
                           prevalence_rendering_min_ages, "_",
                           prevalence_rendering_max_ages, "_tot",
                           sep = "")
  case_age_cols <- paste("cases_",
                         incidence_rendering_min_ages, "_",
                         incidence_rendering_max_ages, "_tot",
                         sep = "")
  clin_case_age_cols <- paste("clin_cases_",
                              clinical_incidence_rendering_min_ages , "_",
                              clinical_incidence_rendering_max_ages, "_tot",
                              sep = "")
  sev_case_age_cols <- paste("sev_cases_",
                             clinical_incidence_rendering_min_ages , "_",
                             clinical_incidence_rendering_max_ages, "_tot",
                             sep = "")
  pfpr_age_cols <- paste("pfpr_",
                         pfpr_min_bounds, "_",
                         pfpr_max_bounds, "_mean",
                         sep = "")
  
  total_age_cols <- (
    length(n_age_cols) + length(detect_age_cols) +
      length(case_age_cols) + length(clin_case_age_cols) +
      length(sev_case_age_cols) + length(pfpr_age_cols)
  )
  
  annual_output <- as.data.frame(
    matrix(rep(0, projection_window_yr * (20 + total_age_cols)), nrow = 9)
  )
  all_summary_names <- c(
    "ISO2", "fs_name_1", "urbanicity", "fs_area", "fs_area_id", "net_name",
    "mass_int_yr", "biennial_reduction", "net_costings", "no_future_nets",
    "routine_baseline", "sample_id", "year_id", "CMC_start", "CMC_end",
    "ann_routine_nets_dist", "ann_camp_nets_dist",
    "adj_ann_routine_nets_dist", "adj_ann_camp_nets_dist", "net_strategy",
    n_age_cols, detect_age_cols, case_age_cols, clin_case_age_cols,
    sev_case_age_cols, pfpr_age_cols
  )
  if (no_future_nets) {net_name <- NA}
  if (no_future_nets | routine_baseline) {mass_int_yr <- NA}
  names(annual_output) <- all_summary_names
  annual_output$ISO2 <- rep(ISO2, projection_window_yr)
  annual_output$fs_name_1 <- rep(fs_name_1, projection_window_yr)
  annual_output$urbanicity <- rep(urbanicity, projection_window_yr)
  annual_output$fs_area <- rep(fs_area, projection_window_yr)
  annual_output$fs_area_id <- rep(fs_area_id, projection_window_yr)
  annual_output$net_name <- rep(net_name, projection_window_yr)
  annual_output$mass_int_yr <- rep(mass_int_yr, projection_window_yr)
  annual_output$biennial_reduction <- rep(biennial_reduction,
                                          projection_window_yr)
  annual_output$net_costings <- rep(net_costings, projection_window_yr)
  annual_output$no_future_nets <- rep(no_future_nets, projection_window_yr)
  annual_output$routine_baseline <- rep(routine_baseline, projection_window_yr)
  annual_output$sample_id <- rep(sid, projection_window_yr)
  annual_output$net_strategy <- rep(net_strategy, projection_window_yr)
  
  new_itn_campaigns_started <- FALSE
  yr <- 0
  
  m_proj_backfilled <- c(rep(-1, proj_camp_1 + N_early - 1), m_proj)
  new_MDCs <- which(m_proj_backfilled == 0)
  new_net_MDCs <- new_MDCs[new_MDCs >= new_net_start_mn]
  first_new_net_MDC <- min(new_net_MDCs)
  if (no_future_nets | routine_baseline) {first_new_net_MDC <- new_net_start_mn}
  for (t in 1:length(output_net_times)) {
    t0 <- output_net_times[t]
    if (t == length(output_net_times)) {
      t1 <- bednet_pars$timesteps
    } else {
      t1 <- output_net_times[t+1] - 1
    }
    if (t >= first_new_net_MDC) {
      # if (no_future_nets | routine_baseline) {
      #   new_itn_campaigns_started <- TRUE
      # } else if (!new_itn_campaigns_started & m_proj[t] == 0) {
      #   new_itn_campaigns_started <- TRUE
      # }
      if (no_future_nets | routine_baseline) {
        if (((t - new_net_start_mn) %% 12 ) == 0) {
          yr <- yr + 1
          annual_output[yr,"year_id"] <- yr
          annual_output[yr, "CMC_start"] <- CMC_Jan2000 + t - 1
          annual_output[yr, "CMC_end"] <- CMC_Jan2000 + t + 10
          annual_output[yr, "ann_camp_nets_dist"] <- 0
          if (no_future_nets) {
            annual_output[yr, "ann_routine_nets_dist"] <- 0
          } else {
            annual_output[yr,
                          "ann_routine_nets_dist"] <- ann_npc_D * sim_population
          }
          if (is.na(annual_output[yr, "ann_routine_nets_dist"])) {
            annual_output[yr, "ann_routine_nets_dist"] <- 0
          } else if (annual_output[yr, "ann_routine_nets_dist"] < 0) {
            annual_output[yr, "ann_routine_nets_dist"] <- 0
          }
          annual_output[yr, "adj_ann_routine_nets_dist"] <- (
            annual_output[yr, "ann_routine_nets_dist"] * 2 / exp(1)
          )
        }
      } else {
        #if ((m_proj[t - new_net_start_mn + 1] %% 12) == 0) {
        if ((m_proj[t - first_new_net_MDC + 1] %% 12) == 0) {
          yr <- yr + 1
          annual_output[yr,"year_id"] <- yr
          annual_output[yr, "CMC_start"] <- CMC_Jan2000 + t - 1
          annual_output[yr, "CMC_end"] <- CMC_Jan2000 + t + 11
          annual_output[yr, "ann_routine_nets_dist"] <- (
            ann_npc_D * sim_population
          )
          annual_output[yr, "adj_ann_routine_nets_dist"] <- (
            annual_output[yr, "ann_routine_nets_dist"] * 2 / exp(1)
          )
          tot_npc_this_yr <- sum(P_npc_dist[t:(t+11)])
          annual_output[yr, "ann_camp_nets_dist"] <- (
            (tot_npc_this_yr - ann_npc_D) * sim_population
          )
          annual_output[yr, "adj_ann_camp_nets_dist"] <- (
            annual_output[yr, "ann_camp_nets_dist"] * 2 / exp(1)
          )
          if (is.na(annual_output[yr, "ann_camp_nets_dist"])) {
            annual_output[yr, "ann_camp_nets_dist"] <- 0
          } else if (annual_output[yr, "ann_camp_nets_dist"] < 0) {
            annual_output[yr, "ann_camp_nets_dist"] <- 0
          }
        }
      }
      if (yr > 0) {
        for (a in 1:length(n_age_cols)) {
          annual_output[yr, n_age_cols[a]] <- (
            annual_output[yr, n_age_cols[a]] + sum(
              output[
                t0:t1,
                paste0(
                  "n_",
                  prevalence_rendering_min_ages[a],
                  "_",
                  prevalence_rendering_max_ages[a]
                )
              ]
            ) / 365
          )
        }
        for (a in 1:length(detect_age_cols)) {
          annual_output[yr, detect_age_cols[a]] <- (
            annual_output[yr, detect_age_cols[a]] + sum(
              output[
                t0:t1,
                paste0(
                  "n_detect_",
                  prevalence_rendering_min_ages[a],
                  "_",
                  prevalence_rendering_max_ages[a]
                )
              ]
            ) / 365
          )
        }
        for (a in 1:length(case_age_cols)) {
          annual_output[yr, case_age_cols[a]] <- (
            annual_output[yr, case_age_cols[a]] + sum(
              output[
                t0:t1,
                paste0(
                  "n_inc_",
                  incidence_rendering_min_ages[a],
                  "_",
                  incidence_rendering_max_ages[a]
                )
              ]
            )
          )
        }
        for (a in 1:length(clin_case_age_cols)) {
          annual_output[yr, clin_case_age_cols[a]] <- (
            annual_output[yr, clin_case_age_cols[a]] + sum(
              output[
                t0:t1,
                paste0(
                  "n_inc_clinical_",
                  clinical_incidence_rendering_min_ages[a],
                  "_",
                  clinical_incidence_rendering_max_ages[a]
                )
              ]
            )
          )
        }
        for (a in 1:length(sev_case_age_cols)) {
          annual_output[yr, sev_case_age_cols[a]] <- (
            annual_output[yr, sev_case_age_cols[a]] + sum(
              output[
                t0:t1,
                paste0(
                  "n_inc_severe_",
                  clinical_incidence_rendering_min_ages[a],
                  "_",
                  clinical_incidence_rendering_max_ages[a]
                )
              ]
            )
          )
        }
        for (tt in t0:t1) {
          n_tot <- 0
          n_detect_tot <- 0
          for (a in 1:length(detect_age_cols)) {
            n_tot <- n_tot + output[
              tt,
              paste0(
                "n_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
            n_detect_tot <- n_detect_tot + output[
              tt,
              paste0(
                "n_detect_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
          }
          annual_output[yr, "pfpr_0_36499_mean"] <- (
            annual_output[yr, "pfpr_0_36499_mean"] + n_detect_tot / n_tot / 365
          )
          n_tot <- 0
          n_detect_tot <- 0
          for (a in 2:3) {
            n_tot <- n_tot + output[
              tt,
              paste0(
                "n_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
            n_detect_tot <- n_detect_tot + output[
              tt,
              paste0(
                "n_detect_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
          }
          annual_output[yr, "pfpr_182_1824_mean"] <- (
            annual_output[yr, "pfpr_182_1824_mean"] + n_detect_tot / n_tot / 365
          )
          n_tot <- 0
          n_detect_tot <- 0
          for (a in 3:4) {
            n_tot <- n_tot + output[
              tt,
              paste0(
                "n_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
            n_detect_tot <- n_detect_tot + output[
              tt,
              paste0(
                "n_detect_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
          }
          annual_output[yr, "pfpr_730_3649_mean"] <- (
            annual_output[yr, "pfpr_730_3649_mean"] + n_detect_tot / n_tot / 365
          )
        }
        
      }
      
    }
  }
  
  annual_output <- annual_output[1:(projection_window_yr-3),]
  
  return(annual_output)
  
}

run_malsim_nets_sequential_v4 <- function(
    dataset,
    areas_included = NULL,
    N_cores = 0,
    areas_per_core = 1,
    sim_population = 1e5,
    ref_CMC = 1500,
    use_hipercow = FALSE,
    bv_beta = NULL,
    bv_gamma = NULL,
    local_sitefiles = FALSE,
    sitefile_folder = "./data_private/sitefiles/"
) {
  
  # Set number of cores
  #total_runs <- N_reps*length(areas_included)
  if (N_cores <= 0) {N_cores <- ceiling(length(areas_included)/areas_per_core)}
  if (use_hipercow) {max_cores <- 32} else {max_cores <- 20}
  if (max_cores > length(areas_included)) {max_cores <- length(areas_included)}
  if (N_cores > max_cores) {N_cores <- max_cores}
  
  # Simulation time
  CMC_sim_start <- CMC_Jan2000
  CMC_sim_end <- CMC_last
  N_CMC_old_nets <- CMC_last - CMC_sim_start + 1
  new_net_start_mn <- N_CMC_old_nets + 1
  N_CMC_sim <- CMC_sim_end - CMC_sim_start + 1
  
  # dataframe for storing output
  output_df <- data.frame(NULL)
  
  # Empty parameter list
  param_list <- list()
  
  # current country
  current_country <- "XXX"
  
  for (i in 1:N_fs_areas) {
    
    if (fs_id_link$fs_area[i] %in% areas_included) {
      
      # Warning for foresite mismatch
      if (fs_id_link$fs_area_id[i] != i) {
        print("Warning: Foresite id mismatch")
      }
      
      # Area-time indices
      area_id <- fs_id_link$new_area_id[i]
      area_time_ref_id <- which(dataset$area_id == area_id &
                                  dataset$CMC == ref_CMC)
      area_time_ids <- which(dataset$area_id == area_id)
      
      # get samples
      invlam_u_mid <- invlam_u[, area_id] %>%
        as.vector %>% unname %>% unlist %>% mean
      lam_u_mid <- 1 / invlam_u_mid
      ret_u_mid <- ret_u[, area_time_ref_id] %>%
        as.vector %>% unname %>% unlist %>% mean
      P_u_mid <- P_u[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      P0_u_mid <- P0_u[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      D_u_mid <- D_u[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      
      invlam_a_mid <- invlam_a[, area_id] %>%
        as.vector %>% unname %>% unlist %>% mean
      lam_a_mid <- 1 / invlam_a_mid
      ret_a_mid <- ret_a[, area_time_ref_id] %>%
        as.vector %>% unname %>% unlist %>% mean
      P_a_mid <- P_a[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      P0_a_mid <- P0_a[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      D_a_mid <- D_a[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      
      # Nets per capita samples
      npc_beta_mid <- bv_beta %>% mean
      npc_gamma_mid <- bv_gamma %>% mean
      
      # Identify last mass campaign
      # Estimated from latest peak in use
      # Previous approach:
      # last_camp_month <- apply(P_u_samples, 1, which.max)
      u0 <- 1
      u1 <- 0
      t <- length(P_u_mid)
      while ((u1 < u0 )&(t > 1)){
        u0 <- P_u_mid[t-1]
        u1 <- P_u_mid[t]
        if (u0 < u1) {last_camp_month <- t}
        t <- t - 1
      }
      
      # Generate ISO code for current admin
      admin_country <- countrycode(fs_id_link$ISO2[i], "iso2c", "iso3c")
      
      adm_site <- site::subset_site(
        site = ctry_site,
        site_filter = data.frame(
          country = countrycode(fs_id_link$ISO2[i], "iso2c", "country.name"),
          iso3c = admin_country,
          name_1 = fs_id_link$fs_name_1[i],
          urban_rural = fs_id_link$urbanicity[i]
        )
      )
      
      N_species <- dim(adm_site$vectors$vector_species)[1]
      
      # Pf EIR
      Pf_eir <- adm_site$eir$eir[1]

      if (Pf_eir > 0) {
        # Create parameter inputs

        site_pars <- site_parameters(
          interventions = adm_site$interventions,
          demography = adm_site$demography,
          vectors = adm_site$vectors$vector_species,
          seasonality = adm_site$seasonality$seasonality_parameters,
          eir = adm_site$eir$eir,
          overrides = list(
            human_population = sim_population,
            individual_mosquitoes = FALSE
          )
        )
        
       
        # Set age groups
        full_age_groups_min <- round(c(0.5, 2, 5) * 365)
        full_age_groups_max <- round(c(2, 5, 10) * 365)
        site_pars$prevalence_rendering_min_ages <- full_age_groups_min
        site_pars$prevalence_rendering_max_ages <- full_age_groups_max
        
        # Combine vector and matrix parameters for parLapply function
        site_pars$P_u_mid <- P_u_mid
        site_pars$P0_u_mid <- P0_u_mid
        site_pars$D_u_mid <- D_u_mid
        site_pars$lam_u_mid <- lam_u_mid
        site_pars$P_a_mid <- P_a_mid
        site_pars$P0_a_mid <- P0_a_mid
        site_pars$D_a_mid <- D_a_mid
        site_pars$lam_a_mid <- lam_a_mid
        site_pars$DOY_1st <- DOY_1st
        site_pars$DOY_mid <- DOY_mid
        
        # Combine single parameters
        site_pars$ISO2 <- fs_id_link$ISO2[i]
        site_pars$fs_area <- fs_id_link$fs_area[i]
        site_pars$fs_name_1 <- fs_id_link$fs_name_1[i]
        site_pars$urbanicity <- fs_id_link$urbanicity[i]
        site_pars$fs_area_id <- fs_id_link$fs_area_id[i]
        site_pars$N_species <- N_species
        site_pars$CMC_first <- CMC_first
        site_pars$CMC_Jan2000 <- CMC_Jan2000
        site_pars$N_CMC <- N_CMC
        site_pars$N_CMC_sim <- N_CMC_sim
        site_pars$sim_population <- sim_population

      

        site_pars$mean_retu <- ret_u_mid
        site_pars$mean_reta <- ret_a_mid
        site_pars$npc_beta <- npc_beta_mid
        site_pars$npc_gamma <- npc_gamma_mid
        
        param_list[[length(param_list) + 1]] <- site_pars
        
        
        
      }
      
      
      pc1 <- round(100 * ii / N_total_its)
      if (pc1 > pc0) {
        pc0 <- pc1
        print(paste(pc0, "% complete", sep = ""))
      }
      
    }
  }
  
  
  
  
  
  
  if (use_hipercow) {
    resources <- hipercow_resources(cores = N_cores)
    if (hiper_debug == TRUE) {
      par_id <- task_create_expr(par_net_region_sequential_v3(param_list[[1]]),
                                 resources = resources)
    } else {
      par_id <- task_create_expr(
        parallel::parLapply(NULL, param_list, par_net_region_sequential_v3),
        parallel = hipercow_parallel("parallel"),
        resources = resources)
    }
    
    
    # par_out <- task_result(par_output)
    # comb_output <- do.call(rbind.data.frame, par_out)
    # 
    # comb_output <- do.call(rbind.data.frame, par_output)
    # output_df <- rbind(output_df, comb_output)
  } else {
    #if (!use_hipercow) {
    cl <- makeCluster(N_cores)
    clusterExport(cl, c(#"param_list",
      "set_bednets",
      "run_simulation",
      "fit_usage_sequential_distributions",
      "population_usage_t"))
    #par_output <- lapply(param_list, par_net_region_sequential3)
    par_output <- parLapply(cl, param_list, par_net_region_sequential_v4)
    comb_output <- do.call(rbind.data.frame, par_output)
    output_df <- rbind(output_df, comb_output)
    stopCluster(cl)
  }
  
  if (use_hipercow) {
    return(par_id)
  } else {
    return(output_df)
  }
}