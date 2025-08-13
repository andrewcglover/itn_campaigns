# malsim.R
# Adapted from old sequential_v3
# 
# hipercow_country_loop <- function(cc,
#                                   mass_int_yr,
#                                   net_type,
#                                   no_future_nets,
#                                   routine_baseline,
#                                   net_costings)
#   #     assign(dynam_id,
#   #            task_create_expr(par_net_region_sequential3(hipercow_params)),
#   #            envir = .GlobalEnv)

par_net_region_sequential_v4 <- function(param_list) {
  
  #library(magrittr)
  
  #N_reps <- length(param_list)
  
  
  # Extract parameters from parameter list
  site_pars <- param_list#[[1]] #Uncomment for testing
  sid <- site_pars$sample_index
  mean_retu <- site_pars$mean_retu
  mean_reta <- site_pars$mean_reta
  npc_beta <- site_pars$npc_beta
  npc_gamma <- site_pars$npc_gamma
  net_type <- site_pars$net_type
  net_name <- site_pars$net_name
  net_strategy <- site_pars$net_strategy
  month_offset <- site_pars$month_offset
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
  projection_window_mn <- site_pars$projection_window_mn
  N_CMC <- site_pars$N_CMC
  N_CMC_sim <- site_pars$N_CMC_sim
  tail_pop <- site_pars$tail_pop
  sim_population <- site_pars$sim_population
  P_a_samples <- site_pars$P_a_samples
  P0_a_samples <- site_pars$P0_a_samples
  D_a_samples <- site_pars$D_a_samples
  lam_a_samples <- site_pars$lam_a_samples
  P_u_samples <- site_pars$P_u_samples
  P0_u_samples <- site_pars$P0_u_samples
  D_u_samples <- site_pars$D_u_samples
  lam_u_samples <- site_pars$lam_u_samples
  dn0_mat <- site_pars$dn0_mat
  rn_mat <- site_pars$rn_mat
  rnm_mat <- site_pars$rnm_mat
  gam_vec <- site_pars$gam_vec
  DOY_1st <- site_pars$DOY_1st
  DOY_mid <- site_pars$DOY_mid
  net_costings <- site_pars$net_costings
  cost_factor <- site_pars$cost_factor
  biennial_reduction <- site_pars$biennial_reduction
  routine_baseline <- site_pars$routine_baseline
  new_net_cost <- site_pars$new_net_cost
  no_future_nets <- site_pars$no_future_nets
  override_cost <- site_pars$override_cost
  override_mdc_only <- site_pars$override_mdc_only
  override_cost_value <- site_pars$override_cost_value
  new_net_start_mn <- site_pars$new_net_start_mn
  prevalence_rendering_min_ages <- site_pars$prevalence_rendering_min_ages
  prevalence_rendering_max_ages <- site_pars$prevalence_rendering_max_ages
  incidence_rendering_min_ages <- site_pars$incidence_rendering_min_ages
  incidence_rendering_max_ages <- site_pars$incidence_rendering_max_ages
  clinical_incidence_rendering_min_ages <- site_pars$clinical_incidence_rendering_min_ages
  clinical_incidence_rendering_max_ages <- site_pars$clinical_incidence_rendering_max_ages
  severe_incidence_rendering_min_ages <- site_pars$severe_incidence_rendering_min_ages
  severe_incidence_rendering_max_ages <- site_pars$severe_incidence_rendering_max_ages
  EIR_fit <- site_pars$EIR_fit
  budget_pc <- site_pars$budget_pc
  
  # if (biennial_reduction & (mass_int_mn < 25)) {
  #   net_strategy <- paste0(net_strategy, "_bien_costed")
  # } else if ((cost_factor < 0.9999) | (cost_factor > 1.0001)) {
  #   net_strategy <- paste0(net_strategy, "_costed")
  # }
  
  if (no_future_nets) {biennial_reduction <- FALSE}
  if (no_future_nets) {net_costings <- FALSE}
  if (no_future_nets) {routine_baseline <- FALSE}
  if (routine_baseline) {biennial_reduction <- FALSE}
  if (routine_baseline) {no_future_nets <- FALSE}
  
  
  mass_int_yr <- round(mass_int_mn / 12)
  
  net_cost_logical <- 1 * net_costings
  if (mass_int_mn != 24) {biennial_reduction <- FALSE}
  biennial_reduction_logical <- 1 * biennial_reduction
  routine_baseline_logical <- 1 * routine_baseline
  
  fs_area_undrscr <- gsub(" ", "_", fs_area)
  
  year = 365
  obs_window = projection_window_mn * 365 / 12
  projection_window_yr <- round(projection_window_mn / 12)
  
  # Monthly rates
  decay_a <- 1 / mean_reta
  decay_u <- 1 / mean_retu
  decay_uga <- decay_u - decay_a
  decay_npc <- decay_a / npc_gamma
  
  # convert retention to days
  mean_retu_dy <- 365 * mean_retu / 12
  mean_reta_dy <- 365 * mean_reta / 12
  
  # Central time point for first regular mass campaign
  proj_camp_1 <- last_camp + mass_int_mn + month_offset
  
  # Define period from first simulated campaign (including projection)
  proj_end <- N_CMC + projection_window_mn
  N_proj <- proj_end - proj_camp_1 + 1
  t_proj <- seq(1, N_proj)
  m_proj <- (t_proj - 1) %% mass_int_mn
  m_long <- seq(1, proj_end) - last_camp
  m_tail <- m_long[last_camp:proj_end]
  
  # Back-fill preliminary period
  N_early <- CMC_first - CMC_Jan2000
  
  # New net times
  new_start <- N_early + N_CMC + 1
  new_end <- new_start + projection_window_mn - 1
  
  # Times for fitting
  times_mn <- seq(1, proj_end)
  times_yr <- rep(seq(0, ceiling(N_CMC_sim / 12)), each=12)
  times_1st_dy <- DOY_1st + (times_yr * year)
  times_mid_dy <- DOY_mid + (times_yr * year)
  input_net_times <- times_mid_dy[1:N_CMC_sim]    # usage for fitting
  output_net_times <- times_1st_dy[1:N_CMC_sim]   # distribution times for netz
  
  #-----------------------------------------------------------------------------
  # usage profiles
  
  # Extract values for selected sample
  P_ui <- P_u_samples[sid,]
  P0_ui <- P0_u_samples[sid,]
  D_ui <- D_u_samples[sid,]
  lambda_ui <- lam_u_samples[sid]
  
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
  P_ai <- P_a_samples[sid,]
  P0_ai <- P0_a_samples[sid,]
  D_ai <- D_a_samples[sid,]
  lambda_ai <- lam_a_samples[sid]
  
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
  bednet_pars <- malariasimulation::set_equilibrium(bednet_pars, init_EIR = EIR_fit)
  output <- malariasimulation::run_simulation(timesteps = bednet_pars$timesteps,
                                              parameters = bednet_pars)
  
  N_timesteps <- bednet_pars$timesteps
  
  #TEST WITH NEW AGE RENDERING OUTPUTS FIRST!!!!
  #save outputs by future ITN years
  
  pfpr_min_bounds <- prevalence_rendering_min_ages[1:3]
  pfpr_max_bounds <- c(tail(prevalence_rendering_max_ages, 1),
                       prevalence_rendering_max_ages[3:4])
  
  n_age_cols <- paste("n_age_",
                      prevalence_rendering_min_ages, "_",
                      prevalence_rendering_max_ages, "_tot",
                      sep = "")
  detect_age_cols <- paste("n_detect_lm_",
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
    matrix(rep(0, projection_window_yr * (21 + total_age_cols)), nrow = 9)
  )
  all_summary_names <- c(
    "ISO2", "fs_name_1", "urbanicity", "fs_area", "fs_area_id", "net_name",
    "mass_int_yr", "biennial_reduction", "net_costings", "no_future_nets",
    "routine_baseline", "sample_id", "year_id", "CMC_start", "CMC_end",
    "ann_routine_nets_dist", "ann_camp_nets_dist",
    "adj_ann_routine_nets_dist", "adj_ann_camp_nets_dist", "net_strategy",
    "budget_pc",
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
  annual_output$budget_pc <- rep(budget_pc, projection_window_yr)
  
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
                  "n_age_",
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
                  "n_detect_lm_",
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
          n_detect_lm_tot <- 0
          for (a in 1:length(detect_age_cols)) {
            n_tot <- n_tot + output[
              tt,
              paste0(
                "n_age_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
            n_detect_lm_tot <- n_detect_lm_tot + output[
              tt,
              paste0(
                "n_detect_lm_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
          }
          annual_output[yr, "pfpr_0_36499_mean"] <- (
            annual_output[yr, "pfpr_0_36499_mean"] +
              n_detect_lm_tot / n_tot / 365
          )
          n_tot <- 0
          n_detect_lm_tot <- 0
          for (a in 2:3) {
            n_tot <- n_tot + output[
              tt,
              paste0(
                "n_age_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
            n_detect_lm_tot <- n_detect_lm_tot + output[
              tt,
              paste0(
                "n_detect_lm_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
          }
          annual_output[yr, "pfpr_182_1824_mean"] <- (
            annual_output[yr, "pfpr_182_1824_mean"] +
              n_detect_lm_tot / n_tot / 365
          )
          n_tot <- 0
          n_detect_lm_tot <- 0
          for (a in 3:4) {
            n_tot <- n_tot + output[
              tt,
              paste0(
                "n_age_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
            n_detect_lm_tot <- n_detect_lm_tot + output[
              tt,
              paste0(
                "n_detect_lm_",
                prevalence_rendering_min_ages[a],
                "_",
                prevalence_rendering_max_ages[a]
              )
            ]
          }
          annual_output[yr, "pfpr_730_3649_mean"] <- (
            annual_output[yr, "pfpr_730_3649_mean"] +
              n_detect_lm_tot / n_tot / 365
          )
        }
        
      }
      
    }
  }
  
  annual_output <- annual_output[1:(projection_window_yr-3),]
  
  return(annual_output)
  
}

run_malsim_nets_sequential_v4 <- function(dataset,
                                          areas_included = NULL,
                                          N_reps = 100,
                                          N_cores = 0,
                                          areas_per_core = 1,
                                          mass_int_yr = c(2,3),
                                          sim_population = 1e5,
                                          ref_CMC = 1500,
                                          only = FALSE,
                                          pbo = FALSE,
                                          pyrrole = FALSE,
                                          net_costings = FALSE,
                                          biennial_reduction = FALSE,
                                          month_default_offset = 0,
                                          rep_offset = 0,
                                          use_hipercow = FALSE,
                                          routine_baseline = FALSE,
                                          no_future_nets = FALSE,
                                          debugging = FALSE,
                                          override_cost = FALSE,
                                          override_mdc_only = FALSE,
                                          override_cost_value = 1,
                                          hiper_debug = FALSE,
                                          bv_beta = NULL,
                                          bv_gamma = NULL,
                                          CMC_last = 1500,
                                          projection_window_yr = 9,
                                          budget_pc = 100) {
  
  # Set number of cores
  # N_cores takes precedence over areas_per core if N_cores set to > 0
  total_runs <- N_reps*length(areas_included)
  if (N_cores <= 0) {N_cores <- ceiling(length(areas_included)/areas_per_core)}
  if (use_hipercow) {max_cores <- 32} else {max_cores <- 18}
  if (max_cores > total_runs) {max_cores <- total_runs}
  if (N_cores > max_cores) {N_cores <- max_cores}
  
  # Projection window
  projection_window_mn <- projection_window_yr * 12
  projection_window_dy <- projection_window_yr * 365
  
  # Simulation time
  CMC_sim_start <- CMC_Jan2000
  CMC_sim_end <- CMC_last + projection_window_mn
  N_CMC_old_nets <- CMC_last - CMC_sim_start + 1
  new_net_start_mn <- N_CMC_old_nets + 1
  N_CMC_sim <- CMC_sim_end - CMC_sim_start + 1
  
  # Number of samples
  N_samples <- dim(P_u)[1]
  
  # Create sample ids
  if (max(long_sample_ids) > N_samples) {
    print("Warning: Some sample ids outwith range")
  }
  rep_ids <- seq(1, N_reps) + rep_offset
  sample_ids <- long_sample_ids[1:N_reps]
  npc_sample_ids <- round(dim(bv_beta)[1] * sample_ids / max(long_sample_ids))
  
  # dataframe for storing output
  output_df <- data.frame(NULL)
  
  # progress indicator
  N_net_types <- only + pbo + pyrrole
  N_int_vals <- length(mass_int_yr)
  N_areas_included <- length(areas_included)
  N_total_its <- N_net_types * N_int_vals * N_areas_included
  pc0 <- 0
  ii <- 0
  
  # Empty parameter list
  param_list <- list()
  
  for (l in 1:3) {
    
    #if ((l==1 & only & !net_costings) | (l==2 & pbo) | (l==3 & pyrrole)) {
    if ((l==1 & only) | (l==2 & pbo) | (l==3 & pyrrole)) {
      
      if (net_costings) {
        if (l==1 & only) {cost_factor <- budget_pc / 100}
        if (l==2 & pbo) {
          cost_factor <- scaled_pbo_nets_equiv_only * budget_pc / 100
        }
        if (l==3 & pyrrole) {
          cost_factor <- scaled_pyrrole_nets_equiv_only * budget_pc / 100
        }
      } else {
        cost_factor <- 1.0
      }
      
      if (l==1 & only) {new_net_cost <- only_total_cost}
      if (l==2 & pbo) {new_net_cost <- pbo_total_cost}
      if (l==3 & pyrrole) {new_net_cost <- pyrrole_total_cost}
      
      old_dat_res <- dat_res_pyr_ll
      if (l==1 & only) {new_dat_res <- dat_res_pyr}
      if (l==2 & pbo) {new_dat_res <- dat_res_pbo}
      if (l==3 & pyrrole) {new_dat_res <- dat_res_pp}
      
      # Correct for precision rounding errors
      old_dat_res$resistance %<>% round(2)
      new_dat_res$resistance %<>% round(2)
      
      if (l==1 & only) {net_name <- "Pyrethroid-only"}
      if (l==2 & pbo) {net_name <- "Pyrethroid-PBO"}
      if (l==3 & pyrrole) {net_name <- "Pyrethroid-Pyrrole"}
      
      for (k in 1:N_int_vals) {
        
        # mass interval
        mass_int_mn <- mass_int_yr[k] * 12
        
        # current country
        current_country <- "XXX"
        
        for (i in 1:N_fs_areas) {
          
          if (fs_id_link$fs_area[i] %in% areas_included) {
            
            # Warning for foresite mismatch
            if (fs_id_link$fs_area_id[i] != i) {
              print("Warning: Foresite id mismatch")
            }
            
            ii <- ii + 1
            
            # Area-time indices
            area_id <- fs_id_link$new_area_id[i]
            area_time_ids <- which(usage_list$a == area_id)
            area_time_ref_id <- area_time_ids[ref_CMC - CMC_first + 1]
            #area_time_ref_id <- which(dataset$area_id == area_id &
            #                            dataset$CMC == ref_CMC)
            #area_time_ids <- which(dataset$area_id == area_id)
            
            # get samples
            invlam_u_samples <- invlam_u[sample_ids, area_id] %>%
              as.vector %>% unname %>% unlist
            lam_u_samples <- 1 / invlam_u_samples
            ret_u_samples <- ret_u[sample_ids, area_time_ref_id] %>%
              as.vector %>% unname %>% unlist
            P_u_samples <- P_u[sample_ids, area_time_ids] %>%
              as.matrix %>% unname
            P0_u_samples <- P0_u[sample_ids, area_time_ids] %>%
              as.matrix %>% unname
            D_u_samples <- D_u[sample_ids, area_time_ids] %>%
              as.matrix %>% unname
            
            invlam_a_samples <- invlam_a[sample_ids, area_id] %>%
              as.vector %>% unname %>% unlist
            lam_a_samples <- 1 / invlam_a_samples
            ret_a_samples <- ret_a[sample_ids, area_time_ref_id] %>%
              as.vector %>% unname %>% unlist
            P_a_samples <- P_a[sample_ids, area_time_ids] %>%
              as.matrix %>% unname
            P0_a_samples <- P0_a[sample_ids, area_time_ids] %>%
              as.matrix %>% unname
            D_a_samples <- D_a[sample_ids, area_time_ids] %>%
              as.matrix %>% unname
            
            # Nets per capita samples
            npc_beta_samples <- bv_beta[npc_sample_ids]
            npc_gamma_samples <- bv_gamma[npc_sample_ids]
            
            # Identify last mass campaign
            # Estimated from latest peak in use
            # Previous approach:
            # last_camp_month <- apply(P_u_samples, 1, which.max)
            u0 <- 1
            u1 <- 0
            t <- dim(P_u_samples)[2]
            while ((u1 < u0 )&(t > 1)){
              u0 <- P_u_samples[1,(t-1)]
              u1 <- P_u_samples[1,t]
              if (u0 < u1) {last_camp_month <- t}
              t <- t - 1
            }
            
            # First regular MDC with new nets
            #first_new_net_month <- last_camp_month + mass_int_mn + 
            
            # Generate ISO code for current admin
            admin_country <- countrycode(fs_id_link$ISO2[i], "iso2c", "iso3c")
            
            # Pull country site
            if (admin_country != current_country) {
              current_country <- admin_country
              if (local_sitefiles) {
                ctry_site <- readRDS(
                  eval(paste0(sitefile_folder, current_country, "_site.rds"))
                )
              } else {
                ctry_site <- site::fetch_site(current_country)
              }
            }
            
            # Isolate single admin site
            adm_site <- site::subset_site(
              site = ctry_site,
              site_filter = data.frame(
                country = countrycode(fs_id_link$ISO2[i], "iso2c", "country.name"),
                iso3c = admin_country,
                name_1 = fs_id_link$fs_name_1[i],
                urban_rural = fs_id_link$urbanicity[i]
              )
            )
            
            # # Pull country site
            # if (admin_country != current_country) {
            #   current_country <- admin_country
            #   ctry_site <- get_site(current_country)
            # }
            # 
            # # Isolate a single site from a country
            # adm_site_index <- which(ctry_site$sites$name_1 == fs_id_link$fs_name_1[i] &
            #                           ctry_site$sites$urban_rural == fs_id_link$urbanicity[i])
            # 
            # # If no foresite file for urban/rural, then revert to other
            # if (identical(adm_site_index, integer(0))) {
            #   if (fs_id_link$urbanicity[i] == "urban") {
            #     adm_site_index <- which(ctry_site$sites$name_1 == fs_id_link$fs_name_1[i] &
            #                               ctry_site$sites$urban_rural == "rural")
            #   } else {
            #     adm_site_index <- which(ctry_site$sites$name_1 == fs_id_link$fs_name_1[i] &
            #                               ctry_site$sites$urban_rural == "urban")
            #   }
            # }
            # 
            # # Check for successful foresite match
            # if (identical(adm_site_index, integer(0))) {
            #   print(paste0("Warning: foresite not linked for admin region ",
            #                fs_id_link$fs_name_1[i], " (index ", i, ") in ",
            #                current_country))
            # }
            # 
            # # Create admin site file
            # adm_site <- site::single_site(ctry_site, adm_site_index)
            
            # Repeat interventions
            adm_site %<>% expand_interventions(expand_year = projection_window_yr,
                                               delay = 0,
                                               counterfactual = FALSE)
            
            # TESTING NEW SITE FILE FORMAT
            # new_demography <- readRDS("NGA_new_demography.rds")
            # new_demography %<>%
            #   dplyr::rename(mortality_rate = adjusted_mortality_rates)
            # 
            # old_adm_site_demography <- adm_site$demography
            # adm_site$demography <- new_demography
            
            # Tail population
            tail_ids <- which(adm_site$demography$year >= 2025 &
                                adm_site$demography$year <= (2025 + projection_window_yr))
            tail_pop <- mean(adm_site$demography$population[tail_ids])
            
            # Pyrethroid resistance
            yearly_res <- adm_site$vectors$pyrethroid_resistance$pyrethroid_resistance
            monthly_res <- rep(yearly_res, each = 12)
            round_monthly_res <- round(monthly_res, 2)
            
            # Name net strategy
            net_strategy <- net_name
            # allow for costed routine
            if (routine_baseline) {
              net_strategy %<>% paste("routine baseline")
            } else {
              net_strategy %<>% paste(mass_int_yr[k], "year interval")
            }
            if (override_cost) {
              if (override_mdc_only) {
                net_strategy %<>% paste("mdc")
              } else {
                net_strategy %<>% paste("all")
              }
              net_strategy %<>% paste("overriden costed")
            } else if (net_costings) {
              if (biennial_reduction & mass_int_yr[k] == 2) {
                net_strategy %<>% paste("biennial and type costed")
              } else {
                net_strategy %<>% paste("type costed")
              }
            } else {
              if (biennial_reduction & mass_int_yr[k] == 2) {
                net_strategy %<>% paste("biennial costed")
              }
            }
            if (budget_pc != 100) {
              net_strategy %<>% paste0(" ", budget_pc, "pc budget")
            }
            
            # Manually change Kedougou EIR = 250 from:
            #Dia, I., et al., Bionomics of Anopheles gambiae Giles, An. arabiensis Patton, An. funestus Giles and An.
            #nili (Theobald) (Diptera: Culicidae) and transmission of Plasmodium falciparum in a Sudano-Guinean
            #zone (Ngari, Senegal). Journal of medical entomology, 2003. 40(3).
            # if (adm_site$eir$name_1[1] == "Kédougou") {
            #   # adm_site$eir$eir[1] <- 250
            #   # Set to highest estimated PfEIR of other admin regions
            #   adm_site$eir$eir[1] <- max(ctry_site$eir$eir[SN_site$eir$spp == "pf"])
            # }
            
            # Pf EIR
            adm_site$eir$eir[1] <- fs_id_link$EIR_fit[i]
            Pf_eir <- adm_site$eir$eir[1]
            
            # Update bednet interventions
            
            if (Pf_eir > 0) {
              # Create parameter inputs
              site_pars <- site::site_parameters(
                interventions = adm_site$interventions,
                demography = adm_site$demography,
                vectors = adm_site$vectors$vector_species,
                seasonality = adm_site$seasonality$seasonality_parameters,
                eir = Pf_eir,
                overrides = list(human_population = sim_population,
                                 individual_mosquitoes = FALSE)
              )
              
              # Set age groups
              full_age_groups_min <- round(c(0, 0.5, 2, seq(5, 65, 5)) * 365)
              full_age_groups_max <- full_age_groups_min[2:length(full_age_groups_min)]
              full_age_groups_max <- c(full_age_groups_max - 1, 36499)
              site_pars$prevalence_rendering_min_ages <- full_age_groups_min
              site_pars$prevalence_rendering_max_ages <- full_age_groups_max
              site_pars$incidence_rendering_min_ages <- full_age_groups_min
              site_pars$incidence_rendering_max_ages <- full_age_groups_max
              site_pars$clinical_incidence_rendering_min_ages <- full_age_groups_min
              site_pars$clinical_incidence_rendering_max_ages <- full_age_groups_max
              site_pars$severe_incidence_rendering_min_ages <- full_age_groups_min
              site_pars$severe_incidence_rendering_max_ages <- full_age_groups_max
              site_pars$age_group_rendering_min_ages <- full_age_groups_min
              site_pars$age_group_rendering_max_ages <- full_age_groups_max
              
              
              # Combine vector and matrix parameters for parLapply function
              site_pars$P_u_samples <- P_u_samples
              site_pars$P0_u_samples <- P0_u_samples
              site_pars$D_u_samples <- D_u_samples
              site_pars$lam_u_samples <- lam_u_samples
              site_pars$P_a_samples <- P_a_samples
              site_pars$P0_a_samples <- P0_a_samples
              site_pars$D_a_samples <- D_a_samples
              site_pars$lam_a_samples <- lam_a_samples
              site_pars$DOY_1st <- DOY_1st
              site_pars$DOY_mid <- DOY_mid
              
              # Combine single parameters
              site_pars$budget_pc <- budget_pc
              site_pars$net_type <- l
              site_pars$net_name <- net_name
              site_pars$net_strategy <- net_strategy
              site_pars$last_camp <- last_camp_month#[j]
              site_pars$mass_int_mn <- mass_int_mn
              site_pars$ISO2 <- fs_id_link$ISO2[i]
              site_pars$fs_area <- fs_id_link$fs_area[i]
              site_pars$fs_name_1 <- fs_id_link$fs_name_1[i]
              site_pars$urbanicity <- fs_id_link$urbanicity[i]
              site_pars$fs_area_id <- fs_id_link$fs_area_id[i]
              site_pars$N_species <- N_species
              site_pars$CMC_first <- CMC_first
              site_pars$CMC_Jan2000 <- CMC_Jan2000
              site_pars$projection_window_mn <- projection_window_mn
              site_pars$N_CMC <- N_CMC
              site_pars$N_CMC_sim <- N_CMC_sim
              site_pars$tail_pop <- tail_pop
              site_pars$sim_population <- sim_population
              site_pars$net_costings <- net_costings
              site_pars$cost_factor <- cost_factor
              site_pars$biennial_reduction <- biennial_reduction
              site_pars$routine_baseline <- routine_baseline
              site_pars$new_net_cost <- new_net_cost
              site_pars$no_future_nets <- no_future_nets
              site_pars$override_cost <- override_cost
              site_pars$override_mdc_only <- override_mdc_only
              site_pars$override_cost_value <- override_cost_value
              site_pars$new_net_start_mn <- new_net_start_mn
              site_pars$EIR_fit <- fs_id_link$EIR_fit[i]
              
              for (j in 1:N_reps) {
                
                jj <- j + rep_offset
                site_pars$sample_index <- jj
                site_pars$month_offset <- long_month_offset[jj]
                site_pars$mean_retu <- ret_u_samples[jj]
                site_pars$mean_reta <- ret_a_samples[jj]
                site_pars$npc_beta <- npc_beta_samples[jj]
                site_pars$npc_gamma <- npc_gamma_samples[jj]
                
                # New Resistance sampling - added on 20/12/24
                N_species <- length(adm_site$vectors$vector_species$species)
                
                dn0_old <- rep(NA, length(monthly_res))
                dn0_vec <- rep(NA, length(monthly_res))
                rn_old <- rep(NA, length(monthly_res))
                rn_vec <- rep(NA, length(monthly_res))
                gam_old <- rep(NA, length(monthly_res))
                gam_vec <- rep(NA, length(monthly_res))
                
                comb_dat_res_sample <- old_dat_res[0,]
                
                res_kk <- old_dat_res[old_dat_res$resistance == round_monthly_res[1],]
                
                comb_dat_res_sample <- old_dat_res[0,]
                
                res_kk <- old_dat_res[old_dat_res$resistance == round_monthly_res[1],]
                
                for (kk in 1:N_CMC_sim) {
                  
                  x <- (jj-1) * N_CMC_sim + kk
                  y <- long_sample_ids[(x-1) %% length(long_sample_ids) + 1]
                  
                  if (kk < N_CMC_old_nets) {
                    
                    if (kk > 1) {
                      if (round_monthly_res[kk] != round_monthly_res[kk-1]) {
                        res_kk <- old_dat_res[old_dat_res$resistance == round_monthly_res[kk],]
                      }
                    }
                    
                    sample_kk <- res_kk[res_kk$draw == y,]
                    
                    if (dim(sample_kk)[1] < 1) {
                      print(paste("Warning:",
                                  round_monthly_res[kk],
                                  "resistance match fail"))
                      
                    }
                    
                  } else {
                    
                    if (kk > 1) {
                      if (round_monthly_res[kk] != round_monthly_res[kk-1]) {
                        res_kk <- new_dat_res[new_dat_res$resistance == round_monthly_res[kk],]
                      }
                    }
                    
                    sample_kk <- res_kk[res_kk$draw == y,]
                    
                    if (dim(sample_kk)[1] < 1) {
                      print(paste("Warning:",
                                  round_monthly_res[kk],
                                  "resistance match fail"))
                    }
                  }
                  comb_dat_res_sample %<>% rbind.data.frame(sample_kk)
                }
                
                dn0_mat <- matrix(rep(comb_dat_res_sample$dn0, N_species),
                                  nrow = N_CMC_sim,
                                  ncol = N_species)
                
                rn_mat <- matrix(rep(comb_dat_res_sample$rn0, N_species),
                                 nrow = N_CMC_sim,
                                 ncol = N_species)
                
                rnm_mat <- matrix(rep(0.24, N_CMC_sim * N_species),
                                  nrow = N_CMC_sim,
                                  ncol = N_species)
                
                gam_vec <- 365 * comb_dat_res_sample$gamman / log(2)
              
                site_pars$dn0_mat <- dn0_mat
                site_pars$rn_mat <- rn_mat
                site_pars$rnm_mat <- rnm_mat
                site_pars$gam_vec <- gam_vec
                
                param_list[[length(param_list) + 1]] <- site_pars
                
              }
              
            }
            
            
            pc1 <- round(100 * ii / N_total_its)
            if (pc1 > pc0) {
              pc0 <- pc1
              print(paste(pc0, "% complete", sep = ""))
            }
            
          }
        }
        
      }
    }
  }
  

  if (use_hipercow) {
    resources <- hipercow_resources(cores = N_cores)
    if (hiper_debug == TRUE) {
      par_id <- task_create_expr(par_net_region_sequential_v4(param_list[[1]]),
                                 resources = resources)
    } else {
      par_id <- task_create_expr(
        parallel::parLapply(NULL, param_list, par_net_region_sequential_v4),
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