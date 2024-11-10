# simulation_repeat.R
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

par_net_region_sequential_costed <- function(param_list) {
  
  #library(magrittr)
  
  #N_reps <- length(param_list)
  
  
  # Extract parameters from parameter list
  site_pars <- param_list#[[1]]
  sid <- site_pars$sample_index
  mean_retention <- site_pars$mean_ret
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
  P_samples <- site_pars$P_samples
  P0_samples <- site_pars$P0_samples
  D_samples <- site_pars$D_samples
  lam_samples <- site_pars$lam_samples
  dn0_mat <- site_pars$dn0_mat
  rn_mat <- site_pars$rn_mat
  rnm_mat <- site_pars$rnm_mat
  gam_vec <- site_pars$gam_vec
  DOY_1st <- site_pars$DOY_1st
  DOY_mid <- site_pars$DOY_mid
  net_costings <- site_pars$net_costings
  cost_factor <- site_pars$cost_factor
  interval_reduction <- site_pars$interval_reduction
  routine_baseline <- site_pars$routine_baseline
  new_net_cost <- site_pars$new_net_cost
  no_future_nets <- site_pars$no_future_nets
  override_cost <- site_pars$override_cost
  override_mdc_only <- site_pars$override_mdc_only
  override_cost_value <- site_pars$override_cost_value
  ref_usage_D <- site_pars$ref_usage_D
  ref_usage_C <- site_pars$ref_usage_C
  ref_usage <- site_pars$ref_usage
  lam_a_samples <- site_pars$lam_a_samples
  ref_access_D <- site_pars$ref_access_D
  ref_access_C <- site_pars$ref_access_C
  ref_access <- site_pars$ref_access
  bv_access <- site_pars$bv_access
  bv_npc <- site_pars$bv_npc
  output_sim <- site_pars$output_sim
  camp_offset_from_peak <- site_pars$camp_offset_from_peak
  proj_camp_override_month <- site_pars$proj_camp_override_month
  N_backfill <- site_pars$N_backfill
  #peak_season_dy <- site_pars$peak_season_dy
  
  # if (interval_reduction & (mass_int_mn < 25)) {
  #   net_strategy <- paste0(net_strategy, "_bien_costed")
  # } else if ((cost_factor < 0.9999) | (cost_factor > 1.0001)) {
  #   net_strategy <- paste0(net_strategy, "_costed")
  # }
  
  net_cost_logical <- 1 * net_costings
  #if (mass_int_mn > 25) {interval_reduction <- FALSE}
  interval_reduction_logical <- 1 * interval_reduction
  routine_baseline_logical <- 1 * routine_baseline
  
  fs_area_undrscr <- gsub(" ", "_", fs_area)
  
  year = 365
  obs_window = 12 * year
  
  # convert retention to days
  mean_retention_dy <- 365 * mean_retention / 12
  
  # Central time point for first regular mass campaign
  if (camp_offset_from_peak) {
  #   peak_month <- round(12*peak_season_dy/year)
  #   camp_month_of_year <- (peak_month + month_offset_from_peak) %% 12
  #   last_camp_month_index <- last_camp %% 12
  #   proj_camp_1 <- 12 * (last_camp %/% 12) + mass_int_mn + camp_month_of_year
  #   proj_camp_1_int <- proj_camp_1 - last_camp
  #   if (proj_camp_1_int < (mass_int_mn - 6)) {proj_camp_1 <- proj_camp_1 + 12}
  #   if (proj_camp_1_int > (mass_int_mn + 6)) {proj_camp_1 <- proj_camp_1 - 12}
    proj_camp_1 <- proj_camp_override_month
  } else {
    proj_camp_1 <- last_camp + mass_int_mn + month_offset
  }
  
  # Define period from first simulated campaign (including projection)
  proj_end <- N_CMC + projection_window_mn
  N_proj <- proj_end - proj_camp_1 + 1
  t_proj <- seq(1, N_proj)
  m_proj <- (t_proj - 1) %% mass_int_mn
  
  # Extract values for selected sample
  P <- P_samples[sid,]
  P0 <- P0_samples[sid,]
  D <- D_samples[sid,]
  lambda <- lam_samples[sid]
  lambda_a <- lam_a_samples[sid]
  Du_ref <- ref_usage_D[sid]
  Cu_ref <- ref_usage_C[sid]
  u_ref <- ref_usage[sid]
  Da_ref <- ref_access_D[sid]
  Ca_ref <- ref_access_C[sid]
  a_ref <- ref_access[sid]
  
  # End values
  P0_end <- tail(P0, n = 1)
  D_end <- tail(D, n = 1)
  C0_end <- P0_end - D_end
  
  # Approx average use given access and usage vs nets per capita conversion
  uga_ref <- min(1.0, u_ref / a_ref)
  bv_use <- bv_access * uga_ref
  saveRDS(bv_access, "val_test.rds")
  npc_ref <- bv_npc[which(abs(u_ref - bv_use) == (min(abs(u_ref - bv_use))))]
  npc_D <- bv_npc[which(abs(D_end - bv_use) == (min(abs(D_end - bv_use))))]
  
  # Costed nets per capita
  #npc_scaled <- (C0_end-(36/mass_int_mn)*exp(-mass_int_mn*lam))/(C0_end-exp(-36*lambda))
  #npc_new <- npc_D + (mass_int_mn/36.0) * (npc_ref-npc_D) #interval adjustment
  npc_costed <- cost_factor * npc_ref #new net cost adjustment
  
  # Calculate new costed use
  avg_use_costed <- bv_use[which(abs(npc_costed-bv_npc) == (min(abs(npc_costed-bv_npc))))]
  avg_C_costed <- avg_use_costed - D_end
  C0_costed <- (avg_C_costed * lambda * 36) / (1 - exp(-lambda * 36))
  if (interval_reduction) {
    C0_new_numer <- C0_costed * mass_int_mn * (1 - exp(- lambda * 36))
    C0_new_denom <- 36 * (1 - exp(- lambda * mass_int_mn))
    C0_new <- C0_new_numer / C0_new_denom
  } else {
    C0_new <- C0_costed
  }
  P0_new <- C0_new + D_end
  if (P0_new > 1) {
    P0_new <- 1
    C0_new <- P0_new - D_end
  }
  saveRDS(P0_new, "P0_new.rds")
  
  # Extend values
  P0_long <- c(P0, rep(P0_new, projection_window_mn))#P0_new, projection_window_mn))
  saveRDS(P0_long, "P0_test.rds")
  D_long <- c(D, rep(D_end, projection_window_mn))
  m_long <- seq(1, proj_end) - last_camp
  m_tail <- m_long[last_camp:proj_end]
  P0_tail <- P0_long[last_camp:proj_end]
  D_tail <- D_long[last_camp:proj_end]
  C0_tail <- P0_tail - D_tail
  decay_tail <- exp(-lambda * m_tail)
  C_tail <- C0_tail * decay_tail
  P_tail <- C_tail + D_tail
  P_long <- c(P[1:(last_camp-1)], P_tail)
  C_long <- P_long - D_long
  
  # P0 and D over projection window
  P0_proj <- tail(P0_long, n = N_proj)
  D_proj <- tail(D_long, n = N_proj)
  C0_proj <- P0_proj - D_proj
  
  # Back-fill first value for early usage
  N_early <- CMC_first - CMC_Jan2000
  D_early <- rep(D[1], N_early)
  P_early <- rep(P[1], N_early)
  C_early <- P_early - D_early
  
  # Calculate routine usage
  D_full <- c(D_early, D_long[1:(proj_camp_1-1)], D_proj)
  
  # Calculate campaign usage
  decay_proj <- exp(-lambda * m_proj)
  C_proj <- C0_proj * decay_proj
  C_pre_proj <- c(C_early, C_long[1:(proj_camp_1-1)])
  C_full <- c(C_early, C_long[1:(proj_camp_1-1)], C_proj)
  
  # Calculate overall usage
  P_full <- C_full + D_full
  
  # Usage with no future mass campaigns
  P_D_proj_only <- c(P_early,P_long)
  
  # Times for fitting
  times_mn <- seq(1, proj_end)
  times_yr <- rep(seq(0, ceiling(N_CMC_sim / 12)), each=12)
  times_1st_dy <- DOY_1st + (times_yr * year)
  times_mid_dy <- DOY_mid + (times_yr * year)
  input_net_times <- times_mid_dy[1:N_CMC_sim]    # usage for fitting
  output_net_times <- times_1st_dy[1:N_CMC_sim]   # distribution times for netz
  
  # Fit nets with no future mass campaigns
  output_nets_no_future_mdc <- fit_usage_sequential(
    target_usage = P_D_proj_only,
    target_usage_timesteps = input_net_times,
    distribution_timesteps = output_net_times,
    mean_retention = mean_retention_dy
  )
  
  saveRDS(P_full, "P_test.rds")
  
  
  # Fit nets with future mass campaigns
  output_nets_future_mdc <- fit_usage_sequential(
    target_usage = P_full,
    target_usage_timesteps = input_net_times,
    distribution_timesteps = output_net_times,
    mean_retention = mean_retention_dy
  )
  

  
  # Future nets distributed through mass campaigns
  future_mdc_nets_only <- output_nets_future_mdc - output_nets_no_future_mdc
  
  # Interval costings
  
  
  # interval adjustment
  # if (interval_reduction) {
  #   future_mdc_nets_only <- future_mdc_nets_only * mass_int_mn / 36
  # }
  
  # Routine baseline
  if (routine_baseline) {
    all_output_nets <- output_nets_no_future_mdc
  } else {
    all_output_nets <- output_nets_no_future_mdc + future_mdc_nets_only
  }
  
  # New net range (month ids)
  new_net_range <- seq(N_CMC_sim - N_proj + 1, N_CMC_sim)
  
  # # net type costing
  # if (net_costings) {
  #   all_output_nets[new_net_range] <- all_output_nets[new_net_range] * cost_factor
  # }
  
  # Combine campaign and routine nets
  # all_output_nets <- output_nets_no_future_mdc + future_mdc_nets_only
  
  # no future nets
  if (no_future_nets) {
    all_output_nets[new_net_range] <- rep(0, length(new_net_range))
    net_strategy <- "no future nets"
  } else {
    #output_new_nets <- all_output_nets[new_net_range]
    #new_net_tail <- tail(output_new_nets, mass_int_mn)
    saveRDS(new_net_range, "test1")
    test2 <- rep(tail(all_output_nets,
                      mass_int_mn),
               length.out = length(new_net_range))
    saveRDS(test2, "test2")
    
    #Repeat final distribution tail
    netend <- tail(all_output_nets, mass_int_mn)
    ordered_netend <- head(c(netend[which(netend==max(netend)):length(netend)],
                             netend[1:which(netend==max(netend))]),
                           length(netend))
    
    all_output_nets[new_net_range] <- rep(ordered_netend,
                                          length.out = length(new_net_range))
  }
  
  # tail nets
  avg_tail_nets <- sum(tail(all_output_nets * tail_pop, n = 12 * 12)) / 12
  avg_pop_adj_tail_nets <- avg_tail_nets * npc_costed #/ 1.8
  tail_net_cost <- avg_pop_adj_tail_nets * new_net_cost
  
  # cost override
  # if (override_cost) {
  #   if (override_mdc_only) {
  #     avg_tail_nets_routine <- sum(tail(output_nets_no_future_mdc * tail_pop, n = 6 * 12)) / 6
  #     avg_pop_adj_tail_nets_routine <- avg_tail_nets_routine / 1.8
  #     tail_net_cost_routine <- avg_pop_adj_tail_nets_routine * new_net_cost
  #     tail_mdc_cost <- tail_net_cost - tail_net_cost_routine
  #     mdc_budget <- max(0, override_cost_value - tail_net_cost_routine)
  #     override_cost_factor <- mdc_budget / tail_mdc_cost
  #     future_mdc_nets_only[new_net_range] <- future_mdc_nets_only[new_net_range] * override_cost_factor
  #     all_output_nets <- output_nets_no_future_mdc + future_mdc_nets_only
  #   } else {
  #     override_cost_factor <- override_cost_value / tail_net_cost
  #     all_output_nets[new_net_range] <- all_output_nets[new_net_range] * override_cost_factor
  #   }
  # }
  
  # Wider future distributions
  # monthly_routine <- fit_usage_sequential(target_usage=rep(D_end, 2),
  #                                         target_usage_timesteps=seq(15,45,30),
  #                                         distribution_timesteps=seq(1,31,30),
  #                                         mean_retention=mean_retention_dy)[2]
  # proj_b <- N_early + proj_end
  # proj_a <- proj_b - N_proj
  # proj_sd <- 365/12
  # smooth_output_nets <- all_output_nets
  # for (i in proj_a:proj_b) {
  #   camp_nets <- all_output_nets[i] - monthly_routine
  #   smooth_output_nets[i] <- smooth_output_nets[i] - camp_nets
  #   camp_weights <- dnorm(output_net_times, output_net_times[i], proj_sd)
  #   camp_sum <- sum(camp_weights)
  #   camp_norm_weights <- camp_weights / camp_sum
  #   smooth_camp_nets <- camp_nets * camp_norm_weights
  #   running_use <- smooth_output_nets[1]
  #   smooth_output_nets[1] <- smooth_output_nets[1]+smooth_camp_nets[1]/(1-running_use)
  #   for (j in 2:proj_b) {
  #     running_use <- running_use * exp(-(output_net_times[j]-output_net_times[j-1])/mean_retention_dy)
  #     smooth_output_nets[j] <- smooth_output_nets[j]+smooth_camp_nets[j]/(1-running_use)
  #     running_use <- running_use + smooth_camp_nets[j]
  #   }
  # }
  
  # set bednets
  bednet_pars <- malariasimulation::set_bednets(site_pars,
                                                timesteps = output_net_times,
                                                coverages = all_output_nets,
                                                retention = mean_retention_dy,
                                                dn0 = dn0_mat,
                                                rn = rn_mat,
                                                rnm = rnm_mat,
                                                gamman = gam_vec)
  
  
  # run simulation
  output <- malariasimulation::run_simulation(timesteps = bednet_pars$timesteps,
                                              parameters = bednet_pars)
  
  N_timesteps <- bednet_pars$timesteps
  
  # 
  # # return avg annual infections over observation window
  # 
  # output_df <- 
  
  # collate model outputs
  S_count <- output$S_count
  A_count <- output$A_count
  D_count <- output$D_count
  U_count <- output$U_count
  Tr_count <- output$Tr_count
  n_total <- S_count + A_count + D_count + U_count + Tr_count
  timestep_yr <- bednet_pars$baseline_year + (output$timestep - 1) / 365
  pfin_all_ages <- output$n_infections / n_total
  pfpr_730_3649 <- output$n_detect_730_3649 / output$n_730_3649
  
  obs_start <- N_timesteps - obs_window
  obs_infections <- sum(output$n_infections[obs_start:N_timesteps])
  annual_infections <- 365 * obs_infections / obs_window
  pred_ann_infect <- tail_pop * annual_infections / sim_population
  
  avg_pfpr <- sum(pfpr_730_3649[obs_start:N_timesteps]) / obs_window
  
  area_net_strategy <- paste(fs_area, net_strategy, sep = " ")
  
  prop_use <- output$n_use_net / n_total
  projection_final_yr <- ceiling(projection_window_mn / 12)
  yr_end <- N_timesteps
  yr_srt <- yr_end - 365 + 1
  for (i in projection_final_yr:1) {
    x_use <- mean(prop_use[yr_srt:yr_end])
    yrs_post_camp <- (i-1) %% (round(mass_int_mn/12)) + 1
    
    if (yrs_post_camp == 1) {
      if(exists("yr1_use")) {yr1_use <- c(yr1_use, x_use)} else {yr1_use <- x_use}
    } else if (yrs_post_camp == 2) {
      if(exists("yr2_use")) {yr2_use <- c(yr2_use, x_use)} else {yr2_use <- x_use}
    } else if (yrs_post_camp == 3) {
      if(exists("yr3_use")) {yr3_use <- c(yr3_use, x_use)} else {yr3_use <- x_use}
    }
    yr_end <- yr_end - 365
    yr_srt <- yr_srt - 365
  }
  avg_yr1_use <- NA
  avg_yr2_use <- NA
  avg_yr3_use <- NA
  if (mass_int_mn >= 12) {
    avg_yr1_use <- mean(yr1_use, na.rm = TRUE)
    if (mass_int_mn >= 24) {
      avg_yr2_use <- mean(yr2_use, na.rm = TRUE)
      if (mass_int_mn >= 36) {
        avg_yr3_use <- mean(yr3_use, na.rm = TRUE)
      }
    }
  }
  
  output_summary_df <- data.frame("fs_area_id" = fs_area_id,
                                  "fs_area" = fs_area,
                                  "ISO2" = ISO2,
                                  "fs_name_1" = fs_name_1,
                                  "urbanicity" = urbanicity,
                                  "pop" = tail_pop,
                                  "net_strategy" = net_strategy,
                                  "net_name" = net_name,
                                  "mass_int" = mass_int_mn/12,
                                  "net_costings" = net_cost_logical,
                                  "interval_costings" = interval_reduction_logical,
                                  "routine_baseline" = routine_baseline_logical,
                                  "sample_index" = sid,
                                  "area_net_strategy" = area_net_strategy,
                                  "annual_infections" = annual_infections,
                                  "pred_ann_infect" = pred_ann_infect,
                                  "avg_pfpr" = avg_pfpr,
                                  "avg_ann_nets_distrib" = avg_tail_nets,
                                  "avg_ann_net_cost" = tail_net_cost,
                                  "avg_yr1_use" = avg_yr1_use,
                                  "avg_yr2_use" = avg_yr2_use,
                                  "avg_yr3_use" = avg_yr3_use
  )
  
  if (output_sim) {
    output_summary_rep_df <- do.call("rbind", replicate(dim(output)[1],
                                                        output_summary_df,
                                                        simplify = FALSE))
    output_df <- cbind.data.frame(output, output_summary_rep_df)
  } else {
    output_df <- output_summary_df
  }
  
  
  # area_net_strategy <- paste(fs_area, net_strategy, sep = " ")
  # 
  # timestamp <- as.numeric(Sys.time())*100000
  # pcname <- Sys.info()[[4]]
  # 
  # if (net_name == "pyrethroid-PBO") {
  #   csvpath <- "./outputs/malsim0_pbo_bicost/"
  # } else if (net_name == "pyrethroid-pyrrole") {
  #   csvpath <- "./outputs/malsim0_pyrrole_bicost/"
  # } else {
  #   csvpath <- "./outputs/malsim0_only_bicost/"
  # }
  # 
  # csvname <- paste0(fs_area_undrscr, "_",
  #                   net_strategy, "_",
  #                   pcname, "_",
  #                   timestamp, ".csv")
  # csvpathname <- paste0(csvpath, csvname)
  # write.table(output_df,
  #             file = csvpathname,
  #             sep = ",",
  #             col.names = TRUE,
  #             row.names = FALSE)
  
  return(output_df)
  
}

run_malsim_nets_sequential_costed <- function(dataset,
                                           areas_included = NULL,
                                           N_reps = 100,
                                           N_cores = 0,
                                           mass_int_yr = c(2,3),
                                           sim_population = 5e4,
                                           ref_CMC = 1476,
                                           only = FALSE,
                                           pbo = FALSE,
                                           pyrrole = FALSE,
                                           net_costings = FALSE,
                                           interval_reduction = FALSE,
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
                                           bv_access = NULL,
                                           bv_npc = NULL,
                                           output_sim = FALSE,
                                           camp_offset_from_peak = TRUE,
                                           month_offset_from_peak = -4,
                                           future_camps_month_sd = 1) {
  
  # Set number of cores
  if (N_cores <= 0) {N_cores <- length(areas_included)}
  if (use_hipercow) {max_cores <- 32} else {max_cores <- 20}
  if (N_cores > max_cores) {N_cores <- max_cores}
  
  # Simulation time
  CMC_sim_start <- CMC_Jan2000
  CMC_sim_end <- CMC_last + 12 * 12 #projection_window_mn
  N_CMC_sim <- CMC_sim_end - CMC_sim_start + 1
  N_old_nets <- CMC_last - CMC_sim_start + 1
  N_backfill <- CMC_first - CMC_sim_start
  
  # Number of samples
  N_samples <- dim(P_u)[1]
  
  # Create sample ids
  if (max(long_sample_ids) > N_samples) {
    print("Warning: Some sample ids outwith range")
  }
  rep_ids <- seq(1, N_reps) + rep_offset
  sample_ids <- long_sample_ids[1:N_reps]
  
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
    
    if ((l==1 & only & !net_costings) | (l==2 & pbo) | (l==3 & pyrrole)) {
      
      #if (net_costings) {
        if (l==1 & only) {cost_factor <- 1.0}
        if (l==2 & pbo) {cost_factor <- scaled_pbo_nets_equiv_only}
        if (l==3 & pyrrole) {cost_factor <- scaled_pyrrole_nets_equiv_only}
      #} else {
        # cost_factor <- 1.0
      #}
      
      if (l==1 & only) {new_net_cost <- only_total_cost}
      if (l==2 & pbo) {new_net_cost <- pbo_total_cost}
      if (l==3 & pyrrole) {new_net_cost <- pyrrole_total_cost}
      
      for (k in 1:N_int_vals) {
        
        # mass interval
        mass_int_mn <- mass_int_yr[k] * 12
        mass_int_dy <- mass_int_yr[k] * 365
        
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
            area_time_ref_id <- which(dataset$area_id == area_id &
                                        dataset$CMC == ref_CMC)
            area_time_ids <- which(dataset$area_id == area_id)
            
            # get samples
            invlam_samples <- invlam_u[sample_ids, area_id] %>%
              as.vector %>% unname %>% unlist
            lam_samples <- 1 / invlam_samples
            ret_ref_samples <- ret_u[sample_ids, area_time_ref_id] %>%
              as.vector %>% unname %>% unlist
            P_samples <- P_u[sample_ids, area_time_ids] %>%
              as.matrix %>% unname
            P0_samples <- P0_u[sample_ids, area_time_ids] %>%
              as.matrix %>% unname
            D_samples <- D_u[sample_ids, area_time_ids] %>%
              as.matrix %>% unname
            
            # get final usage and access samples (set to ref_CMC: 1476 = DEC 22)
            invlam_a_samples <- invlam_a[sample_ids, area_id] %>%
              as.vector %>% unname %>% unlist
            lam_a_samples <- 1 / invlam_a_samples
            P0_ref_samples <- P0_u[sample_ids, area_time_ref_id] %>%
              as.vector %>% unname %>% unlist
            D_ref_samples <- D_u[sample_ids, area_time_ref_id] %>%
              as.vector %>% unname %>% unlist
            C0_ref_samples <- P0_ref_samples - D_ref_samples
            P0_a_ref_samples <- P0_a[sample_ids, area_time_ref_id] %>%
              as.vector %>% unname %>% unlist
            D_a_ref_samples <- D_a[sample_ids, area_time_ref_id] %>%
              as.vector %>% unname %>% unlist
            C0_a_ref_samples <- P0_a_ref_samples - D_a_ref_samples
            
            # avg ref usage and access
            ref_usage_D <- D_ref_samples
            lam_del <- lam_samples * 36
            ref_usage_C <- 1 - exp(-lam_del)
            ref_usage_C <- ref_usage_C * C0_ref_samples / lam_del
            ref_usage <- ref_usage_D + ref_usage_C
            ref_access_D <- D_a_ref_samples
            lam_del_a <- lam_a_samples * 36
            ref_access_C <- 1 - exp(-lam_del_a)
            ref_access_C <- ref_access_C * C0_a_ref_samples / lam_del_a
            ref_access <- ref_access_D + ref_access_C
            
            # Identify month with max predicted usage (estimated last mass campaign)
            # last_camp_month <- apply(P_samples, 1, which.max)
            last_camp_mn <- max(which(net_data$months_post_mdc[area_time_ids] == 0))
            
            # First regular MDC with new nets
            #first_new_net_month <- last_camp_month + mass_int_mn + 
            
            # Generate ISO code for current admin
            admin_country <- countrycode(fs_id_link$ISO2[i], "iso2c", "iso3c")
            
            # Pull country site
            if (admin_country != current_country) {
              current_country <- admin_country
              ctry_site <- get_site(current_country)
            }
            
            # Isolate a single site from a country
            adm_site_index <- which(ctry_site$sites$name_1 == fs_id_link$fs_name_1[i] &
                                      ctry_site$sites$urban_rural == fs_id_link$urbanicity[i])
            
            # If no foresite file for urban/rural, then revert to other
            if (identical(adm_site_index, integer(0))) {
              if (fs_id_link$urbanicity[i] == "urban") {
                adm_site_index <- which(ctry_site$sites$name_1 == fs_id_link$fs_name_1[i] &
                                          ctry_site$sites$urban_rural == "rural")
              } else {
                adm_site_index <- which(ctry_site$sites$name_1 == fs_id_link$fs_name_1[i] &
                                          ctry_site$sites$urban_rural == "urban")
              }
            }
            
            # Check for successful foresite match
            if (identical(adm_site_index, integer(0))) {
              print(paste0("Warning: foresite not linked for admin region ",
                           fs_id_link$fs_name_1[i], " (index ", i, ") in ",
                           current_country))
            }
            
            # Create admin site file
            adm_site <- site::single_site(ctry_site, adm_site_index)
            
            # Repeat interventions
            adm_site %<>% expand_interventions(expand_year = 12,
                                               delay = 0,
                                               counterfactual = FALSE)
            
            # Tail population
            tail_ids <- which(adm_site$population$year >= 2023 &
                                adm_site$population$year <= 2028)
            tail_pop <- mean(adm_site$population$pop[tail_ids])
            
            # Pyrethroid resistance
            yearly_res <- adm_site$pyrethroid_resistance$pyrethroid_resistance
            monthly_res <- rep(yearly_res, each = 12)
            round_monthly_res <- round(monthly_res, 2)
            
            # if (l==1) {old_res <- res_only}
            # if (l==2) {old_res <- res_pbo}
            # if (l==3) {old_res <- res_pyrrole}
            
            # Old pyrethroid-only nets
            old_res <- res_only
            
            # New net types
            if (l==1) {new_res <- res_only}
            if (l==2) {new_res <- res_pbo}
            if (l==3) {new_res <- res_pyrrole}
            
            if (l==1) {net_name <- "pyrethroid-only"}
            if (l==2) {net_name <- "pyrethroid-PBO"}
            if (l==3) {net_name <- "pyrethroid-pyrrole"}
            
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
              if (interval_reduction & mass_int_yr[k] == 2) {
                net_strategy %<>% paste("interval and type costed")
              } else {
                net_strategy %<>% paste("type costed")
              }
            } else {
              if (interval_reduction & mass_int_yr[k] == 2) {
                net_strategy %<>% paste("interval costed")
              }
            }
            
            # Identify resistance id matches (same for old and new net types)
            res_ids <- match(round_monthly_res, old_res$resistance)
            N_species <- length(adm_site$vectors$species)
            
            # Create dn0 matrix
            dn0_old <- old_res$dn0_med[res_ids]
            dn0_vec <- new_res$dn0_med[res_ids]   # initially set for new net
            dn0_vec[1:N_old_nets] <- dn0_old[1:N_old_nets]  # earlier dates to old net
            dn0_vec <- dn0_vec[1:N_CMC_sim]       # restrict to sim length
            dn0_mat <- matrix(rep(dn0_vec, N_species),
                              nrow = N_CMC_sim,
                              ncol = N_species)
            
            
            # Create dn0 matrix
            rn_old <- old_res$rn0_med[res_ids]
            rn_vec <- new_res$rn0_med[res_ids]   # initially set for new net
            rn_vec[1:N_old_nets] <- rn_old[1:N_old_nets]  # earlier dates to old net
            rn_vec <- rn_vec[1:N_CMC_sim]       # restrict to sim length
            rn_mat <- matrix(rep(rn_vec, N_species),
                             nrow = N_CMC_sim,
                             ncol = N_species)
            
            # rn_vec <- old_res$rn0_med[res_ids]
            # rn_vec <- rn_vec[1:N_CMC_sim]
            # rn_mat <<- matrix(rep(rn_vec, N_species),
            #                   nrow = N_CMC_sim,
            #                   ncol = N_species)
            
            # Create rnm matrix
            rnm_mat <- matrix(rep(0.24, N_CMC_sim * N_species),
                              nrow = N_CMC_sim,
                              ncol = N_species)
            
            # Create gamman vector
            gam_old <- 365 * old_res$gamman_med[res_ids] / log(2)
            gam_vec <- 365 * new_res$gamman_med[res_ids] / log(2)
            gam_vec[1:N_old_nets] <- gam_old[1:N_old_nets]
            gam_vec <- gam_vec[1:N_CMC_sim]
            
            # Manually change Kedougou EIR = 250 from:
            #Dia, I., et al., Bionomics of Anopheles gambiae Giles, An. arabiensis Patton, An. funestus Giles and An.
            #nili (Theobald) (Diptera: Culicidae) and transmission of Plasmodium falciparum in a Sudano-Guinean
            #zone (Ngari, Senegal). Journal of medical entomology, 2003. 40(3).
            if (adm_site$eir$name_1[1] == "Kédougou") {
              adm_site$eir$eir[1] <- 250
            }
            
            # Pf EIR
            Pf_eir <- adm_site$eir$eir[1]
            
            # Update bednet interventions
            
            if (Pf_eir > 0) {
              # Create parameter inputs
              site_pars <- site::site_parameters(
                interventions = adm_site$interventions,
                demography = adm_site$demography,
                vectors = adm_site$vectors,
                seasonality = adm_site$seasonality,
                eir = Pf_eir,
                overrides = list(human_population = sim_population,
                                 individual_mosquitoes = FALSE)
              )
              
              # peak season
              peak_season_dy <- peak_season_offset(site_pars)
              
              # Combine vector and matrix parameters for parLapply function
              site_pars$P_samples <- P_samples
              site_pars$P0_samples <- P0_samples
              site_pars$D_samples <- D_samples
              site_pars$lam_samples <- lam_samples
              site_pars$dn0_mat <- dn0_mat
              site_pars$rn_mat <- rn_mat
              site_pars$rnm_mat <- rnm_mat
              site_pars$gam_vec <- gam_vec
              site_pars$DOY_1st <- DOY_1st
              site_pars$DOY_mid <- DOY_mid
              site_pars$ref_usage_D <- ref_usage_D
              site_pars$ref_usage_C <- ref_usage_C
              site_pars$ref_usage <- ref_usage
              site_pars$lam_a_samples <- lam_a_samples
              site_pars$ref_access_D <- ref_access_D
              site_pars$ref_access_C <- ref_access_C
              site_pars$ref_access <- ref_access
              
              # Combine single parameters
              site_pars$net_type <- l
              site_pars$net_name <- net_name
              site_pars$net_strategy <- net_strategy
              site_pars$last_camp <- last_camp_mn#[j]
              site_pars$mass_int_mn <- mass_int_mn
              site_pars$ISO2 <- fs_id_link$ISO2[i]
              site_pars$fs_area <- fs_id_link$fs_area[i]
              site_pars$fs_name_1 <- fs_id_link$fs_name_1[i]
              site_pars$urbanicity <- fs_id_link$urbanicity[i]
              site_pars$fs_area_id <- fs_id_link$fs_area_id[i]
              site_pars$N_species <- N_species
              site_pars$CMC_first <- CMC_first
              site_pars$CMC_Jan2000 <- CMC_Jan2000
              site_pars$projection_window_mn <- 12 * 12 #projection_window_mn
              site_pars$N_CMC <- N_CMC
              site_pars$N_CMC_sim <- N_CMC_sim
              site_pars$tail_pop <- tail_pop
              site_pars$sim_population <- sim_population
              site_pars$net_costings <- net_costings
              site_pars$cost_factor <- cost_factor
              site_pars$interval_reduction <- interval_reduction
              site_pars$routine_baseline <- routine_baseline
              site_pars$new_net_cost <- new_net_cost
              site_pars$no_future_nets <- no_future_nets
              site_pars$override_cost <- override_cost
              site_pars$override_mdc_only <- override_mdc_only
              site_pars$override_cost_value <- override_cost_value
              site_pars$output_sim <- output_sim
              site_pars$N_backfill <- N_backfill
              
              # Combine BV access and npc
              site_pars$bv_access <- bv_access
              site_pars$bv_npc <- bv_npc
              
              #Offset from peak
              site_pars$camp_offset_from_peak <- camp_offset_from_peak
              # site_pars$month_offset_from_peak <- month_offset_from_peak
              
              # site_pars$peak_season_dy <- peak_season_dy
              
              # Set mean day offset
              if (camp_offset_from_peak) {
                mean_day_offset <- peak_season_dy + month_offset_from_peak * 365 / 12
                mean_day_offset <- mean_day_offset %% 365
                last_camp_dy <- (last_camp_mn - 0.5) * 365 / 12
                last_camp_doy <- last_camp_dy %% 365
                last_camp_prev_yrs <- last_camp_dy %/% 365
                proj_camp_override_dy <- 365 * last_camp_prev_yrs + mass_int_dy + last_camp_dy
                proj_camp_override_dy_int <- proj_camp_override_dy - last_camp_dy
                # Adjust mean offset to account for change of year
                if (proj_camp_override_dy_int < (mass_int_dy - 182.5)) {
                  proj_camp_override_dy <- proj_camp_override_dy + 365
                  proj_camp_override_dy_int <- proj_camp_override_dy_int + 365
                } else if (proj_camp_override_dy_int > (mass_int_dy + 182.5)) {
                  proj_camp_override_dy <- proj_camp_override_dy - 365
                  proj_camp_override_dy_int <- proj_camp_override_dy_int - 365
                }
              }
              
              for (j in 1:N_reps) {
                
                jj <- j + rep_offset
                site_pars$sample_index <- jj
                site_pars$mean_ret <- ret_ref_samples[jj]
                
                if (camp_offset_from_peak) {
                  rand_camp_offset <- future_camps_month_sd * rnormvals[jj] * 365 / 12
                  proj_camp_override_dy <- proj_camp_override_dy + rand_camp_offset
                  proj_camp_override_month <- round(proj_camp_override_dy * 12 / 365)
                  site_pars$proj_camp_override_month <- proj_camp_override_month
                  site_pars$month_offset <- 0
                } else {
                  site_pars$proj_camp_override_month <- 0
                  site_pars$month_offset <- long_month_offset[jj]
                }
                
                param_list[[length(param_list) + 1]] <- site_pars
                # 
                # if (use_hipercow) {
                #   dynam_id <- paste("id", i, j, k, l, ii, jj, sep = "_")
                #   hipercow_params <- param_list[[jj]]
                #   if (debugging) {
                #     assign(dynam_id,
                #            task_create_expr(par_net_region_sequential3(hipercow_params)),
                #            envir = .GlobalEnv)
                #   } else {
                #     assign(dynam_id,
                #            task_create_expr(par_net_region_sequential3(hipercow_params)))
                #   }
                # }
                
                
              }
              
            }
            
            
            pc1 <- round(100 * ii / N_total_its)
            if (pc1 > pc0) {
              pc0 <- pc1
              print(paste(pc0, "% complete", sep = ""))
            }
            
            #print(paste("MDC val ", k, " for region ", i, " of ", N_areas, " complete", sep = ""))
          }
        }
        
      }
    }
  }
  
  
  if (use_hipercow) {
    resources <- hipercow_resources(cores = N_cores)
    if (N_total_its == 1) {
      par_id <- task_create_expr(par_net_region_sequential_costed(param_list[[1]]),
                                 resources = resources)
    } else {
      par_id <- task_create_expr(
        parallel::parLapply(NULL, param_list, par_net_region_sequential_costed),
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
      "fit_usage_sequential",
      "population_usage_t"))
    #par_output <- lapply(param_list, par_net_region_sequential3)
    par_output <- parLapply(cl, param_list, par_net_region_sequential_costed)
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