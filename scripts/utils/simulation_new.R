# simulation_repeat.R

par_net_region_sequential_repeat <- function(param_list) {
  
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
  net_cost_strategy_id <- site_pars$net_cost_strategy_id
  cost_strategy <- site_pars$cost_strategy
  cost_factor <- site_pars$cost_factor
  biennial_reduction <- site_pars$biennial_reduction
  
  if (biennial_reduction & (mass_int_mn < 25)) {
    net_strategy <- paste0(net_strategy, "_bien_costed")
  } else if ((cost_factor < 0.9999) | (cost_factor > 1.0001)) {
    net_strategy <- paste0(net_strategy, "_costed")
  }
  
  fs_area_undrscr <- gsub(" ", "_", fs_area)
  
  year = 365
  obs_window = 6 * year
  
  # convert retention to days
  mean_retention_dy <- 365 * mean_retention / 12
  
  # Central time point for first regular mass campaign
  proj_camp_1 <- last_camp + mass_int_mn + month_offset
  
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
  
  # Extend values
  P0_end <- tail(P0, n = 1)
  D_end <- tail(D, n = 1)
  P0_long <- c(P0, rep(P0_end, projection_window_mn))
  D_long <- c(D, rep(D_end, projection_window_mn))
  m_long <- seq(1, proj_end) - last_camp
  m_tail <- m_long[last_camp:proj_end]
  P0_tail <- P0_long[last_camp:proj_end]
  D_tail <- D_long[last_camp:proj_end]
  decay_tail <- exp(-lambda * m_tail)
  P_tail <- P0_tail * decay_tail + (1 - decay_tail) * D_tail
  P_long <- c(P[1:(last_camp-1)], P_tail)
  
  # P0 and D over projection window
  P0_proj <- tail(P0_long, n = N_proj)
  D_proj <- tail(D_long, n = N_proj)
  
  # Generate projected usage
  decay_proj <- exp(-lambda * m_proj)
  P_proj <- P0_proj * decay_proj + (1 - decay_proj) * D_proj
  
  # Generate whole time series of usage
  # Refinements to P_early
  P_early <- rep(P[1], CMC_first - CMC_Jan2000)
  input_net_usage <- c(P_early, P_long[1:(proj_camp_1-1)], P_proj)
  N_input <- length(input_net_usage)
  
  times_mn <- seq(1, proj_end)
  times_yr <- rep(seq(0, ceiling(N_input / 12)), each=12)
  times_1st_dy <- DOY_1st + (times_yr * year)
  times_mid_dy <- DOY_mid + (times_yr * year)
  
  input_net_times <- times_mid_dy[1:N_input]    # usage for fitting
  output_net_times <- times_1st_dy[1:N_input]   # distribution times for netz
  
  # netz fit
  output_nets_distrib <- fit_usage_sequential(target_usage = input_net_usage,
                                              target_usage_timesteps = input_net_times,
                                              distribution_timesteps = output_net_times,
                                              mean_retention = mean_retention_dy)
  output_nets_distrib <- output_nets_distrib * cost_factor
  
  # biennial adjustment
  if (biennial_reduction & (mass_int_mn < 25)) {
    prop_camp_proj <- mean((P_proj - D_proj) / P_proj)
    bien_factor <- (2.0/3.0) / prop_camp_proj
    output_nets_distrib[last_camp:proj_end] <- output_nets_distrib * bien_factor
  }
  
  # tail nets
  avg_tail_nets <- sum(tail(output_nets_distrib * tail_pop, n = 6 * 12)) / 6
  
  # set bednets
  bednet_pars <- malariasimulation::set_bednets(site_pars,
                                                timesteps = output_net_times,
                                                coverages = output_nets_distrib,
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
  annual_infections <- obs_infections / obs_window
  pred_ann_infect <- tail_pop * annual_infections / sim_population
  
  avg_pfpr <- sum(pfpr_730_3649[obs_start:N_timesteps]) / obs_window
  
  area_net_strategy <- paste(fs_area, net_strategy, sep = " ")
  
  output_df <- data.frame("fs_area_id" = fs_area_id,
                          "ISO2" = ISO2,
                          "fs_area" = fs_area,
                          "fs_name_1" = fs_name_1,
                          "urbanicity" = urbanicity,
                          "pop" = tail_pop,
                          "net_strategy" = net_strategy,
                          "net_name" = net_name,
                          "mass_int" = mass_int_mn/12,
                          "sample_index" = sid,
                          "area_sample_id" = 1000 * fs_area_id + sid,
                          "area_net_strategy" = area_net_strategy,
                          "annual_infections" = annual_infections,
                          "pred_ann_infect" = pred_ann_infect,
                          "avg_pfpr" = avg_pfpr,
                          "avg_ann_nets_distrib" = avg_tail_nets
  )

  
  area_net_strategy <- paste(fs_area, net_strategy, sep = " ")
  
  timestamp <- as.numeric(Sys.time())*100000
  pcname <- Sys.info()[[4]]
  
  if (net_name == "pyrethroid-PBO") {
    csvpath <- "./outputs/malsim0_pbo_bicost/"
  } else if (net_name == "pyrethroid-pyrrole") {
    csvpath <- "./outputs/malsim0_pyrrole_bicost/"
  } else {
    csvpath <- "./outputs/malsim0_only_bicost/"
  }
  
  csvname <- paste0(fs_area_undrscr, "_",
                    net_strategy, "_",
                    pcname, "_",
                    timestamp, ".csv")
  csvpathname <- paste0(csvpath, csvname)
  write.table(output_df,
              file = csvpathname,
              sep = ",",
              col.names = TRUE,
              row.names = FALSE)
  
  return(output_df)
  
}

run_malsim_nets_sequential_new <- function(dataset,
                                        areas_included = NULL,
                                        N_reps = 100,
                                        N_cores = 8,
                                        mass_int_yr = c(2,3),
                                        ref_CMC = 1476,
                                        only = TRUE,
                                        pbo = TRUE,
                                        pyrrole = TRUE,
                                        net_costings = TRUE,
                                        biennial_reduction = FALSE,
                                        month_default_offset = 0,
                                        use_hipercow = FALSE,
                                        debugging = FALSE) {
  
  # Simulation time
  CMC_sim_start <- CMC_Jan2000
  CMC_sim_end <- CMC_last + projection_window_mn
  N_CMC_sim <- CMC_sim_end - CMC_sim_start + 1
  
  # Number of samples
  N_samples <- dim(P_u)[1]
  
  # Create sample ids
  if (max(long_sample_ids) > N_samples) {
    print("Warning: Some sample ids outwith range")
  }
  sample_id <- long_sample_ids[1:N_reps]
  
  # dataframe for storing output
  output_df <- data.frame(NULL)
  
  # progress indicator
  N_net_types <- only + pbo + pyrrole
  N_int_vals <- length(mass_int_yr)
  N_areas_included <- length(areas_included)
  N_total_its <- N_net_types * N_int_vals * N_areas_included
  pc0 <- 0
  ii <- 0
  jj <- 1
  
  for (l in 1:3) {
    
    if ((l==1 & only & !net_costings) | (l==2 & pbo) | (l==3 & pyrrole)) {
      
      if (net_costings) {
        if (pbo) {cost_factor <- scaled_pbo_nets_equiv_only}
        if (pyrrole) {cost_factor <- scaled_pbo_nets_equiv_only}
      } else {
        cost_factor <- 1.0
      }
      
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
            area_time_ref_id <- which(dataset$area_id == area_id &
                                        dataset$CMC == ref_CMC)
            area_time_ids <- which(dataset$area_id == area_id)
            
            # get samples
            invlam_samples <- invlam_u[sample_id, area_id] %>%
              as.vector %>% unname %>% unlist
            lam_samples <- 1 / invlam_samples
            ret_ref_samples <- ret_u[sample_id, area_time_ref_id] %>%
              as.vector %>% unname %>% unlist
            P_samples <- P_u[sample_id, area_time_ids] %>%
              as.matrix %>% unname
            P0_samples <- P0_u[sample_id, area_time_ids] %>%
              as.matrix %>% unname
            D_samples <- D_u[sample_id, area_time_ids] %>%
              as.matrix %>% unname
            
            # Identify month with max predicted usage (estimated last mass campaign)
            # last_camp_month <- apply(P_samples, 1, which.max)
            last_camp_month <- max(apply(P_samples, 1, which.max))
            
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
            adm_site %<>% expand_interventions(expand_year = 6,
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
            
            # if (l==1) {pyr_res <- res_only}
            # if (l==2) {pyr_res <- res_pbo}
            # if (l==3) {pyr_res <- res_pyrrole}
            
            old_res <- res_only
            
            if (l==1) {new_res <- res_only}
            if (l==2) {new_res <- res_pbo}
            if (l==3) {new_res <- res_pyrrole}
            
            if (l==1) {net_name <- "pyrethroid-only"}
            if (l==2) {net_name <- "pyrethroid-PBO"}
            if (l==3) {net_name <- "pyrethroid-pyrrole"}
            
            # Name net strategy
            net_strategy <- paste(net_name, mass_int_yr[k], "year interval")
            if (net_costings) {
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
            
            
            res_ids <- match(round_monthly_res, pyr_res$resistance)
            N_species <- length(adm_site$vectors$species)
            
            dn0_vec <- pyr_res$dn0_med[res_ids]
            dn0_vec <- dn0_vec[1:N_CMC_sim]
            dn0_mat <<- matrix(rep(dn0_vec, N_species),
                               nrow = N_CMC_sim,
                               ncol = N_species)
            
            rn_vec <- pyr_res$rn0_med[res_ids]
            rn_vec <- rn_vec[1:N_CMC_sim]
            rn_mat <<- matrix(rep(rn_vec, N_species),
                              nrow = N_CMC_sim,
                              ncol = N_species)
            
            rnm_mat <<- matrix(rep(0.24, N_CMC_sim * N_species),
                               nrow = N_CMC_sim,
                               ncol = N_species)
            
            gam_vec <<- 365 * pyr_res$gamman_med[res_ids] / log(2)
            gam_vec <<- gam_vec[1:N_CMC_sim]
            
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
              
              # Combine parameters for parLapply function
              
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
              
              param_list <- list()
              for (j in 1:N_reps) {
                param_list[[jj]] <- c(site_pars,
                                      "sample_index" = j,
                                      "mean_ret" = ret_ref_samples[j],
                                      "net_type" = l,
                                      "net_name" = net_name,
                                      "net_strategy" = net_strategy,
                                      "month_offset" = long_month_offset[j+month_default_offset],
                                      "last_camp" = last_camp_month,#[j],
                                      #"top_up_int" = top_up_int,
                                      "mass_int_mn" = mass_int_mn,
                                      #"mass_start" = mass_start,
                                      "ISO2" = fs_id_link$ISO2[i],
                                      "fs_area" = fs_id_link$fs_area[i],
                                      "fs_name_1" = fs_id_link$fs_name_1[i],
                                      "urbanicity" = fs_id_link$urbanicity[i],
                                      "fs_area_id" = fs_id_link$fs_area_id[i],
                                      "N_species" = N_species,
                                      "CMC_first" = CMC_first,
                                      "CMC_Jan2000" = CMC_Jan2000,
                                      "projection_window_mn" = projection_window_mn,
                                      "N_CMC" = N_CMC,
                                      "tail_pop" = tail_pop,
                                      "sim_population" = sim_population,
                                      "net_cost_strategy_id" = net_cost_strategy_id,
                                      "cost_strategy" = cost_strategy,
                                      "cost_factor" = cost_factor,
                                      "biennial_reduction" = biennial_reduction)
                
                if (use_hipercow) {
                  dynam_id <- paste("id", i, j, k, l, ii, jj, sep = "_")
                  hipercow_params <- param_list[[jj]]
                  if (debugging) {
                    assign(dynam_id,
                           task_create_expr(par_net_region_sequential3(hipercow_params)),
                           envir = .GlobalEnv)
                  } else {
                    assign(dynam_id,
                           task_create_expr(par_net_region_sequential3(hipercow_params)))
                  }
                }
                
                jj <- jj + 1
                
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
  
  if (!use_hipercow) {
    cl <- makeCluster(N_cores)
    clusterExport(cl, c("param_list",
                        "set_bednets",
                        "run_simulation",
                        "fit_usage_sequential",
                        "population_usage_t"))
    #par_output <- lapply(param_list, par_net_region_sequential3)
    par_output <- parLapply(cl, param_list, par_net_region_sequential3)
    comb_output <- do.call(rbind.data.frame, par_output)
    output_df <- rbind(output_df, comb_output)
    stopCluster(cl)
  }
  
  if (use_hipercow) {
    return(NULL)
  } else {
    return(output_df)
  }
}