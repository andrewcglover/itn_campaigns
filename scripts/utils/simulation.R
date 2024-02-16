# Simulation.R

extend_usage_samples <- function(P_samples,
                                 int_days,
                                 month_offset) {
  
  # campaign interval in months
  int_months <- round(12 * int_days / year)
  
  # time index of last predicted campaign
  last_camp_time_id <- apply(P_samples, 1, which.max)
  
  # 
  
}


par_net_region <- function(param_list) {
  
  #year <- 365
  obs_window <- 6 * year
  
#   param_list[[j]] <- c(site_pars,
#                        "mean_ret" = ret_ref_samples[j],
#                        "P" = P_samples[j,],
#                        "P0" = P0_samples[j],
#                        "D" = D_samples[j],
#                        "net_type" = l,
#                        "month_offset" = month_offset[j],
#                        "top_up_int" = top_up_int,
#                        "mass_int" = mass_int[k],
#                        "mass_start" = mass_start,
#                        "area_id" = area_id)
# }
  
  site_pars <- param_list[[1]]
  # print(param_list[[1]])
  # print(param_list$lambda)
  # print(param_list[[1]]$net_offset)
  # P <- site_pars$P
  # P0 <- site_pars$P0
  # D <- site_pars$D
  # lambda <- site_pars$lambda
  sid <- site_pars$sample_index
  mean_retention <- site_pars$mean_ret
  net_type <- site_pars$net_type
  last_camp <- site_pars$last_camp
  month_offset <- site_pars$month_offset
  mass_int_mn_here <- site_pars$mass_int
  mass_start <- site_pars$mass_start
  N_timesteps <- site_pars$timesteps
  area_id <- site_pars$area_id
  
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
  output_usage <- netz::fit_usage_sequential(target_usage = input_net_usage,
                                             target_usage_timesteps = input_net_times,
                                             distribution_timesteps = output_net_times,
                                             mean_retention = mean_retention_dy)
    
  
  # Final predicted month
  
  # Append usage
  apply(which.max)
  
  
  ### bednet campaigns
  
  bednet_times <- round(seq(1 + net_offset, N_timesteps, top_up_int))
  N_dists <- length(bednet_times)
  mass_times <- seq(mass_start + net_offset, N_timesteps, mass_int)
  mass_ids <- which(bednet_times %in% mass_times)
  
  top_up_cov <- P * D * (1 - exp(-lambda * top_up_int))
  
  no_cov_pre_mass <- (1 - P * D - P * (1-D) * exp(-lambda * mass_int))
  mass_cov <- (P - (1 - no_cov_pre_mass)) / no_cov_pre_mass
  
  bednet_cov <- rep(top_up_cov, length(bednet_times))
  bednet_cov[mass_ids] <- bednet_cov[mass_ids] + mass_cov
  
  bednet_cov[which(bednet_cov > 1)] <- 1
  
  original_bednet_times <- site_pars$bednet_timesteps
  
  bednet_time_id <- rep(NA, N_dists)
  for (i in 1:N_dists) {
    bednet_time_id[i] <- which.min(abs(original_bednet_times - bednet_times[i]))
  }
  
  bednet_pars <- set_bednets(
    site_pars,
    timesteps = bednet_times,
    coverages = bednet_cov,
    retention = mean_retention,
    dn0 = site_pars$bednet_dn0[bednet_time_id,],
    rn = site_pars$bednet_rn[bednet_time_id,],
    rnm = site_pars$bednet_rnm[bednet_time_id,],
    gamman = site_pars$bednet_gamman[bednet_time_id]
  )
  
  # run simulation
  
  output <- run_simulation(timesteps = bednet_pars$timesteps,
                           parameters = bednet_pars)
  
  obs_start <- N_timesteps - obs_window
  obs_infections <- sum(output$n_infections[obs_start:N_timesteps])
  annual_infections <- obs_infections * year / obs_window
  
  # return avg annual infections over observation window
  
  output_df <- data.frame("area_id" = area_id,
                          "mass_int" = mass_int,
                          "annual_avg_infections" = annual_infections)
  
  return(output_df)
  
}

run_malsim_nets <- function(dataset,
                            N_reps = 100,
                            mass_int_yr = c(2,3)*365,
                            ref_CMC = 1476,
                            CMC_sim_start = 1297,
                            CMC_sim_end = 1476,
                            CMC_sim_int = c(1,1,2017),
                            only = TRUE,
                            pbo = TRUE,
                            pyrrole = TRUE,
                            int_end) {
  
  # # Sim year set to 360 to facilate monthly conversion to days
  # sim_year <- 360
  
  # Number of samples
  N_samples <- dim(P_u)[1]
  
  # Monthly offset for future mass campaigns
  month_offset <- sample.int(13, N_reps, replace = TRUE) - 7
  
  # Create sample ids
  sample_id <- sample.int(N_samples, N_reps , replace = TRUE)
  
  # dataframe for storing output
  malsim_out <- output_df
  
  for (l in 1:3) {
    
    if ((l==1 & only) | (l==2 & pbo) | (l==3 & pyrrole)) {
  
      for (k in 1:length(mass_int_yr)) {
        
        # mass interval
        mass_int_mn <- mass_int_yr[k] * 12
        
        # current country
        current_country <- "XXX"
        
        for (i in 1:N_fs_areas) {
          
          # Warning for foresite mismatch
          if (fs_id_link$fs_area_id[i] != i) {print("Warning: Foresite id mismatch")}
          
          # Area-time indices
          area_id <- fs_id_link$new_area_id[i]
          area_time_ref_id <- which(net_data$area_id == area_id &
                                  net_data$CMC == ref_CMC)
          area_time_ids <- which(net_data$area_id == area_id)
          
          # # get samples
          # invlam_samples <- invlam_u[sample_id, area_id] %>%
          #   as.vector %>% unname %>% unlist
          # lam_samples <- 1 / invlam_samples
          # P_samples <<- P_u[sample_id, area_time_ids] %>%
          #   as.matrix %>% unname
          # P0_ref_samples <- P0_u[sample_id, area_time_ref_id] %>%
          #   as.vector %>% unname %>% unlist
          # D_ref_samples <- D_u[sample_id, area_time_ref_id] %>%
          #   as.vector %>% unname %>% unlist
          # ret_ref_samples <- ret_u[sample_id, area_time_ref_id] %>%
          #   as.vector %>% unname %>% unlist
          
          # get samples
          invlam_samples <<- invlam_u[sample_id, area_id] %>%
            as.vector %>% unname %>% unlist
          lam_samples <<- 1 / invlam_samples
          ret_ref_samples <- ret_u[sample_id, area_time_ref_id] %>%
            as.vector %>% unname %>% unlist
          P_samples <<- P_u[sample_id, area_time_ids] %>%
            as.matrix %>% unname
          P0_samples <<- P0_u[sample_id, area_time_ids] %>%
            as.matrix %>% unname
          D_samples <<- D_u[sample_id, area_time_ids] %>%
            as.matrix %>% unname
          
          # Identify month with max predicted usage (estimated last mass campaign)
          last_camp_month <- apply(P_samples, 1, which.max)

          
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
            if (ctry_df$urbanicity[i] == "urban") {
              adm_site_index <- which(ctry_site$sites$name_1 == fs_id_link$fs_name_1[i] &
                                        ctry_site$sites$urban_rural == "rural")
            } else {
              adm_site_index <- which(ctry_site$sites$name_1 == fs_id_link$fs_name_1[i] &
                                        ctry_site$sites$urban_rural == "urban")
            }
          }
          
          if (identical(adm_site_index, integer(0))) {
            print(paste0("Warning: foresite not linked for admin region ",
                         fs_id_link$fs_name_1[i], " (index ", i, ") in ",
                         current_country))
          }
            
          adm_site <- site::single_site(ctry_site, adm_site_index)
          
          # Repeat interventions
          adm_site %<>% expand_interventions(expand_year = 6,
                                             delay = 0,
                                             counterfactual = FALSE)
          
          # Pyrethroid resistance
          yearly_res <- adm_site$pyrethroid_resistance$pyrethroid_resistance
          monthly_res <- rep(yearly_res, each = 12)
          round_monthly_res <- round(monthly_res, 2)

          if (l==1) {pyr_res <- res_only}
          if (l==2) {pyr_res <- res_pbo}
          if (l==3) {pyr_res <- res_pyrrole}
          
          res_ids <- match(round_monthly_res, pyr_res$resistance)
          N_species <- length(adm_site$vectors$species)
          
          
          
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
              overrides = list(human_population = sim_population)
            )
            
            # Combine parameters for parLapply function
            
            param_list <- list()
            for (j in 1:N_reps) {
              param_list[[j]] <- c(site_pars,
                                   # #"P_u" = P_samples[j,],
                                   # "P0" = P0_ref_samples[j],
                                   # "D" = D_ref_samples[j],
                                   # "lambda" = lam_samples[j],
                                   "sample_index" = j,
                                   "mean_ret" = ret_ref_samples[j],
                                   "net_type" = l,
                                   "month_offset" = month_offset[j],
                                   "last_camp" = last_camp_month[j],
                                   "top_up_int" = top_up_int,
                                   "mass_int" = mass_int_mn[k],
                                   "mass_start" = mass_start,
                                   "area_id" = area_id,
                                   "N_species" = N_species)
              # param_list[[j]] <- c(site_pars,
              #                      "P_u" = P_samples[j,],
              #                      "P0" = P0_ref_samples[j],
              #                      "D" = D_ref_samples[j],
              #                      "lambda" = lam_samples[j],
              #                      "mean_ret" = ret_ref_samples[j],
              #                      "net_type" = l,
              #                      "month_offset" = month_offset[j],
              #                      "top_up_int" = top_up_int,
              #                      "mass_int" = mass_int[k],
              #                      "mass_start" = mass_start,
              #                      "area_id" = area_id)
            }
            
            cl <- makeCluster(n_cores)
            clusterExport(cl, c("param_list", "set_bednets", "run_simulation"))
            #par_output <- lapply(param_list, par_net_region)
            par_output <- parLapply(cl, param_list, par_net_region)
            comb_output <- do.call(rbind.data.frame, par_output)
            output_df <- rbind(output_df, comb_output)
            
          }
        }
        
      }
    }
  }
  
}