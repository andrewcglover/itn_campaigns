# malsim.R
# Functions for malaria simulation model including pre-processing

fit_nets <- function(input_usage,
                     input_times,
                     mean_retention,
                     output_times) {
  output_usage <- netz::fit_usage_sequential(target_usage = input_usage,
                                             target_usage_timesteps = input_times,
                                             distribution_timesteps = output_times,
                                             mean_retention = mean_retention)
}

par_net_region <- function(param_list) {
  
  year <- 365
  obs_window <- 6 * year
  
  site_pars <- param_list#[[1]]
  # print(param_list[[1]])
  # print(param_list$lambda)
  # print(param_list[[1]]$net_offset)
  lambda <- site_pars$lambda
  P0 <- site_pars$P0
  D <- site_pars$D
  mean_retention <- site_pars$mean_ret
  net_offset <- site_pars$net_offset
  top_up_int <- site_pars$top_up_int
  mass_int <- site_pars$mass_int
  mass_start <- site_pars$mass_start
  N_timesteps <- site_pars$timesteps
  area_id <- site_pars$area_id
  
  # fit usage
  
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
                            mass_int = c(2,3),
                            ref_CMC = 1476,
                            CMC_sim_start = 1297,
                            CMC_sim_end = 1476,
                            CMC_sim_int = c(1,1,2017),
                            int_end) {
  
  # # Sim year set to 360 to facilate monthly conversion to days
  # sim_year <- 360
  
  # Number of samples
  N_samples <- dim(P_u)[1]
  
  # Day of year for mass campaigns
  net_offset <- sample.int(year, N_reps, replace = TRUE) - 1
  
  # Create sample ids
  sample_id <- sample.int(N_samples, N_reps , replace = TRUE)
  
  # dataframe for storing output
  malsim_out <- output_df
  
  
  for (k in 1:length(mass_int)) {
    
    # current country
    current_country <- NULL
    
    for (i in 1:N_fs_areas) {
      
      # Warning for foresite mismatch
      if (fs_id_link[i] != i) {print("Warning: Foresite id mismatch")}
      
      # Area-time index
      area_id <- fs_id_link$new_area_id[i]
      area_time_id <- which(net_data$area_id == area_id &
                              net_data$CMC == ref_CMC)
      
      # get samples
      # inv_lambda_samples <- growth_fit_samples$inv_lambda[sample_id, area_id]
      ret_samples <- ret_u[sample_id, area_time_id]
      P0_samples <- P0_u[sample_id, area_time_id] %>%
        as.vector %>% unname %>% unlist
      D_samples <- D_u[sample_id, area_time_id] %>%
        as.vector %>% unname %>% unlist
      
      # Generate ISO code for current admin
      admin_country <- countrycode(fs_id_link$ISO2, "iso2c", "iso3c")
      
      # Pull country site
      if (admin_country != current_country) {
        current_country <- admin_country
        ctry_site <- get_site(current_country)
      }
      
      # Isolate a single site from a country
      adm_site_index <- which(ctry_site$sites$name_1 == ctry_df$fs_name_1[i] &
                                ctry_site$sites$urban_rural == ctry_df$urbanicity[i])
      if (identical(adm_site_index, integer(0))) {
        if (ctry_df$urbanicity[i] == "urban") {
          adm_site_index <- which(ctry_site$sites$name_1 == ctry_df$fs_name_1[i] &
                                    ctry_site$sites$urban_rural == "rural")
        } else {
          adm_site_index <- which(ctry_site$sites$name_1 == ctry_df$fs_name_1[i] &
                                    ctry_site$sites$urban_rural == "urban")
        }
      }
      adm_site <- site::single_site(ctry_site, adm_site_index)
      
      adm_site_new <- expand_interventions(adm_site, 6, 0, FALSE)
      
      # Pf EIR
      Pf_eir <- adm_site$eir$eir[1]
      
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
                               "mean_ret" = ret_samples[j],
                               "P0" = P0_samples[j],
                               "D" = D_samples[j],
                               "net_offset" = net_offset[j],
                               "top_up_int" = top_up_int,
                               "mass_int" = mass_int[k],
                               "mass_start" = mass_start,
                               "area_id" = area_id)
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












adm_site_index <- which(ctry_site$sites$name_1 == "Dakar" &
                          ctry_site$sites$urban_rural == "urban")


# Function for running malaria simulation in parallel with bednets
par_net_region <- function(param_list,
                           obs_years = 6,
                           overide_all_net_dists = FALSE) {
  
  year <- 365
  obs_window <- obs_years * year
  
  # Fetch net-specific parameters
  site_pars <- param_list
  lambda <- site_pars$lambda
  P <- site_pars$P
  D <- site_pars$D
  net_offset <- site_pars$net_offset
  top_up_int <- site_pars$top_up_int
  mass_int <- site_pars$mass_int
  mass_start <- site_pars$mass_start
  N_timesteps <- site_pars$timesteps
  area_id <- site_pars$area_id
  
  mean_retention <- 1 / (lambda * (1-D))
  
  ### bednet campaigns
  if (overide_all_net_dists) {
    bednet_times <- round(seq(1 + net_offset, N_timesteps, top_up_int))
  }
  
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