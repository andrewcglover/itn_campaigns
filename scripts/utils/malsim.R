# malsim.R
# Functions for malaria simulation model including pre-processing

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