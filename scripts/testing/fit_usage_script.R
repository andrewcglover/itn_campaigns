#fit_usage_script.R

#-------------------------------------------------------------------------------
# Fit custom usage function

fit_usage_sequential_distributions <- function(
    target_usage = NULL,
    target_usage_timesteps = NULL,
    distribution_timesteps = NULL, 
    distribution_lower = rep(0, length(distribution_timesteps)), 
    distribution_upper = rep(1, length(distribution_timesteps)),
    mean_retention = 365 * 2
){
  loss_rate <- 1 / mean_retention
  distribution <- rep(0, length(distribution_timesteps))
  start_points <- rep(0, length(distribution_timesteps))
  end_points <- rep(0, length(distribution_timesteps))
  for(t in 1:length(distribution_timesteps)){
    # Usage at time point of next distribution
    put <- population_usage_t(distribution_timesteps[t], distribution,
                              distribution_timesteps, mean_retention)
    # Find next target usage
    time_offset <- target_usage_timesteps  - distribution_timesteps[t]
    if(max(time_offset) < 0){
      distribution[t] <- NA
    } else {
      nearest <- min(time_offset[time_offset >= 0])
      index <- which(time_offset == nearest)
      start_points[t] <- target_usage[index] /
        exp(-loss_rate * time_offset[index])
      if(index > 1) {
        last_index <- index - 1
        del_t <- distribution_timesteps[index] -
          distribution_timesteps[last_index]
        end_points[t] <- start_points[t-1] * exp (-loss_rate * del_t)
      } else {
        end_points[t] <- 0
      }
      distribution[t] <- 1 - (1 - start_points[t]) / (1 - put)
      distribution[t] <- min(distribution_upper[t], distribution[t])
      distribution[t] <- max(distribution_lower[t], distribution[t])
    }
  }
  return(list(distribution, start_points, end_points))
}

population_usage_t <- function(
    timesteps = NULL,
    distribution = NULL,
    distribution_timesteps = NULL,
    mean_retention = 365 * 2
){
  loss_rate <- 1 / mean_retention
  
  # Estimate the cumulative usage at distribution time points
  cumulative_usage <- distribution[1]
  if(length(distribution_timesteps) > 1){
    for(t in 2:length(distribution_timesteps)){
      time_offset <- distribution_timesteps[t] - distribution_timesteps[t - 1]
      remaining <- cumulative_usage[t - 1] * exp(-loss_rate * time_offset)
      cumulative_usage[t] <- 1 - (1 - remaining) * (1 - distribution[t])
    }
  }
  
  # Estimate the usage at target time points
  usage <- c()
  for(t in seq_along(timesteps)){
    time_offset <- timesteps[t] - distribution_timesteps
    if(max(time_offset) < 0){
      usage[t] <- 0
    } else {
      nearest <- min(time_offset[time_offset >= 0])
      index <- which(time_offset == nearest)
      usage[t] <- cumulative_usage[index] * exp(-loss_rate * time_offset[index])
    }
  }
  return(usage)
}

#-------------------------------------------------------------------------------
# Input parameters

year <- 365

# Simulation parameters
population <- 1e4
starting_EIR <- 100
N_species <- 1
invlambda <- 1.5 * year       # mean duration of ITN use in days

# target usage series by month
usage <- c(0.126029, 0.124984, 0.123974, 0.122998, 0.122057, 0.121149, 
           0.120273, 0.119429, 0.118616, 0.117833, 0.11708, 0.116356, 0.11566, 
           0.114992, 0.11435, 0.113735, 0.113145, 0.112581, 0.11204, 0.111524, 
           0.11103, 0.110559, 0.11011, 0.109683, 0.109276, 0.108889, 0.108523, 
           0.108175, 0.107847, 0.107536, 0.107244, 0.106968, 0.106709, 0.106467, 
           0.10624, 0.106029, 0.105832, 0.10565, 0.105483, 0.105328, 0.105187, 
           0.105059, 0.104943, 0.10484, 0.104748, 0.104667, 0.104597, 0.104538, 
           0.104489, 0.10445, 0.10442, 0.1044, 0.104389, 0.104386, 0.104391, 
           0.104405, 0.104426, 0.104455, 0.104491, 0.104533, 0.104583, 0.104639, 
           0.104701, 0.104768, 0.177572, 0.481952, 0.659621, 0.717081, 0.699224, 
           0.681871, 0.665009, 0.648625, 0.632707, 0.617243, 0.60222, 0.587628, 
           0.573454, 0.559688, 0.546318, 0.533335, 0.520728, 0.508487, 0.496602, 
           0.485064, 0.473862, 0.462989, 0.452435, 0.442191, 0.432249, 0.4226, 
           0.413238, 0.404152, 0.395337, 0.386785, 0.378487, 0.370438, 0.36263, 
           0.355056, 0.347711, 0.340586, 0.333677, 0.384951, 0.624247, 0.762818, 
           0.803655, 0.78249, 0.761961, 0.742048, 0.722734, 0.704002, 0.685835, 
           0.668218, 0.651133, 0.634567, 0.618504, 0.602929, 0.587829, 0.573188, 
           0.558995, 0.545237, 0.531899, 0.51897, 0.506439, 0.494293, 0.482521, 
           0.471111, 0.460055, 0.449339, 0.438956, 0.428895, 0.419145, 0.409699, 
           0.400547, 0.448574, 0.682008, 0.816561, 0.854661, 0.831434, 0.808928, 
           0.787122, 0.765996, 0.745529, 0.7257, 0.70649, 0.687882, 0.669855, 
           0.652394, 0.63548, 0.619097, 0.603229, 0.58786, 0.572975, 0.558558, 
           0.544597, 0.531076, 0.517982, 0.505302, 0.493023, 0.481133, 0.46962, 
           0.458473, 0.447679, 0.437228, 0.42711, 0.417313, 0.407829, 0.398647, 
           0.389758, 0.381152, 0.372822, 0.364758, 0.422613, 0.6927, 0.848191, 
           0.893142, 0.868302, 0.844255, 0.820976, 0.798441, 0.776627, 0.755511, 
           0.735071, 0.715286, 0.696135, 0.677598, 0.659656, 0.64229, 0.625483, 
           0.609215, 0.59347, 0.578233, 0.563485, 0.549213, 0.5354, 0.522033, 
           0.509097, 0.496579, 0.484464, 0.472741, 0.461397, 0.450419, 0.439797, 
           0.429518, 0.419573, 0.409949, 0.400637, 0.455946, 0.719707, 0.871228, 
           0.914152, 0.888411, 0.863504, 0.839404, 0.816085, 0.793523, 0.771692, 
           0.75057, 0.730133, 0.71036, 0.691229, 0.67272, 0.654812, 0.637487, 
           0.620725, 0.604508, 0.588819, 0.573641, 0.558957, 0.544751, 0.531008, 
           0.517713, 0.504851, 0.492409, 0.480372, 0.487845, 0.515731, 0.502928, 
           0.556105, 0.6335, 0.721516)

# Net parameter values
dn0_val <- 0.340774
rn0_val <- 0.637101
rnm_val <- 0.24
gamman_val <- 2.64 * year / log(2) # Mean duration of insecticidal activity

#-------------------------------------------------------------------------------
# Prepare parameter inputs

# Times
N_months <- length(usage)
output_net_times <- round(seq(0, N_months - 1) * year / 12) + 1 # net dist days
input_net_times <- output_net_times + 15 # target use days
N_timesteps <- ceiling(N_months * year / 12)

# Net parameter vectors
# Note these do not have to be the same value and can be varied appropriately
# for net type and resistance but they must have length of N_months
dn0_vec <- rep(dn0_val, N_months)
rn0_vec <- rep(rn0_val, N_months)
rnm_vec <- rep(rnm_val, N_months)
gam_vec <- rep(gamman_val, N_months)

# Net parameter matrices
# Converts parameter
dn0_mat <- matrix(rep(dn0_vec, N_species),
                  nrow = N_months,
                  ncol = N_species)

rn_mat <- matrix(rep(rn0_vec, N_species),
                 nrow = N_months,
                 ncol = N_species)

rnm_mat <- matrix(rep(rnm_vec, N_species),
                  nrow = N_months,
                  ncol = N_species)

#-------------------------------------------------------------------------------
# Fit custom usage profile

# Fit nets with future mass campaigns
dist_usage <- fit_usage_sequential_distributions(
  target_usage = usage,
  target_usage_timesteps = input_net_times,
  distribution_timesteps = output_net_times,
  mean_retention = invlambda
)

all_output_nets <- dist_usage[[1]]

#-------------------------------------------------------------------------------
# Set malariasimulation parameters

simparams <- malariasimulation::get_parameters(
  list(human_population = population)
)

simparams <- malariasimulation::set_equilibrium(
  parameters = simparams, init_EIR = starting_EIR
)

# set bednets
bednet_pars <- malariasimulation::set_bednets(simparams,
                                              timesteps = output_net_times,
                                              coverages = all_output_nets,
                                              retention = invlambda,
                                              dn0 = dn0_mat,
                                              rn = rn_mat,
                                              rnm = rnm_mat,
                                              gamman = gam_vec)

#-------------------------------------------------------------------------------
# Run simulation

output <- malariasimulation::run_simulation(timesteps = N_timesteps,
                                            parameters = bednet_pars)

#-------------------------------------------------------------------------------
# Plot outputs

# Nets distributed
ylab_str <- "Monthly proportion of population receiving nets"
plot(x = output_net_times / year, y = all_output_nets,
     pch = 20, ylim = c(0,1), xlim = c(0, ceiling(N_timesteps/year)),
     xlab = "Year",
     ylab = ylab_str,
     xaxs = "i", yaxs = "i")
grid(lty = 2, col = "grey80", lwd = 0.5)
axis(side = 1, lty = 1, col = "black", pos = 0)
axis(side = 2, lty = 1, col = "black")

# Proportion using nets
output$prop_use_net <- output$n_use_net / population
plot(x = output$timestep / year, y = output$prop_use_net, type = "l", 
     lwd = 2.5, ylim = c(0,1), xlim = c(0, ceiling(N_timesteps/year)),
     xlab = "Year", ylab = "Proportion of population using ITNs",
     xaxs = "i", yaxs = "i")
lines(x = input_net_times / year, y = usage, type = "l",
      lwd = 2.5, col = "dodgerblue")
grid(lty = 2, col = "grey80", lwd = 0.5)
axis(side = 1, lty = 1, col = "black", pos = 0)
axis(side = 2, lty = 1, col = "black")