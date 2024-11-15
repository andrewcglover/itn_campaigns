# netz_usage_sequential_branch_funs.R

fit_usage_sequential_distributions <- function(
    target_usage,
    target_usage_timesteps,
    distribution_timesteps, 
    distribution_lower = rep(0, length(distribution_timesteps)), 
    distribution_upper = rep(1, length(distribution_timesteps)),
    mean_retention = 365 * 5
){
  loss_rate <- 1 / mean_retention
  distribution <- rep(0, length(distribution_timesteps))
  start_points <- rep(0, length(distribution_timesteps))
  end_points <- rep(0, length(distribution_timesteps))
  for(t in 1:length(distribution_timesteps)){
    # Usage at time point of next distribution
    put <- population_usage_t(distribution_timesteps[t], distribution, distribution_timesteps, mean_retention)
    # Find next target usage
    time_offset <- target_usage_timesteps  - distribution_timesteps[t]
    if(max(time_offset) < 0){
      distribution[t] <- NA
    } else {
      nearest <- min(time_offset[time_offset >= 0])
      index <- which(time_offset == nearest)
      start_points[t] <- target_usage[index] / exp(-loss_rate * time_offset[index])
      if(index > 1) {
        last_index <- index - 1
        del_t <- distribution_timesteps[index] - distribution_timesteps[last_index]
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
    timesteps,
    distribution,
    distribution_timesteps,
    mean_retention = 365 * 5
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