# netz_usage_sequential_branch_funs.R

fit_usage_sequential <- function(
    target_usage,
    target_usage_timesteps,
    distribution_timesteps, 
    distribution_lower = rep(0, length(distribution_timesteps)), 
    distribution_upper = rep(1, length(distribution_timesteps)),
    mean_retention = 365 * 5
){
  loss_rate <- 1 / mean_retention
  distribution <- rep(0, length(distribution_timesteps))
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
      start_point <- target_usage[index] / exp(-loss_rate * time_offset[index])
      distribution[t] <- 1 - (1 - start_point) / (1 - put)
      distribution[t] <- min(distribution_upper[t], distribution[t])
      distribution[t] <- max(distribution_lower[t], distribution[t])
    }
  }
  return(distribution)
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