#use_given_access_test.R

use_given_access_test <- function(data) {
  local_pop_sampled <- NULL
  local_use <- NULL
  local_access <- NULL
  uni_areas_here <- unique(data$area)
  for (i in 1:length(uni_areas_here)) {
    local_data <- data[which(data$area == uni_areas_here[i]),]
    local_pop_sampled[i] <- dim(local_data)[1]
    local_use[i] <- sum(local_data$hml20, na.rm = TRUE)
    local_access[i] <- sum(local_data$access, na.rm = TRUE)
  }
  local_df <- data.frame("N" = local_pop_sampled,
                         "N_used" = local_use,
                         "N_access" = local_access)
  local_df$prop_use <- local_df$N_used / local_df$N
  local_df$prop_access <- local_df$N_access / local_df$N
  local_df$use_given_access <- local_df$prop_use / local_df$prop_access
  return(local_df)
}