# mdc.R
# Functions for estimating mass distribution campaign timings

append_access_meanlife <- function(dataset) {
  area_ids <- (dataset$area_id)
  N_area_ids <- length(area_ids)
  dataset$prior_mean_access_meanlife <- rep(NA, dim(dataset)[1])
  for (i in 1:N_area_ids) {
    ids <- which(dataset$area_id == area_ids[i])
    dataset$prior_mean_access_meanlife[ids] <- prior_mean_access_meanlife[i]
  }
  return(dataset)
}

calculate_net_receipt_weights <- function(dataset) {
  # hv005 = dhs weighting
  # hml4 = months since net obtained
  dhs_weight <- dataset$hv005
  net_age <- dataset$hml4
  growth_rate <- 1 / dataset$prior_mean_access_meanlife
  dataset$receipt_weight <- dhs_weight * exp(net_age * growth_rate)
  return(dataset)
}