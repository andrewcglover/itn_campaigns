#-------------------------------------------------------------------------------
# Usage and access Stan fitting
# Dependencies in usage_access_fitting.R

# Number of individuals for beta-binomial sampling
N_bb <- 100000

net_data$MDC_round <- net_data$MDC_round + 1

# Create lists 
create_usage_access_list(usage = TRUE)
create_usage_access_list(usage = FALSE)

# Adjust round number
# usage_list$rho <- usage_list$rho + 1
# access_list$rho <- access_list$rho + 1

# Run Stan models

# Sys.setenv(MAKEFLAGS = paste0("-j",parallel::detectCores()))
# 
# install.packages(c("StanHeaders","rstan"),type="source")

# usage_access_stan_fit(usage = TRUE)
# usage_access_stan_fit(usage = FALSE)

# > r_tau_orig<-usage_list$r_tau
# > usage_list$r_tau<-usage_list$r_tau/4
# usage_list$r_tau <- 2.0 * usage_list$r_tau / usage_list$r_tau

usage_access_cmdstanr_fit(usage = TRUE)
usage_access_cmdstanr_fit(usage = FALSE)

# running to here 07/02/24

# Append mean parameters and credible intervals to net data
#net_data <- net_data[-c(43:dim(net_data)[2])]
# net_data %<>% append_time_series_fits(cmdstanr = TRUE, access = FALSE)

# Create new index following stan runs
net_data$uastan_id <- seq(1, dim(net_data)[1])

extract_time_series_draws()
net_data %<>%
  append_time_series_stats()

net_data %<>% append_time_series_fits(cmdstanr = TRUE)