# 05_post_use_access_fitting.R
# Code for updating retention estimates, linking data to site files, and
# generating nets per capita curve after use and access timeseries fitting

#-------------------------------------------------------------------------------
# Calculate retention
# Dependencies in retention.R

first_ret_CMC <- date_to_CMC(first_ret_date[1], first_ret_date[2])
last_ret_CMC <- date_to_CMC(last_ret_date[1], last_ret_date[2])

retention_period <- net_data %>%
  fetch_retention_period(CMCa = first_ret_CMC,
                         CMCb = last_ret_CMC)

#-------------------------------------------------------------------------------
# DEPRICATED - Replaced with new site package
# Link data to foresite
# Dependencies in foresite.R

fs_net_data <- net_data %>%
  append_foresite_names(uni_ISO2) %>%
  create_new_foresite_regions(uni_ISO2) %>%
  append_fs_area_names %>%
  append_fs_area_ids

#-------------------------------------------------------------------------------
# Generate nets per capita curve
# Dependencies in npc_stan.R

bv_pred <- stan_npc_fit()
