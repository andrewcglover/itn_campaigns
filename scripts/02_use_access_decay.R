#-------------------------------------------------------------------------------
# Net decay estimation
# Dependencies stored in net_decay.R unless otherwise indicated

# Generate new distribution of nets based on DHS weightings
used_nets_weighted <- net_weighting_fun(access = FALSE) %>%
  filter_weighted_by_net_data
access_nets_weighted <- net_weighting_fun(access = TRUE) %>%
  filter_weighted_by_net_data

# Store original net_data
original_net_data <- net_data

# Subset net data
net_data <- original_net_data %>%
  subset_net_data %>%
  filter_net_by_weighted_data %>%
  create_new_ids 

# Update global variables
update_global_vars_after_new_ids()

# Add new ids to weighted data and remove rows not linked
used_nets_weighted %<>% append_new_ids %>% remove_area_na
access_nets_weighted %<>% append_new_ids %>% remove_area_na

# Create linking data frame
fetch_area_link(net_data)

# Generate and assign country ids
link_country_ids()

# Fetch oldest and youngest nets
fetch_extreme_nets()                          # Function in cleaning.R

# Run Stan
used_decay_fit <- stan_decay_fit(used_nets_weighted, area_link)
used_decay_samples <- rstan::extract(used_decay_fit)
access_decay_fit <- stan_decay_fit(access_nets_weighted, area_link)
access_decay_samples <- rstan::extract(access_decay_fit)
fetch_decay_summary()

# Check where double recording of access is occurring in all_net_data

# Update ids for original individual data set
original_all_net_data <- all_net_data

# Update all net data with updated area ids
all_net_data <- original_all_net_data %>%
  filter_net_by_weighted_data %>%
  append_new_ids %>%
  remove_area_na

# Fetch statistics for informative priors for usage and access
fetch_prior_access_usage_params()

net_data_02 <- net_data
all_net_data_02 <- all_net_data