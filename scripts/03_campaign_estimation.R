# 03_campaign_estimation.R

#-------------------------------------------------------------------------------
# Mass distribution campaigns
# Dependencies in mdc.R unless otherwise indicated

# Append mean net retention by area and calculate receipt weights
all_net_data %<>%
  append_access_meanlife %>%
  calculate_net_receipt_weights

# Append weights to net data totals dataframe
net_data %<>%
  append_total_weights_by_interview_date %>%
  append_weight_window %>%
  append_total_receipt_weights %>%
  append_adj_receipt_weights %>%
  append_reference_nets                       # Function in reference_data.R

net_data_03a <- net_data
all_net_data_03a <- all_net_data

# Combine weight density using weighted avg of total sum of dhs weights
if(urban_split_MDC) {
  net_data %<>% combine_weights("rcpt_grw_w")
  dhs_den <- "rcpt_grw_w"
} else {
  # Combine weight density using weighted avg of total sum of dhs weights
  net_data %<>% combine_weights("rcpt_grw_w")
  dhs_den <- "urb_comb_w"
}

# Normalise densities
columns_to_normalise <- c("ref_nets", dhs_den)
net_data %<>% normalise_area_densities(columns_to_normalise,
                                       norm_over_net_rec_range = FALSE,
                                       time_unit = "years")

# Smooth reference density and identify MDC regions
net_data %<>%
  mode_smoothing("ref_nets") %>%
  identify_antimodes("ref_nets") %>%
  add_antimodes_near_bounds("ref_nets", early_antimode = TRUE) %>%
  add_antimodes_near_bounds("ref_nets", early_antimode = FALSE) %>%
  deselect_adjacent_antimodes("ref_nets")
net_data %<>%
  mode_smoothing("ref_nets_norm") %>%
  identify_antimodes("ref_nets_norm") %>%
  add_antimodes_near_bounds("ref_nets_norm", early_antimode = TRUE) %>%
  add_antimodes_near_bounds("ref_nets_norm", early_antimode = FALSE) %>%
  deselect_adjacent_antimodes("ref_nets_norm")

#additional_early_antimode("ref_nets_norm") %>%

# Fetch mdc period dataframe
net_data %>% fetch_mdc_period_df("antimodes_ref_nets_norm")

# Generate composite density
net_data %<>%
  generate_compostie_density(rec_name = "urb_comb_w_norm",
                             ref_name = "ref_nets_norm",
                             scale_from_means = TRUE,
                             use_predefined_extreme_nets = FALSE) %>%
  overide_comp_density_sections(ref_name = "ref_nets_norm") %>%
  normalise_area_densities("over_comp_nets",
                           norm_over_net_rec_range = FALSE,
                           time_unit = "years") %>%
  mode_smoothing("over_comp_nets_norm")

# Generate mixture densities
#net_data %<>% 


# Estimate MDC timings
N_mdc_uncert_bands <- 3
tau_rank_vals <- c(1, 2, 3)#c(1, 1.5, 2)
net_data %<>%
  estimate_mdc_timings(mdc_bounds_name = "antimodes_ref_nets_norm",
                       density_name = "smth_over_comp_nets_norm",
                       append_uncertainty = TRUE,
                       append_ranked_tau = TRUE)
# net_data %<>%
#   estimate_mdc_timings(mdc_bounds_name = "antimodes_ref_nets_norm",
#                        density_name = "smth_over_comp_nets_norm",
#                        append_uncertainty = TRUE,
#                        uncertainty_bands = N_mdc_uncert_bands)

# Append comparison MDC timings
net_data %<>% append_comparison_mdcs(SN_comparison)

# Estimate uncertainty around MDC timings

#-------------------------------------------------------------------------------
# Plot MDC timings
# Dependencies in plotting.R over_comp_nets_norm

#net_data %>% generate_mdc_plots

#-------------------------------------------------------------------------------
# Number MDC rounds
# Dependencies in mdc_rounds.R

net_data %<>% append_mdc_rounds
unique_areas_included_check()
# generate_MDC_round_matrices(max_tau = 12)
# max_tau value acts as a placeholder for non-observed mass campaigns
# the maximum tau (standard deviation) of observed mass campaigns will equal
# max_tau - 1
matrix_list <- generate_MDC_round_matrices(use_ranked_tau = TRUE, max_tau = 4)
MDC_matrix <- matrix_list[[1]]
MDC_tau_matrix <- matrix_list[[2]]
max_rounds <- dim(MDC_matrix)[2]

#-------------------------------------------------------------------------------
# Plot MDC timings
# Dependencies in plotting.R over_comp_nets_norm

#net_data %>% generate_mdc_plots
