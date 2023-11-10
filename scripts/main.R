# main.R

#-------------------------------------------------------------------------------
# Libraries required

library(rdhs)
library(magrittr)
library(spatstat.utils)
library(colf)
library(geofacet)
library(ggplot2)
library(stringr)
library(countrycode)
library(tidyverse)
#library(foresite)
library(data.table)
library(plyr)
library(stringi)
library(viridis)
library(scales)
library(rstan)
#library(rstanarm)
library(labelled)
#library(rethinking)

#-------------------------------------------------------------------------------
# Variable inputs

# ISO2 codes for included countries
# Enter in alphabetical order of country name, not two character ISO code
# Currently tested for "BF",	"GH",	"MW",	"ML", "MZ", "SN"
# Other countries may require standardise_names to be updated
SSA_ISO2 <- c("BF",	"GH",	"MW",	"ML", "MZ", "SN")

# Time period
first_year <- 2008
final_year <- 2021

# Urban/rural split
urban_split <- TRUE
urban_split_MDC <- FALSE  # Split by urbanicity for mass campaign timings

# Area usage threshold (exclude areas with fewer than this total number of
# individuals recorded using a net)
area_usage_threshold <- 10

# Rules for local regression curve fitting
MDC_min <- first_year
MDC_max <- final_year

prop_max_kde_mdc <- 0.1   # An MDC must be greater than this proportion of the 

min_kde_int_mdc <- 18     # MDCs must have a minimum spacing of 18 months

local_mode_window <- 9    # Number of preceding and subsequent months compared
# for candidate MDC

peak_window_ratio <- 1 # Minimum ratio between candidate MDC mode and mean
# values over preceding and subsquent window

max_modes <- 5            # Maximum MDCs. If <=0, the value will be set to:
# ceiling(total number of months in time series / 36)

#ksmooth_bandwidth <- 12#12
default_bandwidth <- 12
#dhs_bw <- 12    #DHS net kde bandwidth in months
#dst_bw <- 12    #reference MDC kde bandwidth in months

# NB current version only compatible with admin level MDCs (urbanicity split
# is possible)
MDC_kde_national <- FALSE
MDC_kde_global <- FALSE

DHS_for_MDC <- TRUE
AMP_for_MDC <- FALSE

# Maximum default time since last MDC
max_m <- 72

# Additional antimode selection criteria given properties of first antimode
min_antimode_overall_prop <- 0.2
min_antimode_min_ratio <- 2

# Seed value
set.seed(12345)

# (Weighted) DHS density to use
dhs_den <- "rcpt_grw_w"

#-------------------------------------------------------------------------------
# Rules for estimating MDC timing from reference data

use_ref_data_for_MDCs <- TRUE

#-------------------------------------------------------------------------------
# Load function files

file.sources = list.files("./scripts/utils",
                          pattern="*.R$",
                          full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources, source, .GlobalEnv)

#-------------------------------------------------------------------------------
# reference national ITN distributions

national_itn_data <- read.csv("./data/input_itn_distributions.csv")

#-------------------------------------------------------------------------------
# reference SN admin MDCs

SN_comparison <- read.csv("./data/SN_mdc.csv")

#-------------------------------------------------------------------------------
# Package options

# rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
decay_iter <- 400
decay_warmup <- 200
decay_chains <- 4
decay_init_r <- 2           # default value = 2
decay_adapt_delta <- 0.95   # default values = 0.8

# Private function to set rdhs package credentials using set_rdhs_config()
source("./private/rdhs_creds.R")
call_set_rdhs_config()

#-------------------------------------------------------------------------------
# Generate timestamp

timestamp <- format(Sys.time(), "%y%m%d%H%M")

#-------------------------------------------------------------------------------
# Load in reference data
# Dependencies in reference_data.R

fetch_reference_data(national_itn_data)

#-------------------------------------------------------------------------------
# Extract DHS data
# Dependencies in extraction.R

# Extract data
extracted_surveys <- get_net_data(cc = SSA_ISO2, start_year = first_year)

#-------------------------------------------------------------------------------
# Clean DHS
# Dependencies in cleaning.R

# Clean data
all_net_data <- extracted_surveys %>%
  delabel_data %>%
  standardise_names %>%
  remove_unknown_sleep_location %>%
  remove_low_usage %>%
  generate_unique_ids

# Get global variables
fetch_init_global_vars()

# Generate area data frame
fetch_area_df()

# CMC limits for minimum and maximum net receipt dates. By default these are
# equal to the bounds of the DHS surveys called but can be changed.
CMC_net_min <- CMC_first
CMC_net_max <- CMC_last

#-------------------------------------------------------------------------------
# Usage and access
# Dependencies in usage_access.R unless otherwise indicated

all_net_data %<>%
  append_CMC_net_obtained %>%
  simulate_unknown_net_source %>%
  return_all_access

# Remove DHS data prior to start of MDCs (input countries, years and months as
# vectors). remove_pre_mdc_dhs() found in cleaning.R
all_net_data %<>%
  remove_pre_mdc_dhs("GH", date_to_CMC(year = 2010, month = 1))

# Fetch net data from (total values)
fetch_net_data()
#global_camp_nets <- rep(0, N_CMC)

# Append access
net_data %<>% append_usage_access

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

#-------------------------------------------------------------------------------
# Mass distribution campaigns
# Dependencies in mdc.R unless otherwise indicated

# Append area net decay meanlives and calculate receipt weights
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

# Combine weight density using weighted avg of total sum of dhs weights
if(!urban_split_MDC) {
  # Combine weight density using weighted avg of total sum of dhs weights
  net_data %<>% combine_weights(dhs_den)
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
net_data %<>% 
  

# Estimate MDC timings
net_data %<>%
  estimate_mdc_timings(mdc_bounds_name = "antimodes_ref_nets_norm",
                       density_name = "smth_over_comp_nets_norm")

# Append comparison MDC timings
net_data %<>% append_comparison_mdcs(SN_comparison)

#-------------------------------------------------------------------------------
# Plot MDC timings
# Dependencies in plotting.R over_comp_nets_norm

timestamp <- format(Sys.time(), "%y%m%d%H%M")
net_data %>% plot_MDCs(densities = "over_comp_nets_norm",
                       periods_dataset = mdc_period_df,
                       colvals = "springgreen4",
                       cap_extreme = FALSE,
                       plot_step_dens = TRUE,
                       plot_smth_dens = TRUE,
                       plot_modes = FALSE,
                       plot_antimodes = FALSE,
                       plot_mdc_pts = TRUE,
                       plot_vert_periods = TRUE,
                       plot_comparison_mdc = TRUE)

net_data %<>% combine_weights(density_name = "rcpt_dhs_w",
                              out_name_from_input = TRUE)

net_data %<>% normalise_area_densities("comb_rcpt_grw_w",
                                       norm_over_net_rec_range = FALSE,
                                       time_unit = "years")

timestamp <- format(Sys.time(), "%y%m%d%H%M")
net_data %>% plot_MDCs(densities = "ref_nets",
                       periods_dataset = mdc_period_df,
                       colvals = "royalblue",
                       cap_extreme = FALSE,
                       plot_step_dens = TRUE,
                       plot_smth_dens = TRUE,
                       plot_modes = FALSE,
                       plot_antimodes = FALSE,
                       plot_mdc_pts = FALSE,
                       plot_vert_periods = TRUE,
                       plot_comparison_mdc = FALSE)

timestamp <- format(Sys.time(), "%y%m%d%H%M")
net_data %>% plot_MDCs(densities = c("over_comp_nets_norm","ref_nets_norm"),
                       periods_dataset = mdc_period_df,
                       colvals = c("royalblue","darkorange3"),
                       cap_extreme = FALSE,
                       plot_step_dens = TRUE,
                       plot_smth_dens = TRUE,
                       plot_modes = FALSE,
                       plot_antimodes = FALSE,
                       plot_mdc_pts = TRUE,
                       plot_vert_periods = TRUE,
                       plot_comparison_mdc = TRUE,
                       MDC_density_number = 2)

  
  
  
net_data %>% plot_MDCs

#-------------------------------------------------------------------------------