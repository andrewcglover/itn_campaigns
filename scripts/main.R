# main.R

#-------------------------------------------------------------------------------
# Libraries required


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
library(cmdstanr)
#library(rethinking)
library(foresite)
library(rdhs)

#-------------------------------------------------------------------------------
# Variable inputs

# ISO2 codes for included countries
# Enter in alphabetical order of country name, not two character ISO code
# Currently tested for "BF",	"GH",	"MW",	"ML", "MZ", "SN"
# Other countries may require standardise_names to be updated
SSA_ISO2 <- c("BF",	"GH", "MW",	"ML", "MZ", "SN")
#SSA_ISO2 <- c("GH", "MW",	"ML", "MZ", "SN")

# Surveys for removal
corrupted_surveys <- c("GHPR8ADT")

# Time period
first_year <- 2008
final_year <- 2022

# Recorded retention period (enter as vectors of year followed by month)
first_ret_date <- c(2016, 7)
last_ret_date <- c(2022, 6)

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
# malariasimulation parameters

n_cores <- 12

ISO2 <- "GH"
ISO3 <- "GHA"

ref_CMC <- 1453   #SN = 1453 (2021-1)
cal_year <- 2021

sim_population <- 2000

N_reps <- 500
year <- 365

top_up_int <- year / 12
mass_int <- c(2, 3) * year
mass_start <- 5 * year + 1

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
# rstan options

# general options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# net decay model options
decay_iter <- 600
decay_warmup <- 500
decay_chains <- 16
decay_init_r <- 2           # default value = 2
decay_adapt_delta <- 0.95   # default values = 0.8

# usage cmdstanr model options
Ucmd_seed <- 123
Ucmd_init <- 0.5
Ucmd_chains <- 16
Ucmd_parallel_chains <- 16
Ucmd_warmup <- 250
Ucmd_sampling <- 50
Ucmd_refresh <- 25

# access cmdstanr model options
Acmd_seed <- 123
Acmd_init <- 0.5
Acmd_chains <- 16
Acmd_parallel_chains <- 16
Acmd_warmup <- 250
Acmd_sampling <- 50
Acmd_refresh <- 25

#-------------------------------------------------------------------------------
# rdhs options

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
extracted_surveys <- get_net_data(cc = SSA_ISO2,
                                  start_year = first_year,
                                  end_year = final_year)

# Remove any corrupted surveys
retained_surveys <- !(names(extracted_surveys) %in% corrupted_surveys)
extracted_surveys <- extracted_surveys[retained_surveys]

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
# all_net_data %<>%
#   remove_pre_mdc_dhs("GH", date_to_CMC(year = 2010, month = 1))
# Retracted as of 09/01/2024 - function removing data with NA for CMC_net_obtained

# Fetch net data from (total values)
fetch_net_data()
#global_camp_nets <- rep(0, N_CMC)

# Append access, usage and net source information
net_data %<>% append_net_info

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
#net_data %<>% 
  

# Estimate MDC timings
N_mdc_uncert_bands <- 3
tau_rank_vals <- c(1, 1.5, 2)
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
matrix_list <- generate_MDC_round_matrices(use_ranked_tau = TRUE, max_tau = 2)
MDC_matrix <- matrix_list[[1]]
MDC_tau_matrix <- matrix_list[[2]]
max_rounds <- dim(MDC_matrix)[2]

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

usage_access_cmdstanr_fit(usage = TRUE)
usage_access_cmdstanr_fit(usage = FALSE)

# Append mean parameters and credible intervals to net data
net_data <- net_data[-c(43:dim(net_data)[2])]
# net_data %<>% append_time_series_fits(cmdstanr = TRUE, access = FALSE)

extract_time_series_draws()
net_data %<>% append_time_series_stats()

net_data %<>% append_time_series_fits(cmdstanr = TRUE)


#-------------------------------------------------------------------------------
# Calculate retention
# Dependencies in retention.R

first_ret_CMC <- date_to_CMC(first_ret_date[1], first_ret_date[2])
last_ret_CMC <- date_to_CMC(last_ret_date[1], last_ret_date[2])

retention_period <- net_data %>%
  fetch_retention_period(CMCa = first_ret_CMC,
                         CMCb = last_ret_CMC)

#-------------------------------------------------------------------------------
# Link data to foresite
# Dependencies in foresite.R

net_data %<>%
  append_foresite_names(uni_ISO2) %>%
  create_new_foresite_regions(uni_ISO2) %>%
  append_fs_area_names %>%
  append_fs_area_ids
  

#-------------------------------------------------------------------------------
# Malaria Simulation

net_data %<>%
  




#-------------------------------------------------------------------------------
# Usage and access plotting
# Dependencies in usage_access_plotting.R

#net_data %>% plot_usage("BF")

#-------------------------------------------------------------------------------
# Foresite areas
# Dependencies in foresite.R

