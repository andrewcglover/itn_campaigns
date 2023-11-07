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
first_year <- 2009
final_year <- 2022

# Urban/rural split
urban_split <- TRUE
urban_split_MDC <- FALSE  # Split by urbanicity for mass campaign timings

# Area usage threshold (exclude areas with fewer than this total number of
# individuals recorded using a net)
area_usage_threshold <- 10

# Rules for local regression curve fitting
MDC_min <- 2009
MDC_max <- final_year

prop_max_kde_mdc <- 0.1   # An MDC must be greater than this proportion of the 

min_kde_int_mdc <- 18     # MDCs must have a minimum spacing of 18 months

local_mode_window <- 9    # Number of preceding and subsequent months compared
# for candidate MDC

peak_window_ratio <- 1 # Minimum ratio between candidate MDC mode and mean
# values over preceding and subsquent window

max_modes <- 5            # Maximum MDCs. If <=0, the value will be set to:
# ceiling(total number of months in time series / 36)

ksmooth_bandwidth <- 12#12
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
min_first_antimode_overall_prop <- 0.25
min_first_antimode_min_ratio <- 2


# Seed value
set.seed(12345)

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
# Package options

# rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
decay_iter <- 200
decay_warmup <- 100
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
# Dependencies in usage_access.R

all_net_data %<>%
  append_CMC_net_obtained %>%
  simulate_unknown_net_source %>%
  return_all_access

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

# Generate desired weight density
net_den_base <- "rcpt_grw_w"
if(urban_split_MDC) {
  net_den_MDC <- net_den_base
} else {
  # Combine weight density using weighted avg of total sum of dhs weights
  net_data %<>% combine_weights(net_den_base)
  net_den_MDC <- "urb_comb_w"
}

# Normalise densities
columns_to_normalise <- c("ref_nets", "urb_comb_w")
net_data %<>% normalise_area_densities(columns_to_normalise,
                                       norm_over_net_rec_range = TRUE,
                                       time_unit = "years")

# Estimate MDC timings using smoothing method
net_data %<>%
  mode_smoothing(net_density_name = net_den_MDC) %>%
  identify_antimodes(density_name = net_den_MDC) %>%
  deselect_adjacent_antimodes(density_name = net_den_MDC)

net_data %<>%
  mode_smoothing(net_density_name = "ref_nets") %>%
  identify_antimodes(density_name = "ref_nets") %>%
  additional_early_antimode(density_name = "ref_nets") %>%
  deselect_adjacent_antimodes(density_name = "ref_nets")

# Generate composite density
net_data %<>%
  generate_compostie_density(rec_name = "urb_comb_w",
                             ref_name = "ref_nets",
                             scale_from_means = TRUE,
                             use_predefined_extreme_nets = FALSE) %>%
  normalise_area_densities("comp_nets",
                           norm_over_net_rec_range = FALSE,
                           time_unit = "years") %>%
  mode_smoothing(net_density_name = "comp_nets")

# Estimate MDC timings
net_data %<>%
  estimate_mdc_timings(mdc_bounds_name = "antimodes_ref_nets",
                       density_name = "smth_comp_nets")

timestamp <- format(Sys.time(), "%y%m%d%H%M")
net_data %>% plot_MDCs(densities = "comp_nets",
                       colvals = "springgreen4",
                       cap_extreme = FALSE,
                       plot_step_dens = TRUE,
                       plot_smth_dens = TRUE,
                       plot_modes = FALSE,
                       plot_antimodes = FALSE,
                       plot_mdc_pts = TRUE)


#-------------------------------------------------------------------------------
# Plot MDC timings
# Dependencies in plotting.R

timestamp <- format(Sys.time(), "%y%m%d%H%M")

net_data %>% plot_MDCs(densities = "ref_nets",
                       colvals = "darkorange3",
                       cap_extreme = FALSE,
                       plot_step_dens = TRUE,
                       plot_smth_dens = TRUE,
                       plot_modes = TRUE,
                       plot_antimodes = TRUE),
                       plot_vert_periods = TRUE)

  
  
  
net_data %>% plot_MDCs

#-------------------------------------------------------------------------------








binomial_df <- data.frame("ISO2" = campnets_df$ISO2,
                          "ADM1" = campnets_df$ADM1,
                          "area" = campnets_df$area,
                          "CMC" = campnets_df$CMC,
                          #"post_MDC" = campnets_df$months_post_MDC,
                          #"post_prior_MDC" = campnets_df$months_post_prior_MDC,
                          #"rain" = campnets_df$avg_rain,
                          "used" = campnets_df$used,
                          "access" = campnets_df$access,
                          "total" = campnets_df$used + campnets_df$not_used,
                          #"source_rec" = campnets_df$source_rec,
                          #"camp_rec" = campnets_df$camp_rec,
                          "CTRY" = rep(NA, dim(campnets_df)[1]))#,
                          #"MDC_round" = campnets_df$MDC_round)

# binomial_df$MDC_round[which(is.na(binomial_df$MDC_round))] <- 0
# binomial_df$MDC_round <- binomial_df$MDC_round + 1

# binomial_df <- binomial_df[which(!is.na(binomial_df$post_MDC)),]

unique_areas_included <- unique(binomial_df$area)
binomial_df$area_id <- match(binomial_df$area, unique_areas_included)


#-------------------------------------------------------------------------------

#INSERT DECAY RATE ESTIMATION


# Retain areas remaining after generating DHS weighted distributions
binomial_df <- binomial_df[which(binomial_df$area %in% access_nets_weighted$area),]
binomial_df <- binomial_df[which(binomial_df$area %in% used_nets_weighted$area),]
access_nets_weighted <- access_nets_weighted[which(access_nets_weighted$area %in% binomial_df$area),]
used_nets_weighted <- used_nets_weighted[which(used_nets_weighted$area %in% binomial_df$area),]

# Create data frames for linking individual and net data
adm_ind_link <- data.frame("ADM_id" = binomial_df$area_id,
                           "ISO2" = binomial_df$ISO2)
adm_ind_link <- unique(adm_ind_link)

adm_net_link <- data.frame("ADM_id" = binomial_df$area_id,
                           "ISO2" = binomial_df$ISO2)
adm_net_link <- unique(adm_net_link)

# Link individual and net data
for (i in 1:N_ISO2) {
  binomial_df$CTRY[which(binomial_df$ISO2 == uni_ISO2[i])] <- i
  access_nets_weighted$CTRY[which(access_nets_weighted$ISO2 == uni_ISO2[i])] <- i
  used_nets_weighted$CTRY[which(used_nets_weighted$ISO2 == uni_ISO2[i])] <- i
  adm_ind_link$ISO2[which(adm_ind_link$ISO2 == uni_ISO2[i])] <- i
  adm_net_link$ISO2[which(adm_net_link$ISO2 == uni_ISO2[i])] <- i
}
binomial_df$CTRY <- as.integer(binomial_df$CTRY)
access_nets_weighted$CTRY <- as.integer(access_nets_weighted$CTRY)
used_nets_weighted$CTRY <- as.integer(used_nets_weighted$CTRY)
adm_ind_link$CTRY <- as.integer(adm_ind_link$ISO2)
adm_net_link$CTRY <- as.integer(adm_net_link$ISO2)

uni_indiv_areas <- data.frame("area_id" = binomial_df$area_id,
                              "area" = binomial_df$area)
uni_indiv_areas <- uni_indiv_areas[!duplicated(uni_indiv_areas$area_id),]

access_nets_weighted$area_id <- rep(NA, length(access_nets_weighted$area_id))
used_nets_weighted$area_id <- rep(NA, length(used_nets_weighted$area_id))

for (i in 1:dim(uni_indiv_areas)[1]) {
  access_nets_weighted$area_id[which(access_nets_weighted$area == uni_indiv_areas$area[i])] <- uni_indiv_areas$area_id[i]
  used_nets_weighted$area_id[which(used_nets_weighted$area == uni_indiv_areas$area[i])] <- uni_indiv_areas$area_id[i]
}

#N_t <- max_net_obtained - min_net_obtained + 1

# Function to call Stan code to fit hierarchical net decay
access_decay_samples <- stan_decay_fit(access_nets_weighted, adm_net_link)
used_decay_samples <- stan_decay_fit(used_nets_weighted, adm_net_link)


binomial_df_backup <- binomial_df
binomial_df <- binomial_df[which(binomial_df$CMC >= date_to_CMC(2010,1)),]


#-------------------------------------------------------------------------------
# Record CMC nets obtained


N_CMC <- length(CMC_series)

dates_df <- CMC_to_date(CMC_series)
dates_df[which(dates_df[,2] < 10),2] <- (
  paste("0", dates_df[which(dates_df[,2] < 10), 2], sep = ""))
date_series <- as.Date(paste(dates_df[,1],dates_df[,2],"01",sep="-"),
                       format="%Y-%m-%d")

global_camp_nets <- rep(0, N_CMC)

##unique nets data frame

nets_only <- all_net_data[which(!is.na(all_net_data$netid)),]
nets_only <- nets_only[!duplicated(nets_only$netid),]

campnets_df <- data.frame("area" = rep(uni_areas, each = N_CMC),
                          "area_id" = rep(uni_area_ids, each = N_CMC),
                          "ISO2" = rep(areas_df$ISO2, each = N_CMC),
                          "ADM1" = rep(areas_df$ADM1, each = N_CMC),
                          "urbanicity" = rep(areas_df$urbanicity, each = N_CMC),
                          "ISO2" = rep(areas_df$ISO2, each = N_CMC),
                          "CMC" = rep(CMC_series, N_areas),
                          "Date" = rep(date_series, N_areas),
                          "source_rec" = rep(0, N_areas*N_CMC),
                          "camp_rec" = rep(0, N_areas*N_CMC),
                          "camp_nets_w_pseudo" = rep(0, N_areas*N_CMC),
                          "scaled_camp_nets" = rep(0, N_areas*N_CMC),
                          "uni_nets_ided" = rep(0, N_areas*N_CMC),
                          "comb_net_series" = rep(0, N_areas*N_CMC),
                          "smth_dhs" = rep(0, N_areas*N_CMC),
                          "smth_dist" = rep(0, N_areas*N_CMC),
                          "monthly_national" = rep(0, N_areas*N_CMC),
                          "scaled_national" = rep(0, N_areas*N_CMC),
                          "MDC" = rep(FALSE, N_areas*N_CMC),
                          "MDC_comb_series" = rep(NA, N_areas*N_CMC))




#-------------------------------------------------------------------------------

#Loop over admin-1 units to record estimated MDC CMC dates
pc0 <- 0

for (n in 1:N_areas) {
  ccx <- areas_df$ISO2[n]
  adx <- areas_df$ADM1[n]
  urbx <- areas_df$urbanicity[n]
  
  #subset of net data for admin unit
  if (!urban_split_MDC | is.na(urbx)) {
    admin_data <- all_net_data[which(all_net_data$ISO2 == ccx
                                     & all_net_data$ADM1NAME == adx),]
    admin_nets <- nets_only[which(nets_only$ISO2 == ccx
                                  & nets_only$ADM1NAME == adx),]
  } else {
    admin_data <- all_net_data[which(all_net_data$ISO2 == ccx
                                     & all_net_data$ADM1NAME == adx
                                     & all_net_data$urbanicity == urbx),]
    admin_nets <- nets_only[which(nets_only$ISO2 == ccx
                                  & nets_only$ADM1NAME == adx
                                  & nets_only$urbanicity == urbx),]
  }
  
  #national nets
  ccx3 <- countrycode("SN", origin = 'iso2c', destination = 'iso3c')
  annual_dist_nets <- national_itn_data[which(national_itn_data$ISO3 == ccx3),]
  
  min_net_age_ided <- FALSE
  
  for (t in 1:N_CMC) {
    i <- t + (n - 1) * N_CMC
    #040823 changes here
    source_rec_here <- length(which(admin_data$hml22 <= 3
                                    & admin_data$CMC_net_obtained == CMC_series[t]))
    camp_rec_here <- length(which(admin_data$hml22 == 1
                                  & admin_data$CMC_net_obtained == CMC_series[t]))
    source_rec_survey <- length(which(admin_data$hml22 <= 3
                                      & admin_data$hv008 == CMC_series[t]))
    camp_rec_survey <- length(which(admin_data$hml22 == 1
                                    & admin_data$hv008 == CMC_series[t]))
    nets_here <- length(which(admin_data$all_camp == 1
                              & admin_data$CMC_net_obtained == CMC_series[t]))
    uni_nets_ided_here <- length(which(admin_nets$CMC_net_obtained == CMC_series[t]))
    campnets_df$source_rec_survey[i] <- source_rec_survey
    campnets_df$camp_rec_survey[i] <- camp_rec_survey
    campnets_df$camp_nets_w_pseudo[i] <- nets_here
    campnets_df$uni_nets_ided[i] <- uni_nets_ided_here
    global_camp_nets[t] <- global_camp_nets[t] + nets_here
    
    ref_LLIN_id <- which(annual_dist_nets$year == CMC_to_date(CMC_series[t])[[1]])
    if (identical(ref_LLIN_id, integer(0))) {
      campnets_df$monthly_national[i]
    } else {
      campnets_df$monthly_national[i] <- annual_dist_nets$LLIN[ref_LLIN_id] / 12
    }
    
  }
  
  i_1 <- 1 + (n - 1) * N_CMC
  i_n <- n * N_CMC
  
  admin_camp_nets <- campnets_df$uni_nets_ided[i_1:i_n]
  
  # areas_df$min_net_age_rec[n] <- min(which(admin_camp_nets != 0))
  # areas_df$max_net_age_rec[n] <- max(which(admin_camp_nets != 0))
  # admin_camp_nets[1:(min_age_id-1)] <- NA
  # admin_camp_nets[(max_age_id+1):length(admin_camp_nets)] <- NA
  
  scaled_adm_camp_nets <- admin_camp_nets / sum(admin_camp_nets, na.rm = TRUE)
  campnets_df$scaled_camp_nets[i_1:i_n] <- scaled_adm_camp_nets
  adm_dist_nets <- campnets_df$monthly_national[i_1:i_n]
  scaled_adm_dist_nets <- adm_dist_nets / sum(adm_dist_nets, na.rm = TRUE)
  campnets_df$scaled_national[i_1:i_n] <- scaled_adm_dist_nets
  
  min_age_id <- min(which(admin_camp_nets != 0))
  max_age_id <- max(which(admin_camp_nets != 0))
  areas_df$min_net_age_rec[n] <- CMC_series[min_age_id]
  areas_df$max_net_age_rec[n] <- CMC_series[max_age_id]
  
  
  if (!MDC_kde_global & !MDC_kde_national) {
    kde_lt <- ksmth_fun(DHS_for_MDC, AMP_for_MDC,
                        scaled_adm_camp_nets, scaled_adm_dist_nets,
                        CMC_series, CMC_first, CMC_last,
                        min_kde_mode, prop_max_kde_mdc, max_modes,
                        local_mode_window, peak_window_ratio,
                        min_kde_int_mdc, dhs_bw, dst_bw)
    
    comb_net_series <- kde_lt[[1]]
    selected_nodes_id <- kde_lt[[2]]
    
    campnets_df$comb_net_series[i_1:i_n] <- comb_net_series$comb_net_series 
    campnets_df$smth_dhs[i_1:i_n] <- comb_net_series$smth_dhs.y
    campnets_df$smth_dist[i_1:i_n] <- comb_net_series$smth_dist.y
    campnets_df$MDC[i_1:i_n] <- comb_net_series$selected_nodes
    campnets_df$MDC_comb_series[i_1 - 1 + selected_nodes_id] <- (
      campnets_df$comb_net_series[i_1 - 1 + selected_nodes_id])
  }
  
  pc1 <- round(100 * n / N_areas)
  if (pc1 > pc0) {
    pc0 <- pc1
    print(paste(pc0, "% complete", sep = ""))
  }
}

if (!MDC_kde_global & MDC_kde_national) {
  all_country_df <- data.frame("ISO2" = rep(SSA_ISO2, each = N_CMC),
                               "CMC" = rep(CMC_series, length(SSA_ISO2)),
                               "source_rec" = rep(0, length(SSA_ISO2)*N_CMC),
                               "comb_net_series" = rep(0, length(SSA_ISO2)*N_CMC),
                               "smth_dhs" = rep(0, length(SSA_ISO2)*N_CMC),
                               "smth_dist" = rep(0, length(SSA_ISO2)*N_CMC),
                               "monthly_national" = rep(0, length(SSA_ISO2)*N_CMC),
                               "scaled_national" = rep(0, length(SSA_ISO2)*N_CMC),
                               "MDC" = rep(FALSE, length(SSA_ISO2)*N_CMC),
                               "MDC_comb_series" = rep(NA, length(SSA_ISO2)*N_CMC))
  MDC_series <- NULL
  kde_series <- NULL
  s_kde_series <- NULL
  k0 <- 1
  kk0 <- 1
  for (i in 1:N_ISO2) {
    country_df <- campnets_df[which(campnets_df$ISO2==SSA_ISO2[i]),]
    country_camp_nets <- rep(NA, N_CMC)
    country_dist_nets <- rep(NA, N_CMC)
    for (j in 1:N_CMC) {
      country_camp_nets[j] <- sum(country_df$uni_nets_ided[which(country_df$CMC==CMC_series[j])],
                                  na.rm=TRUE)
      country_dist_nets[j] <- sum(country_df$monthly_national[which(country_df$CMC==CMC_series[j])],
                                  na.rm=TRUE)
    }
    scaled_country_camp_nets <- country_camp_nets / sum(country_camp_nets)
    scaled_country_dist_nets <- country_dist_nets / sum(country_dist_nets)
    kde_lt <- ksmth_fun(DHS_for_MDC, AMP_for_MDC,
                        scaled_country_camp_nets, scaled_country_dist_nets,
                        CMC_series, CMC_first, CMC_last,
                        min_kde_mode, prop_max_kde_mdc, max_modes,
                        local_mode_window, peak_window_ratio,
                        min_kde_int_mdc, dhs_bw, dst_bw)
    
    comb_net_series <- kde_lt[[1]]
    selected_nodes_id <- kde_lt[[2]]
    
    country_areas <- length(unique(country_df$area))
    
    ctry_rep_comb_net_series <- rep(comb_net_series$comb_net_series, country_areas)
    ctry_rep_smth_dhs <- rep(comb_net_series$smth_dhs.y, country_areas)
    ctry_rep_smth_dist <- rep(comb_net_series$smth_dist.y, country_areas)
    ctry_rep_MDC <- rep(comb_net_series$selected_nodes, country_areas)
    ctry_rep_MDC_comb_series <- ctry_rep_comb_net_series * ctry_rep_MDC
    ctry_rep_MDC_comb_series[ctry_rep_MDC_comb_series == 0] <- NA
    
    Nk <- length(ctry_rep_comb_net_series)
    km <- k0 + Nk - 1
    kkm <- kk0 + N_CMC - 1
    
    # campnets_df$comb_net_series[k0:km] <- ctry_rep_comb_net_series
    # campnets_df$smth_dhs[k0:km] <- ctry_rep_smth_dhs
    # campnets_df$smth_dist[k0:km] <- ctry_rep_smth_dist
    # campnets_df$MDC[k0:km] <- ctry_rep_MDC
    # campnets_df$MDC_comb_series[k0:km] <- ctry_rep_MDC_comb_series
    
    ids <- which(campnets_df$ISO2 == uni_ISO2[i])
    campnets_df$comb_net_series[ids] <- ctry_rep_comb_net_series
    campnets_df$smth_dhs[ids] <- ctry_rep_smth_dhs
    campnets_df$smth_dist[ids] <- ctry_rep_smth_dist
    campnets_df$MDC[ids] <- ctry_rep_MDC
    campnets_df$MDC_comb_series[ids] <- ctry_rep_MDC_comb_series
    
    ids <- which(all_country_df$ISO2 == uni_ISO2[i])
    all_country_df$comb_net_series[ids] <- comb_net_series$comb_net_series
    all_country_df$smth_dhs[ids] <- comb_net_series$smth_dhs.y
    all_country_df$smth_dist[ids] <- comb_net_series$smth_dist.y
    all_country_df$MDC[ids] <- comb_net_series$selected_nodes
    all_country_df$MDC_comb_series[ids] <- comb_net_series$selected_nodes * comb_net_series$comb_net_series
    all_country_df$MDC_comb_series[all_country_df$MDC_comb_series == 0] <- NA
    
    k0 <- km + 1
    kk0 <- kkm + 1
    
    # for (n in 1:country_areas) {
    #   i_1 <- 1 + (n - 1) * N_CMC
    #   campnets_df$MDC_comb_series[i_1 - 1 + selected_nodes_id] <- (
    #     campnets_df$comb_net_series[i_1 - 1 + selected_nodes_id])
    # }
    
    # kde_df <- kde_lt[[1]]
    # selected_nodes_id <- kde_lt[[2]]
    # MDC_series_i <- rep(kde_df$selected_nodes, length.out=dim(country_df)[1])
    # MDC_series <- c(MDC_series, MDC_series_i)
    # country_areas <- length(unique(country_df$area))
    # kde_series <- c(kde_series, rep(kde_df$kde,country_areas))
    # s_kde_series <- c(s_kde_series, rep(kde_df$scaled_kde,country_areas))
  }
  
  # campnets_df$MDC <- MDC_series
  # campnets_df$kde <- kde_series
  # campnets_df$scaled_kde <- s_kde_series
}


if (MDC_kde_global & !MDC_kde_national) {
  kde_lt <- kde_mdc(global_camp_nets, CMC_series, CMC_first, CMC_last,
                    min_kde_mode, prop_max_kde_mdc, max_modes,
                    local_mode_window, peak_window_ratio, min_kde_int_mdc)
  
  kde_df <- kde_lt[[1]]
  selected_nodes_id <- kde_lt[[2]]
  
  for (n in 1:N_areas) {
    i_1 <- 1 + (n - 1) * N_CMC
    i_n <- n * N_CMC
    campnets_df$kde[i_1:i_n] <- kde_df$kde
    campnets_df$scaled_kde[i_1:i_n] <- kde_df$scaled_kde
    campnets_df$MDC[i_1:i_n] <- kde_df$selected_nodes
    campnets_df$MDC_scaled_kde[i_1 - 1 + selected_nodes_id] <- (
      campnets_df$scaled_kde[i_1 - 1 + selected_nodes_id])
  }
}

# if (!MDC_kde_national) {
#   campnets_df$scaled_kde[which(campnets_df$scaled_kde < 0)] <- 0 
# }

###
campnets_df$ADM_urb <- substr(campnets_df$area,4,nchar(campnets_df$area))
campnets_df$Year <- CMC_to_date(campnets_df$CMC)[1] + CMC_to_date(campnets_df$CMC)[2] / 12


#-------------------------------------------------------------------------------

plot_MDCs()
#plot_MDCs(campnets_df, urban_split_MDC)