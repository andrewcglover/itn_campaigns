# 00_initialisation.R

# Previous libraries:
# devtools::install_github("mrc-ide/malariasimulation@v1.6.0")
#-------------------------------------------------------------------------------
# Libraries required

library(dplyr)
library(magrittr)
library(ggplot2)
library(rdhs)
library(malariasimulation)
library(site)
library(countrycode)
library(rstan)
library(stringi)
library(stringr)
library(spatstat.utils)
library(rstan)

# library(colf)
# library(geofacet)
# library(ggplot2)
# library(stringr)
# library(countrycode)
# library(tidyverse)
# #library(foresite)
# library(data.table)
# library(plyr)
# library(stringi)
# library(viridis)
# library(scales)
# library(rstan)
# #library(rstanarm)
# library(labelled)
# library(cmdstanr)
# #library(rethinking)
# library(foresite)
# library(rdhs)
# library(malariasimulation)
# #library(doSNOW)
# library(parallel)
# library(tictoc)
# library(dplyr)
# #install_github("lmhaile/site")
# library(site)
# #library(devtools)
# #devtools::install_github("mrc-ide/netz@usage_sequential")
# library(netz)
# #library(hipercow)
# library(ggnewscale)
# library(ggrepel)
# #library(grid)
# library(hipercow)
# library(ggh4x)

#-------------------------------------------------------------------------------
# Load function files

file.sources = c(
  list.files(
    "./scripts/utils",
    pattern="*.R$",
    full.names=TRUE,
    ignore.case=TRUE
    ),
  list.files(
    "./scripts/extraction",
    pattern="*.R$",
    full.names=TRUE,
    ignore.case=TRUE
  ),
  list.files(
    "./scripts/use_access_decay",
    pattern="*.R$",
    full.names=TRUE,
    ignore.case=TRUE
  ),
  list.files(
    "./scripts/campaign_estimation",
    pattern="*.R$",
    full.names=TRUE,
    ignore.case=TRUE
  )
)
sapply(file.sources, source, .GlobalEnv)

#-------------------------------------------------------------------------------
# Variable inputs

# ISO2 codes for included countries
# Enter in alphabetical order of country name, not two character ISO code
# Currently tested for "BF",	"GH",	"MW",	"ML", "MZ", "SN"
# Other countries may require standardise_names to be updated
SSA_ISO2 <- c("BF",	"GH", "MW",	"ML", "MZ", "SN")

# Surveys for removal
corrupted_surveys <- NULL #c("GHPR8ADT")

# Time period
first_year <- 2008
final_year <- 2024

# Recorded retention period (enter as vectors of year followed by month)
first_ret_date <- c(2019, 1)
last_ret_date <- c(2024, 12)

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
#dhs_den <- "rcpt_grw_w"

#-------------------------------------------------------------------------------
# Time definitions

year <- 365
DOY_1st <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
DOY_mid <- c(17, 46, 76, 106, 137, 167, 198, 229, 259, 290, 320, 351)
CMC_Jan2000 <- date_to_CMC(2000, 1) # ./utils/formatting.R

#-------------------------------------------------------------------------------
# malariasimulation parameters

# net-types simulated (pyrethroid-)
only <- TRUE
pbo <- TRUE
pyrrole <- TRUE

mass_int_yr <- c(2, 3)
projection_window_yr <- 9

# malsim_cores <- 20

# Reference CMC for plotting and use given access
#ref_CMC <- 1453   #SN = 1453 (2021-1)
ref_CMC <- 1465   #SN = 1453 (2022-1)

sim_population <- 1e5
N_reps <- 100

long_month_offset <- readRDS(
  "./data_public/random_numbers/long_month_offset.Rds"
)
# Generated originally from:
# long_month_offset <- sample.int(13, 10000, replace = TRUE) - 7

#-------------------------------------------------------------------------------
# reference SN admin MDCs

SN_comparison <- read.csv("./data_private/SN_mdc.csv")

#-------------------------------------------------------------------------------
# access vs nets per capita (data from Bertozzi-Villa et al, 2022)

bv_access_npc <- read.csv(
  "./data_public/BertozziVilla2021/fig_4_access_npc.csv"
)

bv_fit <- loess(access_mean ~ percapita_nets_mean,
                data = bv_access_npc,
                span = 0.75)

bv_npc <- seq(min(bv_access_npc$percapita_nets_mean),
              max(bv_access_npc$percapita_nets_mean),
              length.out = 1e4)
bv_access <- predict(bv_fit, bv_npc)

#-------------------------------------------------------------------------------
# rstan options

# general options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# net decay model options
decay_iter <- 4000  # warmup + sampling
decay_warmup <- 2000
decay_chains <- 4
decay_init_r <- 2           # default value = 2
decay_adapt_delta <- 0.999   # default values = 0.8

# usage cmdstanr model options
Ucmd_seed <- 123
Ucmd_init <- 0.1
Ucmd_chains <- 4
Ucmd_parallel_chains <- 4
Ucmd_warmup <- 200
Ucmd_sampling <- 200
Ucmd_refresh <- 50
Ucmc_adapt_delta <- 0.99

# access cmdstanr model options
Acmd_seed <- 123
Acmd_init <- 0.1
Acmd_chains <- 4
Acmd_parallel_chains <- 4
Acmd_warmup <- 200
Acmd_sampling <-200
Acmd_refresh <- 50
Acmc_adapt_delta <- 0.99

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

#national_itn_data <- read.csv("./data_private/MAP_AMP_distrib_2024.csv")
reference_data <- "./data_private/MAP_AMP_distrib_2024.csv" %>%
  read.csv() %>%
  extract_reference_data()

#-------------------------------------------------------------------------------
# Extracted DHS survey path

dhs_surveys_path <- "./data_private/dhs/extracted_surveys_2000_2024.rds"