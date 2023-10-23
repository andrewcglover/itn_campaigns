#regression_model.R


library(rdhs)
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

source("./Scripts/reg_funs.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#SSA_ISO2 <- c("SN")
SSA_ISO2 <- c("BF",	"GH",	"ML",	"MW",	"MZ", "SN")


N_ISO2 <- length(SSA_ISO2)

## set up your credentials
## for writing to a project-specific cache
set_rdhs_config(email = "a.glover18@imperial.ac.uk",
                project = "Optimizing ITN deployment frequency",
                #config_path = "~/.rdhs.json",
                config_path = "rdhs.json",
                cache_path = "LLIN_DHS_data",
                timeout = 240,
                password_prompt = FALSE,
                global = FALSE)

#Time period
first_year <- 2009
final_year <- 2022
CMC_first <- date_to_CMC(first_year, 1)
CMC_last <- date_to_CMC(final_year,12)
CMC_series <- CMC_first:CMC_last

#Urban/rural split
urban_split <- TRUE


# Rules for local regression curve fitting
MDC_min <- 2009
MDC_max <- final_year
CMC_MDC_min <- date_to_CMC(MDC_min, 1)
CMC_MDC_max <- date_to_CMC(MDC_max, 12)

prop_max_kde_mdc <- 0.1   #An MDC must be greater than this proportion of the 
#maximum kde value

# min_kde_mode <- 10        #An MDC must have at least this many nets
#distributed in a month from the kde curve

min_kde_int_mdc <- 18     #MDCs must have a minimum spacing of 18 months

local_mode_window <- 12     #Number of preceding and subsequent months compared
#for candidate MDC

peak_window_ratio <- 1      #Minimum ratio between candidate MDC mode and mean
#values over preceding and subsquent window

max_modes <- ceiling((CMC_MDC_max - CMC_MDC_min) / 36)   #Maximum number of MDCs

dhs_bw <- 12    #DHS net kde bandwidth in months
dst_bw <- 12    #reference MDC kde bandwidth in months

MDC_kde_national <- TRUE
MDC_kde_global <- FALSE

#maximum default time since last MDC
max_m <- 72


set.seed(12345)


national_itn_data <- read.csv("./Data/input_itn_distributions.csv")

