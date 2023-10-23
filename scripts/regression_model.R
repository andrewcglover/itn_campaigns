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

source("./scripts/reg_funs.R")

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

#-------------------------------------------------------------------------------
# reference national ITN distributions

national_itn_data <- read.csv("./data/input_itn_distributions.csv")


#-------------------------------------------------------------------------------
# Get DHS data

#get data
extract_surveys <- get_net_data(cc = SSA_ISO2, start_year = first_year)
all_net_data <- plyr::ldply(extract_surveys, data.frame)
labelled::remove_val_labels(all_net_data)

#Remove special characters and capitalize admin regions
all_net_data$ADM1NAME <- sub("-"," ",all_net_data$ADM1NAME)
all_net_data$ADM1NAME <- stringr::str_to_title(stri_trans_general(all_net_data$ADM1NAME,
                                                                  "Latin-ASCII"))

#Change all nulls to NA
all_net_data$ADM1NAME[which(all_net_data$ADM1NAME == "Null")] <- NA

#remove NAs
all_net_data <- all_net_data[which(all_net_data$ADM1NAME != "NA"),]
all_net_data <- all_net_data[which(!is.na(all_net_data$ADM1NAME)),]

#Correct for alternative spelling
all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Zuguinchor")]<-"Ziguinchor"
all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Maputo Cidade")]<-"Maputo City"
all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Maputo Province")]<-"Maputo"
all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Maputo Provincia")]<-"Maputo"
all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Toumbouctou")]<-"Tombouctou"

all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Boucle De Mouhoun")]<-"Boucle Du Mouhoun"
all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Hauts Basins")]<-"Hauts Bassins"

# all_net_data$ADM1NAME[which(all_net_data$ISO2=="MW" & all_net_data$ADM1NAME=="North")]<-"Northern"
# all_net_data$ADM1NAME[which(all_net_data$ISO2=="MW" & all_net_data$ADM1NAME=="South")]<-"Southern"

all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="North")]<-"Northern"
all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="South")]<-"Southern"






#Record countries and admin 1 locations with DHS survey data
all_net_data$ISO2 <- substr(all_net_data$SurveyId,1,2)
#ADM1 <- all_net_data$ADM1NAME

#Record unique countries and admin units
uni_ISO2 <- unique(all_net_data$ISO2)
uni_ADM1 <- unique(all_net_data$ADM1NAME)
uni_ADM1_ISO2 <- unique(paste(all_net_data$ISO2,all_net_data$ADM1NAME,sep=" "))

#append urbanicity
all_net_data$urbanicity <- rep(NA, length(all_net_data$hv025))

if (urban_split == TRUE) {
  all_net_data$urbanicity[which(all_net_data$hv025 == 1)] <- "urban"
  all_net_data$urbanicity[which(all_net_data$hv025 == 2)] <- "rural"
}

#unique area code
all_net_data$area <- paste(all_net_data$ISO2,
                           all_net_data$ADM1NAME,
                           all_net_data$urbanicity,
                           sep = " ")

uni_areas <- unique(all_net_data$area)
N_areas <- length(uni_areas)

uni_area_ids <- 1:N_areas

all_net_data$area_ID <- match(all_net_data$area, uni_areas)

#areas data frame
matched_area_IDs <- match(uni_area_ids, all_net_data$area_ID)
areas_df <- data.frame("area" = uni_areas,
                       "area_ID" = uni_area_ids,
                       "ISO2" = substr(uni_areas,1,2),
                       "ADM1" = all_net_data$ADM1NAME[matched_area_IDs],
                       "urbanicity" = all_net_data$urbanicity[matched_area_IDs],
                       "min_net_age_rec" = rep(NA, N_areas),
                       "max_net_age_rec" = rep(NA, N_areas))


#Record total entries
dim_net_data <- dim(all_net_data)
N_net_data <- dim_net_data[1]

#Create household ids
all_net_data$hhid <- paste(all_net_data$.id, all_net_data$hv001,
                           all_net_data$hv002, sep = "_")

#Find avg proportion from campaigns over SSA given known source
netsx <- extract_camp_usage(all_net_data)

SSA_camp_prop <- netsx$camp/(netsx$camp+netsx$other)

#-------------------------------------------------------------------------------
# Simulate net source for unknown

unknown_source_id <- which(is.na(all_net_data$hml22) | all_net_data$hml22==9)
N_unknown <- length(unknown_source_id)
rand_vals <- runif(N_unknown, 0, 1)
pseudo_camp <- rep(0, N_unknown)
pseudo_camp[which(rand_vals < SSA_camp_prop)] <- 1


#-------------------------------------------------------------------------------
# Combine with recorded net source data

all_net_data$pseudo_camp <- rep(NA, N_net_data)
all_net_data$pseudo_camp[unknown_source_id] <- pseudo_camp
all_net_data$all_camp <- rep(0, N_net_data)
all_net_data$all_camp[which(all_net_data$pseudo_camp == 1)] <- 1
all_net_data$all_camp[which(all_net_data$hml22 == 1)] <- 1

#-------------------------------------------------------------------------------
# Combine with recorded net source data

#remove NA values from slept there question (hv103)
all_net_data <- all_net_data[which(!is.na(all_net_data$hv103)),]
all_net_data <- return_all_access(all_net_data)

#-------------------------------------------------------------------------------
# Record CMC nets obtained

rec_months_since_obt <- all_net_data$hml4
rec_months_since_obt[which(rec_months_since_obt > 36)] <- NA
all_net_data$CMC_net_obtained <- all_net_data$hv008 - rec_months_since_obt

N_CMC <- length(CMC_series)

dates_df <- CMC_to_date(CMC_series)
dates_df[which(dates_df[,2] < 10),2] <- (
  paste("0", dates_df[which(dates_df[,2] < 10), 2], sep = ""))
date_series <- as.Date(paste(dates_df[,1],dates_df[,2],"01",sep="-"),
                       format="%Y-%m-%d")

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

national_camp_nets <- rep(0, N_CMC)


##unique nets dataframe
all_net_data$netid <- paste(all_net_data$hhid, all_net_data$hmlidx, sep = "_")
all_net_data$netid[which(is.na(all_net_data$hmlidx))] <- NA

nets_only <- all_net_data[which(!is.na(all_net_data$netid)),]
nets_only <- nets_only[!duplicated(nets_only$netid),]

