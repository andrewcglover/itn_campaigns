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

source("./scripts/utils/reg_funs.R")
source("./scripts/utils/cleaning.R")
source("./scripts/utils/plotting.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#SSA_ISO2 <- c("SN")
SSA_ISO2 <- c("BF",	"GH",	"MW",	"ML", "MZ", "SN")


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
urban_split_MDC <- FALSE


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

local_mode_window <- 9     #Number of preceding and subsequent months compared
#for candidate MDC

peak_window_ratio <- 1.05      #Minimum ratio between candidate MDC mode and mean
#values over preceding and subsquent window

max_modes <- ceiling((CMC_MDC_max - CMC_MDC_min) / 36)   #Maximum number of MDCs

dhs_bw <- 12    #DHS net kde bandwidth in months
dst_bw <- 12    #reference MDC kde bandwidth in months

MDC_kde_national <- FALSE
MDC_kde_global <- FALSE

DHS_for_MDC <- TRUE
AMP_for_MDC <- FALSE

#maximum default time since last MDC
max_m <- 72


set.seed(12345)

#-------------------------------------------------------------------------------
# reference national ITN distributions

national_itn_data <- read.csv("./data/input_itn_distributions.csv")


#-------------------------------------------------------------------------------
# Get DHS data

#get data
extracted_surveys <- get_net_data(cc = SSA_ISO2, start_year = first_year)

#clean data
all_net_data <- clean_net_data(extracted_surveys)

#Record unique countries and admin units
uni_ISO2 <- unique(all_net_data$ISO2)
uni_ADM1 <- unique(all_net_data$ADM1NAME)
uni_ADM1_ISO2 <- unique(paste(all_net_data$ISO2,all_net_data$ADM1NAME,sep=" "))

uni_areas <- unique(all_net_data$area)
N_areas <- length(uni_areas)

uni_area_ids <- 1:N_areas

#append area IDs to data frame
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



national_camp_nets <- rep(0, N_CMC)


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
    national_camp_nets[t] <- national_camp_nets[t] + nets_here
    
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
  kde_lt <- kde_mdc(national_camp_nets, CMC_series, CMC_first, CMC_last,
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