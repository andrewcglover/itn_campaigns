#reg_funs.R

### Regression Model Functions ###

#Convert DHS calendar month code (CMC) to year and month
CMC_to_date <- function(CMC){
  yr <- floor((CMC - 1) / 12) + 1900
  mn <- CMC - ((yr-1900) * 12)
  return(data.frame("year" = yr, "month" = mn))
}

#Convert date (year and month) to CMC format
date_to_CMC <- function(year = 2000, month = 1) {
  CMC <- ((year-1900) * 12) + month
  return(CMC)
}


#Function to return net specific DHS data using the rdhs package
get_net_data <- function(cc = "AO", start_year = 2015, end_year = NULL){
  # Selects desired surveys
  survs <- dhs_surveys(countryIds = c(cc), surveyYearStart = start_year,
                       surveyYearEnd = end_year,
                       surveyCharacteristicIds =c(57, 89))
  
  # Selects the desired datasets
  dataset <- dhs_datasets(surveyIds = survs$SurveyId, 
                          fileFormat = "DT", 
                          fileType = "PR")
  
  # Downloads the chosen datasets
  downloads <- get_datasets(dataset$FileName, clear_cache = TRUE)
  
  # Selects the variables that want to investigate
  questions <- search_variables(dataset$FileName, variables = 
                                  c("hv001",   #Cluster number
                                    "hv002",   #Household number
                                    "hv005",   #Household sample weight
                                    "hv008",   #Date of interview (CMC)
                                    "hv008a",  #Date of interview (CDC)
                                    "hv013",   #Total de facto household members
                                    "hv024",   #Region code
                                    "hv025",   #Type of residence (urban/rural)
                                    "hv103",   #Slept last night
                                    "hmlidx",  #Net number
                                    "hml1",    #Number of nets a household owns
                                    "hml4",    #Months ago net obtained
                                    "hml20",   #Slept under a LLIN
                                    "hml21",   #Anyone used
                                    "hml22"    #Source of net
                                  ))
  
  # Extracts the data
  extract <- extract_dhs(questions, add_geo = TRUE)
  
  return (extract)
}

#Extract usage data for campaign-sourced nets
extract_camp_usage <- function(data){
  
  # Percentage people who slept under a LLIN
  all_net <- data[which(data$hml20 == 1),]
  used_weighted <- sum(all_net$hv005/1e6)
  slept_weighted <- sum(data$hv005/1e6)
  
  # Source of net
  all_camp_net <- all_net[which(all_net$hml22 == 1),]
  used_camp_weighted <- sum(all_camp_net$hv005/1e6)
  
  all_other_net <- all_net[which(all_net$hml22 == 0
                                 | all_net$hml22 == 2
                                 | all_net$hml22 == 3),]
  used_other_weighted <- sum(all_other_net$hv005/1e6)
  
  # Return result
  return (data.frame("used" = used_weighted,
                     "slept" = slept_weighted,
                     "camp" = used_camp_weighted,
                     "other" = used_other_weighted))
}



return_all_access <- function(data) {
  
  n_indiv <- length(data$.id)
  
  household <- paste(data$hv002, data$hv001, data$hv024, data$hv008, sep = "_")
  usage <- data$hml20
  slept_there <- data$hv103
  household_nets <- data$hml1
  defacto_members <- data$hv013
  
  access <- rep(NA, n_indiv)
  
  pot_access <- household_nets[1] * 2
  
  if (slept_there[1] == 0) {
    access[1] <- NA
  } else if (pot_access <= 0) {
    access[1] <- 0
  } else {
    access[1] <- 1
    pot_access <- pot_access - 1
  }
  
  pc0 <- 0
  for (i in 2:n_indiv) {
    
    if (household[i] == household[i-1]) {
      pot_access <- household_nets[i] * 2
    }
    
    if (slept_there[i] == 0) {
      access[i] <- NA
    } else if (pot_access <= 0) {
      access[i] <- 0
    } else {
      access[i] <- 1
      pot_access <- pot_access - 1
    }
    
    pc1 <- round(100 * i / n_indiv)
    if (pc1 > pc0) {
      pc0 <- pc1
      print(paste("returning access: ", pc0, "% complete", sep = ""))
    }
    
  }
  
  data_acc <- cbind.data.frame(data, access)
  return(data_acc)
}




net_weighting_fun <- function(all_net_data, CMC_net_min, CMC_net_max,
                              access = FALSE) {
  
  nets <- data.frame("survey_id" = all_net_data$SurveyId,
                     "ISO2" = all_net_data$ISO2,
                     "ADM1" = all_net_data$ADM1NAME,
                     "urbanicity" = all_net_data$urbanicity,
                     "area" = all_net_data$area,
                     "area_id" = all_net_data$area_ID,
                     "cluster" = all_net_data$hv001,
                     "household" = all_net_data$hv002,
                     "household_w" = all_net_data$hv005,
                     "CMC_interview" = all_net_data$hv008,
                     "hh_net_id" = all_net_data$hmlidx,
                     "months_ago_obtained" = all_net_data$hml4,
                     "CMC_obtained" = all_net_data$raw_CMC_net_obtained)
  
  if (access) {
    nets <- data.frame(nets,
                       "anyone_access" = all_net_data$access)
    nets <- nets[which(nets$anyone_access == 1),]
  } else {
    nets <- data.frame(nets,
                       "anyone_used" = all_net_data$hml21)
    nets <- nets[which(nets$anyone_used == 1),]
  }
  nets <- nets[which(!is.na(nets$CMC_obtained)),]
  nets$net_id <- paste(nets$area_id,
                       nets$cluster,
                       nets$household,
                       nets$hh_net_id,
                       sep = ".")
  nets <- nets[!duplicated(nets$net_id),]
  nets <- nets[which(nets$CMC_obtained >= CMC_net_min),]
  nets <- nets[which(nets$CMC_obtained <= CMC_net_max),]
  nets <- nets[which(nets$months_ago_obtained <= 96),]
  nets$CMC_obt_int_area_id <- paste(nets$CMC_obtained,
                                    nets$CMC_interview,
                                    nets$area_id,
                                    sep = ".")
  
  uni_CMC_obt_int_area_id <- unique(nets$CMC_obt_int_area_id)
  N_ucoia <- length(uni_CMC_obt_int_area_id)
  
  nets_totals <- nets[!duplicated(nets$CMC_obt_int_area_id),]
  if (access) {
    nets_totals = subset(nets_totals,
                         select = -c(cluster,
                                     household,
                                     household_w,
                                     hh_net_id,
                                     anyone_access,
                                     net_id))
  } else {
    nets_totals = subset(nets_totals,
                         select = -c(cluster,
                                     household,
                                     household_w,
                                     hh_net_id,
                                     anyone_used,
                                     net_id))
  }
  nets_totals$weighted_counts <- rep(NA, N_ucoia)
  
  for (i in 1:N_ucoia) {
    matched_ids <- which(nets$CMC_obt_int_area_id == nets_totals$CMC_obt_int_area_id[i])
    weighted_counts_here <- round (sum(nets$household_w[matched_ids]) / 1e6)
    nets_totals$weighted_counts[i] <- weighted_counts_here
  }
  
  nets_totals <- nets_totals[which(nets_totals$weighted_counts > 0),]
  
  M <- dim(nets_totals)[1]
  N <- sum(nets_totals$weighted_counts)
  
  nets_weighted <- data.frame("ISO2" = rep(NA, N),
                              "ADM1" = rep(NA, N),
                              "urbanicity" = rep(NA, N),
                              "area" = rep(NA, N),
                              "area_id" = rep(NA, N),
                              "net_obtained" = rep(NA, N),
                              "months_since_obtained" = rep(NA, N),
                              "CTRY" = rep(NA, N))
  
  j <- 1
  for (i in 1:M) {
    counts <- nets_totals$weighted_counts[i]
    k <- j + counts - 1
    nets_weighted$ISO2[j:k] <- nets_totals$ISO2[i]
    nets_weighted$ADM1[j:k] <- nets_totals$ADM1[i]
    nets_weighted$urbanicity[j:k] <- nets_totals$urbanicity[i]
    nets_weighted$area[j:k] <- nets_totals$area[i]
    nets_weighted$area_id[j:k] <- nets_totals$area_id[i]
    nets_weighted$net_obtained[j:k] <- nets_totals$CMC_obtained[i]
    nets_weighted$months_since_obtained[j:k] <- nets_totals$months_ago_obtained[i]
    j <- j + counts
  }
  
  return(nets_weighted)
  
}

stan_decay_fit <- function(nets_weighted, adm_net_link) {
  
  N_a <- length(unique(nets_weighted$area))
  N_c <- N_ISO2
  
  net_decay_dat <- list(N = dim(nets_weighted)[1],
                        N_a = N_a,
                        N_c = N_c,
                        a = nets_weighted$area_id,
                        c = nets_weighted$CTRY,
                        cc = adm_net_link$CTRY,
                        m = nets_weighted$months_since_obtained)
  
  net_decay_fit <- stan('exp_model_d2.stan',
                        data = net_decay_dat,
                        iter = 1000,
                        warmup = 500,
                        chains = 4,
                        #init_r = 1e-2,
                        control = list(adapt_delta = 0.95))
  
  net_decay_samples <- extract(net_decay_fit)
  
  return(net_decay_samples)
  
}



