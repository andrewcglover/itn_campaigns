# dhs.R
# Functions for DHS data extraction and formatting

#-------------------------------------------------------------------------------

# Function to return net specific DHS data using the rdhs package
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

#-------------------------------------------------------------------------------

# Extract usage data for campaign-sourced nets
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

#-------------------------------------------------------------------------------
# Append date net obtained in CMC format
 
append_CMC_net_obtained <- function(dataset) {
  rec_months_since_obt <- dataset$hml4
  rec_months_since_obt[which(rec_months_since_obt > 36)] <- NA
  dataset$CMC_net_obtained <- dataset$hv008 - rec_months_since_obt
  return(dataset)
}

#-------------------------------------------------------------------------------
# Generate pseudo-nets for unknown net source
simulate_unknown_net_source <- function(dataset) {
  
  #Find avg proportion from campaigns over SSA given known source
  netsx <- extract_camp_usage(dataset)
  SSA_camp_prop <- netsx$camp/(netsx$camp+netsx$other)
  
  # Simulate net source for unknown
  unknown_source_id <- which(is.na(dataset$hml22) | dataset$hml22==9)
  N_unknown <- length(unknown_source_id)
  rand_vals <- runif(N_unknown, 0, 1)
  pseudo_camp <- rep(0, N_unknown)
  pseudo_camp[which(rand_vals < SSA_camp_prop)] <- 1
  
  # Combine with recorded net source data and record total entries
  dim_net_data <- dim(dataset)
  N_net_data <- dim_net_data[1]
  dataset$pseudo_camp <- rep(NA, N_net_data)
  dataset$pseudo_camp[unknown_source_id] <- pseudo_camp
  dataset$all_camp <- rep(0, N_net_data)
  dataset$all_camp[which(dataset$pseudo_camp == 1)] <- 1
  dataset$all_camp[which(dataset$hml22 == 1)] <- 1
  
  # Return with simulated source nets
  return(dataset)
}