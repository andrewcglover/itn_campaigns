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

# Function to return access
return_all_access <- function(data) {
  
  # Total number of individuals
  n_indiv <- length(data$.id)
  
  # Temporary vectors for calculating access
  household <- paste(data$hv002, data$hv001, data$hv024, data$hv008, sep = "_")
  usage <- data$hml20
  slept_there <- data$hv103
  household_nets <- data$hml1
  defacto_members <- data$hv013
  access <- rep(NA, n_indiv)
  
  # Rolling value for potential access for each individual. N.B. One net is
  # assumed to provide access for up to two individuals in the same household.
  pot_access <- household_nets[1] * 2
  
  # Loop over all individuals to determine access
  pc0 <- 0    # Progress counter
  for (i in 1:n_indiv) {
    
    # If a subsequent individual is from a different household, reset pot_access
    if (i > 1) {
      if (household[i] != household[i-1]) {
        pot_access <- household_nets[i] * 2
      }
    }
    
    # Determine access
    if (slept_there[i] == 0) {
      access[i] <- NA
    } else if (pot_access <= 0) {
      access[i] <- 0
    } else {
      access[i] <- 1
      pot_access <- pot_access - 1
    }
    
    # Output progress
    pc1 <- round(100 * i / n_indiv)
    if (pc1 > pc0) {
      pc0 <- pc1
      print(paste("returning access: ", pc0, "% complete", sep = ""))
    }
    
  }
  
  # Combine access with original data and return
  data_acc <- cbind.data.frame(data, access)
  return(data_acc)
}