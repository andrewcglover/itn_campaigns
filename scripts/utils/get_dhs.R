# get_dhs.R
# Function for retrieving DHS data

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