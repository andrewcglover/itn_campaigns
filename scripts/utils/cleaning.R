# cleaning.R
# Data cleaning functions

# Remove labels
delabel_data <- function(dataset) {
  all_net_data <- plyr::ldply(dataset, data.frame)
  labelled::remove_val_labels(all_net_data)
}

# Standardise names
standardise_names <- function(dataset) {
  # Remove special characters and capitalize admin regions
  dataset$ADM1NAME <- sub("-"," ",dataset$ADM1NAME)
  dataset$ADM1NAME <- stringr::str_to_title(stri_trans_general(dataset$ADM1NAME,
                                                                    "Latin-ASCII"))
  # Change all nulls to NA
  dataset$ADM1NAME[which(dataset$ADM1NAME == "Null")] <- NA
  
  # Remove NAs
  dataset <- dataset[which(dataset$ADM1NAME != "NA"),]
  dataset <- dataset[which(!is.na(dataset$ADM1NAME)),]
  
  # Correct for alternative spelling
  dataset$ADM1NAME[which(dataset$ADM1NAME=="Zuguinchor")]<-"Ziguinchor"
  dataset$ADM1NAME[which(dataset$ADM1NAME=="Maputo Cidade")]<-"Maputo City"
  dataset$ADM1NAME[which(dataset$ADM1NAME=="Maputo Province")]<-"Maputo"
  dataset$ADM1NAME[which(dataset$ADM1NAME=="Maputo Provincia")]<-"Maputo"
  dataset$ADM1NAME[which(dataset$ADM1NAME=="Toumbouctou")]<-"Tombouctou"
  
  dataset$ADM1NAME[which(dataset$ADM1NAME=="Boucle De Mouhoun")]<-"Boucle Du Mouhoun"
  dataset$ADM1NAME[which(dataset$ADM1NAME=="Hauts Basins")]<-"Hauts Bassins"
  
  dataset$ADM1NAME[which(dataset$ADM1NAME=="North")]<-"Northern"
  dataset$ADM1NAME[which(dataset$ADM1NAME=="South")]<-"Southern"
  
  # Append ISO2 code
  dataset$ISO2 <- substr(dataset$SurveyId,1,2)
  
  #append urbanicity
  dataset$urbanicity <- rep(NA, length(dataset$hv025))
  if (urban_split == TRUE) {
    dataset$urbanicity[which(dataset$hv025 == 1)] <- "urban"
    dataset$urbanicity[which(dataset$hv025 == 2)] <- "rural"
  }
  
  # Return clean named data
  return(dataset)
}

# Remove unknown location of where individual slept last night
remove_unknown_sleep_location <- function(dataset) {
  #remove NA values from slept there question (hv103)
  dataset <- dataset[which(!is.na(dataset$hv103)),]
  return(dataset)
}

# Function to remove areas with low usage
remove_low_usage <- function(dataset) {
  
  # Generate temporary ids
  if (urban_split) {
    temp_ids <- paste0(dataset$ISO2, dataset$ADM1NAME, dataset$hv025)
  } else {
    temp_ids <- paste0(dataset$ISO2, dataset$ADM1NAME)
  }
  id <- unique(temp_ids)
  
  # Sum total recorded usage and remove areas below usage threshold
  for (i in 1:length(id)) {
    total_rec_usage <- sum(dataset$hml20[which(temp_ids == id[i])], na.rm=TRUE)
    if (total_rec_usage < area_usage_threshold) {
      dataset <- dataset[which(temp_ids != id[i]),]
    }
  }
  
  # Return dataset with low usage areas removed
  return(dataset)
}

# Generate unique identifiers
generate_unique_ids <- function(dataset){
  # Unique area code
  dataset$area <- paste(dataset$ISO2, dataset$ADM1NAME, dataset$urbanicity,
                        sep = " ")
  uni_areas <- unique(dataset$area)
  dataset$area_id <- match(dataset$area, uni_areas)
  
  # Create household ids
  dataset$hhid <- paste(dataset$.id, dataset$hv001, dataset$hv002, sep = "_")
  
  # Create net ids
  dataset$netid <- paste(dataset$hhid, dataset$hmlidx, sep = "_")
  dataset$netid[which(is.na(dataset$hmlidx))] <- NA
  
  # Print to state initial data cleaning complete
  print("initial data cleaning: 100% complete")
  
  # Return dataset with unique IDs
  return(dataset)
}

#-------------------------------------------------------------------------------
# Return global variables

fetch_init_global_vars <- function() {
  # Variables from inputs
  
  # Number of countries
  N_ISO2 <<- length(SSA_ISO2)
  
  # Convert dates to DHS calendar month code format
  CMC_first <<- date_to_CMC(first_year, 1)
  CMC_last <<- date_to_CMC(final_year,12)
  CMC_series <<- CMC_first:CMC_last
  CMC_MDC_min <<- date_to_CMC(MDC_min, 1)
  CMC_MDC_max <- date_to_CMC(MDC_max, 12)
  N_CMC <<- length(CMC_series)
  
  # Date series
  dates_df <<- CMC_to_date(CMC_series)
  dates_df[which(dates_df[,2] < 10),2] <<- (
    paste("0", dates_df[which(dates_df[,2] < 10), 2], sep = ""))
  date_series <<- as.Date(paste(dates_df[,1],dates_df[,2],"01",sep="-"),
                         format="%Y-%m-%d")
  
  # Logical to set max number of mass campaigns to default
  if (max_modes <= 0) {
    max_modes <<- ceiling((CMC_MDC_max - CMC_MDC_min) / 36)
  }
  
  # Variables derived after loading DHS surveys
  uni_ISO2 <<- unique(all_net_data$ISO2)
  uni_ADM1 <<- unique(all_net_data$ADM1NAME)
  uni_areas <<- unique(all_net_data$area)
  uni_area_ids <<- unique(all_net_data$area_id)
  uni_ADM1_ISO2 <<- unique(paste(all_net_data$ISO2,all_net_data$ADM1NAME,sep=" "))
  N_areas <<- length(uni_area_ids)
}

# Return area data frame
fetch_area_df <- function() {
  matched_area_ids <<- match(uni_area_ids, all_net_data$area_id)
  areas_df <<- data.frame("area" = uni_areas,
                         "area_ID" = uni_area_ids,
                         "ISO2" = substr(uni_areas,1,2),
                         "ADM1" = all_net_data$ADM1NAME[matched_area_ids],
                         "urbanicity" = all_net_data$urbanicity[matched_area_ids],
                         "min_net_age_rec" = rep(NA, N_areas),
                         "max_net_age_rec" = rep(NA, N_areas))
}

update_global_vars_after_new_ids <- function() {
  
  # Store old global variables
  old_uni_ADM1 <<- uni_ADM1
  old_uni_areas <<- uni_areas
  old_uni_area_ids <<- uni_area_ids
  old_uni_ADM1_ISO2 <<- uni_ADM1_ISO2
  old_N_areas <<- N_areas
  
  # Update with new ids
  uni_ISO2 <<- unique(net_data$ISO2)
  uni_ADM1 <<- unique(net_data$ADM1)
  uni_areas <<- unique(net_data$area)
  uni_area_ids <<- unique(net_data$area_id)
  uni_ADM1_ISO2 <<- unique(paste(net_data$ISO2,net_data$ADM1,sep=" "))
  N_areas <<- length(uni_area_ids)
  
}