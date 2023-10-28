# cleaning.R
# Data cleaning functions

clean_net_data <- function(extract_surveys) {
  
  # Remove labels
  all_net_data <- plyr::ldply(extract_surveys, data.frame)
  labelled::remove_val_labels(all_net_data)
  
  # Remove special characters and capitalize admin regions
  all_net_data$ADM1NAME <- sub("-"," ",all_net_data$ADM1NAME)
  all_net_data$ADM1NAME <- stringr::str_to_title(stri_trans_general(all_net_data$ADM1NAME,
                                                                    "Latin-ASCII"))
  # Change all nulls to NA
  all_net_data$ADM1NAME[which(all_net_data$ADM1NAME == "Null")] <- NA
  
  # Remove NAs
  all_net_data <- all_net_data[which(all_net_data$ADM1NAME != "NA"),]
  all_net_data <- all_net_data[which(!is.na(all_net_data$ADM1NAME)),]
  
  # Correct for alternative spelling
  all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Zuguinchor")]<-"Ziguinchor"
  all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Maputo Cidade")]<-"Maputo City"
  all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Maputo Province")]<-"Maputo"
  all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Maputo Provincia")]<-"Maputo"
  all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Toumbouctou")]<-"Tombouctou"
  
  all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Boucle De Mouhoun")]<-"Boucle Du Mouhoun"
  all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="Hauts Basins")]<-"Hauts Bassins"
  
  all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="North")]<-"Northern"
  all_net_data$ADM1NAME[which(all_net_data$ADM1NAME=="South")]<-"Southern"
  
  #Record countries and admin 1 locations with DHS survey data
  all_net_data$ISO2 <- substr(all_net_data$SurveyId,1,2)

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
  
  #Create household ids
  all_net_data$hhid <- paste(all_net_data$.id, all_net_data$hv001,
                             all_net_data$hv002, sep = "_")
  
  #Create net ids
  all_net_data$netid <- paste(all_net_data$hhid, all_net_data$hmlidx, sep = "_")
  all_net_data$netid[which(is.na(all_net_data$hmlidx))] <- NA
  
  #Print to state initial data cleaning complete
  print("initial data cleaning: 100% complete")
  
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
  
  #Record total entries
  dim_net_data <- dim(all_net_data)
  N_net_data <- dim_net_data[1]
  
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
  
  #return cleaned data frame
  return(all_net_data)
}