#cleaning.R

### Data cleaning functions ###

clean_net_data <- function(extract_surveys) {
  
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
  
  #return cleaned data frame
  return(all_net_data)
}