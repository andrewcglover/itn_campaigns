# foresite.R
# Functions for compatibility with foresite package

# Fuction to append foresite names
append_foresite_names <- function(dataset, cc = NULL) {
  # Note: no naming differences for GH
  
  if (cc %>% is.null) {print("Warning: no countries entered for foresite link")}
  
  dataset$fs_name_1 <- dataset$ADM1
  
  if ("BF" %in% cc) {
    dataset$fs_name_1[which(dataset$fs_name_1 == "Centre Est")] <- "Centre-Est"
    dataset$fs_name_1[which(dataset$fs_name_1 == "Centre Nord")] <- "Centre-Nord"
    dataset$fs_name_1[which(dataset$fs_name_1 == "Centre Ouest")] <- "Centre-Ouest"
    dataset$fs_name_1[which(dataset$fs_name_1 == "Centre Sud")] <- "Centre-Sud"
    dataset$fs_name_1[which(dataset$fs_name_1 == "Hauts Bassins")] <- "Hauts-Bassins"
    dataset$fs_name_1[which(dataset$fs_name_1 == "Plateau Central")] <- "Plateau Central"
    dataset$fs_name_1[which(dataset$fs_name_1 == "Sud Ouest")] <- "Sud-Ouest"
  }
  if ("MW" %in% cc) {
    print(paste0("Warning: DHS treats the 3 regions of Malawi as adm1, while ",
                 "foresite uses the 28 districts"))
  }
  if ("ML" %in% cc) {
    dataset$fs_name_1[which(dataset$fs_name_1 == "Segou")] <- "Ségou"
    dataset$fs_name_1[which(dataset$fs_name_1 == "Tombouctou")] <- "Timbuktu"
  }
  if ("ML" %in% cc) {
    dataset$fs_name_1[which(dataset$fs_name_1 == "Segou")] <- "Ségou"
    dataset$fs_name_1[which(dataset$fs_name_1 == "Tombouctou")] <- "Timbuktu"
  }
  if ("MZ" %in% cc) {
    dataset$fs_name_1[which(dataset$fs_name_1 == "Niassa")] <- "Nassa"
  }
  if ("SN" %in% cc) {
    dataset$fs_name_1[which(dataset$fs_name_1 == "Kedougou")] <- "Kédougou"
    dataset$fs_name_1[which(dataset$fs_name_1 == "Saint Louis")] <- "Saint-Louis"
    dataset$fs_name_1[which(dataset$fs_name_1 == "Sedhiou")] <- "Sédhiou"
    dataset$fs_name_1[which(dataset$fs_name_1 == "Thies")] <- "Thiès"
  }
  
  return(dataset)
  
}

foresite_subADM_rows <- function(dataset, topADMdata = NULL, subADM = NULL) {
  # Number of sub admin units
  Nsub <- length(subADM)
  # Create new entries for sub admin units
  for (i in 1:Nsub) {
    subADMdata <- topADMdata
    subADMdata$fs_name_1 <- rep(subADM[i], dim(topADMdata)[1])
    dataset %<>% rbind(subADMdata)
  }
  return(dataset)
}

create_new_foresite_regions <- function(dataset, cc = NULL) {
  
  if (cc %>% is.null) {print("Warning: no countries entered for foresite link")}
  
  if ("MW" %in% cc) {
    
    # Malawi districts (in foresite format)
    Northern_districts <- c("Chitipa",
                            "Karonga",
                            "Likoma",
                            "Mzimba",
                            "Nkhata Bay",
                            "Rumphi")
    Central_districts <- c("Dedza",
                           "Dowa",
                           "Kasungu",
                           "Lilongwe",
                           "Mchinji",
                           "Nkhotakota",
                           "Ntcheu",
                           "Ntchisi",
                           "Salima")
    Southern_districts <- c("Balaka",
                            "Blantyre",
                            "Chikwawa",
                            "Chiradzulu",
                            "Machinga",
                            "Mangochi",
                            "Mulanje",
                            "Nsanje",
                            "Thyolo",
                            "Phalombe",
                            "Zomba",
                            "Neno")
    
    # Malawi data subsets
    MW_dataset <- dataset %>% filter(ISO2 == "MW")
    N_dataset <- MW_dataset %>% filter(fs_name_1 == "Northern")
    C_dataset <- MW_dataset %>% filter(fs_name_1 == "Central")
    S_dataset <- MW_dataset %>% filter(fs_name_1 == "Southern")
    
    # Remove existing entries for Malawi from main dataset
    dataset %<>% filter(!(ISO2 == "MW"))
    
    # Append district rows
    dataset %<>%
      foresite_subADM_rows(N_dataset, Northern_districts) %>%
      foresite_subADM_rows(C_dataset, Central_districts) %>%
      foresite_subADM_rows(S_dataset, Southern_districts)
    
  }
  
  return(dataset)
  
}

# Create new foresite area names
append_fs_area_names <- function(dataset) {
  dataset$fs_area <- paste(dataset$ISO2,
                           dataset$fs_name_1,
                           dataset$urbanicity,
                           sep = " ")
  return(dataset)
}

# Create new foresite ids
append_fs_area_ids <- function(dataset) {
  unique_fs_areas_included <- unique(dataset$fs_area)
  dataset$fs_area_id <- match(dataset$fs_area, unique_fs_areas_included)
  fs_id_link <<- unique(data.frame("ISO2" = dataset$ISO2,
                                   "ADM1" = dataset$ADM1,
                                   "area" = dataset$area,
                                   "fs_name_1" = dataset$fs_name_1,
                                   "fs_area" = dataset$fs_area,
                                   "CTRY" = dataset$CTRY,
                                   "old_area_id" = dataset$old_area_id,
                                   "new_area_id" = dataset$area_id,
                                   "fs_area_id" = dataset$fs_area_id))
  N_fs_areas <<- max(dataset$fs_area_id)
  return(dataset)
}