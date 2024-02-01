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

append_initial_foresite_area_ids <- function(dataset) {
  dataset$fs_area_id <- dataset$area_id
  return(dataset)
}

create_new_foresite_regions <- function(dataset, cc = NULL) {
  
  if (cc %>% is.null) {print("Warning: no countries entered for foresite link")}
  
  if ("MW" %in% cc) {
    
    # Record maximum area id
    max_orig_area_id <- max(dataset$fs_area_id)
    
    # Malawi data subsets
    MW_dataset <- dataset %>% filter(ISO2 == "MW")
    N_dataset <- MW_dataset %>% filter(fs_name_1 == "Northern")
    C_dataset <- MW_dataset %>% filter(fs_name_1 == "Central")
    S_dataset <- MW_dataset %>% filter(fs_name_1 == "Southern")
    
    # Remove existing entries for Malawi
    
    
    # Northern region
    N_dataset <- MW_dataset %>% filter(fs_name_1 == "Northern")
    dataset %<>%
      mutate(fs_name_1 = ifelse(ISO2 == "MW" & fs_name_1 == "Northern",
                                 "Chitipa",
                                 fs_name_1))
    
    
    
    Chitipa
    
    # Subset datasets for Malawi regions
    MW_Nr_dataset <- MW_dataset %>% filter(area == "MW Northern rural")
    MW_Nu_dataset <- MW_dataset %>% filter(area == "MW Northern urban")
    MW_Cr_dataset <- MW_dataset %>% filter(area == "MW Central rural")
    MW_Cu_dataset <- MW_dataset %>% filter(area == "MW Central urban")
    MW_Sr_dataset <- MW_dataset %>% filter(area == "MW Southern rural")
    MW_Su_dataset <- MW_dataset %>% filter(area == "MW Southern urban")
    
    # Copy new
    
  }
  
}