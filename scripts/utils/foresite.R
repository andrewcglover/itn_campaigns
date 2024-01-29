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