# usage_access.R
# Usage and access functions

#-------------------------------------------------------------------------------
# Function to return access from DHS data
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

#-------------------------------------------------------------------------------
# Return net data frame (total values)
# N.B. net_data replaces campnets_df in earlier iterations

fetch_net_data <- function() {
  net_data <<- data.frame("area" = rep(uni_areas, each = N_CMC),
                          "area_id" = rep(uni_area_ids, each = N_CMC),
                          "ISO2" = rep(areas_df$ISO2, each = N_CMC),
                          "ADM1" = rep(areas_df$ADM1, each = N_CMC),
                          "urbanicity" = rep(areas_df$urbanicity, each = N_CMC),
                          "CMC" = rep(CMC_series, N_areas),
                          "Date" = rep(date_series, N_areas),
                          "used" = rep(NA, N_areas*N_CMC),
                          "not_used" = rep(NA, N_areas*N_CMC),
                          "access" = rep(NA, N_areas*N_CMC),
                          "no_access" = rep(NA, N_areas*N_CMC))
  return(NULL)
}

#-------------------------------------------------------------------------------
# Append usage and access

append_usage_access <- function(dataset) {
  pc0 <- 0
  for (n in 1:N_areas) {
    ccx <- areas_df$ISO2[n]
    adx <- areas_df$ADM1[n]
    urbx <- areas_df$urbanicity[n]
    
    #subset of net data for admin unit
    if (is.na(urbx)) {
      admin_data <- all_net_data[which(all_net_data$ISO2 == ccx
                                       & all_net_data$ADM1NAME == adx),]
    } else {
      admin_data <- all_net_data[which(all_net_data$ISO2 == ccx
                                       & all_net_data$ADM1NAME == adx
                                       & all_net_data$urbanicity == urbx),]
    }
    
    for (t in 1:N_CMC) {
      i <- t + (n - 1) * N_CMC
      
      admin_data_now <- admin_data[which(admin_data$hv008 == CMC_series[t]),]
      
      # Weighted usage and access
      all_net_used <- admin_data_now[which(admin_data_now$hml20 == 1),]
      num_used <- round(sum(all_net_used$hv005/1e6))
      all_net_acc <- admin_data_now[which(admin_data_now$access == 1),]
      num_acc <- round(sum(all_net_acc$hv005/1e6))
      num_all <- round(sum(admin_data_now$hv005/1e6))
      
      dataset$used[i] <- num_used
      dataset$not_used[i] <- num_all - num_used
      
      dataset$access[i] <- num_acc
      dataset$no_access[i] <- num_all - num_acc
      
    }
    
    pc1 <- round(100 * n / N_areas)
    if (pc1 > pc0) {
      pc0 <- pc1
      print(paste0("Appending usage and access: ", pc0, "% complete"))
    }
    
  }
  
  dataset <- data.frame(dataset,
                        "ADM" = dataset$ADM1,
                        "total" = dataset$used + dataset$not_used)
  dataset <- data.frame(dataset,
                        "prop_used" = dataset$used / dataset$total,
                        "prop_access" = dataset$access / dataset$total)
  
  return(dataset)
}