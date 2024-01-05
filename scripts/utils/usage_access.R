# usage_access.R
# Usage and access functions

#-------------------------------------------------------------------------------
# Append date net obtained in CMC format

append_CMC_net_obtained <- function(dataset) {
  rec_months_since_obt <- dataset$hml4
  rec_months_since_obt[which(rec_months_since_obt > 36)] <- NA
  dataset$CMC_net_obtained <- dataset$hv008 - rec_months_since_obt
  return(dataset)
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
  N_at <- N_areas * N_CMC
  net_data <<- data.frame("area" = rep(uni_areas, each = N_CMC),
                          "area_id" = rep(uni_area_ids, each = N_CMC),
                          "ISO2" = rep(areas_df$ISO2, each = N_CMC),
                          "ADM1" = rep(areas_df$ADM1, each = N_CMC),
                          "urbanicity" = rep(areas_df$urbanicity, each = N_CMC),
                          "CMC" = rep(CMC_series, N_areas),
                          "Date" = rep(date_series, N_areas),
                          "used" = rep(NA, N_at),
                          "not_used" = rep(NA, N_at),
                          "access" = rep(NA, N_at),
                          "no_access" = rep(NA, N_at),
                          "source_rec" = rep(0, N_at),
                          "camp_rec" = rep(0, N_at))
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
      
      # Subset of admin data by survey date
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

#-------------------------------------------------------------------------------
# Append usage, access and net source data

append_net_info <- function(dataset) {
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
      
      # Subset of admin data by survey date
      admin_data_now <- admin_data[which(admin_data$hv008 == CMC_series[t]),]
      
      # Subset of admin data by net receipt date
      #admin_data_obt <- admin_data[which(admin_data$CMC_net_obtained == CMC_series[t]),]
      
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
      
      # Net source
      #source_rec_obt <- length(which(admin_data_obt$hml22 <= 3))
      #camp_rec_obt <- length(which(admin_data_obt$hml22 == 1))
      source_rec_survey <- length(which(admin_data_now$hml22 <= 3))
      camp_rec_survey <- length(which(admin_data_now$hml22 == 1))
      #nets_here <- length(which(admin_data_obt$all_camp == 1))
      dataset$source_rec[i] <- source_rec_survey
      dataset$camp_rec[i] <- camp_rec_survey
      #campnets_df$camp_nets_w_pseudo[i] <- nets_here
      
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