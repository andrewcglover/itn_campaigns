# dhs_prevalence.R


fetch_prev_data <- function() {
  N_at <- N_areas * N_CMC
  net_data <- data.frame("area" = rep(uni_areas, each = N_CMC),
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
  return(net_data)
}

#-------------------------------------------------------------------------------
# Append usage, access and net source data

append_prev_info <- function(dataset) {
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
      
      # Subset of 6 - 59 month olds
      admin_kids <- admin_data_now[
        which(admin_data_now$hml16a >= 6 & admin_data_now$hml16a <= 59),
      ]
      
      # Weighted prevlaence
      RDT_all <- admin_kids[
        which(admin_kids$hml35 >= 0 & admin_kids$hml35 <= 1),
      ]
      RDT_positive <- admin_kids[which(admin_kids$hml35 == 1),]
      num_rdt_positive <- sum(RDT_positive$hv005/1e6) 
      num_rdt <- sum(RDT_all$hv005/1e6)
      rdt_prev <- num_rdt_positive / num_rdt
      
      
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