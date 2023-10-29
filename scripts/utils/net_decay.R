# net_decay.R

#-------------------------------------------------------------------------------

# Account for DHS weightings before hierarchical net decay model
# This function creates a new crop of nets 
net_weighting_fun <- function(all_net_data, CMC_net_min, CMC_net_max,
                              access = FALSE) {
  
  # Create net data frame
  nets <- data.frame("survey_id" = all_net_data$SurveyId,
                     "ISO2" = all_net_data$ISO2,
                     "ADM1" = all_net_data$ADM1NAME,
                     "urbanicity" = all_net_data$urbanicity,
                     "area" = all_net_data$area,
                     "area_id" = all_net_data$area_ID,
                     "cluster" = all_net_data$hv001,
                     "household" = all_net_data$hv002,
                     "household_w" = all_net_data$hv005,
                     "CMC_interview" = all_net_data$hv008,
                     "hh_net_id" = all_net_data$hmlidx,
                     "months_ago_obtained" = all_net_data$hml4,
                     "CMC_obtained" = all_net_data$CMC_net_obtained)
                     #"CMC_obtained" = all_net_data$raw_CMC_net_obtained)
  
  # Logical check for access or usage weights
  if (access) {
    nets <- data.frame(nets,
                       "anyone_access" = all_net_data$access)
    nets <- nets[which(nets$anyone_access == 1),]
  } else {
    nets <- data.frame(nets,
                       "anyone_used" = all_net_data$hml21)
    nets <- nets[which(nets$anyone_used == 1),]
  }
  
  # Remove nets received outside of bounds or without CMC recorded
  nets <- nets[which(!is.na(nets$CMC_obtained)),]
  nets$net_id <- paste(nets$area_id,
                       nets$cluster,
                       nets$household,
                       nets$hh_net_id,
                       sep = ".")
  nets <- nets[!duplicated(nets$net_id),]
  nets <- nets[which(nets$CMC_obtained >= CMC_net_min),]
  nets <- nets[which(nets$CMC_obtained <= CMC_net_max),]
  nets <- nets[which(nets$months_ago_obtained <= 96),] #Error codes >=97
  nets$CMC_obt_int_area_id <- paste(nets$CMC_obtained,
                                    nets$CMC_interview,
                                    nets$area_id,
                                    sep = ".")
  
  # Unique identifiers
  uni_CMC_obt_int_area_id <- unique(nets$CMC_obt_int_area_id)
  N_ucoia <- length(uni_CMC_obt_int_area_id)
  
  # Usage vs access subset for unique nets
  nets_totals <- nets[!duplicated(nets$CMC_obt_int_area_id),]
  if (access) {
    nets_totals = subset(nets_totals,
                         select = -c(cluster,
                                     household,
                                     household_w,
                                     hh_net_id,
                                     anyone_access,
                                     net_id))
  } else {
    nets_totals = subset(nets_totals,
                         select = -c(cluster,
                                     household,
                                     household_w,
                                     hh_net_id,
                                     anyone_used,
                                     net_id))
  }
  nets_totals$weighted_counts <- rep(NA, N_ucoia)
  
  # Calculate weighted counts
  for (i in 1:N_ucoia) {
    matched_ids <- which(nets$CMC_obt_int_area_id == nets_totals$CMC_obt_int_area_id[i])
    weighted_counts_here <- round (sum(nets$household_w[matched_ids]) / 1e6)
    nets_totals$weighted_counts[i] <- weighted_counts_here
  }
  
  # Only include non-zero weighted counts
  nets_totals <- nets_totals[which(nets_totals$weighted_counts > 0),]
  
  M <- dim(nets_totals)[1]
  N <- sum(nets_totals$weighted_counts)
  
  nets_weighted <- data.frame("ISO2" = rep(NA, N),
                              "ADM1" = rep(NA, N),
                              "urbanicity" = rep(NA, N),
                              "area" = rep(NA, N),
                              "area_id" = rep(NA, N),
                              "net_obtained" = rep(NA, N),
                              "months_since_obtained" = rep(NA, N),
                              "CTRY" = rep(NA, N))
  
  j <- 1
  for (i in 1:M) {
    counts <- nets_totals$weighted_counts[i]
    k <- j + counts - 1
    nets_weighted$ISO2[j:k] <- nets_totals$ISO2[i]
    nets_weighted$ADM1[j:k] <- nets_totals$ADM1[i]
    nets_weighted$urbanicity[j:k] <- nets_totals$urbanicity[i]
    nets_weighted$area[j:k] <- nets_totals$area[i]
    nets_weighted$area_id[j:k] <- nets_totals$area_id[i]
    nets_weighted$net_obtained[j:k] <- nets_totals$CMC_obtained[i]
    nets_weighted$months_since_obtained[j:k] <- nets_totals$months_ago_obtained[i]
    j <- j + counts
  }
  
  return(nets_weighted)
  
}

stan_decay_fit <- function(nets_weighted, adm_net_link) {
  
  N_a <- length(unique(nets_weighted$area))
  N_c <- N_ISO2
  
  net_decay_dat <- list(N = dim(nets_weighted)[1],
                        N_a = N_a,
                        N_c = N_c,
                        a = nets_weighted$area_id,
                        c = nets_weighted$CTRY,
                        cc = adm_net_link$CTRY,
                        m = nets_weighted$months_since_obtained)
  
  net_decay_fit <- stan('exp_model_d2.stan',
                        data = net_decay_dat,
                        iter = 1000,
                        warmup = 500,
                        chains = 4,
                        #init_r = 1e-2,
                        control = list(adapt_delta = 0.95))
  
  net_decay_samples <- extract(net_decay_fit)
  
  return(net_decay_samples)
  
}