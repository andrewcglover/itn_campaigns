#mdc_rounds.R

# Define last and penultimate rounds
append_mdc_rounds <- function(dataset) {
  for (n in 1:N_areas) {
    most_recent_MDC <- NA
    secondmost_recent_MDC <- NA
    round_num <- NA
    for (t in 1:N_CMC) {
      i <- t + (n - 1) * N_CMC
      if (dataset$mdc[i]) {
        secondmost_recent_MDC <- most_recent_MDC
        most_recent_MDC <- CMC_series[t]
        if (is.na(round_num)) {
          round_num <- 1
        } else {
          round_num <- round_num + 1
        }
      }
      # append round number
      dataset$MDC_round[i] <- round_num
      # append months since previous round
      dataset$months_post_mdc[i] <- min(max_m,
                                        dataset$CMC[i] - most_recent_MDC,
                                        na.rm = TRUE)
      # append months since penultimate round
      dataset$months_post_prior_mdc[i] <- min(max_m,
                                              dataset$CMC[i] - secondmost_recent_MDC,
                                              na.rm = TRUE)
    }
  }
  # Assign NA rounds to zero
  dataset$MDC_round[is.na(dataset$MDC_round)] <- 0
  # Add one to all round ids
  #turned off 030124
  #dataset$MDC_round <- dataset$MDC_round + 1
  
  print(paste0("Assigning rounds complete."))
  return(dataset)
}

# Check unique areas included
unique_areas_included_check <- function() {
  unique_areas_included <<- unique(net_data$area)
}

# Generate MDC round matrices
generate_MDC_round_matrices <- function(max_tau = 12) {
  
  # Generate MDC round data frame
  MDC_rounds <<- data.frame("ISO2" = net_data$ISO2,
                            "ADM1" = net_data$ADM1,
                            "area" = net_data$area,
                            "round" = net_data$MDC_round,
                            "CMC" = net_data$CMC,
                            "tau" = net_data$mdc_tau)
  MDC_rounds <<- MDC_rounds[which(net_data$mdc),]
  MDC_rounds$area_id <<- match(MDC_rounds$area, unique_areas_included)
  max_rounds <<- max(MDC_rounds$round) + 1
  #MDC_rounds <<- data.frame(MDC_rounds,
  #                          "round_tau" = rep(6, dim(MDC_rounds)[1]))
  
  # Generate MDC matrices
  MDC_tau_matrix <<- matrix(-1, nrow = N_areas, ncol = max_rounds)
  MDC_matrix <<- matrix(-9999, nrow = N_areas, ncol = max_rounds)
  MDC_matrix[,1] <<- rep(date_to_CMC(2000,1), N_areas)
  for (i in 1:dim(MDC_rounds)[1]) {
    x <- MDC_rounds$area_id[i]
    y <- MDC_rounds$round[i] + 1
    MDC_matrix[x,y] <<- MDC_rounds$CMC[i]
    MDC_tau_matrix[x,y] <<- MDC_rounds$tau[i]
  }
  
  #MDC_tau_matrix <<- 2 * MDC_tau_matrix / mean(MDC_rounds$round_tau)
  MDC_tau_matrix[MDC_tau_matrix < 0] <<- max_tau
  MDC_tau_matrix[MDC_tau_matrix > max_tau] <<- max_tau
}