# formatting.R
# Functions for formatting date format

# Convert DHS calendar month code (CMC) to year and month
CMC_to_date <- function(CMC){
  yr <- floor((CMC - 1) / 12) + 1900
  mn <- CMC - ((yr-1900) * 12)
  return(data.frame("year" = yr, "month" = mn))
}

# Convert date (year and month) to CMC format
date_to_CMC <- function(year = 2000, month = 1) {
  CMC <- ((year-1900) * 12) + month
  return(CMC)
}

# Identify oldest and youngest nets
fetch_extreme_nets <- function() {
  extreme_nets <<- data.frame("ISO2" = id_link$ISO2,
                             "ADM1" = id_link$ADM1,
                             "area" = id_link$area,
                             "area_id" = id_link$new_area_id,
                             "min_rec" = rep(NA, N_areas),
                             "max_rec" = rep(NA, N_areas))
  for (i in 1:N_areas) {
    obtained_dates <- all_net_data$CMC_net_obtained[all_net_data$area_id == i]
    extreme_nets$min_rec[i] <<- min(obtained_dates, na.rm = TRUE)
    extreme_nets$max_rec[i] <<- max(obtained_dates, na.rm = TRUE)
  }
}