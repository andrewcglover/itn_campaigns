# misc.R
# Micellaneous functions

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