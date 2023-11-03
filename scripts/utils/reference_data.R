#reference_data.R

fetch_reference_data <- function(input_data) {
  reference_data <<- input_data
  reference_data$ISO2 <<- countrycode(reference_data$ISO3,
                                     origin = 'iso3c',
                                     destination = 'iso2c')
  #reference_data$all_nets <<- reference_data$ITN + reference_data$LLIN
  reference_data$monthly_nets <<- reference_data$LLIN / 12
  reference_data$monthly_nets[is.na(reference_data$monthly_nets)] <<- 0
}

append_reference_nets <- function(dataset) {
  
  if (dim(dataset)[1] != N_areas * N_CMC) {
    print("warning: unexpected dataset dimensions")
  }
  
  dataset$ref_nets <- rep(0, dim(dataset)[1])
  
  for (i in 1:dim(dataset)[1]) {
    cx <- dataset$ISO2[i]
    tx <- dataset$CMC[i]
    yx <- CMC_to_date(tx)[[1]]
    id <- which(reference_data$ISO2 == cx & reference_data$year == yx)
    if (!identical(id, integer(0))) {
      dataset$ref_nets[i] <- reference_data$monthly_nets[id]
    }
  }
  
  return(dataset)
  
}