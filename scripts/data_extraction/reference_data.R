#reference_data.R

# Function to extract MAP and AMP reference data

extract_reference_data <- function(input_data) {
  reference_data <- input_data
  reference_data$ISO2 <- countrycode::countrycode(reference_data$ISO3,
                                     origin = 'iso3c',
                                     destination = 'iso2c')
  reference_data$annual_nets <- rowMeans(
    cbind(reference_data$MAP_LLIN, reference_data$AMP_TOT),
    na.rm = TRUE
    )
  reference_data$monthly_nets <- reference_data$annual_nets / 12
  return(reference_data)
}

#-------------------------------------------------------------------------------

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