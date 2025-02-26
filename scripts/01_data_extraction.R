# 01_data_extraction.R
# Script for extracting and cleaning DHS survey data

#-------------------------------------------------------------------------------
# Extract DHS data
# Dependencies in extraction.R
# Surveys were originally extracted from the DHS using the rdhs package
# See 01b_dhs_survey_extraction.R for example code

extracted_surveys <- readRDS(dhs_surveys_path)

#-------------------------------------------------------------------------------
# Clean DHS
# Dependencies in cleaning.R

# Clean and extract usage and access data
all_net_data <- extracted_surveys %>%
  delabel_data() %>%
  standardise_names() %>%
  remove_unknown_sleep_location() %>%
  remove_low_usage() %>%
  generate_unique_ids()

# Get global variables
fetch_init_global_vars()

# Generate area data frame
fetch_area_df()

# CMC limits for minimum and maximum net receipt dates. By default these are
# equal to the bounds of the DHS surveys called but can be changed.
CMC_net_min <- CMC_first
CMC_net_max <- CMC_last

#-------------------------------------------------------------------------------
# Usage and access
# Dependencies in usage_access.R unless otherwise indicated

all_net_data %<>%
  append_CMC_net_obtained() %>%
  simulate_unknown_net_source() %>%
  return_all_access()

# Remove DHS data prior to start of MDCs (input countries, years and months as
# vectors). remove_pre_mdc_dhs() found in cleaning.R
# all_net_data %<>%
#   remove_pre_mdc_dhs("GH", date_to_CMC(year = 2010, month = 1))
# Retracted as of 09/01/2024 - function removing data with NA for CMC_net_obtained

# Fetch net data from (total values)
#fetch_net_data()

# Append access, usage and net source information
net_data <- initialise_net_data() %>%
  append_net_info()













#-------------------------------------------------------------------------------
# Clean DHS data
# Dependencies in ./utils/cleaning.R

extracted_data <- extracted_surveys %>%
  delabel_data %>%
  standardise_names %>%
  generate_unique_ids

#-------------------------------------------------------------------------------
# Extract usage and access data
# Dependencies in ./utils/cleaning.R and ./data_extraction/usage_access.R

# Clean and extract usage and access data
all_net_data <- extracted_data %>%
  remove_unknown_sleep_location %>%
  remove_low_usage %>%
  append_CMC_net_obtained %>%
  simulate_unknown_net_source %>%
  return_all_access

# Append access, usage and net source information
net_data <- initialise_net_data(append_prevalence = TRUE) %>%
  append_net_info()

#-------------------------------------------------------------------------------
# Extract prevalence data
# Dependencies in cleaning.R and prev.R

# Clean and extract usage and access data
all_prev_data <- extracted_data %>%
  delabel_data %>%
  standardise_names %>%

#-------------------------------------------------------------------------------
# Get global variables
fetch_init_global_vars()

# Generate area data frame
fetch_area_df()

# CMC limits for minimum and maximum net receipt dates. By default these are
# equal to the bounds of the DHS surveys called but can be changed.
CMC_net_min <- CMC_first
CMC_net_max <- CMC_last