# 05_post_use_access_fitting.R
# Code for updating retention estimates, linking data to site files, and
# generating nets per capita curve after use and access timeseries fitting
# Bridges appendices 4 and 5 of Glover et al. (2025)
# Functions can be found in ./post_use_access_fitting and ./utils
# Run after 04_use_access_fitting.R

#-------------------------------------------------------------------------------
# Create timestamped data cache for section 05
timestamp_05 <- format(Sys.time(), "%Y%m%d%H%M%S")
cache_05 <- paste0("./data_cache/05_", timestamp_05)
dir.create(cache_05, recursive = TRUE)

#-------------------------------------------------------------------------------
# Calculate retention
# Dependencies in retention.R

first_ret_CMC <- date_to_CMC(first_ret_date[1], first_ret_date[2])
last_ret_CMC <- date_to_CMC(last_ret_date[1], last_ret_date[2])

retention_period <- net_data %>%
  fetch_retention_period(CMCa = first_ret_CMC,
                         CMCb = last_ret_CMC)

#-------------------------------------------------------------------------------
# DEPRICATED - Replaced with new site package
# Link data to foresite
# Dependencies in foresite.R

fs_net_data <- net_data %>%
  append_foresite_names(uni_ISO2) %>%
  create_new_foresite_regions(uni_ISO2) %>%
  append_fs_area_names %>%
  append_fs_area_ids

# Optional: Run for testing
# net_data_pre_fs_link <- net_data

net_data <- fs_net_data

#-------------------------------------------------------------------------------
# Generate nets per capita curve
# Dependencies in npc_stan.R

bv_pred <- stan_npc_fit()

#-------------------------------------------------------------------------------
# Read in net parameters
read_net_params()

#-------------------------------------------------------------------------------
# running 050325_00_05.RData from here (after cache folder initialisation)
#-------------------------------------------------------------------------------
# Reduce memory

# Optional: Run to identify large objects
# sort( sapply(ls(),function(x){object.size(get(x))}))

# Define object names to remove
obj_rm <- c(
  NULL
)

# Define object names to cache and remove
obj_cc <- c(
  "Pbb_u", "Pbb_a", "PC_u", "PC_a", "original_all_net_data",
  "all_net_data_02", "all_net_data_03", "extracted_surveys",
  "MWI_site", "old_net_data", "dataset", "GHA_site", "GH_site",
  "MOZ_site", "SEN_site", "ctry_site", "example_site",
  "BFA_site", "BF_site", "MLI_site"
)

# Define object names to save but keep
obj_sv <- c(
  "bv_gamma", "bv_beta"
  )

# Remove only
handle_objects(obj_rm, remove = TRUE)

# Cache and remove
handle_objects(obj_cc, cache = TRUE, remove = TRUE, cache_path = cache_05)

# Cache but keep
handle_objects(obj_sv, cache = TRUE, remove = FALSE, cache_path = cache_05)

# Force garbage collection and print memory usage after cleanup
cat("Memory usage after cleanup:\n")
gc()

#-------------------------------------------------------------------------------
# Hipercow options

long_sample_ids <- readRDS("./data_public/random_numbers/800_sample_ids.rds")
rnormvals <- readRDS("./data_public/random_numbers/rnormvals.rds")


hipercow::hipercow_configuration()
hipercow::hipercow_init(driver = "dide-windows")
hipercow::windows_authenticate()

hipercow::hipercow_environment_create(
  sources = c(
    "./scripts/post_use_access_fitting/cali_funs.R",
    "./scripts/post_use_access_fitting/EIR_cali.R",
    "./scripts/transmission_model/malsim.R",
    "./scripts/transmission_model/netz_usage_sequential_branch_adapted.R"
  )
)
hipercow::hipercow_provision()
#hipercow_provision(method = "pkgdepends")

# hipercow_provision(
#   requirements = "pkgdepends.txt",
#   repositories = c(
#     mrcide = "https://mrc-ide.r-universe.dev",
#     CRAN = "https://cloud.r-project.org"
#   )
# )

options(hipercow.max_size_local = 1e10)

#-------------------------------------------------------------------------------
# Calibrate EIRs

# Optional: run if testing
# orig_fs_id_link <- fs_id_link

# Optional: debugging code
# debug_dataset <- net_data
# debug_areas <- fs_areas_included
# debug_beta <- bv_beta
# debug_gamma <- bv_gamma
# 
# debug_result <- run_eir_calibration(
#   dataset = debug_dataset,
#   areas_included = debug_areas,
#   bv_beta = debug_beta,
#   bv_gamma = debug_gamma,
#   use_hipercow = FALSE,
#   hiper_debug = TRUE  # Just in case
# )

# Loop over countries
EIR_hipercow <- TRUE # Logical indicator for hipercow
for (i in 1:length(SSA_ISO2)) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  
  EIR_hc_name <- paste("cali", SSA_ISO2[i], "EIR", sep = "_") # Hipercow id name
  
  if (EIR_hipercow) {
    # Calibrate with hipercow
    assign(
      EIR_hc_name,
      net_data %>% run_eir_calibration(areas_included = fs_areas_included,
                                       bv_beta = bv_beta,
                                       bv_gamma = bv_gamma,
                                       use_hipercow = TRUE,
                                       hiper_debug = FALSE,
                                       sim_population = 1e6,
                                       cali_attempts = 30)
    )
    saveRDS(get(EIR_hc_name), paste0(cache_05, "/", EIR_hc_name, "_hc.rds"))
  } else {
    # Calibrate locally
    result_names <- paste("cali", SSA_ISO2[i], "results", sep = "_")
    assign(
      result_names,
      net_data %>% run_eir_calibration(areas_included = fs_areas_included,
                                       bv_beta = bv_beta,
                                       bv_gamma = bv_gamma,
                                       use_hipercow = FALSE)
    )
    saveRDS(get(EIR_hc_name), paste0(EIR_hc_name, ".rds"))
  }
  
}

#-------------------------------------------------------------------------------
# Extract calibration results

if (EIR_hipercow) {
  # All input object names
  all_input_names <- paste("cali", SSA_ISO2, "EIR", sep = "_")
  
  # Check task statuses
  task_statuses <- sapply(all_input_names, function(n) {
    obj <- get(n, envir = .GlobalEnv)
    task_status(obj)
  })
  names(task_statuses) <- SSA_ISO2
  
  # Separate successful tasks and unsuccessful (still running or failed)
  successful_ISOs <- names(task_statuses[task_statuses == "success"])
  unsuccessful_ISOs <- names(task_statuses[task_statuses != "success"])
  
  # Report any not completed
  if (length(unsuccessful_ISOs) > 0) {
    print(paste(
      "The following tasks are not successful (they may still be running):",
      paste(unsuccessful_ISOs, collapse = ", ")
    ))
  }
   
  # Extract, save, and collect results for successful runs
  for (iso in successful_ISOs) {
    input_name <- paste("cali", iso, "EIR", sep = "_")
    output_name <- paste("cali", iso, "results", sep = "_")
    
    assign(
      output_name,
      extract_hipercow_net_runs(get(input_name))
    )
    
    saveRDS(
      get(output_name),
      file = file.path(cache_05, paste0(output_name, ".rds"))
    )
  }
  
  # Combine into single data frame
  if (length(successful_ISOs) > 0) {
    result_names <- paste("cali", successful_ISOs, "results", sep = "_")
    result_list <- mget(result_names, envir = .GlobalEnv)
    EIR_results_df <- do.call(rbind, result_list)
    print("Successfully combined results from completed tasks.")
  } else {
    EIR_results_df <- NULL
    print("No successful results to combine yet.")
  }
  
} else {
  # Get local results
  result_list <- mget(result_names, envir = .GlobalEnv)
  EIR_results_df <- do.call(rbind, result_list)
}

#-------------------------------------------------------------------------------
# Append EIR to fs_id_link

# Get just one EIR_fit per region
EIR_fit_df <- EIR_results_df %>%
  distinct(fs_area_id, CMC, .keep_all = TRUE) %>%
  filter(CMC == CMC_first) %>%
  #select(fs_name_1, EIR_fit)
  select(fs_area_id, EIR_fit)

# Ensure existing EIR_fit is removed BEFORE join
fs_id_link <- fs_id_link %>%
  select(-any_of("EIR_fit")) %>%
  left_join(EIR_fit_df, by = "fs_area_id")

#-------------------------------------------------------------------------------
# Repeat with 8% under-reporting assumption
for (i in 1:length(SSA_ISO2)) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  
  EIR_urep_hc_name <- paste("cali", SSA_ISO2[i], "EIR_urep", sep = "_") # Hipercow id name
  
  if (EIR_hipercow) {
    # Calibrate with hipercow
    assign(
      EIR_urep_hc_name,
      net_data %>% run_eir_calibration(areas_included = fs_areas_included,
                                       bv_beta = bv_beta,
                                       bv_gamma = bv_gamma,
                                       use_hipercow = TRUE,
                                       hiper_debug = FALSE,
                                       sim_population = 1e6,
                                       cali_attempts = 30,
                                       usage_scale_factor = 0.92)
    )
    saveRDS(get(EIR_urep_hc_name), paste0(cache_05, "/", EIR_urep_hc_name, "_hc.rds"))
  } else {
    # Calibrate locally
    urep_result_names <- paste("cali", SSA_ISO2[i], "results", sep = "_")
    assign(
      urep_result_names,
      net_data %>% run_eir_calibration(areas_included = fs_areas_included,
                                       bv_beta = bv_beta,
                                       bv_gamma = bv_gamma,
                                       use_hipercow = FALSE,
                                       usage_scale_factor = 0.92)
    )
    saveRDS(get(EIR_urep_hc_name), paste0(EIR_urep_hc_name, ".rds"))
  }
  
}

# Extract calibration results

if (EIR_hipercow) {
  # All input object names
  all_input_names <- paste("cali", SSA_ISO2, "EIR_urep", sep = "_")
  
  # Check task statuses
  task_statuses <- sapply(all_input_names, function(n) {
    obj <- get(n, envir = .GlobalEnv)
    task_status(obj)
  })
  names(task_statuses) <- SSA_ISO2
  
  # Separate successful tasks and unsuccessful (still running or failed)
  successful_ISOs <- names(task_statuses[task_statuses == "success"])
  unsuccessful_ISOs <- names(task_statuses[task_statuses != "success"])
  
  # Report any not completed
  if (length(unsuccessful_ISOs) > 0) {
    print(paste(
      "The following tasks are not successful (they may still be running):",
      paste(unsuccessful_ISOs, collapse = ", ")
    ))
  }
  
  # Extract, save, and collect results for successful runs
  for (iso in successful_ISOs) {
    input_name <- paste("cali", iso, "EIR_urep", sep = "_")
    output_name <- paste("cali", iso, "results_urep", sep = "_")
    
    assign(
      output_name,
      extract_hipercow_net_runs(get(input_name))
    )
    
    saveRDS(
      get(output_name),
      file = file.path(cache_05, paste0(output_name, ".rds"))
    )
  }
  
  # Combine into single data frame
  if (length(successful_ISOs) > 0) {
    result_names <- paste("cali", successful_ISOs, "results_urep", sep = "_")
    result_list <- mget(result_names, envir = .GlobalEnv)
    EIR_urep_results_df <- do.call(rbind, result_list)
    print("Successfully combined results from completed tasks.")
  } else {
    EIR_urep_results_df <- NULL
    print("No successful results to combine yet.")
  }
  
} else {
  # Get local results
  result_urep_list <- mget(result_names, envir = .GlobalEnv)
  EIR_urep_results_df <- do.call(rbind, result_urep_list)
}

#-------------------------------------------------------------------------------
# Append EIR_urep to fs_id_link

# Get just one EIR_urep_fit per region
EIR_urep_fit_df <- EIR_urep_results_df %>%
  distinct(fs_area_id, CMC, .keep_all = TRUE) %>%
  filter(CMC == CMC_first) %>%
  rename(EIR_urep_fit = EIR_fit) %>%
  #select(fs_name_1, EIR_urep_fit)
  select(fs_area_id, EIR_urep_fit)

# Ensure existing EIR_urep_fit is removed BEFORE join
fs_id_link <- fs_id_link %>%
  select(-any_of("EIR_urep_fit")) %>%
  left_join(EIR_urep_fit_df, by = "fs_area_id")