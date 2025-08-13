# 05_post_use_access_fitting.R
# Code for updating retention estimates, linking data to site files, and
# generating nets per capita curve after use and access timeseries fitting

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

#net_data_pre_fs_link <- net_data

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
# identify large objects
sort( sapply(ls(),function(x){object.size(get(x))}))

# Define object names to remove
obj_rm <- c(
  
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

# # Print memory usage before cleanup
# cat("Memory usage before cleanup:\n")
# gc()
# 
# # Loop through each object to remove
# for (name in obj_rm) {
#   if (exists(name, envir = .GlobalEnv)) {
#     rm(list = name, envir = .GlobalEnv)
#   } else {
#     cat("Object does not exist and was skipped:", name, "\n")
#   }
# }
# 
# # Loop through each object to cache and remove
# for (name in obj_names) {
#   if (exists(name, envir = .GlobalEnv)) {
#     saveRDS(get(name), file = file.path(cache_05, paste0("/", name, ".rds")))
#     rm(list = name, envir = .GlobalEnv)
#   } else {
#     cat("Object does not exist and was skipped:", name, "\n")
#   }
# }
# 
# # Loop through each object to save but retain
# for (name in obj_sv) {
#   if (exists(name, envir = .GlobalEnv)) {
#     saveRDS(get(name), file = file.path(cache_05, paste0("/", name, ".rds")))
#   } else {
#     cat("Object does not exist and was skipped:", name, "\n")
#   }
# }

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
#orig_fs_id_link <- fs_id_link

# EIR_BF_df <- readRDS("EIR_BF_df.rds")
# EIR_GH_df <- readRDS("EIR_GH_df.rds")
# EIR_ML_df <- readRDS("EIR_ML_df.rds")
# EIR_MZ_df <- readRDS("EIR_MZ_df.rds")
# EIR_MW_df <- readRDS("EIR_MW_df.rds")
# EIR_SN_df <- readRDS("EIR_SN_df.rds")
# 
# EIR_GH_df2 <- EIR_GH_df[!duplicated(cbind(EIR_GH_df$CMC,EIR_GH_df$fs_area_id)), ]

debug_dataset <- net_data
debug_areas <- fs_areas_included
debug_beta <- bv_beta
debug_gamma <- bv_gamma

debug_result <- run_eir_calibration(
  dataset = debug_dataset,
  areas_included = debug_areas,
  bv_beta = debug_beta,
  bv_gamma = debug_gamma,
  use_hipercow = FALSE,
  hiper_debug = TRUE  # Just in case
)


for (i in 1:length(SSA_ISO2)) {
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  
  #EIR_df_name <- paste("EIR10", SSA_ISO2[i], "fit", sep = "_")
  #EIR_df_name <- paste("EIR", SSA_ISO2[i], "fit", sep = "_")
  EIR_df_name <- paste("cali", SSA_ISO2[i], "EIR", sep = "_")
  #EIR_df_name <- "testpar10_id"
  
  EIR_hipercow <- TRUE
  
  if (EIR_hipercow) {
    assign(
      EIR_df_name,
      net_data %>% run_eir_calibration(areas_included = fs_areas_included,
                                       bv_beta = bv_beta,
                                       bv_gamma = bv_gamma,
                                       use_hipercow = TRUE,
                                       hiper_debug = FALSE,
                                       sim_population = 1e6,
                                       cali_attempts = 30)
    )
    # saveRDS(get(EIR_df_name), paste0("./data_cache/", EIR_df_name, "_hc.rds"))
    saveRDS(get(EIR_df_name), paste0(cache_05, "/", EIR_df_name, "_hc.rds"))
  } else {
    assign(
      EIR_df_name,
      net_data %>% run_eir_calibration(areas_included = fs_areas_included,
                                       bv_beta = bv_beta,
                                       bv_gamma = bv_gamma,
                                       use_hipercow = FALSE)
    )
    saveRDS(get(EIR_df_name), paste0(EIR_df_name, ".rds"))
  }
  
}

# for (i in 1:length(SSA_ISO2)) {
#   EIR_df_name <- paste("cali", SSA_ISO2[i], "EIR", sep = "_")
#   assign(
#     EIR_df_name,
#     #readRDS(paste0("./data_cache/", EIR_df_name, "_hc.rds"))
#     readRDS(paste0(cache_05, "/", EIR_df_name, "_hc.rds"))
#   )
# }

#-------------------------------------------------------------------------------
# Extract hipercow runs

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

#-------------------------------------------------------------------------------
# Append EIR to fs_id_link

# Get just one EIR_fit per region
EIR_fit_df <- EIR_results_df %>%
  filter(CMC == CMC_first) %>%
  select(fs_area_id, EIR_fit)

# Ensure existing EIR_fit is removed BEFORE join
fs_id_link <- fs_id_link %>%
  select(-any_of("EIR_fit")) %>%
  left_join(EIR_fit_df, by = "fs_area_id")





# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# all_input_names <- paste0("cali_", SSA_ISO2, "_EIR")
# 
# # Check if task_status(...) == "success" for all
# all_success <- all(sapply(all_input_names, function(n) {
#   obj <- get(n, envir = .GlobalEnv)
#   task_status(obj) == "success"
# }))
# 
# # Optional: print result
# if (all_success) {
#   print("all tasks completed successfully.")
#   for (iso in SSA_ISO2) {
#     # Build the input object name (e.g. cali_BF_EIR)
#     input_name <- paste("cali", iso, "EIR", sep = "_")
#     
#     # Build the output object name (e.g. cali_BF_results)
#     EIR_results_name <- paste("cali", iso, "results", sep = "_")
#     
#     # Call the function with the dynamic object and assign result
#     assign(
#       EIR_results_name,
#       extract_hipercow_net_runs(get(input_name))
#     )
#     
#     # Save result to cache
#     saveRDS(
#       get(EIR_results_name),
#       file = file.path(cache_05, paste0(EIR_results_name, ".rds"))
#     )
#   }
#   
#   # Combine into a single data frame
#   EIR_results_list <- paste("cali", SSA_ISO2, "results", sep = "_") %>%
#     mget(envir = .GlobalEnv)
#   EIR_results_df <- do.call(rbind, EIR_results_list)
# } else {
#   print("some tasks have failed or are not yet finished.")
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# BF_test <- extract_hipercow_net_runs(cali_BF_EIR)
# GH_test <- extract_hipercow_net_runs(cali_GH_EIR)
# 
# for (i in 1:length(test)) {
#   #print(paste0("dim(test[[", i, "]] = ", dim(test[[i]])))
#   print(unique(test[[i]]$fs_area))
#   print(dim(test[[i]]))
# }
# 
# for (i in 1:length(test)) {
#   test_area <- unique(test[[i]]$fs_area)
#   test_df <- net_data %>% filter(fs_area == test_area)
#   print(test_area)
#   print(dim(test_df))
# }
# 
# 
# ggplot(data = GH_test %>% dplyr::filter(CMC == CMC_first)) +
#   geom_col(aes(x = fs_area, y = EIR_fit))
# 
# 
# test_df <- extract_hipercow_net_runs(test_EIR_BF_fit)
# 
# ggplot(data = GH_test) +
#   geom_line(aes(x = CMC, y = pfpr_2_10_fit)) +
#   geom_step(aes(x = CMC, y = pfpr_map)) +
#   facet_grid(cols = vars(urbanicity), rows = vars(fs_name_1))
# 
# test <- dplyr::left_join(EIR_MW_df, fs_id_link)
# 
# ggplot(test,
#        aes(x = CMC, y = pfpr_dhs, ymin = pfpr_dhs_lo, ymax = pfpr_dhs_hi,
#            colour = ADM1, shape = urbanicity)) +
#   geom_errorbar() +
#   geom_point() +
#   facet_wrap(vars(ADM1))
# 
# 
# 
# #-------------------------------------------------------------------------------
# # Calibrate EIRs
# #orig_fs_id_link <- fs_id_link
# 
# EIR_df <- net_data%>% run_eir_calibration(
#   areas_included = fs_id_link$fs_area,
#   bv_beta = bv_beta,
#   bv_gamma = bv_gamma
# )
# 
# EIR_BF_df <- EIR_BF_df[!duplicated(cbind(EIR_BF_df$CMC,EIR_BF_df$fs_area_id)), ]
# EIR_ML_df <- EIR_ML_df[!duplicated(cbind(EIR_ML_df$CMC,EIR_ML_df$fs_area_id)), ]
# EIR_MZ_df <- EIR_MZ_df[!duplicated(cbind(EIR_MZ_df$CMC,EIR_MZ_df$fs_area_id)), ]
# EIR_MW_df <- EIR_MW_df[!duplicated(cbind(EIR_MW_df$CMC,EIR_MW_df$fs_area_id)), ]
# 
# #eir_net_data_old <- eir_net_data
# eir_net_data <- dplyr::left_join(
#   net_data,
#   rbind.data.frame(
#     EIR_BF_df,
#     EIR_GH_df
#     #EIR_ML_df,
#     #EIR_MZ_df,
#     #EIR_MW_df
#     ),#, EIR_MZ_df, EIR_SN_df),
#   by = c("ISO2", "fs_area", "fs_name_1", "urbanicity", "fs_area_id", "CMC")
# )
# 
# #EIR_only_df_old <- EIR_only_df
# 
# 
# EIR_only_df <- rbind.data.frame(
#   EIR_BF_df,
#   EIR_GH_df
#   #EIR_ML_df,
#   #EIR_MZ_df,
#   #EIR_MW_df
#   ) %>%
#   dplyr::filter(CMC == CMC_first)
# 
# fs_id_link <- dplyr::left_join(
#   orig_fs_id_link,
#   EIR_only_df,
#   by = c("ISO2", "fs_area", "fs_name_1", "urbanicity", "fs_area_id")
# )
# #EIR_BF_df_old <- EIR_BF_df
# #fs_id_link_old <- fs_id_link
# # Routine ITN channels only
# #bv_beta_backup <- bv_beta
# #bv_gamma_backup <- bv_gamma
# for (i in 3:3) {
#   
#   # Sub-set areas by country
#   fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
#   
#   EIR_df_name <- paste("EIR", SSA_ISO2[i], "df", sep = "_")
#   
#   assign(
#     #paste("EIR", SSA_ISO2[i], "df", sep = "_"),
#     EIR_df_name,
#     net_data %>% run_eir_calibration(areas_included = fs_areas_included,
#                                      bv_beta = bv_beta,
#                                      bv_gamma = bv_gamma)
#   )
#   
#   saveRDS(get(EIR_df_name), paste0(EIR_df_name, ".rds"))
#   
# }
# 
# saveRDS(EIR_GH_df, "EIR_GH_df.rds")
# 
# 
# 
# fs_areas_included <- fs_id_link$fs_area[
#   which(
#     fs_id_link$ISO2 == "BF"
#     )
#   ]
# EIR_id_name <- "EIR_BF_id"
# assign(
#   #paste("EIR", SSA_ISO2[i], "df", sep = "_"),
#   EIR_id_name,
#   net_data %>% run_eir_calibration(areas_included = fs_areas_included,
#                                    bv_beta = bv_beta,
#                                    bv_gamma = bv_gamma,
#                                    use_hipercow = TRUE,
#                                    hiper_debug = TRUE)
# )
# 
# # fs_areas_included <- fs_id_link$fs_area[
# #   which(
# #     fs_id_link$ISO2 == "BF" & fs_id_link$urbanicity == "rural"
# #     )
# #   ]
# # EIR_df_name <- "EIR_BF_rural_id"
# # assign(
# #   #paste("EIR", SSA_ISO2[i], "df", sep = "_"),
# #   EIR_df_name,
# #   net_data %>% run_eir_calibration(areas_included = fs_areas_included,
# #                                    bv_beta = bv_beta,
# #                                    bv_gamma = bv_gamma,
# #                                    use_hipercow = TRUE)
# # )
# # 
# # fs_areas_included <- fs_id_link$fs_area[
# #   which(
# #     fs_id_link$ISO2 == "BF" & fs_id_link$urbanicity == "urban"
# #   )
# # ]
# # EIR_df_name <- "EIR_BF_urban_id"
# # assign(
# #   #paste("EIR", SSA_ISO2[i], "df", sep = "_"),
# #   EIR_df_name,
# #   net_data %>% run_eir_calibration(areas_included = fs_areas_included,
# #                                    bv_beta = bv_beta,
# #                                    bv_gamma = bv_gamma,
# #                                    use_hipercow = TRUE)
# # )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ggplot(
#   data = EIR_only_df,
#   aes(
#     x = fs_name_1,
#     y = EIR_fit,
#     fill = urbanicity
#   )
# ) +
#   geom_col(position = 'dodge')