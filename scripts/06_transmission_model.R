# 06_transmission_model.R

#-------------------------------------------------------------------------------
# Malaria Simulation

#-------------------------------------------------------------------------------
# Create timestamped data cache for section 05

# Create timestamped data cache for section 05
# timestamp_06 <- format(Sys.time(), "%Y%m%d%H%M%S")
cache_06 <- paste0("./data_cache/06_", timestamp_05) # use 05 timestamp
dir.create(cache_06, recursive = TRUE)

#-------------------------------------------------------------------------------
# Reduce memory

# Optional: Run to identify large objects
# sort( sapply(ls(),function(x){object.size(get(x))}))

# Define object names to remove
obj_rm <- c(
  "x", "y", "row_equal", "num_different_rows", "diff_cols", # testing objs
  "param_list", "adm_site", "BF_test", "GH_test"
)

# Define object names to cache and remove
obj_cc <- c(
  "net_data_pre_fs_link", "eir_net_data", "new_net_data", "all_net_data_03a",
  "net_data_03a", "net_data_prev"
)

# Define object names to save but keep
obj_sv <- c(
  "bv_gamma", "bv_beta"
)

# Remove only
handle_objects(obj_rm, remove = TRUE)

# Cache and remove
handle_objects(obj_cc, cache = TRUE, remove = TRUE, cache_path = cache_06)

# Cache but keep
handle_objects(obj_sv, cache = TRUE, remove = FALSE, cache_path = cache_06)

# Force garbage collection and print memory usage after cleanup
cat("Memory usage after cleanup:\n")
gc()

#-------------------------------------------------------------------------------

# Convert projection times to months
mass_int_mn <- mass_int_yr * 12
projection_window_mn <- projection_window_yr * 12
projection_window_dy <- projection_window_yr * 365

fs_areas_included <- unique(fs_id_link$fs_area)
fs_excluded <- c()
fs_areas_included <- fs_areas_included[! fs_areas_included %in% fs_excluded]

# Number of samples
N_samples <- dim(P_u)[1]

# Monthly offset for future mass campaigns
# month_offset <- sample.int(13, N_reps, replace = TRUE) - 7

# Create sample ids
# Uncomment lines below to create sample ids
# long_sample_ids <- sample.int(N_samples, 10000 , replace = TRUE)
# saveRDS(long_sample_ids, "./data_public/random_numbers/800_sample_ids.rds")
long_sample_ids <- readRDS("./data_public/random_numbers/800_sample_ids.rds")
rnormvals <- readRDS("./data_public/random_numbers/rnormvals.rds")


hipercow::hipercow_configuration()
hipercow::hipercow_init(driver = "dide-windows")
hipercow::windows_authenticate()
hipercow_environment_create(
  sources = c(
    #"./scripts/post_use_access_fitting/cali_funs.R",
    #"./scripts/post_use_access_fitting/EIR_cali.R",
    "./scripts/transmission_model/malsim.R",
    "./scripts/transmission_model/netz_usage_sequential_branch_adapted.R"
  )
)
hipercow_provision()
#hipercow_provision(method = "pkgdepends")

options(hipercow.max_size_local = 1e11)

# Previous estimates (TGF Grant Cycle 6)
# only_cost <- 1.95
# pbo_cost <- 2.54
# pyrrole_cost <- 2.56
# dist_cost <- 2.75 # estimated

# New estimates (TGF grant cycle 7)
only_cost <- 1.85
pbo_cost <- 2.14
pyrrole_cost <- 2.56
dist_cost <- 2



only_total_cost <- dist_cost + only_cost
pbo_total_cost <- dist_cost + pbo_cost
pyrrole_total_cost <- dist_cost + pyrrole_cost

scaled_pbo_nets_equiv_only <- only_total_cost / pbo_total_cost
scaled_pyrrole_nets_equiv_only <- only_total_cost / pyrrole_total_cost

# Simulation parameters
mass_int_yr <- c(2,3)     # Mass campaign intervals
sim_reps <- 100
sim_cores <- 16          # Takes precedence over areas per core
sim_areas_per_core <- 1   # Requires sim_cores = 0, otherwise ignored
sim_id <- timestamp_05
# budget_val <- c(100, 75, 50)


# past_sim_df <- data.frame("id.iso" = NULL,
#                              "id.name" = NULL,
#                              "id.string" = NULL)

urep_past_sim_df <- data.frame("id.iso" = NULL,
                          "id.name" = NULL,
                          "id.string" = NULL)


# Baseline sims
results_06 <- paste0("./data_results/urep_06_", timestamp_05) # use 05 timestamp
dir.create(results_06, recursive = TRUE)


results_folder <- results_06



for (i in 2:6){#1:length(SSA_ISO2)) {
  
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  
  fs_areas_included <- fs_id_link$fs_area[
    which(paste0("fs_id_", fs_id_link$fs_area_id) %in% incompletex$fs_area_id)
  ]
  
  for (j in 1:length(fs_areas_included)) {
    
    area_selected <- fs_areas_included[j]
    
    run_malsim_nets_sequential_v5(
      N_reps= sim_reps,
      areas_per_core = sim_areas_per_core,
      N_cores = sim_cores,
      areas_included = area_selected,
      mass_int_yr = 3,
      only = TRUE,
      use_hipercow = TRUE,
      bv_beta = bv_beta,
      bv_gamma = bv_gamma,
      past_sim = TRUE,
      results_path = results_06,
      under_reporting = TRUE
    )
    
    # # Only 2-year
    # area_sim_id <- paste(
    #   "urep_past_sim",
    #   sim_id,
    #   gsub(" ", "_", area_selected),
    #   sep = "_"
    # )
    # assign(
    #   area_sim_id,
    #   #net_data %>%
    #     run_malsim_nets_sequential_v5(
    #       N_reps= sim_reps,
    #       areas_per_core = sim_areas_per_core,
    #       N_cores = sim_cores,
    #       areas_included = area_selected,
    #       mass_int_yr = 2,
    #       only = TRUE,
    #       use_hipercow = TRUE,
    #       bv_beta = bv_beta,
    #       bv_gamma = bv_gamma,
    #       past_sim = TRUE,
    #       results_path = results_06,
    #       under_reporting = TRUE
    #     )
    # )
    # new_past_sim_df <- data.frame("id.iso" = SSA_ISO2[i],
    #                               "id.name" = area_sim_id,
    #                               "id.string" = get(area_sim_id))
    # urep_past_sim_df %<>% rbind.data.frame(new_past_sim_df)
    
  }
  
}

urep_fut_sim_df <- data.frame("id.iso" = NULL,
                               "id.name" = NULL,
                               "id.string" = NULL)


sim_cores <- 16          # Takes precedence over areas per core

for (i in 1:6){#1:length(SSA_ISO2)) {
  
  for (k in 1:2) {
    
    #fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
    
    fs_areas_included <- fs_id_link$fs_area[
      which(paste0("fs_id_", fs_id_link$fs_area_id) %in% incompletex$fs_area_id)
    ]
    
    for (j in 1:length(fs_areas_included)) {
      
      area_selected <- fs_areas_included[j]
      
      #Only routine + campaigns
      run_malsim_nets_sequential_v5(
        N_reps= sim_reps,
        areas_per_core = sim_areas_per_core,
        N_cores = sim_cores,
        areas_included = area_selected,
        mass_int_yr = mass_int_yr[k],
        only = TRUE,
        use_hipercow = TRUE,
        bv_beta = bv_beta,
        bv_gamma = bv_gamma,
        future_sim = TRUE,
        results_path = results_06,
        under_reporting = TRUE,
        check_write = TRUE
      )

      
      # PBO routine + campaigns
      run_malsim_nets_sequential_v5(
        N_reps= sim_reps,
        areas_per_core = sim_areas_per_core,
        N_cores = sim_cores,
        areas_included = area_selected,
        mass_int_yr = mass_int_yr[k],
        pbo = TRUE,
        use_hipercow = TRUE,
        bv_beta = bv_beta,
        bv_gamma = bv_gamma,
        future_sim = TRUE,
        results_path = results_06,
        under_reporting = TRUE
      )

      # Pyrrole routine + campaigns
        run_malsim_nets_sequential_v5(
          N_reps= sim_reps,
          areas_per_core = sim_areas_per_core,
          N_cores = sim_cores,
          areas_included = area_selected,
          mass_int_yr = mass_int_yr[k],
          pyrrole = TRUE,
          use_hipercow = TRUE,
          bv_beta = bv_beta,
          bv_gamma = bv_gamma,
          future_sim = TRUE,
          results_path = results_06,
          under_reporting = TRUE
        )
      
    }
    
  }
  
}

sim_cores <- 16

for (i in 1:6){#1:length(SSA_ISO2)) {
  
  #fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  
  fs_areas_included <- fs_id_link$fs_area[
    which(paste0("fs_id_", fs_id_link$fs_area_id) %in% incompletex$fs_area_id)
  ]
  
  for (j in 1:length(fs_areas_included)) {
    
    area_selected <- fs_areas_included[j]
    
    #No future nets
    run_malsim_nets_sequential_v5(
      N_reps= sim_reps,
      areas_per_core = sim_areas_per_core,
      N_cores = sim_cores,
      areas_included = area_selected,
      mass_int_yr = 3,
      only = TRUE,
      routine_baseline = TRUE,
      no_future_nets = TRUE,
      use_hipercow = TRUE,
      bv_beta = bv_beta,
      bv_gamma = bv_gamma,
      future_sim = TRUE,
      results_path = results_06,
      under_reporting = TRUE
    )
    
    #Only routine
    run_malsim_nets_sequential_v5(
      N_reps= sim_reps,
      areas_per_core = sim_areas_per_core,
      N_cores = sim_cores,
      areas_included = area_selected,
      mass_int_yr = 3,
      only = TRUE,
      routine_baseline = TRUE,
      use_hipercow = TRUE,
      bv_beta = bv_beta,
      bv_gamma = bv_gamma,
      future_sim = TRUE,
      results_path = results_06,
      under_reporting = TRUE
    )


    # PBO routine
    run_malsim_nets_sequential_v5(
      N_reps= sim_reps,
      areas_per_core = sim_areas_per_core,
      N_cores = sim_cores,
      areas_included = area_selected,
      mass_int_yr = 3,
      pbo = TRUE,
      routine_baseline = TRUE,
      use_hipercow = TRUE,
      bv_beta = bv_beta,
      bv_gamma = bv_gamma,
      future_sim = TRUE,
      results_path = results_06,
      under_reporting = TRUE
    )


    # Pyrrole routine
    run_malsim_nets_sequential_v5(
      N_reps= sim_reps,
      areas_per_core = sim_areas_per_core,
      N_cores = sim_cores,
      areas_included = area_selected,
      mass_int_yr = 3,
      pyrrole = TRUE,
      routine_baseline = TRUE,
      use_hipercow = TRUE,
      bv_beta = bv_beta,
      bv_gamma = bv_gamma,
      future_sim = TRUE,
      results_path = results_06,
      under_reporting = TRUE
    )
    
  }
  
}


read_country_simsx <- function(iso2x, base_pathx) {
  net_strategiesx <- c(
    "no_future_nets",
    "Pyrethroid-only_2_year_interval",
    "Pyrethroid-only_3_year_interval",
    "Pyrethroid-PBO_2_year_interval",
    "Pyrethroid-PBO_3_year_interval",
    "Pyrethroid-Pyrrole_2_year_interval",
    "Pyrethroid-Pyrrole_3_year_interval",
    "Pyrethroid-only_routine_baseline",
    "Pyrethroid-PBO_routine_baseline",
    "Pyrethroid-Pyrrole_routine_baseline"
  )
  
  all_rds_filesx <- purrr::map(net_strategiesx, function(strategyx) {
    strat_pathx <- file.path(base_pathx, iso2x, strategyx)
    if (!dir.exists(strat_pathx)) return(NULL)
    
    fs_dirsx <- list.dirs(
      strat_pathx,
      full.names = TRUE,
      recursive = FALSE
    )
    fs_dirsx <- fs_dirsx[grepl("fs_id_\\d+$", basename(fs_dirsx))]
    
    purrr::map(fs_dirsx, function(fs_dirx) {
      list.files(
        fs_dirx,
        pattern = "^sim\\d+\\.rds$",
        full.names = TRUE
      )
    }) %>% unlist()
  }) %>% unlist()
  
  purrr::map_dfr(all_rds_filesx, function(filex) {
    tryCatch(readRDS(filex), error = function(e) NULL)
  })
}

compile_country_sim_data <- function(
    base_path, start_year, iso2, sim_pop = 100000
) {
  message("Starting simulation compilation for: ", iso2)
  
  # Convert to ISO3 and prepare sitefile path
  iso3 <- countrycode::countrycode(
    iso2, origin = "iso2c", destination = "iso3c"
  )
  sitefile_folder <- "./data_private/newsitefiles/"
  
  message("Reading raw simulation data...")
  raw_data <- read_country_simsx(iso2, base_path)
  
  message("Reading site metadata...")
  site_path <- file.path(sitefile_folder, paste0(iso3, "_site.rds"))
  site_data <- readRDS(site_path)
  pop_tbl <- site_data$population$population_total
  
  message("Merging population and region metadata...")
  sim_data <- raw_data %>%
    mutate(calendar_year = start_year + year_id - 1) %>%
    left_join(
      pop_tbl,
      by = c(
        "fs_name_1" = "name_1",
        "urbanicity" = "urban_rural",
        "calendar_year" = "year"
      )
    ) %>%
    left_join(
      fs_id_link %>%
        select(fs_area_id, new_area_id, EIR_urep_fit),
      by = "fs_area_id"
    )
  
  message("Precomputing matrix lookup indices...")
  sim_data <- sim_data %>%
    mutate(
      idx = long_sample_ids[sample_id],
      idy = map_int(new_area_id, ~ which(usage_list$a == .x)[usage_list$N_t])
    )
  
  message("Calculating untreated net use parameters...")
  
  C0_u_mat <- as.matrix(C0_u)
  D_u_mat <- as.matrix(D_u)
  ret_u_mat <- as.matrix(ret_u)
  invlam_u_mat <- as.matrix(invlam_u)
  
  C0_a_mat <- as.matrix(C0_a)
  D_a_mat <- as.matrix(D_a)
  ret_a_mat <- as.matrix(ret_a)
  invlam_a_mat <- as.matrix(invlam_a)
  
  sim_data <- sim_data %>%
    mutate(
      C0_u = C0_u_mat[cbind(idx, idy)],
      D_u  = D_u_mat[cbind(idx, idy)],
      ret_u = ret_u_mat[cbind(idx, idy)],
      invlam_u = invlam_u_mat[cbind(idx, new_area_id)],
      lam_u = 1 / invlam_u,
      scale_factor = 1 - (
        0.08 + rnormvals[sample_id] * 0.05 / 1.96
      ),
      C0_u = C0_u * scale_factor,
      D_u  = D_u * scale_factor,
      mean_u = case_when(
        !is.na(mass_int_yr) ~ D_u + (
          C0_u / (lam_u * (mass_int_yr * 12))
        ) * (1 - exp(-lam_u * (mass_int_yr * 12))),
        is.na(mass_int_yr) & net_strategy == "no future nets" ~ 0,
        is.na(mass_int_yr) & net_strategy != "no future nets" ~ D_u,
        TRUE ~ NA_real_
      )
    )
  
  message("Calculating treated net access parameters...")
  sim_data <- sim_data %>%
    mutate(
      C0_a = C0_a_mat[cbind(idx, idy)],
      D_a  = D_a_mat[cbind(idx, idy)],
      ret_a = ret_a_mat[cbind(idx, idy)],
      invlam_a = invlam_a_mat[cbind(idx, new_area_id)],
      lam_a = 1 / invlam_a,
      mean_a = case_when(
        !is.na(mass_int_yr) ~ D_a + (
          C0_a / (lam_a * (mass_int_yr * 12))
        ) * (1 - exp(-lam_a * (mass_int_yr * 12))),
        is.na(mass_int_yr) & net_strategy == "no future nets" ~ 0,
        is.na(mass_int_yr) & net_strategy != "no future nets" ~ D_a,
        TRUE ~ NA_real_
      )
    )
  
  message("Calculating mean use given access")
  sim_data <- sim_data %>%
    mutate(
      mean_uga = mean_u / mean_a
    )
  
  message("Calculating total nets distributed and campaign quantifier")
  sim_data <- sim_data %>%
    mutate(
      adj_ann_total_nets_dist =
        adj_ann_routine_nets_dist + adj_ann_camp_nets_dist,
      adj_camp_proc = sim_pop / adj_ann_camp_nets_dist
    )
  
  message("Finalising and assigning object to global environment...")
  assign(paste0(iso2, "_sim_data"), sim_data, envir = .GlobalEnv)
  message(paste0("✔ Finished compiling data for ", iso2))
}



# compile_country_sim_data <- function(base_path, start_year, iso2) {
#   message("Starting simulation compilation for: ", iso2)
#   
#   # Convert to ISO3 and prepare sitefile path
#   iso3 <- countrycode::countrycode(
#     iso2, origin = "iso2c", destination = "iso3c"
#   )
#   sitefile_folder <- "./data_private/newsitefiles/"
#   
#   message("Reading raw simulation data...")
#   raw_data <- read_country_simsx(iso2, base_path)
#   
#   message("Reading site metadata...")
#   site_path <- file.path(sitefile_folder, paste0(iso3, "_site.rds"))
#   site_data <- readRDS(site_path)
#   pop_tbl <- site_data$population$population_total
# 
#   message("Merging population and region metadata...")
#   sim_data <- raw_data %>%
#     mutate(calendar_year = start_year + year_id - 1) %>%
#     left_join(
#       pop_tbl,
#       by = c(
#         "fs_name_1" = "name_1",
#         "urbanicity" = "urban_rural",
#         "calendar_year" = "year"
#       )
#     ) %>%
#     left_join(
#       fs_id_link %>%
#         select(fs_area_id, new_area_id, EIR_urep_fit),
#       by = "fs_area_id"
#     )
#   
#   message("Calculating untreated net use parameters...")
#   sim_data <- sim_data %>%
#     mutate(
#       C0_u = pmap_dbl(list(new_area_id, sample_id), ~ {
#         idx <- long_sample_ids[as.integer(..2)]
#         idy <- which(usage_list$a == as.integer(..1))[usage_list$N_t]
#         C0_u[idx, idy] %>%
#           as.matrix() %>%
#           unname()
#       }),
#       D_u = pmap_dbl(list(new_area_id, sample_id), ~ {
#         idx <- long_sample_ids[as.integer(..2)]
#         idy <- which(usage_list$a == as.integer(..1))[usage_list$N_t]
#         D_u[idx, idy] %>%
#           as.matrix() %>%
#           unname()
#       }),
#       invlam_u = pmap_dbl(list(new_area_id, sample_id), ~ {
#         idx <- long_sample_ids[as.integer(..2)]
#         invlam_u[idx, as.integer(..1)] %>%
#           as.vector() %>%
#           unname() %>%
#           unlist()
#       }),
#       lam_u = 1 / invlam_u,
#       scale_factor = 1 - (
#         0.08 + rnormvals[sample_id] * 0.05 / 1.96
#       ),
#       C0_u = C0_u * scale_factor,
#       D_u = D_u * scale_factor,
#       mean_u = case_when(
#         !is.na(mass_int_yr) ~ D_u + (
#           C0_u / (lam_u * (mass_int_yr * 12))
#         ) * (1 - exp(-lam_u * (mass_int_yr * 12))),
#         is.na(mass_int_yr) & net_strategy == "no future nets" ~ 0,
#         is.na(mass_int_yr) & net_strategy != "no future nets" ~ D_u,
#         TRUE ~ NA_real_
#       )
#     )
#   
#   message("Calculating treated net access parameters...")
#   sim_data <- sim_data %>%
#     mutate(
#       C0_a = pmap_dbl(list(new_area_id, sample_id), ~ {
#         idx <- long_sample_ids[as.integer(..2)]
#         idy <- which(usage_list$a == as.integer(..1))[usage_list$N_t]
#         C0_a[idx, idy] %>%
#           as.matrix() %>%
#           unname()
#       }),
#       D_a = pmap_dbl(list(new_area_id, sample_id), ~ {
#         idx <- long_sample_ids[as.integer(..2)]
#         idy <- which(usage_list$a == as.integer(..1))[usage_list$N_t]
#         D_a[idx, idy] %>%
#           as.matrix() %>%
#           unname()
#       }),
#       invlam_a = pmap_dbl(list(new_area_id, sample_id), ~ {
#         idx <- long_sample_ids[as.integer(..2)]
#         invlam_a[idx, as.integer(..1)] %>%
#           as.vector() %>%
#           unname() %>%
#           unlist()
#       }),
#       lam_a = 1 / invlam_a,
#       mean_a = case_when(
#         !is.na(mass_int_yr) ~ D_a + (
#           C0_a / (lam_a * (mass_int_yr * 12))
#         ) * (1 - exp(-lam_a * (mass_int_yr * 12))),
#         is.na(mass_int_yr) & net_strategy == "no future nets" ~ 0,
#         is.na(mass_int_yr) & net_strategy != "no future nets" ~ D_a,
#         TRUE ~ NA_real_
#       ),
#       adj_ann_total_nets_dist =
#         adj_ann_routine_nets_dist + adj_ann_camp_nets_dist
#     )
#   
#   message("Finalising and assigning object to global environment...")
#   assign(paste0(iso2, "_sim_data"), sim_data, envir = .GlobalEnv)
#   message(paste0("✔ Finished compiling data for ", iso2))
# }

append_use_given_access <- function(iso2) {
  obj_name <- paste0(iso2, "_sim_data")
  
  # Check that object exists
  if (!exists(obj_name, envir = .GlobalEnv)) {
    stop(paste("Object", obj_name, "not found in global environment."))
  }
  
  sim_data <- get(obj_name, envir = .GlobalEnv)
  
  # Check if column already exists
  if ("mean_use_given_access" %in% colnames(sim_data)) {
    response <- readline(
      prompt = paste(
        "Column 'mean_use_given_access' already exists.",
        "\nDo you want to overwrite it? (y/n): "
      )
    )
    if (tolower(response) != "y") {
      message("Aborting without changes.")
      return(invisible(NULL))
    }
  }
  
  # Append the column
  sim_data <- sim_data %>%
    mutate(mean_use_given_access = ifelse(
      is.na(mean_a) | mean_a == 0, NA, mean_u / mean_a
    ))
  
  assign(obj_name, sim_data, envir = .GlobalEnv)
  
  message(paste0("Appended 'mean_use_given_access' to ", obj_name))
}


# append_summary_cases <- function(iso2) {
#   obj_name <- paste0(iso2, "_sim_data")
#   
#   if (!exists(obj_name, envir = .GlobalEnv)) {
#     stop(paste("Object", obj_name, "not found in global environment."))
#   }
#   
#   sim_data <- get(obj_name, envir = .GlobalEnv)
#   
#   clin_all <- c(
#     "clin_cases_0_181_tot", "clin_cases_182_729_tot", 
#     "clin_cases_730_1824_tot", "clin_cases_1825_3649_tot", 
#     "clin_cases_3650_5474_tot", "clin_cases_5475_7299_tot", 
#     "clin_cases_7300_9124_tot", "clin_cases_9125_10949_tot", 
#     "clin_cases_10950_12774_tot", "clin_cases_12775_14599_tot", 
#     "clin_cases_14600_16424_tot", "clin_cases_16425_18249_tot", 
#     "clin_cases_18250_20074_tot", "clin_cases_20075_21899_tot", 
#     "clin_cases_21900_23724_tot", "clin_cases_23725_36499_tot"
#   )
#   
#   clin_u5 <- c(
#     "clin_cases_0_181_tot", "clin_cases_182_729_tot", 
#     "clin_cases_730_1824_tot"
#   )
#   
#   clin_5_14 <- c(
#     "clin_cases_1825_3649_tot", "clin_cases_3650_5474_tot"
#   )
#   
#   clin_15plus <- setdiff(clin_all, c(clin_u5, clin_5_14))
#   
#   sev_all <- sub("clin_", "sev_", clin_all)
#   sev_u5 <- sub("clin_", "sev_", clin_u5)
#   sev_5_14 <- sub("clin_", "sev_", clin_5_14)
#   sev_15plus <- sub("clin_", "sev_", clin_15plus)
#   
#   new_cols <- c(
#     "clin_cases_all_ages", "clin_cases_under5", "clin_cases_5_14", 
#     "clin_cases_15plus", "sev_cases_all_ages", "sev_cases_under5", 
#     "sev_cases_5_14", "sev_cases_15plus"
#   )
#   
#   already_exists <- intersect(new_cols, colnames(sim_data))
#   if (length(already_exists) > 0) {
#     response <- readline(prompt = paste(
#       "The following columns already exist:",
#       paste(already_exists, collapse = ", "),
#       "\nDo you want to overwrite them? (y/n): "
#     ))
#     if (tolower(response) != "y") {
#       message("Aborting without changes.")
#       return(invisible(NULL))
#     }
#   }
#   
#   sim_data <- sim_data %>%
#     mutate(
#       clin_cases_all_ages = rowSums(across(all_of(clin_all)), na.rm = TRUE),
#       clin_cases_under5   = rowSums(across(all_of(clin_u5)), na.rm = TRUE),
#       clin_cases_5_14     = rowSums(across(all_of(clin_5_14)), na.rm = TRUE),
#       clin_cases_15plus   = rowSums(across(all_of(clin_15plus)), na.rm = TRUE),
#       sev_cases_all_ages  = rowSums(across(all_of(sev_all)), na.rm = TRUE),
#       sev_cases_under5    = rowSums(across(all_of(sev_u5)), na.rm = TRUE),
#       sev_cases_5_14      = rowSums(across(all_of(sev_5_14)), na.rm = TRUE),
#       sev_cases_15plus    = rowSums(across(all_of(sev_15plus)), na.rm = TRUE)
#     )
#   
#   assign(obj_name, sim_data, envir = .GlobalEnv)
#   message(paste0("Updated ", obj_name, " with extended summary case columns."))
# }

append_summary_cases <- function(iso2) {
  obj_name <- paste0(iso2, "_sim_data")
  
  if (!exists(obj_name, envir = .GlobalEnv)) {
    stop(paste("Object", obj_name, "not found in global environment."))
  }
  
  sim_data <- get(obj_name, envir = .GlobalEnv)
  
  clin_all <- c(
    "clin_cases_0_181_tot", "clin_cases_182_729_tot",
    "clin_cases_730_1824_tot", "clin_cases_1825_3649_tot",
    "clin_cases_3650_5474_tot", "clin_cases_5475_7299_tot",
    "clin_cases_7300_9124_tot", "clin_cases_9125_10949_tot",
    "clin_cases_10950_12774_tot", "clin_cases_12775_14599_tot",
    "clin_cases_14600_16424_tot", "clin_cases_16425_18249_tot",
    "clin_cases_18250_20074_tot", "clin_cases_20075_21899_tot",
    "clin_cases_21900_23724_tot", "clin_cases_23725_36499_tot"
  )
  
  clin_u5 <- c(
    "clin_cases_0_181_tot", "clin_cases_182_729_tot",
    "clin_cases_730_1824_tot"
  )
  
  clin_5_14 <- c("clin_cases_1825_3649_tot", "clin_cases_3650_5474_tot")
  clin_15plus <- setdiff(clin_all, c(clin_u5, clin_5_14))
  
  sev_all <- sub("clin_cases_", "sev_cases_", clin_all)
  sev_u5 <- sub("clin_cases_", "sev_cases_", clin_u5)
  sev_5_14 <- sub("clin_cases_", "sev_cases_", clin_5_14)
  sev_15plus <- sub("clin_cases_", "sev_cases_", clin_15plus)
  
  n_all <- sub("clin_cases_", "n_age_", clin_all)
  n_u5 <- sub("clin_cases_", "n_age_", clin_u5)
  n_5_14 <- sub("clin_cases_", "n_age_", clin_5_14)
  n_15plus <- sub("clin_cases_", "n_age_", clin_15plus)
  
  det_all <- sub("clin_cases_", "n_detect_lm_", clin_all)
  det_u5 <- sub("clin_cases_", "n_detect_lm_", clin_u5)
  det_5_14 <- sub("clin_cases_", "n_detect_lm_", clin_5_14)
  det_15plus <- sub("clin_cases_", "n_detect_lm_", clin_15plus)
  
  tot_all <- sub("clin_cases_", "cases_", clin_all)
  tot_u5 <- sub("clin_cases_", "cases_", clin_u5)
  tot_5_14 <- sub("clin_cases_", "cases_", clin_5_14)
  tot_15plus <- sub("clin_cases_", "cases_", clin_15plus)
  
  new_cols_list <- list(
    clin_cases_all_ages = clin_all,
    clin_cases_under5 = clin_u5,
    clin_cases_5_14 = clin_5_14,
    clin_cases_15plus = clin_15plus,
    sev_cases_all_ages = sev_all,
    sev_cases_under5 = sev_u5,
    sev_cases_5_14 = sev_5_14,
    sev_cases_15plus = sev_15plus,
    n_age_all_ages = n_all,
    n_age_under5 = n_u5,
    n_age_5_14 = n_5_14,
    n_age_15plus = n_15plus,
    n_detect_lm_all_ages = det_all,
    n_detect_lm_under5 = det_u5,
    n_detect_lm_5_14 = det_5_14,
    n_detect_lm_15plus = det_15plus,
    cases_all_ages = tot_all,
    cases_under5 = tot_u5,
    cases_5_14 = tot_5_14,
    cases_15plus = tot_15plus
  )
  
  for (new_col in names(new_cols_list)) {
    if (new_col %in% colnames(sim_data)) {
      message(paste("Skipping existing column:", new_col))
    } else {
      sim_data[[new_col]] <- sim_data %>%
        select(all_of(new_cols_list[[new_col]])) %>%
        rowSums(na.rm = TRUE)
    }
  }
  
  assign(obj_name, sim_data, envir = .GlobalEnv)
  message(paste0(
    "Updated ", obj_name,
    " with additional summary case and population columns."
  ))
}





generate_pop_summary <- function(iso2) {
  obj_name <- paste0(iso2, "_sim_data")
  
  if (!exists(obj_name, envir = .GlobalEnv)) {
    stop("Object ", obj_name, " does not exist.")
  }
  
  sim_data <- get(obj_name, envir = .GlobalEnv)
  
  # Ensure relevant columns exist
  required_cols <- c(
    "fs_area_id", "fs_area", "fs_name_1", "urbanicity",
    "new_area_id", "calendar_year", "pop"
  )
  missing_cols <- setdiff(required_cols, colnames(sim_data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
  }
  
  message("Extracting population summary for: ", iso2)
  
  pop_sum <- sim_data %>%
    select(all_of(required_cols)) %>%
    distinct() %>%
    arrange(fs_area_id, calendar_year)
  
  assign(paste0(iso2, "_pop_sum"), pop_sum, envir = .GlobalEnv)
  message("✔ Population summary assigned to ", iso2, "_pop_sum")
}


generate_pop_summary <- function(iso2) {
  obj_name <- paste0(iso2, "_sim_data")
  
  if (!exists(obj_name, envir = .GlobalEnv)) {
    stop("Object ", obj_name, " does not exist.")
  }
  
  sim_data <- get(obj_name, envir = .GlobalEnv)
  
  required_cols <- c(
    "fs_area_id", "fs_area", "fs_name_1", "urbanicity",
    "new_area_id", "calendar_year", "sample_id",
    "net_strategy", "pop"
  )
  missing_cols <- setdiff(required_cols, colnames(sim_data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
  }
  
  message("Extracting population summary for: ", iso2)
  
  # Check for consistency of population values
  pop_check <- sim_data %>%
    select(fs_area_id, calendar_year, pop, sample_id, net_strategy) %>%
    group_by(fs_area_id, calendar_year) %>%
    summarise(
      n_distinct_pops = n_distinct(pop),
      .groups = "drop"
    ) %>%
    filter(n_distinct_pops > 1)
  
  if (nrow(pop_check) > 0) {
    warning("⚠ Population is not consistent across sample_id or net_strategy!")
    print(head(pop_check, 10))
  } else {
    message("✔ Population is consistent across sample_id and net_strategy.")
  }
  
  # Construct clean summary
  pop_sum <- sim_data %>%
    select(
      fs_area_id, fs_area, fs_name_1, urbanicity,
      new_area_id, calendar_year, pop
    ) %>%
    distinct() %>%
    arrange(fs_area_id, calendar_year)
  
  assign(paste0(iso2, "_pop_sum"), pop_sum, envir = .GlobalEnv)
  message("✔ Population summary assigned to ", iso2, "_pop_sum")
}

generate_country_pop_summary <- function(iso2) {
  pop_obj <- paste0(iso2, "_pop_sum")
  
  if (!exists(pop_obj, envir = .GlobalEnv)) {
    stop("Object ", pop_obj, " not found. Run generate_pop_summary() first.")
  }
  
  pop_data <- get(pop_obj, envir = .GlobalEnv)
  
  # Create year_id column assuming years are consecutive from min
  year_map <- pop_data %>%
    distinct(calendar_year) %>%
    arrange(calendar_year) %>%
    mutate(year_id = row_number())
  
  ctry_pop <- pop_data %>%
    left_join(year_map, by = "calendar_year") %>%
    group_by(year_id, calendar_year) %>%
    summarise(
      pop = sum(pop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year_id)
  
  # Save to global env
  assign(paste0(iso2, "_ctry_pop_sum"), ctry_pop, envir = .GlobalEnv)
  message("✔ Country-wide population summary saved to ", iso2, "_ctry_pop_sum")
}

append_ctry_pop_and_weights <- function(iso2) {
  sim_name <- paste0(iso2, "_sim_data")
  ctry_pop_name <- paste0(iso2, "_ctry_pop_sum")
  
  # Check objects exist
  if (!exists(sim_name, envir = .GlobalEnv)) {
    stop(paste("Object", sim_name, "not found in global environment."))
  }
  if (!exists(ctry_pop_name, envir = .GlobalEnv)) {
    stop(paste("Object", ctry_pop_name, "not found in global environment."))
  }
  
  sim_data <- get(sim_name, envir = .GlobalEnv)
  ctry_pop_data <- get(ctry_pop_name, envir = .GlobalEnv)
  
  # Columns to be added
  new_cols <- c("ctry_pop", "pop_weight")
  already_exists <- intersect(new_cols, names(sim_data))
  
  # Overwrite prompt
  if (length(already_exists) > 0) {
    response <- readline(
      prompt = paste(
        "The following columns already exist:",
        paste(already_exists, collapse = ", "),
        "\nDo you want to overwrite them? (y/n): "
      )
    )
    if (tolower(response) != "y") {
      message("Aborting without changes.")
      return(invisible(NULL))
    }
  }
  
  message("Joining national population and computing weights...")
  
  updated_data <- sim_data %>%
    left_join(
      ctry_pop_data,
      by = c("year_id", "calendar_year")
    ) %>%
    mutate(
      ctry_pop = pop.y,
      pop_weight = pop.x / pop.y
    ) %>%
    select(-pop.y) %>%
    rename(pop = pop.x)
  
  assign(sim_name, updated_data, envir = .GlobalEnv)
  
  message(paste0("Appended ctry_pop and pop_weight to ", sim_name))
}

# append_real_case_counts <- function(iso2, sim_pop = 100000) {
#   obj_name <- paste0(iso2, "_sim_data")
#   
#   if (!exists(obj_name, envir = .GlobalEnv)) {
#     stop(paste("Object", obj_name, "not found in global environment."))
#   }
#   
#   sim_data <- get(obj_name, envir = .GlobalEnv)
#   
#   base_cols <- c(
#     "clin_cases_all_ages", "clin_cases_under5", "sev_cases_all_ages",
#     "sev_cases_under5", "clin_cases_5_14", "clin_cases_15plus",
#     "sev_cases_5_14", "sev_cases_15plus"
#   )
#   
#   real_cols <- paste0(base_cols, "_real")
#   
#   already_exists <- intersect(real_cols, colnames(sim_data))
#   if (length(already_exists) > 0) {
#     response <- readline(prompt = paste(
#       "The following columns already exist:",
#       paste(already_exists, collapse = ", "),
#       "\nDo you want to overwrite them? (y/n): "
#     ))
#     if (tolower(response) != "y") {
#       message("Aborting without changes.")
#       return(invisible(NULL))
#     }
#   }
#   
#   for (i in seq_along(base_cols)) {
#     sim_data[[real_cols[i]]] <- sim_data[[base_cols[i]]] *
#       (sim_data$pop / sim_pop)
#   }
#   
#   assign(obj_name, sim_data, envir = .GlobalEnv)
#   message(paste0(
#     "Appended real case counts to ", obj_name, 
#     " using sim pop = ", sim_pop, "."
#   ))
# }

append_real_case_counts <- function(iso2, sim_pop = 100000) {
  obj_name <- paste0(iso2, "_sim_data")
  
  if (!exists(obj_name, envir = .GlobalEnv)) {
    stop(paste("Object", obj_name, "not found in global environment."))
  }
  
  sim_data <- get(obj_name, envir = .GlobalEnv)
  
  prefixes <- c(
    "clin_cases_", "sev_cases_",
    "n_age_", "n_detect_lm_", "cases_"
  )
  
  # Identify base columns by prefix
  base_cols <- unlist(lapply(prefixes, function(prefix) {
    grep(paste0("^", prefix), names(sim_data), value = TRUE)
  }))
  
  real_cols <- paste0(base_cols, "_real")
  
  # Check for existing real columns
  existing_cols <- intersect(real_cols, names(sim_data))
  new_cols <- setdiff(real_cols, existing_cols)
  
  if (length(existing_cols) > 0) {
    response <- readline(prompt = paste(
      "The following real columns already exist:",
      paste(existing_cols, collapse = ", "),
      "\nDo you want to overwrite them? (y/n): "
    ))
    if (tolower(response) != "y") {
      message("Only generating real columns not already present.")
      base_cols <- base_cols[!real_cols %in% existing_cols]
      real_cols <- real_cols[!real_cols %in% existing_cols]
    }
  }
  
  if (length(base_cols) == 0) {
    message("No new real columns to compute. Exiting.")
    return(invisible(NULL))
  }
  
  message("Appending real columns using sim population = ", sim_pop, "...")
  
  for (i in seq_along(base_cols)) {
    sim_data[[real_cols[i]]] <- sim_data[[base_cols[i]]] *
      (sim_data$pop / sim_pop)
  }
  
  assign(obj_name, sim_data, envir = .GlobalEnv)
  message("✔ Real case-related columns appended to ", obj_name)
}


append_year_weights <- function(iso2) {
  obj_name <- paste0(iso2, "_sim_data")
  
  if (!exists(obj_name, envir = .GlobalEnv)) {
    stop(paste("Object", obj_name, "not found in global environment."))
  }
  
  sim_data <- get(obj_name, envir = .GlobalEnv)
  
  if ("year_weight" %in% names(sim_data)) {
    response <- readline(
      prompt = "Column 'year_weight' already exists. Overwrite? (y/n): "
    )
    if (tolower(response) != "y") {
      message("Aborting: 'year_weight' not overwritten.")
      return(invisible(NULL))
    }
  }
  
  message("➕ Calculating year_weight within each fs_area_id, sample_id, strategy...")
  
  sim_data <- sim_data %>%
    group_by(fs_area_id, sample_id, net_strategy) %>%
    mutate(
      year_weight = pop / sum(pop, na.rm = TRUE)
    ) %>%
    ungroup()
  
  message("✅ year_weight column appended.")
  
  message("🔍 Checking consistency of year_weight across sample_id and strategy...")
  
  weight_check <- sim_data %>%
    group_by(fs_area_id, year_id) %>%
    summarise(
      n_unique_weights = n_distinct(year_weight),
      .groups = "drop"
    ) %>%
    filter(n_unique_weights > 1)
  
  if (nrow(weight_check) > 0) {
    warning("⚠ year_weight is not consistent across sample_id or net_strategy!")
    print(head(weight_check, 10))
  } else {
    message("✔ year_weight is consistent across sample_id and net_strategy.")
  }
  
  assign(obj_name, sim_data, envir = .GlobalEnv)
  invisible(sim_data)
}


generate_annual_summary <- function(iso2) {
  message("Generating annual summary data for ", iso2, "...")
  
  sim_data_name <- paste0(iso2, "_sim_data")
  ann_data_name <- paste0(iso2, "_ann_data")
  
  if (!exists(sim_data_name, envir = .GlobalEnv)) {
    stop("Data object ", sim_data_name,
         " does not exist in the global environment.")
  }
  
  sim_data <- get(sim_data_name, envir = .GlobalEnv)
  
  message("Grouping and summarising across years...")
  
  ann_data <- sim_data %>%
    group_by(sample_id, fs_area_id, net_strategy) %>%
    summarise(
      across(c(ISO2, fs_name_1, urbanicity, fs_area, net_name,
               mass_int_yr, biennial_reduction, net_costings,
               no_future_nets, routine_baseline, CMC_start,
               CMC_end, ann_routine_nets_dist, ann_camp_nets_dist,
               adj_ann_routine_nets_dist, adj_ann_camp_nets_dist,
               budget_pc, country, iso3c, new_area_id, EIR_urep_fit,
               C0_u, D_u, invlam_u, lam_u, scale_factor, mean_u, ret_u,
               C0_a, D_a, invlam_a, lam_a, mean_a, ret_a, mean_uga,
               adj_camp_proc),
             ~ first(.x), .names = "{.col}"),
      across(c(starts_with("n_age_"), starts_with("n_detect_lm_"),
               starts_with("cases_"), starts_with("clin_cases_"),
               starts_with("sev_cases_"), starts_with("pfpr_")),
             ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
      across(c("pop", "par", "par_pf", "par_pv",
               "adj_ann_total_nets_dist", "mean_use_given_access",
               "ctry_pop",
               "pop_weight"
               ),
             ~ mean(.x, na.rm = TRUE), .names = "{.col}")
    ) %>%
    ungroup()
  
  message("Summary complete. Assigning to ", ann_data_name)
  assign(ann_data_name, ann_data, envir = .GlobalEnv)
}


generate_country_annual_summary <- function(iso2, sim_pop = 100000) {
  ann_name <- paste0(iso2, "_ann_data")
  out_name <- paste0(iso2, "_ctry_ann_sum")
  
  if (!exists(ann_name, envir = .GlobalEnv)) {
    stop(paste("Object", ann_name, "not found in global environment."))
  }
  
  ann_data <- get(ann_name, envir = .GlobalEnv)
  
  summary_df <- ann_data %>%
    dplyr::group_by(sample_id, net_strategy) %>%
    dplyr::summarise(
      ISO2 = dplyr::first(ISO2),
      net_name = dplyr::first(net_name),
      mass_int_yr = dplyr::first(mass_int_yr),
      pop_weight_sum = sum(pop_weight, na.rm = TRUE),
      adj_ann_routine_nets_dist = sum(
        adj_ann_routine_nets_dist * pop_weight, na.rm = TRUE
      ),
      adj_camp_nets_dist = sum(
        adj_ann_camp_nets_dist * pop_weight, na.rm = TRUE
      ),
      clin_cases_all_ages = sum(clin_cases_all_ages * pop_weight, na.rm = TRUE),
      clin_cases_under5 = sum(clin_cases_under5 * pop_weight, na.rm = TRUE),
      clin_cases_5_14 = sum(clin_cases_5_14 * pop_weight, na.rm = TRUE),
      clin_cases_15plus = sum(clin_cases_15plus * pop_weight, na.rm = TRUE),
      mean_u = sum(mean_u * pop_weight, na.rm = TRUE),
      mean_a = sum(mean_a * pop_weight, na.rm = TRUE),
      mean_use_given_access = sum(
        mean_use_given_access * pop_weight, na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  summary_df$adj_camp_proc = sim_pop / summary_df$adj_camp_nets_dist
  
  assign(out_name, summary_df, envir = .GlobalEnv)
  message(paste0("Created ", out_name, " with country-level annual summaries."))
}

append_site_file_data <- function(iso2, ref_year = 2024, sim_pop = 100000,
                            data_obj = "_ann_data") {
  message("Starting simulation compilation for: ", iso2)
  
  obj_name <- paste0(iso2, data_obj)
  
  if (!exists(obj_name, envir = .GlobalEnv)) {
    stop(paste("Object", obj_name, "not found in global environment."))
  }
  
  sim_data <- get(obj_name, envir = .GlobalEnv)
  
  if ("MAP_pfpr" %in% names(sim_data)) {
    response <- readline(
      prompt = "MAP prevalence already exists. Overwrite? (y/n): "
    )
    if (tolower(response) != "y") {
      message("Aborting: MAP prevalence not overwritten.")
      return(invisible(NULL))
    }
  }
  
  # Convert to ISO3 and prepare sitefile path
  iso3 <- countrycode::countrycode(
    iso2, origin = "iso2c", destination = "iso3c"
  )
  sitefile_folder <- "./data_private/newsitefiles/"
  
  message("Reading site metadata...")
  site_path <- file.path(sitefile_folder, paste0(iso3, "_site.rds"))
  site_data <- readRDS(site_path)
  
  prev_tbl <- site_data$prevalence %>%
    dplyr::filter(year == ref_year)
  
  res_tbl <- site_data$vectors$pyrethroid_resistance %>%
    dplyr::filter(year == ref_year)
  
  message("Merging population and region metadata...")
  sim_data <- sim_data %>%
    left_join(
      prev_tbl %>%
        select(name_1, urban_rural, pfpr),
      by = c(
        "fs_name_1" = "name_1",
        "urbanicity" = "urban_rural"
      )
    ) %>%
    dplyr::rename("MAP_pfpr" = "pfpr") %>%
    left_join(
      res_tbl %>%
        select(name_1, pyrethroid_resistance),
      by = c(
        "fs_name_1" = "name_1"
      )
    )
  # sim_data$adj_camp_proc <- sim_pop / sim_data$adj_ann_camp_nets_dist
  
  message("Finalising and assigning object to global environment...")
  assign(paste0(iso2, data_obj), sim_data, envir = .GlobalEnv)
  message(paste0("✔ Finished appending MAP prevalence data for ", iso2))
}





run_full_pipeline <- function(iso2, base_path, start_year = 2025) {
  
  compile_country_sim_data(
    base_path = base_path,
    start_year = start_year,
    iso2 = iso2
  )
  
  append_use_given_access(iso2)
  append_summary_cases(iso2)
  generate_pop_summary(iso2)
  generate_country_pop_summary(iso2)
  append_ctry_pop_and_weights(iso2)
  append_real_case_counts(iso2)
  append_year_weights(iso2)
  generate_annual_summary(iso2)
  append_site_file_data(iso2)
  append_site_file_data(iso2, data_obj = "_sim_data")
  #summary_df$adj_camp_proc = sim_pop / summary_df$adj_camp_nets_dist
  #ML_ann_data$adj_camp_proc <- 100000 / ML_ann_data$adj_ann_camp_nets_dist
  generate_country_annual_summary(iso2)
  
}


shorten_sim_data <- function(iso2, shorten_annual_data = FALSE) {
  
  if (shorten_annual_data) {
    sim_data_name <- paste0(iso2, "_ann_data")
    short_sim_data_name <- paste0(iso2, "_short_ann_data")
  } else {
    sim_data_name <- paste0(iso2, "_sim_data")
    short_sim_data_name <- paste0(iso2, "_short_sim_data")
  }
  
  if (!exists(sim_data_name, envir = .GlobalEnv)) {
    stop("Data object ", sim_data_name,
         " does not exist in the global environment.")
  }
  
  sim_data <- get(sim_data_name, envir = .GlobalEnv)
  
  # short_sim_data <- sim_data %>%
  #   dplyr::select(
  #     ISO2, fs_name_1, urbanicity, fs_area, fs_area_id, new_area_id, pop,
  #     pop_weight, D_u, C0_u, lam_u, D_a, C0_a, lam_a, mean_u, mean_a, ret_u,
  #     ret_a, mean_uga, pyrethroid_resistance,
  #     EIR_urep_fit, year_id, sample_id, net_strategy, net_name, mass_int_yr,
  #     no_future_nets, n_age_under5, adj_ann_total_nets_dist,
  #     adj_ann_camp_nets_dist, adj_ann_routine_nets_dist,
  #     cases_all_ages, cases_under5, clin_cases_all_ages, clin_cases_under5,
  #     sev_cases_all_ages, sev_cases_under5, pfpr_0_36499_mean,
  #     pfpr_182_1824_mean, pfpr_730_3649_mean, MAP_pfpr
  #   )
  
  common_cols <- c(
    "ISO2", "fs_name_1", "urbanicity", "fs_area", "fs_area_id", "new_area_id", "pop",
    "pop_weight", "D_u", "C0_u", "lam_u", "D_a", "C0_a", "lam_a", "mean_u", "mean_a", "ret_u",
    "ret_a", "mean_uga", "pyrethroid_resistance", "EIR_urep_fit", "sample_id",
    "net_strategy", "net_name", "mass_int_yr", "no_future_nets", "n_age_under5",
    "adj_ann_total_nets_dist", "adj_ann_camp_nets_dist", "adj_ann_routine_nets_dist",
    "cases_all_ages", "cases_under5", "clin_cases_all_ages", "clin_cases_under5",
    "sev_cases_all_ages", "sev_cases_under5", "pfpr_0_36499_mean",
    "pfpr_182_1824_mean", "pfpr_730_3649_mean", "MAP_pfpr"
  )
  
  # Conditionally add year_id
  col_list <- if (shorten_annual_data) common_cols else c(common_cols, "year_id")
  
  # Select columns
  short_sim_data <- sim_data %>% dplyr::select(all_of(col_list))
  
  
  message("Shortening of sim data complete. Assigning to ", short_sim_data_name)
  assign(short_sim_data_name, short_sim_data, envir = .GlobalEnv)
  
}

if(!exists("all_sim_data")) {all_sim_data <- data.frame(NULL)}
if(!exists("all_ann_data")) {all_ann_data <- data.frame(NULL)}

run_full_pipeline(
  iso2 = "BF",
  base_path = "data_results/urep_06_20250604112237/summary"
)
shorten_sim_data("BF")
shorten_sim_data("BF", shorten_annual_data = TRUE)
all_sim_data <- rbind.data.frame(all_sim_data, BF_short_sim_data)
all_ann_data <- rbind.data.frame(all_ann_data, BF_short_ann_data)
write.csv(BF_sim_data, paste0("BF_sim_data", timestamp_05, ".csv"))
rm(BF_sim_data)
write.csv(BF_ann_data, paste0("BF_ann_data", timestamp_05, ".csv"))
rm(BF_ann_data)
gc()

run_full_pipeline(
  iso2 = "GH",
  base_path = "data_results/urep_06_20250604112237/summary"
)
shorten_sim_data("GH")
shorten_sim_data("GH", shorten_annual_data = TRUE)
all_sim_data <- rbind.data.frame(all_sim_data, GH_short_sim_data)
all_ann_data <- rbind.data.frame(all_ann_data, GH_short_ann_data)
write.csv(GH_sim_data, paste0("GH_sim_data", timestamp_05, ".csv"))
rm(GH_sim_data)
write.csv(GH_ann_data, paste0("GH_ann_data", timestamp_05, ".csv"))
rm(GH_ann_data)
gc()

run_full_pipeline(
  iso2 = "MW",
  base_path = "data_results/urep_06_20250604112237/summary"
)
shorten_sim_data("MW")
shorten_sim_data("MW", shorten_annual_data = TRUE)
all_sim_data <- rbind.data.frame(all_sim_data, MW_short_sim_data)
all_ann_data <- rbind.data.frame(all_ann_data, MW_short_ann_data)
# rm(MW_sim_data)
# rm(MW_ann_data)
gc()

run_full_pipeline(
  iso2 = "ML",
  base_path = "data_results/urep_06_20250604112237/summary"
)
shorten_sim_data("ML")
shorten_sim_data("ML", shorten_annual_data = TRUE)
all_sim_data <- rbind.data.frame(all_sim_data, ML_short_sim_data)
all_ann_data <- rbind.data.frame(all_ann_data, ML_short_ann_data)
write.csv(ML_sim_data, paste0("ML_sim_data", timestamp_05, ".csv"))
rm(ML_sim_data)
write.csv(ML_ann_data, paste0("ML_ann_data", timestamp_05, ".csv"))
rm(ML_ann_data)
gc()

run_full_pipeline(
  iso2 = "MZ",
  base_path = "data_results/urep_06_20250604112237/summary"
)
shorten_sim_data("MZ")
shorten_sim_data("MZ", shorten_annual_data = TRUE)
all_sim_data <- rbind.data.frame(all_sim_data, MZ_short_sim_data)
all_ann_data <- rbind.data.frame(all_ann_data, MZ_short_ann_data)
rm(MZ_sim_data)
rm(MZ_ann_data)
gc()

run_full_pipeline(
  iso2 = "SN",
  base_path = "data_results/urep_06_20250604112237/summary"
)
shorten_sim_data("SN")
shorten_sim_data("SN", shorten_annual_data = TRUE)
all_sim_data <- rbind.data.frame(all_sim_data, SN_short_sim_data)
all_ann_data <- rbind.data.frame(all_ann_data, SN_short_ann_data)
rm(SN_sim_data)
rm(SN_ann_data)
gc()


calculate_cases_per_100k_under5 <- function(
    iso2 = NULL, obj_name_ending = "_short_ann_data", obj_name = NULL
) {
  
  if (!is.null(iso2)) {
    obj_name <- paste0(iso2, obj_name_ending)
  }
  
  # Check that object exists
  if (!exists(obj_name, envir = .GlobalEnv)) {
    stop(paste("Object", obj_name, "not found in global environment."))
  }
  
  obj_data <- get(obj_name, envir = .GlobalEnv)
  
  # Check if column already exists
  if (
    ("cases_under5_p100kU5" %in% colnames(obj_data)) |
    ("clin_cases_under5_p100kU5" %in% colnames(obj_data)) |
    ("sev_cases_under5_p100kU5" %in% colnames(obj_data))
  ) {
    response <- readline(
      prompt = paste(
        "At least one case column per 100k under 5s already exists.",
        "\nDo you want to overwrite it? (y/n): "
      )
    )
    if (tolower(response) != "y") {
      message("Aborting without changes.")
      return(invisible(NULL))
    }
  }
  
  obj_data <- obj_data %>%
    dplyr::mutate(
      cases_under5_p100kU5 = cases_under5 * 100000 / n_age_under5,
      clin_cases_under5_p100kU5 = clin_cases_under5 * 100000 / n_age_under5,
      sev_cases_under5_p100kU5 = sev_cases_under5 * 100000 / n_age_under5
    )
  
  message("Appending cases per 100k under 5s complete. Assigning to ",
          obj_name)
  assign(obj_name, obj_data, envir = .GlobalEnv)
  
}

calculate_cases_per_100k_under5(obj_name = "all_sim_data")
calculate_cases_per_100k_under5(obj_name = "all_ann_data")


create_full_sim_comparison <- function(
    iso2 = NULL, short_sim_data_name = NULL
    ) {
  
  if (is.null(iso2)) {
    data_comparison_name <- paste0(short_sim_data_name, "_comp")
  } else {
    short_sim_data_name <- paste0(iso2, "_short_ann_data")
    data_comparison_name <- paste0(iso2, "_comp_data")
  }
  
  if (!exists(short_sim_data_name, envir = .GlobalEnv)) {
    stop("Data object ", short_sim_data_name,
         " does not exist in the global environment.")
  }
  
  short_sim_data <- get(short_sim_data_name, envir = .GlobalEnv)
  
  unique_strategies <- unique(short_sim_data$net_strategy)
  N_strategies <- length(unique_strategies)
  
  # Initialise dataframe for all combinations
  all_comb_sim_data <- data.frame(NULL)
  
  #print(short_sim_data_name)
  #print(names(short_sim_data))
  
  # Loop over comparator strategies
  for (i in 1:N_strategies) {
    
    current_comp_strategy <- unique_strategies[i]
    message("Processing comparator strategy: ", current_comp_strategy,
            " (", i, "/", N_strategies, ")")
    
    # Filter for comparator
    comp_sim_data <- short_sim_data %>%
      dplyr::filter(net_strategy == unique_strategies[i]) %>%
      dplyr::rename(
        net_strategy_comp = net_strategy,
        net_name_comp = net_name,
        mass_int_yr_comp = mass_int_yr,
        no_future_nets_comp = no_future_nets, 
        n_age_under5_comp = n_age_under5,
        mean_u_comp = mean_u,
        mean_a_comp = mean_a,
        mean_uga_comp = mean_uga,
        ret_u_comp = ret_u,
        ret_a_comp = ret_a,
        adj_ann_total_nets_dist_comp = adj_ann_total_nets_dist,
        cases_all_ages_comp = cases_all_ages,
        cases_under5_comp = cases_under5,
        cases_under5_p100kU5_comp = cases_under5_p100kU5,
        clin_cases_all_ages_comp = clin_cases_all_ages,
        clin_cases_under5_comp = clin_cases_under5,
        clin_cases_under5_p100kU5_comp = clin_cases_under5_p100kU5,
        sev_cases_all_ages_comp = sev_cases_all_ages,
        sev_cases_under5_comp = sev_cases_under5,
        sev_cases_under5_p100kU5_comp = sev_cases_under5_p100kU5,
        pfpr_0_36499_mean_comp = pfpr_0_36499_mean,
        pfpr_182_1824_mean_comp = pfpr_182_1824_mean,
        pfpr_730_3649_mean_comp = pfpr_730_3649_mean
      )
    
    # Loop over intervention strategies
    for (j in 1:N_strategies) {
      
      # Filter for intervention
      int_sim_data <- short_sim_data %>%
        dplyr::filter(net_strategy == unique_strategies[j]) %>%
        dplyr::select(
          net_strategy, net_name, mass_int_yr, no_future_nets, n_age_under5,
          mean_u, mean_a, mean_uga, ret_u, ret_a,# pyrethroid_resistance,
          adj_ann_total_nets_dist, cases_all_ages, cases_under5,
          cases_under5_p100kU5, clin_cases_all_ages, clin_cases_under5,
          clin_cases_under5_p100kU5, sev_cases_all_ages, sev_cases_under5,
          sev_cases_under5_p100kU5, pfpr_0_36499_mean, pfpr_182_1824_mean,
          pfpr_730_3649_mean
        ) %>%
        dplyr::rename(
          net_strategy_int = net_strategy,
          net_name_int = net_name,
          mass_int_yr_int = mass_int_yr,
          no_future_nets_int = no_future_nets, 
          n_age_under5_int = n_age_under5,
          mean_u_int = mean_u,
          mean_a_int = mean_a,
          mean_uga_int = mean_uga,
          ret_u_int = ret_u,
          ret_a_int = ret_a,
          adj_ann_total_nets_dist_int = adj_ann_total_nets_dist,
          cases_all_ages_int = cases_all_ages,
          cases_under5_int = cases_under5,
          cases_under5_p100kU5_int = cases_under5_p100kU5,
          clin_cases_all_ages_int = clin_cases_all_ages,
          clin_cases_under5_int = clin_cases_under5,
          clin_cases_under5_p100kU5_int = clin_cases_under5_p100kU5,
          sev_cases_all_ages_int = sev_cases_all_ages,
          sev_cases_under5_int = sev_cases_under5,
          sev_cases_under5_p100kU5_int = sev_cases_under5_p100kU5,
          pfpr_0_36499_mean_int = pfpr_0_36499_mean,
          pfpr_182_1824_mean_int = pfpr_182_1824_mean,
          pfpr_730_3649_mean_int = pfpr_730_3649_mean
        )
      
      # Combine selected comparator and intervention dataframes
      comb_sim_data <- cbind.data.frame(comp_sim_data, int_sim_data)
      
      # Append to dataframe with for all comparisons
      all_comb_sim_data <- rbind.data.frame(all_comb_sim_data, comb_sim_data)
      
    }
  }
  
  print(names(all_comb_sim_data))
  
  # find differences
  all_comb_sim_data <- all_comb_sim_data %>%
    dplyr::mutate(
      
      adj_ann_total_nets_dist_add =
        adj_ann_total_nets_dist_int - adj_ann_total_nets_dist_comp,
      
      cases_all_ages_avert =
        cases_all_ages_comp - cases_all_ages_int,
      
      cases_under5_avert =
        cases_under5_comp - cases_under5_int,
      
      cases_under5_p100kU5_avert =
        cases_under5_p100kU5_comp - cases_under5_p100kU5_int,
      
      clin_cases_all_ages_avert =
        clin_cases_all_ages_comp - clin_cases_all_ages_int,
      
      clin_cases_under5_avert =
        clin_cases_under5_comp - clin_cases_under5_int,
      
      clin_cases_under5_p100kU5_avert =
        clin_cases_under5_p100kU5_comp - clin_cases_under5_p100kU5_int,
      
      sev_cases_all_ages_avert =
        sev_cases_all_ages_comp - sev_cases_all_ages_int,
      
      sev_cases_under5_avert =
        sev_cases_under5_comp - sev_cases_under5_int,
      
      sev_cases_under5_p100kU5_avert =
        sev_cases_under5_p100kU5_comp - sev_cases_under5_p100kU5_int,
      
      pfpr_0_36499_mean_red = 
        pfpr_0_36499_mean_comp - pfpr_0_36499_mean_int,
      
      pfpr_182_1824_mean_red = 
        pfpr_182_1824_mean_comp - pfpr_182_1824_mean_int,
      
      pfpr_730_3649_mean_red = 
        pfpr_730_3649_mean_comp - pfpr_730_3649_mean_int,
      
    ) %>%
    dplyr::mutate(
      
      adj_ann_total_nets_dist_propinc =
        adj_ann_total_nets_dist_add / adj_ann_total_nets_dist_comp,
      
      cases_all_ages_propdec =
        cases_all_ages_avert / cases_all_ages_comp,
      
      cases_under5_propdec =
        cases_under5_avert / cases_under5_comp,
      
      cases_under5_p100kU5_propdec =
        cases_under5_p100kU5_avert / cases_under5_p100kU5_comp,
      
      clin_cases_all_ages_propdec =
        clin_cases_all_ages_avert / clin_cases_all_ages_comp,
      
      clin_cases_under5_propdec =
        clin_cases_under5_avert / clin_cases_under5_comp,
      
      clin_cases_under5_p100kU5_propdec =
        clin_cases_under5_p100kU5_avert / clin_cases_under5_p100kU5_comp,
      
      sev_cases_all_ages_propdec =
        sev_cases_all_ages_avert / sev_cases_all_ages_comp,
      
      sev_cases_under5_propdec =
        sev_cases_under5_avert / sev_cases_under5_comp,
      
      sev_cases_under5_p100kU5_propdec =
        sev_cases_under5_p100kU5_avert / sev_cases_under5_p100kU5_comp,
    )
  
  message("Comparisons calculated. Assigning to ", data_comparison_name)
  assign(data_comparison_name, all_comb_sim_data, envir = .GlobalEnv)
  
}


#create_full_sim_comparison(short_sim_data_name = "all_sim_data")
create_full_sim_comparison(short_sim_data_name = "all_ann_data")



summarise_full_comparison <- function(iso2 = NULL, all_countries = TRUE,
                                      summarise_annual_data = TRUE,
                                      append_3yr_comparison = TRUE,
                                      append_3yrpbo_comparison = TRUE,
                                      append_3yrpyrrole_comparison = TRUE) {
  
  if (all_countries) {
    if (summarise_annual_data) {
      data_comparison_name <- "all_ann_data_comp"
      sum_comparison_name <- "all_ann_data_sum"
    } else {
      data_comparison_name <- "all_sim_data_comp"
      sum_comparison_name <- "all_sim_data_sum"
    }
  } else {
    data_comparison_name <- paste0(iso2, "_comp_data")
    sum_comparison_name <- paste0(iso2, "_comp_sum")
  }
  
  if (!exists(data_comparison_name, envir = .GlobalEnv)) {
    stop("Data object ", data_comparison_name,
         " does not exist in the global environment.")
  }
  
  all_comb_sim_data <- get(data_comparison_name, envir = .GlobalEnv)
  
  names(all_comb_sim_data)
  
  all_comb_sim_sum <- all_comb_sim_data %>%
    dplyr::group_by(
      ISO2, fs_name_1, urbanicity, fs_area, fs_area_id, new_area_id, pop,
      pop_weight, EIR_urep_fit, pyrethroid_resistance,
      net_strategy_comp, net_name_comp, mass_int_yr_comp,
      no_future_nets_comp, net_strategy_int, net_name_int, mass_int_yr_int,
      no_future_nets_int
    ) %>%
    dplyr::summarise(
      
      
      clin_cases_all_ages_comp_med = median(
        clin_cases_all_ages_comp, na.rm = TRUE
      ),
      clin_cases_all_ages_comp_lo = quantile(
        clin_cases_all_ages_comp, 0.025, na.rm = TRUE
      ),
      clin_cases_all_ages_comp_hi = quantile(
        clin_cases_all_ages_comp, 0.975, na.rm = TRUE
      ),
      
      clin_cases_under5_p100kU5_comp_med = median(
        clin_cases_under5_p100kU5_comp, na.rm = TRUE
      ),
      clin_cases_under5_p100kU5_comp_lo = quantile(
        clin_cases_under5_p100kU5_comp, 0.025, na.rm = TRUE
      ),
      clin_cases_under5_p100kU5_comp_hi = quantile(
        clin_cases_under5_p100kU5_comp, 0.975, na.rm = TRUE
      ),
      
      pfpr_182_1824_mean_comp_med = median(
        pfpr_182_1824_mean_comp, na.rm = TRUE
      ),
      pfpr_182_1824_mean_comp_lo = quantile(
        pfpr_182_1824_mean_comp, 0.025, na.rm = TRUE
      ),
      pfpr_182_1824_mean_comp_hi = quantile(
        pfpr_182_1824_mean_comp, 0.975, na.rm = TRUE
      ),
      
      pfpr_730_3649_mean_comp_med = median(
        pfpr_730_3649_mean_comp, na.rm = TRUE
      ),
      pfpr_730_3649_mean_comp_lo = quantile(
        pfpr_730_3649_mean_comp, 0.025, na.rm = TRUE
      ),
      pfpr_730_3649_mean_comp_hi = quantile(
        pfpr_730_3649_mean_comp, 0.975, na.rm = TRUE
      ),
      
      pfpr_0_36499_mean_comp_med = median(
        pfpr_0_36499_mean_comp, na.rm = TRUE
      ),
      pfpr_0_36499_mean_comp_lo = quantile(
        pfpr_0_36499_mean_comp, 0.025, na.rm = TRUE
      ),
      pfpr_0_36499_mean_comp_hi = quantile(
        pfpr_0_36499_mean_comp, 0.975, na.rm = TRUE
      ),
      
      mean_u_int_med = median(
        mean_u_int, na.rm = TRUE
      ),
      mean_u_int_lo = quantile(
        mean_u_int, 0.025, na.rm = TRUE
      ),
      mean_u_int_hi = quantile(
        mean_u_int, 0.975, na.rm = TRUE
      ),
      
      mean_a_int_med = median(
        mean_a_int, na.rm = TRUE
      ),
      mean_a_int_lo = quantile(
        mean_a_int, 0.025, na.rm = TRUE
      ),
      mean_a_int_hi = quantile(
        mean_a_int, 0.975, na.rm = TRUE
      ),
      
      mean_uga_int_med = median(
        mean_uga_int, na.rm = TRUE
      ),
      mean_uga_int_lo = quantile(
        mean_uga_int, 0.025, na.rm = TRUE
      ),
      mean_uga_int_hi = quantile(
        mean_uga_int, 0.975, na.rm = TRUE
      ),
      
      ret_u_int_med = median(
        ret_u_int, na.rm = TRUE
      ),
      ret_u_int_lo = quantile(
        ret_u_int, 0.025, na.rm = TRUE
      ),
      ret_u_int_hi = quantile(
        ret_u_int, 0.975, na.rm = TRUE
      ),
      
      ret_a_int_med = median(
        ret_a_int, na.rm = TRUE
      ),
      ret_a_int_lo = quantile(
        ret_a_int, 0.025, na.rm = TRUE
      ),
      ret_a_int_hi = quantile(
        ret_a_int, 0.975, na.rm = TRUE
      ),
      
      adj_ann_total_nets_dist_int_med = median(
        adj_ann_total_nets_dist_int, na.rm = TRUE
      ),
      adj_ann_total_nets_dist_int_lo = quantile(
        adj_ann_total_nets_dist_int, 0.025, na.rm = TRUE
      ),
      adj_ann_total_nets_dist_int_hi = quantile(
        adj_ann_total_nets_dist_int, 0.975, na.rm = TRUE
      ),
      
      adj_ann_total_nets_dist_add_med = median(
        adj_ann_total_nets_dist_add, na.rm = TRUE
      ),
      adj_ann_total_nets_dist_add_lo = quantile(
        adj_ann_total_nets_dist_add, 0.025, na.rm = TRUE
      ),
      adj_ann_total_nets_dist_add_hi = quantile(
        adj_ann_total_nets_dist_add, 0.975, na.rm = TRUE
      ),
      
      adj_ann_total_nets_dist_propinc_med = median(
        adj_ann_total_nets_dist_propinc, na.rm = TRUE
      ),
      adj_ann_total_nets_dist_propinc_lo = quantile(
        adj_ann_total_nets_dist_propinc, 0.025, na.rm = TRUE
      ),
      adj_ann_total_nets_dist_propinc_hi = quantile(
        adj_ann_total_nets_dist_propinc, 0.975, na.rm = TRUE
      ),
      
      cases_all_ages_avert_med = median(
        cases_all_ages_avert, na.rm = TRUE
      ),
      cases_all_ages_avert_lo = quantile(
        cases_all_ages_avert, 0.025, na.rm = TRUE
      ),
      cases_all_ages_avert_hi = quantile(
        cases_all_ages_avert, 0.975, na.rm = TRUE
      ),
      
      cases_all_ages_propdec_med = median(
        cases_all_ages_propdec, na.rm = TRUE
      ),
      cases_all_ages_propdec_lo = quantile(
        cases_all_ages_propdec, 0.025, na.rm = TRUE
      ),
      cases_all_ages_propdec_hi = quantile(
        cases_all_ages_propdec, 0.975, na.rm = TRUE
      ),
      
      cases_under5_p100kU5_avert_med = median(
        cases_under5_p100kU5_avert, na.rm = TRUE
      ),
      cases_under5_p100kU5_avert_lo = quantile(
        cases_under5_p100kU5_avert, 0.025, na.rm = TRUE
      ),
      cases_under5_p100kU5_avert_hi = quantile(
        cases_under5_p100kU5_avert, 0.975, na.rm = TRUE
      ),
      
      cases_under5_p100kU5_propdec_med = median(
        cases_under5_p100kU5_propdec, na.rm = TRUE
      ),
      cases_under5_p100kU5_propdec_lo = quantile(
        cases_under5_p100kU5_propdec, 0.025, na.rm = TRUE
      ),
      cases_under5_p100kU5_propdec_hi = quantile(
        cases_under5_p100kU5_propdec, 0.975, na.rm = TRUE
      ),
      
      clin_cases_all_ages_avert_med = median(
        clin_cases_all_ages_avert, na.rm = TRUE
      ),
      clin_cases_all_ages_avert_lo = quantile(
        clin_cases_all_ages_avert, 0.025, na.rm = TRUE
      ),
      clin_cases_all_ages_avert_hi = quantile(
        clin_cases_all_ages_avert, 0.975, na.rm = TRUE
      ),
      
      clin_cases_all_ages_propdec_med = median(
        clin_cases_all_ages_propdec, na.rm = TRUE
      ),
      clin_cases_all_ages_propdec_lo = quantile(
        clin_cases_all_ages_propdec, 0.025, na.rm = TRUE
      ),
      clin_cases_all_ages_propdec_hi = quantile(
        clin_cases_all_ages_propdec, 0.975, na.rm = TRUE
      ),
      
      clin_cases_under5_p100kU5_avert_med = median(
        clin_cases_under5_p100kU5_avert, na.rm = TRUE
      ),
      clin_cases_under5_p100kU5_avert_lo = quantile(
        clin_cases_under5_p100kU5_avert, 0.025, na.rm = TRUE
      ),
      clin_cases_under5_p100kU5_avert_hi = quantile(
        clin_cases_under5_p100kU5_avert, 0.975, na.rm = TRUE
      ),
      
      clin_cases_under5_p100kU5_propdec_med = median(
        clin_cases_under5_p100kU5_propdec, na.rm = TRUE
      ),
      clin_cases_under5_p100kU5_propdec_lo = quantile(
        clin_cases_under5_p100kU5_propdec, 0.025, na.rm = TRUE
      ),
      clin_cases_under5_p100kU5_propdec_hi = quantile(
        clin_cases_under5_p100kU5_propdec, 0.975, na.rm = TRUE
      ),
      
      sev_cases_all_ages_avert_med = median(
        sev_cases_all_ages_avert, na.rm = TRUE
      ),
      sev_cases_all_ages_avert_lo = quantile(
        sev_cases_all_ages_avert, 0.025, na.rm = TRUE
      ),
      sev_cases_all_ages_avert_hi = quantile(
        sev_cases_all_ages_avert, 0.975, na.rm = TRUE
      ),
      
      sev_cases_all_ages_propdec_med = median(
        sev_cases_all_ages_propdec, na.rm = TRUE
      ),
      sev_cases_all_ages_propdec_lo = quantile(
        sev_cases_all_ages_propdec, 0.025, na.rm = TRUE
      ),
      sev_cases_all_ages_propdec_hi = quantile(
        sev_cases_all_ages_propdec, 0.975, na.rm = TRUE
      ),

      sev_cases_under5_p100kU5_avert_med = median(
        sev_cases_under5_p100kU5_avert, na.rm = TRUE
      ),
      sev_cases_under5_p100kU5_avert_lo = quantile(
        sev_cases_under5_p100kU5_avert, 0.025, na.rm = TRUE
      ),
      sev_cases_under5_p100kU5_avert_hi = quantile(
        sev_cases_under5_p100kU5_avert, 0.975, na.rm = TRUE
      ),
      
      sev_cases_under5_p100kU5_propdec_med = median(
        sev_cases_under5_p100kU5_propdec, na.rm = TRUE
      ),
      sev_cases_under5_p100kU5_propdec_lo = quantile(
        sev_cases_under5_p100kU5_propdec, 0.025, na.rm = TRUE
      ),
      sev_cases_under5_p100kU5_propdec_hi = quantile(
        sev_cases_under5_p100kU5_propdec, 0.975, na.rm = TRUE
      ),
      
      pfpr_0_36499_mean_int_med = median(
        pfpr_0_36499_mean_int, na.rm = TRUE
      ),
      pfpr_0_36499_mean_int_lo = quantile(
        pfpr_0_36499_mean_int, 0.025, na.rm = TRUE
      ),
      pfpr_0_36499_mean_int_hi = quantile(
        pfpr_0_36499_mean_int, 0.975, na.rm = TRUE
      ),
      
      pfpr_182_1824_mean_int_med = median(
        pfpr_182_1824_mean_int, na.rm = TRUE
      ),
      pfpr_182_1824_mean_int_lo = quantile(
        pfpr_182_1824_mean_int, 0.025, na.rm = TRUE
      ),
      pfpr_182_1824_mean_int_hi = quantile(
        pfpr_182_1824_mean_int, 0.975, na.rm = TRUE
      ),
      
      pfpr_730_3649_mean_int_med = median(
        pfpr_730_3649_mean_int, na.rm = TRUE
      ),
      pfpr_730_3649_mean_int_lo = quantile(
        pfpr_730_3649_mean_int, 0.025, na.rm = TRUE
      ),
      pfpr_730_3649_mean_int_hi = quantile(
        pfpr_730_3649_mean_int, 0.975, na.rm = TRUE
      ),
      
      pfpr_0_36499_mean_red_med = median(
        pfpr_0_36499_mean_red, na.rm = TRUE
      ),
      pfpr_0_36499_mean_red_lo = quantile(
        pfpr_0_36499_mean_red, 0.025, na.rm = TRUE
      ),
      pfpr_0_36499_mean_red_hi = quantile(
        pfpr_0_36499_mean_red, 0.975, na.rm = TRUE
      ),
      
      pfpr_182_1824_mean_red_med = median(
        pfpr_182_1824_mean_red, na.rm = TRUE
      ),
      pfpr_182_1824_mean_red_lo = quantile(
        pfpr_182_1824_mean_red, 0.025, na.rm = TRUE
      ),
      pfpr_182_1824_mean_red_hi = quantile(
        pfpr_182_1824_mean_red, 0.975, na.rm = TRUE
      ),
      
      pfpr_730_3649_mean_red_med = median(
        pfpr_730_3649_mean_red, na.rm = TRUE
      ),
      pfpr_730_3649_mean_red_lo = quantile(
        pfpr_730_3649_mean_red, 0.025, na.rm = TRUE
      ),
      pfpr_730_3649_mean_red_hi = quantile(
        pfpr_730_3649_mean_red, 0.975, na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  message("Medians and credible intervals calculated.")
  
  if (append_3yr_comparison) {
    
    comparison_3yr <- all_comb_sim_sum %>%
      dplyr::filter(
        net_strategy_int == "Pyrethroid-only 3-year campaigns uncosted",
        net_strategy_comp == "Pyrethroid-only 3-year campaigns uncosted"
      ) %>%
      dplyr::select(
        fs_area,
        mean_u_int_med, mean_u_int_lo, mean_u_int_hi,
        mean_a_int_med, mean_a_int_lo, mean_a_int_hi,
        mean_uga_int_med, mean_uga_int_lo, mean_uga_int_hi,
        ret_u_int_med, ret_u_int_lo, ret_u_int_hi,
        ret_a_int_med, ret_a_int_lo, ret_a_int_hi
      ) %>%
      dplyr::distinct() %>%
      dplyr::rename(
        mean_u_3yr_med = mean_u_int_med,
        mean_u_3yr_lo  = mean_u_int_lo,
        mean_u_3yr_hi  = mean_u_int_hi,
        mean_a_3yr_med = mean_a_int_med,
        mean_a_3yr_lo  = mean_a_int_lo,
        mean_a_3yr_hi  = mean_a_int_hi,
        mean_uga_3yr_med = mean_uga_int_med,
        mean_uga_3yr_lo  = mean_uga_int_lo,
        mean_uga_3yr_hi  = mean_uga_int_hi,
        ret_u_3yr_med = ret_u_int_med,
        ret_u_3yr_lo  = ret_u_int_lo,
        ret_u_3yr_hi  = ret_u_int_hi,
        ret_a_3yr_med = ret_a_int_med,
        ret_a_3yr_lo  = ret_a_int_lo,
        ret_a_3yr_hi  = ret_a_int_hi
      )
    
    all_comb_sim_sum <- all_comb_sim_sum %>%
      dplyr::left_join(comparison_3yr, by = "fs_area")
    
  }
  
  if (append_3yrpbo_comparison) {
    
    comparison_3yrpbo <- all_comb_sim_sum %>%
      dplyr::filter(
        net_strategy_comp == "Pyrethroid-PBO 3-year campaigns uncosted",
        net_strategy_int == "Pyrethroid-PBO 3-year campaigns uncosted"
      ) %>%
      dplyr::select(
        fs_area,
        clin_cases_all_ages_comp_med,
        clin_cases_all_ages_comp_lo,
        clin_cases_all_ages_comp_hi,
        clin_cases_under5_p100kU5_comp_med,
        clin_cases_under5_p100kU5_comp_lo,
        clin_cases_under5_p100kU5_comp_hi,
        pfpr_182_1824_mean_comp_med,
        pfpr_182_1824_mean_comp_lo,
        pfpr_182_1824_mean_comp_hi,
        pfpr_730_3649_mean_comp_med,
        pfpr_730_3649_mean_comp_lo,
        pfpr_730_3649_mean_comp_hi,
        pfpr_0_36499_mean_comp_med,
        pfpr_0_36499_mean_comp_lo,
        pfpr_0_36499_mean_comp_hi
      ) %>%
      dplyr::distinct() %>%
      dplyr::rename(
        clin_cases_all_ages_3yrpbo_med = clin_cases_all_ages_comp_med,
        clin_cases_all_ages_3yrpbo_lo  = clin_cases_all_ages_comp_lo,
        clin_cases_all_ages_3yrpbo_hi  = clin_cases_all_ages_comp_hi,
        clin_cases_under5_p100kU5_3yrpbo_med = clin_cases_under5_p100kU5_comp_med,
        clin_cases_under5_p100kU5_3yrpbo_lo  = clin_cases_under5_p100kU5_comp_lo,
        clin_cases_under5_p100kU5_3yrpbo_hi  = clin_cases_under5_p100kU5_comp_hi,
        pfpr_182_1824_mean_3yrpbo_med = pfpr_182_1824_mean_comp_med,
        pfpr_182_1824_mean_3yrpbo_lo  = pfpr_182_1824_mean_comp_lo,
        pfpr_182_1824_mean_3yrpbo_hi  = pfpr_182_1824_mean_comp_hi,
        pfpr_730_3649_mean_3yrpbo_med = pfpr_730_3649_mean_comp_med,
        pfpr_730_3649_mean_3yrpbo_lo  = pfpr_730_3649_mean_comp_lo,
        pfpr_730_3649_mean_3yrpbo_hi  = pfpr_730_3649_mean_comp_hi,
        pfpr_0_36499_mean_3yrpbo_med = pfpr_0_36499_mean_comp_med,
        pfpr_0_36499_mean_3yrpbo_lo  = pfpr_0_36499_mean_comp_lo,
        pfpr_0_36499_mean_3yrpbo_hi  = pfpr_0_36499_mean_comp_hi
      )
    
    all_comb_sim_sum <- all_comb_sim_sum %>%
      dplyr::left_join(comparison_3yrpbo, by = "fs_area")
    
  }
  
  if (append_3yrpyrrole_comparison) {
    
    comparison_3yrpyrrole <- all_comb_sim_sum %>%
      dplyr::filter(
        net_strategy_comp == "Pyrethroid-Pyrrole 3-year campaigns uncosted",
        net_strategy_int == "Pyrethroid-Pyrrole 3-year campaigns uncosted"
      ) %>%
      dplyr::select(
        fs_area,
        clin_cases_all_ages_comp_med,
        clin_cases_all_ages_comp_lo,
        clin_cases_all_ages_comp_hi,
        clin_cases_under5_p100kU5_comp_med,
        clin_cases_under5_p100kU5_comp_lo,
        clin_cases_under5_p100kU5_comp_hi,
        pfpr_182_1824_mean_comp_med,
        pfpr_182_1824_mean_comp_lo,
        pfpr_182_1824_mean_comp_hi,
        pfpr_730_3649_mean_comp_med,
        pfpr_730_3649_mean_comp_lo,
        pfpr_730_3649_mean_comp_hi,
        pfpr_0_36499_mean_comp_med,
        pfpr_0_36499_mean_comp_lo,
        pfpr_0_36499_mean_comp_hi
      ) %>%
      dplyr::distinct() %>%
      dplyr::rename(
        clin_cases_all_ages_3yrpyrrole_med = clin_cases_all_ages_comp_med,
        clin_cases_all_ages_3yrpyrrole_lo  = clin_cases_all_ages_comp_lo,
        clin_cases_all_ages_3yrpyrrole_hi  = clin_cases_all_ages_comp_hi,
        clin_cases_under5_p100kU5_3yrpyrrole_med = clin_cases_under5_p100kU5_comp_med,
        clin_cases_under5_p100kU5_3yrpyrrole_lo  = clin_cases_under5_p100kU5_comp_lo,
        clin_cases_under5_p100kU5_3yrpyrrole_hi  = clin_cases_under5_p100kU5_comp_hi,
        pfpr_182_1824_mean_3yrpyrrole_med = pfpr_182_1824_mean_comp_med,
        pfpr_182_1824_mean_3yrpyrrole_lo  = pfpr_182_1824_mean_comp_lo,
        pfpr_182_1824_mean_3yrpyrrole_hi  = pfpr_182_1824_mean_comp_hi,
        pfpr_730_3649_mean_3yrpyrrole_med = pfpr_730_3649_mean_comp_med,
        pfpr_730_3649_mean_3yrpyrrole_lo  = pfpr_730_3649_mean_comp_lo,
        pfpr_730_3649_mean_3yrpyrrole_hi  = pfpr_730_3649_mean_comp_hi,
        pfpr_0_36499_mean_3yrpyrrole_med = pfpr_0_36499_mean_comp_med,
        pfpr_0_36499_mean_3yrpyrrole_lo  = pfpr_0_36499_mean_comp_lo,
        pfpr_0_36499_mean_3yrpyrrole_hi  = pfpr_0_36499_mean_comp_hi
      )
    
    all_comb_sim_sum <- all_comb_sim_sum %>%
      dplyr::left_join(comparison_3yrpyrrole, by = "fs_area")
    
  }
  
  message("Comparison summary calculated. Assigning to ", sum_comparison_name)
  assign(sum_comparison_name, all_comb_sim_sum, envir = .GlobalEnv)
  
}

summarise_full_comparison()

write.csv(all_ann_data_sum, "all_ann_data_sum.csv", row.names = FALSE)

all_ann_data_sum <- all_ann_data_sum %>%
  dplyr::mutate(
    country_name = countrycode::countrycode(ISO2, "iso2c", "country.name")
  )





all_ann_data_sum_no_MW <- subset(all_ann_data_sum, ISO2 != "MW")


# eLife


clin_plt_width <- 11
clin_plt_height <- 8.5




cfp_plt_1 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpyrrole_med",
  x_lo = "clin_cases_all_ages_3yrpyrrole_lo",
  x_hi = "clin_cases_all_ages_3yrpyrrole_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-Chlorfenapyr campaigns"),
  x_breaks = seq(0,1000,100),
  x_limits = c(0, 850),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-175, 550),
  xsf = 1/100,
  reference_lines = c(0,10,20,30,40,50,60,75,100,150,300)
)
ggsave(filename = "clin_cfp_cases.pdf", plot = cfp_plt_1, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

cfp_plt_1_BF <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpyrrole_med",
  x_lo = "clin_cases_all_ages_3yrpyrrole_lo",
  x_hi = "clin_cases_all_ages_3yrpyrrole_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-Chlorfenapyr campaigns"),
  alpha_val = 0.6,
  x_breaks = seq(0,1000,100),
  x_limits = c(200, 750),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 550),
  xsf = 1/100,
  reference_lines = c(0,10,20,30,40,50,60,75,100,150,250,500),
  iso2 = "BF"
)
ggsave(filename = "clin_cfp_cases_change_BF.pdf", plot = cfp_plt_1_BF,
       width = 8, height = clin_plt_height, units = "in",
       device = cairo_pdf)

cfp_plt_1_GH <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpyrrole_med",
  x_lo = "clin_cases_all_ages_3yrpyrrole_lo",
  x_hi = "clin_cases_all_ages_3yrpyrrole_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-Chlorfenapyr campaigns"),
  alpha_val = 0.6,
  x_breaks = seq(0,1000,100),
  x_limits = c(150, 700),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 550),
  xsf = 1/100,
  reference_lines = c(0,10,20,30,40,50,60,75,100,150),
  iso2 = "GH"
)
ggsave(filename = "clin_cfp_cases_change_GH.pdf", plot = cfp_plt_1_GH,
       width = 8, height = clin_plt_height, units = "in",
       device = cairo_pdf)

cfp_plt_1_ML <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpyrrole_med",
  x_lo = "clin_cases_all_ages_3yrpyrrole_lo",
  x_hi = "clin_cases_all_ages_3yrpyrrole_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-Chlorfenapyr campaigns"),
  alpha_val = 0.6,
  x_breaks = seq(0,1000,100),
  x_limits = c(0, 700),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 550),
  xsf = 1/100,
  reference_lines = c(0,10,20,30,40,50,60,75,100,150,250,500),
  iso2 = "ML"
)
ggsave(filename = "clin_cfp_cases_change_ML.pdf", plot = cfp_plt_1_ML,
       width = 9, height = clin_plt_height, units = "in",
       device = cairo_pdf)

cfp_plt_1_MZ <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpyrrole_med",
  x_lo = "clin_cases_all_ages_3yrpyrrole_lo",
  x_hi = "clin_cases_all_ages_3yrpyrrole_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-Chlorfenapyr campaigns"),
  alpha_val = 0.6,
  x_breaks = seq(0,1000,100),
  x_limits = c(0, 700),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 550),
  xsf = 1/100,
  reference_lines = c(0,10,20,30,40,50,60,75,100,150,250,500),
  iso2 = "MZ"
)
ggsave(filename = "clin_cfp_cases_change_MZ.pdf", plot = cfp_plt_1_MZ,
       width = 9, height = clin_plt_height, units = "in",
       device = cairo_pdf)

cfp_plt_1_SN <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpyrrole_med",
  x_lo = "clin_cases_all_ages_3yrpyrrole_lo",
  x_hi = "clin_cases_all_ages_3yrpyrrole_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-Chlorfenapyr campaigns"),
  alpha_val = 0.6,
  x_breaks = seq(0,1000,50),
  x_limits = c(0, 225),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 425),
  xsf = 1/100,
  reference_lines = c(0,50,100,150,200,300,500,1000),
  iso2 = "SN"
)
ggsave(filename = "clin_cfp_cases_change_SN.pdf", plot = cfp_plt_1_SN,
       width = clin_plt_width, height = clin_plt_height, units = "in",
       device = cairo_pdf)

cfp_plt_2 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "no future nets",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpyrrole_med",
  x_lo = "clin_cases_all_ages_3yrpyrrole_lo",
  x_hi = "clin_cases_all_ages_3yrpyrrole_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-Chlorfenapyr campaigns"),
  xsf = 1/100,
  x_breaks = seq(0,1000,100),
  x_limits = c(0, 850),
  y_breaks = seq(0,1000,100),
  y_limits = c(0, 750),
  flip_y = FALSE
)
ggsave(filename = "clin_cfp_cases_avert.pdf", plot = cfp_plt_2, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

cfp_plt_3 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "pfpr_0_36499_mean_3yrpyrrole_med",
  x_lo = "pfpr_0_36499_mean_3yrpyrrole_lo",
  x_hi = "pfpr_0_36499_mean_3yrpyrrole_hi",
  x_label = expression("Mean " * italic(Pf) * "PR with " *
                         "triennial Pyrethroid-Chlorfenapyr campaigns (%)"),
  x_breaks = seq(0,100,5),
  x_limits = c(0,37),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 550),
  xsf = 100
)
ggsave(filename = "clin_cfp_pfpr.pdf", plot = cfp_plt_3, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

cfp_plt_4 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "mean_u_3yr_med",
  x_lo = "mean_u_3yr_lo",
  x_hi = "mean_u_3yr_hi",
  x_label = paste("Mean use under a triennial strategy (%)"),
  x_breaks = seq(0,100,5),
  x_limits = c(23, 62.8),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 550),
  xsf = 100
)
ggsave(filename = "clin_cfp_meanu.pdf", plot = cfp_plt_4, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

cfp_plt_5 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "mean_a_3yr_med",
  x_lo = "mean_a_3yr_lo",
  x_hi = "mean_a_3yr_hi",
  x_label = paste("Mean access under a triennial strategy (%)"),
  x_breaks = seq(0,100,5),
  x_limits = c(35,77.5),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 550),
  xsf = 100
)
ggsave(filename = "clin_cfp_meana.pdf", plot = cfp_plt_5, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

cfp_plt_6 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "mean_uga_3yr_med",
  x_lo = "mean_uga_3yr_lo",
  x_hi = "mean_uga_3yr_hi",
  x_label = paste("Mean use given access under a triennial strategy (%)"),
  x_breaks = seq(0,100,5),
  x_limits = c(40, 95),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 550),
  xsf = 100
)
ggsave(filename = "clin_cfp_meanuga.pdf", plot = cfp_plt_6, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

cfp_plt_7 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "ret_u_int_med",
  x_lo = "ret_u_int_lo",
  x_hi = "ret_u_int_hi",
  x_label = paste("Mean duration of use (months)"),
  x_breaks = seq(0,100,3),
  x_limits = c(9, 42),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 550)
  )
ggsave(filename = "clin_cfp_retu.pdf", plot = cfp_plt_7, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

cfp_plt_8 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "ret_a_int_med",
  x_lo = "ret_a_int_lo",
  x_hi = "ret_a_int_hi",
  x_label = paste("Mean duration of access (months)"),
  x_breaks = seq(0,100,3),
  x_limits = c(12, 45),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 550)
)
ggsave(filename = "clin_cfp_reta.pdf", plot = cfp_plt_8, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

cfp_plt_9 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-Pyrrole 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "pyrethroid_resistance",
  x_label = paste("Pyrethroid resistance (%)"),
  x_breaks = seq(0,100,5),
  x_limits = c(25, 95),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-200, 550),
  xsf = 100
)
ggsave(filename = "clin_cfp_res.pdf", plot = cfp_plt_9, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)



clin_plt_width <- 9
clin_plt_height <- 7

pbo_plt_1 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpbo_med",
  x_lo = "clin_cases_all_ages_3yrpbo_lo",
  x_hi = "clin_cases_all_ages_3yrpbo_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-PBO campaigns"),
  x_breaks = seq(0,1000,100),
  x_limits = c(0, 850),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500),
  xsf = 1/100,
  reference_lines = c(0,10,20,30,40,50,60,75,100,150,300)
)
ggsave(filename = "clin_pbo_cases.pdf", plot = pbo_plt_1, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

pbo_plt_1_BF <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpbo_med",
  x_lo = "clin_cases_all_ages_3yrpbo_lo",
  x_hi = "clin_cases_all_ages_3yrpbo_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-PBO campaigns"),
  alpha_val = 0.6,
  x_breaks = seq(0,1000,100),
  x_limits = c(0, 850),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500),
  xsf = 1/100,
  reference_lines = c(0,10,20,30,40,50,60,75,100,150,300),
  iso2 = "BF",
  roman_labels = TRUE
)
# ggsave(filename = "clin_pbo_cases_change_BF.pdf", plot = pbo_plt_1_BF,
#        width = 8, height = clin_plt_height, units = "in",
#        device = cairo_pdf)
BF_quad_cases <- combine_plots_vertical(BF_quad, pbo_plt_1_BF,
                                        legend_from_first = TRUE,
                                        plot_heights = c(0.6,1))
ggsave(filename = "BF_quad_pbo_cases_change.pdf", plot = BF_quad_cases,
       width = 9, height = 12, units = "in",
       device = cairo_pdf)

pbo_plt_1_GH <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpbo_med",
  x_lo = "clin_cases_all_ages_3yrpbo_lo",
  x_hi = "clin_cases_all_ages_3yrpbo_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-PBO campaigns"),
  alpha_val = 0.6,
  x_breaks = seq(0,1000,100),
  x_limits = c(0, 850),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500),
  xsf = 1/100,
  reference_lines = c(0,10,20,30,40,50,60,75,100,150,300),
  iso2 = "GH",
  roman_labels = TRUE
)
# ggsave(filename = "clin_pbo_cases_change_GH.pdf", plot = pbo_plt_1_GH,
#        width = clin_plt_width, height = clin_plt_height, units = "in",
#        device = cairo_pdf)
GH_quad_cases <- combine_plots_vertical(GH_quad, pbo_plt_1_GH,
                                        legend_from_first = TRUE,
                                        plot_heights = c(0.6,1))
ggsave(filename = "GH_quad_pbo_cases_change.pdf", plot = GH_quad_cases,
       width = 9, height = 12, units = "in",
       device = cairo_pdf)

pbo_plt_1_ML <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpbo_med",
  x_lo = "clin_cases_all_ages_3yrpbo_lo",
  x_hi = "clin_cases_all_ages_3yrpbo_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-PBO campaigns"),
  alpha_val = 0.6,
  x_breaks = seq(0,1000,100),
  x_limits = c(0, 850),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500),
  xsf = 1/100,
  reference_lines = c(0,10,20,30,40,50,60,75,100,150,300),
  iso2 = "ML",
  roman_labels = TRUE
)
# ggsave(filename = "clin_pbo_cases_change_ML.pdf", plot = pbo_plt_1_ML,
#        width = clin_plt_width, height = clin_plt_height, units = "in",
#        device = cairo_pdf)
ML_quad_cases <- combine_plots_vertical(ML_quad, pbo_plt_1_ML,
                                        legend_from_first = TRUE,
                                        plot_heights = c(0.6,1))
ggsave(filename = "ML_quad_pbo_cases_change.pdf", plot = ML_quad_cases,
       width = 9, height = 12, units = "in",
       device = cairo_pdf)

pbo_plt_1_MZ <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpbo_med",
  x_lo = "clin_cases_all_ages_3yrpbo_lo",
  x_hi = "clin_cases_all_ages_3yrpbo_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-PBO campaigns"),
  alpha_val = 0.6,
  x_breaks = seq(0,1000,100),
  x_limits = c(0, 850),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500),
  xsf = 1/100,
  reference_lines = c(0,10,20,30,40,50,60,75,100,150,300),
  iso2 = "MZ",
  roman_labels = TRUE
)
# ggsave(filename = "clin_pbo_cases_change_MZ.pdf", plot = pbo_plt_1_MZ,
#        width = clin_plt_width, height = clin_plt_height, units = "in",
#        device = cairo_pdf)
MZ_quad_cases <- combine_plots_vertical(MZ_quad, pbo_plt_1_MZ,
                                        legend_from_first = TRUE,
                                        plot_heights = c(0.6,1))
ggsave(filename = "MZ_quad_pbo_cases_change.pdf", plot = MZ_quad_cases,
       width = 9, height = 12, units = "in",
       device = cairo_pdf)

pbo_plt_1_SN <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpbo_med",
  x_lo = "clin_cases_all_ages_3yrpbo_lo",
  x_hi = "clin_cases_all_ages_3yrpbo_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-PBO campaigns"),
  alpha_val = 0.6,
  x_breaks = seq(0,1000,50),
  x_limits = c(0, 250),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-250, 450),
  xsf = 1/100,
  reference_lines = c(0,50,100,150,200,300,500),
  iso2 = "SN",
  roman_labels = TRUE
)
# ggsave(filename = "clin_pbo_cases_change_SN.pdf", plot = pbo_plt_1_SN,
#        width = clin_plt_width, height = clin_plt_height, units = "in",
#        device = cairo_pdf)
SN_quad_cases <- combine_plots_vertical(SN_quad, pbo_plt_1_SN,
                                        legend_from_first = TRUE,
                                        plot_heights = c(0.6,1))
ggsave(filename = "SN_quad_pbo_cases_change.pdf", plot = SN_quad_cases,
       width = 9, height = 12, units = "in",
       device = cairo_pdf)


pbo_plt_2 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "no future nets",
  age_group_label = "all-age",
  x = "clin_cases_all_ages_3yrpbo_med",
  x_lo = "clin_cases_all_ages_3yrpbo_lo",
  x_hi = "clin_cases_all_ages_3yrpbo_hi",
  x_label = paste("Mean annual clinical cases per 1,000 people with triennial",
                  "Pyrethroid-PBO campaigns"),
  xsf = 1/100,
  x_breaks = seq(0,1000,100),
  x_limits = c(0, 850),
  y_breaks = seq(0,1000,100),
  y_limits = c(0, 750),
  flip_y = FALSE
)
ggsave(filename = "clin_pbo_cases_avert.pdf", plot = pbo_plt_2, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

pbo_plt_3 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "pfpr_0_36499_mean_3yrpbo_med",
  x_lo = "pfpr_0_36499_mean_3yrpbo_lo",
  x_hi = "pfpr_0_36499_mean_3yrpbo_hi",
  x_label = expression("Mean P" * italic(f) * "PR with " *
                         "triennial Pyrethroid-PBO campaigns (%)"),
  x_breaks = seq(0,100,5),
  x_limits = c(0,37),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500),
  xsf = 100
)
ggsave(filename = "clin_pbo_pfpr.pdf", plot = pbo_plt_3, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

pbo_plt_4 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "mean_u_3yr_med",
  x_lo = "mean_u_3yr_lo",
  x_hi = "mean_u_3yr_hi",
  x_label = paste("Mean use under a triennial strategy (%)"),
  x_breaks = seq(0,100,5),
  x_limits = c(23, 62.8),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500),
  xsf = 100
)
ggsave(filename = "clin_pbo_meanu.pdf", plot = pbo_plt_4, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

pbo_plt_5 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "mean_a_3yr_med",
  x_lo = "mean_a_3yr_lo",
  x_hi = "mean_a_3yr_hi",
  x_label = paste("Mean access under a triennial strategy (%)"),
  x_breaks = seq(0,100,5),
  x_limits = c(35,77.5),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500),
  xsf = 100
)
ggsave(filename = "clin_pbo_meana.pdf", plot = pbo_plt_5, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

pbo_plt_6 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "mean_uga_3yr_med",
  x_lo = "mean_uga_3yr_lo",
  x_hi = "mean_uga_3yr_hi",
  x_label = paste("Mean use given access under a triennial strategy (%)"),
  x_breaks = seq(0,100,5),
  x_limits = c(40, 95),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500),
  xsf = 100
)
ggsave(filename = "clin_pbo_meanuga.pdf", plot = pbo_plt_6, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

pbo_plt_7 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "ret_u_int_med",
  x_lo = "ret_u_int_lo",
  x_hi = "ret_u_int_hi",
  x_label = paste("Mean duration of use (months)"),
  x_breaks = seq(0,100,3),
  x_limits = c(9, 42),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500)
)
ggsave(filename = "clin_pbo_retu.pdf", plot = pbo_plt_7, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

pbo_plt_8 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "ret_a_int_med",
  x_lo = "ret_a_int_lo",
  x_hi = "ret_a_int_hi",
  x_label = paste("Mean duration of access (months)"),
  x_breaks = seq(0,100,3),
  x_limits = c(12, 45),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500)
)
ggsave(filename = "clin_pbo_reta.pdf", plot = pbo_plt_8, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)

pbo_plt_9 <- all_ann_data_sum_no_MW %>% make_clin_cases_plot(
  comparison_name = "Pyrethroid-PBO 3-year campaigns uncosted",
  age_group_label = "all-age",
  x = "pyrethroid_resistance",
  x_label = paste("Pyrethroid resistance (%)"),
  x_breaks = seq(0,100,10),
  x_limits = c(25, 95),
  y_breaks = seq(-1000,1000,100),
  y_limits = c(-325, 500),
  xsf = 100
)
ggsave(filename = "clin_pbo_res.pdf", plot = pbo_plt_9, width = clin_plt_width,
       height = clin_plt_height, units = "in", device = cairo_pdf)



combine_plots <- function(...,
                          legend_from = 1,
                          override_legends = character(0),
                          plot_widths = NULL,
                          plot_heights = NULL) {
  plots <- list(...)
  n_plots <- length(plots)
  
  if (legend_from < 1 || legend_from > n_plots) {
    stop("legend_from must be between 1 and the number of plots")
  }
  
  all_aes <- c("colour", "fill", "shape", "linetype", "size", "alpha", "stroke")
  guides_to_remove <- setdiff(all_aes, override_legends)
  
  # Remove legends from all plots except the one selected
  plots <- lapply(seq_along(plots), function(i) {
    p <- plots[[i]]
    if (i != legend_from) {
      for (aesthetic in guides_to_remove) {
        p <- p + do.call(guides, setNames(list("none"), aesthetic))
      }
    }
    p
  })
  
  wrap_plots(plots,
             nrow = 1,
             guides = "collect",
             widths = plot_widths,
             heights = plot_heights) +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")
}


ret_quad_case_wd <- 10
ret_quad_case_ht <- 14

BF_ret_quad <- combine_plots(BF_retention, BF_quad, legend_from = 2)
ggsave(filename = "BF_ret_quad.pdf", plot = BF_ret_quad, width = 10,
       height = 5, units = "in", device = cairo_pdf)
BF_ret_quad_cases <- combine_plots_vertical(BF_ret_quad, pbo_plt_1_BF,
                                            plot_heights = c(0.5,1))
ggsave(filename = "BF_ret_quad_cases.pdf", plot = BF_ret_quad_cases,
       width = ret_quad_case_wd, height = ret_quad_case_ht,
       units = "in", device = cairo_pdf)

GH_ret_quad <- combine_plots(GH_retention, GH_quad, legend_from = 2)
ggsave(filename = "GH_ret_quad.pdf", plot = GH_ret_quad, width = 10,
       height = 5, units = "in", device = cairo_pdf)
GH_ret_quad_cases <- combine_plots_vertical(GH_ret_quad, pbo_plt_1_GH,
                                            plot_heights = c(0.5,1))
ggsave(filename = "GH_ret_quad_cases.pdf", plot = GH_ret_quad_cases,
       width = ret_quad_case_wd*1.1, height = ret_quad_case_ht*1.1,
       units = "in", device = cairo_pdf)

MW_ret_quad <- combine_plots(MW_retention, MW_quad, legend_from = 2)
ggsave(filename = "MW_ret_quad.pdf", plot = MW_ret_quad, width = 10,
       height = 5, units = "in", device = cairo_pdf)
# MW_ret_quad_cases <- combine_plots_vertical(MW_ret_quad, pbo_plt_1_MW,
#                                             plot_heights = c(0.5,1))
# ggsave(filename = "MW_ret_quad_cases.pdf", plot = MW_ret_quad_cases, width = 9,
#        height = 12, units = "in", device = cairo_pdf)

ML_ret_quad <- combine_plots(ML_retention, ML_quad, legend_from = 2)
ggsave(filename = "ML_ret_quad.pdf", plot = ML_ret_quad, width = 10,
       height = 5, units = "in", device = cairo_pdf)
ML_ret_quad_cases <- combine_plots_vertical(ML_ret_quad, pbo_plt_1_ML,
                                            plot_heights = c(0.5,1))
ggsave(filename = "ML_ret_quad_cases.pdf", plot = ML_ret_quad_cases,
       width = ret_quad_case_wd, height = ret_quad_case_ht,
       units = "in", device = cairo_pdf)

MZ_ret_quad <- combine_plots(MZ_retention, MZ_quad, legend_from = 2)
ggsave(filename = "MZ_ret_quad.pdf", plot = MZ_ret_quad, width = 10,
       height = 5, units = "in", device = cairo_pdf)
MZ_ret_quad_cases <- combine_plots_vertical(MZ_ret_quad, pbo_plt_1_MZ,
                                            plot_heights = c(0.5,1))
ggsave(filename = "MZ_ret_quad_cases.pdf", plot = MZ_ret_quad_cases,
       width = ret_quad_case_wd, height = ret_quad_case_ht,
       units = "in", device = cairo_pdf)

SN_ret_quad <- combine_plots(SN_retention, SN_quad, legend_from = 2)
ggsave(filename = "SN_ret_quad.pdf", plot = SN_ret_quad, width = 10,
       height = 5, units = "in", device = cairo_pdf)
SN_ret_quad_cases <- combine_plots_vertical(SN_ret_quad, pbo_plt_1_SN,
                                            plot_heights = c(0.6,1))
ggsave(filename = "SN_ret_quad_cases.pdf", plot = SN_ret_quad_cases,
       width = ret_quad_case_wd, height = ret_quad_case_ht,
       units = "in", device = cairo_pdf)


library(patchwork)

combine_plots <- function(...,
                          legend_from_first = TRUE,
                          override_legends = character(0),
                          plot_widths = NULL,
                          plot_heights = NULL) {
  plots <- list(...)
  
  all_aes <- c("colour", "fill", "shape", "linetype", "size", "alpha", "stroke")
  guides_to_remove <- setdiff(all_aes, override_legends)
  
  if (legend_from_first && length(plots) > 1) {
    plots[-1] <- lapply(plots[-1], function(p) {
      for (aesthetic in guides_to_remove) {
        p <- p + do.call(guides, setNames(list("none"), aesthetic))
      }
      p
    })
  }
  
  wrap_plots(plots,
             nrow = 1,
             guides = "collect",
             widths = plot_widths,
             heights = plot_heights) +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")
}



combine_plots_vertical <- function(...,
                                   legend_from_first = TRUE,
                                   override_legends = character(0),
                                   plot_widths = NULL,
                                   plot_heights = NULL) {
  plots <- list(...)
  
  all_aes <- c("colour", "fill", "shape", "linetype", "size", "alpha", "stroke")
  guides_to_remove <- setdiff(all_aes, override_legends)
  
  if (legend_from_first && length(plots) > 1) {
    plots[-1] <- lapply(plots[-1], function(p) {
      # Remove guides by explicitly setting them to "none"
      for (aesthetic in guides_to_remove) {
        p <- p + do.call(guides, setNames(list("none"), aesthetic))
      }
      p
    })
  }
  
  wrap_plots(plots,
             ncol = 1,
             guides = "collect",
             widths = plot_widths,
             heights = plot_heights) +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "right")
}



# e.g.
# eir_age_plt <- combine_plots_vertical(plt_eir_all_age, plt_eir_u5)
# ggsave(filename = "eir_age_plt.pdf", plot = eir_age_plt, width = 9,
#        height = 11, units = "in", device = cairo_pdf)



#-------------------------------------------------------------------------------

library(dplyr)
library(tidyr)

# Variables to process
vars <- c("ret_u", "ret_a", "mean_u", "mean_a", "mean_uga", "D_u",
          "C0_u", "P0_u", "prop_routine_dist")

# Region ID columns
region_vars <- c("ISO2", "fs_name_1", "urbanicity", "fs_area",
                 "fs_area_id", "new_area_id", "pop")

# Filtered dataset
filtered_data <- all_ann_data %>%
  filter(net_strategy == "Pyrethroid-only 3-year campaigns uncosted") %>%
  mutate(prop_routine_dist = adj_ann_routine_nets_dist / adj_ann_total_nets_dist,
         P0_u = D_u + C0_u)

# =============================
# 1. REGION-LEVEL (UNWEIGHTED)
# =============================

# 1. Compute one value per region × sample_id × variable
region_sample_summary <- vars %>%
  lapply(function(v) {
    filtered_data %>%
      select(all_of(region_vars), sample_id, !!sym(v)) %>%
      rename(value = !!sym(v)) %>%
      mutate(variable = v)
  }) %>%
  bind_rows()

# 2. Summarise across samples to get mean and CrI
region_summary <- region_sample_summary %>%
  group_by(across(all_of(region_vars)), variable) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    lwr = quantile(value, 0.025, na.rm = TRUE),
    upr = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )



# =============================
# 2. COUNTRY-LEVEL (WEIGHTED)
# =============================

# Compute weighted mean per ISO2 and sample_id for each variable
country_sample_summary <- vars %>%
  lapply(function(v) {
    filtered_data %>%
      group_by(ISO2, sample_id) %>%
      summarise(
        value = sum(.data[[v]] * pop, na.rm = TRUE) / sum(pop, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(variable = v)
  }) %>%
  bind_rows()

# Then summarise across samples
country_weighted_summary <- country_sample_summary %>%
  group_by(ISO2, variable) %>%
  summarise(
    mean_weighted = mean(value),
    lwr = quantile(value, 0.025),
    upr = quantile(value, 0.975),
    .groups = "drop"
  )


# =============================
# 3. OVERALL (WEIGHTED)
# =============================

# Compute weighted mean across all rows for each sample_id
overall_sample_summary <- vars %>%
  lapply(function(v) {
    filtered_data %>%
      group_by(sample_id) %>%
      summarise(
        value = sum(.data[[v]] * pop, na.rm = TRUE) / sum(pop, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(variable = v)
  }) %>%
  bind_rows()

# Summarise across the 100 samples
overall_weighted_summary <- overall_sample_summary %>%
  group_by(variable) %>%
  summarise(
    mean_weighted = mean(value),
    lwr = quantile(value, 0.025),
    upr = quantile(value, 0.975),
    .groups = "drop"
  )


# Filtered dataset
filtered_data2 <- all_ann_data %>%
  filter(net_strategy == "Pyrethroid-only 2-year campaigns uncosted") %>%
  mutate(prop_routine_dist = adj_ann_routine_nets_dist / adj_ann_total_nets_dist)


# =============================
# 1. REGION-LEVEL (UNWEIGHTED)
# =============================

# 1. Compute one value per region × sample_id × variable
region_sample_summary2 <- vars %>%
  lapply(function(v) {
    filtered_data2 %>%
      select(all_of(region_vars), sample_id, !!sym(v)) %>%
      rename(value = !!sym(v)) %>%
      mutate(variable = v)
  }) %>%
  bind_rows()

# 2. Summarise across samples to get mean and CrI
region_summary2 <- region_sample_summary2 %>%
  group_by(across(all_of(region_vars)), variable) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    lwr = quantile(value, 0.025, na.rm = TRUE),
    upr = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )



# =============================
# 2. COUNTRY-LEVEL (WEIGHTED)
# =============================

# Compute weighted mean per ISO2 and sample_id for each variable
country_sample_summary2 <- vars %>%
  lapply(function(v) {
    filtered_data2 %>%
      group_by(ISO2, sample_id) %>%
      summarise(
        value = sum(.data[[v]] * pop, na.rm = TRUE) / sum(pop, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(variable = v)
  }) %>%
  bind_rows()

# Then summarise across samples
country_weighted_summary2 <- country_sample_summary2 %>%
  group_by(ISO2, variable) %>%
  summarise(
    mean_weighted = mean(value),
    lwr = quantile(value, 0.025),
    upr = quantile(value, 0.975),
    .groups = "drop"
  )


# =============================
# 3. OVERALL (WEIGHTED)
# =============================

# Compute weighted mean across all rows for each sample_id
overall_sample_summary2 <- vars %>%
  lapply(function(v) {
    filtered_data2 %>%
      group_by(sample_id) %>%
      summarise(
        value = sum(.data[[v]] * pop, na.rm = TRUE) / sum(pop, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(variable = v)
  }) %>%
  bind_rows()

# Summarise across the 100 samples
overall_weighted_summary2 <- overall_sample_summary2 %>%
  group_by(variable) %>%
  summarise(
    mean_weighted = mean(value),
    lwr = quantile(value, 0.025),
    upr = quantile(value, 0.975),
    .groups = "drop"
  )


mean_uga_region_summary <- region_summary %>%
  filter(variable == "mean_uga")
dim(mean_uga_region_summary)
sum(mean_uga_region_summary$mean_val > 0.8)

mean_uga_region_summary_noMW <- region_summary %>%
  filter(variable == "mean_uga", ISO2 != "MW")
dim(mean_uga_region_summary_noMW)
sum(mean_uga_region_summary$mean_val > 0.8)

mean_uga_region_summary2 <- region_summary2 %>%
  filter(variable == "mean_uga")
dim(mean_uga_region_summary2)
sum(mean_uga_region_summary2$mean_val > 0.8)

mean_uga_region_summary2_noMW <- region_summary2 %>%
  filter(variable == "mean_uga", ISO2 != "MW")
dim(mean_uga_region_summary2_noMW)
sum(mean_uga_region_summary2$mean_val > 0.8)

# relative increase in use

rel_inc_u_df <- filtered_data %>%
  dplyr::select(fs_area_id, new_area_id, sample_id, pop, mean_u) %>%
  dplyr::rename(mean_u_3yr = mean_u) %>%
  dplyr::left_join(filtered_data2 %>%
                     dplyr::select(fs_area_id, new_area_id, sample_id, pop, mean_u)) %>%
  dplyr::rename(mean_u_2yr = mean_u) %>%
  dplyr::mutate(rel_inc_u = (mean_u_2yr - mean_u_3yr) / mean_u_3yr )
rel_inc_u_global <- rel_inc_u_df %>%
  filter(!is.na(pop), !is.na(rel_inc_u)) %>%
  group_by(sample_id) %>%
  summarise(
    weighted_mean_rel_inc_u = if (sum(pop) > 0) {
      sum(rel_inc_u * pop) / sum(pop)
    } else {
      NA_real_
    },
    .groups = "drop"
  )
median(rel_inc_u_global$weighted_mean_rel_inc_u)
quantile(rel_inc_u_global$weighted_mean_rel_inc_u, probs = c(0.025, 0.975))

rel_inc_a_df <- filtered_data %>%
  dplyr::select(fs_area_id, new_area_id, sample_id, pop, mean_a) %>%
  dplyr::rename(mean_a_3yr = mean_a) %>%
  dplyr::left_join(filtered_data2 %>%
                     dplyr::select(fs_area_id, new_area_id, sample_id, pop, mean_a)) %>%
  dplyr::rename(mean_a_2yr = mean_a) %>%
  dplyr::mutate(rel_inc_a = (mean_a_2yr - mean_a_3yr) / mean_a_3yr )
rel_inc_a_global <- rel_inc_a_df %>%
  filter(!is.na(pop), !is.na(rel_inc_a)) %>%
  group_by(sample_id) %>%
  summarise(
    weighted_mean_rel_inc_a = if (sum(pop) > 0) {
      sum(rel_inc_a * pop) / sum(pop)
    } else {
      NA_real_
    },
    .groups = "drop"
  )
median(rel_inc_a_global$weighted_mean_rel_inc_a)
quantile(rel_inc_a_global$weighted_mean_rel_inc_a, probs = c(0.025, 0.975))

rel_inc_uga_df <- filtered_data %>%
  dplyr::select(fs_area_id, new_area_id, sample_id, pop, mean_uga) %>%
  dplyr::rename(mean_uga_3yr = mean_uga) %>%
  dplyr::left_join(filtered_data2 %>%
                     dplyr::select(fs_area_id, new_area_id, sample_id, pop, mean_uga)) %>%
  dplyr::rename(mean_uga_2yr = mean_uga) %>%
  dplyr::mutate(rel_inc_uga = (mean_uga_2yr - mean_uga_3yr) / mean_uga_3yr )
rel_inc_uga_global <- rel_inc_uga_df %>%
  filter(!is.na(pop), !is.na(rel_inc_uga)) %>%
  group_by(sample_id) %>%
  summarise(
    weighted_mean_rel_inc_uga = if (sum(pop) > 0) {
      sum(rel_inc_uga * pop) / sum(pop)
    } else {
      NA_real_
    },
    .groups = "drop"
  )
median(rel_inc_uga_global$weighted_mean_rel_inc_uga)
quantile(rel_inc_uga_global$weighted_mean_rel_inc_uga, probs = c(0.025, 0.975))

#-------------------------------------------------------------------------------
# Proportion continuous

library(dplyr)
library(tidyr)
library(ggplot2)
library(countrycode)

# Prepare the data
plot_data <- region_summary %>%
  distinct(new_area_id, variable, mean_val, lwr, upr, .keep_all = TRUE) %>%
  filter(variable %in% c("D_u", "mean_u")) %>%
  select(ISO2, fs_area_id, urbanicity, pop, variable, mean_val, lwr, upr) %>%
  pivot_wider(
    names_from = variable,
    values_from = c(mean_val, lwr, upr),
    names_glue = "{.value}_{variable}"
  ) %>%
  mutate(country = countrycode(ISO2, origin = "iso2c", destination = "country.name"))

# Global regression
global_model <- lm(mean_val_D_u ~ 0 + mean_val_mean_u, data = plot_data, weights = pop)
global_slope <- coef(global_model)[[1]]
global_slope_bounds <- confint(global_model, level = 0.95)

# Country-specific slopes
country_slopes <- plot_data %>%
  # group_by(country) %>%
  # summarise(
  #   slope = coef(lm(mean_val_D_u ~ 0 + mean_val_mean_u, weights = pop))[1],
  #   .groups = "drop"
  # ) %>%
  group_by(country) %>%
  summarise(
    model = list(lm(mean_val_D_u ~ 0 + mean_val_mean_u, weights = pop)),
    .groups = "drop"
  ) %>%
  mutate(
    slope = purrr::map_dbl(model, ~ coef(.x)[[1]]),
    ci = purrr::map(model, ~ confint(.x)[1, ]),
    ci_lower = purrr::map_dbl(ci, 1),
    ci_upper = purrr::map_dbl(ci, 2)
  ) %>%
  select(-model, -ci) %>%
  mutate(
    ratio = slope / global_slope,
    label = paste0(
      "continuous = ",
      round(slope * 100, 1),
      "% × overall (",
      round(ci_lower * 100, 1),
      ", ",
      round(ci_upper * 100, 1),
      ")"
      ),
    angle = atan(slope) * 180 / pi  # angle in degrees
  )

# Label positions for regression line annotations (scaled)
label_positions <- country_slopes %>%
  mutate(
    x = 75,  # previously 0.7
    y = slope * x + if_else(country == "Burkina Faso", -1.2, 1.2),  # offset in percent units
    hjust = 0.5
  )

# Global regression annotation
global_label <- data.frame(
  country = "Overall",
  x = 75,
  y = global_slope * 75 + 1.2,
  label = paste0("continuous = ", round(global_slope * 100, 1), "% × overall (",
                 round(global_slope_bounds[1] * 100, 1), ", ",
                 round(global_slope_bounds[2] * 100, 1), ")"),
  angle = atan(global_slope) * 180 / pi
)

# Routine ITN labels (per country, in %)
routine_labels <- country_weighted_summary %>%
  filter(variable == "prop_routine_dist") %>%
  mutate(
    label = paste0(
      round(mean_weighted * 100, 1),
      "% (", round(lwr * 100, 1), ", ", round(upr * 100, 1), ")"
    )
  ) %>%
  left_join(plot_data %>% distinct(ISO2, country), by = "ISO2") %>%
  arrange(country) %>%
  mutate(
    x = 22,
    y = 58 - 2 * row_number()  # stacked vertically in percent space
  )

# Overall routine ITN label
overall_label <- overall_weighted_summary %>%
  filter(variable == "prop_routine_dist") %>%
  mutate(
    label = paste0(
      round(mean_weighted * 100, 1),
      "% (", round(lwr * 100, 1), ", ", round(upr * 100, 1), ")"
    ),
    x = 22,
    y = 58 - 2 * (nrow(routine_labels) + 1)
  )

# Heading label
routine_heading <- data.frame(
  x = 22,
  y = 58,
  label = "Continuous ITNs distributed"
)

# Plot
cont_prop_plt <- ggplot(
  plot_data,
  aes(
    x = mean_val_mean_u * 100,
    y = mean_val_D_u * 100,
    color = country
  )
) +
  geom_point(aes(size = pop, shape = urbanicity), alpha = 0.5) +
  geom_errorbarh(aes(xmin = lwr_mean_u * 100, xmax = upr_mean_u * 100), size = 0.5, alpha = 0.5) +
  geom_errorbar(aes(ymin = lwr_D_u * 100, ymax = upr_D_u * 100), size = 0.5, alpha = 0.5) +
  
  # Country regression lines
  geom_smooth(
    method = "lm",
    formula = y ~ 0 + x,
    aes(group = country, fill = country, weight = pop),
    linetype = "dashed",
    se = TRUE,
    alpha = 0.15,
    fullrange = TRUE
  ) +
  
  # Global regression line
  geom_smooth(
    method = "lm",
    formula = y ~ 0 + x,
    aes(weight = pop),
    linetype = "dashed",
    color = "black",
    fill = "black",
    se = TRUE,
    alpha = 0.15,
    fullrange = TRUE
  ) +
  
  # Country labels for regression lines
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = label, angle = angle, color = country),
    inherit.aes = FALSE,
    size = 2,
    show.legend = FALSE
  ) +
  
  # Global label
  geom_text(
    data = global_label,
    aes(x = x, y = y, label = label, angle = angle),
    color = "black",
    inherit.aes = FALSE,
    size = 2,
    show.legend = FALSE
  ) +
  
  # # Country routine ITN % labels
  # geom_text(
  #   data = routine_labels,
  #   aes(x = x, y = y, label = label, color = country),
  #   inherit.aes = FALSE,
  #   hjust = 0,
  #   size = 3.2,
  #   show.legend = FALSE
  # ) +
  # 
  # # Overall label
  # geom_text(
  #   data = overall_label,
  #   aes(x = x, y = y, label = label),
  #   inherit.aes = FALSE,
  #   hjust = 0,
  #   size = 3.2,
  #   color = "black",
  #   show.legend = FALSE
  # ) +
  # 
  # # Heading
  # geom_text(
  #   data = routine_heading,
  #   aes(x = x, y = y, label = label),
  #   inherit.aes = FALSE,
  #   hjust = 0,
  #   size = 3.4,
  #   fontface = "bold",
  #   color = "black"
  # ) +
  # 
  # Axis scales: now in percentages
  scale_x_continuous(
    breaks = seq(0,100,10),
    limits = c(20, 90)
  ) +
  scale_y_continuous(
    breaks = seq(0,100,10),
    limits = c(0, 60)
  ) +
  
  
  # Axis and legend labels
  labs(
    x = "Mean annual use of any ITN (%)",
    y = "Use of continuously-distributed ITNs (%)",
    color = "Country",
    fill = "Country",
    shape = "Urbanicity",
    size = "Population"
  ) +
  coord_fixed() +
  theme_bw()

# Save to PDF
ggsave("cont_prop_plt_v3.pdf", cont_prop_plt, width = 8, height = 6)








# Prepare the data
plot_start_data <- region_summary %>%
  distinct(new_area_id, variable, mean_val, lwr, upr, .keep_all = TRUE) %>%
  filter(variable %in% c("D_u", "P0_u")) %>%
  select(ISO2, fs_area_id, urbanicity, pop, variable, mean_val, lwr, upr) %>%
  pivot_wider(
    names_from = variable,
    values_from = c(mean_val, lwr, upr),
    names_glue = "{.value}_{variable}"
  ) %>%
  mutate(country = countrycode(ISO2, origin = "iso2c", destination = "country.name"))

# Global regression
global_model <- lm(mean_val_D_u ~ 0 + mean_val_P0_u, data = plot_start_data, weights = pop)
global_slope <- coef(global_model)[[1]]
global_slope_bounds <- confint(global_model, level = 0.95)

# Country-specific slopes
country_slopes <- plot_start_data %>%
  # group_by(country) %>%
  # summarise(
  #   slope = coef(lm(mean_val_D_u ~ 0 + mean_val_mean_u, weights = pop))[1],
  #   .groups = "drop"
  # ) %>%
  group_by(country) %>%
  summarise(
    model = list(lm(mean_val_D_u ~ 0 + mean_val_P0_u, weights = pop)),
    .groups = "drop"
  ) %>%
  mutate(
    slope = purrr::map_dbl(model, ~ coef(.x)[[1]]),
    ci = purrr::map(model, ~ confint(.x)[1, ]),
    ci_lower = purrr::map_dbl(ci, 1),
    ci_upper = purrr::map_dbl(ci, 2)
  ) %>%
  select(-model, -ci) %>%
  mutate(
    ratio = slope / global_slope,
    label = paste0(
      "continuous = ",
      round(slope * 100, 1),
      "% × overall (",
      round(ci_lower * 100, 1),
      ", ",
      round(ci_upper * 100, 1),
      ")"
    ),
    angle = atan(slope) * 180 / pi  # angle in degrees
  )

# Label positions for regression line annotations (scaled)
label_positions <- country_slopes %>%
  mutate(
    x = 125,  # previously 0.7
    y = slope * x + if_else(
      country == "Senegal",
      -1.2,
      ifelse(
        country == "Mozambique",
        2.4,
        ifelse(
          country == "Burkina Faso",
          1.6,
          1.2
          )
        )
      ),  # offset in percent units
    hjust = 0.5
  )

# Global regression annotation
global_label <- data.frame(
  country = "Overall",
  x = 125,
  y = global_slope * 125 + 1.2,
  label = paste0("continuous = ", round(global_slope * 100, 1), "% × overall (",
                 round(global_slope_bounds[1] * 100, 1), ", ",
                 round(global_slope_bounds[2] * 100, 1), ")"),
  angle = atan(global_slope) * 180 / pi
)

# Routine ITN labels (per country, in %)
routine_labels <- country_weighted_summary %>%
  filter(variable == "prop_routine_dist") %>%
  mutate(
    label = paste0(
      round(mean_weighted * 100, 1),
      "% (", round(lwr * 100, 1), ", ", round(upr * 100, 1), ")"
    )
  ) %>%
  left_join(plot_data %>% distinct(ISO2, country), by = "ISO2") %>%
  arrange(country) %>%
  mutate(
    x = 22,
    y = 58 - 2 * row_number()  # stacked vertically in percent space
  )

# Overall routine ITN label
overall_label <- overall_weighted_summary %>%
  filter(variable == "prop_routine_dist") %>%
  mutate(
    label = paste0(
      round(mean_weighted * 100, 1),
      "% (", round(lwr * 100, 1), ", ", round(upr * 100, 1), ")"
    ),
    x = 22,
    y = 58 - 2 * (nrow(routine_labels) + 1)
  )

# Heading label
routine_heading <- data.frame(
  x = 22,
  y = 58,
  label = "Continuous ITNs distributed"
)

# Plot
cont_prop_start_plt <- ggplot(
  plot_start_data,
  aes(
    x = mean_val_P0_u * 100,
    y = mean_val_D_u * 100,
    color = country
  )
) +
  geom_point(aes(size = pop, shape = urbanicity), alpha = 0.5) +
  geom_errorbarh(aes(xmin = lwr_P0_u * 100, xmax = upr_P0_u * 100), size = 0.5, alpha = 0.5) +
  geom_errorbar(aes(ymin = lwr_D_u * 100, ymax = upr_D_u * 100), size = 0.5, alpha = 0.5) +
  
  # Country regression lines
  geom_smooth(
    method = "lm",
    formula = y ~ 0 + x,
    aes(group = country, fill = country, weight = pop),
    linetype = "dashed",
    se = TRUE,
    alpha = 0.15,
    fullrange = TRUE
  ) +
  
  # Global regression line
  geom_smooth(
    method = "lm",
    formula = y ~ 0 + x,
    aes(weight = pop),
    linetype = "dashed",
    color = "black",
    fill = "black",
    se = TRUE,
    alpha = 0.15,
    fullrange = TRUE
  ) +
  
  # Country labels for regression lines
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = label, angle = angle, color = country),
    inherit.aes = FALSE,
    size = 2,
    show.legend = FALSE
  ) +
  
  # Global label
  geom_text(
    data = global_label,
    aes(x = x, y = y, label = label, angle = angle),
    color = "black",
    inherit.aes = FALSE,
    size = 2,
    show.legend = FALSE
  ) +
  
  # # Country routine ITN % labels
  # geom_text(
  #   data = routine_labels,
  #   aes(x = x, y = y, label = label, color = country),
  #   inherit.aes = FALSE,
  #   hjust = 0,
  #   size = 3.2,
  #   show.legend = FALSE
  # ) +
  # 
  # # Overall label
  # geom_text(
  #   data = overall_label,
  #   aes(x = x, y = y, label = label),
  #   inherit.aes = FALSE,
  #   hjust = 0,
  #   size = 3.2,
  #   color = "black",
  #   show.legend = FALSE
  # ) +
  # 
  # # Heading
  # geom_text(
  #   data = routine_heading,
  #   aes(x = x, y = y, label = label),
  #   inherit.aes = FALSE,
  #   hjust = 0,
  #   size = 3.4,
  #   fontface = "bold",
  #   color = "black"
  # ) +
  
  # Axis scales: now in percentages
  scale_x_continuous(
    breaks = seq(0,100,10),
    limits = c(40, 140)
  ) +
  scale_y_continuous(
    breaks = seq(0,100,10),
    limits = c(0, 60)
  ) +
  
  
  # Axis and legend labels
  labs(
    x = "Use of any ITN immediately after a campaign (%)",
    y = "Use of continuously-distributed ITNs (%)",
    color = "Country",
    fill = "Country",
    shape = "Urbanicity",
    size = "Population"
  ) +
  coord_fixed() +
  theme_bw()

# Save to PDF
ggsave("cont_prop_plt_v4.pdf", cont_prop_start_plt, width = 8, height = 6)


prop_cont_comb <- combine_plots_vertical(cont_prop_start_plt,
                                         cont_prop_plt)
ggsave(filename = "prop_cont_comb.pdf", plot = prop_cont_comb,
       width = 8, height = 10,
       units = "in", device = cairo_pdf)













summarize_and_plot_ctry_ann <- function(ctry_ann_sum_df, sf = 1/100) {
  df <- ctry_ann_sum_df
  
  df$clin_cases_all_ages <- df$clin_cases_all_ages * sf
  df$clin_cases_under5 <- df$clin_cases_under5 * sf
  df$clin_cases_5_14  <- df$clin_cases_5_14  * sf
  df$clin_cases_15plus <- df$clin_cases_15plus * sf
  df$routine <- df$adj_ann_routine_nets_dist * sf
  df$campaign <- df$adj_camp_nets_dist * sf
  
  df <- df %>%
    mutate(net_strategy = str_remove(net_strategy, " uncosted$")) %>%
    mutate(
      facet_group = case_when(
        net_strategy == "no future nets"             ~ "No ITNs",
        is.na(mass_int_yr)                           ~ "Routine distribution only",
        mass_int_yr == 2                              ~ "2-year campaigns + routine",
        mass_int_yr == 3                              ~ "3-year campaigns + routine",
        TRUE                                          ~ NA_character_
      ),
      facet_group = factor(facet_group, levels = c(
        "No ITNs",
        "Routine distribution only",
        "2-year campaigns + routine",
        "3-year campaigns + routine"
      ))
    )
  
  annot_df <- df %>%
    filter(net_name == "Pyrethroid-only") %>%
    group_by(facet_group) %>%
    summarise(
      mean_u_med  = median(mean_u, na.rm = TRUE),
      mean_u_lb   = quantile(mean_u, 0.025, na.rm = TRUE),
      mean_u_ub   = quantile(mean_u, 0.975, na.rm = TRUE),
      mean_a_med  = median(mean_a, na.rm = TRUE),
      mean_a_lb   = quantile(mean_a, 0.025, na.rm = TRUE),
      mean_a_ub   = quantile(mean_a, 0.975, na.rm = TRUE),
      uga_med     = median(mean_u / mean_a, na.rm = TRUE),
      uga_lb      = quantile(mean_u / mean_a, 0.025, na.rm = TRUE),
      uga_ub      = quantile(mean_u / mean_a, 0.975, na.rm = TRUE),
      routine_med = median(routine, na.rm = TRUE),
      camp_med = median(campaign, na.rm = TRUE),
      proc_med = median(adj_camp_proc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = if_else(
        facet_group == "Routine distribution only",
        sprintf(
          "Use: %d%% [%d, %d]\nAccess: %d%% [%d, %d]\nUse | Access: %d%% [%d, %d]\nAnnual routine ITNs: %d\u2030\n\n",
          round(mean_u_med * 100), round(mean_u_lb * 100), round(mean_u_ub * 100),
          round(mean_a_med * 100), round(mean_a_lb * 100), round(mean_a_ub * 100),
          round(uga_med * 100),    round(uga_lb * 100),    round(uga_ub * 100),
          round(routine_med)
        ),
        sprintf(
          "Use: %d%% [%d, %d]\nAccess: %d%% [%d, %d]\nUse | Access: %d%% [%d, %d]\nAnnual routine ITNs: %d\u2030\nCampaign ITNs: %d\u2030\nCampaign quantifier: %.1f",
          round(mean_u_med * 100), round(mean_u_lb * 100), round(mean_u_ub * 100),
          round(mean_a_med * 100), round(mean_a_lb * 100), round(mean_a_ub * 100),
          round(uga_med * 100),    round(uga_lb * 100),    round(uga_ub * 100),
          round(routine_med), round(camp_med), proc_med
        )
      )
    )
  
  
  y_max <- max(df$clin_cases_all_ages, na.rm = TRUE)
  
  return(
    ggplot(df, aes(x = net_strategy, y = clin_cases_all_ages, fill = net_name)) +
      geom_violin(color = NA, alpha = 0.3, trim = FALSE, scale = "width") +
      stat_summary(
        fun = median,
        geom = "point",
        shape = 21,
        size = 2,
        aes(color = net_name),
        show.legend = FALSE
      ) +
      stat_summary(
        fun.data = function(x) {
          data.frame(
            y = median(x),
            ymin = quantile(x, 0.025),
            ymax = quantile(x, 0.975)
          )
        },
        geom = "errorbar",
        width = 0.2,
        aes(color = net_name),
        show.legend = FALSE
      ) +
      scale_fill_discrete(na.translate = TRUE, labels = function(x) ifelse(is.na(x), "none", x)) +
      scale_color_discrete(na.translate = TRUE, labels = function(x) ifelse(is.na(x), "none", x)) +
      labs(
        #title = "All-Age Clinical Cases by Net Strategy",
        x = NULL,
        y = "All-age annual clinical cases per 1,000",
        fill = "Net Type"
      ) +
      facet_grid(. ~ facet_group, scales = "free_x", space = "free_x", switch = "x") +
      geom_text(
        data = annot_df,
        aes(
          x = 0.5,
          y = y_max * 0.95,
          label = label
        ),
        inherit.aes = FALSE,
        hjust = 0,
        size = 3.5
      ) +
      theme_bw() +
      theme(
        panel.spacing = unit(0, "pt"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
      )
  )
}



ML_plt <- summarize_and_plot_ctry_ann(ML_ctry_ann_sum)

ggsave("ML_plot.pdf", plot = ML_plt, width = 9, height = 6,
       units = "in", device = cairo_pdf)




plot_regional_violin_pages <- function(ann_data,
                                       urban_type = "urban",
                                       n_regions_per_page = 4,
                                       output_prefix = "region_plot",
                                       sf = 1 / 100,
                                       yr_sf = 6,
                                       fixed_y = TRUE,
                                       ymin = 0,
                                       ymax = 1000,
                                       show_legend = TRUE,
                                       y_tick_interval = 100,
                                       stretch_last_page = FALSE) {
  require(dplyr)
  require(ggplot2)
  require(stringr)
  require(ggtext)
  require(countrycode)
  
  # Get ISO2 and full country name
  iso2 <- unique(ann_data$ISO2)
  if (length(iso2) != 1) stop("More than one ISO2 code in dataset.")
  country_name <- countrycode(iso2, origin = "iso2c", destination = "country.name")
  
  # Step 1: Filter by urbanicity
  df <- ann_data %>%
    filter(urbanicity == urban_type)
  
  # Step 2: Compute MAP_pfpr per region and facet labels (with markdown)
  regions_ordered <- df %>%
    group_by(fs_area_id, fs_name_1, urbanicity) %>%
    summarise(MAP_pfpr = first(MAP_pfpr), .groups = "drop") %>%
    arrange(MAP_pfpr) %>%
    mutate(facet_row = sprintf(
      "**%s**<br>2024 *Pf*PR<sub>2–10</sub> = %.1f%%",
      fs_name_1, MAP_pfpr * 100
    ))
  
  # Step 3: Split fs_area_ids into batches
  region_batches <- split(regions_ordered$fs_area_id,
                          ceiling(seq_along(regions_ordered$fs_area_id) / n_regions_per_page))
  
  # Step 4: Loop through batches
  for (i in seq_along(region_batches)) {
    batch_ids <- region_batches[[i]]
    
    batch_df <- df %>%
      filter(fs_area_id %in% batch_ids) %>%
      left_join(regions_ordered[, c("fs_area_id", "facet_row", "MAP_pfpr")], by = "fs_area_id") %>%
      mutate(
        facet_row = factor(
          facet_row,
          levels = unique(regions_ordered$facet_row[regions_ordered$fs_area_id %in% batch_ids])
        ),
        net_strategy = str_remove(net_strategy, " uncosted$"),
        clin_cases_all_ages = clin_cases_all_ages * sf,
        routine = adj_ann_routine_nets_dist * sf,
        campaign = adj_ann_camp_nets_dist * sf,
        tot = adj_ann_total_nets_dist * sf * yr_sf
      ) %>%
      mutate(
        facet_group = case_when(
          net_strategy == "no future nets"             ~ "No ITNs",
          is.na(mass_int_yr)                           ~ "Routine distribution only",
          mass_int_yr == 3                             ~ "3-year campaigns + routine",
          mass_int_yr == 2                             ~ "2-year campaigns + routine",
          TRUE                                         ~ NA_character_
        ),
        facet_group = factor(facet_group, levels = c(
          "No ITNs",
          "Routine distribution only",
          "3-year campaigns + routine",
          "2-year campaigns + routine"
        ))
      )
    
    # Step 5: Generate annotation labels
    annot_df <- batch_df %>%
      filter(net_name == "Pyrethroid-only") %>%
      group_by(facet_row, facet_group) %>%
      summarise(
        mean_u_med  = median(mean_u, na.rm = TRUE),
        mean_u_lb   = quantile(mean_u, 0.025, na.rm = TRUE),
        mean_u_ub   = quantile(mean_u, 0.975, na.rm = TRUE),
        mean_a_med  = median(mean_a, na.rm = TRUE),
        mean_a_lb   = quantile(mean_a, 0.025, na.rm = TRUE),
        mean_a_ub   = quantile(mean_a, 0.975, na.rm = TRUE),
        uga_med     = median(mean_u / mean_a, na.rm = TRUE),
        uga_lb      = quantile(mean_u / mean_a, 0.025, na.rm = TRUE),
        uga_ub      = quantile(mean_u / mean_a, 0.975, na.rm = TRUE),
        routine_med = median(routine, na.rm = TRUE),
        camp_med    = median(campaign, na.rm = TRUE),
        proc_med    = median(adj_camp_proc, na.rm = TRUE),
        tot_med     = median(tot, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        label = if_else(
          facet_group == "Routine distribution only",
          sprintf(
            "\n\n\nUse: %d%% [%d, %d]\nAccess: %d%% [%d, %d]\nUse | Access: %d%% [%d, %d]\nITNs dist over %d yrs (\u2030): %d\n\n",
            round(mean_u_med * 100), round(mean_u_lb * 100), round(mean_u_ub * 100),
            round(mean_a_med * 100), round(mean_a_lb * 100), round(mean_a_ub * 100),
            round(uga_med * 100),    round(uga_lb * 100),    round(uga_ub * 100),
            yr_sf, round(tot_med)
          ),
          sprintf(
            "\n\nUse: %d%% [%d, %d]\nAccess: %d%% [%d, %d]\nUse | Access: %d%% [%d, %d]\nITNs dist over %d yrs (\u2030): %d\nCampaign quantifier: %.1f",
            round(mean_u_med * 100), round(mean_u_lb * 100), round(mean_u_ub * 100),
            round(mean_a_med * 100), round(mean_a_lb * 100), round(mean_a_ub * 100),
            round(uga_med * 100),    round(uga_lb * 100),    round(uga_ub * 100),
            yr_sf, round(tot_med), proc_med
          )
        )
      )
    
    # Step 6: Create the plot
    p <- ggplot(batch_df, aes(x = net_strategy, y = clin_cases_all_ages, fill = net_name)) +
      geom_violin(color = NA, alpha = 0.3, trim = FALSE, scale = "width") +
      stat_summary(
        fun = median,
        geom = "point",
        shape = 21,
        size = 2,
        aes(color = net_name),
        show.legend = FALSE
      ) +
      stat_summary(
        fun.data = function(x) {
          data.frame(
            y = median(x),
            ymin = quantile(x, 0.025),
            ymax = quantile(x, 0.975)
          )
        },
        geom = "errorbar",
        width = 0.2,
        aes(color = net_name),
        show.legend = FALSE
      ) +
      geom_text(
        data = annot_df,
        aes(x = 0.5, y = ymax * 0.95, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        size = 3
      ) +
      scale_fill_discrete(na.translate = TRUE, labels = function(x) ifelse(is.na(x), "none", x)) +
      scale_color_discrete(na.translate = TRUE, guide = "none") +
      guides(fill = guide_legend(title = "Net Type")) +
      scale_y_continuous(breaks = seq(ymin, ymax, by = y_tick_interval)) +
      labs(
        title = paste(country_name, "-", urban_type, "areas", i),
        x = NULL,
        y = "All-age mean annual clinical cases per 1,000",
        fill = "Net Type"
      ) +
      facet_grid(
        rows = vars(facet_row),
        cols = vars(facet_group),
        scales = "free_x",
        space = "free_x",
        switch = "x",
        labeller = label_value
      ) +
      theme_bw() +
      theme(
        panel.spacing = unit(0, "pt"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y = ggtext::element_markdown(size = 11),
        strip.text.x = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal"
      )
    
    if (fixed_y) {
      p <- p + coord_cartesian(ylim = c(ymin, ymax))
    }
    
    # Adjust height based on padding toggle
    n_actual_rows <- length(unique(batch_df$facet_row))
    height_per_row <- 11.5 / n_regions_per_page
    plot_height <- if (stretch_last_page || n_actual_rows == n_regions_per_page) {
      11.5
    } else {
      height_per_row * n_actual_rows + 1
    }
    
    filename <- sprintf("%s_%s_%s_page%02d.pdf", iso2, output_prefix, urban_type, i)
    ggsave(filename, plot = p, width = 8.5, height = plot_height, units = "in", device = cairo_pdf)
  }
  
  # Step 7: Combine individual plots into one multi-page PDF
  combined_file <- sprintf("%s_%s_%s_ALL_PAGES.pdf", iso2, output_prefix, urban_type)
  individual_files <- sprintf("%s_%s_%s_page%02d.pdf", iso2, output_prefix, urban_type, seq_along(region_batches))
  
  if (requireNamespace("pdftools", quietly = TRUE)) {
    pdftools::pdf_combine(input = individual_files, output = combined_file)
  } else {
    warning("The 'pdftools' package is not installed. Cannot combine individual PDFs into a multi-page file.")
  }
}


plot_regional_violin_pages(BF_ann_data, urban_type = "urban", ymax = 1200)
plot_regional_violin_pages(BF_ann_data, urban_type = "rural", ymax = 1200)
plot_regional_violin_pages(ML_ann_data, urban_type = "urban", ymax = 1200)
plot_regional_violin_pages(ML_ann_data, urban_type = "rural", ymax = 1200)

plot_regional_violin_pages(SN_ann_data, urban_type = "urban", ymax = 400)




create_short_df <- function(XX) {
  # Construct full data frame name (e.g., "BF_ann_data")
  full_df_name <- paste0(XX, "_ann_data")
  
  # Get the actual data frame from the global environment
  df <- get(full_df_name, envir = .GlobalEnv)
  
  # Define columns to extract
  cols_to_keep <- c(
    "country", "ISO2", "fs_name_1", "urbanicity", "fs_area", "fs_area_id", 
    "new_area_id", "EIR_urep_fit", "pop", "MAP_pfpr", "net_strategy", 
    "net_name", "mass_int_yr", "no_future_nets", "sample_id", 
    "adj_ann_routine_nets_dist", "adj_ann_camp_nets_dist", 
    "adj_ann_total_nets_dist", "adj_camp_proc", "mean_u", "mean_a", 
    "mean_use_given_access", "n_age_all_ages", "n_age_under5", 
    "n_age_5_14", "n_age_15plus", "cases_all_ages", "cases_under5", 
    "cases_5_14", "cases_15plus", "clin_cases_all_ages", 
    "clin_cases_under5", "clin_cases_5_14", "clin_cases_15plus", 
    "sev_cases_all_ages", "sev_cases_under5", "sev_cases_5_14", 
    "sev_cases_15plus", "pfpr_0_36499_mean", "pfpr_182_1824_mean", 
    "pfpr_730_3649_mean"
  )
  
  # Subset
  short_df <- df[, cols_to_keep, drop = FALSE]
  
  # Define new names
  new_names <- c(
    "country", "ISO2", "fs_name_1", "urbanicity", "fs_area", "fs_area_id", 
    "new_area_id", "EIR_urep_fit", "pop", "pfpr_MAP_2024", "net_strategy", 
    "net_name", "mass_int_yr", "no_future_nets", "sample_id", 
    "avg_ann_routine_nets_dist", "yr1_campaign_nets_dist", 
    "avg_ann_total_nets_dist", "camp_quantifier", "mean_use", "mean_access", 
    "mean_use_given_access", "n_age_all_ages", "n_age_under5", 
    "n_age_5_14", "n_age_15plus", "cases_all_ages", "cases_under5", 
    "cases_5_14", "cases_15plus", "clin_cases_all_ages", 
    "clin_cases_under5", "clin_cases_5_14", "clin_cases_15plus", 
    "sev_cases_all_ages", "sev_cases_under5", "sev_cases_5_14", 
    "sev_cases_15plus", "pfpr_all_ages", "pfpr_6_59mo", 
    "pfpr_2_10"
  )
  colnames(short_df) <- new_names
  
  # Assign to global environment
  assign(paste0(XX, "_ann_data_short"), short_df, envir = .GlobalEnv)
}


create_short_df("BF")
create_short_df("ML")



summarise_short_df <- function(XX) {
  # Load short data
  df_name <- paste0(XX, "_ann_data_short")
  df <- get(df_name, envir = .GlobalEnv)
  
  # Drop sample_id
  df <- df %>% select(-sample_id)
  
  # Variables to summarise
  vars_to_summarise <- c(
    "avg_ann_routine_nets_dist", "yr1_campaign_nets_dist", 
    "avg_ann_total_nets_dist", "camp_quantifier", "mean_use", "mean_access", 
    "mean_use_given_access", "n_age_all_ages", "n_age_under5", 
    "n_age_5_14", "n_age_15plus", "cases_all_ages", "cases_under5", 
    "cases_5_14", "cases_15plus", "clin_cases_all_ages", 
    "clin_cases_under5", "clin_cases_5_14", "clin_cases_15plus", 
    "sev_cases_all_ages", "sev_cases_under5", "sev_cases_5_14", 
    "sev_cases_15plus", "pfpr_all_ages", "pfpr_6_59mo", 
    "pfpr_2_10"
  )
  
  # Grouping variables
  group_vars <- c("fs_area_id", "net_strategy")
  
  # Metadata vars = everything else (excluding group vars + summary vars)
  meta_vars <- setdiff(names(df), c(vars_to_summarise, group_vars))
  
  # Summarise numeric variables
  df_sum <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      across(all_of(vars_to_summarise), list(
        lo = ~quantile(.x, 0.025, na.rm = TRUE),
        med = ~median(.x, na.rm = TRUE),
        hi = ~quantile(.x, 0.975, na.rm = TRUE)
      )),
      .groups = "drop"
    )
  
  # Get representative metadata per group
  meta_unique <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(across(all_of(meta_vars), ~first(.x)), .groups = "drop")
  
  # Join and save
  final_df <- left_join(meta_unique, df_sum, by = group_vars)
  assign(paste0(XX, "_ann_data_sum"), final_df, envir = .GlobalEnv)
}



summarise_short_df("BF")
summarise_short_df("ML")


create_vshort_sum_df <- function(XX) {
  # Construct input and output names
  input_name <- paste0(XX, "_ann_data_sum")
  output_name <- paste0(XX, "_ann_data_vshort_sum")
  
  # Load input data
  df <- get(input_name, envir = .GlobalEnv)
  
  # Define the column order
  cols_to_keep <- c(
    "fs_area_id", "net_strategy", "country", "ISO2", "fs_name_1", "urbanicity", "fs_area", 
    "new_area_id", "EIR_urep_fit", "pop", "pfpr_MAP_2024", "net_name", "mass_int_yr", "no_future_nets",
    
    "avg_ann_routine_nets_dist_lo", "avg_ann_routine_nets_dist_med", "avg_ann_routine_nets_dist_hi",
    "yr1_campaign_nets_dist_lo", "yr1_campaign_nets_dist_med", "yr1_campaign_nets_dist_hi",
    "avg_ann_total_nets_dist_lo", "avg_ann_total_nets_dist_med", "avg_ann_total_nets_dist_hi",
    "camp_quantifier_lo", "camp_quantifier_med", "camp_quantifier_hi",
    "mean_use_lo", "mean_use_med", "mean_use_hi",
    "mean_access_lo", "mean_access_med", "mean_access_hi",
    "mean_use_given_access_lo", "mean_use_given_access_med", "mean_use_given_access_hi",
    
    "n_age_under5_lo", "n_age_under5_med", "n_age_under5_hi",
    
    "cases_all_ages_lo", "cases_all_ages_med", "cases_all_ages_hi",
    "cases_under5_lo", "cases_under5_med", "cases_under5_hi",
    
    "clin_cases_all_ages_lo", "clin_cases_all_ages_med", "clin_cases_all_ages_hi",
    "clin_cases_under5_lo", "clin_cases_under5_med", "clin_cases_under5_hi",
    
    "sev_cases_all_ages_lo", "sev_cases_all_ages_med", "sev_cases_all_ages_hi",
    "sev_cases_under5_lo", "sev_cases_under5_med", "sev_cases_under5_hi",
    
    "pfpr_all_ages_lo", "pfpr_all_ages_med", "pfpr_all_ages_hi",
    "pfpr_6_59mo_lo", "pfpr_6_59mo_med", "pfpr_6_59mo_hi",
    "pfpr_2_10_lo", "pfpr_2_10_med", "pfpr_2_10_hi"
  )
  
  # Subset and reorder
  df_vshort <- df[, cols_to_keep, drop = FALSE]
  
  # Save to global environment
  assign(output_name, df_vshort, envir = .GlobalEnv)
}

create_vshort_sum_df("BF")
create_vshort_sum_df("ML")

write.csv(BF_ann_data_vshort_sum, "BF_ann_data_vshort_sum.csv", row.names = FALSE)
write.csv(ML_ann_data_vshort_sum, "ML_ann_data_vshort_sum.csv", row.names = FALSE)

write.csv(data.frame(Variable = names(BF_ann_data_vshort_sum)), 
          "Description_XX_ann_data_vshort_sum.csv", 
          row.names = FALSE)




generate_multiannual_summary <- function(iso2) {
  message("Generating annual summary data for ", iso2, "...")
  
  sim_data_name <- paste0(iso2, "_sim_data")
  ann_data_name <- paste0(iso2, "_multiann_data")
  
  if (!exists(sim_data_name, envir = .GlobalEnv)) {
    stop("Data object ", sim_data_name,
         " does not exist in the global environment.")
  }
  
  sim_data <- get(sim_data_name, envir = .GlobalEnv)
  
  message("Grouping and summarising across years...")
  
  ann_data <- sim_data %>%
    group_by(sample_id, fs_area_id, net_strategy, year_id) %>%
    summarise(
      across(c(ISO2, fs_name_1, urbanicity, fs_area, net_name,
               mass_int_yr, biennial_reduction, net_costings,
               no_future_nets, routine_baseline, CMC_start,
               CMC_end, ann_routine_nets_dist, ann_camp_nets_dist,
               adj_ann_routine_nets_dist, adj_ann_camp_nets_dist,
               budget_pc, country, iso3c, new_area_id, EIR_urep_fit,
               C0_u, D_u, invlam_u, lam_u, scale_factor, mean_u,
               C0_a, D_a, invlam_a, lam_a, mean_a),
             ~ first(.x), .names = "{.col}"),
      across(c(starts_with("n_age_"), starts_with("n_detect_lm_"),
               starts_with("cases_"), starts_with("clin_cases_"),
               starts_with("sev_cases_"), starts_with("pfpr_")),
             ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
      across(c("pop", "par", "par_pf", "par_pv",
               "adj_ann_total_nets_dist", "mean_use_given_access",
               "ctry_pop",
               "pop_weight"
      ),
      ~ mean(.x, na.rm = TRUE), .names = "{.col}")
    ) %>%
    ungroup()
  
  message("Summary complete. Assigning to ", ann_data_name)
  assign(ann_data_name, ann_data, envir = .GlobalEnv)
}



generate_multiannual_summary("BF")
append_site_file_data("BF", data_obj = "_multiann_data")


batch_df <- BF_multiann_data %>% filter(urbanicity == "rural")

batch_extra_df <- batch_df %>% dplyr::filter(year_id == 6)
batch_extra_df$year_id <- batch_extra_df$year_id + 1
batch_df <- rbind.data.frame(batch_df, batch_extra_df)

sf <- 1/100

batch_df <- batch_df %>%
  mutate(
    net_strategy = str_remove(net_strategy, " uncosted$"),
    clin_cases_all_ages = clin_cases_all_ages * sf,
    routine = adj_ann_routine_nets_dist * sf,
    campaign = adj_ann_camp_nets_dist * sf,
    tot = adj_ann_total_nets_dist * sf
  ) %>%
  mutate(
    facet_group = case_when(
      net_strategy == "no future nets"             ~ "No ITNs",
      is.na(mass_int_yr)                           ~ "Routine distribution only",
      mass_int_yr == 3                             ~ "3-year campaigns + routine",
      mass_int_yr == 2                             ~ "2-year campaigns + routine",
      TRUE                                         ~ NA_character_
    ),
    facet_group = factor(facet_group, levels = c(
      "No ITNs",
      "Routine distribution only",
      "3-year campaigns + routine",
      "2-year campaigns + routine"
    ))
  ) %>%
  group_by(year_id, fs_area, urbanicity, net_name, mass_int_yr) %>%
  summarise(
    clin_cases_all_ages_med = median(
      clin_cases_all_ages, na.rm = TRUE
    ),
    clin_cases_all_ages_lo = quantile(
      clin_cases_all_ages, 0.025, na.rm = TRUE
    ),
    clin_cases_all_ages_hi = quantile(
      clin_cases_all_ages, 0.975, na.rm = TRUE
    ),
    clin_cases_under5_med = median(
      clin_cases_under5, na.rm = TRUE
    ),
    clin_cases_under5_lo = quantile(
      clin_cases_under5, 0.025, na.rm = TRUE
    ),
    clin_cases_under5_hi = quantile(
      clin_cases_under5, 0.975, na.rm = TRUE
    ),
    cases_all_ages_med = median(
      cases_all_ages, na.rm = TRUE
    ),
    cases_all_ages_lo = quantile(
      cases_all_ages, 0.025, na.rm = TRUE
    ),
    cases_all_ages_hi = quantile(
      cases_all_ages, 0.975, na.rm = TRUE
    ),
    cases_under5_med = median(
      cases_under5, na.rm = TRUE
    ),
    cases_under5_lo = quantile(
      cases_under5, 0.025, na.rm = TRUE
    ),
    cases_under5_hi = quantile(
      cases_under5, 0.975, na.rm = TRUE
    ),
    .groups = "drop"
  )


p <- ggplot(batch_df %>% filter(mass_int_yr %in% c(2,3)),
            aes(x = year_id - 1,
                y = clin_cases_under5_med,
                colour = net_name,
                linetype = as.factor(mass_int_yr)
            )
) +
  geom_step() +
  theme_bw() +
  facet_wrap(facets = vars(fs_area), ncol = 5)
  
print(p)
  
  
  
  geom_step(color = NA, alpha = 0.3, trim = FALSE, scale = "width") +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 2,
    aes(color = net_name),
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data = function(x) {
      data.frame(
        y = median(x),
        ymin = quantile(x, 0.025),
        ymax = quantile(x, 0.975)
      )
    },
    geom = "errorbar",
    width = 0.2,
    aes(color = net_name),
    show.legend = FALSE
  ) +
  geom_text(
    data = annot_df,
    aes(x = 0.5, y = ymax * 0.95, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3
  ) +
  scale_fill_discrete(na.translate = TRUE, labels = function(x) ifelse(is.na(x), "none", x)) +
  scale_color_discrete(na.translate = TRUE, guide = "none") +
  guides(fill = guide_legend(title = "Net Type")) +
  scale_y_continuous(breaks = seq(ymin, ymax, by = y_tick_interval)) +
  labs(
    title = paste(country_name, "-", urban_type, "areas", i),
    x = NULL,
    y = "All-age mean annual clinical cases per 1,000",
    fill = "Net Type"
  ) +
  facet_grid(
    rows = vars(facet_row),
    cols = vars(facet_group),
    scales = "free_x",
    space = "free_x",
    switch = "x",
    labeller = label_value
  ) +
  theme_bw() +
  theme(
    panel.spacing = unit(0, "pt"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y = ggtext::element_markdown(size = 11),
    strip.text.x = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )



















plot_regional_violin_pages <- function(ann_data,
                                       urban_type = "urban",
                                       n_regions_per_page = 4,
                                       output_prefix = "region_plot",
                                       sf = 1 / 100,
                                       yr_sf = 6,
                                       fixed_y = TRUE,
                                       ymin = 0,
                                       ymax = 1000,
                                       show_legend = TRUE,
                                       y_tick_interval = 100,
                                       stretch_last_page = FALSE) {
  require(dplyr)
  require(ggplot2)
  require(stringr)
  require(ggtext)
  require(countrycode)
  
  # Get ISO2 and full country name
  iso2 <- unique(ann_data$ISO2)
  if (length(iso2) != 1) stop("More than one ISO2 code in dataset.")
  country_name <- countrycode(iso2, origin = "iso2c", destination = "country.name")
  
  # Step 1: Filter by urbanicity
  df <- ann_data %>%
    filter(urbanicity == urban_type)
  
  # Step 2: Compute MAP_pfpr per region and facet labels (with markdown)
  regions_ordered <- df %>%
    group_by(fs_area_id, fs_name_1, urbanicity) %>%
    summarise(MAP_pfpr = first(MAP_pfpr), .groups = "drop") %>%
    arrange(MAP_pfpr) %>%
    mutate(facet_row = sprintf(
      "**%s**<br>2024 *Pf*PR<sub>2–10</sub> = %.1f%%",
      fs_name_1, MAP_pfpr * 100
    ))
  
  # Step 3: Split fs_area_ids into batches
  region_batches <- split(regions_ordered$fs_area_id,
                          ceiling(seq_along(regions_ordered$fs_area_id) / n_regions_per_page))
  
  # Step 4: Loop through batches
  for (i in seq_along(region_batches)) {
    batch_ids <- region_batches[[i]]
    
    batch_df <- df %>%
      filter(fs_area_id %in% batch_ids) %>%
      left_join(regions_ordered[, c("fs_area_id", "facet_row", "MAP_pfpr")], by = "fs_area_id") %>%
      mutate(
        facet_row = factor(
          facet_row,
          levels = unique(regions_ordered$facet_row[regions_ordered$fs_area_id %in% batch_ids])
        ),
        net_strategy = str_remove(net_strategy, " uncosted$"),
        clin_cases_all_ages = clin_cases_all_ages * sf,
        routine = adj_ann_routine_nets_dist * sf,
        campaign = adj_ann_camp_nets_dist * sf,
        tot = adj_ann_total_nets_dist * sf * yr_sf
      ) %>%
      mutate(
        facet_group = case_when(
          net_strategy == "no future nets"             ~ "No ITNs",
          is.na(mass_int_yr)                           ~ "Routine distribution only",
          mass_int_yr == 3                             ~ "3-year campaigns + routine",
          mass_int_yr == 2                             ~ "2-year campaigns + routine",
          TRUE                                         ~ NA_character_
        ),
        facet_group = factor(facet_group, levels = c(
          "No ITNs",
          "Routine distribution only",
          "3-year campaigns + routine",
          "2-year campaigns + routine"
        ))
      ) %>%
      group_by(year_id) %>%
      summarise(
        clin_case_ages_all_ages_med = median(
          clin_case_ages_all_ages, na.rm = TRUE
        ),
        clin_case_ages_all_ages_lo = quantile(
          clin_case_ages_all_ages, 0.025, na.rm = TRUE
        ),
        clin_case_ages_all_ages_hi = quantile(
          clin_case_ages_all_ages, 0.975, na.rm = TRUE
        ),
        clin_case_ages_under5_med = median(
          clin_case_ages_under5, na.rm = TRUE
        ),
        clin_case_ages_all_ages_lo = quantile(
          clin_case_ages_under5, 0.025, na.rm = TRUE
        ),
        clin_case_ages_all_ages_hi = quantile(
          clin_case_ages_under5, 0.975, na.rm = TRUE
        ),
        .groups = "drop"
      )
    
    # # Step 5: Generate annotation labels
    # annot_df <- batch_df %>%
    #   filter(net_name == "Pyrethroid-only") %>%
    #   group_by(facet_row, facet_group) %>%
      # summarise(
      #   mean_u_med  = median(mean_u, na.rm = TRUE),
      #   mean_u_lb   = quantile(mean_u, 0.025, na.rm = TRUE),
      #   mean_u_ub   = quantile(mean_u, 0.975, na.rm = TRUE),
      #   mean_a_med  = median(mean_a, na.rm = TRUE),
      #   mean_a_lb   = quantile(mean_a, 0.025, na.rm = TRUE),
      #   mean_a_ub   = quantile(mean_a, 0.975, na.rm = TRUE),
      #   uga_med     = median(mean_u / mean_a, na.rm = TRUE),
      #   uga_lb      = quantile(mean_u / mean_a, 0.025, na.rm = TRUE),
      #   uga_ub      = quantile(mean_u / mean_a, 0.975, na.rm = TRUE),
      #   routine_med = median(routine, na.rm = TRUE),
      #   camp_med    = median(campaign, na.rm = TRUE),
      #   proc_med    = median(adj_camp_proc, na.rm = TRUE),
      #   tot_med     = median(tot, na.rm = TRUE),
      #   .groups = "drop"
      # ) %>%
    #   mutate(
    #     label = if_else(
    #       facet_group == "Routine distribution only",
    #       sprintf(
    #         "\n\n\nUse: %d%% [%d, %d]\nAccess: %d%% [%d, %d]\nUse | Access: %d%% [%d, %d]\nITNs dist over %d yrs (\u2030): %d\n\n",
    #         round(mean_u_med * 100), round(mean_u_lb * 100), round(mean_u_ub * 100),
    #         round(mean_a_med * 100), round(mean_a_lb * 100), round(mean_a_ub * 100),
    #         round(uga_med * 100),    round(uga_lb * 100),    round(uga_ub * 100),
    #         yr_sf, round(tot_med)
    #       ),
    #       sprintf(
    #         "\n\nUse: %d%% [%d, %d]\nAccess: %d%% [%d, %d]\nUse | Access: %d%% [%d, %d]\nITNs dist over %d yrs (\u2030): %d\nCampaign quantifier: %.1f",
    #         round(mean_u_med * 100), round(mean_u_lb * 100), round(mean_u_ub * 100),
    #         round(mean_a_med * 100), round(mean_a_lb * 100), round(mean_a_ub * 100),
    #         round(uga_med * 100),    round(uga_lb * 100),    round(uga_ub * 100),
    #         yr_sf, round(tot_med), proc_med
    #       )
    #     )
    #   )
    
    # Step 6: Create the plot
    p <- ggplot(batch_df,
                aes(x = year_id - 1,
                    y = clin_cases_all_ages_med,
                    fill = net_name)
                ) +
      geom_violin(color = NA, alpha = 0.3, trim = FALSE, scale = "width") +
      stat_summary(
        fun = median,
        geom = "point",
        shape = 21,
        size = 2,
        aes(color = net_name),
        show.legend = FALSE
      ) +
      stat_summary(
        fun.data = function(x) {
          data.frame(
            y = median(x),
            ymin = quantile(x, 0.025),
            ymax = quantile(x, 0.975)
          )
        },
        geom = "errorbar",
        width = 0.2,
        aes(color = net_name),
        show.legend = FALSE
      ) +
      geom_text(
        data = annot_df,
        aes(x = 0.5, y = ymax * 0.95, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        size = 3
      ) +
      scale_fill_discrete(na.translate = TRUE, labels = function(x) ifelse(is.na(x), "none", x)) +
      scale_color_discrete(na.translate = TRUE, guide = "none") +
      guides(fill = guide_legend(title = "Net Type")) +
      scale_y_continuous(breaks = seq(ymin, ymax, by = y_tick_interval)) +
      labs(
        title = paste(country_name, "-", urban_type, "areas", i),
        x = NULL,
        y = "All-age mean annual clinical cases per 1,000",
        fill = "Net Type"
      ) +
      facet_grid(
        rows = vars(facet_row),
        cols = vars(facet_group),
        scales = "free_x",
        space = "free_x",
        switch = "x",
        labeller = label_value
      ) +
      theme_bw() +
      theme(
        panel.spacing = unit(0, "pt"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y = ggtext::element_markdown(size = 11),
        strip.text.x = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal"
      )
    
    if (fixed_y) {
      p <- p + coord_cartesian(ylim = c(ymin, ymax))
    }
    
    # Adjust height based on padding toggle
    n_actual_rows <- length(unique(batch_df$facet_row))
    height_per_row <- 11.5 / n_regions_per_page
    plot_height <- if (stretch_last_page || n_actual_rows == n_regions_per_page) {
      11.5
    } else {
      height_per_row * n_actual_rows + 1
    }
    
    filename <- sprintf("%s_%s_%s_page%02d.pdf", iso2, output_prefix, urban_type, i)
    ggsave(filename, plot = p, width = 8.5, height = plot_height, units = "in", device = cairo_pdf)
  }
  
  # Step 7: Combine individual plots into one multi-page PDF
  combined_file <- sprintf("%s_%s_%s_ALL_PAGES.pdf", iso2, output_prefix, urban_type)
  individual_files <- sprintf("%s_%s_%s_page%02d.pdf", iso2, output_prefix, urban_type, seq_along(region_batches))
  
  if (requireNamespace("pdftools", quietly = TRUE)) {
    pdftools::pdf_combine(input = individual_files, output = combined_file)
  } else {
    warning("The 'pdftools' package is not installed. Cannot combine individual PDFs into a multi-page file.")
  }
}

















































# 
# 
# # Example usage:
# compile_country_sim_data(
#   base_path = "data_results/urep_06_20250604112237/summary",
#   start_year = 2025,
#   iso2 = "BF"
# )
# 
# append_use_given_access("BF")
# 
# append_summary_cases("BF")
# 
# generate_pop_summary("BF")
# 
# generate_country_pop_summary("BF")
# 
# append_ctry_pop_and_weights("BF")
# 
# append_real_case_counts("BF")               # Uses default sim_pop = 100000
# 
# append_year_weights("BF")
# 
# 
# generate_annual_summary("BF")
# 









generate_sim_summary <- function(iso2_code) {
  # Dynamically get the data object (e.g. ML_sim_data)
  sim_data_name <- paste0(iso2_code, "_sim_data")
  sim_data <- get(sim_data_name, envir = .GlobalEnv)
  
  # Summarise annual mean per region/sample (unweighted)
  sim_sum_unweighted <- sim_data %>%
    group_by(
      ISO2, fs_name_1, urbanicity, fs_area, fs_area_id,
      net_name, mass_int_yr, biennial_reduction, net_costings,
      no_future_nets, routine_baseline, sample_id, net_strategy, budget_pc
    ) %>%
    summarise(
      across(
        .cols = matches("(_tot$|_mean$)"),
        .fns  = mean,
        na.rm = TRUE
      ),
      adj_ann_routine_nets_dist = mean(adj_ann_routine_nets_dist, na.rm = TRUE),
      adj_ann_camp_nets_dist    = mean(adj_ann_camp_nets_dist, na.rm = TRUE),
      pop = mean(pop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      adj_ann_total_nets_dist = adj_ann_routine_nets_dist + adj_ann_camp_nets_dist
    )
  
  # Compute population weights for each region/sample/year
  weights_df <- sim_data %>%
    group_by(sample_id, fs_area_id) %>%
    summarise(
      mean_region_pop = mean(pop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(sample_id) %>%
    mutate(
      country_pop = sum(mean_region_pop, na.rm = TRUE),
      pop_weight  = mean_region_pop / country_pop
    ) %>%
    select(sample_id, fs_area_id, pop_weight)
  
  # Join weights and compute weighted means
  sim_sum_weighted <- sim_data %>%
    left_join(weights_df, by = c("sample_id", "fs_area_id")) %>%
    group_by(
      ISO2, fs_name_1, urbanicity, fs_area, fs_area_id,
      net_name, mass_int_yr, biennial_reduction, net_costings,
      no_future_nets, routine_baseline, sample_id, net_strategy, budget_pc
    ) %>%
    summarise(
      weighted_EIR_urep_fit              = sum(EIR_urep_fit * pop_weight, na.rm = TRUE),
      weighted_mean_u                    = sum(mean_u * pop_weight, na.rm = TRUE),
      weighted_mean_a                    = sum(mean_a * pop_weight, na.rm = TRUE),
      weighted_adj_ann_routine_nets_dist = sum(adj_ann_routine_nets_dist * pop_weight, na.rm = TRUE),
      weighted_adj_ann_camp_nets_dist    = sum(adj_ann_camp_nets_dist * pop_weight, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      weighted_adj_ann_total_nets_dist = weighted_adj_ann_routine_nets_dist + weighted_adj_ann_camp_nets_dist
    )
  
  # Final merge
  sim_sum <- sim_sum_unweighted %>%
    left_join(sim_sum_weighted,
              by = c(
                "ISO2", "fs_name_1", "urbanicity", "fs_area", "fs_area_id",
                "net_name", "mass_int_yr", "biennial_reduction", "net_costings",
                "no_future_nets", "routine_baseline", "sample_id", "net_strategy", "budget_pc"
              ))
  
  # Assign dynamically
  assign(paste0(iso2_code, "_sim_sum"), sim_sum, envir = .GlobalEnv)
}



generate_sim_summary("BF")
# → Creates BF_sim_sum in your global environment








ml_data_sum_pop <- ml_data_pop %>%
  group_by(
    ISO2, fs_name_1, urbanicity, fs_area, fs_area_id,
    net_name, mass_int_yr, biennial_reduction, net_costings,
    no_future_nets, routine_baseline, sample_id, net_strategy, budget_pc
  ) %>%
  summarise(
    across(
      .cols = matches("(_tot$|_mean$)"),  # existing epidemiological measures
      .fns = mean, na.rm = TRUE
    ),
    EIR_urep_fit = mean(EIR_urep_fit, na.rm = TRUE),
    mean_u       = mean(mean_u, na.rm = TRUE),
    mean_a       = mean(mean_a, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    clin_cases_all_ages = rowSums(across(63:78), na.rm = TRUE),
    clin_cases_under5   = rowSums(across(63:65), na.rm = TRUE)
  )



ml_data_realpop <- ml_data_pop %>%
  mutate(across(
    .cols = 22:101,
    .fns = ~ (.x / 100000) * pop
  )) %>%
  rename_with(
    .cols = 22:101,
    .fn = ~ str_replace(.x, "_tot$", "_pop")
  )

ml_data_sum_realpop <- ml_data_realpop %>%
  group_by(
    ISO2, fs_name_1, urbanicity, fs_area, fs_area_id,
    net_name, mass_int_yr, biennial_reduction, net_costings,
    no_future_nets, routine_baseline, sample_id, net_strategy, budget_pc
  ) %>%
  summarise(
    across(
      .cols = matches("(_pop$|_mean$)"),
      .fns = mean,
      na.rm = TRUE
    ),
    pop = sum(pop, na.rm = TRUE),  # keep total pop for weighting
    .groups = "drop"
  ) %>%
  left_join(
    ml_data_realpop %>%
      group_by(
        ISO2, fs_name_1, urbanicity, fs_area, fs_area_id,
        net_name, mass_int_yr, biennial_reduction, net_costings,
        no_future_nets, routine_baseline, sample_id, net_strategy, budget_pc
      ) %>%
      summarise(
        EIR_urep_fit = sum(EIR_urep_fit * pop, na.rm = TRUE) / sum(pop, na.rm = TRUE),
        mean_u       = sum(mean_u * pop, na.rm = TRUE) / sum(pop, na.rm = TRUE),
        mean_a       = sum(mean_a * pop, na.rm = TRUE) / sum(pop, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c(
      "ISO2", "fs_name_1", "urbanicity", "fs_area", "fs_area_id",
      "net_name", "mass_int_yr", "biennial_reduction", "net_costings",
      "no_future_nets", "routine_baseline", "sample_id", "net_strategy", "budget_pc"
    )
  ) %>%
  mutate(
    clin_cases_all_ages = rowSums(across(63:78), na.rm = TRUE),
    clin_cases_under5   = rowSums(across(63:65), na.rm = TRUE)
  )

                                  





# 
# ml_data_realpop <- ml_data_pop %>%
#   mutate(across(
#     .cols = 22:101,
#     .fns = ~ (.x / 100000) * pop
#   )) %>%
#   rename_with(
#     .cols = 22:101,
#     .fn = ~ str_replace(.x, "_tot$", "_pop")
#   )
# 
# ml_data_sum_realpop <- ml_data_realpop %>%
#   group_by(
#     ISO2, fs_name_1, urbanicity, fs_area, fs_area_id,
#     net_name, mass_int_yr, biennial_reduction, net_costings,
#     no_future_nets, routine_baseline, sample_id, net_strategy, budget_pc
#   ) %>%
#   summarise(
#     across(
#       .cols = matches("(_tot$|_mean$)"),  # existing epidemiological measures
#       .fns = mean, na.rm = TRUE
#     ),
#     EIR_urep_fit = mean(EIR_urep_fit, na.rm = TRUE),
#     mean_u       = mean(mean_u, na.rm = TRUE),
#     mean_a       = mean(mean_a, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     clin_cases_all_ages = rowSums(across(63:78), na.rm = TRUE),
#     clin_cases_under5   = rowSums(across(63:65), na.rm = TRUE)
#   )



# 
# ml_ctrysum_real_pop <- ml_data_sum_realpop %>%
#   group_by(ISO2, net_strategy, net_name, mass_int_yr, sample_id) %>%  # keep net_name here
#   summarise(across(
#     .cols = c(
#       starts_with("n_age_"),
#       starts_with("n_detect_lm_"),
#       starts_with("cases_"),
#       starts_with("clin_cases_"),
#       starts_with("sev_cases_"),
#       matches("pfpr_"),
#       clin_cases_all_ages,
#       clin_cases_under5
#     ),
#     .fns = sum,
#     na.rm = TRUE
#   ), .groups = "drop")

ml_ctrysum_real_pop <- ml_data_sum_realpop %>%
  group_by(ISO2, net_strategy, net_name, mass_int_yr, sample_id) %>%
  summarise(
    across(
      .cols = c(
        starts_with("n_age_"),
        starts_with("n_detect_lm_"),
        starts_with("cases_"),
        starts_with("clin_cases_"),
        starts_with("sev_cases_"),
        matches("pfpr_"),
        clin_cases_all_ages,
        clin_cases_under5
      ),
      .fns = sum,
      na.rm = TRUE
    ),
    EIR_urep_fit = first(EIR_urep_fit),
    mean_u       = first(mean_u),
    mean_a       = first(mean_a),
    .groups = "drop"
  )



ml_ctrysum_real_pop <- ml_ctrysum_real_pop %>%
  mutate(net_strategy = str_remove(net_strategy, " uncosted$"))

ml_ctrysum_real_pop <- ml_ctrysum_real_pop %>%
  mutate(
    facet_group = case_when(
      net_strategy == "no future nets"             ~ "No ITNs",
      is.na(mass_int_yr)                           ~ "Routine distribution only",
      mass_int_yr == 2                              ~ "2-year campaigns + routine",
      mass_int_yr == 3                              ~ "3-year campaigns + routine",
      TRUE                                          ~ NA_character_
    ),
    facet_group = factor(facet_group, levels = c(
      "No ITNs",
      "Routine distribution only",
      "2-year campaigns + routine",
      "3-year campaigns + routine"
    ))
  )

annot_df <- ml_ctrysum_real_pop %>%
  filter(net_name == "Pyrethroid-only") %>%
  group_by(facet_group) %>%
  summarise(
    mean_u_med  = median(mean_u, na.rm = TRUE),
    mean_u_lb   = quantile(mean_u, 0.025, na.rm = TRUE),
    mean_u_ub   = quantile(mean_u, 0.975, na.rm = TRUE),
    mean_a_med  = median(mean_a, na.rm = TRUE),
    mean_a_lb   = quantile(mean_a, 0.025, na.rm = TRUE),
    mean_a_ub   = quantile(mean_a, 0.975, na.rm = TRUE),
    uga_med     = median(mean_u / mean_a, na.rm = TRUE),
    uga_lb      = quantile(mean_u / mean_a, 0.025, na.rm = TRUE),
    uga_ub      = quantile(mean_u / mean_a, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf(
      "Use: %d%% [%d, %d]\nAccess: %d%% [%d, %d]\nUse | Access: %d%% [%d, %d]",
      round(mean_u_med * 100), round(mean_u_lb * 100), round(mean_u_ub * 100),
      round(mean_a_med * 100), round(mean_a_lb * 100), round(mean_a_ub * 100),
      round(uga_med * 100),    round(uga_lb * 100),    round(uga_ub * 100)
    )
  )


# Determine max y for positioning annotations
y_max <- max(ml_ctrysum_real_pop$clin_cases_all_ages, na.rm = TRUE)

# Final plot
ggplot(ml_ctrysum_real_pop, aes(x = net_strategy, y = clin_cases_all_ages, fill = net_name)) +
  geom_violin(color = NA, alpha = 0.3, trim = FALSE, scale = "width") +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 2,
    aes(color = net_name),
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data = function(x) {
      data.frame(
        y = median(x),
        ymin = quantile(x, 0.025),
        ymax = quantile(x, 0.975)
      )
    },
    geom = "errorbar",
    width = 0.2,
    aes(color = net_name),
    show.legend = FALSE
  ) +
  scale_fill_discrete(na.translate = TRUE, labels = function(x) ifelse(is.na(x), "none", x)) +
  scale_color_discrete(na.translate = TRUE, labels = function(x) ifelse(is.na(x), "none", x)) +
  labs(
    title = "All-Age Clinical Cases by Net Strategy",
    x = NULL,
    y = "Annual Clinical Cases (All Ages)",
    fill = "Net Type"
  ) +
  facet_grid(. ~ facet_group, scales = "free_x", space = "free_x", switch = "x") +
  geom_text(
    data = annot_df,
    aes(
      x = 0.5,  # shift slightly right
      y = y_max * 0.95,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,  # left-justified
    size = 3.5
  ) +
  theme_bw() +
  theme(
    panel.spacing = unit(0, "pt"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )













# Summarise mean_u by strategy and facet
summary_mean_u <- ml_ctrysum_real_pop %>%
  group_by(net_strategy, facet_group) %>%
  summarise(
    median_u = median(mean_u, na.rm = TRUE),
    lower_u = quantile(mean_u, 0.025, na.rm = TRUE),
    upper_u = quantile(mean_u, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


# Final plot
ggplot(ml_ctrysum_real_pop, aes(x = net_strategy, y = clin_cases_all_ages, fill = net_name)) +
  geom_violin(color = NA, alpha = 0.3, trim = FALSE, scale = "width") +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 2,
    aes(color = net_name),
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data = function(x) {
      data.frame(
        y = median(x),
        ymin = quantile(x, 0.025),
        ymax = quantile(x, 0.975)
      )
    },
    geom = "errorbar",
    width = 0.2,
    aes(color = net_name),
    show.legend = FALSE
  ) +
  scale_fill_discrete(na.translate = TRUE, labels = function(x) ifelse(is.na(x), "none", x)) +
  scale_color_discrete(na.translate = TRUE, labels = function(x) ifelse(is.na(x), "none", x)) +
  labs(
    title = "All-Age Clinical Cases by Net Strategy",
    x = "Net Strategy",
    y = "Annual Clinical Cases (All Ages)",
    fill = "Net Type"
  ) +
  # facet_grid(. ~ facet_group, scales = "free_x", space = "free_x") +
  facet_grid(. ~ facet_group, scales = "free_x", space = "free_x", switch = "x") +
  theme_bw() +
  theme(
    panel.spacing = unit(0, "pt"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )








current_country <- "MLI"

ctry_site <- readRDS(
  eval(paste0(sitefile_folder, current_country, "_site.rds"))
)
























# Define the base directory
base_dir <- "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary/Pyrethroid-only_routine_baseline/BF"

# List all fs_id_* directories
fs_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

# Function to read all sim*.rds files in numeric order
read_sim_files <- function(folder) {
  sim_files <- list.files(folder, pattern = "^sim\\d+\\.rds$", full.names = TRUE)
  
  # Sort files by numeric simulation number
  sim_files <- sim_files[order(as.integer(gsub("^.*sim(\\d+)\\.rds$", "\\1", sim_files)))]
  
  purrr::map_dfr(sim_files, readRDS, .id = "file_id")
}

# Read and combine all data
BFtest_combined_data <- purrr::map_dfr(fs_dirs, read_sim_files, .id = "fs_id")

# Add fs_id and sim number columns
# BFtest_combined_data <- dplyr::mutate(
#   BFtest_combined_data,
#   fs_id = gsub(".*/(fs_id_\\d+).*", "\\1", fs_id),
#   sim = as.integer(gsub("^sim(\\d+)\\.rds$", "\\1", basename(file_id)))
# )

# Write to CSV
#write.csv(BFonly3_combined_data, "M:/Andrew/Github/itn_campaigns/BF_annual_urep_only3.csv", row.names = FALSE)














# Simulation parameters
mass_int_yr <- c(2,3)     # Mass campaign intervals
sim_reps <- 100
sim_cores <- 16          # Takes precedence over areas per core
sim_areas_per_core <- 1   # Requires sim_cores = 0, otherwise ignored
sim_id <- timestamp_05













for (i in 6:6){#1:length(SSA_ISO2)) {
  
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  
  for (j in 1:length(fs_areas_included)) {
    
    area_selected <- fs_areas_included[j]
    
    # Only 2-year
    area_sim_id <- paste(
      "past_sim",
      sim_id,
      gsub(" ", "_", area_selected),
      sep = "_"
    )
    assign(
      area_sim_id,
      net_data %>%
        run_malsim_nets_sequential_v5(
          N_reps= sim_reps,
          areas_per_core = sim_areas_per_core,
          N_cores = sim_cores,
          areas_included = area_selected,
          mass_int_yr = 2,
          only = TRUE,
          use_hipercow = TRUE,
          bv_beta = bv_beta,
          bv_gamma = bv_gamma,
          past_sim = TRUE,
          results_path = results_06
        )
    )
    new_past_sim_df <- data.frame("id.iso" = SSA_ISO2[i],
                                  "id.name" = area_sim_id,
                                  "id.string" = get(area_sim_id))
    past_sim_df %<>% rbind.data.frame(new_past_sim_df)
    
  }
  
}



  
saveRDS(past_sim_df, file = file.path(results_06, paste0("past_sim_df.rds")))


#hipercow_id_df <- readRDS("BF_costed_03APR_hipercow_id_df.rds")

#BF_short_costed <- extract_hipercow_net_runs(hipercow_id_df$id.string)


























for (i in 1:length(SSA_ISO2)) {
  
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  fs_area_EIRs <- fs_id_link$EIR_fit[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  
  
  
  hc_sim_id <- paste(
    "sim",
    sim_id,
    SSA_ISO2[i],
    "test",
    sep = "_"
  )
  assign(
    hc_sim_id,#!!!!!
    net_data %>%
      run_malsim_nets_sequential_v5(
        N_reps= sim_reps,
        areas_per_core = sim_areas_per_core,
        areas_included = fs_areas_included,
        mass_int_yr = 3,
        only = TRUE,
        use_hipercow = TRUE,
        bv_beta = bv_beta,
        bv_gamma = bv_gamma,
        baseline_sim = TRUE,
        results_folder = results_06
      )
  )
  
}


for (k in 1:length(budget_val)) {
  
  for (i in 1:1) {
    
    fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
    
    
    for (j in 1:length(fs_areas_included)) {
      
      area_selected <- fs_areas_included[j]
      
      # Only 2-year
      area_sim_id <- paste(
        "sim",
        sim_id,
        gsub(" ", "_", area_selected),
        2,
        budget_val[k],
        "bud_only",
        sep = "_"
      )
      assign(
        area_sim_id,
        eir_net_data %>%
          run_malsim_nets_sequential_v4(
            N_reps= sim_reps,
            areas_per_core = sim_areas_per_core,
            N_cores = sim_cores,
            areas_included = area_selected,
            mass_int_yr = 2,
            only = TRUE,
            use_hipercow = TRUE,
            bv_beta = bv_beta,
            bv_gamma = bv_gamma,
            budget_pc = budget_val[k],
            net_costings = TRUE,
            biennial_reduction = TRUE
          )
      )
      new_hipercow_id_df <- data.frame("id.iso" = SSA_ISO2[i],
                                       "id.budget" = budget_val[k],
                                       "id.name" = area_sim_id,
                                       "id.string" = get(area_sim_id))
      hipercow_id_df %<>% rbind.data.frame(new_hipercow_id_df)
      
      # PBO 2-year
      area_sim_id <- paste(
        "sim",
        sim_id,
        gsub(" ", "_", area_selected),
        2,
        budget_val[k],
        "bud_pbo",
        sep = "_"
      )
      assign(
        area_sim_id,
        eir_net_data %>%
          run_malsim_nets_sequential_v4(
            N_reps= sim_reps,
            areas_per_core = sim_areas_per_core,
            N_cores = sim_cores,
            areas_included = area_selected,
            mass_int_yr = 2,
            pbo = TRUE,
            use_hipercow = TRUE,
            bv_beta = bv_beta,
            bv_gamma = bv_gamma,
            budget_pc = budget_val[k],
            net_costings = TRUE,
            biennial_reduction = TRUE
          )
      )
      new_hipercow_id_df <- data.frame("id.iso" = SSA_ISO2[i],
                                       "id.budget" = budget_val[k],
                                       "id.name" = area_sim_id,
                                       "id.string" = get(area_sim_id))
      hipercow_id_df %<>% rbind.data.frame(new_hipercow_id_df)
      
      # Pyrrole 2-year
      area_sim_id <- paste(
        "sim",
        sim_id,
        gsub(" ", "_", area_selected),
        2,
        budget_val[k],
        "bud_pyrrole",
        sep = "_"
      )
      assign(
        area_sim_id,
        eir_net_data %>%
          run_malsim_nets_sequential_v4(
            N_reps= sim_reps,
            areas_per_core = sim_areas_per_core,
            N_cores = sim_cores,
            areas_included = area_selected,
            mass_int_yr = 2,
            pyrrole = TRUE,
            use_hipercow = TRUE,
            bv_beta = bv_beta,
            bv_gamma = bv_gamma,
            budget_pc = budget_val[k],
            net_costings = TRUE,
            biennial_reduction = TRUE
          )
      )
      new_hipercow_id_df <- data.frame("id.iso" = SSA_ISO2[i],
                                       "id.budget" = budget_val[k],
                                       "id.name" = area_sim_id,
                                       "id.string" = get(area_sim_id))
      hipercow_id_df %<>% rbind.data.frame(new_hipercow_id_df)
      
      # only 3-year
      area_sim_id <- paste(
        "sim",
        sim_id,
        gsub(" ", "_", area_selected),
        3,
        budget_val[k],
        "bud_only",
        sep = "_"
      )
      assign(
        area_sim_id,
        eir_net_data %>%
          run_malsim_nets_sequential_v4(
            N_reps= sim_reps,
            areas_per_core = sim_areas_per_core,
            N_cores = sim_cores,
            areas_included = area_selected,
            mass_int_yr = 3,
            only = TRUE,
            use_hipercow = TRUE,
            bv_beta = bv_beta,
            bv_gamma = bv_gamma,
            budget_pc = budget_val[k],
            net_costings = TRUE
          )
      )
      new_hipercow_id_df <- data.frame("id.iso" = SSA_ISO2[i],
                                       "id.budget" = budget_val[k],
                                       "id.name" = area_sim_id,
                                       "id.string" = get(area_sim_id))
      hipercow_id_df %<>% rbind.data.frame(new_hipercow_id_df)
      
      # PBO 3-year
      area_sim_id <- paste(
        "sim",
        sim_id,
        gsub(" ", "_", area_selected),
        3,
        budget_val[k],
        "bud_pbo",
        sep = "_"
      )
      assign(
        area_sim_id,
        eir_net_data %>%
          run_malsim_nets_sequential_v4(
            N_reps= sim_reps,
            areas_per_core = sim_areas_per_core,
            N_cores = sim_cores,
            areas_included = area_selected,
            mass_int_yr = 3,
            pbo = TRUE,
            use_hipercow = TRUE,
            bv_beta = bv_beta,
            bv_gamma = bv_gamma,
            budget_pc = budget_val[k],
            net_costings = TRUE
          )
      )
      new_hipercow_id_df <- data.frame("id.iso" = SSA_ISO2[i],
                                       "id.budget" = budget_val[k],
                                       "id.name" = area_sim_id,
                                       "id.string" = get(area_sim_id))
      hipercow_id_df %<>% rbind.data.frame(new_hipercow_id_df)
      
      # Pyrrole 3-year
      area_sim_id <- paste(
        "sim",
        sim_id,
        gsub(" ", "_", area_selected),
        3,
        budget_val[k],
        "bud_pyrrole",
        sep = "_"
      )
      assign(
        area_sim_id,
        eir_net_data %>%
          run_malsim_nets_sequential_v4(
            N_reps= sim_reps,
            areas_per_core = sim_areas_per_core,
            N_cores = sim_cores,
            areas_included = area_selected,
            mass_int_yr = 3,
            pyrrole = TRUE,
            use_hipercow = TRUE,
            bv_beta = bv_beta,
            bv_gamma = bv_gamma,
            budget_pc = budget_val[k],
            net_costings = TRUE
          )
      )
      new_hipercow_id_df <- data.frame("id.iso" = SSA_ISO2[i],
                                       "id.budget" = budget_val[k],
                                       "id.name" = area_sim_id,
                                       "id.string" = get(area_sim_id))
      hipercow_id_df %<>% rbind.data.frame(new_hipercow_id_df)
      
    }
    
  }
  
}

saveRDS(hipercow_id_df, "BF_costed_03APR_hipercow_id_df.rds")

hipercow_id_df <- readRDS("BF_costed_03APR_hipercow_id_df.rds")

BF_short_costed <- extract_hipercow_net_runs(hipercow_id_df$id.string)













# Simulation parameters
mass_int_yr <- c(2,3)     # Mass campaign intervals
sim_reps <- 100
sim_cores <- 32          # Takes precedence over areas per core
sim_areas_per_core <- 1   # Requires sim_cores = 0, otherwise ignored
sim_id <- "MAR06b"

# Routine ITN channels only
# for (i in 1:N_ISO2) {
for (i in 5:5) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  # fs_excluded <- c("BF Hauts-Bassins rural",
  #                  "BF Hauts-Bassins urban")
  # fs_areas_included <- fs_areas_included[! fs_areas_included %in% fs_excluded]
  #fs_areas_included <- "BF Boucle du Mouhoun rural"
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_only", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             only = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pbo", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             pbo = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pyrrole", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             pyrrole = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
}


for (i in 1:1) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  # fs_excluded <- c("BF Hauts-Bassins rural",
  #                  "BF Hauts-Bassins urban")
  fs_areas_included <- fs_areas_included[! fs_areas_included %in% fs_excluded]
  
  # No future
  assign(paste("sim", sim_id, SSA_ISO2[i], "0", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             only = TRUE,
             routine_baseline = TRUE,
             no_future_nets = TRUE,
             use_hipercow = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  # Mass campaigns
  for (j in 1:length(mass_int_yr)) {
    assign(paste("sim", sim_id, SSA_ISO2[i], "only", mass_int_yr[j], sep = "_"), eir_net_data %>%
             run_malsim_nets_sequential_v4(
               N_reps= sim_reps,
               areas_per_core = sim_areas_per_core,
               N_cores = sim_cores,
               areas_included = fs_areas_included,
               mass_int_yr = mass_int_yr[j],
               only = TRUE,
               use_hipercow = TRUE,
               bv_beta = bv_beta,
               bv_gamma = bv_gamma
             )
    )
    assign(paste("sim", sim_id, SSA_ISO2[i], "pbo", mass_int_yr[j], sep = "_"), eir_net_data %>%
             run_malsim_nets_sequential_v4(
               N_reps= sim_reps,
               areas_per_core = sim_areas_per_core,
               N_cores = sim_cores,
               areas_included = fs_areas_included,
               mass_int_yr = mass_int_yr[j],
               pbo = TRUE,
               use_hipercow = TRUE,
               bv_beta = bv_beta,
               bv_gamma = bv_gamma
             )
    )
    assign(
      paste("sim", sim_id, SSA_ISO2[i], "pyrrole", mass_int_yr[j], sep = "_"),
      eir_net_data %>%
        run_malsim_nets_sequential_v4(
          N_reps= sim_reps,
          areas_per_core = sim_areas_per_core,
          N_cores = sim_cores,
          areas_included = fs_areas_included,
          mass_int_yr = mass_int_yr[j],
          pyrrole = TRUE,
          use_hipercow = TRUE,
          bv_beta = bv_beta,
          bv_gamma = bv_gamma
        )
    )
  }
}



# Simulation parameters
mass_int_yr <- c(2,3)     # Mass campaign intervals
sim_reps <- 100
sim_cores <- 0          # Takes precedence over areas per core
sim_areas_per_core <- 1   # Requires sim_cores = 0, otherwise ignored
sim_id <- "MAR06c"

# Routine ITN channels only
# for (i in 1:N_ISO2) {
for (i in 6:6) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  # fs_excluded <- c("BF Hauts-Bassins rural",
  #                  "BF Hauts-Bassins urban")
  # fs_areas_included <- fs_areas_included[! fs_areas_included %in% fs_excluded]
  #fs_areas_included <- "BF Boucle du Mouhoun rural"
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_only", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             only = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pbo", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             pbo = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pyrrole", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             pyrrole = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
}


for (i in 6:6) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  # fs_excluded <- c("BF Hauts-Bassins rural",
  #                  "BF Hauts-Bassins urban")
  fs_areas_included <- fs_areas_included[! fs_areas_included %in% fs_excluded]
  
  # No future
  assign(paste("sim", sim_id, SSA_ISO2[i], "0", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             only = TRUE,
             routine_baseline = TRUE,
             no_future_nets = TRUE,
             use_hipercow = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  # Mass campaigns
  for (j in 1:length(mass_int_yr)) {
    assign(paste("sim", sim_id, SSA_ISO2[i], "only", mass_int_yr[j], sep = "_"), eir_net_data %>%
             run_malsim_nets_sequential_v4(
               N_reps= sim_reps,
               areas_per_core = sim_areas_per_core,
               N_cores = sim_cores,
               areas_included = fs_areas_included,
               mass_int_yr = mass_int_yr[j],
               only = TRUE,
               use_hipercow = TRUE,
               bv_beta = bv_beta,
               bv_gamma = bv_gamma
             )
    )
    assign(paste("sim", sim_id, SSA_ISO2[i], "pbo", mass_int_yr[j], sep = "_"), eir_net_data %>%
             run_malsim_nets_sequential_v4(
               N_reps= sim_reps,
               areas_per_core = sim_areas_per_core,
               N_cores = sim_cores,
               areas_included = fs_areas_included,
               mass_int_yr = mass_int_yr[j],
               pbo = TRUE,
               use_hipercow = TRUE,
               bv_beta = bv_beta,
               bv_gamma = bv_gamma
             )
    )
    assign(
      paste("sim", sim_id, SSA_ISO2[i], "pyrrole", mass_int_yr[j], sep = "_"),
      eir_net_data %>%
        run_malsim_nets_sequential_v4(
          N_reps= sim_reps,
          areas_per_core = sim_areas_per_core,
          N_cores = sim_cores,
          areas_included = fs_areas_included,
          mass_int_yr = mass_int_yr[j],
          pyrrole = TRUE,
          use_hipercow = TRUE,
          bv_beta = bv_beta,
          bv_gamma = bv_gamma
        )
    )
  }
}












# Simulation parameters
mass_int_yr <- c(2,3)     # Mass campaign intervals
sim_reps <- 10
sim_cores <- 0           # Takes precedence over areas per core
sim_areas_per_core <- 1   # Requires sim_cores = 0, otherwise ignored
sim_id <- "MAR05_10runs"

# Routine ITN channels only
# for (i in 1:N_ISO2) {
for (i in 1:1) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  # fs_excluded <- c("BF Hauts-Bassins rural",
  #                  "BF Hauts-Bassins urban")
  # fs_areas_included <- fs_areas_included[! fs_areas_included %in% fs_excluded]
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_only", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             only = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pbo", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             pbo = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pyrrole", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             pyrrole = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
}


for (i in 1:1) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  # fs_excluded <- c("BF Hauts-Bassins rural",
  #                  "BF Hauts-Bassins urban")
  fs_areas_included <- fs_areas_included[! fs_areas_included %in% fs_excluded]
  
  # No future
  assign(paste("sim", sim_id, SSA_ISO2[i], "0", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             only = TRUE,
             routine_baseline = TRUE,
             no_future_nets = TRUE,
             use_hipercow = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  # Mass campaigns
  for (j in 1:length(mass_int_yr)) {
    assign(paste("sim", sim_id, SSA_ISO2[i], "only", mass_int_yr[j], sep = "_"), eir_net_data %>%
             run_malsim_nets_sequential_v4(
               N_reps= sim_reps,
               areas_per_core = sim_areas_per_core,
               N_cores = sim_cores,
               areas_included = fs_areas_included,
               mass_int_yr = mass_int_yr[j],
               only = TRUE,
               use_hipercow = TRUE,
               bv_beta = bv_beta,
               bv_gamma = bv_gamma
             )
    )
    assign(paste("sim", sim_id, SSA_ISO2[i], "pbo", mass_int_yr[j], sep = "_"), eir_net_data %>%
             run_malsim_nets_sequential_v4(
               N_reps= sim_reps,
               areas_per_core = sim_areas_per_core,
               N_cores = sim_cores,
               areas_included = fs_areas_included,
               mass_int_yr = mass_int_yr[j],
               pbo = TRUE,
               use_hipercow = TRUE,
               bv_beta = bv_beta,
               bv_gamma = bv_gamma
             )
    )
    assign(paste("sim", sim_id, SSA_ISO2[i], "pyrrole", mass_int_yr[j], sep = "_"), eir_net_data %>%
             run_malsim_nets_sequential_v4(
               N_reps= sim_reps,
               areas_per_core = sim_areas_per_core,
               N_cores = sim_cores,
               areas_included = fs_areas_included,
               mass_int_yr = mass_int_yr[j],
               pyrrole = TRUE,
               use_hipercow = TRUE,
               bv_beta = bv_beta,
               bv_gamma = bv_gamma
             )
    )
  }
}






# Simulation parameters
mass_int_yr <- c(2,3)     # Mass campaign intervals
sim_reps <- 100
sim_cores <- 32          # Takes precedence over areas per core
sim_areas_per_core <- 1   # Requires sim_cores = 0, otherwise ignored
sim_id <- "MAR06c_rural"

# Routine ITN channels only
# for (i in 1:N_ISO2) {
for (i in 2:2) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[
    which(fs_id_link$ISO2 == SSA_ISO2[i] & fs_id_link$urbanicity == "rural")
  ]
  # fs_excluded <- c("BF Hauts-Bassins rural",
  #                  "BF Hauts-Bassins urban")
  # fs_areas_included <- fs_areas_included[! fs_areas_included %in% fs_excluded]
  #fs_areas_included <- "BF Boucle du Mouhoun rural"
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_only", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             only = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pbo", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             pbo = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pyrrole", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             pyrrole = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
}


for (i in 2:2) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[
    which(fs_id_link$ISO2 == SSA_ISO2[i] & fs_id_link$urbanicity == "rural")
  ]
  
  # No future
  assign(paste("sim", sim_id, SSA_ISO2[i], "0", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             only = TRUE,
             routine_baseline = TRUE,
             no_future_nets = TRUE,
             use_hipercow = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  # Mass campaigns
  for (j in 1:length(mass_int_yr)) {
    assign(paste("sim", sim_id, SSA_ISO2[i], "only", mass_int_yr[j], sep = "_"), eir_net_data %>%
             run_malsim_nets_sequential_v4(
               N_reps= sim_reps,
               areas_per_core = sim_areas_per_core,
               N_cores = sim_cores,
               areas_included = fs_areas_included,
               mass_int_yr = mass_int_yr[j],
               only = TRUE,
               use_hipercow = TRUE,
               bv_beta = bv_beta,
               bv_gamma = bv_gamma
             )
    )
    assign(paste("sim", sim_id, SSA_ISO2[i], "pbo", mass_int_yr[j], sep = "_"), eir_net_data %>%
             run_malsim_nets_sequential_v4(
               N_reps= sim_reps,
               areas_per_core = sim_areas_per_core,
               N_cores = sim_cores,
               areas_included = fs_areas_included,
               mass_int_yr = mass_int_yr[j],
               pbo = TRUE,
               use_hipercow = TRUE,
               bv_beta = bv_beta,
               bv_gamma = bv_gamma
             )
    )
    assign(
      paste("sim", sim_id, SSA_ISO2[i], "pyrrole", mass_int_yr[j], sep = "_"),
      eir_net_data %>%
        run_malsim_nets_sequential_v4(
          N_reps= sim_reps,
          areas_per_core = sim_areas_per_core,
          N_cores = sim_cores,
          areas_included = fs_areas_included,
          mass_int_yr = mass_int_yr[j],
          pyrrole = TRUE,
          use_hipercow = TRUE,
          bv_beta = bv_beta,
          bv_gamma = bv_gamma
        )
    )
  }
}





# Simulation parameters
mass_int_yr <- c(2,3)     # Mass campaign intervals
sim_reps <- 100
sim_cores <- 32          # Takes precedence over areas per core
sim_areas_per_core <- 1   # Requires sim_cores = 0, otherwise ignored
sim_id <- "MAR06c_urban"

# Routine ITN channels only
# for (i in 1:N_ISO2) {
for (i in 2:2) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[
    which(fs_id_link$ISO2 == SSA_ISO2[i] & fs_id_link$urbanicity == "urban")
  ]
  # fs_excluded <- c("BF Hauts-Bassins rural",
  #                  "BF Hauts-Bassins urban")
  # fs_areas_included <- fs_areas_included[! fs_areas_included %in% fs_excluded]
  #fs_areas_included <- "BF Boucle du Mouhoun rural"
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_only", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             only = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pbo", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             pbo = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pyrrole", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             pyrrole = TRUE,
             routine_baseline = TRUE,
             use_hipercow = TRUE,
             #hiper_debug = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
}


for (i in 2:2) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[
    which(fs_id_link$ISO2 == SSA_ISO2[i] & fs_id_link$urbanicity == "urban")
  ]
  
  # No future
  assign(paste("sim", sim_id, SSA_ISO2[i], "0", sep = "_"), eir_net_data %>%
           run_malsim_nets_sequential_v4(
             N_reps= sim_reps,
             areas_per_core = sim_areas_per_core,
             N_cores = sim_cores,
             areas_included = fs_areas_included,
             mass_int_yr = 3,
             only = TRUE,
             routine_baseline = TRUE,
             no_future_nets = TRUE,
             use_hipercow = TRUE,
             bv_beta = bv_beta,
             bv_gamma = bv_gamma
           )
  )
  
  # Mass campaigns
  for (j in 1:length(mass_int_yr)) {
    assign(paste("sim", sim_id, SSA_ISO2[i], "only", mass_int_yr[j], sep = "_"), eir_net_data %>%
             run_malsim_nets_sequential_v4(
               N_reps= sim_reps,
               areas_per_core = sim_areas_per_core,
               N_cores = sim_cores,
               areas_included = fs_areas_included,
               mass_int_yr = mass_int_yr[j],
               only = TRUE,
               use_hipercow = TRUE,
               bv_beta = bv_beta,
               bv_gamma = bv_gamma
             )
    )
    assign(paste("sim", sim_id, SSA_ISO2[i], "pbo", mass_int_yr[j], sep = "_"), eir_net_data %>%
             run_malsim_nets_sequential_v4(
               N_reps= sim_reps,
               areas_per_core = sim_areas_per_core,
               N_cores = sim_cores,
               areas_included = fs_areas_included,
               mass_int_yr = mass_int_yr[j],
               pbo = TRUE,
               use_hipercow = TRUE,
               bv_beta = bv_beta,
               bv_gamma = bv_gamma
             )
    )
    assign(
      paste("sim", sim_id, SSA_ISO2[i], "pyrrole", mass_int_yr[j], sep = "_"),
      eir_net_data %>%
        run_malsim_nets_sequential_v4(
          N_reps= sim_reps,
          areas_per_core = sim_areas_per_core,
          N_cores = sim_cores,
          areas_included = fs_areas_included,
          mass_int_yr = mass_int_yr[j],
          pyrrole = TRUE,
          use_hipercow = TRUE,
          bv_beta = bv_beta,
          bv_gamma = bv_gamma
        )
    )
  }
}




extract_hipercow_net_runs <- function(ids) {
  sim_data <- NULL
  for (i in 1:length(ids)) {
    scenario_data <- ids[i] %>% task_result %>% do.call(rbind.data.frame, .)
    sim_data %<>% rbind.data.frame(scenario_data)
  }
  return(sim_data)
}


BF_sim_MAR06 <- extract_hipercow_net_runs(
  c(
    "92c29643bf0ce367f6f56cfb9360e9c9",
    "4ad03b2be73d3d5945e04ca534994134",
    "f5e4a8d397c9f6d47be9302101a246da",
    "63ed2918f836ddb54eaf5876a8de3769",
    "f4b4f5336effca8b9587dbd0829b52d6",
    "1c12dd42a25ce852f56dc15366f4f1ba",
    "260063ce61e8aa85b072c25be7c8c07d",
    "9c40ac6dc5560f5ca67430c741197896",
    "7e399f783d77ba4b596a1b947cac4ce5",
    "54fe76f5bb4f72d0d120ec2cae06f8a9"
  )
)
BF_sim_MAR06$adj_ann_camp_nets_dist[BF_sim_MAR06$ann_camp_nets_dist == 0] <- 0
write.csv(BF_sim_MAR06, "BF_sim_02APR.csv")



ML_sim_MAR06 <- extract_hipercow_net_runs(
  c(
    "ff28596b4cd8286d0ecd91ed24abca27",
    "f51ec8d06ec60c7456ed805bd7f6cebe",
    "16466979f4301623cd6d32ad7f6fc662",
    "1866519d8c503d8358743828f656e891",
    "57b703b5b344770431ab036cd1f82ef2",
    "93563db696992a33c040e515db0cd386",
    "d10571627cfd4347ac03f8703b40ca85",
    "c037008ea3ed8309c8df1532a80c3e5b",
    "2ed59b4080a09fd4b94ba878fd737620",
    "c24a9303bacbabe97dad8572117d37ce"
  )
)
ML_sim_MAR06$adj_ann_camp_nets_dist[ML_sim_MAR06$ann_camp_nets_dist == 0] <- 0
write.csv(ML_sim_MAR06, "ML_sim_02APR.csv")



SN_sim_MAR06 <- extract_hipercow_net_runs(
  c(
    "d797d5780ee6e39fa3fc9b24d2093a7c",
    "a1daa1d054859fca0d96cd644ce478d3",
    "aed9a5ca9fd6075e6dccf6992f63f95c",
    "649e127e2943bb3716f37bd92f313ad6",
    "cc8937c79a81082868473086575ecca2",
    "c93138d5777598b8fdeb33a393b09d0f",
    "f1bd025691b86576c5496f4f80b50414",
    "6dcb1de996d2418776ba474cc7b3cca8",
    "35d115896c8723af7f394b93c8cc5a4b",
    "7e950dab27253686746638792a875838"
  )
)
SN_sim_MAR06$adj_ann_camp_nets_dist[SN_sim_MAR06$ann_camp_nets_dist == 0] <- 0
write.csv(SN_sim_MAR06, "SN_sim_02APR.csv")





# sim_MAR06_BF_0,
# sim_MAR06_BF_routine_only,
# sim_MAR06_BF_routine_pbo,
# sim_MAR06_BF_routine_pyrrole,
# sim_MAR06_BF_only_2,
# sim_MAR06_BF_pbo_2,
# sim_MAR06_BF_pyrrole_2,
# sim_MAR06_BF_only_3,
# sim_MAR06_BF_pbo_3,
# sim_MAR06_BF_pyrrole_3#,


sim_MAR05_10runs_uncosted_data <- extract_hipercow_net_runs(
  c(
    sim_MAR05_10runs_BF_0,
    sim_MAR05_10runs_BF_routine_only,
    sim_MAR05_10runs_BF_routine_pbo,
    sim_MAR05_10runs_BF_routine_pyrrole,
    sim_MAR05_10runs_BF_only_2,
    sim_MAR05_10runs_BF_pbo_2,
    sim_MAR05_10runs_BF_pyrrole_2,
    sim_MAR05_10runs_BF_only_3,
    sim_MAR05_10runs_BF_pbo_3,
    sim_MAR05_10runs_BF_pyrrole_3#,
    # sim_31JAN24_BF_only_costed_2,
    # sim_31JAN24_BF_pbo_costed_2,
    # sim_31JAN24_BF_pyrrole_costed_2,
    # sim_31JAN24_BF_pbo_costed_3,
    # sim_31JAN24_BF_pyrrole_costed_3
  )
)

write.csv(sim_MAR05_10runs_uncosted_data, "sim_MAR05_10runs_uncosted_data.csv")


save(sim_12FEB24b_BF_0,
     sim_12FEB24b_BF_routine_only,
     sim_12FEB24b_BF_routine_pbo,
     sim_12FEB24b_BF_routine_pyrrole,
     sim_12FEB24b_BF_only_2,
     sim_12FEB24b_BF_pbo_2,
     sim_12FEB24b_BF_pyrrole_2,
     sim_12FEB24b_BF_only_3,
     sim_12FEB24b_BF_pbo_3,
     sim_12FEB24b_BF_pyrrole_3,
     file="12FEB24b_hipercow_ids.RData")

load("./data/12FEB24b_hipercow_ids.RData")


write.csv(sim_12FEB24b_uncosted_data, "sim_12FEB24b_uncosted_data.csv")



# 
# for (i in 1:1) {
#   
#   # Sub-set areas by country
#   fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
#   # fs_excluded <- c("BF Hauts-Bassins rural",
#   #                  "BF Hauts-Bassins urban")
#   fs_areas_included <- fs_areas_included[! fs_areas_included %in% fs_excluded]
#   
#   # Mass campaigns
#   for (j in 1:length(mass_int_yr)) {
#     assign(paste("sim", sim_id, SSA_ISO2[i], "only_costed", mass_int_yr[j], sep = "_"), net_data %>%
#              run_malsim_nets_sequential_v4(
#                N_reps= 100,
#                #areas_per_core = 1,
#                N_cores = 32,
#                areas_included = fs_areas_included,
#                mass_int_yr = mass_int_yr[j],
#                only = TRUE,
#                use_hipercow = TRUE,
#                bv_beta = bv_beta,
#                bv_gamma = bv_gamma,
#                biennial_reduction = TRUE
#                #net_costings = TRUE
#              )
#     )
#     assign(paste("sim", sim_id, SSA_ISO2[i], "pbo_costed", mass_int_yr[j], sep = "_"), net_data %>%
#              run_malsim_nets_sequential_v4(
#                N_reps = 100,
#                #areas_per_core = 1,
#                N_cores = 32,
#                areas_included = fs_areas_included,
#                mass_int_yr = mass_int_yr[j],
#                pbo = TRUE,
#                use_hipercow = TRUE,
#                bv_beta = bv_beta,
#                bv_gamma = bv_gamma,
#                biennial_reduction = TRUE,
#                net_costings = TRUE
#              )
#     )
#     assign(paste("sim", sim_id, SSA_ISO2[i], "pyrrole_costed", mass_int_yr[j], sep = "_"), net_data %>%
#              run_malsim_nets_sequential_v4(
#                N_reps = 100,
#                #areas_per_core = 1,
#                N_cores = 32,
#                areas_included = fs_areas_included,
#                mass_int_yr = mass_int_yr[j],
#                pyrrole = TRUE,
#                use_hipercow = TRUE,
#                bv_beta = bv_beta,
#                bv_gamma = bv_gamma,
#                biennial_reduction = TRUE,
#                net_costings = TRUE
#              )
#     )
#   }
# }
# 
