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

options(hipercow.max_size_local = 1e10)

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
sim_cores <- 25          # Takes precedence over areas per core
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


sim_cores <- 2          # Takes precedence over areas per core

for (i in 5:6){#1:length(SSA_ISO2)) {
  
  for (k in 1:2) {
    
    fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
    
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

sim_cores <- 2

for (i in 5:6){#1:length(SSA_ISO2)) {
  
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  
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
    
    fs_dirsx <- list.dirs(strat_pathx, full.names = TRUE, recursive = FALSE)
    fs_dirsx <- fs_dirsx[grepl("fs_id_\\d+$", basename(fs_dirsx))]
    
    purrr::map(fs_dirsx, function(fs_dirx) {
      list.files(fs_dirx, pattern = "^sim\\d+\\.rds$", full.names = TRUE)
    }) %>% unlist()
  }) %>% unlist()
  
  # Safely read and bind all .rds files
  combinedx <- purrr::map_dfr(all_rds_filesx, function(filex) {
    tryCatch(readRDS(filex), error = function(e) NULL)
  })
  
  return(combinedx)
}


merge_populationx <- function(sim_datax, iso3x, start_yearx, sitefile_folderx) {
  
  # Use fixed site file folder
  sitefile_folderx <- "./data_private/newsitefiles/"
  
  # Load site metadata
  site_pathx <- file.path(sitefile_folderx, paste0(iso3x, "_site.rds"))
  ctry_sitex <- readRDS(site_pathx)
  
  # # Load site metadata
  # site_pathx <- file.path(sitefile_folderx, paste0(iso3x, "_site.rds"))
  # ctry_sitex <- readRDS(site_pathx)
  
  # Extract population table
  popx <- ctry_sitex$population$population_total
  
  # Derive calendar year in simulation data
  sim_datax <- sim_datax %>%
    mutate(calendar_year = start_yearx + year_id - 1)
  
  # Join on fs_name_1 <-> name_1, urbanicity <-> urban_rural, calendar_year <-> year
  sim_with_popx <- sim_datax %>%
    left_join(popx,
              by = c(
                "fs_name_1" = "name_1",
                "urbanicity" = "urban_rural",
                "calendar_year" = "year"
              ))
  
  return(sim_with_popx)
}


# Example usage:
base_path <- "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary"
ML_data <- read_country_simsx("ML", base_path)

ml_data_pop <- merge_populationx(
  sim_datax = ML_data,
  iso3x = "MLI",
  start_yearx = 2025
) %>%
  left_join(fs_id_link %>% select(fs_area_id, new_area_id, EIR_urep_fit), by = "fs_area_id")


ml_data_pop <- ml_data_pop %>%
  mutate(
    # --- Calculate usage values for untreated nets ---
    C0_u = pmap_dbl(
      list(new_area_id, sample_id),
      ~ {
        area <- as.integer(..1)
        samp_id <- as.integer(..2)
        sample_index <- long_sample_ids[samp_id]
        C0_u[sample_index, usage_list$N_t] %>%
          as.matrix() %>%
          unname()
      }
    ),
    
    D_u = pmap_dbl(
      list(new_area_id, sample_id),
      ~ {
        area <- as.integer(..1)
        samp_id <- as.integer(..2)
        sample_index <- long_sample_ids[samp_id]
        D_u[sample_index, usage_list$N_t] %>%
          as.matrix() %>%
          unname()
      }
    ),
    
    invlam_u = pmap_dbl(
      list(new_area_id, sample_id),
      ~ {
        area <- as.integer(..1)
        samp_id <- as.integer(..2)
        sample_index <- long_sample_ids[samp_id]
        invlam_u[sample_index, area] %>%
          as.vector() %>%
          unname() %>%
          unlist()
      }
    ),
    
    lam_u = 1 / invlam_u
  ) %>%
  # --- Apply scaling factor to untreated usage ---
  mutate(
    scale_factor = 1 - (0.08 + rnormvals[sample_id] * 0.05 / 1.96),
    C0_u = C0_u * scale_factor,
    D_u  = D_u  * scale_factor
  ) %>%
  # --- Compute mean untreated usage ---
  mutate(
    mean_u = case_when(
      !is.na(mass_int_yr) ~ D_u + (C0_u / (lam_u * (mass_int_yr * 12))) * (1 - exp(-lam_u * (mass_int_yr * 12))),
      is.na(mass_int_yr) & net_strategy == "no future nets" ~ 0,
      is.na(mass_int_yr) & net_strategy != "no future nets" ~ D_u,
      TRUE ~ NA_real_
    )
  ) %>%
  # --- Now compute values for treated nets (_a) ---
  mutate(
    C0_a = pmap_dbl(
      list(new_area_id, sample_id),
      ~ {
        area <- as.integer(..1)
        samp_id <- as.integer(..2)
        sample_index <- long_sample_ids[samp_id]
        C0_a[sample_index, usage_list$N_t] %>%
          as.matrix() %>%
          unname()
      }
    ),
    
    D_a = pmap_dbl(
      list(new_area_id, sample_id),
      ~ {
        area <- as.integer(..1)
        samp_id <- as.integer(..2)
        sample_index <- long_sample_ids[samp_id]
        D_a[sample_index, usage_list$N_t] %>%
          as.matrix() %>%
          unname()
      }
    ),
    
    invlam_a = pmap_dbl(
      list(new_area_id, sample_id),
      ~ {
        area <- as.integer(..1)
        samp_id <- as.integer(..2)
        sample_index <- long_sample_ids[samp_id]
        invlam_a[sample_index, area] %>%
          as.vector() %>%
          unname() %>%
          unlist()
      }
    ),
    
    lam_a = 1 / invlam_a,
    
    mean_a = case_when(
      !is.na(mass_int_yr) ~ D_a + (C0_a / (lam_a * (mass_int_yr * 12))) * (1 - exp(-lam_a * (mass_int_yr * 12))),
      is.na(mass_int_yr) & net_strategy == "no future nets" ~ 0,
      is.na(mass_int_yr) & net_strategy != "no future nets" ~ D_a,
      TRUE ~ NA_real_
    )
  )




# 
# 
# ml_data_sum_pop <- ml_data_pop %>%
#   group_by(ISO2, fs_name_1, urbanicity, fs_area, fs_area_id,
#            net_name, mass_int_yr, biennial_reduction, net_costings,
#            no_future_nets, routine_baseline, sample_id, net_strategy, budget_pc) %>%
#   summarise(across(
#     matches("(_tot$|_mean$)"),  # selects all renamed *_pop and pfpr_*_mean columns
#     mean, na.rm = TRUE
#   ), .groups = "drop") %>%
#   mutate(
#     clin_cases_all_ages = rowSums(across(63:78), na.rm = TRUE),
#     clin_cases_under5   = rowSums(across(63:65), na.rm = TRUE)
#   )

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
