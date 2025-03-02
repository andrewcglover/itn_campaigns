# 06_transmission_model.R

#-------------------------------------------------------------------------------
# Malaria Simulation

# Load net resistance data

#if (only) {res_only <- read.csv("./data/pyrethroid_only_nets.csv")}
#if (pbo) {res_pbo <- read.csv("./data/pyrethroid_pbo_nets.csv")}
#if (pyrrole) {res_pyrrole <- read.csv("./data/pyrethroid_pyrrole_nets.csv")}

# Convert projection times to months
mass_int_mn <- mass_int_yr * 12
projection_window_mn <- projection_window_yr * 12
projection_window_dy <- projection_window_yr * 365

fs_areas_included <- unique(fs_id_link$fs_area)
# fs_excluded <- c(fs_id_link$fs_area[fs_id_link$ISO2=="GH"],
#                  fs_id_link$fs_area[fs_id_link$ISO2=="MW"],
#                  "BF Hauts-Bassins rural",
#                  "BF Hauts-Bassins urban",
#                  "BF Plateau Central rural",
#                  "BF Plateau Central urban")
fs_excluded <- c("GH Bono rural",
                 "GH Bono urban",
                 "GH Bono East rural",
                 "GH Bono East urban",
                 "GH Ahafo rural",
                 "GH Ahafo urban",
                 "GH Savannah rural",
                 "GH Savannah urban",
                 "GH North East rural",
                 "GH North East urban",
                 "GH Oti rural",
                 "GH Oti urban"
                 #"BF Hauts-Bassins rural",
                 #"BF Hauts-Bassins urban"
)
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
hipercow_environment_create(sources = c("./scripts/utils/simulation_v3.R",
                                        "./scripts/utils/netz_usage_sequential_branch_adapted.R"))
hipercow_provision()

options(hipercow.max_size_local = 1e10)

only_cost <- 1.95
pbo_cost <- 2.54
pyrrole_cost <- 2.56
dist_cost <- 2.75

only_total_cost <- dist_cost + only_cost
pbo_total_cost <- dist_cost + pbo_cost
pyrrole_total_cost <- dist_cost + pyrrole_cost

scaled_pbo_nets_equiv_only <- only_total_cost / pbo_total_cost
scaled_pyrrole_nets_equiv_only <- only_total_cost / pyrrole_total_cost

# Simulation parameters
mass_int_yr <- c(2,3)     # Mass campaign intervals
sim_reps <- 100
sim_cores <- 32           # Takes precedence over areas per core
sim_areas_per_core <- 1   # Requires sim_cores = 0, otherwise ignored

# Routine ITN channels only
for (i in 1:N_ISO2) {
  
  # Sub-set areas by country
  fs_areas_included <- fs_id_link$fs_area[which(fs_id_link$ISO2 == SSA_ISO2[i])]
  # fs_excluded <- c("BF Hauts-Bassins rural",
  #                  "BF Hauts-Bassins urban")
  # fs_areas_included <- fs_areas_included[! fs_areas_included %in% fs_excluded]
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_only", sep = "_"), net_data %>%
           run_malsim_nets_sequential_v3(
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
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pbo", sep = "_"), net_data %>%
           run_malsim_nets_sequential_v3(
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
  
  assign(paste("sim", sim_id, SSA_ISO2[i], "routine_pyrrole", sep = "_"), net_data %>%
           run_malsim_nets_sequential_v3(
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