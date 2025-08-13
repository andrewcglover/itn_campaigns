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
# fs_areas_included <- c("SN Dakar urban",
#                        "SN Sédhiou rural",
#                        "SN Kolda rural")

# Number of samples
N_samples <- dim(P_u)[1]

# Monthly offset for future mass campaigns
# month_offset <- sample.int(13, N_reps, replace = TRUE) - 7

# Create sample ids
#long_sample_ids <- sample.int(N_samples, 10000 , replace = TRUE)
#saveRDS(long_sample_ids, "./data/800_sample_ids.rds")
long_sample_ids <- readRDS("./data/800_sample_ids.rds")
rnormvals <- readRDS("./data/rnormvals.rds")


hipercow::hipercow_configuration()
hipercow::hipercow_init(driver = "dide-windows")
hipercow::windows_authenticate()
# hipercow_environment_create(sources = c("./scripts/utils/simulation_new.R",
#                                         "./scripts/utils/simulation_costed.R",
#                                         #"./scripts/utils/simulation.R",
#                                         #"./scripts/utils/simulation2.R",
#                                         "./scripts/utils/netz_usage_sequential_branch_funs.R"))
#hipercow_provision(method="pkgdepends",refs=c("mrc-ide/malariasimulation mrc-ide/netz@usage_sequential"))
# hipercow_environment_create(sources = c("./scripts/utils/simulation_npc.R",
#                                         "./scripts/utils/netz_usage_sequential_branch_adapted.R"))
hipercow_environment_create(sources = c("./scripts/utils/simulation_v3.R",
                                        "./scripts/utils/netz_usage_sequential_branch_adapted.R"))
hipercow_provision()
#a<-as.numeric(Sys.time())*100000

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

