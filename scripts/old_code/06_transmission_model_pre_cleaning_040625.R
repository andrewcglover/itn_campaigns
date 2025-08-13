# 06_transmission_model.R

#-------------------------------------------------------------------------------
# Malaria Simulation

# Load net resistance data

#if (only) {res_only <- read.csv("./data/pyrethroid_only_nets.csv")}
#if (pbo) {res_pbo <- read.csv("./data/pyrethroid_pbo_nets.csv")}
#if (pyrrole) {res_pyrrole <- read.csv("./data/pyrethroid_pyrrole_nets.csv")}

# Create timestamped data cache for section 05
timestamp_05 <- format(Sys.time(), "%Y%m%d%H%M%S")
cache_05 <- paste0("./data_cache/05_", timestamp_05)
dir.create(cache_05, recursive = TRUE)

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
# fs_excluded <- c("GH Bono rural",
#                  "GH Bono urban",
#                  "GH Bono East rural",
#                  "GH Bono East urban",
#                  "GH Ahafo rural",
#                  "GH Ahafo urban",
#                  "GH Savannah rural",
#                  "GH Savannah urban",
#                  "GH North East rural",
#                  "GH North East urban",
#                  "GH Oti rural",
#                  "GH Oti urban"
#                  #"BF Hauts-Bassins rural",
#                  #"BF Hauts-Bassins urban"
# )
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
    "./scripts/post_use_access_fitting/cali_funs.R",
    "./scripts/post_use_access_fitting/EIR_cali.R",
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
sim_reps <- 20
sim_cores <- 2          # Takes precedence over areas per core
sim_areas_per_core <- 1   # Requires sim_cores = 0, otherwise ignored
sim_id <- "03APRcost"
budget_val <- c(100, 75, 50)


hipercow_id_df <- data.frame("id.iso" = NULL,
                             "id.budget" = NULL,
                             "id.name" = NULL,
                             "id.string" = NULL)

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
