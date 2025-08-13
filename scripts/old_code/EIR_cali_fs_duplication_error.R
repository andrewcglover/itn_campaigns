# EIR_cali.R

err_par_EIR_fit <- function(param_list) {
  
  #library(magrittr)
  
  #N_reps <- length(param_list)
  
  #browser()
  
  # Extract parameters from parameter list
  site_pars <- param_list#[[1]] #Uncomment for testing
  mean_retu <- site_pars$mean_retu
  mean_reta <- site_pars$mean_reta
  npc_beta <- site_pars$npc_beta
  npc_gamma <- site_pars$npc_gamma
  net_type <- site_pars$net_type
  net_name <- site_pars$net_name
  net_strategy <- site_pars$net_strategy
  last_camp <- site_pars$last_camp
  mass_int_mn <- site_pars$mass_int_mn
  fs_area <- site_pars$fs_area
  ISO2 <- site_pars$ISO2
  fs_name_1 <- site_pars$fs_name_1
  urbanicity <- site_pars$urbanicity
  fs_area_id <- site_pars$fs_area_id
  N_species <- site_pars$N_species
  CMC_first <- site_pars$CMC_first
  CMC_Jan2000 <- site_pars$CMC_Jan2000
  N_CMC <- site_pars$N_CMC
  N_CMC_sim <- site_pars$N_CMC_sim
  sim_population <- site_pars$sim_population
  P_a_mid <- site_pars$P_a_mid
  P0_a_mid <- site_pars$P0_a_mid
  D_a_mid <- site_pars$D_a_mid
  lam_a_mid <- site_pars$lam_a_mid
  P_u_mid <- site_pars$P_u_mid
  P0_u_mid <- site_pars$P0_u_mid
  D_u_mid <- site_pars$D_u_mid
  lam_u_mid <- site_pars$lam_u_mid
  dn0_mat <- site_pars$dn0_mat
  rn_mat <- site_pars$rn_mat
  rnm_mat <- site_pars$rnm_mat
  gam_vec <- site_pars$gam_vec
  DOY_1st <- site_pars$DOY_1st
  DOY_mid <- site_pars$DOY_mid
  prevalence_rendering_min_ages <- site_pars$prevalence_rendering_min_ages
  prevalence_rendering_max_ages <- site_pars$prevalence_rendering_max_ages
  pfpr_full <- site_pars$pfpr_full
  pfpr_alpha <- site_pars$pfpr_alpha
  pfpr_beta <- site_pars$pfpr_beta
  pfpr_var <- site_pars$pfpr_var
  pfpr_lo <- site_pars$pfpr_lo
  pfpr_hi <- site_pars$pfpr_hi
  target <- site_pars$target
  use_existing_eir <- site_pars$use_existing_eir
  existing_eir <- site_pars$existing_eir
  use_hipercow <- site_pars$use_hipercow
  cali_attempts <- site_pars$cali_attempts

  fs_area_undrscr <- gsub(" ", "_", fs_area)
  
  year = 365

  # Monthly rates
  decay_a <- 1 / mean_reta
  decay_u <- 1 / mean_retu
  decay_uga <- decay_u - decay_a
  decay_npc <- decay_a / npc_gamma
  
  # convert retention to days
  mean_retu_dy <- 365 * mean_retu / 12
  mean_reta_dy <- 365 * mean_reta / 12
  
  # Central time point for first regular mass campaign

  # Back-fill preliminary period
  N_early <- CMC_first - CMC_Jan2000
  
  # Times for fitting
  times_mn <- seq(1, N_CMC)
  times_yr <- rep(seq(0, ceiling(N_CMC_sim / 12)), each=12)
  times_1st_dy <- DOY_1st + (times_yr * year)
  times_mid_dy <- DOY_mid + (times_yr * year)
  input_net_times <- times_mid_dy[1:N_CMC_sim]    # usage for fitting
  output_net_times <- times_1st_dy[1:N_CMC_sim]   # distribution times for netz
  
  #-----------------------------------------------------------------------------
  # usage profiles
  
  # Extract values for selected sample
  P_ui <- P_u_mid
  P0_ui <- P0_u_mid
  D_ui <- D_u_mid
  lambda_ui <- lam_u_mid
  
  # Backfill routine nets
  P_u_full <- c(rep(D_ui[1], N_early), P_ui)
  
  # Fit nets for calibration
  cali_dist_list <- fit_usage_sequential_distributions(
    target_usage = P_u_full,
    target_usage_timesteps = input_net_times,
    distribution_timesteps = output_net_times,
    mean_retention = mean_retu_dy
  )
  
  # Extract use "coverages"
  cali_dist <- cali_dist_list[[1]]
  
  dn0_mat <-  site_pars$bednet_dn0[
    rep(seq_len(nrow(site_pars$bednet_dn0)), each = 12),
  ]
  rn_mat <-  site_pars$bednet_rn[
    rep(seq_len(nrow(site_pars$bednet_rn)), each = 12),
  ]
  rnm_mat <-  site_pars$bednet_rnm[
    rep(seq_len(nrow(site_pars$bednet_rnm)), each = 12),
  ]
  gamman_vec <- rep(site_pars$bednet_gamman, each = 12)
  
  # set bednets
  bednet_pars <- malariasimulation::set_bednets(site_pars,
                                                timesteps = output_net_times,
                                                coverages = cali_dist,
                                                retention = mean_retu_dy,
                                                dn0 = dn0_mat,
                                                rn = rn_mat,
                                                rnm = rnm_mat,
                                                gamman = gamman_vec)
  
  # Prevalence target
  inc_months <- which(!is.nan(pfpr_full))
  target_prev <- pfpr_full[inc_months]
  target_lo <- pfpr_lo[inc_months]
  target_hi <- pfpr_hi[inc_months]
  target_var <- pfpr_var[inc_months]
  weights <- 1 / target_var
  weights <- weights / sum(weights)
  weighted_target <- target_prev * weights
  inc_months_with_backfill <- inc_months + N_early
  
  #adm_site$prevalence$year
  #target <- adm_site$prevalence$pfpr
  
  # Cali summary function
  # weighted_filtered_annual_pfpr_summary <- function(x, w = weights){
  #   month <- ceiling(x$timestep / (365/12))
  #   pfpr_6_59_mo <- (x$n_detect_lm_182_730 + x$n_detect_lm_730_1825) /
  #     (x$n_age_182_730 + x$n_age_730_1825)
  #   filtered_months <- which(month %in% inc_months_with_backfill)
  #   pfpr_6_59_mo <- pfpr_6_59_mo[filtered_months]
  #   month <- month[filtered_months]
  #   tapply(pfpr_6_59_mo, month, mean) * w
  # }
  
  filtered_annual_pfpr_summary <- function(x){
    year <- ceiling(x$timestep / 365)
    pfpr_2_10 <- (x$n_detect_lm_730_1824 + x$n_detect_lm_1825_3649) /
      (x$n_age_730_1824 + x$n_age_1825_3649)
    tapply(pfpr_2_10, year, mean)
  }
  
  # debug_file <- paste0("debug2_output_", Sys.getpid(), ".txt")
  # writeLines(capture.output({
  #   str(use_hipercow)
  # }), debug_file)
  
  # Calibrate EIR
  if (use_existing_eir) {
    EIR_out <- existing_eir
  } else {
    set.seed(123)
    if (use_hipercow) {
      # EIR_out <- cali_calibrate(
      #   parameters = bednet_pars,
      #   target = target,
      #   summary_function = filtered_annual_pfpr_summary,
      #   eq_prevalence = target[1],
      #   max_attempts = 30,
      #   use_pfpr_6_59_mo = FALSE
      # )
      # EIR_out <- cali::calibrate(
      #   parameters = bednet_pars,
      #   target = target,
      #   summary_function = filtered_annual_pfpr_summary,
      #   eq_prevalence = target[1],
      #   max_attempts = 30
      # )
      EIR_out <- tryCatch({
        cali::calibrate(
          parameters = bednet_pars,
          target = target,
          summary_function = filtered_annual_pfpr_summary,
          eq_prevalence = target[1],
          max_attempts = cali_attempts
        )
      }, error = function(e) {
        log_file <- paste0("EIR_cali_error_", Sys.getpid(), ".txt")
        writeLines(
          paste0("Error in cali::calibrate(): ", conditionMessage(e)),
          con = log_file
        )
        return(NA)  # Or return a sentinel value
      })
    } else {
      EIR_out <- cali::calibrate(
        parameters = bednet_pars,
        target = target,
        summary_function = filtered_annual_pfpr_summary,
        eq_prevalence = target[1],
        max_attempts = cali_attempts,
        use_pfpr_6_59_mo = FALSE
      )
    }
  }
  
  debug_file <- paste0("EIR1_output_", Sys.getpid(), ".txt")
  writeLines(capture.output({
    str(EIR_out)
  }), debug_file)


  
  # Malariasimulation test run
  bednet_pars$human_population <- sim_population
  bednet_pars <- malariasimulation::set_equilibrium(bednet_pars, init_EIR = EIR_out)
  raw <- malariasimulation::run_simulation(bednet_pars$timesteps, parameters = bednet_pars)
  #raw$pfpr <- raw$n_detect_lm_182_1825  / raw$n_age_182_1825
  raw$pfpr_2_10 <- (raw$n_detect_lm_730_1824 + raw$n_detect_lm_1825_3649) /
    (raw$n_age_730_1824 + raw$n_age_1825_3649)
  raw$pfpr_6_59_mo <- (raw$n_detect_lm_182_729 + raw$n_detect_lm_730_1824) /
    (raw$n_age_182_729 + raw$n_age_730_1824)
  
  
  raw_months <- floor(raw$timestep / (365/12))
  fit_CMC <- unique(raw_months) + CMC_Jan2000
  fit_pfpr_2_10 <- tapply(raw$pfpr_2_10, raw_months, mean)
  fit_pfpr_6_59_mo <- tapply(raw$pfpr_6_59_mo, raw_months, mean)
  fit_CMC <- fit_CMC[(N_early+1):(length(fit_CMC)-1)]
  fit_pfpr_2_10 <- fit_pfpr_2_10[(N_early+1):(length(fit_pfpr_2_10)-1)]
  fit_pfpr_6_59_mo <- fit_pfpr_6_59_mo[(N_early+1):(length(fit_pfpr_6_59_mo)-1)]
  
  #pfpr_map <- rep(NA,length(target)*12)
  # for (i in 1:length(target)) {
  #   pfpr_map[(i-1)*12+6] <- target[i]
  # }
  pfpr_map <- rep(target, each = 12)
  pfpr_map <- pfpr_map[(N_early+1):(length(pfpr_map))]
  
  fs_area <- site_pars$fs_area
  ISO2 <- site_pars$ISO2
  fs_name_1 <- site_pars$fs_name_1
  urbanicity <- site_pars$urbanicity
  fs_area_id <- site_pars$fs_area_id
  
  debug_file <- paste0("EIR2_output_", Sys.getpid(), ".txt")
  writeLines(capture.output({
    str(EIR_out)
  }), debug_file)
  
  
  eir_output <- data.frame(
    "fs_area" = rep(fs_area, length(fit_CMC)),
    "ISO2" = rep(ISO2, length(fit_CMC)),
    "fs_name_1" = rep(fs_name_1, length(fit_CMC)),
    "urbanicity" = rep(urbanicity, length(fit_CMC)),
    "fs_area_id" = rep(fs_area_id, length(fit_CMC)),
    "CMC" = fit_CMC,
    "pfpr_2_10_fit" = fit_pfpr_2_10,
    "pfpr_6_59_mo_fit" = fit_pfpr_6_59_mo,
    "pfpr_map" = pfpr_map,
    "pfpr_dhs" = pfpr_full,
    "pfpr_dhs_lo" = pfpr_lo,
    "pfpr_dhs_hi" = pfpr_hi,
    "EIR_fit" = rep(EIR_out, length(fit_CMC))
    )
# 
#   ggplot(
#     data = fit_df,
#     aes(
#       x = CMC,
#       y = pfpr_2_10_fit,
#     )
#   ) +
#     geom_line(
#       col = "purple"
#     ) +
#     geom_step(
#       aes(y = pfpr_map),
#       position = position_nudge(x = -0.5),
#       col = "dodgerblue"
#     ) +
#     ylim(0, 1) +
#     ylab(expression(P*italic(f)*Pr["6-59mo"])) +
#     xlab("Time") +
#     theme_bw()
# 
#   ggplot(
#     data = fit_df,
#     aes(
#       x = CMC,
#       y = pfpr_dhs,
#       ymin = pfpr_dhs_lo,
#       ymax = pfpr_dhs_hi
#     )
#   ) +
#     geom_line(
#       aes(y = pfpr_6_59_mo_fit),
#       col = "purple"
#     ) +
#     geom_point(
#       col = "dodgerblue"
#     ) +
#     geom_errorbar(
#       col = "dodgerblue"
#     ) +
#     ylim(0, 1) +
#     ylab(expression(P*italic(f)*Pr["6-59mo"])) +
#     xlab("Time") +
#     theme_bw()

  # 
  # 
  # raw_sub <- raw[((N_early*365)+1):(dim(raw)[1]-1),]
  # 
  # # 
  # ggplot() +
  #   geom_line(
  #     data = raw_sub,
  #     aes(
  #       x = timestep,
  #       y = pfpr,
  #     ),
  #     col = "deeppink",
  #     linewidth = 0.5
  #   ) +
  #   geom_point(
  #     aes(
  #       x = inc_months * 365 / 12,
  #       y = target
  #     ),
  #     col = "dodgerblue",
  #     size = 2
  #   ) +
  #   geom_errorbar(
  #     aes(
  #       x =  inc_months * 365 / 12,
  #       y = target,
  #       ymin = target_lo,
  #       ymax = target_hi
  #     ),
  #     col = "dodgerblue"
  #   ) +
  #   ylim(0, 1) +
  #   ylab(expression(P*italic(f)*Pr["6-59mo"])) +
  #   xlab("Time") +
  #   theme_bw()
  
  
  return(eir_output)
  
}

err_run_eir_calibration <- function(
    dataset,
    areas_included = NULL,
    N_cores = 0,
    areas_per_core = 1,
    sim_population = 1e5,
    ref_CMC = 1500,
    use_hipercow = FALSE,
    bv_beta = NULL,
    bv_gamma = NULL,
    local_sitefiles = TRUE,
    sitefile_folder = "./data_private/newsitefiles/",
    use_existing_eir = FALSE,
    cap_dataset_prior_to_expansion = TRUE,
    run_matches = TRUE,
    run_mismatches = TRUE,
    hiper_debug = FALSE,
    cali_attempts = 10
) {
  
  # Restrict dataset prior to exansion to align with use-access Stan fitting
  #if (cap_dataset_prior_to_expansion) {
  #  dataset <- dataset[1:(N_areas*N_CMC),]
  #}
  
  # Set number of cores
  #total_runs <- N_reps*length(areas_included)
  if (N_cores <= 0) {N_cores <- ceiling(length(areas_included)/areas_per_core)}
  if (use_hipercow) {max_cores <- 32} else {max_cores <- 18}
  if (max_cores > length(areas_included)) {max_cores <- length(areas_included)}
  if (N_cores > max_cores) {N_cores <- max_cores}
  
  # Simulation time
  CMC_sim_start <- CMC_Jan2000
  CMC_sim_end <- CMC_last
  N_CMC_old_nets <- CMC_last - CMC_sim_start + 1
  new_net_start_mn <- N_CMC_old_nets + 1
  N_CMC_sim <- CMC_sim_end - CMC_sim_start + 1
  
  # dataframe for storing output
  output_df <- data.frame(NULL)
  
  # Empty parameter list
  param_list <- list()
  
  # current country
  current_country <- "XXX"
  
  pc0 <- 0
  pc1 <- 0
  for (i in 1:N_fs_areas) {
    
    if (fs_id_link$fs_area[i] %in% areas_included) {
      
      print(paste("area index", i, "selected:"))
      print(paste(fs_id_link$fs_area[i]))
      
      # Warning for foresite mismatch
      if (fs_id_link$fs_area_id[i] != i) {
        print("Warning: Foresite id mismatch")
      }
      
      # Area-time indices
      area_id <- fs_id_link$new_area_id[i]
      area_time_ids <- which(usage_list$a == area_id)
      area_time_ref_id <- area_time_ids[ref_CMC - CMC_first + 1]
      #area_time_ref_id <- which(dataset$area_id == area_id &
      #                            dataset$CMC == ref_CMC)
      #area_time_ids <- which(dataset$area_id == area_id)
      
      # get samples
      invlam_u_mid <- invlam_u[, area_id] %>%
        as.vector %>% unname %>% unlist %>% mean
      lam_u_mid <- 1 / invlam_u_mid
      ret_u_mid <- ret_u[, area_time_ref_id] %>%
        as.vector %>% unname %>% unlist %>% mean
      P_u_mid <- P_u[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      P0_u_mid <- P0_u[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      D_u_mid <- D_u[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      
      invlam_a_mid <- invlam_a[, area_id] %>%
        as.vector %>% unname %>% unlist %>% mean
      lam_a_mid <- 1 / invlam_a_mid
      ret_a_mid <- ret_a[, area_time_ref_id] %>%
        as.vector %>% unname %>% unlist %>% mean
      P_a_mid <- P_a[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      P0_a_mid <- P0_a[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      D_a_mid <- D_a[, area_time_ids] %>%
        as.matrix %>% unname %>% colMeans
      
      # Nets per capita samples
      npc_beta_mid <- bv_beta %>% mean
      npc_gamma_mid <- bv_gamma %>% mean
      
      # Identify last mass campaign
      # Estimated from latest peak in use
      # Previous approach:
      # last_camp_month <- apply(P_u_samples, 1, which.max)
      u0 <- 1
      u1 <- 0
      t <- length(P_u_mid)
      while ((u1 < u0 )&(t > 1)){
        u0 <- P_u_mid[t-1]
        u1 <- P_u_mid[t]
        if (u0 < u1) {last_camp_month <- t}
        t <- t - 1
      }
      
      # Generate ISO code for current admin
      admin_country <- countrycode(fs_id_link$ISO2[i], "iso2c", "iso3c")
      
      # Pull country site
      if (admin_country != current_country) {
        current_country <- admin_country
        if (local_sitefiles) {
          ctry_site <- readRDS(
            eval(paste0(sitefile_folder, current_country, "_site.rds"))
            )
        } else {
          ctry_site <- site::fetch_site(current_country)
        }
      }
      
      # Isolate single admin site
      adm_site_matches <- ctry_site$sites[
        which(ctry_site$sites$name_1 == fs_id_link$fs_name_1[i]),
      ]
      if (fs_id_link$urbanicity[i] %in% adm_site_matches$urban_rural) {
        urb_select <- fs_id_link$urbanicity[i]
      } else if (fs_id_link$urbanicity[i] == "urban") {
        urb_select <- "rural"
      } else {
        urb_select <- "urban"
      }
      adm_site <- site::subset_site(
        site = ctry_site,
        site_filter = data.frame(
          country = countrycode(fs_id_link$ISO2[i], "iso2c", "country.name"),
          iso3c = admin_country,
          name_1 = fs_id_link$fs_name_1[i],
          urban_rural = urb_select
        )
      )
      
      N_species <- dim(adm_site$vectors$vector_species)[1]
      
      # Prevalence
      area_id_selected <- area_id
      adm_dataset <- dataset %>% dplyr::filter(area_id == area_id_selected)
      pfpr_full <- adm_dataset$rdt_prev
      pfpr_alpha <- adm_dataset$rdt_pos + 0.5
      pfpr_beta <- adm_dataset$rdt_neg + 0.5
      pfpr_var <- (pfpr_alpha * pfpr_beta) /
        ((pfpr_alpha + pfpr_beta)^2 * (pfpr_alpha + pfpr_beta + 1))
      pfpr_lo <- adm_dataset$lo_prev
      pfpr_hi <- adm_dataset$hi_prev
      # pfpr_full <- dataset$rdt_prev[area_time_ids]
      # pfpr_alpha <- dataset$rdt_pos[area_time_ids] + 0.5
      # pfpr_beta <- dataset$rdt_neg[area_time_ids] + 0.5
      # pfpr_var <- (pfpr_alpha * pfpr_beta) /
      #   ((pfpr_alpha + pfpr_beta)^2 * (pfpr_alpha + pfpr_beta + 1))
      # pfpr_lo <- dataset$lo_prev[area_time_ids]
      # pfpr_hi <- dataset$hi_prev[area_time_ids]
      
      target <- adm_site$prevalence$pfpr

        # Create parameter inputs

        site_pars <- site_parameters(
          interventions = adm_site$interventions,
          demography = adm_site$demography,
          vectors = adm_site$vectors$vector_species,
          seasonality = adm_site$seasonality$seasonality_parameters,
          eir = adm_site$eir$eir[1],
          overrides = list(
            human_population = sim_population,
            individual_mosquitoes = FALSE
          )
        )
        
       
        # Set age groups
        full_age_groups_min <- round(c(0.5, 2, 5) * 365)
        full_age_groups_max <- round(c(2, 5, 10) * 365) - 1
        site_pars$prevalence_rendering_min_ages <- full_age_groups_min
        site_pars$prevalence_rendering_max_ages <- full_age_groups_max
        
        # Combine vector and matrix parameters for parLapply function
        site_pars$P_u_mid <- P_u_mid
        site_pars$P0_u_mid <- P0_u_mid
        site_pars$D_u_mid <- D_u_mid
        site_pars$lam_u_mid <- lam_u_mid
        site_pars$P_a_mid <- P_a_mid
        site_pars$P0_a_mid <- P0_a_mid
        site_pars$D_a_mid <- D_a_mid
        site_pars$lam_a_mid <- lam_a_mid
        site_pars$DOY_1st <- DOY_1st
        site_pars$DOY_mid <- DOY_mid
        
        # Combine single parameters
        site_pars$ISO2 <- fs_id_link$ISO2[i]
        site_pars$fs_area <- fs_id_link$fs_area[i]
        site_pars$fs_name_1 <- fs_id_link$fs_name_1[i]
        site_pars$urbanicity <- fs_id_link$urbanicity[i]
        site_pars$fs_area_id <- fs_id_link$fs_area_id[i]
        site_pars$N_species <- N_species
        site_pars$CMC_first <- CMC_first
        site_pars$CMC_Jan2000 <- CMC_Jan2000
        site_pars$N_CMC <- N_CMC
        site_pars$N_CMC_sim <- N_CMC_sim
        site_pars$sim_population <- sim_population

        # Net parameters
        site_pars$mean_retu <- ret_u_mid
        site_pars$mean_reta <- ret_a_mid
        site_pars$npc_beta <- npc_beta_mid
        site_pars$npc_gamma <- npc_gamma_mid
        
        # Prev parameters
        site_pars$pfpr_full <- pfpr_full
        site_pars$pfpr_alpha <- pfpr_alpha
        site_pars$pfpr_beta <- pfpr_beta
        site_pars$pfpr_var <- pfpr_var
        site_pars$pfpr_lo <- pfpr_lo
        site_pars$pfpr_hi <- pfpr_hi
        site_pars$target <- target
        
        # Existing EIR
        site_pars$use_existing_eir <- use_existing_eir
        if(use_existing_eir) {
          site_pars$existing_eir <- fs_id_link$EIR_fit[i]
        } else {
          site_pars$existing_eir <- NULL
        }
        
        # Hipercow indicator
        site_pars$use_hipercow <- use_hipercow
        
        # Max cali attempts
        site_pars$cali_attempts <- cali_attempts
        
        # Add site parameters to list
        if (run_matches &
            (fs_id_link$urbanicity[i] %in% adm_site_matches$urban_rural)) {
          param_list[[length(param_list) + 1]] <- site_pars
        }
        if (run_mismatches &
            !(fs_id_link$urbanicity[i] %in% adm_site_matches$urban_rural)) {
          param_list[[length(param_list) + 1]] <- site_pars
        }
      
      # 
      pc1 <- round(100 * i / N_fs_areas)
      if (pc1 > pc0) {
       pc0 <- pc1
       print(paste(pc0, "% complete", sep = ""))
      }
      
    }
  }
  
  
  # writeLines(capture.output(str(param_list[[1]]$use_hipercow)), "my_variable.txt")
  
  
  
  if (use_hipercow) {
    resources <- hipercow_resources(cores = N_cores)
    if (hiper_debug == TRUE) {
      # par_id <- task_create_expr(par_EIR_fit(param_list[[1]]),
      #                            resources = resources)
      par_id <- task_create_expr(par_EIR_fit(param_list[[1]]),
                                 resources = hipercow_resources(cores = 1))
    } else {
      par_id <- task_create_expr(
        parallel::parLapply(NULL, param_list, par_EIR_fit),
        parallel = hipercow_parallel("parallel"),
        resources = resources)
    }
    
    
    # par_out <- task_result(par_output)
    # comb_output <- do.call(rbind.data.frame, par_out)
    # 
    # comb_output <- do.call(rbind.data.frame, par_output)
    # output_df <- rbind(output_df, comb_output)
  } else {
    #if (!use_hipercow) {
    cl <- makeCluster(N_cores)
    clusterExport(cl, c(#"param_list",
      "set_bednets",
      "run_simulation",
      "set_equilibrium",
      "fit_usage_sequential_distributions",
      "calibrate",
      "population_usage_t"))
    #par_output <- lapply(param_list, par_net_region_sequential3)
    par_output <- parLapply(cl, param_list, par_EIR_fit)
    comb_output <- do.call(rbind.data.frame, par_output)
    output_df <- rbind(output_df, comb_output)
    output_df <- output_df[!duplicated(cbind(output_df$CMC,
                                             output_df$fs_area_id)), ]
    stopCluster(cl)
  }
  
  if (use_hipercow) {
    return(par_id)
  } else {
    return(output_df)
  }
}