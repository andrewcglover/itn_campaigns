shorten_sim_data <- function(iso2) {
  
  sim_data_name <- paste0(iso2, "_ann_data")
  short_sim_data_name <- paste0(iso2, "_short_ann_data")
  
  if (!exists(sim_data_name, envir = .GlobalEnv)) {
    stop("Data object ", sim_data_name,
         " does not exist in the global environment.")
  }
  
  sim_data <- get(sim_data_name, envir = .GlobalEnv)
  
  short_sim_data <- sim_data %>%
    dplyr::select(
      ISO2, fs_name_1, urbanicity, fs_area, fs_area_id, new_area_id, pop,
      pop_weight, sample_id, net_strategy, net_name, mass_int_yr,
      no_future_nets, n_age_under5, adj_ann_total_nets_dist,
      cases_all_ages, cases_under5, clin_cases_all_ages, clin_cases_under5,
      sev_cases_all_ages, sev_cases_under5
    )
  
  message("Shortening of sim data complete. Assigning to ", short_sim_data_name)
  assign(short_sim_data_name, short_sim_data, envir = .GlobalEnv)
  
}

calculate_cases_per_100k_under5 <- function(
    iso2, obj_name_ending = "_short_ann_data"
) {
  
  obj_name <- paste0(iso2, obj_name_ending)
  
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

create_full_sim_comparison <- function(iso2) {
  
  short_sim_data_name <- paste0(iso2, "_short_ann_data")
  data_comparison_name <- paste0(iso2, "_comp_data")
  
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
        adj_ann_total_nets_dist_comp = adj_ann_total_nets_dist,
        cases_all_ages_comp = cases_all_ages,
        cases_under5_comp = cases_under5,
        cases_under5_p100kU5_comp = cases_under5_p100kU5,
        clin_cases_all_ages_comp = clin_cases_all_ages,
        clin_cases_under5_comp = clin_cases_under5,
        clin_cases_under5_p100kU5_comp = clin_cases_under5_p100kU5,
        sev_cases_all_ages_comp = sev_cases_all_ages,
        sev_cases_under5_comp = sev_cases_under5,
        sev_cases_under5_p100kU5_comp = sev_cases_under5_p100kU5
      )
    
    # Loop over intervention strategies
    for (j in 1:N_strategies) {
      
      # Filter for intervention
      int_sim_data <- short_sim_data %>%
        dplyr::filter(net_strategy == unique_strategies[j]) %>%
        dplyr::select(
          net_strategy, net_name, mass_int_yr, no_future_nets, n_age_under5,
          adj_ann_total_nets_dist, cases_all_ages, cases_under5,
          cases_under5_p100kU5, clin_cases_all_ages, clin_cases_under5,
          clin_cases_under5_p100kU5, sev_cases_all_ages, sev_cases_under5,
          sev_cases_under5_p100kU5
        ) %>%
        dplyr::rename(
          net_strategy_int = net_strategy,
          net_name_int = net_name,
          mass_int_yr_int = mass_int_yr,
          no_future_nets_int = no_future_nets, 
          n_age_under5_int = n_age_under5,
          adj_ann_total_nets_dist_int = adj_ann_total_nets_dist,
          cases_all_ages_int = cases_all_ages,
          cases_under5_int = cases_under5,
          cases_under5_p100kU5_int = cases_under5_p100kU5,
          clin_cases_all_ages_int = clin_cases_all_ages,
          clin_cases_under5_int = clin_cases_under5,
          clin_cases_under5_p100kU5_int = clin_cases_under5_p100kU5,
          sev_cases_all_ages_int = sev_cases_all_ages,
          sev_cases_under5_int = sev_cases_under5,
          sev_cases_under5_p100kU5_int = sev_cases_under5_p100kU5
        )
      
      # Combine selected comparator and intervention dataframes
      comb_sim_data <- cbind.data.frame(comp_sim_data, int_sim_data)
      
      # Append to dataframe with for all comparisons
      all_comb_sim_data <- rbind.data.frame(all_comb_sim_data, comb_sim_data)
      
    }
  }
  
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
        sev_cases_under5_p100kU5_comp - sev_cases_under5_p100kU5_int
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

create_ctry_sim_comparison <- function(iso2) {
  
  data_comparison_name <- paste0(iso2, "_comp_data")
  ctry_comparison_name <- paste0(iso2, "_ctry_comp_data")
  
  if (!exists(data_comparison_name, envir = .GlobalEnv)) {
    stop("Data object ", data_comparison_name,
         " does not exist in the global environment.")
  }
  
  all_comb_sim_data <- get(data_comparison_name, envir = .GlobalEnv)
  
  ctry_comb_data <- all_comb_sim_data %>%
    dplyr::mutate(
      
      n_age_under5_comp_wt = n_age_under5_comp * pop_weight,
      
      adj_ann_total_nets_dist_comp_wt = adj_ann_total_nets_dist_comp * pop_weight, 
      
      cases_all_ages_comp_wt = cases_all_ages_comp * pop_weight, 
      cases_under5_comp_wt = cases_under5_comp * pop_weight,
      cases_under5_p100kU5_comp_wt = cases_under5_p100kU5_comp * pop_weight,
      
      clin_cases_all_ages_comp_wt = clin_cases_all_ages_comp * pop_weight,
      clin_cases_under5_comp_wt = clin_cases_under5_comp * pop_weight,
      clin_cases_under5_p100kU5_comp_wt = clin_cases_under5_p100kU5_comp * pop_weight,
      
      sev_cases_all_ages_comp_wt = sev_cases_all_ages_comp * pop_weight, 
      sev_cases_under5_comp_wt = sev_cases_under5_comp * pop_weight,
      sev_cases_under5_p100kU5_comp_wt = sev_cases_under5_p100kU5_comp * pop_weight, 
      
      adj_ann_total_nets_dist_int_wt = adj_ann_total_nets_dist_int * pop_weight,
      
      cases_all_ages_int_wt = cases_all_ages_int * pop_weight, 
      cases_under5_int_wt = cases_under5_int * pop_weight,
      cases_under5_p100kU5_int_wt = cases_under5_p100kU5_int * pop_weight, 
      
      clin_cases_all_ages_int_wt = clin_cases_all_ages_int * pop_weight,        
      clin_cases_under5_int_wt = clin_cases_under5_int * pop_weight,
      clin_cases_under5_p100kU5_int_wt = clin_cases_under5_p100kU5_int * pop_weight,
      
      sev_cases_all_ages_int_wt = sev_cases_all_ages_int * pop_weight, 
      sev_cases_under5_int_wt = sev_cases_under5_int * pop_weight,
      sev_cases_under5_p100kU5_int_wt = sev_cases_under5_p100kU5_int * pop_weight, 
      
      adj_ann_total_nets_dist_add_wt = adj_ann_total_nets_dist_add * pop_weight, 
      
      cases_all_ages_avert_wt = cases_all_ages_avert * pop_weight, 
      cases_under5_avert_wt = cases_under5_avert * pop_weight, 
      cases_under5_p100kU5_avert_wt = cases_under5_p100kU5_avert * pop_weight,
      
      clin_cases_all_ages_avert_wt = clin_cases_all_ages_avert * pop_weight, 
      clin_cases_under5_avert_wt = clin_cases_under5_avert * pop_weight,
      clin_cases_under5_p100kU5_avert_wt = clin_cases_under5_p100kU5_avert * pop_weight,
      
      sev_cases_all_ages_avert_wt = sev_cases_all_ages_avert * pop_weight, 
      sev_cases_under5_avert_wt = sev_cases_under5_avert * pop_weight,
      sev_cases_under5_p100kU5_avert_wt = sev_cases_under5_p100kU5_avert * pop_weight,
      
    ) %>%
    dplyr::group_by(
      ISO2, sample_id, net_strategy_comp, net_name_comp, mass_int_yr_comp,
      no_future_nets_comp, net_strategy_int, net_name_int, mass_int_yr_int,
      no_future_nets_int
    ) %>%
    dplyr::summarise(
      
      pop_nat = sum(pop),
      
      n_age_under5_comp_nat = sum(n_age_under5_comp_wt),
      
      adj_ann_total_nets_dist_comp_nat = sum(adj_ann_total_nets_dist_comp_wt),
      
      cases_all_ages_comp_nat = sum(cases_all_ages_comp_wt), 
      cases_under5_comp_nat = sum(cases_under5_comp_wt), 
      cases_under5_p100kU5_comp_nat = sum(cases_under5_p100kU5_comp_wt), 
      
      clin_cases_all_ages_comp_nat = sum(clin_cases_under5_comp_wt),
      clin_cases_under5_comp_nat = sum(clin_cases_under5_comp_wt),
      clin_cases_under5_p100kU5_comp_nat = sum(clin_cases_under5_p100kU5_comp_wt), 
      
      sev_cases_all_ages_comp_nat = sum(sev_cases_all_ages_comp_wt), 
      sev_cases_under5_comp_nat = sum(sev_cases_under5_comp_wt),
      sev_cases_under5_p100kU5_comp_nat = sum(sev_cases_under5_p100kU5_comp_wt), 
      
      adj_ann_total_nets_dist_int_nat = sum(adj_ann_total_nets_dist_int_wt), 
      
      cases_all_ages_int_nat = sum(cases_all_ages_int_wt), 
      cases_under5_int_nat = sum(cases_under5_int_wt),
      cases_under5_p100kU5_int_nat = sum(cases_under5_p100kU5_int_wt),
      
      clin_cases_all_ages_int_nat = sum(clin_cases_all_ages_int_wt),        
      clin_cases_under5_int_nat = sum(clin_cases_under5_int_wt),
      clin_cases_under5_p100kU5_int_nat = sum(clin_cases_under5_p100kU5_int_wt), 
      
      sev_cases_all_ages_int_nat = sum(sev_cases_all_ages_int_wt), 
      sev_cases_under5_int_nat = sum(sev_cases_under5_int_wt), 
      sev_cases_under5_p100kU5_int_nat = sum(sev_cases_under5_p100kU5_int_wt), 
      
      adj_ann_total_nets_dist_add_nat = sum(adj_ann_total_nets_dist_add_wt),  
      
      cases_all_ages_avert_nat = sum(cases_all_ages_avert_wt), 
      cases_under5_avert_nat = sum(cases_under5_avert_wt),
      cases_under5_p100kU5_avert_nat = sum(cases_under5_p100kU5_avert_wt),
      
      clin_cases_all_ages_avert_nat = sum(clin_cases_all_ages_avert_wt), 
      clin_cases_under5_avert_nat = sum(clin_cases_under5_avert_wt),
      clin_cases_under5_p100kU5_avert_nat = sum(clin_cases_under5_p100kU5_avert_wt),
      
      sev_cases_all_ages_avert_nat = sum(sev_cases_all_ages_avert_wt), 
      sev_cases_under5_avert_nat = sum(sev_cases_under5_avert_wt),
      sev_cases_under5_p100kU5_avert_nat = sum(sev_cases_under5_p100kU5_avert_wt),
      
    ) %>%
    dplyr::mutate(
      
      adj_ann_total_nets_dist_propinc_nat =
        adj_ann_total_nets_dist_add_nat / adj_ann_total_nets_dist_comp_nat,
      
      cases_all_ages_propdec_nat =
        cases_all_ages_avert_nat / cases_all_ages_comp_nat,
      cases_under5_p100kU5_propdec_nat =
        cases_under5_p100kU5_avert_nat / cases_under5_p100kU5_comp_nat,
      
      clin_cases_all_ages_propdec_nat =
        clin_cases_all_ages_avert_nat / clin_cases_all_ages_comp_nat,
      clin_cases_under5_p100kU5_propdec_nat =
        clin_cases_under5_p100kU5_avert_nat / clin_cases_under5_p100kU5_comp_nat,
      
      sev_cases_all_ages_propdec_nat =
        sev_cases_all_ages_avert_nat / sev_cases_all_ages_comp_nat,
      sev_cases_under5_p100kU5_propdec_nat =
        sev_cases_under5_p100kU5_avert_nat / sev_cases_under5_p100kU5_comp_nat,
      
    )
    
  message("Country values calculated. Assigning to ", ctry_comparison_name)
  assign(ctry_comparison_name, ctry_comb_data, envir = .GlobalEnv)
  
}

summarise_full_comparison <- function(iso2) {
  
  data_comparison_name <- paste0(iso2, "_comp_data")
  sum_comparison_name <- paste0(iso2, "_comp_sum")
  
  if (!exists(data_comparison_name, envir = .GlobalEnv)) {
    stop("Data object ", data_comparison_name,
         " does not exist in the global environment.")
  }
  
  all_comb_sim_data <- get(data_comparison_name, envir = .GlobalEnv)
  
  all_comb_sim_sum <- all_comb_sim_data %>%
    dplyr::group_by(
      ISO2, fs_name_1, urbanicity, fs_area, fs_area_id, new_area_id, pop,
      pop_weight, net_strategy_comp, net_name_comp, mass_int_yr_comp,
      no_future_nets_comp, net_strategy_int, net_name_int, mass_int_yr_int,
      no_future_nets_int
    ) %>%
    dplyr::summarise(
      
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
      
      # cases_under5_avert_med = median(
      #   cases_under5_avert, na.rm = TRUE
      # ),
      # cases_under5_avert_lo = quantile(
      #   cases_under5_avert, 0.025, na.rm = TRUE
      # ),
      # cases_under5_avert_hi = quantile(
      #   cases_under5_avert, 0.975, na.rm = TRUE
      # ),
      
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
      
      # clin_cases_under5_avert_med = median(
      #   clin_cases_under5_avert, na.rm = TRUE
      # ),
      # clin_cases_under5_avert_lo = quantile(
      #   clin_cases_under5_avert, 0.025, na.rm = TRUE
      # ),
      # clin_cases_under5_avert_hi = quantile(
      #   clin_cases_under5_avert, 0.975, na.rm = TRUE
      # ),
      
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
      
      # sev_cases_under5_avert_med = median(
      #   sev_cases_under5_avert, na.rm = TRUE
      # ),
      # sev_cases_under5_avert_lo = quantile(
      #   sev_cases_under5_avert, 0.025, na.rm = TRUE
      # ),
      # sev_cases_under5_avert_hi = quantile(
      #   sev_cases_under5_avert, 0.975, na.rm = TRUE
      # ),
      
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
      
      .groups = "drop"
    )
  
  message("Comparison summary calculated. Assigning to ", sum_comparison_name)
  assign(sum_comparison_name, all_comb_sim_sum, envir = .GlobalEnv)

}

summarise_ctry_comparison <- function(iso2) {
  
  data_comparison_name <- paste0(iso2, "_ctry_comp_data")
  sum_comparison_name <- paste0(iso2, "_ctry_comp_sum")
  
  if (!exists(data_comparison_name, envir = .GlobalEnv)) {
    stop("Data object ", data_comparison_name,
         " does not exist in the global environment.")
  }
  
  all_comb_sim_data <- get(data_comparison_name, envir = .GlobalEnv)
  
  all_comb_sim_sum <- all_comb_sim_data %>%
    dplyr::group_by(
      ISO2, net_strategy_comp, net_name_comp, mass_int_yr_comp,
      no_future_nets_comp, net_strategy_int, net_name_int, mass_int_yr_int,
      no_future_nets_int
    ) %>%
    dplyr::summarise(
      
      adj_ann_total_nets_dist_add_med = median(
        adj_ann_total_nets_dist_add_nat, na.rm = TRUE
      ),
      adj_ann_total_nets_dist_add_lo = quantile(
        adj_ann_total_nets_dist_add_nat, 0.025, na.rm = TRUE
      ),
      adj_ann_total_nets_dist_add_hi = quantile(
        adj_ann_total_nets_dist_add_nat, 0.975, na.rm = TRUE
      ),
      
      adj_ann_total_nets_dist_propinc_med = median(
        adj_ann_total_nets_dist_propinc_nat, na.rm = TRUE
      ),
      adj_ann_total_nets_dist_propinc_lo = quantile(
        adj_ann_total_nets_dist_propinc_nat, 0.025, na.rm = TRUE
      ),
      adj_ann_total_nets_dist_propinc_hi = quantile(
        adj_ann_total_nets_dist_propinc_nat, 0.975, na.rm = TRUE
      ),
      
      cases_all_ages_avert_med = median(
        cases_all_ages_avert_nat, na.rm = TRUE
      ),
      cases_all_ages_avert_lo = quantile(
        cases_all_ages_avert_nat, 0.025, na.rm = TRUE
      ),
      cases_all_ages_avert_hi = quantile(
        cases_all_ages_avert_nat, 0.975, na.rm = TRUE
      ),
      
      cases_all_ages_propdec_med = median(
        cases_all_ages_propdec_nat, na.rm = TRUE
      ),
      cases_all_ages_propdec_lo = quantile(
        cases_all_ages_propdec_nat, 0.025, na.rm = TRUE
      ),
      cases_all_ages_propdec_hi = quantile(
        cases_all_ages_propdec_nat, 0.975, na.rm = TRUE
      ),
      
      # cases_under5_avert_med = median(
      #   cases_under5_avert_nat, na.rm = TRUE
      # ),
      # cases_under5_avert_lo = quantile(
      #   cases_under5_avert_nat, 0.025, na.rm = TRUE
      # ),
      # cases_under5_avert_hi = quantile(
      #   cases_under5_avert_nat, 0.975, na.rm = TRUE
      # ),
      
      cases_under5_p100kU5_avert_med = median(
        cases_under5_p100kU5_avert_nat, na.rm = TRUE
      ),
      cases_under5_p100kU5_avert_lo = quantile(
        cases_under5_p100kU5_avert_nat, 0.025, na.rm = TRUE
      ),
      cases_under5_p100kU5_avert_hi = quantile(
        cases_under5_p100kU5_avert_nat, 0.975, na.rm = TRUE
      ),
      
      cases_under5_p100kU5_propdec_med = median(
        cases_under5_p100kU5_propdec_nat, na.rm = TRUE
      ),
      cases_under5_p100kU5_propdec_lo = quantile(
        cases_under5_p100kU5_propdec_nat, 0.025, na.rm = TRUE
      ),
      cases_under5_p100kU5_propdec_hi = quantile(
        cases_under5_p100kU5_propdec_nat, 0.975, na.rm = TRUE
      ),
      
      clin_cases_all_ages_avert_med = median(
        clin_cases_all_ages_avert_nat, na.rm = TRUE
      ),
      clin_cases_all_ages_avert_lo = quantile(
        clin_cases_all_ages_avert_nat, 0.025, na.rm = TRUE
      ),
      clin_cases_all_ages_avert_hi = quantile(
        clin_cases_all_ages_avert_nat, 0.975, na.rm = TRUE
      ),
      
      clin_cases_all_ages_propdec_med = median(
        clin_cases_all_ages_propdec_nat, na.rm = TRUE
      ),
      clin_cases_all_ages_propdec_lo = quantile(
        clin_cases_all_ages_propdec_nat, 0.025, na.rm = TRUE
      ),
      clin_cases_all_ages_propdec_hi = quantile(
        clin_cases_all_ages_propdec_nat, 0.975, na.rm = TRUE
      ),
      
      # clin_cases_under5_avert_med = median(
      #   clin_cases_under5_avert_nat, na.rm = TRUE
      # ),
      # clin_cases_under5_avert_lo = quantile(
      #   clin_cases_under5_avert_nat, 0.025, na.rm = TRUE
      # ),
      # clin_cases_under5_avert_hi = quantile(
      #   clin_cases_under5_avert_nat, 0.975, na.rm = TRUE
      # ),
      
      clin_cases_under5_p100kU5_avert_med = median(
        clin_cases_under5_p100kU5_avert_nat, na.rm = TRUE
      ),
      clin_cases_under5_p100kU5_avert_lo = quantile(
        clin_cases_under5_p100kU5_avert_nat, 0.025, na.rm = TRUE
      ),
      clin_cases_under5_p100kU5_avert_hi = quantile(
        clin_cases_under5_p100kU5_avert_nat, 0.975, na.rm = TRUE
      ),
      
      clin_cases_under5_p100kU5_propdec_med = median(
        clin_cases_under5_p100kU5_propdec_nat, na.rm = TRUE
      ),
      clin_cases_under5_p100kU5_propdec_lo = quantile(
        clin_cases_under5_p100kU5_propdec_nat, 0.025, na.rm = TRUE
      ),
      clin_cases_under5_p100kU5_propdec_hi = quantile(
        clin_cases_under5_p100kU5_propdec_nat, 0.975, na.rm = TRUE
      ),
      
      sev_cases_all_ages_avert_med = median(
        sev_cases_all_ages_avert_nat, na.rm = TRUE
      ),
      sev_cases_all_ages_avert_lo = quantile(
        sev_cases_all_ages_avert_nat, 0.025, na.rm = TRUE
      ),
      sev_cases_all_ages_avert_hi = quantile(
        sev_cases_all_ages_avert_nat, 0.975, na.rm = TRUE
      ),
      
      sev_cases_all_ages_propdec_med = median(
        sev_cases_all_ages_propdec_nat, na.rm = TRUE
      ),
      sev_cases_all_ages_propdec_lo = quantile(
        sev_cases_all_ages_propdec_nat, 0.025, na.rm = TRUE
      ),
      sev_cases_all_ages_propdec_hi = quantile(
        sev_cases_all_ages_propdec_nat, 0.975, na.rm = TRUE
      ),
      
      # sev_cases_under5_avert_med = median(
      #   sev_cases_under5_avert_nat, na.rm = TRUE
      # ),
      # sev_cases_under5_avert_lo = quantile(
      #   sev_cases_under5_avert_nat, 0.025, na.rm = TRUE
      # ),
      # sev_cases_under5_avert_hi = quantile(
      #   sev_cases_under5_avert_nat, 0.975, na.rm = TRUE
      # ),
      
      sev_cases_under5_p100kU5_avert_med = median(
        sev_cases_under5_p100kU5_avert_nat, na.rm = TRUE
      ),
      sev_cases_under5_p100kU5_avert_lo = quantile(
        sev_cases_under5_p100kU5_avert_nat, 0.025, na.rm = TRUE
      ),
      sev_cases_under5_p100kU5_avert_hi = quantile(
        sev_cases_under5_p100kU5_avert_nat, 0.975, na.rm = TRUE
      ),
      
      sev_cases_under5_p100kU5_propdec_med = median(
        sev_cases_under5_p100kU5_propdec_nat, na.rm = TRUE
      ),
      sev_cases_under5_p100kU5_propdec_lo = quantile(
        sev_cases_under5_p100kU5_propdec_nat, 0.025, na.rm = TRUE
      ),
      sev_cases_under5_p100kU5_propdec_hi = quantile(
        sev_cases_under5_p100kU5_propdec_nat, 0.975, na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  message("Comparison summary calculated. Assigning to ", sum_comparison_name)
  assign(sum_comparison_name, all_comb_sim_sum, envir = .GlobalEnv)
  
}

shorten_sim_data("BF")
calculate_cases_per_100k_under5("BF")
create_full_sim_comparison("BF")
create_ctry_sim_comparison("BF")
summarise_full_comparison("BF")
summarise_ctry_comparison("BF")

write.csv(BF_comp_sum, "BF_comp_sum.csv", row.names = FALSE)
write.csv(BF_ctry_comp_sum, "BF_ctry_comp_sum.csv", row.names = FALSE)


shorten_sim_data("ML")
calculate_cases_per_100k_under5("ML")
create_full_sim_comparison("ML")
create_ctry_sim_comparison("ML")
summarise_full_comparison("ML")
summarise_ctry_comparison("ML")

write.csv(ML_comp_sum, "ML_comp_sum.csv", row.names = FALSE)
write.csv(ML_ctry_comp_sum, "ML_ctry_comp_sum.csv", row.names = FALSE)
