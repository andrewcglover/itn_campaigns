#plotting.R

# Function to identify 
plot_elements <- function(N_dens = 1,
                          colvals = "black",
                          plot_step_dens = FALSE,
                          plot_smth_dens = FALSE,
                          plot_modes = FALSE,
                          plot_antimodes = FALSE,
                          plot_mdc_pts = FALSE,
                          plot_vert_periods = FALSE,
                          plot_comparison_mdc = FALSE,
                          vert_dataset = NULL) {
  list(
    if (plot_step_dens) 
      geom_step(aes(y = den_1), alpha = 0.6, color = colvals[1], size = 0.5),
      if (N_dens == 2)
        geom_step(aes(y = den_2), alpha = 0.6, color = colvals[2], size = 0.5),
    if (plot_smth_dens) 
      geom_area(aes(y = smth_1), alpha = 0.3, color = NA, fill = colvals[1]),
      if (N_dens == 2)
        geom_area(aes(y = smth_2), alpha = 0.3, color = NA, fill = colvals[2]),
    if (plot_modes) 
      geom_point(aes(y = modes_val_1), alpha = 0.6, color = colvals[1], size = 2, shape = 15),
      if (N_dens == 2)
        geom_point(aes(y = modes_val_2), alpha = 0.6, color = colvals[2], size = 2, shape = 15),
    if (plot_antimodes) 
      geom_point(aes(y = antimodes_val_1), alpha = 0.6, color = colvals[1], size = 2, shape = 15),
      if (N_dens == 2)
        geom_point(aes(y = antimodes_val_2), alpha = 0.6, color = colvals[2], size = 2, shape = 15),
    if (plot_comparison_mdc) 
      geom_point(aes(y = cmdc_val), alpha = 0.6, color = "darkorchid2", size = 3, shape = 15),
    if (plot_mdc_pts) 
      geom_point(aes(y = mdc_val), alpha = 0.6, color = "black", size = 3, shape = 17),
    if (plot_vert_periods) 
      geom_vline(aes(xintercept = CMC_mdc_bounds),data = vert_dataset)
  )
}

# Function to generate plots
generate_MDC_plots <- function(dataset,
                               densities,
                               periods_dataset,
                               cap_extreme,
                               N_dens,
                               colvals,
                               plot_step_dens,
                               plot_smth_dens,
                               plot_modes,
                               plot_antimodes,
                               plot_mdc_pts,
                               plot_vert_periods,
                               plot_comparison_mdc,
                               MDC_density_number,
                               folderpath,
                               fileprefix,
                               timestamp,
                               multi_plt = FALSE) {
                   
  
  # Initialise values
  plt_df <- NULL
  xlbs_ids <- seq(date_to_CMC(2010,1), CMC_last,24)
  xlbs <- CMC_to_date(xlbs_ids)[,1]
  if (MDC_kde_national) {
    first_month_national <- rep(0, N_ISO2)
    last_month_national <- rep(0, N_ISO2)
    ylim_global <- 0
  } else {
    ylim_national <- rep(0, N_ISO2)
  }
  
  # Trim dataset
  plt_df <- NULL
  plt_vert <- NULL
  if (MDC_kde_national) {
    # Changes required before using previous implementation
  } else {
    
    for (i in 1:N_areas) {
      if (cap_extreme) {
        t0 <- extreme_nets$min_rec[i]
        tm <- extreme_nets$max_rec[i]
        plt_df <- rbind.data.frame(plt_df,
                                   dataset[dataset$area_id == i &
                                             !(dataset$CMC < t0
                                               | dataset$CMC > tm),])
        if (plot_vert_periods) {
          plt_vert <- rbind.data.frame(plt_vert,
                                       periods_dataset[periods_dataset$area_id == i &
                                                 !(periods_dataset$CMC_mdc_bounds < t0
                                                   | periods_dataset$CMC_mdc_bounds > tm),])
        }
      } else {
        plt_df <- rbind.data.frame(plt_df, dataset[dataset$area_id == i,])
        if (plot_vert_periods) {
          plt_vert <- rbind.data.frame(plt_vert,
                                       periods_dataset[periods_dataset$area_id == i,])
        }
      }
    }
  }
  
  # Prepare data frame for plotting
  plt_df$countryname <- countrycode(plt_df$ISO2,
                                    origin = 'iso2c',
                                    destination = 'country.name')
  if (plot_vert_periods) {
    plt_vert$countryname <- countrycode(plt_vert$ISO2,
                                        origin = 'iso2c',
                                        destination = 'country.name')
  }
  
  # Plotting Loop
  if (MDC_kde_national) {pages <- 1} else {pages <- 1:N_ISO2}
  for (i in pages) {
    
    if (MDC_kde_national) {
      plt_this <- plt_df
      plt_this2 <- plt_vert
    } else {
      plt_this <- plt_df[which(plt_df$ISO2 == uni_ISO2[i]),]
      plt_this2 <- plt_vert[which(plt_vert$ISO2 == uni_ISO2[i]),]
    }
    
    # Title
    if (!MDC_kde_national) {
      ttl_nm <- countrycode(SSA_ISO2[i],
                            origin = 'iso2c',
                            destination = 'country.name')
    }
    
    # Y axis aesthetics
    if (MDC_kde_national) {
      ylim <- ylim_global
    } else {
      ylim <- ylim_national[i]
    }
    ylbs_ids <- seq(0, ylim, 0.2)
    ylbs <- as.character(ylbs_ids)
    ylbs[length(ylbs)] <- paste0(">",ylbs[length(ylbs)])
    
    # Plotting object
    base_plt <- ggplot(plt_this, aes(x = CMC))
    plt_obj <- base_plt +
      print(plot_elements(N_dens,
                          colvals,
                          plot_step_dens,
                          plot_smth_dens,
                          plot_modes,
                          plot_antimodes,
                          plot_mdc_pts,
                          plot_vert_periods,
                          plot_comparison_mdc,
                          plt_this2)) +
      scale_x_continuous(breaks = xlbs_ids, labels = xlbs) +
      ylab("Density") + 
      xlab("Year") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
    
    # Plot facets
    if (MDC_kde_national) {
      facet_nets <- plt_obj + facet_wrap(~countryname,  ncol=1, labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free_y", strip.position = "right")
    } else {
      if (urban_split_MDC) {
        facet_nets <- plt_obj + facet_grid(ADM1 ~ urbanicity, labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free_y") + ggtitle(ttl_nm)# + geom_blank(data = areas_plt_df, aes(y = ylim))
      } else {
        facet_nets <- plt_obj + facet_wrap(.~ADM1,  ncol=1, labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free_y", strip.position = "right") + theme(strip.text.y = element_text(size = 7)) + ggtitle(ttl_nm)
      }
    }
    
    # Save plot
    if (multi_plt & !MDC_kde_national) {
      filename_ctry <- paste0(folderpath, fileprefix, timestamp, "_", uni_ISO2[i], ".pdf")
      # pdf(filename_ctry, width = 7.8, height = 11.2, paper="a4")
      pdf(filename_ctry, width = 5, height = 15)#, paper="a4")
      print(facet_nets)
      dev.off()
    } else {
      print(facet_nets)
    }
    
  }
}


plot_MDCs <- function(dataset,
                      densities = NULL,
                      periods_dataset = NULL,
                      colvals = NULL,
                      cap_extreme = FALSE,
                      plot_step_dens = FALSE,
                      plot_smth_dens = FALSE,
                      plot_modes = FALSE,
                      plot_antimodes = FALSE,
                      plot_mdc_pts = FALSE,
                      plot_vert_periods = FALSE,
                      plot_comparison_mdc = FALSE,
                      MDC_density_number = 1) {
  
  # Number of densities
  N_dens <- length(densities)
  
  # Warning - only a maximum of two density types will be plotted
  if (N_dens > 2) {print("warning: maximum of two densities expected")}
  
  # Warning - current version not tested for synchronised national MDCs
  if (MDC_kde_national) {print("warning: not tested for MDC_kde_national")}
  
  # Return error if modes, antimodes or MDCs plotted without smooth density
  if ((plot_modes|plot_antimodes|plot_mdc_pts) & !plot_smth_dens) {
    print("warning: smooth densities expected for mode, antimode or MDC plots")
  }
  
  # Create new data frames of densities for plotting
  N_tot <- dim(dataset)[1]
  if (plot_step_dens) {step_dens <- data.frame(NULL)}
  if (plot_smth_dens) {smth_dens <- data.frame(NULL)}
  if (plot_modes) {mode_dens <- data.frame(NULL)}
  if (plot_modes) {mode_val_dens <- data.frame(NULL)}
  if (plot_antimodes) {antimode_dens <- data.frame(NULL)}
  if (plot_antimodes) {antimode_val_dens <- data.frame(NULL)}
  for (i in 1:N_dens) {
    if (plot_step_dens) {
      step_dens[1:N_tot, paste0("den_",i)] <- dataset[, densities[i]]
    }
    if (plot_smth_dens) {
        smth_name <- paste0("smth_", densities[i])
        smth_dens[1:N_tot, paste0("smth_",i)] <- dataset[, smth_name]
      if (plot_modes) {
        mode_name <- paste0("modes_",densities[i])
        mode_dens[1:N_tot, paste0("modes_",i)] <- dataset[, mode_name]
        mode_val_name <- paste0("modes_val_", densities[i])
        mode_val_dens[1:N_tot, paste0("modes_val_",i)] <- smth_dens[1:N_tot, paste0("smth_",i)] *
          mode_dens[1:N_tot, paste0("modes_",i)]
      }
      if (plot_antimodes) {
        antimode_name <- paste0("antimodes_", densities[i])
        antimode_dens[1:N_tot, paste0("antimodes_",i)] <- dataset[, antimode_name]
        antimode_val_name <- paste0("antimodes_val_", densities[i])
        antimode_val_dens[1:N_tot, paste0("antimodes_val_",i)] <- smth_dens[1:N_tot, paste0("smth_",i)] *
          antimode_dens[1:N_tot, paste0("antimodes_",i)]
      }
    }
  }
  
  # Calculate MDC point values to dataset
  smth_name <- paste0("smth_", densities[MDC_density_number])
  if (plot_mdc_pts) {
    # smth_name <- paste0("smth_", densities[MDC_density_number])
    MDC_val_den <- data.frame("mdc_val" = dataset[, smth_name] * dataset$mdc)
  }
  if (plot_comparison_mdc) {
    # smth_name <- paste0("smth_", densities[MDC_density_number])
    cMDC_val_den <- data.frame("cmdc_val" = dataset[, smth_name] * dataset$cmdc)
  }
  
  # Change zero values to NA
  if (plot_modes) {mode_val_dens[mode_val_dens == 0] <- NA}
  if (plot_antimodes) {antimode_val_dens[antimode_val_dens == 0] <- NA}
  if (plot_mdc_pts) {MDC_val_den[MDC_val_den == 0] <- NA}
  if (plot_comparison_mdc) {cMDC_val_den[cMDC_val_den == 0] <- NA}
  
  # Combine new data frames with original dataset
  if (plot_step_dens) {dataset <- cbind.data.frame(dataset, step_dens)}
  if (plot_smth_dens) {dataset <- cbind.data.frame(dataset, smth_dens)}
  if (plot_modes) {dataset <- cbind.data.frame(dataset, mode_val_dens)}
  if (plot_antimodes) {dataset <- cbind.data.frame(dataset, antimode_val_dens)}
  if (plot_mdc_pts) {dataset <- cbind.data.frame(dataset, MDC_val_den)}
  if (plot_comparison_mdc) {dataset <- cbind.data.frame(dataset, cMDC_val_den)}
  
  #Bound time series by first and last net recorded
  if (!urban_split_MDC) {
    dataset <- dataset[which(dataset$urbanicity=="urban"),]
  }
  plt_ids <- NULL
  for (i in 1:N_areas) {
    if (cap_extreme) {
      plt_ids <- c(plt_ids,
                   which((dataset$area_id == i &
                            !((dataset$CMC < extreme_nets$min_rec[i]) |
                                (dataset$CMC > extreme_nets$max_rec[i])))))
    } else {
      plt_ids <- c(plt_ids, which(dataset$area_id == i))
    }
  }
  dataset <- dataset[plt_ids,]
  
  #filename
  folderpath <- "./outputs/mdc_timings/"
  if (MDC_kde_national) {
    fileprefix <- "MDC_timings_ADMsyn_"
  } else {
    fileprefix <- "MDC_timings_ADMvar_"
  }
  
  #generate plots
  if (MDC_kde_national) {
    filename_all <- paste0(folderpath, fileprefix, timestamp, ".pdf")
    pdf(filename_all, width = 7.8, height = 11.2, paper="a4")
    generate_MDC_plots(dataset,
                       densities,
                       periods_dataset,
                       cap_extreme,
                       N_dens,
                       colvals,
                       plot_step_dens,
                       plot_smth_dens,
                       plot_modes,
                       plot_antimodes,
                       plot_mdc_pts,
                       plot_vert_periods,
                       plot_comparison_mdc,
                       MDC_density_number,
                       folderpath,
                       fileprefix,
                       timestamp,
                       multi_plt = FALSE)
    dev.off()
  } else {
    filename_all <- paste0(folderpath, fileprefix, timestamp, "_all.pdf")
    #single pdf plot
    #pdf(filename_all, width = 7.8, height = 11.2, paper="a4")
    pdf(filename_all, width = 5.5, height = 10)
    generate_MDC_plots(dataset,
                       densities,
                       periods_dataset,
                       cap_extreme,
                       N_dens,
                       colvals,
                       plot_step_dens,
                       plot_smth_dens,
                       plot_modes,
                       plot_antimodes,
                       plot_mdc_pts,
                       plot_vert_periods,
                       plot_comparison_mdc,
                       MDC_density_number,
                       folderpath,
                       fileprefix,
                       timestamp,
                       multi_plt = FALSE)
    dev.off()
    #multiple pdf plots
    generate_MDC_plots(dataset,
                       densities,
                       periods_dataset,
                       cap_extreme,
                       N_dens,
                       colvals,
                       plot_step_dens,
                       plot_smth_dens,
                       plot_modes,
                       plot_antimodes,
                       plot_mdc_pts,
                       plot_vert_periods,
                       plot_comparison_mdc,
                       MDC_density_number,
                       folderpath,
                       fileprefix,
                       timestamp,
                       multi_plt = TRUE)
  }
  
  return(NULL)
}

#-------------------------------------------------------------------------------

# Example net ages
plot_net_ages <- function(area_id_input = 1, use_access = TRUE) {
  
  if (use_access) {
    dataset <- access_nets_weighted
    meanlife <- prior_mean_access_meanlife[area_id_input]
    meanlife_samples <- (access_decay_samples$inv_lambda[,area_id_input])
  } else {
    dataset <- used_nets_weighted
    meanlife <- prior_mean_used_meanlife[area_id_input]
    meanlife_samples <- (used_decay_samples$inv_lambda[,area_id_input])
  }
  
  area_ids <- which(dataset$area_id == area_id_input)
  area_data <- dataset[area_ids,]
  area_ages <- area_data$months_since_obtained
  
  ages <- seq(0,37)
  N_ages <- length(ages)
  raw_den <- rep(0, N_ages)
  fit_den <- rep(0, N_ages)

  for (i in 1:N_ages) {
    raw_den[i] <- sum(area_ages == ages[i])
  }
  
  norm_den <- raw_den / sum(raw_den)
  
  # Fit distribution
  meanlife_CrI <- quantile(meanlife_samples, c(0.025, 0.975))
  meanlife_LB <- meanlife_CrI[[1]]
  meanlife_UB <- meanlife_CrI[[2]]
  lambda_LB <- 1 / meanlife_UB
  lambda_UB <- 1 / meanlife_LB
  lambda <- 1 / meanlife
  
  fit_LB <- lambda_UB * exp(-lambda_UB*ages)
  fit_UB <- lambda_LB * exp(-lambda_LB*ages)
  fit_den <- lambda * exp(-lambda*ages)
  
  plt_data <- data.frame("Age" = ages,
                         "Density" = norm_den,
                         "Fit" = fit_den,
                         "Fit_LB" = fit_LB,
                         "Fit_UB" = fit_UB)
  
  ggplot(plt_data, aes(x = Age)) +
    geom_step(aes(y = Density), alpha = 0.8, color = "royalblue4", size = 1) +
    geom_line(aes(y = Fit), alpha = 0.8, color = "maroon", size = 1) +
    geom_ribbon(aes(ymin=Fit_LB, ymax=Fit_UB), fill = "maroon", alpha=0.2) +
    ylim(0,0.21)

  #print(plt)
  #show(plt)
  
  
}

# Example net ages
plot_country_net_ages <- function(cc = "SN", use_access = TRUE) {
  
  country_id <- unique(area_link$CTRY[which(area_link$ISO2==cc)])

  if (use_access) {
    dataset <- access_nets_weighted
    mu_meanlife_samples <- (access_decay_samples$mu_c[,country_id])
    sigma_meanlife_samples <- (access_decay_samples$sigma_c[,country_id])
  } else {
    dataset <- used_nets_weighted
    mu_meanlife_samples <- (used_decay_samples$mu_c[,country_id])
    sigma_meanlife_samples <- (used_decay_samples$sigma_c[,country_id])
  }
  meanlife <- mean(mu_meanlife_samples)
  
  N_samples <- length(mu_meanlife_samples)
  meanlife_samples <- rep(NA, N_samples)
  for (i in 1:N_samples) {
    meanlife_samples[i] <- rnorm(1, mean = mu_meanlife_samples[i],
                                 sd = sigma_meanlife_samples[i])
  }

  ctry_ids <- which(dataset$CTRY == country_id)
  ctry_data <- dataset[ctry_ids,]
  ctry_ages <- ctry_data$months_since_obtained
  
  ages <- seq(0,37)
  N_ages <- length(ages)
  raw_den <- rep(0, N_ages)
  fit_den <- rep(0, N_ages)
  
  for (i in 1:N_ages) {
    raw_den[i] <- sum(ctry_ages == ages[i])
  }
  
  norm_den <- raw_den / sum(raw_den)
  
  # Fit distribution
  # meanlife_CrI <- quantile(meanlife_samples, c(0.025, 0.975))
  # meanlife_LB <- meanlife_CrI[[1]]
  # meanlife_UB <- meanlife_CrI[[2]]
  # lambda_LB <- 1 / meanlife_UB
  # lambda_UB <- 1 / meanlife_LB
  lambda <- 1 / meanlife
  lambda_samples <- 1 / meanlife_samples
  
  fit_samples <- matrix(nrow = N_ages, ncol = N_samples)
  for (j in 1:N_samples) {
    fit_samples[,j] <- lambda_samples[j] * exp(-lambda_samples[j]*ages)
  }
  fit_LB <- rep(NA, N_ages)
  fit_UB <- rep(NA, N_ages)
  fit_den <- rep(NA, N_ages)
  for (i in 1:N_ages) {
    fits <- fit_samples[i,which(!is.na(fit_samples[i,]))]
    fit_den[i] <- mean(fits)
    CrI <- quantile(fits, c(0.025, 0.5, 0.975))
    fit_LB[i] <- CrI[[1]]
    fit_den[i] <- CrI[[2]]
    fit_UB[i] <- CrI[[3]]
  }
  
  plt_data <- data.frame("Age" = ages,
                         "Density" = norm_den,
                         "Fit" = fit_den,
                         "Fit_LB" = fit_LB,
                         "Fit_UB" = fit_UB)
  
  ggplot(plt_data, aes(x = Age)) +
    geom_step(aes(y = Density), alpha = 0.8, color = "royalblue4", size = 1) +
    geom_line(aes(y = Fit), alpha = 0.8, color = "maroon", size = 1) +
    geom_ribbon(aes(ymin=Fit_LB, ymax=Fit_UB), fill = "maroon", alpha=0.2) +
    ylim(0,0.21)
  
}


# Example net ages
plot_universal_net_ages <- function(use_access = TRUE) {
  
  if (use_access) {
    dataset <- access_nets_weighted
    mu_samples <- (access_decay_samples$mu_u)
    sigma_samples <- (access_decay_samples$sigma_u)
    tau_samples <- (access_decay_samples$tau_u)
    rho_samples <- (access_decay_samples$rho_u)
  } else {
    dataset <- used_nets_weighted
    mu_samples <- (used_decay_samples$mu_u)
    sigma_samples <- (used_decay_samples$sigma_u)
    tau_samples <- (used_decay_samples$tau_u)
    rho_samples <- (used_decay_samples$rho_u)
  }
  meanlife <- mean(mu_samples)
  
  N_samples <- length(mu_samples)
  mu_c_samples <- rep(NA, N_samples)
  sigma_c_samples <- rep(NA, N_samples)
  meanlife_samples <- rep(NA, N_samples)
  for (i in 1:N_samples) {
    mu_c_samples[i] <- rnorm(1, mean = mu_samples[i],
                             sd = sigma_samples[i])
    sigma_c_samples[i] <- rnorm(1, mean = tau_samples[i],
                                sd = rho_samples[i])
    sigma_c_samples[sigma_c_samples<0] <- NA
    meanlife_samples[i] <- rnorm(1, mean = mu_c_samples[i],
                                 sd = sigma_c_samples[i])
  }

  univ_ages <- dataset$months_since_obtained
  
  ages <- seq(0,37)
  N_ages <- length(ages)
  raw_den <- rep(0, N_ages)
  fit_den <- rep(0, N_ages)
  
  for (i in 1:N_ages) {
    raw_den[i] <- sum(univ_ages == ages[i])
  }
  
  norm_den <- raw_den / sum(raw_den)
  
  # Fit distribution
  # meanlife_CrI <- quantile(meanlife_samples, c(0.025, 0.975))
  # meanlife_LB <- meanlife_CrI[[1]]
  # meanlife_UB <- meanlife_CrI[[2]]
  # lambda_LB <- 1 / meanlife_UB
  # lambda_UB <- 1 / meanlife_LB
  lambda <- 1 / meanlife
  lambda_samples <- 1 / meanlife_samples
  
  fit_samples <- matrix(nrow = N_ages, ncol = N_samples)
  for (j in 1:N_samples) {
    fit_samples[,j] <- lambda_samples[j] * exp(-lambda_samples[j]*ages)
  }
  fit_LB <- rep(NA, N_ages)
  fit_UB <- rep(NA, N_ages)
  fit_den <- rep(NA, N_ages)
  for (i in 1:N_ages) {
    fits <- fit_samples[i,which(!is.na(fit_samples[i,]))]
    fit_den[i] <- mean(fits)
    CrI <- quantile(fits, c(0.025, 0.5, 0.975))
    fit_LB[i] <- CrI[[1]]
    fit_den[i] <- CrI[[2]]
    fit_UB[i] <- CrI[[3]]
  }

  plt_data <- data.frame("Age" = ages,
                         "Density" = norm_den,
                         "Fit" = fit_den,
                         "Fit_LB" = fit_LB,
                         "Fit_UB" = fit_UB)
  
  ggplot(plt_data, aes(x = Age)) +
    geom_step(aes(y = Density), alpha = 0.8, color = "royalblue4", size = 1) +
    geom_line(aes(y = Fit), alpha = 0.8, color = "maroon", size = 1) +
    geom_ribbon(aes(ymin=Fit_LB, ymax=Fit_UB), fill = "maroon", alpha=0.2) +
    ylim(0,0.21)
  
}
  