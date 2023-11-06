#plotting.R

# Function to identify 
plot_elements <- function(N_dens = 1,
                          colvals = "black",
                          plot_step_dens = FALSE,
                          plot_smth_dens = FALSE,
                          plot_modes = FALSE,
                          plot_antimodes = FALSE,
                          plot_mdc_pts = FALSE) {
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
      geom_point(aes(y = modes_val_1), alpha = 0.6, color = colvals[1], size = 2),
      if (N_dens == 2)
        geom_point(aes(y = modes_val_2), alpha = 0.6, color = colvals[2], size = 2),
    if (plot_antimodes) 
      geom_point(aes(y = antimodes_val_1), alpha = 0.6, color = colvals[1], size = 2, shape = 15),
      if (N_dens == 2)
        geom_point(aes(y = antimodes_val_2), alpha = 0.6, color = colvals[2], size = 2, shape = 15),
    if (plot_mdc_pts) 
      geom_point(aes(y = mdc_val), alpha = 0.6, color = "black", size = 2, shape = 17)
  )
}

# Function to generate plots
generate_MDC_plots <- function(dataset,
                               densities,
                               cap_extreme,
                               N_dens,
                               colvals,
                               plot_step_dens,
                               plot_smth_dens,
                               plot_modes,
                               plot_antimodes,
                               plot_mdc_pts,
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
      } else {
        plt_df <- rbind.data.frame(plt_df, dataset[dataset$area_id == i,])
      }
    }
  }
  
  # Prepare data frame for plotting
  plt_df$countryname <- countrycode(plt_df$ISO2,
                                    origin = 'iso2c',
                                    destination = 'country.name')
  
  # Plotting Loop
  if (MDC_kde_national) {pages <- 1} else {pages <- 1:N_ISO2}
  for (i in pages) {
    
    if (MDC_kde_national) {
      plt_this <- plt_df
    } else {
      plt_this <- plt_df[which(plt_df$ISO2 == uni_ISO2[i]),]
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
                          plot_mdc_pts)) +
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
      pdf(filename_ctry, width = 7.8, height = 11.2, paper="a4")
      print(facet_nets)
      dev.off()
    } else {
      print(facet_nets)
    }
    
  }
}


plot_MDCs <- function(dataset,
                      densities = NULL,
                      colvals = NULL,
                      cap_extreme = FALSE,
                      plot_step_dens = FALSE,
                      plot_smth_dens = FALSE,
                      plot_modes = FALSE,
                      plot_antimodes = FALSE,
                      plot_mdc_pts = FALSE,
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
  if (plot_mdc_pts) {
    smth_name <- paste0("smth_", densities[MDC_density_number])
    MDC_val_den <- data.frame("mdc_val" = dataset[, smth_name] * dataset$mdc)
  }
  
  # Change zero values to NA
  if (plot_modes) {mode_val_dens[mode_val_dens == 0] <- NA}
  if (plot_antimodes) {antimode_val_dens[antimode_val_dens == 0] <- NA}
  if (plot_mdc_pts) {MDC_val_den[MDC_val_den == 0] <- NA}
  
  # Combine new data frames with original dataset
  if (plot_step_dens) {dataset <- cbind.data.frame(dataset, step_dens)}
  if (plot_smth_dens) {dataset <- cbind.data.frame(dataset, smth_dens)}
  if (plot_modes) {dataset <- cbind.data.frame(dataset, mode_val_dens)}
  if (plot_antimodes) {dataset <- cbind.data.frame(dataset, antimode_val_dens)}
  if (plot_mdc_pts) {dataset <- cbind.data.frame(dataset, MDC_val_den)}
  
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
                       cap_extreme,
                       N_dens,
                       colvals,
                       plot_step_dens,
                       plot_smth_dens,
                       plot_modes,
                       plot_antimodes,
                       plot_mdc_pts,
                       MDC_density_number,
                       folderpath,
                       fileprefix,
                       timestamp,
                       multi_plt = FALSE)
    dev.off()
  } else {
    filename_all <- paste0(folderpath, fileprefix, timestamp, "_all.pdf")
    #single pdf plot
    pdf(filename_all, width = 7.8, height = 11.2, paper="a4")
    generate_MDC_plots(dataset,
                       densities,
                       cap_extreme,
                       N_dens,
                       colvals,
                       plot_step_dens,
                       plot_smth_dens,
                       plot_modes,
                       plot_antimodes,
                       plot_mdc_pts,
                       MDC_density_number,
                       folderpath,
                       fileprefix,
                       timestamp,
                       multi_plt = FALSE)
    dev.off()
    #multiple pdf plots
    generate_MDC_plots(dataset,
                       densities,
                       cap_extreme,
                       N_dens,
                       colvals,
                       plot_step_dens,
                       plot_smth_dens,
                       plot_modes,
                       plot_antimodes,
                       plot_mdc_pts,
                       MDC_density_number,
                       folderpath,
                       fileprefix,
                       timestamp,
                       multi_plt = TRUE)
  }
  
  return(NULL)
}