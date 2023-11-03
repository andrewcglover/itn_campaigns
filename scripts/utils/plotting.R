#plotting.R


# Function to generate plots
generate_MDC_plots <- function(dataset,
                               net_density_name,
                               smth_density_name,
                               ref_density_name,
                               folderpath,
                               fileprefix,
                               timestamp,
                               multi_plt) {
  
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
  trimmed_dataset <- NULL
  if (MDC_kde_national) {
    # Changes required before using previous implementation
  } else {
    for (i in 1:N_areas) {
      t0 <- extreme_nets$min_rec[i]
      tm <- extreme_nets$max_rec[i]
      trimmed_dataset <- rbind.data.frame(trimmed_dataset,
                                          dataset[dataset$area_id == i &
                                                    !(dataset$CMC < t0
                                                      | dataset$CMC > tm),])
    }
  }
  
  # Prepare data frame for plotting
  trimmed_dataset$countryname <- countrycode(dataset$ISO2,
                                             origin = 'iso2c',
                                             destination = 'country.name')
  plt_df <- data.frame("ISO2" = trimmed_dataset$ISO2,
                       "Admin" = trimmed_dataset$ADM1,
                       "Urbanicity" = trimmed_dataset$urbanicity,
                       "Country" = trimmed_dataset$countryname,
                       "CMC" = trimmed_dataset$CMC,
                       "MDC" = trimmed_dataset$mdc_points,
                       "Net_density" = trimmed_dataset[, net_density_name],
                       "Smooth_density" = trimmed_dataset[, smth_density_name],
                       "Ref_density" = trimmed_dataset[, ref_density_name])
  
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
    net_plt <- ggplot() +
      geom_step(data = plt_this,
                aes(x = CMC, y = Ref_density),
                alpha = 1, color = "darkorange", size = 0.5) +
      geom_step(data = plt_this,
                aes(x = CMC, y = Net_density),
                alpha = 1, color = "royalblue3", size = 0.5) +
      geom_area(data = plt_this,
                aes(x = CMC, y = Smooth_density),
                color = NA, alpha = 0.5, fill = "royalblue1") +
      # geom_path(data = plt_df,
      #           aes(x = CMC, y = smth_dhs_norm),
      #           color = "black") +
      geom_point(data = plt_this,
                 aes(x = CMC, y = MDC),
                 color = "black", size = 2) +
      #scale_y_continuous(limits = c(0, 0.8)) +
      scale_x_continuous(breaks = xlbs_ids, labels = xlbs) +
      ylab("Density") + 
      xlab("Year") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
    
    # Plot facets
    if (MDC_kde_national) {
      facet_nets <- net_plt + facet_wrap(~Country,  ncol=1, labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free_y", strip.position = "right")
    } else {
      if (urban_split_MDC) {
        facet_nets <- net_plt + facet_grid(Admin ~ Urbanicity, labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free_y") + ggtitle(ttl_nm)# + geom_blank(data = areas_plt_df, aes(y = ylim))
      } else {
        facet_nets <- net_plt + facet_wrap(~Admin,  ncol=1, labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free_y", strip.position = "right") + theme(strip.text.y = element_text(size = 7)) + ggtitle(ttl_nm)
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
                      net_density_name = "norm_fit_nets",
                      smth_density_name = "norm_smth_nets",
                      ref_density_name = "norm_ref_nets") {
  
  if (MDC_kde_national) {print("warning: not tested for MDC_kde_national")}
  
  # Append MDC points for plotting
  dataset$mdc_points <- dataset[, smth_density_name] * dataset$mdc
  dataset$mdc_points[dataset$mdc_points == 0] <- NA
  
  #Bound time series by first and last net recorded
  if (urban_split_MDC) {
    dataset_filter <- dataset
  } else {
    dataset_filter <- dataset[which(dataset$urbanicity=="urban"),]
  }
  
  plt_ids <- NULL
  for (i in 1:N_areas) {
    plt_ids <- c(plt_ids,
                 which((dataset_filter$area_id == i &
                          !((dataset_filter$CMC < extreme_nets$min_rec[i]) |
                              (dataset_filter$CMC > extreme_nets$max_rec[i])))))
  }
  dataset_filter <- dataset_filter[plt_ids,]
  
  
  #filename
  folderpath <- "./outputs/mdc_timings/"
  timestamp <- format(Sys.time(), "%y%m%d%H%M")
  if (MDC_kde_national) {
    fileprefix <- "MDC_timings_ADMsyn_"
  } else {
    fileprefix <- "MDC_timings_ADMvar_"
  }
  
  #generate plots
  if (MDC_kde_national) {
    filename_all <- paste0(folderpath, fileprefix, timestamp, ".pdf")
    pdf(filename_all, width = 7.8, height = 11.2, paper="a4")
    generate_MDC_plots(dataset_filter,
                       net_density_name,
                       smth_density_name,
                       ref_density_name,
                       folderpath,
                       fileprefix,
                       timestamp,
                       multi_plt = FALSE)
    dev.off()
  } else {
    filename_all <- paste0(folderpath, fileprefix, timestamp, "_all.pdf")
    #single pdf plot
    pdf(filename_all, width = 7.8, height = 11.2, paper="a4")
    generate_MDC_plots(dataset_filter,
                       net_density_name,
                       smth_density_name,
                       ref_density_name,
                       folderpath,
                       fileprefix,
                       timestamp,
                       multi_plt = FALSE)
    dev.off()
    #multiple pdf plots
    generate_MDC_plots(dataset_filter,
                       net_density_name,
                       smth_density_name,
                       ref_density_name,
                       folderpath,
                       fileprefix,
                       timestamp,
                       multi_plt = TRUE)
  }
  
  return(NULL)
}