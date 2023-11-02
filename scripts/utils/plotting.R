#plotting.R

# Function for normalising densities
normalise_MDC_density <- function(region_df, first_month, last_month) {
  
  # Normalisation factor
  series_length <- last_month - first_month
  norm_fac <- series_length / 12
  
  # Normalise smooth DHS density
  region_df$smth_dhs_norm <- region_df$smth_dhs * norm_fac
  
  # Y limit for plotting (density ceiling rounded to nearest 0.2)
  adm_ylim <- ceiling(max(region_df$smth_dhs_norm) * 5) / 5
  
  # Normalise net density
  region_df$nets_norm <- region_df$scaled_camp_nets * norm_fac

  # Return region data frame and y limit
  outlist <- list(region_df, adm_ylim)
  return(outlist)
  
}

# Function to generate plots
generate_MDC_plots <- function(dataset_filter,
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
  
  # Normalisation loop over countries
  for (i in 1:N_ISO2) {

    # Select country values
    country_df <- dataset_filter[which(dataset_filter$ISO2 == uni_ISO2[i]),]
    if (!urban_split_MDC) {
      country_df <- country_df[which(country_df$urbanicity == "urban"),]
    }
    
    # Normalise for plotting
    if (MDC_kde_national) {
      country_df_no_reps <- country_df
      area1 <- which(country_df_no_reps$area_id == min(unique(country_df$area_id)))
      country_df_no_reps <- country_df_no_reps[area1,]
      first_month_national[i] <- min(areas_df$min_net_age_rec[which(areas_df$ISO2 == uni_ISO2[i])])
      last_month_national[i] <- max(areas_df$max_net_age_rec[which(areas_df$ISO2 == uni_ISO2[i])])
      norm_list <- normalise_MDC_density(country_df_no_reps,
                                         first_month_national[i],
                                         last_month_national[i])
      norm_ctry_df <- norm_list[[1]]
      ctry_ylim <- norm_list[[2]]
      # Set y limit and value after final recorded net density to zero
      #norm_ctry_df$nets_norm[which(norm_ctry_df$CMC == min(CMC_last, last_month_national[i] + 1))] <- 0
      norm_ctry_df$nets_norm[which(norm_ctry_df$CMC == first_month_national[i])] <- 0
      norm_ctry_df$nets_norm[which(norm_ctry_df$CMC == last_month_national[i])] <- 0
      if (ctry_ylim > ylim_global) {ylim_global <- ctry_ylim}
      # Combine with plotting data frame
      plt_df <- rbind.data.frame(plt_df, norm_ctry_df)
    } else {
      #plt_df <- NULL
      ctry_area_ids <- unique(country_df$area_id)
      for (id in ctry_area_ids) {
        adm_df <- country_df[which(country_df$area_id == id),]
        first_month <- areas_df$min_net_age_rec[which(areas_df$area_ID == id)]
        last_month <- areas_df$max_net_age_rec[which(areas_df$area_ID == id)]
        norm_list <- normalise_MDC_density(adm_df, first_month, last_month)
        norm_adm_df <- norm_list[[1]]
        adm_ylim <- norm_list[[2]]
        # Set y limit and value after final recorded net density to zero
        #norm_adm_df$nets_norm[which(norm_adm_df$CMC == min(CMC_last, last_month + 1))] <- 0
        norm_adm_df$nets_norm[which(norm_adm_df$CMC == first_month)] <- 0
        norm_adm_df$nets_norm[which(norm_adm_df$CMC == last_month)] <- 0
        if (adm_ylim > ylim_national[i]) {ylim_national[i] <- adm_ylim}
        # Combine with plotting data frame
        plt_df <- rbind.data.frame(plt_df, norm_adm_df)
      }
    }
  }
  
  # Bound net density for plotting
  if (MDC_kde_national) {
    plt_df$nets_norm[which(plt_df$nets_norm > ylim_global)] <- ylim_global
  } else {
    for (i in 1:N_ISO2) {
      plt_df$nets_norm[which(plt_df$ISO2 == uni_ISO2[i] &
                               plt_df$nets_norm > ylim_national[i])] <- ylim_national[i]
      
    }
  }
  
  # MDC points
  plt_df$MDC_norm <- plt_df$smth_dhs_norm * plt_df$MDC
  plt_df$MDC_norm[which(plt_df$MDC_norm == 0)] <- NA
  
  # Prepare data frame for plotting
  plt_df$countryname <- countrycode(plt_df$ISO2,
                                    origin = 'iso2c',
                                    destination = 'country.name')
  plt_this <- plt_df
  
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
                aes(x = CMC, y = nets_norm),
                alpha = 1, color = "#458AFC", size = 0.5) +
      geom_area(data = plt_this,
                aes(x = CMC, y = smth_dhs_norm),
                color = NA, alpha = 0.5, fill = "#6AA1FD") +
      # geom_path(data = plt_df,
      #           aes(x = CMC, y = smth_dhs_norm),
      #           color = "black") +
      geom_point(data = plt_this,
                 aes(x = CMC, y = MDC_norm),
                 color = "black", size = 2) +
      #scale_y_continuous(limits = c(0, 0.8)) +
      scale_x_continuous(breaks = xlbs_ids, labels = xlbs) +
      ylab("Density") + 
      xlab("Year") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
    
    # Plot facets
    if (MDC_kde_national) {
      facet_nets <- net_plt + facet_wrap(~countryname,  ncol=1, labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free_y", strip.position = "right")
    } else {
      if (urban_split_MDC) {
        facet_nets <- net_plt + facet_grid(ADM1 ~ urbanicity, labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free_y") + ggtitle(ttl_nm)# + geom_blank(data = areas_plt_df, aes(y = ylim))
      } else {
        facet_nets <- net_plt + facet_wrap(~ADM1,  ncol=1, labeller = label_wrap_gen(width = 2, multi_line = TRUE), scales = "free_y", strip.position = "right") + theme(strip.text.y = element_text(size = 7)) + ggtitle(ttl_nm)
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


plot_MDCs <- function(dataset, net_density_name, ref_density_name) {
  
  #Bound time series by first and last net recorded
  if (urban_split_MDC) {
    dataset_filter <- dataset
  } else {
    dataset_filter <- dataset[which(dataset$urbanicity=="urban"),]
  }
  for (n in 1:N_areas) {
    plt_ids <- which(!(dataset_filter$area_id == areas_df$area_ID[n] &
                         ((dataset_filter$CMC < areas_df$min_net_age_rec[n]) |
                            (dataset_filter$CMC > areas_df$max_net_age_rec[n]))))
    dataset_filter <- dataset_filter[plt_ids,]
  }
  
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
                       folderpath,
                       fileprefix,
                       timestamp,
                       multi_plt = FALSE)
    dev.off()
    #multiple pdf plots
    generate_MDC_plots(dataset_filter,
                       folderpath,
                       fileprefix,
                       timestamp,
                       multi_plt = TRUE)
  }
  
  return(NULL)
}



# 
# 
# 
# 
# all_country_df <- NULL
# 
# for (i in 1:N_ISO2) {
#   all_country_df <- rbind.data.frame(all_country_df,dataset_filter[1:length(CMC_series),])
# }
# all_country_df$missing_dhs <- all_country_df$scaled_camp_nets * norm_fac
# all_country_df$missing_dhs[which(all_country_df$missing_dhs < y_lim)] <- NA
# all_country_df$missing_dhs[which(all_country_df$missing_dhs > y_lim)] <- y_lim
# 
# all_country_df$scaled_camp_capped <- all_country_df$scaled_camp_nets * norm_fac
# all_country_df$scaled_camp_capped[which(all_country_df$scaled_camp_capped > y_lim)] <- y_lim
# 
# all_country_df$smth_dhs_capped <- all_country_df$smth_dhs * norm_fac
# all_country_df$smth_dhs_capped[which(all_country_df$smth_dhs_capped > y_lim)] <- y_lim
# 
# net_plt <- ggplot() +
#   # geom_step(data = all_country_df,
#   #           aes(x = CMC, y = scaled_camp_capped),
#   #           alpha = 1, color = "#458AFC", size = 0.5) +
#   geom_step(data = all_country_df,
#             aes(x = CMC, y = scaled_national * norm_fac),
#             alpha = 1, color = "#F66C19", size = 0.5) +
#   # geom_area(data = all_country_df,
#   #           aes(x = CMC, y = smth_dhs_capped),
#   #           color = NA, alpha = 0.5, fill = "#6AA1FD") +
#   geom_area(data = all_country_df,
#             aes(x = CMC, y = smth_dist * norm_fac),
#             color = NA, alpha = 0.5, fill = "#F88947") +
#   geom_path(data = all_country_df,
#             aes(x = CMC, y = comb_net_series * norm_fac),
#             color = "#30123B") +
#   geom_point(data = all_country_df,
#              aes(x = CMC, y = MDC_comb_series * norm_fac),
#              color = "#30123B", size = 2) +
#   scale_y_continuous(limits = c(0, y_lim)) +
#   scale_x_continuous(breaks = ylbs_ids, labels = ylbs) +
#   ylab("Density") + 
#   xlab("Year") +
#   theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
# 
# net_plt + facet_wrap(~ISO2, ncol = 1)
# 
# facet_nets <- net_plt + facet_grid(ISO2, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + ggtitle(ttl_nm)
