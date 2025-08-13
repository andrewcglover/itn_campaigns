#eir_plotting.R

plot_eir_series <- function(ISO2 = "SN",
                            plt_map_prop_step = FALSE,
                            plt_map_prop_point = FALSE,
                            plt_dhs_prop = FALSE,
                            plt_map_eir = FALSE,
                            plt_dhs_eir = FALSE,
                            sf = 100,
                            ylab_str = NULL,
                            save_str = NULL,
                            ymax = 1) {
  
  ctry_data <- eir_net_data[which(eir_net_data$ISO2 == ISO2),]
  
  ctry_data$month <- CMC_to_date(ctry_data$CMC)[2]
  
  # if(plt_map_prop_point) {
  #   ctry_data$month <- CMC_to_date(ctry_data$CMC)[2]
  #   ctry_data$pfpr_map_pts <- ctry_data$pfpr_map
  #   ctry_data$pfpr_map_pts[ctry_data$month != 6] <- NA
  # }
  
  # ctry_data$uga <- ctry_data$used / ctry_data$access
  # ctry_data$uga[which(is.nan(ctry_data$uga))] <- NA
  # ctry_data$prop_used[which(is.nan(ctry_data$prop_used))] <- NA
  # ctry_data$prop_access <- ctry_data$access / ctry_data$total
  # ctry_data$prop_access[which(is.nan(ctry_data$prop_access))] <- NA
  
  time_series <- CMC_to_date(CMC_series)
  jan_ids <- which(time_series[,2] == 1)
  xvals <- CMC_series[jan_ids]
  xlbs <- paste0("01/",substr(time_series[jan_ids,1],3,4))
  
  shp_size <- 1
  
  ctry_data$pfpr_2_10_fit[ctry_data$pfpr_2_10_fit > ymax] <- ymax
  ctry_data$pfpr_6_59_mo_fit[ctry_data$pfpr_6_59_mo_fit > ymax] <- ymax
  ctry_data$pfpr_dhs[ctry_data$pfpr_dhs > ymax] <- ymax
  ctry_data$pfpr_dhs_lo[ctry_data$pfpr_dhs_lo > ymax] <- ymax
  ctry_data$pfpr_dhs_hi[ctry_data$pfpr_dhs_hi > ymax] <- ymax
  ctry_data$pfpr_map[ctry_data$pfpr_map > ymax] <- ymax
  
  tplt <- ggplot(data = ctry_data,
                 aes(x = CMC))
  if(plt_map_prop_step) {
    tplt <- tplt +
      geom_step(
        aes(y = pfpr_map * sf, colour = ADM1),
        position = position_nudge(x = -0.5),
        linetype = "dotted"
      )
  }
  if(plt_map_prop_point) {
    tplt <- tplt +
      geom_point(data = subset(ctry_data, month == 7 & urbanicity == "rural"),
                 aes(y = pfpr_map * sf, colour = ADM1),
                 #position = position_nudge(x = 1),
                 shape = 1,
                 size = shp_size
      ) +
      geom_point(data = subset(ctry_data, month == 7 & urbanicity == "urban"),
                 aes(y = pfpr_map * sf, colour = ADM1),
                 #position = position_nudge(x = 1),
                 shape = 2,
                 size = shp_size
      )
  }
  if(plt_dhs_prop) {
    tplt <- tplt +
      geom_linerange(data = subset(ctry_data, !is.na(pfpr_dhs)),
                     aes(y = pfpr_dhs * sf, ymin = pfpr_dhs_lo * sf,
                         ymax = pfpr_dhs_hi * sf,
                         color = ADM1),
                     alpha = 0.5,
                     lwd = 0.35) +
      geom_point(data = subset(ctry_data, !is.na(pfpr_dhs)),
                 aes(y = pfpr_dhs * sf,
                     color = ADM1, shape = urbanicity),
                 alpha = 0.5,
                 size = shp_size)
    
    # geom_point(data = subset(ctry_data, !is.na(pfpr_dhs)),
    #            aes(y = prop_used * sf, color = ADM1, shape = urbanicity),
    #            alpha = 0.5) +
    # geom_errorbar(data = subset(ctry_data, !is.na(pfpr_dhs)),
    #               aes(y = prop_used * sf, color = ADM1, shape = urbanicity),
    #               alpha = 0.5)
  }
  if(plt_map_eir) {
    tplt <- tplt +
      geom_path(aes(y = pfpr_2_10_fit * sf, colour = ADM1),
                alpha = 1)
  }
  if(plt_dhs_eir) {
    tplt <- tplt +
      geom_path(aes(y = pfpr_6_59_mo_fit * sf, colour = ADM1),
                alpha = 0.5)
  }
  tplt <- tplt +
    # scale_fill_viridis(
    #   discrete = TRUE,
    #   option = "H",
    #   guide = guide_legend(title = "Region")#,
    #   #labels = "3 year interval"#only_label_vals
    # ) +
    # scale_color_viridis(
    #   discrete = TRUE,
    #   option = "H",
    #   guide = guide_legend(title = "People\nsurveyed")#,
    #   #labels = "3 year interval"#only_label_vals
  # ) +
  ylab(ylab_str) + 
    xlab("Year") +
    #scale_y_continuous(trans = 'log10') +
    scale_y_continuous(breaks = seq(0,1,0.2)*100,limits = c(0, ymax)*100) +
    #scale_y_continuous(breaks = seq(0,1,0.2)*100,limits = c(0, 1)*100) +
    #scale_x_continuous(breaks = ylbs_ids, labels = ylbs, limits = c(date_to_CMC(2014,1), date_to_CMC(2021,1))) +
    #scale_x_continuous(breaks = ylbs_ids, labels = ylbs, limits = c(date_to_CMC(2013,6), date_to_CMC(2021,6))) +
    #scale_x_continuous(breaks = ylbs_ids, labels = ylbs, limits = c(date_to_CMC(2010,6), date_to_CMC(2020,6))) +
    scale_x_continuous(breaks = xvals, labels = xlbs) +#, limits = c(date_to_CMC(2010,1), date_to_CMC(2021,6))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
          axis.title.x = element_blank(),
          strip.text.y = element_text(angle = 0)) +
    theme(
      #panel.background = element_rect(fill = "transparent",
      #                                colour = NA_character_), # necessary to avoid drawing panel outline
      #panel.grid.major = element_blank(), # get rid of major grid
      #panel.grid.minor = element_blank(), # get rid of minor grid
      plot.background = element_rect(fill = "transparent",
                                     colour = NA_character_), # necessary to avoid drawing plot outline
      legend.background = element_rect(fill = "transparent"),
      legend.box.background = element_rect(fill = "transparent"),
      legend.key = element_rect(fill = "transparent")
    ) +
    guides(colour = "none",
           fill = "none",
           shape = "none",
           size = "none")
  
  
  #tplt + facet_wrap(~ADM1, nrow = 2)
  tplt + facet_grid(ADM1~urbanicity)
  #tplt <- tplt + facet_grid(ADM1~urbanicity)
  
  #tplt <- tplt + scale_y_cut(breaks=c(10), which=c(1, 2), scales=c(3, 0.5))
  #tplt + scale_y_cut(breaks=c(10), which=c(1, 2), scales=c(3, 0.5))
  
  
  ggsave(paste0(ISO2, "_", save_str, ".pdf"), bg = "transparent",
         w = 8, h = 10, dpi = 450)
  
  #print(tplt)
  
}



plot_eir_series(ISO2 = "ML",
                ymax = 1,
                plt_map_prop_step = FALSE,
                plt_map_prop_point = TRUE,
                plt_dhs_prop = TRUE,
                plt_map_eir = TRUE,
                plt_dhs_eir = TRUE,
                ylab_str = expression(italic(Pf)*Pr),
                save_str = "eir_MAR25")