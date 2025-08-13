library(cowplot)

#eir_plotting.R

plot_eir_time_cow <- function(ISO2 = "SN",
                              plt_map_prop_step = FALSE,
                              plt_map_prop_point = FALSE,
                              plt_dhs_prop = FALSE,
                              plt_map_eir = FALSE,
                              plt_dhs_eir = FALSE,
                              plt_usage_prop = FALSE,
                              plt_theoretical_max_usage = FALSE,
                              plt_theoretical_min_usage = FALSE,
                              plt_usage = FALSE,
                              plt_usage_p0 = FALSE,
                              plt_usage_d = FALSE,
                              plt_usage_bb = FALSE,
                              plt_access_prop = FALSE,
                              plt_access = FALSE,
                              plt_access_d = FALSE,
                              plt_access_bb = FALSE,
                              plt_uga_prop = FALSE,
                              plt_uga = FALSE,
                              CrI50_alpha = 0.2,
                              CrI95_alpha = 0.2,
                              sf = 100,
                              ylab_str1 = NULL,
                              ylab_str2 = NULL,
                              save_str = NULL,
                              ymax2 = 1,
                              ybreaks2 = 0.2,
                              yr0 = 2005,
                              mn0 = 1,
                              yr1 = 2024,
                              mn1 = 12,
                              urban_rural = "rural") {
  
  # Usage or access timeseries plot
  
  ctry_data <- net_data[which(net_data$ISO2 == ISO2),]
  ctry_data$uga <- ctry_data$used / ctry_data$access
  ctry_data$uga[which(is.nan(ctry_data$uga))] <- NA
  ctry_data$prop_used[which(is.nan(ctry_data$prop_used))] <- NA
  ctry_data$prop_access <- ctry_data$access / ctry_data$total
  ctry_data$prop_access[which(is.nan(ctry_data$prop_access))] <- NA
  
  if (plt_uga) {
    ctry_data$uga[ctry_data$uga > 1] <- 1
    ctry_data$P_condu_mean[ctry_data$P_condu_mean > 1] <- 1
    ctry_data$P_condu_LB1[ctry_data$P_condu_LB1 > 1] <- 1
    ctry_data$P_condu_UB1[ctry_data$P_condu_UB1 > 1] <- 1
  }
  
  
  if (plt_theoretical_max_usage) {
    ctry_data$theoret_max_u <- (ctry_data$P0_u_mean-ctry_data$D_u_mean) * exp(
      - ctry_data$months_post_mdc / ctry_data$invlam_u_mean
    ) + ctry_data$D_u_mean
  }
  
  time_series <- CMC_to_date(CMC_series)
  jan_ids <- which(time_series[,2] == 1)
  xvals <- CMC_series[jan_ids]
  xlbs <- paste0("01/",substr(time_series[jan_ids,1],3,4))
  
  ctry_data %<>% filter(urbanicity == urban_rural)
  
  plt1 <- ggplot(data = ctry_data,
                 aes(x = CMC, fill = ADM1))
  if(plt_usage_prop) {
    plt1 <- plt1 +
      geom_count(data = subset(ctry_data, !is.na(prop_used)),
                 aes(y = prop_used * sf, size = total, color = ADM1, shape = urbanicity),
                 alpha = 0.5)
  }
  if(plt_usage_bb) {
    plt1 <- plt1 +
      geom_ribbon(aes(ymin = Pbb_u_LB1 * sf, ymax = Pbb_u_UB1 * sf), alpha = CrI50_alpha) +
      geom_ribbon(aes(ymin = Pbb_u_LB2 * sf, ymax = Pbb_u_UB2 * sf), alpha = CrI50_alpha) +
      geom_ribbon(aes(ymin = Pbb_u_LB3 * sf, ymax = Pbb_u_UB3 * sf), alpha = CrI95_alpha) +
      geom_path(aes(y = Pbb_u_mean * sf, colour = ADM1))
  }
  if(plt_usage) {
    plt1 <- plt1 +
      geom_ribbon(aes(ymin = P_u_LB1*  sf, ymax = P_u_UB1 * sf), alpha = 0.4) +
      geom_path(aes(y = P_u_mean * sf, colour = ADM1), alpha = 1)
  }
  if (plt_theoretical_max_usage) {
    plt1 <- plt1 +
      geom_path(aes(y = theoret_max_u * sf, colour = ADM1), alpha = 0.65)#, linetype = "dotted")
  }
  if(plt_usage_p0) {
    plt1 <- plt1 +
      geom_path(aes(y = P0_u_mean * sf, colour = ADM1), linetype = "dotted")
  }
  if(plt_usage_d) {
    plt1 <- plt1 +
      geom_path(aes(y = D_u_mean * sf, colour = ADM1), linetype = "dashed")
  }
  if(plt_access_prop) {
    plt1 <- plt1 +
      geom_count(data = subset(ctry_data, !is.na(prop_used)),
                 aes(y = prop_used * sf, size = total, color = ADM1, shape = urbanicity),
                 alpha = 0.5)
  }
  if(plt_access_bb) {
    plt1 <- plt1 +
      geom_ribbon(aes(ymin = Pbb_a_LB1 * sf, ymax = Pbb_a_UB1 * sf), alpha = CrI50_alpha) +
      geom_ribbon(aes(ymin = Pbb_a_LB2 * sf, ymax = Pbb_a_UB2 * sf), alpha = CrI50_alpha) +
      geom_ribbon(aes(ymin = Pbb_a_LB3 * sf, ymax = Pbb_a_UB3 * sf), alpha = CrI95_alpha) +
      geom_path(aes(y = Pbb_a_mean * sf, colour = ADM1))
  }
  if(plt_access) {
    plt1 <- plt1 +
      #geom_ribbon(aes(ymin = P_a_LB1*  sf, ymax = P_a_UB1 * sf), alpha = 0.4) +
      geom_path(aes(y = P_a_mean * sf, colour = ADM1), linetype = "dotted")
  }
  if(plt_access_d) {
    plt1 <- plt1 +
      geom_path(aes(y = D_a_mean * sf, colour = ADM1), linetype = "dotted")
  }
  if(plt_uga_prop) {
    plt1 <- plt1 +
      geom_count(data = subset(ctry_data, !is.na(uga)),
                 aes(y = uga * sf, size = total, color = ADM1, shape = urbanicity),
                 alpha = 0.5)
  }
  if(plt_uga) {
    plt1 <- plt1 +
      geom_ribbon(aes(ymin = P_condu_LB1 * sf, ymax = P_condu_UB1 * sf), alpha = 0.4) +
      geom_path(aes(y = P_condu_mean * sf, colour = ADM1))
  }
  
  plt1 <- plt1 +
    ylab(ylab_str1) + 
    xlab("Year") +
    scale_y_continuous(breaks = seq(0,1,0.2)*100,limits = c(0, 1)*100) +
    scale_x_continuous(
      breaks = xvals, labels = xlbs, limits = c(
        date_to_CMC(yr0,mn0), date_to_CMC(yr1,mn1)
      )
    ) +    theme_bw() +
    theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
          axis.title.x = element_blank(),
          strip.text.y = element_text(angle = 0)) +
    theme(
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
  
  #plt1 + facet_wrap(~ADM1, nrow = 2)
  plt1 <- plt1 +
    # facet_wrap(~ADM1, ncol = 1) +
    # theme(strip.background = element_blank(),
    #       strip.text.y = element_blank()) +
    facet_grid(ADM1~urbanicity) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          strip.text.y = element_blank())
  
  # EIR plot
  
  ctry_data <- eir_net_data[which(eir_net_data$ISO2 == ISO2),]
  
  ctry_data$month <- CMC_to_date(ctry_data$CMC)[2]
  
  time_series <- CMC_to_date(CMC_series)
  jan_ids <- which(time_series[,2] == 1)
  xvals <- CMC_series[jan_ids]
  xlbs <- paste0("01/",substr(time_series[jan_ids,1],3,4))
  
  shp_size <- 1
  
  ctry_data$pfpr_2_10_fit[ctry_data$pfpr_2_10_fit > ymax2] <- ymax2
  ctry_data$pfpr_6_59_mo_fit[ctry_data$pfpr_6_59_mo_fit > ymax2] <- ymax2
  ctry_data$pfpr_dhs[ctry_data$pfpr_dhs > ymax2] <- ymax2
  ctry_data$pfpr_dhs_lo[ctry_data$pfpr_dhs_lo > ymax2] <- ymax2
  ctry_data$pfpr_dhs_hi[ctry_data$pfpr_dhs_hi > ymax2] <- ymax2
  ctry_data$pfpr_map[ctry_data$pfpr_map > ymax2] <- ymax2
  
  ctry_data %<>% filter(urbanicity == urban_rural)
  
  plt2 <- ggplot(data = ctry_data,
                 aes(x = CMC))
  if(plt_map_prop_step) {
    plt2 <- plt2 +
      geom_step(
        aes(y = pfpr_map * sf, colour = ADM1),
        position = position_nudge(x = -0.5),
        linetype = "dotted"
      )
  }
  if(plt_map_prop_point) {
    plt2 <- plt2 +
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
    plt2 <- plt2 +
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
  }
  if(plt_map_eir) {
    plt2 <- plt2 +
      geom_path(aes(y = pfpr_2_10_fit * sf, colour = ADM1),
                alpha = 1)
  }
  if(plt_dhs_eir) {
    plt2 <- plt2 +
      geom_path(aes(y = pfpr_6_59_mo_fit * sf, colour = ADM1),
                alpha = 0.5)
  }
  
  plt2 <- plt2 +
    
    ylab(ylab_str2) + 
    xlab("Year") +
    #scale_y_continuous(trans = 'log10') +
    scale_y_continuous(breaks = seq(0,ymax2,ybreaks2)*100,limits = c(0, ymax2)*100) +
    scale_x_continuous(
      breaks = xvals, labels = xlbs, limits = c(
        date_to_CMC(yr0,mn0), date_to_CMC(yr1,mn1)
      )
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
          axis.title.x = element_blank(),
          strip.text.y = element_text(angle = 0)) +
    theme(
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
  
  plt2 <- plt2 +
    facet_grid(ADM1~urbanicity) +
    theme(strip.background.x = element_blank(),
          strip.text.x = element_blank())
  
  
  plot_grid(plt1, plt2, labels = c('A', 'B'), rel_widths = c(1, 1.24))
  
  ggsave(paste0(ISO2, "_", save_str, ".pdf"), bg = "transparent",
         w = 8, h = 10, dpi = 450)
  
  #print(tplt)
  
}





plot_eir_time_cow(ISO2 = "SN",
                  ymax2 = 0.3,
                  ybreaks2 = 0.1,
                  yr0 = 2008,
                  mn0 = 1,
                  yr1 = 2024,
                  mn1 = 12,
                  plt_map_prop_step = FALSE,
                  plt_map_prop_point = TRUE,
                  plt_dhs_prop = TRUE,
                  plt_map_eir = TRUE,
                  plt_dhs_eir = TRUE,
                  ylab_str2 = expression(italic(Pf)*PR~"(%)"),
                  plt_usage_prop = TRUE,
                  plt_usage_bb = TRUE,
                  plt_usage_d = TRUE,
                  plt_access = TRUE,
                  ylab_str1 = "Probability (%)",
                  save_str = "full_use_pfpr")
