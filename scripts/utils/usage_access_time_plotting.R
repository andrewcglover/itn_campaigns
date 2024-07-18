#usage_access_time_plotting.R

plot_timeseries <- function(ISO2 = "SN",
                            plt_usage_prop = FALSE,
                            plt_usage = FALSE,
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
                            sf = 100) {
  
  ctry_data <- net_data[which(net_data$ISO2 == ISO2),]
  ctry_data$uga <- ctry_data$used / ctry_data$access
  ctry_data$uga[which(is.nan(ctry_data$uga))] <- NA
  ctry_data$prop_used[which(is.nan(ctry_data$prop_used))] <- NA
  ctry_data$prop_access <- ctry_data$access / ctry_data$total
  ctry_data$prop_access[which(is.nan(ctry_data$prop_access))] <- NA
  
  ctry_data$uga[ctry_data$uga > 1] <- 1
  ctry_data$P_condu_mean[ctry_data$P_condu_mean > 1] <- 1
  ctry_data$P_condu_LB1[ctry_data$P_condu_LB1 > 1] <- 1
  ctry_data$P_condu_UB1[ctry_data$P_condu_UB1 > 1] <- 1
  # ctry_data$Pbb_condu_mean[ctry_data$Pbb_condu_mean > 1] <- 1
  # ctry_data$Pbb_condu_LB1[ctry_data$Pbb_condu_LB1 > 1] <- 1
  # ctry_data$Pbb_condu_UB1[ctry_data$Pbb_condu_UB1 > 1] <- 1
  # ctry_data$Pbb_condu_LB2[ctry_data$Pbb_condu_LB2 > 1] <- 1
  # ctry_data$Pbb_condu_UB2[ctry_data$Pbb_condu_UB2 > 1] <- 1
  # ctry_data$Pbb_condu_LB3[ctry_data$Pbb_condu_LB3 > 1] <- 1
  # ctry_data$Pbb_condu_UB3[ctry_data$Pbb_condu_UB3 > 1] <- 1
  
  time_series <- CMC_to_date(CMC_series)
  jan_ids <- which(time_series[,2] == 1)
  xvals <- CMC_series[jan_ids]
  xlbs <- paste0("01/",substr(time_series[jan_ids,1],3,4))
  
  
  tplt <- ggplot(data = ctry_data,
                 aes(x = CMC, fill = ADM1))
  if(plt_usage_prop) {
    tplt <- tplt +
      geom_count(data = subset(ctry_data, !is.na(prop_used)),
                 aes(y = prop_used * sf, size = total, color = ADM1, shape = urbanicity),
                 alpha = 0.5)
  }
  if(plt_usage_bb) {
    tplt <- tplt +
      geom_ribbon(aes(ymin = Pbb_u_LB1 * sf, ymax = Pbb_u_UB1 * sf), alpha = CrI50_alpha) +
      geom_ribbon(aes(ymin = Pbb_u_LB2 * sf, ymax = Pbb_u_UB2 * sf), alpha = CrI50_alpha) +
      geom_ribbon(aes(ymin = Pbb_u_LB3 * sf, ymax = Pbb_u_UB3 * sf), alpha = CrI95_alpha) +
      geom_path(aes(y = Pbb_u_mean * sf, colour = ADM1))
      #geom_path(aes(y = Pbb_u_LB1, colour = ADM1), linetype = "dotted") +
      #geom_path(aes(y = Pbb_u_UB1, colour = ADM1), linetype = "dotted")
  }
  if(plt_usage) {
    tplt <- tplt +
      geom_ribbon(aes(ymin = P_u_LB1*  sf, ymax = P_u_UB1 * sf), alpha = 0.4) +
      geom_path(aes(y = P_u_mean * sf, colour = ADM1))
  }
  if(plt_usage_d) {
    tplt <- tplt +
      #geom_ribbon(aes(ymin = P_u_LB1, ymax = P_u_UB1), alpha = CrI_alpha) +
      geom_path(aes(y = D_u_mean * sf, colour = ADM1), linetype = "dashed")
  }
  if(plt_access_prop) {
    tplt <- tplt +
      geom_count(data = subset(ctry_data, !is.na(prop_used)),
                 aes(y = prop_used * sf, size = total, color = ADM1, shape = urbanicity),
                 alpha = 0.5)
  }
  if(plt_access_bb) {
    tplt <- tplt +
      geom_ribbon(aes(ymin = Pbb_a_LB1 * sf, ymax = Pbb_a_UB1 * sf), alpha = CrI50_alpha) +
      geom_ribbon(aes(ymin = Pbb_a_LB2 * sf, ymax = Pbb_a_UB2 * sf), alpha = CrI50_alpha) +
      geom_ribbon(aes(ymin = Pbb_a_LB3 * sf, ymax = Pbb_a_UB3 * sf), alpha = CrI95_alpha) +
      geom_path(aes(y = Pbb_a_mean * sf, colour = ADM1))
  }
  if(plt_access) {
    tplt <- tplt +
      geom_ribbon(aes(ymin = P_a_LB1*  sf, ymax = P_a_UB1 * sf), alpha = 0.4) +
      geom_path(aes(y = P_a_mean * sf, colour = ADM1))
  }
  if(plt_access_d) {
    tplt <- tplt +
      geom_path(aes(y = D_a_mean * sf, colour = ADM1), linetype = "dashed")
  }
  if(plt_uga_prop) {
    tplt <- tplt +
      geom_count(data = subset(ctry_data, !is.na(uga)),
                 aes(y = uga * sf, size = total, color = ADM1, shape = urbanicity),
                 alpha = 0.5)
  }
  if(plt_uga) {
    tplt <- tplt +
      geom_ribbon(aes(ymin = P_condu_LB1 * sf, ymax = P_condu_UB1 * sf), alpha = 0.4) +
      #geom_ribbon(aes(ymin = Pbb_condu_LB1 * sf, ymax = Pbb_condu_UB1 * sf), alpha = 0.2) +
      #geom_ribbon(aes(ymin = Pbb_condu_LB2 * sf, ymax = Pbb_condu_UB2 * sf), alpha = 0.2) +
      #geom_ribbon(aes(ymin = Pbb_condu_LB3 * sf, ymax = Pbb_condu_UB3 * sf), alpha = 0.2) +
      geom_path(aes(y = P_condu_mean * sf, colour = ADM1))
  }
    # geom_ribbon(data = country_df,
    #             aes(x = CMC, ymax = d_UB, ymin = d_LB, fill = urbanicity),
    #             alpha = 0.3) +
    # geom_path(data = country_df,
    #           aes(x = CMC, y = mean_d, color = urbanicity),
    #           linetype = "dashed") +
    # geom_ribbon(data = country_df,
    #             aes(x = CMC, ymax = p_UB, ymin = p_LB, fill = urbanicity),
    #             alpha = 0.2) +
    # geom_path(data = country_df,
    #           aes(x = CMC, y = mean_p, color = urbanicity)) +
    # geom_count(data = campnets_df_plt,
    #            aes(x = CMC, y = prop_used, size = total, color = urbanicity),
    #            alpha = 0.2, stroke = NA) +
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
    ylab("Usage (%)") + 
    xlab("Year") +
    scale_y_continuous(breaks = seq(0,1,0.2)*100,limits = c(0, 1)*100) +
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
  
  ggsave(paste0(ISO2,"_usage_bb_alpha_final.pdf"), bg = "transparent",
         w = 8, h = 10, dpi = 450)
  
  #print(tplt)
  
}