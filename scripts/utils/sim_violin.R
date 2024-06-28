# net_cases_plot_test.R

sim_violin_plot <- function(sim_data,
                         fs_areas_included = NULL,
                         plotting_var = NULL,
                         costed = FALSE,
                         costed_and_uncosted = FALSE,
                         baseline_dist = "pyrethroid-only routine baseline") {
  
  # if (plot_cases_averted & !(plotting_var == "var_empty")) {
  #   print("Warning: cases averted plot priotised over entered variable")
  # }
  # 
  # if (plotting_var == "cases_averted") {
  #   
  # }
  
  for (i in 1:length(fs_areas_included)) {
    
    single_adm <- sim_data[which(sim_data$fs_area == fs_areas_included[i]),]
    
    # if (plotting_var == "cases_averted") {
    #   single_adm$baseline_infections <- rep(baseline_infecions,
    #                                         length.out = dim(single_adm)[1])
    #   single_adm$cases_averted <- single_adm$baseline_infections - single_adm$annual_infections
    # }
    
    if (plotting_var == "cases_averted") {
      #single_adm$var_den <- single_adm$baseline_cases - single_adm$pred_ann_infect
      single_adm$var_den <- single_adm$cases_averted
      yax_val <- "Annual cases averted"
      ylim_vals <- c(-0e5, 1.6e6)
    } else if (plotting_var == "cases_averted_per_USD") {
      single_adm$var_den <- single_adm$cases_averted_per_USD
      yax_val <- "Annual cases averted per USD"
      ylim_vals <- c(0, 1.5)
    } else if (plotting_var == "cases_averted_pp") {
      single_adm$var_den <- 10000 * (single_adm$baseline_cases - single_adm$pred_ann_infect) / single_adm$pop
      yax_val <- "Annual cases averted per 10,000"
      #ylim_vals <- c(-500, 1500)
      ylim_vals <- c(0, 10000)
      #ylim_vals <- c(-2500, 17500)
    } else if (plotting_var == "avg_pfpr") {
      single_adm$var_den <- single_adm$avg_pfpr * 10000
      yax_val <- "2-10 y/o prevalence per 10,000"
      ylim_vals <- c(0, 1300)
    } else if (plotting_var == "ann_incidence") {
      single_adm$var_den <- single_adm$ann_incidence * 10000
      yax_val <- "Incidence per 10,000 for all ages"
      ylim_vals <- c(0, 1000)
    } else if (plotting_var == "add_cases_averted") {
      single_adm$var_den <- single_adm$add_cases_averted
      yax_val <- "Additional annual cases averted"
      ylim_vals <- c(-0e5, 1.6e6)
    }
    
    
    sim_sum <- single_adm %>%
      group_by(net_strategy) %>%
      dplyr::summarise(var_mid = mean(var_den),
                       var_lo = quantile(var_den, 0.025),
                       var_hi = quantile(var_den, 0.975),
                       var_max = max(var_den))
    
    cost_sum <- single_adm %>%
      group_by(net_strategy) %>%
      dplyr::summarise(mean_add_cost = mean(add_cost))
    sim_sum$add_mil_cost <- paste0("$",round(cost_sum$mean_add_cost/1e6,2),"M")
    
    only_label_vals <- c("2 year interval", "3 year interval", "Continous only")
    label_vals <- c("2 year interval", "3 year interval", "Continuous only")
    
    if (costed) {costed_str <- "_costed"} else {costed_str <- "_uncosted"}
    if (costed_and_uncosted) {costed_str <- "_costed_and_uncosted"}
    if (costed_and_uncosted) {label_vals <- c("2 year interval",
                                              "2 year interval (costed)",
                                              "3 year interval",
                                              "3 year interval (costed)")}
    if (costed_and_uncosted) {only_label_vals <- c("2 year interval",
                                              "2 year interval (costed)",
                                              "3 year interval",
                                              "Routine only")}
    
    filename <- paste0(plotting_var, costed_str, i, ".png")
    
    single_adm_no_base <- single_adm %>%
      #filter(net_strategy != baseline_dist) %>%
      filter(net_strategy != "no future nets")
    sim_sum_no_base <- sim_sum %>%
      #filter(net_strategy != baseline_dist) %>%
      filter(net_strategy != "no future nets")
    
    # routine_df <- single_adm %>% filter(net_strategy == baseline_dist)
    # routine_sumdf <- sim_sum %>% filter(net_strategy == baseline_dist)
    # routine_df$net_strategy <- rep("a pyrethroid-only rouine")
    # routine_sumdf$net_strategy <- rep("a pyrethroid-only rouine")
    
    
    only_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-only")
    pbo_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-PBO")
    #pbo_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-PBO")
    
    pyrrole_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-pyrrole")
    
    only_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-only", net_strategy))
    pbo_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-PBO", net_strategy))
    pyrrole_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-pyrrole", net_strategy))
    
    cost_text_size <- 4
    
    #single_adm_no_base %>%
      ggplot() +
      #   geom_violin(data = routine_df,
      #               aes(x = net_strategy,
      #                   y = var_den,
      #                   fill = net_strategy),
      #               alpha = 0.4,
      #               color = NA) +
      #   geom_errorbar(data = routine_sumdf,
      #                 aes(x = net_strategy,
      #                     ymin = var_lo,
      #                     ymax = var_hi,
      #                     color = net_strategy),
      #                 alpha = 0.8) +
      #   geom_point(data = routine_sumdf,
      #              aes(x = net_strategy,
      #                  y = var_mid,
      #                  color = net_strategy),
      #              alpha = 1) +
      #   # geom_text(data = only_sumdf %>% filter (net_strategy == baseline_dist),
      #   #           aes(x = net_strategy,
      #   #               y = var_max,
      #   #               label = add_mil_cost),
      #   #           vjust = -2,
      #   #           size = cost_text_size) +
      #   #   geom_text(data = only_sumdf %>% filter (net_strategy != baseline_dist),
      #   #             aes(x = net_strategy,
      #   #                 y = var_max,
      #   #                 label = add_mil_cost),
      #   #             vjust = -2,
      # #             size = cost_text_size) +
      # scale_fill_viridis(
      #   alpha = 1,
      #   begin = 0,
      #   end = 0,
      #   direction = 1,
      #   discrete = TRUE,
      #   option = "A",
      #   guide = guide_legend(title = "Pyrethroid-only"),
      #   labels = "Routine only"
      # ) +
      #   scale_color_viridis(
      #     alpha = 1,
      #     begin = 0,
      #     end = 0,
      #     direction = 1,
      #     discrete = TRUE,
      #     option = "A",
      #     guide = guide_legend(title = "Pyrethroid-only"),
      #     labels = "Routine only"
      #   ) +
      #   new_scale_fill() +
      #   new_scale_colour() +
        geom_violin(data = only_df,
                    aes(x = net_strategy,
                        y = var_den,
                        fill = net_strategy),
                    alpha = 0.4,
                    color = NA) +
        geom_errorbar(data = only_sumdf,
                      aes(x = net_strategy,
                          ymin = var_lo,
                          ymax = var_hi,
                          color = net_strategy),
                      alpha = 0.8) +
        geom_point(data = only_sumdf,
                   aes(x = net_strategy,
                       y = var_mid,
                       color = net_strategy),
                   alpha = 1) +
        geom_text(data = only_sumdf,
                  aes(x = net_strategy,
                      y = var_max,
                      label = add_mil_cost),
                  vjust = -1,
                  size = cost_text_size) +
        # geom_text(data = only_sumdf %>% filter (net_strategy == baseline_dist),
        #           aes(x = net_strategy,
        #               y = var_max,
        #               label = add_mil_cost),
        #           vjust = -2,
        #           size = cost_text_size) +
        #   geom_text(data = only_sumdf %>% filter (net_strategy != baseline_dist),
        #             aes(x = net_strategy,
        #                 y = var_max,
        #                 label = add_mil_cost),
        #             vjust = -2,
        #             size = cost_text_size) +
      scale_fill_viridis(
        alpha = 1,
        #begin = 0.8,
        #end = 0.95,
        begin = 0.875,
        end = 0.875,
        direction = -1,
        discrete = TRUE,
        option = "H",
        guide = guide_legend(title = "Pyrethroid-only"),
        labels = "3 year interval"#only_label_vals
        ) +
        scale_color_viridis(
          alpha = 1,
          #begin = 0.8,
          #end = 0.95,
          begin = 0.875,
          end = 0.875,
          direction = -1,
          discrete = TRUE,
          option = "H",
          guide = guide_legend(title = "Pyrethroid-only"),
          labels = "3 year interval"#only_label_vals
        ) +
        new_scale_fill() +
        new_scale_colour() +
        geom_violin(data = pbo_df,
                    aes(x = net_strategy,
                        y = var_den,
                        fill = net_strategy),
                    alpha = 0.4,
                    color = NA) +
        geom_errorbar(data = pbo_sumdf,
                      aes(x = net_strategy,
                          ymin = var_lo,
                          ymax = var_hi,
                          color = net_strategy),
                      alpha = 0.8) +
        geom_point(data = pbo_sumdf,
                   aes(x = net_strategy,
                       y = var_mid,
                       color = net_strategy),
                   alpha = 1) +
        geom_text(data = pbo_sumdf,
                  aes(x = net_strategy,
                      y = var_max,
                      label = add_mil_cost),
                  vjust = -1,
                  size = cost_text_size) +
      scale_fill_viridis(
        alpha = 1,
        begin = 0.45,
        end = 0.75,
        direction = 1,
        discrete = TRUE,
        option = "D",
        guide = guide_legend(title = "Pyrethroid-PBO"),
        labels = label_vals
      ) +
        scale_color_viridis(
          alpha = 1,
          begin = 0.45,
          end = 0.75,
          direction = 1,
          discrete = TRUE,
          option = "D",
          guide = guide_legend(title = "Pyrethroid-PBO"),
          labels = label_vals
        ) +
        new_scale_fill() +
        new_scale_colour() +
        geom_violin(data = pyrrole_df,
                    aes(x = net_strategy,
                        y = var_den,
                        fill = net_strategy),
                    alpha = 0.4,
                    color = NA) +
        geom_errorbar(data = pyrrole_sumdf,
                      aes(x = net_strategy,
                          ymin = var_lo,
                          ymax = var_hi,
                          color = net_strategy),
                      alpha = 0.8) +
        geom_point(data = pyrrole_sumdf,
                   aes(x = net_strategy,
                       y = var_mid,
                       color = net_strategy),
                   alpha = 1) +
        geom_text(data = pyrrole_sumdf,
                  aes(x = net_strategy,
                      y = var_max,
                      label = add_mil_cost),
                  vjust = -1,
                  size = cost_text_size) +
        scale_fill_viridis(
          alpha = 1,
          #begin = 0.05,
          begin = 0.2,
          end = 0.2,
          direction = 1,
          discrete = TRUE,
          option = "H",
          guide = guide_legend(title = "Pyrethroid-pyrrole"),
          labels = "2 year interval"#label_vals
        ) +
        scale_color_viridis(
          alpha = 1,
          #begin = 0.05,
          begin = 0.2,
          end = 0.2,
          direction = 1,
          discrete = TRUE,
          option = "H",
          guide = guide_legend(title = "Pyrethroid-pyrrole"),
          labels = "2 year interval"#label_vals
        ) +
        # guides(fill = guide_legend(title = "Net Strategy")) +
        # guides(colour = guide_legend(title = "Net Strategy")) +
        # labs(colour = "Pyrethroid-PBO",
        #      y = "add_cases_averted_per_usd") +
        ylab(yax_val) +
        xlab(NULL) +
        scale_y_continuous(limits = ylim_vals,
                           labels = label_comma()) + #,
                           #breaks = round(seq(ylim_vals[1], ylim_vals[2], by = 0.25),2)) +
        scale_x_discrete(breaks = NULL) +
        ggtitle(fs_areas_included[i])
    
      # ggsave(filename, bg = "white",
      #        w = 10, h = 4, dpi = 450)
      
      ggsave(filename, bg = "white",
             w = 4, h = 3, dpi = 450)
      
      
      #print to screen
      print(fs_areas_included[i])

      print(paste("only_var_lo =", only_sumdf$var_lo))
      print(paste("only_var_mid =", only_sumdf$var_mid))
      print(paste("only_var_hi =", only_sumdf$var_hi))
      # print(paste("only_pc_lo =", (only_sumdf$var_lo[1] - only_sumdf$var_lo[2])/only_sumdf$var_lo[2]))
      # print(paste("only_pc_mid =", (only_sumdf$var_mid[1] - only_sumdf$var_mid[2])/only_sumdf$var_mid[2]))
      # print(paste("only_pc_hi =", (only_sumdf$var_hi[1] - only_sumdf$var_hi[2])/only_sumdf$var_hi[2]))
      # 
      # print(paste("pbo_var_lo =", pbo_sumdf$var_lo))
      # print(paste("pbo_var_mid =", pbo_sumdf$var_mid))
      # print(paste("pbo_var_hi =", pbo_sumdf$var_hi))
      # print(paste("pbo_pc_lo =", (pbo_sumdf$var_lo[1] - pbo_sumdf$var_lo[2])/pbo_sumdf$var_lo[2]))
      # print(paste("pbo_pc_mid =", (pbo_sumdf$var_mid[1] - pbo_sumdf$var_mid[2])/pbo_sumdf$var_mid[2]))
      # print(paste("pbo_pc_hi =", (pbo_sumdf$var_hi[1] - pbo_sumdf$var_hi[2])/pbo_sumdf$var_hi[2]))
      # 
      print(paste("pyrrole_var_lo =", pyrrole_sumdf$var_lo))
      print(paste("pyrrole_var_mid =", pyrrole_sumdf$var_mid))
      print(paste("pyrrole_var_hi =", pyrrole_sumdf$var_hi))
      # print(paste("pyrrole_pc_lo =", (pyrrole_sumdf$var_lo[1] - pyrrole_sumdf$var_lo[2])/pyrrole_sumdf$var_lo[2]))
      # print(paste("pyrrole_pc_mid =", (pyrrole_sumdf$var_mid[1] - pyrrole_sumdf$var_mid[2])/pyrrole_sumdf$var_mid[2]))
      # print(paste("pyrrole_pc_hi =", (pyrrole_sumdf$var_hi[1] - pyrrole_sumdf$var_hi[2])/pyrrole_sumdf$var_hi[2]))
      # 
      
    # ggsave(filename, bg = "white",
    #        w = 7, h = 5, dpi = 450)
    
    
    
  }
  
}




sim_violin_plot_no_costings <- function(sim_data,
                            fs_areas_included = NULL,
                            plotting_var = NULL,
                            costed = FALSE,
                            costed_and_uncosted = FALSE,
                            baseline_dist = "pyrethroid-only routine baseline") {
    
    if (plotting_var == "cases_averted") {
      #single_adm$var_den <- single_adm$baseline_cases - single_adm$pred_ann_infect
      single_adm$var_den <- single_adm$cases_averted
      yax_val <- "Annual cases averted"
      ylim_vals <- c(-0e5, 1.6e6)
    } else if (plotting_var == "cases_averted_per_USD") {
      single_adm$var_den <- single_adm$cases_averted_per_USD
      yax_val <- "Annual cases averted per USD"
      ylim_vals <- c(0, 1.5)
    } else if (plotting_var == "cases_averted_pp") {
      single_adm$var_den <- 10000 * (single_adm$baseline_cases - single_adm$pred_ann_infect) / single_adm$pop
      yax_val <- "Annual cases averted per 10,000"
      #ylim_vals <- c(-500, 1500)
      ylim_vals <- c(0, 10000)
      #ylim_vals <- c(-2500, 17500)
    } else if (plotting_var == "avg_pfpr") {
      single_adm$var_den <- single_adm$avg_pfpr * 10000
      yax_val <- "2-10 y/o prevalence per 10,000"
      ylim_vals <- c(0, 1300)
    } else if (plotting_var == "ann_incidence") {
      single_adm$var_den <- single_adm$ann_incidence * 10000
      yax_val <- "Incidence per 10,000 for all ages"
      ylim_vals <- c(0, 1000)
    } else if (plotting_var == "add_cases_averted") {
      sim_data$var_den <- sim_data$add_cases_averted
      yax_val <- "Additional annual cases averted"
      ylim_vals <- c(-0e5, 1.6e6)
    }
    
    
    sim_data_sum <- sim_data %>%
      group_by(net_strategy, fs_area) %>%
      dplyr::summarise(var_mid = mean(var_den),
                       var_lo = quantile(var_den, 0.025),
                       var_hi = quantile(var_den, 0.975),
                       var_max = max(var_den))
    
    cost_sum <- single_adm %>%
      group_by(net_strategy) %>%
      dplyr::summarise(mean_add_cost = mean(add_cost))
    sim_sum$add_mil_cost <- paste0("$",round(cost_sum$mean_add_cost/1e6,2),"M")
    
    only_label_vals <- c("2 year interval", "3 year interval", "Continous only")
    label_vals <- c("2 year interval", "3 year interval", "Continuous only")
    
    if (costed) {costed_str <- "_costed"} else {costed_str <- "_uncosted"}
    if (costed_and_uncosted) {costed_str <- "_costed_and_uncosted"}
    if (costed_and_uncosted) {label_vals <- c("2 year interval",
                                              "2 year interval (costed)",
                                              "3 year interval",
                                              "3 year interval (costed)")}
    if (costed_and_uncosted) {only_label_vals <- c("2 year interval",
                                                   "2 year interval (costed)",
                                                   "3 year interval",
                                                   "Routine only")}
    
    filename <- paste0(plotting_var, costed_str, i, ".png")
    
    single_adm_no_base <- single_adm %>%
      #filter(net_strategy != baseline_dist) %>%
      filter(net_strategy != "no future nets")
    sim_sum_no_base <- sim_sum %>%
      #filter(net_strategy != baseline_dist) %>%
      filter(net_strategy != "no future nets")
    
    # routine_df <- single_adm %>% filter(net_strategy == baseline_dist)
    # routine_sumdf <- sim_sum %>% filter(net_strategy == baseline_dist)
    # routine_df$net_strategy <- rep("a pyrethroid-only rouine")
    # routine_sumdf$net_strategy <- rep("a pyrethroid-only rouine")
    
    
    only_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-only")
    pbo_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-PBO")
    #pbo_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-PBO")
    
    pyrrole_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-pyrrole")
    
    only_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-only", net_strategy))
    pbo_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-PBO", net_strategy))
    pyrrole_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-pyrrole", net_strategy))
    
    cost_text_size <- 4
    
    #single_adm_no_base %>%
    ggplot() +
      #   geom_violin(data = routine_df,
      #               aes(x = net_strategy,
      #                   y = var_den,
      #                   fill = net_strategy),
      #               alpha = 0.4,
      #               color = NA) +
      #   geom_errorbar(data = routine_sumdf,
      #                 aes(x = net_strategy,
      #                     ymin = var_lo,
      #                     ymax = var_hi,
      #                     color = net_strategy),
    #                 alpha = 0.8) +
    #   geom_point(data = routine_sumdf,
    #              aes(x = net_strategy,
    #                  y = var_mid,
    #                  color = net_strategy),
    #              alpha = 1) +
    #   # geom_text(data = only_sumdf %>% filter (net_strategy == baseline_dist),
    #   #           aes(x = net_strategy,
    #   #               y = var_max,
    #   #               label = add_mil_cost),
    #   #           vjust = -2,
    #   #           size = cost_text_size) +
    #   #   geom_text(data = only_sumdf %>% filter (net_strategy != baseline_dist),
    #   #             aes(x = net_strategy,
    #   #                 y = var_max,
    #   #                 label = add_mil_cost),
    #   #             vjust = -2,
    # #             size = cost_text_size) +
    # scale_fill_viridis(
    #   alpha = 1,
    #   begin = 0,
    #   end = 0,
    #   direction = 1,
    #   discrete = TRUE,
    #   option = "A",
    #   guide = guide_legend(title = "Pyrethroid-only"),
    #   labels = "Routine only"
    # ) +
    #   scale_color_viridis(
    #     alpha = 1,
    #     begin = 0,
    #     end = 0,
    #     direction = 1,
    #     discrete = TRUE,
    #     option = "A",
    #     guide = guide_legend(title = "Pyrethroid-only"),
    #     labels = "Routine only"
    #   ) +
    #   new_scale_fill() +
    #   new_scale_colour() +
    geom_violin(data = only_df,
                aes(x = net_strategy,
                    y = var_den,
                    fill = net_strategy),
                alpha = 0.4,
                color = NA) +
      geom_errorbar(data = only_sumdf,
                    aes(x = net_strategy,
                        ymin = var_lo,
                        ymax = var_hi,
                        color = net_strategy),
                    alpha = 0.8) +
      geom_point(data = only_sumdf,
                 aes(x = net_strategy,
                     y = var_mid,
                     color = net_strategy),
                 alpha = 1) +
      geom_text(data = only_sumdf,
                aes(x = net_strategy,
                    y = var_max,
                    label = add_mil_cost),
                vjust = -1,
                size = cost_text_size) +
      # geom_text(data = only_sumdf %>% filter (net_strategy == baseline_dist),
      #           aes(x = net_strategy,
      #               y = var_max,
      #               label = add_mil_cost),
      #           vjust = -2,
      #           size = cost_text_size) +
      #   geom_text(data = only_sumdf %>% filter (net_strategy != baseline_dist),
      #             aes(x = net_strategy,
      #                 y = var_max,
      #                 label = add_mil_cost),
      #             vjust = -2,
    #             size = cost_text_size) +
    scale_fill_viridis(
      alpha = 1,
      #begin = 0.8,
      #end = 0.95,
      begin = 0.875,
      end = 0.875,
      direction = -1,
      discrete = TRUE,
      option = "H",
      guide = guide_legend(title = "Pyrethroid-only"),
      labels = "3 year interval"#only_label_vals
    ) +
      scale_color_viridis(
        alpha = 1,
        #begin = 0.8,
        #end = 0.95,
        begin = 0.875,
        end = 0.875,
        direction = -1,
        discrete = TRUE,
        option = "H",
        guide = guide_legend(title = "Pyrethroid-only"),
        labels = "3 year interval"#only_label_vals
      ) +
      new_scale_fill() +
      new_scale_colour() +
      geom_violin(data = pbo_df,
                  aes(x = net_strategy,
                      y = var_den,
                      fill = net_strategy),
                  alpha = 0.4,
                  color = NA) +
      geom_errorbar(data = pbo_sumdf,
                    aes(x = net_strategy,
                        ymin = var_lo,
                        ymax = var_hi,
                        color = net_strategy),
                    alpha = 0.8) +
      geom_point(data = pbo_sumdf,
                 aes(x = net_strategy,
                     y = var_mid,
                     color = net_strategy),
                 alpha = 1) +
      geom_text(data = pbo_sumdf,
                aes(x = net_strategy,
                    y = var_max,
                    label = add_mil_cost),
                vjust = -1,
                size = cost_text_size) +
      scale_fill_viridis(
        alpha = 1,
        begin = 0.45,
        end = 0.75,
        direction = 1,
        discrete = TRUE,
        option = "D",
        guide = guide_legend(title = "Pyrethroid-PBO"),
        labels = label_vals
      ) +
      scale_color_viridis(
        alpha = 1,
        begin = 0.45,
        end = 0.75,
        direction = 1,
        discrete = TRUE,
        option = "D",
        guide = guide_legend(title = "Pyrethroid-PBO"),
        labels = label_vals
      ) +
      new_scale_fill() +
      new_scale_colour() +
      geom_violin(data = pyrrole_df,
                  aes(x = net_strategy,
                      y = var_den,
                      fill = net_strategy),
                  alpha = 0.4,
                  color = NA) +
      geom_errorbar(data = pyrrole_sumdf,
                    aes(x = net_strategy,
                        ymin = var_lo,
                        ymax = var_hi,
                        color = net_strategy),
                    alpha = 0.8) +
      geom_point(data = pyrrole_sumdf,
                 aes(x = net_strategy,
                     y = var_mid,
                     color = net_strategy),
                 alpha = 1) +
      geom_text(data = pyrrole_sumdf,
                aes(x = net_strategy,
                    y = var_max,
                    label = add_mil_cost),
                vjust = -1,
                size = cost_text_size) +
      scale_fill_viridis(
        alpha = 1,
        #begin = 0.05,
        begin = 0.2,
        end = 0.2,
        direction = 1,
        discrete = TRUE,
        option = "H",
        guide = guide_legend(title = "Pyrethroid-pyrrole"),
        labels = "2 year interval"#label_vals
      ) +
      scale_color_viridis(
        alpha = 1,
        #begin = 0.05,
        begin = 0.2,
        end = 0.2,
        direction = 1,
        discrete = TRUE,
        option = "H",
        guide = guide_legend(title = "Pyrethroid-pyrrole"),
        labels = "2 year interval"#label_vals
      ) +
      # guides(fill = guide_legend(title = "Net Strategy")) +
      # guides(colour = guide_legend(title = "Net Strategy")) +
      # labs(colour = "Pyrethroid-PBO",
      #      y = "add_cases_averted_per_usd") +
      ylab(yax_val) +
      xlab(NULL) +
      scale_y_continuous(limits = ylim_vals,
                         labels = label_comma()) + #,
      #breaks = round(seq(ylim_vals[1], ylim_vals[2], by = 0.25),2)) +
      scale_x_discrete(breaks = NULL) +
      ggtitle(fs_areas_included[i])
    
    # ggsave(filename, bg = "white",
    #        w = 10, h = 4, dpi = 450)
    
    ggsave(filename, bg = "white",
           w = 4, h = 3, dpi = 450)
    
    
    #print to screen
    print(fs_areas_included[i])
    
    print(paste("only_var_lo =", only_sumdf$var_lo))
    print(paste("only_var_mid =", only_sumdf$var_mid))
    print(paste("only_var_hi =", only_sumdf$var_hi))
    # print(paste("only_pc_lo =", (only_sumdf$var_lo[1] - only_sumdf$var_lo[2])/only_sumdf$var_lo[2]))
    # print(paste("only_pc_mid =", (only_sumdf$var_mid[1] - only_sumdf$var_mid[2])/only_sumdf$var_mid[2]))
    # print(paste("only_pc_hi =", (only_sumdf$var_hi[1] - only_sumdf$var_hi[2])/only_sumdf$var_hi[2]))
    # 
    # print(paste("pbo_var_lo =", pbo_sumdf$var_lo))
    # print(paste("pbo_var_mid =", pbo_sumdf$var_mid))
    # print(paste("pbo_var_hi =", pbo_sumdf$var_hi))
    # print(paste("pbo_pc_lo =", (pbo_sumdf$var_lo[1] - pbo_sumdf$var_lo[2])/pbo_sumdf$var_lo[2]))
    # print(paste("pbo_pc_mid =", (pbo_sumdf$var_mid[1] - pbo_sumdf$var_mid[2])/pbo_sumdf$var_mid[2]))
    # print(paste("pbo_pc_hi =", (pbo_sumdf$var_hi[1] - pbo_sumdf$var_hi[2])/pbo_sumdf$var_hi[2]))
    # 
    print(paste("pyrrole_var_lo =", pyrrole_sumdf$var_lo))
    print(paste("pyrrole_var_mid =", pyrrole_sumdf$var_mid))
    print(paste("pyrrole_var_hi =", pyrrole_sumdf$var_hi))
    # print(paste("pyrrole_pc_lo =", (pyrrole_sumdf$var_lo[1] - pyrrole_sumdf$var_lo[2])/pyrrole_sumdf$var_lo[2]))
    # print(paste("pyrrole_pc_mid =", (pyrrole_sumdf$var_mid[1] - pyrrole_sumdf$var_mid[2])/pyrrole_sumdf$var_mid[2]))
    # print(paste("pyrrole_pc_hi =", (pyrrole_sumdf$var_hi[1] - pyrrole_sumdf$var_hi[2])/pyrrole_sumdf$var_hi[2]))
    # 
    
    # ggsave(filename, bg = "white",
    #        w = 7, h = 5, dpi = 450)
    
    
    
#  }
  
}

sim_violin_plot_cost_usage <- function(sim_data,
                                       sim_sum,
                                       country = "SN",
                                       plotting_var = "avert",
                                       costed = FALSE,
                                       avg_yr1_usage = FALSE,
                                       per_xpop = 1,
                                       baseline_dist = "pyrethroid-only routine baseline") {
  
  # sim_data %<>% filter(ISO2 == country & fs_name_1 == "Dakar" & urbanicity == "urban")
  # sim_sum %<>% filter(ISO2 == country & fs_name_1 == "Dakar" & urbanicity == "urban")
  
  sim_data %<>% filter(ISO2 == country)
  sim_sum %<>% filter(ISO2 == country)
  
  if (plotting_var == "avert") {
    sim_data$var_den <- sim_data$cases_averted_per_pop * per_xpop
    sim_sum$var_lo <- sim_sum$LB_avert_percap * per_xpop
    sim_sum$var_mid <- sim_sum$mean_avert_percap * per_xpop
    sim_sum$var_hi <- sim_sum$UB_avert_percap * per_xpop
    yax_val <- "Annual cases averted"
  } else if (plotting_var == "add_avert") {
    sim_data$var_den <- sim_data$add_cases_averted_per_pop * per_xpop
    sim_sum$var_lo <- sim_sum$LB_add_avert_percap * per_xpop
    sim_sum$var_mid <- sim_sum$mean_add_avert_percap * per_xpop
    sim_sum$var_hi <- sim_sum$UB_add_avert_percap * per_xpop
    yax_val <- "Additional annual cases averted"
  }
  
  if (per_xpop == 1) {
    yax_val %<>% paste("per capita")
  } else {
    yax_val %<>% paste("per", per_xpop)
  }
  
  sim_sum$mil_cost <- paste0("$",round(sim_sum$mean_avg_ann_net_cost/1e6,2),"M")
  
  only_label_vals <- c("2 year interval", "3 year interval")
  label_vals <- c("2 year interval", "3 year interval")
  
  text_size <- 3
  text_angle <- 90
  
  sim_sum$cost_ypos <- sim_sum$var_hi + per_xpop / 5
  sim_sum$use_ypos <- sim_sum$var_lo - per_xpop / 8
  
  sim_sum$mean_use <- paste0(round(sim_sum$mean_avg_yr1_use, 2) * 100, "%")
  
  only_data <- sim_data %>% filter(net_name == "pyrethroid-only")
  pbo_data <- sim_data %>% filter(net_name == "pyrethroid-PBO")
  pyrrole_data <- sim_data %>% filter(net_name == "pyrethroid-pyrrole")
  
  only_sum <- sim_sum %>% filter(net_name == "pyrethroid-only")
  pbo_sum <- sim_sum %>% filter(net_name == "pyrethroid-PBO")
  pyrrole_sum <- sim_sum %>% filter(net_name == "pyrethroid-pyrrole")
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  N_a <- length(unique(sim_sum$fs_name_1))
  reg_strip <- strip_themed(background_x = elem_list_rect(fill = gg_color_hue(N_a)))
  
  plt <- ggplot() +
    geom_violin(data = only_data,
                aes(x = net_strategy,
                    y = var_den,
                    fill = net_strategy),
                alpha = 0.4,
                color = NA) +
    geom_errorbar(data = only_sum,
                  aes(x = net_strategy,
                      ymin = var_lo,
                      ymax = var_hi,
                      color = net_strategy),
                  alpha = 0.8) +
    geom_point(data = only_sum,
               aes(x = net_strategy,
                   y = var_mid,
                   color = net_strategy),
               alpha = 1) +
    geom_text(data = only_sum,
              aes(x = net_strategy,
                  y = cost_ypos,
                  label = mil_cost),
              vjust = 0.5,
              size = text_size,
              angle = text_angle) +
    geom_text(data = only_sum,
              aes(x = net_strategy,
                  y = use_ypos,
                  label = mean_use),
              vjust = 0.5,
              size = text_size,
              angle = text_angle) +
    scale_fill_viridis(
      alpha = 1,
      begin = 0.8,
      end = 0.95,
      # begin = 0.875,
      # end = 0.875,
      direction = -1,
      discrete = TRUE,
      option = "H",
      guide = guide_legend(title = "Pyrethroid-only:",
                           order = 1),
      labels = only_label_vals
    ) +
    scale_color_viridis(
      alpha = 1,
      begin = 0.8,
      end = 0.95,
      #begin = 0.875,
      #end = 0.875,
      direction = -1,
      discrete = TRUE,
      option = "H",
      guide = guide_legend(title = "Pyrethroid-only:",
                           order = 1),
      labels = only_label_vals
    ) +
    new_scale_fill() +
    new_scale_colour() +
    geom_violin(data = pbo_data,
                aes(x = net_strategy,
                    y = var_den,
                    fill = net_strategy),
                alpha = 0.4,
                color = NA) +
    geom_errorbar(data = pbo_sum,
                  aes(x = net_strategy,
                      ymin = var_lo,
                      ymax = var_hi,
                      color = net_strategy),
                  alpha = 0.8) +
    geom_point(data = pbo_sum,
               aes(x = net_strategy,
                   y = var_mid,
                   color = net_strategy),
               alpha = 1) +
    geom_text(data = pbo_sum,
              aes(x = net_strategy,
                  y = cost_ypos,
                  label = mil_cost),
              vjust = 0.5,
              size = text_size,
              angle = text_angle) +
    geom_text(data = pbo_sum,
              aes(x = net_strategy,
                  y = use_ypos,
                  label = mean_use),
              vjust = 0.5,
              size = text_size,
              angle = text_angle) +
    scale_fill_viridis(
      alpha = 1,
      begin = 0.45,
      end = 0.75,
      direction = 1,
      discrete = TRUE,
      option = "D",
      guide = guide_legend(title = "Pyrethroid-PBO:",
                           order = 2),
      labels = label_vals
    ) +
    scale_color_viridis(
      alpha = 1,
      begin = 0.45,
      end = 0.75,
      direction = 1,
      discrete = TRUE,
      option = "D",
      guide = guide_legend(title = "Pyrethroid-PBO:",
                           order = 2),
      labels = label_vals
    ) +
    new_scale_fill() +
    new_scale_colour() +
    geom_violin(data = pyrrole_data,
                aes(x = net_strategy,
                    y = var_den,
                    fill = net_strategy),
                alpha = 0.4,
                color = NA) +
    geom_errorbar(data = pyrrole_sum,
                  aes(x = net_strategy,
                      ymin = var_lo,
                      ymax = var_hi,
                      color = net_strategy),
                  alpha = 0.8) +
    geom_point(data = pyrrole_sum,
               aes(x = net_strategy,
                   y = var_mid,
                   color = net_strategy),
               alpha = 1) +
    geom_text(data = pyrrole_sum,
              aes(x = net_strategy,
                  y = cost_ypos,
                  label = mil_cost),
              vjust = 0.5,
              size = text_size,
              angle = text_angle) +
    geom_text(data = pyrrole_sum,
              aes(x = net_strategy,
                  y = use_ypos,
                  label = mean_use),
              vjust = 0.5,
              size = text_size,
              angle = text_angle) +
    scale_fill_viridis(
      alpha = 1,
      begin = 0.05,
      #begin = 0.2,
      end = 0.2,
      direction = 1,
      discrete = TRUE,
      option = "H",
      guide = guide_legend(title = "Pyrethroid-pyrrole:",
                           order = 3),
      labels = label_vals
    ) +
    scale_color_viridis(
      alpha = 1,
      begin = 0.05,
      #begin = 0.2,
      end = 0.2,
      direction = 1,
      discrete = TRUE,
      option = "H",
      guide = guide_legend(title = "Pyrethroid-pyrrole:",
                           order = 3),
      labels = label_vals
    ) +
    # guides(fill = guide_legend(title = "Net Strategy")) +
    # guides(colour = guide_legend(title = "Net Strategy")) +
    # labs(colour = "Pyrethroid-PBO",
    #      y = "add_cases_averted_per_usd") +
    ylab(yax_val) +
    xlab(NULL) +
    scale_y_continuous(limits = c(-0.2,1.7)) +
    scale_x_discrete(breaks = NULL) +
    theme(legend.position="bottom") +
    facet_grid2(urbanicity ~ fs_name_1, strip = reg_strip)
  
  ggsave("SN_avert_cost_yr1use_coltest.png", bg = "white",
         w = 16, h = 6, dpi = 450)
  
  
  
  # scale_y_continuous(limits = ylim_vals,
  #                    labels = label_comma()) + #,
  #breaks = round(seq(ylim_vals[1], ylim_vals[2], by = 0.25),2)) +
  # scale_x_discrete(breaks = NULL) +
  # ggtitle(fs_areas_included[i])
  
  # ggsave(filename, bg = "white",
  #        w = 10, h = 4, dpi = 450)
  
  # ggsave(filename, bg = "white",
  #        w = 4, h = 3, dpi = 450)
  
  print(plt)
  
  
  
  
}

