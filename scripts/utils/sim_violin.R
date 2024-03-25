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
      single_adm$var_den <- single_adm$baseline_cases - single_adm$pred_ann_infect
      yax_val <- "Annual cases averted"
      ylim_vals <- c(-0e5, 6.5e5)
    } else if (plotting_var == "add_cases_averted_per_usd") {
      single_adm$var_den <- single_adm$add_cases_averted_per_usd
      yax_val <- "Annual cases averted per USD"
      ylim_vals <- c(0, 1.75)
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
    
    only_label_vals <- c("2 year interval", "3 year interval")
    label_vals <- c("2 year interval", "3 year interval")
    
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
      filter(net_strategy != baseline_dist) %>%
      filter(net_strategy != "no future nets")
    sim_sum_no_base <- sim_sum %>%
      filter(net_strategy != baseline_dist) %>%
      filter(net_strategy != "no future nets")
    
    routine_df <- single_adm %>% filter(net_strategy == baseline_dist)
    routine_sumdf <- sim_sum %>% filter(net_strategy == baseline_dist)
    routine_df$net_strategy <- rep("a pyrethroid-only rouine")
    routine_sumdf$net_strategy <- rep("a pyrethroid-only rouine")
    
    
    only_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-only")
    pbo_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-PBO")
    #pbo_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-PBO")
    
    pyrrole_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-pyrrole")
    
    only_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-only", net_strategy))
    pbo_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-PBO", net_strategy))
    pyrrole_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-pyrrole", net_strategy))
    
    cost_text_size <- 3
    
    #single_adm_no_base %>%
      ggplot() +
        geom_violin(data = routine_df,
                    aes(x = net_strategy,
                        y = var_den,
                        fill = net_strategy),
                    alpha = 0.4,
                    color = NA) +
        geom_errorbar(data = routine_sumdf,
                      aes(x = net_strategy,
                          ymin = var_lo,
                          ymax = var_hi,
                          color = net_strategy),
                      alpha = 0.8) +
        geom_point(data = routine_sumdf,
                   aes(x = net_strategy,
                       y = var_mid,
                       color = net_strategy),
                   alpha = 1) +
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
        begin = 0,
        end = 0,
        direction = 1,
        discrete = TRUE,
        option = "A",
        guide = guide_legend(title = "Pyrethroid-only"),
        labels = "Routine only"
      ) +
        scale_color_viridis(
          alpha = 1,
          begin = 0,
          end = 0,
          direction = 1,
          discrete = TRUE,
          option = "A",
          guide = guide_legend(title = "Pyrethroid-only"),
          labels = "Routine only"
        ) +
        new_scale_fill() +
        new_scale_colour() +
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
        begin = 0.8,
        end = 0.95,
        direction = -1,
        discrete = TRUE,
        option = "H",
        guide = guide_legend(title = "Pyrethroid-only"),
        labels = only_label_vals
        ) +
        scale_color_viridis(
          alpha = 1,
          begin = 0.8,
          end = 0.95,
          direction = -1,
          discrete = TRUE,
          option = "H",
          guide = guide_legend(title = "Pyrethroid-only"),
          labels = only_label_vals
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
        # geom_text(data = pbo_sumdf,
        #           aes(x = net_strategy,
        #               y = var_max,
        #               label = add_mil_cost),
        #           vjust = -1,
        #           size = cost_text_size) +
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
        # geom_text(data = pyrrole_sumdf,
        #           aes(x = net_strategy,
        #               y = var_max,
        #               label = add_mil_cost),
        #           vjust = -1,
        #           size = cost_text_size) +
        scale_fill_viridis(
          alpha = 1,
          begin = 0.05,
          end = 0.2,
          direction = 1,
          discrete = TRUE,
          option = "H",
          guide = guide_legend(title = "Pyrethroid-pyrrole"),
          labels = label_vals
        ) +
        scale_color_viridis(
          alpha = 1,
          begin = 0.05,
          end = 0.2,
          direction = 1,
          discrete = TRUE,
          option = "H",
          guide = guide_legend(title = "Pyrethroid-pyrrole"),
          labels = label_vals
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
    
      ggsave(filename, bg = "white",
             w = 8, h = 3, dpi = 450)
      
    # ggsave(filename, bg = "white",
    #        w = 7, h = 5, dpi = 450)
    
    
    
  }
  
}


