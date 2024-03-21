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
      basline_ids <- which(single_adm$net_strategy == baseline_dist)
      baseline_infecions <- single_adm$pred_ann_infect[basline_ids]
      single_adm$baseline_infections <- rep(baseline_infecions,
                                            length.out = dim(single_adm)[1])
      single_adm$var_den <- single_adm$baseline_infections - single_adm$pred_ann_infect
      sim_sum <- single_adm %>%
        group_by(net_strategy) %>%
        dplyr::summarise(var_mid = mean(var_den),
                         var_lo = quantile(var_den, 0.025),
                         var_hi = quantile(var_den, 0.975),
                         var_max = max(var_den))
      yax_val <- "Annual cases averted"
      ylim_vals <- c(-0.5e5, 1.5e6)
    } else if (plotting_var == "add_cases_averted_per_usd") {
      # basline_ids <- which(single_adm$net_strategy == baseline_dist)
      # baseline_infecions <- single_adm$pred_ann_infect[basline_ids]
      # single_adm$baseline_infections <- rep(baseline_infecions,
      #                                       length.out = dim(single_adm)[1])
      # single_adm$cases_averted <- single_adm$baseline_infections - single_adm$pred_ann_infect
      sim_sum <- single_adm %>%
        group_by(net_strategy) %>%
        dplyr::summarise(var_mid = mean(add_cases_averted_per_usd),
                         var_lo = quantile(add_cases_averted_per_usd, 0.025),
                         var_hi = quantile(add_cases_averted_per_usd, 0.975))
      yax_val <- "Cases averted vs triennial pyrethroid-only strategy"
    } else if (plotting_var == "cases_averted_pp") {
      basline_ids <- which(single_adm$net_strategy == baseline_dist)
      baseline_infecions <- single_adm$pred_ann_infect[basline_ids]
      single_adm$baseline_infections <- rep(baseline_infecions,
                                            length.out = dim(single_adm)[1])
      single_adm$cases_averted_pp <- 10000 * (single_adm$baseline_infections - single_adm$pred_ann_infect) / single_adm$pop
      sim_sum <- single_adm %>%
        group_by(net_strategy) %>%
        dplyr::summarise(var_mid = mean(cases_averted_pp),
                         var_lo = quantile(cases_averted_pp, 0.025),
                         var_hi = quantile(cases_averted_pp, 0.975))
      yax_val <- "Cases averted per 10,000 vs triennial pyrethroid-only strategy"
    } else if (plotting_var == "pred_ann_infect") {
      sim_sum <- single_adm %>%
        group_by(net_strategy) %>%
        dplyr::summarise(var_mid = mean(pred_ann_infect),
                         var_lo = quantile(pred_ann_infect, 0.025),
                         var_hi = quantile(pred_ann_infect, 0.975))
      yax_val <- "Mean annual cases"
    } else if (plotting_var == "avg_pfpr") {
      single_adm$avg_pfpr_scaled <- single_adm$avg_pfpr * 10000
      sim_sum <- single_adm %>%
        group_by(net_strategy) %>%
        dplyr::summarise(var_mid = mean(avg_pfpr_scaled),
                         var_lo = quantile(avg_pfpr_scaled, 0.025),
                         var_hi = quantile(avg_pfpr_scaled, 0.975))
      yax_val <- "Prevalence per 10,000 in 2-10 year olds"
    } else if (plotting_var == "ann_incidence") {
      single_adm$ann_inc_scaled <- single_adm$ann_incidence * 10000
      sim_sum <- single_adm %>%
        group_by(net_strategy) %>%
        dplyr::summarise(var_mid = mean(ann_inc_scaled),
                         var_lo = quantile(ann_inc_scaled, 0.025),
                         var_hi = quantile(ann_inc_scaled, 0.975))
      yax_val <- "Incidence per 10,000 for all ages"
    } else if (plotting_var == "avg_ann_nets_distrib") {
      sim_sum <- single_adm %>%
        group_by(net_strategy) %>%
        dplyr::summarise(var_mid = mean(avg_ann_nets_distrib),
                         var_lo = quantile(avg_ann_nets_distrib, 0.025),
                         var_hi = quantile(avg_ann_nets_distrib, 0.975))
      yax_val <- "Mean total nets distributed per year"
    } else if (plotting_var == "ann_nets_pp") {
      single_adm$ann_nets_pp <- single_adm$avg_ann_nets_distrib / single_adm$pop
      sim_sum <- single_adm %>%
        group_by(net_strategy) %>%
        dplyr::summarise(var_mid = mean(ann_nets_pp),
                         var_lo = quantile(ann_nets_pp, 0.025),
                         var_hi = quantile(ann_nets_pp, 0.975))
      yax_val <- "Mean total nets distributed per person per year"
    } else if (plotting_var == "avert_per_net") {
      basline_ids <- which(single_adm$net_strategy == baseline_dist)
      baseline_infecions <- single_adm$pred_ann_infect[basline_ids]
      single_adm$baseline_infections <- rep(baseline_infecions,
                                            length.out = dim(single_adm)[1])
      single_adm$cases_averted <- single_adm$baseline_infections - single_adm$pred_ann_infect
      single_adm$avert_per_net <- 100000 * single_adm$cases_averted / single_adm$avg_ann_nets_distrib
      sim_sum <- single_adm %>%
        group_by(net_strategy) %>%
        dplyr::summarise(var_mid = mean(avert_per_net),
                         var_lo = quantile(avert_per_net, 0.025),
                         var_hi = quantile(avert_per_net, 0.975))
      yax_val <- "Cases averted vs baseline per 100,000 nets"
    }
    
    cost_sum <- single_adm %>%
      group_by(net_strategy) %>%
      dplyr::summarise(mean_add_cost = mean(add_cost))
    sim_sum$add_mil_cost <- paste0("$",round(cost_sum$mean_add_cost/1e6,2),"M")
    
    label_vals <- c("2 year interval", "3 year interval")
    
    if (costed) {costed_str <- "_costed"} else {costed_str <- "_uncosted"}
    if (costed_and_uncosted) {costed_str <- "_costed_and_uncosted"}
    if (costed_and_uncosted) {label_vals <- c("2 year interval",
                                              "2 year interval (costed)",
                                              "3 year interval",
                                              "3 year interval (costed)")}
    
    filename <- paste0(plotting_var, costed_str, i, ".png")
    
    single_adm_no_base <- single_adm %>% filter(net_strategy != baseline_dist)
    sim_sum_no_base <- sim_sum %>% filter(net_strategy != baseline_dist)
    
    routine_df <- single_adm %>% filter(net_strategy == baseline_dist)
    routine_sumdf <- single_adm %>% filter(net_strategy == baseline_dist)
    
    only_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-only")
    pbo_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-PBO")
    pyrrole_df <- single_adm_no_base %>% filter(net_name == "pyrethroid-pyrrole")
    
    only_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-only", net_strategy))
    pbo_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-PBO", net_strategy))
    pyrrole_sumdf <- sim_sum_no_base %>% filter(grepl("pyrethroid-pyrrole", net_strategy))
    
    cost_text_size <- 2.5
    
    #single_adm_no_base %>%
      ggplot() +
      geom_violin(data = routine_df,
                  aes(x = net_strategy,
                      y = pred_ann_infect,
                      fill = net_strategy),
                  alpha = 0.4,
                  color = NA) +
      # geom_errorbar(data = routine_sumdf,
      #               aes(x = net_strategy,
      #                   ymin = var_lo,
      #                   ymax = var_hi,
      #                   color = net_strategy),
      #               alpha = 0.8) +
      # geom_point(data = routine_sumdf,
      #            aes(x = net_strategy,
      #                y = var_mid,
      #                color = net_strategy),
      #            alpha = 1) +
      # geom_text(data = routine_sumdf,
      #           aes(x = net_strategy,
      #               y = var_max,
      #               label = add_mil_cost),
      #           vjust = -1,
      #           size = cost_text_size) +
      scale_fill_viridis(
        alpha = 1,
        begin = 0,
        end = 0,
        direction = 1,
        discrete = TRUE,
        option = "A",
        guide = guide_legend(title = "Routine-only"),
        labels = NULL
      ) +
      scale_color_viridis(
        alpha = 1,
        begin = 0,
        end = 0,
        direction = 1,
        discrete = TRUE,
        option = "A",
        guide = guide_legend(title = "Routine-only"),
        labels = NULL
      ) +
      new_scale_fill() +
      new_scale_colour() +
      geom_violin(data = only_df,
                  aes(x = net_strategy,
                      y = var_den,
                      fill = net_strategy),
                  alpha = 0.4,
                  color = NA) +
      geom_violin(data = only_df %>% filter (net_strategy == "pyrethroid-only 3 year interval"),
                  aes(x = net_strategy,
                      y = var_den,
                      fill = net_strategy),
                  alpha = 0.2,
                  show.legend = FALSE) +
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
    scale_fill_viridis(
      alpha = 1,
      begin = 0.8,
      end = 0.95,
      direction = -1,
      discrete = TRUE,
      option = "H",
      guide = guide_legend(title = "Pyrethroid-only"),
      labels = label_vals
      ) +
      scale_color_viridis(
        alpha = 1,
        begin = 0.8,
        end = 0.95,
        direction = -1,
        discrete = TRUE,
        option = "H",
        guide = guide_legend(title = "Pyrethroid-only"),
        labels = label_vals
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
                         labels = label_comma()) +
      scale_x_discrete(breaks = NULL) +
    ggtitle(fs_areas_included[i])
    
    ggsave(filename, bg = "white",
           w = 7, h = 5, dpi = 450)
    
    
    
  }
  
}


