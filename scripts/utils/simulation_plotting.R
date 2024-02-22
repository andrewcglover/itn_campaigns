# net_cases_plot_test.R

plot_sim_bars <- function(sim_data,
                          fs_areas_included = NULL,
                          plotting_var = NULL,
                          #plot_cases_averted = FALSE,
                          baseline_dist = "pyrethroid-only 3 year interval") {

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
      single_adm$cases_averted <- single_adm$baseline_infections - single_adm$pred_ann_infect
      sim_sum <- single_adm %>%
        group_by(net_strategy) %>%
        dplyr::summarise(var_mid = mean(cases_averted),
                         var_lo = quantile(cases_averted, 0.025),
                         var_hi = quantile(cases_averted, 0.975))
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

    filename <- paste0(plotting_var, i, ".png")
    
    ggplot(sim_sum) +
      geom_col(aes(x = net_strategy,
                   y = var_mid,
                   fill = net_strategy),
               alpha = 0.4) +
      geom_errorbar(aes(x = net_strategy,
                        ymin = var_lo,
                        ymax = var_hi,
                        color = net_strategy),
                    alpha = 0.9) +
      ylab(yax_val) + 
      xlab(NULL) +
      #scale_y_continuous(limits = c(0, 0.6)) +
      scale_x_discrete(breaks = NULL) +
      #guides(fill = guide_colorsteps(title = "Distribution strategy"))
      # guides(fill = guide_legend(title = "Distribution strategy"))
      guides(fill = guide_legend(title = "Net Strategy")) +
      guides(color = guide_legend(title = "Net Strategy")) +
      ggtitle(fs_areas_included[i])
    
    ggsave(filename, bg = "white",
           w = 5, h = 4, dpi = 450)
    
    
    
  }

}


