cases_averted_scatter <- function(sim_sum,
                                  country = NULL,
                                  rm.country = NULL,
                                  var_name = "pfeir",
                                  per_xpop = 1000,
                                  only3_comparison = FALSE){
  
  sim_sum$country <- countrycode(sim_sum$ISO2,"iso2c","country.name")
  
  if (!is.null(country)) {sim_sum %<>% filter(ISO2 == country)}
  if (!is.null(rm.country)) {sim_sum %<>% filter(ISO2 != rm.country)}
  
  if (only3_comparison) {
    sim_sum %<>% filter(net_strategy != "pyrethroid-only 3 year interval")
    only2_col <- 0.95
    only3_col <- 0.8
    sim_sum %<>% cbind(
      "mean_cases_avert_plt" = sim_sum$mean_add_avert_percap * per_xpop,
      "LB_cases_avert_plt" = sim_sum$LB_add_avert_percap * per_xpop,
      "UB_cases_avert_plt" = sim_sum$UB_add_avert_percap * per_xpop
      )
  } else {
    only2_col <- 0.95
    only3_col <- 0.95
    sim_sum %<>% cbind(
      "mean_cases_avert_plt" = sim_sum$mean_avert_percap * per_xpop,
      "LB_cases_avert_plt" = sim_sum$LB_avert_percap * per_xpop,
      "UB_cases_avert_plt" = sim_sum$UB_avert_percap * per_xpop
    )
  }
  
  if (var_name == "pfeir") {sim_sum$var_val <- sim_sum$pfeir}
  if (var_name == "pyrethroid_resistance") {sim_sum$var_val <- sim_sum$pyrethroid_resistance}
  if (var_name == "uret") {sim_sum$var_val <- sim_sum$mean_retu}

  
  # sim_sum$var_LB <- sim_sum$var_val
  # sim_sum$var_UB <- sim_sum$var_val
  
  #sim_sum$LB_avert_percap[which(sim_sum$LB_avert_percap < 0.01)] <- 0.01
  

  
  pos <- position_dodge(width = 0)
  alpha_val <- 0.5
  size_val <- 0.4
  label_vals <- c("2 year interval", "3 year interval")
  if (var_name == "uret") {
    plt <- ggplot() +
      geom_errorbar(
        data = sim_sum %>% filter(net_name == "pyrethroid-only"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          ymin = LB_cases_avert_plt,
          ymax = UB_cases_avert_plt,
          color = net_strategy,
          shape = urbanicity
        ),
        position = pos,
        alpha = alpha_val,
        size = size_val,
        stroke=NA
      ) +
      geom_errorbarh(
        data = sim_sum %>% filter(net_name == "pyrethroid-only"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          xmin = LB_retu,
          xmax = UB_retu,
          color = net_strategy,
          shape = urbanicity
        ),
        position = pos,
        alpha = alpha_val,
        size = size_val,
        stroke=NA
      ) +
      scale_colour_viridis(
        alpha = 1,
        begin = 0.8,#only3_col,
        end = 0.95,#only2_col,
        direction = -1,
        discrete = TRUE,
        option = "H",
        guide = guide_legend(title = "Pyrethroid-only", order = 1),
        labels = label_vals
      ) +
      new_scale_colour() +
      geom_errorbar(
        data = sim_sum %>% filter(net_name == "pyrethroid-PBO"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          ymin = LB_cases_avert_plt,
          ymax = UB_cases_avert_plt,
          color = net_strategy,
          shape = urbanicity
        ),
        position = pos,
        alpha = alpha_val,
        size = size_val,
        stroke=NA
      ) +
      geom_errorbarh(
        data = sim_sum %>% filter(net_name == "pyrethroid-PBO"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          xmin = LB_retu,
          xmax = UB_retu,
          color = net_strategy,
          shape = urbanicity
        ),
        position = pos,
        alpha = alpha_val,
        size = size_val,
        stroke=NA
      ) +
      scale_colour_viridis(
        alpha = 1,
        begin = 0.45,
        end = 0.75,
        direction = 1,
        discrete = TRUE,
        option = "D",
        guide = guide_legend(title = "Pyrethroid-PBO", order = 2),
        labels = label_vals
      ) +
      new_scale_colour() +
      geom_errorbar(
        data = sim_sum %>% filter(net_name == "pyrethroid-pyrrole"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          ymin = LB_cases_avert_plt,
          ymax = UB_cases_avert_plt,
          color = net_strategy,
          shape = urbanicity
        ),
        position = pos,
        alpha = alpha_val,
        size = size_val,
        stroke=NA
      ) +
      geom_errorbarh(
        data = sim_sum %>% filter(net_name == "pyrethroid-pyrrole"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          xmin = LB_retu,
          xmax = UB_retu,
          color = net_strategy,
          shape = urbanicity
        ),
        position = pos,
        alpha = alpha_val,
        size = size_val,
        stroke=NA
      ) +
      scale_colour_viridis(
        alpha = 1,
        begin = 0.05,
        end = 0.2,
        direction = 1,
        discrete = TRUE,
        option = "H",
        guide = guide_legend(title = "Pyrethroid-pyrrole", order = 3),
        labels = label_vals
      )
  } else {
    plt <- ggplot() +
      geom_pointrange(
        data = sim_sum %>% filter(net_name == "pyrethroid-only"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          ymin = LB_cases_avert_plt,
          ymax = UB_cases_avert_plt,
          color = net_strategy,
          shape = urbanicity
        ),
        position = pos,
        alpha = alpha_val,
        size = size_val,
        stroke=NA
      ) +
      scale_colour_viridis(
        alpha = 1,
        begin = 0.8,#only3_col,
        end = 0.95,#only2_col,
        direction = -1,
        discrete = TRUE,
        option = "H",
        guide = guide_legend(title = "Pyrethroid-only", order = 1),
        labels = label_vals
      ) +
      new_scale_colour() +
      geom_pointrange(
        data = sim_sum %>% filter(net_name == "pyrethroid-PBO"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          ymin = LB_cases_avert_plt,
          ymax = UB_cases_avert_plt,
          color = net_strategy,
          shape = urbanicity
        ),
        position = pos,
        alpha = alpha_val,
        size = size_val,
        stroke=NA
      ) +
      scale_colour_viridis(
        alpha = 1,
        begin = 0.45,
        end = 0.75,
        direction = 1,
        discrete = TRUE,
        option = "D",
        guide = guide_legend(title = "Pyrethroid-PBO", order = 2),
        labels = label_vals
        ) +
      new_scale_colour() +
      geom_pointrange(
        data = sim_sum %>% filter(net_name == "pyrethroid-pyrrole"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          ymin = LB_cases_avert_plt,
          ymax = UB_cases_avert_plt,
          color = net_strategy,
          shape = urbanicity
        ),
        position = pos,
        alpha = alpha_val,
        size = size_val,
        stroke=NA
      ) +
      scale_colour_viridis(
        alpha = 1,
        begin = 0.05,
        end = 0.2,
        direction = 1,
        discrete = TRUE,
        option = "H",
        guide = guide_legend(title = "Pyrethroid-pyrrole", order = 3),
        labels = label_vals
        )
  }
  
  if (var_name == "pfeir") {plt <- plt + scale_x_continuous(trans='log10')}
  if (var_name == "pfeir") {plt <- plt + xlab("PfEIR")}
  if (var_name == "pyrethroid_resistance") {plt <- plt + xlab("Pyrethroid resistance")}
  if (var_name == "uret") {plt <- plt + xlab("Mean used net retention (months)")}
  if (var_name == "uga") {plt <- plt + xlab("Mean usage given access")}
  if (var_name == "pyrethroid_resistance") {plt <- plt + xlab("Pyrethroid resistance")}
  if (per_xpop == 1) {per_xpop <- "capita"}
  ylab_base <- paste0("annual cases averted per ", per_xpop)
  if (only3_comparison) {plt <- plt + ylab(paste0("Additional mean ", ylab_base))}
  if (!only3_comparison) {plt <- plt + ylab(paste0("Mean ", ylab_base))}
  
  plt <- plt + facet_grid(cols = vars(net_name),
                          rows = vars(country),
                          scales="free_y")

  print(plt)
  
  if(only3_comparison) {add_str <- "add"} else {add_str <- ""}
  #ggsave(paste0("allctry_", add_str, "caseavrt_vs_",var_name,".png"), bg = "white",
  #       w = 6, h = 8, dpi = 450)
  
}