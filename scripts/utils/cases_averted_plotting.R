cases_averted_scatter <- function(sim_sum,
                                  country = NULL,
                                  var_name = "pfeir",
                                  per_xpop = 1000,
                                  only3_comparison = FALSE){
  
  sim_sum$country <- countrycode(sim_sum$ISO2,"iso2c","country.name")
  
  if (!is.null(country)) {sim_sum %<>% filter(ISO2 == country)}
  
  if (only3_comparison) {
    sim_sum %<>% filter(net_strategy != "pyrethroid-only 3 year interval")
    only2_col <- 0.95
    only3_col <- 0.8
  } else {
    only2_col <- 0.95
    only3_col <- 0.95
  }
  
  if (var_name == "pfeir") {sim_sum$var_val <- sim_sum$pfeir}
  if (var_name == "pyrethroid_resistance") {sim_sum$var_val <- sim_sum$pyrethroid_resistance}
  
  #sim_sum$LB_avert_percap[which(sim_sum$LB_avert_percap < 0.01)] <- 0.01
  

  
  pos <- position_dodge(width = 0)
  alpha_val <- 0.6
  size_val <- 0.6
  label_vals <- c("2 year interval", "3 year interval")
  plt <- ggplot() +
    geom_pointrange(
      data = sim_sum %>% filter(net_name == "pyrethroid-only"),
      aes(
        x = var_val,
        y = mean_avert_percap * per_xpop,
        ymin = LB_avert_percap * per_xpop,
        ymax = UB_avert_percap * per_xpop,
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
      begin = only3_col,
      end = only2_col,
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
        y = mean_avert_percap * per_xpop,
        ymin = LB_avert_percap * per_xpop,
        ymax = UB_avert_percap * per_xpop,
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
        y = mean_avert_percap * per_xpop,
        ymin = LB_avert_percap * per_xpop,
        ymax = UB_avert_percap * per_xpop,
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
  
  if (var_name == "pfeir") {plt <- plt + scale_x_continuous(trans='log10')}
  if (var_name == "pfeir") {plt <- plt + xlab("PfEIR")}
  if (var_name == "pyrethroid_resistance") {plt <- plt + xlab("Pyrethroid resistance")}
  ylab_base <- paste0("annual cases averted per ", per_xpop)
  if (only3_comparison) {plt <- plt + ylab(paste0("Additional mean ", ylab_base))}
  if (!only3_comparison) {plt <- plt + ylab(paste0("Mean ", ylab_base))}
  
  plt <- plt + facet_grid(cols = vars(net_name),
                          rows = vars(country),
                          scales="free_y")

  print(plt)
  
}