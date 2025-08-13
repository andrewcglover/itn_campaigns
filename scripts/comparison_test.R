comparison_name <- "no future nets"
plt_eir_all_age <- eir_df %>%
  dplyr::filter(
    net_strategy_comp == comparison_name
  ) %>%
  dplyr::filter(age_group == "all-age") %>%
  ggplot(
    aes(
      x = EIR_urep_fit,
      y = med/100,
      ymin = lo/100,
      ymax = hi/100,
      colour = country_name,
      size = adj_ann_total_nets_dist_int_med / 1e5,
      fill = country_name,
      shape = urbanicity
    )
  ) +
  geom_errorbar(size = errorbar_size, alpha = alpha_val) +
  #geom_point(shape = 16, alpha = 0.4) +
  geom_point(alpha = alpha_val) +
  scale_size(range = size_range) +
  #scale_x_log10() +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000),
    minor_breaks =c(
      seq(1, 10, by = 1),
      seq(10, 100, by = 10),
      seq(100, 1000, by = 100)
    )
  ) +
  labs(
    x = expression("Baseline " * italic("Pf") * "EIR"),
    y = "Annual mean clinical cases averted per 1,000 (under 5's)",
    colour = "Country",
    fill = "Country",
    size = "Annual mean\nITNs distributed\nper person",
    shape = "Urbanicity"
  ) +
  theme_bw() +
  theme(
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10)  # optional: tweak facet label text
  ) +
  facet_grid(rows = vars(net_name_int), cols = vars(dist_strategy_int))


comparison_name <- "no future nets"
plt_eir_u5 <- eir_df %>%
  dplyr::filter(
    net_strategy_comp == comparison_name
  ) %>%
  dplyr::filter(age_group == "under 5") %>%
  ggplot(
    aes(
      x = EIR_urep_fit,
      y = med/100,
      ymin = lo/100,
      ymax = hi/100,
      colour = country_name,
      size = adj_ann_total_nets_dist_int_med / 1e5,
      fill = country_name,
      shape = urbanicity
    )
  ) +
  geom_errorbar(size = errorbar_size, alpha = alpha_val) +
  #geom_point(shape = 16, alpha = 0.4) +
  geom_point(alpha = alpha_val) +
  scale_size(range = size_range) +
  #scale_x_log10() +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000),
    minor_breaks =c(
      seq(1, 10, by = 1),
      seq(10, 100, by = 10),
      seq(100, 1000, by = 100)
    )
  ) +
  labs(
    x = expression("Baseline " * italic("Pf") * "EIR"),
    y = "Annual mean clinical cases averted per 1,000 (under 5's)",
    colour = "Country",
    fill = "Country",
    size = "Annual mean\nITNs distributed\nper person",
    shape = "Urbanicity"
  ) +
  theme_bw() +
  theme(
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10)  # optional: tweak facet label text
  ) +
  facet_grid(rows = vars(net_name_int), cols = vars(dist_strategy_int))

comparison_name <- "Pyrethroid-only 3-year campaigns uncosted"
plt_eir_u5_only <- eir_df %>%
  dplyr::filter(
    net_strategy_comp == comparison_name
  ) %>%
  dplyr::filter(age_group == "under 5") %>%
  ggplot(
    aes(
      x = EIR_urep_fit,
      y = med/100,
      ymin = lo/100,
      ymax = hi/100,
      colour = country_name,
      size = adj_ann_total_nets_dist_int_med / 1e5,
      fill = country_name,
      shape = urbanicity
    )
  ) +
  geom_errorbar(size = errorbar_size, alpha = alpha_val) +
  #geom_point(shape = 16, alpha = 0.4) +
  geom_point(alpha = alpha_val) +
  scale_size(range = size_range) +
  #scale_x_log10() +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000),
    minor_breaks =c(
      seq(1, 10, by = 1),
      seq(10, 100, by = 10),
      seq(100, 1000, by = 100)
    )
  ) +
  labs(
    x = expression("Baseline " * italic("Pf") * "EIR"),
    y = "Annual mean clinical cases averted per 1,000 (under 5's)",
    colour = "Country",
    fill = "Country",
    size = "Annual mean\nITNs distributed\nper person",
    shape = "Urbanicity"
  ) +
  theme_bw() +
  theme(
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10)  # optional: tweak facet label text
  ) +
  facet_grid(rows = vars(net_name_int), cols = vars(dist_strategy_int))


comparison_name <- "Pyrethroid-PBO 3-year campaigns uncosted"
plt_eir_u5_pbo <- eir_df %>%
  dplyr::filter(
    net_strategy_comp == comparison_name
  ) %>%
  dplyr::filter(age_group == "under 5") %>%
  ggplot(
    aes(
      x = EIR_urep_fit,
      y = med/100,
      ymin = lo/100,
      ymax = hi/100,
      colour = country_name,
      size = adj_ann_total_nets_dist_int_med / 1e5,
      fill = country_name,
      shape = urbanicity
    )
  ) +
  geom_errorbar(size = errorbar_size, alpha = alpha_val) +
  #geom_point(shape = 16, alpha = 0.4) +
  geom_point(alpha = alpha_val) +
  scale_size(range = size_range) +
  #scale_x_log10() +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000),
    minor_breaks =c(
      seq(1, 10, by = 1),
      seq(10, 100, by = 10),
      seq(100, 1000, by = 100)
    )
  ) +
  labs(
    x = expression("Baseline " * italic("Pf") * "EIR"),
    y = "Annual mean clinical cases averted per 1,000 (under 5's)",
    colour = "Country",
    fill = "Country",
    size = "Annual mean\nITNs distributed\nper person",
    shape = "Urbanicity"
  ) +
  theme_bw() +
  theme(
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10)  # optional: tweak facet label text
  ) +
  facet_grid(rows = vars(net_name_int), cols = vars(dist_strategy_int))


comparison_name <- "Pyrethroid-Pyrrole 3-year campaigns uncosted"
plt_eir_u5_pyrrole <- eir_df %>%
  dplyr::filter(
    net_strategy_comp == comparison_name
  ) %>%
  dplyr::filter(age_group == "under 5") %>%
  ggplot(
    aes(
      x = EIR_urep_fit,
      y = med/100,
      ymin = lo/100,
      ymax = hi/100,
      colour = country_name,
      size = adj_ann_total_nets_dist_int_med / 1e5,
      fill = country_name,
      shape = urbanicity
    )
  ) +
  geom_errorbar(size = errorbar_size, alpha = alpha_val) +
  #geom_point(shape = 16, alpha = 0.4) +
  geom_point(alpha = alpha_val) +
  scale_size(range = size_range) +
  #scale_x_log10() +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000),
    minor_breaks =c(
      seq(1, 10, by = 1),
      seq(10, 100, by = 10),
      seq(100, 1000, by = 100)
    )
  ) +
  labs(
    x = expression("Baseline " * italic("Pf") * "EIR"),
    y = "Annual mean clinical cases averted per 1,000 (under 5's)",
    colour = "Country",
    fill = "Country",
    size = "Annual mean\nITNs distributed\nper person",
    shape = "Urbanicity"
  ) +
  theme_bw() +
  theme(
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10)  # optional: tweak facet label text
  ) +
  facet_grid(rows = vars(net_name_int), cols = vars(dist_strategy_int))

