generate_retention_plt <- function(
    iso2 = NULL,
    force_dhs_adm = FALSE,
    plot_rank = TRUE
) {
  
  # Country-level summary for horizontal reference lines
  country_df <- country_weighted_summary %>%
    dplyr::filter(
      ISO2 == iso2,
      variable %in% c("ret_u", "ret_a")
    )
  
  # Region-level data
  regional_df <- region_summary %>%
    dplyr::filter(
      ISO2 == iso2,
      variable %in% c("ret_u", "ret_a")
    )
  
  # Optionally rename regions to DHS ADM labels
  if (force_dhs_adm) {
    for (i in 1:nrow(regional_df)) {
      regional_df$fs_name_1[i] <- fs_id_link$ADM1[
        which(fs_id_link$fs_area_id == regional_df$fs_area_id[i])
      ]
    }
    regional_df <- regional_df %>%
      dplyr::distinct(fs_name_1, urbanicity, mean_val, .keep_all = TRUE)
  }
  
  # Rank regions
  regional_df <- regional_df %>%
    dplyr::mutate(region_rank = match(fs_name_1, sort(unique(fs_name_1))))
  
  # Choose x-axis: rank or region name
  if (plot_rank) {
    regional_df <- regional_df %>%
      dplyr::mutate(x_region = as.factor(region_rank))
  } else {
    regional_df <- regional_df %>%
      dplyr::mutate(x_region = fs_name_1)
  }
  
  # Assign fill only for ret_a, NA for ret_u
  regional_df <- regional_df %>%
    dplyr::mutate(
      fill_val = ifelse(variable == "ret_a", fs_name_1, NA_character_),
      fill_val = factor(fill_val, levels = sort(unique(na.omit(fs_name_1))))
    )
  
  # Y-axis limits
  ret_axis_lims <- c(
    floor(min(regional_df$lwr, na.rm = TRUE)),
    ceiling(max(regional_df$upr, na.rm = TRUE))
  )
  
  # Compute fill colors, if needed
  fill_levels <- levels(regional_df$fill_val)
  fill_colors <- if (length(fill_levels) > 0) {
    scales::hue_pal()(length(fill_levels))
  } else {
    NULL
  }
  
  # Generate plot
  ggplot(
    regional_df,
    aes(
      x = x_region,
      y = mean_val,
      ymin = lwr,
      ymax = upr,
      shape = urbanicity,
      colour = x_region,
      fill = fill_val
    )
  ) +
    geom_pointrange(
      position = position_dodge(width = 0.5),
      size = 0.5,
      stroke = 1
    ) +
    geom_hline(
      data = country_df,
      aes(yintercept = mean_weighted, linetype = variable),
      colour = "black",
      show.legend = FALSE
    ) +
    scale_shape_manual(values = c("urban" = 24, "rural" = 21)) +
    scale_fill_manual(
      values = fill_colors,
      na.value = "transparent"
    ) +
    scale_linetype_manual(
      values = c("ret_a" = "solid", "ret_u" = "dashed")
    ) +
    scale_y_continuous(
      breaks = seq(0, 48, 3),
      limits = ret_axis_lims
    ) +
    labs(
      x = "Region",
      y = "Mean duration of use / retention (months)"
    ) +
    theme_bw() +
    theme(
      #axis.title.x = element_blank(),
      legend.position = "none"
    )
}


BF_retention <- generate_retention_plt("BF")
GH_retention <- generate_retention_plt("GH")
ML_retention <- generate_retention_plt("ML")
MW_retention <- generate_retention_plt("MW", force_dhs_adm = TRUE)
MZ_retention <- generate_retention_plt("MZ")
SN_retention <- generate_retention_plt("SN")

