

make_clin_cases_plot <- function(
    all_data,
    comparison_name,
    age_group_label,
    x = "EIR_urep_fit",
    x_lo = NULL,
    x_hi = NULL,
    size_range = c(0.2, 3),
    alpha_val = 0.3,
    errorbar_size = 0.5,
    x_breaks = NULL,
    x_label = NULL,
    x_limits = NULL,
    y_limits = NULL,
    y_breaks = NULL,
    log_x = FALSE,
    xsf = 1,
    flip_y = TRUE,
    rename_pyrrole = TRUE,
    roman_labels = FALSE,
    reference_lines = c(12.5, 25, 50, 100),
    iso2 = NULL,
    size_by_nets_dist = FALSE
) {
  
  # Helper function to clip a segment to x/y limits
  clip_segment <- function(m, sign, x_min, x_max, y_min, y_max) {
    # Start with full range
    x1 <- x_min
    x2 <- x_max
    y1 <- sign * m * x1
    y2 <- sign * m * x2
    
    # If y2 is outside bounds, clip x2
    if (y2 > y_max) {
      x2 <- y_max / (sign * m)
      y2 <- y_max
    } else if (y2 < y_min) {
      x2 <- y_min / (sign * m)
      y2 <- y_min
    }
    
    data.frame(x1, x2, y1, y2, slope = m, sign = sign)
  }
  
  if (grepl("^clin_cases_all_ages_3yr(pbo|pyrrole)_med$", x)) {
    xtick_angle <- 90
    xtick_vjust <- 0.5
  } else {
    xtick_angle <- 0
    xtick_vjust <- 0
  }
  
  if (!is.null(iso2)) {all_data %<>% dplyr::filter(ISO2 == iso2)}
  
  df <- all_data %>%
    dplyr::filter(net_strategy_int != "no future nets") %>%
    dplyr::select(
      country_name, fs_name_1, urbanicity, fs_area, fs_area_id,
      EIR_urep_fit, pop, pop_weight, net_strategy_comp, net_strategy_int,
      net_name_int, mass_int_yr_int, pyrethroid_resistance,
      ret_u_int_med, ret_u_int_lo, ret_u_int_hi, ret_a_int_med,
      ret_a_int_lo, ret_a_int_hi, mean_u_3yr_med, mean_u_3yr_lo,
      mean_u_3yr_hi, mean_a_3yr_med, mean_a_3yr_lo, mean_a_3yr_hi,
      mean_uga_3yr_med, mean_uga_3yr_lo, mean_uga_3yr_hi,
      adj_ann_total_nets_dist_int_med,
      clin_cases_all_ages_avert_med, clin_cases_all_ages_avert_lo,
      clin_cases_all_ages_avert_hi, clin_cases_under5_p100kU5_avert_med,
      clin_cases_under5_p100kU5_avert_lo, clin_cases_under5_p100kU5_avert_hi,
      clin_cases_all_ages_3yrpbo_med, clin_cases_all_ages_3yrpbo_lo,
      clin_cases_all_ages_3yrpbo_hi, clin_cases_under5_p100kU5_3yrpbo_med,
      clin_cases_under5_p100kU5_3yrpbo_lo,
      clin_cases_under5_p100kU5_3yrpbo_hi,
      pfpr_182_1824_mean_3yrpbo_med, pfpr_182_1824_mean_3yrpbo_lo,
      pfpr_182_1824_mean_3yrpbo_hi, pfpr_730_3649_mean_3yrpbo_med,
      pfpr_730_3649_mean_3yrpbo_lo, pfpr_730_3649_mean_3yrpbo_hi,
      pfpr_0_36499_mean_3yrpbo_med, pfpr_0_36499_mean_3yrpbo_lo,
      pfpr_0_36499_mean_3yrpbo_hi,
      clin_cases_all_ages_3yrpyrrole_med,
      clin_cases_all_ages_3yrpyrrole_lo,
      clin_cases_all_ages_3yrpyrrole_hi,
      clin_cases_under5_p100kU5_3yrpyrrole_med,
      clin_cases_under5_p100kU5_3yrpyrrole_lo,
      clin_cases_under5_p100kU5_3yrpyrrole_hi,
      pfpr_182_1824_mean_3yrpyrrole_med, pfpr_182_1824_mean_3yrpyrrole_lo,
      pfpr_182_1824_mean_3yrpyrrole_hi, pfpr_730_3649_mean_3yrpyrrole_med,
      pfpr_730_3649_mean_3yrpyrrole_lo, pfpr_730_3649_mean_3yrpyrrole_hi,
      pfpr_0_36499_mean_3yrpyrrole_med, pfpr_0_36499_mean_3yrpyrrole_lo,
      pfpr_0_36499_mean_3yrpyrrole_hi
    ) %>%
    # pivot_longer(
    #   cols = starts_with("clin_cases"),
    #   names_to = c("age_group", "stat"),
    #   names_pattern = "clin_cases_(.*)_avert_(.*)",
    #   values_to = "value"
    # ) %>%
    pivot_longer(
      cols = starts_with("clin_cases") & 
        !matches("3yrpbo") & 
        !matches("3yrpyrrole"),
      names_to = c("age_group", "stat"),
      names_pattern = "clin_cases_(.*)_avert_(.*)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    dplyr::mutate(
      age_group = recode(
        age_group,
        "all_ages" = "all-age",
        "under5_p100kU5" = "under 5"
      ),
      age_group = factor(age_group, levels = c("under 5", "all-age")),
      dist_strategy_int = case_when(
        is.na(mass_int_yr_int) ~ "Continuous only",
        mass_int_yr_int == 2 ~ "2-year campaigns\nwith continuous",
        mass_int_yr_int == 3 ~ "3-year campaigns\nwith continuous"
      ),
      dist_strategy_int = factor(
        dist_strategy_int,
        levels = c(
          "Continuous only",
          "3-year campaigns\nwith continuous",
          "2-year campaigns\nwith continuous"
        )
      )
    )
  
  if (rename_pyrrole) {
    df <- df %>%
      dplyr::mutate(
        net_name_int = dplyr::recode(
          net_name_int,
          "Pyrethroid-Pyrrole" = "Pyrethroid-Chlorfenapyr"
        ),
        net_name_int = factor(
          net_name_int,
          levels = c(
            "Pyrethroid-only",
            "Pyrethroid-PBO",
            "Pyrethroid-Chlorfenapyr"
          )
        )
      )
  } else {
    df <- df %>%
      dplyr::mutate(
        net_name_int = factor(
          net_name_int,
          levels = c(
            "Pyrethroid-only",
            "Pyrethroid-PBO",
            "Pyrethroid-Pyrrole"
          )
        )
      )
  }
  
  
  facet_labels <- df %>%
    dplyr::filter(
      net_strategy_comp == comparison_name,
      age_group == age_group_label
    ) %>%
    dplyr::distinct(net_name_int, dist_strategy_int) %>%
    dplyr::arrange(net_name_int, dist_strategy_int) %>%  # ← Use factor levels
    dplyr::mutate(label = if (roman_labels) {
      tolower(as.character(utils::as.roman(row_number())))
    } else {
      LETTERS[row_number()]
    })
  
  df <- df %>%
    dplyr::left_join(facet_labels, by = c("net_name_int", "dist_strategy_int"))
  
  # Define the variable name used for x
  x_var <- x  # This is already a function argument
  
  # Compute leftward label positions based on panel-specific min(x)
  facet_labels_nudged <- df %>%
    dplyr::filter(
      net_strategy_comp == comparison_name,
      age_group == age_group_label
    ) %>%
    dplyr::group_by(net_name_int, dist_strategy_int) %>%
    dplyr::summarise(
      #x = x_breaks[1] * 1, #min(.data[[x_var]], na.rm = TRUE) * xsf * 0.95,  # nudge right
      y = Inf,
      .groups = "drop"
    ) %>%
    dplyr::mutate(x = x_limits[1]) %>%
    dplyr::left_join(facet_labels, by = c("net_name_int", "dist_strategy_int"))
    
  
  p_df <- df %>%
    dplyr::filter(
      net_strategy_comp == comparison_name,
      age_group == age_group_label
    ) 
  
  if (size_by_nets_dist) {
    p <- p_df %>%
      ggplot(aes(
        x = .data[[x]] * xsf,
        y = if (flip_y) -med / 100 else med / 100,
        ymin = if (flip_y) -hi / 100 else lo / 100,
        ymax = if (flip_y) -lo / 100 else hi / 100,
        colour = if (is.null(iso2)) country_name else fs_name_1,
        fill = if (is.null(iso2)) country_name else fs_name_1,
        shape = urbanicity,
        size = adj_ann_total_nets_dist_int_med / 1e5
      ))
  } else {
    p <- p_df %>%
      ggplot(aes(
        x = .data[[x]] * xsf,
        y = if (flip_y) -med / 100 else med / 100,
        ymin = if (flip_y) -hi / 100 else lo / 100,
        ymax = if (flip_y) -lo / 100 else hi / 100,
        colour = if (is.null(iso2)) country_name else fs_name_1,
        fill = if (is.null(iso2)) country_name else fs_name_1,
        shape = urbanicity
      ))
  }
    
  
  if (
    (comparison_name == "Pyrethroid-Pyrrole 3-year campaigns uncosted" &&
     age_group_label == "all-age" &&
     x == "clin_cases_all_ages_3yrpyrrole_med") |
    (comparison_name == "Pyrethroid-PBO 3-year campaigns uncosted" &&
     age_group_label == "all-age" &&
     x == "clin_cases_all_ages_3yrpbo_med")
  ) {
    # Define slopes and labels
    # 1. Set x/y limits
    x_min <- if (!is.null(x_limits)) x_limits[1] else 0
    x_max <- if (!is.null(x_limits)) x_limits[2] else max(df[[x]] * xsf, na.rm = TRUE)
    
    y_min <- if (!is.null(y_limits)) y_limits[1] else min((if (flip_y) -df$hi else df$lo) / 100, na.rm = TRUE)
    y_max <- if (!is.null(y_limits)) y_limits[2] else max((if (flip_y) -df$lo else df$hi) / 100, na.rm = TRUE)
    
    # 2. Define slopes
    slopes_df <- tibble::tibble(
      pct = reference_lines,
      slope = pct / 100,
      label = dplyr::case_when(
        pct == 0       ~ "0%",
        pct > 100      ~ paste0("+ ", pct, "%"),
        TRUE           ~ paste0("+/- ", pct, "%")
      ),
      show_negative = pct <= 100 & pct != 0
    )
    
    
    # 3. Use clip_segment and build segments
    segments <- purrr::map_dfr(seq_len(nrow(slopes_df)), function(i) {
      s <- slopes_df[i, ]
      pos <- cbind(
        clip_segment(s$slope, 1, x_min, x_max, y_min, y_max),
        label = s$label,
        show_in_legend = TRUE
      )
      if (s$show_negative) {
        neg <- cbind(
          clip_segment(s$slope, -1, x_min, x_max, y_min, y_max),
          label = s$label,
          show_in_legend = FALSE
        )
        rbind(pos, neg)
      } else {
        pos
      }
    })
    segments$.label <- segments$label
    

    
    ordered_labels <- c("+/- 12.5%", "+/- 25%", "+/- 50%", "+/- 100%")
    
    p <- p + geom_segment(
      data = segments,
      aes(
        x = x1, xend = x2,
        y = y1, yend = y2,
        linetype = .label,
        group = .label
      ),
      colour = "grey75",
      linewidth = 0.4,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
    
    # # Define label positions
    # annotation_positions <- segments %>%
    #   dplyr::filter(sign == 1) %>%  # just top-right versions
    #   dplyr::group_by(label) %>%
    #   dplyr::summarise(
    #     x = max(x2, na.rm = TRUE),
    #     y = max(y2, na.rm = TRUE),
    #     .groups = "drop"
    #   )
    
    # Calculate max limits
    x_max <- if (!is.null(x_limits)) x_limits[2] else max(df[[x]] * xsf, na.rm = TRUE)
    y_max <- if (!is.null(y_limits)) y_limits[2] else max((if (flip_y) -df$lo else df$hi) / 100, na.rm = TRUE)
    
    # Define label coordinates for positive slopes
    # annotation_positions <- segments %>%
    #   dplyr::filter(sign == 1) %>%
    #   dplyr::distinct(label, slope) %>%
    #   dplyr::mutate(
    #     # Try placing at x = max(x), compute y = m * x
    #     x = x_max,
    #     y = slope * x,
    #     # If y is too high, compute x = y_max / slope
    #     y = ifelse(y > y_max, y_max, y),
    #     x = ifelse(y == y_max, y_max / slope, x),
    #     # Nudge slightly left for visibility
    #     #x = x * 0.99,
    #     x = x * 0.97,      # ← smaller = more left
    #     #y = y * 1.03,        # ← larger = more up
    #     angle = atan(slope) * 180 / pi
    #   )
    
    annotation_positions <- segments %>%
      dplyr::filter(sign == 1) %>%
      dplyr::distinct(label, slope) %>%
      dplyr::mutate(
        x = x_max,
        y = slope * x,
        y_clipped = y > y_max,
        y = ifelse(y_clipped, y_max - 0.1, y),
        x = ifelse(y_clipped, y_max / slope, x),
        x = x * 0.95,  # leftward nudge
        angle = if (!is.null(iso2)) {
          if (iso2 != "SN") atan(slope) * 180 / pi else 0
        } else {
          atan(slope) * 180 / pi
        }
      )
    
    
    p <- p + geom_text(
      data = annotation_positions,
      aes(x = x, y = y, label = label, angle = angle),
      inherit.aes = FALSE,
      size = 2.5,
      hjust = 1,
      vjust = 0,
      colour = "grey40"
    )
    

    # p <- p + geom_text(
    #   data = annotation_positions,
    #   aes(x = x, y = y, label = label, angle = angle),
    #   inherit.aes = FALSE,
    #   size = 2.5,
    #   hjust = 1.02,
    #   vjust = 0.2,
    #   colour = "grey70"
    # )
    
    
    
    
    # p <- p + scale_linetype_manual(
    #   name = "% change in\nclinical cases",
    #   values = c(
    #     "+/- 12.5%" = "dotdash",
    #     "+/- 25%" = "dotted",
    #     "+/- 50%" = "dashed",
    #     "+/- 100%" = "solid"
    #   ),
    #   breaks = c("+/- 12.5%", "+/- 25%", "+/- 50%", "+/- 100%"),
    #   guide = guide_legend(order = 1)  # Optional: move it to the top
    # )
    
    
  }
  
  
  
  if (!is.null(x_lo) && !is.null(x_hi)) {
    p <- p + geom_errorbarh(aes(
      xmin = .data[[x_lo]] * xsf,
      xmax = .data[[x_hi]] * xsf
    ), size = errorbar_size, alpha = alpha_val)
  }
  
  p <- p +
    geom_errorbar(size = errorbar_size, alpha = alpha_val) +
    geom_point(alpha = alpha_val) +
    scale_size_continuous(range = size_range)
  
  p <- p +
    geom_text(
      data = facet_labels_nudged,
      aes(x = x, y = y, label = label),
      hjust = 0,
      vjust = 2,
      inherit.aes = FALSE,
      fontface = "plain",
      size = 3
    )
  
  
  if (log_x) {
    p <- p + scale_x_log10(
      breaks = x_breaks %||% 10^(0:5),
      minor_breaks = c(
        1:9, unlist(lapply(1:4, function(i) seq(10^i, 10^(i + 1), 10^i)))
      ),
      limits = x_limits
    )
  } else {
    p <- p + scale_x_continuous(
      breaks = x_breaks,
      limits = x_limits
    )
  }
  
  p <- p + scale_y_continuous(
    breaks = y_breaks,
    limits = y_limits
  )
  
  if (flip_y) {y_prefix <- "Change in "} else {y_prefix <- "Additional "}
  if (flip_y) {y_cases_type <- "cases"} else {y_cases_type <- "cases averted"}
  
  p <- p + labs(
    x = x_label %||% x,
    y = if (comparison_name == "no future nets") {
      if (age_group_label == "all-age") {
        paste("Mean annual clinical", y_cases_type, "per 1,000 people")
      } else {
        paste0("Mean annual clinical ", y_cases_type, " per 1,000 (",
               age_group_label, "s)")
      }
    } else {
      if (age_group_label == "all-age") {
        paste(y_prefix, "mean annual clinical", y_cases_type, "per 1,000 people")
      } else {
        paste0(y_prefix,
               "mean annual clinical cases averted per 1,000 (",
               age_group_label, "s)")
      }
    },
    colour = if (is.null(iso2)) "Country" else "Region",
    fill = if (is.null(iso2)) "Country" else "Region",
    size = "Mean annual\nITNs distributed\nper person",
    shape = "Urbanicity"
  ) +
    theme_bw() +
    theme(
      panel.spacing = unit(0, "lines"),
      strip.background = element_blank(),
      strip.text = element_text(size = 10),
      axis.text.x = element_text(angle = xtick_angle, vjust = xtick_vjust)
    ) +
    facet_grid(
      rows = vars(net_name_int),
      cols = vars(dist_strategy_int)
    )
  
  if (
    (comparison_name == "Pyrethroid-Pyrrole 3-year campaigns uncosted" &&
     age_group_label == "all-age" &&
     x == "clin_cases_all_ages_3yrpyrrole_med") |
    (comparison_name == "Pyrethroid-PBO 3-year campaigns uncosted" &&
     age_group_label == "all-age" &&
     x == "clin_cases_all_ages_3yrpbo_med")
  ) {
    if (!is.null(iso2)) {if (iso2 != "SN") p <- p + coord_fixed(ratio = 1)}
  }
  
  return(p)

}
