uga_ret_quadrant_plt <- function(iso2 = "SN",
                                 per_xpop = 1,
                                 plot_quadrants = TRUE,
                                 annot_labels = NULL,
                                 rural_labels = NULL,
                                 urban_labels = NULL,
                                 force_dhs_adm = FALSE,
                                 horiz_offset = 0,
                                 flip_xy = TRUE,
                                 retu_limits = TRUE) {
  
  sim_sum_long <- region_summary %>%
    dplyr::filter(ISO2 == iso2) %>%
    dplyr::filter(
      variable == "ret_a" | variable == "mean_uga" | variable == "ret_u")
  
  sim_sum <- sim_sum_long %>%
    pivot_wider(
      id_cols = c(ISO2, fs_name_1, urbanicity, fs_area, fs_area_id, new_area_id, pop),
      names_from = variable,
      values_from = c(mean_val, lwr, upr),
      names_glue = "{variable}_{.value}"
    ) %>%
    rename(
      uga_med  = mean_uga_mean_val,
      uga_lo   = mean_uga_lwr,
      uga_hi   = mean_uga_upr,
      reta_med = ret_a_mean_val,
      reta_lo  = ret_a_lwr,
      reta_hi  = ret_a_upr,
      retu_med = ret_u_mean_val,
      retu_lo  = ret_u_lwr,
      retu_hi  = ret_u_upr
    )
  
  if (force_dhs_adm) {
    for (i in 1:nrow(sim_sum)) {
      sim_sum$fs_name_1[i] <- fs_id_link$ADM1[which(fs_id_link$fs_area_id == sim_sum$fs_area_id[i])]
    }
    sim_sum %<>% dplyr::distinct(fs_name_1, urbanicity, .keep_all = TRUE)
  }
  
  sim_sum %<>%
    group_by(fs_name_1) %>%
    mutate(group_id = cur_group_id())
  
  if (is.null(annot_labels)) {
    sim_sum$inc_group_id <- sim_sum$group_id
    sim_sum$inc_group_id[sim_sum$urbanicity == "rural" & sim_sum$fs_name_1 %in% rural_labels] <- ""
    sim_sum$inc_group_id[sim_sum$urbanicity == "urban" & sim_sum$fs_name_1 %in% urban_labels] <- ""
    
    sim_sum$exc_group_id <- sim_sum$group_id
    sim_sum$exc_group_id[sim_sum$urbanicity == "rural" & !(sim_sum$fs_name_1 %in% rural_labels)] <- ""
    sim_sum$exc_group_id[sim_sum$urbanicity == "urban" & !(sim_sum$fs_name_1 %in% urban_labels)] <- ""
  } else {
    annot_labels <- annot_labels[annot_labels[,1] == country,]
    N_annot <- nrow(annot_labels)
    if (is.null(N_annot)) {
      annot_labels <- t(as.matrix(annot_labels))
      N_annot <- 1
    }
    annot_ids <- rep(NA, N_annot)
    for (i in 1:N_annot) {
      annot_ids[i] <- which(
        sim_sum$ISO2 == annot_labels[i,1] &
          sim_sum$mass_int == annot_labels[i,2] &
          sim_sum$net_name == annot_labels[i,3] &
          sim_sum$urbanicity == annot_labels[i,4] &
          sim_sum$fs_name_1 == annot_labels[i,5]
      )
    }
    sim_sum$inc_group_id <- sim_sum$group_id
    sim_sum$inc_group_id[annot_ids] <- ""
    
    sim_sum$exc_group_id <- rep("", length(sim_sum$group_id))
    sim_sum$exc_group_id[annot_ids] <- sim_sum$group_id[annot_ids]
  }
  
  # X-axis (uret)
  sim_sum$xmean <- sim_sum$reta_med 
  sim_sum$xLB <- sim_sum$reta_lo
  sim_sum$xUB <- sim_sum$reta_hi
  xround <- 1.5
  xbreaks <- seq(0, 48, 3)
  if (retu_limits) {
    # horiz_min <- floor(min(sim_sum$ret_u_lwr)/xround)*xround
    # horiz_max <- ceiling(max(sim_sum$ret_u_upr)/xround)*xround
    horiz_min <- floor(min(sim_sum$retu_lo))
    horiz_max <- ceiling(max(sim_sum$reta_hi))
  } else {
    horiz_min <- floor(min(sim_sum$xLB)/xround)*xround
    horiz_max <- ceiling(max(sim_sum$xUB)/xround)*xround
  }
  print(min(sim_sum$retu_lo))
  print(max(sim_sum$reta_hi))
  print(min(horiz_min))
  print(max(horiz_max))
  horiz_mid <- round(mean(sim_sum$xmean)/xround)*xround
  horiz_delta <- max(horiz_mid - horiz_min, horiz_max - horiz_mid)
  #horiz_min <- horiz_mid - horiz_delta
  #horiz_max <- horiz_mid + horiz_delta
  horiz_high <- horiz_max - 0.1 * horiz_delta
  horiz_low <- horiz_min + 0.1 * horiz_delta
  
  # Y-axis (uga)
  sim_sum$ymean <- sim_sum$uga_med * 100
  sim_sum$yLB <- sim_sum$uga_lo * 100
  sim_sum$yUB <- sim_sum$uga_hi * 100
  sim_sum$yUB[sim_sum$yUB > 100] <- 100
  vert_min <- floor(min(sim_sum$yLB) / 2.5) * 2.5
  vert_max <- ceiling(max(sim_sum$yUB) / 2.5) * 2.5
  vert_mid <- round(mean(sim_sum$ymean) / 2.5) * 2.5
  vert_delta <- max(vert_mid - vert_min, vert_max - vert_mid)
  vert_min <- vert_mid - vert_delta
  vert_max <- vert_mid + vert_delta
  vert_high <- vert_max
  vert_low <- vert_min
  ybreaks <- seq(0, 100, 10)
  
  #horiz_mean <- sum(sim_sum$xmean * sim_sum$pop) / sum(sim_sum$pop)
  #vert_mean <- sum(sim_sum$ymean * sim_sum$pop) / sum(sim_sum$pop)
  horiz_mean <- country_weighted_summary$mean_weighted[
    which(
      country_weighted_summary$ISO2 == iso2 &
        country_weighted_summary$variable == "ret_a"
    )
  ]
  vert_mean <- country_weighted_summary$mean_weighted[
    which(
      country_weighted_summary$ISO2 == iso2 &
        country_weighted_summary$variable == "mean_uga"
    )
  ] * 100
  
  sim_sum$numbered_areas <- paste0(sim_sum$fs_name_1, " (", sim_sum$group_id, ")")
  pos <- if (force_dhs_adm) {
    position_dodge(width = 0)
  } else {
    position_dodge(width = 0.25)
  }
  
  plt <- ggplot(sim_sum,
                aes(x = xmean,
                    y = ymean,
                    shape = urbanicity,
                    colour = numbered_areas,
                    label = as.factor(inc_group_id),
                    group = interaction(urbanicity, numbered_areas)))
  
  if (plot_quadrants) {
    plt <- plt + 
      geom_hline(yintercept = vert_mean, colour = "black") +
      geom_vline(xintercept = horiz_mean, colour = "black")
  }
  
  plt <- plt +
    #geom_pointrange(aes(xmin = xLB, xmax = xUB), size = 1.3, position = pos, alpha = 0.4) +
    #geom_pointrange(aes(ymin = yLB, ymax = yUB), size = 1.3, position = pos, alpha = 0.4) +
    geom_errorbarh(aes(xmin = xLB, xmax = xUB), height = 0, position = pos,
                   alpha = 0.6) +
    geom_errorbar(aes(ymin = yLB, ymax = yUB), width = 0, position = pos,
                  alpha = 0.6) +
    geom_point(position = pos, size = 4, alpha = 0.6) +
    scale_x_continuous("Mean ITN retention time (months)", breaks = xbreaks, labels = xbreaks, limits = c(horiz_min, horiz_max)) +
    scale_y_continuous("Mean use given access (%)", breaks = ybreaks, labels = ybreaks, limits = c(vert_min, vert_max)) +
    geom_text(colour = 'black', size = 2.5, show.legend = FALSE, position = pos)
  
  if (!(is.null(rural_labels) & is.null(urban_labels) & is.null(annot_labels))) {
    plt <- plt +
      geom_text_repel(aes(label = exc_group_id), colour = 'black', size = 2.5, position = pos,
                      segment.colour = NA, max.overlaps = 30)
  }
  
  if (plot_quadrants) {
    if (flip_xy) {
      plt <- plt +
        annotate("text",
                 y = vert_high - (vert_high - vert_low) * horiz_offset,
                 x = horiz_high, label = "Quadrant 1") +
        annotate("text",
                 y = vert_low + (vert_high - vert_low) * horiz_offset,
                 x = horiz_high, label = "Quadrant 2") +
        annotate("text",
                 y = vert_low + (vert_high - vert_low) * horiz_offset,
                 x = horiz_low,  label = "Quadrant 3") +
        annotate("text",
                 y = vert_high - (vert_high - vert_low) * horiz_offset,
                 x = horiz_low,  label = "Quadrant 4")
    } else {
      plt <- plt +
        annotate("text",
                 x = horiz_high - (horiz_high - horiz_low) * horiz_offset,
                 y = vert_high, label = "Quadrant 1") +
        annotate("text",
                 x = horiz_low + (horiz_high - horiz_low) * horiz_offset,
                 y = vert_high, label = "Quadrant 2") +
        annotate("text",
                 x = horiz_low + (horiz_high - horiz_low) * horiz_offset,
                 y = vert_low,  label = "Quadrant 3") +
        annotate("text",
                 x = horiz_high - (horiz_high - horiz_low) * horiz_offset,
                 y = vert_low,  label = "Quadrant 4")
    }
    
  }
  
  plt <- plt +
    guides(colour = guide_legend(title = "Region"),
           shape = guide_legend(title = "Urbanicity"),
           label = "none") +
    theme_bw() +
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", color = "transparent"),
      legend.box.background = element_rect(fill = "transparent", color = "transparent"),
      legend.key = element_rect(fill = "transparent"),
      legend.text = element_text(size = 8)
    )
  
  if (flip_xy) {
    plt <- plt + coord_flip()
  }
  

  plt <- plt +
    scale_x_continuous("Mean ITN retention time (months)", breaks = xbreaks, labels = xbreaks, limits = c(horiz_min, horiz_max))

  ggsave(paste(iso2, "uga", "uret", "quadrantfinal.pdf", sep = "_"),
         plot = plt, bg = "transparent", w = 8, h = 6, dpi = 450)
  
  return(plt)
}

library(ggrepel)

BF_quad <- uga_ret_quadrant_plt(iso2 = "BF",
                                horiz_offset = 0.15)

GH_quad <- uga_ret_quadrant_plt(
  iso2 = "GH",
  horiz_offset = 0.1,
  rural_labels = c(
    "Northern", "Savannah", "Western", "Oti", "Bono East", "Ahafo"
  ),
  urban_labels = c(
    "North East", "Savannah", "Western", "Oti", "Bono East", "Ahafo"
  )
)

ML_quad <- uga_ret_quadrant_plt(iso2 = "ML",
                                horiz_offset = 0.1,
                                plot_quadrants = TRUE)

MW_quad <- uga_ret_quadrant_plt(iso2 = "MW",
                                horiz_offset = 0.1, force_dhs_adm = TRUE)

MZ_quad <- uga_ret_quadrant_plt(iso2 = "MZ",
                                horiz_offset = 0.1)

SN_quad <- uga_ret_quadrant_plt(iso2 = "SN",
                                horiz_offset = 0.1,
                                urban_labels = "Fatick",
                     rural_labels = "Kaffrine")



