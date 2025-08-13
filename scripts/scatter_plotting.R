cases_averted_scatter <- function(sim_sum,
                                  country = NULL,
                                  rm.country = NULL,
                                  var_name = "pfeir",
                                  per_xpop = 1,
                                  only3_comparison = FALSE){
  
  sim_sum$country <- countrycode(sim_sum$ISO2,"iso2c","country.name")
  
  if (!is.null(country)) {sim_sum %<>% filter(ISO2 == country)}
  if (!is.null(rm.country)) {sim_sum %<>% filter(ISO2 != rm.country)}
  
  sim_sum$net_strategy <- paste(sim_sum$net_name,
                                sim_sum$mass_int,
                                "year interval")
  
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
  if (var_name == "uret") {sim_sum$var_LB <- sim_sum$LB_retu}
  if (var_name == "uret") {sim_sum$var_UB <- sim_sum$UB_retu}
  if (var_name == "aret") {sim_sum$var_val <- sim_sum$mean_reta}
  if (var_name == "aret") {sim_sum$var_LB <- sim_sum$LB_reta}
  if (var_name == "aret") {sim_sum$var_UB <- sim_sum$UB_reta}
  if (var_name == "uga") {sim_sum$var_val <- sim_sum$mean_uga}
  if (var_name == "uga") {sim_sum$var_LB <- sim_sum$LB_uga}
  if (var_name == "uga") {sim_sum$var_UB <- sim_sum$UB_uga}

  
  # sim_sum$var_LB <- sim_sum$var_val
  # sim_sum$var_UB <- sim_sum$var_val
  
  #sim_sum$LB_avert_percap[which(sim_sum$LB_avert_percap < 0.01)] <- 0.01
  

  
  #pos <- position_dodge(width = 0.5)
  alpha_val <- 0.5
  size_val <- 0.4
  label_vals <- c("2 year interval", "3 year interval")
  if (var_name == "uret" || var_name == "aret" || var_name == "uga") {
    plt <- ggplot() +
      geom_linerange(
        data = sim_sum %>% filter(net_name == "pyrethroid-only"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          xmin = var_LB,
          xmax = var_UB,
          color = net_strategy,
          shape = urbanicity
        ),
        #position = pos,
        alpha = alpha_val,
        #size = size_val,
        stroke=NA
      ) +
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
        #position = pos,
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
        guide = guide_legend(title = "Pyrethroid-only",
                             order = 1,
                             ncol = 1,
                             title.position = "top"
                             ),
        labels = label_vals
      ) +
      new_scale_colour() +
      geom_linerange(
        data = sim_sum %>% filter(net_name == "pyrethroid-PBO"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          xmin = var_LB,
          xmax = var_UB,
          color = net_strategy,
          shape = urbanicity
        ),
        #position = pos,
        alpha = alpha_val,
        size = size_val,
        stroke=NA
      ) +
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
        #position = pos,
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
        guide = guide_legend(title = "Pyrethroid-PBO",
                             order = 2,
                             ncol = 1,
                             title.position = "top"
                             ),
        labels = label_vals
      ) +
      new_scale_colour() +
      geom_linerange(
        data = sim_sum %>% filter(net_name == "pyrethroid-pyrrole"),
        aes(
          x = var_val,
          y = mean_cases_avert_plt,
          xmin = var_LB,
          xmax = var_UB,
          color = net_strategy,
          shape = urbanicity
        ),
        #position = pos,
        alpha = alpha_val,
        size = size_val,
        stroke=NA
      ) +
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
        #position = pos,
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
        guide = guide_legend(title = "Pyrethroid-pyrrole",
                             order = 3,
                             ncol = 1,
                             title.position = "top"
                             ),
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
        #position = pos,
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
        guide = guide_legend(title = "Pyrethroid-only",
                             order = 1,
                             ncol = 1,
                             title.position = "top"
                             ),
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
        #position = pos,
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
        guide = guide_legend(title = "Pyrethroid-PBO",
                             order = 2,
                             ncol = 1,
                             title.position = "top"
                             ),
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
        #position = pos,
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
        guide = guide_legend(title = "Pyrethroid-pyrrole",
                             order = 3,
                             ncol = 1,
                             title.position = "top"
                             ),
        labels = label_vals
        )
  }
  
  log_breaks <- 10^(-10:10)
  log_minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
  
  plt <- plt +
    scale_y_continuous(limits = c(0, NA))
  
  #if (var_name == "pfeir") {plt <- plt + scale_x_continuous(trans='log10')}
  if (var_name == "pfeir") {plt <- plt + scale_x_log10(breaks = log_breaks,
                                                       minor_breaks = log_minor_breaks)}
  if (var_name == "pfeir") {plt <- plt + xlab(expression(~italic(Pf)~'EIR'))}
  if (var_name == "pyrethroid_resistance") {plt <- plt + xlab("Pyrethroid resistance")}
  if (var_name == "uret") {plt <- plt + xlab("Mean duration of use (months)")}
  if (var_name == "aret") {plt <- plt + xlab("Mean duration of access (months)")}
  if (var_name == "uga") {plt <- plt + xlab("Mean usage given access")}
  if (var_name == "pyrethroid_resistance") {plt <- plt + xlab("Pyrethroid resistance")}
  if (per_xpop == 1) {per_xpop <- "capita"}
  ylab_base <- paste0("annual cases averted per ", per_xpop)
  if (only3_comparison) {plt <- plt + ylab(paste0("Additional mean ", ylab_base))}
  if (!only3_comparison) {plt <- plt + ylab(paste0("Mean ", ylab_base))}
  
  plt <- plt +
    theme_bw() +
    theme(
      #panel.background = element_rect(fill = "transparent",
      #                                colour = NA_character_), # necessary to avoid drawing panel outline
      #panel.grid.major = element_blank(), # get rid of major grid
      #panel.grid.minor = element_blank(), # get rid of minor grid
      plot.background = element_rect(fill = "transparent",
                                     colour = NA_character_), # necessary to avoid drawing plot outline
      #legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
      # legend.background = element_rect(fill = "transparent"),
      #legend.background = element_rect(fill = "transparent",colour = "black", linetype='solid'),
      #legend.box.background = element_rect(fill = "transparent"),
      legend.key = element_rect(fill = "transparent")
    ) +
    facet_grid(cols = vars(country),
               rows = vars(net_name),
               #scales="free_y"
    ) +
    theme(legend.position="bottom") +
    guides(shape = guide_legend(title = "Urbanicity",
                                order = 4,
                                ncol = 1,
                                title.position = "top"))
    # facet_grid(cols = vars(net_name),
    #                       rows = vars(country),
    #                       #scales="free_y"
    #            )

  print(plt)
  
  if(only3_comparison) {add_str <- "add"} else {add_str <- ""}
  # ggsave(paste0("allctry_", add_str, "caseavrt_vs_",var_name,"_horiz.pdf"), bg = "transparent",
  #        w = 6, h = 8, dpi = 450)
  ggsave(paste0("allctry_", add_str, "caseavrt_vs_",var_name,"_horiz_18NOV24.pdf"), bg = "transparent",
         w = 8, h = 5.5, dpi = 450)
  
}

quadrant_mean_retention <- function(sim_sum,
                                    country = "SN",
                                    xvar = "uret",
                                    yvar = "uga",
                                    per_xpop = 1,
                                    plot_quadrants = TRUE,
                                    facets_on = FALSE,
                                    annot_labels = NULL,
                                    rural_labels = NULL,
                                    urban_labels = NULL,
                                    force_dhs_adm = FALSE){
  
  sim_sum %<>% filter(ISO2 == country)
  if (facets_on) {
    if (yvar == "add_avert") {
      sim_sum %<>% dplyr::filter(net_strategy == "uncosted" & !(mass_int == 3 & net_name == "pyrethroid-only"))
    }
  } else {
    sim_sum %<>% dplyr::filter(net_strategy == "uncosted" & mass_int == 3 & net_name == "pyrethroid-only")
  }
  
  sim_sum$mass_int_name <- paste0(sim_sum$mass_int, "-year interval")
  
  if (force_dhs_adm) {
    for (i in 1:dim(sim_sum)[1]) {
      sim_sum$fs_name_1[i] <- fs_id_link$ADM1[which(
        fs_id_link$fs_area_id == sim_sum$fs_area_id[i])]
    }
  }
  
  # Fill in missing EIR for Kedougou, Senegal
  if (country == "SN") {
    sim_sum$pfeir[sim_sum$pfeir == 0] <- max(sim_sum$pfeir)
  }
  
  sim_sum %<>%
    dplyr::group_by(fs_name_1) %>%
    dplyr::mutate(group_id = dplyr::cur_group_id())
  
  if (is.null(annot_labels)) {
    
    sim_sum$inc_group_id <- sim_sum$group_id
    sim_sum$inc_group_id[which((sim_sum$urbanicity == "rural") &
                                 (sim_sum$fs_name_1 %in% rural_labels))] <- ""
    sim_sum$inc_group_id[which((sim_sum$urbanicity == "urban") &
                                 (sim_sum$fs_name_1 %in% urban_labels))] <- ""
    
    sim_sum$exc_group_id <- sim_sum$group_id
    sim_sum$exc_group_id[which((sim_sum$urbanicity == "rural") &
                                 !(sim_sum$fs_name_1 %in% rural_labels))] <- ""
    sim_sum$exc_group_id[which((sim_sum$urbanicity == "urban") &
                                 !(sim_sum$fs_name_1 %in% urban_labels))] <- ""
    
  } else {
    
    annot_labels <- annot_labels[annot_labels[,1] == country,]
    N_annot <- dim(annot_labels)[1]
    if (is.null(N_annot)) {
      annot_labels %<>% as.matrix %>% t
      N_annot <- 1
    }
    annot_ids <- rep(NA, N_annot)
    
    for (i in 1:N_annot) {
      #print(annot_labels[i,])
      annot_ids[i] <- which(
        (sim_sum$ISO2 == annot_labels[i,1]) &
          (sim_sum$mass_int == annot_labels[i,2]) &
          (sim_sum$net_name == annot_labels[i,3]) &
          (sim_sum$urbanicity == annot_labels[i,4]) &
          (sim_sum$fs_name_1 == annot_labels[i,5])
      )
    }
    
    sim_sum$inc_group_id <- sim_sum$group_id
    sim_sum$inc_group_id[annot_ids] <- ""
    
    sim_sum$exc_group_id <- rep("", length(sim_sum$group_id))
    sim_sum$exc_group_id[annot_ids] <- sim_sum$group_id[annot_ids]
    
  }
  
  #test <<- sim_sum$exc_group_id
  
  if (xvar == "uret") {
    sim_sum$xmean <- sim_sum$mean_retu
    sim_sum$xLB <- sim_sum$LB_retu
    sim_sum$xUB <- sim_sum$UB_retu
    xround <- 1.5
    xbreaks <- seq(0,48,3)
    horiz_min <- floor(min(sim_sum$xLB)/xround)*xround
    horiz_max <- ceiling(max(sim_sum$xUB)/xround)*xround
    horiz_mid <- round(mean(sim_sum$xmean)/xround)*xround
    horiz_delta <- max(c(horiz_mid - horiz_min, horiz_max - horiz_mid))
    horiz_min <- horiz_mid - horiz_delta
    horiz_max <- horiz_mid + horiz_delta
    horiz_high <- horiz_max - 0.1 * horiz_delta
    horiz_low <- horiz_min + 0.1 * horiz_delta
  } else if (xvar == "aret") {
    sim_sum$xmean <- sim_sum$mean_reta
    sim_sum$xLB <- sim_sum$LB_reta
    sim_sum$xUB <- sim_sum$UB_reta
    xround <- 1.5
    xbreaks <- seq(0,48,3)
    horiz_min <- floor(min(sim_sum$xLB)/xround)*xround
    horiz_max <- ceiling(max(sim_sum$xUB)/xround)*xround
    horiz_mid <- round(mean(sim_sum$xmean)/xround)*xround
    horiz_delta <- max(c(horiz_mid - horiz_min, horiz_max - horiz_mid))
    horiz_min <- horiz_mid - horiz_delta
    horiz_max <- horiz_mid + horiz_delta
    horiz_high <- horiz_max - 0.1 * horiz_delta
    horiz_low <- horiz_min + 0.1 * horiz_delta
  } else if (xvar == "pfeir") {
    sim_sum$xmean <- sim_sum$pfeir
  }
  
  if (yvar == "uga") {
    sim_sum$ymean <- sim_sum$mean_uga*100
    sim_sum$yLB <- sim_sum$LB_uga*100
    sim_sum$yUB <- sim_sum$UB_uga*100
    sim_sum$yUB[sim_sum$yUB > 100] <- 100
    vert_min <- floor(min(sim_sum$yLB)/2.5)*2.5
    vert_max <- ceiling(max(sim_sum$yUB)/2.5)*2.5
    vert_mid <- round(mean(sim_sum$ymean)/2.5)*2.5
    vert_delta <- max(c(vert_mid - vert_min, vert_max - vert_mid))
    vert_min <- vert_mid - vert_delta
    vert_max <- vert_mid + vert_delta
    vert_high <- vert_max - 0 * vert_delta
    vert_low <- vert_min + 0 * vert_delta
    ybreaks <- seq(0,100,10)
  } else if (yvar == "add_avert") {
    sim_sum$ymean <- sim_sum$mean_add_avert_percap*per_xpop
    sim_sum$yLB <- sim_sum$LB_add_avert_percap*per_xpop
    sim_sum$yUB <- sim_sum$UB_add_avert_percap*per_xpop
  } else if (yvar == "avert") {
    sim_sum$ymean <- sim_sum$mean_avert_percap*per_xpop
    sim_sum$yLB <- sim_sum$LB_avert_percap*per_xpop
    sim_sum$yUB <- sim_sum$UB_avert_percap*per_xpop
  }
  
  horiz_mean <- sum(sim_sum$xmean * sim_sum$pop / sum(sim_sum$pop))
  vert_mean <- sum(sim_sum$ymean * sim_sum$pop / sum(sim_sum$pop))
  
  sim_sum$numbered_areas <- paste0(sim_sum$fs_name_1," (",sim_sum$group_id,")")

  alpha_val <- 0.8
  log_breaks <- 10^(-10:10)
  log_minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
  
  if (force_dhs_adm) {
    pos <- position_dodge(width = 0)
  } else {
    pos <- position_dodge(width = 0.5)
  }
  #pos <- position_dodge2(width = 0.3)
  
  plt <- ggplot(data = sim_sum,
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
  if (xvar == "uret") {
    plt <- plt +
      geom_pointrange(aes(xmin = xLB, xmax = xUB),
                      size = 1.3,
                      position = pos,
                      alpha = 0.4) +
      scale_x_continuous("Mean duration of usage (months)",
                         breaks = xbreaks,
                         labels = xbreaks,
                         limits = c(horiz_min, horiz_max))
  } else if (xvar == "aret") {
    plt <- plt +
      geom_pointrange(aes(xmin = xLB, xmax = xUB),
                      size = 1.3,
                      position = pos,
                      alpha = 0.4) +
      scale_x_continuous("Mean duration of access (months)",
                         breaks = xbreaks,
                         labels = xbreaks,
                         limits = c(horiz_min, horiz_max))
  } else if (xvar == "pfeir") {
    plt <- plt +
      #scale_x_continuous("PfEIR",
      #                   trans = "log10")
      scale_x_log10(expression(~italic(Pf)~'EIR'),
                    breaks = log_breaks,
                    minor_breaks = log_minor_breaks)
  }
  if (yvar == "uga") {
    plt <- plt +
      geom_pointrange(aes(ymin = yLB, ymax = yUB),
                      size = 1.3,
                      position = pos,
                      alpha = 0.4) +
      scale_y_continuous("Use given access (%)",
                         breaks = ybreaks,
                         labels = ybreaks,
                         limits = c(vert_min, vert_max))
  } else if (yvar == "add_avert") {
    if (per_xpop == 1) {per_char <- "capita"} else {per_char <- per_xpop}
    ytitle <- paste("Additional mean annual cases averted per", per_char)
    plt <- plt +
      geom_pointrange(aes(ymin = yLB, ymax = yUB),
                      size = 1.3,
                      position = pos,
                      alpha = 0.6) +
      scale_y_continuous(ytitle,
                         #expand = c(0, 0),
                         limits = c(0, NA))
  } else if (yvar == "avert") {
    if (per_xpop == 1) {per_char <- "capita"} else {per_char <- per_xpop}
    ytitle <- paste("Annual mean cases averted per", per_char)
    plt <- plt +
      geom_pointrange(aes(ymin = yLB, ymax = yUB),
                      size = 1.3,
                      position = pos,
                      alpha = 0.6) +
      scale_y_continuous(ytitle,
                         #expand = c(0, 0),
                         limits = c(0, NA))
  }
  plt <- plt +
    geom_text(colour = 'black',
              size = 2.5,
              show.legend = FALSE,
              position = pos)
  if (!(is.null(rural_labels) & is.null(urban_labels) & is.null(annot_labels))) {
    plt <- plt +
      geom_text_repel(aes(label = exc_group_id),
                      colour = 'black',
                      size = 2.5,
                      position = pos,
                      segment.colour = NA,
                      max.overlaps = 30)
  }
  if (plot_quadrants) {
    plt <- plt +
      annotate("text",
               x = horiz_high,
               y = vert_high,
               label = "Quadrant 1") +
      annotate("text",
               x = horiz_low,
               y = vert_high,
               label = "Quadrant 2") +
      annotate("text",
               x = horiz_low,
               y = vert_low,
               label = "Quadrant 3") +
      annotate("text",
               x = horiz_high,
               y = vert_low,
               label = "Quadrant 4")
  }
  plt <- plt +
    theme_bw() +
    theme(
      #panel.background = element_rect(fill = "transparent",
      #                                colour = NA_character_), # necessary to avoid drawing panel outline
      #panel.grid.major = element_blank(), # get rid of major grid
      #panel.grid.minor = element_blank(), # get rid of minor grid
      plot.background = element_rect(fill = "transparent",
                                     colour = NA_character_), # necessary to avoid drawing plot outline
      #legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
       legend.background = element_rect(fill = "transparent", color = "transparent"),
      #legend.background = element_rect(fill = "transparent",colour = "black", linetype='solid'),
      legend.box.background = element_rect(fill = "transparent", color = "transparent"),
      legend.key = element_rect(fill = "transparent")
    )
  
  if (facets_on) {
    plt <- plt +
      guides(colour = guide_legend(title = "Region",
                                   ncol = 5,
                                   title.position = "top"),
             shape = guide_legend(title = "Urbanicity",
                                  title.position = "top"),
             label = "none") +
      theme(legend.position="bottom",) +
      theme(text=element_text(size=15),
            legend.text=element_text(size=12)) +
      facet_grid(cols = vars(mass_int_name), rows = vars(net_name))
      #facet_grid(cols = vars(net_name), rows = vars(mass_int_name))
    
    # ggsave(paste(country,yvar,xvar,"facet_long_18NOV24.pdf", sep = "_"),
    #        plot = plt, bg = "white", w = 15, h = 10, dpi = 450)
    ggsave(paste(country,yvar,xvar,"facet_long_18NOV24.pdf", sep = "_"),
          plot = plt, bg = "transparent", w = 10, h = 12.5, dpi = 450)
  } else {
    plt <- plt +
      guides(colour = guide_legend(title = "Region"),
             shape = guide_legend(title = "Urbanicity"),
             label = "none")
    
    ggsave(paste(country,yvar,xvar,"quadrantfinal.pdf", sep = "_"),
           plot = plt, bg = "transparent", w = 8, h = 6, dpi = 450)
  }
    
  
  #print(plt)
  
  
  
}

cases_averted_cost_plot <- function(sim_sum,
                                    country = "SN",
                                    #xvar = "uret",
                                    yvar = "avert",
                                    strategy = "pyrethroid-only 3 year interval",
                                    per_xpop = 1,
                                    plot_quadrants = TRUE,
                                    facets_on = TRUE,
                                    annot_labels = NULL,
                                    rural_labels = NULL,
                                    urban_labels = NULL){
  
  sim_sum %<>% filter(ISO2 == country)
  if (facets_on) {
    if (yvar == "add_avert") {
      sim_sum %<>% filter(net_strategy != "pyrethroid-only 3 year interval")
    }
  } else {
    sim_sum %<>% filter(net_strategy == strategy)
  }
  
  sim_sum$mass_int_name <- paste0(sim_sum$mass_int, "-year interval")
  
  if (yvar == "add_avert") {
    sim_sum$ymean <- sim_sum$mean_add_avert_percap*per_xpop
    sim_sum$yLB <- sim_sum$LB_add_avert_percap*per_xpop
    sim_sum$yUB <- sim_sum$UB_add_avert_percap*per_xpop
  } else if (yvar == "avert") {
    sim_sum$ymean <- sim_sum$mean_avert_percap*per_xpop
    sim_sum$yLB <- sim_sum$LB_avert_percap*per_xpop
    sim_sum$yUB <- sim_sum$UB_avert_percap*per_xpop
  }
  
  plt <- ggplot(data = sim_sum,
                aes(x = fs_name_1,
                    y = ymean,
                    shape = urbanicity,
                    colour = fs_name_1))
  if (yvar == "add_avert") {
    if (per_xpop == 1) {per_char <- "capita"} else {per_char <- per_xpop}
    ytitle <- paste("Additional annual cases averted per", per_char)
    plt <- plt +
      geom_pointrange(aes(ymin = yLB, ymax = yUB),
                      size = 1.3,
                      position = pos,
                      alpha = 0.6) +
      scale_y_continuous(ytitle)
  } else if (yvar == "avert") {
    if (per_xpop == 1) {per_char <- "capita"} else {per_char <- per_xpop}
    ytitle <- paste("Annual cases averted per", per_char)
    plt <- plt +
      geom_pointrange(aes(ymin = yLB, ymax = yUB),
                      size = 1.3,
                      position = pos,
                      alpha = 0.6) +
      scale_y_continuous(ytitle)
  }
  
  if (facets_on) {
    plt <- plt +
      theme(text=element_text(size=20)) +
      facet_grid(cols = vars(net_name), rows = vars(mass_int_name))
    
    # ggsave(paste(country,yvar,xvar,"facet.png", sep = "_"), bg = "white",
    #        w = 15, h = 10, dpi = 450)
  } else {
    # ggsave(paste(country,yvar,xvar,"quadrant.png", sep = "_"), bg = "white",
    #        w = 8, h = 6, dpi = 450)
  }
  
  
  print(plt)
  
}
