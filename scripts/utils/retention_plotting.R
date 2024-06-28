# Retention_plotting.R

create_ret_density <- function(ref_CMC = 1476,
                               access = TRUE,
                               usage = TRUE) {
  
  # Prepare retention data frame
  N_samples <- dim(ret_a)[1]
  ret_df <- fs_id_link[rep(row.names(fs_id_link), each = N_samples),]
  N_area_samples <- dim(ret_df)[1]
  # if (access) {ret_df$ret_a_samples <- rep(NA, N_area_samples)}
  # if (usage) {ret_df$ret_u_samples <- rep(NA, N_area_samples)}
  
  # Loop over all foresite areas
  all_ret_a_samples <- c()
  all_ret_u_samples <- c()
  all_invlam_a_samples <- c()
  all_invlam_u_samples <- c()
  
  for (i in 1:N_fs_areas) {
      
      # Area-time indices
      area_id <- fs_id_link$new_area_id[i]
      area_time_ref_id <- which(net_data$area_id == area_id &
                                  net_data$CMC == ref_CMC)
      area_time_ids <- which(net_data$area_id == area_id)
      
      # get access samples
      if (access) {
        # ret_a_samples <- ret_a[, area_time_ref_id] %>%
        #   as.vector %>% unname %>% unlist
        area_ret_a_samples <- ret_a[, area_time_ref_id] %>%
          as.vector %>% unname %>% unlist
        area_invlam_a_samples <- invlam_a[, area_id] %>%
          as.vector %>% unname %>% unlist
        all_ret_a_samples <- c(all_ret_a_samples, area_ret_a_samples)
        all_invlam_a_samples <- c(all_invlam_a_samples, area_invlam_a_samples)
      }
      
      # get usage samples
      if (usage) {
        area_ret_u_samples <- ret_u[, area_time_ref_id] %>%
          as.vector %>% unname %>% unlist
        area_invlam_u_samples <- invlam_u[, area_id] %>%
          as.vector %>% unname %>% unlist
        all_ret_u_samples <- c(all_ret_u_samples, area_ret_u_samples)
        all_invlam_u_samples <- c(all_invlam_u_samples, area_invlam_u_samples)
      }
      #print(paste0("i =", i, "; rolling_access_samples = ", length(rolling_access_samples)))
  }
  
  if (access) {ret_df$ret_a_samples <- all_ret_a_samples}
  if (access) {ret_df$invlam_a_samples <- all_invlam_a_samples}
  if (usage) {ret_df$ret_u_samples <- all_ret_u_samples}
  if (usage) {ret_df$invlam_u_samples <- all_invlam_u_samples}
  
  # Return retention df samples
  return(ret_df)
}

plot_ctry_retention <- function(ISO2 = "SN",
                                plotting_var = "ret_u") {
  
  ctry_ret_data <- ret_df[which(ret_df$ISO2 == ISO2),]
  
  if (plotting_var == "ret_u") {
    ctry_ret_data$ret_samples <- ctry_ret_data$ret_u_samples/12
  } else if (plotting_var == "ret_a") {
    ctry_ret_data$ret_samples <- ctry_ret_data$ret_a_samples/12
  } else if (plotting_var == "invlam_u") {
    ctry_ret_data$ret_samples <- ctry_ret_data$invlam_u_samples/12
  } else if (plotting_var == "invlam_a") {
    ctry_ret_data$ret_samples <- ctry_ret_data$invlam_a_samples/12
  }
    
  ctry_ret_sum <- ctry_ret_data %>%
    dplyr::group_by(ADM1, urbanicity) %>%
    dplyr::summarise(ret_mean = mean(ret_samples),
                     ret_LB = quantile(ret_samples, probs = 0.025, na.rm = TRUE),
                     ret_UB = quantile(ret_samples, probs = 0.975, na.rm = TRUE))

  
  pos <- position_dodge(0.5)
  ylablims <- c(0,4)
  ylabticks <- seq(0,4,0.5)
  ggplot() +
    geom_violin(data = ctry_ret_data,
                aes(x = ADM1,
                    y = ret_samples,
                    fill = urbanicity),
                alpha = 0.4,
                color = NA,
                position = pos) +
  geom_errorbar(data = ctry_ret_sum,
                aes(x = ADM1,
                    ymin = ret_LB,
                    ymax = ret_UB,
                    color = urbanicity),
                alpha = 0.6,
                position = pos) +
    geom_point(data = ctry_ret_sum,
               aes(x = ADM1,
                   y = ret_mean,
                   color = urbanicity),
               alpha = 0.6,
               shape=16,
               position = pos) +
    guides(fill = guide_legend(title = "Urbanicity")) +
    guides(color = guide_legend(title = "Urbanicity")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    ylab("Mean retention (years)") +
    xlab(NULL) +
    scale_y_continuous(limits = ylablims,
                       breaks = ylabticks,
                       labels = ylabticks) +
  coord_flip()
    #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))
  
  ggsave(paste0(ISO2,"_",plotting_var,"_post2.png"), bg = "white",
         w = 8, h = 6, dpi = 450)
  
  
}






plot_ctry_access_usage_retention <- function(ISO2 = "SN",
                                             use_invlam = FALSE) {
  
  if (use_invlam) {ret_id <- "ua_invlam"} else {ret_id <- "ua_ret"} 
  
  ctry_ret_data <- ret_df[which(ret_df$ISO2 == ISO2),]
  # N0 <- dim(ctry_ret_data)[1]
  # 
  # if (use_invlam) {
  #   ret_ua_samples <- c(ctry_ret_data$invlam_u_samples/12,
  #                       ctry_ret_data$invlam_a_samples/12)
  # } else {
  #   ret_ua_samples <- c(ctry_ret_data$ret_U_samples/12,
  #                       ctry_ret_data$ret_a_samples/12)
  # }
  # 
  # ctry_ret_data %<>% rbind(ctry_ret_data)
  # ctry_ret_data$u
  # 
  if (use_invlam) {
    ctry_ret_data$ret_u_samples <- ctry_ret_data$invlam_u_samples/12
    ctry_ret_data$ret_a_samples <- ctry_ret_data$invlam_a_samples/12
  } else {
    ctry_ret_data$ret_u_samples <- ctry_ret_data$ret_u_samples/12
    ctry_ret_data$ret_a_samples <- ctry_ret_data$ret_a_samples/12
  }
  
  ctry_ret_sum <- ctry_ret_data %>%
    dplyr::group_by(ADM1, urbanicity) %>%
    dplyr::summarise(ret_u_mean = mean(ret_u_samples),
                     ret_u_LB = quantile(ret_u_samples,
                                         probs = 0.025,
                                         na.rm = TRUE),
                     ret_u_UB = quantile(ret_u_samples,
                                         probs = 0.975,
                                         na.rm = TRUE),
                     ret_a_mean = mean(ret_a_samples),
                     ret_a_LB = quantile(ret_a_samples,
                                         probs = 0.025,
                                         na.rm = TRUE),
                     ret_a_UB = quantile(ret_a_samples,
                                         probs = 0.975,
                                         na.rm = TRUE))
  
  # ctry_ret_sum <- ctry_ret_data %>%
  #   dplyr::group_by(ADM1, urbanicity) %>%
  #   dplyr::summarise(ret_a_mean = mean(ret_a_samples),
  #                    ret_a_LB = quantile(ret_a_samples,
  #                                        probs = 0.025,
  #                                        na.rm = TRUE),
  #                    ret_a_UB = quantile(ret_a_samples,
  #                                        probs = 0.975,
  #                                        na.rm = TRUE))
  
  pos <- position_dodge(0.5)
  ylablims <- c(0,4)
  ylabticks <- seq(0,4,0.5)
  ggplot() +
    # geom_violin(data = ctry_ret_data,
    #             aes(x = ADM1,
    #                 y = ret_u_samples,
    #                 fill = urbanicity),
    #             alpha = 0.4,
    #             color = NA,
    #             position = pos) +
    geom_errorbar(data = ctry_ret_sum,
                  aes(x = ADM1,
                      ymin = ret_u_LB,
                      ymax = ret_u_UB,
                      color = urbanicity),
                  alpha = 0.6,
                  width=0,
                  position = pos) +
    geom_point(data = ctry_ret_sum,
               aes(x = ADM1,
                   y = ret_u_mean,
                   color = urbanicity,
                   shape = urbanicity),
               alpha = 0.6,
               shape=16,
               position = pos) +
    scale_fill_viridis(
      alpha = 1,
      begin = 0.875,
      end = 0.975,
      direction = -1,
      discrete = TRUE,
      option = "H",
      guide = guide_legend(title = "Used nets")#,
      #labels = "3 year interval"#only_label_vals
    ) +
    scale_color_viridis(
      #alpha = 1,
      begin = 0.875,
      end = 0.975,
      direction = -1,
      discrete = TRUE,
      option = "H",
      guide = guide_legend(title = "Used nets")#,
      #labels = "3 year interval"#only_label_vals
    ) +
    new_scale_fill() +
    new_scale_colour() +
    # geom_violin(data = ctry_ret_data,
    #             aes(x = ADM1,
    #                 y = ret_a_samples,
    #                 fill = urbanicity),
    #             alpha = 0.4,
    #             color = NA,
    #             position = pos) +
    geom_errorbar(data = ctry_ret_sum,
                  aes(x = ADM1,
                      ymin = ret_a_LB,
                      ymax = ret_a_UB,
                      color = urbanicity),
                  alpha = 0.6,
                  width=0,
                  position = pos) +
    geom_point(data = ctry_ret_sum,
               aes(x = ADM1,
                   y = ret_a_mean,
                   color = urbanicity,
                   shape = urbanicity),
               alpha = 0.6,
               shape=16,
               position = pos) +
    scale_fill_viridis(
      #alpha = 1,
      begin = 0.025,
      end = 0.125,
      direction = 1,
      discrete = TRUE,
      option = "H",
      guide = guide_legend(title = "Accessible\nnets")
    ) +
    scale_color_viridis(
      #alpha = 1,
      begin = 0.025,
      end = 0.125,
      direction = 1,
      discrete = TRUE,
      option = "H",
      guide = guide_legend(title = "Accessible\nnets")
    ) +
    #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    ylab("Mean retention (years)") +
    xlab(NULL) +
    scale_y_continuous(limits = ylablims,
                       breaks = ylabticks,
                       labels = ylabticks) +
    guides(shape = "none") +
    #coord_flip()
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  ggsave(paste0(ISO2,"_",ret_id,"_post3.png"), bg = "white",
          w = 8, h = 6, dpi = 450)
  
  
}


plot_ctry_access_usage_retention2 <- function(ISO2 = "SN",
                                             use_invlam = FALSE) {
  
  if (use_invlam) {ret_id <- "ua_invlam"} else {ret_id <- "ua_ret"} 
  
  ctry_ret_data <- ret_df[which(ret_df$ISO2 == ISO2),]

  if (use_invlam) {
    ctry_ret_data$ret_u_samples <- ctry_ret_data$invlam_u_samples/12
    ctry_ret_data$ret_a_samples <- ctry_ret_data$invlam_a_samples/12
  } else {
    ctry_ret_data$ret_u_samples <- ctry_ret_data$ret_u_samples/12
    ctry_ret_data$ret_a_samples <- ctry_ret_data$ret_a_samples/12
  }
  
  ctry_ret_sum <- ctry_ret_data %>%
    dplyr::group_by(ADM1, urbanicity) %>%
    dplyr::summarise(ret_u_mean = mean(ret_u_samples),
                     ret_u_LB = quantile(ret_u_samples,
                                         probs = 0.025,
                                         na.rm = TRUE),
                     ret_u_UB = quantile(ret_u_samples,
                                         probs = 0.975,
                                         na.rm = TRUE),
                     ret_a_mean = mean(ret_a_samples),
                     ret_a_LB = quantile(ret_a_samples,
                                         probs = 0.025,
                                         na.rm = TRUE),
                     ret_a_UB = quantile(ret_a_samples,
                                         probs = 0.975,
                                         na.rm = TRUE))
  
  N <- dim(ctry_ret_sum)[1]
  
  ret_mean <- c(ctry_ret_sum$ret_u_mean, ctry_ret_sum$ret_a_mean)
  ret_LB <- c(ctry_ret_sum$ret_u_LB, ctry_ret_sum$ret_a_LB)
  ret_UB <- c(ctry_ret_sum$ret_u_UB, ctry_ret_sum$ret_a_UB)
  ctry_ret_sum %<>% rbind(ctry_ret_sum)
  ctry_ret_sum$ret_mean <- ret_mean
  ctry_ret_sum$ret_LB <- ret_LB
  ctry_ret_sum$ret_UB <- ret_UB
  ctry_ret_sum$usage_access <- c(rep("usage", N), rep("access", N))
  
  pos <- position_dodge(0.5)
  ylablims <- c(0,4)
  ylabticks <- seq(0,4,0.5)
  ggplot() +
    geom_errorbar(data = ctry_ret_sum,
                  aes(x = ADM1,
                      ymin = ret_LB,
                      ymax = ret_UB,
                      color = ADM1,
                      alpha = usage_access,
                      group = urbanicity
                      ),
                  #alpha = 1,
                  width=0,
                  position = pos
                  ) +
    geom_point(data = ctry_ret_sum,
               aes(x = ADM1,
                   y = ret_mean,
                   color = ADM1,
                   shape = urbanicity,
                   alpha = usage_access,
                   group = urbanicity
                   ),
               #alpha = 1,
               position = pos
               ) +
    # new_scale_fill() +
    # new_scale_colour() +
    # geom_errorbar(data = ctry_ret_sum,
    #               aes(x = ADM1,
    #                   ymin = ret_a_LB,
    #                   ymax = ret_a_UB,
    #                   color = ADM1,
    #                   group = ADM1),
    #               alpha = 0.5,
    #               #width=0,
    #               position = pos) +
    # geom_point(data = ctry_ret_sum,
    #            aes(x = ADM1,
    #                y = ret_a_mean,
    #                color = ADM1,
    #                shape = urbanicity),
    #            alpha = 0.5,
    #            position = pos) +
    ylab("Mean retention (years)") +
    xlab(NULL) +
    scale_y_continuous(limits = ylablims,
                       breaks = ylabticks,
                       labels = ylabticks) +
    scale_alpha_manual(values = c(0.4, 1)) +
    guides(shape = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  ggsave(paste0(ISO2,"_",ret_id,"_postcol.png"), bg = "white",
         w = 8, h = 6, dpi = 450)
  
  
}

plot_ctry_access_usage_retention3 <- function(ISO2 = "SN",
                                              use_invlam = FALSE,
                                              use_sim_sum = TRUE) {
  
  if (use_invlam) {ret_id <- "ua_invlam"} else {ret_id <- "ua_ret"} 
  
  if (use_sim_sum) {
    ctry_ret_sum <- sim_sum[which(sim_sum$ISO2 == ISO2),]
    ret_mean <- c(ctry_ret_sum$mean_retu, ctry_ret_sum$mean_reta)
    ret_LB <- c(ctry_ret_sum$LB_retu, ctry_ret_sum$LB_reta)
    ret_UB <- c(ctry_ret_sum$UB_retu, ctry_ret_sum$UB_reta)
    ctry_ret_sum$ADM1 <- ctry_ret_sum$fs_name_1
  } else {
    ctry_ret_data <- ret_df[which(ret_df$ISO2 == ISO2),]
    
    if (use_invlam) {
      ctry_ret_data$ret_u_samples <- ctry_ret_data$invlam_u_samples/12
      ctry_ret_data$ret_a_samples <- ctry_ret_data$invlam_a_samples/12
    } else {
      ctry_ret_data$ret_u_samples <- ctry_ret_data$ret_u_samples/12
      ctry_ret_data$ret_a_samples <- ctry_ret_data$ret_a_samples/12
    }
    
    ctry_ret_sum <- ctry_ret_data %>%
      dplyr::group_by(ADM1, urbanicity) %>%
      dplyr::summarise(ret_u_mean = mean(ret_u_samples),
                       ret_u_LB = quantile(ret_u_samples,
                                           probs = 0.025,
                                           na.rm = TRUE),
                       ret_u_UB = quantile(ret_u_samples,
                                           probs = 0.975,
                                           na.rm = TRUE),
                       ret_a_mean = mean(ret_a_samples),
                       ret_a_LB = quantile(ret_a_samples,
                                           probs = 0.025,
                                           na.rm = TRUE),
                       ret_a_UB = quantile(ret_a_samples,
                                           probs = 0.975,
                                           na.rm = TRUE))
    
    ret_mean <- c(ctry_ret_sum$ret_u_mean, ctry_ret_sum$ret_a_mean)
    ret_LB <- c(ctry_ret_sum$ret_u_LB, ctry_ret_sum$ret_a_LB)
    ret_UB <- c(ctry_ret_sum$ret_u_UB, ctry_ret_sum$ret_a_UB)
  }
  
  N <- dim(ctry_ret_sum)[1]
  
  ctry_ret_sum %<>% rbind(ctry_ret_sum)
  ctry_ret_sum$ret_mean <- ret_mean
  ctry_ret_sum$ret_LB <- ret_LB
  ctry_ret_sum$ret_UB <- ret_UB
  ctry_ret_sum$usage_access <- c(rep("Used nets", N), rep("Accessible nets", N))
  
  ctry_ret_sum$ret_UB[ctry_ret_sum$ret_UB > 42] <- 42
  
  pos <- position_dodge(0.5)
  ylablims <- c(1,3.5)*12
  ylabticks <- seq(0,4,0.25)*12
  ggplot() +
    geom_hline(yintercept = seq(12,36,12), colour = "grey80") +
    geom_errorbar(data = ctry_ret_sum,
                  aes(x = ADM1,
                      ymin = ret_LB,
                      ymax = ret_UB,
                      color = ADM1,
                      #alpha = usage_access,
                      group = urbanicity
                  ),
                  #alpha = 1,
                  width=0,
                  position = pos
    ) +
    geom_point(data = ctry_ret_sum,
               aes(x = ADM1,
                   y = ret_mean,
                   color = ADM1,
                   shape = urbanicity,
                   #alpha = usage_access,
                   group = urbanicity
               ),
               #alpha = 1,
               position = pos
    ) +
    # new_scale_fill() +
    # new_scale_colour() +
    # geom_errorbar(data = ctry_ret_sum,
    #               aes(x = ADM1,
    #                   ymin = ret_a_LB,
    #                   ymax = ret_a_UB,
    #                   color = ADM1,
    #                   group = ADM1),
    #               alpha = 0.5,
    #               #width=0,
    #               position = pos) +
  # geom_point(data = ctry_ret_sum,
  #            aes(x = ADM1,
  #                y = ret_a_mean,
  #                color = ADM1,
  #                shape = urbanicity),
  #            alpha = 0.5,
  #            position = pos) +
  ylab("Mean duration (months)") +
    xlab(NULL) +
    scale_y_continuous(limits = ylablims,
                       breaks = ylabticks,
                       labels = ylabticks) +
    scale_alpha_manual(values = c(0.4, 1)) +
    guides(colour = "none",
           shape = guide_legend("Urbanicity")) +
    labs(shape = "Urbanicity") +
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1),
          #legend.position = "bottom"
          legend.position = c(.97, .97),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6)
          ) +
    facet_grid(cols = vars(usage_access))
  
  ggsave(paste0(ISO2,"_",ret_id,"_postfacetest.png"), bg = "white",
         w = 8, h = 6, dpi = 450)
  
  
}

