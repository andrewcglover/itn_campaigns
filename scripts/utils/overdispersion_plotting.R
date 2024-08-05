# # overdispersion_plotting.R
# 
# repeated_overdispersion <- function(overdisp_df,
#                                     usage = FALSE,
#                                     access = FALSE,
#                                     LB = 0.025,
#                                     UB = 0.975,
#                                     mean_p = 0.5,
#                                     p_vals = seq(0.01,0.99,0.01)
# ) {
#   
#   if (usage & access) {
#     print("Warning: plot either usage or access, not both")
#     return(NULL)
#   }
#   
#   if (!usage & !access) {
#     print("Warning: neither usage or access called")
#     return(NULL)
#   }
#   
#   if (usage) {overdisp <- overdisp_u}
#   if (access) {overdisp <- overdisp_a}
#   
#   overdisp %<>% as.matrix
#   
#   N_p <- length(p_vals)
#   N_a <- dim(overdisp_df)[1]
#   N_r <- dim(overdisp)[1]
#   
#   print(N_r)
#   
#   overdisp_df_rep <- overdisp_df[rep(seq_len(nrow(overdisp_df)), each = N_p), ]
#   
#   overdisp_df_rep$LB_cumprobu <- rep(NA, N_a * N_p)
#   overdisp_df_rep$mid_cumprobu <- rep(NA, N_a * N_p)
#   overdisp_df_rep$mean_cumprobu <- rep(NA, N_a * N_p)
#   overdisp_df_rep$UB_cumprobu <- rep(NA, N_a * N_p)
#   overdisp_df_rep$p_vals <- rep(NA, N_a * N_p)
#   
#   l <- 1
#   for (i in 1:N_a) {
#     print(i)
#     for (j in 1:N_p) {
#       u_draws <- rep(NA, N_r)
#       for (k in 1:N_r) {
#         alpha0 <- overdisp[k,i]
#         alpha <- mean_p * alpha0
#         beta <- (1 - mean_p) * alpha0
#         #print(alpha)
#         u_draws[k] <- dbeta(p_vals[j], alpha, beta)
#       }
#       overdisp_df_rep$LB_cumprobu[l] <- quantile(u_draws, probs = LB)
#       overdisp_df_rep$mid_cumprobu[l] <- quantile(u_draws, probs = 0.5)
#       overdisp_df_rep$UB_cumprobu[l] <- quantile(u_draws, probs = UB)
#       overdisp_df_rep$mean_cumprobu[l] <- mean(u_draws)
#       overdisp_df_rep$p_vals[l] <- p_vals[j]
#       l <- l + 1
#     }
#   }
#   
#   return(overdisp_df_rep)
# }
# 
# plot_repeated_overdispersion <- function(overdisp_df_rep,
#                                  country = "SN",
#                                 mean_p = 0.5) {
#   
#   cc_overdisp_df_rep <- overdisp_df_rep %<>% filter(ISO2 == country)
#   
#   plt <- ggplot(data = cc_overdisp_df_rep,
#                 aes(x = p_vals,
#                     y = mid_cumprobu,
#                     ymin = LB_cumprobu,
#                     ymax = UB_cumprobu,
#                     colour = fs_name_1,
#                     fill = fs_name_1))
#   if (mean_p == 0.5) {plt <- plt + geom_abline(intercept = 0, slope = 1)}
#   if (mean_p == 0.5) {plt <- plt + geom_hline(yintercept = mean_p)}
#   plt <- plt +
#     geom_vline(xintercept = mean_p) +
#     geom_ribbon(alpha = 0.4, colour = NA) +
#     geom_path() +
#     facet_grid(urbanicity ~ fs_name_1)
#     
#   mean_pc <- mean_p * 100
#   
#   ggsave(paste0(ISO2,"_",mean_pc,"pc_overdisp_denprob_test_narrow.png"), bg = "white",
#          w = 16, h = 4, dpi = 450)
#     
#     print(plt)
#     
# 
#     
# }
# 
# 
# 
# 
# repeated_mean_cumulative_overdispersion <- function(
#     overdisp_df,
#     usage = FALSE,
#     access = FALSE,
#     LB = 0.025,
#     UB = 0.975,
#     default_p = 0.5,
#     use_yr1_usage = TRUE,
#     p_vals = seq(0.1,1,0.1)
# ) {
#   
#   if (!usage & !access) {
#     print("Warning: neither usage or access called")
#     return(NULL)
#   }
#   
#   diff_p <- p_vals[2] - p_vals[1]
#   
#   if (usage) {overdisp_u_mid <- overdisp_u %>% apply(2, median)}
#   if (access) {overdisp_a_mid <- overdisp_a %>% apply(2, median)}
#   
#   N_p <- length(p_vals)
#   N_a <- dim(overdisp_df)[1]
#   
#   overdisp_cum_df_rep <- overdisp_df[rep(seq_len(nrow(overdisp_df)), each = N_p), ]
#   
#   overdisp_cum_df_rep$upper_p_vals <- rep(NA, N_a * N_p)
#   if (usage) {overdisp_cum_df_rep$beta_u_den <- rep(NA, N_a * N_p)}
#   if (access) {overdisp_cum_df_rep$beta_a_den <- rep(NA, N_a * N_p)}
#   
#   k <- 1
#   for (i in 1:N_a) {
#     if (use_yr1_usage) {
#       if (overdisp_df$ISO2[i] %in% costed_sim_sum$ISO2 &
#           overdisp_df$fs_name_1[i] %in% costed_sim_sum$fs_name_1 &
#           overdisp_df$urbanicity[i] %in% costed_sim_sum$urbanicity) {
#         id <- which(costed_sim_sum$ISO2 == overdisp_df$ISO2[i] &
#                       costed_sim_sum$fs_name_1 == overdisp_df$fs_name_1[i] &
#                       costed_sim_sum$urbanicity == overdisp_df$urbanicity[i] &
#                       costed_sim_sum$net_name == "pyrethroid-only" &
#                       costed_sim_sum$mass_int == 3)
#         mean_p <- costed_sim_sum$mean_avg_yr1_use[id]
#         print(id)
#       } else {
#         mean_p <- default_p
#       }
#     }
#     for (j in 1:N_p) {
#       if (usage) {
#         alpha_u <- overdisp_u_mid[i] * mean_p
#         beta_u <- overdisp_u_mid[i] * (1 - mean_p)
#         pnew <- pbeta(p_vals[j], alpha_u, beta_u)
#         if (j > 1) {
#           pold <- pbeta(p_vals[j-1], alpha_u, beta_u)
#           overdisp_cum_df_rep$beta_u_den[k] <- pnew - pold
#         } else {
#           overdisp_cum_df_rep$beta_u_den[k] <- pnew
#         }
#       }
#       if (access) {
#         alpha_a <- overdisp_a_mid[i] * mean_p
#         beta_a <- overdisp_a_mid[i] * (1 - mean_p)
#         pnew <- pbeta(p_vals[j], alpha_a, beta_a)
#         if (j > 1) {
#           pold <- pbeta(p_vals[j-1], alpha_a, beta_a)
#           overdisp_cum_df_rep$beta_a_den[k] <- pnew - pold
#         } else {
#           overdisp_cum_df_rep$beta_a_den[k] <- pnew
#         }
#       }
#       overdisp_cum_df_rep$upper_p_vals[k] <- p_vals[j]
#       k <- k + 1
#     }
#   }
#   
#   overdisp_cum_df_rep$mid_p_vals <- overdisp_cum_df_rep$upper_p_vals - diff_p/2
#   
#   return(overdisp_cum_df_rep)
# }


repeated_cumulative_overdispersion <- function(
    overdisp_df,
    usage = FALSE,
    access = FALSE,
    LB = 0.025,
    UB = 0.975,
    default_p = 0.5,
    use_yr1_usage = TRUE,
    mean_only = FALSE,
    p_vals = seq(0.1,1,0.1)
) {
  
  if (!usage & !access) {
    print("Warning: neither usage or access called")
    return(NULL)
  }
  
  diff_p <- p_vals[2] - p_vals[1]
  
  if (usage) {overdisp_u_mid <<- overdisp_u %>% apply(2, mean)}
  if (access) {overdisp_a_mid <<- overdisp_a %>% apply(2, mean)}
  
  N_p <- length(p_vals)
  N_a <- dim(overdisp_df)[1]
  
  overdisp_cum_df_rep <- overdisp_df[rep(seq_len(nrow(overdisp_df)), each = N_p), ]
  
  overdisp_cum_df_rep$upper_p_vals <- rep(NA, N_a * N_p)
  if (usage) {overdisp_cum_df_rep$beta_u_den <- rep(NA, N_a * N_p)}
  if (access) {overdisp_cum_df_rep$beta_a_den <- rep(NA, N_a * N_p)}
  
  sim_sum_here <<- costed_sim_sum %>%
    left_join(sim_sum %>%
                filter(net_strategy == "pyrethroid-only 3 year interval") %>%
                ungroup() %>%
                select(fs_area_id, LB_uga, mean_uga, UB_uga),
              by = "fs_area_id")

  
  k <- 1
  for (i in 1:N_a) {
    if (use_yr1_usage) {
      if (overdisp_df$ISO2[i] %in% sim_sum_here$ISO2 &
          overdisp_df$fs_name_1[i] %in% sim_sum_here$fs_name_1 &
          overdisp_df$urbanicity[i] %in% sim_sum_here$urbanicity) {
        id <- which(sim_sum_here$ISO2 == overdisp_df$ISO2[i] &
                      sim_sum_here$fs_name_1 == overdisp_df$fs_name_1[i] &
                      sim_sum_here$urbanicity == overdisp_df$urbanicity[i] &
                      sim_sum_here$net_name == "pyrethroid-only" &
                      sim_sum_here$mass_int == 3)
        mean_p <- overdisp_df$ref_yr1_yr3_u[i]
        mean_pa <- overdisp_df$ref_yr1_yr3_a[i]
        print(id)
      } else {
        mean_p <- default_p
        mean_pa <- default_p
      }
    }
    if (i==132) {
      print(paste(overdisp_df$fs_name_1[i], overdisp_df$urbanicity[i],
                  "; usage = ", mean_p, "overdisp = ", overdisp_u_mid[i],
                  "; access = ", mean_pa, "overdisp = ", overdisp_a_mid[i]))
    }
    for (j in 1:N_p) {
      if (usage) {
        alpha_u <<- overdisp_u_mid[i] * mean_p
        beta_u <<- overdisp_u_mid[i] * (1 - mean_p)
        pnew_u <<- pbeta(p_vals[j], alpha_u, beta_u)
        if (j > 1) {
          pold_u <<- pbeta(p_vals[j-1], alpha_u, beta_u)
          overdisp_cum_df_rep$beta_u_den[k] <- pnew_u - pold_u
          if (i==130) {print(paste(p_vals[j],p_vals[j-1],pnew,pold,
                                   alpha_u, beta_u,
                                   mean_p, overdisp_u_mid[i],
                                   overdisp_cum_df_rep$beta_u_den[k],
                                   overdisp_cum_df_rep$fs_area[k]))}
        } else {
          overdisp_cum_df_rep$beta_u_den[k] <- pnew_u
        }
        overdisp_cum_df_rep$alpha_u[k] <- alpha_u
        overdisp_cum_df_rep$beta_u[k] <- beta_u
      }
      if (access) {
        alpha_a <- overdisp_a_mid[i] * mean_pa
        beta_a <- overdisp_a_mid[i] * (1 - mean_pa)
        pnew <- pbeta(p_vals[j], alpha_a, beta_a)
        if (j > 1) {
          pold <- pbeta(p_vals[j-1], alpha_a, beta_a)
          overdisp_cum_df_rep$beta_a_den[k] <- (pnew - pold)
        } else {
          overdisp_cum_df_rep$beta_a_den[k] <- pnew
        }
        overdisp_cum_df_rep$alpha_a[k] <- alpha_a
        overdisp_cum_df_rep$beta_a[k] <- beta_a
      }
      overdisp_cum_df_rep$upper_p_vals[k] <- p_vals[j]
      k <- k + 1
    }
  }
  
  overdisp_cum_df_rep$mid_p_vals <- overdisp_cum_df_rep$upper_p_vals - diff_p/2
  
  return(overdisp_cum_df_rep)
}

plot_mean_overdispersion <- function(overdisp_df_rep,
                                     country = "SN",
                                     filter_ADM = FALSE,
                                     fs_names_included = NULL,
                                     manual_cols = FALSE,
                                     plot_usage = FALSE,
                                     plot_access = FALSE,
                                     plot_rural = TRUE,
                                     plot_urban = TRUE,
                                     mean_p = 0.5,
                                     LB1 = 0.025,
                                     UB1 = 0.975,
                                     LB2 = 0.25,
                                     UB2 = 0.75) {
  
  # Filter country
  cc_overdisp_df_rep <- overdisp_df_rep %>% filter(ISO2 == country)
  
  # Filter urbanicity
  if (plot_rural & !plot_urban) {
    cc_overdisp_df_rep %<>% filter(urbanicity == "rural")
  } else if (!plot_rural & plot_urban) {
    cc_overdisp_df_rep %<>% filter(urbanicity == "urban")
  }
  
  # Filter countries selected
  if (filter_ADM) {
    cc_overdisp_df_rep %<>% filter(fs_name_1 %in% fs_names_included)
  }

  # Create plotting dataframe
  if (!plot_usage & !plot_access) {
    print("Warning: Neither usage and access or selected. Defaulting to usage")
    cc_overdisp_df_rep$beta_den <- cc_overdisp_df_rep$beta_u_den
    cc_overdisp_df_rep$mean_p <- cc_overdisp_df_rep$ref_yr1_yr3_u
    cc_overdisp_df_rep$LB1_p <- qbeta(LB1,
                                      cc_overdisp_df_rep$alpha_u,
                                      cc_overdisp_df_rep$beta_u)
    cc_overdisp_df_rep$UB1_p <- qbeta(UB1,
                                      cc_overdisp_df_rep$alpha_u,
                                      cc_overdisp_df_rep$beta_u)
    cc_overdisp_df_rep$LB2_p <- qbeta(LB2,
                                      cc_overdisp_df_rep$alpha_u,
                                      cc_overdisp_df_rep$beta_u)
    cc_overdisp_df_rep$UB2_p <- qbeta(UB2,
                                      cc_overdisp_df_rep$alpha_u,
                                      cc_overdisp_df_rep$beta_u)
  } else if (plot_usage & !plot_access) {
    cc_overdisp_df_rep$beta_den <- cc_overdisp_df_rep$beta_u_den
    cc_overdisp_df_rep$mean_p <- cc_overdisp_df_rep$ref_yr1_yr3_u
    cc_overdisp_df_rep$LB1_p <- qbeta(LB1,
                                      cc_overdisp_df_rep$alpha_u,
                                      cc_overdisp_df_rep$beta_u)
    cc_overdisp_df_rep$UB1_p <- qbeta(UB1,
                                      cc_overdisp_df_rep$alpha_u,
                                      cc_overdisp_df_rep$beta_u)
    cc_overdisp_df_rep$LB2_p <- qbeta(LB2,
                                      cc_overdisp_df_rep$alpha_u,
                                      cc_overdisp_df_rep$beta_u)
    cc_overdisp_df_rep$UB2_p <- qbeta(UB2,
                                      cc_overdisp_df_rep$alpha_u,
                                      cc_overdisp_df_rep$beta_u)
  } else if (!plot_usage & plot_access) {
    cc_overdisp_df_rep$beta_den <- cc_overdisp_df_rep$beta_a_den
    cc_overdisp_df_rep$mean_p <- cc_overdisp_df_rep$ref_yr1_yr3_a
    cc_overdisp_df_rep$LB1_p <- qbeta(LB1,
                                      cc_overdisp_df_rep$alpha_a,
                                      cc_overdisp_df_rep$beta_a)
    cc_overdisp_df_rep$UB1_p <- qbeta(UB1,
                                      cc_overdisp_df_rep$alpha_a,
                                      cc_overdisp_df_rep$beta_a)
    cc_overdisp_df_rep$LB2_p <- qbeta(LB2,
                                      cc_overdisp_df_rep$alpha_a,
                                      cc_overdisp_df_rep$beta_a)
    cc_overdisp_df_rep$UB2_p <- qbeta(UB2,
                                      cc_overdisp_df_rep$alpha_a,
                                      cc_overdisp_df_rep$beta_a)
  } else if (plot_usage & plot_access) {
    if (plot_rural & plot_urban) {
      Print("Warning: Usage and access selected without urban/rural filter")
    }
    N <- dim(cc_overdisp_df_rep)[1]
    beta_u_den <- cc_overdisp_df_rep$beta_u_den
    beta_a_den <- cc_overdisp_df_rep$beta_a_den
    mean_u <- cc_overdisp_df_rep$ref_yr1_yr3_u
    mean_a <- cc_overdisp_df_rep$ref_yr1_yr3_a
    LB1_u <- qbeta(LB1, cc_overdisp_df_rep$alpha_u, cc_overdisp_df_rep$beta_u)
    UB1_u <- qbeta(UB1, cc_overdisp_df_rep$alpha_u, cc_overdisp_df_rep$beta_u)
    LB2_u <- qbeta(LB2, cc_overdisp_df_rep$alpha_u, cc_overdisp_df_rep$beta_u)
    UB2_u <- qbeta(UB2, cc_overdisp_df_rep$alpha_u, cc_overdisp_df_rep$beta_u)
    LB1_a <- qbeta(LB1, cc_overdisp_df_rep$alpha_a, cc_overdisp_df_rep$beta_a)
    UB1_a <- qbeta(UB1, cc_overdisp_df_rep$alpha_a, cc_overdisp_df_rep$beta_a)
    LB2_a <- qbeta(LB2, cc_overdisp_df_rep$alpha_a, cc_overdisp_df_rep$beta_a)
    UB2_a <- qbeta(UB2, cc_overdisp_df_rep$alpha_a, cc_overdisp_df_rep$beta_a)
    cc_overdisp_df_rep %<>% rbind(cc_overdisp_df_rep)
    cc_overdisp_df_rep$beta_den <- c(beta_u_den, beta_a_den)
    cc_overdisp_df_rep$mean_p <- c(mean_u, mean_a)
    cc_overdisp_df_rep$LB1_p <- c(LB1_u, LB1_a)
    cc_overdisp_df_rep$UB1_p <- c(UB1_u, UB1_a)
    cc_overdisp_df_rep$LB2_p <- c(LB2_u, LB2_a)
    cc_overdisp_df_rep$UB2_p <- c(UB2_u, UB2_a)
    cc_overdisp_df_rep$usage_access <- c(rep("Usage", N), rep("Access", N))
    cc_overdisp_df_rep$rural_labs <- c(rep("Rural\nThiès", N/2),
                                         rep("Rural\nZiguinchor", N/2),
                                         rep("Rural\nThiès", N/2),
                                         rep("Rural\nZiguinchor", N/2))
  }
  

  cc_overdisp_df_rep <<- cc_overdisp_df_rep
  
  xbreaks <- seq(0,1,0.1)
  xlabels <- seq(0,100,10)
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  if (manual_cols) {
    all_cols <- gg_color_hue(14)
    man_cols <- all_cols[13:14] 
  }
  
  wdth <- 0.5
  
  plt <- ggplot(data = cc_overdisp_df_rep,
                aes(x = mid_p_vals,
                    y = beta_den * 100,
                    fill = fs_name_1,
                    colour = fs_name_1
                    )) +
    geom_vline(aes(xintercept = mean_p), size = wdth) +
    geom_vline(aes(xintercept = LB2_p), linetype = "dashed", size = wdth) +
    geom_vline(aes(xintercept = UB2_p), linetype = "dashed", size = wdth) +
    #geom_vline(aes(xintercept = LB2_p), linetype = "dashed", size = wdth) +
    #geom_vline(aes(xintercept = UB2_p), linetype = "dashed", size = wdth) +
    geom_col(alpha = 0.6, width = 0.1, linewidth = 0.25) +
    ylab("Proportion of the population (%)") +
    scale_x_continuous(breaks = xbreaks, labels = xlabels) +
    guides(fill = "none", colour = "none") +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 0),
      plot.background = element_rect(fill = "transparent",
                                     colour = NA_character_), # necessary to avoid drawing plot outline
      legend.background = element_rect(fill = "transparent"),
      legend.box.background = element_rect(fill = "transparent"),
      legend.key = element_rect(fill = "transparent")
    )
  
  if (!plot_usage & !plot_access) {
    plt <- plt + xlab("Probability of usage (%)")
  } else if (plot_usage & !plot_access) {
    plt <- plt + xlab("Probability of usage (%)")
  } else if (!plot_usage & plot_access) {
    plt <- plt + xlab("Probability of access (%)")
  } else if (plot_usage & plot_access) {
    plt <- plt + xlab("Probability of access or usage (%)")
  }
  
  
    
  
  if (plot_usage & plot_access) {
    plt <- plt + facet_grid(rural_labs ~ usage_access)
  } else {
    plt <- plt + facet_grid(fs_name_1 ~ urbanicity)
    
  }

  if (manual_cols) {
    plt <- plt +
      scale_fill_manual(values=man_cols) +
      scale_colour_manual(values=man_cols)
  }
  
  mean_pc <- mean_p * 100
  
  if (manual_cols) {
    ggsave(paste0(country,"_subset_yr1_3_overdisp_50CrI.pdf"), bg = "transparent",
           w = 6, h = 3, dpi = 450)
  } else if (plot_usage) {
    ggsave(paste0(country,"_yr1_3_overdisp_50CrI_usage.pdf"), bg = "transparent",
           w = 8, h = 10, dpi = 450)
  } else if (plot_access) {
    ggsave(paste0(country,"_yr1_3_overdisp_50CrI_access.pdf"), bg = "transparent",
           w = 8, h = 10, dpi = 450)
  }

  
  print(plt)
  
  
  
}



