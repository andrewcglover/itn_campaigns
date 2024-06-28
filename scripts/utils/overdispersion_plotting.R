# overdispersion_plotting.R

repeated_overdispersion <- function(overdisp_df,
                                    usage = FALSE,
                                    access = FALSE,
                                    LB = 0.025,
                                    UB = 0.975,
                                    mean_p = 0.5,
                                    p_vals = seq(0.01,0.99,0.01)
) {
  
  if (usage & access) {
    print("Warning: plot either usage or access, not both")
    return(NULL)
  }
  
  if (!usage & !access) {
    print("Warning: neither usage or access called")
    return(NULL)
  }
  
  if (usage) {overdisp <- overdisp_u}
  if (access) {overdisp <- overdisp_a}
  
  overdisp %<>% as.matrix
  
  N_p <- length(p_vals)
  N_a <- dim(overdisp_df)[1]
  N_r <- dim(overdisp)[1]
  
  print(N_r)
  
  overdisp_df_rep <- overdisp_df[rep(seq_len(nrow(overdisp_df)), each = N_p), ]
  
  overdisp_df_rep$LB_cumprobu <- rep(NA, N_a * N_p)
  overdisp_df_rep$mid_cumprobu <- rep(NA, N_a * N_p)
  overdisp_df_rep$mean_cumprobu <- rep(NA, N_a * N_p)
  overdisp_df_rep$UB_cumprobu <- rep(NA, N_a * N_p)
  overdisp_df_rep$p_vals <- rep(NA, N_a * N_p)
  
  l <- 1
  for (i in 1:N_a) {
    print(i)
    for (j in 1:N_p) {
      u_draws <- rep(NA, N_r)
      for (k in 1:N_r) {
        alpha0 <- overdisp[k,i]
        alpha <- mean_p * alpha0
        beta <- (1 - mean_p) * alpha0
        #print(alpha)
        u_draws[k] <- dbeta(p_vals[j], alpha, beta)
      }
      overdisp_df_rep$LB_cumprobu[l] <- quantile(u_draws, probs = LB)
      overdisp_df_rep$mid_cumprobu[l] <- quantile(u_draws, probs = 0.5)
      overdisp_df_rep$UB_cumprobu[l] <- quantile(u_draws, probs = UB)
      overdisp_df_rep$mean_cumprobu[l] <- mean(u_draws)
      overdisp_df_rep$p_vals[l] <- p_vals[j]
      l <- l + 1
    }
  }
  
  return(overdisp_df_rep)
}

plot_overdispersion <- function(overdisp_df_rep,
                                 country = "SN",
                                mean_p = 0.5) {
  
  cc_overdisp_df_rep <- overdisp_df_rep %<>% filter(ISO2 == country)
  
  plt <- ggplot(data = cc_overdisp_df_rep,
                aes(x = p_vals,
                    y = mid_cumprobu,
                    ymin = LB_cumprobu,
                    ymax = UB_cumprobu,
                    colour = fs_name_1,
                    fill = fs_name_1))
  if (mean_p == 0.5) {plt <- plt + geom_abline(intercept = 0, slope = 1)}
  if (mean_p == 0.5) {plt <- plt + geom_hline(yintercept = mean_p)}
  plt <- plt +
    geom_vline(xintercept = mean_p) +
    geom_ribbon(alpha = 0.4, colour = NA) +
    geom_path() +
    facet_grid(urbanicity ~ fs_name_1)
    
  mean_pc <- mean_p * 100
  
  ggsave(paste0(ISO2,"_",mean_pc,"pc_overdisp_denprob_test_narrow.png"), bg = "white",
         w = 16, h = 4, dpi = 450)
    
    print(plt)
    

    
}