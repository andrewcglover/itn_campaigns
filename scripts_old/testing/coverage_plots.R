# coverage_plots.R

cov_df <- fs_id_link

cov_df$D_a_lo <- rep(NA, N_fs_areas)
cov_df$D_a_mid <- rep(NA, N_fs_areas)
cov_df$D_a_hi <- rep(NA, N_fs_areas)
cov_df$mean_a_lo <- rep(NA, N_fs_areas)
cov_df$mean_a_mid <- rep(NA, N_fs_areas)
cov_df$mean_a_hi <- rep(NA, N_fs_areas)

CMC_date <- date_to_CMC(2022,12)
del_t <- 36

LB <- 0.025
mid <- 0.5
UB <- 0.975
scale <- 100

for (i in 1:N_fs_areas) {
  ii <- fs_id_link$new_area_id[i]
  j <- which(net_data$area_id == ii & net_data$CMC == CMC_date)
  D_ai <- D_a[,j] %>% as.vector %>% unname %>% unlist
  P0_ai <- P0_a[,j] %>% as.vector %>% unname %>% unlist
  invlam_ai <- invlam_a[,ii] %>% as.vector %>% unname %>% unlist
  lam_ai <- 1 / invlam_ai
  
  C0_ai <- P0_ai - D_ai
  mean_C_ai <- (C0_ai / (lam_ai * del_t)) * (1 - exp(-lam_ai * del_t))
  mean_P_ai <- mean_C_ai + D_ai
  
  cov_df$D_a_lo[i] <- quantile(D_ai, probs = LB) * scale
  cov_df$D_a_mid[i] <- quantile(D_ai, probs = mid) * scale
  cov_df$D_a_hi[i] <- quantile(D_ai, probs = UB) * scale
  cov_df$mean_a_lo[i] <- quantile(mean_P_ai, probs = LB) * scale
  cov_df$mean_a_mid[i] <- quantile(mean_P_ai, probs = mid) * scale
  cov_df$mean_a_hi[i] <- quantile(mean_P_ai, probs = UB) * scale
}

cov_df$Country <- countrycode(cov_df$ISO2, "iso2c", "country.name")

ggplot(data = cov_df,
       aes(x = D_a_mid,
           xmin = D_a_lo,
           xmax = D_a_hi,
           y = mean_a_mid,
           ymin = mean_a_lo,
           ymax = mean_a_hi,
           colour = Country)) +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point() +
  xlab("Access to routine ITNs (%)") +
  ylab("Mean access to any ITN (%)")