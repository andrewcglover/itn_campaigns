N_fs_areas <- dim(fs_id_link)[1]

net_retention_df <- fs_id_link

ref_CMC = 1476

for (i in 1:N_fs_areas) {
  area_id <- fs_id_link$new_area_id[i]
  area_time_ref_id <- which(net_data$area_id == area_id &
                              net_data$CMC == ref_CMC)
  ret_ref_samples <- ret_u[, area_time_ref_id] %>%
    as.vector %>% unname %>% unlist
  net_retention_df$mean_ret[i] <- mean(ret_ref_samples)
  net_retention_df$LB95_ret[i] <- quantile(ret_ref_samples,
                                           probs = 0.025,
                                           names = FALSE)
  net_retention_df$UB95_ret[i] <- quantile(ret_ref_samples,
                                           probs = 0.975,
                                           names = FALSE)
}