#appending_malsimsum.R

append_previous_estimates <- function(dataset = sim_sum,
                                      retention_usage = TRUE,
                                      retention_access = TRUE,
                                      use_given_access = TRUE,
                                      ref_access = TRUE,
                                      ref_yr1_to_yr3_access = TRUE,
                                      ref_usage = TRUE,
                                      ref_yr1_to_yr3_usage = TRUE,
                                      LB_val = 0.025,
                                      UB_val = 0.975,
                                      ref_CMC = 1476,
                                      ref_CMC_low = 1405,
                                      ref_CMC_high = 1476) {
  
  N <- dim(dataset)[1]
  
  if (retention_usage) {
    dataset$LB_retu = rep(NA, N)
    dataset$mean_uretu = rep(NA, N)
    dataset$UB_retu = rep(NA, N)
  }
  
  if (retention_access) {
    dataset$LB_reta = rep(NA, N)
    dataset$mean_reta = rep(NA, N)
    dataset$UB_reta = rep(NA, N)
  }
  
  if(use_given_access)
  {
    dataset$LB_uga = rep(NA, N)
    dataset$mean_uga = rep(NA, N)
    dataset$UB_uga = rep(NA, N)
    if (!exists("P_uga")) {P_uga <<- P_u / P_a}
  }  
  
  for (i in 1:N){
    
    fs_area_id <- dataset$fs_area_id[i]
    fs_id_match <- which(fs_id_link$fs_area_id == fs_area_id)
    new_area_id <- fs_id_link$new_area_id[fs_id_match]
    area_time_ref_low_id <- which(net_data$area_id == new_area_id &
                                net_data$CMC == ref_CMC_low)
    area_time_ref_high_id <- which(net_data$area_id == new_area_id &
                                    net_data$CMC == ref_CMC_high)
    area_time_ref_id <- which(net_data$area_id == new_area_id &
                                    net_data$CMC == ref_CMC)
    
    if (retention_usage) {
      retu_samples <- ret_u[, area_time_ref_low_id:area_time_ref_high_id]
      retu_samples %<>% rowMeans
      dataset$LB_retu[i] <- quantile(retu_samples, LB_val, na.rm = TRUE)# %>% unname
      dataset$mean_retu[i] <- mean(retu_samples, na.rm = TRUE)
      dataset$UB_retu[i] <- quantile(retu_samples, UB_val, na.rm = TRUE)# %>% unname
    }
    
    if (retention_access) {
      reta_samples <- ret_a[, area_time_ref_low_id:area_time_ref_high_id]
      reta_samples %<>% rowMeans
      dataset$LB_reta[i] <- quantile(reta_samples, LB_val, na.rm = TRUE)# %>% unname
      dataset$mean_reta[i] <- mean(reta_samples, na.rm = TRUE)
      dataset$UB_reta[i] <- quantile(reta_samples, UB_val, na.rm = TRUE)# %>% unname
    }
    
    if (use_given_access) {
      uga_samples <- P_uga[, area_time_ref_low_id:area_time_ref_high_id]
      uga_samples %<>% rowMeans
      dataset$LB_uga[i] <- quantile(uga_samples, LB_val, na.rm = TRUE) #%>% unname
      dataset$mean_uga[i] <- mean(uga_samples, na.rm = TRUE)
      dataset$UB_uga[i] <- quantile(uga_samples, UB_val, na.rm = TRUE) #%>% unname
    }
    
    if (ref_access) {
      lambda_a_samples <- lamrep_a[, area_time_ref_id]
      dataset$lambda_a_ref_LB[i] <- lambda_a_samples %>% quantile(LB_val, na.rm = TRUE)
      dataset$lambda_a_ref_mean[i] <- lambda_a_samples %>% mean(na.rm = TRUE)
      dataset$lambda_a_ref_UB[i] <- lambda_a_samples %>% mean(UB_val, na.rm = TRUE)
      P0_a_samples <- P0_a[, area_time_ref_id] %>% unlist
      #P0_a %<>% rowMeans
      dataset$P0_a_ref_LB[i] <- P0_a_samples %>% quantile(LB_val, na.rm = TRUE)
      dataset$P0_a_ref_mean[i] <- P0_a_samples %>% mean(na.rm = TRUE)
      dataset$P0_a_ref_UB[i] <- P0_a_samples %>% quantile(UB_val, na.rm = TRUE)
      D_a_samples <- D_a[, area_time_ref_id] %>% unlist
      #D_a_samples %<>% rowMeans
      dataset$D_a_ref_LB[i] <- D_a_samples %>% quantile(LB_val, na.rm = TRUE)
      dataset$D_a_ref_mean[i] <- D_a_samples %>% mean(na.rm = TRUE)
      dataset$D_a_ref_UB[i] <- D_a_samples %>% quantile(UB_val, na.rm = TRUE)
    }
    
    if (ref_yr1_to_yr3_access) {
      d <- dataset$D_a_ref_mean[i]
      p0 <- dataset$P0_a_ref_mean[i]
      lam <- dataset$lambda_a_ref_mean[i]
      t <- 1 / ((36*lam^2)/(-36*lam+exp(36*lam)-1) + lam)
      a <- (p0-d) * exp(-lam*t) + d
      dataset$ref_yr1_yr3_a[i] <- a
    }
    
    if (ref_usage) {
      lambda_u_samples <- lamrep_u[, area_time_ref_id]
      dataset$lambda_u_ref_LB[i] <- lambda_u_samples %>% quantile(LB_val, na.rm = TRUE)
      dataset$lambda_u_ref_mean[i] <- lambda_u_samples %>% mean(na.rm = TRUE)
      dataset$lambda_u_ref_UB[i] <- lambda_u_samples %>% mean(UB_val, na.rm = TRUE)
      P0_u_samples <- P0_u[, area_time_ref_id] %>% unlist
      #P0_u %<>% rowMeans
      dataset$P0_u_ref_LB[i] <- P0_u_samples %>% quantile(LB_val, na.rm = TRUE)
      dataset$P0_u_ref_mean[i] <- P0_u_samples %>% mean(na.rm = TRUE)
      dataset$P0_u_ref_UB[i] <- P0_u_samples %>% quantile(UB_val, na.rm = TRUE)
      D_u_samples <- D_u[, area_time_ref_id] %>% unlist
      #D_u_samples %<>% rowMeans
      dataset$D_u_ref_LB[i] <- D_u_samples %>% quantile(LB_val, na.rm = TRUE)
      dataset$D_u_ref_mean[i] <- D_u_samples %>% mean(na.rm = TRUE)
      dataset$D_u_ref_UB[i] <- D_u_samples %>% quantile(UB_val, na.rm = TRUE)
    }
    
    if (ref_yr1_to_yr3_usage) {
      d <- dataset$D_u_ref_mean[i]
      p0 <- dataset$P0_u_ref_mean[i]
      lam <- dataset$lambda_u_ref_mean[i]
      t <- 1 / ((36*lam^2)/(-36*lam+exp(36*lam)-1) + lam)
      u <- (p0-d) * exp(-lam*t) + d
      dataset$ref_yr1_yr3_u[i] <- u
    }

  }
  
  return(dataset)
  
}