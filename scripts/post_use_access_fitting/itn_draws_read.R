read_net_params <- function() {
  
  net_params_folder <- "data_private/net_params"
  
  dat_res_pyr_ll <<- readRDS(
    paste(
      ".",
      net_params_folder,
      "pyrethroid_uncertainty 2",
      sep = "/"
    )
  ) %>%
    dplyr::rename(rn0 = rn_pyr, gamman = mean_duration) %>% 
    dplyr::mutate(bioassay_mortality = 1 - resistance) %>%
    dplyr::select(dn0, rn0, gamman, resistance, bioassay_mortality) %>%
    cbind.data.frame("draw" = rep(1:1000, each = 101)) #%>%
  #cbind.data.frame("res_index" = rep(seq(1,101), 1000)) %>%
  #cbind.data.frame("index" = seq(1, 101000))
  
  dat_res_pbo_ll <<- readRDS(
    paste(
      ".",
      net_params_folder,
      "pbo_uncertainty_using_pyrethroid_dn0_for_mn_durability 2",
      sep = "/"
    )
  ) %>%
    dplyr::rename(rn0 = rn_pbo, gamman = mn_dur, resistance = resistance.x) %>%
    dplyr::mutate(bioassay_mortality = 1 - resistance) %>%
    dplyr::select(dn0, rn0, gamman, resistance, bioassay_mortality) %>%
    cbind.data.frame("draw" = rep(1:1000, each = 101)) #%>%
  #cbind.data.frame("res_index" = rep(seq(1,101), 1000)) %>%
  #cbind.data.frame("index" = seq(1, 101000))
  
  dat_res_pp_ll <<- readRDS(
    paste(
      ".",
      net_params_folder,
      "pyrrole_uncertainty_using_pyrethroid_dn0_for_mn_durability 2",
      sep = "/"
    )
  ) %>%
    dplyr::rename(rn0 = rn_pbo, gamman = mn_dur, resistance = resistance.x) %>%
    dplyr::mutate(bioassay_mortality = 1 - resistance) %>%
    dplyr::select(dn0, rn0, gamman, resistance, bioassay_mortality) %>%
    cbind.data.frame("draw" = rep(1:1000, each = 101)) #%>%
  # cbind.data.frame("res_index" = rep(seq(1,101), 1000)) %>%
  # cbind.data.frame("index" = seq(1, 101000))
  
  
  dat_res_pyr <<- readRDS(
    paste(
      ".",
      net_params_folder,
      "pyrethroid_binomial_uncertainty 1",
      sep = "/"
    )
  ) %>%
    dplyr::rename(rn0 = rn_pyr, gamman = mean_duration) %>% 
    dplyr::mutate(bioassay_mortality = 1 - resistance)  %>%
    dplyr::select(dn0, rn0, gamman, resistance, bioassay_mortality) %>%
    cbind.data.frame("draw" = rep(1:1000, each = 101)) #%>%
  # cbind.data.frame("res_index" = rep(seq(1,101), 1000)) %>%
  # cbind.data.frame("index" = seq(1, 101000))
  
  dat_res_pbo <<- readRDS(
    paste(
      ".",
      net_params_folder,
      "pbo_binomial_uncertainty_using_pyrethroid_dn0_for_mn_durability 1",
      sep = "/"
    )
  ) %>%
    dplyr::rename(rn0 = rn_pbo, gamman = mn_dur, resistance = resistance.x) %>%
    dplyr::mutate(bioassay_mortality = 1 - resistance) %>%
    dplyr::select(dn0, rn0, gamman, resistance, bioassay_mortality) %>%
    cbind.data.frame("draw" = rep(1:1000, each = 101)) #%>%
  # cbind.data.frame("res_index" = rep(seq(1,101), 1000)) %>%
  # cbind.data.frame("index" = seq(1, 101000))
  
  dat_res_pp <<- readRDS(
    paste(
      ".",
      net_params_folder,
      "pyrrole_binomial_uncertainty_using_pyrethroid_dn0_for_mn_durability 1",
      sep = "/"
    )
  ) %>%
    dplyr::rename(rn0 = rn_pbo, gamman = mn_dur, resistance = resistance.x) %>%
    dplyr::mutate(bioassay_mortality = 1 - resistance) %>%
    dplyr::select(dn0, rn0, gamman, resistance, bioassay_mortality) %>%
    cbind.data.frame("draw" = rep(1:1000, each = 101)) #%>%
  # cbind.data.frame("res_index" = rep(seq(1,101), 1000)) %>%
  # cbind.data.frame("index" = seq(1, 101000))
  
  dat_res <<- rbind.data.frame(
    cbind.data.frame(
      dat_res_pyr_ll,
      "itn_name" = rep("pyrethroid-only", dim(dat_res_pyr_ll)[1]),
      "rds_label" = rep(2, dim(dat_res_pyr_ll)[1])
    ),
    cbind.data.frame(
      dat_res_pbo_ll,
      "itn_name" = rep("pyrethroid-PBO", dim(dat_res_pyr_ll)[1]),
      "rds_label" = rep(2, dim(dat_res_pyr_ll)[1])
    ),
    cbind.data.frame(
      dat_res_pp_ll,
      "itn_name" = rep("pyrethroid-pyrrole", dim(dat_res_pyr_ll)[1]),
      "rds_label" = rep(2, dim(dat_res_pyr_ll)[1])
    ),
    cbind.data.frame(
      dat_res_pyr,
      "itn_name" = rep("pyrethroid-only", dim(dat_res_pyr_ll)[1]),
      "rds_label" = rep(1, dim(dat_res_pyr_ll)[1])
    ),
    cbind.data.frame(
      dat_res_pbo,
      "itn_name" = rep("pyrethroid-PBO", dim(dat_res_pyr_ll)[1]),
      "rds_label" = rep(1, dim(dat_res_pyr_ll)[1])
    ),
    cbind.data.frame(
      dat_res_pp,
      "itn_name" = rep("pyrethroid-pyrrole", dim(dat_res_pyr_ll)[1]),
      "rds_label" = rep(1, dim(dat_res_pyr_ll)[1])
    )
  )
  
  dat_res_filter <<- dat_res %>% dplyr::filter(resistance %in% seq(0,1,0.05))  
}