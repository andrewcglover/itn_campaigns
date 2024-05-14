N_fs_areas <- dim(fs_id_link)[1]

summary_df <- fs_id_link

ref_year <- 2022
ref_y0 <- 2020
ref_y1 <- 2022
ref_CMC <- date_to_CMC(ref_year, 6)
ref_CMC_series <- seq(date_to_CMC(ref_y0, 1), date_to_CMC(ref_y1, 12))

for (i in 1:N_fs_areas) {
  area_id <- fs_id_link$new_area_id[i]
  area_time_ref_id <- which(net_data$area_id == area_id &
                              net_data$CMC == ref_CMC)
  area_series_ref_id <- which(net_data$area_id == area_id &
                              net_data$CMC %in% ref_CMC_series)
  
  # get samples
  ret_ref_samples <- ret_u[, area_time_ref_id]
  avg_u_samples <- P_u[, area_time_ref_id] %>% rowMeans()
  avg_a_samples <- P_a[, area_time_ref_id] %>% rowMeans()
  avg_uga_samples <- avg_u_samples / avg_a_samples
  avg_uga_samples[which(uga_samples > 1)] <- 1
  # P0u_samples <- P0_u[, area_time_ref_id] %>%
  #   as.vector %>% unname %>% unlist
  # Pu_samples <- P_u[, area_time_ref_id] %>%
  #   as.vector %>% unname %>% unlist
  # Pa_samples <- P_a[, area_time_ref_id] %>%
  #   as.vector %>% unname %>% unlist
  # uga_samples <- Pu_samples / Pa_samples
  # uga_samples[which(uga_samples > 1)] <- 1
  
  summary_df$mean_ret[i] <- mean(ret_ref_samples)
  summary_df$LB95_ret[i] <- quantile(ret_ref_samples,
                                           probs = 0.025,
                                           names = FALSE)
  summary_df$UB95_ret[i] <- quantile(ret_ref_samples,
                                           probs = 0.975,
                                           names = FALSE)
  
  summary_df$mean_avgu[i] <- mean(avg_u_samples)
  summary_df$LB95_avgu[i] <- quantile(avg_u_samples,
                                     probs = 0.025,
                                     names = FALSE)
  summary_df$UB95_avgu[i] <- quantile(avg_u_samples,
                                     probs = 0.975,
                                     names = FALSE)
  
  summary_df$mean_avga[i] <- mean(avg_a_samples)
  summary_df$LB95_avga[i] <- quantile(avg_a_samples,
                                     probs = 0.025,
                                     names = FALSE)
  summary_df$UB95_avga[i] <- quantile(avg_a_samples,
                                     probs = 0.975,
                                     names = FALSE)
  
  summary_df$mean_avguga[i] <- mean(avg_uga_samples)
  summary_df$LB95_avguga[i] <- quantile(avg_uga_samples,
                                      probs = 0.025,
                                      names = FALSE)
  summary_df$UB95_avguga[i] <- quantile(avg_uga_samples,
                                      probs = 0.975,
                                      names = FALSE)
  
  # summary_df$mean_u0[i] <- mean(P0u_samples)
  # summary_df$LB95_u0[i] <- quantile(P0u_samples,
  #                                  probs = 0.025,
  #                                  names = FALSE)
  # summary_df$UB95_u0[i] <- quantile(P0u_samples,
  #                                  probs = 0.975,
  #                                  names = FALSE)
  # 
  # summary_df$mean_u[i] <- mean(Pu_samples)
  # summary_df$LB95_u[i] <- quantile(Pu_samples,
  #                                    probs = 0.025,
  #                                    names = FALSE)
  # summary_df$UB95_u[i] <- quantile(Pu_samples,
  #                                    probs = 0.975,
  #                                    names = FALSE)
  # 
  # summary_df$mean_a[i] <- mean(Pa_samples)
  # summary_df$LB95_a[i] <- quantile(Pa_samples,
  #                                    probs = 0.025,
  #                                    names = FALSE)
  # summary_df$UB95_a[i] <- quantile(Pa_samples,
  #                                    probs = 0.975,
  #                                    names = FALSE)
  # 
  # summary_df$mean_uga[i] <- mean(uga_samples)
  # summary_df$LB95_uga[i] <- quantile(uga_samples,
  #                                    probs = 0.025,
  #                                    names = FALSE)
  # summary_df$UB95_uga[i] <- quantile(uga_samples,
  #                                    probs = 0.975,
  #                                    names = FALSE)
}

ctry_df <- summary_df %>% filter(ISO2 == "SN")
ctry_df$fs_name_1[which(ctry_df$fs_name_1 == "Hauts-Bassins")] <- "Haut-Bassins"
ctry_df$fs_name_1[which(ctry_df$fs_name_1 == "Plateau Central")] <- "Plateau-Central"

fs_ctry <- get_site("SEN")

ctry_df$eir <- rep(NA, dim(ctry_df)[1])
for (i in 1:dim(ctry_df)[1]) {
  # Isolate a single site from a country
  adm_site_index <- which(fs_ctry$eir$name_1 == ctry_df$fs_name_1[i] &
                            fs_ctry$eir$urban_rural == ctry_df$urbanicity[i] &
                            fs_ctry$eir$spp == "pf")
  
  # If no foresite file for urban/rural, then revert to other
  if (identical(adm_site_index, integer(0))) {
    if (ctry_df$urbanicity[i] == "urban") {
      adm_site_index <- which(fs_ctry$eir$name_1 == ctry_df$fs_name_1[i] &
                                fs_ctry$eir$urban_rural == "rural" &
                                fs_ctry$eir$spp == "pf")
    } else {
      adm_site_index <- which(fs_ctry$eir$name_1 == ctry_df$fs_name_1[i] &
                                fs_ctry$eir$urban_rural == "urban" &
                                fs_ctry$eir$spp == "pf")
    }
  }
  ctry_df$eir[i] <- fs_ctry$eir$eir[adm_site_index]
}

for (i in 1:dim(ctry_df)[1]) {
  # Isolate a single site from a country
  adm_site_index <- which(fs_ctry$pyrethroid_resistance$name_1 == ctry_df$fs_name_1[i] &
                            fs_ctry$pyrethroid_resistance$urban_rural == ctry_df$urbanicity[i] &
                            fs_ctry$pyrethroid_resistance$year == ref_year)
  
  # If no foresite file for urban/rural, then revert to other
  if (identical(adm_site_index, integer(0))) {
    if (ctry_df$urbanicity[i] == "urban") {
      adm_site_index <- which(fs_ctry$pyrethroid_resistance$name_1 == ctry_df$fs_name_1[i] &
                                fs_ctry$pyrethroid_resistance$urban_rural == "rural" &
                                fs_ctry$pyrethroid_resistance$year == ref_year)
    } else {
      adm_site_index <- which(fs_ctry$pyrethroid_resistance$name_1 == ctry_df$fs_name_1[i] &
                                fs_ctry$pyrethroid_resistance$urban_rural == "urban" &
                                fs_ctry$pyrethroid_resistance$year == ref_year)
    }
  }
  ctry_df$pyr_res[i] <- fs_ctry$pyrethroid_resistance$pyrethroid_resistance[adm_site_index]
}

bi_df <- pyrrole2
tri_df <- pyrrole3

bi_df <- readRDS("BFpyrrole2.rds")
tri_df <- readRDS("BFpyrrole3.rds")

avert_df <- bi_df
avert_df$tri_to_bi_avert <- tri_df$pred_ann_infect - bi_df$pred_ann_infect
avert_df$tri_to_bi_avert_pcap <- avert_df$tri_to_bi_avert / avert_df$pop
#ctry_df$pyrrole_tri_to_bi_avert <- rep(NA, dim(ctry_df)[1])
for (i in 1:dim(ctry_df)[1]){
  adm_sample_ids <- which(avert_df$fs_area == ctry_df$fs_area[i])
  averted_samples <- avert_df$tri_to_bi_avert_pcap[adm_sample_ids]
  
  ctry_df$mean_3to2av[i] <- mean(averted_samples,
                                 na.rm = TRUE)
  ctry_df$LB95_3to2av[i] <- quantile(averted_samples,
                                     probs = 0.025,
                                     names = FALSE,
                                     na.rm = TRUE)
  ctry_df$UB95_3to2av[i] <- quantile(averted_samples,
                                     probs = 0.975,
                                     names = FALSE,
                                     na.rm = TRUE)
}

ctry_df$eir[which(ctry_df$eir < 1)] <- 1
ctry_df$LB95_avguga[which(ctry_df$LB95_avguga > 1)] <- 1
ctry_df$mean_avguga[which(ctry_df$mean_avguga > 1)] <- 1
ctry_df$UB95_avguga[which(ctry_df$UB95_avguga > 1)] <- 1

ctry_df$avert_uncert <- ctry_df$UB95_3to2av - ctry_df$LB95_3to2av
ctry_df$avert_uncert_flip <- max(ctry_df$avert_uncert, na.rm = TRUE) - ctry_df$avert_uncert

#SEN
uga_breaks <- seq(0.2, 1, 0.2) * 100
ret_breaks <- seq(12, 24, 3)
avt_breaks <- seq(0, 0.4, 0.05)
eir_breaks <- c(1,2,5,10,20,50,100,200,400)
uga_lim <- c(0.15, 1.05) * 100
ret_lim <- c(11, 25)
avt_lim <- c(0, 0.4)
#eir_lim <- 

#BFA
uga_breaks <- seq(0.7, 1, 0.1) * 100
ret_breaks <- seq(24, 39, 3)
eir_breaks <- c(1,2,5,10,20,50,100,200,500)
uga_lim <- c(0.65, 1.05) * 100
ret_lim <- c(23, 40)
#eir_lim <- 

bar_alpha = 0.6
size_lim = c(0, 1.0)
# avt_breaks <- seq(0, 2.5, 0.25)
# avt_lim = c(0, 2.5)


ggplot(data = ctry_df, aes(x = mean_ret,
                           xmin = LB95_ret,
                           xmax = UB95_ret,
                           y = mean_avguga * 100,
                           ymin = LB95_avguga * 100,
                           ymax = UB95_avguga * 100,
                           colour = fs_name_1,
                           shape = urbanicity)) +
  geom_errorbar(alpha = bar_alpha) +
  geom_errorbarh(alpha = bar_alpha) +
  geom_point() +
  scale_size(limits = c(0, 0.5)) +
  paletteer::scale_colour_paletteer_d("ggsci::category20_d3") +
  scale_x_continuous(breaks = ret_breaks,
                     limits = ret_lim) +
  scale_y_continuous(breaks = uga_breaks,
                     limits = uga_lim) +
  xlab("Mean retention (months)") +
  ylab("Usage given access (%)") +
  guides(shape = "none",
         size = "none",
         colour = "none") +
  theme_bw(base_size = 16)
#ggsave(paste0("SN_ADM_urleg",".png"), bg = "white",
ggsave(paste0("SN_ADM_uga_vs_ret",".png"), bg = "white",
       w = 3.75, h = 3.75, dpi = 400)

ggplot(data = ctry_df, aes(x = mean_ret,
                           xmin = LB95_ret,
                           xmax = UB95_ret,
                           y = eir,
                           colour = fs_name_1,
                           shape = urbanicity)) +
  geom_errorbarh(alpha = bar_alpha) +
  geom_point() +
  scale_size(limits = c(0, 0.5)) +
  paletteer::scale_colour_paletteer_d("ggsci::category20_d3") +
  scale_x_continuous(breaks = ret_breaks,
                     limits = ret_lim) +
  scale_y_continuous(breaks = eir_breaks,
                     trans = "log") +
  xlab("Mean retention (months)") +
  ylab("EIR") +
  guides(shape = "none",
         size = "none",
         colour = "none") +
  theme_bw(base_size = 16)
ggsave(paste0("SN_ADM_eir_vs_ret",".png"), bg = "white",
       w = 3.75, h = 3.75, dpi = 400)


ggplot(data = ctry_df, aes(y = eir,
                           x = mean_avguga * 100,
                           xmin = LB95_avguga * 100,
                           xmax = UB95_avguga * 100,
                           colour = fs_name_1,
                           shape = urbanicity)) +
  geom_errorbarh(alpha = bar_alpha) +
  geom_point() +
  scale_size(limits = c(0, 0.5)) +
  paletteer::scale_colour_paletteer_d("ggsci::category20_d3") +
  scale_y_continuous(breaks = eir_breaks,
                     trans = "log") +
  scale_x_continuous(breaks = uga_breaks,
                  limits = uga_lim) +
  # scale_x_reverse(breaks = uga_breaks,
  #                 limits = rev(uga_lim)) +
  ylab("EIR") +
  xlab("Usage given access (%)") +
  guides(shape = "none",
         size = "none",
         colour = "none") +
  theme_bw(base_size = 16)
ggsave(paste0("SN_ADM_eir_vs_uga",".png"), bg = "white",
       w = 3.75, h = 3.75, dpi = 400)



ctry_df_filter <- ctry_df %>% filter(!is.na(mean_3to2av))
  
ggplot(data = ctry_df_filter, aes(x = mean_ret,
                             xmin = LB95_ret,
                             xmax = UB95_ret,
                             y = mean_avguga * 100,
                             ymin = LB95_avguga * 100,
                             ymax = UB95_avguga * 100,
                             colour = mean_3to2av,
                             shape = urbanicity)) +
  geom_errorbar(alpha = bar_alpha) +
  geom_errorbarh(alpha = bar_alpha) +
  geom_point(aes(size = avert_uncert_flip)) +
  scale_size(limits = size_lim) +
  # scale_colour_viridis(option = "viridis",
  #                      breaks = avt_breaks,
  #                      labels = avt_breaks,
  #                      limits = avt_lim
  #                      ) +
  paletteer::scale_colour_paletteer_c("pals::ocean.haline",
                                      breaks = avt_breaks,
                                      labels = avt_breaks,
                                      limits = avt_lim
                                      ) +
  scale_x_continuous(breaks = ret_breaks,
                     limits = ret_lim) +
  scale_y_continuous(breaks = uga_breaks,
                     limits = uga_lim) +
  xlab("Mean retention (months)") +
  ylab("Usage given access (%)") +
  guides(shape = "none",
         size = "none",
         colour = "none"
         ) +
  theme_bw(base_size = 16)
ggsave(paste0("SN_3to2CFPav_uga_vs_ret_ocha",".png"), bg = "white",
       w = 3.75, h = 3.75, dpi = 400)


ggplot(data = ctry_df_filter, aes(x = mean_ret,
                           xmin = LB95_ret,
                           xmax = UB95_ret,
                           y = eir,
                           colour = mean_3to2av,
                           shape = urbanicity)) +
  geom_errorbarh(alpha = bar_alpha) +
  geom_point(aes(size = avert_uncert_flip)) +
  scale_size(limits = size_lim) +
  # scale_colour_viridis(option = "turbo",
  #                      breaks = avt_breaks,
  #                      labels = avt_breaks,
  #                      limits = avt_lim) +
  paletteer::scale_colour_paletteer_c("pals::ocean.haline",
                                      breaks = avt_breaks,
                                      labels = avt_breaks,
                                      limits = avt_lim
  ) +
  scale_x_continuous(breaks = ret_breaks,
                     limits = ret_lim) +
  scale_y_continuous(breaks = eir_breaks,
                     trans = "log") +
  xlab("Mean retention (months)") +
  ylab("EIR") +
  guides(shape = "none",
         size = "none",
         colour = "none") +
  theme_bw(base_size = 16)
ggsave(paste0("BF_3to2CFPav_eir_vs_ret_ocha",".png"), bg = "white",
       w = 3.75, h = 3.75, dpi = 400)


ggplot(data = ctry_df_filter, aes(y = eir,
                           x = mean_avguga * 100,
                           xmin = LB95_avguga * 100,
                           xmax = UB95_avguga * 100,
                           colour = mean_3to2av,
                           shape = urbanicity)) +
  geom_errorbarh(alpha = bar_alpha) +
  geom_point(aes(size = avert_uncert_flip)) +
  scale_size(limits = size_lim) +
  # scale_colour_viridis(option = "turbo",
  #                      breaks = avt_breaks,
  #                      labels = avt_breaks,
  #                      limits = avt_lim) +
  paletteer::scale_colour_paletteer_c("pals::ocean.haline",
                                      breaks = avt_breaks,
                                      labels = avt_breaks,
                                      limits = avt_lim
  ) +
  scale_y_continuous(breaks = eir_breaks,
                     trans = "log") +
  # scale_x_reverse(breaks = uga_breaks,
  #                 limits = rev(uga_lim)) +
  scale_x_continuous(breaks = uga_breaks,
                  limits = uga_lim) +
  ylab("EIR") +
  xlab("Usage given access (%)") +
  guides(shape = "none",
         size = "none",
         colour = "none") +
  # guides(shape = "none",
  #        size = "none",
  #        colour = guide_colorbar(title = "",
  #                              #title.position="top",
  #                              barwidth = 18,
  #                              barheight = 0.5)) +
  theme_bw(base_size = 16) 
  # theme(legend.position="bottom",
  #       legend.title.align=0.5)
ggsave(paste0("SN_3to2CFPav_eir_vs_uga_ocha",".png"), bg = "white",
       w = 3.75, h = 3.75, dpi = 400)
# ggsave(paste0("leg_3to2CFP_ocha",".png"), bg = "white",
#        w = 5, h = 5, dpi = 400)






ggplot(data = ctry_df, aes(x = mean_ret,
                           xmin = LB95_ret,
                           xmax = UB95_ret,
                           y = mean_3to2av,
                           ymin = LB95_3to2av,
                           ymax = UB95_3to2av,
                           colour = fs_name_1,
                           shape = urbanicity)) +
  geom_point() +
  geom_errorbar() +
  geom_errorbarh() +
  paletteer::scale_colour_paletteer_d("ggsci::category20_d3") +
  theme_bw()

ggplot(data = ctry_df, aes(x = mean_avguga,
                           xmin = LB95_avguga,
                           xmax = UB95_avguga,
                           y = mean_3to2av,
                           ymin = LB95_3to2av,
                           ymax = UB95_3to2av,
                           colour = fs_name_1,
                           shape = urbanicity)) +
  geom_point() +
  geom_errorbar() +
  geom_errorbarh() +
  scale_x_continuous(trans = "log") +
  paletteer::scale_colour_paletteer_d("ggsci::category20_d3") +
  theme_bw()

ggplot(data = ctry_df, aes(x = eir,
                           y = mean_3to2av,
                           ymin = LB95_3to2av,
                           ymax = UB95_3to2av,
                           colour = fs_name_1,
                           shape = urbanicity)) +
  geom_errorbar(position=position_dodge(width=0.1),
                alpha = 0.8) +
  geom_point(aes(size = pyr_res),
             position=position_dodge(width=0.1),
             alpha = 0.8) +
  scale_x_continuous(trans = "log") +
  paletteer::scale_colour_paletteer_d("ggsci::category20_d3") +
  theme_bw()

ggplot(data = ctry_df, aes(x = pyr_res,
                           y = mean_3to2av,
                           ymin = LB95_3to2av,
                           ymax = UB95_3to2av,
                           colour = fs_name_1,
                           shape = urbanicity)) +
  geom_point(aes(size = eir),
             position=position_dodge(width=0.01)) +
  geom_errorbar(position=position_dodge(width=0.01)) +
  paletteer::scale_colour_paletteer_d("ggsci::category20_d3") +
  theme_bw()



ggplot(data = ctry_df) +
  geom_point(aes(x = mean_ret,
                 y = mean_uga,
                 colour = fs_name_1,
                 shape = urbanicity)) +
  geom_errorbar(aes(x = mean_ret,
                    ymin = LB95_uga,
                    ymax = UB95_uga,
                    colour = fs_name_1)) +
  geom_errorbarh(aes(y = mean_uga,
                     xmin = LB95_ret,
                     xmax = UB95_ret,
                     colour = fs_name_1)) +
  paletteer::scale_colour_paletteer_d("ggsci::category20_d3") +
  theme_bw()


ggplot(data = ctry_df, aes(x = pfpr,
                           y = mean_uga,
                           colour = fs_name_1,
                           shape = urbanicity)) +
  geom_point(position=position_dodge(width=10)) +
  # geom_errorbar(aes(x = mean_ret,
  #                   ymin = LB95_uga,
  #                   ymax = UB95_uga,
  #                   colour = fs_name_1)) +
  geom_errorbar(aes(ymin = LB95_uga,
                    ymax = UB95_uga),
                position=position_dodge(width=10)) +
  paletteer::scale_colour_paletteer_d("ggsci::category20_d3") +
  theme_bw()