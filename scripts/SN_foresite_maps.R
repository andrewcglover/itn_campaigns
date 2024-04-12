fs_SEN <- get_site("SEN")

SN_year <- 2019

SN_prev <- fs_SEN$prevalence %>% filter(year == SN_year)
SN_pop <- fs_SEN$population %>% filter(year == SN_year)
SN_res <- fs_SEN$pyrethroid_resistance %>% filter(year == SN_year)

uni_name_1 <- unique(SN_prev$name_1)
N_n1 <- length(uni_name_1)

SN_fs_df <- data.frame("name_1" = uni_name_1,
                       "pfpr" = rep(0, N_n1),
                       "res" = rep(0, N_n1))


for (i in 1:N_n1) {
  
  SN_fs_df$name_1[i] <- uni_name_1[i]
  
  urb_id <- which(SN_prev$name_1 == uni_name_1[i] & SN_prev$urban_rural == "urban")
  urb_pop <- SN_pop$pop[urb_id]
  urb_pfpr <- SN_prev$pfpr[urb_id]
  urb_res <- SN_res$pyrethroid_resistance[urb_id]
  if (identical(numeric(0), urb_pop)) {urb_pop <- 0}
  if (identical(numeric(0), urb_pfpr)) {urb_pfpr <- 0}
  if (identical(numeric(0), urb_res)) {urb_res <- 0}
  
  rur_id <- which(SN_prev$name_1 == uni_name_1[i] & SN_prev$urban_rural == "rural")
  rur_pop <- SN_pop$pop[rur_id]
  rur_pfpr <- SN_prev$pfpr[rur_id]
  rur_res <- SN_res$pyrethroid_resistance[rur_id]
  if (identical(numeric(0), urb_pop)) {urb_pop <- 0}
  if (identical(numeric(0), rur_pfpr)) {rur_pfpr <- 0}
  if (identical(numeric(0), rur_res)) {rur_res <- 0}
  
  adm_pop <- urb_pop + rur_pop
  urb_weight <- urb_pop / adm_pop
  rur_weight <- rur_pop / adm_pop
  
  SN_fs_df$pfpr[i] <- urb_pfpr * urb_weight + rur_pfpr * rur_weight
  SN_fs_df$res[i] <- rur_res * urb_weight + rur_res * rur_weight
  
}

# get shapefiles
adm1.shp <- raster::getData("GADM", country = "SEN", level = 1)
adm1.shp.f <- sf::st_as_sf(adm1.shp, region = "NAME_1") #fortify
#adm1.shp.f <- tidy(adm1.shp, region = "NAME_1") #fortify

adm1.shp.f %<>%
  dplyr::rename(name_1 = NAME_1)

SN_fs_shapes <- merge(adm1.shp.f, as_tibble(SN_fs_df), by = "name_1")

ggplot() +
  geom_sf(data = SN_fs_shapes,
          aes(group = name_1,
              fill = pfpr * 100)) +
  scale_fill_viridis(option = "turbo",
                     #trans = "log",
                     begin = 0.4,
                     end = 0,
                     direction = 1,
                     #breaks = cases_avert_per_1000_breaks,
                     #labels = cases_avert_per_1000_breaks,
                     limits = c(0,8)
  ) +
  #guides(fill = guide_colorbar(title = bquote(PfPR[2-10]))) +
  guides(fill = "none") +
  theme_map()
ggsave(paste0("SN_pfpr",".png"), bg = "white",
       w = 5, h = 5, dpi = 1000)

ggplot() +
  geom_sf(data = SN_fs_shapes,
          aes(group = name_1,
              fill = res * 100)) +
  scale_fill_viridis(option = "turbo",
                     #trans = "log",
                     begin = 1,
                     end = 0.6,
                     direction = -1,
                     #breaks = cases_avert_per_1000_breaks,
                     #labels = cases_avert_per_1000_breaks,
                     limits = c(10,70)
  ) +
  #guides(fill = guide_colorbar(title = "Pyrethroid\nresistance (%)")) +
  guides(fill = "none") +
  theme_map()
ggsave(paste0("SN_res",".png"), bg = "white",
       w = 5, h = 5, dpi = 1000)

