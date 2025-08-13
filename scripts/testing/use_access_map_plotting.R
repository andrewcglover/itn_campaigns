library(cartogram)
library(tidyverse)
library(ggpattern)
library(sf)
library(purrr)
library(cowplot)
library(tmap)

net_data_end <- net_data %>%
  dplyr::filter(CMC == CMC_last)

site_list <- list()
for (i in 1:N_ISO2) {
  site_list[[i]] <- fetch_site(
    iso3c = countrycode(SSA_ISO2[i], "iso2c", "iso3c")
  )
}



ua_end_stats <- data.frame(
  "fs_area_id" = net_data_end$fs_area_id,
  "D_a" = net_data_end$D_a_mean,
  "C0_a" = net_data_end$P0_a_mean - net_data_end$D_a_mean,
  "lambda_a" = 1 / net_data_end$invlam_a_mean,
  "D_u" = net_data_end$D_u_mean,
  "C0_u" = net_data_end$P0_u_mean - net_data_end$D_u_mean,
  "lambda_u" = 1 / net_data_end$invlam_u_mean
)

ua_end <- left_join(fs_id_link, ua_end_stats, by = "fs_area_id")

ISO_pops <- NULL
for (i in 1:N_ISO2) {
  pop_last <- site_list[[i]]$population$population_total %>%
    filter(year == unname(unlist(CMC_to_date(CMC_last)[1])))
  ISO_pops %<>% rbind(pop_last)
}
ISO_pops %<>% rename(c("urbanicity" = "urban_rural", "fs_name_1" = "name_1"))
ISO_pops$ISO2 <- countrycode(ISO_pops$iso3c, "iso3c", "iso2c")

ua_end %<>% left_join(ISO_pops, by = c("ISO2", "fs_name_1", "urbanicity"))

ua_end %<>% group_by(ISO2, fs_name_1) %>%
  mutate(
    adm_pop = sum(pop)
  )

ua_end %<>% rbind.data.frame(ua_end)

avg_use_access <- function(D, C0, lambda, del_t) {
  avg_val <- D + ( (C0 / (lambda * del_t) ) * (1 - exp(-lambda * del_t) ) )
}

time_above_threshold <- function(D, C0, lambda, del_t, thres){
  tau <- -log( (thres - D) / C0) / lambda
}

avg_use_given_access <- function(
    D_u, C0_u, lambda_u, D_a, C0_a, lambda_a, del_t, nsteps = NULL, cap1 = TRUE
) {
  nregions <- length(D_u)
  avg_uga <- rep(NA, nregions)
  for (j in 1:nregions) {
    if (is.null(nsteps)) {nsteps <- del_t[j]}
    del_s <- del_t[j] / nsteps
    cum_uga <- 0
    Cold_a <- C0_a[j]
    Cold_u <- C0_u[j]
    for (i in 1:nsteps) {
      uga_step <- (
        avg_use_access(D = D_u[j], C0 = Cold_u, lambda = lambda_u[j], del_s) /
          avg_use_access(D = D_a[j], C0 = Cold_a, lambda = lambda_a[j], del_s)
      )
      if (cap1) {if (uga_step > 1) {uga_step <- 1} }
      cum_uga <- cum_uga + uga_step
      Cold_a <- Cold_a * exp(-lambda_a[j] * del_s)
      Cold_u <- Cold_u * exp(-lambda_u[j] * del_s)
    }
    avg_uga[j] <- cum_uga / nsteps
  }
  return(avg_uga)
}

ua_end$interval_yr <- c(rep(2, dim(ua_end)[1]/2), rep(3, dim(ua_end)[1]/2))
ua_end$interval_mn <- ua_end$interval_yr * 12

ua_end$avg_a <- avg_use_access(
  D = ua_end$D_a,
  C0 = ua_end$C0_a,
  lambda = ua_end$lambda_a,
  del_t = ua_end$interval_mn
)

ua_end$avg_u <- avg_use_access(
  D = ua_end$D_u,
  C0 = ua_end$C0_u,
  lambda = ua_end$lambda_u,
  del_t = ua_end$interval_mn
)

ua_end$avg_uga <- avg_use_given_access(
  D_u = ua_end$D_u,
  C0_u = ua_end$C0_u,
  lambda_u = ua_end$lambda_u,
  D_a = ua_end$D_a,
  C0_a = ua_end$C0_a,
  lambda_a = ua_end$lambda_a,
  del_t = ua_end$interval_mn
)

ua_end$over80time_a <- time_above_threshold(
  D = ua_end$D_a,
  C0 = ua_end$C0_a,
  lambda = ua_end$lambda_a,
  del_t = ua_end$interval_mn,
  thres = 0.8
)
ua_end$over80time_a[ua_end$over80time_a < 0] <- 0
ua_end$over80time_a <- ifelse(
  is.nan(ua_end$over80time_a),
  ua_end$interval_mn,
  ua_end$over80time_a
  )

ua_end$over80time_u <- time_above_threshold(
  D = ua_end$D_u,
  C0 = ua_end$C0_u,
  lambda = ua_end$lambda_u,
  del_t = ua_end$interval_mn,
  thres = 0.8
)
ua_end$over80time_u[ua_end$over80time_u < 0] <- 0
ua_end$over80time_u <- ifelse(
  is.nan(ua_end$over80time_u),
  ua_end$interval_mn,
  ua_end$over80time_u
)

ua_end$propover80_a <- ua_end$over80time_a / ua_end$interval_mn
ua_end$propover80_u <- ua_end$over80time_u / ua_end$interval_mn

# ISO_sitefile_shapes <- NULL
# for (i in 1:N_ISO2) {
#   ISO_sitefile_shapes %<>% rbind(site_list[[i]]$shape$level_1)
# }

adm1.shp <- geodata::gadm(country = SSA_ISO2, level = 1, 
                             path = tempdir(),
                             version = "4.0")
adm1.shp.f <- sf::st_as_sf(adm1.shp, region = "NAME_1")
# adm1.shp.f <- rename(adm1.shp.f,c("GADM4pt0_geom" = "geometry",
# "country" = "COUNTRY", "name_1" = "NAME_1"))
adm1.shp.f <- rename(adm1.shp.f,c("geom" = "geometry",
                                  "country" = "COUNTRY", "fs_name_1" = "NAME_1"))
adm1.shp.f <- as_tibble(adm1.shp.f)
# adm1.shp.f %<>% dplyr::select(country, name_1, GADM4pt0_geom)
adm1.shp.f %<>% dplyr::select(country, fs_name_1, geom)
adm1.shp.f$ISO2 <- countrycode(adm1.shp.f$country,"country.name","iso2c")

adm1.shp.f %<>%st_as_sfc(crs = st_crs(4326))

#ISO_shapes <- merge(ISO_sitefile_shapes, adm1.shp.f)
# ISO_shapes <- left_join(ISO_sitefile_shapes, adm1.shp.f)
# 
# 
# ISO_shapes$site_geom <- ISO_shapes$geom
# 
# ISO_shapes$geom <- ifelse(
#   ISO_shapes$iso3c == "GHA",
#   ISO_shapes$GADM4pt0_geom,
#   ISO_shapes$geom
# ) %>%
#   st_as_sfc(crs = st_crs(4326))
# 
# ISO_shapes$fs_name_1 <- ISO_shapes$name_1
# ISO_shapes$ISO2 <- countrycode(ISO_shapes$iso3c, "iso3c", "iso2c")

ua_tibble <- #as_tibble(
  merge(
    ua_end,
    adm1.shp.f,#ISO_shapes,
    by = c("ISO2", "fs_name_1")
    #)
  )

ua_tibble$geom <- ifelse(
  st_geometry_type(ua_tibble$geom) == "POLYGON",
  st_cast(ua_tibble$geom, "MULTIPOLYGON"),
  ua_tibble$geom
) %>%
  st_as_sfc(crs = st_crs(4326))

ua_sf_shape <- ua_tibble %>%
  mutate(
    geometry = map(geom,
                   ~ do.call(rbind, .) %>% # make each list a matrix
                     list() %>% # st_polygon() requires a list
                     st_multipolygon()
    )
  ) %>% 
  st_as_sf()


ua_sf_shape$country <- countrycode(ua_sf_shape$ISO2,"iso2c","country.name")

#cols4all::c4a_gui()
tm_min <- 0.2
tm_mid <- 0.8
tm_max <- 1
base_ncols <- 500
lo_ncols <- round((tm_mid-tm_min)*2*base_ncols)
hi_ncols <- round((tm_max-tm_mid)*2*base_ncols)
#lo_pal <- cols4all::c4a("seaborn.rocket", n = lo_ncols, reverse = FALSE)
#hi_pal <- cols4all::c4a("seaborn.mako", n = hi_ncols, reverse = TRUE)
lo_pal <- cols4all::c4a("carto.pink_yl", n = lo_ncols, reverse = TRUE)
#lo_pal <- cols4all::c4a("carto.sunset", n = lo_ncols, reverse = TRUE)
hi_pal <- cols4all::c4a("carto.blu_yl", n = hi_ncols, reverse = FALSE)

lo_pal <-  cols4all::c4a("matplotlib.spectral", n = lo_ncols*2, reverse = FALSE)
hi_pal <-  cols4all::c4a("matplotlib.spectral", n = hi_ncols*2, reverse = FALSE)
lo_pal <- lo_pal[1:lo_ncols]
hi_pal <- tail(hi_pal,hi_ncols)
custom_pal <- c(lo_pal, hi_pal)
time_pal <- cols4all::c4a("seaborn.rocket", n = base_ncols, reverse = FALSE)


lo_pal <-  cols4all::c4a("met.homer1", n = lo_ncols*2, reverse = FALSE)
hi_pal <-  cols4all::c4a("met.homer1", n = hi_ncols*2, reverse = FALSE)
lo_pal <- lo_pal[1:lo_ncols]
hi_pal <- tail(hi_pal,hi_ncols)
custom_pal <- c(lo_pal, hi_pal)
time_pal <- cols4all::c4a("carto.sunset_dark", n = base_ncols, reverse = FALSE)


lo_pal <-  viridis::turbo(n = lo_ncols*2 , direction = -1)
hi_pal <-  viridis::turbo(n = hi_ncols*2, direction  = -1)
lo_pal <- lo_pal[1:lo_ncols]
hi_pal <- tail(hi_pal, hi_ncols)

custom_pal <- c(lo_pal, hi_pal)
time_pal <- viridis::turbo(n = base_ncols*2, direction = -1)
time_pal <- time_pal[1:round(base_ncols*0.8)]
# time_pal <- viridis::turbo(n = 24, direction = -1)
# time_pal <- time_pal[1:8]



lo_pal <-  cols4all::c4a("matplotlib.spectral", n = lo_ncols*2, reverse = FALSE)
hi_pal <-  cols4all::c4a("matplotlib.spectral", n = hi_ncols*2, reverse = FALSE)
lo_pal <- lo_pal[1:lo_ncols]
hi_pal <- tail(hi_pal, hi_ncols)
custom_pal <- c(lo_pal, hi_pal)
time_pal <- cols4all::c4a("seaborn.rocket", n = base_ncols, reverse = FALSE)


if(is.null(ua_sf_shape$mn80_u)) {
  ua_sf_shape$mn80_u <- floor(ua_sf_shape$over80time_u)
}

tm_plt <- 
  ua_sf_shape %>%
  dplyr::filter(interval_yr == 2) %>%
  dplyr::filter(urbanicity == "urban") %>%
  # st_transform("ESRI:102022") %>%
  # cartogram_cont(weight = "pop") %>%
  tm_shape() +
  tm_polygons(
    fill = c("avg_u", "avg_a", "avg_uga", "mn80_u"),
    fill.scale = 
      list(tm_scale_continuous(values = custom_pal,
                               limits = c(tm_min, tm_max)),
           tm_scale_continuous(values = custom_pal,
                               limits = c(tm_min, tm_max)),
           tm_scale_continuous(values = custom_pal,
                               limits = c(tm_min, tm_max)),
           # tm_scale_continuous(values = time_pal,
           #                     limits = c(0, 8))
           # tm_scale_intervals(n = 7,
           #                    style = "cat",
           #                    #breaks = seq(0,8),
           #                    values = time_pal),
           tm_scale_discrete(values = time_pal)
      ),
    fill.legend = list(tm_legend(""),
                       tm_legend(""),
                       tm_legend(""),
                       tm_legend(""))
  ) +
  tm_facets_grid(columns = "country") + 
  tm_layout(frame = FALSE,
            legend.frame = FALSE,
            panel.labels = list(
              c("Mean use", "Mean access", "Mean use given access",
                "Months with use above 80%"),
              c("Burkina Faso", "Ghana", "Malawi", "Mali", "Mozambique", "Senegal"))
  )

  tmap_save(tm_plt, "use_access_maps_urban_2.pdf", width = 6, height = 6)

#c4a_gui()
time_pal <- cols4all::c4a("tol.incandescent", n = hi_ncols, reverse = TRUE)
tm_plt <- ua_sf_shape %>%
  dplyr::filter(interval_yr == 3) %>%
  dplyr::filter(urbanicity == "urban") %>%
  #st_transform("ESRI:102022") %>%
  #cartogram_cont(weight = "pop") %>%
  tm_shape() +
  tm_polygons(
    fill = c("over80time_u"),
    fill.scale = 
      list(
           tm_scale_continuous(values = time_pal,
                               limits = c(0, 8))
      ),
    fill.legend = list(tm_legend(""))
  ) +
  tm_facets_grid(columns = "country") + 
  tm_layout(frame = FALSE,
            legend.frame = FALSE,
            panel.labels = list(
              c("Time over 80"),
              c("Burkina Faso", "Ghana", "Malawi", "Mali", "Mozambique", "Senegal"))
  )



all_int_ua_sf_shape <- ua_sf_shape %>% 
  dplyr::filter(interval_yr == 3) %>%
  dplyr::select(ISO2, fs_name_1, ADM1, area, urbanicity, fs_area, country, pop,
                adm_pop, over80time_u, over80time_a, geometry, geom)

int2_ua_sf_shape <- ua_sf_shape %>% 
  dplyr::filter(interval_yr == 2) %>%
  dplyr::select(ISO2, fs_name_1, ADM1, area, urbanicity, fs_area, country, pop,
                adm_pop, over80time_u, over80time_a, geometry, geom,
                avg_u, avg_a, avg_uga)

int3_ua_sf_shape <- ua_sf_shape %>% 
  dplyr::filter(interval_yr == 3) %>%
  dplyr::select(ISO2, fs_name_1, ADM1, area, urbanicity, fs_area, country, pop,
                adm_pop, over80time_u, over80time_a, geometry, geom,
                avg_u, avg_a, avg_uga)

all_int_ua_sf_shape$avg_u2 <- int2_ua_sf_shape$avg_u
all_int_ua_sf_shape$avg_a2 <- int2_ua_sf_shape$avg_a
all_int_ua_sf_shape$avg_uga2 <- int2_ua_sf_shape$avg_uga
all_int_ua_sf_shape$avg_u3 <- int3_ua_sf_shape$avg_u
all_int_ua_sf_shape$avg_a3 <- int3_ua_sf_shape$avg_a
all_int_ua_sf_shape$avg_uga3 <- int3_ua_sf_shape$avg_uga

#all_int_ua_sf_shape <- all_int_ua_sf_shape[1:(dim(ua_sf_shape)[1]/2),]
# all_int_ua_sf_shape$avg_u2 <- ua_sf_shape$avg_u[ua_sf_shape$interval_yr == 2]
# all_int_ua_sf_shape$avg_a2 <- ua_sf_shape$avg_a[ua_sf_shape$interval_yr == 2]
# all_int_ua_sf_shape$avg_uga2 <- ua_sf_shape$avg_uga[ua_sf_shape$interval_yr == 2]
# all_int_ua_sf_shape$avg_u3 <- ua_sf_shape$avg_u[ua_sf_shape$interval_yr == 3]
# all_int_ua_sf_shape$avg_a3 <- ua_sf_shape$avg_a[ua_sf_shape$interval_yr == 3]
# all_int_ua_sf_shape$avg_uga3 <- ua_sf_shape$avg_uga[ua_sf_shape$interval_yr == 3]


#cols4all::c4a_gui()
tm_min <- 0.2
tm_mid <- 0.8
tm_max <- 1
base_ncols <- 500
lo_ncols <- round((tm_mid-tm_min)*2*base_ncols)
hi_ncols <- round((tm_max-tm_mid)*2*base_ncols)
#lo_pal <- cols4all::c4a("seaborn.rocket", n = lo_ncols, reverse = FALSE)
#hi_pal <- cols4all::c4a("seaborn.mako", n = hi_ncols, reverse = TRUE)
#lo_pal <- cols4all::c4a("carto.pink_yl", n = lo_ncols, reverse = TRUE)
lo_pal <- cols4all::c4a("carto.teal_rose", n = lo_ncols, reverse = TRUE)
hi_pal <- cols4all::c4a("carto.teal_grn", n = hi_ncols, reverse = TRUE)
time_pal <- cols4all::c4a("carto.temps", n = hi_ncols, reverse = TRUE)
custom_pal <- c(lo_pal, hi_pal)
#tm_plt <- 
  all_int_ua_sf_shape %>%
  dplyr::filter(urbanicity == "rural") %>%
  #st_transform("ESRI:102022") %>%
  #cartogram_cont(weight = "pop") %>%
  tm_shape() +
  tm_polygons(
    fill = c("avg_u3", "avg_u2", "avg_a3", "avg_a2", "avg_uga3", "avg_uga2", "over80time_u"),
    fill.scale = 
      list(tm_scale_continuous(values = lo_pal,
                               limits = c(tm_min, tm_mid)),
           tm_scale_continuous(values = lo_pal,
                               limits = c(tm_min, tm_mid)),
           tm_scale_continuous(values = lo_pal,
                               limits = c(tm_min, tm_mid)),
           tm_scale_continuous(values = lo_pal,
                               limits = c(tm_min, tm_mid)),
           tm_scale_continuous(values = lo_pal,
                               limits = c(tm_min, tm_mid)),
           tm_scale_continuous(values = tail(custom_pal,n=length(custom_pal)/2),
                               limits = c(0.6, tm_max)),
           tm_scale_continuous(values = time_pal,
                               limits = c(0, 8))
      ),
    fill.legend = list(tm_legend(""),
                       tm_legend(""),
                       tm_legend(""),
                       tm_legend(""),
                       tm_legend(""),
                       tm_legend(""),
                       tm_legend(""))
  ) +
  tm_facets_grid(columns = "country") + 
  tm_layout(frame = FALSE,
            legend.frame = FALSE,
            panel.labels = list(
              c("3-year mean use", "2-year mean use",
                "3-year mean access", "2-year mean access",
                "3-year mean use given access", "2-year mean use given access",
                "Months use above 80%"),
              c("Burkina Faso", "Ghana", "Malawi", "Mali", "Mozambique", "Senegal"))
  )













ggplot() +
  geom_sf(
    data = ua_sf_shape %>%
      filter(interval_yr == 2) %>%
      filter(ISO2 == "SN") %>%
      filter(urbanicity == "urban") #%>%
      #st_transform("ESRI:102022") %>%
    #cartogram_cont(weight = "pop")
    ,
    aes(
      group = name_1,
      fill = avg_u
    )
  )

custom_pal <- pals::ocean.curl(100)


ua_sf_shape %>%
  dplyr::filter(interval_yr == 2) %>%
  dplyr::filter(urbanicity == "rural") %>%
#st_transform("ESRI:102022") %>%
#cartogram_cont(weight = "pop") %>%
tm_shape() +
  tm_polygons(fill = 'avg_u',
              #title = 'PKO targeting\nevents (logged)',
              style = 'cont') +                  # continuous variable
  tm_facets('ISO2',nrow = 1) +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "bottom",
            legend.position = c(0.25, 0.25))
  tm_layout(legend.outside.position =  "bottom", # legend outside below
            legend.position = c(.8, 1.1)
            )        # manually position legend

  ua_sf_shape$country <- countrycode(ua_sf_shape$ISO2,"iso2c","country.name")
#c4a_gui()
tm_min <- 0.2
tm_mid <- 0.8
tm_max <- 1
base_ncols <- 500
lo_ncols <- round((tm_mid-tm_min)*2*base_ncols)
hi_ncols <- round((tm_max-tm_mid)*2*base_ncols)
# lo_pal <- rev(pals::brewer.purples(lo_ncols))
# hi_pal <- pals::brewer.bugn(hi_ncols)
#lo_pal <- pals::ocean.dense(lo_ncols)
#hi_pal <- rev(pals::ocean.algae(hi_ncols))
#lo_pal <- cols4all::c4a("carto.magenta", n = lo_ncols, reverse = TRUE)
#hi_pal <- cols4all::c4a("carto.mint", n = hi_ncols, reverse = FALSE)
lo_pal <- cols4all::c4a("seaborn.rocket", n = lo_ncols, reverse = FALSE)
hi_pal <- cols4all::c4a("seaborn.mako", n = hi_ncols, reverse = TRUE)
custom_pal <- c(lo_pal, hi_pal)
time_pal <- cols4all::c4a("tol.incandescent", n = hi_ncols, reverse = TRUE)
# base_left_pal <- pals::ocean.curl(left_ncols)
# base_right_pal <- pals::ocean.curl(right_ncols)
# custom_pal <- rev(c(base_left_pal[1:(left_ncols/2)],
#                     base_right_pal[(right_ncols/2):right_ncols]))
tm_plt <- ua_sf_shape %>%
    dplyr::filter(interval_yr == 3) %>%
    dplyr::filter(urbanicity == "urban") %>%
    #st_transform("ESRI:102022") %>%
    #cartogram_cont(weight = "pop") %>%
    tm_shape() +
    tm_polygons(
      fill = c("avg_u", "avg_a", "avg_uga", "over80time_u"),
      #style = 'cont',
      fill.scale = 
        list(tm_scale_continuous(values = custom_pal,
                                 limits = c(tm_min, tm_max)),
             tm_scale_continuous(values = custom_pal,
                                 limits = c(tm_min, tm_max)),
             tm_scale_continuous(values = custom_pal,
                                 limits = c(tm_min, tm_max)),
             tm_scale_continuous(values = time_pal,
                                 limits = c(0, 8))
        ),
      # fill.legend = list(tm_legend("mean\nuse"),
      #                    tm_legend("mean\naccess"),
      #                    tm_legend("mean use\ngiven access"))
      fill.legend = list(tm_legend(""),
                         tm_legend(""),
                         tm_legend(""),
                         tm_legend(""))
      ) +
        # list(tm_scale_intervals(values = "-brewer.reds"),
        #      #tm_scale_intervals(values = "brewer.purples"),
        #      tm_scale_intervals(values = "brewer.blues")
        #      )) +
    tm_facets_grid(columns = "country") + 
  tm_layout(frame = FALSE,
            legend.frame = FALSE,
            #legend.text.size = 0.5,
            panel.labels = list(
              c("mean use", "mean access", "mean use given access"),
              c("Burkina Faso", "Ghana", "Malawi", "Mali", "Mozambique", "Senegal"))
            )

tmap_save(tm_plt, "tmap.pdf", width = 8, height = 6)



g <- purrr::map(ua_sf_shape$ISO2,
                function(x) {
                  ggplot() +
                    geom_sf(data = filter(ua_sf_shape, ISO2 == x)) +
                    guides(fill = FALSE) +
                    ggtitle(x)
                })

g2 <- cowplot::plot_grid(plotlist = g)



## create maps in separate plots, force common scale between them
maps <- purrr::map(.x = SSA_ISO2, 
            .f = function(x) ua_sf_shape %>% 
              dplyr::filter(ISO2 == x) %>% 
              #st_join(acled) %>% 
              #group_by(NAME_0, NAME_1, NAME_2) %>% 
              #summarize(attacks = log1p(sum(!is.na(event_id_cnty)))) %>% 
              ggplot(aes(fill = avg_u)) +
              geom_sf(lwd = NA) +
              scale_fill_continuous(name = 'PKO targeting\nevents (logged)') +
              theme_rw() +
              theme(axis.text = element_blank(),
                    axis.ticks = element_blank()))







ggsave(paste0("test",".png"), bg = "white",
       w = 8, h = 2, dpi = 250)










ua_tibble %<>% filter(ISO2 != "GH")



ISO_shapes <- merge(
  x = ISO_sitefile_shapes,
  y = adm1.shp.f[ , c("country", "name_1", "GADM4pt0_geom")],
  by = c("country", "name_1", all.x=TRUE)
)

ISO_shapes$GADM4pt0_geom <- left_join()

ISO_shapes$geom <- ifelse(ISO_shapes$iso3c == "GHA",
                          )

ISO_shapes$fs_name_1 <- ISO_shapes$name_1

ua_tibble <- as_tibble(
  merge(
    ua_end,
    ISO_shapes,
    by = "fs_name_1"
  )
)

ua_tibble$geom <- ifelse(
  st_geometry_type(ua_tibble$geom) == "POLYGON",
  st_cast(ua_tibble$geom, "MULTIPOLYGON"),
  ua_tibble$geom
)


ua_sf_shape <- ua_tibble %>%
  mutate(
    geometry = map(geom,
                   ~ do.call(rbind, .) %>% # make each list a matrix
                     list() %>% # st_polygon() requires a list
                     st_multipolygon()
    )
  ) %>% 
  st_as_sf()


GH.adm1.shp <- geodata::gadm(country = "GH", level = 1, 
                          path = tempdir(),
                          version = "4.0")
GH.adm1.shp.f <- sf::st_as_sf(GH.adm1.shp, region = "NAME_1")


ua_sf_shape %>% 
  st_transform("ESRI:102022") %>%
  group_by(ISO2)
  tm_shape() +
  tm_polygons(col = 'avg_u',
              title = 'PKO targeting\nevents (logged)',
              style = 'cont') +                  # continuous variable
  tm_facets('') +
  tm_layout(legend.outside.position =  "bottom", # legend outside below
            legend.position = c(.8, 1.1))        # manually position legend




shape_test <- ua_sf_shape %>% filter(ISO2 == "BF")
shape_test <- shape_test %>% st_set_crs(4979)

baseline_carto_shapes <- st_transform(shape_test) %>%
  cartogram_cont(weight = "pfpr")

carto_test <- ua_sf_shape %>%
  filter(ISO2 == "BF" & !is.na(EIR_fit)) %>%
  cartogram_cont(weight = "fs_area_id")

  #cartogram_cont(weight = "EIR_fit", itermax = 500, maxSizeError = 1.1)

ggplot() +
  geom_sf(
    data = ua_sf_shape %>%
      filter(ISO2 == "MW") #%>%
      #cartogram_cont(weight = "fs_area_id")
    ,
    aes(
      group = name_1,
        fill = avg_u
      )
  )



ua_df <- fs_id_link
ua_df <- left_join(
  fs_id_link
)



site_file$metadata$boundary

BF_site$metadata$boundary

rur_prev24 <- BF_site$prevalence %>% dplyr::filter(year == 2024 & urban_rural == "rural")

tib_shpe <- as_tibble(merge(rur_prev24, BF_site$shape$level_1, by = "name_1"))

sf_shape <- tib_shpe %>%
  mutate(geometry = map(geom,
                        ~ do.call(rbind, .) %>% # make each list a matrix
                          list() %>% # st_polygon() requires a list
                          st_polygon()
  )
  ) %>% 
  st_as_sf()

baseline_carto_shapes <- st_transform(sf_shape, "ESRI:102022") %>%
  cartogram_cont(weight = "pfpr")

ggplot() +
  geom_sf_pattern(data = baseline_carto_shapes,
                  aes(group = name_1,
                      fill = pfpr))