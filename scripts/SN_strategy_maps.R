library(tidyverse)
library(raster)
library(sf)
library(ggplot2)
library(viridis)
library(cartogram)
library(ggthemes)

# Read in SN strategies
only0 <- readRDS("SNonly0.rds")
only3 <- readRDS("SNonly3.rds")
pyrroleD <- readRDS("SNpyrroleD.rds")
pyrrole2 <- readRDS("SNpyrrole2.rds")

# Fetch_baselines
#record_baselines()

# Append cases averted and costs
only0 %<>% append_cases_averted(baseline_df = only0)
only3 %<>% append_cases_averted(baseline_df = only0)
pyrroleD %<>% append_cases_averted(baseline_df = only0)
pyrrole2 %<>% append_cases_averted(baseline_df = only0)

# Urban-rural weightings
append_admin_pops_cases_costs <- function(net_df = NULL) {
  reps <- max(net_df$sample_index)
  admin_pops <- net_df %>%
    group_by(fs_name_1) %>%
    dplyr::summarise(admin_pop = sum(pop)/reps,
                     admin_cases = sum(pred_ann_infect)/reps,
                     admin_cost = sum(avg_ann_net_cost)/reps,
                     admin_cases_per_USD = sum(pred_ann_infect)/sum(avg_ann_net_cost))
  net_df %<>% left_join(admin_pops)
}
only0 %<>% append_admin_pops_cases_costs()
only3 %<>% append_admin_pops_cases_costs()
pyrroleD %<>% append_admin_pops_cases_costs()
pyrrole2 %<>% append_admin_pops_cases_costs()

append_adm_weighted_data <- function(net_df = NULL) {
  net_df$urbanicity_weighting <- net_df$pop / net_df$admin_pop
  net_df$adm_cases <- 
  net_df
}
only0 %<>% append_adm_weighted_data()
only3 %<>% append_adm_weighted_data()
pyrroleD %<>% append_adm_weighted_data()
pyrrole2 %<>% append_adm_weighted_data()



# Admin cases


# # Summarise
# summarise_admin_costs_cases <- function(net_df = NULL){
#   net_df %>%
#     group_by(fs_area) %>%
#     dplyr::summarise(fs_name_1,
#                      pop,
#                      urbanicity_weighting,
#                      cases_avert_mean = mean(cases_averted),
#                      cases_avert_95LB = quantile(cases_averted,
#                                                  probs = 0.025,
#                                                  names = FALSE),
#                      cases_avert_95UB = quantile(cases_averted,
#                                                  probs = 0.975,
#                                                  names = FALSE),
#                      cases_avert_per_cap_mean = mean(cases_averted_per_capita),
#                      cases_avert_per_cap_95LB = quantile(cases_averted_per_capita,
#                                                          probs = 0.025,
#                                                          names = FALSE),
#                      cases_avert_per_cap_95UB = quantile(cases_averted_per_capita,
#                                                          probs = 0.975,
#                                                          names = FALSE),
#                      cost_mean = mean(avg_ann_net_cost),
#                      cost_95LB = quantile(avg_ann_net_cost,
#                                           probs = 0.025,
#                                           names = FALSE),
#                      cost_95UB = quantile(avg_ann_net_cost,
#                                           probs = 0.975,
#                                           names = FALSE),
#                      cases_avert_per_USD_mean = mean(cases_averted_per_USD),
#                      cases_avert_per_USD_95LB = quantile(cases_averted_per_USD,
#                                                          probs = 0.025,
#                                                          names = FALSE),
#                      cases_avert_per_USD_95UB = quantile(cases_averted_per_USD,
#                                                          probs = 0.975,
#                                                          names = FALSE)
#     ) %>%
#     unique()
# }

# Summarise
summarise_admin_costs_cases <- function(net_df = NULL){
  net_df %<>%
    group_by(fs_area) %>%
    dplyr::summarise(fs_name_1,
                     admin_pop,
                     admin_cases,
                     admin_cases_per_cap = admin_cases / admin_pop,
                     admin_cost,
                     admin_cases_per_USD,
                     # cases_avert_mean = mean(cases_averted),
                     # cases_avert_per_cap_mean = mean(cases_averted_per_capita),
                     # cost_mean = mean(avg_ann_net_cost),
                     # cases_avert_per_USD_mean = mean(cases_averted_per_USD)
    ) %>%
    unique()
  #net_df %>% select(!fs_area)
}

pyrroleD_admur_sum <- pyrroleD %>% summarise_admin_costs_cases# %>% as.data.frame()
pyrroleD_admur_sum <- pyrroleD_admur_sum[,2:dim(pyrroleD_admur_sum)[2]] %>% unique()
pyrrole2_admur_sum <- summarise_admin_costs_cases(pyrrole2)
only3_admur_sum <- summarise_admin_costs_cases(only3)

pyrrole2_adm_sum <- 

# get shapefiles
adm1.shp <- raster::getData("GADM", country = "SEN", level = 1)
adm1.shp.f <- sf::st_as_sf(adm1.shp, region = "NAME_1") #fortify
#adm1.shp.f <- tidy(adm1.shp, region = "NAME_1") #fortify

adm1.shp.f %<>%
  dplyr::rename(fs_name_1 = NAME_1)

cases_averted_lgd <- "Cases averted\nper 1,000"

pyrrole2_shapes <- merge(adm1.shp.f, as_tibble(pyrroleD_admur_sum), by = "fs_name_1")
pyrrole2_carto_shapes <- st_transform(pyrrole2_shapes, "ESRI:102022") %>%
  cartogram_cont(weight = "admin_pop")


pyrrole2_cartog <- cartogram_ncont(st_transform(pyrrole2_shapes), weight = "pop")

test <- merge(st_transform(adm1.shp), pyrrole2_sum, by = "fs_name_1")

ggplot() +
  geom_sf(data = pyrrole2_carto_shapes,
          aes(group = fs_name_1,
              fill = admin_cases_per_cap * 1e3)) +
  # geom_sf(data = pyrrole2_carto_shapes,
  #         aes(group = fs_name_1,
  #             fill = cases_avert_mean)) +
  scale_fill_viridis(option = "mako") +
  guides(fill = guide_colorbar(title = "")) +
  theme_map()
# ggsave(paste0("pboprevred",timestamp,".png"), bg = "white",
#        w = 15, h = 5, dpi = 450)


# scale_fill_viridis(begin = .5,
#                    end = 1,
#                    direction = -1,
#                    option = "turbo") +

# Generate adm-urbanicity sum
adm_urbanicity_sum pyrrole

# 
only3_country_cost <- sum(only3_sum$cost_mean)
pyrroleD_country_cost <- sum(pyrroleD$cost_mean)

priority_ids <- order(pyrrole2_sum$cases_avert_per_cap_mean,
                                     decreasing = TRUE)

pyrroleD2mix_cost <- pyrroleD_country_cost
for 








total_test <- only3 %>%
  group_by(sample_index) %>%
  dplyr::summarise(total_cost = sum(avg_ann_net_cost))

mtcars %>%
  group_by(cyl) %>%
  dplyr::summarise(mean = mean(disp), n = n())