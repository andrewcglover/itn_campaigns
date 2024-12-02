library(tidyverse)
library(raster)
library(sf)
library(ggplot2)
library(viridis)
library(cartogram)
library(ggthemes)
library(cowplot)
library(grid)
library(gridExtra)
library(paletteer)
library(geodata)
library(ggpattern)

# Select data
uncosted_data <- sim_18NOV24_data_uncosted

uncosted_data$strategy_id <- paste(uncosted_data$net_strategy,
                                   uncosted_data$net_name,
                                   uncosted_data$mass_int)
uncosted_data$ann_infect_percap <- uncosted_data$annual_infections / 5e4
uncosted_data$avg_ann_tot_cost <- uncosted_data$avg_ann_npc_cost * uncosted_data$pop
#uncosted_data$cases_avert_per_usd <- uncosted_data$ann_infect_percap / uncosted_data$avg_ann_npc_cost

# Country specific
ISO2_map <- "SN"
net_type <- "pyrethroid-only"
ctry_uncosted <- uncosted_data %>% filter(ISO2 == ISO2_map)
no_future_ann_infect_percap <- ctry_uncosted$ann_infect_percap[ctry_uncosted$routine_baseline == 1]
ctry_uncosted$no_future_ann_infect_percap <- rep(no_future_ann_infect_percap,
                                                 length.out = dim(ctry_uncosted)[1])
ctry_uncosted$add_avert_percap <- ctry_uncosted$no_future_ann_infect_percap - ctry_uncosted$ann_infect_percap
ctry_uncosted$cases_avert_per_usd <- ctry_uncosted$add_avert_percap / ctry_uncosted$avg_ann_npc_cost


urban_pop <- ctry_uncosted %>%
  filter(urbanicity == "urban") %>%
  group_by(fs_name_1) %>%
  dplyr::summarise(urban_pop = mean(pop)) %>%
  as.data.frame()

rural_pop <- ctry_uncosted %>%
  filter(urbanicity == "rural") %>%
  group_by(fs_name_1) %>%
  dplyr::summarise(rural_pop = mean(pop)) %>%
  as.data.frame()

pop_summary <- urban_pop
pop_summary$rural_pop <- rural_pop$rural_pop
pop_summary$all_pop <- pop_summary$urban_pop + pop_summary$rural_pop
pop_summary$urban_weight <- pop_summary$urban_pop / pop_summary$all_pop
pop_summary$rural_weight <- pop_summary$rural_pop / pop_summary$all_pop

N_regions <- dim(pop_summary)[1]

ctry_summary <- ctry_uncosted %>%
  dplyr::group_by(fs_area,
                  fs_name_1,
                  urbanicity,
                  strategy_id,
                  pop,
                  net_name,
                  net_strategy,
                  mass_int) %>%
  dplyr::summarise(lo_cost = quantile(avg_ann_tot_cost,
                                      probs = 0.025),
                   mid_cost = quantile(avg_ann_tot_cost,
                                       probs = 0.5),
                   hi_cost = quantile(avg_ann_tot_cost,
                                      probs = 0.975),
                   lo_avertperusd = quantile(cases_avert_per_usd,
                                             probs = 0.025),
                   mid_avertperusd = quantile(cases_avert_per_usd,
                                              probs = 0.5),
                   hi_avertperusd = quantile(cases_avert_per_usd,
                                             probs = 0.975),
                   lo_pfpr = quantile(avg_pfpr,
                                      probs = 0.025),
                   mid_pfpr = quantile(avg_pfpr,
                                       probs = 0.5),
                   hi_pfpr = quantile(avg_pfpr,
                                      probs = 0.975)) %>%
  as.data.frame()

baseline_summary <- ctry_summary %>%
  filter(strategy_id == "uncosted pyrethroid-only 3")

baseline_lo_cost <- sum(baseline_summary$lo_cost)
baseline_mid_cost <- sum(baseline_summary$mid_cost)
baseline_hi_cost <- sum(baseline_summary$hi_cost)

net_summary <- ctry_summary %>%
  filter(net_name == net_type, net_strategy == "uncosted")

int2_summary <- net_summary %>%
  filter(mass_int == 2) %>%
  arrange(desc(mid_avertperusd))

int4_summary <- net_summary %>%
  filter(mass_int == 4) %>%
  arrange(desc(mid_avertperusd))

running_cost <- sum(int4_summary$mid_cost)

urban_deprioritised <- c()
rural_deprioritised <- c()
for (i in 1:dim(int2_summary)[1]) {
  candidate_cost <- running_cost - int4_summary$mid_cost[i] + int2_summary$mid_cost[i]
  if (candidate_cost <= baseline_hi_cost) {
    running_cost <- candidate_cost
  } else {
    if(int2_summary$urbanicity[i] == "urban") {
      urban_deprioritised <- c(urban_deprioritised, int2_summary$fs_name_1[i])
    } else {
      rural_deprioritised <- c(rural_deprioritised, int2_summary$fs_name_1[i])
    }
  }
}
  

baseline_adm <- pop_summary
baseline_adm$mid_pfpr <- rep(NA, N_regions)
for (i in 1:N_regions) {
  urban_id <- which((baseline_summary$fs_name_1 == baseline_adm$fs_name_1[i]) &
      (baseline_summary$urbanicity == "urban"))
  rural_id <- which((baseline_summary$fs_name_1 == baseline_adm$fs_name_1[i]) &
                      (baseline_summary$urbanicity == "rural"))
  baseline_adm$mid_pfpr[i] <- (baseline_summary$mid_pfpr[urban_id] * baseline_adm$urban_weight[i]) + (baseline_summary$mid_pfpr[rural_id] * baseline_adm$rural_weight[i])
}




# get shapefiles
#adm1.shp <- raster::getData("GADM", country = "SEN", level = 1)
adm1.shp <- geodata::gadm(country = "SEN", level = 1, 
                          path = tempdir(),
                          version = "latest")
adm1.shp.f <- sf::st_as_sf(adm1.shp, region = "NAME_1") #fortify
#adm1.shp.f <- tidy(adm1.shp, region = "NAME_1") #fortify

adm1.shp.f %<>%
  dplyr::rename(fs_name_1 = NAME_1)

cases_averted_lgd <- "Cases averted\nper 1,000"

baseline_shapes <- merge(adm1.shp.f, as_tibble(baseline_adm), by = "fs_name_1")
baseline_carto_shapes <- st_transform(baseline_shapes, "ESRI:102022") %>%
  cartogram_cont(weight = "all_pop")

ggplot() +
  geom_sf(data = baseline_carto_shapes,
          aes(group = fs_name_1,
              fill = mid_pfpr*100,
              alpha = mid_pfpr*100),
  ) +
  scale_fill_viridis(option = "turbo",
                     begin = 0,
                     end = 1,
                     direction = 1,
                     #breaks = cases_avert_per_1000_breaks,
                     #labels = cases_avert_per_1000_breaks,
                     limits = c(0,60)
  ) +
  guides(fill = guide_colorbar(title = expression(P*italic(f)*PR[2-10]*" (%)"),
                                title.position="top",
                                barwidth = 8,
                                barheight = 0.5)) +
  #guides(fill = "none") +
  theme_map() +
  theme(legend.position="bottom",
        legend.title.align=0.5)















cost_mid_df <- ctry_uncosted %>%
  dplyr::group_by(fs_area, strategy_id) %>%
  dplyr::summarise(mid_cost = quantile(avg_ann_tot_cost, probs = 0.5))
cost_hi_df <- ctry_uncosted %>%
  dplyr::group_by(fs_area, strategy_id) %>%
  dplyr::summarise(mid_cost = quantile(avg_ann_tot_cost, probs = 0.975))







baseline_strategy <- ctry_uncosted %>%
  filter(strategy_id == "uncosted pyrethroid-only 3")


baseline_cost_lo_df <- baseline_strategy %>%
  dplyr::group_by(fs_area) %>%
  dplyr::summarise(mid_cost = quantile(avg_ann_tot_cost, probs = 0.025))
baseline_cost_lo <- sum(baseline_cost_lo_df[,2])

baseline_cost_mid_df <- baseline_strategy %>%
  dplyr::group_by(fs_area) %>%
  dplyr::summarise(mid_cost = quantile(avg_ann_tot_cost, probs = 0.5))
baseline_cost_mid <- sum(baseline_cost_mid_df[,2])

baseline_cost_hi_df <- baseline_strategy %>%
  dplyr::group_by(fs_area) %>%
  dplyr::summarise(mid_cost = quantile(avg_ann_tot_cost, probs = 0.975))
baseline_cost_hi <- sum(baseline_cost_hi_df[,2])


ctry_uncosted %>%
  dplyr::group_by(fs_area,
                  strategy_id,
                  net_strategy,
                  net_name,
                  mass_int) %>%
  dplyr::summarise(mean_cost = quantile(avg_ann_tot_cost, probs = 0.5))

ctry_only_2 <- ctry_uncosted %>%
  filter(strategy_id == "uncosted pyrethroid-only 2")
ctry_only_4 <- ctry_uncosted %>%
  filter(strategy_id == "uncosted pyrethroid-only 4")
ctry_pbo_2 <- ctry_uncosted %>%
  filter(strategy_id == "uncosted pyrethroid-PBO 2")
ctry_pbo_4 <- ctry_uncosted %>%
  filter(strategy_id == "uncosted pyrethroid-PBO 4")
ctry_pyrrole_2 <- ctry_uncosted %>%
  filter(strategy_id == "uncosted pyrethroid-pyrrole 2")
ctry_pyrrole_4 <- ctry_uncosted %>%
  filter(strategy_id == "uncosted pyrethroid-pyrrole 4")




dakar_test <- ctry_uncosted[ctry_uncosted$fs_area=="SN Dakar urban",]
ggplot(data = dakar_test,
       aes(x = avg_ann_tot_cost,
           colour = net_name,
           fill = net_name,
           linetype = mass_int)) +
  geom_density(alpha = 0.2)


# Read in SN strategies
only0 <- readRDS("SNonly0.rds")
only3 <- readRDS("SNonly3.rds")
pyrrole3 <- readRDS("SNpyrrole3.rds")
pyrrole3c <- readRDS("SNpyrrole3c.rds")
pyrrole2 <- readRDS("SNpyrrole2.rds")
pyrrole2c <- readRDS("SNpyrrole2c.rds")
pyrroleD <- readRDS("SNpyrroleD.rds")

# Deprioritised admins
deprioritised2_adms <- c("Fatick",
                        "Saint-Louis",
                        "Thiès",
                        "Dakar")
deprioritised3_adms <- c("Dakar")

# Mixed strategy
pyrrole2m <- pyrrole2
deprioritised2_IDs <-which(pyrrole2m$fs_name_1 %in% deprioritised2_adms)
pyrrole2m[deprioritised2_IDs,] <- pyrroleD[deprioritised2_IDs,]
pyrrole3m <- pyrrole3
deprioritised3_IDs <-which(pyrrole3m$fs_name_1 %in% deprioritised3_adms)
pyrrole3m[deprioritised3_IDs,] <- pyrroleD[deprioritised3_IDs,]

# Fetch_baselines
#record_baselines()

# Append cases averted and costs
only0 %<>% append_cases_averted(baseline_df = only0)
only3 %<>% append_cases_averted(baseline_df = only0)
pyrrole3 %<>% append_cases_averted(baseline_df = only0)
pyrrole3c %<>% append_cases_averted(baseline_df = only0)
pyrrole2 %<>% append_cases_averted(baseline_df = only0)
pyrrole2c %<>% append_cases_averted(baseline_df = only0)
pyrroleD %<>% append_cases_averted(baseline_df = only0)
pyrrole2m %<>% append_cases_averted(baseline_df = only0)
pyrrole3m %<>% append_cases_averted(baseline_df = only0)

# Append default baseline cases
only3$default_baseline_cases <- only0$pred_ann_infect
pyrrole3$default_baseline_cases <- only0$pred_ann_infect
pyrrole3c$default_baseline_cases <- only0$pred_ann_infect
pyrrole2$default_baseline_cases <- only0$pred_ann_infect
pyrrole2c$default_baseline_cases <- only0$pred_ann_infect
pyrroleD$default_baseline_cases <- only0$pred_ann_infect
pyrrole2m$default_baseline_cases <- only0$pred_ann_infect
pyrrole3m$default_baseline_cases <- only0$pred_ann_infect

# Default annual averted
only3$default_cases_averted <- only3$default_baseline_cases - only3$pred_ann_infect
pyrrole3$default_cases_averted <- pyrrole3$default_baseline_cases - pyrrole3$pred_ann_infect
pyrrole3c$default_cases_averted <- pyrrole3c$default_baseline_cases - pyrrole3c$pred_ann_infect
pyrrole2$default_cases_averted <- pyrrole2$default_baseline_cases - pyrrole2$pred_ann_infect
pyrrole2c$default_cases_averted <- pyrrole2c$default_baseline_cases - pyrrole2c$pred_ann_infect
pyrroleD$default_cases_averted <- pyrroleD$default_baseline_cases - pyrroleD$pred_ann_infect
pyrrole2m$default_cases_averted <- pyrrole2m$default_baseline_cases - pyrrole2m$pred_ann_infect
pyrrole3m$default_cases_averted <- pyrrole3m$default_baseline_cases - pyrrole3m$pred_ann_infect

# Country cases averted
all_cases_averted <- function(net_df) {
  N_admins <- length(unique(net_df$fs_area))
  N_reps <- max(net_df$sample_index)
  all_cases_averted_samples <- c()
  for (i in 1:N_reps) {
    ca_count <- 0
    for (j in 1:N_admins) {
      ca_count <- ca_count + net_df$default_cases_averted[i+(j-1)*N_reps]
    }
    all_cases_averted_samples[i] <- ca_count
  }
  return(all_cases_averted_samples)
}
all_cost <- function(net_df) {
  N_admins <- length(unique(net_df$fs_area))
  N_reps <- max(net_df$sample_index)
  var_samples <- c()
  for (i in 1:N_reps) {
    var_count <- 0
    for (j in 1:N_admins) {
      var_count <- var_count + net_df$avg_ann_net_cost[i+(j-1)*N_reps]
    }
    var_samples[i] <- var_count
  }
  return(var_samples)
}

all_ca <- pyrrole3m %>% all_cases_averted()
mean(all_ca)
quantile(all_ca, probs = c(0.025, 0.975), names = FALSE)

all_ca <- pyrrole3m %>% all_cost()
mean(all_ca)
quantile(all_ca, probs = c(0.025, 0.975), names = FALSE)

# Urban-rural weightings
append_admin_pops_cases_costs <- function(net_df = NULL) {
  reps <- max(net_df$sample_index)
  admin_pops <- net_df %>%
    group_by(fs_name_1) %>%
    dplyr::summarise(admin_pop = sum(pop)/reps,
                     admin_baseline_cases = sum(default_baseline_cases)/reps,
                     admin_cases = sum(pred_ann_infect)/reps,
                     admin_cases_avert = sum(default_cases_averted)/reps,
                     admin_cases_avert_mean = mean(default_cases_averted),
                     admin_cost = sum(avg_ann_net_cost)/reps,
                     admin_cases_avert_per_USD = sum(default_cases_averted)/sum(avg_ann_net_cost))
  net_df %<>% left_join(admin_pops)
}
#only0 %<>% append_admin_pops_cases_costs()
only3 %<>% append_admin_pops_cases_costs()
pyrrole3 %<>% append_admin_pops_cases_costs()
pyrrole3c %<>% append_admin_pops_cases_costs()
pyrrole2 %<>% append_admin_pops_cases_costs()
pyrrole2c %<>% append_admin_pops_cases_costs()
pyrroleD %<>% append_admin_pops_cases_costs()
pyrrole2m %<>% append_admin_pops_cases_costs()
pyrrole3m %<>% append_admin_pops_cases_costs()

# Append cases averted CrI
# append_cases_averted_CrI <- function(net_df = NULL){
#   reps <- max(net_df$sample_index)
#   fs_name_1s_here <- unique(net_df$fs_name_1)
#   new_net_df <- data.frame()
#   for (i in 1:length(fs_name_1s_here)) {
#     admin_df <- net_df %>% filter(fs_name_1 == fs_name_1s_here[i])
#     cases_averted_here <- c()
#     for (j in 1:reps) {
#       ca1 <- admin_df$default_cases_averted[j]# * admin_df$pop[j] / admin_df$admin_pop[j]
#       ca2 <- admin_df$default_cases_averted[j+reps]# * admin_df$pop[j+reps] / admin_df$admin_pop[j+reps]
#       cases_averted_here[j] <- ca1 + ca2
#     }
#     caLB <- quantile(cases_averted_here, probs = 0.025, names = FALSE)
#     caUB <- quantile(cases_averted_here, probs = 0.975, names = FALSE)
#     admin_df$cases_averted_LB <- rep(caLB, reps * 2)
#     admin_df$cases_averted_UB <- rep(caUB, reps * 2)
#     new_net_df <- rbind.data.frame(new_net_df, admin_df)
#   }
#   return(new_net_df)
# }

# Summarise
summarise_admin_costs_cases <- function(net_df = NULL){
  net_df %<>%
    group_by(fs_area) %>%
    dplyr::summarise(fs_name_1,
                     admin_pop,
                     admin_cases,
                     admin_cases_avert,
                     admin_cases_avert_per_cap = admin_cases_avert / admin_pop,
                     admin_cost,
                     admin_cases_avert_per_USD,
                     # cases_avert_mean = mean(cases_averted),
                     # cases_avert_per_cap_mean = mean(cases_averted_per_capita),
                     # cost_mean = mean(avg_ann_net_cost),
                     # cases_avert_per_USD_mean = mean(cases_averted_per_USD)
    ) %>%
    unique()
  #net_df %>% select(!fs_area)
}

only3_admur_sum <- only3 %>% summarise_admin_costs_cases
only3_admur_sum <- only3_admur_sum[,2:dim(only3_admur_sum)[2]] %>% unique()
pyrrole3_admur_sum <- pyrrole3 %>% summarise_admin_costs_cases
pyrrole3_admur_sum <- pyrrole3_admur_sum[,2:dim(pyrrole3_admur_sum)[2]] %>% unique()
pyrrole3c_admur_sum <- pyrrole3c %>% summarise_admin_costs_cases
pyrrole3c_admur_sum <- pyrrole3c_admur_sum[,2:dim(pyrrole3c_admur_sum)[2]] %>% unique()
pyrrole2_admur_sum <- pyrrole2 %>% summarise_admin_costs_cases
pyrrole2_admur_sum <- pyrrole2_admur_sum[,2:dim(pyrrole2_admur_sum)[2]] %>% unique()
pyrrole2c_admur_sum <- pyrrole2c %>% summarise_admin_costs_cases
pyrrole2c_admur_sum <- pyrrole2c_admur_sum[,2:dim(pyrrole2c_admur_sum)[2]] %>% unique()
pyrroleD_admur_sum <- pyrroleD %>% summarise_admin_costs_cases
pyrroleD_admur_sum <- pyrroleD_admur_sum[,2:dim(pyrroleD_admur_sum)[2]] %>% unique()
pyrrole2m_admur_sum <- pyrrole2m %>% summarise_admin_costs_cases
pyrrole2m_admur_sum <- pyrrole2m_admur_sum[,2:dim(pyrrole2m_admur_sum)[2]] %>% unique()
pyrrole3m_admur_sum <- pyrrole3m %>% summarise_admin_costs_cases
pyrrole3m_admur_sum <- pyrrole3m_admur_sum[,2:dim(pyrrole3m_admur_sum)[2]] %>% unique()
#pyrrole2_admur_sum <- summarise_admin_costs_cases(pyrrole2)
#only3_admur_sum <- summarise_admin_costs_cases(only3)

depri2_admur_sum <- pyrrole2m_admur_sum %>% filter(fs_name_1 %in% deprioritised2_adms)
depri3_admur_sum <- pyrrole3m_admur_sum %>% filter(fs_name_1 %in% deprioritised3_adms)

only3_admur_sum$add_admin_cases_avert_per_cap <- only3_admur_sum$admin_cases_avert_per_cap - only3_admur_sum$admin_cases_avert_per_cap
pyrrole2_admur_sum$add_admin_cases_avert_per_cap <- pyrrole2_admur_sum$admin_cases_avert_per_cap - only3_admur_sum$admin_cases_avert_per_cap
pyrrole2c_admur_sum$add_admin_cases_avert_per_cap <- pyrrole2c_admur_sum$admin_cases_avert_per_cap - only3_admur_sum$admin_cases_avert_per_cap
pyrrole2m_admur_sum$add_admin_cases_avert_per_cap <- pyrrole2m_admur_sum$admin_cases_avert_per_cap - only3_admur_sum$admin_cases_avert_per_cap
pyrrole3_admur_sum$add_admin_cases_avert_per_cap <- pyrrole3_admur_sum$admin_cases_avert_per_cap - only3_admur_sum$admin_cases_avert_per_cap
pyrrole3c_admur_sum$add_admin_cases_avert_per_cap <- pyrrole3c_admur_sum$admin_cases_avert_per_cap - only3_admur_sum$admin_cases_avert_per_cap
pyrrole3m_admur_sum$add_admin_cases_avert_per_cap <- pyrrole3m_admur_sum$admin_cases_avert_per_cap - only3_admur_sum$admin_cases_avert_per_cap



only3_admur_sum$add_admin_cost <- only3_admur_sum$admin_cost - only3_admur_sum$admin_cost
pyrrole2_admur_sum$add_admin_cost <- pyrrole2_admur_sum$admin_cost - only3_admur_sum$admin_cost
pyrrole2c_admur_sum$add_admin_cost <- pyrrole2c_admur_sum$admin_cost - only3_admur_sum$admin_cost
pyrrole2m_admur_sum$add_admin_cost <- pyrrole2m_admur_sum$admin_cost - only3_admur_sum$admin_cost
pyrrole3_admur_sum$add_admin_cost <- pyrrole3_admur_sum$admin_cost - only3_admur_sum$admin_cost
pyrrole3c_admur_sum$add_admin_cost <- pyrrole3c_admur_sum$admin_cost - only3_admur_sum$admin_cost
pyrrole3m_admur_sum$add_admin_cost <- pyrrole3m_admur_sum$admin_cost - only3_admur_sum$admin_cost

# get shapefiles
adm1.shp <- raster::getData("GADM", country = "SEN", level = 1)
adm1.shp.f <- sf::st_as_sf(adm1.shp, region = "NAME_1") #fortify
#adm1.shp.f <- tidy(adm1.shp, region = "NAME_1") #fortify

adm1.shp.f %<>%
  dplyr::rename(fs_name_1 = NAME_1)

cases_averted_lgd <- "Cases averted\nper 1,000"

only3_shapes <- merge(adm1.shp.f, as_tibble(only3_admur_sum), by = "fs_name_1")
only3_carto_shapes <- st_transform(only3_shapes, "ESRI:102022") %>%
  cartogram_cont(weight = "admin_pop")
# pyrrole3_shapes <- merge(adm1.shp.f, as_tibble(pyrrole3_admur_sum), by = "fs_name_1")
# pyrrole3_carto_shapes <- st_transform(pyrrole3_shapes, "ESRI:102022") %>%
#   cartogram_cont(weight = "admin_pop")
pyrrole2_shapes <- merge(adm1.shp.f, as_tibble(pyrrole2_admur_sum), by = "fs_name_1")
pyrrole2_carto_shapes <- st_transform(pyrrole2_shapes, "ESRI:102022") %>%
  cartogram_cont(weight = "admin_pop")
pyrrole2c_shapes <- merge(adm1.shp.f, as_tibble(pyrrole2c_admur_sum), by = "fs_name_1")
pyrrole2c_carto_shapes <- st_transform(pyrrole2c_shapes, "ESRI:102022") %>%
  cartogram_cont(weight = "admin_pop")
pyrroleD_shapes <- merge(adm1.shp.f, as_tibble(pyrroleD_admur_sum), by = "fs_name_1")
pyrroleD_carto_shapes <- st_transform(pyrroleD_shapes, "ESRI:102022") %>%
  cartogram_cont(weight = "admin_pop")
pyrrole2m_shapes <- merge(adm1.shp.f, as_tibble(pyrrole2m_admur_sum), by = "fs_name_1")
pyrrole2m_carto_shapes <- st_transform(pyrrole2m_shapes, "ESRI:102022") %>%
  cartogram_cont(weight = "admin_pop")

pyrrole3_shapes <- merge(adm1.shp.f, as_tibble(pyrrole3_admur_sum), by = "fs_name_1")
pyrrole3_carto_shapes <- st_transform(pyrrole3_shapes, "ESRI:102022") %>%
  cartogram_cont(weight = "admin_pop")
pyrrole3c_shapes <- merge(adm1.shp.f, as_tibble(pyrrole3c_admur_sum), by = "fs_name_1")
pyrrole3c_carto_shapes <- st_transform(pyrrole3c_shapes, "ESRI:102022") %>%
  cartogram_cont(weight = "admin_pop")
pyrrole3m_shapes <- merge(adm1.shp.f, as_tibble(pyrrole3m_admur_sum), by = "fs_name_1")
pyrrole3m_carto_shapes <- st_transform(pyrrole3m_shapes, "ESRI:102022") %>%
  cartogram_cont(weight = "admin_pop")



depri2_carto_shapes <- pyrrole2m_carto_shapes %>% filter(fs_name_1 %in% deprioritised2_adms)
depri3_carto_shapes <- pyrrole3m_carto_shapes %>% filter(fs_name_1 %in% deprioritised3_adms)


cases_avert_per_1000_breaks = seq(0,1.6,0.4)#c(100,1000,10000)

ggplot() +
  geom_sf(data = pyrrole3m_carto_shapes,
          aes(group = fs_name_1,
              fill = admin_cases_avert_per_cap),
          ) +
  scale_fill_viridis(option = "turbo",
                     begin = 0.4,
                     end = 1,
                     direction = -1,
                     breaks = cases_avert_per_1000_breaks,
                     labels = cases_avert_per_1000_breaks,
                     limits = c(0,1.6)
  ) +
  geom_sf(data = depri3_carto_shapes,
          fill = "transparent",
          color = "gray20",
          linewidth = 1.5
  )+
  # guides(fill = guide_colorbar(title = "Annualised cases\naverted per capita",
  #                              title.position="top",
  #                              barwidth = 8,
  #                              barheight = 0.5)) +
  guides(fill = "none") +
  theme_map()
  # theme(legend.position="bottom",
  #       legend.title.align=0.5)
ggsave(paste0("pyrrole3m_SN_casesavertpercap",".png"), bg = "white",
       w = 5, h = 5, dpi = 1000)

cost_breaks = seq(0,8,2)
ggplot() +
  geom_sf(data = pyrrole3m_carto_shapes,
          aes(group = fs_name_1,
              fill = admin_cost/1e6)) +
  scale_fill_viridis(option = "rocket",
                     #trans = "log",
                     direction = -1,
                     breaks = cost_breaks,
                     labels = cost_breaks,
                     limits = c(0,8)
  ) +
  geom_sf(data = depri3_carto_shapes,
          fill = "transparent",
          color = "gray20",
          linewidth = 1.5
  )+
  # scale_size_identity() +
  # guides(fill = guide_colorbar(title = "Annualised cost\n($M USD)",
  #                              title.position="top",
  #                              barwidth = 8,
  #                              barheight = 0.5)) +
  guides(fill = "none") +
  theme_map()
  # theme(legend.position="bottom",
  #       legend.title.align=0.5)
ggsave(paste0("pyrrole3m_SN_cost",".png"), bg = "white",
       w = 5, h = 5, dpi = 1000)



pyrrole2m_carto_shapes$admin_cases_avert_per_USD[pyrrole2m_carto_shapes$admin_cases_avert_per_USD > 1.2] <- 1.2
pyrrole2c_carto_shapes$admin_cases_avert_per_USD[pyrrole2c_carto_shapes$admin_cases_avert_per_USD > 1.2] <- 1.2

cases_averted_per_USD_breaks = seq(0,1.2,0.4)
cases_averted_per_USD_brlabs = c("0","0.4","0.8",">1.2")
ggplot() +
  geom_sf(data = pyrrole3_carto_shapes,
          aes(group = fs_name_1,
              fill = admin_cases_avert_per_USD)) +
  scale_fill_viridis(option = "viridis",
                     #trans = "log",
                     breaks = cases_averted_per_USD_breaks,
                     labels = cases_averted_per_USD_brlabs,
                     limits = c(0,1.2)
  ) +
  # geom_sf(data = depri_carto_shapes,
  #         fill = "transparent",
  #         color = "gray20",
  #         linewidth = 1.5
  # )+
  # scale_size_identity() +
  #guides(fill = guide_colorbar(title = "Cases averted\nper $USD")) +
  # guides(fill = guide_colorbar(title = "Annualised cases\naverted per $USD",
  #                              title.position="top",
  #                              barwidth = 8,
  #                              barheight = 0.5)) +
  guides(fill = "none") +
  theme_map()
  # theme(legend.position="bottom",
  #       legend.title.align=0.5)
ggsave(paste0("pyrrole3_SN_perusd",".png"), bg = "white",
       w = 5, h = 5, dpi = 1000)


ggplot() +
  geom_sf(data = pyrrole3m_carto_shapes,
          aes(group = fs_name_1,
              fill = add_admin_cases_avert_per_cap),
  ) +
  paletteer::scale_fill_paletteer_c("pals::ocean.curl",
                                    limits = c(-0.8,0.8),
                                    direction = -1) +
  geom_sf(data = depri3_carto_shapes,
          fill = "transparent",
          color = "gray20",
          linewidth = 1.5
  ) +
  # scale_size_identity() +
  # guides(fill = guide_colorbar(title = "Additional ann.cases\naverted per capita",
  #                              title.position="top",
  #                              barwidth = 8,
  #                              barheight = 0.5)) +
  guides(fill = "none") +
  theme_map()
  # theme(legend.position="bottom",
  #       legend.title.align=0.5)
ggsave(paste0("pyrrole3m_SN_diffcasesavertpercap",".png"), bg = "white",
       w = 5, h = 5, dpi = 1000)


ggplot() +
  geom_sf(data = pyrrole2m_carto_shapes,
          aes(group = fs_name_1,
              fill = add_admin_cost / 1e6),
  ) +
  paletteer::scale_fill_paletteer_c("pals::ocean.balance",
                                    limits = c(-4,4),
                                    direction = 1) +
  geom_sf(data = depri_carto_shapes,
          fill = "transparent",
          color = "gray20",
          linewidth = 1.5
  ) +
  scale_size_identity() +
  guides(fill = guide_colorbar(title = "Additional ann.\ncost ($M USD)",
                               title.position="top",
                               barwidth = 8,
                               barheight = 0.5)) +
  #guides(fill = "none") +
  theme_map() +
  theme(legend.position="bottom",
        legend.title.align=0.5)
ggsave(paste0("legSN_addcost",".png"), bg = "white",
       w = 5, h = 5, dpi = 1000)


paletteer::scale_fill_paletteer_c("pals::ocean.curl",
                                  limits = c(-120,120),
                                  direction = -1) +



# 
# append_adm_weighted_data <- function(net_df = NULL) {
#   net_df$urbanicity_weighting <- net_df$pop / net_df$admin_pop
#   net_df$adm_cases <- 
#   net_df
# }
# only0 %<>% append_adm_weighted_data()
# only3 %<>% append_adm_weighted_data()
# pyrrole3 %<>% append_adm_weighted_data()
# pyrrole2 %<>% append_adm_weighted_data()



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

pyrrole3_admur_sum <- pyrrole3 %>% summarise_admin_costs_cases# %>% as.data.frame()
pyrrole3_admur_sum <- pyrrole3_admur_sum[,2:dim(pyrrole3_admur_sum)[2]] %>% unique()
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

pyrrole2_shapes <- merge(adm1.shp.f, as_tibble(pyrrole3_admur_sum), by = "fs_name_1")
pyrrole2_carto_shapes <- st_transform(pyrrole2_shapes, "ESRI:102022") %>%
  cartogram_cont(weight = "admin_pop")


pyrrole2_cartog <- cartogram_ncont(st_transform(pyrrole2_shapes), weight = "pop")

#test <- merge(st_transform(adm1.shp), pyrrole2_sum, by = "fs_name_1")
my_breaks = c(0.1,1,10,100,1000)#c(2,5,10,20,50,100,200,500,1000)
ggplot() +
  geom_sf(data = only3_carto_shapes,
          aes(group = fs_name_1,
              fill = admin_cases_per_cap * 1e3)) +
  # geom_sf(data = pyrrole2_carto_shapes,
  #         aes(group = fs_name_1,
  #             fill = cases_avert_mean)) +
  scale_fill_viridis(option = "mako",
                     trans = "log",
                     breaks = my_breaks,
                     labels = my_breaks,
                     limits = c(0.1,5000)
                     ) +
  guides(fill = guide_colorbar(title = "")) +
  #guides(fill = "none") +
  theme_map()
ggsave(paste0("only3SN_casesper1000",timestamp,".png"), bg = "white",
       w = 5, h = 5, dpi = 1000)


# scale_fill_viridis(begin = .5,
#                    end = 1,
#                    direction = -1,
#                    option = "turbo") +

# Generate adm-urbanicity sum
adm_urbanicity_sum pyrrole

# 
only3_country_cost <- sum(only3_sum$cost_mean)
pyrrole3_country_cost <- sum(pyrrole3$cost_mean)

priority_ids <- order(pyrrole2_sum$cases_avert_per_cap_mean,
                                     decreasing = TRUE)

pyrrole32mix_cost <- pyrrole3_country_cost
for 








total_test <- only3 %>%
  group_by(sample_index) %>%
  dplyr::summarise(total_cost = sum(avg_ann_net_cost))

mtcars %>%
  group_by(cyl) %>%
  dplyr::summarise(mean = mean(disp), n = n())