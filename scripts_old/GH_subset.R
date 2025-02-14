#GH_subset.R

GH_used_nets <- net_data %>%
  filter(ISO2 == "GH") %>%
  select(ISO2, ADM1, urbanicity, old_area_id, area_id, CMC, used, total,
         Pbb_u_mean, Pbb_u_LB1, Pbb_u_UB1,
         Pbb_u_LB2, Pbb_u_UB2, Pbb_u_LB3, Pbb_u_UB3,
         P_u_mean, P_u_LB1, P_u_UB1)

GH_used_nets_former_ADM1 <- GH_used_nets

replace_regions <- function(dataset,
                            recorded_region,
                            updated_region,
                            usage = TRUE,
                            access = FALSE,
                            total = TRUE) {
  for (i in 1:2) {
    if (i == 1){
      recorded_ids <- which(dataset$ADM1 == recorded_region &
                              dataset$urbanicity == "urban")
      updated_ids <- which(dataset$ADM1 == updated_region &
                             dataset$urbanicity == "urban")
    } else {
      recorded_ids <- which(dataset$ADM1 == recorded_region &
                              dataset$urbanicity == "rural")
      updated_ids <- which(dataset$ADM1 == updated_region &
                             dataset$urbanicity == "rural")
    }
    if (usage) {
      dataset$used[updated_ids] <- dataset$used[updated_ids] + dataset$used[recorded_ids]
    }
    if (access) {
      dataset$access[updated_ids] <- dataset$access[updated_ids] + dataset$access[recorded_ids]
    }
    if (total) {
      dataset$total[updated_ids] <- dataset$total[updated_ids] + dataset$total[recorded_ids]
    }
  }
  dataset %<>% filter(ADM1 != recorded_region)
  return(dataset)
}

GH_used_nets_former_ADM1 %<>% replace_regions("Bono", "Brong-Ahafo")
GH_used_nets_former_ADM1 %<>% replace_regions("Bono East", "Brong-Ahafo")
GH_used_nets_former_ADM1 %<>% replace_regions("Ahafo", "Brong-Ahafo")
GH_used_nets_former_ADM1 %<>% replace_regions("Savannah", "Northern")
GH_used_nets_former_ADM1 %<>% replace_regions("North East", "Northern")
GH_used_nets_former_ADM1 %<>% replace_regions("Oti", "Volta")
GH_used_nets_former_ADM1 %<>% replace_regions("Western North", "Western")

GH_used_nets_former_ADM1$year = CMC_to_date(GH_used_nets_former_ADM1$CMC)[1]
GH_used_nets_former_ADM1$month = CMC_to_date(GH_used_nets_former_ADM1$CMC)[2]