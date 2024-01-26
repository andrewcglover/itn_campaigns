# retention.R

append_usage_ret_cols <- function(col_names) {
  col_names <- c(col_names, "invlam_u_mean", "invlam_u_LB1", "invlam_u_UB1",
                "ret_u_mean", "ret_u_LB1", "ret_u_UB1")
  return(col_names)
}

append_access_ret_cols <- function(col_names) {
  ret_cols <- c(col_names, "invlam_u_mean", "invlam_u_LB1", "invlam_u_UB1",
                "ret_u_mean", "ret_u_LB1", "ret_u_UB1")
  return(col_names)
}

annual_ret_adjust <- function(ret_dataset) {
  Nretcol <- dim(ret_dataset)[2]
  ret_dataset[8:Nretcol] %<>% divide_by(12)
  return(ret_dataset)
}

fetch_retention_final <- function(dataset,
                                  usage = TRUE,
                                  access = TRUE,
                                  annual_adjustment = TRUE) {
  
  # Declare column names for selection
  ret_cols <- c("ISO2", "ADM1", "area", "urbanicity", "area_id")
  if (usage) {ret_cols %<>% append_usage_ret_cols}
  if (access) {ret_cols %<>% append_access_ret_cols}
  
  # Filter for final retention values
  retention_final <- dataset %>% filter(CMC == CMC_last)
  retention_final <- retention_final[ret_cols]
  
  # Annual adjustment
  if (annual_adjustment) {retention_final %<>% annual_ret_adjust}
  
  # Return final retention dataframe
  return(retention_final)
}

fetch_retention_period <- function(dataset, usage = TRUE, access = TRUE,
                                   CMCa = CMC_first, CMCb = CMC_last) {
  
  # Declare column names for selection
  ret_cols <- c("ISO2", "ADM1", "area", "urbanicity", "area_id")
  if (usage) {ret_cols %<>% append_usage_ret_cols}
  if (access) {ret_cols %<>% append_access_ret_cols}
  
  # Filter for final retention values
  retention_period <- dataset %>% filter(CMC >= CMCa & CMC <= CMCb)
  retention_period <- retention_period[ret_cols]
  
  # Annual adjustment
  if (annual_adjustment) {retention_period %<>% annual_ret_adjust}
  
  # Calculate mean retention over period
  if (usage) {
    if (access) {
      final_retention %<>%
        group_by(ISO2,
                 ADM1,
                 area,
                 urbanicity,
                 old_area_id,
                 area_id,
                 invlam_u_mean,
                 invlam_u_LB1,
                 invlam_u_UB1,
                 invlam_a_mean,
                 invlam_a_LB1,
                 invlam_a_UB1) %>%
        summarise_at(vars(ret_u_mean,
                          ret_u_LB1,
                          ret_u_UB1,
                          ret_a_mean,
                          ret_a_LB1,
                          ret_a_UB1),
                     list(avg = mean))
    } else {
      final_retention %<>%
        group_by(ISO2,
                 ADM1,
                 area,
                 urbanicity,
                 old_area_id,
                 area_id,
                 invlam_u_mean,
                 invlam_u_LB1,
                 invlam_u_UB1) %>%
        summarise_at(vars(ret_u_mean,
                          ret_u_LB1,
                          ret_u_UB1),
                     list(avg = mean))
    }
  } else if (access) {
    final_retention %<>%
      group_by(ISO2,
               ADM1,
               area,
               urbanicity,
               old_area_id,
               area_id,
               invlam_a_mean,
               invlam_a_LB1,
               invlam_a_UB1) %>%
      summarise_at(vars(ret_a_mean,
                        ret_a_LB1,
                        ret_a_UB1),
                   list(avg = mean))
  } else {
    print(paste0("Warning: Neither usage or access requested for retention ",
                 "reference period."))
  }
  
  # Return final retention dataframe
  return(final_retention)
}