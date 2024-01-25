# retention.R

fetch_retention_final <- function(dataset, usage = TRUE, access = TRUE) {
  
  # Declare column names for selection
  ret_cols <- c("ISO2", "ADM1", "area", "urbanicity", "area_id")
  if (usage) {
    ret_cols <- c(ret_cols, "invlam_u_mean", "invlam_u_LB1", "invlam_u_UB1",
                  "ret_u_mean", "ret_u_LB1", "ret_u_UB1")
  }
  if (access) {
    ret_cols <- c(ret_cols, "invlam_a_mean", "invlam_a_LB1", "invlam_a_UB1",
                  "ret_a_mean", "ret_a_LB1", "ret_a_UB1")
  }
  
  # Filter for final retention values
  final_retention <- dataset %>% filter(CMC == CMC_last)
  final_retention <- final_retention[ret_cols]
  
  # Return final retention dataframe
  return(final_retention)
}

fetch_retention_period <- function(dataset, usage = TRUE, access = TRUE,
                                   CMCa = CMC_first, CMC_b = CMC_last) {
  
  # Declare column names for selection
  ret_cols <- c("ISO2", "ADM1", "area", "urbanicity", "area_id")
  if (usage) {
    ret_cols <- c(ret_cols, "invlam_u_mean", "invlam_u_LB1", "invlam_u_UB1",
                  "ret_u_mean", "ret_u_LB1", "ret_u_UB1")
  }
  if (access) {
    ret_cols <- c(ret_cols, "invlam_a_mean", "invlam_a_LB1", "invlam_a_UB1",
                  "ret_a_mean", "ret_a_LB1", "ret_a_UB1")
  }
  
  # Filter for final retention values
  final_retention <- dataset %>% filter(CMC == CMC_last)
  final_retention <- final_retention[ret_cols]
  
  # Return final retention dataframe
  return(final_retention)
}