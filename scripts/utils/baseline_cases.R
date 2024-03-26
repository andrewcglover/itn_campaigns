# baseline_cases.R

# Fetch baseline cases and costs
record_baselines <- function(no_net_baseline = TRUE,
                             only_baseline = TRUE,
                             pbo_baseline = TRUE,
                             pyrrole_baseline = TRUE) {
  # Baseline cases
  if (no_net_baseline) {no_net_cases <<- only0$pred_ann_infect}
  if (only_baseline) {onlyD_cases <<- pboD$pred_ann_infect}
  if (pbo_baseline) {pboD_cases <<- pboD$pred_ann_infect}
  if (pyrrole_baseline) {pyrroleD_cases <<- pyrroleD$pred_ann_infect}
  # Baseline cases
  if (no_net_baseline) {no_net_cost <<- only0$avg_ann_net_cost}
  if (only_baseline) {onlyD_cost <<- pboD$avg_ann_net_cost}
  if (pbo_baseline) {pboD_cost <<- pboD$avg_ann_net_cost}
  if (pyrrole_baseline) {pyrroleD_cost <<- pyrroleD$avg_ann_net_cost}
  # NULL return
  return(NULL)
}

# Append cases averted
append_cases_averted <- function(net_df = NULL,
                                 baseline_df = only0
                                 ) {
  net_df$cases_baseline <- baseline_df$pred_ann_infect
  net_df$cost_baseline <- baseline_df$avg_ann_net_cost
  net_df$cases_averted <-net_df$ cases_baseline - net_df$pred_ann_infect
  net_df$cases_averted_per_capita <- net_df$cases_averted /  net_df$pop
  net_df$add_cost <- net_df$avg_ann_net_cost - net_df$cost_baseline
  net_df$cases_averted_per_USD <- net_df$cases_averted / net_df$add_cost
  net_df$cases_averted_per_USD[is.na(net_df$cases_averted_per_USD)] <- 0
  # if (cases_total) {net_df$add_cases_averted <- cases_averted}
  # if (cases_per_capita) {net_df$cases_averted_per_capita <- cases_averted_per_capita}
  # if (cases_per_USD) {net_df$add_cases_averted_per_USD <- cases_averted_per_USD}
  return(net_df)
}

# # Append cases averted
# append_cases_averted <- function(net_df = NULL,
#                                  baseline_df = NULL,
#                                  baseline_cases = TRUE,
#                                  baseline_cost = TRUE,
#                                  cases_total = TRUE,
#                                  cases_per_capita = TRUE,
#                                  cases_per_USD = TRUE) {
#   cases_baseline <- baseline_df$pred_ann_infect
#   cost_baseline <- baseline_df$avg_ann_net_cost
#   cases_averted <- cases_baseline - net_df$pred_ann_infect
#   cases_averted_per_capita <- cases_averted /  net_df$pop
#   add_cost <- net_df$avg_ann_net_cost - cost_baseline
#   cases_averted_per_USD <- cases_averted / add_cost
#   cases_averted_per_USD[is.na(cases_averted_per_USD)] <- 0
#   if (cases_total) {net_df$add_cases_averted <- cases_averted}
#   if (cases_per_capita) {net_df$cases_averted_per_capita <- cases_averted_per_capita}
#   if (cases_per_USD) {net_df$add_cases_averted_per_USD <- cases_averted_per_USD}
#   return(net_df)
# }