# hipercow.R
# Helper functions for Hipercow

extract_hipercow_net_runs <- function(ids) {
  sim_data <- NULL
  for (i in 1:length(ids)) {
    scenario_data <- ids[i] %>% task_result %>% do.call(rbind.data.frame, .)
    sim_data %<>% rbind.data.frame(scenario_data)
  }
  return(sim_data)
}
