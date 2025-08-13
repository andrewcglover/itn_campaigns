library(tidyverse)

# Base path to the root of the country folders
base_pathx <- "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary"

# ISO2 country codes
countriesx <- c("BF", "GH", "ML", "MZ", "MW", "SN")

# Target net strategy folders
net_strategiesx <- c(
  "no_future_nets",
  "Pyrethroid-only_2_year_interval",
  "Pyrethroid-only_3_year_interval",
  "Pyrethroid-PBO_2_year_interval",
  "Pyrethroid-PBO_3_year_interval",
  "Pyrethroid-Pyrrole_2_year_interval",
  "Pyrethroid-Pyrrole_3_year_interval",
  "Pyrethroid-only_routine_baseline",
  "Pyrethroid-PBO_routine_baseline",
  "Pyrethroid-Pyrrole_routine_baseline"
)

# Function to count sim files in a given fs_id folder
count_sim_filesx <- function(pathx) {
  filesx <- list.files(pathx, pattern = "^sim\\d+\\.rds$", full.names = TRUE)
  length(filesx)
}

# Build the data frame
resultsx <- purrr::map_dfr(countriesx, function(iso2x) {
  purrr::map_dfr(net_strategiesx, function(strategyx) {
    strat_pathx <- file.path(base_pathx, iso2x, strategyx)
    if (!dir.exists(strat_pathx)) return(NULL)
    
    fs_dirsx <- list.dirs(strat_pathx, full.names = TRUE, recursive = FALSE)
    fs_dirsx <- fs_dirsx[grepl("fs_id_\\d+$", basename(fs_dirsx))]
    
    purrr::map_dfr(fs_dirsx, function(fs_dirx) {
      tibble(
        ISO2 = iso2x,
        net_strategy = strategyx,
        fs_area_id = basename(fs_dirx),
        tot_sims = count_sim_filesx(fs_dirx)
      )
    })
  })
})

incompletex <- resultsx %>% filter(tot_sims != 100)



