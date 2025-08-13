# adjustsite.R

#-------------------------------------------------------------------------------
# Functions from ForesiteExplorer.R

# Helper Functions
copy_and_prepare_rows <- function(row, num_copies, start_year) {
  new_rows <- replicate(num_copies, row, simplify = FALSE)
  for (i in seq_along(new_rows)) {
    new_rows[[i]]$year <- start_year + i
  }
  do.call(rbind, new_rows)
}

insert_new_rows <- function(site_data, new_rows, insert_after_idx) {
  site_data$interventions <- rbind(
    site_data$interventions[1:insert_after_idx, ],
    new_rows,
    site_data$interventions[(insert_after_idx + 1):nrow(site_data$interventions), ]
  )
  site_data
}

adjust_column_types <- function(site_data, updater_row_types) {
  for (j in seq_along(updater_row_types)) {
    if (typeof(site_data$interventions[[j]]) != updater_row_types[j]) {
      site_data$interventions[[j]] <- as(site_data$interventions[[j]], updater_row_types[j])
    }
  }
  site_data
}

adjust_row_indices <- function(ridx, num_added_rows, length_ri) {
  ridx + (num_added_rows) * 0:(length_ri - 1)
}

counterfactual_replacement <- function(interventions, counterfactual) {
  if (!counterfactual) {
    return(interventions)
  }
  
  interventions %>%
    dplyr::mutate(
      dn0 = ifelse(counterfactual, 0.541979954, dn0),
      rn0 = ifelse(counterfactual, 0.456350279, rn0),
      gamman = ifelse(counterfactual, 2.64, gamman)
    )
}

# Main Function
expand_interventions <- function(site_data, expand_year, delay, counterfactual) {
  max_year <- max(site_data$interventions$year)
  ridx <- which(site_data$interventions$year == max_year)
  updateridx <- adjust_row_indices(ridx, expand_year, length(ridx))
  
  for (i in seq_along(updateridx)) {
    updater_row <- site_data$interventions[updateridx[i], ]
    updater_row_types <- sapply(updater_row, class)
    new_rows <- copy_and_prepare_rows(updater_row, expand_year, max_year)
    colnames(new_rows) <- colnames(site_data$interventions)
    site_data <- insert_new_rows(site_data, new_rows, updateridx[i])
    site_data <- adjust_column_types(site_data, updater_row_types)
  }
  
  site_data$interventions <- site_data$interventions[1:(max(updateridx) + expand_year), ]
  
  new_max_year <- max(site_data$interventions$year) - expand_year - 1
  new_ridx <- which(site_data$interventions$year == new_max_year)
  total_inserted_rows <- 0
  
  for (idx in new_ridx) {
    adjusted_idx <- idx + total_inserted_rows
    delay_row <- site_data$interventions[adjusted_idx, ]
    new_rows <- copy_and_prepare_rows(delay_row, delay, delay_row$year)
    site_data <- insert_new_rows(site_data, new_rows, adjusted_idx)
    total_inserted_rows <- total_inserted_rows + delay
    
    update_start <- adjusted_idx + delay + 1
    update_end <- min(update_start + expand_year, nrow(site_data$interventions))
    if (update_start <= update_end) {
      site_data$interventions[update_start:update_end, "year"] <- 
        site_data$interventions[update_start:update_end, "year"] + delay
    }
  }
  
  site_data$interventions <- counterfactual_replacement(site_data$interventions, counterfactual)
  
  # Save interventions data
  # save_interventions_data(site_data)
  
  return(site_data)
  
}

#-------------------------------------------------------------------------------
# Function to rebase site files

trim_years_keep_actual <- function(adm_site, years_to_trim) {
  time_series_components <- c(
    "interventions", "cases_deaths", "prevalence",
    "population", "population_by_age", "demography",
    "vectors", "seasonality", "eir"
  )
  
  data_frames_to_process <- list()
  
  for (comp in time_series_components) {
    if (is.data.frame(adm_site[[comp]])) {
      data_frames_to_process[[comp]] <- adm_site[[comp]]
    } else if (is.list(adm_site[[comp]])) {
      for (subname in names(adm_site[[comp]])) {
        if (is.data.frame(adm_site[[comp]][[subname]]) &&
            "year" %in% names(adm_site[[comp]][[subname]])) {
          key <- paste(comp, subname, sep = "$")
          data_frames_to_process[[key]] <- adm_site[[comp]][[subname]]
        }
      }
    }
  }
  
  min_years <- sapply(
    data_frames_to_process,
    function(df) min(df$year, na.rm = TRUE)
  )
  global_min_year <- min(min_years)
  cutoff_year <- global_min_year + years_to_trim
  
  for (key in names(data_frames_to_process)) {
    df <- data_frames_to_process[[key]]
    if (!"year" %in% names(df)) next
    df <- df[df$year >= cutoff_year, , drop = FALSE]
    
    if (grepl("\\$", key)) {
      parts <- strsplit(key, "\\$", fixed = TRUE)[[1]]
      adm_site[[parts[1]]][[parts[2]]] <- df
    } else {
      adm_site[[key]] <- df
    }
  }
  
  return(adm_site)
}
