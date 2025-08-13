# Reformating code for removing the 7th row from som summary outputs

# Define all base paths
base_paths_x <- c(
  "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary/BF/Pyrethroid-only_2_year_interval",
  "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary/BF/Pyrethroid-only_3_year_interval",
  "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary/BF/Pyrethroid-PBO_2_year_interval",
  "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary/BF/Pyrethroid-PBO_3_year_interval",
  "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary/BF/Pyrethroid-Pyrrole_2_year_interval",
  "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary/BF/Pyrethroid-Pyrrole_3_year_interval"
)

# Total number of files
total_files_x <- length(base_paths_x) * 26 * 100
counter_x <- 0  # progress tracker

# Loop through each base path
for (base_path_x in base_paths_x) {
  for (fs_id_x in 1:26) {
    fs_folder_x <- file.path(base_path_x, paste0("fs_id_", fs_id_x))
    
    for (sim_id_x in 1:100) {
      sim_file_x <- file.path(fs_folder_x, paste0("sim", sim_id_x, ".rds"))
      
      # Read the .rds file
      sim_data_x <- readRDS(sim_file_x)
      
      # Remove the 7th row
      sim_data_x <- sim_data_x[-7, ]
      
      # Overwrite the file
      saveRDS(sim_data_x, sim_file_x)
      
      # Update and print progress
      counter_x <- counter_x + 1
      if (counter_x %% 100 == 0 || counter_x == total_files_x) {
        percent_x <- round((counter_x / total_files_x) * 100, 1)
        cat(paste0("Progress: ", percent_x, "% complete (", counter_x, "/", total_files_x, ")\n"))
      }
    }
  }
}

# Define file path to inspect
check_file_x <- "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary/BF/Pyrethroid-only_3_year_interval/fs_id_1/sim1.rds"

# Read the file
sim_check_x <- readRDS(check_file_x)

# View the first 13 columns
print(sim_check_x[, 1:13])

# Define paths to the new and old (backup) files
new_file_x <- "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary/BF/Pyrethroid-only_3_year_interval/fs_id_1/sim1.rds"
old_file_x <- "M:/Andrew/Github/itn_campaigns/data_results/urep_06_20250604112237/summary/BF/backup7lines/Pyrethroid-only_3_year_interval/fs_id_1/sim1.rds"

# Load both versions
new_sim_x <- readRDS(new_file_x)
old_sim_x <- readRDS(old_file_x)

# Compare the first 6 rows and first 13 columns
identical_subset_x <- identical(new_sim_x, old_sim_x[1:6,])

# Print result
if (identical_subset_x) {
  cat("Y: The first 6 rows are identical between the new and old files.\n")
} else {
  cat("N: The first 6 rows differ between the new and old files.\n")
}