#cache.R
# functions for caching and removing objs

handle_objects <- function(names, cache = FALSE, remove = FALSE, cache_path = NULL) {
  for (name in names) {
    if (exists(name, envir = .GlobalEnv)) {
      if (cache && !is.null(cache_path)) {
        saveRDS(get(name, envir = .GlobalEnv), 
                file = file.path(cache_path, paste0(name, ".rds")))
      }
      if (remove) {
        rm(list = name, envir = .GlobalEnv)
      }
    } else {
      cat("Object does not exist and was skipped:", name, "\n")
    }
  }
}
