library(ggplot2)

# Some useful functions for converting between dates and calendar month codes
date_to_CMC <- function(year = 2000, month = 1) {
  CMC <- ((year-1900) * 12) + month
  return(CMC)
}
CMC_to_date <- function(CMC){
  yr <- floor((CMC - 1) / 12) + 1900
  mn <- CMC - ((yr-1900) * 12)
  return(data.frame("year" = yr, "month" = mn))
}

# Read in samples from usage model (all samples across all regions)
usage_draws_df <- readRDS("usage_draws_df.csv")

# Define when you want to simulate mass campaigns
# N.B. the time window is between Jan 2005 and Dec 2024 inclusive
campaign_yrs <- c(2006, 2009, 2012, 2015, 2018, 2021)
campaign_mns <- c(1,    3,    5,    2,    6,    8   )

# Append years and dates using calendar month codes
cmc_dates <- CMC_to_date(usage_draws_df$CMC)
usage_draws_df$year <- cmc_dates$year
usage_draws_df$month <- cmc_dates$month
usage_draws_df$decimal_year <- usage_draws_df$year + (usage_draws_df$month - 1) / 12

# Append mean duration of use
usage_draws_df$ret_u <- with(usage_draws_df, {P0 * invlam / C0})

# Calculate campaign CMCs
campaign_cmc <- mapply(date_to_CMC, year = campaign_yrs, month = campaign_mns)

# Calculate the months since the last campaign
index <- findInterval(usage_draws_df$CMC, campaign_cmc)
no_prior_match <- index == 0
index[no_prior_match] <- NA
last_cmc <- campaign_cmc[index]
months_since_last <- usage_draws_df$CMC - last_cmc
months_since_last[no_prior_match] <- 240
usage_draws_df$m <- months_since_last

# Calculate prob of using a campaign (C) or any net (P) given campaign timings
usage_draws_df$C <- with(usage_draws_df, {C0 * exp(-m/invlam)})
usage_draws_df$P <- with(usage_draws_df, {C + D})

# Generate a sample of parameter draws
n_sample_ids <- 1000
unique_ids <- unique(usage_draws_df$sample_id)
set.seed(123)
sampled_ids <- sample(unique_ids, n_sample_ids)
usage_draws_sample <- usage_draws_df[usage_draws_df$sample_id %in% sampled_ids,]

# Plot overall ITN use for the given sample
ggplot(usage_draws_sample, aes(x = decimal_year, y = P, group = sample_id)) +
  geom_line(alpha = 0.1) +
  labs(x = "Year", y = "ITN use") +
  theme_minimal()