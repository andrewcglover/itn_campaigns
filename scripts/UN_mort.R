# UN_mort.R
# UN WPP mortality rate calculation

library(peeps)
library(dplyr)
library(magrittr)

# Parameters
min_year <- 2000
max_year <- 2030
sim_pop <- 10000
sim_iso3c <- "NGA"

# Read in
orig_mort_data <- read.csv("NGA_mort_2000_2100.csv")
orig_pop_data <- read.csv("NGA_pop_2000_2100.csv")

# Filter for years of interest
mort_data <- orig_mort_data |>
  dplyr::filter(Time <= min_year) |>
  dplyr::filter(Time <= max_year)
pop_data <- orig_pop_data |>
  dplyr::filter(Time <= min_year) |>
  dplyr::filter(Time <= max_year)

# Unique regions, years and ages
UN_regions <- unique(mort_data$Location)
UN_iso2c <- unique(mort_data$Iso2)
UN_iso3c <- unique(mort_data$Iso3)
N_regions <- length(UN_region)

UN_years <- unique(mort_data$Time)
N_years <- length(UN_years)

UN_age_lower <- unique(mort_data$AgeStart)
UN_age_upper <- UN_age_lower + 1
N_ages <- length(UN_age_lower)
UN_age_upper[N_ages] <- 200

N_total_UN <- N_regions * N_years * N_ages

# Generate data frame of UN data
UN_WPP_df <- data.frame(
  "region" = rep(NA, N_total_UN),
  "iso3c" = rep(NA, N_total_UN),
  "year" = rep(NA, N_total_UN),
  "age_lower" = rep(NA, N_total_UN),
  "age_upper" = rep(NA, N_total_UN),
  "qx" = rep(NA, N_total_UN),
  "pop" = rep(NA, N_total_UN),
  "p" = rep(NA, N_total_UN)
)

# Populate UN data frame
x <- 1
for (i in 1:N_regions) {
  for(j in 1:N_years) {
    for(k in 1:N_ages) {
      UN_WPP_df$region[x] <- UN_regions[i]
      UN_WPP_df$iso3c[x] <- UN_iso3c[i]
      UN_WPP_df$year[x] <- UN_years[j]
      UN_WPP_df$age_lower[x] <- UN_age_lower[k]
      UN_WPP_df$age_upper[x] <- UN_age_upper[k]
      if ((mort_data$Location[x] == UN_regions[i]) &
          (mort_data$Time[x] == UN_years[j]) &
          (mort_data$AgeStart[x] == UN_age_lower[k])) {
        UN_WPP_df$qx[x] <- mort_data$Value[x]
      } else {
        print(paste0("Warning: mortality data mismatch at id = ", x))
      }
      if ((pop_data$Location[x] == UN_regions[i]) &
          (pop_data$Time[x] == UN_years[j]) &
          (pop_data$AgeStart[x] == UN_age_lower[k])) {
        UN_WPP_df$pop[x] <- pop_data$Value[x]
      } else {
        print(paste0("Warning: population data mismatch at id = ", x))
      }
      x <- x + 1
    }
  }
}

# Calculate proportions in each age group
UN_WPP_df <- UN_WPP_df |>
  dplyr::group_by(region, year) |>
  dplyr::mutate(total_pop = sum(pop))
UN_WPP_df$p <- UN_WPP_df$pop / UN_WPP_df$total_pop

# Filter country of interest
country_df <- UN_WPP_df |>
  dplyr::filter(iso3c == sim_iso3c)

#-------------------------------------------------------------------------------
# Model 1 - generalisable code for calculating & simulating new mortality rates

# Objective function from peeps
optim_objective <- function (mortality_rates, target_age_distribution) 
{
  sum(
    (peeps::equilibrium_age_distribution(
      mortality_rates
      ) - target_age_distribution)^2
    )
}

# Calculate new mortality rates
mortality_rates <- NULL
for (j in 1:N_years) {
  j0 <- 1 + (j - 1) * N_ages
  j1 <- j * N_ages
  year_df <- country_df[j0:j1,]
  if (!(all(year_df$year == UN_years[j]))) {print("Warning: year mismatch")}
  # adapted from peeps::estimate_mortality_rates
  opt <- stats::optim(par = year_df$qx, fn = optim_objective, 
                      method = "L-BFGS-B", lower = year_df$qx/2, 
                      upper = 0.99, target_age_distribution = year_df$p)
  mortality_rates <- c(mortality_rates, opt$par)
  print(j)
}
country_df$mortality_rates <- mortality_rates
daily_mortality_rates <- 1 - (1 - mortality_rates)^(1/365)

# Create parameter list for malariasimulation
params <- malariasimulation::get_parameters(
  overrides = list(
    human_population = sim_pop
  )
) |>
  malariasimulation::set_demography(
    agegroups = UN_age_upper * 365,
    timesteps = (UN_years - UN_years[1]) * 365,
    deathrates = matrix(
      daily_mortality_rates, nrow = N_years, byrow = TRUE
    )
  )

# Set rendering age groups
params$age_group_rendering_min_ages = UN_age_lower * 365
params$age_group_rendering_max_ages = UN_age_upper * 365

# Run simulation
sim_output <- malariasimulation::run_simulation(timesteps = 365 * N_years,
                                                parameters = params)

# Select the age-outputs
age_output <- sim_output |> 
  dplyr::select(paste0("n_age_", params$age_group_rendering_min_ages,
                "_", params$age_group_rendering_max_ages))

# Estimate proportions in each age group by year
age_proportion <- NULL
for (j in 1:31) {
  j0 <- 1 + (j - 1) * 365
  j1 <- j * 365
  age_year <- age_output[j0:j1,]
  age_proportion <- c(
    age_proportion,
    colMeans(age_year, na.rm = TRUE) / sum(colMeans(age_year, na.rm = TRUE))
  )
}


#-------------------------------------------------------------------------------
# Model 2 - Simulating mortality rates from new site file for Nigeria

# From new NGA site file
daily_mortality_rates2 <- 1 - (1 - new_demography$mortality_rate)^(1/365)
daily_mortality_rates2 <- daily_mortality_rates2[
  1:length(daily_mortality_rates)
  ]
params2 <- malariasimulation::get_parameters(
  overrides = list(
    human_population = sim_pop
  )
) |>
  malariasimulation::set_demography(
    agegroups = UN_age_upper * 365,
    timesteps = (UN_years - UN_years[1]) * 365,
    deathrates = matrix(
      daily_mortality_rates2, nrow = N_years, byrow = TRUE
    )
  )

# Set rendering age groups
params2$age_group_rendering_min_ages = UN_age_lower * 365
params2$age_group_rendering_max_ages = UN_age_upper * 365

# Run simulation
sim_output2 <- malariasimulation::run_simulation(timesteps = 365 * N_years,
                                                 parameters = params2)

# Select the age-outputs
age_output2 <- sim_output2 |> 
  dplyr::select(paste0("n_age_", params2$age_group_rendering_min_ages,
                       "_", params2$age_group_rendering_max_ages))

# Estimate proportions in each age group by year
age_proportion2 <- NULL
for (j in 1:N_years) {
  j0 <- 1 + (j - 1) * 365
  j1 <- j * 365
  age_year <- age_output2[j0:j1,]
  age_proportion2 <- c(
    age_proportion2,
    colMeans(age_year, na.rm = TRUE) / sum(colMeans(age_year, na.rm = TRUE))
  )
}

#-------------------------------------------------------------------------------
# Plot comparison


plot_data <- data.frame(
  year = rep(rep(seq(2000,2030), each = N_ages), 3),
  age = rep(UN_age_lower, N_years * 3),
  Output = rep(c("UNWPP", "Model 1", "Model 2"), each = 31 * N_ages),
  p = c(country_filter$p, age_proportion, age_proportion2)
)

# Plot the comparison of target and modelled age_distribution
ggplot(plot_data, aes(x = age, y = p, colour = Output)) +
  geom_path(alpha = 0.6) +
  xlab("Age") +
  ylab("Propotion of population") +
  facet_wrap(vars(year), ncol = 8) +
  theme_bw()