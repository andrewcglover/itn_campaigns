# retention histogram

uret_hist <- region_summary %>%
  dplyr::filter(variable == "ret_u") %>%
  ggplot(aes(x = mean_val, y = ..density.., weight = pop, fill = ISO2)) +
  geom_histogram(alpha = 0.5) +
  theme_bw()
uret_hist

uret_hist <- region_summary %>%
  dplyr::filter(variable == "ret_u") %>%
  dplyr::filter(!is.nan(pop)) %>%
  ggplot(aes(x = mean_val, weight = pop, fill = ISO2, colour = ISO2)) +
  geom_density(alpha = 0.1) +
  theme_bw()
uret_hist

ret_hist <- region_summary %>%
  dplyr::filter(variable %in% c("ret_u", "ret_a")) %>%
  dplyr::mutate(country = countrycode(ISO2, "iso2c", "country.name")) %>%
  dplyr::mutate(
    var_name = ifelse(variable == "ret_a","Retention time", "Duration of use")
  ) %>%
  ggplot(aes(x = mean_val, y = ..density.., weight = pop, fill = country)) +
  geom_histogram(alpha = 0.5, breaks = seq(12, 42, 1)) +
  scale_x_continuous(breaks = seq(12, 42, 6), limits = c(12, 42)) +
  theme_bw() +
  facet_grid(rows = vars(var_name), cols = vars(country))
ret_hist


ret_hist <- region_summary %>%
  dplyr::filter(variable %in% c("ret_u", "ret_a")) %>%
  dplyr::mutate(country = countrycode(ISO2, "iso2c", "country.name")) %>%
  dplyr::mutate(
    var_name = ifelse(variable == "ret_a","Retention time", "Duration of use")
  ) %>%
  ggplot(aes(y = mean_val, x = ..density.., weight = pop, fill = country)) +
  geom_histogram(alpha = 0.5, breaks = seq(12, 42, 1)) +
  scale_y_continuous(breaks = seq(12, 42, 6), limits = c(12, 42)) +
  theme_bw() +
  facet_grid(rows = vars(var_name), cols = vars(country))
ret_hist

library(ggridges)

# Scale pop down & uncount
ret_expanded <- region_summary %>%
  filter(variable %in% c("ret_u", "ret_a")) %>%
  mutate(
    pop_scaled = round(pop / 10000)   # scale weights down by factor of 10,000
  ) %>%
  filter(pop_scaled > 0) %>%          # remove rows that would be dropped (pop=0)
  uncount(weights = pop_scaled)       # expands the data

# Plot ridge plot
ret_ridge <- ret_expanded %>%
  mutate(
    country = countrycode(ISO2, "iso2c", "country.name"),
    var_name = ifelse(variable == "ret_a", "Retention time", "Duration of use")
  ) %>%
  ggplot(aes(x = mean_val, y = country, fill = country)) +
  stat_density_ridges(alpha = 0.5, scale = 1, rel_min_height = 0.01, colour = NA) +
  scale_x_continuous(breaks = seq(12, 42, 6), limits = c(12, 42)) +
  theme_bw()

# --- prepare summary data with country names to match y-axis ---
summary_lines <- country_weighted_summary %>%
  filter(variable == "ret_u") %>%
  mutate(country = countrycode(ISO2, "iso2c", "country.name"))

# --- add reference lines and CI shading ---
ret_ridge <- ret_ridge +
  geom_vline(
    data = summary_lines,
    aes(xintercept = mean_weighted, colour = country),
    linetype = "dashed"
  )
ret_ridge







library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(countrycode)

# 1️⃣ Scale pop down & uncount to avoid exploding the dataset
ret_expanded <- region_summary %>%
  filter(variable %in% c("ret_u", "ret_a")) %>%
  mutate(
    pop_scaled = round(pop / 10000)    # scale weights down by factor of 10,000
  ) %>%
  filter(pop_scaled > 0) %>%           # remove rows where scaled weight is 0
  uncount(weights = pop_scaled)        # replicate rows by scaled weight

# 2️⃣ Prepare ridge plot data
country_levels <- sort(unique(countrycode(ret_expanded$ISO2, "iso2c", "country.name")))

ret_expanded <- ret_expanded %>%
  mutate(
    country = countrycode(ISO2, "iso2c", "country.name"),
    var_name = ifelse(variable == "ret_a", "Retention time", "Duration of use")
  )

# 3️⃣ Prepare summary data for reference lines
summary_lines <- country_weighted_summary %>%
  filter(variable %in% c("ret_u", "ret_a")) %>%
  mutate(
    country = countrycode(ISO2, "iso2c", "country.name"),
    y_pos = match(country, country_levels),   # numeric y position for each country
    var_name = ifelse(variable == "ret_a", "Retention time", "Duration of use")
  ) %>%
dplyr::mutate(y_pos = (max(y_pos) + 1) - y_pos)

# 4️⃣ Plot with stat_density_ridges and country-coloured reference lines
ret_ridge <- ggplot(ret_expanded, aes(x = mean_val, y = country, fill = country)) +
  # --- Ridges ---
  stat_density_ridges(alpha = 0.5, scale = 1, rel_min_height = 0.01, colour = NA) +
  
  # --- Short dashed reference lines per ridge (colored by country) ---
  geom_segment(
    data = summary_lines,
    aes(
      x = mean_weighted, xend = mean_weighted,
      y = y_pos,   # start slightly below the ridge
      yend = y_pos + 1,
      colour = country   # 👈 map colour to country
    ),
    linetype = "dashed",
    size = 0.7,
    inherit.aes = FALSE
  ) +
  
  # --- Scales & theme ---
  scale_x_continuous(breaks = seq(12, 42, 6), limits = c(12, 42)) +
  scale_y_discrete(limits = country_levels) +
  theme_bw()

ret_ridge <- ret_ridge + facet_grid(cols = vars(var_name))

ret_ridge




library(ggplot2)
library(ggridges)


ret_hist <- ggplot(ret_expanded, aes(
  x = mean_val,
  y = country,
  fill = country
)) +
  
  # Global summary
  geom_segment(
    data = overall_summary %<>%
      dplyr::filter(variable %in% c("ret_u", "ret_a")) %>%
      dplyr::mutate(
        var_name = ifelse(variable == "ret_a", "Retention time", "Duration of use")
      ),
    aes(
      x = mean_val,
      y = 0.4,
      yend = 7
    ),
    inherit.aes = FALSE,
    colour = "black",
    linetype = "dashed",
    size = 0.4
  ) +
  
 
  
  # country-separated “histogram ridges”
  stat_binline(
    binwidth = 1,                     # like breaks seq(12,42,1)
    boundary = 11.5,
    alpha = 0.6,
    scale = 0.8,                        # keep each country's bins on its own level
    draw_baseline = FALSE,
    colour = NA
  ) +
  
  # dashed mean lines for each country (drawn short so they only cross that country’s row)
  geom_segment(
    data = summary_lines,
    aes(
      x = mean_weighted, xend = mean_weighted,
      y = y_pos,   # start a little below the strip
      yend = y_pos + 1,
      colour = country
    ),
    inherit.aes = FALSE,
    linetype = "dashed",
    size = 0.4
  ) +
  
  # X-axis & faceting
  scale_x_continuous(breaks = seq(12, 42, 6), limits = c(11.5, 42.5)) +
  scale_y_discrete(limits = rev(country_levels)) +
  facet_grid(cols = vars(var_name)) +

# ✅ Axis and legend labels
labs(
  x = "months",  # <-- change this text
  y = NULL                      # removes the y-axis title
) +
  
  # ✅ Theme tweaks
  theme_bw() +
  theme(
    legend.position = "none",     # remove the legend entirely
    strip.background = element_blank(),
    axis.title.y = element_blank() # hides the y-axis title just in case
  )

ret_hist

ggsave("country_ret.pdf", plot = ret_hist, bg = "transparent",
       w = 6, h = 3, dpi = 450)
