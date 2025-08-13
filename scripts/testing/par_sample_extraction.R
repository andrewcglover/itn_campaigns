# Load packages
library(ggplot2)
library(bayestestR)  # for median, CCI, and HDI

# usage_draws <- usage_fit_raw$draws(format = "draws_df")
# Pbb_u <<- usage_draws %>% select(starts_with("u_tilde[")) %>% divide_by(N_bb)
# P_u <<- usage_draws %>% select(starts_with("P["))
# P0_u <<- usage_draws %>% select(starts_with("P0["))
# D_u <<- usage_draws %>% select(starts_with("D["))
# C_u <<- usage_draws %>% select(starts_with("C["))
# C0_u <<- usage_draws %>% select(starts_with("C0["))
# invlam_u <<- usage_draws %>% select(starts_with("inv_lambda["))
# overdisp_u <<- usage_draws %>% select(starts_with("alpha0["))

t_vec <- usage_list$t
t_mean <- mean(t_vec)
t_sd <- sd(t_vec)
t_hat_vec <- (t_vec - t_mean) / t_sd

N_draws <- dim(P0_u)[1]
N_area_times <- dim(P0_u)[2]

P0_u_vec <- as.vector(t(as.matrix(P0_u)))
C0_u_vec <- as.vector(t(as.matrix(C0_u)))
D_u_vec <- as.vector(t(as.matrix(D_u)))
invlam_u_vec <- rep(as.vector(t(as.matrix(invlam_u))), each = usage_list$N_t)

usage_draws_df <- data.frame(
  "sample_id" = rep(seq(1:(N_draws*usage_list$N_a)), each = usage_list$N_t),
  "CMC" = rep(t_vec, N_draws),
  "P0" = P0_u_vec,
  "C0" = C0_u_vec,
  "D" = D_u_vec,
  "invlam" = invlam_u_vec)
saveRDS(usage_draws_df, "usage_draws_df.rds")

usage_draws_df <- readRDS("usage_draws_df.csv")

campaign_yrs <- c(2006, 2009, 2012, 2015, 2018, 2021)
campaign_mns <- c(1,    3,    5,    2,    6,    8   )

campaign_cmc <- mapply(date_to_CMC, year = campaign_yrs, month = campaign_mns)

index <- findInterval(usage_draws_df$CMC, campaign_cmc)
no_prior_match <- index == 0
index[no_prior_match] <- NA
last_cmc <- campaign_cmc[index]
months_since_last <- usage_draws_df$CMC - last_cmc
months_since_last[no_prior_match] <- 240
usage_draws_df$m <- months_since_last

usage_draws_df$C <- with(usage_draws_df, {C0 * exp(-m/invlam)})
usage_draws_df$P <- with(usage_draws_df, {C + D})

n_sample_ids <- 100
unique_ids <- unique(usage_draws_df$sample_id)
set.seed(123)
sampled_ids <- sample(unique_ids, n_sample_ids)
usage_draws_sample <- usage_draws_df[usage_draws_df$sample_id %in% sampled_ids,]

ggplot(usage_draws_sample, aes(x = CMC, y = P, group = sample_id)) +
  geom_line(alpha = 0.2) +
  labs(x = "Calendar Month Code", y = "ITN use") +
  theme_minimal()



beta0 <- usage_draws %>% dplyr::select(starts_with("beta_0["))
betat <- usage_draws %>% dplyr::select(starts_with("beta_t["))


use_draws_df <- data.frame("CMC" = usage_list$t)


ref_ret_u <- ret_u[, which(usage_list$t == ref_CMC)]

ref_ret_u <- ret_u[, which(usage_list$t == ref_CMC)]

ref_ret_u_vec <- as.vector(as.matrix(ref_ret_u))

# Convert to data frame
ref_ret_u_df <- data.frame(ret_u = ref_ret_u_vec)

# Compute summary stats
med_ret_u <- median(ref_ret_u_vec)
cci_ret_u <- ci(ref_ret_u_vec, method = "ETI", ci = 0.95)
hdi_ret_u <- hdi(ref_ret_u_vec, ci = 0.95)

# Plot
ggplot(ref_ret_u_df, aes(x = ret_u)) +
  geom_density(fill = "skyblue", alpha = 0.4) +
  
  # Median line
  geom_vline(xintercept = med_ret_u, color = "black", linetype = "dashed", size = 1) +
  
  # Central credible interval (ETI)
  geom_segment(aes(x = cci_ret_u$CI_low, xend = cci_ret_u$CI_high, y = 0, yend = 0),
               color = "blue", size = 2) +
  
  # Highest density interval (HDI)
  geom_segment(aes(x = hdi_ret_u$CI_low, xend = hdi_ret_u$CI_high, y = -0.005, yend = -0.005),
               color = "red", size = 2) +
  
  labs(title = "Density of ref_ret_u_vec with Median, CCI, and HDI",
       x = "ret_u", y = "Density") +
  theme_minimal()




# Load ggplot2
library(ggplot2)

ref_ret_u <- ret_u[, which(usage_list$t == ref_CMC)]

ref_ret_u_vec <- as.vector(as.matrix(ref_ret_u))

# Convert to data frame
ref_ret_u_df <- data.frame(ret_u = ref_ret_u_vec)

# Compute median and 95% central credible interval (CCI)
med_ret_u <- median(ref_ret_u_vec)
cci_ret_u <- quantile(ref_ret_u_vec, probs = c(0.025, 0.975))

# Summary row for plotting
summary_df <- data.frame(
  y = -0.005,
  x = med_ret_u,
  xmin = cci_ret_u[1],
  xmax = cci_ret_u[2]
)

# Round values for annotation
label_x <- round(summary_df$x, 1)
label_xmin <- round(summary_df$xmin, 1)
label_xmax <- round(summary_df$xmax, 1)
label_y <- summary_df$y - 0.002  # adjust vertical position slightly

# Plot
ggplot() +
  geom_density(data = ref_ret_u_df, aes(x = ret_u), fill = "skyblue", colour = NA, alpha = 0.6) +
  
  # Error bar (horizontal line) — all aesthetics defined here
  geom_errorbarh(data = summary_df, aes(xmin = xmin, xmax = xmax, y = y), height = 0, color = "grey20", size = 1.5) +
  
  # Median point
  geom_point(data = summary_df, aes(x = x, y = y), size = 3, color = "grey20") +
  
  # Annotations
  annotate("text", x = label_x, y = label_y, label = paste("Median:", label_x), vjust = 0.8, size = 3.5) +
  annotate("text", x = label_xmin, y = label_y, label = paste("Lower:", label_xmin), vjust = 0, hjust = 1.1, size = 3.5) +
  annotate("text", x = label_xmax, y = label_y, label = paste("Upper:", label_xmax), vjust = 0, hjust = -0.1, size = 3.5) +
  
  
  labs(#title = "Density with Median (dot) and 95% CCI (line)",
       x = "Mean duration of use (months)", y = "Density") +
  theme_minimal()

ggsave("ref_ret_u_density_plot.png", width = 9, height = 6, dpi = 300,
       units = "in", bg = "white")

saveRDS(ref_ret_u_vec, file = "ref_ret_u_vec.rds")

# Load required libraries
library(ggplot2)
library(sn)

# 1. Simulate or input your vector
ref_ret_u <- ret_u[, which(usage_list$t == ref_CMC)]
ref_ret_u_vec <- as.vector(as.matrix(ref_ret_u))

# 2. Convert to data frame
ref_ret_u_df <- data.frame(ret_u = ref_ret_u_vec)

# 3. Summary statistics
med_ret_u <- median(ref_ret_u_vec)
cci_ret_u <- quantile(ref_ret_u_vec, probs = c(0.025, 0.975))

# 4. Prepare for point + line summary
summary_df <- data.frame(
  y = -0.005,
  x = med_ret_u,
  xmin = cci_ret_u[1],
  xmax = cci_ret_u[2]
)

# 5. Fit skew-normal using known quantiles
objective <- function(params) {
  xi <- params[1]
  omega <- params[2]
  alpha <- params[3]
  q <- qsn(c(0.025, 0.5, 0.975), xi, omega, alpha)
  sum((q - c(cci_ret_u[1], med_ret_u, cci_ret_u[2]))^2)
}

fit <- optim(par = c(med_ret_u, sd(ref_ret_u_vec), 0),
             fn = objective,
             method = "L-BFGS-B",
             lower = c(-Inf, 1e-6, -20), upper = c(Inf, Inf, 20))

fitted_params <- fit$par
names(fitted_params) <- c("xi", "omega", "alpha")

# 6. Generate skew-normal density curve
x_vals <- seq(min(ref_ret_u_vec), max(ref_ret_u_vec), length.out = 500)
skew_dens <- dsn(x_vals, xi = fitted_params["xi"], omega = fitted_params["omega"], alpha = fitted_params["alpha"])
skew_df <- data.frame(x = x_vals, y = skew_dens)

# 7. Plot
ggplot() +
  # Original density
  geom_density(data = ref_ret_u_df, aes(x = ret_u), fill = "skyblue", colour = NA, alpha = 0.4) +
  
  # Skew-normal fitted curve
  geom_line(data = skew_df, aes(x = x, y = y), color = "red", size = 1) +
  
  # CCI line
  geom_errorbarh(data = summary_df, aes(xmin = xmin, xmax = xmax, y = y), height = 0, color = "grey20", size = 1.5) +
  
  # Median dot
  geom_point(data = summary_df, aes(x = x, y = y), size = 3, color = "grey20") +
  
  # Annotations above the line
  annotate("text", x = round(summary_df$x, 1), y = summary_df$y + 0.003,
           label = paste("Median:", round(summary_df$x, 1)), vjust = 0, size = 3.5) +
  annotate("text", x = round(summary_df$xmin, 1), y = summary_df$y + 0.003,
           label = paste("Lower:", round(summary_df$xmin, 1)), vjust = 0, hjust = 1.1, size = 3.5) +
  annotate("text", x = round(summary_df$xmax, 1), y = summary_df$y + 0.003,
           label = paste("Upper:", round(summary_df$xmax, 1)), vjust = 0, hjust = -0.1, size = 3.5) +
  
  labs(title = "Density with Median, 95% CCI, and Fitted Skew-Normal",
       x = "ret_u", y = "Density") +
  theme_minimal()




# Load required libraries
library(ggplot2)
library(MASS)
library(sn)

# 1. Extract data
ref_ret_u <- ret_u[, which(usage_list$t == ref_CMC)]
ref_ret_u_vec <- as.vector(as.matrix(ref_ret_u))
ref_ret_u_df <- data.frame(ret_u = ref_ret_u_vec)

# 2. Set up colors and names
dist_names <- c("Normal", "Skew-Normal", "Gamma", "Log-Normal", "Weibull")
dist_colors <- palette()[1:5]
names(dist_colors) <- dist_names

# 3. X range for densities
x_vals <- seq(min(ref_ret_u_vec), max(ref_ret_u_vec), length.out = 500)

# 4. Create containers
densities <- list()
quantiles <- list()
logliks <- numeric(5)
params <- list()

# 5. Normal
mu <- mean(ref_ret_u_vec)
sd_val <- sd(ref_ret_u_vec)
densities[["Normal"]] <- dnorm(x_vals, mean = mu, sd = sd_val)
quantiles[["Normal"]] <- qnorm(c(0.025, 0.5, 0.975), mean = mu, sd = sd_val)
logliks[1] <- sum(dnorm(ref_ret_u_vec, mean = mu, sd = sd_val, log = TRUE))
params[["Normal"]] <- sprintf("mu=%.2f, sd=%.2f", mu, sd_val)

# 6. Skew-Normal
skew_fit <- selm(ref_ret_u_vec ~ 1, family = "SN")
skew_params <- coef(skew_fit, "DP")
xi <- skew_params["xi"]
omega <- skew_params["omega"]
alpha <- skew_params["alpha"]
densities[["Skew-Normal"]] <- dsn(x_vals, xi = xi, omega = omega, alpha = alpha)
quantiles[["Skew-Normal"]] <- qsn(c(0.025, 0.5, 0.975), xi = xi, omega = omega, alpha = alpha)
logliks[2] <- sum(dsn(ref_ret_u_vec, xi = xi, omega = omega, alpha = alpha, log = TRUE))
params[["Skew-Normal"]] <- sprintf("xi=%.2f, omega=%.2f, alpha=%.2f", xi, omega, alpha)

# 7. Gamma
gamma_fit <- fitdistr(ref_ret_u_vec, densfun = "gamma")
shape_g <- gamma_fit$estimate["shape"]
rate_g <- gamma_fit$estimate["rate"]
densities[["Gamma"]] <- dgamma(x_vals, shape = shape_g, rate = rate_g)
quantiles[["Gamma"]] <- qgamma(c(0.025, 0.5, 0.975), shape = shape_g, rate = rate_g)
logliks[3] <- sum(dgamma(ref_ret_u_vec, shape = shape_g, rate = rate_g, log = TRUE))
params[["Gamma"]] <- sprintf("shape=%.2f, rate=%.2f", shape_g, rate_g)

# 8. Log-Normal
lnorm_fit <- fitdistr(ref_ret_u_vec, "log-normal")
meanlog <- lnorm_fit$estimate["meanlog"]
sdlog <- lnorm_fit$estimate["sdlog"]
densities[["Log-Normal"]] <- dlnorm(x_vals, meanlog = meanlog, sdlog = sdlog)
quantiles[["Log-Normal"]] <- qlnorm(c(0.025, 0.5, 0.975), meanlog = meanlog, sdlog = sdlog)
logliks[4] <- sum(dlnorm(ref_ret_u_vec, meanlog = meanlog, sdlog = sdlog, log = TRUE))
params[["Log-Normal"]] <- sprintf("meanlog=%.2f, sdlog=%.2f", meanlog, sdlog)

# 9. Weibull
weib_fit <- fitdistr(ref_ret_u_vec, "weibull")
shape_w <- weib_fit$estimate["shape"]
scale_w <- weib_fit$estimate["scale"]
densities[["Weibull"]] <- dweibull(x_vals, shape = shape_w, scale = scale_w)
quantiles[["Weibull"]] <- qweibull(c(0.025, 0.5, 0.975), shape = shape_w, scale = scale_w)
logliks[5] <- sum(dweibull(ref_ret_u_vec, shape = shape_w, scale = scale_w, log = TRUE))
params[["Weibull"]] <- sprintf("shape=%.2f, scale=%.2f", shape_w, scale_w)

# 10. Density data frame
dens_df <- do.call(rbind, lapply(seq_along(dist_names), function(i) {
  data.frame(
    x = x_vals,
    y = densities[[dist_names[i]]],
    distribution = factor(dist_names[i], levels = dist_names)
  )
}))

# 11. Summary data frame for intervals
summary_df <- do.call(rbind, lapply(seq_along(dist_names), function(i) {
  q <- quantiles[[dist_names[i]]]
  data.frame(
    x = q[2], xmin = q[1], xmax = q[3],
    y = -0.01 - 0.005 * i,
    distribution = factor(dist_names[i], levels = dist_names)
  )
}))

# 12. Annotations
annotation_df <- do.call(rbind, lapply(seq_along(dist_names), function(i) {
  q <- quantiles[[dist_names[i]]]
  y_pos <- -0.009 - 0.005 * i
  data.frame(
    x = c(q[1], q[2], q[3]),
    y = rep(y_pos, 3),
    label = c("Lower: ", "Median: ", "Upper: "),
    value = round(q, 2),
    distribution = factor(dist_names[i], levels = dist_names)
  )
}))

# 13. Raw data summary
raw_q <- quantile(ref_ret_u_vec, c(0.025, 0.5, 0.975))
raw_y <- -0.0075  # slightly closer to density curve
raw_summary <- data.frame(
  x = raw_q[2], xmin = raw_q[1], xmax = raw_q[3],
  y = raw_y
)
raw_annot <- data.frame(
  x = c(raw_q[1], raw_q[2], raw_q[3]),
  y = rep(raw_y + 0.001, 3),
  label = c("Lower: ", "Median: ", "Upper: "),
  value = round(raw_q, 2)
)

# 14. Plot
plot_final <- ggplot() +
  geom_density(data = ref_ret_u_df, aes(x = ret_u), fill = "grey60", colour = NA, alpha = 0.5) +
  geom_line(data = dens_df, aes(x = x, y = y, color = distribution), size = 1) +
  geom_errorbarh(data = raw_summary, aes(xmin = xmin, xmax = xmax, y = y), height = 0, color = "grey40", size = 1) +
  geom_point(data = raw_summary, aes(x = x, y = y), color = "grey40", size = 3) +
  geom_errorbarh(data = summary_df, aes(xmin = xmin, xmax = xmax, y = y, color = distribution), height = 0) +
  geom_point(data = summary_df, aes(x = x, y = y, color = distribution), size = 3) +
  geom_text(data = annotation_df, aes(x = x, y = y + 0.001, label = paste0(label, value), color = distribution), size = 2, vjust = 0) +
  geom_text(data = raw_annot, aes(x = x, y = y, label = paste0(label, value)), color = "grey40", size = 2, vjust = 0) +
  scale_color_manual(values = dist_colors, breaks = dist_names, labels = sapply(seq_along(dist_names), function(i) {
    paste0(dist_names[i], " (", params[[dist_names[i]]], ", LL=", round(logliks[i], 1), ")")
  })) +
  labs(x = "Mean duration of use (months)", y = "Density", color = "Distribution (Params, LL)") +
  theme_minimal()

# 15. Save

ggsave("ref_ret_u_fit_comparison_final.png", plot = plot_final, width = 9, height = 7.5, dpi = 300, units = "in", bg = "white")