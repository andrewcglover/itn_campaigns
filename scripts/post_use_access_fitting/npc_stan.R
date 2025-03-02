#npc_stan.R

stan_npc_fit <- function() {
  bv_N <- dim(bv_access_npc)[1]
  bv_mu_a <- bv_access_npc$access_mean
  bv_sigma_a <- (bv_access_npc$access_upper - bv_access_npc$access_lower) / 3.92
  bv_mu_n <- bv_access_npc$percapita_nets_mean
  bv_sigma_n <- (bv_access_npc$percapita_nets_upper - bv_access_npc$percapita_nets_lower) / 3.92
  
  bv_dat <- list(N = bv_N,
                 mu_a = bv_mu_a,
                 sigma_a = bv_sigma_a,
                 mu_n = bv_mu_n,
                 sigma_n = bv_sigma_n)
  
  bv_iter <- 4000
  bv_warm <- 2000
  bv_chains <- 4
  bv_stan_fit <- stan('./scripts/post_use_access_fitting/bv_stan_test.stan',
                      data = bv_dat,
                      iter = bv_iter,
                      warmup = bv_warm,
                      chains = bv_chains)#,
  #init_r = decay_init_r,
  #control = list(adapt_delta = decay_adapt_delta))
  
  bv_stan_samples <- rstan::extract(bv_stan_fit)
  
  bv_beta <<- bv_stan_samples$beta
  bv_gamma <<- bv_stan_samples$gamma
  bv_n_pred <- seq(0,1,0.001)
  bv_a_pred <- matrix(nrow = length(bv_beta), ncol = length(bv_n_pred))
  for (i in 1:length(bv_beta)) {
    bv_a_pred[i,] <- bv_beta[i] * bv_n_pred ^ bv_gamma[i]
  }
  bv_a_LB <- rep(NA, length(bv_n_pred))
  bv_a_mid <- rep(NA, length(bv_n_pred))
  bv_a_UB <- rep(NA, length(bv_n_pred))
  bv_a_lin <- rep(NA, length(bv_n_pred))
  for (j in 1:length(bv_n_pred)) {
    bv_a_LB[j] <- min(1, quantile(bv_a_pred[,j], 0.025))
    bv_a_mid[j] <- quantile(bv_a_pred[,j], 0.5)
    bv_a_UB[j] <- min(1, quantile(bv_a_pred[,j], 0.975))
    bv_a_lin[j] <- bv_n_pred[j] * 1.8
  }
  
  bv_pred <- data.frame("a_LB" = bv_a_LB,
                        "a_mid" = bv_a_mid,
                        "a_UB" = bv_a_UB,
                        "a_lin" = bv_a_lin,
                        "n" = bv_n_pred)
  
  return(bv_pred)
}

# 
# 
# ggplot() +
#   geom_errorbar(data = bv_access_npc,
#                 aes(x = percapita_nets_mean,
#                     ymin = access_lower,
#                     ymax = access_upper),
#                 width = 0,
#                 colour = "grey70",
#                 alpha = 0.5) +
#   geom_errorbarh(data = bv_access_npc,
#                  aes(xmin = percapita_nets_lower,
#                      xmax = percapita_nets_upper,
#                      y = access_mean),
#                  height = 0,
#                  colour = "grey70",
#                  alpha = 0.5) +
#   geom_point(data = bv_access_npc,
#              aes(x = percapita_nets_mean,
#                  y = access_mean),
#              ) +
#   geom_path(data = bv_pred,
#             aes(x = n,
#                 y = a_lin),
#             linewidth = 1) +
#   geom_smooth(data = bv_access_npc,
#               aes(x = percapita_nets_mean,
#                    y = access_mean),
#               colour = "cornflowerblue",
#               se = FALSE,
#               linewidth = 1) +
#   geom_ribbon(data = bv_pred,
#             aes(x = n,
#                 ymin = a_LB,
#                 ymax = a_UB),
#             colour = NA,
#             fill = "firebrick",
#             alpha = 0.5) +
#   geom_path(data = bv_pred,
#             aes(x = n,
#                 y = a_mid),
#             colour = "firebrick",
#             linewidth = 1) +
#   scale_x_continuous(breaks = seq(0,1,0.1),
#                      limits = c(0,1)) + 
#   scale_y_continuous(breaks = seq(0,1,0.1),
#                      limits = c(0,1)) +
#   xlab("nets per capita") +
#   ylab("access") +
#   theme_bw()
# 
# hist(bv_stan_samples$beta)
# 
# plot(bv_n, bv_a, xlim = c(0,1), ylim = c(0,1))
# for (i in seq(1, (bv_iter - bv_warm) * bv_chains, 50)) {
#   gam <- bv_stan_samples$gamma[i]
#   bet <- bv_stan_samples$beta[i]
#   nn <- seq(0,1,0.001)
#   aa <- bet * nn ^ gam
#   lines(nn,aa)
# }
# 
# ggplot(data = )