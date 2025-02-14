# bv_a <- bv_access_npc$access_mean
# bv_n <- bv_access_npc$percapita_nets_mean
# 
# log_a <- log(bv_a)
# log_n <- log(bv_n)
# 
# bv_log_fit <- lm(log_a ~ log_n -1)
# 
# #bv_beta <- exp(bv_log_fit$coefficients[[1]])
# bv_gamma <- bv_log_fit$coefficients[[1]]
# 
# plot(bv_n, bv_a, xlim = c(0,1), ylim = c(0,1))
# nn <- seq(0, 1, 0.01)
# #aa <- 1 * nn ^ bv_gamma
# aa <- nn ^ 0.57
# lines(nn, aa)
# lines(nn,1.8*nn)
