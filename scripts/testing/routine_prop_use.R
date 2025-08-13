library(geomtextpath)

avg_use_access <- function(D, C0, lambda, del_t) {
  avg_val <- D + ( (C0 / (lambda * del_t) ) * (1 - exp(-lambda * del_t) ) )
}

end_use_access <- function(D, C0, lambda, del_t) {
  end_val <- D + ( C0 * exp(-lambda * del_t) )
}

sf <- 100

ref_data <- net_data %>% filter(CMC == ref_CMC)

interval_mn <- 36

ref_data$avg_u_mean <- avg_use_access(
  D = ref_data$D_u_mean,
  C0 = ref_data$P0_u_mean - ref_data$D_u_mean,
  lambda = 1/ref_data$invlam_u_mean,
  del_t = rep(interval_mn, dim(ref_data)[1])
)

ref_data$avg_u_LB <- avg_use_access(
  D = ref_data$D_u_LB1,
  C0 = ref_data$P0_u_LB1 - ref_data$D_u_LB1,
  lambda = 1/ref_data$invlam_u_LB1,
  del_t = rep(interval_mn, dim(ref_data)[1])
)

ref_data$avg_u_UB <- avg_use_access(
  D = ref_data$D_u_UB1,
  C0 = ref_data$P0_u_UB1 - ref_data$D_u_UB1,
  lambda = 1/ref_data$invlam_u_UB1,
  del_t = rep(interval_mn, dim(ref_data)[1])
)

ref_data$end_u_mean <- end_use_access(
  D = ref_data$D_u_mean,
  C0 = ref_data$P0_u_mean - ref_data$D_u_mean,
  lambda = 1/ref_data$invlam_u_mean,
  del_t = rep(interval_mn, dim(ref_data)[1])
)

ref_data$end_u_LB <- end_use_access(
  D = ref_data$D_u_LB1,
  C0 = ref_data$P0_u_LB1 - ref_data$D_u_LB1,
  lambda = 1/ref_data$invlam_u_LB1,
  del_t = rep(interval_mn, dim(ref_data)[1])
)

ref_data$end_u_UB <- end_use_access(
  D = ref_data$D_u_UB1,
  C0 = ref_data$P0_u_UB1 - ref_data$D_u_UB1,
  lambda = 1/ref_data$invlam_u_UB1,
  del_t = rep(interval_mn, dim(ref_data)[1])
)

ref_data$Country <- ref_data$ISO2 %>% countrycode("iso2c", "country.name")

use.lm <- lm(D_u_mean ~ 0 + avg_u_mean, ref_data)

use.lm.df <- data.frame("ISO2" = SSA_ISO2,
                        "Country" = countrycode(SSA_ISO2, "iso2c", "country.name"),
                        "slope" = rep(NA, length(SSA_ISO2)))

for (i in 1:length(SSA_ISO2)) {
  ctry_data <- ref_data %>% filter(ISO2 == SSA_ISO2[i])
  ctry.lm <- lm(D_u_mean ~ 0 + avg_u_mean, ctry_data)
  use.lm.df$slope[i] <- coef(ctry.lm)[["avg_u_mean"]]
}
                        

  
  
  lm(D_u_mean ~ 0 + avg_u_mean, ref_data)


ggplot(data = ref_data,
       aes(x = avg_u_mean * sf,
           y = D_u_mean * sf,
           colour = Country)) +
  geom_pointrange(aes(ymin = D_u_LB1 * sf, ymax = D_u_UB1 * sf)) +
  geom_errorbarh(aes(xmin = avg_u_LB * sf, xmax = avg_u_UB * sf, height = 0)) +
  #geom_abline(slope = coef(use.lm)[["avg_u_mean"]]) +
  #geom_textabline(slope = coef(use.lm)[["avg_u_mean"]], label = coef(use.lm)[["avg_u_mean"]],
  #                gap = FALSE, offset = unit(0.2, "lines"), text_only = TRUE) +
  geom_textabline(aes(slope = slope, intercept = 0,
                  label = paste("routine =", round(slope,3), "* overall"),
                  color = Country),
                  data = use.lm.df, hjust = 0.9,
                  size = 3) + #, #vjust = -0.2) +
  geom_textabline(slope = coef(use.lm)[["avg_u_mean"]],
                  label = paste("routine =", round(coef(use.lm)[["avg_u_mean"]],3), "* overall"),
                  color = "black", hjust = 0.9, size = 3) +
  scale_x_continuous(breaks = seq(0,1,0.1)*sf,limits = c(0, 1)*sf) +
  scale_y_continuous(breaks = seq(0,1,0.1)*sf,limits = c(0, 0.5)*sf) +
  theme_bw() +
  xlab("Mean overall ITN use with three-year campaigns (%)") +
  ylab("Use of routine ITNs (%)")
  #labs(colour = "Country")
  



use.lm0 <- lm(D_u_mean ~ 0 + P0_u_mean, ref_data)


ggplot(data = ref_data,
       aes(x = P0_u_mean * sf,
           y = D_u_mean * sf,
           colour = Country)) +
  geom_pointrange(aes(ymin = D_u_LB1 * sf, ymax = D_u_UB1 * sf)) +
  geom_errorbarh(aes(xmin = P0_u_LB1 * sf, xmax = P0_u_UB1 * sf, height = 0)) +
  geom_abline(slope = coef(use.lm0)[["P0_u_mean"]]) +
  scale_x_continuous(breaks = seq(0,1,0.1)*sf,limits = c(0, 1)*sf) +
  scale_y_continuous(breaks = seq(0,1,0.1)*sf,limits = c(0, 0.5)*sf) +
  theme_bw() +
  xlab("Mean overall ITN use with three-year campaigns (100%)") +
  ylab("Use of routine ITNs (100%)")
#labs(colour = "Country")
  
  #geom_abline(slope = coef(use.lm)[["avg_u_mean"]], 
  #            intercept = coef(use.lm)[["(Intercept)"]])


use.lm.end <- lm(D_u_mean ~ 0 + end_u_mean, ref_data)


ggplot(data = ref_data,
       aes(x = end_u_mean * sf,
           y = D_u_mean * sf,
           colour = Country)) +
  geom_pointrange(aes(ymin = D_u_LB1 * sf, ymax = D_u_UB1 * sf)) +
  geom_errorbarh(aes(xmin = end_u_LB * sf, xmax = end_u_UB * sf, height = 0)) +
  geom_abline(slope = coef(use.lm.end)[["end_u_mean"]]) +
  scale_x_continuous(breaks = seq(0,1,0.1)*sf,limits = c(0, 0.8)*sf) +
  scale_y_continuous(breaks = seq(0,1,0.1)*sf,limits = c(0, 0.5)*sf) +
  theme_bw() +
  xlab("Mean overall ITN use with three-year campaigns (100%)") +
  ylab("Use of routine ITNs (100%)")
#labs(colour = "Country")
