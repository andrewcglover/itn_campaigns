#routine_vs_mass_coverage.R

# oo <- options(repos = "https://cran.r-project.org/")
# utils::install.packages("Matrix")
# utils::install.packages("lme4")
# options(oo)

library(lme4)

mass_camp_df <- net_data %>% filter(months_post_mdc == 0)

final_mass_camp_df <- mass_camp_df %>% group_by(area) %>%
  filter(CMC == max(CMC))

final_mass_camp_df$Urbanicity <- final_mass_camp_df$urbanicity
final_mass_camp_df$Country <- countrycode(final_mass_camp_df$ISO2,
                                             "iso2c","country.name")


# cov_mod <- lmer(D_u_mean ~ P_u_mean + (P_u_mean | ISO2),
#                 data = final_mass_camp_df)
cov_mod <- lmer(D_u_mean ~ 0 + P0_u_mean + (P0_u_mean | ISO2),
                data = final_mass_camp_df)

cov_mod_test <- lm(D_u_mean ~ P_u_mean, data = final_mass_camp_df)

Number_of_boots <- 1000

# Extract the fixed effect coefficients.
FE_df <- fixef(cov_mod) %>% 
  t() %>%
  as.data.frame()

cov_df <- data.frame("P_u_mean" = seq(round(min(final_mass_camp_df$P0_u_mean),2),
                                      round(max(final_mass_camp_df$P0_u_mean),2),
                                      0.01))
cov_df$D_u_mean = cov_df$P_u_mean * FE_df[[1]]

#predict(cov_mod, cov_df)



alpha_val <- 0.5
ggplot() +
  geom_errorbar(data = final_mass_camp_df,
                aes(x = P0_u_mean, ymin = D_u_LB1, ymax = D_u_UB1,
                    colour = Country),
                alpha = alpha_val) +
  geom_errorbarh(data = final_mass_camp_df,
                 aes(y = D_u_mean, xmin = P0_u_LB1, xmax = P0_u_UB1,
                     colour = Country),
                 alpha = alpha_val) +
  geom_point(data = final_mass_camp_df,
             aes(y = D_u_mean,
                 x = P0_u_mean,
                 colour = Country,
                 shape = Urbanicity),
             alpha = 0.8) +
  geom_path(data = cov_df, aes(x = P_u_mean, y = D_u_mean)) +
  theme_bw() +
  ylab("Proportion using routine ITNs") +
  xlab("Proportion using any ITN") +
  ggtitle("Usage following the last mass campaign prior to 2023")