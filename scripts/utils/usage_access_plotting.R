# usage_access_plotting.R

plot_usage <- function(dataset, cc = "SN") {
  
  # Data subset for country
  cc_dataset <- dataset[which(dataset$ISO2 == cc),]
  
  # Plot
  
  ylbs_ids <- seq(CMC_first,CMC_last,24)
  ylbs <- CMC_to_date(ylbs_ids)[,1]
  
  baseplt <- ggplot() +
    geom_ribbon(data = cc_dataset,
                aes(x = CMC, ymax = D_u_UB, ymin = D_u_LB, fill = urbanicity),
                alpha = 0.3) +
    geom_path(data = cc_dataset,
              aes(x = CMC, y = D_u_mean, color = urbanicity),
              linetype = "dashed") +
    # geom_ribbon(data = country_df,
    #             aes(x = CMC, ymax = p0_UB, ymin = p0_LB, fill = urbanicity),
    #             alpha = 0.3) +
    # geom_path(data = country_df,
    #           aes(x = CMC, y = mean_p0, color = urbanicity),
  #           linetype = "dotted") +
  geom_ribbon(data = cc_dataset,
              aes(x = CMC, ymax = P_u_UB, ymin = P_u_LB, fill = urbanicity),
              alpha = 0.2) +
  geom_path(data = cc_dataset,
            aes(x = CMC, y = P_u_mean, color = urbanicity)) +
  geom_count(data = cc_dataset,
             aes(x = CMC, y = prop_used, size = total, color = urbanicity),
             alpha = 0.4, stroke = NA) +
    ylab("Surveyed usage") + 
    xlab("Year") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = ylbs_ids, labels = ylbs, limits = c(date_to_CMC(2014,1), date_to_CMC(2021,1))) +
    #scale_x_continuous(breaks = ylbs_ids, labels = ylbs, limits = c(date_to_CMC(2013,6), date_to_CMC(2021,6))) +
    #scale_x_continuous(breaks = ylbs_ids, labels = ylbs, limits = c(date_to_CMC(2010,6), date_to_CMC(2020,6))) +
    scale_x_continuous(breaks = ylbs_ids, labels = ylbs, limits = c(date_to_CMC(2010,1), date_to_CMC(2021,6))) +
    theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1))
  
  baseplt + facet_wrap(~ADM1, nrow = 2)
}