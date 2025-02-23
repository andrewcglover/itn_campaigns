#editable function for generating mdc plots

generate_mdc_plots <- function(net_data) {
  timestamp <- format(Sys.time(), "%y%m%d%H%M")
  net_data %>% plot_MDCs(densities = "over_comp_nets_norm",
                         periods_dataset = mdc_period_df,
                         colvals = "springgreen4",
                         cap_extreme = FALSE,
                         plot_step_dens = TRUE,
                         plot_smth_dens = TRUE,
                         plot_modes = FALSE,
                         plot_antimodes = FALSE,
                         plot_mdc_pts = TRUE,
                         plot_vert_periods = TRUE,
                         plot_comparison_mdc = TRUE)
  
  net_data %<>% combine_weights(density_name = "rcpt_dhs_w",
                                out_name_from_input = TRUE)
  
  net_data %<>% normalise_area_densities("comb_rcpt_grw_w",
                                         norm_over_net_rec_range = FALSE,
                                         time_unit = "years")
  
  timestamp <- format(Sys.time(), "%y%m%d%H%M")
  net_data %>% plot_MDCs(densities = "ref_nets",
                         periods_dataset = mdc_period_df,
                         colvals = "royalblue",
                         cap_extreme = FALSE,
                         plot_step_dens = TRUE,
                         plot_smth_dens = TRUE,
                         plot_modes = FALSE,
                         plot_antimodes = FALSE,
                         plot_mdc_pts = FALSE,
                         plot_vert_periods = TRUE,
                         plot_comparison_mdc = FALSE)
  
  timestamp <- format(Sys.time(), "%y%m%d%H%M")
  net_data %>% plot_MDCs(densities = c("over_comp_nets_norm","ref_nets_norm"),
                         periods_dataset = mdc_period_df,
                         colvals = c("royalblue","darkorange3"),
                         cap_extreme = FALSE,
                         plot_step_dens = TRUE,
                         plot_smth_dens = TRUE,
                         plot_modes = FALSE,
                         plot_antimodes = FALSE,
                         plot_mdc_pts = TRUE,
                         plot_vert_periods = TRUE,
                         plot_comparison_mdc = TRUE,
                         MDC_density_number = 2)
}