
extracted_surveys <- get_net_data(cc = SSA_ISO2,
                                  start_year = 2000,
                                  end_year = 2024)

SSA_ISO2_m_CI <- all_SSA_ISO2[all_SSA_ISO2 != "CI"]


# Extract data
extracted_surveys_2008_2022 <- get_net_data(cc = SSA_ISO2,
                                            start_year = 2008,
                                            end_year = 2022)

saveRDS(extracted_surveys_2008_2022, "extracted_surveys_2008_2022.rds")

extracted_surveys_2008_2024 <- get_net_data(cc = SSA_ISO2,
                                            start_year = 2008,
                                            end_year = 2024)
saveRDS(extracted_surveys_2008_2024, "extracted_surveys_2008_2024.rds")

extracted_surveys_SSA_2008_2024 <- get_net_data(cc = all_SSA_ISO2,
                                                start_year = 2008,
                                                end_year = 2024)
saveRDS(extracted_surveys_SSA_2008_2024, "extracted_surveys_SSA_2008_2024.rds")


# Extract data
extracted_surveys_2000_2022 <- get_net_data(cc = SSA_ISO2,
                                            start_year = 2000,
                                            end_year = 2022)

saveRDS(extracted_surveys_2000_2022, "extracted_surveys_2000_2022.rds")

extracted_surveys_2000_2024 <- get_net_data(cc = SSA_ISO2,
                                            start_year = 2000,
                                            end_year = 2024)
saveRDS(extracted_surveys_2000_2024, "extracted_surveys_2000_2024.rds")

extracted_surveys_SSA_2000_2024 <- get_net_data(cc = SSA_ISO2_m_CI,
                                                start_year = 2000,
                                                end_year = 2024)
saveRDS(extracted_surveys_SSA_2000_2024, "extracted_surveys_SSA_2000_2024.rds")

# N_test <- dim(test)[1]
# test_sum <- test %>%
#   dplyr::group_by(hml16) %>%
#   dplyr::summarise(prev = sum(hml35==1) / (sum(hml35<=1)))

# Extract data
extracted_surveys <- get_net_data(cc = SSA_ISO2,
                                  start_year = first_year,
                                  end_year = final_year)

# Remove any corrupted surveys
retained_surveys <- !(names(extracted_surveys) %in% corrupted_surveys)
extracted_surveys <- extracted_surveys[retained_surveys]