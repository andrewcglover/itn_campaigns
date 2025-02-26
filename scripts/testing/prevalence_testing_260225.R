test_df <- data.frame(
  "ISO2" = net_data$ISO2,
  "admin" = net_data$ADM1,
  "urbanicity" = net_data$urbanicity,
  "CMC" = net_data$CMC,
  "used" = net_data$used,
  "not_used" = (net_data$used / net_data$prop_used) * (1 - net_data$prop_used),
  "prop_used" = net_data$prop_used
)

test_df$lo_prop_used <- qbeta(0.025, test_df$used + 0.5, test_df$not_used + 0.5)
test_df$hi_prop_used <- qbeta(0.975, test_df$used + 0.5, test_df$not_used + 0.5)

ggplot(
  data = test_df %>%
    dplyr::filter(ISO2 == "ML") %>%
    dplyr::filter(!(is.na(prop_used))),
  aes(
    x = CMC,
    y = prop_used,
    ymin = lo_prop_used,
    ymax = hi_prop_used,
    colour = admin
  )
) +
  geom_point() +
  geom_errorbar() +
  facet_grid(
    rows = vars(admin),
    cols = vars(urbanicity)
  )



test_df <- data.frame(
  "ISO2" = net_data$ISO2,
  "admin" = net_data$ADM1,
  "urbanicity" = net_data$urbanicity,
  "CMC" = net_data$CMC,
  "used" = net_data$used,
  "not_used" = (net_data$used / net_data$prop_used) * (1 - net_data$prop_used),
  "prop_used" = net_data$prop_used,
  "rdt_prev" = net_data$rdt_prev,
  "lo_prev" = net_data$lo_prev,
  "hi_prev" = net_data$hi_prev
)

ggplot(
  data = test_df %>%
    dplyr::filter(ISO2 == "BF") %>%
    dplyr::filter(!(is.na(rdt_prev))),
  aes(
    x = CMC,
    y = rdt_prev,
    ymin = lo_prev,
    ymax = hi_prev,
    colour = admin
  )
) +
  geom_point() +
  geom_errorbar() +
  facet_grid(
    rows = vars(admin),
    cols = vars(urbanicity)
  )

saveRDS(net_data,"net_data_250225.rds")

net_data_prev <- readRDS("net_data_250225.rds")

net_prev_only <- net_data_prev %>%
  dplyr::select(
    area, CMC, rdt_pos, rdt_neg, rdt_num, rdt_prev, lo_prev, hi_prev
  )

new_net_data <- dplyr::left_join(
  net_data,
  net_prev_only,
  by = c("area", "CMC")
)

old_net_data <- net_data
net_data <- new_net_data