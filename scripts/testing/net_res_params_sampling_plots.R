plot(dat_res_pbo_ll$resistance[1:101],dat_res_pbo_ll$dn0[1:101])
plot(dat_res_pbo_ll$resistance[1:101],dat_res_pbo_ll$rn0[1:101])
plot(dat_res_pbo_ll$resistance[1:101],dat_res_pbo_ll$gamman[1:101])




test_draws <- 4

res_draw_test <- data.frame("Month" = rep(seq(1, N_CMC_sim), test_draws * 2),
                            "Samping" = c(rep("sequential",
                                              test_draws * N_CMC_sim),
                                          rep("random",
                                              test_draws * N_CMC_sim)),
                            "dn0" = rep(0, N_CMC_sim * test_draws * 2),
                            "rn0" = rep(0, N_CMC_sim * test_draws * 2),
                            "gamman" = rep(0, N_CMC_sim * test_draws * 2))


for (jj in 1:test_draws) {
  
  old_dat_res_sample <- old_dat_res[0,]
  new_dat_res_sample <- new_dat_res[0,]
  
  for (kk in 1:N_CMC_sim) {
    
    x <- (jj-1) * N_CMC_sim + kk
    y <- jj
    
    old_sample_kk <- old_dat_res %>%
      dplyr::filter(resistance == round_monthly_res[kk]) %>%
      dplyr::filter(draw == y)
    if (dim(old_sample_kk)[1] < 1) {
      print(paste("Warning:",
                  round_monthly_res[kk],
                  "resistance match fail"))
    }
    old_dat_res_sample %<>% rbind.data.frame(old_sample_kk)
    
    new_sample_kk <- new_dat_res %>%
      dplyr::filter(resistance == round_monthly_res[kk]) %>%
      dplyr::filter(draw == y)
    if (dim(new_sample_kk)[1] < 1) {
      print(paste("Warning:",
                  round_monthly_res[kk],
                  "resistance match fail"))
    }
    new_dat_res_sample %<>% rbind.data.frame(new_sample_kk)
    
  }
  
  comb_dat_res_sample <- new_dat_res_sample
  comb_dat_res_sample[1:N_CMC_old_nets,] <- old_dat_res_sample[1:N_CMC_old_nets,]
  
  i0 <- 1 + (jj - 1) * N_CMC_sim
  i1 <- i0 + N_CMC_sim - 1
  
  res_draw_test$dn0[i0:i1] <- comb_dat_res_sample$dn0
  res_draw_test$rn0[i0:i1] <- comb_dat_res_sample$rn0
  res_draw_test$gamman[i0:i1] <- comb_dat_res_sample$gamman
  
  
}



for (jj in 1:test_draws) {
  
  old_dat_res_sample <- old_dat_res[0,]
  new_dat_res_sample <- new_dat_res[0,]
  
  for (kk in 1:N_CMC_sim) {
    
    x <- (jj-1) * N_CMC_sim + kk
    y <- long_sample_ids[(x-1) %% length(long_sample_ids) + 1]

    old_sample_kk <- old_dat_res %>%
      dplyr::filter(resistance == round_monthly_res[kk]) %>%
      dplyr::filter(draw == y)
    if (dim(old_sample_kk)[1] < 1) {
      print(paste("Warning:",
                  round_monthly_res[kk],
                  "resistance match fail"))
    }
    old_dat_res_sample %<>% rbind.data.frame(old_sample_kk)
    
    new_sample_kk <- new_dat_res %>%
      dplyr::filter(resistance == round_monthly_res[kk]) %>%
      dplyr::filter(draw == y)
    if (dim(new_sample_kk)[1] < 1) {
      print(paste("Warning:",
                  round_monthly_res[kk],
                  "resistance match fail"))
    }
    new_dat_res_sample %<>% rbind.data.frame(new_sample_kk)
    
  }
  
  comb_dat_res_sample <- new_dat_res_sample
  comb_dat_res_sample[1:N_CMC_old_nets,] <- old_dat_res_sample[1:N_CMC_old_nets,]
  
  i0 <- 1 + ((jj + test_draws) - 1) * N_CMC_sim
  i1 <- i0 + N_CMC_sim - 1
  
  res_draw_test$dn0[i0:i1] <- comb_dat_res_sample$dn0
  res_draw_test$rn0[i0:i1] <- comb_dat_res_sample$rn0
  res_draw_test$gamman[i0:i1] <- comb_dat_res_sample$gamman
  
}

res_draw_test_flat <- data.frame(
  "Month" = rep(seq(1, N_CMC_sim), test_draws * 6),
  "Sampling" = rep(c(rep("sequential", test_draws * N_CMC_sim),
                    rep("random", test_draws * N_CMC_sim)), 3),
  "Sample" = rep(rep(seq(1, test_draws), each = N_CMC_sim), 2),
  "Variable" = rep(NA, N_CMC_sim * test_draws * 6),
  "Value" = rep(NA, N_CMC_sim * test_draws * 6)
)

i0 <- 1
i1 <- i0 + N_CMC_sim * test_draws * 2 - 1
res_draw_test_flat$Variable[i0:i1] <- "dn0"
res_draw_test_flat$Value[i0:i1] <- res_draw_test$dn0
i0 <- i1 + 1
i1 <- i0 + N_CMC_sim * test_draws * 2 - 1
res_draw_test_flat$Variable[i0:i1] <- "rn0"
res_draw_test_flat$Value[i0:i1] <- res_draw_test$rn0
i0 <- i1 + 1
i1 <- i0 + N_CMC_sim * test_draws * 2 - 1
res_draw_test_flat$Variable[i0:i1] <- "gamman"
res_draw_test_flat$Value[i0:i1] <- res_draw_test$gamman



ggplot(res_draw_test_flat) +
  geom_point(aes(x = Month,
                 y = Value,
                 colour = as.factor(Sample),
                 fill = as.factor(Sample)),
             alpha = 0.2,
             stroke = NA,
             #color = "transparent",
             size = 2,
             shape = 21) +
  facet_grid(rows = vars(Variable),
             cols = vars(Sampling),
             scales = "free_y")


# Timing comparison
tic()
for (jj in 1:test_draws) {
  
  old_dat_res_sample <- old_dat_res[0,]
  new_dat_res_sample <- new_dat_res[0,]
  
  for (kk in 1:N_CMC_sim) {
    
    x <- (jj-1) * N_CMC_sim + kk
    y <- jj
    
    old_sample_kk <- old_dat_res %>%
      dplyr::filter(resistance == round_monthly_res[kk]) %>%
      dplyr::filter(draw == y)
    if (dim(old_sample_kk)[1] < 1) {
      print(paste("Warning:",
                  round_monthly_res[kk],
                  "resistance match fail"))
    }
    old_dat_res_sample %<>% rbind.data.frame(old_sample_kk)
    
    new_sample_kk <- new_dat_res %>%
      dplyr::filter(resistance == round_monthly_res[kk]) %>%
      dplyr::filter(draw == y)
    if (dim(new_sample_kk)[1] < 1) {
      print(paste("Warning:",
                  round_monthly_res[kk],
                  "resistance match fail"))
    }
    new_dat_res_sample %<>% rbind.data.frame(new_sample_kk)
    
  }
  
  comb_dat_res_sample <- new_dat_res_sample
  comb_dat_res_sample[1:N_CMC_old_nets,] <- old_dat_res_sample[1:N_CMC_old_nets,]
  
}
toc()

tic()
for (jj in 1:test_draws) {
  
  comb_dat_res_sample <- old_dat_res[0,]

  for (kk in 1:N_CMC_sim) {
    
    x <- (jj-1) * N_CMC_sim + kk
    y <- jj
    
    if (kk < N_CMC_old_nets) {
      sample_kk <- old_dat_res %>%
        dplyr::filter(resistance == round_monthly_res[kk]) %>%
        dplyr::filter(draw == y)
      if (dim(old_sample_kk)[1] < 1) {
        print(paste("Warning:",
                    round_monthly_res[kk],
                    "resistance match fail"))
      }
    } else {
      sample_kk <- new_dat_res %>%
        dplyr::filter(resistance == round_monthly_res[kk]) %>%
        dplyr::filter(draw == y)
      if (dim(new_sample_kk)[1] < 1) {
        print(paste("Warning:",
                    round_monthly_res[kk],
                    "resistance match fail"))
      }
    }
    comb_dat_res_sample %<>% rbind.data.frame(sample_kk)
  }
    
}
toc()


tic()
for (jj in 1:test_draws) {
  
  comb_dat_res_sample <- old_dat_res[0,]
  
  res_kk <- old_dat_res %>%
    dplyr::filter(resistance == round_monthly_res[1])
  
  for (kk in 1:N_CMC_sim) {
    
    x <- (jj-1) * N_CMC_sim + kk
    y <- jj
    
    if (kk < N_CMC_old_nets) {
      
      if (kk > 1) {
        if (round_monthly_res[kk] != round_monthly_res[kk-1]) {
          res_kk <- old_dat_res %>%
            dplyr::filter(resistance == round_monthly_res[kk])
        }
      }
      
      sample_kk <- res_kk %>%
        dplyr::filter(draw == y)
      
      if (dim(sample_kk)[1] < 1) {
        print(paste("Warning:",
                    round_monthly_res[kk],
                    "resistance match fail"))
        
      }
      
    } else {
      
      if (kk > 1) {
        if (round_monthly_res[kk] != round_monthly_res[kk-1]) {
          res_kk <- new_dat_res %>%
            dplyr::filter(resistance == round_monthly_res[kk])
        }
      }
      
      sample_kk <- res_kk %>%
        dplyr::filter(draw == y)
      
      if (dim(sample_kk)[1] < 1) {
        print(paste("Warning:",
                    round_monthly_res[kk],
                    "resistance match fail"))
      }
    }
    comb_dat_res_sample %<>% rbind.data.frame(sample_kk)
  }
  
}
toc()


tic()
for (jj in 1:test_draws) {
  
  comb_dat_res_sample <- old_dat_res[0,]
  
  res_kk <- old_dat_res[old_dat_res$resistance == round_monthly_res[1],]
  
  for (kk in 1:N_CMC_sim) {
    
    x <- (jj-1) * N_CMC_sim + kk
    y <- long_sample_ids[(x-1) %% length(long_sample_ids) + 1]
    
    if (kk < N_CMC_old_nets) {
      
      if (kk > 1) {
        if (round_monthly_res[kk] != round_monthly_res[kk-1]) {
          res_kk <- old_dat_res[old_dat_res$resistance == round_monthly_res[kk],]
        }
      }
      
      sample_kk <- res_kk[res_kk$draw == y,]
      
      if (dim(sample_kk)[1] < 1) {
        print(paste("Warning:",
                    round_monthly_res[kk],
                    "resistance match fail"))
        
      }
      
    } else {
      
        if (kk > 1) {
          if (round_monthly_res[kk] != round_monthly_res[kk-1]) {
            res_kk <- new_dat_res[new_dat_res$resistance == round_monthly_res[kk],]
          }
        }
        
        sample_kk <- res_kk[res_kk$draw == y,]
      
      if (dim(sample_kk)[1] < 1) {
        print(paste("Warning:",
                    round_monthly_res[kk],
                    "resistance match fail"))
      }
    }
    comb_dat_res_sample %<>% rbind.data.frame(sample_kk)
  }
  
  i0 <- 1 + ((jj + test_draws) - 1) * N_CMC_sim
  i1 <- i0 + N_CMC_sim - 1
  
  res_draw_test$dn0[i0:i1] <- comb_dat_res_sample$dn0
  res_draw_test$rn0[i0:i1] <- comb_dat_res_sample$rn0
  res_draw_test$gamman[i0:i1] <- comb_dat_res_sample$gamman
  
}
toc()