bv_access_npc <- read.csv("./data/fig_4_access_npc.csv")

lo <- loess(access_mean ~ percapita_nets_mean,
                          data = bv_access_npc,
            span = 0.75)

lo.x <- seq(min(bv_access_npc$percapita_nets_mean),
            max(bv_access_npc$percapita_nets_mean),
            length.out = 50)
lo.y <- predict(lo, lo.x)

slp.y <- seq(min(bv_access_npc$access_mean),
             max(bv_access_npc$access_mean),
             length.out = 50)
slp.x <- slp.y / 1.8

xval <- approx(x = lo$fitted, y = lo$x, xout = 0.6)$y

plot(bv_access_npc$percapita_nets_mean,
     bv_access_npc$access_mean,
     pch = 16,
     col = "grey50",
     )
lines(slp.x, slp.y, lwd = 2)
lines(lo.x, lo.y, col = "cornflowerblue", lwd = 4)
points(xval, 0.6, col = "blue", lwd = 5)