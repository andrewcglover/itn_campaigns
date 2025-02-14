library(ggplot2)

exp_loss <- function(u0, lambda, t) {
  u <- u0 * exp(-lambda * t)
}

calc_tau <- function(meanlife, kappa = 20) {
  tau <- meanlife / sqrt(1 - kappa / (kappa - log(0.5)))
}

smooth_comp_loss <- function(v0, kappa = 20, tau, t) {
  v <- v0 * exp(kappa - kappa / (1 - (t / tau)^2))
}

year <- 365
t0 <- 0
t1 <- 4 * year
u0 <- seq(0.5, 1, 0.1)
meanlife <- seq(1, 3, 0.5) * year
t <- seq(t0, t1)
Nt <- length(t)
Nu <- length(u0)
Nl <- length(meanlife)

df <- data.frame(
  "Function" = c(rep("Exponential", Nt * Nu * Nl),
                 rep("Smooth-compact", Nt * Nu *Nl)),
  "Time" = rep(t, Nu * Nl * 2),
  "Initial_nets" = rep(rep(u0, each = Nt), Nl * 2),
  "Central_estimate" = rep(rep(meanlife, each = Nt * Nu), 2),
  "Nets" = rep(NA, Nt * Nu * Nl * 2)
)

id <- 1
for (i in 1:2) {
  for (j in 1:Nl) {
    for (k in 1:Nu) {
      for (l in 1:Nt) {
        if (i == 1) {
          df$Nets[id] <- exp_loss(u0[k], lambda = 1/meanlife[j], t = t[l])
        } else {
          tau <- calc_tau(meanlife = meanlife[j])
          v0 <- u0[k] * 2 / exp(1)
          df$Nets[id] <- smooth_comp_loss(v0 = v0, tau = tau, t = t[l])
        }
        id <- id + 1
      }
    }
  }
}

ggplot(data = df,
       aes(x = Time / year,
           y = Nets,
           colour = Function),) +
  geom_vline(aes(xintercept = Central_estimate / year)) +
  geom_vline(aes(xintercept = Central_estimate * 0.433609 / year),
             linetype = "dashed") +
  geom_hline(aes(yintercept = Initial_nets)) +
  geom_hline(aes(yintercept = Initial_nets * 2 / exp(1)),
             linetype = "dashed") +
  geom_line() +
  facet_grid(cols = vars(Central_estimate / year),
             rows = vars(Initial_nets)) +
  theme_bw()
