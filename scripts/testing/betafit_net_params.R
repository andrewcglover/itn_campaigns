library(extraDistr)

#
# Logistic transformation of the Beta CDF.
#
f.beta <- function(alpha, beta, x, lower=0, upper=1) {
  p <- pbeta((x-lower)/(upper-lower), alpha, beta)
  log(p/(1-p))
}
#
# Sums of squares.
#
delta.beta <- function(fit, actual) sum((fit-actual)^2)
#
# The objective function handles the transformed parameters `theta` and
# uses `f.beta` and `delta` to fit the values and measure their discrepancies.
#
objective.beta <- function(theta, x, prob, ...) {
  ab <- exp(theta) # Parameters are the *logs* of alpha and beta
  fit <- f.beta(ab[1], ab[2], x, ...)
  return (delta.beta(fit, prob))
}
#


x.prob <- c(0.1, 0.5, 0.9)
x.p <- (function(p) log(p/(1-p)))(x.prob)
x <- c(0.126132416,	0.191397819, 0.256531177)
start <- log(c(x[2],1-x[2])*1e2) 

sol <- nlm(objective.beta, start, x=x, prob=x.p, lower=0, upper=1,
           typsize=c(1,1), fscale=1e-12, steptol = 1e-12, gradtol=1e-12)

parms <- exp(sol$estimate) 

curve(pbeta(x, parms[1], parms[2]), n=1001, col="Red")
points(x, x.prob)



x.prob <- c(0.1, 0.5, 0.9)
x.p <- (function(p) log(p/(1-p)))(x.prob)
x <- c(0.699339199,	0.736874872,	0.76083029)
start <- log(c(x[2],1-x[2])*1e2) 

sol <- nlm(objective.beta, start, x=x, prob=x.p, lower=0, upper=1,
           typsize=c(1,1), fscale=1e-12, steptol = 1e-12, gradtol=1e-12)

parms <- exp(sol$estimate) 

curve(pbeta(x, parms[1], parms[2]), n=1001, col="Blue")
points(x, x.prob)



x.prob <- c(0.1, 0.5, 0.9)
x.p <- (function(p) log(p/(1-p)))(x.prob)
y <- c(1.641616652, 1.96516878, 2.281355244)
x <- y / (1 + y)
start <- log(c(x[2],1-x[2])*1e2) 

sol <- nlm(objective.beta, start, x=x, prob=x.p, lower=0, upper=1,
           typsize=c(1,1), fscale=1e-12, steptol = 1e-12, gradtol=1e-12)

parms <- exp(sol$estimate) 

# curve(pbeta(x, parms[1], parms[2]), n=1001, col="Red")
# points(x, x.prob)

plot(seq(0.1,3,0.1),pbetapr(seq(0.1,3,0.1), parms[1], parms[2]), col="Green", type = "l")

points(y, x.prob)









# Solve two problems.
#
par(mfrow=c(1,2))
alpha <- 15; beta <- 22 # The true parameters
for (x in list(c(1e-3, 2e-3), c(1/3, 2/3))) {
  x.p <- f.beta(alpha, beta, x)        # The correct values of the p_i
  start <- log(c(1e1, 1e1))            # A good guess is useful here
  sol <- nlm(objective.beta, start, x=x, prob=x.p, lower=0, upper=1,
             typsize=c(1,1), fscale=1e-12, gradtol=1e-12)
  parms <- exp(sol$estimate)           # Estimates of alpha and beta
  #
  # Display the actual and estimated values.
  #
  print(rbind(Actual=c(alpha=alpha, beta=beta), Fit=parms))
  #
  # Plot the true and estimated CDFs.
  #      
  curve(pbeta(x, alpha, beta), 0, 1, n=1001, lwd=2)
  curve(pbeta(x, parms[1], parms[2]), n=1001, add=TRUE, col="Red")
  points(x, pbeta(x, alpha, beta))
}