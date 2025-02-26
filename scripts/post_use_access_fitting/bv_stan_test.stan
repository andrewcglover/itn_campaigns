//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] mu_a;
  vector[N] sigma_a;
  vector[N] mu_n;
  vector[N] sigma_n;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  // vector<lower = 0, upper = 1>[N] a;
  // vector<lower = 0, upper = 1>[N] n;
  // real<lower = 0> beta;
  // real<lower = 0, upper = 1> gamma;
  // real<lower = 0> sigma;
  // vector<lower = 0, upper = 1>[N] a;
  vector<lower = 0>[N] a;
  vector<lower = 0>[N] n;
  real<lower = 0> beta;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  beta ~ normal(1, 0.5);
  gamma ~ normal(0.5, 0.5);
  sigma ~ normal(0.5, 0.5);
  a ~ normal(mu_a, sigma_a);
  n ~ normal(mu_n, sigma_n);
  a ~ normal(beta * n ^ gamma, sigma);
}
