// hier_net_decay.stan
// Cloned from exp_model_d2.stan on 28/10/23

data {
  int<lower = 1> N;                      // number of observations
  int<lower = 1> N_a;                    // number of admins
  int<lower = 1> N_c;                    // number of countries
  int<lower = 1, upper = N_a> a[N];      // admin ids
  int<lower = 1, upper = N_c> c[N];      // country ids
  int<lower = 1, upper = N_c> cc[N_a];   // country id by admin
  int<lower = 0> m[N];                   // months since net obtained 
}

parameters {
  real<lower = 0> mu_u;
  real<lower = 0> sigma_u;
  real<lower = 0> tau_u;
  real<lower = 0> rho_u;
  vector<lower = 0>[N_c] mu_c;
  vector<lower = 0>[N_c] sigma_c;
  vector<lower = 0>[N_a] inv_lambda;
}

model {
  // hyperhyperpriors
  mu_u ~ normal(24,12);
  sigma_u ~ normal(1,1);
  tau_u ~ normal(1,1);
  rho_u ~ normal(1,1);
  // hyperpriors
  mu_c[] ~ normal(mu_u, sigma_u);
  sigma_c[] ~ normal(tau_u, rho_u);
  // prior
  inv_lambda[] ~ normal(mu_c[cc], sigma_c[cc]);
  // likelihood
  for (i in 1:N) {
    if (m[i] > 36) {
      target += exponential_lccdf(36.5 | 1.0 / inv_lambda[a[i]]);
    } else {
      target += exponential_lpdf(m[i] | 1.0 / inv_lambda[a[i]]);
    }
  }
}