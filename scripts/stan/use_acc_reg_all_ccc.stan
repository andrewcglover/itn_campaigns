// use_acc_reg_all_cc.stan
// Usage/access regression model for all countries
// Originally from log_growth_MDC_decay_beta_bin_MDC21_q.stan

functions {
  
  // discretisation function
  vector disc(real LB, real UB, int steps) {
    vector[steps] x;
    real interval = (UB - LB) / steps;
    x[1] = LB + interval / 2.0;
    for (i in 2:steps) {
      x[i] = x[i-1] + interval;
    }
    return x;
  }
  
  // count number of negative values
  int num_negative(vector x) {
    int N = rows(x);
    int y[N];
    for (i in 1:N) {
      if(x[i] < 0) {
        y[i] = 1;
      } else {
        y[i] = 0;
      }
    }
    int t = sum(y);
    return t;
  }
  
  // calculate proportion of campaign nets m months after a mass campaign
  vector q_fun(int N, int N_disc, vector q0, matrix m_p, matrix m, vector inv_lambda, int[] a) {
    vector[N] q = rep_vector(0,N);
    for (i in 1:N) {
      for (j in 1:N_disc) {
        real m_decay = q0[a[i]] * exp(-m[j,i] / inv_lambda[a[i]]);
        q[i] = q[i] + m_p[j,i] * m_decay / (m_decay + (1-q0[a[i]]));
      }
    }
    return q;
  }
}

data {
  int<lower = 0> N_c;               //number of countries
  int<lower = 0> N_a;               //number of admin units across all countries
  int<lower = 0> N_t;               //number of time points
  int<lower = 0, upper = N_a> cc[N_a]; //country id by admin
  vector<lower = 0>[N_a] mu_n;      //mean net life
  vector<lower = 0>[N_a] sigma_n;   //sd net life
  int<lower = 0> N;                 //number of trials/time points over admins
  int<lower = 1, upper = N_c> c[N]; //countries
  int<lower = 1, upper = N_a> a[N]; //admin one units
  vector<lower = 0>[N] t;           //time (Calendar Month Code)
  int<lower = 0> u[N];              //number used for each trial
  int<lower = 0> n[N];              //total for each trial
  int<lower = 0> s[N];              //number with source of net recorded
  int<lower = 0> z[N];              //number with campaign as source
  int<lower = 0> max_rho;
  matrix[N_a, max_rho] r_meas;
  matrix[N_a, max_rho] r_tau;
  int<lower = 0> rho[N];
  real<lower = 0> disc_rnge;
  int<lower = 1> N_disc;
  real<lower = 0> max_m;
  real no_round;
}

transformed data {
  real<lower = 0> t_mean = mean(t);
  real<lower = 0> t_sd = sd(t);
  real hlf_rnge = disc_rnge / 2.0;
  real max_m_hat = max_m / t_sd;
  vector[N] t_hat = (t - t_mean) / t_sd;
  matrix[N_a, max_rho] r_meas_hat = (r_meas - t_mean) / t_sd;
  matrix[N_a, max_rho] r_tau_hat = r_tau / t_sd;
  matrix[max_rho, N_disc] r_hat[N_a];
  matrix[N_disc, N] m_hat;
  matrix[N_disc, N] m_hat_p;
  matrix[N_a, max_rho] r_hat_LB = r_meas_hat - r_tau_hat * hlf_rnge;
  matrix[N_a, max_rho] r_hat_UB = r_meas_hat + r_tau_hat * hlf_rnge;
  for (i in 1:N_a) {
    for (j in 1:max_rho) {
      if (r_meas[i,j] < 0) {
        r_hat[i,j,] = rep_row_vector(no_round, N_disc);
      } else {
        r_hat[i,j,] = to_row_vector(disc(r_hat_LB[i,j], r_hat_UB[i,j], N_disc));
      }
    }
  }
  vector[N_disc] SN_z = disc(-hlf_rnge, hlf_rnge, N_disc);
  vector[N_disc] SN_den;
  for (i in 1:N_disc) {
    SN_den[i] = exp(std_normal_lpdf(SN_z[i]));
  }
  real norm_cnst = sum(SN_den);
  vector[N_disc] prpr_SN_den = SN_den / norm_cnst;
  for (j in 1:N) {
    int j1;
    int j2;
    if (rho[j] == max_rho) {
      j1 = rho[j];
      j2 = rho[j] - 1;
    } else if (r_meas[a[j],rho[j]+1] < 0) {
      j1 = rho[j];
      j2 = rho[j] - 1;
    } else if ((t[j] - r_meas[a[j],rho[j]]) < (r_meas[a[j],rho[j]+1] - t[j])) {
      j1 = rho[j];
      j2 = rho[j] - 1;
    } else {
      j1 = rho[j] + 1;
      j2 = rho[j];
    }
    m_hat[,j] = t_hat[j] - to_vector(r_hat[a[j],j1,]);
    m_hat_p[,j] = prpr_SN_den;
    int N_neg = num_negative(m_hat[,j]);
    int N_pos = N_disc - N_neg;
    int N_posp1 = N_pos + 1;
    real pos_p_sum = 0;
        if (N_pos > 0) {
          pos_p_sum = sum(m_hat_p[1:N_pos,j]);
        }
    if (N_neg > 0) {
      //print("N_neg = ", N_neg, "; N_pos = ", N_pos, "; N_posp1 = ", N_posp1, "; N_disc = ", N_disc, "; j = ", j, "; j2 = ", j2);
      if (j2 > 0) {
        m_hat[N_posp1:N_disc,j] = t_hat[j] - to_vector(disc(r_hat_LB[a[j],j2], r_hat_UB[a[j],j2], N_neg));
        vector[N_neg] SN_z_prev = disc(-hlf_rnge, hlf_rnge, N_neg);
        vector[N_neg] SN_den_prev;
        for (k in 1:N_neg) {
          SN_den_prev[k] = exp(std_normal_lpdf(SN_z_prev[k]));
        }
        real norm_cnst_prev = sum(SN_den_prev) / (1.0 - pos_p_sum);
        m_hat_p[N_posp1:N_disc,j] = SN_den_prev / norm_cnst_prev;
      } else {
        m_hat[N_posp1:N_disc,j] = rep_vector(max_m_hat, N_neg);
        m_hat_p[N_posp1:N_disc,j] = rep_vector((1.0 - pos_p_sum / (1.0 * N_neg)), N_neg);
      }
    }
  }
}

parameters {
  vector<lower = 0>[N_a] inv_lambda;
  vector<lower = -1e3, upper = 16>[N_a] beta_0;
  vector<lower = 0, upper = 16>[N_a] beta_t;
  vector<lower = 0, upper = 1>[N_a] k;
  vector<lower = 0, upper = 1>[N_a] q0;
  vector<lower = 0>[N_a] alpha0;
}

transformed parameters {
  vector[N] inv_lambda_hat = inv_lambda / t_sd;
  vector[N] beta_hat = beta_0[a] + beta_t[a] .* t_hat;
  vector[N] P0 = k[a] .* inv_logit(beta_hat);
  vector[N] C0 = q0[a] .* P0;
  vector[N] D = P0 - C0;
  vector[N] C = rep_vector(0, N);
  for (i in 1:N) {
    for (j in 1:N_disc) {
      C[i] = C[i] + C0[i] * m_hat_p[j,i] * exp(-m_hat[j,i] / inv_lambda_hat[a[i]]);
    }
  }
  vector[N] P = C + D;
}

model {
  inv_lambda ~ normal(mu_n, sigma_n);
  beta_0[] ~ normal(0, 0.1);
  beta_t[] ~ normal(t_sd/t_mean, 0.1);
  k[] ~ beta(0.5, 0.5);
  q0[] ~ beta(0.5, 0.5);
  z ~ binomial(s, q_fun(N, N_disc, q0, m_hat_p, m_hat, inv_lambda_hat, a));
  alpha0[] ~ exponential(1e-3);
  u ~ binomial(n, P);
  u ~ beta_binomial(n, alpha0[a] .* P, alpha0[a] .* (1.0 - P));
}

generated quantities{
  int<lower = 0> n_tilde[N];
  n_tilde = rep_array(10000, N);
  int<lower = 0> u_tilde[N];
  u_tilde = beta_binomial_rng(n_tilde, alpha0[a] .* P, alpha0[a] .* (1.0 - P));
}
