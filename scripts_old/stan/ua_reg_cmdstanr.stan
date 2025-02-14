// use_acc_reg_all_cc.stan
// Usage/access regression model for all countries
// Originally from log_growth_MDC_decay_beta_bin_MDC21_q.stan

functions {
  
  //discretisation function
  vector disc(real LB, real UB, int steps) {
    vector[steps] x;
    real interval = (UB - LB) / steps;
    x[1] = LB + interval / 2.0;
    for (i in 2:steps) {
      x[i] = x[i-1] + interval;
    }
    return x;
  }
  
  //count number of negative values
  int num_negative(vector x) {
    int N = rows(x);
    //int y[N];
    array[N] int y;
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
  
  //calculate proportion of campaign nets m months after a mass campaign
  vector q_fun(int N, int N_disc, vector q0, matrix m_p, matrix m, vector inv_lambda, array[] int a) {
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
  array[N_a] int<lower = 0, upper = N_a> cc; //country id by admin
  vector<lower = 0>[N_a] mu_n;      //mean net life
  vector<lower = 0>[N_a] sigma_n;   //sd net life
  int<lower = 0> N;                 //number of trials/time points over admins
  array[N] int<lower = 1, upper = N_c> c; //countries
  array[N] int<lower = 1, upper = N_a> a; //admin one units
  vector<lower = 0>[N] t;           //time (Calendar Month Code)
  array[N] int<lower = 0> u;              //number used (or with access) for each trial
  array[N] int<lower = 0> n;              //total for each trial
  array[N] int<lower = 0> s;              //number with source of net recorded
  array[N] int<lower = 0> z;              //number with campaign as source
  int<lower = 0> max_rho;           //maximum number of mass campaign rounds
  matrix[N_a, max_rho] r_meas;      //campaign round timing central estimates
  matrix[N_a, max_rho] r_tau;       //campaian round timing uncertainty
  array[N] int<lower = 0> rho;            //most recent campaign round index
  real<lower = 0> disc_rnge;        //discretisation range
  int<lower = 1> N_disc;            //number of discretisation values
  real<lower = 0> max_m;            //maximum number of months since a campaign
  real no_round;                    //no previous campaign round indicator
  int<lower = 1> N_bb;              //number of beta-binomial samples
  
}

transformed data {
  
  //calculate half discretisation range
  real hlf_rnge = disc_rnge / 2.0;
  
  //normalise time values
  real<lower = 0> t_mean = mean(t);
  real<lower = 0> t_sd = sd(t);
  vector[N] t_hat = (t - t_mean) / t_sd;
  
  //normalise months since campaign and prepare discretised grids for uncertainty
  real max_m_hat = max_m / t_sd;
  matrix[N_disc, N] m_hat;            //normalised m over discretised range
  matrix[N_disc, N] m_hat_p;          //associated probabilities
  
  //normalise timings of mass campaigns
  matrix[N_a, max_rho] r_meas_hat = (r_meas - t_mean) / t_sd;
  matrix[N_a, max_rho] r_tau_hat = r_tau / t_sd;
  //matrix[max_rho, N_disc] r_hat[N_a];
  array[N_a] matrix[max_rho, N_disc] r_hat; //array of matrices for round timing
                                      //x=admin, y=round, z=time discretisation
  matrix[N_a, max_rho] r_hat_LB = r_meas_hat - r_tau_hat * hlf_rnge;
  matrix[N_a, max_rho] r_hat_UB = r_meas_hat + r_tau_hat * hlf_rnge;
  
  //create discretised grid for round timings around central estimates
  for (i in 1:N_a) {
    for (j in 1:max_rho) {
      if (r_meas[i,j] < 0) {
        r_hat[i,j,] = rep_row_vector(no_round, N_disc);
      } else {
        r_hat[i,j,] = to_row_vector(disc(r_hat_LB[i,j], r_hat_UB[i,j], N_disc));
      }
    }
  }
  
  //generate a discretised grid of standard normal densities
  vector[N_disc] SN_z = disc(-hlf_rnge, hlf_rnge, N_disc);
  vector[N_disc] SN_den;
  for (i in 1:N_disc) {
    SN_den[i] = exp(std_normal_lpdf(SN_z[i]));
  }
  real norm_cnst = sum(SN_den);
  vector[N_disc] prpr_SN_den = SN_den / norm_cnst; //proper std norm density
  
  //calculate the probabilities for the number of months since the last campaign
  //the probabilities for the closest round are considered first
  //the probabiliites for the penultimate round to this are then considered
  //loop over all trials/time points over admin levels
  for (j in 1:N) {
    int j1;   //closest round
    int j2;   //penultimate round to closest
    if (rho[j] == max_rho) {
      //if the preceding round is the maximum permitted
      j1 = rho[j];
      j2 = rho[j] - 1;
    } else if (r_meas[a[j], rho[j]+1] < 0) {
      //if preceding round is flagged as no round
      j1 = rho[j];
      j2 = rho[j] - 1;
    } else if ((t[j] - r_meas[a[j],rho[j]]) < (r_meas[a[j],rho[j]+1] - t[j])) {
      //if preceding round is closer than succeding round
      j1 = rho[j];
      j2 = rho[j] - 1;
    } else {
      //if succeding round is closer than preceding round
      j1 = rho[j] + 1;
      j2 = rho[j];
    }
    //discretised grid of months since closest campaign
    m_hat[,j] = t_hat[j] - to_vector(r_hat[a[j],j1,]);
    m_hat_p[,j] = prpr_SN_den;   //set probabilites to proper discrete std normal
    //identify the number of negative months since closest campaign
    //(i.e. probabilities for the closest campaign to be succeding)
    int N_neg = num_negative(m_hat[,j]);
    int N_pos = N_disc - N_neg;
    int N_posp1 = N_pos + 1;
    real pos_p_sum = 0.0;
    if (N_pos > 0) {
      pos_p_sum = sum(m_hat_p[1:N_pos,j]);
    }
    //if some months were negative
    if (N_neg > 0) {
      //if the penultimate campaign to the closest one exists
      if (j2 > 0) {
        //calculate probabilities for the months since the penultimate campaign
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
        m_hat_p[N_posp1:N_disc,j] = rep_vector(((1.0 - pos_p_sum) / (1.0 * N_neg)), N_neg);
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
  vector[N_a] inv_lambda_hat = inv_lambda / t_sd;
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
  
  //strongly-informative prior
  inv_lambda ~ normal(mu_n, sigma_n);
  
  //weakly-informative priors
  beta_0[] ~ normal(0, 0.1);
  beta_t[] ~ normal(t_sd/t_mean, 0.1);
  alpha0[] ~ exponential(1e-3);
  
  //non-informative (Jeffrey's) priors
  k[] ~ beta(0.5, 0.5);
  q0[] ~ beta(0.5, 0.5);
  
  //joint likelihood
  z ~ binomial(s, q_fun(N, N_disc, q0, m_hat_p, m_hat, inv_lambda_hat, a));
  u ~ binomial(n, P);
  u ~ beta_binomial(n, alpha0[a] .* P, alpha0[a] .* (1.0 - P));
  
}

generated quantities{
  array[N] int<lower = 0> n_tilde;
  //n_tilde = rep_array(10000, N);
  n_tilde = rep_array(N_bb, N);
  array[N] int<lower = 0> u_tilde;
  u_tilde = beta_binomial_rng(n_tilde, alpha0[a] .* P, alpha0[a] .* (1.0 - P));
}