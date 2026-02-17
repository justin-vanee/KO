data {
  int<lower=1> N;                 // observations
  int<lower=0> N_p;               // observations (predicted)
  int<lower=1> P;                 // fixed-effect columns
  int<lower=1> G;                 // number of genotypes
  int<lower=1> T;                 // number of traits
  int<lower=0> T_a;               // number of traits observed in all environments 
  int<lower=1> E;                 // number of environments
  int<lower=1> q;                 // number of factors 

  // Observations 
  array[N] real y;     

  // indices (obs)
  array[N] int<lower=1, upper=G> g_idx;
  array[N] int<lower=1, upper=T> t_idx;
  array[N] int<lower=1, upper=E> e_idx;
  
  // indices (pred)
  array[N_p] int<lower=1, upper=G> g_idx_p;
  array[N_p] int<lower=1, upper=T> t_idx_p;
  array[N_p] int<lower=1, upper=E> e_idx_p;

  // Design Matrices 
  matrix[E, P] X;         
}

parameters {
  // Loading matrices for each environment
  matrix[T, q] W;  

  // Fixed effects
  matrix[q, P] eta;

  // Error variance trait x environment (latent)
  matrix[T_a, E] log_sigma2; // different priors for traits observed in all versus some environments         
  matrix<lower=0>[T-T_a, E] sigma2_s; 
  vector<lower=0>[T_a] sigma2_sigma2;
  vector[T_a] mu_sigma2;
}

transformed parameters {
  // different priors for traits observed in all versus some environments   
  matrix<lower=0>[T, E] sigma2;
  for(t in 1:T){
    if(t <= T_a){
      // hierichical prior for traits observed in many environments 
      sigma2[t,] = exp(log_sigma2[t,]); 
    } else {
      // half-cauchy(0,1) prior for traits observed in 1 to a few environments 
      sigma2[t,] = sigma2_s[t-T_a,];
    }
  }
    
  // Environmental Associations 
  matrix[P, T] beta;         
  for (t in 1:T) {
    for (p in 1:P) {
      beta[p, t] = dot_product(eta[, p], W[t,]);
    }
  }
}

model {
  // Loading matrix (unidentifiable)
  to_vector(W) ~ normal(0, 1);
  
  // Regression coefficients (unidentifiable)
  to_vector(eta) ~ normal(0, 1);

  // Error and hierarchical variances 
  to_vector(sigma2_s) ~ gamma(1, 1);
  sigma2_sigma2 ~ gamma(1, 1);
  mu_sigma2 ~ normal(0, 1);
  for (t in 1:T_a)
    log_sigma2[t,] ~ normal(mu_sigma2[t], sqrt(sigma2_sigma2[t]));

  // Likelihood (observations are independent normals conditional on u)
  for (i in 1:N) {
    int t = t_idx[i];
    int e = e_idx[i];
    int g = g_idx[i];
    y[i] ~ normal(X[e,] * beta[, t], sqrt(sigma2[t, e]));
  }
}

generated quantities {
  // Posterior predictive draws (for requested N_p locations)
  vector[N_p] y_pred;
  for (i in 1:N_p) {
    int g = g_idx_p[i];
    int t = t_idx_p[i];
    int e = e_idx_p[i];
    y_pred[i] = normal_rng(X[e,] * beta[, t], sqrt(sigma2[t, e]));
  }
}