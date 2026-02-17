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

  // Genomic relationship matrix (must be PSD, aligned to genotype levels)
  matrix[G, G] K;
}

transformed data {
  // stabilized Cholesky of K for non-centered parameterization
  matrix[G, G] K_jitter = K;
  for (i in 1:G) K_jitter[i, i] += 1e-8;
  matrix[G, G] L_K = cholesky_decompose(K_jitter);
}

parameters {
  // Loadings / rotations for each trait (lower-dimensional representation)
  matrix[T, q] W;  

  // Latent-variable regression coefficients (raw)
  matrix[q, P] eta;
  
  // Standard normals for genomic random effects (non-centered)
  array[E] matrix[G, T] z_u;  
   
  // principal component genetic covariances (variances)
  matrix<lower=0>[q, E] zeta2;

  // Error variance trait x environment (latent)
  matrix[T_a, E] log_sigma2; // hierarchical for traits observed in all envs
  matrix<lower=0>[T-T_a, E] sigma2_s; 
  vector<lower=0>[T_a] sigma2_sigma2;
  vector[T_a] mu_sigma2;
}

transformed parameters {
  // different priors for traits observed in all versus some environments   
  matrix<lower=0>[T, E] sigma2;
  for (t in 1:T) {
    if (t <= T_a) {
      sigma2[t,] = exp(log_sigma2[t,]); 
    } else {
      sigma2[t,] = sigma2_s[t - T_a,];
    }
  }

  // Implied genetic variances (per trait × env)
  matrix[T, E] var_g;
  for (t in 1:T) {
    for (e in 1:E) {
      var_g[t, e] = dot_product(zeta2[, e], square(W[t,])); // zeta2 is variance per PC
    }
  }

  // Environmental Associations (map from latent PC space to observed predictors)
  matrix[P, T] beta;
  for (t in 1:T) {
    for (p in 1:P) {
      beta[p, t] = dot_product(eta[, p], W[t,]);
    }
  }

  // Genomic random effects (non-centered -> u)
  array[E] matrix[G, T] u;
  for (e in 1:E) {
    for (t in 1:T) {
      // fixed effect contribution + genetic random effect (non-centered)
      // note: L_K * z_u[e][, t] is G-vector, scaled by sqrt(var_g[t,e])
      u[e][, t] = rep_vector( dot_product(X[e,], beta[, t]), G ) + (L_K * z_u[e][, t]) * sqrt(var_g[t, e]);
    }
  }
}

model {
  // Loading matrix (unidentifiable)
  to_vector(W) ~ normal(0, 1);
  
  // Regression coefficients (unidentifiable)
  to_vector(eta) ~ normal(0, 1);

  // Standard normal random variables 
  for (e in 1:E)
    to_vector(z_u[e]) ~ normal(0, 1);         

  // Latent variable variances (unidentifiable)
  to_vector(zeta2) ~ gamma(1.0 / q, 1);
  
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
    y[i] ~ normal(u[e][g, t], sqrt(sigma2[t, e]));
  }
}

generated quantities {
  // Heritability (per trait × env)
  matrix[T, E] H2;
  for (t in 1:T) {
    for (e in 1:E) {
      H2[t, e] = var_g[t, e] / (var_g[t, e] + sigma2[t, e]);
    }
  }

  // Posterior predictive draws (for requested N_p locations)
  vector[N_p] y_pred;
  for (i in 1:N_p) {
    int g = g_idx_p[i];
    int t = t_idx_p[i];
    int e = e_idx_p[i];
    y_pred[i] = normal_rng(u[e][g, t], sqrt(sigma2[t, e]));
  }
}

