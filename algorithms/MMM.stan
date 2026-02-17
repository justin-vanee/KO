functions {
  // Kronecker product for Cholesky factors (lower-triangular)
  matrix kron_chol(matrix L_A, matrix L_B) {
    int n_A = rows(L_A);
    int n_B = rows(L_B);
    int n_Z = n_A * n_B;
    matrix[n_Z, n_Z] Z;
    
    Z = rep_matrix(0, n_Z, n_Z);
    
    for (i in 1:n_A) {
      for (j in 1:i) {  // lower triangle only
        for (p in 1:n_B) {
          for (q in 1:p) {
            int row = (i - 1) * n_B + p;
            int col = (j - 1) * n_B + q;
            Z[row, col] = L_A[i,j] * L_B[p,q];
          }
        }
      }
    }
    return Z;
  }
  
}

data {
  
  int<lower=1> N;       // observations
  int<lower=1> N_p;     // predicted observations
  int<lower=1> P;       // fixed-effect columns
  int<lower=1> G;       // genotypes
  int<lower=1> T;       // traits
  int<lower=0> T_a;     // traits observed in all environments 
  int<lower=1> E;       // environments

  array[N] real y;       // response
  matrix[E, P] X;        // fixed design matrix

  array[N] int<lower=1, upper=G> g_idx;
  array[N] int<lower=1, upper=T> t_idx;
  array[N] int<lower=1, upper=E> e_idx;

  array[N_p] int<lower=1, upper=G> g_idx_p;
  array[N_p] int<lower=1, upper=T> t_idx_p;
  array[N_p] int<lower=1, upper=E> e_idx_p;

  matrix[G, G] K;        // genomic relationship matrix
  
}

transformed data { 
    matrix[G, G] L_K = cholesky_decompose(K);
}

parameters {
  // Environmental Associations 
  matrix[P, T] beta;

  cholesky_factor_corr[T] L_Omega;   // trait correlation
  matrix<lower=0>[T, E] zeta2;       // env-specific scaling

  // Standard normal latent variables for each environment
  array[E] vector[G*T] z_u;

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
  
  // scaled covariance matrices for each environment 
  array[E] matrix[T, T] L_Sigma;
  array[E] matrix[T, T] Sigma_trait;
  for(e in 1:E){
    L_Sigma[e] = diag_pre_multiply(sqrt(zeta2[,e]), L_Omega);
    Sigma_trait[e] = multiply_lower_tri_self_transpose(L_Sigma[e]);
  }
  
  // Latent (TxG) random effects per environment
  array[E] matrix[G, T] u;
  for (e in 1:E) {
    matrix[G*T, G*T] L_kron = kron_chol(L_Sigma[e], L_K);
    u[e] = to_matrix(L_kron * z_u[e], G, T);
  }
}

model {
  // Priors
  to_vector(beta) ~ normal(0, 1);
  to_vector(zeta2) ~ gamma(1, 1);
  L_Omega ~ lkj_corr_cholesky(2);

  // Standard normals for Kronecker transformation
  for (e in 1:E) {
    z_u[e] ~ normal(0, 1);
  }
  
  // Error variance for traits only observed in a few environments 
  to_vector(sigma2_s) ~ gamma(1, 1);
  
  // Hierichical variance (variance of error variances)
  sigma2_sigma2 ~ gamma(1, 1);
  
  // Weakly informative priors for zeta2 and sigma2 
  mu_sigma2 ~ normal(0, 1); // results in a roughly uniform prior heritability  
  for(t in 1:T_a)
    log_sigma2[t,] ~ normal(mu_sigma2[t], sqrt(sigma2_sigma2[t]));

  // Likelihood
  for (n in 1:N) {
    int t = t_idx[n];
    int e = e_idx[n];
    int g = g_idx[n];
    real mu = dot_product(X[e,], beta[, t]) + u[e][g, t];
    y[n] ~ normal(mu, sqrt(sigma2[t,e]));
  }
  
}

generated quantities {
  // Heritability
  matrix[T, E] H2;
  for (t in 1:T) {
    for (e in 1:E) {
      real var_g = zeta2[t, e];
      real var_e = sigma2[t, e];
      H2[t,e] = var_g / (var_g + var_e);
    }
  }

  // Posterior predictive
  vector[N_p] y_pred;
  for (i in 1:N_p) {
    int t = t_idx_p[i];
    int e = e_idx_p[i];
    int g = g_idx_p[i];
    real mu = dot_product(X[e,], beta[, t]) + u[e][g, t];
    y_pred[i] = normal_rng(mu, sqrt(sigma2[t,e]));
  }
  
}
