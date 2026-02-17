###
### Helper functions 
###

simulate_h2_prior <- function(n = 10000, q = 4, shape = 1, rate = 1, hierarchical = TRUE) {
  # Genetic variance as sum of q gamma-chi-squared mixtures
  g <- numeric(n)
  for (i in seq_len(q)) {
    g <- g + rgamma(n, shape = 1 / q, rate = rate) * rnorm(n)^2
  }
  
  # Environmental variance
  if (hierarchical) {
    mean <- rnorm(n, sd = 1)
    var <- rgamma(n, shape, rate)
    e <- exp(rnorm(n, mean, sqrt(var)))
  } else {
    e <- rgamma(n, shape = shape, rate = rate)
  }
  
  # Heritability
  h2 <- g / (g + e)
  return(h2)
}

### This function produces a binary matrix based on which traits where observed in which environments 
get_env_trait_matrix <- function(df) {
  # Get unique environments and traits
  envs   <- sort(unique(df$e))
  traits <- sort(unique(df$t))
  
  # Initialize E x T matrix with zeros
  mat <- matrix(0, nrow = length(envs), ncol = length(traits),
                dimnames = list(envs, traits))
  
  # Loop through each environment and trait
  for (e in envs) {
    for (t in traits) {
      vals <- df$value[df$e == e & df$t == t]
      if (any(!is.na(vals) & !is.nan(vals))) {
        mat[as.character(e), as.character(t)] <- 1
      }
    }
  }
  
  return(t(mat))
}

rotation_mat <- function(theta){
  return(matrix(
    c(cos(theta), -sin(theta), sin(theta), cos(theta)),
    2,
    2
  ))
}

tmp_func <- function(theta, W){
  return((W[1,]%*%rotation_mat(theta))[1])
}

heritability_func <- function(zeta2, sigma2){
  return(zeta2 / (sigma2 + zeta2))
}

simulate_data <- 
  function(
    G = 20, T = 6, E = 5, P = 1, # Dimensions 
    target_var = 0.9, percent_remove = 0.0, trait_cor = 0.95, # Hyperparameters  
    dir # Directory for kinship matrix 
    ) {

  # --- Parameters ---
  Sigma <- matrix(trait_cor, T, T); diag(Sigma) <- 1
  zeta2 <- matrix(rgamma(T * E, 1, 1), T, E)
  sigma2 <- matrix(rgamma(T * E, 1, 1), T, E)
  beta <- matrix(rnorm(T * (P + 1)), P + 1, T)
  
  # --- Kinship matrix ---
  # Load the example DO mouse data from the qtl2 package
  dirpath <- system.file("extdata", "iron.zip", package = "qtl2")
  cross <- read_cross2(dirpath)
  
  # Calculate genotype probabilities with an error probability of 0.002
  probs <- calc_genoprob(cross, map = cross$gmap, error_prob = 0.002)
  
  # Calculate the kinship matrix using the overall method
  kinship <- calc_kinship(probs)
  idx <- sample(1:ncol(kinship), G, replace = FALSE)
  K <- cov2cor(kinship[idx, idx])
  
  # --- Covariates ---
  X <- cbind(1, matrix(rnorm(P * E), E, P))
  
  # --- Simulate trait data ---
  Y <- array(NA, dim = c(G, T, E))
  for (i in 1:E) {
    Sigma_trait <- diag(sqrt(zeta2[, i])) %*% Sigma %*% diag(sqrt(zeta2[, i]))
    u <- rmvn(1, rep(0, G * T), Sigma_trait %x% K) %>%
      matrix(G, T)
    for (j in 1:T) {
      mu <- X[i, ] %*% beta[, j]
      Y[, j, i] <- rnorm(G, rep(mu, G) + u[, j], sqrt(sigma2[j, i]))
    }
  }
  for (j in 1:T) Y[, j, ] <- scale(c(Y[, j, ])) %>% matrix(G, E)
  Y_long <- do.call(rbind, map(1:E, ~ Y[, , .x]))

  # --- Heritability ---
  H2_true <- heritability_func(zeta2, sigma2)
  
  # --- Prepare for Stan ---
  data_stan <- Y_long %>%
    as.data.frame() %>%
    mutate(g = rep(1:G, E),
           e = rep(1:E, each = G)) %>%
    pivot_longer(-c(g, e), names_to = "t", values_to = "value") %>%
    mutate(t = as.integer(as.factor(t)),
           truth = value) %>%
    arrange(e, t)
  
  # --- Remove observations for out-of-sample ---
  n_remove <- floor(percent_remove * G * E)
  GxE <- expand.grid(g = 1:G, e = 1:E)
  remove_combos <- GxE %>%
    slice_sample(n = n_remove)
  # Set value = NA for all rows matching those combinations
  data_stan$value <- ifelse(
    paste(data_stan$g, data_stan$e) %in% paste(remove_combos$g, remove_combos$e),
    NA,
    data_stan$value
  )
  y_true <-
    data_stan %>%
    filter(is.na(value)) %>%
    pluck("truth")
  NAs <- which(is.na(data_stan$value))
  
  # --- Principal components ---
  Y_cor <- cor(Y_long, use = "pairwise.complete.obs")
  vals <- eigen(Y_cor)$values
  q <- which(cumsum(vals / sum(vals)) >= target_var)[1]
  
  if(length(NAs) > 0){
    
    # --- Return data list ---
    list(
      N = prod(G, E, T) - length(NAs),
      N_p = length(NAs),
      P = ncol(X),
      G = G,
      E = E,
      T = T,
      T_a = T, 
      q = q,
      g_idx = data_stan$g[-NAs],
      e_idx = data_stan$e[-NAs],
      t_idx = data_stan$t[-NAs],
      y = data_stan$value[-NAs],
      Y = map(1:E, ~ Y[, , .x]),
      X = X,
      K = K,
      g_idx_p = data_stan$g[NAs],
      e_idx_p = data_stan$e[NAs],
      t_idx_p = data_stan$t[NAs],
      sigma2 = sigma2,
      zeta2 = zeta2,
      beta_true = beta,
      H2_true = H2_true,
      y_true = y_true
    )
    
  } else {
    
    # --- Return data list ---
    list(
      N = prod(G, E, T) - length(NAs),
      N_p = length(NAs),
      P = ncol(X),
      G = G,
      E = E,
      T = T,
      T_a = T, 
      q = q,
      g_idx = data_stan$g,
      e_idx = data_stan$e,
      t_idx = data_stan$t,
      y = data_stan$value,
      Y = map(1:E, ~ Y[, , .x]),
      X = X,
      K = K,
      g_idx_p = data_stan$g[NAs],
      e_idx_p = data_stan$e[NAs],
      t_idx_p = data_stan$t[NAs],
      sigma2 = sigma2,
      zeta2 = zeta2,
      beta_true = beta,
      H2_true = H2_true,
      y_true = y_true
    )
    
  }
}

fit_model <- function(
    data,
    dir, 
    chains, iter_warmup, iter_sampling, refresh,
    model = c("KO", "MMM", "UO", "CO", "CCO"),
    pathfinder = TRUE
) {
  model <- match.arg(model)
  
  # Construct path using paste0()
  stan_file <- file.path(dir, "algorithms", paste0(model, ".stan"))
  mod <- cmdstan_model(stan_file)
  
  # Default pathfinder init functions
  init_fun <- NULL
  if (model == "KO") {
    init_fun <- function() list(W = matrix(rnorm(data$T * data$q, 0, 0.1), data$T, data$q))
  } else if (model == "MMM") {
    init_fun <- function() list(
      beta = matrix(0, nrow = data$P, ncol = data$T),
      L_Omega = diag(1, data$T),
      tau_trait = rep(1, data$T),
      sigma_g_env = rep(1, data$E),
      u = array(0, dim = c(data$E, data$G, data$T)),
      u_mat = array(0, dim = c(data$G * data$T, data$E)),
      sigma = matrix(1, nrow = data$T, ncol = data$E),
      sigma_sigma = rep(0.1, data$T),
      mu_sigma = rep(0, data$T)
    )
  }
  
  # --- Fitting logic ---
  if (pathfinder) {
    pathfinder_fit <- mod$pathfinder(
      data = data,
      num_paths = 1,
      init = init_fun
    )
    
    fit <- mod$sample(
      data = data,
      chains = chains,
      parallel_chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      init = pathfinder_fit$draws(format = "list"),
      refresh = refresh
    )
  } else {
    fit <- mod$sample(
      data = data,
      chains = chains,
      parallel_chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      refresh = refresh
    )
  }
  
  return(fit)
}

extract_performance <- function(fit, tib, data, H2 = TRUE, beta = TRUE) {
  
  tib$time <- fit$time()$total
  
  if (H2) {
    
    ### Extract Performance (Heritability)
    heritability_results <-
      fit$summary("H2") %>%
      mutate(truth = c(data$H2_true),
             coverage = truth > q5 & truth < q95,
             bias = abs((mean - truth)))
    
    tib$coverage <- mean(heritability_results$coverage)
    tib$bias <- mean(heritability_results$bias)
    tib$sd <- mean(heritability_results$sd)
    tib$cor <- cor(heritability_results$mean, heritability_results$truth)
    tib$ess <- mean(heritability_results$ess_bulk)
    
  } else {
    
    tib$coverage <- NA
    tib$bias <- NA
    tib$sd <- NA
    tib$cor <- NA
    tib$ess <- NA
    
  }
  
  if (beta) {
    
    ### Extract Performance (environmental associations)
    beta_results <-
      fit$summary("beta") %>%
      mutate(truth = c(data$beta_true),
             intercept = rep(c(TRUE, rep(FALSE, data$P-1)), data$T),
             coverage = truth > q5 & truth < q95,
             bias = abs((mean - truth))) %>%
      filter(!intercept)
    
    tib$coverage_b <- mean(beta_results$coverage)
    tib$bias_b <- mean(beta_results$bias)
    tib$sd_b <- mean(beta_results$sd)
    tib$cor_b <- cor(beta_results$mean, beta_results$truth)
    tib$ess_b <- mean(beta_results$ess_bulk)
    
  } else {
    
    tib$coverage_b <- NA
    tib$bias_b <- NA
    tib$sd_b <- NA
    tib$cor_b <- NA
    tib$ess_b <- NA
    
  }
  
  ### Extract Performance (Prediction)
  prediction_results <-
    fit$summary("y_pred") %>%
    mutate(truth = data$y_true,
           coverage = truth > q5 & truth < q95,
           bias = abs((mean - truth)))
  
  tib$coverage_Y <- mean(prediction_results$coverage)
  tib$bias_Y <- mean(prediction_results$bias)
  tib$sd_Y <- mean(prediction_results$sd)
  tib$cor_Y <- cor(prediction_results$mean, prediction_results$truth)
  tib$ess_Y <- mean(prediction_results$ess_bulk)
  
  return(tib)
}

extract_performance_CV <- function(fit, mod, test) {
  
  prediction_results <-
    fit$summary("y_pred") %>%
    cbind(test) %>%
    mutate(coverage = value > q5 & value < q95,
           bias = abs(mean - value),
           time = fit$time()$total,
           model = mod)
  
  return(prediction_results)
}


simulation_study <- function(data, models = c("KO", "MMM", "CO", "UO", "CCO"),
                             chains = 1, parallel_chains = 1,
                             iter_warmup = 200, iter_sampling = 1000,
                             refresh = 100) {
  
  # Initialize empty tibbles for each model
  metrics_list <- lapply(models, function(m) tibble(model = m))
  names(metrics_list) <- models

  
  # Loop over models
  for (model_name in models) {
    fit <- fit_model(
      data = data,
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      refresh = refresh,
      model = model_name
    )
    
    # Special case for CO
    if (model_name == "CO") {
      metrics_list[[model_name]] <- extract_performance(
        fit = fit,
        tib = metrics_list[[model_name]],
        beta = TRUE,
        H2 = FALSE,
        data = data
      )
    } else if(model_name == "UO"){
      metrics_list[[model_name]] <- extract_performance(
        fit = fit,
        tib = metrics_list[[model_name]],
        beta = FALSE,
        H2 = TRUE,
        data = data
      )
    } else {
      metrics_list[[model_name]] <- extract_performance(
        fit = fit,
        tib = metrics_list[[model_name]],
        data = data
      )
    }
    
    # Add dataset dimensions
    metrics_list[[model_name]] <- metrics_list[[model_name]] %>%
      mutate(T = data$T, G = data$G, E = data$E)
  }
  
  # Combine all metrics into a single tibble
  metrics_tib <- bind_rows(metrics_list)
  return(metrics_tib)
}


