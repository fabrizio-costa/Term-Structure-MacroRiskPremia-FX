source("code/matrix.notation.R")
source("code/posterior.R")

# ---------------------------- (1) Gibb Sampler  ----------------------------

gibbs_sampler1 <- function(g_t, v_t, S_bar, num_iter = 1000, burn_in = 200, seed = 42) {
  set.seed(seed)
  t <- length(g_t)
  K <- ncol(v_t)  # v_t is T x K matrix
  
  # initialize storage matrices for samples
  sigma_wg2_samples <- numeric(num_iter)
  rho_g_samples <- matrix(0, nrow = num_iter, ncol = S_bar + 2)
  eta_g_samples <- matrix(0, nrow = num_iter, ncol = K)
  
  # random initial values
  rho_g <- rep(0.1, S_bar + 2)             # initial rho_g
  eta_g <- rep(1 / sqrt(K), K)             # initial eta_g (unit norm)
  
  # mean of v_t
  mu_v <- colMeans(v_t)  # K-dimensional vector
  
  # Gibbs sampling loop
  for (iter in 1:num_iter) {
    # Step 1: Sample sigma_wg²
    sigma_wg2 <- get.sigma.post(g_t, rho_g, eta_g, v_t, S_bar)
    
    # Step 2: Construct V_rho
    V_rho <- get.V_rho(v_t, eta_g, mu_v, S_bar)
    
    # Step 3: Sample rho_g (adjusted version)
    rho_g <- get.rho_g.adjusted(G = get.G(g_t, S_bar), sigma_wg2, V_rho, S_bar)
    
    # Step 4: Construct V_eta
    V_eta <- get.V_eta(v_t, rho_g, mu_v, S_bar)
    
    # Step 5: Sample eta_g (already normalized inside)
    eta_g <- get.eta_g(get.G_bar(g_t, S_bar), sigma_wg2, rho_g, V_eta, shrink_eta = FALSE)
    
    # Store samples
    sigma_wg2_samples[iter] <- sigma_wg2
    rho_g_samples[iter, ] <- as.numeric(rho_g)
    eta_g_samples[iter, ] <- as.numeric(eta_g)
    
    if (iter %% 100 == 0) cat("Iteration", iter, "\n")
  }
  
  # Remove burn-in
  if (burn_in > 0 && burn_in < num_iter) {
    sigma_wg2_samples <- sigma_wg2_samples[(burn_in + 1):num_iter]
    rho_g_samples <- rho_g_samples[(burn_in + 1):num_iter, ]
    eta_g_samples <- eta_g_samples[(burn_in + 1):num_iter, ]
  }
  
  return(list(sigma_wg2 = sigma_wg2_samples,
              rho_g = rho_g_samples,
              eta_g = eta_g_samples))
}


# ---------------------------- (2) Gibb Sampler  ----------------------------

gibbs_sampler2 <- function(R, v_mat, num_iter = 1000, burnin = 200, seed = 42) {
  set.seed(seed)
  
  t <- nrow(R)
  N <- ncol(R)
  K <- ncol(v_mat)
  
  # storage for ALL draws
  B_r_samples <- array(NA, dim = c(num_iter, K + 1, N))
  Sigma_wr_samples <- array(NA, dim = c(num_iter, N, N))
  
  # initialization
  mu_r <- colMeans(R)
  beta_v <- matrix(0, nrow = K, ncol = N)
  
  for (iter in 1:num_iter) {
    Sigma_wr <- get.Sigma_wr(R, v_mat, mu_r, beta_v)
    B_r <- get.B_r.post(R, v_mat, Sigma_wr)
    
    # # DEBUG dimensions
    # if (iter == 1) {
    #   cat("Iteration 1:\n")
    #   cat("- B_r dim: "); print(dim(B_r))
    #   cat("- Sigma_wr dim: "); print(dim(Sigma_wr))
    # }
    
    mu_r <- B_r[1, ]
    beta_v <- B_r[-1, , drop = FALSE]
    
    B_r_samples[iter, , ] <- B_r
    Sigma_wr_samples[iter, , ] <- Sigma_wr
    
    if (iter %% 100 == 0) cat("Iteration", iter, "complete\n")
  }
  
  # cat("Final B_r_samples dim: "); print(dim(B_r_samples))
  # cat("Final Sigma_wr_samples dim: "); print(dim(Sigma_wr_samples))
  # cat("Burn-in:", burnin, " Num_iter:", num_iter, "\n")
  # 
  # --- Burn-in trimming --- (robusto per num_iter = 1)
  if (burnin >= num_iter) {
    warning("⚠️ burn-in >= num_iter: skipping trimming, returning all draws.")
    keep_idx <- 1:num_iter
  } else {
    keep_idx <- (burnin + 1):num_iter
  }
  
  B_r_samples <- B_r_samples[keep_idx, , , drop = FALSE]
  Sigma_wr_samples <- Sigma_wr_samples[keep_idx, , , drop = FALSE]
  
  # cat("Post burn-in B_r_samples dim: "); print(dim(B_r_samples))
  # cat("Post burn-in Sigma_wr_samples dim: "); print(dim(Sigma_wr_samples))
  # 
  # Also return separated mu_r and beta_v
  mu_r_samples <- B_r_samples[, 1, , drop = FALSE]      # draws x N
  beta_v_samples <- B_r_samples[, -1, , drop = FALSE]   # draws x K x N
  
  # cat("mu_r_samples dim: "); print(dim(mu_r_samples))
  # cat("beta_v_samples dim: "); print(dim(beta_v_samples))
  
  return(list(
    B_r = B_r_samples,
    mu_r = mu_r_samples,
    beta_v = beta_v_samples,
    Sigma_wr = Sigma_wr_samples
  ))
}

gibbs_sampler2 <- function(R, v_mat, num_iter = 1000, burnin = 200, seed = 42) {
  set.seed(seed)
  
  t <- nrow(R)
  N <- ncol(R)
  K <- ncol(v_mat)
  
  # storage for ALL draws
  B_r_samples <- array(NA, dim = c(num_iter, K + 1, N))
  Sigma_wr_samples <- array(NA, dim = c(num_iter, N, N))
  
  # initialization
  mu_r <- colMeans(R)
  beta_v <- matrix(0, nrow = K, ncol = N)
  
  for (iter in 1:num_iter) {
    Sigma_wr <- get.Sigma_wr(R, v_mat, mu_r, beta_v)
    B_r <- get.B_r.post(R, v_mat, Sigma_wr)
    
    mu_r <- B_r[1, ]
    beta_v <- B_r[-1, , drop = FALSE]
    
    # --- Enforce positive sign convention on each factor ---
    for (k in 1:K) {
      mean_sign <- sign(mean(beta_v[k, ]))
      if (mean_sign == -1) {
        beta_v[k, ] <- -beta_v[k, ]
        v_mat[, k] <- -v_mat[, k]  # flip latent factor accordingly
      }
    }
    
    # Update B_r with sign-adjusted beta_v
    B_r[-1, ] <- beta_v
    
    B_r_samples[iter, , ] <- B_r
    Sigma_wr_samples[iter, , ] <- Sigma_wr
    
    if (iter %% 100 == 0) cat("Iteration", iter, "complete\n")
  }
  
  # --- Burn-in trimming ---
  keep_idx <- if (burnin >= num_iter) 1:num_iter else (burnin + 1):num_iter
  
  B_r_samples <- B_r_samples[keep_idx, , , drop = FALSE]
  Sigma_wr_samples <- Sigma_wr_samples[keep_idx, , , drop = FALSE]
  
  mu_r_samples <- B_r_samples[, 1, , drop = FALSE]      # draws x N
  beta_v_samples <- B_r_samples[, -1, , drop = FALSE]   # draws x K x N
  
  return(list(
    B_r = B_r_samples,
    mu_r = mu_r_samples,
    beta_v = beta_v_samples,
    Sigma_wr = Sigma_wr_samples
  ))
}


# ---------------------------- (3) Gibb Sampler  ----------------------------

gibbs_sampler3 <- function(R, mu_r, beta_v, Sigma_wr, shrinkage = FALSE, 
                           num_iter = 1000, burnin = 200, seed = 42) {
  set.seed(seed)
  
  t <- nrow(R)           # number of time points
  N <- ncol(R)           # number of assets
  K <- nrow(beta_v)      # number of latent factors (from beta_v: K x N)
  
  # Initialize latent factors: T x K
  v_t <- matrix(rnorm(t * K), nrow = t, ncol = K)
  
  # Initial values for mu_v and Sigma_v
  mu_v <- rep(0, K)
  Sigma_v <- diag(K)
  
  # Storage for samples
  mu_v_samples <- matrix(NA, nrow = num_iter, ncol = K)
  Sigma_v_samples <- array(NA, dim = c(num_iter, K, K))
  v_t_samples <- array(NA, dim = c(num_iter, t, K))  # each row = time, col = factor
  
  for (iter in 1:num_iter) {
    # Step 1: Sample v_t at each time point
    for (tt in 1:t) {  # avoid overwriting `t`
      v_t[tt, ] <- get.v_t(
        r_t = R[tt, ],
        mu_r = mu_r,
        beta_v = beta_v,
        Sigma_wr = Sigma_wr,
        mu_v = mu_v,
        Sigma_v = Sigma_v  
      )
    }

    # Insert diagnostic lines here — before computing Sigma_v
    # cat(sprintf("[Iter %d] sd(v_t): %s\n", iter, paste(round(apply(v_t, 2, sd), 3), collapse = ", ")))
    # cat(sprintf("[Iter %d] max(abs(mean(v_t) - mu_v)): %.3e\n", iter, max(abs(colMeans(v_t) - mu_v))))

    # Step 2: Sample Sigma_v | v_t
    Sigma_v <- get.Sigma_v(v_t, shrinkage = shrinkage, verbose = TRUE, iter = iter)
    
    # Step 3: Sample mu_v | v_t, Sigma_v
    mu_v <- get.mu_v(v_t, Sigma_v, verbose = TRUE, iter = iter)
    
    # if (iter %% 50 == 0) {
    #   cat("Iteration", iter, "\n")
    #   cat("Sigma_v cond number:", kappa(Sigma_v), "\n")
    #   cat("mu_v:", round(mu_v, 8), "\n")
    #   cat("v_t SDs:", round(apply(v_t, 2, sd), 8), "\n")
    # }
    
    # Store draws
    mu_v_samples[iter, ] <- mu_v
    Sigma_v_samples[iter, , ] <- Sigma_v
    v_t_samples[iter, , ] <- v_t  # now: T x K
    
    if (iter %% 100 == 0) cat("Iteration", iter, "\n")
  }
  
  # Drop burn-in
  keep_idx <- (burnin + 1):num_iter
  
  mu_v_out <- mu_v_samples[keep_idx, , drop = FALSE]
  Sigma_v_out <- Sigma_v_samples[keep_idx, , , drop = FALSE]
  v_t_out <- v_t_samples[keep_idx, , , drop = FALSE]
  
  return(list(
    mu_v = mu_v_out,              # always (num_samples x K)
    Sigma_v = Sigma_v_out,        # always (num_samples x K x K)
    v_t = v_t_out                 # always (num_samples x T x K)
  ))
}


# ---------------------------- (4) Gibb Sampler  Full ----------------------------

gibbs_sampler_full <- function(g_t, R, S_bar, K, num_iter = 1000, burnin = 200, seed = 42) {
  set.seed(seed)
  
  t <- nrow(R)
  N <- ncol(R)
  
  # Initial values
  mu_v <- rep(0, K)
  Sigma_v <- diag(seq(1, 0.5, length.out = K))
  v_t <- MASS::mvrnorm(n = t, mu = mu_v, Sigma = Sigma_v)
  
  eta_g <- rnorm(K); eta_g <- eta_g / sqrt(sum(eta_g^2))
  tau <- 0.2
  rho_g <- rnorm(S_bar + 2, mean = 0, sd = tau)
  sigma_wg2 <- 1 / rgamma(1, shape = 2, rate = 1)
  if (!is.finite(sigma_wg2)) sigma_wg2 <- 1.0
  
  mu_r <- rnorm(N, mean = 0, sd = 0.1)
  beta_v <- matrix(rnorm(K * N, mean = 0, sd = 0.1), nrow = K)
  V_r <- get.V_r(v_t)
  B_r <- rbind(mu_r, beta_v)
  residuals <- R - V_r %*% B_r
  Sigma_wr <- diag(apply(residuals, 2, var) + 1e-4)
  
  # Storage
  samples <- list(
    sigma_wg2 = numeric(num_iter),
    rho_g = matrix(NA, num_iter, S_bar + 2),
    eta_g = matrix(NA, num_iter, K),
    mu_r = matrix(NA, num_iter, N),
    beta_v = array(NA, c(num_iter, K, N)),
    B_r = array(NA, c(num_iter, K + 1, N)),
    Sigma_wr = array(NA, c(num_iter, N, N)),
    mu_v = matrix(NA, num_iter, K),
    Sigma_v = array(NA, c(num_iter, K, K)),
    v_t = array(NA, c(num_iter, t, K)),
    lambda_v = matrix(NA, num_iter, K),
    lambda_g = matrix(NA, num_iter, S_bar + 1)
  )
  
  # Flags: set to TRUE to activate freezing after burn-in
  freeze_eta_g <- FALSE
  freeze_lambda_v <- FALSE
  freeze_rho_g <- FALSE
  
  for (iter in 1:num_iter) {
    
    ## Step 1: g_t block
    step1 <- gibbs_sampler1(g_t, v_t, S_bar = S_bar, num_iter = 1, burn_in = 0)
    sigma_wg2 <- step1$sigma_wg2[1]
    rho_g <- step1$rho_g[1, ]
    eta_g <- step1$eta_g[1, ]
    
    ## Step 2: r_t block
    step2 <- gibbs_sampler2(R, v_t, num_iter = 1, burnin = 0)
    if (any(!is.finite(step2$B_r))) stop(sprintf("Invalid B_r at iter %d", iter))
    
    mu_r <- drop(step2$mu_r[1, , ])
    beta_v <- drop(step2$beta_v[1, , ])
    B_r <- drop(step2$B_r[1, , ])
    Sigma_wr <- drop(step2$Sigma_wr[1, , ])
    
    ## Step 3: v_t block
    step3 <- gibbs_sampler3(R, mu_r, beta_v, Sigma_wr, shrinkage = FALSE, num_iter = 1, burnin = 0)
    mu_v <- step3$mu_v[1, ]
    Sigma_v <- step3$Sigma_v[1, , ]
    v_t <- step3$v_t[1, , ]
    
    ## freeze values after burn-in
    if (iter == burnin + 1) {
      freeze_idx <- floor(burnin * (2/3)) + 1  # start of final third
      if (freeze_eta_g) eta_g_fixed <- colMeans(samples$eta_g[freeze_idx:burnin, , drop = FALSE])
      if (freeze_lambda_v) lambda_v_fixed <- colMeans(samples$lambda_v[freeze_idx:burnin, , drop = FALSE])
      if (freeze_rho_g) rho_g_fixed <- apply(samples$rho_g[freeze_idx:burnin, , drop = FALSE], 2, median)
    }
    
    if (iter > burnin) {
      if (freeze_eta_g) eta_g <- eta_g_fixed
      if (freeze_lambda_v) lambda_v_iter <- lambda_v_fixed else lambda_v_iter <- get.lambda_v(beta_v, mu_r, Sigma_wr)
      if (freeze_rho_g) rho_g <- rho_g_fixed
    } else {
      lambda_v_iter <- get.lambda_v(beta_v, mu_r, Sigma_wr)
    }
    
    ## Store draws
    samples$sigma_wg2[iter] <- sigma_wg2
    samples$rho_g[iter, ] <- rho_g
    samples$eta_g[iter, ] <- eta_g
    samples$mu_r[iter, ] <- mu_r
    samples$beta_v[iter, , ] <- beta_v
    samples$B_r[iter, , ] <- B_r
    samples$Sigma_wr[iter, , ] <- Sigma_wr
    samples$mu_v[iter, ] <- mu_v
    samples$Sigma_v[iter, , ] <- Sigma_v
    samples$v_t[iter, , ] <- v_t
    samples$lambda_v[iter, ] <- lambda_v_iter
    samples$lambda_g[iter, ] <- get.lambda_g(lambda_v_iter, eta_g, rho_g, S_bar)
    
    ## Optional diagnostics
    if (iter %% 50 == 0) {
      cat(sprintf("[Iter %d] lambda_g:\n", iter))
      cat(sprintf("  %s\n", paste(round(samples$lambda_g[iter, ], 3), collapse = "  ")))
    }
    if (iter %% 100 == 0) cat(sprintf("Iteration %d completed\n", iter))
  }
  
  # Drop burn-in
  keep <- (burnin + 1):num_iter
  
  return(list(
    sigma_wg2 = samples$sigma_wg2[keep],
    rho_g = samples$rho_g[keep, ],
    eta_g = samples$eta_g[keep, ],
    mu_r = samples$mu_r[keep, ],
    beta_v = samples$beta_v[keep, , ],
    B_r = samples$B_r[keep, , ],
    Sigma_wr = samples$Sigma_wr[keep, , ],
    mu_v = samples$mu_v[keep, ],
    Sigma_v = samples$Sigma_v[keep, , ],
    v_t = samples$v_t[keep, , ],
    lambda_v = samples$lambda_v[keep, ],
    lambda_g = samples$lambda_g[keep, ]
  ))
}


