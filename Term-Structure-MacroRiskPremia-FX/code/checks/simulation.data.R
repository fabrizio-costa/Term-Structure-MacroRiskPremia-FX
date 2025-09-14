rm(list = ls())
if (!require(xts)) install.packages("xts"); library(xts)
if (!require(forecast)) install.packages("forecast"); library(forecast)

# Importing the posterior and matrix notation
source("code/matrix.notation.R")
source("code/posterior.R")

# ---------------------------- (1) Simulate data  ----------------------------
simulate_data <- function(t = 200, K = 3, S_bar = 2,
                          rho_g_true = c(0.7, -0.4, 0.2, 0.1),
                          sigma2_wg_true = 2,
                          eta_g_true = NULL) {
  
  # Step 1: simulate latent factors v_t as T x K
  v_t <- matrix(rnorm(t * K), nrow = t, ncol = K)
  mu_v <- colMeans(v_t)  # 1 x K
  
  # Step 2: use provided eta_g or simulate a random unit vector
  if (is.null(eta_g_true)) {
    eta_g_true <- rnorm(K)
  }
  eta_g_true <- eta_g_true / sqrt(sum(eta_g_true^2))  # normalize
  
  # Step 3: compute f_t = (v_t - mu_v) %*% eta_g_true
  v_centered <- v_t - matrix(mu_v, nrow = t, ncol = K, byrow = TRUE)  # T x K
  f_t <- as.numeric(v_centered %*% eta_g_true)  # T x 1 vector
  
  # Step 4: simulate g_t using MA(S_bar) convolution
  g_t <- numeric(t)
  for (tt in (S_bar + 1):t) {
    lag_sum <- sum(sapply(0:S_bar, function(s) rho_g_true[s + 2] * f_t[tt - s]))
    g_t[tt] <- rho_g_true[1] + lag_sum + rnorm(1, sd = sqrt(sigma2_wg_true))
  }
  
  # Return simulated data
  return(list(
    g_t = g_t,
    v_t = v_t,
    eta_g_true = eta_g_true,
    rho_g_true = rho_g_true,
    sigma_wg2_true = sigma2_wg_true
  ))
}

# ---------------------------- (2) Simulate data  ----------------------------

simulate_data2 <- function(time, K, N, seed = 123) {
  set.seed(seed)
  
  mu_r_true <- rnorm(N, mean = 0.01, sd = 0.02)
  beta_v_true <- matrix(rnorm(K * N, mean = 0.1, sd = 0.2), nrow = K)
  Sigma_wr_true <- diag(runif(N, 0.1, 2))
  
  v_mat <- matrix(rnorm(time * K), nrow = time)
  V_r <- get.V_r(v_mat)
  B_r_true <- get.B_r(mu_r_true, beta_v_true)
  
  R <- V_r %*% B_r_true + matrix(rnorm(time * N), time, N) %*% chol(Sigma_wr_true)
  
  list(R = R, v_mat = v_mat, mu_r_true = mu_r_true,
       beta_v_true = beta_v_true, Sigma_wr_true = Sigma_wr_true)
}

# ---------------------------- (3) Simulate data  ----------------------------
simulate_data3 <- function(t = 500, K = 3, N = 10, seed = 123) {
  set.seed(seed)
  
  mu_r_true <- rnorm(N, mean = 0.01, sd = 0.02)
  beta_v_true <- matrix(rnorm(K * N), nrow = K, ncol = N)
  
  Sigma_wr_true <- diag(runif(N, 0.1, 0.5))
  mu_v_true <- rnorm(K, 0, 0.5)
  
  A <- matrix(rnorm(K * K), K)
  Sigma_v_true <- crossprod(A)
  
  # Add ridge to ensure Sigma_v_true is well-conditioned
  ridge <- 1e-3
  Sigma_v_true <- Sigma_v_true + ridge * diag(K)
  
  # Force minimum eigenvalue
  eig <- eigen(Sigma_v_true)
  min_eig <- 0.05
  eig$values[eig$values < min_eig] <- min_eig
  Sigma_v_true <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  
  # Check conditioning
  # cat("Sigma_v_true eigenvalues:\n")
  # print(round(eigen(Sigma_v_true)$values, 6))
  # cat("Condition number of Sigma_v_true:", kappa(Sigma_v_true), "\n")
  
  v_t_true <- MASS::mvrnorm(n = t, mu = mu_v_true, Sigma = Sigma_v_true)
  
  # Check variability of simulated latent factors
  # cat("v_t_true column-wise SD:\n")
  # print(round(apply(v_t_true, 2, sd), 6))
  
  R <- matrix(0, nrow = t, ncol = N)
  for (tt in 1:t) {
    mean_r_t <- mu_r_true + t(beta_v_true) %*% v_t_true[tt, ]
    R[tt, ] <- MASS::mvrnorm(n = 1, mu = as.numeric(mean_r_t), Sigma = Sigma_wr_true)
  }
  
  return(list(
    R = R,
    v_t_true = v_t_true,
    mu_r_true = mu_r_true,
    beta_v_true = beta_v_true,
    Sigma_wr_true = Sigma_wr_true,
    mu_v_true = mu_v_true,
    Sigma_v_true = Sigma_v_true
  ))
}


# ---------------------------- (4) Simulate data  ----------------------------
simulate_data4 <- function(K = 3, N = 10, S_bar = 2, tau = 0.5, seed = 42) {
  set.seed(seed)
  
  # Step 1: Simulate beta_v (K x N)
  beta_v_true <- matrix(rnorm(K * N, mean = 0.1, sd = 0.2), nrow = K, ncol = N)
  
  # Step 2: Simulate mu_r (N-vector)
  mu_r_true <- rnorm(N, mean = 0.01, sd = 0.02)
  
  # Step 3: Simulate Sigma_wr (N x N diagonal)
  Sigma_wr_true <- diag(runif(N, min = 0.1, max = 0.3))
  
  # Step 4: Simulate eta_g (K-vector, normalized)
  eta_g_true <- rnorm(K)
  eta_g_true <- eta_g_true / sqrt(sum(eta_g_true^2))  # unit norm
  
  # Step 5: Simulate rho_g ~ N(0, tau^2)
  rho_g_true <- rnorm(S_bar + 2, mean = 0, sd = tau)
  
  # Step 6: Compute lambda_v (true)
  BtB <- beta_v_true %*% t(beta_v_true)
  ridge <- 1e-6
  if (rcond(BtB) < 1e-12) BtB <- BtB + ridge * diag(K)
  
  mu_r_tilde <- mu_r_true + 0.5 * diag(Sigma_wr_true)
  lambda_v_true <- solve(BtB, beta_v_true %*% mu_r_tilde)  # K x 1
  
  # Step 7: Compute lambda_g^S for S = 0,...,S_bar
  lambda_g_true <- numeric(S_bar + 1)
  rho_terms <- rho_g_true[2:(S_bar + 2)]  # exclude intercept
  eta_lambda <- sum(eta_g_true * lambda_v_true)
  
  for (S in 0:S_bar) {
    sum_rho <- 0
    for (tau_i in 0:S) {
      sum_rho <- sum_rho + sum(rho_terms[1:(tau_i + 1)])
    }
    lambda_g_true[S + 1] <- (sum_rho / (1 + S)) * eta_lambda
  }
  
  return(list(
    K = K,
    N = N,
    S_bar = S_bar,
    tau = tau,
    beta_v_true = beta_v_true,
    mu_r_true = mu_r_true,
    Sigma_wr_true = Sigma_wr_true,
    eta_g_true = eta_g_true,
    rho_g_true = rho_g_true,
    lambda_v_true = as.numeric(lambda_v_true),
    lambda_g_true = lambda_g_true
  ))
}




# ---------------------------- (5) Simulate Data Full  ----------------------------

# log returns and log changes in macroeconomic variables
simulate_data_full <- function(t = 300, K = 3, N = 10, S_bar = 2, seed = 123) {
  set.seed(seed)
  
  # ---------- (1) Parameters ----------
  mu_v_true <- rnorm(K, mean = 0, sd = 0.5)
  
  # Sigma_v_true: positive-definite, well-conditioned
  eig_vals <- seq(1, 0.1, length.out = K)
  Q <- qr.Q(qr(matrix(rnorm(K^2), K)))
  Sigma_v_true <- Q %*% diag(eig_vals) %*% t(Q)
  
  eta_g_true <- rnorm(K)
  eta_g_true <- eta_g_true / sqrt(sum(eta_g_true^2))  # normalize
  
  mu_g_true <- 0.3
  tau <- 0.2
  rho_g_true <- c(mu_g_true, rnorm(S_bar + 1, mean = 0, sd = tau))
  sigma2_wg_true <- 0.2
  
  mu_r_true <- rnorm(N, mean = 0.01, sd = 0.02)
  beta_v_true <- matrix(rnorm(K * N, mean = 0.1, sd = 0.2), nrow = K, ncol = N)
  Sigma_wr_true <- diag(runif(N, min = 0.1, max = 0.3))
  
  # ---------- (2) Simulate latent factors ----------
  v_t_true <- MASS::mvrnorm(n = t, mu = mu_v_true, Sigma = Sigma_v_true)
  
  # ---------- (3) Compute f_t and g_t ----------
  v_centered <- v_t_true - matrix(mu_v_true, nrow = t, ncol = K, byrow = TRUE)
  f_t <- as.numeric(v_centered %*% eta_g_true)
  
  # Pad f_t to safely simulate full g_t with MA(S_bar)
  f_t_padded <- c(rep(mean(f_t), S_bar), f_t)
  g_t <- numeric(t)
  for (tt in 1:t) {
    lag_sum <- sum(sapply(0:S_bar, function(s) rho_g_true[s + 2] * f_t_padded[tt + S_bar - s]))
    g_t[tt] <- rho_g_true[1] + lag_sum + rnorm(1, sd = sqrt(sigma2_wg_true))
  }
  
  # ---------- (4) Simulate R ----------
  R <- matrix(0, nrow = t, ncol = N)
  for (tt in 1:t) {
    mu_r_t <- mu_r_true + t(beta_v_true) %*% v_t_true[tt, ]
    R[tt, ] <- MASS::mvrnorm(n = 1, mu = as.numeric(mu_r_t), Sigma = Sigma_wr_true)
  }
  
  # ---------- (5) Compute lambda_v and lambda_g ----------
  mu_r_tilde <- mu_r_true + 0.5 * diag(Sigma_wr_true)
  BtB <- beta_v_true %*% t(beta_v_true)
  if (rcond(BtB) < 1e-12) BtB <- BtB + 1e-6 * diag(K)
  lambda_v_true <- solve(BtB, beta_v_true %*% mu_r_tilde)
  
  lambda_g_true <- numeric(S_bar + 1)
  for (S in 0:S_bar) {
    sum_rho <- sum(sapply(0:S, function(tau) sum(rho_g_true[2:(tau + 2)])))
    lambda_g_true[S + 1] <- (sum_rho / (1 + S)) * sum(eta_g_true * lambda_v_true)
  }
  
  # ---------- (6) Return ----------
  return(list(
    g_t = g_t,
    R = R,
    v_t_true = v_t_true,
    mu_v_true = mu_v_true,
    Sigma_v_true = Sigma_v_true,
    eta_g_true = eta_g_true,
    rho_g_true = rho_g_true,
    sigma2_wg_true = sigma2_wg_true,
    mu_r_true = mu_r_true,
    beta_v_true = beta_v_true,
    Sigma_wr_true = Sigma_wr_true,
    lambda_v_true = as.numeric(lambda_v_true),
    lambda_g_true = lambda_g_true,
    S_bar = S_bar
  ))
}
