if (!require(MASS)) install.packages("MASS"); library(MASS)
if (!require(MCMCpack)) install.packages("MCMCpack"); library(MCMCpack)

# importing the matrix notations V_rho, V_eta, G and G_bar
source("code/matrix.notation.R")

# ---------------------------------- (1) ----------------------------------

# Function to sample from inverse-gamma distribution
rinvgamma <- function(n, shape_p, rate_p) {
  return(1 / stats::rgamma(n, shape = shape_p, rate = rate_p))
}

# get.sigma.post
# ----------------------------------
# Input:
# - g_t: numeric vector of length T-S (the factor time series)
# - rho_g: numeric vector (coefficients for lag polynomial)
# - eta_g: numeric vector of length K (direction in factor space)
# - v_t: T x K matrix of latent factors (columns = time)
# - S_bar: integer, number of lags in the MA representation
#
# Output:
# - sigma_wg2: scalar, posterior draw for the idiosyncratic variance of g_t

# Function to compute posterior draw of sigma^2_wg
get.sigma.post <- function(g_t, rho_g, eta_g, v_t, S_bar) {
  t <- length(g_t)             # g_t is a vector of length T
  K <- ncol(v_t)               # v_t is T x K (rows = time, cols = factors)
  
  G <- get.G(g_t, S_bar)       # (T - S_bar) x 1
  mu_v <- colMeans(v_t)        # 1 x K
  V_rho <- get.V_rho(v_t, eta_g, mu_v, S_bar)  # (T - S_bar) x (S_bar + 2)
  
  residuals <- G - V_rho %*% rho_g             # (T - S_bar) x 1
  df <- (t - S_bar) / 2
  rate_param <- as.numeric(t(residuals) %*% residuals) / 2  # scalar
  
  sigma_wg2 <- rinvgamma(1, shape = df, rate = rate_param)
  return(sigma_wg2)
}


# get.rho_g
# ----------------------------------
# Input:
# - G: vector of length T - S_bar (dependent variable in g_t equation)
# - sigma_wg2: scalar variance
# - V_rho: (T - S_bar) x (S_bar + 2) matrix of lagged factor loadings
#
# Output:
# - rho_g: vector of coefficients (length S_bar + 2)

get.rho_g <- function(G, sigma_wg2, V_rho, tau2_rho = 1) {
  XtX <- t(V_rho) %*% V_rho 
  
  # Conditional Ridge 
  if (rcond(XtX) < 1e-12) {
    XtX <- XtX + 1e-6 * diag(ncol(XtX))
  }
  
  Sigma_rho <- solve(XtX)
  mean_rho <- Sigma_rho %*% t(V_rho) %*% G
  rho_g <- MASS::mvrnorm(1, mean_rho, Sigma_rho * sigma_wg2) 
  return(rho_g)
}


# get.rho_g.adjusted
# ----------------------------------
# Input:
# - G: vector of length T - S_bar (dependent variable in g_t equation)
# - sigma_wg2: scalar variance (not directly used here, robust inference instead)
# - V_rho: (T - S_bar) x (S_bar + 2) matrix of lagged factor loadings
# - S_bar: integer (used for Newey-West truncation lag)
#
# Output:
# - rho_g: vector of coefficients (length S_bar + 2) with NW-adjusted covariance

get.rho_g.adjusted <- function(G, sigma_wg2, V_rho, S_bar) {
  T_eff <- nrow(V_rho)
  p <- ncol(V_rho)
  L <- S_bar
  
  # Conditional ridge
  XtX <- t(V_rho) %*% V_rho
  if (rcond(XtX) < 1e-12) {
    XtX <- XtX + 1e-6 * diag(p)
  }
  XtX_inv <- solve(XtX)
  
  rho_hat <- XtX_inv %*% t(V_rho) %*% G
  w_hat <- G - V_rho %*% rho_hat
  
  S_hat <- matrix(0, p, p)
  for (t in 1:T_eff) {
    v_t <- V_rho[t, , drop = FALSE]
    S_hat <- S_hat + t(v_t) %*% v_t * as.numeric(w_hat[t]^2)
  }
  S_hat <- S_hat / T_eff
  
  for (l in 1:L) {
    Gamma_l <- matrix(0, p, p)
    for (t in (l + 1):T_eff) {
      v_t    <- V_rho[t, , drop = FALSE]
      v_lag  <- V_rho[t - l, , drop = FALSE]
      wt     <- w_hat[t]
      wt_lag <- w_hat[t - l]
      Gamma_l <- Gamma_l + t(v_t) %*% v_lag * as.numeric(wt * wt_lag)
    }
    Gamma_l <- Gamma_l / T_eff
    weight <- 1 - l / (L + 1)
    S_hat <- S_hat + weight * (Gamma_l + t(Gamma_l))
  }
  
  Sigma_rho <- XtX_inv %*% (T_eff * S_hat) %*% XtX_inv
  rho_g <- MASS::mvrnorm(1, mu = as.numeric(rho_hat), Sigma = Sigma_rho)
  return(rho_g)
}

# get.eta_g
# ----------------------------------
# Input:
# - G_bar: (T - S_bar) x 1 vector (rotated g_t observations)
# - sigma_wg2: scalar variance
# - rho_g: vector of coefficients (length S_bar + 2)
# - V_eta: (T - S_bar) x K matrix for eta estimation
#
# Output:
# - eta_sample: normalized K-length vector (||eta|| = 1)

get.eta_g <- function(G_bar, sigma_wg2, rho_g, V_eta, shrink_eta = FALSE) {
  XtX <- t(V_eta) %*% V_eta
  
  if (rcond(XtX) < 1e-12) {
    XtX <- XtX + 1e-6 * diag(ncol(XtX))
  }
  
  tau_eta2 = 0.01
  Sigma_eta <- solve(XtX)
  
  mean_eta <- Sigma_eta %*% t(V_eta) %*% G_bar
  
  # Optional shrinkage toward prior mean
  if (shrink_eta) {
    K <- ncol(V_eta)
    eta_prior <- rep(1 / sqrt(K), K)  # simple prior direction
    Sigma_prior <- diag(tau_eta2, K)  # prior variance
    
    precision_post <- solve(Sigma_eta * sigma_wg2)
    precision_prior <- solve(Sigma_prior)
    
    Sigma_combined <- solve(precision_post + precision_prior)
    mean_combined <- Sigma_combined %*% (precision_post %*% mean_eta + precision_prior %*% eta_prior)
    
    eta_sample <- MASS::mvrnorm(1, mean_combined, Sigma_combined)
  } else {
    eta_sample <- MASS::mvrnorm(1, mean_eta, Sigma_eta * sigma_wg2)
  }
  
  # Normalize
  eta_sample <- eta_sample / sqrt(sum(eta_sample^2))
  return(as.numeric(eta_sample))
}


# ---------------------------------- (2) ----------------------------------

# sample_Sigma_wr
# ----------------------------------
# Input:
# - R: T x N matrix of asset returns
# - v_mat: T x K matrix of latent factors
# - mu_r: N-length numeric vector
# - beta_v: K x N matrix
#
# Output:
# - Sigma_wr: N x N covariance matrix sampled from inverse-Wishart

get.Sigma_wr <- function(R, v_mat, mu_r, beta_v) {
  stopifnot(is.matrix(R), is.matrix(v_mat), is.matrix(beta_v))
  
  t <- nrow(R)
  V_r <- get.V_r(v_mat)                # T x (1 + K)
  B_r <- get.B_r(mu_r, beta_v)         # (1 + K) x N
  
  residuals <- R - V_r %*% B_r         # T x N
  S <- t(residuals) %*% residuals      # N x N
  
  # Soft protection only when needed
  eigvals <- eigen(S, only.values = TRUE)$values
  if (min(eigvals) < 1e-8) {
    warning("S is nearly singular. Smallest eigenvalue = ", signif(min(eigvals), 3))
    S <- S + 1e-6 * diag(ncol(S))
  }
  
  Sigma_wr <- MCMCpack::riwish(t, S)
  return(Sigma_wr)
}

# Posterior for B_r = [mu_r; beta_v]
# --------------------------------------
# Input:
# - R: T x N matrix of returns
# - v_mat: T x K matrix of latent factors
# - Sigma_wr: N x N diagonal covariance matrix
# Output:
# - B_r: (1 + K) x N matrix

# sample_B_r
# ----------------------------------
# Input:
# - R: T x N matrix of asset returns
# - v_mat: T x K matrix of latent factors
# - Sigma_wr: N x N covariance matrix (from inverse-Wishart)
#
# Output:
# - B_r: (1 + K) x N matrix of coefficients sampled from matrix-normal

get.B_r.post <- function(R, v_mat, Sigma_wr, tau2 = 0.01) {
  stopifnot(is.matrix(R), is.matrix(v_mat), is.matrix(Sigma_wr))
  
  V_r <- get.V_r(v_mat)                    # T x (1 + K)
  K_plus_1 <- ncol(V_r)
  N <- ncol(R)
  
  # Ridge-regularized V'V
  VtV <- crossprod(V_r) + (1 / tau2) * diag(K_plus_1)
  VtV_inv <- solve(VtV)
  VtR <- crossprod(V_r, R)                 # (1 + K) x N
  
  M_post <- VtV_inv %*% VtR                # Posterior mean
  
  B_r <- matrix(NA, nrow = K_plus_1, ncol = N)
  for (j in 1:N) {
    B_r[, j] <- MASS::mvrnorm(
      n = 1,
      mu = M_post[, j],
      Sigma = VtV_inv * Sigma_wr[j, j]
    )
  }
  
  return(B_r)
}


# ---------------------------------- (3) ----------------------------------

# Get.v_t
# ----------------------------------
# Input:
# - r_t: N-dimensional numeric vector (excess returns at time t)
# - mu_r: N-dimensional numeric vector (asset-specific intercepts)
# - beta_v: K x N matrix of factor loadings
# - Sigma_wr: N x N covariance matrix of returns
# - mu_v: K-dimensional numeric vector (prior mean of latent factors)
# - Sigma_v: K x K prior covariance matrix for latent factors
#
# Output:
# - v_t: K-dimensional numeric vector (sample from posterior of latent factor at time t)

get.v_t <- function(r_t, mu_r, beta_v, Sigma_wr, mu_v, Sigma_v) {
  stopifnot(length(r_t) == length(mu_r))
  stopifnot(ncol(Sigma_wr) == length(r_t))
  stopifnot(nrow(beta_v) == length(mu_v))
  
  mu_v <- matrix(mu_v, ncol = 1)  # ensure column vector
  ridge <- 1e-6
  K <- length(mu_v)
  
  # Regularize Sigma_wr if ill-conditioned
  if (rcond(Sigma_wr) < 1e-12) {
    Sigma_wr <- Sigma_wr + ridge * diag(nrow(Sigma_wr))
  }
  Sigma_wr_inv <- solve(Sigma_wr)
  
  # Regularize Sigma_v if ill-conditioned
  if (rcond(Sigma_v) < 1e-12) {
    Sigma_v <- Sigma_v + ridge * diag(nrow(Sigma_v))
  }
  Sigma_v_inv <- solve(Sigma_v)
  
  # Posterior covariance and mean
  precision_post <- beta_v %*% Sigma_wr_inv %*% t(beta_v) + Sigma_v_inv
  if (rcond(precision_post) < 1e-12) {
    precision_post <- precision_post + ridge * diag(K)
  }
  Sigma_post <- solve(precision_post)
  
  eig_vals <- eigen(Sigma_post, only.values = TRUE)$values
  min_eig <- min(eig_vals)
  
  if (min_eig <= 0 || any(is.nan(eig_vals))) {
    cat("Sigma_post not positive definite!\n")
    cat("Min eigenvalue:", min_eig, "\n")
    cat("Condition number:", kappa(Sigma_post), "\n")
    stop("Aborting: Sigma_post not positive definite")
  }
  
  
  mean_post <- Sigma_post %*% (
    beta_v %*% Sigma_wr_inv %*% (r_t - mu_r) +
      Sigma_v_inv %*% mu_v
  )
  
  # Final ridge if still problematic
  if (!isSymmetric(Sigma_post)) {
    Sigma_post <- (Sigma_post + t(Sigma_post)) / 2
  }
  if (any(eigen(Sigma_post, only.values = TRUE)$values <= 0)) {
    cat("Sigma_post not PD → applying extra ridge.\n")
    Sigma_post <- Sigma_post + 1e-4 * diag(K)
  }
  
  # Safe sampling
  v_t <- tryCatch(
    MASS::mvrnorm(1, mu = as.numeric(mean_post), Sigma = Sigma_post),
    error = function(e) {
      cat("MASS::mvrnorm failed — using fallback mean draw\n")
      return(as.numeric(mean_post))  # fallback: just return the mean
    }
  )
  
  return(v_t)
}


# sample Sigma_v 
# ------------------------------------------
# Inputs:
# - v_t: T x K matrix of latent factors (each column is v_t at time t)
#
# Output:
# - Sigma_v: K x K sampled covariance matrix

get.Sigma_v <- function(v_t, shrinkage = FALSE, verbose = FALSE, iter = NULL) {
  t <- nrow(v_t)
  K <- ncol(v_t)
  ridge <- 1e-4
  shrink_level <- 0.3
  
  v_bar <- colMeans(v_t)
  centered_v <- sweep(v_t, 2, v_bar, "-")
  empirical <- t(centered_v) %*% centered_v  # empirical covariance * (T - 1)
  
  # Apply shrinkage directly if requested
  scale_matrix <- if (shrinkage) {
    (shrink_level * empirical + (1 - shrink_level) * diag(K))
  } else {
    empirical
  }
  
  # Ridge stabilization if needed
  cond_number <- kappa(scale_matrix)
  if (cond_number > 1e8) {
    if (verbose && !is.null(iter)) {
      cat(sprintf("[%04d] Sigma_v cond > 1e8 — adding ridge\n", iter))
    }
    scale_matrix <- scale_matrix + ridge * diag(K)
  }
  
  # Inverse-Wishart draw
  Sigma_v <- tryCatch({
    MCMCpack::riwish(t - 1, scale_matrix)
  }, error = function(e) {
    cat("Error in Sigma_v draw. Iter:", iter, "\n")
    print(scale_matrix)
    stop(e)
  })
  
  # Validate output
  eig_vals <- eigen(Sigma_v, only.values = TRUE)$values
  if (any(is.nan(eig_vals)) || any(eig_vals <= 0)) {
    stop("Drawn Sigma_v is not positive definite")
  }
  
  return(Sigma_v)
}

# Sample mu_v 
# ----------------------------------------
# Inputs:
# - v_t: T x K matrix of latent factors (each column is v_t at time t)
# - Sigma_v: K x K covariance matrix (current draw or prior)
#
# Output:
# - mu_v: K-dimensional vector (posterior draw of the factor mean)

get.mu_v <- function(v_t, Sigma_v, verbose = FALSE, iter = NULL) {
  T <- nrow(v_t)
  K <- ncol(v_t)
  ridge <- 1e-6
  
  v_bar <- colMeans(v_t)
  
  # Stabilize Sigma_v if needed
  if (rcond(Sigma_v) < 1e-12) {
    if (verbose && !is.null(iter)) {
      cat(sprintf("[%04d] Sigma_v ill-conditioned — adding ridge\n", iter))
    }
    Sigma_v <- Sigma_v + ridge * diag(K)
  }
  
  Sigma_post <- Sigma_v / T
  
  # Ensure symmetry
  if (!isSymmetric(Sigma_post)) {
    Sigma_post <- (Sigma_post + t(Sigma_post)) / 2
  }
  
  # Check PD
  eig_vals <- eigen(Sigma_post, only.values = TRUE)$values
  if (any(eig_vals <= 0 | is.nan(eig_vals))) {
    cat(sprintf("[%04d] Sigma_post for mu_v not PD — fallback to mean\n", iter))
    return(as.numeric(v_bar))
  }
  
  # Sample safely
  mu_v <- tryCatch(
    MASS::mvrnorm(1, mu = v_bar, Sigma = Sigma_post),
    error = function(e) {
      cat(sprintf("[%04d] mvrnorm failed for mu_v — using mean\n", iter))
      return(as.numeric(v_bar))
    }
  )
  
  return(mu_v)
}

# ---------------------------------- (4) ----------------------------------

# get.lambda_v
# ----------------------------------
# Input:
# - beta_v: K x N matrix of factor loadings
# - mu_r: N-dimensional vector (intercepts from return equation)
# - Sigma_wr: N x N covariance matrix of asset returns
#
# Output:
# - lambda_v: K-dimensional vector of risk premia on latent factors

get.lambda_v <- function(beta_v, mu_r, Sigma_wr) {
  stopifnot(is.matrix(beta_v), length(mu_r) == ncol(beta_v))
  
  K <- nrow(beta_v)
  N <- ncol(beta_v)
  
  # Step 1: Compute tilde mu_r = mu_r + 0.5 * diag(Sigma_wr)
  mu_r_tilde <- mu_r + 0.5 * diag(Sigma_wr)  # N × 1
  
  # Step 2: Compute BtB = beta_v %*% t(beta_v) → K × K
  BtB <- beta_v %*% t(beta_v)
  ridge <- 1e-6
  if (rcond(BtB) < 1e-12) {
    BtB <- BtB + ridge * diag(K)
  }
  
  # Step 3: Compute lambda_v = (BtB)^{-1} * beta_v * mu_r_tilde
  lambda_v <- solve(BtB, beta_v %*% mu_r_tilde)  # (K × N) × (N × 1) = K × 1
  return(as.numeric(lambda_v))  # convert to numeric vector
}

# get.lambda_g
# ----------------------------------
# Input:
# - lambda_v: K-dimensional vector (from get.lambda_v)
# - eta_g: K-dimensional normalized vector
# - rho_g: (S_bar + 2)-length vector (intercept + rho_0 to rho_Sbar)
# - S_bar: integer (maximum lag)
#
# Output:
# - lambda_g: (S_bar + 1)-length numeric vector with lambda_g^S for S = 0,...,S_bar

get.lambda_g <- function(lambda_v, eta_g, rho_g, S_bar) {
  stopifnot(length(lambda_v) == length(eta_g))
  
  rho_terms <- rho_g[2:(S_bar + 2)]  # extract rho_0 to rho_Sbar (exclude intercept)
  
  # Precompute ηᵗ λ
  eta_lambda <- sum(eta_g * lambda_v)
  
  # Compute λ_g^S for S = 0,...,S_bar
  lambda_g <- numeric(S_bar + 1)
  for (S in 0:S_bar) {
    sum_rho <- 0
    for (tau in 0:S) {
      sum_rho <- sum_rho + sum(rho_terms[1:(tau + 1)])
    }
    lambda_g[S + 1] <- (sum_rho / (1 + S)) * eta_lambda
  }
  
  return(lambda_g)
}



