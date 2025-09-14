# Note: I assume the inputs are the following:

# (1)
# - v_t a T x K matrix, where K is the number of latent factors and T the time step
# - mu_v is a 1 x K vector/matrix
# - eta_g is a K x 1 vector/matrix
# - rho_s is a (S_bar + 2) x 1 vector(matrix)

get.V_rho <- function(v_t, eta_g, mu_v, S_bar) {
  T <- nrow(v_t)  # number of time periods
  K <- ncol(v_t)  # number of latent factors
  V_rho <- matrix(1, nrow = T - S_bar, ncol = S_bar + 2)  # First column = intercept
  
  for (s in 0:S_bar) {
    # Select v_{t-s} for t = (S_bar+1):T, i.e., v_{S_bar+1 - s} to v_{T - s}
    v_lagged <- v_t[(S_bar + 1 - s):(T - s), , drop = FALSE]  # (T - S_bar) x K
    centered_v <- v_lagged - matrix(mu_v, nrow = T - S_bar, ncol = K, byrow = TRUE)
    
    # Project each row onto eta_g: result is a vector of length (T - S_bar)
    V_rho[, s + 2] <- centered_v %*% eta_g
  }
  
  return(V_rho)  # (T - S_bar) x (S_bar + 2)
}


# Function to compute V_eta
get.V_eta <- function(v_t, rho_s, mu_v, S_bar) {
  T <- nrow(v_t)       # time dimension
  K <- ncol(v_t)       # number of factors
  V_eta <- matrix(0, nrow = T - S_bar, ncol = K)  # output: (T - S_bar) x K
  
  for (s in 0:S_bar) {
    # v_{t-s} for t = S_bar+1,...,T → rows (S_bar+1 - s):(T - s)
    v_lagged <- v_t[(S_bar + 1 - s):(T - s), , drop = FALSE]  # (T - S_bar) x K
    centered_v <- v_lagged - matrix(mu_v, nrow = T - S_bar, ncol = K, byrow = TRUE)
    
    # Add lag with weight rho_s[s + 1]
    # Assumes: rho_s = c(mu_g, rho_0, ..., rho_S_bar)
    V_eta <- V_eta + rho_s[s + 2] * centered_v  # skip intercept and rho_0
  }
  
  return(V_eta)  # (T - S_bar) x K
}

# Function to compute G
get.G <- function(g_t, S_bar) {
  t <- length(g_t)
  return(matrix(g_t[(S_bar + 1):t], ncol = 1, nrow = t - S_bar))
}

# Function to compute G_bar
get.G_bar <- function(g_t, S_bar) {
  t <- length(g_t)
  mu_g <- mean(g_t)
  return(matrix(g_t[(S_bar + 1):t] - mu_g, ncol = 1, nrow = t - S_bar))
}

# (2)

# Input:
# - r_mat: a T x N matrix of excess returns
#          (rows = time, columns = assets)
# Output:
# - R: same matrix (identity function, included for symmetry)

get.R <- function(r_mat) {
  stopifnot(is.matrix(r_mat))
  return(r_mat)
}


# Input:
# - v_mat: a T x K matrix of latent factors
#          (rows = time, columns = latent factors)
# Output:
# - V_r: a T x (1 + K) matrix where the first column is 1 (intercept)
#        and the remaining are the centered latent factors

get.V_r <- function(v_mat) {
  stopifnot(is.matrix(v_mat))
  
  mu_v <- colMeans(v_mat)                    # Mean of each factor
  V_centered <- sweep(v_mat, 2, mu_v, "-")   # Center each column
  V_r <- cbind(1, V_centered)                # Add intercept column
  return(V_r)
}

# Input:
# - mu_r: an N-dimensional numeric vector of asset-specific intercepts
# - beta_v: a K x N matrix of factor loadings
# Output:
# - B_r: a (1 + K) x N matrix stacking mu_r and beta_v

get.B_r <- function(mu_r, beta_v) {
  stopifnot(length(mu_r) == ncol(beta_v))
  
  B_r <- rbind(mu_r, beta_v)  # (1 + K) x N
  return(B_r)
}

