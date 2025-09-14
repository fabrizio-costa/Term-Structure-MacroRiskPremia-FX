rm(list = ls())
if (dev.cur() != 1) dev.off()

# Load required code
source("code/posterior.R")
source("code/checks/simulation.data.R")

# ---------------------------- (1) Simulate minimal data for Step 4 ----------------------------
set.seed(42)
K <- 5
N <- 50
S_bar <- 25

sim_data <- simulate_data4(K = K, N = N, S_bar = S_bar)

beta_v_true <- sim_data$beta_v_true
mu_r_true <- sim_data$mu_r_true
Sigma_wr_true <- sim_data$Sigma_wr_true
eta_g_true <- sim_data$eta_g_true
rho_g_true <- sim_data$rho_g_true

lambda_v_true <- sim_data$lambda_v_true
lambda_g_true <- sim_data$lambda_g_true

# ---------------------------- (2) Generate posterior draws ----------------------------

B <- 1000
lambda_v_draws <- matrix(NA, nrow = B, ncol = K)
lambda_g_draws <- matrix(NA, nrow = B, ncol = S_bar + 1)

for (b in 1:B) {
  mu_r_b <- MASS::mvrnorm(1, mu = mu_r_true, Sigma = diag(0.01, N))
  Sigma_wr_b <- Sigma_wr_true
  
  lambda_v_b <- get.lambda_v(beta_v_true, mu_r_b, Sigma_wr_b)
  lambda_g_b <- get.lambda_g(lambda_v_b, eta_g_true, rho_g_true, S_bar)
  
  lambda_v_draws[b, ] <- lambda_v_b
  lambda_g_draws[b, ] <- lambda_g_b
}

# ---------------------------- (3) Plot histograms to PDF ----------------------------
pdf("code/checks/plot/4posterior.checks.pdf", width = 11, height = 8.5)

## ---- lambda_v histograms ----
par(mfrow = c(1, K), mar = c(3, 2.5, 3, 1))  # Adjust margins
for (k in 1:K) {
  hist(lambda_v_draws[, k], breaks = 40, col = "lightgray",
       main = bquote(lambda[v][.(k)]),
       xlab = "", ylab = "")
  abline(v = lambda_v_true[k], col = "red", lwd = 2)
}

## ---- lambda_g^S histograms across multiple pages ----
for (s in 1:(S_bar + 1)) {
  par(mfrow = c(1, 1), mar = c(4, 4, 4, 2))  # Reset for single plot per page
  hist(lambda_g_draws[, s], breaks = 40, col = "lightgray",
       main = bquote(lambda[g]^.(s - 1)),
       xlab = "", ylab = "")
  abline(v = lambda_g_true[s], col = "red", lwd = 2)
}

dev.off()
