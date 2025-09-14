rm(list = ls())
if (dev.cur() != 1) dev.off()

source("code/checks/simulation.data.R")
source("code/gibbs.sampler.R")

# Inputs
t <- 250     # time periods
K <- 4       # number of latent factors
N <- 50      # number of assets
num_iter <- 5000
burnin <- 1000
seed <- 42

# Simulate
data_sim <- simulate_data2(t, K, N, seed)
gibbs_res <- gibbs_sampler2(data_sim$R, data_sim$v_mat, num_iter = num_iter, burnin = burnin, seed = 123)

B_samples <- gibbs_res$B_r
Sigma_samples <- gibbs_res$Sigma_wr
mu_post <- drop(gibbs_res$mu_r)             # Drop dimensione [iter x 1 x N] -> [iter x N]
beta_post <- gibbs_res$beta_v               # [iter x K x N]

# Grid function
make_grid <- function(n) {
  rows <- ceiling(sqrt(n))
  cols <- ceiling(n / rows)
  return(c(rows, cols))
}

pdf("code/checks/plot/2check.gibbs.sampler.pdf", width = 11, height = 8.5)
par_default <- par(no.readonly = TRUE)

# ------------------------ TRACE PLOTS ------------------------
layout_dims <- make_grid(N)
par(mfrow = layout_dims, mar = c(2.5, 2.5, 2, 1))
for (j in 1:N) {
  diag_vals <- Sigma_samples[, j, j]
  plot(diag_vals, type = "l", main = bquote(Sigma[wr]^2[.(j)]), xlab = "Iter", ylab = "")
  abline(h = data_sim$Sigma_wr_true[j, j], col = "red", lwd = 2)
}

layout_dims <- make_grid(N)
par(mfrow = layout_dims, mar = c(2.5, 2.5, 2, 1))
for (j in 1:N) {
  plot(mu_post[, j], type = "l", main = bquote(mu[r][.(j)]), xlab = "Iter", ylab = "")
  abline(h = data_sim$mu_r_true[j], col = "red", lwd = 2)
}

for (k in 1:K) {
  layout_dims <- make_grid(N)
  par(mfrow = layout_dims, mar = c(2.5, 2.5, 2, 1))
  for (j in 1:N) {
    plot(beta_post[, k, j], type = "l", main = bquote(beta[.(k)][.(j)]), xlab = "Iter", ylab = "")
    abline(h = data_sim$beta_v_true[k, j], col = "red", lwd = 2)
  }
}

# ------------------------ DENSITY PLOTS ------------------------
layout_dims <- make_grid(N)
par(mfrow = layout_dims, mar = c(2.5, 2.5, 2, 1))
for (j in 1:N) {
  diag_vals <- Sigma_samples[, j, j]
  plot(density(diag_vals), main = bquote("Density of " ~ Sigma[wr]^2[.(j)]), xlab = "")
  abline(v = data_sim$Sigma_wr_true[j, j], col = "red", lwd = 2)
}

layout_dims <- make_grid(N)
par(mfrow = layout_dims, mar = c(2.5, 2.5, 2, 1))
for (j in 1:N) {
  plot(density(mu_post[, j]), main = bquote("Density of " ~ mu[r][.(j)]), xlab = "")
  abline(v = data_sim$mu_r_true[j], col = "red", lwd = 2)
}

for (k in 1:K) {
  layout_dims <- make_grid(N)
  par(mfrow = layout_dims, mar = c(2.5, 2.5, 2, 1))
  for (j in 1:N) {
    plot(density(beta_post[, k, j]), main = bquote("Density of " ~ beta[.(k)][.(j)]), xlab = "")
    abline(v = data_sim$beta_v_true[k, j], col = "red", lwd = 2)
  }
}

# ------------------------ ACF PLOTS ------------------------
layout_dims <- make_grid(N)
par(mfrow = layout_dims, mar = c(2.5, 2.5, 2, 1))
for (j in 1:N) {
  diag_vals <- Sigma_samples[, j, j]
  acf(diag_vals, main = bquote("ACF of " ~ Sigma[wr]^2[.(j)]))
}

layout_dims <- make_grid(N)
par(mfrow = layout_dims, mar = c(2.5, 2.5, 2, 1))
for (j in 1:N) {
  acf(mu_post[, j], main = bquote("ACF of " ~ mu[r][.(j)]))
}

for (k in 1:K) {
  layout_dims <- make_grid(N)
  par(mfrow = layout_dims, mar = c(2.5, 2.5, 2, 1))
  for (j in 1:N) {
    acf(beta_post[, k, j], main = bquote("ACF of " ~ beta[.(k)][.(j)]))
  }
}

# Close PDF
par(par_default)
dev.off()
