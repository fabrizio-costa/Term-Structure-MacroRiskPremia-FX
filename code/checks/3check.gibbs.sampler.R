rm(list = ls())
if (dev.cur() != 1) dev.off()

source("code/checks/simulation.data.R")
source("code/gibbs.sampler.R")

# Settings
time <- 250
K <- 4
N <- 50
num_iter <- 5000
burnin <- 1000
seed <- 42

# Simulate data and run Gibbs sampler
data_sim <- simulate_data3(t = time, K = K, N = N, seed = seed)
gibbs_res <- gibbs_sampler3(
  R = data_sim$R,
  mu_r = data_sim$mu_r_true,
  beta_v = data_sim$beta_v_true,
  Sigma_wr = data_sim$Sigma_wr_true,
  num_iter = num_iter,
  burnin = burnin,
  seed = seed
)

# Extract post-burn-in draws
mu_samples <- gibbs_res$mu_v
Sigma_samples <- gibbs_res$Sigma_v
v_t_samples <- gibbs_res$v_t  # [iter, time, K]

# Helper: Layout
make_grid <- function(n) {
  rows <- ceiling(sqrt(n))
  cols <- ceiling(n / rows)
  return(c(rows, cols))
}

pdf("code/checks/plot/3check.gibbs.sampler.pdf", width = 11, height = 8.5)
par_default <- par(no.readonly = TRUE)

# --------- TRACE PLOTS ----------
# mu_v
layout_dims <- make_grid(K)
par(mfrow = layout_dims)
for (k in 1:K) {
  plot(mu_samples[, k], type = "l",
       main = bquote(mu[v][.(k)]), xlab = "Iter", ylab = "")
  abline(h = data_sim$mu_v_true[k], col = "red", lwd = 2)
}

# diag elements of Sigma_v
layout_dims <- make_grid(K)
par(mfrow = layout_dims)
for (k in 1:K) {
  plot(Sigma_samples[, k, k], type = "l",
       main = bquote(Sigma[v][.(k) * "," * .(k)]), xlab = "Iter", ylab = "")
  abline(h = data_sim$Sigma_v_true[k, k], col = "red", lwd = 2)
}

# v_{t, k}: select specific time points
selected_times <- c(1, floor(time / 2), time)
layout_dims <- make_grid(length(selected_times) * K)
par(mfrow = layout_dims)

for (t in selected_times) {
  for (k in 1:K) {
    v_series <- v_t_samples[, t, k]  
    plot(v_series, type = "l",
         main = bquote("Trace of " ~ v[.(k) * "," * .(t)]),
         xlab = "Iter", ylab = "")
    abline(h = data_sim$v_t_true[t, k], col = "red", lwd = 2)  
  }
}

# --------- DENSITY PLOTS ----------
layout_dims <- make_grid(K)
par(mfrow = layout_dims)
for (k in 1:K) {
  plot(density(mu_samples[, k]), main = bquote("Density of " ~ mu[v][.(k)]), xlab = "")
  abline(v = data_sim$mu_v_true[k], col = "red", lwd = 2)
}

layout_dims <- make_grid(K)
par(mfrow = layout_dims)
for (k in 1:K) {
  plot(density(Sigma_samples[, k, k]), main = bquote("Density of " ~ Sigma[v][.(k) * "," * .(k)]), xlab = "")
  abline(v = data_sim$Sigma_v_true[k, k], col = "red", lwd = 2)
}

layout_dims <- make_grid(length(selected_times) * K)
par(mfrow = layout_dims)
for (t in selected_times) {
  for (k in 1:K) {
    v_series <- v_t_samples[, t, k]  
    v_true_val <- data_sim$v_t_true[t, k]  
    plot(density(v_series),
         main = bquote("Density of " ~ v[.(k) * "," * .(t)]),
         xlab = "")
    abline(v = v_true_val, col = "red", lwd = 2)
  }
}

# --------- ACF PLOTS ----------
layout_dims <- make_grid(K)
par(mfrow = layout_dims)
for (k in 1:K) {
  acf(mu_samples[, k], main = bquote("ACF of " ~ mu[v][.(k)]))
}

layout_dims <- make_grid(K)
par(mfrow = layout_dims)
for (k in 1:K) {
  acf(Sigma_samples[, k, k], main = bquote("ACF of " ~ Sigma[v][.(k) * "," * .(k)]))
}

layout_dims <- make_grid(length(selected_times) * K)
par(mfrow = layout_dims)
for (t in selected_times) {
  for (k in 1:K) {
    v_series <- v_t_samples[, t, k]  
    acf(v_series, main = bquote("ACF of " ~ v[.(k) * "," * .(t)]))
  }
}

# closing pdf
par(par_default)
dev.off()
