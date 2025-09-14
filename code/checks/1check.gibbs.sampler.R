rm(list = ls())
if (dev.cur() != 1) dev.off()

# Load required sources
source("code/checks/simulation.data.R")
source("code/gibbs.sampler.R")

# ---------------------------- Simulate Data ----------------------------
set.seed(42)
t <- 250
k <- 4
s_bar <- 24
tau <- 2
rho_g.sim <- rnorm(s_bar + 2, mean = 0, sd = tau)
sigma2.sim <- 0.2
eta_g.sim <- rnorm(k)
eta_g.sim <- eta_g.sim / sqrt(sum(eta_g.sim^2))

sim_data <- simulate_data(t = t, K = k, S_bar = s_bar, 
                          rho_g_true = rho_g.sim, 
                          sigma2_wg_true = sigma2.sim, 
                          eta_g_true = eta_g.sim)

g_t <- sim_data$g_t
v_t <- sim_data$v_t

# Run Gibbs sampler
sigma.rho.eta.post <- gibbs_sampler1(g_t, v_t, s_bar, num_iter = 5000, burn_in = 1000)
sigma.post <- sigma.rho.eta.post$sigma_wg2
rho.post   <- sigma.rho.eta.post$rho_g
eta.post   <- sigma.rho.eta.post$eta_g

# Sign alignment
eta_mean <- colMeans(eta.post)
sign_correction <- sign(sum(eta_mean * eta_g.sim))
eta.post <- eta.post * sign_correction
rho.post[, -1] <- rho.post[, -1] * sign_correction

# ------------------------- Begin PDF -------------------------
pdf("code/checks/plot/1check.gibbs.sampler.pdf", width = 11, height = 8.5)
par_default <- par(no.readonly = TRUE)

# ------------------------ TRACE PLOTS ------------------------
layout(matrix(1:16, ncol = 4, byrow = TRUE))

# sigma^2 trace
# plot(sigma.post, type = "l", main = expression(Trace~plot~of~sigma^2),
#      xlab = "Iteration", ylab = "")
# abline(h = sigma2.sim, col = "red", lwd = 2)

# rho trace
for (j in 1:ncol(rho.post)) {
  plot(rho.post[, j], type = "l",
       main = bquote(Trace~plot~of~rho[.(j - 1)]),
       xlab = "Iteration", ylab = "")
  abline(h = rho_g.sim[j], col = "red", lwd = 2)
}

# mean rho trace
# rho_mean <- rowMeans(rho.post)
# plot(rho_mean, type = "l",
#      main = expression("Trace plot of mean of " * rho),
#      xlab = "Iteration", ylab = expression(bar(rho)))
# abline(h = mean(rho_g.sim), col = "red", lwd = 2)

# eta trace
for (j in 1:ncol(eta.post)) {
  plot(eta.post[, j], type = "l",
       main = bquote(Trace~plot~of~eta[.(j)]),
       xlab = "Iteration", ylab = "")
  abline(h = eta_g.sim[j], col = "red", lwd = 2)
}

# mean eta trace
# eta_mean <- rowMeans(eta.post)
# plot(eta_mean, type = "l",
#      main = expression("Trace plot of mean of " * eta),
#      xlab = "Iteration", ylab = expression(bar(eta)))
# abline(h = mean(eta_g.sim), col = "red", lwd = 2)

# ------------------------ ACF PLOTS ------------------------
layout(matrix(1:16, ncol = 4, byrow = TRUE))

# sigma^2 ACF
# acf(sigma.post, main = expression(ACF~of~sigma^2))

# rho ACF
for (j in 1:ncol(rho.post)) {
  acf(rho.post[, j], main = bquote("ACF of " ~ rho[.(j - 1)]))
}

# eta ACF
for (j in 1:ncol(eta.post)) {
  acf(eta.post[, j], main = bquote("ACF of " ~ eta[.(j)]))
}

# ------------------------ DENSITY PLOTS ------------------------
layout(matrix(1:16, ncol = 4, byrow = TRUE))

# sigma^2 density
# plot(density(sigma.post), main = expression(Density~of~sigma^2), xlab = "")
# abline(v = sigma2.sim, col = "red", lwd = 2)

# rho density
for (j in 1:ncol(rho.post)) {
  plot(density(rho.post[, j]), main = bquote(Density~of~rho[.(j - 1)]), xlab = "")
  abline(v = rho_g.sim[j], col = "red", lwd = 2)
}

# eta density
for (j in 1:ncol(eta.post)) {
  plot(density(eta.post[, j]), main = bquote(Density~of~eta[.(j)]), xlab = "")
  abline(v = eta_g.sim[j], col = "red", lwd = 2)
}

# ------------------------ Close Device ------------------------
par(par_default)
dev.off()
