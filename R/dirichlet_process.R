# DIRICHLET PROCESS FUNCTIONS

# Update k-sigma and k-beta
update_DP_k <- function(k, eps, zetas, p, c_DP, d_DP, alpha) {
  ind <- zetas == k
  ind_pos <- ind & (unname(eps) >= 0)
  ind_neg <- ind & (unname(eps) < 0)
  r <- sum(ind)

  sigma_k <- rigamma(
    1,
    c_DP + r,
    d_DP + p * sum(eps[ind_pos]) + (1-p) * sum(-eps[ind_neg])
  )

  beta_k <- rbeta(n = 1, shape1 = 1 + r, shape2 = alpha + sum(zetas > k))

  return(list(
    sigma_k = sigma_k,
    beta_k = beta_k
  ))
}


# Update Dirichlet Process
update_DP_ks <- function(eps, zetas, p, c_DP, d_DP, alpha, N) {
  return(lapply(
    1:N,
    update_DP_k,
    eps = eps,
    zetas = zetas,
    p = p,
    c_DP = c_DP,
    d_DP = d_DP,
    alpha = alpha
  ))
}

update_sigmas <- function(DP_ks) unlist(
  lapply(DP_ks, function(DP_k) DP_k$sigma_k)
)

update_betas <- function(DP_ks) unlist(
  lapply(DP_ks, function(DP_k) DP_k$beta_k)
)

turn_betas_into_pis <- function(betas) {
  pis_hat <- betas * c(1, cumprod(1-betas)[1:(length(betas)-1)])
  pis <- pis_hat / sum(pis_hat)
  return(pis)
}
