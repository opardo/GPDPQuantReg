# CLASSES FUNCTIONS

update_class_i <- function(i, xis, zetas, pis, eps, sigmas, N, p) {
  u <- runif(n = 1, min = 0, max = xis[zetas[i]])
  if(sum(is.na((xis[1:N] > u) * pis / xis[1:N] * dalp(eps[i], sigma = sigmas, p = p))) > 0) browser()
  zeta_i <- sample(
    1,
    x = 1:N,
    prob = (xis[1:N] > u) *
      pis / xis[1:N] *
      dalp(eps[i], sigma = sigmas, p = p)
  )
  N_i <- max((1:length(xis))[xis > u])
  return(list(
    zeta_i = zeta_i,
    N_i = N_i
  ))
}

update_classes <- function(xis, zetas, pis, eps, sigmas, N, p, m) {
  return(lapply(
    1:m,
    update_class_i,
    xis = xis,
    zetas = zetas,
    pis = pis,
    eps = eps,
    sigmas = sigmas,
    N = N,
    p = p
  ))
}

update_zetas <- function(classes) unlist(
  lapply(classes, function(class) class$zeta_i)
)

update_N <- function(classes) max(unlist(
  lapply(classes, function(class) class$N_i)
))
