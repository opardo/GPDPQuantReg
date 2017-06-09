# ASYMMETRIC LAPLACE FUNCTIONS

# Density
dalp <- function(x, sigma, p) {
  kernel <- ifelse (
    x >= 0,
    exp(- p * x / sigma),
    exp((1 - p) * x / sigma)
  )
  d <- (p * (1 - p) / sigma) * kernel
  return(d)
}

# Random one sigma
ralp <- function(n, sigma = 1, p = 0.5) {
  x <- seq(-10, 10, 0.05)
  prob <- dalp(x, sigma, p)
  sample_alp <- sample(
    x = x,
    size = n,
    prob = prob
  )
  return(sample_alp)
}

# Random different sigmas vector
random_asymmetric_error <- function(sigmas, p) {
  unique_sigmas <- unique(sigmas)
  error <- vector(length = length(sigmas))

  for (sigma in unique_sigmas) {
    dim <- sum(sigmas == sigma)
    simulations <- ralp(dim, sigma = sigma, p)
    error[sigmas == sigma] <- simulations
  }

  return(error)
}
