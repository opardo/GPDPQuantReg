# GAUSSIAN PROCESS FUNCTIONS

update_b <- function(m, p, sigmas, zetas) {
  return(ifelse(runif(m) > p, p / sigmas[zetas], -(1-p) / sigmas[zetas]))
}

update_f <- function(Y, M_X, K_XX, lambda, b) {
  f <- tmvtnorm::rtmvnorm(
    n = 1,
    mean = as.vector(M_X + (lambda * K_XX) %*% b),
    sigma = lambda * K_XX,
    lower = ifelse(b >= 0, -Inf, Y),
    upper = ifelse(b >= 0, Y, Inf),
    algorithm = "gibbs"
  )[1,]
  if(any(is.nan(f)) > 0 | any(is.infinite(f)) > 0) {
    f <- NULL
  }
  return(f)
}

update_lambda <- function(f, M_X, b, K_XX, K_XX_inv, n, c_lambda, d_lambda) {

  bar_F <- 1/2 * as.numeric(t(f-as.vector(M_X)) %*% K_XX_inv %*% (f-as.vector(M_X)))
  bar_B <- 1/2 * as.numeric(t(b) %*% K_XX %*% b)
  bar_c <- c_lambda + n/2
  bar_d <- d_lambda + bar_F

  lambdas <- seq(0.005, 5, 0.005)
  fx <- function(lambda) {
    return(lambda ^ (-(bar_c + 1)) * exp(- bar_d / lambda - bar_B * lambda))
  }
  prob <- fx(lambdas)

  if (sum(prob > 0) > 0) {
    lambda <- sample(x = lambdas[prob > 0], size = 1, prob = prob[prob > 0])
  } else {
    lambda <- pscl::rigamma(1, c_lambda, d_lambda)
  }
  return(lambda)
}