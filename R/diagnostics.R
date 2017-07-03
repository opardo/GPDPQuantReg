#' Diagnose GPDPQuantReg
#'
#' Run diagnostics for the GPDPQuantReg's Monte Carlo Markov Chains
#'
#' @param GPDP_MCMC A GPDP_MCMC object
#' @return plots traces, autocorrelations, crosscorrelations and ergodicity
#' @author Omar Pardo (omarpardog@gmail.com)
#' @examples
#' m <- 35
#' x <- sort(sample(seq(-15, 15, 0.005), m))
#' f_x <- function(x) return((1/40) * x^2 - (1/20) * x - 2)
#' data <- data.frame(x = x, y = f_x(x) + rnorm(m, 0, 1))
#' GPDP_MCMC <- GPDPQuantReg(y ~ x, data, p = 0.250)
#' diagnose(GPDP_MCMC)
#' @export diagnose.GPDP_MCMC
diagnose.GPDP_MCMC <- function(GPDP_MCMC) {
  rand_obs <- sample(size = 1, x = 1:length(GPDP_MCMC$parameters$f[[1]]))

  zeta <- matrix(
    unlist(GPDP_MCMC$parameters$zetas),
    nrow = length(GPDP_MCMC$parameters$zetas),
    byrow = T
  )[ , rand_obs]

  sigma <- unlist(lapply(
    1: length(zeta),
    function(i) GPDP_MCMC$parameters$sigmas[[i]][zeta[i]]
  ))

  lambda <- unlist(GPDP_MCMC$parameters$lambda)

  f <- matrix(
    unlist(GPDP_MCMC$parameters$f),
    nrow = length(GPDP_MCMC$parameters$f),
    byrow = T
  )[, rand_obs]

  MCMC_test <- mcmc(cbind(
    sigma,
    lambda,
    f
  ))

  diagnostic_trace(MCMC_test)
  diagnostic_autocorrelation(MCMC_test)
  diagnostic_crosscorrelation(MCMC_test)
  diagnostic_ergodicity(MCMC_test)
}

diagnostic_autocorrelation <- function(MCMC_test) {
  return(autocorr.plot(MCMC_test))
}

diagnostic_trace <- function(MCMC_test) {
  return(traceplot(MCMC_test))
}

diagnostic_ergodicity <- function(MCMC_test) {
  it_means <- sweep(apply(MCMC_test, 2, cumsum), 1, 1:nrow(MCMC_test), "/")
  it_means <- it_means[, which(colnames(it_means) %in% c("sigma", "N", "lambda", "f"))]
  melted_it_means <- melt(it_means)
  colnames(melted_it_means) <- c("iterations", "parameter", "value")
  plot(xyplot(
    value ~ iterations | parameter,
    type = "l",
    data = melted_it_means,
    scales = list(y = list(relation="free")),
    main="Ergodicity"
  ))
  # return(ergodicity_plot)
}

diagnostic_crosscorrelation <- function(MCMC_test) {
  return(crosscorr.plot(MCMC_test))
}
