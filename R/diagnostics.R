diagnostic_GPDPQuantReg <- function(GPDP_MCMC) {
  rand_obs <- sample(size = 1, x = 1:length(GPDP_MCMC$parameters$f[[1]]))
  
  zetas <- matrix(
    unlist(GPDP_MCMC$parameters$zetas),
    nrow = length(GPDP_MCMC$parameters$zetas),
    byrow = T
  )
  zeta <- zetas[ , rand_obs]
  
  sigma <- unlist(lapply(
    1: length(zeta), 
    function(i) GPDP_MCMC$parameters$sigmas[[i]][zeta[i]]
  ))
  
  N <- unlist(GPDP_MCMC$parameters$N)
  
  lambda <- unlist(GPDP_MCMC$parameters$lambda)
  
  fs <- matrix(
    unlist(GPDP_MCMC$parameters$f),
    nrow = length(GPDP_MCMC$parameters$f),
    byrow = T
  )
  f <- fs[, rand_obs]
  
  MCMC_test <- mcmc(cbind(
    sigma,
    zeta,
    N,
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
    main="Ergodicity", 
  ))
  # return(ergodicity_plot)
}

diagnostic_crosscorrelation <- function(MCMC_test) {
  return(crosscorr.plot(MCMC_test))
}
