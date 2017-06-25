predict_GPDP <- function(GPDP_MCMC, Xp, credibility) {
  
  # Load needed metadata
  X <- GPDP_MCMC$metadata$X
  Y <- GPDP_MCMC$metadata$Y
  MCMC_length <- GPDP_MCMC$metadata$mcit / GPDP_MCMC$metadata$thin
  
  # Get scale parameters
  X_scaled_means <- as.numeric(attr(scale(X), "scaled:center"))
  X_scaled_sigmas <- as.numeric(attr(scale(X), "scaled:scale"))
  y_scaled_mean <- attr(scale(Y), "scaled:center")
  y_scaled_sigma <- attr(scale(Y), "scaled:scale")
  
  # Get scaled f & lambda posteriors from the fitted model
  fs <- lapply(
    1:MCMC_length,
    function(i) (GPDP_MCMC$parameters$f[[i]] - y_scaled_mean) /  y_scaled_sigma
  )
  lambdas <- unlist(GPDP_MCMC$parameters$lambda) / y_scaled_sigma
  fitted_parameters <- lapply(
    1:MCMC_length,
    function(i) list(f = fs[[i]], lambda = lambdas[[i]])
  )
  
  # Scale predictive data
  scaled_Xp <- matrix(nrow = nrow(Xp), ncol = ncol(Xp))
  for(i in 1:ncol(Xp)) {
    scaled_Xp[,i] <- (Xp[,i] - X_scaled_means[i]) / X_scaled_sigmas[i]
  }
  
  # Get means and correlations for the Gaussian Process
  # original data
  M_X <- M(X)
  K_XX_inv <- solve(get_K_XX(X))
  # predictive data
  M_Xp <- M(scaled_Xp)
  K_XpXp <- get_K_XX(scaled_Xp)
  # both
  K_XpX <- get_K_XpX(scaled_Xp, X)
  
  # Calculate fixed conditional parameters
  Sigma <- K_XpXp - K_XpX %*% K_XX_inv %*% t(K_XpX)
  Mu_aux <- K_XpX %*% K_XX_inv
  
  # Calculate predictive data's posterior
  fp_list <- lapply(
    parameters,
    FUN = simulate_fp,
    Sigma = Sigma, 
    M_Xp = M_Xp, 
    M_X = M_X, 
    Mu_aux = Mu_aux
  )
  
  fp_matrix <- data.frame(matrix(
    unlist(fp_list),
    nrow = length(fp_list),
    byrow = T
  ))
  
  # Calculate mean and median for each x
  fp_mean <- unname(apply(fp_matrix, 2, mean))
  fp_median <- unname(apply(fp_matrix, 2, get_quantile, p = 0.5))
  
  # Calculate credible interval
  p_lower <- (1 - credibility) / 2
  p_upper <- 1 - p_lower
  fp_lower <- unname(apply(fp_matrix, 2, get_quantile, p_lower))
  fp_upper <- unname(apply(fp_matrix, 2, get_quantile, p_upper))
  
  # Organize prediction
  results <- data.frame(
    mean = fp_mean,
    median = fp_median,
    lower = fp_lower,
    upper = fp_upper
  )
  prediction <- cbind(Xp, results)
  
  return(prediction)
}

simulate_fp <- function(params, Sigma, M_Xp, M_X, Mu_aux) {
  f <- params$f
  lambda <- params$lambda
  Mu <- M_Xf + Mu_aux %*% (f - M_X)
  scaled_fp <- as.vector(rmvnorm(1, mean = Mu, sigma = lambda * Sigma))
  fp <- (scaled_fp * y_scaled_sigma) + y_scaled_mean
  return(ff)
}

