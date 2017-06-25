fit_GPDPQuantReg <- function(
  formula,
  data,
  p = 0.5,
  c_DP = 2,
  d_DP = 1,
  c_lambda = 2,
  d_lambda = 0.5,
  alpha = sqrt(length(Y)),
  M = zero_function,
  mcit = 3e4,
  burn = 1e4,
  thin = 10
){
  
  # Load X and Y from data
  formula <- update(formula, ~ . + 0)
  mf <- model.frame(formula = formula, data = data)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- model.response(data = mf)
  
  # Save used variables in formula
  y_name <- toString(formula[2])
  x_names <- colnames(X)
  formula <- as.formula(paste0(y_name, " ~ ", paste(x_names, collapse= "+"), "+ 0"))
  
  # Repositories
  GPDP_MCMC <- list()
  
  GPDP_MCMC$parameters <- list()
  GPDP_MCMC$parameters$sigmas <- list()
  GPDP_MCMC$parameters$pis <- list()
  GPDP_MCMC$parameters$zetas <- list()
  GPDP_MCMC$parameters$N <- list()
  GPDP_MCMC$parameters$b <- list()
  GPDP_MCMC$parameters$lambda <- list()
  GPDP_MCMC$parameters$f <- list()
  GPDP_MCMC$parameters$alternative <- list()
  
  GPDP_MCMC$a_priori <- list()
  GPDP_MCMC$a_priori$c_DP <- c_DP 
  GPDP_MCMC$a_priori$d_DP <- d_DP
  GPDP_MCMC$a_priori$c_lambda <- c_lambda
  GPDP_MCMC$a_priori$d_lambda <- d_lambda
  GPDP_MCMC$a_priori$alpha <- alpha
  GPDP_MCMC$a_priori$M <- M
  
  GPDP_MCMC$metadata <- list()
  GPDP_MCMC$metadata$formula <- formula
  GPDP_MCMC$metadata$data <- data
  GPDP_MCMC$metadata$p <- p
  GPDP_MCMC$metadata$mcit <- mcit
  GPDP_MCMC$metadata$burn <- burn
  GPDP_MCMC$metadata$thin <- thin
  
  cont <- 1
  
  # Scale data
  y_scaled_mean <- attr(scale(Y), "scaled:center")
  y_scaled_sigma <- attr(scale(Y), "scaled:scale")
  
  X <- as.matrix(scale(X))
  Y <- as.vector(scale(Y))

  # Fixed values
  xis <- (0.5)*(0.5)^(0:1e3)
  m <- dim(X)[1]
  n <- dim(X)[2]
  K_XX <- get_K_XX(X)
  K_XX_inv <- solve(K_XX)
  M_X <- M(X)

  # Initial values
  N <- 50
  zetas <- sample(1:N, size = length(Y), replace = T)
  f <- Y
  eps <- Y - f

  # Gibbs sampler
  cat(sprintf("There are %12d iterations.\n", mcit + burn))
  for(i in 1:(mcit+burn)){

    # Update Dirichlet Process
    DP_ks <- update_DP_ks(eps, zetas, p, c_DP, d_DP, alpha, N)
    sigmas <- update_sigmas(DP_ks)
    betas <- update_betas(DP_ks)
    pis <- turn_betas_into_pis(betas)

    # Update each observation's class
    classes <- update_classes(xis, zetas, pis, eps, sigmas, N, p, m)
    zetas <- update_zetas(classes)

    # Update Y >= f or Y < f
    b <- update_b(m, p, sigmas, zetas)

    # Update Gaussian Process
    try_f <- NULL
    global_chances <- 1

    while (global_chances <= 3 & is.null(try_f) ) {

      try_lambda <- update_lambda(f, M_X, b, K_XX, K_XX_inv, n, c_lambda, d_lambda)

      f_chances <- 1
      while(f_chances <= 3 & is.null(try_f)) {
        try_f <- update_f(Y, M_X, K_XX, try_lambda, b)
        f_chances <- f_chances + 1
      }

      global_chances <- global_chances + 1
    }

    if(is.null(try_f)) {
      f <- Y - random_asymmetric_error(sigmas[zetas], p)
      lambda <- update_lambda(f, M_X, b, K_XX, K_XX_inv, n, c_lambda, d_lambda)
      alternative <- 1
    } else {
      lambda <- try_lambda
      f <- try_f
      alternative <- 0
    }

    eps <- Y-f
    
    # Update number of classes truncation
    N <- update_N(classes)

    # Aux to delete burning simulations
    if(i > burn && (i - burn) %% thin == 0){
      GPDP_MCMC$parameters$sigmas[[cont]] <- sigmas * y_scaled_sigma
      GPDP_MCMC$parameters$pis[[cont]] <- pis
      GPDP_MCMC$parameters$zetas[[cont]] <- zetas
      GPDP_MCMC$parameters$N[[cont]] <- N
      GPDP_MCMC$parameters$b[[cont]] <- b
      GPDP_MCMC$parameters$lambda[[cont]] <- lambda * y_scaled_sigma
      GPDP_MCMC$parameters$f[[cont]] <- f * y_scaled_sigma + y_scaled_mean
      GPDP_MCMC$parameters$alternative[[cont]] <- alternative
      cont <- cont + 1
    }

    if(i %% 1000 == 0){
      cat(sprintf("%12d iterations.\n", i))
    }
  }

  return(GPDP_MCMC)
}