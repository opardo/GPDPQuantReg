# TOP LEVEL FUNCTION

MCMC_GPDPQuantReg <- function(
  X,
  Y,
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

  # Repositories
  output <- list()
  output$sigmas <- list()
  output$pis <- list()
  output$zetas <- list()
  output$N <- list()
  output$b <- list()
  output$lambda <- list()
  output$f <- list()
  output$bad <- list()
  cont <- 1

  # Scale data
  scaled_mean <- attr(scale(Y), "scaled:center")
  scaled_sigma <- attr(scale(Y), "scaled:scale")
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
  cat(sprintf("Son %6d iteraciones.\n", mcit + burn))
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
    global_chances <- 0

    while (global_chances < 3 & is.null(try_f) ) {

      try_lambda <- update_lambda(f, M_X, b, K_XX, K_XX_inv, n, c_lambda, d_lambda)

      f_chances <- 0
      while(f_chances < 4 & is.null(try_f)) {
        try_f <- update_f(Y, M_X, K_XX, try_lambda, b)
        f_chances <- f_chances + 1
      }

      global_chances <- global_chances + 1
    }

    if(is.null(try_f)) {
      f <- Y - random_asymmetric_error(sigmas[zetas], p)
      lambda <- update_lambda(f, M_X, b, K_XX, K_XX_inv, n, c_lambda, d_lambda)
      bad <- 1
    } else {
      lambda <- try_lambda
      f <- try_f
      bad <- 0
    }

    eps <- Y-f

    # Aux to delete burning simulations
    if(i > burn && (i - burn) %% thin == 0){
      output$sigmas[[cont]] <- sigmas
      output$pis[[cont]] <- pis
      output$zetas[[cont]] <- zetas
      output$N[[cont]] <- N
      output$b[[cont]] <- b
      output$lambda[[cont]] <- lambda
      output$f[[cont]] <- f
      output$bad[[cont]] <- bad
      cont <- cont + 1
    }

    if(i %% 1000 == 0){
      cat(sprintf("Van %6d iteraciones.\n", i))
    }

    # Update number of classes truncation
    N <- update_N(classes)
  }

  # Unscale parameters
  output$sigmas <- lapply(output$sigmas, function(sigmas) sigmas * scaled_sigma)
  output$f <- lapply(output$f, function(f) f * scaled_sigma + scaled_mean)
  output$lambda <- lapply(output$lambda, function(lambda) lambda * scaled_sigma)

  return(output)
}