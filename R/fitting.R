#' Fits GPDPQuantReg
#'
#' @description Runs a MCMC algorithm to fit a quantile regression model, using Gaussian and
#' Dirichlet Processes. Returns a GPDP_MCMC object with thinned chains.
#' Scales data in the process so the induced correlation makes sense, and
#' unscales at the end. Take it into account for the initial parameters.
#'
#' @usage GPDPQuantReg(
#'   formula, data, p = 0.5, c_DP = 2, d_DP = 1, c_lambda = 2,
#'   d_lambda = 0.5, alpha = sqrt(nrow(data)), M = zero_function,
#'   mcit = 3e4, burn = 1e4, thin = 10
#' )
#' @param formula formula object with the dependent and independent variables
#' @param data data frame
#' @param p real between (0,1), probability corresponding to the estimated quantile
#' @param c_DP real > 0, shape parameter for the DP's base distribution
#' @param d_DP real > 0, scale parameter fot the DP's base distribution
#' @param c_lambda real > 0, shape parameter for the GP's lambda distribution
#' @param d_lambda real > 0, scale parameter fot the GP's lambda distribution
#' @param alpha real > 0, DP's concentration parameter
#' @param M function, a priori estimation for the final function
#' @param mcit integer > 0, number of MCMC algorithm's valid chains
#' @param burn integer > 0, number of MCMC algorithm's first burned chains
#' @param thin integer > 0 and < mcit, MCMC algorithm's thinning
#' @return GPDP_MCMC object, with the MCMC algorithm's chains
#' @author Omar Pardo (omarpardog@gmail.com)
#' @examples
#' m <- 35
#' x <- sort(sample(seq(-15, 15, 0.005), m))
#' f_x <- function(x) return((1/40) * x^2 - (1/20) * x - 2)
#' data <- data.frame(x = x, y = f_x(x) + rnorm(m, 0, 1))
#' GPDP_MCMC <- GPDPQuantReg(y ~ x, data, p = 0.250)
#' @export
GPDPQuantReg <- function(
  formula,
  data,
  p = 0.5,
  c_DP = 2,
  d_DP = 1,
  c_lambda = 2,
  d_lambda = 0.5,
  alpha = sqrt(nrow(data)),
  M = zero_function,
  mcit = 3e4,
  burn = 1e4,
  thin = 10
){

  ptm <- proc.time()

  # Load X and Y from data
  formula <- delete_intercept(formula)
  mf <- model.frame(formula = formula, data = data)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- model.response(data = mf)

  # Save used variables in formula
  y_name <- toString(formula[2])
  x_names <- colnames(X)
  formula <- as.formula(paste0(y_name, " ~ ", paste(x_names, collapse= "+"), "+ 0"))

  # Build GPDP_MCMC S3 object
  GPDP_MCMC <- build_GPDP_MCMC(
    formula,
    data,
    p,
    mcit,
    burn,
    thin,
    c_DP,
    d_DP,
    c_lambda,
    d_lambda,
    alpha,
    M
  )

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
  cat(sprintf("Are total %6d iterations.\n", mcit + burn))
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
      cat(sprintf("iteration: %6d \n", i))
    }
  }

  GPDP_MCMC$metadata$time <- proc.time() - ptm

  return(GPDP_MCMC)
}
