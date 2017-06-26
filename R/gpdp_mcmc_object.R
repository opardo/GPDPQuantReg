build_GPDP_MCMC <- function(
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
){

  GPDP_MCMC <- structure(list(), class = "GPDP_MCMC")

  GPDP_MCMC$metadata <- list()
  GPDP_MCMC$metadata$formula <- formula
  GPDP_MCMC$metadata$data <- data
  GPDP_MCMC$metadata$p <- p
  GPDP_MCMC$metadata$mcit <- mcit
  GPDP_MCMC$metadata$burn <- burn
  GPDP_MCMC$metadata$thin <- thin

  GPDP_MCMC$a_priori <- list()
  GPDP_MCMC$a_priori$c_DP <- c_DP
  GPDP_MCMC$a_priori$d_DP <- d_DP
  GPDP_MCMC$a_priori$c_lambda <- c_lambda
  GPDP_MCMC$a_priori$d_lambda <- d_lambda
  GPDP_MCMC$a_priori$alpha <- alpha
  GPDP_MCMC$a_priori$M <- M

  GPDP_MCMC$parameters <- list()
  GPDP_MCMC$parameters$sigmas <- list()
  GPDP_MCMC$parameters$pis <- list()
  GPDP_MCMC$parameters$zetas <- list()
  GPDP_MCMC$parameters$N <- list()
  GPDP_MCMC$parameters$b <- list()
  GPDP_MCMC$parameters$lambda <- list()
  GPDP_MCMC$parameters$f <- list()
  GPDP_MCMC$parameters$alternative <- list()

  return(GPDP_MCMC)
}


#Dispatchers. Please register your functions here.
dispatcher_creator <- function(function_name){
  return(
    function (x, ...) {
      UseMethod(function_name, x)
    }
  )
}

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
#' @export
diagnose <- dispatcher_creator("diagnose")
