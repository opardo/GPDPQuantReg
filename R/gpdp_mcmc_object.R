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

diagnose <- dispatcher_creator("diagnose")
predict_GPDPQuantReg <- dispatcher_creator("predict_GPDPQuantReg")
# GPDPQuantReg <- dispatcher_creator("GPDPQuantReg")
