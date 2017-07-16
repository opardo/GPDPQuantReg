get_quant_p_from_Y_function <- function(formula, data, p){
  mf <- model.frame(formula = formula, data = data)
  Y <- scale(model.response(data = mf))
  Q <- as.numeric(quantile(Y, p))
  quant_p <- function(X) rep(Q, nrow(X))
  return(quant_p)
}

default_DP <- function(p){
  return(p * (1-p) / sqrt(2 * (1 - 2 * p * (1-p))))
}