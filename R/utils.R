# Correlation matrix
get_K_XX <- function(X) {
  m <- dim(X)[1]
  K <- matrix(nrow = m, ncol = m)
  for(i in 1:m) {
    for(j in 1:m){
      K[i,j] <- exp(-as.numeric(dist(rbind(X[i,],X[j,]))))
    }
  }
  return(K)
}

# Correlation between predictive and original data
get_K_XpX <- function(Xp, X) {
  r <- dim(Xp)[1]
  m <- dim(X)[1]
  K <- matrix(nrow = r, ncol = m)
  for(i in 1:r) {
    for(j in 1:m){
      K[i,j] <- exp(-as.numeric(dist(rbind(Xp[i,], X[j,]))))
    }
  }
  return(K)
}

get_quantile <- function(values, p) as.numeric(quantile(values, p))

zero_function <- function(x) return(0 * x)

delete_intercept <- function(formula) {
  non_intercept_in_formula <- grepl("+0", gsub(" ", "", toString(formula[3])))
  if(!non_intercept_in_formula) {
    formula <- as.formula(paste(
      toString(formula[2]),
      "~",
      toString(formula[3]),
      "+0"
    ))  
  }
  return(formula)
}