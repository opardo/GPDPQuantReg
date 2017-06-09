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

zero_function <- function(x) return(0 * x)