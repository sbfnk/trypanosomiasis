foifm <- function(mixing, b, gamma, mu)
  {
    beta <- matrix(nrow == nrow(mixing), ncol == nrow(mixing))
    for (i in 1:nrow(mixing)) {
      for (j in 1:nrow(mixing)) {
        beta[i,j]=b[mixing[i,j]]
      }
    }

    foi <- function(l) {
      r <- (l - (beta %*% (l/(l+mu+gamma))))
      as.vector(r)
    }

    
  }
