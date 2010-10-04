pfm <- function(mixing, b, gamma, mu, density = FALSE, N = NA)
  {
    beta <- matrix(NA, nrow(mixing), ncol(mixing))

    solvfun <- function(prev) {
      for (i in 1:nrow(mixing)) {
        for (j in 1:nrow(mixing)) {
          beta[i,j]=b[mixing[i,j]]
        }
      }
      if (density) {
        lambda <- beta %*% (prev*N)
      } else {
        lambda <- beta %*% prev
      }
      
      r <- 1/(1+(gamma+mu)/lambda) - prev
      as.vector(r)
    }

    prev0 <- rep(0.5,nrow(mixing))

    ans <- dfsane(par = prev0, fn = solvfun, control = list(trace = FALSE))

    ans$par
  }
