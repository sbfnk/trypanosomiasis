foifm <- function(mixing, b, gamma, mu)
  {
    beta <- matrix(NA, nrow(mixing), ncol(mixing))
    for (i in 1:nrow(mixing)) {
      for (j in 1:nrow(mixing)) {
        beta[i,j]=b[mixing[i,j]]
      }
    }

    foi <- function(l) {
      r <- (l - (beta %*% (l/(l+mu+gamma))))
      as.vector(r)
    }

    l0 <- rep(1,nrow(mixing))

    ans <- dfsane(par = l0, fn = foi, control = list(trace = FALSE))

    ans$par
  }
