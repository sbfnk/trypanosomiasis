mffoi <- function(lambda, mixing, gamma, mu)
  {
    beta <- matrix(NA, nrow(mixing), ncol(mixing))

    foi <- function(b) {
      for (i in 1:nrow(mixing)) {
        for (j in 1:nrow(mixing)) {
          beta[i,j]=b[mixing[i,j]]
        }
      }
      r <- (lambda - (beta %*% (lambda/(lambda+mu+gamma))))
      as.vector(r)
    }

    b0 <- rep(1,max(mixing))

    ans <- dfsane(par = b0, fn = foi, control = list(trace = FALSE))

    for (i in 1:nrow(mixing)) {
      for (j in 1:nrow(mixing)) {
        beta[i,j]=ans$par[mixing[i,j]]
      }
    }

    beta
  }
