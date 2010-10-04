prand <- function(a, gamma, mu, density = FALSE, N = NA)
  {
    beta <- matrix(NA, length(a), length(a))

    solvfun <- function(prev) {
      for (i in 1:length(a)) {
        for (j in 1:length(a)) {
          beta[i,j]=a[i]*a[j]
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

    prev0 <- rep(0.5,length(a))

    ans <- dfsane(par = prev0, fn = solvfun, control = list(trace = FALSE))
    ans$par
  }
