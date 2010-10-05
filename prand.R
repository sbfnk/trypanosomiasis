prand <- function(mixing=NA, pars, gamma, mu, density = FALSE, N = NA)
  {
    beta <- matrix(NA, length(pars), length(pars))

    solvfun <- function(prev) {
      for (i in 1:length(pars)) {
        for (j in 1:length(pars)) {
          beta[i,j]=pars[i]*pars[j]
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

    prev0 <- rep(0.5,length(pars))

    ans <- dfsane(par = prev0, fn = solvfun, control = list(trace = FALSE))
    ans$par
  }
