pfm <- function(pars, gamma, mu, mixing_structure = NA, density = FALSE, N = NA)
  {
    beta <- mixing(pars, mixing_structure)

    solvfun <- function(prev) {
      if (density) {
        lambda <- beta %*% (prev*N)
      } else {
        lambda <- beta %*% prev
      }
      
      r <- 1/(1+(gamma+mu)/lambda) - prev
      as.vector(r)
    }

    prev0 <- rep(0.5,nrow(beta))

    ans <- dfsane(par = prev0, fn = solvfun, control = list(trace = FALSE))

    ans$par
  }
