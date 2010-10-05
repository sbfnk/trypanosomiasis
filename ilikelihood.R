ilikelihood <- function(prev, mu, gamma, M, N)
  {
    sum(M[prev>0]*log(prev[prev>0]) +
        (N[prev>0]-M[prev>0])*log(1-prev[prev>0]))
  }
