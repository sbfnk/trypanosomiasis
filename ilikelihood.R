ilikelihood <- function(prev, mu, gamma, M, N)
  {
    sum(M*log(prev) +
        (N-M)*log(1-prev))
  }
