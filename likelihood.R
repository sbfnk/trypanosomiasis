likelihood <- function(lambda, mu, gamma, M, N)
  {
    sum(M*log(lambda/(lambda+mu+gamma)) +
        (N-M)*log((mu+gamma)/(lambda+mu+gamma)))
  }
