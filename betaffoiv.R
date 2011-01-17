betaffoiv <- function(rlambda, vlambda, theta = NULL, b = NULL,
                     rabundance = NULL, vdensity = NULL, area_convert = NULL)
  {
    foi <- function(beta) {
      factor <- beta
      if (!is.null(theta)) {
        factor <- factor * lambda
      }
      if (!is.null(b)) {
        factor <- factor * b
      }
      if (!is.null(rabundance)) {
        factor <- factor / b
      }
      if (!is.null(vdensity)) {
        factor <- factor * vdensity
      }
      if (!is.null(area_convert)) {
        factor <- factor * area_convert
      }
      r <- beta * (t(beta) %*% (prev*N)) /
        (mu_vector + (t(beta) %*% (prev*N)))
      as.vector(r)
    }

    b0 <- rep(1,length(theta))

    ans <- dfsane(par = b0, fn = foi, control = list(trace = FALSE))

    ans
  }
