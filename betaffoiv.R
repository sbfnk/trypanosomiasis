betaffoiv <- function(rlambda, vdensity, rabundance, factor = 1,
                      vprev = NULL, rprev = NULL, vmu = NULL)
  {
    foi <- function(beta) {
      r <- beta * factor
      if (is.null(vprev)) {
        # nonlinear equation
        r <- r * vdensity / rabundance * sum(beta * factor * rprev) /
          sum(vmu + beta * factor * rprev)
      } else {
        # linear equation
        r <- r * vdensity / rabundance * vprev
      }
      r <- r - rlambda
      as.vector(r)
    }

    b0 <- rep(1,length(rlambda))

    ans <- dfsane(par = b0, fn = foi, control = list(trace = FALSE))

    ans
  }
