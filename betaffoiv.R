betaffoiv <- function(rlambda, vdensity, rabundance, theta,
                      biting_rate, area_convert,
                      vprev = NULL, rprev = NULL, vmu = NULL)
  {
    ## cat ("rlambda: ", rlambda, "\n")
    ## cat ("vdensity: ", vdensity, "\n")
    ## cat ("rabundance: ", rabundance, "\n")
    ## cat ("theta: ", theta, "\n")
    ## cat ("biting_rate: ", biting_rate, "\n")
    ## cat ("area_convert: ", area_convert, "\n")
    ## cat ("vprev: ", vprev, "\n")
    ## cat ("rprev: ", rprev, "\n")
    ## cat ("vmu: ", vmu, "\n")

    foi <- function(beta) {
      r <- beta
      if (is.null(vprev)) {
        # nonlinear equation
        r <- r * area_convert * vdensity / rabundance *
          biting_rate*sum(beta * theta * rprev) /
          (vmu + biting_rate * sum(beta * theta * rprev))
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
