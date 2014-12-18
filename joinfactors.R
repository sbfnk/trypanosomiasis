joinfactors <- function(theta = NULL, b = NULL, rabundance = NULL,
                        vdensity = NULL, area_convert = NULL)
  {
    factor <- 1
    if (!is.null(theta)) {
      factor <- factor * theta
    }
    if (!is.null(b)) {
      factor <- factor * b
    }
    if (!is.null(area_convert)) {
      factor <- factor * area_convert
    }
    factor
  }
