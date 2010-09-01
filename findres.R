findres <- function(NGM, projection, depth = 1, start = 1)
  {
    res <- FALSE
    if (depth == 1) {
      for (i in c(start:ncol(NGM))) {
        unit=diag(1,38,38)
        proj <- projection
        proj[i,i] = 1
        U <- proj %*% NGM
        Q <- (unit-proj) %*% NGM
        u_sr <- max(abs(eigen(U)$values))
        q_sr <- max(abs(eigen(Q)$values))
        if (u_sr > 1 && q_sr < 1) {
          species <- c()
          for (i in c(1:ncol(NGM))) {
            if (proj[i,i] == 1) {
              species <- c(species, i)
            }
          }
          cat ("Maintenance: ", species, "\n")
          res <- TRUE
        }
      }
    } else {
      for (i in c(start:(ncol(NGM)-depth+1))) {
        proj <- projection
        proj[i,i] = 1
        res <- findres(NGM, proj, depth - 1, i + 1)
      }
    }
    res
  }
