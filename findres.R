findres <- function(NGM, projection = matrix(0,nrow(NGM), ncol(NGM)), depth = 1, start = 1)
  {
    res <- FALSE
    if (depth == 1) {
      cat ("R0 = ", max(abs(eigen(NGM)$values)), "\n")
      for (i in c(start:ncol(NGM))) {
        unit=diag(1,nrow(NGM),ncol(NGM))
        proj <- projection
        proj[i,i] = 1
        U <- proj %*% NGM
        Q <- (unit-proj) %*% NGM
        u_sr <- max(abs(eigen(U)$values))
        q_sr <- max(abs(eigen(Q)$values))
        pvec <- c()
        for (j in 1:nrow(proj))  {
          if (proj[j,j] == 1) {
            pvec <- c(pvec, j)
          }
        }
        cat ("pvec = ", pvec, ", usr = ", u_sr, ", q_sr = ", q_sr, "\n")
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
        if (findres(NGM, proj, depth - 1, i + 1)) {
          res <- TRUE
        }
      }
    }
    res
  }
