findres <- function(NGM, projection = matrix(0,nrow(NGM), ncol(NGM)), depth = 1, start = 1)
  {
    res <- FALSE
    species_list <- c()
    u <- c()
    q <- c()
    if (depth == 1) {
      R0 <- max(abs(eigen(NGM)$values))
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
          species_list <- c(species_list, species)
          u <- c(u, u_sr)
          q <- c(q, q_sr)
          res <- TRUE
        }
      }
    } else {
      for (i in c(start:(ncol(NGM)-depth+1))) {
        proj <- projection
        proj[i,i] = 1
        search_result <- findres(NGM, proj, depth - 1, i + 1)
        R0 <- search_result$R0
        if (search_result$res) {
          species_list <- c(species_list, as.vector(search_result$combos))
          u <- c(u, search_result$u)
          q <- c(q, search_result$q)
          res <- TRUE
        }
      }
    }
    if (length(species_list) > 0) {
      combos <- matrix(species_list, ncol=depth, byrow=T)
    } else {
      combos <- NA
    }
    return(list(res = res, combos = combos, u = u, q = q, R0 = R0))
  }
