components <- function(NGM)
  {
    humans <- c(1)
    domestic <- c(2:5)
    wild <- c(6:38)

    proj_labels <- c("humans only", "domestic only", "wild only",
                     "humans and domestic", "domestic and wild",
                     "humans and wild", "all")

    projections <- matrix(0, nrow=(3+3+1),
                          ncol=length(humans)+length(domestic)+length(wild))
    
    projections[1,humans] <- 1
    projections[1,domestic] <- 0
    projections[1,wild] <- 0

    projections[2,humans] <- 0
    projections[2,domestic] <- 1
    projections[2,wild] <- 0

    projections[3,humans] <- 0
    projections[3,domestic] <- 0
    projections[3,wild] <- 1

    projections[4,humans] <- 1
    projections[4,domestic] <- 1
    projections[4,wild] <- 0

    projections[5,humans] <- 0
    projections[5,domestic] <- 1
    projections[5,wild] <- 1

    projections[6,humans] <- 1
    projections[6,domestic] <- 0
    projections[6,wild] <- 1

    projections[7,humans] <- 1
    projections[7,domestic] <- 1
    projections[7,wild] <- 1

    unit=diag(1,nrow(NGM),ncol(NGM))

    for (i in 1:nrow(projections)) {
      proj <- matrix(0,nrow(NGM), ncol(NGM))
      for (j in 1:ncol(projections)) {
        proj[j,j] <- projections[i,j]
      }

      U <- proj %*% NGM
      Q <- (unit-proj) %*% NGM
      u_sr <- max(abs(eigen(U)$values))
      q_sr <- max(abs(eigen(Q)$values))
      cat (proj_labels[i], ": rho(u) =", u_sr, ", rho(q) =", q_sr, "\n")
    }
  }

