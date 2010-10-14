mixing <- function(pars, mixing_structure=NA)
  {
    if (is.na(mixing_structure)) {
      # random
      
      beta <- matrix(NA, length(pars), length(pars))
      beta <- beta[pars>0,][,pars>0]

      b <- pars[pars > 0]
      
      for (i in 1:nrow(beta)) {
        for (j in 1:ncol(beta)) {
          beta[i,j]=b[i]*b[j]
        }
      }
    } else {
      #structured
      beta <- matrix(NA, nrow(mixing_structure), ncol(mixing_structure))

      for (i in 1:nrow(mixing_structure)) {
        for (j in 1:nrow(mixing_structure)) {
	  if (mixing_structure[i,j] > 0) {
            beta[i,j]=pars[mixing_structure[i,j]]
	  } else {
            beta[i,j]=0
	  }
        }
      }
    }
    
    beta

  }
