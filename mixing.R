mixing <- function(pars, mixing_structure=NA)
  {
    if (is.na(mixing_structure)) {
      # random
      beta <- matrix(NA, length(pars), length(pars))
      
      for (i in 1:length(pars)) {
        for (j in 1:length(pars)) {
          beta[i,j]=pars[i]*pars[j]
        }
      }
    } else {
      #structured
      beta <- matrix(NA, nrow(mixing_structure), ncol(mixing_structure))

      for (i in 1:nrow(mixing_structure)) {
        for (j in 1:nrow(mixing_structure)) {
          beta[i,j]=pars[mixing_structure[i,j]]
        }
      }
    }
    
    beta

  }
