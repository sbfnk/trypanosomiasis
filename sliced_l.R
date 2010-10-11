sliced_l <- function(pars, gamma, mu, M, mixing_structure = NA,
                       density = FALSE, N = rep(1, length(M)),
                       param=1, width=0.2)
{
  library('BB')
  source('pfm.R')
  source('mixing.R')
  source('ilikelihood.R')
  
  prev <- pfm(b, gamma=gamma, mu=mu, mixing_structure=mixing_structure,
              density=density, N=N)
  slices = seq(pars[param]*(1-width/2), pars[param]*(1+width/2),
              pars[param]*width/10)
  slice_matrix <- matrix(rep(b, length(slices)), ncol=length(slices))
  slice_matrix[param,] <- slices
  l_slices = rep(0, length(slices))
  for (i in 1:length(slices)) {
    l_slices[i] <- ilikelihood(pfm(slice_matrix[,i], gamma=gamma, mu=mu,
                                   mixing_structure=mixing_structure,
                                   density=density, N=N),
                               mu=mu, gamma=gamma, M=M, N=N)
  }
  return(list(x = slices, y = l_slices))
}
