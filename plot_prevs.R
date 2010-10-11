plot_prevs <- function(pars, gamma, mu, M, mixing_structure = NA,
                       density = FALSE, N = rep(1, length(M)),
                       labels = rep("", length(M)))
{
  library('BB')
  source('pfm.R')
  source('mixing.R')
  
  par(mar = c(12, 3, 1, 1), las=2)
  plot(pfm(b, gamma=gamma, mu=mu, mixing_structure=mixing_structure,
           density=density, N=N)*N,
       col="red", pch=19, xlab="", ylab="# infected",axes=FALSE, ylim=c(0,120))
  par(new=T)
  plot(M, col="blue", axes=F, pch="+", ylab="", xlab="")
  axis(1,labels=labels,at=1:38,lwd.ticks=0)
  axis(2,at=20*0:120)
}
