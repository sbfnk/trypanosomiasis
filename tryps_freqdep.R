#!/usr/bin/Rscript

# reservoirs.R
# Work out reservoirs from mixing matrix and prevalence data
# Invoke % Rscript tryps_freqdep.R 

library('getopt');

opt = getopt(c(
  'verbose', 'v', 2, "integer",
  'help', 'h', 0, "logical",
  'count', 'c', 1, "integer",
  'mean', 'm', 1, "double",
  'sd', 's', 1, "double"
));


mixing <- read.csv(file="mixing.csv", head=FALSE, sep=",")



# tryps <- read.csv(file="epidemic.csv",head=TRUE,sep=",")
# l=(tryps$pos_tbng + tryps$pos_tbg)/(tryps$N-(tryps$pos_tbng + tryps$pos_tbg)) * (tryps$mortality + tryps$rec_rate)
#b=l/sqrt(sum(l*l/(l+tryps$mortality+tryps$rec_rate)))
#a=l2/(b*sum(b*l2/(l2+tryps$mortality+tryps$rec_rate)))
#a=l2/l
#a[is.nan(a)] <- 0
#R = ((a*b) %o% b) %*% diag(1/(tryps$mortality+tryps$rec_rate))
#R0 = max(Re(eigen(R)$values))
#cat ("R0 = ", R0, "\n")
#unit=diag(1,38,38)
#for (i in c(1:38)) {
#  proj=matrix(0,38,38)
#  proj[i,i]=1
#  U=(unit-proj) %*% R
#  V=proj %*% R
#  usr=max(abs(eigen(U)$values))
#  vsr=max(abs(eigen(V)$values))
#  cat ("i = ", i, ", usr = ", usr, ", vsr = ", vsr, "\n")
#}
