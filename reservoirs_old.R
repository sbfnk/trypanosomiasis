#!/usr/bin/Rscript

# reservoirs.R
# Work out reservoirs from mixing matrix and prevalence data
# Invoke % Rscript reservoirs.R 

library('getopt');

opt = getopt(c(
  'mixing', 'm', 1, "character"
));

# estimate mixing matrix from data -- either assume random mixing or get 
# structure of the mixing matrix from csv file
stopifnot ( !is.null(opt$mixing) );
  
if (opt$mixing == "random") {
  
} else {
  mixing <- read.csv(file=opt$mixing, head=FALSE, sep=",")
}



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
