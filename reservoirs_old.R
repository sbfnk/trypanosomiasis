#!/usr/bin/Rscript

# reservoirs.R
# Work out reservoirs from mixing matrix and prevalence data
# Invoke % Rscript reservoirs.R 

library('getopt');

opt = getopt(c(
  'mixing', 'm', 1, "character",
  'data', 'd', 1, "character"
));

# read data
stopifnot ( !is.null(opt$mixing) );
data <- read.csv(file=opt$data, head=TRUE, sep=",")

# estimate mixing matrix from data -- either assume random mixing or get 
# structure of the mixing matrix from csv file
stopifnot ( !is.null(opt$mixing) );

M <- data$pos_tbg + data$pos_tbg
N <- data$N
mu <- data$mortality
gamma <- data$rec_rate

if (opt$mixing == "random") {
  lambda <- M/(N-M)*(mu+gamma)
  b <- lambda/sqrt(sum(lambda*lambda/(lambda+mu+gamma)))
  beta <- b %o% b
} else {
  mixing <- read.csv(file=opt$mixing, head=FALSE, sep=",")
  if (max(mixing) == nrow(data)) {
    # number of parameters equal to number of data rows -> estimate beta
    lambda <- M/(N-M)*(mu+gamma)
    beta <- mffoi(lambda, mixing, gamma, mu)
  } else {
    # number of parameters smaller than number of data rows -> estimate lambda
    # this is more complicated -- will have to do some kind of sophisticated
    # random walk
    #lambda <- foifm(mixing, 
  }
}

R <- beta %*% diag(1/(mu+lambda))
R0 <- max(Re(eigen(R)$values))

cat ("R0 = ", R0, "\n")
unit=diag(1,38,38)
for (i in c(1:38)) {
  proj=matrix(0,38,38)
  proj[i,i]=1
  U=(unit-proj) %*% R
  V=proj %*% R
  usr=max(abs(eigen(U)$values))
  vsr=max(abs(eigen(V)$values))
  cat ("i = ", i, ", usr = ", usr, ", vsr = ", vsr, "\n")
}
