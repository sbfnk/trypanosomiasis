#!/usr/bin/Rscript

# reservoirs.R
# Work out reservoirs from mixing matrix and prevalence data
# Invoke % Rscript reservoirs.R 

library('getopt');

source('mffoi.R');
source('findres.R');

opt = getopt(c(
  'mixing', 'm', 1, "character",
  'data', 'd', 1, "character",
  'density', 'y', 0, "logical"
));

# read data
stopifnot ( !is.null(opt$data) );
data <- read.csv(file=opt$data, head=TRUE, sep=",")

# estimate mixing matrix from data -- either assume random mixing or get 
# structure of the mixing matrix from csv file
stopifnot ( !is.null(opt$mixing) );

M <- data$pos_tbg + data$pos_tbng
N <- data$N
mu <- data$mortality
gamma <- data$rec_rate

if (opt$mixing == "random") {
  lambda <- M/(N-M)*(mu+gamma)
  if (!is.null(opt$density)) {
    b <- lambda/sqrt(sum(N*lambda*lambda/(lambda+mu+gamma)))
  } else {
    b <- lambda/sqrt(sum(lambda*lambda/(lambda+mu+gamma)))
  }
  beta <- b %o% b
} else {
  mixing <- read.csv(file=opt$mixing, head=FALSE, sep=",")
  if (max(mixing) == nrow(data)) {
    # number of parameters equal to number of data rows -> estimate beta
    lambda <- M/(N-M)*(mu+gamma)
    if (!is.null(opt$density)) {
      beta <- mffoi(lambda, mixing, gamma, mu, TRUE, N)
    } else {
      beta <- mffoi(lambda, mixing, gamma, mu)
    }
  } else {
    # number of parameters smaller than number of data rows -> estimate lambda
    # this is more complicated -- will have to do some kind of sophisticated
    # random walk
    #lambda <- foifm(mixing, 
  }
}

if (!is.null(opt$density)) {
  beta <- t(beta * N)
}

R <- beta %*% diag(1/(mu+lambda))
R0 <- max(Re(eigen(R)$values))

cat ("R0 = ", R0, "\n")

found <- FALSE
depth <- 1
while (!found) {
  found <- findres(R, matrix(0,38,38), depth)
  depth <- depth + 1
}


unit=diag(1,38,38)
for (i in c(1:38)) {
  proj=matrix(0,38,38)
  proj[i,i]=1
   U=proj %*% R
   Q=(unit-proj) %*% R
   u_sr=max(abs(eigen(U)$values))
   q_sr=max(abs(eigen(Q)$values))
   cat ("i = ", sprintf("%2d", i), ", u_sr = ", sprintf("%.2f", u_sr),
        ", q_sr = ", sprintf("%.2f", q_sr), "\n")
}
