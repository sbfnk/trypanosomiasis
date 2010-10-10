#!/usr/bin/Rscript

# reservoirs.R
# Work out reservoirs from mixing matrix and prevalence data
# Invoke % Rscript reservoirs.R 

library('getopt');
library('BB');
library('mnormt');

source('mffoi.R');
source('findres.R');
source('foifm.R');
source('pfm.R');
source('mixing.R');
source('likelihood.R');
source('ilikelihood.R');

opt = getopt(c(
  'mixing', 'm', 1, "character",
  'data', 'd', 1, "character",
  'density', 'y', 0, "logical",
  'iterations', 'i', 1, "integer"
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

ndata <- nrow(data)

cat ("Target prev=", M/N, "\n")

if (opt$mixing == "random") {
  mixing_matrix <- NA
  b <- rep(1, length(N))
} else {
  mixing_matrix <- read.csv(file=opt$mixing, head=FALSE, sep=",")
  b <- rep(1, max(mixing_matrix))
}

if (is.null(opt$iterations)) {
  iter <- 1000
} else {
  iter <- opt$iterations
}
  
l <- -Inf
  
for (i in 1:iter) {
  # propose update
  saveb <- b
  savel <- l
  
  mult <- FALSE

#  nr <- sample(length(b), 1)
#  nr <- 1
  nr <- ceiling(abs(rmnorm(1, 0, 2)))

  el <- sample(length(b), nr)
#  b[el] <- rmnorm(1, b[el], diag(b[el]/100, nrow=nr, ncol=nr))
#  b[b<0] <- 0 
#  r <- sample(9*4, nr, replace = TRUE)
  for (j in 1:nr) {
#    if (r[j] > 18) {
#      r[j] <- r[j] - 18 
#      if (r[j] > 9) {
#        r[j] <- r[j] - 9
#        b[el[j]] <- b[el[j]] * r[j]
#      } else {
#        b[el[j]] <- b[el[j]] / r[j]
#      }
#    } else {
      nonzero <- TRUE
      while (nonzero == TRUE) {
        b[el[j]] <- rnorm(1, mean=b[el[j]], sd=b[el[j]]/10)
#        cat(j, " ", b[el[j]], "\n")
        if (b[el[j]]>0) { nonzero <- FALSE }
      }
    }
#  }
      
  prev <- pfm(pars=b, gamma=gamma, mu=mu, mixing_structure=mixing_matrix, 
              density=!is.null(opt$density), N=N)
  cat ("i=", i, ", testb=", b, "\n")
  cat ("i=", i, ", prev=",prev, "\n")
  l <- ilikelihood(prev, mu, gamma, M, N)

  cat ("i=", i, ", savel=", savel, ", l=", l, ", savel-l=", savel-l, "\n")
#  accept <- min(c(1, exp(-(savel-l))))
#  cat ("i=", i, ", accept=",accept,"savel=",savel,"l=",l,"\n")
#  if (!is.nan(accept) && runif(1) < accept) {
   if (!is.nan(l) && l > savel) {
    # accept
    cat ("i=", i, ", accepted\n")
    savel <- l
    saveb <- b
  } else {
    l <- savel
    b <- saveb
  }
  cat("i=", i, ", b=", saveb,"\n")
}
