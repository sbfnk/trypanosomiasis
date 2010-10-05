#!/usr/bin/Rscript

# reservoirs.R
# Work out reservoirs from mixing matrix and prevalence data
# Invoke % Rscript reservoirs.R 

library('getopt');
library('BB');

source('mffoi.R');
source('findres.R');
source('foifm.R');
source('pfm.R');
source('prand.R');
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
  mixing <- NA
  b <- rep(1, length(N))
  fn <- prand
} else {
  mixing <- read.csv(file=opt$mixing, head=FALSE, sep=",")
  b <- rep(1, max(mixing))
  fn <- pfm
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

  nr <- sample(length(b), 1)
  el <- sample(length(b), nr)
  r <- sample(9*2, nr, replace = TRUE)
  for (j in 1:nr) {
    if (r[j] > 9) {
      mult <- TRUE
      r[j] <- r[j] - 9
    }
    
    if (mult) {
      b[el[j]] <- b[el[j]] * r[j]
    } else {
      b[el[j]] <- b[el[j]] / r[j]
    }
  }
      
  prev <- fn(mixing=mixing, pars=b, gamma=gamma, mu=mu,
            density=!is.null(opt$density), N=N)
  cat ("prev=",prev,"\n")
  l <- ilikelihood(prev, mu, gamma, M, N)

  cat ("savel=", savel, ", l=", l, ", savel-l=", savel-l, "\n")
  accept <- min(c(1, exp(-(savel-l))))
  cat ("accept=",accept,"savel=",savel,"l=",l,"\n")
  if (!is.nan(accept) && runif(1) < accept) {
    # accept
    #        cat ("accepted\n")
    savel <- l
    saveb <- b
  } else {
    l <- savel
    b <- saveb
  }
  cat("b=", saveb,"\n")
}
