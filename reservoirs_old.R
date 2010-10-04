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
source('likelihood.R');
source('ilikelihood.R');

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

ndata <- nrow(data)

if (opt$mixing == "random") {
  lambda <- M/(N-M)*(mu+gamma )
  if (!is.null(opt$density)) {
    b <- lambda/sqrt(sum(N*lambda*lambda/(lambda+mu+gamma)))
  } else {
    b <- lambda/sqrt(sum(lambda*lambda/(lambda+mu+gamma)))
  }
  beta <- b %o% b
} else {
  mixing <- read.csv(file=opt$mixing, head=FALSE, sep=",")
  if (max(mixing) == ndata) {
    # number of parameters equal to number of data rows -> estimate beta
    lambda <- M/(N-M)*(mu+gamma)
    if (!is.null(opt$density)) {
      beta <- mffoi(lambda, mixing, gamma, mu, TRUE, N)
    } else {
      beta <- mffoi(lambda, mixing, gamma, mu)
    }
  } else {
    # number of parameters smaller than number of data rows -> estimate lambda

    stepsize <- 0.1
    b <- rep(1, max(mixing))
    l <- -Inf

    for (i in 1:1000) {
      # propose update
      saveb <- b
      savel <- l

      mult <- FALSE
      
      r <- sample(max(mixing)*10*2, 1)
      r <- r - max(mixing)*10
      if (r > 0) {
        mult <- TRUE
      } else {
        r <- -r
      }
      factor <- r %% 10 + 1
      el <- r %/% 10 + 1

      cat ("el=", el, "\n")
      if (mult) {
        b[el] <- b[el] * factor
      } else {
        b[el] <- b[el] / factor
      }
      
      ##   el <- el - max(mixing)
      ##   b[el] = b[el] - stepsize
      ## } else {
      ##   b[el] = b[el] + stepsize
      ## }

      prev <- pfm(mixing, b, gamma, mu, !is.null(opt$density), N)
      l <- ilikelihood(prev, mu, gamma, M, N)

      accept <- min(c(1, exp(-(savel-l))))
      cat ("accept=",accept,"savel=",savel,"l=",l,"\n")
      if (runif(1) < accept) {
        # accept
        savel <- l
        saveb <- b
      } else {
        l <- savel
        b <- saveb
      }
      cat(saveb,"\n")
    }
  }
}
