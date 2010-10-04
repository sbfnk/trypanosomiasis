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

cat ("Target prev=", M/N, "\n")

if (opt$mixing == "random") {

  a <- rep(1, length(N))
  
  l <- -Inf
  
  for (i in 1:1000) {
    # propose update
    savea <- a
    savel <- l
    
    mult <- FALSE
    
    r <- sample(length(a)*9*2-2, 1)
    r <- r - length(a)*9
    if (r > 0) {
      mult <- TRUE
    } else {
      r <- -r
    }
    factor <- r %% 9 + 2
    el <- r %/% 9 + 1
    
    if (mult) {
      a[el] <- a[el] * factor
    } else {
      a[el] <- a[el] / factor
    }
      
    prev <- prand(a, gamma, mu, !is.null(opt$density), N)
    cat ("prev=",prev,"\n")
    l <- ilikelihood(prev, mu, gamma, M, N)

    accept <- min(c(1, exp(-(savel-l))))
    cat ("accept=",accept,"savel=",savel,"l=",l,"\n")
    if (runif(1) < accept) {
        # accept
        savel <- l
        savea <- a
      } else {
        l <- savel
        a <- savea
      }
      cat("a=", savea,"\n")
    }
  }
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
      
      r <- sample(max(mixing)*9*2-2, 1)
      r <- r - max(mixing)*9
      if (r > 0) {
        mult <- TRUE
      } else {
        r <- -r
      }
      factor <- r %% 9 + 2
      el <- r %/% 9 + 1

      if (mult) {
        b[el] <- b[el] * factor
      } else {
        b[el] <- b[el] / factor
      }
      
#      cat ("el=", el, ", factor=", factor, ", mult=", mult, "\n")
#      cat(b, "\n")
      ##   el <- el - max(mixing)
      ##   b[el] = b[el] - stepsize
      ## } else {
      ##   b[el] = b[el] + stepsize
      ## }

      prev <- pfm(mixing, b, gamma, mu, !is.null(opt$density), N)
      cat ("prev=",prev,"\n")
      l <- ilikelihood(prev, mu, gamma, M, N)

#      cat ("savel=", savel, ", l=", l, ", savel-l=", savel-l, "\n")
      accept <- min(c(1, exp(-(savel-l))))
      cat ("accept=",accept,"savel=",savel,"l=",l,"\n")
      if (runif(1) < accept) {
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
  }
}
