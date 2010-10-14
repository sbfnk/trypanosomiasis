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
  'iterations', 'i', 1, "integer",
  'ignorezeroes', 'z', 0, "logical",
  'gambiense', 'g', 0, "logical",
  'nongambiense', 'n', 0, "logical",
  'convergence', 'c', 1, "integer"
));

# read data
stopifnot ( !is.null(opt$data) );
data <- read.csv(file=opt$data, head=TRUE, sep=",")

# estimate mixing matrix from data -- either assume random mixing or get 
# structure of the mixing matrix from csv file
stopifnot ( !is.null(opt$mixing) );

if (is.null(opt$gambiense)) {
   if (is.null(opt$nongambiense)) {
     M <- data$pos_tbg + data$pos_tbng
   } else {
     M <- data$pos_tbng
   }
} else {
  M <- data$pos_tbg
}


N <- data$N
mu <- data$mortality
gamma <- data$rec_rate

if (!is.null(opt$ignorezeroes)) {
  N <- N[M>0]
  mu <- mu[M>0]
  gamma <- gamma[M>0]
}

if (opt$mixing == "random") {
  mixing_matrix <- NA
  mixing_els <- 1:length(M)
  if (!is.null(opt$ignorezeroes)) {
    mixing_els <- mixing_els[M>0]
  }
  b <- rep(0, length(M))
} else {
  mixing_matrix <- as.matrix(read.csv(file=opt$mixing, head=FALSE, sep=","))
  if (!is.null(opt$ignorezeroes)) {
    mixing_matrix <- mixing_matrix[M>0,][,M>0]
  }
  mixing_els <- unique(as.vector(mixing_matrix))
  mixing_els <- mixing_els[mixing_els != 0]
  cat ("mixing els ", mixing_els, "\n")
  b <- rep(0, max(mixing_matrix))
}

b[mixing_els] <- 0.1
cat ("Species: ", rownames(data)[M>0], "\n")

if (!is.null(opt$ignorezeroes)) {
  M <- M[M>0]
}

cat ("Target prev=", M/N, "\n")

if (is.null(opt$iterations)) {
  iter <- 1000
} else {
  iter <- opt$iterations
}
  
l <- -Inf

## cat ("gamma=", gamma, "\n")
## cat ("mu=", mu, "\n")
## cat ("N=", N, "\n")
## cat ("mixing_structure=", mixing_matrix, "\n")

maxl <- -Inf
maxb <- b
convcount <- 0

if (!is.null(opt$convergence)) {
  iter <- Inf
}

i <- 0

while (i < iter && (is.null(opt$convergence) || convcount < opt$convergence)) {
  i <- i + 1
  # propose update
  saveb <- b
  savel <- l

  convcount <- convcount + 1
  
  mult <- FALSE

  nr <- ceiling(abs(rmnorm(1, 0, 2)))

  el <- sample(x=mixing_els, size=nr)
  for (j in 1:nr) {
    nonzero <- TRUE
    while (nonzero == TRUE) {
      test <- rnorm(1, mean=b[el[j]], sd=b[el[j]]/10)
      if (test>0) { 
        nonzero <- FALSE
        b[el[j]] <- test 
      }
    }
  }
  
  prev <- pfm(pars=b, gamma=gamma, mu=mu, mixing_structure=mixing_matrix, 
              density=!is.null(opt$density), N=N)
  cat ("i=", i, ", testb=", b, "\n")
  cat ("i=", i, ", prev=",prev, "\n")
  l <- ilikelihood(prev, mu, gamma, M, N)

  cat ("i=", i, ", savel=", savel, ", l=", l, ", savel-l=", savel-l, ", maxl=",
       maxl, "\n")
  accept <- min(c(1, exp(-(savel-l))))
  if (!is.nan(accept) && runif(1) < accept) {
    # accept
    cat ("i=", i, ", accepted\n")
    savel <- l
    saveb <- b
    if (l > maxl) {
      maxl <- l
      maxb <- b
      cat ("i=", i, ", reset_convcount=", convcount, "\n")
      convcount <- 0
    }
  } else {
    cat ("i=", i, ", convcount=", convcount, "\n")
    l <- savel
    b <- saveb
  }
  cat("i=", i, ", b=", b,"\n")
  cat("i=", i, ", maxb=", maxb,"\n")
  if (!is.null(opt$convergence) && convcount == opt$convergence) {
    cat ("Converged\n")
    cat("i=", i, ", maxl=", maxl,"\n")
    cat("i=", i, ", b=", maxb,"\n")
    break
  }
}
