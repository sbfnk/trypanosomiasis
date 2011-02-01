#!/usr/bin/Rscript

# Work out reservoirs from mixing matrix and prevalence data by matrix
# inversion (after conversation with Nishiura)
# Invoke % Rscript reservoirs.R 

library('getopt');
library('BB');
library('mnormt');

source('joinfactors.R');
source('betaffoiv.R');

opt = getopt(c(
  'data', 'd', 1, "character",
  'vector', 'v', 1, "character",
  'params', 'p', 1, "character",
  'iterations', 'i', 1, "integer",
  'ignorezeroes', 'z', 0, "logical",
  'gambiense', 'g', 0, "logical",
  'nongambiense', 'n', 0, "logical",
  'convergence', 'c', 1, "integer",
  'area_convert', 'a', 0, "logical",
  'vector_prevalence', 'l', 0, "logical"
))

# read data
stopifnot ( !is.null(opt$data) );
stopifnot ( !is.null(opt$vector) );
data <- read.csv(file=opt$data, head=TRUE, sep=",")
vector <- read.csv(file=opt$vector, head=TRUE, sep=",")
params <- read.csv(file=opt$params, head=TRUE, sep=",")

if (is.null(opt$gambiense)) {
   if (is.null(opt$nongambiense)) {
     rM <- data$pos_tbg + data$pos_tbng
     vM <- vector$pos_tbg + vector$pos_tbng
   } else {
     rM <- data$pos_tbng
     vM <- vector$pos_tbng
   }
} else {
  rM <- data$pos_tbg
  vM <- vector$pos_tbg
}


rN <- data$N
vN <- vector$N
rprev <- rM/rN
if (!is.null(opt$vector_prevalence)) {
  vprev <- vM/vN
} else {
  vprev <- NULL
}
rmu <- data$mortality
vmu <- vector$mortality
rgamma <- data$rec_rate
vgamma <- vector$rec_rate
rlambda <- rprev/(1-rprev)*(rmu + rgamma)
vlambda <- vprev/(1-vprev)*(vmu + vgamma)
theta <- data$theta
rabundance <- data$abundance
vdensity <- vector$density

biting_rate <- vector$biting_rate
if (!is.null(opt$area_convert)) {
  area_convert <- params$area_convert
} else {
  area_convert <- NULL
}

if (!is.null(opt$ignorezeroes)) {
  rN <- rN[rM>0]
  rmu <- rmu[rM>0]
  rgamma <- rgamma[rM>0]
  vN <- vN[vM>0]
  vmu <- vmu[vM>0]
  vgamma <- vgamma[vM>0]
}

if (!is.null(opt$ignorezeroes)) {
  rM <- rM[rM>0]
  vM <- vM[vM>0]
}

factor <- joinfactors(theta, biting_rate, area_convert)
res <- betaffoiv(rlambda, vdensity, rabundance, factor, vprev, rprev, vmu)
beta <- res$par

NGM <- matrix(0,length(rgamma)+length(vgamma), length(rgamma)+length(vgamma))

NGM[1:length(vgamma),(length(vgamma)+1):(length(vgamma)+length(rgamma))] <-
  1/vmu %*% t(beta) * factor * vdensity/rabundance
for (i in 1:length(vgamma)) {
  NGM[(length(vgamma)+1):(length(vgamma)+length(rgamma)), i] <-
    beta/(rgamma+rmu)
}

# one vector case

R0 <- sqrt(sum(NGM[2:(length(rgamma)+1), 1] * NGM[1, 2:(length(rgamma)+1)]))
cat ("R0=", R0, "\n");
for (i in 1:length(rgamma)) {
  cat (i, ": ", sqrt(NGM[i+1, 1]*NGM[1, i+1]), "\n")
}
