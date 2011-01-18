#!/usr/bin/Rscript

# reservoirs_invertmatrix.R
# Work out reservoirs from mixing matrix and prevalence data by matrix
# inversion (after conversation with Nishiura)
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
  'data', 'd', 1, "character",
  'vector', 'v', 1, "character",
  'params', 'p', 1, "character",
  'iterations', 'i', 1, "integer",
  'ignorezeroes', 'z', 0, "logical",
  'gambiense', 'g', 0, "logical",
  'nongambiense', 'n', 0, "logical",
  'convergence', 'c', 1, "integer",
  'area_convert', 'a', 0, "logical"
))

# read data
stopifnot ( !is.null(opt$data) );
stopifnot ( !is.null(opt$vector) );
data <- read.csv(file=opt$data, head=TRUE, sep=",")
vector <- read.csv(file=opt$vector, head=TRUE, sep=",")
params <- read.csv(file=opt$params, head=TRUE, sep=",")

# estimate mixing matrix from data -- either assume random mixing or get 
# structure of the mixing matrix from csv file
stopifnot ( !is.null(opt$mixing) );

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
vprev <- vM/vN
rmu <- data$mortality
vmu <- vector$mortality
rgamma <- data$rec_rate
vgamma <- vector$rec_rate
rlambda <- rprev/(1-rprev)*(rmu + rgamma)
vlambda <- vprev/(1-vprev)*(vmu + vgamma)
rtheta <- data$theta
rabundance <- data$abundance
vdensity <- vector$density

biting_rate <- vector$biting_rate
if (!is.null(opt$area_convert) {
  area_convert <- params$area_convert
}

if (!is.null(opt$ignorezeroes)) {
  rN <- rN[rM>0]
  rmu <- rmu[rM>0]
  rgamma <- rgamma[rM>0]
  vN <- vN[vM>0]
  vmu <- vmu[vM>0]
  vgamma <- vgamma[vM>0]
}

cat ("Species: ", rownames(data)[rM>0], "\n")
cat ("Vectors: ", rownames(vector)[vM>0], "\n")

if (!is.null(opt$ignorezeroes)) {
  rM <- rM[rM>0]
  vM <- vM[vM>0]
}

beta <- betaffoiv(rlambda, vlambda, theta, biting_rate, rabundance, vdensity,
                  area_convert, vprev)

NGM <- matrix(0,length(theta), length(theta))
  NGM[1,2:length(theta)] <- beta / (mu + gamma)
  NGM[2:length(theta), 1] <- sum(beta/N)
}
