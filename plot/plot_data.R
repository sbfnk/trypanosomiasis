#!/usr/bin/Rscript

# Work out reservoirs from mixing matrix and prevalence data by matrix
# inversion (after conversation with Nishiura)
# Invoke % Rscript reservoirs.R 

library('getopt');
library('ggplot2');

opt = getopt(c(
  'file', 'f', 1, "character",
  'data', 'd', 1, "character",
  'vector', 'v', 1, "character",
  'params', 'p', 1, "character",
  'area_convert', 'a', 0, "logical",
  'alpha', 'l', 0, "logical"
))

# read data
stopifnot ( !is.null(opt$data) );
stopifnot ( !is.null(opt$vector) );
data <- read.csv(file=opt$data, head=TRUE, sep=",")
vector <- read.csv(file=opt$vector, head=TRUE, sep=",")
params <- read.csv(file=opt$params, head=TRUE, sep=",")

offset <- 13
if (!is.null(opt$alpha)) {
  offset <- offset + 1
}

rmu <- data$mortality
vmu <- vector$mortality
rgamma <- data$rec_rate
theta <- data$theta
rabundance <- data$abundance
vdensity <- vector$density
biting_rate <- vector$biting_rate

if (!is.null(opt$area_convert)) {
  area_convert <- params$area_convert
} else {
  area_convert <- 1
}


res_data <- matrix(scan(opt$file), byrow=TRUE, ncol=2*offset-1)

pdf(paste("human_data.pdf"))
human_data <- sqrt(res_data[,offset+1])
hist(main="Human", breaks=100, freq=F, x=human_data, xlab="R0", ylab="Likelihood")
lines(density(human_data), col="red")
CIrand <- quantile(x=human_data, probs=c(0.025, 0.975))
abline(v=CIrand, col="blue", lwd=2)
abline(v=1, col="green", lwd=2)
mtext(paste("p(R0>1) =", sum(human_data>1)/1000000))
dev.off()

pdf(paste("species_data.pdf"))
for (i in 1:12)  {
  species_data<- sqrt(res_data[,i+offset])
  hist(main=data$name[i], breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
  lines(density(species_data), col="red")
  CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
  abline(v=CIrand, col="blue", lwd=2)
  abline(v=1, col="green", lwd=2)
  mtext(paste("p(R0>1) =", sum(species_data>1)/1000000))
}
dev.off()

pdf(paste("domestic.pdf"))
species_data <- rep(0, 1000000)
for (i in 2:4) {
  species_data <- species_data +
    res_data[,i+offset]
}
species_data <- sqrt(species_data)
hist(main="Domestic cycle", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
lines(density(species_data), col="red")
CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
abline(v=CIrand, col="blue", lwd=2)
abline(v=1, col="green", lwd=2)
mtext(paste("p(R0>1) =", sum(species_data>1)/1000000))
dev.off()

pdf(paste("wildlife.pdf"))
species_data <- rep(0, 1000000)
for (i in 5:12) {
  species_data <- species_data +
    res_data[,i+offset]
}
species_data <- sqrt(species_data)
hist(main="Wildlife cycle", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
lines(density(species_data), col="red")
CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
abline(v=CIrand, col="blue", lwd=2)
abline(v=1, col="green", lwd=2)
mtext(paste("p(R0>1) =", sum(species_data>1)/1000000))
dev.off()

pdf(paste("dom_wild.pdf"))
species_data <- rep(0, 1000000)
for (i in 2:12) {
  species_data <- species_data +
    res_data[,i+offset]
}
species_data <- sqrt(species_data)
hist(main="Domestic+Wildlife", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
lines(density(species_data), col="red")
CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
abline(v=CIrand, col="blue", lwd=2)
abline(v=1, col="green", lwd=2)
mtext(paste("p(R0>1) =", sum(species_data>1)/1000000))
dev.off()

## pdf(paste("human_vs_wildlife.pdf"))
## df=data.frame(humans=human_data, animals=species_data)
## qplot(data=df, humans, animals, geom="point")
## dev.off()

pdf(paste("r0.pdf"))
species_data <- rep(0, 1000000)
for (i in 1:12) {
  species_data <- species_data +
    res_data[,i+offset]
}
species_data <- sqrt(species_data)
hist(main="R0", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
lines(density(species_data), col="red")
CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
abline(v=CIrand, col="blue", lwd=2)
abline(v=1, col="green", lwd=2)
mtext(paste("p(R0>1) =", sum(species_data>1)/1000000))
dev.off()

if (!is.null(opt$alpha)) {
  pdf(paste("vprev.pdf"))
  hist(main="vprev", breaks=100, freq=F, x=res_data[,offset], xlab="vprev", ylab="Likelihood")
  lines(density(res_data[,offset]), col="red")
  CIrand <- quantile(x=res_data[,offset], probs=c(0.025, 0.975))
  abline(v=CIrand, col="blue", lwd=2)
  dev.off()
  pdf(paste("alpha.pdf"))
  hist(main="alpha", breaks=100, freq=F, x=res_data[,2*offset-1], xlab="alpha", ylab="Likelihood")
  lines(density(res_data[,2*offset-1]), col="red")
  CIrand <- quantile(x=res_data[,14], probs=c(0.025, 0.975))
  dev.off()
}
