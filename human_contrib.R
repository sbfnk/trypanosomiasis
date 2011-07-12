#!/usr/bin/Rscript

# Work out the human contribution

library('getopt');

opt = getopt(c(
  'file', 'f', 1, "character",
  'offset', 'o', 1, "integer"
))

if (is.null(opt$offset)) {
  offset <- 27
} else {
  offset <- opt$offset
}

tab <- read.table(opt$file, header=T, sep=",")
nam <- names(tab)
res_data <- as.matrix(tab) 

species_data<-res_data[,1+offset]
h <- hist(main=nam[i+offset], breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
cat(mean(species_data), ",", CIrand[1], ",", CIrand[2])
