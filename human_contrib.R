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

human_data<-res_data[,1+offset]
h <- hist(main=nam[1+offset], breaks=100, freq=F, x=human_data, xlab="R0", ylab="Likelihood")
CIhuman <- quantile(x=human_data, probs=c(0.025, 0.975))

animal_data <- rep(0, nrow(res_data))
for (i in 2:12) {
  animal_data <- animal_data +
    res_data[,i+offset]^2
}
h <- hist(main=nam[1+offset], breaks=100, freq=F, x=animal_data, xlab="R0", ylab="Likelihood")
CIanimal <- quantile(x=animal_data, probs=c(0.025, 0.975))

cat(mean(human_data), ",", CIhuman[1], ",", CIhuman[2], ",", mean(animal_data), ",", CIanimal[1], ",", CIanimal[2], sep="")
