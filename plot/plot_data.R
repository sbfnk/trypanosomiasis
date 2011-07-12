#!/usr/bin/Rscript

# Work out reservoirs from mixing matrix and prevalence data by matrix
# inversion (after conversation with Nishiura)
# Invoke % Rscript reservoirs.R 

library('getopt');
library('ggplot2');

opt = getopt(c(
  'file', 'f', 1, "character",
  'offset', 'o', 1, "integer",
  'alpha', 'a', 0, "logical"
))

if (is.null(opt$offset)) {
  offset <- 27
} else {
  offset <- opt$offset
}

tab <- read.table(opt$file, header=T, sep=",")
nam <- names(tab)
res_data <- as.matrix(tab) 
#res_data <- matrix(scan(opt$file), byrow=TRUE, ncol=2*(offset-1))

pdf(paste("species_data.pdf"))
for (i in 1:12)  {
  species_data<-res_data[,i+offset]
  h <- hist(main=nam[i+offset], breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
  lines(density(species_data), col="red")
  CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
  abline(v=CIrand, col="blue", lwd=2)
  abline(v=1, col="green", lwd=2)
  mtext(paste("p(R0>1) =", sum(species_data>1)/nrow(res_data)))
  if (i == 1) {
    cat("humans",h$mids[h$intensities==max(h$intensities)],CIrand,"\n");
  }
}
dev.off()

pdf(paste("domestic.pdf"))
if (is.null(opt$noplot)) {
  species_data <- rep(0, nrow(res_data))
  for (i in 2:4) {
    species_data <- species_data +
      res_data[,i+offset]^2
  }
  species_data <- sqrt(species_data)
  hist(main="Domestic cycle", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
  lines(density(species_data), col="red")
  CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
  abline(v=CIrand, col="blue", lwd=2)
  abline(v=1, col="green", lwd=2)
  mtext(paste("p(R0>1) =", sum(species_data>1)/nrow(res_data)))
  dev.off()

  pdf(paste("wildlife.pdf"))
  species_data <- rep(0, nrow(res_data))
  for (i in 5:12) {
    species_data <- species_data +
      res_data[,i+offset]^2
  }
  species_data <- sqrt(species_data)
  hist(main="Wildlife cycle", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
  lines(density(species_data), col="red")
  CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
  abline(v=CIrand, col="blue", lwd=2)
  abline(v=1, col="green", lwd=2)
  mtext(paste("p(R0>1) =", sum(species_data>1)/nrow(res_data)))
}
dev.off()

pdf(paste("dom_wild.pdf"))
species_data <- rep(0, nrow(res_data))
for (i in 2:12) {
  species_data <- species_data +
    res_data[,i+offset]^2
}
species_data <- sqrt(species_data)
hist(main="Domestic+Wildlife", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
lines(density(species_data), col="red")
CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
abline(v=CIrand, col="blue", lwd=2)
abline(v=1, col="green", lwd=2)
mtext(paste("p(R0>1) =", sum(species_data>1)/nrow(res_data)))
cat("animals",h$mids[h$intensities==max(h$intensities)],CIrand,"\n");
dev.off()

pdf(paste("human_vs_wildlife.pdf"))
human_data<-res_data[,1+offset]
df=data.frame(humans=human_data, animals=species_data)
qplot(data=df, humans, animals, geom="point")
dev.off()

pdf(paste("r0.pdf"))
if (is.null(opt$noplot)) {
  species_data <- rep(0, nrow(res_data))
  for (i in 1:12) {
    species_data <- species_data +
      res_data[,i+offset]^2
  }
  species_data <- sqrt(species_data)
  hist(main="R0", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
  lines(density(species_data), col="red")
  CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
  abline(v=CIrand, col="blue", lwd=2)
  abline(v=1, col="green", lwd=2)
  mtext(paste("p(R0>1) =", sum(species_data>1)/nrow(res_data)))

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
  dev.off();
}
