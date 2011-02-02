pdf(paste("species_data.pdf"))
for (i in 1:12)  {
  species_data<- sqrt(res[,13+i]*res[,13+i]*theta[i]*biting_rate*area_convert*vdensity/rabundance[i]/vmu*theta[i]*biting_rate/(rgamma[i]+rmu[i]))
  hist(main=data$name[i], breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
  lines(density(species_data), col="red")
  CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
  abline(v=CIrand, col="blue", lwd=2)
  abline(v=1, col="green", lwd=2)
}
dev.off()

pdf(paste("domestic.pdf"))
species_data <- rep(0, 1000000)
for (i in 2:4) {
  species_data <- species_data +
    res[,i+13]*res[,i+13]*theta[i]*biting_rate*area_convert*vdensity/rabundance[i]/vmu*theta[i]*biting_rate/(rgamma[i]+rmu[i])
}
species_data <- sqrt(species_data)
hist(main="Domestic cycle", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
lines(density(species_data), col="red")
CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
abline(v=CIrand, col="blue", lwd=2)
abline(v=1, col="green", lwd=2)
dev.off()

pdf(paste("wildlife.pdf"))
species_data <- rep(0, 1000000)
for (i in 5:12) {
  species_data <- species_data +
    res[,i+13]*res[,i+13]*theta[i]*biting_rate*area_convert*vdensity/rabundance[i]/vmu*theta[i]*biting_rate/(rgamma[i]+rmu[i])
}
species_data <- sqrt(species_data)
hist(main="Wildlife cycle", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
lines(density(species_data), col="red")
CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
abline(v=CIrand, col="blue", lwd=2)
abline(v=1, col="green", lwd=2)
dev.off()

pdf(paste("dom_wild.pdf"))
species_data <- rep(0, 1000000)
for (i in 2:12) {
  species_data <- species_data +
    res[,i+13]*res[,i+13]*theta[i]*biting_rate*area_convert*vdensity/rabundance[i]/vmu*theta[i]*biting_rate/(rgamma[i]+rmu[i])
}
species_data <- sqrt(species_data)
hist(main="Domestic+Wildlife", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
lines(density(species_data), col="red")
CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
abline(v=CIrand, col="blue", lwd=2)
abline(v=1, col="green", lwd=2)
dev.off()

pdf(paste("r0.pdf"))
species_data <- rep(0, 1000000)
for (i in 1:12) {
  species_data <- species_data +
    res[,i+13]*res[,i+13]*theta[i]*biting_rate*area_convert*vdensity/rabundance[i]/vmu*theta[i]*biting_rate/(rgamma[i]+rmu[i])
}
species_data <- sqrt(species_data)
hist(main="R0", breaks=100, freq=F, x=species_data, xlab="R0", ylab="Likelihood")
lines(density(species_data), col="red")
CIrand <- quantile(x=species_data, probs=c(0.025, 0.975))
abline(v=CIrand, col="blue", lwd=2)
abline(v=1, col="green", lwd=2)
dev.off()
