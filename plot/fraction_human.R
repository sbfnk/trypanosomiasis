t <- read.table('human_frac.dat', sep=",")
pdf("human_frac.pdf")
matplot(t[,1]/100,t[,2:7], type="l", lty=c(1,2,2,1,2,2), col=c(1,1,1,2,2,2),
        lwd=1.5, xlab="Fraction of humans involved", ylab="Contribution to R0")
legend("topleft", c("Humans", "Animals"), col=c(1,2), lwd=2, bty="n")
grid()
dev.off()

t <- read.table('asymptomatic_frac.dat', sep=",")
pdf("asymptotic_frac.pdf")
matplot(t[,1]/100,t[,2:7], type="l", lty=c(1,2,2,1,2,2), col=c(1,1,1,2,2,2),
        lwd=1.5, xlab="Fraction of humans involved", ylab="Contribution to R0")
legend("topleft", c("Humans", "Animals"), col=c(1,2), lwd=2, bty="n")
grid()
