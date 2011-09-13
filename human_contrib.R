confint <- function(x) {
  x <- na.omit(x)
  quant <- quantile(x, prob=c(0.025,0.975), names=F)
  mean <- mean(x)
  data.frame(y=mean, ymin=quant[1], ymax=quant[2])
}

species_data <- data.frame(db[hosts], db[domains[-1]], N=renv$Human_N,
  xi=renv$G._palpalis_gambiense_xi, groups=renv$groups, habitat=renv$habitat)

ggplot(subset(species_data, groups=="random"), aes(x=N/3540,y=Human))+
  stat_summary(fun.data="confint", geom="smooth", colour="red",
               alpha=0.4)+ theme_bw()+
  scale_x_continuous("Fraction of human population exposed")+
  scale_y_continuous(substitute(R[0]*" in system of humans and vector"))+
  geom_abline(intercept=1,slope=0) 
