confint <- function(x) {
  x <- na.omit(x)
  quant <- quantile(x, prob=c(0.025,0.975), names=F)
  mean <- mean(x)
  data.frame(y=mean, ymin=quant[1], ymax=quant[2])
}

pdf('human_contrib.pdf')
ggplot(subset(species_data, groups=="random"), aes(x=N/3540,y=Human))+
  stat_summary(fun.data="confint", geom="smooth", colour="red",
               alpha=0.4)+ theme_bw()+
  scale_x_continuous("Fraction of human population exposed")+
  scale_y_continuous(substitute(R[0]*" in system of humans and vector"))+
  geom_abline(intercept=1,slope=0) 
dev.off()

human_animal_random <- melt(subset(species_data, groups=="random"),
  measure.vars=c("Human", "domestic.wildlife"), id.vars=c("N"),
  variable_name="contrib") 

pdf('human_contrib_animal.pdf')
ggplot(human_animal_random, aes(x=N/3540,y=value, linetype=contrib,
                                 color=contrib))+
  stat_summary(fun.data="confint", geom="smooth", alpha=0.4, lwd=2)+
  theme_bw(20)+
  scale_x_continuous("Fraction of human population exposed")+
  scale_y_continuous(substitute(R[0]*" in system of humans/animals and vector"),
                     limits=c(0,1.4))+
  geom_abline(intercept=1,slope=0,lwd=1.25)+
  opts(legend.position="none", axis.title.x = theme_text(size=20,vjust = -0.5),
       axis.title.y = theme_text(size=20,angle=90, vjust = 0.3))+
  geom_vline(xintercept=0.36061, lwd=1.25, linetype=2)+
  geom_vline(xintercept=0.32455, lwd=1.25, linetype=3)+
  scale_color_brewer(palette="Set1")
dev.off()



pdf('human_contrib_hum_domwild.pdf')
ggplot(subset(species_data, groups=="hum_domwild" & N>3000 &
  habitat=="none"), aes(x=xi,y=Human))+
  stat_summary(fun.data="confint", geom="smooth", colour="red",
               alpha=0.4)+ theme_bw()+
  scale_x_continuous("Rate of host switching between humans and animals", limits=c(0,10))+
  scale_y_continuous(substitute(R[0]*" in system of humans and vector"))+
  geom_abline(intercept=1,slope=0) 
dev.off()

human_animal_hum_domwild <-
  melt(subset(species_data, groups=="hum_domwild" & N>3000 & habitat=="none"),
  measure.vars=c("Human", "domestic.wildlife"), id.vars=c("xi"),
  variable_name="contrib")

pdf('human_contrib_switching_animal.pdf')
ggplot(human_animal_hum_domwild, aes(x=xi,y=value, linetype=contrib,
                                     color=contrib))+ 
  stat_summary(fun.data="confint", geom="smooth", alpha=0.4, lwd=2)+
  theme_bw(20)+
  scale_x_log10("Rate of host switching between humans and animals")+
  scale_y_continuous(substitute(R[0]*" in system of humans/animal and vector"),
                     limits=c(0,1.2))+
  geom_abline(intercept=1,slope=0, lwd=1.25)+
  opts(legend.position="none", axis.title.x = theme_text(size=20,vjust = -0.5),
       axis.title.y = theme_text(size=20,angle=90, vjust = 0.3))+
  scale_color_brewer(palette="Set1")
dev.off()
