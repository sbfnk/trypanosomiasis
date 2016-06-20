library(ggplot2)
library(grid)

rt_humdom <- read.csv('reintroduction_times_humdom.csv', header=T, sep=',')
rt_hum_dom <- read.csv('reintroduction_times_hum_dom.csv', header=T, sep=',')

rt_humdom$time <- as.numeric(as.character(rt_humdom$time))
rt_hum_dom$time <- as.numeric(as.character(rt_humdom$time))

rt_humdom <- rt_humdom[rt_humdom$time!=Inf,]
rt_hum_dom <- rt_hum_dom[rt_hum_dom$time!=Inf,]

rt_humdom$xi <- as.factor(rt_humdom$xi)
rt_humdom$xi <- factor(rt_humdom$xi, levels = rev(levels(rt_humdom$xi)))
rt_hum_dom$xi <- as.factor(rt_hum_dom$xi)
rt_hum_dom$xi <- factor(rt_hum_dom$xi, levels = rev(levels(rt_hum_dom$xi)))

postscript(
    "reintroduction_times_humdom.eps", onefile=F, horizontal=F,
    paper="special", width=3.27, height=3.27,
    family = c("arial.afm.gz", "arialbd.afm.gz", "ariali.afm.gz",
    "arialbi.afm.gz") 
    )
ggplot(rt_humdom, aes(x=time,color=xi))+
  geom_line(size=1.25, stat="density")+
  scale_x_log10("Time until reintroduction (in yr)",
                limits=c(0.01,700),
                breaks=c(0.1,1,10,100),
                labels=c("0.1","1","10","100"))+
  theme_bw(12)+
  theme(
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.x = element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10),
        plot.margin=unit(c(2,2,2,2), "points"),
        legend.position=c(0.82,0.75),
        legend.key=element_rect(color=0)
#        legend.title=theme_text(hjust=3)
       )+
#  scale_color_brewer(substitute(paste(fb*yr^{-1}*lb, sep=""),
#    list(fb="Switch rate\n (in ", lb=")")), palette="Set1")+
  scale_color_brewer("Switch rate\n(per yr)",
                     palette="Set1")+
  scale_y_continuous("Probability density", limits=c(0,1.65))
dev.off()

postscript(
    "reintroduction_times_hum_dom.eps", onefile=F, horizontal=F,
    paper="special", width=3.27, height=3.27,
    family = c("arial.afm.gz", "arialbd.afm.gz", "ariali.afm.gz",
    "arialbi.afm.gz") 
    )
ggplot(rt_hum_dom, aes(x=time,color=xi))+
  geom_line(size=1.25, stat="density")+
  scale_x_log10("Time until reintroduction (in yr)",
                limits=c(0.01,500),
                breaks=c(0.1,1,10,100),
                labels=c("0.1","1","10","100"))+
  theme_bw(12)+
  ggtitle("Probability density")+
  theme(
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.x = element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10),
        plot.margin=unit(c(2,2,2,2), "points"),
        legend.position=c(0.82,0.76)
       )+
  scale_color_brewer("Switch rate",palette="Set1")+
  scale_y_continuous(limits=c(0,1.65))
dev.off()
