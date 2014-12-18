library(ggplot2)

data1 <- read.csv('reintroduction_times_0_001.csv', header=F, sep=',')
data2 <- read.csv('reintroduction_times_0_01.csv', header=F, sep=',')
data3 <- read.csv('reintroduction_times_0_1.csv', header=F, sep=',')
data4 <- read.csv('reintroduction_times_1.csv', header=F, sep=',')

data <- rbind(data1, data2, data3, data4)

rtimes <- data[data$V2!=Inf,]

names(rtimes) <- c("xi", "years")
rtimes$xi <- as.factor(rtimes$xi)

pdf('reintroduction_times.pdf')
ggplot(rtimes,
  aes(x=years,fill=xi))+geom_density(alpha=0.2)+scale_x_log10(limits=c(0.01,100),
  breaks=c(0.1,1,10,50,100),
  labels=c("0.1","1","10","50","100"))+theme_bw()+opts(panel.grid.major=theme_blank(),
  panel.grid.minor=theme_blank(), title)
dev.off()

pdf('reintroduction_times2.pdf')
ggplot(rtimes,
  aes(x=years,color=xi))+geom_density(size=2)+scale_x_log10(limits=c(0.01,100),
  breaks=c(0.1,1,10,50,100),
  labels=c("0.1","1","10","50","100"))+theme_bw()+opts(panel.grid.major=theme_blank(),
  panel.grid.minor=theme_blank(), title)
dev.off()
