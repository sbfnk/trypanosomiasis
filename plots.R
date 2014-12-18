library(ggplot2)
library(gdata)

data <- read.csv("all.csv", header=T, sep=",")

melted <-
  melt(data,id.var=c("xi", "params", "groups", "human_frac", "habitat_mixing",
  "Human_prev", "Sheep_prev", "Goat_prev", "Pig_prev",
  "White.eyelid.mangabey_prev", "G..white.nosed.monkey_prev",
  "Blackstriped.duiker_prev", "Blue.duiker_prev",
  "Brush.tailed.porcupine_prev", "Giant.rat_prev", "Small.spotted.genet_prev",
  "Two.spotted.palm.civet_prev", "G._palpalis_gambiense_prev"),
  variable_name="contrib")

melted$contrib_type="single"
melted[melted$contrib=="Human",]$contrib_type="human"
melted[melted$contrib=="Domestic.animals",]$contrib_type="group"
melted[melted$contrib=="Wildlife.animals",]$contrib_type="group"
melted[melted$contrib=="All.animals",]$contrib_type="groups"
melted[melted$contrib=="R0",]$contrib_type="R0"

pdf("human_fraction.pdf")
ggplot(drop.levels(subset(melted, params=="pigsrecover" & groups=="single" &
  habitat_mixing=="none" & xi>9000 &
  (contrib_type=="groups" | contrib_type=="human" | contrib_type == "R0")),reorder=F), aes(x=human_frac/100, y=value, color=contrib)) +
  geom_line(lwd=1.5)+scale_y_continuous("Contribution", limits=c(0,1.2))+scale_x_continuous("Fraction
  of humans involved")+scale_colour_discrete("")
ggplot(drop.levels(subset(melted, params=="pigsrecover" & groups=="random_prev" &
  (contrib_type=="groups" | contrib_type=="human" | contrib_type == "R0")),reorder=F), aes(x=human_frac/100, y=value, color=contrib)) +
  geom_line(lwd=1.5)+scale_y_continuous("Contribution", limits=c(0,1.2))+scale_x_continuous("Fraction
  asymptomatic")+scale_colour_discrete("")
dev.off()


