library(ggplot2)
library(reshape)
library(data.table)
library(grid)

species_labels=c("Human", "Goat", "Pig", "Sheep", "Blackstriped duiker",
  "Blue duiker", "Brush-tailed porcupine", "Giant rat",
  "G. white-nosed monkey", "Small-spotted genet", "Two-spotted palm civet",
  "White-eyelid mangabey")
species_labels_r0 <- c(species_labels, expression(R[0]))
domain_labels=c("Human", "Human+domestic", "Human+wildlife", "Domestic",
  "Domestic+wildlife", "Wildlife", expression(R[0]))

quantiles <- melt(subset(species_data, groups=="random" & N>3000),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)
(CIR0 <- quantile(quantiles$R0, probs=c(0.025, 0.975)))
pdf("stage2_species_contributions.pdf")
#pdf("stage2_species_contributions.pdf")
ggplot(quantiles, aes(species, value))+
  geom_boxplot()+
  theme_bw(20)+
  scale_y_continuous(substitute("          Contribution to "*R[0]))+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()
r0_quantiles <- quantiles[quantiles$species=="Human",]
r0_quantiles$species <- "R0"
r0_quantiles$value <- r0_quantiles$R0
quantiles_r0 <- rbind(quantiles,r0_quantiles)
pdf("stage2_species_contributions_r0.pdf")
ggplot(quantiles_r0, aes(species, value))+
  geom_boxplot()+
  theme_bw(20)+
  scale_y_continuous(substitute("          Contribution to "*R[0]), limits=c(0,2.7))+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
       axis.title.x=element_blank())+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles_r0$species),labels=species_labels_r0)
dev.off()
rm(quantiles)

group_quantiles <- melt(subset(species_data, groups=="random" & N>3000),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
(CIHuman <- quantile(group_quantiles$Human, probs=c(0.025, 0.975)))
pdf("stage2_group_contributions.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw(20)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]),
  limits=c(0,4))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()
rm(group_quantiles)

quantiles <- melt(subset(species_data, groups=="humdom_wild" & N>3000 &
  habitat=="none" & xi < 0.01),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)
pdf("stage2_species_contributions_humdom_wild.pdf")
ggplot(quantiles, aes(species,
  value))+geom_boxplot()+theme_bw(20)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()
rm(quantiles)

group_quantiles <- melt(subset(species_data, groups=="humdom_wild" & N>3000 &
  habitat=="none" & xi < 0.01),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
pdf("stage2_group_contributions_humdom_wild.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw(20)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()
rm(group_quantiles)

quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>3000 &
  habitat=="none" & xi < 0.01),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)
pdf("stage2_species_contributions_hum_dom_wild.pdf")
ggplot(quantiles, aes(species,
  value))+geom_boxplot()+theme_bw(20)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()
rm(quantiles)

group_quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>3000 &
  habitat=="none" & xi < 0.01),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
pdf("stage2_group_contributions_hum_dom_wild.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw(20)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()
rm(group_quantiles)

quantiles <- melt(subset(species_data, groups=="hum_domwild" & N>3000 &
  habitat=="none" & xi < 0.01),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)
pdf("stage2_species_contributions_hum_domwild.pdf")
ggplot(quantiles, aes(species,
  value))+geom_boxplot()+theme_bw(20)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()
rm(quantiles)

group_quantiles <- melt(subset(species_data, groups=="hum_domwild" & N>3000 &
  habitat=="none" & xi < 0.01),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
pdf("stage2_group_contributions_hum_domwild.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw(20)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()
rm(group_quantiles)

quantiles <- melt(subset(species_data, groups=="single" & N>3000 &
  habitat=="binary" & xi > 10),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)
pdf("species_contributions_binary.pdf")
ggplot(quantiles, aes(species,
  value))+geom_boxplot()+theme_bw(20)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()
rm(quantiles)

group_quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>3000 &
  habitat=="binary" & xi > 10),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
pdf("group_contributions_binary.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw(20)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()
rm(group_quantiles)

quantiles <- melt(subset(species_data, groups=="single" & N>3000 &
  habitat=="fractional" & xi > 10),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)
pdf("stage2_species_contributions_fractional.pdf")
ggplot(quantiles, aes(species,
  value))+geom_boxplot()+theme_bw(20)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()
rm(quantiles)

group_quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>3000 &
  habitat=="fractional" & xi > 10),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
pdf("stage2_group_contributions_fractional.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw(20)+scale_y_continuous(substitute("          Contribution to "*R[0]))+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()
rm(group_quantiles)

## group_quantiles <- melt(subset(habitat_data_xitau, groups=="habitat" &  
##   habitat=="fractional"),
##   id=c("N","tau","groups","habitat",hosts[-1]),
##   variable_name="grouped")
group_quantiles <- melt(subset(habitat_data_new, groups=="habitat" &  
  habitat=="fractional" & xi==48.2578),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
(CIwild <- quantile(subset(group_quantiles, grouped=="wildlife")$value, probs=c(0.025, 0.53)))
postscript(
    "group_contributions_fractional.eps", onefile=F, horizontal=F, paper="special",
    width=3.27, height=3.27,
    family = c("arial.afm.gz", "arialbd.afm.gz", "ariali.afm.gz",
    "arialbi.afm.gz") 
    )
## pdf("stage2_group_contributions_single_fractional.pdf",
##     width=2.72,height=4.54,pointsize=8,
##     family = c("arial.afm", "arialbd.afm", "ariali.afm", "arialbi.afm"))
ggplot(group_quantiles, aes(grouped, value))+
  geom_boxplot(size=0.2, outlier.size=0.75)+
  theme_bw(12)+
  theme(
        axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=10),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=10),
        plot.title=element_text(size=12),
        plot.margin=unit(c(2,2,2,2), "points")
        )+
  scale_y_continuous(substitute("   Contribution to "*R[0]),
                     limits=c(0,3.5))+
  geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),
                   labels=domain_labels) 
dev.off()
rm(group_quantiles)

group_quantiles <- melt(subset(species_data, groups=="random"),
  id=c("N","xi","groups","habitat",hosts),
  variable_name="grouped")
quantiles <- melt(subset(species_data, groups=="random"), 
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles <- quantiles[!is.nan(quantiles$value),]
table(quantiles[quantiles$value>2,]$species) /
  table(quantiles$species)
group_plot <- ggplot(group_quantiles, aes(grouped, value))+
  geom_boxplot(size=0.2, outlier.size=0.75)+
  theme_bw(12)+
  theme(
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10),
        axis.title.y=element_blank(),
#        axis.ticks.y=element_blank(),
        plot.title=element_text(size=12),
        plot.margin=unit(c(0,0,0,0), "cm")
        )+
  ggtitle("Groups of species")+
  scale_y_continuous(substitute("Contribution to "*R[0]),
                     limits=c(0,3.5))+
  geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),
                   labels=domain_labels[-1])
single_plot <- ggplot(quantiles, aes(species, value))+
  geom_boxplot(size=0.2, outlier.size=0.75)+
  theme_bw(12)+
  scale_y_continuous(substitute("Contribution to "*R[0]),
                     limits=c(0,2))+
  theme(
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10),
        plot.title=element_text(size=12),
        plot.margin=unit(c(0,0,0,0), "cm")
        )+
  ggtitle("Individual species")+
  geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),
                   labels=species_labels)
vp1 <- viewport(
                width=(unit(0.59, "npc")-unit(2, "points")),
                height=(unit(1, "npc")-unit(4, "points")),
                x=(unit(0.295, "npc") + unit(1, "points")),
                y=0.5
                )
vp2 <- viewport(
                width=(unit(0.39, "npc")-unit(2, "points")),
                height=(unit(0.91, "npc")-unit(4, "points")),
                x=(unit(0.805, "npc") - unit(1, "points")),
                y=(unit(0.545, "npc"))
                )

postscript(
    "contributions_random_mixing.eps", onefile=F, horizontal=F,
    paper="special", 
    width=6.83, height=3.27,
    family = c("arial.afm.gz", "arialbd.afm.gz", "ariali.afm.gz",
    "arialbi.afm.gz") 
    )
print(single_plot, vp = vp1)
print(group_plot, vp = vp2)
grid.text("A", vp = vp1, x = 0,
          y = (unit(1, "npc") - unit(2, "points")),
          just=c("left", "top"),
          gp=gpar(fontface="bold", fontsize=12))
grid.text("B", vp = vp2, x = 0,
          y = (unit(1, "npc") - unit(2, "points")),
          just=c("left", "top"),
          gp=gpar(fontface="bold", fontsize=12))
dev.off()
rm(group_quantiles)
rm(quantiles)
