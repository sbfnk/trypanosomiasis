library(ggplot2)

species_labels=c("Human", "Goat", "Pig", "Sheep", "Blackstriped duiker",
  "Blue duiker", "Brush-tailed porcupine", "Giant rat",
  "G. white-nosed monkey", "Small-spotted genet", "Two-spotted palm civet",
  "White-eyelid mangabey")
domain_labels=c("Human", "Human+domestic", "Human+wildlife", "Domestic",
  "Domestic+wildlife", "Wildlife", expression(R[0]))

quantiles <- melt(subset(species_data, groups=="random" & N>1113 & N<1545),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)
(CIR0 <- quantile(quantiles$R0, probs=c(0.025, 0.975)))
svg("species_contributions.svg")
#pdf("species_contributions.pdf")
ggplot(quantiles, aes(species, value))+
  geom_boxplot()+
  theme_bw(20)+
  scale_y_continuous(substitute("          Contribution to "*R[0]))+
  opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=theme_blank())+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()


## quantiles <- melt(subset(species_data, groups=="random" & N>3000),
##   id=c("N","xi","groups","habitat",domains[-1]),
##   variable_name="species")
## quantiles$species <- factor(quantiles$species, levels=hosts)
## r0_quantiles <- quantiles[quantiles$species=="Human",]
## r0_quantiles$species <- "R0"
## r0_quantiles$value <- r0_quantiles$R0
## quantiles_r0 <- rbind(quantiles,r0_quantiles)
## (CIR0 <- quantile(quantiles$R0, probs=c(0.025, 0.975)))
## species_labels_r0 <- c("Human", "Goat", "Pig", "Sheep", "Blackstriped duiker",
##   "Blue duiker", "Brush-tailed porcupine", "Giant rat",
##   "G. white-nosed monkey", "Small-spotted genet", "Two-spotted palm civet",
##   "White-eyelid mangabey", expression(R[0]))
## svg("species_contributions_r0.svg")
## ggplot(quantiles_r0, aes(species, value))+
##   geom_boxplot()+
##   theme_bw(20)+
##   scale_y_continuous(substitute("          Contribution to "*R[0]), limits=c(0,2.7))+
##   opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
##        axis.title.x=theme_blank())+geom_hline(yintercept=1)+
##   scale_x_discrete(breaks=levels(quantiles_r0$species),labels=species_labels_r0)
## dev.off()

## group_quantiles <- melt(subset(species_data_nh, groups=="random" & N>3000),
##   id=c("N","xi","groups","habitat",hosts[-1]),
##   variable_name="grouped")
## (CIHuman <- quantile(group_quantiles$Human, probs=c(0.025, 0.975)))
## svg("group_contributions.svg")
## ggplot(group_quantiles, aes(grouped,
##   value))+geom_boxplot()+theme_bw(20)+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
##   axis.title.x=theme_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]),
##   limits=c(0,4))+geom_hline(yintercept=1)+
##   scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
## dev.off()

## quantiles <- melt(subset(species_data, groups=="humdom_wild" & N>1113 & N<1545 &
##   habitat=="none" & xi < 0.01),
##   id=c("N","xi","groups","habitat",domains[-1]),
##   variable_name="species")
## quantiles$species <- factor(quantiles$species, levels=hosts)

## svg("species_contributions_humdom_wild.svg")
## ggplot(quantiles, aes(species,
##   value))+geom_boxplot()+theme_bw(20)+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
##   axis.title.x=theme_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
##   scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
## dev.off()

## group_quantiles <- melt(subset(species_data, groups=="humdom_wild" & N>1113 & N<1545 &
##   habitat=="none" & xi < 0.01),
##   id=c("N","xi","groups","habitat",hosts[-1]),
##   variable_name="grouped")
## svg("group_contributions_humdom_wild.svg")
## ggplot(group_quantiles, aes(grouped,
##   value))+geom_boxplot()+theme_bw(20)+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
##   axis.title.x=theme_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
##   scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
## dev.off()

## quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>1113 & N<1545 &
##   habitat=="none" & xi < 0.01),
##   id=c("N","xi","groups","habitat",domains[-1]),
##   variable_name="species")
## quantiles$species <- factor(quantiles$species, levels=hosts)

## svg("species_contributions_hum_dom_wild.svg")
## ggplot(quantiles, aes(species,
##   value))+geom_boxplot()+theme_bw(20)+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
##   axis.title.x=theme_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
##   scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
## dev.off()

## group_quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>1113 & N<1545 &
##   habitat=="none" & xi < 0.01),
##   id=c("N","xi","groups","habitat",hosts[-1]),
##   variable_name="grouped")
## svg("group_contributions_hum_dom_wild.svg")
## ggplot(group_quantiles, aes(grouped,
##   value))+geom_boxplot()+theme_bw(20)+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
##   axis.title.x=theme_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
##   scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
## dev.off()

## quantiles <- melt(subset(species_data, groups=="hum_domwild" & N>1113 & N<1545 &
##   habitat=="none" & xi < 0.01),
##   id=c("N","xi","groups","habitat",domains[-1]),
##   variable_name="species")
## quantiles$species <- factor(quantiles$species, levels=hosts)

## svg("species_contributions_hum_domwild.svg")
## ggplot(quantiles, aes(species,
##   value))+geom_boxplot()+theme_bw(20)+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
##   axis.title.x=theme_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
##   scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
## dev.off()

## group_quantiles <- melt(subset(species_data, groups=="hum_domwild" & N>1113 & N<1545 &
##   habitat=="none" & xi < 0.01),
##   id=c("N","xi","groups","habitat",hosts[-1]),
##   variable_name="grouped")

## svg("group_contributions_hum_domwild.svg")
## ggplot(group_quantiles, aes(grouped,
##   value))+geom_boxplot()+theme_bw(20)+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
##   axis.title.x=theme_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
##   scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
## dev.off()

## quantiles <- melt(subset(species_data, groups=="single" & N>1113 & N<1545 &
##   habitat=="binary" & xi > 10),
##   id=c("N","xi","groups","habitat",domains[-1]),
##   variable_name="species")
## quantiles$species <- factor(quantiles$species, levels=hosts)

## svg("species_contributions_binary.svg")
## ggplot(quantiles, aes(species,
##   value))+geom_boxplot()+theme_bw(20)+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
##   axis.title.x=theme_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
##   scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
## dev.off()

## group_quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>1113 & N<1545 &
##   habitat=="binary" & xi > 10),
##   id=c("N","xi","groups","habitat",hosts[-1]),
##   variable_name="grouped")

## svg("group_contributions_binary.svg")
## ggplot(group_quantiles, aes(grouped,
##   value))+geom_boxplot()+theme_bw(20)+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
##   axis.title.x=theme_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
##   scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
## dev.off()

## quantiles <- melt(subset(species_data, groups=="single" & N>1113 & N<1545 &
##   habitat=="fractional" & xi > 10),
##   id=c("N","xi","groups","habitat",domains[-1]),
##   variable_name="species")
## quantiles$species <- factor(quantiles$species, levels=hosts)

## svg("species_contributions_fractional.svg")
## ggplot(quantiles, aes(species,
##   value))+geom_boxplot()+theme_bw(20)+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
##   axis.title.x=theme_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]))+geom_hline(yintercept=1)+
##   scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
## dev.off()

group_quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>1113 & N<1545 &
  habitat=="fractional" & xi > 10),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")

svg("group_contributions_fractional.svg")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw(20)+scale_y_continuous(substitute("          Contribution to "*R[0]))+
  opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=theme_blank())+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()

group_quantiles <- melt(subset(species_data, groups=="single" & N>3000 & 
  habitat=="fractional" & xi > 10),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
pdf("group_contributions_single_fractional.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw(20)+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=theme_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]),
  limits=c(0,4))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()



