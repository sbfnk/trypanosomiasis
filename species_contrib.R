library(ggplot2)
library(reshape)
library(data.table)

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
pdf("stage2_species_contributions_binary.pdf")
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
pdf("stage2_group_contributions_binary.pdf")
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

group_quantiles <- melt(subset(species_data, groups=="single" & N>3000 & 
  habitat=="fractional" & xi > 10),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
(CIwild <- quantile(subset(group_quantiles, grouped=="domestic.wildlife")$value, probs=c(0.025, 0.31)))
svg("stage2_group_contributions_single_fractional.svg")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw(20)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16),
  axis.title.x=element_blank())+scale_y_continuous(substitute("          Contribution to "*R[0]),
  limits=c(0,4))+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()
rm(group_quantiles)
