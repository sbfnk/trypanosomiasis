library(ggplot2)

species_labels=c("Human", "Goat", "Pig", "Sheep", "Blackstriped duiker",
  "Blue duiker", "Brush-tailed porcupine", "Giant rat",
  "G. white-nosed monkey", "Small-spotted genet", "Two-spotted palm civet",
  "White-eyelid mangabey")
domain_labels=c("Human", "Human+domestic", "Human+wildlife", "Domestic",
  "Domestic+wildlife", "Wildlife", expression(R[0]))

quantiles <- melt(subset(species_data, groups=="random" & N>3000),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)
(CIR0 <- quantile(quantiles$R0, probs=c(0.025, 0.975)))
pdf("species_contributions.pdf")
ggplot(quantiles, aes(species,
  value))+geom_boxplot()+theme_bw()+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1),
  axis.title.x=theme_blank())+scale_y_continuous(title="contribution")+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()

group_quantiles <- melt(subset(species_data, groups=="random" & N>3000),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
(CIHuman <- quantile(group_quantiles$Human, probs=c(0.025, 0.975)))
pdf("group_contributions.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw()+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1),
  axis.title.x=theme_blank())+scale_y_continuous(title="contribution")+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()

quantiles <- melt(subset(species_data, groups=="humdom_wild" & N>3000 &
  habitat=="none" & xi < 0.01),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)

pdf("species_contributions_humdom_wild.pdf")
ggplot(quantiles, aes(species,
  value))+geom_boxplot()+theme_bw()+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1),
  axis.title.x=theme_blank())+scale_y_continuous(title="contribution")+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()

group_quantiles <- melt(subset(species_data, groups=="humdom_wild" & N>3000 &
  habitat=="none" & xi < 0.01),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
pdf("group_contributions_humdom_wild.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw()+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1),
  axis.title.x=theme_blank())+scale_y_continuous(title="contribution")+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()

quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>3000 &
  habitat=="none" & xi < 0.01),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)

pdf("species_contributions_hum_dom_wild.pdf")
ggplot(quantiles, aes(species,
  value))+geom_boxplot()+theme_bw()+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1),
  axis.title.x=theme_blank())+scale_y_continuous(title="contribution")+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()

group_quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>3000 &
  habitat=="none" & xi < 0.01),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")
pdf("group_contributions_hum_dom_wild.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw()+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1),
  axis.title.x=theme_blank())+scale_y_continuous(title="contribution")+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()

quantiles <- melt(subset(species_data, groups=="single" & N>3000 &
  habitat=="binary" & xi > 10),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)

pdf("species_contributions_binary.pdf")
ggplot(quantiles, aes(species,
  value))+geom_boxplot()+theme_bw()+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1),
  axis.title.x=theme_blank())+scale_y_continuous(title="contribution")+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()

group_quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>3000 &
  habitat=="binary" & xi > 10),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")

pdf("group_contributions_binary.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw()+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1),
  axis.title.x=theme_blank())+scale_y_continuous(title="contribution")+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()

quantiles <- melt(subset(species_data, groups=="single" & N>3000 &
  habitat=="fractional" & xi > 10),
  id=c("N","xi","groups","habitat",domains[-1]),
  variable_name="species")
quantiles$species <- factor(quantiles$species, levels=hosts)

pdf("species_contributions_fractional.pdf")
ggplot(quantiles, aes(species,
  value))+geom_boxplot()+theme_bw()+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1),
  axis.title.x=theme_blank())+scale_y_continuous(title="contribution")+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(quantiles$species),labels=species_labels)
dev.off()

group_quantiles <- melt(subset(species_data, groups=="hum_dom_wild" & N>3000 &
  habitat=="fractional" & xi > 10),
  id=c("N","xi","groups","habitat",hosts[-1]),
  variable_name="grouped")

pdf("group_contributions_fractional.pdf")
ggplot(group_quantiles, aes(grouped,
  value))+geom_boxplot()+theme_bw()+opts(axis.text.x=theme_text(angle=45,hjust=1,vjust=1),
  axis.title.x=theme_blank())+scale_y_continuous(title="contribution")+geom_hline(yintercept=1)+
  scale_x_discrete(breaks=levels(group_quantiles$grouped),labels=domain_labels)
dev.off()



