library('ggplot2')
library('cowplot')
library('scales')
library('tidyverse')
library('magrittr')
library('stringr')

lambda <- read_csv("Passive/resultslambda.csv") %>%
  mutate(project=ifelse(project == "Arua", "Arua-Yumbe", project))

project_n <- lambda %>%
  count(project) %>%
  rename(project_n=n)

lambda_hist <- lambda %>%
  mutate(lambda_ML_range=cut(lambda_ML, seq(0, max(lambda_ML)), include.lowest=TRUE)) %>%
  mutate(lower.age.limit=as.integer(sub("^.([0-9]*),.*$","\\1", lambda_ML_range))) %>%
  mutate(lambda_ML_range=
           factor(paste(lower.age.limit, lower.age.limit+0.9, sep="-"),
                  levels=paste(seq_len(max(lower.age.limit)+1)-1,
                               seq_len(max(lower.age.limit)+1)-0.1,
                               sep="-"))) %>%
  count(project, lambda_ML_range) %>%
  left_join(project_n) %>%
  mutate(prop=n/project_n) %>%
  mutate(lambda_ML_range)

p <- ggplot(lambda_hist, aes(x=lambda_ML_range, y=prop, fill=project)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_brewer("Project", palette="Dark2") +
  scale_y_continuous("Proportion of villages", labels=percent) +
  scale_x_discrete(expression(atop(Incidence~(lambda[i])~"in cases",
                                   per~"1,000"~susceptible~"person-months")),
                   drop=FALSE) +
  theme(legend.justification=c(1, 1), legend.position=c(1, 1),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("Passive/foi_dist.pdf", p)

lambda %<>%
  arrange(project, lambda_ML) %>%
  mutate(id=1:n())

p <- ggplot(lambda, aes(x=id)) +
  geom_point(aes(y=lambda_median), shape=4, color="black", alpha=0.5) +
  geom_point(aes(y=lambda_ML), shape=21, fill="grey", color="black") +
  geom_errorbar(aes(ymin=lambda_lower, ymax=lambda_upper), color="black", alpha=0.25) +
  coord_cartesian(ylim=c(0, 15)) +
  facet_grid(~project, scales="free_x", space="free_x") +
  scale_x_continuous("Village") +
  scale_y_continuous(
    expression(atop(Incidence~(lambda[i])~"in cases",
                    per~"1,000"~susceptible~"person-months"))) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("Passive/foi_villages.pdf", p, width=8.2, height=3.8)
