
setwd("~/Documents/important_docs/draft_expev/results_summary")

library(xlsx)
library(ggplot2)
library(dplyr)

table_competis <- read.xlsx("competition_capsules.xlsx", sheetIndex = 3, header = TRUE)
table_points <- read.xlsx("competition_capsules.xlsx", sheetIndex = 4, header = TRUE)

table_competis %>%
  ggplot(aes(x = as.factor(Sample), y = Median-1, color = Condition, fill = Condition)) +
  geom_bar( aes(x=as.factor(Sample), y=Median-1), stat="identity", alpha=0.2) +
  geom_errorbar( aes(x=as.factor(Sample), ymin=Median-1-Sterr, ymax=Median-1+Sterr), width=0, alpha=0.9, size=1.3, col = "darkgrey")+
  theme_bw(base_size = 20) +
  geom_hline(yintercept = 0) +
  ylim(-0.5, 0.9)+
  ylab("Relative fitness (Non-capsulated/Capsulated)") +
  xlab("Clones") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.background = element_blank(), panel.grid = element_blank(), 
          #aspect.ratio = 1, 
          strip.text.x = element_blank(),
          legend.position = "bottom",
          strip.background = element_blank())

table_points %>%
  ggplot(aes(x = as.factor(Sample), y = Rel..growth-1, color = Condition, fill = Condition)) +
  geom_jitter(size = 4)+
  theme_bw(base_size = 20) +
  geom_hline(yintercept = 0) +
  ylim(-0.5, 0.9)+
  ylab("Relative fitness (Non-capsulated/Capsulated)") +
  xlab("Clones") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        strip.text.x = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank())
