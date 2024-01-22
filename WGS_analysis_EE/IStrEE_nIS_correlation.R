
# Definitive script for correlation of initial number of IS1 in the genome
# and the number of IS1 movements when the plasmid pOXA-48 is present during
# the experimental evolution of clinical enterobacteria

setwd("/home/jorge/Documents/important_docs/draft_expev/figures_draft_expev/summary_parsed_vcall_w_plasmids/")

library(ggplot2)
library(xlsx)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(lme4)

# Import data

#merged_summary_vcall_w_plasmids <- read.xlsx("merged_summary_vcall_w_plasmids.xlsx", sheetName = "ISs_per_replicate")
merged_summary_vcall_w_plasmids <- read.xlsx("merged_summary_vcall_w_plasmids.xlsx", sheetName = "merged_summary")
initial_n_IS1 <- read.xlsx("/home/jorge/Documents/important_docs/draft_expev/results_summary/tracking_ISs/IS1_elements_genomes.xlsx", sheetIndex = 1)

# Recalculate number of IS jumps per replicate as we only want IS1 movements in this correlation

is1_table <- subset(merged_summary_vcall_w_plasmids, grepl("IS1\\b|IS1-like\\b", Annotation2))
#is1_table <- is1_table %>% filter(Type == "Chromosome-Plasmid")
is1_table <- rbind(is1_table,subset(merged_summary_vcall_w_plasmids, grepl("IS1\\b|IS1-like\\b", Annotation1)))

summary_table_is1 <- aggregate(is1_table$IS.mediated, list(is1_table$Strain, is1_table$Replicate, is1_table$Condition), FUN = length)
colnames(summary_table_is1) <- c("Strain", "Replicate", "Condition", "IS1_movs")

# Sum IS1 movs per strain

summary_table_is1 <- aggregate(summary_table_is1$IS1_movs, list(summary_table_is1$Strain, summary_table_is1$Condition), FUN = sum)

colnames(summary_table_is1) <- c("Strain", "Condition", "IS1_movs")

#table_corr_nIS_movs_IS$Species <- c(rep("E. coli", 12), rep("C. freundii", 6), rep("K. pneumoniae", 21))

is1_table <- is1_table %>%
  select(Species, Strain, Condition, Replicate)

table_corr_nIS_movs_IS <- merge(summary_table_is1, initial_n_IS1, by = "Strain")

# Append samples with no movements at all

table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("C309", "No pOXA-48", 0, "E. coli", 26, 0, 0, 26))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("C309", "pOXA-48", 0, "E. coli", 26, 0, 0, 26))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("C309", "pOXA-48+AMC", 0, "E. coli", 26, 0, 0, 26))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("CF12", "No pOXA-48", 0, "C. freundii", 0, 1, 1, 1))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("CF13", "No pOXA-48", 0, "C. freundii", 1, 2, 1, 2))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("K091", "No pOXA-48", 0, "K. pneumoniae", 0, 3, 3, 3))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("K163", "No pOXA-48", 0, "K. pneumoniae", 0, 5, 5, 5))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("K209", "No pOXA-48", 0, "K. pneumoniae", 0, 1, 1, 1))


table_corr_nIS_movs_IS$IS1_movs <- as.numeric(table_corr_nIS_movs_IS$IS1_movs)
table_corr_nIS_movs_IS$n_IS1_chromosome <- as.numeric(table_corr_nIS_movs_IS$n_IS1_chromosome)
table_corr_nIS_movs_IS$n_IS1_total_genome <- as.numeric(table_corr_nIS_movs_IS$n_IS1_total_genome)

table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("C309", "No pOXA-48", 0, "E. coli", 26, 0, 0, 26))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("C309", "pOXA-48", 0, "E. coli", 26, 0, 0, 26))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("C309", "pOXA-48+AMC", 0, "E. coli", 26, 0, 0, 26))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("CF12", "No pOXA-48", 0, "C. freundii", 0, 1, 1, 1))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("CF13", "No pOXA-48", 0, "C. freundii", 1, 2, 1, 2))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("K091", "No pOXA-48", 0, "K. pneumoniae", 0, 3, 3, 3))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("K163", "No pOXA-48", 0, "K. pneumoniae", 0, 5, 5, 5))
table_corr_nIS_movs_IS <- rbind(table_corr_nIS_movs_IS, c("K209", "No pOXA-48", 0, "K. pneumoniae", 0, 1, 1, 1))

# Introduce by hand replicates which show 0 movements
#write.xlsx(table_corr_nIS_movs_IS, "table_corr_nIS_movs_IS.xlsx")
# Reload table
table_corr_nIS_movs_IS_full <- read.xlsx("/home/jorge/Documents/important_docs/draft_expev/results_summary/tracking_ISs/table_corr_nIS_movs_IS.xlsx", sheetIndex = 1)

table_corr_nIS_movs_IS_full_by_strain <- aggregate(table_corr_nIS_movs_IS_full$IS1_movs, by = list(table_corr_nIS_movs_IS_full$Strain, table_corr_nIS_movs_IS_full$Condition), FUN = sum)
colnames(table_corr_nIS_movs_IS_full_by_strain) <- c("Strain", "Condition", "IS1_movs")
table_corr_nIS_movs_IS_full_by_strain <- merge(table_corr_nIS_movs_IS_full_by_strain, initial_n_IS1, by = "Strain")

### Calculate increment of IS1 movements for both pOXA-48 and pOXA-48+AMC vs No pOXA-48

table_corr_nIS_movs_IS_full_by_strain_nopoxa <- table_corr_nIS_movs_IS_full_by_strain %>%
  filter(Condition == "No pOXA-48")

table_corr_nIS_movs_IS_full_by_strain_poxa <- table_corr_nIS_movs_IS_full_by_strain %>%
  filter(Condition == "pOXA-48")

table_corr_nIS_movs_IS_full_by_strain_poxa_AMC <- table_corr_nIS_movs_IS_full_by_strain %>%
  filter(Condition == "pOXA-48+AMC")

table_corr_nIS_movs_IS_full_by_strain_poxa$inc_movs_IS1 <- table_corr_nIS_movs_IS_full_by_strain_poxa$IS1_movs - table_corr_nIS_movs_IS_full_by_strain_nopoxa$IS1_movs
table_corr_nIS_movs_IS_full_by_strain_poxa_AMC$inc_movs_IS1 <- table_corr_nIS_movs_IS_full_by_strain_poxa_AMC$IS1_movs - table_corr_nIS_movs_IS_full_by_strain_nopoxa$IS1_movs

table_corr_w_increments <- rbind(table_corr_nIS_movs_IS_full_by_strain_poxa, table_corr_nIS_movs_IS_full_by_strain_poxa_AMC)

### PLOTS

# Plot the increment for all the strains when the plasmid is present independently of the AMC presence
mypal <- c("#b8d7e7ff", "#daa3b8ff", "#a58ac8ff")

C <- table_corr_w_increments %>%
  #filter(Strain != "K153") %>%
  ggplot(aes(x = n_IS1_total_genome, y = inc_movs_IS1)) +
  geom_point(aes(col = Species, shape = Condition), alpha = 0.7, size = 4) +
  scale_shape_manual(values = c(16,1)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.1)+
  #facet_wrap(~Condition) +
  stat_cor(method = "pearson")+
  ylab("Increment of IS1-mediated NJ") +
  xlab("Initial number of IS1 in the genome") +
  scale_color_manual(values = mypal) +
  theme_bw(base_size = 14)+
  ylim(-6, 20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        legend.position = "bottom",
        strip.background = element_blank());C

