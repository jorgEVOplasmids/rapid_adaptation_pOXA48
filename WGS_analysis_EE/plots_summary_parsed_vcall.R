
setwd("/home/jorge/Documents/important_docs/draft_expev/figures_draft_expev/summary_parsed_vcall_w_plasmids")

library(ggplot2)
library(xlsx)
library(dplyr)
library(janitor)
library(ggpubr)
library(wesanderson)
library(colorBlindness)
library(circlize)
library(BioCircos)
library(car)

# Define custom palette

#mypal <- c("#B30000", "#004E80", "#EEB609")

#mypal <- c("#737373", "red3", "#0066CC")

# Import data of SNPs, INDELs and NJs merged (sheet3)

merged_vcall <- read.xlsx("merged_summary_vcall_w_plasmids.xlsx", 3, header = TRUE)

# Plot frequencies (Supplementary Fig 4)
# 3 Species
freq_plot <- merged_vcall %>%
  filter(Event != "INDEL") %>%
  ggplot(aes(x = factor(Event, levels = c("SNP", "INDEL", "NJ")), y = Freq, group = Event, col = Condition, fill = Condition))+
  geom_boxplot(alpha = 0.2)+
  geom_jitter(size = 2, width = 0.1)+
  facet_wrap(~factor(Species, levels = c("E. coli", "C. freundii", "K. pneumoniae"))+Condition, nrow = 3)+
  xlab("Mutational event")+
  ylab("Frequency (%)")+
  theme_bw(base_size = 28)+
  ylim(0,120)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  #scale_color_manual(values = mypal) +
  #scale_fill_manual(values = mypal) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1)); freq_plot

cvdPlot(freq_plot)

# From sheet 3, get number of events in each condition

tabyl(merged_vcall %>% filter(Species == "E. coli"), Type, Condition)
tabyl(merged_vcall %>% filter(Species == "K. pneumoniae"), Type, Condition)
tabyl(merged_vcall %>% filter(Species == "C. freundii"), Type, Condition)

# Also, get the number of IS mediated NJ per replicate

merged_vcall$sprep <- as.vector(interaction(merged_vcall$Strain, merged_vcall$Replicate))
# E COLI
# IS mediated NJ
tabyl(merged_vcall %>% filter(Species == "E. coli", Condition == "No pOXA-48", IS.mediated != "-"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "E. coli", Condition == "pOXA-48", IS.mediated != "-"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "E. coli", Condition == "pOXA-48+AMC", IS.mediated != "-"), IS.mediated, sprep)
# Chr-Plasmid mediated NJ
tabyl(merged_vcall %>% filter(Species == "E. coli", Condition == "No pOXA-48", IS.mediated != "-", Type == "Chromosome-Plasmid"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "E. coli", Condition == "pOXA-48", IS.mediated != "-", Type == "Chromosome-Plasmid"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "E. coli", Condition == "pOXA-48+AMC", IS.mediated != "-", Type == "Chromosome-Plasmid"), IS.mediated, sprep)
# Chr-chr mediated NJ
tabyl(merged_vcall %>% filter(Species == "E. coli", Condition == "No pOXA-48", IS.mediated != "-", Type == "Chromosome-Chromosome"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "E. coli", Condition == "pOXA-48", IS.mediated != "-", Type == "Chromosome-Chromosome"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "E. coli", Condition == "pOXA-48+AMC", IS.mediated != "-", Type == "Chromosome-Chromosome"), IS.mediated, sprep)
# K PNEUMONIAE
# IS mediated NJ
tabyl(merged_vcall %>% filter(Species == "K. pneumoniae", Condition == "No pOXA-48", IS.mediated != "-"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "K. pneumoniae", Condition == "pOXA-48", IS.mediated != "-"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "K. pneumoniae", Condition == "pOXA-48+AMC", IS.mediated != "-"), IS.mediated, sprep)
# Chr-Plasmid mediated NJ
tabyl(merged_vcall %>% filter(Species == "K. pneumoniae", Condition == "No pOXA-48", IS.mediated != "-", Type == "Chromosome-Plasmid"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "K. pneumoniae", Condition == "pOXA-48", IS.mediated != "-", Type == "Chromosome-Plasmid"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "K. pneumoniae", Condition == "pOXA-48+AMC", IS.mediated != "-", Type == "Chromosome-Plasmid"), IS.mediated, sprep)
# Chr-chr mediated NJ
tabyl(merged_vcall %>% filter(Species == "K. pneumoniae", Condition == "No pOXA-48", IS.mediated != "-", Type == "Chromosome-Chromosome"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "K. pneumoniae", Condition == "pOXA-48", IS.mediated != "-", Type == "Chromosome-Chromosome"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "K. pneumoniae", Condition == "pOXA-48+AMC", IS.mediated != "-", Type == "Chromosome-Chromosome"), IS.mediated, sprep)
# C FREUNDII
# IS mediated NJ
tabyl(merged_vcall %>% filter(Species == "C. freundii", Condition == "No pOXA-48", IS.mediated != "-"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "C. freundii", Condition == "pOXA-48", IS.mediated != "-"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "C. freundii", Condition == "pOXA-48+AMC", IS.mediated != "-"), IS.mediated, sprep)
# Chr-Plasmid mediared NJ
tabyl(merged_vcall %>% filter(Species == "C. freundii", Condition == "No pOXA-48", IS.mediated != "-", Type == "Chromosome-Plasmid"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "C. freundii", Condition == "pOXA-48", IS.mediated != "-", Type == "Chromosome-Plasmid"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "C. freundii", Condition == "pOXA-48+AMC", IS.mediated != "-", Type == "Chromosome-Plasmid"), IS.mediated, sprep)
# Chr-chr mediated NJ
tabyl(merged_vcall %>% filter(Species == "C. freundii", Condition == "No pOXA-48", IS.mediated != "-", Type == "Chromosome-Chromosome"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "C. freundii", Condition == "pOXA-48", IS.mediated != "-", Type == "Chromosome-Chromosome"), IS.mediated, sprep)
tabyl(merged_vcall %>% filter(Species == "C. freundii", Condition == "pOXA-48+AMC", IS.mediated != "-", Type == "Chromosome-Chromosome"), IS.mediated, sprep)

# Import occurrences sheet in which data from previous code was stored

occurrences_vcall <- read.xlsx("merged_summary_vcall_w_plasmids.xlsx", 4, header = TRUE)



occurrences_vcall %>%
  #filter(Species == "E. coli") %>%
  ggplot(aes(y = Number, x= factor(Event, levels = c("SNP", "NJ")), group = Event, col = Type, fill = Type))+
  geom_bar(stat="identity", alpha = 0.8)+
  #geom_jitter(size = 2, width = 0.1)+
  facet_wrap(~factor(Species, levels = c("E. coli", "C. freundii", "K. pneumoniae"))+Condition, nrow = 3)+
  xlab("Mutational event")+
  ylab("Number of occurrences")+
  theme_bw(base_size = 28)+
  #ylim(0,110)+
  scale_color_brewer(palette = "Paired")+
  scale_fill_brewer(palette = "Paired")+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# STATS for frequency of occurrences

frequencies_vcall_coli <- merged_vcall %>%
  filter(Species == "E. coli", Event != "INDEL")

frequencies_vcall_coli$g <- interaction(frequencies_vcall_coli$Event, frequencies_vcall_coli$Condition)

# Check normality assumption

model  <- lm(Freq ~ g, data = frequencies_vcall_coli)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

bartlett.test(Freq ~ g, data = frequencies_vcall_coli)

# No normal data -> kruskal wallis test

freqcall.kt <- kruskal.test(Freq ~ g, data = frequencies_vcall_coli); freqcall.kt
summary(freqcall.kt)

# Post Hoc test

pairwise.wilcox.test(x = frequencies_vcall_coli$Freq, g = frequencies_vcall_coli$g, data = frequencies_vcall_coli, p.adjust.method = "holm")

nvcall.tukey <- TukeyHSD(freqcall.kt)
summary(nvcall.tukey)

# Represent the type of SNP or NJ depending on the experimental condition for each species

table_mainfig2 <- read.xlsx("table_mainfig2.xlsx", sheetIndex = 1, header = TRUE) # Summarized from complete summary table of EE

# Main Fig 2

A <- table_mainfig2 %>%
  filter(Type != "Pseudogene", Type != "Plasmid-Plasmid", Species == "E. coli", Event != "NJ") %>%
  ggplot(aes(y = Number, x= factor(Event, levels = c("SNP", "NJ", "IS-mediated NJ")), group = Event, col = Type, fill = Type))+
  geom_bar(stat="identity", alpha = 0.8)+
  #geom_jitter(size = 2, width = 0.1)+
  facet_wrap(~factor(Species, levels = c("E. coli", "C. freundii", "K. pneumoniae"))+Condition, nrow = 1)+
  xlab("Mutational event")+
  ylab("Number of occurrences")+
  theme_bw(base_size = 20)+
  #ylim(0,110)+
  #geom_vline() +
  scale_fill_manual(values = c("Intergenic" = "#e8fbfc",
                               "Non-synonymous" = "#0E7175",
                               "Synonymous" = "#8dedf1",
                               "Non-IS mediated" ="#fff2e6",
                               "Other"= "#febc80",
                               "IS1"="#FD7901"))+
  scale_color_manual(values = c(rep("black",7)))+
  #scale_fill_manual(values = mypal)+
  #scale_color_manual(values = mypal)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        #legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1.25),
        strip.background = element_blank());A

B <- table_mainfig2 %>%
  filter(Type != "Pseudogene", Type != "Plasmid-Plasmid", Species == "K. pneumoniae", Event != "NJ") %>%
  ggplot(aes(y = Number, x= factor(Event, levels = c("SNP", "NJ", "IS-mediated NJ")), group = Event, col = Type, fill = Type))+
  geom_bar(stat="identity", alpha = 0.8)+
  #geom_jitter(size = 2, width = 0.1)+
  facet_wrap(~factor(Species, levels = c("E. coli", "C. freundii", "K. pneumoniae"))+Condition, nrow = 1)+
  xlab("Mutational event")+
  ylab("Number of occurrences")+
  theme_bw(base_size = 20)+
  #ylim(0,110)+
  #geom_vline() +
  scale_fill_manual(values = c("Intergenic" = "#e8fbfc",
                               "Non-synonymous" = "#0E7175",
                               "Synonymous" = "#8dedf1",
                               "Non-IS mediated" ="#fff2e6",
                               "Other"= "#febc80",
                               "IS1"="#FD7901"))+
  scale_color_manual(values = c(rep("black",7)))+
  #scale_fill_manual(values = mypal)+
  #scale_color_manual(values = mypal)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        #legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1.25),
        strip.background = element_blank());B

C <- table_mainfig2 %>%
  filter(Type != "Pseudogene", Type != "Plasmid-Plasmid", Species == "C. freundii", Event != "NJ") %>%
  ggplot(aes(y = Number, x= factor(Event, levels = c("SNP", "NJ", "IS-mediated NJ")), group = Event, col = Type, fill = Type))+
  geom_bar(stat="identity", alpha = 0.8)+
  #geom_jitter(size = 2, width = 0.1)+
  facet_wrap(~factor(Species, levels = c("E. coli", "C. freundii", "K. pneumoniae"))+Condition, nrow = 1)+
  xlab("Mutational event")+
  ylab("Number of occurrences")+
  theme_bw(base_size = 20)+
  #ylim(0,110)+
  #geom_vline() +
  scale_fill_manual(values = c("Intergenic" = "#e8fbfc",
                               "Non-synonymous" = "#0E7175",
                               "Synonymous" = "#8dedf1",
                               "Non-IS mediated" ="#fff2e6",
                               "Other"= "#febc80",
                               "IS1"="#FD7901"))+
  scale_color_manual(values = c(rep("black",7)))+
  #scale_fill_manual(values = mypal)+
  #scale_color_manual(values = mypal)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        #legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1.25),
        strip.background = element_blank());C

ggarrange(A + rremove("ylab") + rremove("xlab"), B + rremove("ylab") + rremove("xlab"), C + rremove("ylab") + rremove("xlab"),
          common.legend = TRUE, nrow = 3, legend = "bottom", align = "hv")

############################################################################################

### REPRESENT IS MEDIATED NJ

# Supplementary Figure 5 

IS_complete_table <- read.xlsx("merged_summary_vcall_w_plasmids.xlsx", 5, header = TRUE)
# The 3 species
IS_complete_table %>%
  ggplot(aes(x = Condition, y = IS_mediated_rearrangements, fill = Condition, col = Condition))+
  geom_boxplot(alpha = 0.2)+
  #geom_jitter()+
  geom_point(size = 2, position = position_jitter(w = 0.3, h = 0), shape = 1, stroke = 1)+
  facet_wrap(~factor(Species, levels = c("E. coli", "K. pneumoniae", "C. freundii")), nrow = 1) +
  xlab("Condition")+
  ylab("IS-mediated NJ")+
  theme_bw(base_size = 26)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  ylim(-0.5, 11)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        #strip.text.x = element_blank(),
        #axis.text.y = element_blank(),
        legend.position = "none",
        strip.text = element_text(face = "italic"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

### STATS IS movements

IS_coli_table <- IS_complete_table %>% filter(Species == "E. coli")

# Check normality assumption

model  <- lm(IS_mediated_rearrangements ~ Condition, data = IS_coli_table)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(IS_mediated_rearrangements ~ Condition, data = IS_coli_table)

# No normal data -> kruskal wallis test

iscoli.kt <- kruskal.test(IS_mediated_rearrangements ~ Condition, data = IS_coli_table); iscoli.kt
summary(iscoli.kt)

# Post Hoc test

pairwise.wilcox.test(x = IS_coli_table$IS_mediated_rearrangements, g = IS_coli_table$Condition, data = IS_coli_table, p.adjust.method = "holm")

###

IS_cf_table <- IS_complete_table %>% filter(Species == "C. freundii")

# Check normality assumption

model  <- lm(IS_mediated_rearrangements ~ Condition, data = IS_cf_table)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(IS_mediated_rearrangements ~ Condition, data = IS_cf_table)

# No normal data -> kruskal wallis test

iscf.kt <- kruskal.test(IS_mediated_rearrangements ~ Condition, data = IS_cf_table); iscf.kt
summary(iscf.kt)

# Post Hoc test

pairwise.wilcox.test(x = IS_cf_table$IS_mediated_rearrangements, g = IS_cf_table$Condition, data = IS_cf_table, p.adjust.method = "holm")

###

IS_kpn_table <- IS_complete_table %>% filter(Species == "K. pneumoniae")

# Check normality assumption

model  <- lm(IS_mediated_rearrangements ~ Condition, data = IS_kpn_table)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

bartlett.test(IS_mediated_rearrangements ~ Condition, data = IS_kpn_table)
leveneTest(IS_mediated_rearrangements ~ Condition, data = IS_kpn_table)

# No normal data -> kruskal wallis test
# anova to check significance with a more powerfull test
iskpn.aov <- aov(IS_mediated_rearrangements ~ Condition, data = IS_kpn_table); iskpn.aov
summary(iskpn.aov)
TukeyHSD(iskpn.aov, "Condition")

# kruskal test to obtain comparable values with the other species
iskpn.kt <- kruskal.test(IS_mediated_rearrangements ~ Condition, data = IS_kpn_table); iskpn.kt
summary(iskpn.kt)

# Post Hoc test

pairwise.wilcox.test(x = IS_kpn_table$IS_mediated_rearrangements, g = IS_kpn_table$Condition, data = IS_kpn_table, p.adjust.method = "holm")

######################################################################################################

### CHROMOSOME CIRCAPLOT

## For main figure 3 A-B

# Import replicate names for the strain from the parsed xlsx file

xlsx_parsed <- read.xlsx("/home/jorge/Documents/TFM/clinical_samples/R_analysis/xlsx_breseqs/parsed_mutations/110822/parsed_K209.xlsx", sheetIndex = 2, header = FALSE)
xlsx_parsed <- xlsx_parsed %>% select(X9,X10,X11,X13,X14,X15,X16,X17,X18)

repl_name <- as.vector(t(xlsx_parsed[1,]))
repl_number <- c(1:3, 1:3, 1:3)
repl_condition <- c(rep("pOXA-48",3), rep("No pOXA-48",3), rep("pOXA-48+AMC",3))
repl_keydf <- data.frame(cond_key = interaction(repl_condition,repl_number), Replicate = repl_name, Rep_num = repl_number)

# Import each strain mutations for circular plotting

strain_summary_circa <- merged_vcall %>%
  filter(Strain == "K209", Type != "Plasmid-Plasmid", Seq.ID1 == "Chromosome") %>%
  select(Position1, Replicate, Type, Condition, Seq.ID1, Gene1, Annotation1, IS.mediated)

#strain_summary_circa$Seq.ID <- "Chromosome"
strain_summary_circa$end <- strain_summary_circa$Position+1
strain_summary_circa <- merge(strain_summary_circa, repl_keydf, by = "Replicate")
strain_summary_circa <- strain_summary_circa %>% relocate(Seq.ID1, Position1, end, Rep_num, Type, Condition, IS.mediated)

# Define chromosome for plotting

len_chr <- read.xlsx("/home/jorge/Documents/TFM/clinical_samples/R_analysis/xlsx_breseqs/parsed_mutations/110822/parsed_K209.xlsx", sheetIndex = 1, header = TRUE)[1,2]

chr_data <- data.frame("Chromosome", 0, len_chr)
colnames(chr_data) <- c("Region", "Begin", "End")

# Define lines for ploting each replicate

repl_lines <- rbind(chr_data, chr_data, chr_data)
repl_lines$Rep_num <- c(1,2,3)

# Begin plot with circlize
circos.clear()
circos.par(start.degree = 90)
circos.genomicInitialize(chr_data)
#unique_labels <- strain_summary_circa[!duplicated(strain_summary_circa$Gene), ]
#circos.labels(sectors = c(rep("Chromosome", length(unique_labels$Position))), x = unique_labels$Position, labels = unique_labels$Gene, side = "outside", cex = 0.6)
circos.genomicTrack(data = strain_summary_circa %>% filter(Condition == "pOXA-48+AMC"),
                    cell.padding = c(0.05, 1.00, 0.05, 1.00),
                    numeric.column = c("Rep_num"),
                    panel.fun = function(region, value,...) {
                      #circos.genomicPoints(region, value,...)
                      #circos.genomicLines(region, value, type = "segment",...)
                    })
#circos.rect(track.index = 2,col = "#4cae4acc")
y <- seq(1,3, by = 1)
circos.segments(0,y,len_chr,y, track.index = 2, lwd = 3.5, col = "#4DAF4A")
# Draw non-synonymus SNPs as orange points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "non-synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "non-synonymous") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 21, cex = 1.5, bg = "black")
# Draw synonymous SNPs as empty points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "synonymous") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 1, cex = 1.5, lwd = 2)
# Draw Intergenic SNPs as grey points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "intergenic") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "intergenic") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 21, cex = 1.5, bg = "grey")
# Draw Chr-Chr NJ as empty diamond
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Chromosome") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Chromosome") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 2, cex = 1.5, lwd = 5)
# Draw Chr-Plasmid NJ as a star
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Plasmid") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Plasmid") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 5, cex = 1.5, lwd = 5)
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Chromosome", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Chromosome", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 24, cex = 1.5, lwd = 5, bg = "orange")
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 23, cex = 1.5, lwd = 5, bg = "orange")
#circos.genomicLines(region = repl_lines %>% select(Begin, End), value = repl_lines %>% select(Rep_num),
#numeric.column = c("Rep_num"), track.index = 2, type = "segment")
circos.genomicTrack(data = strain_summary_circa %>% filter(Condition == "pOXA-48"),
                    cell.padding = c(0.05, 1.00, 0.05, 1.00),
                    numeric.column = c("Rep_num"),
                    panel.fun = function(region, value,...) {
                      #circos.genomicPoints(region, value,...)
                      #circos.genomicLines(region, value, type = "segment",...)
                    })
#circos.rect(track.index = 2,col = "#4cae4acc")
y <- seq(1,3, by = 1)
circos.segments(0,y,len_chr,y, track.index = 3, lwd = 3.5, col = "#377EB8")
# Draw non-synonymus SNPs as empty points
circos.genomicPoints(region = strain_summary_circa %>% 
                       filter(Condition == "pOXA-48", Type == "non-synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "non-synonymous") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 21, cex = 1.5, bg = "black")
# Draw synonymous SNPs as empty points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "synonymous") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 1, cex = 1.5, lwd = 2)
# Draw Intergenic SNPs as grey points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "intergenic") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "intergenic") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 21, cex = 1.5, bg = "grey")
# Draw Chr-Chr NJ as diamond squares
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Chromosome") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Chromosome") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 2, cex = 1.5, lwd = 5)
# Draw Chr-Plasmid NJ as a star
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Plasmid") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Plasmid") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 5, cex = 1.5, lwd = 5)
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Chromosome", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Chromosome", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 24, cex = 1.5, lwd = 5, bg = "orange")
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 23, cex = 1.5, lwd = 5, bg = "orange")
circos.genomicTrack(data = strain_summary_circa %>% filter(Condition == "No pOXA-48"),
                    cell.padding = c(0.05, 1.00, 0.05, 1.00),
                    numeric.column = c("Rep_num"),
                    panel.fun = function(region, value,...) {
                      #circos.genomicPoints(region, value,...)
                      #circos.genomicLines(region, value, type = "segment",...)
                    })

#circos.rect(track.index = 2,col = "#4cae4acc")
y <- seq(1,3, by = 1)
circos.segments(0,y,len_chr,y, track.index = 4, lwd = 3.5, col = "#E41A1C")
# Draw non-synonymus SNPs black points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "non-synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "non-synonymous") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 21, cex = 1.5, bg = "black")
# Draw synonymous SNPs as empty points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "synonymous") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 1, cex = 1.5, lwd = 2)
# Draw Intergenic SNPs as grey points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "intergenic") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "intergenic") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 21, cex = 1.5, bg = "grey")
# Draw Chr-Chr NJ as crossed square
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Chromosome") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Chromosome") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 2, cex = 1.5, lwd = 5)
# Draw Chr-Plasmid NJ as a star
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Plasmid") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Plasmid") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 5, cex = 1.5, lwd = 5)
# Draw IS-mediated events as orange diamonds
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Chomosome", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Chomosome", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 24, cex = 1.5, lwd = 5, bg = "orange")
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 23, cex = 1.5, lwd = 5, bg = "orange")
circos.clear()

# Labels with gene names
circos.par(start.degree = 90)
circos.genomicInitialize(chr_data, plotType = NULL)
unique_labels <- strain_summary_circa[!duplicated(strain_summary_circa$Gene), ]
circos.labels(sectors = c(rep("Chromosome", length(unique_labels$Position1))), x = unique_labels$Position1, labels = unique_labels$Gene1, side = "outside", cex = 0.6)
circos.clear()
# Labels with annotation of genes
circos.par(start.degree = 90)
circos.genomicInitialize(chr_data, plotType = NULL)
unique_labels <- strain_summary_circa[!duplicated(strain_summary_circa$Gene), ]
circos.labels(sectors = c(rep("Chromosome", length(unique_labels$Position1))), x = unique_labels$Position1, labels = unique_labels$Annotation1, side = "outside", cex = 0.6)
circos.clear()

###########################################################################

# PLASMIDS CIRCAPLOT

# Import replicate names for the strain from the parsed xlsx file

# Import replicate names for the strain from the parsed xlsx file

xlsx_parsed <- read.xlsx("/home/jorge/Documents/TFM/clinical_samples/R_analysis/xlsx_breseqs/parsed_mutations/110822/parsed_K163.xlsx", sheetIndex = 2, header = FALSE)
xlsx_parsed <- xlsx_parsed %>% select(X9,X10,X11,X13,X14,X15,X16,X17,X18)

repl_name <- as.vector(t(xlsx_parsed[1,]))
repl_number <- c(1:3, 1:3, 1:3)
repl_condition <- c(rep("pOXA-48",3), rep("No pOXA-48",3), rep("pOXA-48+AMC",3))
repl_keydf <- data.frame(cond_key = interaction(repl_condition,repl_number), Replicate = repl_name, Rep_num = repl_number)

plasmid_name <- "ColRNAI_1"

# Import each strain mutations for circular plotting
# Get mutations in the plasmid
strain_summary_circa <- merged_vcall %>%
  filter(Strain == "K163", Seq.ID1 == plasmid_name) %>%
  select(Position1, Replicate, Type, Condition, Seq.ID1, Gene1, Annotation1, IS.mediated)
# Get NJs in the plasmid if it appears in the position 2 of the NJ event
strain_summary_circa_NJs <- merged_vcall %>%
  filter(Strain == "K163", Seq.ID2 == plasmid_name) %>%
  select(Position2, Replicate, Type, Condition, Seq.ID2, Gene2, Annotation2, IS.mediated)

colnames(strain_summary_circa_NJs) <- colnames(strain_summary_circa)

strain_summary_circa <- rbind(strain_summary_circa, strain_summary_circa_NJs)

#strain_summary_circa$Seq.ID <- "Chromosome"
strain_summary_circa$Position1 <- as.numeric(strain_summary_circa$Position1)
strain_summary_circa$end <- strain_summary_circa$Position1+1
strain_summary_circa <- merge(strain_summary_circa, repl_keydf, by = "Replicate")
strain_summary_circa <- strain_summary_circa %>% relocate(Seq.ID1, Position1, end, Rep_num, Type, Condition, IS.mediated)

# Define chromosome for plotting

len_chr <- read.xlsx("/home/jorge/Documents/TFM/clinical_samples/R_analysis/xlsx_breseqs/parsed_mutations/110822/parsed_K163.xlsx", sheetIndex = 1, header = TRUE)[4,2]

chr_data <- data.frame(plasmid_name, 0, len_chr)
colnames(chr_data) <- c("Region", "Begin", "End")

# Define lines for ploting each replicate

repl_lines <- rbind(chr_data, chr_data, chr_data)
repl_lines$Rep_num <- c(1,2,3)

# Begin plot with circlize
circos.clear()
circos.par(start.degree = 90)
circos.genomicInitialize(chr_data)
#unique_labels <- strain_summary_circa[!duplicated(strain_summary_circa$Gene), ]
#circos.labels(sectors = c(rep("Chromosome", length(unique_labels$Position))), x = unique_labels$Position, labels = unique_labels$Gene, side = "outside", cex = 0.6)
circos.genomicTrack(data = strain_summary_circa,
                    cell.padding = c(0.05, 1.00, 0.05, 1.00),
                    numeric.column = c("Rep_num"),
                    panel.fun = function(region, value,...) {
                      #circos.genomicPoints(region, value,...)
                      #circos.genomicLines(region, value, type = "segment",...)
                    })
#circos.rect(track.index = 2,col = "#4cae4acc")
y <- seq(1,3, by = 1)
circos.segments(0,y,len_chr,y, track.index = 2, lwd = 3.5, col = "#4DAF4A")
# Draw non-synonymus SNPs as orange points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "non-synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "non-synonymous") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 21, cex = 1.5, bg = "black")
# Draw synonymous SNPs as empty points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "synonymous") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 1, cex = 1.5, lwd = 2)
# Draw Intergenic SNPs as grey points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "intergenic") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "intergenic") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 21, cex = 1.5, bg = "grey")
# Draw Chr-Chr NJ as empty diamond
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Plasmid-Plasmid") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Plasmid-Plasmid") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 2, cex = 1.5, lwd = 5)
# Draw Chr-Plasmid NJ as a star
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Plasmid") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Plasmid") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 5, cex = 1.5, lwd = 5)
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Plasmid-Plasmid", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Plasmid-Plasmid", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 24, cex = 1.5, lwd = 5, bg = "orange")
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48+AMC", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 2,
                     pch = 23, cex = 1.5, lwd = 5, bg = "orange")
#circos.genomicLines(region = repl_lines %>% select(Begin, End), value = repl_lines %>% select(Rep_num),
#numeric.column = c("Rep_num"), track.index = 2, type = "segment")
circos.genomicTrack(data = strain_summary_circa,
                    cell.padding = c(0.05, 1.00, 0.05, 1.00),
                    numeric.column = c("Rep_num"),
                    panel.fun = function(region, value,...) {
                      #circos.genomicPoints(region, value,...)
                      #circos.genomicLines(region, value, type = "segment",...)
                    })
#circos.rect(track.index = 2,col = "#4cae4acc")
y <- seq(1,3, by = 1)
circos.segments(0,y,len_chr,y, track.index = 3, lwd = 3.5, col = "#377EB8")
# Draw non-synonymus SNPs as empty points
circos.genomicPoints(region = strain_summary_circa %>% 
                       filter(Condition == "pOXA-48", Type == "non-synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "non-synonymous") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 21, cex = 1.5, bg = "black")
# Draw synonymous SNPs as empty points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "synonymous") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 1, cex = 1.5, lwd = 2)
# Draw Intergenic SNPs as grey points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "intergenic") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "intergenic") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 21, cex = 1.5, bg = "grey")
# Draw Chr-Chr NJ as diamond squares
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Plasmid-Plasmid") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Plasmid-Plasmid") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 2, cex = 1.5, lwd = 5)
# Draw Chr-Plasmid NJ as a star
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Plasmid") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Plasmid") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 5, cex = 1.5, lwd = 5)
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Plasmid-Plasmid", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Plasmid-Plasmid", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 24, cex = 1.5, lwd = 5, bg = "orange")
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "pOXA-48", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 3,
                     pch = 23, cex = 1.5, lwd = 5, bg = "orange")
circos.genomicTrack(data = strain_summary_circa,
                    cell.padding = c(0.05, 1.00, 0.05, 1.00),
                    numeric.column = c("Rep_num"),
                    panel.fun = function(region, value,...) {
                      #circos.genomicPoints(region, value,...)
                      #circos.genomicLines(region, value, type = "segment",...)
                    })

#circos.rect(track.index = 2,col = "#4cae4acc")
y <- seq(1,3, by = 1)
circos.segments(0,y,len_chr,y, track.index = 4, lwd = 3.5, col = "#E41A1C")
# Draw non-synonymus SNPs black points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "non-synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "non-synonymous") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 21, cex = 1.5, bg = "black")
# Draw synonymous SNPs as empty points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "synonymous") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "synonymous") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 1, cex = 1.5, lwd = 2)
# Draw Intergenic SNPs as grey points
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "intergenic") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "intergenic") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 21, cex = 1.5, bg = "grey")
# Draw Chr-Chr NJ as crossed square
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Plasmid-Plasmid") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Plasmid-Plasmid") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 2, cex = 1.5, lwd = 5)
# Draw Chr-Plasmid NJ as a star
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Plasmid") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Plasmid") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 5, cex = 1.5, lwd = 5)
# Draw IS-mediated events as orange diamonds
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Chomosome", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Chomosome", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 24, cex = 1.5, lwd = 5, bg = "orange")
circos.genomicPoints(region = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Position1, end),
                     value = strain_summary_circa %>%
                       filter(Condition == "No pOXA-48", Type == "Chromosome-Plasmid", IS.mediated == "yes") %>%
                       select(Rep_num),
                     track.index = 4,
                     pch = 23, cex = 1.5, lwd = 5, bg = "orange")
circos.clear()

# Labels with gene names
circos.par(start.degree = 90)
circos.genomicInitialize(chr_data, plotType = NULL)
unique_labels <- strain_summary_circa[!duplicated(strain_summary_circa$Gene), ]
circos.labels(sectors = c(rep(plasmid_name, length(unique_labels$Position1))), x = unique_labels$Position1, labels = unique_labels$Gene1, side = "outside", cex = 0.6)
circos.clear()
# Labels with annotation of genes
circos.par(start.degree = 90)
circos.genomicInitialize(chr_data, plotType = NULL)
unique_labels <- strain_summary_circa[!duplicated(strain_summary_circa$Gene), ]
circos.labels(sectors = c(rep(plasmid_name, length(unique_labels$Position1))), x = unique_labels$Position1, labels = unique_labels$Annotation1, side = "outside", cex = 0.6)
circos.clear()
