
setwd("/WGS_analysis_in_vivo/input")

library(ggplot2)
library(xlsx)
library(ggpubr)
library(dplyr)
library(colorBlindness)
library(car)

summary_patients_table <- read.xlsx("Supplementary_Table_6_summary_vcall_patients_expev.xlsx", 3, header = TRUE)
summary_patients_table$SNP_rate <- as.numeric(summary_patients_table$SNP_rate)
summary_patients_table$IS1_rate <- as.numeric(summary_patients_table$IS1_rate)
summary_patients_table$SNPs_INDELS <- as.numeric(summary_patients_table$SNPs_INDELS)
summary_patients_table$SNP_per_day <- as.numeric(summary_patients_table$SNP_per_day)
summary_patients_table$IS1_per_day <- as.numeric(summary_patients_table$IS1_per_day)
summary_patients_table$Adaptive_changes_per_day <- as.numeric(summary_patients_table$Adaptive_changes_per_day)
summary_patients_table$Adaptive_changes <- as.numeric(summary_patients_table$Adaptive_changes)

# Main Fig 5B

mypal <- c("#b8d7e7ff", "#daa3b8ff", "#a58ac8ff")

A <- summary_patients_table %>%
  filter(Species != "C. freundii") %>%
  ggplot(aes(x = Days_between_isolation, y = SNPs_INDELS, col = Species))+
  geom_point(aes(col = Species, stroke = 1.5, size = nIS1_genome), pch = 1)+
  facet_wrap(~factor(Species, levels = c("E. coli", "K. pneumoniae")), scales = "free_x")+
  xlab("Days")+
  ylab("SNP")+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1)+
  stat_cor(method = "pearson")+  theme_bw(base_size = 20)+
  ylim(0,9)+
  scale_y_continuous(breaks = integer_breaks())+
  scale_color_manual(values = mypal) +
  scale_fill_manual(values = mypal) +
  scale_shape_manual(values=c(1, 13)) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        #legend.position = "none",
        axis.title.x=element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank()); A

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

B <- summary_patients_table %>%
  filter(Strain != "K174", Strain != "K281") %>%
  ggplot(aes(x = nIS1_genome, y = IS1_per_day))+
  geom_point(aes(col = Species))+
  #facet_wrap(~factor(Species, levels = c("E. coli", "K. pneumoniae")), scales = "free_x")+
  #xlab("IS1 in the genome")+
  #ylab("IS1 per SNP")+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1)+
  stat_cor(method = "pearson")+  theme_bw(base_size = 20)+
  scale_color_manual(values = mypal) +
  scale_fill_manual(values = mypal) +
  scale_shape_manual(values=c(1, 13)) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank()); B

plist <- list(A,B)

ggarrange(plotlist = plist, nrow = 2, common.legend = TRUE, align = "hv")

library(mgcv)

summary_patients_table <- summary_patients_table %>%
  filter(Strain != "C184", Strain != "K281")

summary_patients_table$Species <- as.factor(summary_patients_table$Species)

gamk3 <- gam(IS1_per_day ~ s(nIS1_genome, k = 6) + Species, data = summary_patients_table, method = "REML")

summary(gamk3)

plot(gamk3, all.terms = TRUE, pages = 1, residuals = TRUE, pch = 1, cex = 1, shade = TRUE)

gam.check(gamk3)

