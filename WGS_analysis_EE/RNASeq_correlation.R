# Building table
corr_table_fc <- as.data.frame(matrix(nrow=10,ncol=3))
colnames(corr_table_fc) <- c("Strain", "log2FC", "numIS1")
corr_table_fc$Strain <- c("C325", "CF13", "EC10", "H53", "K091", "K141", "K147", "K209", "K249", "K275")
corr_table_fc$Species <- c("E. coli", "C. freundii", "E. coli", rep("K. pneumoniae", 7))

# IS1 numbers are obtained from annotated IS1s, multiplying by estimated PCN
corr_table_fc$numIS1 <- c(12, 5, 2, 8, 6, 66, 7, 22, 8, 8)

# Log2FC are obtained from the DESeq2 raw output
corr_table_fc$log2FC <- c(0.336947239785708, 0.809171651123054, 3.22684752401752, 0.761776685244606, 1.65416578660558, -0.298940862287361, 0.334224576505695, -0.235702742630926, 0.975609714589933, 1.04471306249639)

# Normality test
shapiro.test(corr_table_fc$log2FC)
shapiro.test(corr_table_fc$numIS1)

# Correlation plot
ggscatter(corr_table_fc, x = "numIS1", y = "log2FC",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "number of IS1 copies", ylab = "log2FC",
          add.params = list(color = "black", fill = "lightgray"),
          size = 4, cor.coef.size = 6, color="Strain", alpha=0.5)

mypal <- c("#b8d7e7ff", "#daa3b8ff", "#a58ac8ff")

D <- corr_table_fc %>%
  #filter(Strain != "K153") %>%
  ggplot(aes(x = numIS1, y = log2FC)) +
  geom_point(aes(col = Species), size = 4) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.1)+
  #facet_wrap(~Condition) +
  stat_cor(method = "spearman")+
  ylab("Increment of IS1 expression (log2FC)") +
  xlab("Number of IS1 in the genome") +
  scale_color_manual(values = mypal) +
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1,
        legend.position = "bottom",
        strip.background = element_blank());D
