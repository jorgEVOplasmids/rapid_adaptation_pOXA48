
library(magicaxis)
library(rsalvador)
library(ggplot2)
library(scales)
library(tidyverse)
library(xlsx)
library(ggpubr)
library(car)
library(scales)

setwd("~/Downloads")

Mutants<-read.xlsx("traslucidas_merged.xlsx", sheetIndex = 1) 
Viables<-read.xlsx("viables_merged.xlsx", sheetIndex = 1) 

x=4

mu=seq(1,ncol(Mutants),1)
mut_rat=seq(1,ncol(Viables),1)
int95=matrix(NA,nrow=ncol(Mutants),ncol=2)
int_Mut_rat95=matrix(NA,nrow=ncol(Mutants),ncol=2)
int84=matrix(NA,nrow=ncol(Mutants),ncol=2)
int_Mut_rat84=matrix(NA,nrow=ncol(Mutants),ncol=2)
pvalue=matrix(NA, nrow=ncol(Mutants),ncol=ncol(Mutants))
test_statistic=matrix(NA, nrow=ncol(Mutants),ncol=ncol(Mutants))

Mutants<-Mutants/100000 # RSalvador troubleshooting
Viables<-Viables/100000

for (i in 1:ncol(Mutants)){
  
  mu[i]<-newton.LD(Mutants[,i][complete.cases(Mutants[,i])], show.iter = T)
  mut_rat[i]=mu[i]/(mean(Viables[,i][complete.cases(Viables[,i])]))
  int95[i,]<-confint.LD(Mutants[,i][complete.cases(Mutants[,i])])
  int_Mut_rat95[i,]<-int95[i,]/(mean(Viables[,i][complete.cases(Viables[,i])]))
  
  int84[i,]<-confint.LD(Mutants[,i][complete.cases(Mutants[,i])],alpha = .16)
  int_Mut_rat84[i,]<-int84[i,]/(mean(Viables[,i][complete.cases(Viables[,i])]))
  
}

pval<-matrix(ncol=8, nrow=77)
stat=seq(1,ncol(Mutants),1)
z=1
j=1

for(x in seq(1, length(Mutants), by=2)){
  for (h in 1:length(Mutants)) {
    
    print(paste0('WT=',colnames(Mutants[x])))
    print(paste0('h=',h))
    print(paste0('Z=',z))
    Mutants1<-Mutants[,x][complete.cases(Mutants[,x])]
    Viables1<-mean(Viables[,x][complete.cases(Viables[,x])])
    Mutants2<-Mutants[,h][complete.cases(Mutants[,h])]
    Viables2<-mean(Viables[,h][complete.cases(Viables[,h])])
    R=Viables2/Viables1
    
    lorr<- LRT.LD(Mutants1,Mutants2,R=R )
    
    pval[j,h]=lorr[2]
    stat[z]=lorr[1]
    z=z+1
  }
  j=j+1
}

colnames(pval)<-colnames(Mutants[1:ncol(Mutants)])
rownames(pval)<- colnames(Mutants[ seq(1, length(Mutants), by=2)])

Result<-data.frame(mu, mut_rat,int_Mut_rat95, int_Mut_rat84, colSums(!is.na(Mutants)),row.names = as.character(colnames(Mutants)))
colnames(Result)<-c("Mutations per Culture", "Mutation Rate","95% IC Lower Limit","95% IC Upper Limit","84% IC Lower Limit","84% IC Upper Limit",
                    "Number of replicates")
print(Result)

Result$Genotype<- rep(c("No pOXA-48", "pOXA-48", "pOXA-48 + pIS", "pIS"),2)
Result$Strain<-rownames(Result)
Result$Condition <- c(rep("LB+ARA", 4), rep("LB", 4))

### Import directly from xlsx with info

mutation_rates_K25 <- write.xlsx(Result,"mut_rate_K25_exp_merged.xlsx")

##########################    Import results from excel      ##############################################

Result <- read.xlsx("mut_rate_K25_exp_merged.xlsx", sheetIndex = 1)

#Result$log_pr <- log10(Result$Mutation.Rate)

# Plot all the controls

pIS_results <- Result %>%
  ggplot(aes(x = factor(NA., levels = c("K25", "K25p", "K25p_pIS", "K25_pIS", "K25_ARA", "K25p_ARA", "K25p_pIS_ARA", "K25_pIS_ARA")), y = `Mutation.Rate`, col = Genotype, fill = Genotype, shape = Condition)) +
  geom_point(size = 5) +
  theme_bw(base_size = 14) +
  #geom_bar(stat="identity", alpha = 0.2)+
  xlab("Genotype") +
  ylab("log(Phenotype rate)") +
  #facet_wrap(~Assay, scales = "free_x") +
  scale_y_continuous(trans='log10') +
  geom_vline(xintercept = 4.5, linewidth = 1, color = "darkgrey", linetype = "dashed") +
  #coord_flip()+
  geom_errorbar(aes(ymin=`Mutation.Rate`-X95..IC.Lower.Limit, ymax=`Mutation.Rate`+X95..IC.Upper.Limit), width=0,
                position=position_dodge(.9), col = "darkgrey") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1,
        strip.text.x = element_blank(),
        #legend.position = "none",
        strip.background = element_blank()); pIS_results

Result_mainfig <- Result %>%
  filter(Genotype != "pIS", NA. != "K25_ARA" , NA. != "K25p_ARA")

Result_mainfig$Assay <- c("pIS", "No pIS", "No pIS", "pIS")

pIS_results <- Result_mainfig %>%
  ggplot(aes(x = factor(NA., levels = c("K25", "K25p", "K25p_pIS", "K25p_pIS_ARA")), y = `Mutation.Rate`, col = Genotype, fill = Genotype, shape = Condition)) +
  geom_point(size = 5) +
  theme_bw(base_size = 14) +
  #geom_bar(stat="identity", alpha = 0.2)+
  xlab("Genotype") +
  ylab("log(Phenotype rate)") +
  facet_wrap(~Assay, scales = "free_x") +
  scale_y_continuous(trans='log10') +
  #coord_flip()+
  geom_errorbar(aes(ymin=`Mutation.Rate`-X95..IC.Lower.Limit, ymax=`Mutation.Rate`+X95..IC.Upper.Limit), width=0,
                position=position_dodge(.9), col = "darkgrey") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1,
        strip.text.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank()); pIS_results

# Barplot

starting_point <- 0.0000001
Result_mainfig$Mutation.Rate <- Result_mainfig$Mutation.Rate - starting_point

dfPlot <- Result_mainfig %>% 
  mutate(ymin = -log10(X95..IC.Lower.Limit),
         ymax = -log10(X95..IC.Upper.Limit),
         ymean = -log10(Mutation.Rate))

pIS_results <- Result_mainfig %>%
  ggplot(aes(x = factor(Strain, levels = c("K25", "K25p", "K25p_pIS", "K25p_pIS_ARA")), y = `Mutation.Rate`, col = Genotype, fill = Genotype, shape = Condition, ymin=10E-7, ymax=1)) +
  #geom_point(size = 5) +
  geom_rect(position=position_dodge(.8))+
  theme_bw(base_size = 14) +
  geom_bar(stat="identity", alpha = 0.2) +
  xlab("Genotype") +
  ylab("log(Phenotype rate)") +
  facet_wrap(~Assay, scales = "free_x") +
  scale_y_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))+
  geom_errorbar(aes(ymin=`Mutation.Rate`-X95..IC.Lower.Limit, ymax=`Mutation.Rate`+X95..IC.Upper.Limit), width=0,
                position=position_dodge(.9), col = "darkgrey") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1,
        strip.text.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank()); pIS_results
