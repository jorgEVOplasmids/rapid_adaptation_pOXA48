
setwd("~/Documents/TFM/Curvas")

library(tidyverse)
library(readxl)
library(lattice)
library(deSolve)
library(growthrates)
library(dplyr)
library(gridExtra)
library(caTools)
library(tidyr)
library(flux)
library(ggrepel)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(car)
#library(multcomp)

path_to_txt= "txt/"
path_to_Growthrates="GrowthRates_results/"
path_to_output="Output/"
Rosett_Alfonso <- read.delim("Rosett/Rosett_Alfonso")

file.list <- list.files(path =path_to_txt, full.names = F)
df.list <- lapply(paste0(path_to_txt, file.list), 
                  function(x)read.delim(x, header=T, nrows=133, dec=","))

attr(df.list, "names") <- file.list
df <- bind_rows(df.list, .id = "id") %>% mutate(Time=rep(seq(0, (133*10)-10, 10),35))


df_test<-df %>% 
  select( -`T..Optical.Density.600`) %>% gather(-Time,-id,  key = Well, value = OD ) %>% 
  separate(id, into=c("Plate", "Day"), remove = F) 

write.table(df_test, file=paste0(path_to_output, "curves_alf"))


curve_data <- df_test %>% 
  left_join(Rosett_Alfonso %>% mutate(Plate=as.character(Plate), 
                                      Day=as.character(Day),
                                      Well=as.character(Well))) %>% 
  mutate(Day=as.numeric(Day))

curve_data <- curve_data %>% 
  filter(
    Sample!="-") %>% 
  mutate(OD=as.numeric(OD))

# Keep out non-growing samples and contaminated

failed_K153 <- curve_data %>%
  filter(Sample == "K153" & Day == "1" & Replicate == "2" & Plasmid == "YES" & Antibiotic == "NO")

anti_join(curve_data, failed_K153) %>% 
  filter(Sample=="K153") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K153")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# Keep out non-growing samples and contaminated

failed_K163 <- curve_data %>%
  filter(Sample == "K163" & Replicate == "2" & Plasmid == "NO" & Antibiotic == "NO")

anti_join(curve_data, failed_K163) %>% 
  filter(Sample=="K163") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K163")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# Keep out non-growing samples and contaminated

failed_K25 <- rbind(curve_data %>%
                      filter(Sample == "K25" & Replicate == "4" & Plasmid == "NO" & Antibiotic == "NO"), curve_data %>%
                      filter(Sample == "K25" & Replicate == "6" & Plasmid == "YES" & Antibiotic == "NO"), curve_data %>%
                      filter(Sample == "K25" & Replicate == "5" & Plasmid == "YES" & Antibiotic == "YES"))

anti_join(curve_data, failed_K25) %>% 
  filter(Sample=="K25") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K25")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="C011") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C011")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# Keep out non-growing samples and contaminated

failed_C021 <- rbind(curve_data %>%
                       filter(Sample == "C021" & Replicate == "6" & Plasmid == "NO" & Antibiotic == "NO"), curve_data %>%
                       filter(Sample == "C021" & Replicate == "3" & Plasmid == "YES" & Antibiotic == "NO"), curve_data %>%
                       filter(Sample == "C021" & Replicate == "5" & Plasmid == "YES" & Antibiotic == "YES"))

anti_join(curve_data, failed_C021) %>% 
  filter(Sample=="C021") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C021")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# Keep out non-growing samples and contaminated

failed_C324 <- rbind(curve_data %>%
                       filter(Sample == "C324" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "4"), curve_data %>%
                       filter(Sample == "C324" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "6"), curve_data %>%
                       filter(Sample == "C324" & Plasmid == "YES" & Antibiotic == "NO" & Replicate == "3" & Day == "3"), curve_data %>%
                       filter(Sample == "C324" & Plasmid == "YES" & Antibiotic == "NO" & Replicate == "4" & Day == "3"),curve_data %>%
                       filter(Sample == "C324" & Plasmid == "YES" & Antibiotic == "NO" & Replicate == "5" & Day == "3"))

anti_join(curve_data, failed_C324) %>% 
  filter(Sample=="C324") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C324")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="C232") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C232")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# Keep out non-growing samples and contaminated

failed_C286 <- rbind(curve_data %>%
                       filter(Sample == "C286" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "2"), curve_data %>%
                       filter(Sample == "C286" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "5"), curve_data %>%
                       filter(Sample == "C286" & Plasmid == "YES" & Antibiotic == "NO" & Replicate == "2"), curve_data %>%
                       filter(Sample == "C286" & Plasmid == "YES" & Antibiotic == "NO" & Replicate == "4"), curve_data %>%
                       filter(Sample == "C286" & Plasmid == "YES" & Antibiotic == "NO" & Replicate == "5"), curve_data %>%
                       filter(Sample == "C286" & Plasmid == "YES" & Antibiotic == "YES" & Replicate == "5"), curve_data %>%
                       filter(Sample == "C286" & Plasmid == "YES" & Antibiotic == "NO" & Replicate == "6" & Day == "1"))

anti_join(curve_data, failed_C286) %>% 
  filter(Sample=="C286") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C286")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# Keep out non-growing samples and contaminated

failed_C309 <- rbind(curve_data %>%
                       filter(Sample == "C309" & Plasmid == "YES" & Antibiotic == "YES" & Replicate == "4"), curve_data %>%
                       filter(Sample == "C309" & Plasmid == "YES" & Antibiotic == "YES" & Replicate == "5"), curve_data %>%
                       filter(Sample == "C309" & Plasmid == "YES" & Antibiotic == "YES" & Replicate == "6"), curve_data %>%
                       filter(Sample == "C309" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "1" & Day == "1"), curve_data %>%
                       filter(Sample == "C309" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "4" & Day == "1"), curve_data %>%
                       filter(Sample == "C309" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "5" & Day == "1"))

anti_join(curve_data, failed_C309) %>% 
  filter(Sample=="C309") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C309")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# Keep out non-growing samples and contaminated

failed_K091 <- rbind(curve_data %>%
                       filter(Sample == "K091" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "1"), curve_data %>%
                       filter(Sample == "K091" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "4"), curve_data %>%
                       filter(Sample == "K091" & Plasmid == "YES" & Antibiotic == "NO" & Replicate == "1"), curve_data %>%
                       filter(Sample == "K091" & Plasmid == "YES" & Antibiotic == "YES" & Replicate == "3"))

anti_join(curve_data, failed_K091) %>% 
  filter(Sample=="K091") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K091")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="K112") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K112")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="K131") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K131")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# Keep out non-growing samples and contaminated

failed_K209 <- rbind(curve_data %>%
                       filter(Sample == "K209" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "1"), curve_data %>%
                       filter(Sample == "K209" & Plasmid == "YES" & Antibiotic == "NO" & Day == "3"))

anti_join(curve_data, failed_K209) %>% 
  filter(Sample=="K209") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K209")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# Keep out non-growing samples and contaminated

failed_H53 <- curve_data %>%
  filter(Sample == "H53" & Plasmid == "YES" & Antibiotic == "NO" & Replicate == "2")

anti_join(curve_data, failed_H53) %>% 
  filter(Sample=="H53") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="H53")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# Keep out non-growing samples and contaminated

failed_K147 <- rbind(curve_data %>%
                       filter(Sample == "K147" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "1"), curve_data %>%
                       filter(Sample == "K147" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "2"))

anti_join(curve_data, failed_K147) %>% 
  filter(Sample=="K147") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K147")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

failed_CF12 <- rbind(curve_data %>%
                       filter(Sample == "CF12" & Plasmid == "YES" & Antibiotic == "NO" & Replicate == "5"), curve_data %>%
                       filter(Sample == "CF12" & Plasmid == "YES" & Antibiotic == "YES" & Replicate == "6"))

anti_join(curve_data, failed_CF12) %>% 
  filter(Sample=="CF12") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="CF12")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

failed_CF13 <- rbind(curve_data %>%
                       filter(Sample == "CF13" & Plasmid == "NO" & Antibiotic == "NO" & Replicate == "4"), curve_data %>%
                       filter(Sample == "CF13" & Plasmid == "YES" & Antibiotic == "NO" & Replicate == "6"))


anti_join(curve_data, failed_CF13) %>% 
  filter(Sample=="CF13") %>%
  ggplot(aes(y=OD, x=Time/60, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line(size = 0.65)+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="CF13")+
  xlab("Time(h)")+
  ylab("OD600")+
  theme_bw(base_size = 16)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))
###########
curve_data %>% 
  filter(Project!="RIBO") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day~Sample)+
  labs(title="ALL")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Project=="RIBO") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  facet_grid(~Plasmid~Antibiotic~Day~Sample)+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))


if (file.exists(paste0(path_to_Growthrates, "GrowthRates_results"))){
  Growthrate_results<-read.table((paste0(path_to_Growthrates, "GrowthRates_results")), header=T) %>% filter(r2>0.95)
  
} else{
  manysplits<- all_easylinear(OD~Time | Plasmid + Sample +  Antibiotic + Replicate + Day + Project + Species,
                              data=anti_join(curve_data, curve_data %>% filter(Sample=="CF13" & Day==13 & Antibiotic=="NO" & Plasmid=="NO")))
  write.table(results(manysplits), paste0(path_to_Growthrates, "GrowthRates_results"))
  
  Growthrate_results<-results(manysplits) %>% filter(r2>0.95)
}

contaminations <- rbind(failed_C021, failed_C286, failed_C309, failed_C324,
                        failed_H53, failed_K091, failed_K147, failed_K153,
                        failed_K163, failed_K209, failed_K25)

curve_data_clean <- anti_join(curve_data, contaminations)

# Supplementary Fig 1

curve_data_clean %>%
  filter(Project!="RIBO", Sample != "C011", Sample != "C232", Sample != "Delta6", Sample != "MG1656", Sample != "K112", Sample != "K131") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  facet_grid(~Plasmid~Antibiotic~Day~Sample)+
  theme_bw()+
  xlab("Time (min)") +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

data_analysed<- curve_data_clean %>% group_by(Plasmid,Sample, Day,  Replicate, Antibiotic, Project, Species) %>% 
  group_modify(~ as.data.frame(flux::auc(.x$Time, .x$OD))) %>%
  mutate(AUC=`flux::auc(.x$Time, .x$OD)`) %>% 
  select(-`flux::auc(.x$Time, .x$OD)`)%>% 
  ungroup() 

data_analysed_2<- data_analysed %>% 
  left_join(curve_data_clean %>% 
              group_by(Plasmid, Sample, Day, Antibiotic, Project, Replicate, Species) %>% 
              summarise(ODmax=max(OD, na.rm = T)))

data<-data_analysed_2 %>% 
  mutate(Replicate_cntrl=Replicate, 
         Replicate=as.numeric(Replicate)-1) %>% 
  left_join(Growthrate_results)

data<-data %>% 
  filter(Day==1) %>% # Reference day
  group_by(Plasmid, Sample, Day, Replicate, Antibiotic, Project, Species) %>% 
  summarise(Vmax_day3=mumax, 
            lag_day3=lag, 
            AUC_day3=AUC, 
            ODmax_day3=ODmax) %>% 
  left_join(data, by=c("Plasmid", "Sample", "Replicate", "Antibiotic", "Project", "Species")) %>% 
  select(-Day.x) %>% 
  mutate(Day=Day.y) %>% 
  mutate(Treatment=ifelse(Plasmid=="YES" & Antibiotic=="YES", "Plasmid+Ab",
                          ifelse(Plasmid=="YES" & Antibiotic=="NO", "Plasmid", 
                                 ifelse(Plasmid=="NO" & Antibiotic=="NO", "Control", "problems"))))

data$Day <- as.factor(data$Day)

# Plot Fig 1

data %>%
  filter(ODmax>0.2, Day == 1 | Day == 15 | Day ==3, Sample != "C011", Sample != "C232", Sample != "Delta6", Sample != "MG1656", Sample != "K112", Sample != "K131") %>%
  ggplot(aes(y=AUC, x=Day, color=Treatment, group=interaction(Treatment))) +
  #geom_jitter(aes(y=AUC, x=Day, color=Treatment, fill = Treatment, group=interaction(Treatment,Day)), width = 0.1, size=1)+
  geom_boxplot(aes(y=AUC, x=Day, color=Treatment, fill = Treatment, group=interaction(Treatment,Day)), alpha = 0.2, size = 0.75)+
  geom_point(position=position_jitterdodge(jitter.width = 0.1), size=1)+
  #geom_line(size = 1.5)+
  facet_wrap(~factor(Species, levels = c("Ec", "Kp", "Cf")), scales = "free_x", nrow = 1)+
  #scale_y_continuous(limits=c(600,1600))+
  #scale_x_continuous("Day", labels = as.numeric(data$Day), breaks = data$Day)+
  scale_colour_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  ylab("Bacterial growth (AUC)")+
  xlab("Day")+
  theme_bw(base_size = 28)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

### Stats for growth days 1 to 15

data_for_stats <- data %>%
  filter(ODmax>0.2, Day == 1 | Day == 15, Sample != "C011", Sample != "C232", Sample != "Delta6", Sample != "MG1656", Sample != "K112", Sample != "K131")

model <- lm(AUC ~ Day + Plasmid + Antibiotic, data = data_for_stats)
summary(model)
ggqqplot(residuals(model))
shapiro.test(residuals(model))
# Compute Kruskal Wallis
growth.kt <- kruskal.test(AUC ~ interaction(Day, Plasmid, Antibiotic), data = data_for_stats); growth.kt
summary(growth.kt)

growth.kt.day <- wilcox.test(AUC ~ Day, data = data_for_stats); growth.kt.day
growth.kt.plasmid <- wilcox.test(AUC ~ Plasmid, data = data_for_stats); growth.kt.plasmid
growth.kt.ab <- wilcox.test(AUC ~ Antibiotic, data = data_for_stats); growth.kt.ab

# Per species

data_for_stats_Ec <- data %>%
  filter(ODmax>0.2, Species == "Ec", Day == 1 | Day == 15, Sample != "C011", Sample != "C232", Sample != "Delta6", Sample != "MG1656", Sample != "K112", Sample != "K131")

model <- lm(AUC ~ Day + Plasmid + Antibiotic, data = data_for_stats_Ec)
ggqqplot(residuals(model))
shapiro.test(residuals(model))
# Compute ANOVA
growth.aov <- aov(AUC ~ Day + Plasmid + Antibiotic, data = data_for_stats_Ec); growth.aov
summary(growth.aov)

data_for_stats_Kpn <- data %>%
  filter(ODmax>0.2, Species == "Kp", Day == 1 | Day == 15, Sample != "C011", Sample != "C232", Sample != "Delta6", Sample != "MG1656", Sample != "K112", Sample != "K131")

model <- lm(AUC ~ Day + Plasmid + Antibiotic, data = data_for_stats_Kpn)
ggqqplot(residuals(model))
shapiro.test(residuals(model))
# Check homoscedasticity
leveneTest(AUC ~ Day + Plasmid + Antibiotic, data = data_for_stats_Kpn)
# Compute Kruskal Wallis
growth.kt <- kruskal.test(AUC ~ interaction(Day, Plasmid, Antibiotic), data = data_for_stats_Kpn); growth.kt
summary(growth.kt)

growth.kt.day <- wilcox.test(AUC ~ Day, data = data_for_stats_Kpn); growth.kt.day
growth.kt.plasmid <- wilcox.test(AUC ~ Plasmid, data = data_for_stats_Kpn); growth.kt.plasmid
growth.kt.ab <- wilcox.test(AUC ~ Antibiotic, data = data_for_stats_Kpn); growth.kt.ab

data_for_stats_Cf <- data %>%
  filter(ODmax>0.2, Species == "Cf", Day == 1 | Day == 15, Sample != "C011", Sample != "C232", Sample != "Delta6", Sample != "MG1656", Sample != "K112", Sample != "K131")

model <- lm(AUC ~ Day + Plasmid + Antibiotic, data = data_for_stats_Cf)
ggqqplot(residuals(model))
shapiro.test(residuals(model))
# Check homoscedasticity
leveneTest(AUC ~ Day + Plasmid + Antibiotic, data = data_for_stats_Cf)
# Compute Kruskal Wallis
growth.aov <- aov(AUC ~ Day + Plasmid + Antibiotic, data = data_for_stats_Cf); growth.aov
summary(growth.aov)
