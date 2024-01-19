setwd("/phylogeny_and_phenotype/curves_DFE")

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
library(lubridate)
library(FSA)

### PLATE 1

# load data for OD as reported by the plate reader (only OD table)

od_data <- read.table("OD_table_plate_1_day_1.txt", header = FALSE)
od_data <- t(od_data)
colnames(od_data) <- od_data[1,]
od_data <- od_data[-1,]
#od_data
# Reshape table
# remove temperature column

od_data <- as.data.frame(subset(od_data, select = -c(2)))

# convert strings to numeric

od_data <- as.data.frame(lapply(od_data, function(x) {gsub(",", ".", x)}))
od_data$Time <- period_to_seconds(hms(od_data$Time))
od_data <- od_data %>% mutate_if(is.character, as.numeric)

# define dataframe containing 3 cols: Time (transformed to secs), OD, Well

od_data_reshaped_plate_1 <- pivot_longer(od_data, cols = 2:97, names_to = "Well", values_to = "OD")

# load sample information
# load sample id

sample_table <- as.data.frame(read_xlsx("metadata_plate_1.xlsx", sheet = 1))

sample_info <- pivot_longer(sample_table, cols = 2:13, names_to = "well_col", values_to = "Sample")
sample_info$Well <- gsub(" ", "", paste(sample_info$...1,sample_info$well_col))
sample_info <- as.data.frame(sample_info)
sample_info <- subset(sample_info, select = -c(...1, well_col))

# load condition

cond_table <- as.data.frame(read_xlsx("metadata_plate_1.xlsx", sheet = 2))

cond_info <- pivot_longer(cond_table, cols = 2:13, names_to = "well_col", values_to = "Condition")
cond_info$Well <- gsub(" ", "", paste(cond_info$...1,cond_info$well_col))
cond_info <- as.data.frame(cond_info)
cond_info <- subset(cond_info, select = -c(...1, well_col))

# load replicate

repl_table <- as.data.frame(read_xlsx("metadata_plate_1.xlsx", sheet = 3))

repl_info <- pivot_longer(repl_table, cols = 2:13, names_to = "well_col", values_to = "Replicate")
repl_info$Well <- gsub(" ", "", paste(repl_info$...1,repl_info$well_col))
repl_info <- as.data.frame(repl_info)
repl_info <- subset(repl_info, select = -c(...1, well_col))

# load strain

str_table <- as.data.frame(read_xlsx("metadata_plate_1.xlsx", sheet = 4))

str_info <- pivot_longer(str_table, cols = 2:13, names_to = "well_col", values_to = "Strain")
str_info$Well <- gsub(" ", "", paste(str_info$...1,str_info$well_col))
str_info <- as.data.frame(str_info)
str_info <- subset(str_info, select = -c(...1, well_col))

# include everything in the same dataframe and change OD format to numeric

od_data_reshaped_plate_1 <- inner_join(od_data_reshaped_plate_1, sample_info, by = "Well")
od_data_reshaped_plate_1 <- inner_join(od_data_reshaped_plate_1, repl_info, by = "Well")
od_data_reshaped_plate_1 <- inner_join(od_data_reshaped_plate_1, cond_info, by = "Well")
od_data_reshaped_plate_1 <- inner_join(od_data_reshaped_plate_1, str_info, by = "Well")
od_data_reshaped_plate_1$OD <- as.numeric(gsub(",", ".", od_data_reshaped_plate_1$OD))

### PLATE 2

# load data for OD as reported by the plate reader (only OD table)

od_data <- read.table("OD_table_plate_2_day_1.txt", header = FALSE)
od_data <- t(od_data)
colnames(od_data) <- od_data[1,]
od_data <- od_data[-1,]
#od_data
# Reshape table
# remove temperature column

od_data <- as.data.frame(subset(od_data, select = -c(2)))

# convert strings to numeric

od_data <- as.data.frame(lapply(od_data, function(x) {gsub(",", ".", x)}))
od_data$Time <- period_to_seconds(hms(od_data$Time))
od_data <- od_data %>% mutate_if(is.character, as.numeric)

# define dataframe containing 3 cols: Time (transformed to secs), OD, Well

od_data_reshaped_plate_2 <- pivot_longer(od_data, cols = 2:97, names_to = "Well", values_to = "OD")

# load sample information
# load sample id

sample_table <- as.data.frame(read_xlsx("metadata_plate_2.xlsx", sheet = 1))

sample_info <- pivot_longer(sample_table, cols = 2:13, names_to = "well_col", values_to = "Sample")
sample_info$Well <- gsub(" ", "", paste(sample_info$...1,sample_info$well_col))
sample_info <- as.data.frame(sample_info)
sample_info <- subset(sample_info, select = -c(...1, well_col))

# load condition

cond_table <- as.data.frame(read_xlsx("metadata_plate_2.xlsx", sheet = 2))

cond_info <- pivot_longer(cond_table, cols = 2:13, names_to = "well_col", values_to = "Condition")
cond_info$Well <- gsub(" ", "", paste(cond_info$...1,cond_info$well_col))
cond_info <- as.data.frame(cond_info)
cond_info <- subset(cond_info, select = -c(...1, well_col))

# load replicate

repl_table <- as.data.frame(read_xlsx("metadata_plate_2.xlsx", sheet = 3))

repl_info <- pivot_longer(repl_table, cols = 2:13, names_to = "well_col", values_to = "Replicate")
repl_info$Well <- gsub(" ", "", paste(repl_info$...1,repl_info$well_col))
repl_info <- as.data.frame(repl_info)
repl_info <- subset(repl_info, select = -c(...1, well_col))

# load strain

str_table <- as.data.frame(read_xlsx("metadata_plate_2.xlsx", sheet = 4))

str_info <- pivot_longer(str_table, cols = 2:13, names_to = "well_col", values_to = "Strain")
str_info$Well <- gsub(" ", "", paste(str_info$...1,str_info$well_col))
str_info <- as.data.frame(str_info)
str_info <- subset(str_info, select = -c(...1, well_col))

# include everything in the same dataframe and change OD format to numeric

od_data_reshaped_plate_2 <- inner_join(od_data_reshaped_plate_2, sample_info, by = "Well")
od_data_reshaped_plate_2 <- inner_join(od_data_reshaped_plate_2, repl_info, by = "Well")
od_data_reshaped_plate_2 <- inner_join(od_data_reshaped_plate_2, cond_info, by = "Well")
od_data_reshaped_plate_2 <- inner_join(od_data_reshaped_plate_2, str_info, by = "Well")
od_data_reshaped_plate_2$OD <- as.numeric(gsub(",", ".", od_data_reshaped_plate_2$OD))

### PLATE 3

# load data for OD as reported by the plate reader (only OD table)

od_data <- read.table("OD_table_plate_3_day_1.txt", header = FALSE)
od_data <- t(od_data)
colnames(od_data) <- od_data[1,]
od_data <- od_data[-1,]
#od_data
# Reshape table
# remove temperature column

od_data <- as.data.frame(subset(od_data, select = -c(2)))

# convert strings to numeric

od_data <- as.data.frame(lapply(od_data, function(x) {gsub(",", ".", x)}))
od_data$Time <- period_to_seconds(hms(od_data$Time))
od_data <- od_data %>% mutate_if(is.character, as.numeric)

# define dataframe containing 3 cols: Time (transformed to secs), OD, Well

od_data_reshaped_plate_3 <- pivot_longer(od_data, cols = 2:97, names_to = "Well", values_to = "OD")

# load sample information
# load sample id

sample_table <- as.data.frame(read_xlsx("metadata_plate_3.xlsx", sheet = 1))

sample_info <- pivot_longer(sample_table, cols = 2:13, names_to = "well_col", values_to = "Sample")
sample_info$Well <- gsub(" ", "", paste(sample_info$...1,sample_info$well_col))
sample_info <- as.data.frame(sample_info)
sample_info <- subset(sample_info, select = -c(...1, well_col))

# load condition

cond_table <- as.data.frame(read_xlsx("metadata_plate_3.xlsx", sheet = 2))

cond_info <- pivot_longer(cond_table, cols = 2:13, names_to = "well_col", values_to = "Condition")
cond_info$Well <- gsub(" ", "", paste(cond_info$...1,cond_info$well_col))
cond_info <- as.data.frame(cond_info)
cond_info <- subset(cond_info, select = -c(...1, well_col))

# load replicate

repl_table <- as.data.frame(read_xlsx("metadata_plate_3.xlsx", sheet = 3))

repl_info <- pivot_longer(repl_table, cols = 2:13, names_to = "well_col", values_to = "Replicate")
repl_info$Well <- gsub(" ", "", paste(repl_info$...1,repl_info$well_col))
repl_info <- as.data.frame(repl_info)
repl_info <- subset(repl_info, select = -c(...1, well_col))

# load strain

str_table <- as.data.frame(read_xlsx("metadata_plate_3.xlsx", sheet = 4))

str_info <- pivot_longer(str_table, cols = 2:13, names_to = "well_col", values_to = "Strain")
str_info$Well <- gsub(" ", "", paste(str_info$...1,str_info$well_col))
str_info <- as.data.frame(str_info)
str_info <- subset(str_info, select = -c(...1, well_col))

# include everything in the same dataframe and change OD format to numeric

od_data_reshaped_plate_3 <- inner_join(od_data_reshaped_plate_3, sample_info, by = "Well")
od_data_reshaped_plate_3 <- inner_join(od_data_reshaped_plate_3, repl_info, by = "Well")
od_data_reshaped_plate_3 <- inner_join(od_data_reshaped_plate_3, cond_info, by = "Well")
od_data_reshaped_plate_3 <- inner_join(od_data_reshaped_plate_3, str_info, by = "Well")
od_data_reshaped_plate_3$OD <- as.numeric(gsub(",", ".", od_data_reshaped_plate_3$OD))

### PLATE 5

# load data for OD as reported by the plate reader (only OD table)

od_data <- read.table("OD_table_plate_5_day_1.txt", header = FALSE)
od_data <- t(od_data)
colnames(od_data) <- od_data[1,]
od_data <- od_data[-1,]
#od_data
# Reshape table
# remove temperature column

od_data <- as.data.frame(subset(od_data, select = -c(2)))

# convert strings to numeric

od_data <- as.data.frame(lapply(od_data, function(x) {gsub(",", ".", x)}))
od_data$Time <- period_to_seconds(hms(od_data$Time))
od_data <- od_data %>% mutate_if(is.character, as.numeric)

# define dataframe containing 3 cols: Time (transformed to secs), OD, Well

od_data_reshaped_plate_5 <- pivot_longer(od_data, cols = 2:97, names_to = "Well", values_to = "OD")

# load sample information
# load sample id

sample_table <- as.data.frame(read_xlsx("metadata_plate_5.xlsx", sheet = 1))

sample_info <- pivot_longer(sample_table, cols = 2:13, names_to = "well_col", values_to = "Sample")
sample_info$Well <- gsub(" ", "", paste(sample_info$...1,sample_info$well_col))
sample_info <- as.data.frame(sample_info)
sample_info <- subset(sample_info, select = -c(...1, well_col))

# load condition

cond_table <- as.data.frame(read_xlsx("metadata_plate_5.xlsx", sheet = 2))

cond_info <- pivot_longer(cond_table, cols = 2:13, names_to = "well_col", values_to = "Condition")
cond_info$Well <- gsub(" ", "", paste(cond_info$...1,cond_info$well_col))
cond_info <- as.data.frame(cond_info)
cond_info <- subset(cond_info, select = -c(...1, well_col))

# load replicate

repl_table <- as.data.frame(read_xlsx("metadata_plate_5.xlsx", sheet = 3))

repl_info <- pivot_longer(repl_table, cols = 2:13, names_to = "well_col", values_to = "Replicate")
repl_info$Well <- gsub(" ", "", paste(repl_info$...1,repl_info$well_col))
repl_info <- as.data.frame(repl_info)
repl_info <- subset(repl_info, select = -c(...1, well_col))

# load strain

str_table <- as.data.frame(read_xlsx("metadata_plate_5.xlsx", sheet = 4))

str_info <- pivot_longer(str_table, cols = 2:13, names_to = "well_col", values_to = "Strain")
str_info$Well <- gsub(" ", "", paste(str_info$...1,str_info$well_col))
str_info <- as.data.frame(str_info)
str_info <- subset(str_info, select = -c(...1, well_col))

# include everything in the same dataframe and change OD format to numeric

od_data_reshaped_plate_5 <- inner_join(od_data_reshaped_plate_5, sample_info, by = "Well")
od_data_reshaped_plate_5 <- inner_join(od_data_reshaped_plate_5, repl_info, by = "Well")
od_data_reshaped_plate_5 <- inner_join(od_data_reshaped_plate_5, cond_info, by = "Well")
od_data_reshaped_plate_5 <- inner_join(od_data_reshaped_plate_5, str_info, by = "Well")
od_data_reshaped_plate_5$OD <- as.numeric(gsub(",", ".", od_data_reshaped_plate_5$OD))

### PLATE 6

# load data for OD as reported by the plate reader (only OD table)

od_data <- read.table("OD_table_plate_6_day_1.txt", header = FALSE)
od_data <- t(od_data)
colnames(od_data) <- od_data[1,]
od_data <- od_data[-1,]
#od_data
# Reshape table
# remove temperature column

od_data <- as.data.frame(subset(od_data, select = -c(2)))

# convert strings to numeric

od_data <- as.data.frame(lapply(od_data, function(x) {gsub(",", ".", x)}))
od_data$Time <- period_to_seconds(hms(od_data$Time))
od_data <- od_data %>% mutate_if(is.character, as.numeric)

# define dataframe containing 3 cols: Time (transformed to secs), OD, Well

od_data_reshaped_plate_6 <- pivot_longer(od_data, cols = 2:97, names_to = "Well", values_to = "OD")

# load sample information
# load sample id

sample_table <- as.data.frame(read_xlsx("metadata_plate_6.xlsx", sheet = 1))

sample_info <- pivot_longer(sample_table, cols = 2:13, names_to = "well_col", values_to = "Sample")
sample_info$Well <- gsub(" ", "", paste(sample_info$...1,sample_info$well_col))
sample_info <- as.data.frame(sample_info)
sample_info <- subset(sample_info, select = -c(...1, well_col))

# load condition

cond_table <- as.data.frame(read_xlsx("metadata_plate_6.xlsx", sheet = 2))

cond_info <- pivot_longer(cond_table, cols = 2:13, names_to = "well_col", values_to = "Condition")
cond_info$Well <- gsub(" ", "", paste(cond_info$...1,cond_info$well_col))
cond_info <- as.data.frame(cond_info)
cond_info <- subset(cond_info, select = -c(...1, well_col))

# load replicate

repl_table <- as.data.frame(read_xlsx("metadata_plate_6.xlsx", sheet = 3))

repl_info <- pivot_longer(repl_table, cols = 2:13, names_to = "well_col", values_to = "Replicate")
repl_info$Well <- gsub(" ", "", paste(repl_info$...1,repl_info$well_col))
repl_info <- as.data.frame(repl_info)
repl_info <- subset(repl_info, select = -c(...1, well_col))

# load strain

str_table <- as.data.frame(read_xlsx("metadata_plate_6.xlsx", sheet = 4))

str_info <- pivot_longer(str_table, cols = 2:13, names_to = "well_col", values_to = "Strain")
str_info$Well <- gsub(" ", "", paste(str_info$...1,str_info$well_col))
str_info <- as.data.frame(str_info)
str_info <- subset(str_info, select = -c(...1, well_col))

# include everything in the same dataframe and change OD format to numeric

od_data_reshaped_plate_6 <- inner_join(od_data_reshaped_plate_6, sample_info, by = "Well")
od_data_reshaped_plate_6 <- inner_join(od_data_reshaped_plate_6, repl_info, by = "Well")
od_data_reshaped_plate_6 <- inner_join(od_data_reshaped_plate_6, cond_info, by = "Well")
od_data_reshaped_plate_6 <- inner_join(od_data_reshaped_plate_6, str_info, by = "Well")
od_data_reshaped_plate_6$OD <- as.numeric(gsub(",", ".", od_data_reshaped_plate_6$OD))

### STATS MATRICES OF EACH PLATE

### PLATE 1

stats_curves_plate_1 <- read.table("matrix_stats_plate_1_day_1.txt", header = TRUE)
# Import AUC
stats_curves_plate_1 <- as.data.frame(lapply(stats_curves_plate_1, function(x) {gsub(",", ".", x)}))
AUC <- as.data.frame(cbind(stats_curves_plate_1$Well, stats_curves_plate_1$Integral.Integral.OpticalDensity.600.))
colnames(AUC) <- c("Well", "AUC")

od_data_reshaped_plate_1 <- inner_join(od_data_reshaped_plate_1, AUC, by = "Well")
od_data_reshaped_plate_1$AUC <- as.numeric(od_data_reshaped_plate_1$AUC)

# Import maxOD
maxOD <- as.data.frame(cbind(stats_curves_plate_1$Well, stats_curves_plate_1$maxOD.MeanMaxOD.OpticalDensity.600.))
colnames(maxOD) <- c("Well", "maxOD")

od_data_reshaped_plate_1 <- inner_join(od_data_reshaped_plate_1, maxOD, by = "Well")
od_data_reshaped_plate_1$maxOD <- as.numeric(od_data_reshaped_plate_1$maxOD)

# Import mumax
mumax <- as.data.frame(cbind(stats_curves_plate_1$Well, stats_curves_plate_1$CustomStatistics.MaxV.OpticalDensity.600.))
colnames(mumax) <- c("Well", "mumax")

od_data_reshaped_plate_1 <- inner_join(od_data_reshaped_plate_1, mumax, by = "Well")
od_data_reshaped_plate_1$mumax <- as.numeric(od_data_reshaped_plate_1$mumax)

# Import lag time and reformat to integers

lag <- as.data.frame(cbind(stats_curves_plate_1$Well, stats_curves_plate_1$CustomStatistics.Lagtime.OpticalDensity.600.))
colnames(lag) <- c("Well", "lag")

lag$lag <- period_to_seconds(hms(lag$lag))

od_data_reshaped_plate_1 <- inner_join(od_data_reshaped_plate_1, lag, by = "Well")
od_data_reshaped_plate_1$lag <- as.numeric(od_data_reshaped_plate_1$lag)

### PLATE 2

stats_curves_plate_2 <- read.table("matrix_stats_plate_2_day_1.txt", header = TRUE)
# Import AUC
stats_curves_plate_2 <- as.data.frame(lapply(stats_curves_plate_2, function(x) {gsub(",", ".", x)}))
AUC <- as.data.frame(cbind(stats_curves_plate_2$Well, stats_curves_plate_2$Integral.Integral.OpticalDensity.600.))
colnames(AUC) <- c("Well", "AUC")

od_data_reshaped_plate_2 <- inner_join(od_data_reshaped_plate_2, AUC, by = "Well")
od_data_reshaped_plate_2$AUC <- as.numeric(od_data_reshaped_plate_2$AUC)

# Import maxOD
maxOD <- as.data.frame(cbind(stats_curves_plate_2$Well, stats_curves_plate_2$maxOD.MeanMaxOD.OpticalDensity.600.))
colnames(maxOD) <- c("Well", "maxOD")

od_data_reshaped_plate_2 <- inner_join(od_data_reshaped_plate_2, maxOD, by = "Well")
od_data_reshaped_plate_2$maxOD <- as.numeric(od_data_reshaped_plate_2$maxOD)

# Import mumax
mumax <- as.data.frame(cbind(stats_curves_plate_2$Well, stats_curves_plate_2$CustomStatistics.MaxV.OpticalDensity.600.))
colnames(mumax) <- c("Well", "mumax")

od_data_reshaped_plate_2 <- inner_join(od_data_reshaped_plate_2, mumax, by = "Well")
od_data_reshaped_plate_2$mumax <- as.numeric(od_data_reshaped_plate_2$mumax)

# Import lag time and reformat to integers

lag <- as.data.frame(cbind(stats_curves_plate_2$Well, stats_curves_plate_2$CustomStatistics.Lagtime.OpticalDensity.600.))
colnames(lag) <- c("Well", "lag")

lag$lag <- period_to_seconds(hms(lag$lag))

od_data_reshaped_plate_2 <- inner_join(od_data_reshaped_plate_2, lag, by = "Well")
od_data_reshaped_plate_2$lag <- as.numeric(od_data_reshaped_plate_2$lag)

### PLATE 3

stats_curves_plate_3 <- read.table("matrix_stats_plate_3_day_1.txt", header = TRUE)
# Import AUC
stats_curves_plate_3 <- as.data.frame(lapply(stats_curves_plate_3, function(x) {gsub(",", ".", x)}))
AUC <- as.data.frame(cbind(stats_curves_plate_3$Well, stats_curves_plate_3$Integral.Integral.OpticalDensity.600.))
colnames(AUC) <- c("Well", "AUC")

od_data_reshaped_plate_3 <- inner_join(od_data_reshaped_plate_3, AUC, by = "Well")
od_data_reshaped_plate_3$AUC <- as.numeric(od_data_reshaped_plate_3$AUC)

# Import maxOD
maxOD <- as.data.frame(cbind(stats_curves_plate_3$Well, stats_curves_plate_3$maxOD.MeanMaxOD.OpticalDensity.600.))
colnames(maxOD) <- c("Well", "maxOD")

od_data_reshaped_plate_3 <- inner_join(od_data_reshaped_plate_3, maxOD, by = "Well")
od_data_reshaped_plate_3$maxOD <- as.numeric(od_data_reshaped_plate_3$maxOD)

# Import mumax
mumax <- as.data.frame(cbind(stats_curves_plate_3$Well, stats_curves_plate_3$CustomStatistics.MaxV.OpticalDensity.600.))
colnames(mumax) <- c("Well", "mumax")

od_data_reshaped_plate_3 <- inner_join(od_data_reshaped_plate_3, mumax, by = "Well")
od_data_reshaped_plate_3$mumax <- as.numeric(od_data_reshaped_plate_3$mumax)

# Import lag time and reformat to integers

lag <- as.data.frame(cbind(stats_curves_plate_3$Well, stats_curves_plate_3$CustomStatistics.Lagtime.OpticalDensity.600.))
colnames(lag) <- c("Well", "lag")

lag$lag <- period_to_seconds(hms(lag$lag))

od_data_reshaped_plate_3 <- inner_join(od_data_reshaped_plate_3, lag, by = "Well")
od_data_reshaped_plate_3$lag <- as.numeric(od_data_reshaped_plate_3$lag)

### PLATE 5

stats_curves_plate_5 <- read.table("matrix_stats_plate_5_day_1.txt", header = TRUE)
# Import AUC
stats_curves_plate_5 <- as.data.frame(lapply(stats_curves_plate_5, function(x) {gsub(",", ".", x)}))
AUC <- as.data.frame(cbind(stats_curves_plate_5$Well, stats_curves_plate_5$Integral.Integral.OpticalDensity.600.))
colnames(AUC) <- c("Well", "AUC")

od_data_reshaped_plate_5 <- inner_join(od_data_reshaped_plate_5, AUC, by = "Well")
od_data_reshaped_plate_5$AUC <- as.numeric(od_data_reshaped_plate_5$AUC)

# Import maxOD
maxOD <- as.data.frame(cbind(stats_curves_plate_5$Well, stats_curves_plate_5$maxOD.MeanMaxOD.OpticalDensity.600.))
colnames(maxOD) <- c("Well", "maxOD")

od_data_reshaped_plate_5 <- inner_join(od_data_reshaped_plate_5, maxOD, by = "Well")
od_data_reshaped_plate_5$maxOD <- as.numeric(od_data_reshaped_plate_5$maxOD)

# Import mumax
mumax <- as.data.frame(cbind(stats_curves_plate_5$Well, stats_curves_plate_5$CustomStatistics.MaxV.OpticalDensity.600.))
colnames(mumax) <- c("Well", "mumax")

od_data_reshaped_plate_5 <- inner_join(od_data_reshaped_plate_5, mumax, by = "Well")
od_data_reshaped_plate_5$mumax <- as.numeric(od_data_reshaped_plate_5$mumax)

# Import lag time and reformat to integers

lag <- as.data.frame(cbind(stats_curves_plate_5$Well, stats_curves_plate_5$CustomStatistics.Lagtime.OpticalDensity.600.))
colnames(lag) <- c("Well", "lag")

lag$lag <- period_to_seconds(hms(lag$lag))

od_data_reshaped_plate_5 <- inner_join(od_data_reshaped_plate_5, lag, by = "Well")
od_data_reshaped_plate_5$lag <- as.numeric(od_data_reshaped_plate_5$lag)

### PLATE 5

stats_curves_plate_6 <- read.table("matrix_stats_plate_6_day_1.txt", header = TRUE)
# Import AUC
stats_curves_plate_6 <- as.data.frame(lapply(stats_curves_plate_6, function(x) {gsub(",", ".", x)}))
AUC <- as.data.frame(cbind(stats_curves_plate_6$Well, stats_curves_plate_6$Integral.Integral.OpticalDensity.600.))
colnames(AUC) <- c("Well", "AUC")

od_data_reshaped_plate_6 <- inner_join(od_data_reshaped_plate_6, AUC, by = "Well")
od_data_reshaped_plate_6$AUC <- as.numeric(od_data_reshaped_plate_6$AUC)

# Import maxOD
maxOD <- as.data.frame(cbind(stats_curves_plate_6$Well, stats_curves_plate_6$maxOD.MeanMaxOD.OpticalDensity.600.))
colnames(maxOD) <- c("Well", "maxOD")

od_data_reshaped_plate_6 <- inner_join(od_data_reshaped_plate_6, maxOD, by = "Well")
od_data_reshaped_plate_6$maxOD <- as.numeric(od_data_reshaped_plate_6$maxOD)

# Import mumax
mumax <- as.data.frame(cbind(stats_curves_plate_6$Well, stats_curves_plate_6$CustomStatistics.MaxV.OpticalDensity.600.))
colnames(mumax) <- c("Well", "mumax")

od_data_reshaped_plate_6 <- inner_join(od_data_reshaped_plate_6, mumax, by = "Well")
od_data_reshaped_plate_6$mumax <- as.numeric(od_data_reshaped_plate_6$mumax)

# Import lag time and reformat to integers

lag <- as.data.frame(cbind(stats_curves_plate_6$Well, stats_curves_plate_6$CustomStatistics.Lagtime.OpticalDensity.600.))
colnames(lag) <- c("Well", "lag")

lag$lag <- period_to_seconds(hms(lag$lag))

od_data_reshaped_plate_6 <- inner_join(od_data_reshaped_plate_6, lag, by = "Well")
od_data_reshaped_plate_6$lag <- as.numeric(od_data_reshaped_plate_6$lag)

# Merge everything in the same table

od_data_day_1 <- rbind(od_data_reshaped_plate_1, od_data_reshaped_plate_2,
                       od_data_reshaped_plate_3, od_data_reshaped_plate_5,
                       od_data_reshaped_plate_6)

od_data_day_1 <- as.data.frame(od_data_day_1)

od_data_day_1$Time <- as.numeric(od_data_day_1$Time)
od_data_day_1$OD <- as.numeric(od_data_day_1$OD)
od_data_day_1$Replicate <- as.numeric(od_data_day_1$Replicate)


od_data_day_1 %>%
  filter(Time == 86884 | Time == 72560) %>% # get only the last datapoint for each replicate
  ggplot(aes(y=AUC/60, x=Strain, fill = Replicate, group=Replicate)) +
  geom_boxplot(aes(y=AUC/60, x=Strain, fill = Replicate, group=Replicate, outlier.shape = NULL, alpha = 0.2))+
  geom_point(aes(y=AUC/60, x=Strain, fill = Replicate, group=Replicate), position = position_jitter(w = 0.1, h = 0), size = 3)+
  ylab("AUC")+
  xlab("Sample") +
  scale_x_discrete(limits = c("d6", "P1", "P2", "P3", "P4", "P5", "MG1655"))+
  theme_bw(base_size = 28) +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

# Filter useless data and non used strains

od_data_day_1 <- od_data_day_1 %>%
  filter(Sample != "-", Replicate != "-", Condition != "-",
         Strain != "-", Strain != "d6",
         Strain != "C011", Strain != "C232",
         Strain != "K131", Strain != "K112") 

od_data_day_1 %>%
  filter(Sample != "-", Replicate != "-") +
  ggplot(aes(y=OD, x=Time, fill = Replicate, group=Replicate)) +
  geom_line() +
  facet_wrap(~Strain)

od_data_day_1 %>%
  filter(Sample != "-", Replicate != "-") +
  ggplot(aes(y = AUC, x = Strain, group = Replicate)) +
  geom_boxplot()

od_data_day_1_def <- od_data_day_1 %>% filter(Time == 86884 | Time == 72560)
od_data_day_1_def$AUC <- od_data_day_1_def$AUC/60

#write.xlsx(od_data_day_1_def, "od_data_day_1_def.xlsx")

costs_data <- read_xlsx("od_data_day_1_def.xlsx", sheet = 2)
costs_data <- as.data.frame(costs_data)

costs_data %>%
  ggplot(aes(y = Median_cost-1, x = factor(Strain, levels = c("K209", "C286", "CF13", "C324", "K153", "C309", "K147", "CF12", "K163", "K25", "H53", "C021", "K091")), col = Species, fill = Species)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x=Strain, ymin=Median_cost-1-Stderr_cost, ymax=Median_cost-1+Stderr_cost), width=0, alpha=0.4, size=1, col = "black")+
  ylab("Relative fitness (w)") +
  xlab("Strain") +
  theme_bw(base_size = 24)+
  ylim(-0.3,0.3)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

### Load data for stats

costs_data_stats <- read_xlsx("od_data_day_1_def.xlsx", sheet = 1)
costs_data_stats <- as.data.frame(costs_data_stats)

### Check normality

model <- lm(AUC ~ Condition + Strain, data = costs_data_stats)

ggqqplot(residuals(model))

shapiro.test(residuals(model))

leveneTest(AUC ~ Condition, data = costs_data_stats)

stat.test <- costs_data_stats %>%
  group_by(Condition) %>%
  pairwise_t_test(
    Median_AUC ~ Strain, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) %>%
  select(-df, -statistic, -p) # Remove details
stat.test

table_pvalues <- data.frame(Strain=character(0), pvalue=numeric(0))

table_pvalues <- table_pvalues %>% add_row (Strain= "C021",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "C021" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "C021" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "C286",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "C286" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "C286" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "C309",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "C309" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "C309" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "C324",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "C324" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "C324" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "CF12",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "CF12" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "CF12" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "CF13",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "CF13" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "CF13" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "H53",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "H53" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "H53" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "K091",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "K091" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "K091" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "K147",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "K147" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "K147" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "K153",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "K153" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "K153" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "K163",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "K163" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "K163" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "K209",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "K209" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "K209" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "K25",
                                            pvalue= t.test(costs_data_stats[costs_data_stats$Strain == "K25" & costs_data_stats$Condition == "pOXA-48", ]$AUC,
                                                           costs_data_stats[costs_data_stats$Strain == "K25" & costs_data_stats$Condition == "No pOXA-48", ]$AUC,
                                                           paired=TRUE)$p.value)


table_pvalues$p.adj.fdr <- p.adjust(table_pvalues$pvalue, method="bonferroni", n=13)


#### Plot histogram of cost distribution to compare with Aida's paper

costs_data %>%
  ggplot(aes(x = Median_cost)) +
  geom_histogram() +
  ylab("Number of strains") +
  xlab("Relative fitness") +
  theme_bw(base_size = 24)+
  ylim(-0.3,0.3)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))


# Plot barplots with points

costs_data <- read_xlsx("od_data_day_1_def.xlsx", sheet = 1)
df.summary <- read_xlsx("od_data_day_1_def.xlsx", sheet = 2)
costs_data <- as.data.frame(costs_data)
df.summary <- as.data.frame(df.summary)

costs_data <- costs_data %>%
  filter(Condition == "pOXA-48")

mypal <- c("#b8d7e7ff", "#daa3b8ff", "#a58ac8ff")

costs_data%>%
  ggplot(aes(y = Relative_AUC-1, x = factor(Strain, levels = c("K209", "C286", "CF13", "C309", "K153", "CF12", "K147", "C324", "K25", "H53", "K163", "C021", "K091")), col = Species, fill = Species)) +
  geom_bar(stat = "identity", data = df.summary, size = 0) +
  geom_jitter(size = 3, width = 0.1, shape = 1, stroke = 1) +
  geom_errorbar(aes(x=Strain, ymin=Relative_AUC-1-Stderr_rel_AUC, ymax=Relative_AUC-1+Stderr_rel_AUC), width=0, alpha=0.6, size=1, col = "darkgrey", data = df.summary)+
  ylab("Relative AUC") +
  scale_color_manual(values = mypal) +
  scale_fill_manual(values = mypal) +
  xlab("Strain") +
  theme_bw(base_size = 28)+
  ylim(-0.6,0.2)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))
