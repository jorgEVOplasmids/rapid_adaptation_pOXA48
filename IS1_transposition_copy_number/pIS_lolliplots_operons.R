
### R script for representing the mutations on the klebsiella capsule operon and CF rpoS and nlpD operon

setwd("~/workie/figures_draft_expev/summary_parsed_vcall_w_plasmids/capsule_mutations")

library(trackViewer)
library(readxl)
library(dplyr)
library(GenomicRanges)

# Import data

muts_cap <- read_xlsx("capsule_operon_K25_pIS.xlsx", 1)
coords_cap <- read_xlsx("capsule_operon_K25_pIS.xlsx", 2)

# K25 vs K25p

muts_strain <- muts_cap %>%
  filter(Condition == "LB", Genotype == "K25" | Genotype == "K25p")

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(Type)

muts_condition <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(Genotype)

is_annotation <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(`IS mediated`)

second_seq <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(`Seq ID2`)

coords_cap$width <- coords_cap$end - coords_cap$begin

coords_cap_K25 <- coords_cap %>%
  filter(strain == "K25")

score <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(Freq)

mutations_gr <- GRanges("chr1", IRanges(position_muts, width = 1, names = muts_strain$names))
features <- GRanges("chr1", IRanges(as.vector(coords_cap_K25$begin), width = as.vector(coords_cap_K25$width)))
features$height <- c(0.05, 0.05 )

color <- replace(muts_annotation, muts_annotation == "Chromosome-Plasmid", "orange")
color <- replace(color, color == "non-synonymous", "black")
color <- replace(color, color == "synonymous", "white")
color <- replace(color, color == "intergenic", "grey")

second_seq<- replace(second_seq, second_seq == "IncL/M(pOXA-48)_1_pOXA-48", "dodgerblue3")
second_seq<- replace(second_seq, second_seq == "-", "black")
second_seq<- replace(second_seq, second_seq == "IncFII_1_pKP91", "darkmagenta")
side <- replace(muts_condition, muts_condition == "K25", "top")
side <- replace(side, side != "top", "bottom")
mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
#mutations_gr$score <- score
mutations_gr$label.parameter.rot <- 45

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(1710000, 1737000)))


# K25p_pIS LB vs K25p_pIS_ARA

muts_strain <- muts_cap %>%
  filter(Genotype == "K25p_pIS")

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(Type)

muts_condition <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(Condition)

is_annotation <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(`IS mediated`)

second_seq <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(`Seq ID2`)

coords_cap$width <- coords_cap$end - coords_cap$begin

coords_cap_K25 <- coords_cap %>%
  filter(strain == "K25")

score <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(Freq)

mutations_gr <- GRanges("chr1", IRanges(position_muts, width = 1, names = muts_strain$names))
features <- GRanges("chr1", IRanges(as.vector(coords_cap_K25$begin), width = as.vector(coords_cap_K25$width)))
features$height <- c(0.05, 0.05 )

color <- replace(muts_annotation, muts_annotation == "Chromosome-Plasmid", "orange")
color <- replace(color, color == "non-synonymous", "black")
color <- replace(color, color == "synonymous", "white")
color <- replace(color, color == "intergenic", "grey")

second_seq<- replace(second_seq, second_seq == "IncL/M(pOXA-48)_1_pOXA-48", "dodgerblue3")
second_seq<- replace(second_seq, second_seq == "-", "black")
second_seq<- replace(second_seq, second_seq == "IncFII_1_pKP91", "darkmagenta")
second_seq<- replace(second_seq, second_seq == "Chromosome-Chromosome", "black")
second_seq<- replace(second_seq, second_seq == "Chromosome", "black")
side <- replace(muts_condition, muts_condition == "LB+ARA", "top")
side <- replace(side, side != "top", "bottom")
color<- replace(color, color == "Chromosome-Chromosome", "black")
mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
#mutations_gr$score <- score
mutations_gr$label.parameter.rot <- 45

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(1710000, 1737000)))


# Merged repressed (K25 + K25p_pIS_ARA) vs induced (K25p + K25p_pIS)


muts_strain_expev <- muts_cap %>%
  filter(Condition == "LB", Genotype == "K25" | Genotype == "K25p")

muts_pIS <- muts_cap %>%
  filter(Genotype == "K25p_pIS")

muts_strain <- rbind(muts_strain_expev, muts_pIS)

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(Type)

muts_condition <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(Genotype)

is_annotation <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(`IS mediated`)

second_seq <- muts_strain %>%
  filter(Strain == "K25") %>%
  pull(`Seq ID2`)

coords_cap$width <- coords_cap$end - coords_cap$begin

coords_cap_K25 <- coords_cap %>%
  filter(strain == "K25")

score <- muts_strain %>%
  filter(Strain == "K25") %>%
  select(Genotype, Condition)

score <- as.data.frame(score)

merged_condition <- apply(score, 1, function(row) paste(row, collapse = ""))

merged_score <- replace(merged_condition, merged_condition == "K25p_pISLB+ARA", 30)
merged_score <- replace(merged_score, merged_score == "K25LB", 1)
merged_score <- replace(merged_score, merged_score == "K25p_pISLB", 30)
merged_score <- replace(merged_score, merged_score == "K25pLB", 1)
merged_score <- as.numeric(merged_score)


mutations_gr <- GRanges("chr1", IRanges(position_muts, width = 1, names = muts_strain$names))
features <- GRanges("chr1", IRanges(as.vector(coords_cap_K25$begin), width = as.vector(coords_cap_K25$width)))
features$height <- c(0.05, 0.05 )

color <- replace(muts_annotation, muts_annotation == "Chromosome-Plasmid", "orange")
color <- replace(color, color == "Chromosome-Chromosome", "black")
color <- replace(color, color == "non-synonymous", "black")
color <- replace(color, color == "synonymous", "white")
color <- replace(color, color == "intergenic", "grey")

second_seq<- replace(second_seq, second_seq == "IncL/M(pOXA-48)_1_pOXA-48", "dodgerblue3")
second_seq<- replace(second_seq, second_seq == "Chromosome", "black")
second_seq<- replace(second_seq, second_seq == "-", "black")
second_seq<- replace(second_seq, second_seq == "IncFII_1_pKP91", "darkmagenta")

merged_condition <- apply(score, 1, function(row) paste(row, collapse = ""))

side <- replace(merged_condition, merged_condition == "K25p_pISLB+ARA", "top")
side <- replace(side, side == "K25LB", "top")
side <- replace(side, side != "top", "bottom")

mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
mutations_gr$score <- merged_score
mutations_gr$label.parameter.rot <- 45

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(1710000, 1737000)))


