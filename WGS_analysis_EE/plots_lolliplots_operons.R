
### R script for representing the mutations on the klebsiella capsule operon and citrobacter rpoS and nlpD operon

# Replace working directory by directory with xlsx files of mutations placed in the operons

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



# K163

muts_strain <- muts_cap %>%
  filter(Strain == "K163")

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_cap %>%
  filter(Strain == "K163") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_cap %>%
  filter(Strain == "K163") %>%
  pull(Type)

muts_condition <- muts_cap %>%
  filter(Strain == "K163") %>%
  pull(Condition)

is_annotation <- muts_cap %>%
  filter(Strain == "K163") %>%
  pull(`IS mediated`)

second_seq <- muts_cap %>%
  filter(Strain == "K163") %>%
  pull(`Seq ID2`)

coords_cap$width <- coords_cap$end - coords_cap$begin

coords_cap_K163 <- coords_cap %>%
  filter(strain == "K163")

score <- muts_cap %>%
  filter(Strain == "K163") %>%
  pull(Freq)

mutations_gr <- GRanges("chr1", IRanges(position_muts, width = 1, names = muts_strain$names))
features <- GRanges("chr1", IRanges(as.vector(coords_cap_K163$begin), width = as.vector(coords_cap_K163$width)))
features$height <- c(0.05, 0.05 )

color <- replace(muts_annotation, muts_annotation == "Chromosome-Plasmid", "orange")
color <- replace(color, color == "non-synonymous", "black")
color <- replace(color, color == "synonymous", "white")
color <- replace(color, color == "intergenic", "grey")

second_seq<- replace(second_seq, second_seq == "IncL/M(pOXA-48)_1_pOXA-48", "dodgerblue3")
second_seq<- replace(second_seq, second_seq == "-", "black")
second_seq<- replace(second_seq, second_seq == "IncFIB(K)_1_Kpn3/IncFII_1_pKP91", "darkmagenta")
side <- replace(muts_condition, muts_condition == "No pOXA-48", "top")
side <- replace(side, side != "top", "bottom")
mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
mutations_gr$score <- score
mutations_gr$label.parameter.rot <- 45

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(1761000, 1784000)), legend = "legend")

# K25

muts_strain <- muts_cap %>%
  filter(Strain == "K25")

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_cap %>%
  filter(Strain == "K25") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_cap %>%
  filter(Strain == "K25") %>%
  pull(Type)

muts_condition <- muts_cap %>%
  filter(Strain == "K25") %>%
  pull(Condition)

is_annotation <- muts_cap %>%
  filter(Strain == "K25") %>%
  pull(`IS mediated`)

second_seq <- muts_cap %>%
  filter(Strain == "K25") %>%
  pull(`Seq ID2`)

coords_cap$width <- coords_cap$end - coords_cap$begin

coords_cap_K25 <- coords_cap %>%
  filter(strain == "K25")

score <- muts_cap %>%
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
second_seq<- replace(second_seq, second_seq == "IncFIB(K)_1_Kpn3/IncFII_1_pKP91", "darkmagenta")
side <- replace(muts_condition, muts_condition == "No pOXA-48", "top")
side <- replace(side, side != "top", "bottom")
mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
mutations_gr$score <- score
mutations_gr$label.parameter.rot <- 45

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(1710000, 1737000)))

# H53

muts_strain <- muts_cap %>%
  filter(Strain == "H53")

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_cap %>%
  filter(Strain == "H53") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_cap %>%
  filter(Strain == "H53") %>%
  pull(Type)

muts_condition <- muts_cap %>%
  filter(Strain == "H53") %>%
  pull(Condition)

is_annotation <- muts_cap %>%
  filter(Strain == "H53") %>%
  pull(`IS mediated`)

second_seq <- muts_cap %>%
  filter(Strain == "H53") %>%
  pull(`Seq ID2`)

coords_cap$width <- coords_cap$end - coords_cap$begin

coords_cap_H53 <- coords_cap %>%
  filter(strain == "H53")

score <- muts_cap %>%
  filter(Strain == "H53") %>%
  pull(Freq)

mutations_gr <- GRanges("chr1", IRanges(position_muts, width = 1, names = muts_strain$names))
features <- GRanges("chr1", IRanges(as.vector(coords_cap_H53$begin), width = as.vector(coords_cap_H53$width)))
features$height <- c(0.05, 0.05 )

color <- replace(muts_annotation, muts_annotation == "Chromosome-Plasmid", "orange")
color <- replace(color, color == "non-synonymous", "black")
color <- replace(color, color == "Chromosome-Chromosome", "orange")
color <- replace(color, color == "synonymous", "white")
color <- replace(color, color == "intergenic", "grey")

second_seq<- replace(second_seq, second_seq == "IncL/M(pOXA-48)_1_pOXA-48", "dodgerblue3")
second_seq<- replace(second_seq, second_seq == "-", "black")
second_seq<- replace(second_seq, second_seq == "Chromosome", "black")
second_seq<- replace(second_seq, second_seq == "IncFIB(K)_1_Kpn3", "darkmagenta")
side <- replace(muts_condition, muts_condition == "No pOXA-48", "top")
side <- replace(side, side != "top", "bottom")
mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
mutations_gr$score <- score
mutations_gr$label.parameter.rot <- 45

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(1698000, 1721000)))

# K091

muts_strain <- muts_cap %>%
  filter(Strain == "K091")

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_cap %>%
  filter(Strain == "K091") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_cap %>%
  filter(Strain == "K091") %>%
  pull(Type)

muts_condition <- muts_cap %>%
  filter(Strain == "K091") %>%
  pull(Condition)

is_annotation <- muts_cap %>%
  filter(Strain == "K091") %>%
  pull(`IS mediated`)

second_seq <- muts_cap %>%
  filter(Strain == "K091") %>%
  pull(`Seq ID2`)

coords_cap$width <- coords_cap$end - coords_cap$begin

coords_cap_K091 <- coords_cap %>%
  filter(strain == "K091")

score <- muts_cap %>%
  filter(Strain == "K091") %>%
  pull(Freq)

mutations_gr <- GRanges("chr1", IRanges(position_muts, width = 1, names = muts_strain$names))
features <- GRanges("chr1", IRanges(as.vector(coords_cap_K091$begin), width = as.vector(coords_cap_K091$width)))
features$height <- c(0.05, 0.05 )

color <- replace(muts_annotation, muts_annotation == "Chromosome-Plasmid", "orange")
color <- replace(color, color == "non-synonymous", "black")
color <- replace(color, color == "Chromosome-Chromosome", "orange")
color <- replace(color, color == "synonymous", "white")
color <- replace(color, color == "intergenic", "grey")

second_seq<- replace(second_seq, second_seq == "IncL/M(pOXA-48)_1_pOXA-48", "dodgerblue3")
second_seq<- replace(second_seq, second_seq == "-", "black")
second_seq<- replace(second_seq, second_seq == "Chromosome", "black")
second_seq<- replace(second_seq, second_seq == "IncFII_1_pKP91/IncFIA(HI1)_1_HI1", "darkmagenta")
side <- replace(muts_condition, muts_condition == "No pOXA-48", "top")
side <- replace(side, side != "top", "bottom")
mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
mutations_gr$score <- score
mutations_gr$label.parameter.rot <- 45

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(1704000, 1717000)), legend = "legend")

# K147

muts_strain <- muts_cap %>%
  filter(Strain == "K147")

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_cap %>%
  filter(Strain == "K147") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_cap %>%
  filter(Strain == "K147") %>%
  pull(Type)

muts_condition <- muts_cap %>%
  filter(Strain == "K147") %>%
  pull(Condition)

is_annotation <- muts_cap %>%
  filter(Strain == "K147") %>%
  pull(`IS mediated`)

second_seq <- muts_cap %>%
  filter(Strain == "K147") %>%
  pull(`Seq ID2`)

coords_cap$width <- coords_cap$end - coords_cap$begin

coords_cap_K147 <- coords_cap %>%
  filter(strain == "K147")

score <- muts_cap %>%
  filter(Strain == "K147") %>%
  pull(Freq)

mutations_gr <- GRanges("chr1", IRanges(position_muts, width = 1, names = muts_strain$names))
features <- GRanges("chr1", IRanges(as.vector(coords_cap_K147$begin), width = as.vector(coords_cap_K147$width)))
features$height <- c(0.05, 0.05 )

color <- replace(muts_annotation, muts_annotation == "Chromosome-Plasmid", "orange")
color <- replace(color, color == "non-synonymous", "black")
color <- replace(color, color == "Chromosome-Chromosome", "orange")
color <- replace(color, color == "synonymous", "white")
color <- replace(color, color == "intergenic", "grey")

second_seq<- replace(second_seq, second_seq == "IncL/M(pOXA-48)_1_pOXA-48", "dodgerblue3")
second_seq<- replace(second_seq, second_seq == "-", "black")
second_seq<- replace(second_seq, second_seq == "Chromosome", "black")
second_seq<- replace(second_seq, second_seq == "IncFIB(K)_1_Kpn3/IncFII_1_pKP91", "darkmagenta")
side <- replace(muts_condition, muts_condition == "No pOXA-48", "top")
side <- replace(side, side != "top", "bottom")
mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
mutations_gr$score <- score
mutations_gr$label.parameter.rot <- 45

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(1704000, 1729000)))

# K209

muts_strain <- muts_cap %>%
  filter(Strain == "K209")

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_cap %>%
  filter(Strain == "K209") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_cap %>%
  filter(Strain == "K209") %>%
  pull(Type)

muts_condition <- muts_cap %>%
  filter(Strain == "K209") %>%
  pull(Condition)

is_annotation <- muts_cap %>%
  filter(Strain == "K209") %>%
  pull(`IS mediated`)

second_seq <- muts_cap %>%
  filter(Strain == "K209") %>%
  pull(`Seq ID2`)

coords_cap$width <- coords_cap$end - coords_cap$begin

coords_cap_K209 <- coords_cap %>%
  filter(strain == "K209")

score <- muts_cap %>%
  filter(Strain == "K209") %>%
  pull(Freq)

mutations_gr <- GRanges("chr1", IRanges(position_muts, width = 1, names = muts_strain$names))
features <- GRanges("chr1", IRanges(as.vector(coords_cap_K209$begin), width = as.vector(coords_cap_K209$width)))
features$height <- c(0.05, 0.05 )

color <- replace(muts_annotation, muts_annotation == "Chromosome-Plasmid", "orange")
color <- replace(color, color == "non-synonymous", "black")
color <- replace(color, color == "Chromosome-Chromosome", "orange")
color <- replace(color, color == "synonymous", "white")
color <- replace(color, color == "intergenic", "grey")

second_seq<- replace(second_seq, second_seq == "IncL/M(pOXA-48)_1_pOXA-48", "dodgerblue3")
second_seq<- replace(second_seq, second_seq == "-", "black")
second_seq<- replace(second_seq, second_seq == "Chromosome", "black")
second_seq<- replace(second_seq, second_seq == "IncHI1B_1_pNDM-MAR", "darkmagenta")
second_seq<- replace(second_seq, second_seq == "IncFIA(HI1)_1_HI1/IncFII_1_pKP91", "darkmagenta")
side <- replace(muts_condition, muts_condition == "No pOXA-48", "top")
side <- replace(side, side != "top", "bottom")
mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
mutations_gr$score <- score
mutations_gr$label.parameter.rot <- 45

par(mar=c(0.5, 0,0,.5))

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(1704000, 1729000)), rescale = FALSE)


# CF operons

muts_cf <- read_xlsx("CF_operons_for_lolliplot.xlsx", 1)
coords_cf <- read_xlsx("CF_operons_for_lolliplot.xlsx", 2)

# CF12

muts_strain <- muts_cf %>%
  filter(Strain == "CF12")

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_cf %>%
  filter(Strain == "CF12") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_cf %>%
  filter(Strain == "CF12") %>%
  pull(Type)

muts_condition <- muts_cf %>%
  filter(Strain == "CF12") %>%
  pull(Condition)

is_annotation <- muts_cf %>%
  filter(Strain == "CF12") %>%
  pull(`IS mediated`)

second_seq <- muts_cf %>%
  filter(Strain == "CF12") %>%
  pull(`Seq ID2`)

coords_cf$width <- coords_cf$end - coords_cf$begin

coords_cf_CF12 <- coords_cf %>%
  filter(strain == "CF12")

score <- muts_cf %>%
  filter(Strain == "CF12") %>%
  pull(Freq)

mutations_gr <- GRanges("chr1", IRanges(position_muts, width = 1, names = muts_strain$names))
features <- GRanges("chr1", IRanges(as.vector(coords_cf_CF12$begin), width = as.vector(coords_cf_CF12$width)))
features$height <- c(0.05, 0.05 )

color <- replace(muts_annotation, muts_annotation == "Chromosome-Plasmid", "orange")
color <- replace(color, color == "Chromosome-Chromosome", "darkgreen")
color <- replace(color, color == "non-synonymous", "black")
color <- replace(color, color == "synonymous", "white")
color <- replace(color, color == "intergenic", "grey")

second_seq<- replace(second_seq, second_seq == "IncL/M(pOXA-48)_1_pOXA-48", "dodgerblue3")
second_seq<- replace(second_seq, second_seq == "-", "black")
second_seq<- replace(second_seq, second_seq == "Chromosome", "black")
#second_seq<- replace(second_seq, second_seq == "IncFIB(K)_1_Kpn3/IncFII_1_pKP91", "darkmagenta")
side <- replace(muts_condition, muts_condition == "No pOXA-48", "top")
side <- replace(side, side != "top", "bottom")
mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
mutations_gr$score <- score
mutations_gr$label.parameter.rot <- 45

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(3895000, 3897500)), legend = "legend")

# CF13

muts_strain <- muts_cf %>%
  filter(Strain == "CF13")

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_cf %>%
  filter(Strain == "CF13") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_cf %>%
  filter(Strain == "CF13") %>%
  pull(Type)

muts_condition <- muts_cf %>%
  filter(Strain == "CF13") %>%
  pull(Condition)

is_annotation <- muts_cf %>%
  filter(Strain == "CF13") %>%
  pull(`IS mediated`)

second_seq <- muts_cf %>%
  filter(Strain == "CF13") %>%
  pull(`Seq ID2`)

coords_cf$width <- coords_cf$end - coords_cf$begin

coords_cf_CF13 <- coords_cf %>%
  filter(strain == "CF13")

score <- muts_cf %>%
  filter(Strain == "CF13") %>%
  pull(Freq)

mutations_gr <- GRanges("chr1", IRanges(position_muts, width = 1, names = muts_strain$names))
features <- GRanges("chr1", IRanges(as.vector(coords_cf_CF13$begin), width = as.vector(coords_cf_CF13$width)))
features$height <- c(0.05, 0.05 )

color <- replace(muts_annotation, muts_annotation == "Chromosome-Plasmid", "orange")
color <- replace(color, color == "Chromosome-Chromosome", "orange")
color <- replace(color, color == "non-synonymous", "black")
color <- replace(color, color == "synonymous", "white")
color <- replace(color, color == "intergenic", "grey")

second_seq<- replace(second_seq, second_seq == "IncL/M(pOXA-48)_1_pOXA-48", "dodgerblue3")
second_seq<- replace(second_seq, second_seq == "-", "black")
second_seq<- replace(second_seq, second_seq == "Chromosome", "black")
#second_seq<- replace(second_seq, second_seq == "IncFIB(K)_1_Kpn3/IncFII_1_pKP91", "darkmagenta")
side <- replace(muts_condition, muts_condition == "No pOXA-48", "top")
side <- replace(side, side != "top", "bottom")
mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
mutations_gr$score <- score
mutations_gr$label.parameter.rot <- 45

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(3929300, 3931900)), legend = "legend")

# C021 operon

muts_cf <- read_xlsx("C021_operon.xlsx", 1)
coords_cf <- read_xlsx("C021_operon.xlsx", 2)

# C021

muts_strain <- muts_cf %>%
  filter(Strain == "C021")

muts_strain <- muts_strain %>%
  mutate(names = if_else(Type == "Chromosome-Plasmid", substr(Annotation2,1,22), paste(Gene1,"(",Position2,")")))

position_muts <- muts_cf %>%
  filter(Strain == "C021") %>%
  pull(Position1)

position_muts <- as.numeric(as.vector(position_muts))

muts_annotation <- muts_cf %>%
  filter(Strain == "C021") %>%
  pull(Type)

muts_condition <- muts_cf %>%
  filter(Strain == "C021") %>%
  pull(Condition)

is_annotation <- muts_cf %>%
  filter(Strain == "C021") %>%
  pull(`IS mediated`)

second_seq <- muts_cf %>%
  filter(Strain == "C021") %>%
  pull(`Seq ID2`)

coords_cf$width <- coords_cf$end - coords_cf$begin

coords_cf_C021 <- coords_cf %>%
  filter(strain == "C021")

score <- muts_cf %>%
  filter(Strain == "C021") %>%
  pull(Freq)

mutations_gr <- GRanges("chr1", IRanges(position_muts, width = 1, names = muts_strain$names))
features <- GRanges("chr1", IRanges(as.vector(coords_cf_C021$begin), width = as.vector(coords_cf_C021$width)))
features$height <- c(0.05, 0.05 )

color <- replace(muts_annotation, muts_annotation == "Chromosome-Plasmid", "orange")
color <- replace(color, color == "Chromosome-Chromosome", "darkgreen")
color <- replace(color, color == "non-synonymous", "black")
color <- replace(color, color == "synonymous", "white")
color <- replace(color, color == "intergenic", "grey")

second_seq<- replace(second_seq, second_seq == "IncL/M(pOXA-48)_1_pOXA-48", "dodgerblue3")
second_seq<- replace(second_seq, second_seq == "-", "black")
second_seq<- replace(second_seq, second_seq == "Chromosome", "black")
#second_seq<- replace(second_seq, second_seq == "IncFIB(K)_1_Kpn3/IncFII_1_pKP91", "darkmagenta")
side <- replace(muts_condition, muts_condition == "No pOXA-48", "top")
side <- replace(side, side != "top", "bottom")
mutations_gr$color <- color
mutations_gr$SNPsideID <- side
mutations_gr$border <- second_seq
mutations_gr$score <- score
mutations_gr$label.parameter.rot <- 45

lolliplot(mutations_gr, features, ranges = GRanges("chr1", IRanges(1092000, 1094500)), legend = "legend")

