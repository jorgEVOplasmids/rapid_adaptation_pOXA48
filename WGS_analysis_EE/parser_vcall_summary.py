#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 15:13:21 2023

@author: Jorge Sastre Domínguez

Script for summarising the results obtained with the variant calling of the pOXA-48 experimental evolution.

The aim of this code is to summarize the variants identified using breseq after applying the filters indicated 
in the parser_EE.py.

First, this script takes the summary tables (parsed mutations and parsed NJ sheets) and identifies each type of mutation
for each condition, with the intention of getting at the end the number of mutations per each condition depending
on the strain/species and also their frequency.

Appart from summarizing the number and frequency of each SNP or NJ, it also summarizes the type of event, i.e., Non-synonymous,
Synonymous, Intergenic, etc.


"""

import pandas as pd

# LOAD DATAFRAME

strain = "CF13"
species = "C. freundii"

xlsx = pd.ExcelFile("/home/jorge/Documents/TFM/clinical_samples/R_analysis/xlsx_breseqs/parsed_mutations/110822/parsed_"+strain+".xlsx")

parsed_muts = pd.read_excel(xlsx, "parsed mutations")
parsed_NJs = pd.read_excel(xlsx, "parsed NJ")

# Declare empty columns for later dataframe

seq_id_mutated = []
annotation_mutation = []
position_mutated = []
type_mutation = []
event_snp = []
gene_mutated = []
description_mutated = []
replicate_mutated = []
strain_mutated = []
species_mutated = []
condition_mutated = []
freq_mutation = []

condition_reference = ["pOXA-48-Anc", "pOXA-48", "pOXA-48", "pOXA-48", 
                       "No pOXA-48-Anc", "No pOXA-48", "No pOXA-48", 
                       "No pOXA-48", "pOXA-48+AMC", "pOXA-48+AMC","pOXA-48+AMC"]

# Iterate through SNPs

for index, row in parsed_muts.iterrows():
    # Get position of the mutation in each row in numeric object type
    if ":" not in row[2]:
        position = int(row[2].replace(",","")) # position of the mutation
    else:
        position = int(row[2].split(":")[0].replace(",","")) # position of the indel
    
    # Get only those mutations which are not in the ancestors (avoid possible mutations fixed at day 1)
    if row[7] == 0 and row[11] == 0:
        for i in range(7,18): # Iterate through replicates columns in each row
            if row[i] > 15:
                position_mutated.append(position)
                if row[1] == 1:
                    seq_id_mutated.append("Chromosome")
                else:
                    seq_id_mutated.append("Plasmid_"+str(row[1]))
                condition_mutated.append(condition_reference[i-7])
                replicate_mutated.append(parsed_muts.columns[i])
                strain_mutated.append(strain)
                species_mutated.append(species)
                gene_mutated.append(row[5])
                description_mutated.append(row[6])
                freq_mutation.append(row[i])
                annotation_mutation.append(row[3])
                #print(position, condition_reference[i-7], parsed_muts.columns[i], row[5])
                annot = row[4].split("(")
                if annot[0][0] == annot[0][-1]:
                    type_mut = "synonymous"
                if "intergenic" in annot[0]:
                    type_mut = "intergenic"
                if "pseudogene" in annot[0]:
                    type_mut = "pseudogene"
                if "pseudogene" not in annot[0] and "intergenic" not in annot[0] and annot[0][0] != annot[0][-1]:
                    type_mut = "non-synonymous"
                type_mutation.append(type_mut)
                #print(annot)
                if "Δ" in row[3] or "+" in row[3]:
                    event_snp.append("INDEL")
                else:
                    event_snp.append("SNP")

# Create table from the columns 

dict_muts = {"Seq ID": seq_id_mutated, "Event": event_snp, "Type": type_mutation, "Position": position_mutated, "Mutation":annotation_mutation,
             "Gene": gene_mutated, "Annotation": description_mutated, "Replicate": replicate_mutated, "Strain": strain_mutated,
             "Species": species_mutated, "Condition": condition_mutated,"Freq": freq_mutation}
    
df_muts_parsed = pd.DataFrame(dict_muts)

#print(df_muts_parsed)

seq_id1_NJ = []
seq_id2_NJ = []
position1_NJ = []
position2_NJ = []
type_NJ = []
event_NJ = []
gene1_NJ = []
gene2_NJ = []
description1_NJ = []
description2_NJ = []
replicate_NJ = []
strain_NJ = []
species_NJ = []
condition_NJ = []
freq_NJ = []
is_mediated_NJ = []

# Create condition to identify IS mediated NJ events

IS_gene = ["family transposase", "DGQHR"]

# Iterate through NJs

for index, row in parsed_NJs.iterrows():
    contig_1 = row[0]
    contig_2 = row[3]
    position_1 = int(str(row[1]).replace(" ","").replace("=", ""))
    position_2 = int(str(row[4]).replace(" ","").replace("=", ""))
    if row[15] == 0 and row[19] == 0:
        for i in range(15,26):
            if row[i] > 15:
                position1_NJ.append(row[1])
                position2_NJ.append(row[4])
                if row[0] == 1:
                    contig_1 = "Chromosome"
                else:
                    contig_1 = "Plasmid_"+str(row[0])
                if row[3] == 1:
                    contig_2 = "Chromosome"
                else:
                    contig_2 = "Plasmid_"+str(row[3])
                seq_id1_NJ.append(contig_1)
                seq_id2_NJ.append(contig_2)
                type_NJ.append(contig_1+"-"+contig_2)
                event_NJ.append("NJ")
                gene1_NJ.append(row[10])
                description1_NJ.append(row[11])
                gene2_NJ.append(row[13])
                description2_NJ.append(row[14])
                replicate_NJ.append(parsed_NJs.columns[i])
                strain_NJ.append(strain)
                species_NJ.append(species)
                condition_NJ.append(condition_reference[i-15])
                freq_NJ.append(row[i])
                if "family transposase" in row[11] or "family transposase" in row[14] or "DGQHR" in row[11] or "DGQHR" in row[14]:
                    is_mediated_NJ.append("yes")
                else:
                    is_mediated_NJ.append("no")
                    

dict_NJs = {"Seq ID1": seq_id1_NJ, "Seq ID2": seq_id2_NJ, "Event": event_NJ, "Type": type_NJ, 
            "Position1": position1_NJ, "Position2": position2_NJ, "Gene1": gene1_NJ, "Gene2": gene2_NJ,
        "Annotation1": description1_NJ, "Annotation2": description2_NJ, "Replicate": replicate_NJ, 
        "Strain": strain_NJ, "Species": species_NJ, "Condition": condition_NJ,"Freq": freq_NJ, "IS mediated": is_mediated_NJ} 
    
df_NJs_parsed = pd.DataFrame(dict_NJs)

#print(df_NJs_parsed)
outname = "/home/jorge/Documents/important_docs/draft_expev/figures_draft_expev/summary_parsed_vcall_w_plasmids/summary_parsed_"+strain+".xlsx"

writer = pd.ExcelWriter(outname, engine='xlsxwriter')
df_muts_parsed.to_excel(writer, sheet_name='parsed muts table', index=False)
df_NJs_parsed.to_excel(writer, sheet_name='parsed NJs table', index=False)
writer.save()

