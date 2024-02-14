#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import pandas as pd
import pprint
import openpyxl

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""
Created on Thu Aug 11 14:53:53 2022

@author: Jorge Sastre Domínguez

Script to reanalyze experimental evolution data with restricting filters. This script takes as input the parsed xlsx output by 
breseq_index_parser.py and gives as output an xlsx file per strain with the filtered mutations applying diverse filters. 
The input xlsx needs to be ordered in a specific shape, as explained here:
    
    - Sheet2 is "mutations" spreadsheet, containing the SNPs and called JC/MC by breseq. The first 7 columns must be: 
        Evidence, Seq ID, Position, Mutation, Annotation, Gene and Description. The following 11 columns containing the
        expev replicates of the strain must have the following order: ancestor with pOXA-48, 3x evolved with pOXA-48 
        and no Ab, ancestor without pOXA-48, 3x evolved without pOXA-48, 3x evolved with pOXA-48 and Ab. Last column indicates
        the length of the assembled reference sequence and can be ommited.
    
    - Sheet7 must be "JC full" spreadsheet, and contain the New Junction evidences called by breseq. In this case the first 
    15 columns must be: Seq ID 1, Position 1, Reads (Cov) 1, Seq ID 2, Position 2, Reads (Cov) 2, Reads (Cov), Score, Skew,
    Annotation 1, Gene 1, Product 1, Annotation 2, Gene 2, Product 2. The following columns of the expev replicates must have
    the same order as the mutations sheet.
    
The parameters employed to filter out events are:
    
    - Mutations present in ancestors (at nt level) are filtered out.
    
    - Genes mutated more than twice in the same replicate are discarded due to possible false positive calls.
    
    - Parallel mutations, i.e. mutations present in more than one replicate per strain, are filtered out if the frequency is less than 15%.
    
    - Non-parallel mutations are filtered out if the frequency is less than 15%.
    
    - Mutations present in regions close to intergenic mutations (+-10 bp) are filtered out.
    
    - NJs present in ancestors and in surrounding regions close to intergenic mutations are filtered out. Freq filter is again > 15%.
    
For the output mutations, the contig into which the event is located is taken into account, as well as the position onto 
NJ evidences (left/right and strand).

The output xlsx file is structured in multiple sheets:
    
    - Sheet1 ("Contig info"): contains the information of the ID and length of both chromosome and plasmids.
    
    - Sheet2 ("Parsed mutations"): SNPs and INDELs which pass the filters applyed for non-parallel mutations. They are annotated at bp level.
    
    - Sheet3 ("Parallel mutations"): SNPs and INDELs which pass the filters applyed for parallel mutations. These are summarized per gene.
    
    - Sheet4 ("Ancestors ompensatory candidates"): SNPs and INDELs present in ancestors which can have compensated the plasmid cost from the beginning
    of the experimental evolution. All high frequency mutations in both ancestors and present in samples evolved with pOXA-48 are included.
    
    - Sheet5 ("Parsed NJ"): NJ evidences which pass the filters applied for non-parallel NJ. They are annotated at bp level.
    
    - Sheet6 ("Summary NJ"): Summary of NJ depending on the contigs. Rearrangements within the chromosome or each plasmid, as well as junctions of
    the chromosome with pOXA-48 or other plasmids are included. A summary of NJ caused by Insertion Sequences is also included in this sheet.
    
    - Sheet7 ("NJ per gene 1"): Parallel evolution caused by NJ in genes annotated as the extreme 1 of the NJ.
    
    - Sheet8 ("NJ per gene 2"): Parallel evolution caused by NJ in genes annotated as the extreme 2 of the NJ.
    
""")

parser.add_argument('Directory', action="store",
                    help="Master directory containing parsed xlsx output by breseq_index_parser.py")
parser.add_argument('Strain', action="store",
                    help="Strain name of the xlsx file to apply the filters")

args = parser.parse_args()

DirName = args.Directory
DirValidation = os.path.isdir(DirName)
StrainName = args.Strain

if DirValidation not in ['True', True]:
    print("Error. Please, enter a valid master directory to parse. Try --help to see the parser's description. Exiting!")
    exit(1)

FileNames = os.listdir(DirName)

####################################################### DECLARE GLOBAL VARIABLES #######################################################

####################################################### MUTATIONS SPREADSHEET #######################################################

# Diccionary that will contain SNPs present in ancestors from "mutations" sheet
# Key = chromosome/plasmid; Value = list of positions (int) mutated in the chromosome/plasmid in the ancestors
anc_muts = {}

# Dictionary that will contain n of times that X gene is mutated per strain from "mutations" sheet
# Key = gene; Value = list of replicates in which the gene is mutated
genes_per_rep = {}

# Dictionary of positions in which intergenic mutations fall in "mutations" sheet
# Key = replicate; Value = dictionary of surrounding positions (int) to intergenig mutations per chromosome/plasmid
interg_pos_dic = {}

# Dictionary of positions mutated in samples and not present in ancestors
# Will be useful to parse by window length regions which contain more than 2 mutations
pos_mut = {}

# Dictionary with mutations present in overmutated regions
# Key = replicate; Value = dictionary of positions (str) mutated per chr/ plasmid in windows of 400 bp
pos_overmut = {}

# Initialize empty lists what will be converted in the parsed data frames later
parsed_muts = []
parsed_ordered_muts = []
per_gene_muts = []
parallel_clean = []
parallel_parsed = []
compensatory_candidates = []

####################################################### NJ SPREADSHEET #######################################################

# Dictionary that will contain NJ positions present in ancestors from "NJ" sheet
# Key = chromosome/plasmid, Value = list of positions (int) targeted by NJ in the chromosome/plasmid in the ancestors
anc_nj = {}

# Dictionary that will contain NJ events depending on the EE condition
# Key = condition, Value = list of NJ movements chr-chr, chr_oxa, chr_pl, pl_pl, oxa_oxa, chr_isoxa

summary_NJ = {}

# Initialize empty lists that will be converted in the parsed data frames later

parsed_NJ = []
summary_NJ_df = []
parsed_NJ_ordered = []
per_gene_NJ = []
per_gene_NJ_2 = []

# Open breseq parsed xlsx input file for the specified strain
# xlsx files obtained with breseq_index_parser.py are placed in /input/xlsx_raw

for i, element in enumerate(FileNames):
    if StrainName in FileNames[i]:
        DirName=str(DirName)
        print("Parsing "+DirName+FileNames[i])
        xlsx = pd.ExcelFile(DirName+FileNames[i])
        contig_info = pd.read_excel(xlsx, "contig info")
        muts = pd.read_excel(xlsx, "mutations")
        NJ = pd.read_excel(xlsx, "JC full")
        #print(len(muts))
        #print(len(NJ))
        
        ####################################################### PARSE MUTATIONS SPREADSHEET #######################################################
        
        ####################################################### Get dicts and lists for discarding false or noisy mutations #######################################################
        
        print("Getting filters for parsing mutations sheet of "+StrainName)
        
        for index, row in muts.iterrows():
            # Get position of the mutation in each row in numeric object type
            if ":" not in row[2]:
                position = int(row[2].replace(",","")) # position of the mutation
            else:
                position = int(row[2].split(":")[0].replace(",","")) # position of the indel
            
            # Get dictionary of positions mutated in ancestors
            if row[7] != 0 or row[11] != 0:
                if row[1] not in anc_muts.keys():
                    anc_muts[row[1]] = []
                    anc_muts[row[1]].append(position)
                else:
                    anc_muts[row[1]].append(position)
                    
            # Get dictionary with counts of mutations per gene in each replicate
            
            for i in range(7,18): # Iterate through replicates columns in each row
                if row[i] != 0:
                    if row[5] not in genes_per_rep:
                        genes_per_rep[row[5]] = []
                        genes_per_rep[row[5]].append(muts.columns.values.tolist()[i])
                    else:
                        genes_per_rep[row[5]].append(muts.columns.values.tolist()[i])
                        
            # Get regions surrounding to intergenic mutations
            # To deal with numbers, str positions are changed to int
            
            for i in range(7,18): # Iterate through replicates columns in each row
                if row[i] != 0:
                    repl_name = muts.columns.values.tolist()[i] # name of the replicate mutated
                    
                    # Get dictionary with regions surrounding intergenic mutations
                    
                    if "Δ" not in row[3] and "intergenic" in row[4]:
                        interg_pos = position
                        if repl_name not in interg_pos_dic.keys():
                            interg_pos_dic[repl_name] = {} # declare empty nested dictionary
                            if row[1] not in interg_pos_dic[repl_name].keys(): # if contig not in replicate nested dictionary
                                interg_pos_dic[repl_name][row[1]] = []
                                for n in range(interg_pos-10, interg_pos+11):
                                    interg_pos_dic[repl_name][row[1]].append(n)
                            else:                    
                                for n in range(interg_pos-10, interg_pos+11):
                                    interg_pos_dic[repl_name][row[1]].append(n)
                        else:
                            if row[1] not in interg_pos_dic[repl_name].keys():
                                interg_pos_dic[repl_name][row[1]] = []
                                for n in range(interg_pos-10, interg_pos+11):
                                    interg_pos_dic[repl_name][row[1]].append(n)
                            else:                    
                                for n in range(interg_pos-10, interg_pos+11):
                                    interg_pos_dic[repl_name][row[1]].append(n)
                        
                    # Get positions mutated in each contig per replicate
                    
                    if repl_name not in pos_mut.keys():
                        pos_mut[repl_name] = {} # declare empty nested dictionary
                        if row[1] not in pos_mut[repl_name].keys(): # if contig not in replicate nested dictionary
                            pos_mut[repl_name][row[1]] = []
                            pos_mut[repl_name][row[1]].append(position)
                        else:
                            pos_mut[repl_name][row[1]].append(position)
                    else:
                        if row[1] not in pos_mut[repl_name].keys(): # if contig not in replicate nested dictionary
                            pos_mut[repl_name][row[1]] = []
                            pos_mut[repl_name][row[1]].append(position)
                        else:
                            pos_mut[repl_name][row[1]].append(position)
                    
                # Remove duplicates from intergenic dictionary
                
                repl_name = muts.columns.values.tolist()[i]
                if repl_name in interg_pos_dic.keys():
                    if row[1]in interg_pos_dic[repl_name].keys():
                        interg_pos_dic[repl_name][row[1]]=list(set(interg_pos_dic[repl_name][row[1]]))
                        
        # Get positions of mutations contained in windows of 100 bp mutated more than twice                   
        
        for replicate in pos_mut.keys():
            for contig in pos_mut[replicate].keys(): # iterate theough all the replicates and contigs positions mutated
                for p_mutated in pos_mut[replicate][contig]:    
                    counter = 0
                    for p in range(p_mutated-200, p_mutated+201): # count mutations in window of 100 bp
                        if p in pos_mut[replicate][contig]:
                            counter += 1
                    if counter > 2:
                        
                        if replicate not in pos_overmut.keys():
                            pos_overmut[replicate] = {}
                            if contig not in pos_overmut[replicate].keys():
                                pos_overmut[replicate][contig] = []
                                pos_overmut[replicate][contig].append(p_mutated)
                            else:
                                pos_overmut[replicate][contig].append(p_mutated)
                        else:
                            if contig not in pos_overmut[replicate].keys():
                                pos_overmut[replicate][contig] = []
                                pos_overmut[replicate][contig].append(p_mutated)
                            else:
                                pos_overmut[replicate][contig].append(p_mutated)
        
        # Remove positions mutated in enriched windows dictionary if frequency is higher than 95% so that they won't be later filtered
        for index, row in muts.iterrows():
            if ":" not in row[2]:
                position = int(row[2].replace(",",""))
            else:
                position = int(row[2].split(":")[0].replace(",",""))
            
            for i in range(7,18):
                repl_name = muts.columns.values.tolist()[i]
                if row[i] >= 95:
                    if repl_name in pos_overmut.keys() and row[1] in pos_overmut[repl_name]:
                        if position in pos_overmut[repl_name][row[1]]:
                            pos_overmut[repl_name][row[1]].remove(position)
        
        ####################################################### Get mutations which pass the filters #######################################################
        
        print("Filtering mutations of "+StrainName+" which fulfill conditions")
        
        for index, row in muts.iterrows():
            counter_parallel = 0
            if ":" not in row[2]:
                position = int(row[2].replace(",",""))
            else:
                position = int(row[2].split(":")[0].replace(",",""))
            for i in range(7,18):
                if row[i] >=  15: # Filter out mutations with lower frequencies than 15%
                    counter_parallel += 1
                    repl_name = muts.columns.values.tolist()[i]
                    # Declare empty list of ancestral mutations in contigs for avoiding key errors
                    if row[1] not in anc_muts.keys():
                        anc_muts[row[1]] = []
                    # Declare empty nested dictionaries for replicates with no intergenic mutations for avoiding key errors
                    if repl_name not in interg_pos_dic.keys():
                        interg_pos_dic[repl_name] = {}
                    # Declare empty lists for contigs with no intergenic mutations for avoiding key errors
                    if row[1] not in interg_pos_dic[repl_name].keys():
                            interg_pos_dic[repl_name][row[1]] = []
                    # Delete intergenic mutated positions itself for later parsing
                    if "Δ" not in row[3] and "intergenic" in row[4]:
                        interg_pos = position
                        interg_pos_dic[repl_name][row[1]].remove(interg_pos)
                    
                    # Declare empty nested dictionaries and lists in overenriched regions dictionaries for avoiding key errors
                    if repl_name not in pos_overmut.keys():
                            pos_overmut[repl_name] = {}
                    # Declare empty lists for contigs with no overenriched regions for avoiding key errors
                    if row[1] not in pos_overmut[repl_name].keys():
                        pos_overmut[repl_name][row[1]] = []
                        
                    # Filter out positions mutated in ancestors, in overenriched regions, surrounding intergenic mutations
                    if position not in anc_muts[row[1]] and position not in pos_overmut[repl_name][row[1]] and position not in interg_pos_dic[repl_name][row[1]]:
                        # Filter out mutations in genes mutated more than twice in the same replicate in low frequencies
                        if genes_per_rep[row[5]].count(repl_name) >=2 and row[i] >= 50: # if freq is really high, the mutation is included despite the gene being highly targeted
                            parsed_muts.append(row)
                        elif genes_per_rep[row[5]].count(repl_name) <= 2: # if the gene is targeted twice or less in the same replicate, it can be included
                            parsed_muts.append(row)
                    
                    # Include parallel mutations if they are not present in the ancestors. The rest of parameters are not taken into account for these
                    if position not in anc_muts[row[1]] and counter_parallel >= 2:
                        parsed_muts.append(row)
        
        ####################################################### Parallel evolution analysis (SNPs) #######################################################
        print("Summarizing mutations of "+StrainName+" per gene")
        # Reorder dataframe by gene name
        muts_ordered = muts.sort_values(by=["Gene"])
        
        # If mutations not in ancestors neither in overenriched windows of the replicate, keep them in the ordered dataframe
        for index, row in muts_ordered.iterrows():
            if ":" not in row[2]:
                position = int(row[2].replace(",",""))
            else:
                position = int(row[2].split(":")[0].replace(",",""))
            for i in range(7,18):
                if row[i] != 0:
                    repl_name = muts.columns.values.tolist()[i]
                    # Declare empty list of ancestral mutations in contigs for avoiding key errors
                    if row[1] not in anc_muts.keys():
                        anc_muts[row[1]] = []
                     # Declare empty lists for contigs with no overenriched regions for avoiding key errors
                    if row[1] not in pos_overmut[repl_name].keys():
                        pos_overmut[repl_name][row[1]] = []
                    if position not in anc_muts[row[1]] and position not in pos_overmut[repl_name][row[1]]:
                        parsed_ordered_muts.append(row)
        # Convert to dataframe ordered and already filtered
        parsed_ordered_muts_df = pd.DataFrame(parsed_ordered_muts)
        parsed_ordered_muts_df_def = parsed_ordered_muts_df.drop_duplicates(keep = "first")
        
        # Work on the filtered and ordered df to sum frequencies affecting the same genes
        genes_mutated_list = parsed_ordered_muts_df_def["Gene"].tolist()
        multitargeted_genes = []
        for gene in genes_mutated_list: # keep genes mutated more than once that passed the filters
            if genes_mutated_list.count(gene) >= 2:
                multitargeted_genes.append(gene)
        multitargeted_genes = list(set(multitargeted_genes))
        # Get header for parallel spreadsheet
        header_parallel = parsed_ordered_muts_df_def.columns.values.tolist()
        header_parallel.remove("Evidence")
        header_parallel.remove("Position")
        header_parallel.remove("Mutation")
        header_parallel.remove("Annotation")
        per_gene_muts.append(header_parallel)
        # Get parallel mutation row for each gene
        for gene in multitargeted_genes:
            parallel_row = []
            gene_df = parsed_ordered_muts_df_def[parsed_ordered_muts_df_def["Gene"] == gene]
            parallel_row.append(gene_df.iloc[0]["Seq ID"])
            parallel_row.append(gene_df.iloc[0]["Gene"])
            parallel_row.append(gene_df.iloc[0]["Description"])
            for i in range(7, 18):
                repl_name = repl_name = parsed_ordered_muts_df_def.columns.values.tolist()[i]
                parallel_row.append(gene_df[repl_name].sum())
            parallel_row.append(gene_df.iloc[0]["Length"])
            per_gene_muts.append(parallel_row)

            #print(gene_df)
        #for index, row in parsed_ordered_muts_df_def.iterrows():
        per_gene_muts_df = pd.DataFrame(per_gene_muts)
        per_gene_muts_to_parse = per_gene_muts_df.rename(columns=per_gene_muts_df.iloc[0]).drop(per_gene_muts_df.index[0])
        
        # Remove false positives due to hypermutator strains
        for index, row in per_gene_muts_to_parse.iterrows():
            c = 0
            for i in range(3, 14):
                if row[i] != 0:
                    c += 1
            if c >= 2:
                parallel_clean.append(row)
                
        parallel_clean_df = pd.DataFrame(parallel_clean)
        # Filter by frequency
        for index, row in parallel_clean_df.iterrows():
            for i in range(3, 14):
                if row[i] >= 15:
                    parallel_parsed.append(row)
        
        ####################################################### Compensatory candidates analysis #######################################################
        
        print("Annotating possible compensatory candidates present in "+StrainName+" ancestors")
        
        for index, row in muts.iterrows():
            if ":" not in row[2]:
                position = int(row[2].replace(",",""))
            else:
                position = int(row[2].split(":")[0].replace(",",""))
            if position in anc_muts[row[1]]:
                for i in range(8,11):
                    repl_name = parsed_ordered_muts_df_def.columns.values.tolist()[i]
                    if row[i] >= 33 and position not in pos_overmut[repl_name][row[1]] and genes_per_rep[row[5]].count(repl_name) <= 2:
                        compensatory_candidates.append(row)
                for i in range(15, 18):
                    repl_name = repl_name = parsed_ordered_muts_df_def.columns.values.tolist()[i]
                    if row[i] >= 33 and position not in pos_overmut[repl_name][row[1]] and genes_per_rep[row[5]].count(repl_name) <= 2:
                        compensatory_candidates.append(row)
        
        ####################################################### PARSE NEW JUNCTIONS SPREADSHEET #######################################################
        
        ####################################################### Get dicts and lists for discarding false or noisy NJ #######################################################

        print("Getting filters for parsing NJ sheet of "+StrainName)
        
        # Get contig information to identify events into the chromosome, between plasmids and into plasmids
        
        plasmids_id = []
        for index, row in contig_info.iterrows():
            if row[1] > 4000000:
                chr_id = row[0]
            elif row[1] > 63000 and row[1] < 66000:
                poxa_id = row[0]
            else:
                plasmids_id.append(row[0])
        
        #print(chr_id, poxa_id, plasmids_id)
        
        # Get positions in which the ancestors show NJ
        
        for index, row in NJ.iterrows():
            contig_1 = row[0]
            contig_2 = row[3]
            position_1 = int(str(row[1]).replace(" ","").replace("=", ""))
            position_2 = int(str(row[4]).replace(" ","").replace("=", ""))
            #print(contig_1, position_1, contig_2, position_2)
            
            if row[15] != 0 or row[19] != 0:
                if contig_1 not in anc_nj.keys():
                    anc_nj[contig_1] = []
                    anc_nj[contig_1].append(position_1)
                else:
                    anc_nj[contig_1].append(position_1)
                if contig_2 not in anc_nj.keys():
                    anc_nj[contig_2] = []
                    anc_nj[contig_2].append(position_2)
                else:
                    anc_nj[contig_2].append(position_2)

        ####################################################### Filter out NJ not passing filters #######################################################
        
        print("Filtering NJ of "+StrainName+" which fulfill conditions")
            
        # Fuse NJ events repeated due to duplications caused after IS insertions
        # Consider contigs, positions, and extremes in each 2 consecutive rows mutated in the same replicate
        # NJs present in ancestors are already filtered out
        
        for index, row in NJ.iterrows():
            contig_1 = row[0]
            contig_2 = row[3]
            position_1 = int(str(row[1]).replace(" ","").replace("=", ""))
            position_2 = int(str(row[4]).replace(" ","").replace("=", ""))
            coverage_NJ = int(row[6].split(" ")[0])
            
            # Declare empty lists for avoiding key errors
            if contig_1 not in anc_nj.keys():
                anc_nj[contig_1] = []
            if contig_2 not in anc_nj.keys():
                anc_nj[contig_2] = []
            
            for i in range(15,26): # parse by replicate
                repl_name = NJ.columns.values.tolist()[i]
                if row[i] != 0 and position_1 != 1:
                    #print(repl_name, contig_1, contig_2)
                    # Declare empty lists for contigs with no intergenic mutations for avoiding key errors
                    if contig_1 not in interg_pos_dic[repl_name].keys():
                            interg_pos_dic[repl_name][contig_1] = []
                    if contig_2 not in interg_pos_dic[repl_name].keys():
                            interg_pos_dic[repl_name][contig_2] = []
                    if position_1 not in anc_nj[contig_1] and position_2 not in anc_nj[contig_2] and position_1 not in interg_pos_dic[repl_name][contig_1] and position_2 not in interg_pos_dic[repl_name][contig_2] and coverage_NJ >= 10: # filter out those NJ in ancestors and surrounding intergenic SNPs
                        #print(NJ.iloc[index])
                        #print(NJ.iloc[index-1])
                        # Merge rows which indicate same NJ but repeated due to IS duplication event after integration
                        if NJ.iloc[index]["Seq ID 1"] == NJ.iloc[index-1]["Seq ID 1"] and NJ.iloc[index]["Seq ID 2"] == NJ.iloc[index-1]["Seq ID 2"] and NJ.iloc[index]["Product 1"] == NJ.iloc[index-1]["Product 1"] and "IS" in NJ.iloc[index]["Product 2"] or "DGQHR" in NJ.iloc[index]["Product 2"]:
                            dif_targets = int(str(NJ.iloc[index]["Position 1"]).replace(" ","").replace("=", "")) - int(str(NJ.iloc[index-1]["Position 1"]).replace(" ","").replace("=", ""))
                            if  dif_targets < 15 and dif_targets > 1:
                                row[i] = (NJ.iloc[index][i] + NJ.iloc[index-1][i])/2 # Freq of NJ is assigned to the mean of the annotated event duplicated
                                if row[i] >= 15 and coverage_NJ >= 10: # Filter by frequency
                                    parsed_NJ.append(row)
                            elif row[i] >= 15 and coverage_NJ >= 10:
                                parsed_NJ.append(row) # Add NJ similarly annotated but far away in the same gene
                        elif row[i] >= 15 and coverage_NJ >= 10:
                            parsed_NJ.append(row) # Add NJ which are not duplicated
        
        parsed_NJ_df = pd.DataFrame(parsed_NJ)
        parsed_NJ_df_def = parsed_NJ_df.drop_duplicates(keep = "first")
        parsed_NJ_df_def = parsed_NJ_df_def.reset_index()
        parsed_NJ_df_def.drop("index", axis = 1, inplace = True)
        #print(parsed_NJ_df_def)
        #print(parsed_NJ_df_def.drop("index", axis = 1))
        
        # Finish polishing NJ dataframe (finish removing rows that indicate the same NJ)
        rows_to_drop_NJ = []
        for i, row in parsed_NJ_df_def.iterrows():
            #print(i, row)
            contig_1 = row[0]
            contig_2 = row[3]
            position_1 = int(str(row[1]).replace(" ","").replace("=", ""))
            position_2 = int(str(row[4]).replace(" ","").replace("=", ""))
            if i > 0:
                #print(i)
                if parsed_NJ_df_def.iloc[i]["Seq ID 1"] == parsed_NJ_df_def.iloc[i-1]["Seq ID 1"] and parsed_NJ_df_def.iloc[i]["Seq ID 2"] == parsed_NJ_df_def.iloc[i-1]["Seq ID 2"] and parsed_NJ_df_def.iloc[i]["Product 1"] == parsed_NJ_df_def.iloc[i-1]["Product 1"]:
                    dif_targets = int(str(parsed_NJ_df_def.iloc[i]["Position 1"]).replace(" ","").replace("=", "")) - int(str(parsed_NJ_df_def.iloc[i-1]["Position 1"]).replace(" ","").replace("=", ""))
                    if  dif_targets < 15 and dif_targets > 1:
                        rows_to_drop_NJ.append(i-1)
        #print(rows_to_drop_NJ)
        parsed_NJ_df_def.drop(rows_to_drop_NJ, inplace = True)
        
        ####################################################### Get Rearrangements between chr-pl, pl-pl, differentiating pOXA/Ab and EE condition #######################################################
        # Iterate through filtered NJ dataframe
        for index, row in parsed_NJ_df_def.iterrows():
            contig_1 = row[0]
            contig_2 = row[3]
            position_1 = int(str(row[1]).replace(" ","").replace("=", ""))
            position_2 = int(str(row[4]).replace(" ","").replace("=", ""))
            # Consider the condition of each replicate to sum the events
            for i in range(15, 26):
                if row[i] != 0:
                    if i <= 18:
                        condition = "pOXA-48"
                    elif i >= 19 and i <= 22:
                        condition = "No pOXA-48"
                    elif i >= 23:
                        condition = "pOXA-48+Ab"
                    if condition not in summary_NJ.keys(): # Initialize condition row
                        summary_NJ[condition] = [0,0,0,0,0,0,0]
                    # Add events for each of them depending on condition
                    if contig_1 == chr_id and contig_2 == chr_id and i <= 18:
                        summary_NJ["pOXA-48"][0] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i <= 18:
                            summary_NJ["pOXA-48"][6] += 1
                    elif contig_1 == chr_id and contig_2 == poxa_id and i <= 18:
                        summary_NJ["pOXA-48"][1] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i <= 18:
                            summary_NJ["pOXA-48"][5] += 1
                            summary_NJ["pOXA-48"][6] += 1
                    elif contig_1 == chr_id and contig_2 in plasmids_id and i <= 18:
                        summary_NJ["pOXA-48"][2] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i <= 18:
                            summary_NJ["pOXA-48"][6] += 1
                    elif contig_1 in plasmids_id and contig_2 in plasmids_id and i <= 18:
                        summary_NJ["pOXA-48"][3] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i <= 18:
                            summary_NJ["pOXA-48"][6] += 1
                    elif contig_1 == poxa_id and contig_2 == poxa_id and i <= 18:
                        summary_NJ["pOXA-48"][4] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i <= 18:
                            summary_NJ["pOXA-48"][6] += 1
                        
                    if contig_1 == chr_id and contig_2 == chr_id and i >= 19 and i <= 22:
                        summary_NJ["No pOXA-48"][0] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i >= 19 and i <= 22:
                            summary_NJ["No pOXA-48"][6] += 1
                    elif contig_1 == chr_id and contig_2 == poxa_id and i >= 19 and i <= 22:
                        summary_NJ["No pOXA-48"][1] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i >= 19 and i <= 22:
                            summary_NJ["No pOXA-48"][5] += 1
                            summary_NJ["No pOXA-48"][6] += 1
                    elif contig_1 == chr_id and contig_2 in plasmids_id and i >= 19 and i <= 22:
                        summary_NJ["No pOXA-48"][2] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i >= 19 and i <= 22:
                            summary_NJ["No pOXA-48"][6] += 1
                    elif contig_1 in plasmids_id and contig_2 in plasmids_id and i >= 19 and i <= 22:
                        summary_NJ["No pOXA-48"][3] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i >= 19 and i <= 22:
                            summary_NJ["No pOXA-48"][6] += 1
                    elif contig_1 == poxa_id and contig_2 == poxa_id and i >= 19 and i <= 22:
                        summary_NJ["No pOXA-48"][4] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i >= 19 and i <= 22:
                            summary_NJ["No pOXA-48"][6] += 1
                        
                    if contig_1 == chr_id and contig_2 == chr_id and i >= 23:
                        summary_NJ["pOXA-48+Ab"][0] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i >= 23:
                            summary_NJ["pOXA-48+Ab"][6] += 1
                    elif contig_1 == chr_id and contig_2 == poxa_id and i >= 23:
                        summary_NJ["pOXA-48+Ab"][1] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i >= 23:
                            summary_NJ["pOXA-48+Ab"][5] += 1
                            summary_NJ["pOXA-48+Ab"][6] += 1
                    elif contig_1 == chr_id and contig_2 in plasmids_id and i >= 23:
                        summary_NJ["pOXA-48+Ab"][2] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i >= 23:
                            summary_NJ["pOXA-48+Ab"][6] += 1
                    elif contig_1 in plasmids_id and contig_2 in plasmids_id and i >= 23:
                        summary_NJ["pOXA-48+Ab"][3] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i >= 23:
                            summary_NJ["pOXA-48+Ab"][6] += 1
                    elif contig_1 == poxa_id and contig_2 == poxa_id and i >= 23:
                        summary_NJ["pOXA-48+Ab"][4] += 1
                        if "IS" in row[14] or "DGQHR" in row[14] and i >= 23:
                            summary_NJ["pOXA-48+Ab"][6] += 1
                        
        #pprint.pprint(summary_NJ)
        
        # Organize summary NJ dataframe
        
        header_summary_NJ = ["Condition", "Chr-chr", "Chr-pOXA", "Chr-pl", "Pl-pl", "pOXA-pOXA", "Chr-pOXA(IS)", "Total IS"]
        summary_NJ_df.append(header_summary_NJ)
        for cond in summary_NJ.keys():
            row_sum_NJ = [cond]
            row_sum_NJ = row_sum_NJ + summary_NJ[cond]
            summary_NJ_df.append(row_sum_NJ)
            
        summary_NJ_df_def = pd.DataFrame(summary_NJ_df)
        summary_NJ_df_def = summary_NJ_df_def.rename(columns=summary_NJ_df_def.iloc[0]).drop(summary_NJ_df_def.index[0])
        
        #print(summary_NJ_df_def)
        
        ####################################################### Parallel evolution analysis (NJ) #######################################################
        print("Summarizing NJ of "+StrainName+" per gene")
        # Order NJ dataframe by 1st gene
        NJ_ordered = NJ.sort_values(by = ["Gene 1"])
        
        # If NJ not in ancestors and not in surrounding intergenic mutations keep them in the ordered dataframe
        for index, row in NJ_ordered.iterrows():
            contig_1 = row[0]
            contig_2 = row[3]
            position_1 = int(str(row[1]).replace(" ","").replace("=", ""))
            position_2 = int(str(row[4]).replace(" ","").replace("=", ""))
            coverage_NJ = int(row[6].split(" ")[0])
            
            for i in range(15,26): # parse by replicate
                repl_name = NJ_ordered.columns.values.tolist()[i]
                if row[i] != 0 and position_1 not in anc_nj[contig_1] and position_1 not in interg_pos_dic[repl_name][contig_1]:
                    # Merge rows which indicate same NJ but repeated due to IS duplication event after integration
                    if NJ_ordered.iloc[index]["Seq ID 1"] == NJ_ordered.iloc[index-1]["Seq ID 1"] and NJ_ordered.iloc[index]["Seq ID 2"] == NJ_ordered.iloc[index-1]["Seq ID 2"] and NJ_ordered.iloc[index]["Product 1"] == NJ_ordered.iloc[index-1]["Product 1"] and "IS" in NJ_ordered.iloc[index]["Product 2"] or "DGQHR" in NJ_ordered.iloc[index]["Product 2"]:
                        dif_targets = int(str(NJ_ordered.iloc[index]["Position 1"]).replace(" ","").replace("=", "")) - int(str(NJ_ordered.iloc[index-1]["Position 1"]).replace(" ","").replace("=", ""))
                        if  dif_targets < 15 and dif_targets > 1:
                            row[i] = (NJ_ordered.iloc[index][i] + NJ_ordered.iloc[index-1][i])/2
                            if row[i] != 0 and coverage_NJ >= 10: # Filter by frequency
                                    parsed_NJ_ordered.append(row)
                        elif row[i] != 0 and coverage_NJ >= 10:
                            parsed_NJ_ordered.append(row) # Add NJ similarly annotated but far away in the same gene
                    elif row[i] != 0 and coverage_NJ >= 10:
                        parsed_NJ_ordered.append(row)
                if row[i] != 0 and position_2 not in anc_nj[contig_2] and position_2 not in interg_pos_dic[repl_name][contig_2]:
                    # Merge rows which indicate same NJ but repeated due to IS duplication event after integration
                    if NJ_ordered.iloc[index]["Seq ID 1"] == NJ_ordered.iloc[index-1]["Seq ID 1"] and NJ_ordered.iloc[index]["Seq ID 2"] == NJ_ordered.iloc[index-1]["Seq ID 2"] and NJ_ordered.iloc[index]["Product 1"] == NJ_ordered.iloc[index-1]["Product 1"] and "IS" in NJ_ordered.iloc[index]["Product 2"] or "DGQHR" in NJ_ordered.iloc[index]["Product 2"]:
                        dif_targets = int(str(NJ_ordered.iloc[index]["Position 1"]).replace(" ","").replace("=", "")) - int(str(NJ_ordered.iloc[index-1]["Position 1"]).replace(" ","").replace("=", ""))
                        if  dif_targets < 15 and dif_targets > 1:
                            row[i] = (NJ_ordered.iloc[index][i] + NJ_ordered.iloc[index-1][i])/2
                            if row[i] != 0 and coverage_NJ >= 10: # Filter by frequency
                                    parsed_NJ_ordered.append(row)
                        elif row[i] != 0 and coverage_NJ >= 10:
                            parsed_NJ_ordered.append(row) # Add NJ similarly annotated but far away in the same gene
                    elif row[i] != 0 and coverage_NJ >= 10:
                        parsed_NJ_ordered.append(row)
        
        # Convert to dataframe ordered and already filtered
        parsed_NJ_ordered_df = pd.DataFrame(parsed_NJ_ordered)
        parsed_NJ_ordered_df_def = parsed_NJ_ordered_df.drop_duplicates(keep = "first")
        
        # Work on the filtered and ordered df to sum frequencies affecting the same genes
        genes_NJ_list = parsed_NJ_ordered_df_def["Gene 1"].tolist()
        multitargeted_by_NJ_genes = []
        for gene in genes_NJ_list: # keep genes targeted more than once that passed the filters
            if genes_NJ_list.count(gene) >= 2:
                multitargeted_by_NJ_genes.append(gene)
        multitargeted_by_NJ_genes = list(set(multitargeted_by_NJ_genes))
        #print(multitargeted_by_NJ_genes)
        # Get header for parallel NJ dataframe
        header_parallel_NJ = parsed_NJ_ordered_df_def.columns.values.tolist()
        header_parallel_NJ[0] = "Seq ID"
        header_parallel_NJ.remove("Position 1")
        header_parallel_NJ.remove("Reads (Cov) 1")
        header_parallel_NJ.remove("Seq ID 2")
        header_parallel_NJ.remove("Position 2")
        header_parallel_NJ.remove("Reads (Cov) 2")
        header_parallel_NJ.remove("Reads (Cov)" )
        header_parallel_NJ.remove("Score")
        header_parallel_NJ.remove("Skew")
        header_parallel_NJ.remove("Annotation 1")
        header_parallel_NJ.remove("Annotation 2")
        header_parallel_NJ.remove("Gene 2")
        header_parallel_NJ.remove("Product 2")
        header_parallel_NJ.remove("Length 1")
        header_parallel_NJ.remove("Length 2")
        per_gene_NJ.append(header_parallel_NJ)
        
        # Get a row for each gene targeted by NJ at position 1
        for gene in multitargeted_by_NJ_genes:
            parallel_row = []
            gene_df = parsed_NJ_ordered_df_def[parsed_NJ_ordered_df_def["Gene 1"] == gene]
            parallel_row.append(gene_df.iloc[0]["Seq ID 1"])
            parallel_row.append(gene_df.iloc[0]["Gene 1"])
            parallel_row.append(gene_df.iloc[0]["Product 1"])
            for i in range(15,26):
                repl_name = parsed_NJ_ordered_df_def.columns.values.tolist()[i]
                parallel_row.append(gene_df[repl_name].sum())
            per_gene_NJ.append(parallel_row)
            
        per_gene_NJ_df = pd.DataFrame(per_gene_NJ)
        per_gene_NJ_to_parse = per_gene_NJ_df.rename(columns=per_gene_NJ_df.iloc[0]).drop(per_gene_NJ_df.index[0])
        
        # Do the same for genes called as the second sequence of the NJ
        # Work on the filtered and ordered df to sum frequencies affecting the same genes
        genes_NJ_list_2 = parsed_NJ_ordered_df_def["Gene 2"].tolist()
        multitargeted_by_NJ_genes_2 = []
        for gene in genes_NJ_list_2: # keep genes targeted more than once that passed the filters
            if genes_NJ_list_2.count(gene) >= 2:
                multitargeted_by_NJ_genes_2.append(gene)
        multitargeted_by_NJ_genes_2 = list(set(multitargeted_by_NJ_genes_2))
        #print(multitargeted_by_NJ_genes_2)
        # Get header for parallel NJ dataframe
        header_parallel_NJ_2 = parsed_NJ_ordered_df_def.columns.values.tolist()
        header_parallel_NJ_2[0] = "Seq ID"
        header_parallel_NJ_2.remove("Position 1")
        header_parallel_NJ_2.remove("Reads (Cov) 1")
        header_parallel_NJ_2.remove("Seq ID 2")
        header_parallel_NJ_2.remove("Position 2")
        header_parallel_NJ_2.remove("Reads (Cov) 2")
        header_parallel_NJ_2.remove("Reads (Cov)" )
        header_parallel_NJ_2.remove("Score")
        header_parallel_NJ_2.remove("Skew")
        header_parallel_NJ_2.remove("Annotation 1")
        header_parallel_NJ_2.remove("Annotation 2")
        header_parallel_NJ_2.remove("Gene 1")
        header_parallel_NJ_2.remove("Product 1")
        header_parallel_NJ_2.remove("Length 1")
        header_parallel_NJ_2.remove("Length 2")
        per_gene_NJ_2.append(header_parallel_NJ_2)
        
        # Get a row for each gene targeted by NJ at position 2
        for gene in multitargeted_by_NJ_genes_2:
            parallel_row = []
            gene_df = parsed_NJ_ordered_df_def[parsed_NJ_ordered_df_def["Gene 2"] == gene]
            parallel_row.append(gene_df.iloc[0]["Seq ID 2"])
            parallel_row.append(gene_df.iloc[0]["Gene 2"])
            parallel_row.append(gene_df.iloc[0]["Product 2"])
            for i in range(15,26):
                repl_name = parsed_NJ_ordered_df_def.columns.values.tolist()[i]
                parallel_row.append(gene_df[repl_name].sum())
            per_gene_NJ_2.append(parallel_row)
        
        per_gene_NJ_df_2 = pd.DataFrame(per_gene_NJ_2)
        per_gene_NJ_to_parse_2 = per_gene_NJ_df_2.rename(columns=per_gene_NJ_df_2.iloc[0]).drop(per_gene_NJ_df_2.index[0])

#print(anc_muts)
#print(genes_per_rep)
#pprint.pprint(interg_pos_dic)
#pprint.pprint(pos_mut)
#pprint.pprint(pos_overmut)
#print(muts_ordered)
#print(parsed_ordered_muts)

parsed_muts_df = pd.DataFrame(parsed_muts)
parsed_muts_df_def = parsed_muts_df.drop_duplicates(keep = "first")

#print(parsed_muts_df_def)

#parsed_ordered_muts_df = pd.DataFrame(parsed_ordered_muts)
#parsed_ordered_muts_df_def = parsed_ordered_muts_df.drop_duplicates(keep = "first")

#print(parsed_ordered_muts_df_def)
#print(parsed_ordered_muts_df_def[parsed_ordered_muts_df_def["Gene"] == "rpoS→"])

#print(header_parallel)
#print(multitargeted_genes)

#per_gene_muts_df = pd.DataFrame(per_gene_muts)
#per_gene_muts_df_def = per_gene_muts_df.rename(columns=per_gene_muts_df.iloc[0]).drop(per_gene_muts_df.index[0])


#print(per_gene_muts_to_parse)
#print(parallel_clean_df)

parallel_parsed_df = pd.DataFrame(parallel_parsed)
parallel_parsed_df_def = parallel_parsed_df.drop_duplicates(keep = "first")

compensatory_candidates_df = pd.DataFrame(compensatory_candidates)
compensatory_candidates_df_def = compensatory_candidates_df.drop_duplicates(keep = "first")

#print(compensatory_candidates_df_def)
#pprint.pprint(anc_nj)

#print(summary_NJ_df_def)

#print(parsed_NJ_ordered_df)
#print(per_gene_NJ_to_parse)

####################################################### EXPORT TO EXCEL OUTPUT #######################################################

outname = DirName+"parsed_"+StrainName+".xlsx"

writer = pd.ExcelWriter(outname, engine='xlsxwriter')
contig_info.to_excel(writer, sheet_name='contig info', index=False)
parsed_muts_df_def.to_excel(writer, sheet_name = "parsed mutations", index = False)
parallel_parsed_df_def.to_excel(writer, sheet_name = "parallel mutations", index = False)
compensatory_candidates_df_def.to_excel(writer, sheet_name = "compensatory candidates", index = False)
parsed_NJ_df_def.to_excel(writer, sheet_name = "parsed NJ", index = False)
summary_NJ_df_def.to_excel(writer, sheet_name = "summary NJ", index = False)
per_gene_NJ_to_parse.to_excel(writer, sheet_name = "NJ per gene 1", index = False)
per_gene_NJ_to_parse_2.to_excel(writer, sheet_name = "NJ per gene 2", index = False)
writer.save()

print("Output saved to "+outname)
