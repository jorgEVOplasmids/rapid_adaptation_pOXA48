#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import pandas as pd
import pprint
import openpyxl

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""
Created on Thu Oct 20 10:22:49 2022

@author: Jorge Sastre Domínguez

This script takes the output of breseq_index_parser.py in clonal mode and reorganizes the New Junctions sheet to remove duplicated called junctions
to report them at gene level. It also takes the raw output of breseq_index_parser.py to take into account mutations present in ancestral populations
and filters them out.

""")

parser.add_argument('Directory', action="store",
                    help="Master directory containing parsed xlsx output by breseq_index_parser.py")
parser.add_argument('Directory_ancestor', action="store",
                    help="Master directory containing parsed xlsx output by breseq_index_parser.py for populations containing ancestral samples")
parser.add_argument('Strain', action="store",
                    help="Strain name of the xlsx file to apply the filters")

args = parser.parse_args()

DirName = args.Directory
DirAncName = args.Directory_ancestor
DirValidation = os.path.isdir(DirName)
DirAncValidation = os.path.isdir(DirAncName)
StrainName = args.Strain

if DirValidation not in ['True', True]:
    print("Error. Please, enter a valid master directory to parse. Try --help to see the parser's description. Exiting!")
    exit(1)
    
if DirAncValidation not in ['True', True]:
    print("Error. Please, enter a valid master directory to parse. Try --help to see the parser's description. Exiting!")
    exit(1)

FileNames = os.listdir(DirName)
FileAncNames = os.listdir(DirAncName)

# Initialize empty variables that will be converted in the parsed data frames later

anc_muts = {}
anc_nj = {}

parsed_muts = []
parsed_ordered_muts = []

range_clones_muts = []
range_clones_NJs = []
parsed_NJ = []
summary_NJ_df = []
parsed_NJ_ordered = []
per_gene_NJ = []
per_gene_NJ_2 = []

######################################################### GET MUTATIONS PRESENT IN ANCESTORS ##############################################
# Open breseq raw xlsx file for getting ancestors mutations of the specified strain

print("Getting filters for parsing mutations and NJ sheet of "+StrainName)

for i, element in enumerate(FileAncNames):
    if StrainName in FileAncNames[i]:
        DirAncName=str(DirAncName)
        print("Parsing "+DirAncName+FileAncNames[i])
        xlsx_anc = pd.ExcelFile(DirAncName+FileAncNames[i])
        contig_info_anc = pd.read_excel(xlsx_anc, "contig info")
        muts_anc = pd.read_excel(xlsx_anc, "mutations")
        NJ_anc = pd.read_excel(xlsx_anc, "JC full")
        MC_anc = pd.read_excel(xlsx_anc, "MC full")
        
        # Get positions of SNPs present in ancestors to filter them out from the clones
        
        for index, row in muts_anc.iterrows():
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
                    
        # Get positions in which the ancestors show NJ
        
        for index, row in NJ_anc.iterrows():
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

######################################################### PARSE CLONES XLSX FILES ##############################################
# Open breseq xlsx input file for the specified strain clones

for i, element in enumerate(FileNames):
    if StrainName in FileNames[i]:
        DirName=str(DirName)
        print("Parsing "+DirName+FileNames[i])
        xlsx = pd.ExcelFile(DirName+FileNames[i])
        contig_info = pd.read_excel(xlsx, "contig info")
        muts = pd.read_excel(xlsx, "mutations")
        NJ = pd.read_excel(xlsx, "JC full")
        MC = pd.read_excel(xlsx, "MC full")
        #print(len(muts))
        #print(len(NJ))
        
        print("Filtering SNPs present in ancestral populations")
        
        # Get columns of replicates of clones and the range to parse SNPs
        
        for column in muts.columns.values.tolist():
            #print(column)
            if StrainName in column:
                #print(NJ.columns.get_loc(column))
                range_clones_muts.append(muts.columns.get_loc(column))
        
        for index, row in muts.iterrows():
            if ":" not in row[2]:
                position = int(row[2].replace(",",""))
            else:
                position = int(row[2].split(":")[0].replace(",",""))
            for i in range(range_clones_muts[0],range_clones_muts[-1] + 1):
                repl_name = muts.columns.values.tolist()[i]
                # Declare empty list of ancestral mutations in contigs for avoiding key errors
                if row[1] not in anc_muts.keys():
                    anc_muts[row[1]] = []
                    
                # Filter out positions mutated in ancestors                    
                if position not in anc_muts[row[1]]:
                    parsed_muts.append(row)
                    
        print("Rearranging NJ sheet of "+StrainName+ " clones")
        
        # Get columns of replicates of clones and the range to parse NJs
        
        for column in NJ.columns.values.tolist():
            #print(column)
            if StrainName in column:
                #print(NJ.columns.get_loc(column))
                range_clones_NJs.append(NJ.columns.get_loc(column))
        
        #print(range_clones[0],range_clones[-1])
        
        # Get contig information to identify events into the chromosome, between plasmids and into plasmids
        
        plasmids_id = []
        for index, row in contig_info.iterrows():
            if row[1] > 4000000:
                chr_id = row[0]
            elif row[1] > 63000 and row[1] < 66000:
                poxa_id = row[0]
            else:
                plasmids_id.append(row[0])
                
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
            
            for i in range(range_clones_NJs[0],range_clones_NJs[-1]+1): # parse by replicate
                repl_name = NJ.columns.values.tolist()[i]
                if row[i] != 0 and position_1 != 1:
                    #print(repl_name, contig_1, contig_2)
                    if position_1 not in anc_nj[contig_1] and position_2 not in anc_nj[contig_2] and coverage_NJ >= 10: # filter out those NJ in ancestors and surrounding intergenic SNPs
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
            
            for i in range(range_clones_NJs[0],range_clones_NJs[-1]+1): # parse by replicate
                repl_name = NJ_ordered.columns.values.tolist()[i]
                if row[i] != 0 and position_1 not in anc_nj[contig_1]:
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
                if row[i] != 0 and position_2 not in anc_nj[contig_2]:
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
            for i in range(range_clones_NJs[0],range_clones_NJs[-1]+1):
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
            for i in range(range_clones_NJs[0],range_clones_NJs[-1]+1):
                repl_name = parsed_NJ_ordered_df_def.columns.values.tolist()[i]
                parallel_row.append(gene_df[repl_name].sum())
            per_gene_NJ_2.append(parallel_row)
        
        per_gene_NJ_df_2 = pd.DataFrame(per_gene_NJ_2)
        per_gene_NJ_to_parse_2 = per_gene_NJ_df_2.rename(columns=per_gene_NJ_df_2.iloc[0]).drop(per_gene_NJ_df_2.index[0])
        
        
#print(anc_muts)
#print(anc_nj)

parsed_muts_df = pd.DataFrame(parsed_muts)
parsed_muts_df_def = parsed_muts_df.drop_duplicates(keep = "first")

####################################################### EXPORT TO EXCEL OUTPUT #######################################################

outname = DirName+"parsed_"+StrainName+"_clones.xlsx"

writer = pd.ExcelWriter(outname, engine='xlsxwriter')
contig_info.to_excel(writer, sheet_name='contig info', index=False)
parsed_muts_df_def.to_excel(writer, sheet_name = "parsed mutations", index = False)
parsed_NJ_df_def.to_excel(writer, sheet_name = "parsed NJ", index = False)
per_gene_NJ_to_parse.to_excel(writer, sheet_name = "NJ per gene 1", index = False)
per_gene_NJ_to_parse_2.to_excel(writer, sheet_name = "NJ per gene 2", index = False)
writer.save()

print("Output saved to "+outname)