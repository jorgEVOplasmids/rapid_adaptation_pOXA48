#! /usr/bin/python3
# -*- coding: utf-8 -*-

from bs4 import BeautifulSoup
import openpyxl
import pandas as pd
import os, argparse, sys
from functools import reduce


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''=========================================== breseq index parser ============================================

This parser filters the index.html files created by breseq in either consensus or polymorphism mode using
the same reference sequence. It generates an excel table summarizing the information of all given breseq
output directories, including the following table sheets:
- Contig information
- Predicted mutations: comparison of presence/absence (consensus) or frequency (-p) of mutations.
- Mutations info: comparison of alternative allele depth and total read depth of all mutations.
- Unassigned missing coverage evidence: comparison of presence/absence.
- Full MC table with all relevant information.
- Unassigned new junction evidence: comparison of presence/absence.
- Full JC table with all relevant information.

===========================================================================================================''')

parser.add_argument('DIRECTORY', action="store", help="Master directory containing all breseq output directories to parse")
parser.add_argument('-p', action="store_true", help="Include this flag if you ran breseq in polymorphism mode", default=False)

args = parser.parse_args()

DirName = args.DIRECTORY
DirValidation = os.path.isdir(DirName)
if DirValidation not in ['True', True]:
	print("Error. Please, enter a valid master directory to parse. Try --help to see the parser's description. Exiting!")
	exit(1)
DirectoryNames = os.listdir(DirName)
DirectoryNames.sort()

mutListdf = []
mcfullListdf = []
mcListdf = []
jcfullListdf = []
jcListdf = []

print("== breseq index parser ==")


### ITERATING OVER THE INDEX.HTML FILES IN THE MASTER DIRECTORY

for i, element in enumerate(DirectoryNames):
	SubDirName=DirectoryNames[i]
	SubDirName=str(SubDirName)
	print("Parsing "+SubDirName)
	file = DirName+'/'+SubDirName+'/output/index.html'
	frames = []
	
	# Column names for the dataframes
	clonal_col_names=['Evidence', 'Seq ID', 'Position', 'Mutation', 'Annotation', 'Gene', 'Description', SubDirName]
	poly_col_names=['Evidence', 'Seq ID', 'Position', 'Mutation', 'Annotation', 'Gene', 'Description', SubDirName]
	mcfull_names=['Seq ID', 'Start', 'End', 'Size', '<-Reads', 'Reads->', 'Gene', 'Description', SubDirName]
	mc_names=['Seq ID', 'Gene', 'Description', SubDirName]
	njefull_names=['Seq ID 1', 'Position 1', 'Reads (Cov) 1', 'Seq ID 2', 'Position 2', 'Reads (Cov) 2', 'Reads (Cov)', 'Score', 'Skew', 'Annotation 1', 'Gene 1', 'Product 1', 'Annotation 2', 'Gene 2', 'Product 2', SubDirName]
	nje_names=['Seq ID 1', 'Seq ID 2', 'Annotation 1', 'Gene 1', 'Product 1', 'Annotation 2', 'Gene 2', 'Product 2', SubDirName]
	
	# Initializing daframes with column names	
	mcfull_df = pd.DataFrame(columns = mcfull_names)
	njefull_df = pd.DataFrame(columns = njefull_names)
	mc_df = pd.DataFrame(columns = mc_names)
	nje_df = pd.DataFrame(columns = nje_names)
	
	if args.p in ["True", True]:
		snp_df = pd.DataFrame(columns = poly_col_names)
		poly_snp_df = pd.DataFrame(columns = poly_col_names)
	
	if args.p in ["False", False]:
		snp_df = pd.DataFrame(columns = clonal_col_names)
		poly_snp_df = pd.DataFrame(columns = clonal_col_names)
	
	
	if os.path.isfile(file) == True:
		try:
			with open(file, 'rb') as file:
				
				
				### CREATING THE MUTATIONS DATAFRAMES
				
				soup = BeautifulSoup(file, 'lxml')
				normal_table = soup.find_all(attrs={'class':'normal_table_row'})
				poly_table = soup.find_all(attrs={'class':'polymorphism_table_row'})
				
				# Appending rows of mutations with 100% frequency
				for rows in normal_table:
					line = [cel.get_text(strip=True) for cel in rows.find_all('td')]
					line = [el.replace('\xa0',' ') for el in line]
					line = [el.replace(' genes',' genes: ') for el in line]
					line = [el.replace('%','') for el in line]
										
					if args.p in ["True", True]:
						order = [0,1,2,3,5,6,7,4]  # Reorder to get freq below strain name
						line = [line[i] for i in order]
						snp_df = snp_df.append(pd.Series(line, index=poly_col_names), ignore_index=True)
					if args.p in ["False", False]:
						line.append("1")  # Presence of mutation
						snp_df=snp_df = snp_df.append(pd.Series(line, index=clonal_col_names), ignore_index=True)
				
				# Appending rows of mutations in polymorphism mode				
				for poly_rows in poly_table:
					poly_line = [cel.get_text(strip=True) for cel in poly_rows.find_all('td')]
					poly_line = [el.replace('\xa0', ' ') for el in poly_line]
					poly_line = [el.replace(' genes ', ' genes: ') for el in poly_line]
					poly_line = [el.replace('%', '') for el in poly_line]
										
					if args.p in ["True", True]:
						order = [0,1,2,3,5,6,7,4]  # Reorder to get freq below strain name
						poly_line = [poly_line[i] for i in order]
						poly_snp_df = poly_snp_df.append(pd.Series(poly_line, index=poly_col_names), ignore_index=True)
					if args.p in ["False", False]:
						poly_snp_df= poly_snp_df.append(pd.Series(poly_line, index=clonal_col_names), ignore_index=True)
				
				
				### CREATING DE MC UNASSIGNED EVIDENCE DATAFRAMES
				
				# Delimiting the MC table
				alph_soup = str(soup)
				begin_umc = alph_soup.find('<tr><th align="left" class="missing_coverage_header_row" colspan="11">Unassigned missing coverage evidence</th></tr>')
				end_umc = alph_soup.find('<th align="left" class="new_junction_header_row" colspan="12">Unassigned new junction evidence</th>')
				chopped_soup = alph_soup[begin_umc:end_umc]
				soup_mc = BeautifulSoup(chopped_soup, 'lxml')
				
				if soup_mc:  # If there are MC items
					for rows in soup_mc:
						line = [cel.get_text(strip=True) for cel in rows.find_all('td')]
						line = [el.replace('\xa0', ' ') for el in line]
						line = [line[x:x+11] for x in range(0, len(line), 11)]  # Break the list in chunks of 11 size to get a list of lists (rows)
						
						for el in line:
							el = el[3:]  # Remove * and ? fields
							el.append("1")
							# Adding rows for the full MC table
							mcfull_df = mcfull_df.append(pd.Series(el, index=mcfull_names), ignore_index=True)
							# Adding rows for the comparison MC table
							del el[1:6]
							mc_df = mc_df.append(pd.Series(el, index=mc_names), ignore_index=True)
				
				
				### CREATING DE JC UNASSIGNED EVIDENCE DATAFRAMES
				
				# Finding the JC table
				begin_nje = alph_soup.find('<!-- Side 1 Item Lines for New Junction -->')
				chopped_soup2 = alph_soup[begin_nje:]
				soup_nje = BeautifulSoup(chopped_soup2, 'lxml')
				
				if soup_nje:  # If there are JC items
					# Read rows into a list of lists. Both sides of a JC are consecutive in the list
					tab_data = [[celldata.get_text(strip=True) for celldata in rowdata.find_all(['td'])] for rowdata in soup_nje.find_all('tr')]
					tab_data = [[el.replace('\xa0', ' ') for el in el1] for el1 in tab_data]
					tab_data = [[el.replace('%', '') for el in el1] for el1 in tab_data]
					tab_data = [[x for x in el1 if "*" not in x and "?" not in x] for el1 in tab_data]
					new_tab_data = []
					
					# Concatenate boths sides of the JC into a single list (inside the list) as pairs				
					for i in range(0, len(tab_data),2):
						nj_row = tab_data[i]+tab_data[i+1]
						order = [0,1,2,10,11,12,3,4,5,7,8,9,13,14,15,6]  # Reorder the row; freq below strain name
						nj_row = [nj_row[i] for i in order]
						# Adding rows for the full JC table
						njefull_df = njefull_df.append(pd.Series(nj_row, index=njefull_names), ignore_index=True)
						# Adding rows for the comparison JC table
						del nj_row[1:3]
						del nj_row[2:7]
						nje_df = nje_df.append(pd.Series(nj_row, index=nje_names), ignore_index=True)
			
			
			### SAVING THE DATAFRAMES IN A LIST OF DATAFRAMES
			
			# Mutations and full MC and JC dataframes
			frames = frames+[snp_df, poly_snp_df]
			all_snps = pd.concat(frames)
			mutListdf.append(all_snps)
			mcfullListdf.append(mcfull_df)
			jcfullListdf.append(njefull_df)
			# Simplified MC and JC dataframes for comparison tables. Sometimes a region can have multiple MC or JC afecting the same genes;
			# to make the comparison table only the Seq ID and annotations are used, which can produce duplicate rows.
			# Duplicate rows are removed to construct the comparison table (note that it is a presence/absence table with "1" for "present").
			mc_df.drop_duplicates(inplace=True)
			mcListdf.append(mc_df)
			nje_df[SubDirName] = "1"
			nje_df.drop_duplicates(inplace=True)
			jcListdf.append(nje_df)
		
		
		except Exception as e:
			print("Could not parse directory "+SubDirName+ ". Did you correctly call the -p flag? If you didn't check the python script. Error description:")
			print(e)
	
	else:
		print("The index.html file for " +SubDirName + " could not be found")


### ITERATING OVER THE FASTA.FAI IN THE MASTER DIRECTORY

fai_list = []

for i, element in enumerate(DirectoryNames):
	SubDirName=DirectoryNames[i]
	SubDirName=str(SubDirName)
	fai = DirName+'/'+SubDirName+'/data/reference.fasta.fai'
	
	if os.path.isfile(fai) == True:
		try:
			with open(fai, 'r') as fai:
				for line in fai:
					line = line.rstrip()
					line = line.split('\t')
					line = line[0:2]
					if line not in fai_list:
						fai_list.append(line)
		
		except Exception as e:
			print("Could not parse the referece.fasta.fai file in "+SubDirName+ ". Error description:")
			print(e)
	else:
		print("The reference.fasta.fai file for " +SubDirName + " could not be found. You can create it with samtools faidx.")

fai_df = pd.DataFrame(fai_list, columns = ['Seq ID', 'Length'])


### MERGING THE DATAFRAMES

mutdf_merged = reduce(lambda left,right: pd.merge(left,right,on=['Evidence', 'Seq ID', 'Position', 'Mutation', 'Annotation', 'Gene', 'Description'], how='outer'), mutListdf).fillna('0')
mutdf_len_merged = pd.merge(mutdf_merged, fai_df, on='Seq ID')

mcfulldf_merged = reduce(lambda left,right: pd.merge(left,right,on=['Seq ID', 'Start', 'End', 'Size', '<-Reads', 'Reads->', 'Gene', 'Description'], how='outer'), mcfullListdf).fillna('0')
mcfulldf_len_merged = pd.merge(mcfulldf_merged, fai_df, on='Seq ID')

mcdf_merged = reduce(lambda left,right: pd.merge(left,right,on=['Seq ID', 'Gene', 'Description'], how='outer'), mcListdf).fillna('0')
mcdf_len_merged = pd.merge(mcdf_merged, fai_df, on='Seq ID')

jcfulldf_merged = reduce(lambda left,right: pd.merge(left,right,on=['Seq ID 1', 'Position 1', 'Reads (Cov) 1', 'Seq ID 2', 'Position 2', 'Reads (Cov) 2', 'Reads (Cov)', 'Score', 'Skew', 'Annotation 1', 'Gene 1', 'Product 1', 'Annotation 2', 'Gene 2', 'Product 2'], how='outer'), jcfullListdf).fillna('0')
fai_df.columns = ['Seq ID 1', 'Length 1']
jcfulldf_len_merged = pd.merge(jcfulldf_merged, fai_df, on=['Seq ID 1'], how='left')
fai_df.columns = ['Seq ID 2', 'Length 2']
jcfulldf_len_merged = pd.merge(jcfulldf_len_merged, fai_df, on=['Seq ID 2'], how='left')
fai_df.columns = ['Seq ID', 'Length']

jcdf_merged = reduce(lambda left,right: pd.merge(left,right,on=['Seq ID 1', 'Seq ID 2', 'Annotation 1', 'Gene 1', 'Product 1', 'Annotation 2', 'Gene 2', 'Product 2'], how='outer'), jcListdf).fillna('0')
fai_df.columns = ['Seq ID 1', 'Length 1']
jcdf_len_merged = pd.merge(jcdf_merged, fai_df, on=['Seq ID 1'], how='left')
fai_df.columns = ['Seq ID 2', 'Length 2']
jcdf_len_merged = pd.merge(jcdf_len_merged, fai_df, on=['Seq ID 2'], how='left')
fai_df.columns = ['Seq ID', 'Length']


### ITERATING OVER THE VCF FILE IN THE MASTER DIRECTORY

vcfListdf = []
for i, element in enumerate(DirectoryNames):
	SubDirName=DirectoryNames[i]
	SubDirName=str(SubDirName)
	vcf = DirName+'/'+SubDirName+'/data/output.vcf'
	vcf_list = []
	
	if os.path.isfile(vcf) == True:
		try:
			with open(vcf, 'r') as vcf:
				for line in vcf:
					if not line.startswith('#'):
						line = line.rstrip()
						line = line.split('\t')
						if "AD" in line[7]:
							info = line[7].split(';')
							AD = info[1].split('=')
							AD = AD[1]
							DP = info[2].split('=')
							DP = DP[1]
							newline = [line[0], line[1], line[3], line[4], AD, DP]
							vcf_list.append(newline)
		
		except Exception as e:
			print("Could not parse the VCF file in "+SubDirName+ ". Error description:")
			print(e)
	else:
		print("The VCF file for " +SubDirName + " could not be found")
	
	vcf_df = pd.DataFrame(vcf_list, columns = ['Seq ID', 'Position', 'Ref', 'Alt', 'Allele Depth '+SubDirName, 'Total Depth '+SubDirName])
	vcfListdf.append(vcf_df)

vcfdf_merged = reduce(lambda left,right: pd.merge(left,right,on=['Seq ID', 'Position', 'Ref', 'Alt'], how='outer'), vcfListdf).fillna('0')


### EXPORTING TO EXCEL

newDirName = DirName.split("/")
if DirName[-1] == "/":
	newDirName = newDirName[-2]
else:
	newDirName = newDirName[-1]
suff = ""
if args.p:
	suff = "_p"
outname = "summary_table_"+newDirName+suff+".xlsx"

writer = pd.ExcelWriter(outname, engine='xlsxwriter')
fai_df.to_excel(writer, sheet_name='contig info', index=False)
mutdf_len_merged.to_excel(writer, sheet_name='mutations', index=False)
vcfdf_merged.to_excel(writer, sheet_name='mutations info', index=False)
mcdf_len_merged.to_excel(writer, sheet_name='missing coverage', index=False)
mcfulldf_len_merged.to_excel(writer, sheet_name='MC full', index=False)
jcdf_len_merged.to_excel(writer, sheet_name='new junctions', index=False)
jcfulldf_len_merged.to_excel(writer, sheet_name='JC full', index=False)
writer.save()

print("Done! Output can be found in ./"+outname)

