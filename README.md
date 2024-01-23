# Plasmid-encoded insertion sequences promote rapid adaptation in clinical enterobacteria

This repository contains all the code used in the work **Plasmid-encoded insertion sequences promote rapid adaptation in clinical enterobacteria**. The code developed for each of the analysis is organized in order for the sections in the results part of the main text of the paper, following the order in which the data was analyzed. For more detail regarding questions different to the code for analyzing and/or plotting consider consulting the Supplementary Material/Methods section of the article. The version of each software is specified in the Methods section of the article.

## Phenotypic and phylogenetic analyses

#### Phylogenetic tree

For the analysis of the phylogeny of the strains included in our work we used mashtree as shown in [/phylogeny_and_phenotype/phylogeny_script.sh](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/phylogeny_and_phenotype/phylogeny_script.sh) using the closed reference genomes of each strain uploaded at *link_to_SRA_sequences*. We represented the phylogenetic tree as shown in **Figure 1A** using the online tool of iTOL at https://itol.embl.de/. 

#### pOXA-48 fitness effects

To study the fitness effects of pOXA-48 in the bacterial hosts chosen for the experimental evolution assay, we measured the OD of 6 replicates with and without pOXA-48 for each strain during 24 hours (/phylogeny_and_phenotype/curves_DFE) and calculated the growth of pOXA-48-carrying samples relative to pOXA-48-free ones. The code for the analysis, statistics and plotting of this data (**Figure 1B**) can be found at [/phylogeny_and_phenotype/curves_DFE.R](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/phylogeny_and_phenotype/curves_DFE.R).

#### Bacterial growth during experimental evolution

To assess the phenotipic changes during the experimental evolution assay, we measured the bacterial growth of each replicate population each 2 days of the experiment. The raw data for this analysis can be found at /phylogeny_and_phenotype/Curvas and the script for the plotting (**Figure 1D**) and statistical analyses can be checked at [/phylogeny_and_phenotype/curves_non_cont.R](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/phylogeny_and_phenotype/curves_non_cont.R).

## Genomic analyses of the Experimental Evolution

#### Reference genomes

To perform the downstream analysis steps we firstly assembled the reference genomes by hybrid assembly using Unicycler and combining short and long reads. In the cases in which the genomes couldn't be assembled using Unicycler, we employed Flye and polished the assemblies using Medaka and several rounds of Pilon. The final assembled closed genomes were annotated using the NCBI Prokaryotic Genome Annotation Pipeline (PGAP). All these steps are summarized in [/WGS_analysis_EE/assembly_and_annotation.sh](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/assembly_and_annotation.sh).

#### Data processing and variant calling

The first step in our genomic analyses were the usual quality control and trimming of the raw reads from the sequencing. With the clean reads, we performed the variant calling of the samples against their reference genome using breseq. All these steps and the commands used are summarized in /WGS_analysis_EE/breseq_pipeline.sh. The variant calling results reported by breseq as an html table for each replicate were merged into a xlsx table per strain using [/WGS_analysis_EE/breseq_index_parser.py](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/breseq_index_parser.py) (modified from [this script](https://github.com/sirmicrobe/LabScripts/blob/master/BreseqCat3.py)). The dependencies needed for its execution can be found at /WGS_analysis_EE/environment_pybreseqparser.yml. Then, we performed the filtering of the variant calling results as described in the Methods section using the code in [/WGS_analysis_EE/parser_EE.py](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/parser_EE.py) (for populations) and [/WGS_analysis_EE/parser_clones_EE.py](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/parser_clones_EE.py) (for clones) and summarized each table using [/WGS_analysis_EE/parser_vcall_summary.py](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/parser_vcall_summary.py).Finally we merged the data for all the strains into a unique table to ease the downstream analyses. To show the results as in **Figure 2**, we summarized the events depending on their type and plotted using the first part of the code in [/WGS_analysis_EE/plots_summary_parsed_vcall.R](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/plots_summary_parsed_vcall.R).

We also summarized the data per strain (**Supplementary Figures 7-19**) and represented the location of the mutations in the different replicons of each genome (**Figure 3AB**) using the circlize package in the code developed in the second part of the code in [/WGS_analysis_EE/plots_summary_parsed_vcall.R](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/plots_summary_parsed_vcall.R). We also intended to show different pathways of adaptation at the main targets found during our experiment, for which we used the lolliplot function from the TrackViewer package as shown in [/WGS_analysis_EE/plots_lolliplots_operons.R](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/plots_lolliplots_operons.R), to represent **Figure 3CD** and lolliplots of the supplementary panels per strain.

#### Tracking IS1 elemments

To check whether we could confirm the transposition of IS1 elements from pOXA-48 into other regions of the genomes, we studied if there were more identical copies of these in other regions of the genome. To do so, we performed the steps specified in [/WGS_analysis_EE/tracking_IS1.sh](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/tracking_IS1.sh) and determined those strains in which we could track the transposition (i.e. those in which the BLAST hits only showed 100% of identity and length for the pOXA-48 copies).

#### Long-read analysis

To confirm the structural variants detected by short read data as analyzed in the previous sections, we decided to sequence one evolved clone for each _K. pneumoniae_ population which showed high frequency NJ events of interest, as well as one evolved clone for each _E. coli_ population as control. The code for the variant calling analysis of this data is summarized in [/WGS_analysis_EE/long_read_analysis.sh](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/long_read_analysis_EE.sh). The code to performe closed hybrid assemblies of these samples can also be found at [/WGS_analysis_EE/long_read_analysis.sh](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/long_read_analysis_EE.sh).

## Analysis of IS1 transposition depending on genomic IS1 copy number

To test the hypothesis of the genomic IS1 copy number regulating the transposition rate in the strains under study, we analyzed multiple sources of data. We first analyzed the EE variant calling results as well as transcriptomic data. Then, we also complemented our observations with a fluctuation assay testing transposition rate by checking non capsulated mutants. All these analyses shown in **Figure 4** are depicted in the following subsections.

#### Correlation of IS1 transposition in EE & genomic copy number

The first observation we could make about the IS1 copy number regulating the transposition during evolution was the correlation between the increment in IS1 transposition in pOXA-48-carrying populations (w & w/o AMC) relative to pOXA-48-free populations with the genomic IS1 copy number for each of the strains (**Supplementary Figure 21A**). The code for this analysis is found in [/IS1_transposition_copy_number/IStrEE_nIS_correlation.R](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/IS1_transposition_copy_number/IStrEE_nIS_correlation.R).

#### Transcriptomic analyses to support the correlation

To further support the results we got in the previous section, we analyzed the change in the expression of IS1 elements with the entrance of pOXA-48 in a subset of clinical strains. Then we correlated these results with the genomic IS1 copy number in the strains (**Supplementary Figure 21B**). The code for these analyses can be found in [/IS1_transposition_copy_number/RNASeq_correlation.R](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/IS1_transposition_copy_number/RNASeq_correlation.R).

#### Fluctuation assay testing transposition rate

In order to analyze the increment in mutants of interest (non capsulated) in the different conditions of this assay, we used the RSalvador package as determined in [/IS1_transposition_copy_number/RSalvador_script.R](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/IS1_transposition_copy_number/RSalvador_script.R). We checked the mutation profile in the non capsulated mutants for all the conditions by sequencing them and analyzing as described in the **variant calling** subsection of the **Genomic analyses during the EE** section. All these results regarding the fluctuation assay of non capsulated mutants testing IS1 transposition rate can be found in **Figure 4D-E** and **Supplementary Figure 20**.

## Genomic analyses of _in vivo_ evolution

In order to try to extrapolate our EE results beyond the conditions tested, we intended to analyze within-patient evolution of pOXA-48 carrying isolates. The code employed for this, which results can be found in **Figure 5** is specified in the following subsections.

#### Variant calling of _in vivo_ evolved clones

To analyze the mutations present in _in vivo_ evolved clones we performed the same analysis as the one explained for the EE variant calling in /WGS_analysis_EE/breseq_pipeline.sh. We closed the genomes of one of the isolates per lineage doing hybrid assembly and annotated them using PGAP as previously explained in [/WGS_analysis_EE/assembly_and_annotation.sh](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/WGS_analysis_EE/assembly_and_annotation.sh).  

#### GAM model of IS1 transposition _in vivo_

To test the main forces driving the increase in IS1 transposition (species or genomic IS1 copy number), we built a Generalized Additive Model (GAM) using the mgcv package in R. The code for the development of the model can be found in /WGS_analysis_in_vivo/GAM_model_patients.R. The results of this analysis can be found at **Supplementary Figure 22**.

## Additional analyses

In addition to the analyses shown in the main text, we performed complementary analyses to check the sequence of pOXA-48 in each of the strains, calculate the plasmid copy number in each of the samples, calculate the genomic IS1 copy number in the genome of each strain, etc.

#### Calculation of plasmid copy number (PCN) of pOXA-48 and all the other plasmids for each strain

To correct the number of genomic IS1 copies in each strain, we needed to take into account the possible elements encoded in plasmids. Then, the copy number of these would influence the genomic IS1 copy number. The code employed for this can be found in [/additional_analyses/PCN_analysis.sh](https://github.com/jorgEVOplasmids/rapid_adaptation_pOXA48/blob/main/additional_analyses/PCN_analysis.sh).

#### Calculation of genomic IS1 copy number

To determine the genomic IS1 copy number of the strains we used the annotation GBK files from PGAP output for each of the strains (both in the EE and in the _in vivo_ analyses). We considered those elements encoded in plasmids and multiplied them by the corresponding PCN.

#### Checking pOXA-48 sequence

To analyze whether the strains contained the most common variant of pOXA-48 (K8) we performed the variant calling of the strains at day 1 of the EE against a common reference of pOXA-48 K8. The code for these were the same as for the variant calling explained for the EE.
