# Plasmid-encoded insertion sequences promote rapid adaptation in clinical enterobacteria

This repository contains all the code used in the work **Plasmid-encoded insertion sequences promote rapid adaptation in clinical enterobacteria**. The code developed for each of the analysis is organized in order for the sections in the results part of the main text of the paper, following the order in which the data was analyzed. For more detail regarding questions different to the code for analyzing and/or plotting consider consulting the Supplementary Material/Methods section of the article. The version of each software is specified in the Methods section of the article.

## Phenotypic and phylogenetic analyses

#### Phylogenetic tree

For the analysis of the phylogeny of the strains included in our work we used mashtree as shown in *link_to_script* using the closed reference genomes of each strain *link_OR_reference_to_sequences*. We represented the phylogenetic tree as shown in **Figure 1A** using the online tool of iTOL at https://itol.embl.de/. 

#### pOXA-48 fitness effects

To study the fitness effects of pOXA-48 in the bacterial hosts chosen for the experimental evolution assay, we measured the OD of 6 replicates with and without pOXA-48 for each strain during 24 hours (*link_to_table_of_data*) and calculated the growth of pOXA-48-carrying samples relative to pOXA-48-free ones. The code for the analysis, statistics and plotting of this data (**Figure 1B**) can be found at *link_to_script*.

#### Bacterial growth during experimental evolution

To assess the phenotipic changes during the experimental evolution assay, we measured the bacterial growth of each replicate population each 2 days of the experiment. The raw data for this analysis can be found at *link_to_table_of_data* and the script for the plotting (**Figure 1D**) and statistical analyses can be checked at *link_to_script*.

## Genomic analyses of the Experimental Evolution

#### Data processing and variant calling

The first step in our genomic analyses were the usual quality control and trimming of the raw reads from the sequencing. With the clean reads, we performed the variant calling of the samples against their reference genome using breseq. All these steps and the commands used are summarized in *link_to_script*. The variant calling results reported by breseq as an html table for each replicate were merged into a xlsx table per strain using *link_to_script* (modified from *ref_to_original_script_dev*). Then, we performed the filtering of the variant calling results as described in the Methods section using the code in *link_to_script* and merged the data for all the strains into a unique table to ease the downstream analyses with *link_to_script*. To show the results as in **Figure 2**, we summarized the events depending on their type and plotted using the first part of the code in *link_to_script*.

We also summarized the data per strain (**Supplementary Figures 7-19**) and represented the location of the mutations in the different replicons of each genome (**Figure 3AB**) using the circlize package in the code developed in *link_to_script*. We also intended to show different pathways of adaptation at the main targets found during our experiment, for which we used the lolliplot function from the TrackViewer package as shown in *link_to_script*, to represent **Figure 3CD**.

#### Tracking IS1 elemments

To check whether we could confirm the transposition of IS1 elements from pOXA-48 into other regions of the genomes, we studied if there were more identical copies of these in other regions of the genome. To do so, we performed the steps specified in *link_to_script* and determined those strains in which we could track the transposition (i.e. those in which the BLAST hits only showed 100% of identity and length for the pOXA-48 copies).

#### Long-read analysis

To confirm the structural variants detected by short read data as analyzed in the previous sections, we decided to sequence one evolved clone for each _K. pneumoniae_ population which showed high frequency NJ events of interest, as well as one evolved clone for each _E. coli_ population as control. The code for the variant calling analysis of this data is summarized in *link_to_script*. The code to performe closed hybrid assemblies of these samples can be found at *link_to_script*.

## Transcriptomic analyses

## Genomic analyses of _in vivo_ evolution

#### Variant calling of _in vivo_ evolved clones

#### GAM model of IS1 transposition _in vivo_

## Additional analyses

In addition to the analyses shown in the main text, we performed complementary analyses to check the sequence of pOXA-48 in each of the strains, calculate the plasmid copy number in each of the samples, etc.


