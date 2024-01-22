#!/bin/bash

### Extract the IS1 sequence of each pOXA-48

samtools faidx C021p.fasta contig_4:41726-42423 > C021_pOXA_IS1.fasta
samtools faidx C286p.fasta 5:32412-33109 > C286_pOXA_IS1.fasta
samtools faidx C309p.fasta contig_3:7871-8568 > C309_pOXA_IS1.fasta
samtools faidx C324p.fasta 4:32412-33109 > C324_pOXA_IS1.fasta
samtools faidx CF12.fasta 2:30502-31199 > CF12_pOXA_IS1.fasta
samtools faidx CF13.fasta 2:30502-31199 > CF13_pOXA_IS1.fasta
samtools faidx H53.fasta 2:30502-31199 > H53_pOXA_IS1.fasta
samtools faidx K091p.fasta 3:32412-33109 > K091_pOXA_IS1.fasta
samtools faidx K147.fasta 3:32412-33109 > K147_pOXA_IS1.fasta
samtools faidx K153_pOXA48.fasta contig_2:4580-5277 > K153_pOXA_IS1.fasta
samtools faidx K163.fasta 3:32412-33109 > K163_pOXA_IS1.fasta
samtools faidx K209p.fasta 4:32412-33109 > K209_pOXA_IS1.fasta
samtools faidx K25.fasta 3:32412-33109 > K25_pOXA_IS1.fasta

### Build BLAST database for each reference genome

makeblastdb -in C021p.fasta -dbtype nucl
makeblastdb -in C286p.fasta -dbtype nucl
makeblastdb -in C309p.fasta -dbtype nucl
makeblastdb -in C324p.fasta -dbtype nucl
makeblastdb -in CF12.fasta -dbtype nucl
makeblastdb -in CF13.fasta -dbtype nucl
makeblastdb -in H53.fasta -dbtype nucl
makeblastdb -in K091p.fasta -dbtype nucl
makeblastdb -in K147.fasta -dbtype nucl
makeblastdb -in K153_pOXA48.fasta -dbtype nucl
makeblastdb -in K163.fasta -dbtype nucl
makeblastdb -in K209p.fasta -dbtype nucl
makeblastdb -in K25.fasta -dbtype nucl

### Align IS1 element from pOXA-48 against its reference genome

blastn -db C021p.fasta -query C021_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_C021.tab
blastn -db C286p.fasta -query C286_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_C286.tab
blastn -db C309p.fasta -query C309_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_C309.tab
blastn -db C324p.fasta -query C324_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_C324.tab
blastn -db CF12.fasta -query CF12_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_CF12.tab
blastn -db CF13.fasta -query CF13_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_CF13.tab
blastn -db H53.fasta -query H53_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_H53.tab
blastn -db K091p.fasta -query K091_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_K091.tab
blastn -db K147.fasta -query K147_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_K147.tab
blastn -db K153_pOXA48.fasta -query K153_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_K153.tab
blastn -db K163.fasta -query K163_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_K163.tab
blastn -db K209p.fasta -query K209_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_K209.tab
blastn -db K25.fasta -query K25_pOXA_IS1.fasta -outfmt 6 > IS1_blast_results_K25.tab

