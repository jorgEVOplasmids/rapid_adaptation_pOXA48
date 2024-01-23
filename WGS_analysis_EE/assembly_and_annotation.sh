#!/bin/bash

# Close the genomes by hybrid assembly using Unicycler

unicycler -l $longreads -1 $forwardshortreads -2 $reverseshortreads -o $output_ref_genome_folder -t 16 --racon_path $path_to_racon

# If Unicycler failed, try Flye

flye --nano-raw $longreads --threads 16 --plasmids --out-dir $output_ref_genome_folder

# Perform polishing rounds combining consensus sequence of Medaka and polishing with Pilon

medaka_consensus -i $longreads -d $output_ref_genome_folder/assembly.fasta -o $output_ref_genome_folder/medaka

bwa index $output_ref_genome_folder/medaka/consensus.fasta

bwa mem -t 16 $output_ref_genome_folder/medaka/consensus.fasta $forwardshortreads $reverseshortreads | samtools sort -o $output_ref_genome_folder/medaka/alignments.bam

samtools index $output_ref_genome_folder/medaka/alignments.bam

java -Xmx16G -jar pilon-1.24.jar --changes --genome $output_ref_genome_folder/medaka/consensus.fasta --frags $output_ref_genome_folder/medaka/alignments.bam --output $strain.1 --outdir $output_ref_genome_folder/pilon/

# If changes reported by pilon are 0, save assembly. Otherwise, repeat polishing rounds by using the following lines

bwa index $output_ref_genome_folder/pilon/$strain.1.fasta

bwa mem -t 16 $output_ref_genome_folder/pilon/$strain.1.fasta $forwardshortreads $reverseshortreads | samtools sort -o $output_ref_genome_folder/medaka_2/alignments.bam

samtools index $output_ref_genome_folder/medaka_2/alignments.bam

java -Xmx16G -jar pilon-1.24.jar --changes --genome $output_ref_genome_folder/pilon/$strain.1.fasta --frags $output_ref_genome_folder/medaka_2/alignments.bam --output $strain.2 --outdir $output_ref_genome_folder/pilon_2/

# Once the genome is closed, annotate using PGAP

# cd to folder with executable .py
# input.yaml contains info of input genome fasta file path and metadata as explained in https://github.com/ncbi/pgap/wiki/Input-Files

./pgap.py -r -o $output_annotations/$strain $output_ref_genome_folder/input.yaml

