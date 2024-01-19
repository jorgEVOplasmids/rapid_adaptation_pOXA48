#!/bin/bash

# Check ref genome file names and read sample names for execution

### PIPELINE COMMANDS FOR VARIANT CALLING IN EXPERIMENTAL EVOLUTION GENOMIC DATA ###

for folder in ~/data_EE/raw_data/*
do
	strain=data/$( basename $folder )
	
	## Quality control check and deduplication of raw data
	echo "=== Beginning QC and deduplication of "$strain" ==="
	
	for reads1 in ~/data_EE/raw_data/*_1.fastq.gz
	do
		reads2=${reads1%%_1.fastq.gz}"_2.fastq.gz"
		
		fastqc $reads1 -t 2 -o ~/data_EE/raw_data/fastqc_out
		fastqc $reads2 -t 2 -o ~/data_EE/raw_data/fastqc_out
	done
	
	multiqc ~/data_EE/raw_data/fastqc_out -o ~/data_EE/raw_data/fastqc_out/multiqc_out
	
	## Trim sequences and quality control check of trimmed reads
	
	for reads1 in ~/data_EE/raw_data/$strain/*_1.fastq.gz
	do
		reads2=${reads1%%_1.fastq.gz}"_2.fastq.gz"
		trim_galore -q 20 --length 50 --nextera --cores 2 -o ~/data_EE/trimmed_data/$strain --fastqc_args \"-t 2 -o ~/data_EE/fastqc_out/$strain/trimmed_data\" --paired $reads1 $reads2
	done

	multiqc ~/data_EE/fastqc_out/$strain/trimmed_data -o ~/data_EE/fastqc_out/$strain/trimmed_data/multiqc_out
	
	## De novo assembly in contigs
	
	for reads1 in ./reads_trimmed/$strain/*val_1.fq.gz
	do
		reads2=${reads1%%val_1.fq.gz}"val_2.fq.gz"
		sample=$( basename $reads1 )
		sample=$( echo ${sample%%_val_1.fq.gz} )
		sample=$( echo ${sample} | sed 's/^[0-9]_//g' )
		# echo $reads1 $reads2 $sample
		echo spades.py --isolate --cov-cutoff auto -1 $reads1 -2 $reads2 -o ./assembly_spades/$strain/$sample/
	done
	
	## Map clean reads vs reference and sort with samtools
	
	for reads1 in ~/data_EE/trimmed_data/$strain/*R1_001_val_1.fq.gz
	do
		reads2=${reads1%%R1_001_val_1.fq.gz}"R2_001_val_2.fq.gz"
		bwa mem -t 5 $strain/reference/MG1655_ptrna67.fasta $reads1 $reads2 > $strain/mapping/map_clumped_complete/"aln_"$strain".sam"
	done
		
	for file in ~/data_EE/mapping/$strain/*.sam
	do
		filename=$(basename -s .sam $file)
		#echo $filename
		echo samtools view -b $file \> ~/data_EE/mapping/$strain/$filename".bam"
	done

	for file in ~/data_EE/mapping/$strain/*.bam
	do

		filename=$(basename -s .bam $file)
		echo samtools sort $file -o ~/data_EE/mapping/$strain/$filename".sorted.bam"
	done
	
	## Variant calling using breseq and snippy
	
	for reads1 in ~/data_EE/raw_data/$strain/*1_001.fastq.gz
	do

		reads2=${reads1%%1_val_1.fq.gz}"2_val_2.fq.gz"
		strain=$( echo ${reads1%%_R1_001_val_1.fq.gz} | cut -d '/' -f 2 )
		breseq -r ~/data_EE/reference/MG1655_complete.gbk -j 5 -n $strain -o ~/data_EE/breseqs/$strain -p $reads1 $reads2
		snippy --minfrac 0.05 --report --outdir ~/data_EE/snippy/$strain --ref ~/data_EE/reference/ref.gbk --R1 $reads1 --R2 $reads2
	done
done
