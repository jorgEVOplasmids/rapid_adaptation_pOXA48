#!/usr/bin/env bash

data_dir=$1 # Where fast5 folder is placed
barcode=$2 # Native Barcoding kit employed EXP-NBD104 or EXP-NBD114

### Basecalling, QC and Demultiplexing of raw data

# Basecalling with guppy

guppy_basecaller -i $data_dir/fast5 -s $data_dir/fastq --flowcell FLO-MIN106 --kit SQK-LSK109 -x "cuda:0"

# QC of reads with minIONQC

Rscript ~/Documents/programs/minionQC/MinIONQC.R -i $data_dir -o $data_dir/minionQC -p 8

# Demultiplexing with guppy

guppy_barcoder -i $data_dir/fastq/pass/ -s $data_dir/fastq/pass_demult --barcode_kits $barcode -t 8

# Compress the fastq files to fastq.gz

for dir in $data_dir/fastq/pass_demult/barcode*
do
	barc=$( basename $dir )
	#echo $dir
	gzip $dir/*.fastq
	cat $dir/*.fastq.gz > $dir/merged_$barc.fastq.gz
done


### Variant calling using sniffles

# Map long reads against ref genome to create bam file

minimap2 -a $reffile.fasta $longreadfiles.fastq.gz > $mappingfile.sam

# Samtools for reformatting the mapping file

samtools view -b $mappingfile.sam > $mappingfile.bam
samtools sort $mappingfile.bam -o $mappingfile.sorted.bam --threads 16

# Run sniffles for variant calling

sniffles -i $mappingfile.sorted.bam -v $sniffles.vcf --reference $reffile.fasta

