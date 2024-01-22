#!/bin/bash

### Iterate through trimmed reads folder with samples to analyze

for folder in reads_trimmed/*/
do
	strain=$( basename $folder )
	# mkdir plasmid_copy_number/$strain
	
	# Indicate chromosome and plasmid contigs using anvio reformat of fasta headers to ease the identification of chromosome and plasmids
	headers=$(<anvio/$strain/$strain-contigs-reformat.txt)
	chrcontig=$( echo "$headers" | grep -E "c_000000000001" | tr "\t" " " | cut -d " " -f 1 )
	contigs=$( echo "$headers" | grep -E "c_00000000000[^1]" | tr "\t" " " | cut -d " " -f 1 )
	
	# Save results file header in output file
	echo -e 'Strain\tSample\tp_contig\tchr_mean\tchr_median\tchr_std\tp_mean\tp_median\tp_std\tmean_ratio\tmedian_ratio' >> plasmid_copy_number/$strain-pcn.tsv
	
	# Calculate statistics for each contig in each sample of the strain
	for sortedbam in mapping/$strain/*.sorted.bam
	do
		samplename=$( basename $sortedbam )
		samplename=${samplename::-11}
		
		# Calculate chromosome statistics
		chrstats=$( samtools depth -r $chrcontig -a $sortedbam | datamash -R 2 mean 3 median 3 sstdev 3 )
		chrstats=${chrstats//,/.} # to change decimal separator
		chrmean=$( echo $chrstats | awk '{print $1}' )
		chrmedian=$( echo $chrstats | awk '{print $2}' )
		chrsd=$( echo $chrstats | awk '{print $3}' )
		for pcontig in $contigs
		do
			# Calculate plasmid statistics
			pstats=$( samtools depth -r $pcontig -a $sortedbam | datamash -R 2 mean 3 median 3 sstdev 3 )
			pstats=${pstats//,/.}
			pmean=$( echo $pstats | awk '{print $1}' )
			pmedian=$( echo $pstats | awk '{print $2}' )
			psd=$( echo $pstats | awk '{print $3}' )
			
			# Calculate the ratio of coverage plasmid/chromosome
			ratiomean=$(echo "scale=4; x=($pmean/$chrmean); if(x<1) print 0; x" | bc)
			ratiomedian=$(echo "scale=4; x=($pmedian/$chrmedian); if(x<1) print 0; x" | bc )
			
			# Save results in output tsv file
			echo -e $strain'\t'$samplename'\t'$pcontig'\t'$chrmean'\t'$chrmedian'\t'$chrsd'\t'$pmean'\t'$pmedian'\t'$psd'\t'$ratiomean'\t'$ratiomedian >> plasmid_copy_number/$strain-pcn.tsv
		done
	done
done
