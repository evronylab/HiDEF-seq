#!/bin/bash
#Usage: Process_IlluminaWGS_for_HiDEF-part2.sh outputbasename reference.fasta sample_R1_RGXXXXX.fastq.gz sample_R2_RGXXXXX.fastq.gz
#	 - Run this script for each read group

#Align read group
readgroup=`echo $3 | sed 's/.*_RG//;s/\.fastq\.gz//'`

sbatch -t 3400 -N 1 -n 12 --mem-per-cpu=3000 --wrap "bwa mem -R '@RG\\tID:'$readgroup'\\tSM:'$1'\\tLB:'$1'\\tPL:ILLUMINA' -t 8 $2 $3 $4 | samtools view -@4 -b - > $1.$readgroup.bam; rm $3 $4"
