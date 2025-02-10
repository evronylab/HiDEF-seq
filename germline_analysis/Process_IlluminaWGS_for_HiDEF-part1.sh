#!/bin/bash
#Usage: Process_IlluminaWGS_for_HiDEF-part1.sh sample.R1.fastq.gz sample.R2.fastq.gz

## For each fastq, extract individual fastq files for each read group from main fastq file.
for i in $@; do
sbatch -c 4 -t 32:00:00 --wrap "zcat $i | awk -v sample=`basename $i .fastq.gz` -F : '{if(NR%4==1){RG=\$3\".\"\$4} print \$0 | \"gzip >>\"sample\"_RG\"RG\".fastq.gz\"}'"
done
