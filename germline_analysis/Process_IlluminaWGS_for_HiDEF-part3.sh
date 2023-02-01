#!/bin/bash
#Usage: Process_IlluminaWGS_for_HiDEF-part3.sh outputbasename reference.fasta reference.chrsizes.bed [picard.jar path] [gatk path] [deepvariant sif path] [bedGraphToBigWig path]
#	- Merges all outputbasename*.bam files in the current directory, and continues with all downstream processing
# - reference.chrsizes.bed = 3-column tab-separated chr 0 chrsize

reffasta=`readlink -f $2`
refchrsizes=`readlink -f $3`
picard=`readlink -f $4`
gatk=`readlink -f $5`
deepvariant=`readlink -f $6`
bedGraphToBigWig=`readlink -f $7`

## i. Make CRAM file
#Combine read group BAM files.
mergejobid=$(sbatch -c 8 -t 1200 --parsable --wrap "samtools merge -@8 -O BAM $1.merged.bam `ls ${1}*.bam | xargs`")

#Sort by query name
sortqueryjobid=$(sbatch -t 1200 --parsable --dependency=afterok:$mergejobid --wrap "java -jar $picard SortSam SORT_ORDER=queryname INPUT=$1.merged.bam OUTPUT=$1.querysorted.bam R=$reffasta; rm $1.merged.bam")

#Mark duplicates
markdupjobid=$(sbatch -t 1200 --mem=64G --parsable --dependency=afterok:$sortqueryjobid --wrap "java -jar $picard MarkDuplicates I=$1.querysorted.bam O=$1.markdup.bam M=$1.markdup.metrics OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500; rm $1.querysorted.bam")

#Sort by coordinates.
sortcoordsjobid=$(sbatch -t 1200 --parsable --dependency=afterok:$markdupjobid --wrap "java -jar $picard SortSam CREATE_INDEX=TRUE SORT_ORDER=coordinate INPUT=$1.markdup.bam OUTPUT=$1.markdup.sorted.bam R=$reffasta; rm $1.markdup.bam")

#Convert to CRAM and index
cramjobid=$(sbatch -t 1200 --parsable --dependency=afterok:$sortcoordsjobid --wrap "samtools view -O CRAM -T $reffasta $1.markdup.sorted.bam > $1.cram; samtools index $1.cram; rm $1.*.bam* $1.*.bai")

## ii. Call variants with GATK
chroms="`seq 1 22 | xargs printf "chr%s\n" | xargs` chrX chrM"
numchroms=`echo $chroms | tr ' ' '\n' | wc -l`

gatk1jobid=$(sbatch -a 1-$numchroms -t 1660 -N 1 -c 4 --mem=24G --parsable --dependency=afterok:$cramjobid --wrap 'chrom=$(echo `seq 1 22 | xargs printf "chr%s\n" | xargs` chrX chrM | cut -f $SLURM_ARRAY_TASK_ID -d " "); '$gatk' --java-options "-Xmx20g" HaplotypeCaller -L $chrom -I '$1.cram' -O '$1'.GATK.g.vcf.gz.tmp.$chrom.vcf.gz -R '$reffasta' -ERC GVCF -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation')

gatk2jobid=$(sbatch -t 600 --mem=8G --parsable --dependency=afterok:$gatk1jobid --wrap "java -jar $picard MergeVcfs `echo $chroms | xargs printf -- "I=$1.GATK.g.vcf.gz.tmp.%s.vcf.gz\n" | xargs` O=$1.GATK.g.vcf.gz; rm $1.GATK.g.vcf.gz.tmp.*")

#Joint variant calling: Single sample GVCF -> Combined multi-sample GVCF -> Genotyping -> Final output VCF
#gatk_variant-calling_step2.sh [output VCF] [Reference fasta] [sample1 .g.vcf.gz] [sample2 .g.vcf.gz] ... [sampleN .g.vcf.gz]
gatk3jobid=$(sbatch -t 3600 --mem=32G --parsable --dependency=afterok:$gatk2jobid --wrap "$gatk CombineGVCFs -R $reffasta -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation --variant $1.GATK.g.vcf.gz -O $1.tmp.cohort.g.vcf.gz")
gatk4jobid=$(sbatch -t 3600 --mem=32G --parsable --dependency=afterok:$gatk3jobid --wrap "$gatk GenotypeGVCFs -R $reffasta -V $1.tmp.cohort.g.vcf.gz -O $1.GATK.vcf.gz; rm $1.tmp.cohort.g.vcf.gz*")

## iii. Call variants with DeepVariant script
sbatch -t 1440 -N 1 -n 16 --mem-per-cpu=4G --dependency=afterok:$cramjobid --wrap "singularity run -B `pwd`:`pwd` -B `dirname $reffasta`:`dirname $reffasta` -B /usr/lib/locale/:/usr/lib/locale/ $deepvariant /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$reffasta --reads=`readlink -f $1.cram` --regions=$refchrsizes --output_vcf=$1.deepvariant.vcf.gz --output_gvcf=$1.deepvariant.g.vcf.gz --num_shards=16"

##  iv. Make bigwig of total read coverage of germline reference
# Using the same types of filters used later for pileup filtering. This will ensure that the use of this bigwig for calculating the fraction of the genome that was filtered and properly calculating the mutation rate is correct, rather than simply using genomeCov from bedtools. Note using samtools mpileup here instead of bcftools mpileup, because we need an output that can be re-formatted into bedgraph and then into BigWig.
bigwig1jobid=$(sbatch -t 900 --parsable --dependency=afterok:$cramjobid --wrap "samtools mpileup -A -B -d 999999 --ff 1024 $1.cram -f $reffasta 2>/dev/null | awk '{print \$1 \"\t\" \$2-1 \"\t\" \$2 \"\t\" \$4}' > $1.coverage.bg")

#Sort bedgraph
bigwig2jobid=$(sbatch -t 900 --parsable --dependency=afterok:$bigwig1jobid --wrap "sort -k1,1 -k2,2n $1.coverage.bg > $1.coverage.sorted.bg; rm $1.coverage.bg")

#Bedgraph to bigwig
bigwig3jobid=$(sbatch -t 600 --parsable --dependency=afterok:$bigwig2jobid --wrap "tempfile=$(mktemp); cut -f 1,3 $refchrsizes > \$tempfile; $bedGraphToBigWig $1.coverage.sorted.bg \$tempfile $1.coverage.bw; rm $1.coverage.sorted.bg; rm \$tempfile")
