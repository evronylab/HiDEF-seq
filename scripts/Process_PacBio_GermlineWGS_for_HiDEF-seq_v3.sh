#!/bin/bash
#Usage: Process_PacBio_GermlineWGS_for_HiDEF-seq_v3.sh [input_bam] [output_basename] [reference.fasta] [reference.mmi] [hidef-seq sif path] [clair3 sif path] [deepvariant sif path]
input_bam=`readlink -f $1`
output_basename=`readlink -f $2`
reffasta=`readlink -f $3`
refmmi=`readlink -f $4`
hidefseq_sif=`readlink -f $5`
clair3_sif=`readlink -f $6`
deepvariant_sif=`readlink -f $7`

## i. Align with pbmm2
alignjobid=$(sbatch -c 32 -t 3600 --mem 96G --parsable --wrap "singularity exec -B `dirname $refmmi`:`dirname $refmmi` -B `dirname $input_bam`:`dirname $input_bam` -B `dirname $output_basename`:`dirname $output_basename` $hidefseq_sif bash -c 'source /hidef/miniconda3/etc/profile.d/conda.sh; conda activate /hidef/bin/pbconda; if [[ -s $output_basename.bam.bai ]]; then exit 0; fi; pbmm2 align -j 23 -J 8 -m 4G --preset CCS --sort $refmmi $input_bam $output_basename.bam'")

## ii. Call variants with Clair3
sbatch -t 2880 -c 32 --mem 64G --dependency=afterok:$alignjobid --wrap "if [[ -s $output_basename.clair3.vcf.gz.tbi ]]; then exit 0; fi; singularity run -B `pwd`:`pwd` -B `dirname $reffasta`:`dirname $reffasta` -B /usr/lib/locale/:/usr/lib/locale/ -B `dirname $output_basename`:`dirname $output_basename` $clair3_sif /opt/bin/run_clair3.sh --bam_fn=$output_basename.bam --ref_fn=$reffasta --threads=32 --platform=\"hifi\" --model_path=\"/opt/models/hifi_revio\" --include_all_ctgs --output=${output_basename}.clair3 --remove_intermediate_dir; mv ${output_basename}.clair3/merge_output.vcf.gz ${output_basename}.clair3.vcf.gz; mv ${output_basename}.clair3/merge_output.vcf.gz.tbi ${output_basename}.clair3.vcf.gz.tbi"

## iii. Call variants with DeepVariant
sbatch -t 3600 -c 16 --mem 64G --dependency=afterok:$alignjobid --wrap "if [[ -s $output_basename.deepvariant.vcf.gz.tbi ]]; then exit 0; fi; singularity run -B `pwd`:`pwd` -B `dirname $reffasta`:`dirname $reffasta` -B /usr/lib/locale/:/usr/lib/locale/ -B `dirname $output_basename`:`dirname $output_basename` $deepvariant_sif /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=$reffasta --reads=$output_basename.bam --output_vcf=$output_basename.deepvariant.vcf.gz --num_shards=16"
