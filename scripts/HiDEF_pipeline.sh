#!/bin/bash
#Usage: HiDEF_pipeline.sh [config.yaml] [step to run | all]
#  [step to run | all]: Specify one of the following steps:
#   process_subreads
#   make_bamfilezmw_all
#   qc_bamfilezmw_all
#   mutation_filtering
#   print_mutations
#   mutation_frequencies
#   mutation_signatures
#   gcloud_upload
#   all: Run all steps
#
#  Dispatches to SLURM the following scripts, with appropriate dependencies: process_subreads.sh, make_bamfilezmw_all.R, qc_bamfilezmw_all.R, mutation_filtering.R, print_mutations.R, mutation_frequencies.R, and copies output files to GCloud.
#  Prerequisites: SLURM, and see individual pipeline scripts

if [ ! -f $1 ]; then
    echo "Config YAML file not found!"
    exit 0
fi

export configyaml=$1
export scriptdir=`grep "HiDEFpipeline: " $configyaml | cut -d " " -f 2`
export process_subreads_output_path=`grep "process_subreads_output_path: " $configyaml | cut -d " " -f 2`
export analysisoutput_path=`grep "analysisoutput_path: " $configyaml | cut -d " " -f 2`
export gcloudoutput_path=`grep "gcloudoutput_path: " $configyaml | cut -d " " -f 2`
export gsutil_bin=`grep "gsutil_bin: " $configyaml | cut -d " " -f 2`
export slurm_himemjob_options=`grep "slurm_himemjob_options: " $configyaml | cut -d " " -f 2-`
export slurm_add_options=`grep "slurm_add_options: " $configyaml | cut -d " " -f 2-`

if [[ $2 == "process_subreads" || $2 == "all" ]]; then
  #First step of pipeline (process_subreads.sh) is a rapid script that does not need to be submitted to SLURM, but still runs within the singularity container.
  $scriptdir/process_subreads.sh $configyaml
  dependency="--dependency=afterok:`cat $process_subreads_output_path/process_subreads.DONE`"
fi

if [[ $2 == "make_bamfilezmw_all" || $2 == "all" ]]; then
  cmd=$(mktemp ./.XXXXXX)
  echo "$scriptdir/make_bamfilezmw_all.R $configyaml" > $cmd
  chmod +x $cmd
  make_bamfilezmw_all_job=$(sbatch $slurm_add_options $dependency -t 2000 --mem=128G -c 4 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")
  
  rm -f $process_subreads_output_path/process_subreads.DONE
  dependency="--dependency=afterok:$make_bamfilezmw_all_job"
fi

if [[ $2 == "qc_bamfilezmw_all" || $2 == "all" ]]; then
  cmd=$(mktemp ./.XXXXXX)
  echo "$scriptdir/qc_bamfilezmw_all.R $configyaml" > $cmd
  chmod +x $cmd
  qc_bamfilezmw_all_job=$(sbatch $slurm_add_options $dependency -t 240 --mem=96G -c 4 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")
fi

if [[ $2 == "mutation_filtering" || $2 == "all" ]]; then
  cmd=$(mktemp ./.XXXXXX)
  echo "$scriptdir/mutation_filtering.R $configyaml" > $cmd
  chmod +x $cmd
  mutation_filtering_job=$(sbatch $slurm_add_options $dependency -t 2400 $slurm_himemjob_options -c 4 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")
  dependency="--dependency=afterok:$mutation_filtering_job"
fi

if [[ $2 == "print_mutations" || $2 == "all" ]]; then
  cmd=$(mktemp ./.XXXXXX)
  echo "$scriptdir/print_mutations.R $configyaml" > $cmd
  chmod +x $cmd
  print_mutations_job=$(sbatch $slurm_add_options $dependency -t 720 $slurm_himemjob_options -c 1 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")
  dependency="--dependency=afterok:$print_mutations_job"
fi

if [[ $2 == "mutation_frequencies" || $2 == "all" ]]; then
  cmd=$(mktemp ./.XXXXXX)
  echo "$scriptdir/mutation_frequencies.R $configyaml" > $cmd
  chmod +x $cmd
  mutation_frequencies_job=$(sbatch $slurm_add_options $dependency -t 2400 $slurm_himemjob_options -c 4 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")
  dependency="--dependency=afterok:$mutation_frequencies_job"
fi

if [[ $2 == "mutation_signatures" || $2 == "all" ]]; then
  cmd=$(mktemp ./.XXXXXX)
  echo "$scriptdir/mutation_signatures.R $configyaml" > $cmd
  chmod +x $cmd
  mutation_signatures_job=$(sbatch $slurm_add_options $dependency -t 600 --mem=64G -c 4 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")
  dependency="--dependency=afterok:$mutation_signatures_job"
fi

if [[ $2 == "gcloud_upload" || $2 == "all" ]]; then
#GCloud upload operates on SLURM, but not within the singularity container, since GCloud cannot be easily configured for user to operate inside the singularity container.
  if [ ! -z "$gcloudoutput_path" ]; then gcloud_job=$(sbatch $slurm_add_options $dependency -t 600 --parsable --wrap "$gsutil_bin cp $process_subreads_output_path/*.final.bam* $gcloudoutput_path; $gsutil_bin cp $analysisoutput_path/*.bam* $gcloudoutput_path"); fi
fi
