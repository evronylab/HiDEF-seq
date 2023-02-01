#!/bin/bash

##USER CONFIGURATION
####################
#Paths to HPC singularity SLURM config and hidef-seq .sif singularity image
export SINGULARITY_BIN=[Full path to singularity binary on HPC SLURM cluster]
export HIDEF_IMAGE=[Full path to hidef-seq singularity image file (.sif) on HPC SLURM cluster]

#General file systems and required slurm paths to mount to singularity container
export GENERAL_FILE_SYSTEMS=[PATH1],[PATH2],...
export SLURM_REQUIRED_PATHS=$(which sbatch),$(which sacct)

#SLURM bin folder
export SLURM_BIN_FOLDER=$(dirname $(which sbatch))

##SCRIPT
########
args=''
for i in "$@"; do 
  i="${i//\\/\\\\}"
  args="$args \"${i//\"/\\\"}\""
done

if [ "$args" == "" ]; then args="/bin/bash"; fi

#Make empty directory to bind to home, to exclude local HOME directory. Also copy gcloud credential folder.
mkdir -p $HOME/.tmphome
mkdir -p $HOME/.tmphome/.config
cp -r $HOME/.config/gcloud $HOME/.tmphome/.config/

# paths to pipeline scripts inside container and slurm bin folder
export SINGULARITYENV_PREPEND_PATH=/hidef/scripts:$SLURM_BIN_FOLDER

if [ "$NESTED_SINGULARITY_RUN" == "YES" ]; then

  export SINGULARITYENV_NESTED_SINGULARITY_RUN=$NESTED_SINGULARITY_RUN
  export SINGULARITYENV_HIDEF_SINGULARITY_WRAPPER=$HIDEF_SINGULARITY_WRAPPER

  $SINGULARITY_BIN \
    exec --home $HOME/.tmphome $HIDEF_IMAGE \
    /bin/bash -c "eval $args"

elif [ "$NESTED_SINGULARITY_RUN" == "" ]; then
  export NESTED_SINGULARITY_RUN="YES"

  #bind paths
  export SINGULARITY_BINDPATH=$GENERAL_FILE_SYSTEMS,$SLURM_REQUIRED_PATHS,$PWD

  # path to this singularity wrapper
  export HIDEF_SINGULARITY_WRAPPER=$(realpath $0)

  # user permissions file for slurm
  passwd=$PWD/.passwd-$USER
  if [ ! -e $passwd ]; then
      touch $passwd
      chmod 600 $passwd
      getent passwd slurm $USER > $passwd
  fi
  export SINGULARITY_BINDPATH=${passwd}:/etc/passwd:ro,$SINGULARITY_BINDPATH

  export SINGULARITYENV_NESTED_SINGULARITY_RUN=$NESTED_SINGULARITY_RUN
  export SINGULARITYENV_HIDEF_SINGULARITY_WRAPPER=$HIDEF_SINGULARITY_WRAPPER

  $SINGULARITY_BIN \
    exec --home $HOME/.tmphome -e $HIDEF_IMAGE \
    /bin/bash -c "eval $args"

fi
