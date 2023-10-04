## Configuration of the HiDEF-seq pipeline for the SLURM cluster
HiDEF-seq and its wrapper script that runs the pipeline needs to be configured for the SLURM cluster, per these instructions.

### A. Download HiDEF-seq docker/singularity image:
- Download with singularity: ```singularity pull docker://gevrony/hidef-seq:latest```

### B. Configure HiDEF-seq wrapper script: 
- A wrapper script (run-hidef-seq.bash) is used to run the HiDEF-seq singularity image, because this enables submission of SLURM jobs from within the singularity image.
- The wrapper script needs to be configured for each SLURM computing cluster, per the below instructions.
1. Download the [run-HiDEF-seq.bash wrapper script template](https://github.com/gevro/HiDEF-seq/raw/main/scripts/run-HiDEF-seq.bash)
2. Make the wrapper script executable: ```chmod +x run-hidef-seq.bash```
3. Configure SINGULARITY_BIN, HIDEF_IMAGE, and GENERAL_FILE_SYSTEMS parameters, as follows::
    - export SINGULARITY_BIN=[Full path to singularity binary on HPC SLURM cluster]
    - export HIDEF_IMAGE=[Full path to hidef-seq singularity image file (.sif) on HPC SLURM cluster]
    - export GENERAL_FILE_SYSTEMS=[PATH1],[PATH2],...
      - This is a comma-separated list of the mounted volumes and/or folders made available to the singularity image while it is running. This should include only the volumes and/or folders in which the pipeline's input files are located, as specified in the config [YAML].
      - Never mount any volume or folder that begins with /usr or /bin. This will cause conflicts with required folders inside the singularity image.
      - Note that top-level volumes and folders will also include all nested child folders. For example, /scratch will include /scratch/input/bam
    
4. Configure SLURM_REQUIRED_PATHS and SLURM_BIN_FOLDER so that SLURM will work from inside the singularity image.

    **== Background ==**
      
      i. export SLURM_REQUIRED_PATHS=[PATH1],[PATH2],...  
      - This is a comma-separated list of all the paths required by SLURM that were not already included in GENERAL_FILE_SYSTEMS. These are made available to the singularity image while it is running.
      - You should specify the most specific (lower level) paths in order to avoid conflicts between these paths and paths inside the singularity image. For example, if /var/etc/file1 is required for the SLURM configuration, specify /var/etc/file1 rather than /var/etc. However, if the entire /var/etc folder is required, then specify /var/etc.
      - If there is a higher-level folder already specified in GENERAL_FILE_SYSTEMS that includes a folder in SLURM_REQUIRED_PATHS at any lower level of the folder hierarchy, then it does not need to be specified in SLURM_REQUIRED_PATHS. For example, /cm/share/apps/slurm does not need to be specified in SLURM_REQUIRED_PATHS, if /cm is already specified in GENERAL_FILE_SYSTEMS.
      - SLURM_REQUIRED_PATHS must include the path of $SLURM_BIN_FOLDER (see below) or a higher-level folder that includes $SLURM_BIN_FOLDER.

      ii. export SLURM_BIN_FOLDER=[Full path to folder containing both sbatch and sacct bin files of SLURM]
      - The specified folder(s) is added to the $PATH variable when the singularity container runs.
      
      **== Configuration process ==**
    
      The below process for configuring SLURM_REQUIRED_PATHS and SLURM_BIN_FOLDER should work for most HPCs. Most of the process involves iteratively running the sacct command from within the singularity image to identify missing folders/files that need to be added to SLURM_REQUIRED_PATHS, until the sacct command works from within the singularity image.
    
      i. SLURM_REQUIRED_PATHS and SLURM_BIN_FOLDER should start off set to the below in run-hidef-seq.bash (this is the default setting in the run-hidef-seq.bash template):
        ```
        export SLURM_REQUIRED_PATHS=$(which sbatch),$(which sacct)
        export SLURM_BIN_FOLDER=$(dirname $(which sbatch))
        ```
          
      - Modify SLURM_REQUIRED_PATHS and SLURM_BIN_FOLDER to sbatch and sacct binary file locations if sbatch and sacct are not found in $PATH (e.g. if the commands 'which sbatch' or 'which sacct' do not work in your HPC login node). This should be a very uncommon situation, and indicates an issue with your HPC SLURM configuration (contact your HPC admins for assistance).
      - If more than one folder needs to be specified for SLURM_BIN_FOLDER (i.e. if sbatch and sacct are in different folders), they should be separated by : (colon symbol).
    
    ii. On your SLURM login node, find the location of SLURM's configuration file by running: scontrol show conf | grep SLURM_CONF
      - Example output: ```SLURM_CONF = /opt/slurm/etc/slurm.conf```
      - Next, check if the slurm.conf file exists:
          ls [Full path to slurm.conf file from above]
      - ONLY IF the file exists (i.e. you do not get a 'No such file or directory' error), add the full path of the SLURM configuration file (e.g. /opt/slurm/etc/slurm.conf) to SLURM_REQUIRED_PATHS:
        e.g. ```export SLURM_REQUIRED_PATHS=$(which sbatch),$(which sacct),/opt/slurm/etc/slurm.conf```
    
    iii. Run: [full path to run-HiDEF-seq.bash] sacct
      - You will likely get this error:
          ```
          sacct: error while loading shared libraries: libslurmfull.so: cannot open shared object file: No such file or directory
          ```
      - This indicates that you need to find the location of libslurmfull.so

    iv. Run: ```ldd $(which sacct) | grep libslurmfull.so```
      - Example output: ```libslurmfull.so => /opt/slurm/lib64/slurm/libslurmfull.so (0x00007fdf7373b000)```
      - Add the full path of the FOLDER containing libslurmfull.so (e.g. /opt/slurm/lib64/slurm) to SLURM_REQUIRED_PATHS:
        e.g. ```export SLURM_REQUIRED_PATHS=$(which sbatch),$(which sacct),/opt/slurm/etc/slurm.conf,/opt/slurm/lib64/slurm```
        
      - Important: For libslurmfull.so, you must add the FOLDER containing libslurmfull.so, not just /opt/slurm/lib64/slurm/libslurmfull.so
      - Note 1: This will likely be the only required path in SLURM_REQUIRED_PATHS that is a folder instead of a single file
      - Note 2: You may also see additional output from the above command, which you can ignore:
          ```
          ldd: ./alias: No such file or directory
          ldd: ./sacct='sacct: No such file or directory
          ldd: ./-X': No such file or directory
          ```
      
    v. Again, run: ```[full path to run-HiDEF-seq.bash] sacct```
      - You will likely get an error similar to this:
        ```
        sacct: error: plugin_load_from_file: dlopen(/opt/slurm/lib64/slurm/auth_munge.so): libmunge.so.2: cannot open shared object file: No such file or directory
        sacct: error: Couldn't load specified plugin name for auth/munge: Dlopen of plugin file failed
        sacct: error: cannot create auth context for auth/munge
        sacct: error: slurm_send_node_msg: auth_g_create: REQUEST_PERSIST_INIT has authentication error
        sacct: error: slurm_persist_conn_open: failed to send persistent connection init message to slurm-1:6819
        sacct: error: Sending PersistInit msg: Protocol authentication error
        sacct: error: Problem talking to the database: Protocol authentication error
        ```
        
      - This indicates that you need to find the location of libmunge.so.2
      
    vi. Run: ```ldd [full path to auth_munge.so from the output of the prior step] | grep libmunge.so.2```
      - Example output: ```libmunge.so.2 => /usr/lib64/libmunge.so.2 (0x00007f27e660a000)```
      - Add the full path of the file libmunge.so.2 (e.g. /usr/lib64/libmunge.so.2) to SLURM_REQUIRED_PATHS:
        e.g. ```export SLURM_REQUIRED_PATHS=$(which sbatch),$(which sacct),/opt/slurm/etc/slurm.conf,/opt/slurm/lib64/slurm,/usr/lib64/libmunge.so.2```
        
    vii. Again, run: ```[full path to run-HiDEF-seq.bash] sacct```
      - You will likely get an error similar to this:
        ```
        sacct: error: If munged is up, restart with --num-threads=10
        sacct: error: Munge encode failed: Failed to access "/var/run/munge/munge.socket.2": No such file or directory
        sacct: error: slurm_send_node_msg: auth_g_create: REQUEST_PERSIST_INIT has authentication error
        sacct: error: slurm_persist_conn_open: failed to send persistent connection init message to slurm-1:6819
        sacct: error: Sending PersistInit msg: Protocol authentication error
        sacct: error: Problem talking to the database: Protocol authentication error
        ```
        
      - Add the full path of the file munge.socket.2 shown in the error message (e.g. /var/run/munge/munge.socket.2) to SLURM_REQUIRED_PATHS:
        e.g. ```export SLURM_REQUIRED_PATHS=$(which sbatch),$(which sacct),/opt/slurm/etc/slurm.conf,/opt/slurm/lib64/slurm,/usr/lib64/libmunge.so.2,/var/run/munge/munge.socket.2```
        
    viii. Again, run: ```[full path to run-HiDEF-seq.bash] sacct```
      - Now the sacct command is working (example output):
        ```
        JobID           JobName  Partition    Account  AllocCPUS      State ExitCode 
        ------------ ---------- ---------- ---------- ---------- ---------- -------- 
        15095929           wrap         cl      users          4 CANCELLED+      0:0 
        15095929.ba+      batch                 users          4  CANCELLED     0:15 
        15095929.ex+     extern                 users          4  COMPLETED      0:0 
        ...
        ```
      
    ix. Run one more final test, to submit a job via sbatch via the singularity wrapper script: ```[full path to run-HiDEF-seq.bash] sbatch -t 1 --wrap "ls"```
      - Example output: ```"Submitted batch job 16190753"```
      - Confirm that the sbatch job completed successfully
      
    x. If the above process does not work, continue to run ```[full path to run-HiDEF-seq.bash] sacct```, iteratively, and add files to SLURM_REQUIRED_PATHS that are missing per the error messages.
      - Your HPC admins should be able to advise you on any other files required by SLURM to run that need to be made available to the singularity image via the SLURM_REQUIRED_PATHS configuration.
