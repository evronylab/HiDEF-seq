# HiDEF-seq
HiDEF-seq is a single-molecule sequencing method with single-molecule accuracy for single base substitutions, when present in either one or both strands of DNA. This repository contains scripts and a pipeline for analysis of HiDEF-seq data. The protocol for creating and sequencing HiDEF-seq libraries is detailed in the HiDEF-seq manuscript.

Before starting HiDEF-seq analysis, every sample requires germline sequencing data that must be processed prior to running the HiDEF-seq pipeline.

Then, for HiDEF-seq data, the fastest and most straightforward way to analyze it is using the pre-configured docker image. However, individual scripts for each step of the analysis are also provided.


### Outline
- [Computing environment](#computing-environment)
- [Reference genome](#reference-genome)
- [Germline sequencing data processing](#germline-sequencing-data-processing)
  - [Script requirements](#a-script-requirements)
  - [Illumina germline sequencing data processing](#b-illumina-germline-sequencing-data-processing)
  - [PacBio germline sequencing data processing](#c-pacbio-germline-sequencing-data-processing)
- [Configure HiDEF-seq wrapper script for SLURM cluster](#configure-hidef-seq-for-the-slurm-cluster)
- [Run HiDEF-seq Pipeline](#run-hidef-seq-pipeline)
- [Outputs](#outputs)
- [Citation](#citation)

## Computing environment
The analysis requires a SLURM computing cluster.

If a SLURM cluster is not available to you locally, you may setup a virtual SLURM cluster on Google Cloud (https://github.com/SchedMD/slurm-gcp) or on AWS (https://docs.aws.amazon.com/parallelcluster/latest/ug/what-is-aws-parallelcluster.html). Note: this incurs cloud-associated costs.

The cluster must also have Singularity version >= 3.7.1.

## Reference genome
The analyses are performed using the CHM13 T2T v1.0 reference genome (https://github.com/marbl/CHM13), which improves analysis accuracy.

### A. Script requirements
- [bwa](https://github.com/lh3/bwa) in the system PATH 
- [samtools](http://www.htslib.org/download/) in the system PATH

### B. Preparing reference genome files:
1. Download the [CHM13 v1.0 reference genome](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.0.fasta.gz)
2. Extract with gunzip: ```gunzip chm13.draft_v1.0.fasta.gz```
3. Make an 'fai' index file: ```samtools faidx chm13.draft_v1.0.fasta```
4. Make a BED file of chromosome sizes: ```awk '{print $1 "\t0\t" $2}' chm13.draft_v1.0.fasta.fai > chm13.draft_v1.0.chrsize.bed```
5. Make a bwa index: ```bwa index chm13.draft_v1.0.fasta```
6. Download additional required references files:
    - [clip_all.CHM13_v1.bw](https://storage.googleapis.com/hidef-seq-references/clip_all.CHM13_v1.bw)
    - [orphan_cov_all.CHM13_v1.bw](https://storage.googleapis.com/hidef-seq-references/orphan_cov_all.CHM13_v1.bw)
    - [prop_cov_all.CHM13_v1.bw](https://storage.googleapis.com/hidef-seq-references/prop_cov_all.CHM13_v1.bw)

## Germline sequencing data processing
Prior to starting the HiDEF-seq pipeline, germline sequencing data (either Illumina OR PacBio), must be processed for each sample.

The scripts for performing this are in the germline_analysis folder, and instructions for running them are below.

### A. Script requirements
- [bwa](https://github.com/lh3/bwa) in the system PATH 
- [samtools](http://www.htslib.org/download/) in the system PATH
- singularity (for running DeepVariant)
- [GATK](https://github.com/broadinstitute/gatk/releases)
- DeepVariant: ```singularity pull docker://google/deepvariant:latest```. This makes an SIF file ('singularity image file')
- [Picard .jar file](https://github.com/broadinstitute/picard/releases/tag/2.27.5)
- [bedGraphToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/)

### B. Illumina germline sequencing data processing
1. Download the [Illumina germline analysis scripts](https://github.com/evronylab/HiDEF-seq/tree/main/germline_analysis)

2. For each sample's R1 and R2 sequencing data, extract the reads of each read group into a separate FASTQ file:

    ```Process_IlluminaWGS_for_HiDEF-part1.sh [samplename.R1.fastq.gz] [samplename.R2.fastq.gz]```

3. For each extracted read group FASTQ file, align the reads to the genome with bwa:

    ```Process_IlluminaWGS_for_HiDEF-part2.sh [samplename] [chm13.draft_v1.0.fasta] [samplename_R1_RGXXXXX.fastq.gz] [samplename_R2_RGXXXXX.fastq.gz]```
      [samplename]: sample name; should be the same for all read groups that are extracted from the same sample
      XXXXX = Read group ID

    - Run this script for each read group

4. For each sample, Merge BAM files of all the sample's read groups, mark duplicates, sort, convert to CRAM, call variants with GATK and DeepVariant, and make a bigWig file of read coverage across the genome.

    ```Process_IlluminaWGS_for_HiDEF-part3.sh [samplename] [reference.fasta] [chm13.draft_v1.0.chrsize.bed] [picard.jar path] [gatk path] [deepvariant sif path] [bedGraphToBigWig path]```

### C. PacBio germline sequencing data processing
1. Install these two pbbioconda tools: pbmerge, pbmm2 (https://github.com/PacificBiosciences/pbbioconda)
2. Index reference genome: ```pbmm2 index [chm13.draft_v1.0.fasta] [chm13.draft_v1.0.mmi]
3. Generate HiFi reads (at least 3 SMRTcells per sample)
4. Merge HiFi reads from all SMRTcells of each sample with 'pbmerge'
5. Align reads to the genome: ```pbmm2 align --log-level INFO -j [cpus] [chm13.draft_v1.0.mmi] [samplename.merged.hifi.bam] [samplename.merged.hifi.aligned.bam] --preset CCS --sort```
6. Call variants with DeepVariant (enter below configuration parameters before running):
  ```
  REFFASTA=[chm13.draft_v1.0.fasta]
  REFCHRSIZES=[chm13.draft_v1.0.chrsize.bed]
  INPUTDIR=[path of input BAM files]
  OUTPUTDIR=[path for output BAM files]
  INPUTFILE=[input BAM]
  DEEPVARIANT=[path to deepvariant sif file]
  sbatch -c 16 -t 1200 --wrap "singularity run -B /usr/lib/locale/:/usr/lib/locale/,$REFDIR:$REFDIR,$INPUTDIR:$INPUTDIR,$OUTPUTDIR:$OUTPUTDIR $DEEPVARIANT /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO \
  --ref=$REFFASTA --reads=${INPUTDIR}/${INPUTFILE} \
  --regions=$REFCHRSIZES \
  --output_vcf=${OUTPUTDIR}/`basename ${INPUTFILE} .bam`.deepvariant.vcf.gz \
  --output_gvcf=${OUTPUTDIR}/`basename ${INPUTFILE} .bam`.deepvariant.g.vcf.gz \
  --num_shards=16"
  ```
7. Make a bigWig file of read coverage across the genome (enter below configuration parameters before running):
  ```
  HIFIALIGNEDBAM=[samplename.merged.hifi.aligned.bam]
  SAMPLENAME=`basename HIFIALIGNEDBAM .merged.hifi.aligned.bam`
  REFFASTA=[chm13.draft_v1.0.fasta]
  REFCHRSIZES=[chm13.draft_v1.0.chrsize.bed]
  bedGraphToBigWig=[bedGraphToBigWig path]
  
  #mpileup
  bigwig1jobid=$(sbatch -t 900 --parsable --wrap "samtools mpileup -A -B -d 999999 --ff 1024 $HIFIALIGNEDBAM -f $REFFASTA 2>/dev/null | awk '{print \$1 \"\t\" \$2-1 \"\t\" \$2 \"\t\" \$4}' > $SAMPLENAME.coverage.bg")

  #Sort bedgraph
  bigwig2jobid=$(sbatch -t 900 --parsable --dependency=afterok:$bigwig1jobid --wrap "sort -k1,1 -k2,2n $SAMPLENAME.coverage.bg > $SAMPLENAME.coverage.sorted.bg; rm $SAMPLENAME.coverage.bg")

  #Bedgraph to bigwig
  bigwig3jobid=$(sbatch -t 600 --parsable --dependency=afterok:$bigwig2jobid --wrap "tempfile=$(mktemp); cut -f 1,3 $REFCHRSIZES > \$tempfile; $bedGraphToBigWig $SAMPLENAME.coverage.sorted.bg \$tempfile $SAMPLENAME.coverage.bw; rm $SAMPLENAME.coverage.sorted.bg; rm \$tempfile")
  ```

## Configure HiDEF-seq for the SLURM cluster
Install and configure HiDEF-seq for the SLURM cluster, per these [instructions](https://github.com/evronylab/HiDEF-seq/blob/main/docs/slurm_configuration.md).

This only needs to be performed once.

## Run HiDEF-seq Pipeline
### A. Subreads BAM file input:
- The sequencing data input is a subreads.bam file from a single run of the PacBio Sequel instrument, and its associated 'pbi' index file.
- More than one sample can be multiplexed in each sequencing run, so that the subreads BAM file may contain more than one sample.

### B. YAML configuration file: This file contains all the parameters for the pipeline.
- Slightly different configurations are used for samples with Illumina versus PacBio germline sequencing data, and for chr1-22,X vs chrM analysis.
- Template YAML configuration files for each of these 4 scenarios, and documentation of each parameter, are provided [here](config_templates).
1. Create one chr1-22X and one chrM configuration file *for each* PacBio HiDEF-seq sequencing run (i.e. for each subreads BAM file).
 
### C. Run pipeline:
1. chr1-22X analysis: ```[Full path to run-HiDEF-seq.bash] HiDEF_pipeline.sh [chr1-22X YAML config file] all```
2. After chr1-22X analysis outputs the RDS file specified by [bamfilezmw_all_filename] in the YAML config file (i.e., after the process_subreads command completes), you can start the chrM analysis: ```[Full path to run-HiDEF-seq.bash] HiDEF_pipeline.sh [chrM YAML config file] all_except_process_subreads_and_gcloud_upload```

- Note, the following commands are available to run:
  ```
  [Full path to run-HiDEF-seq.bash] HiDEF_pipeline.sh [YAML config file] [command]
    [command]:
      process_subreads
      make_bamfilezmw_all
      qc_bamfilezmw_all
      mutation_filtering
      print_mutations
      mutation_frequencies
      gcloud_upload
      all: Run all steps
      all_except_process_subreads_and_gcloud_upload: Run all steps except process_subreads and gcloud_upload (used for chrM analysis)
  ```

## Outputs
Note: [items in brackets] refer to parameters defined in the YAML config file.

### A. Primary data processing:
In the [process_subreads_output_path] directory:
1. [ccs_BAM_prefix].[samplename].ccs.demux.[barcodename].postccsfilter.aligned.final.bam (and .bai/.pbi indexes): aligned BAM, after creating consensus sequence, basic molecule filters, demultiplexing, and alignment.
2. [bamfilezmw_all_filename] (file name is typically configured to end with .bamfilezmw_all.RDS): Aligned consensus read data processed into a format for analysis in R.
3. make_bamfilezmw_all.output.[timestamp].txt: Log of the script making [bamfilezmw_all_filename]

In the [process_subreads_output_path]/logs directory:

4. subreads.zmwcount.txt: Number of ZMWs (i.e. molecules) in raw data
5. subreads.cxfiltered.zmwcount.txt: Number of ZMWs after filtering ZMWs with 'cx' tag != 3 (ZMWs for which an adaptor was not detected on both ends)
6. [ccs_BAM_prefix].ccsreport.chunk#.txt: ccs (consensus sequence building tool) report for each chunk, plain text format
7. [ccs_BAM_prefix].ccsreport.chunk#.json: ccs (consensus sequence building tool) report for each chunk, json format
8. [ccs_BAM_prefix].ccsmetrics.chunk#.json: ccs (consensus sequence building tool) metrics for each chunk, json format
9. ccs.zmwcount.txt: Number of ZMWs after ccs
10. [ccs_BAM_prefix].ccs.demux.lima.report/counts/summary: metrics from lima demultiplexing
11. [barcodename].lima.zmwcount.txt: Number of ZMWs for each barcode after lima demultiplexing
12. [ccs_BAM_prefix].[samplename].ccs.demux.[barcodename].bam.ccsfilterstats: Statistics of filters applied after ccs
13. [ccs_BAM_prefix].[samplename].ccs.demux.[barcodename].postccsfilter.aligned.bam.pbmm2filterstats: Statistics of filters applied after pbmm2 alignment

### B. QC after primary data processing:
In the [analysisoutput_path] directory:
1. [samplename].[genomereference].[tag].hist.pdf: histograms of ZMW tag values
2. [samplename].[genomereference].[tag].fwdvsrev.scatter.pdf: Scatter plots between the two strands of a molecule
3. [samplename].[genomereference].[tag1].[tag2].scatter.pdf: Scatter plots of all possible pairings of [tag1] versus [tag2]
  
    Where [tag] = 
    - rq: Predicted average read accuracy (see [ccs documentation](https://ccs.how/faq/bam-output.html))
    - ec: Effective coverage per strand (see [ccs documentation](https://ccs.how/faq/bam-output.html))
    - np: Number of full-length passes per strand (see [ccs documentation](https://ccs.how/faq/bam-output.html))
    - qwidth: The length of the ZMW sequence
    - isize: The length of the sequence alignment in the reference space (the sum of the M/D/=/X CIGAR lengths)
    - mapq: mapping quality
  
4. [samplename].[genomereference].[VCF-file-name].[VCF tag].hist.pdf: histogram of tag values for each germline VCF

    Where [VCF tag] = 
    - Depth: Total read depth at locus
    - VAF: Variant allele fraction
    - GQ: Genotype quality
  
### C. Call filtering:
In the [analysisoutput_path] directory:
1. [analysisoutput_basename].RDS: Post-call filtering data in format for analysis in R.
2. mutation_filtering.output.[timestamp].txt: Log of the script making [analysisoutput_basename].RDS
3. print_mutations.output.txt and print_mutations.output.RDS: List of ssDNA mismatches and dsDNA mutations after call filtering, in txt and in RDS format suitable for analysis in R
    - [Reference for output fields](docs/print_mutations.md)
4. [samplename].[genomereference].[threshold-index].subreads.aligned.bam (and .bai/.pbi indexes): subreads with post-filtering calls aligned to the genome, for each [threshold-index] (i.e. index # defining a set of thresholds) defined in the YAML config file.

### D. Call burdens (i.e. mutation frequencies):
In the [analysisoutput_path] directory:
1. mutation_frequencies.RDS: Post-mutation frequencies script data in format for analysis in R
2. mutation_frequencies.output.[timestamp].txt: Log of the script making mutation_frequencies.RDS
3. mutation_frequencies.table.output.[timestamp].txt: All burden-related results in table format
    - [Reference for output fields](docs/mutation_frequencies.table.md)

### E. Mutational signatures:
  - Note that non-collapsed ssDNA spectra/signatures are plotted with 192-trinucleotide contexts with the label 'Transcribed' signifying central pyrimidine and 'Untranscribed' signifying central purine contexts (reverse complement of annotated central pyrimidine context). dsDNA spectra/signatures are plotted with the standard 96-trinucleotide context. Collapsed ssDNA spectrua/signatures are plotted with 96-trinucleotide contexts by summing each central pyrimidine context with its reverse complement central purine context.
  - Note, sigfit signature analysis below is automated output that may work with a broader range of sample types. Analyses in the paper utilized different settings. Settings can be adjusted by modifying the script.

In the [analysisoutput_path] directory:
1. [analysisoutput_basename].[samplename].[genomereference].[threshold-index].subreads.aligned.bam (and .bai/.pbi indexes): subreads of molecules in which ssDNA and dsDNA calls were made.
2. ssDNA/dsDNA.trinuc.counts.uncorrected.txt: ssDNA and dsDNA raw call counts for each trinucleotide context 
3. ssDNA/dsDNA.trinuc.counts.corrected.txt: ssDNA and dsDNA call counts for each trinucleotide context, corrected for the trinucleotide content of the genome (using the GRCh37 reference to allow comparison to COSMIC signatures) relative to the that of interrogated read bases.
4. sigfit.ssDNA/dsDNA.counts.uncorrected.pdf: ssDNA and dsDNA spectra of raw call counts
5. sigfit.ssDNA/dsDNA.probability.uncorrected.pdf: ssDNA and dsDNA spectra of fractional contribution of each trinucleotide context (i.e. raw counts / total counts)
6. sigfit.ssDNA/dsDNA.probability.corrected.pdf: ssDNA and dsDNA spectra of fractional contribution of each trinucleotide context, corrected for the trinucleotide content of the genome relative to the that of interrogated read bases.
7. sigfit.ssDNAcollapsed.counts.uncorrected.pdf: ssDNA raw call counts, collapsed to 96-trinucleotide context.
8. sigfit.ssDNAcollapsed.probability.uncorrected.pdf: ssDNA spectrum of fractional contribution of each trinucleotide context, collapsed to 96-trinucleotide context.
9. sigfit.ssDNAcollapsed.probability.corrected.pdf: ssDNA spectrum of fractional contribution of each trinucleotide context, collapsed to 96-trinucleotide context, corrected for the trinucleotide content of the genome relative to the that of interrogated read bases.
10. [samplename].[genomereference].[threshold-index].sigfit/sigfit_thresh.dsDNA/ssDNAcollapsed_Catalogues_[timestamp].pdf
11. [samplename].[genomereference].[threshold-index].sigfit/sigfit_thresh.dsDNA/ssDNAcollapsed_Signatures_[timestamp].pdf
12. [samplename].[genomereference].[threshold-index].sigfit/sigfit_thresh.dsDNA/ssDNAcollapsed_Exposures_[timestamp].pdf
13. [samplename].[genomereference].[threshold-index].sigfit/sigfit_thresh.dsDNA/ssDNAcollapsed_Reconstructions_[timestamp].pdf
    - sigfit: fitting to COSMIC SBS signatures
    - sigfit_thresh: re-fitting only to COSMIC SBS signatures with exposures whose lower Bayesian HPD interval is > 0.01, plus SBS1/5/18/40 that are common in healthy tissues.
    For each:
      - Catalogues: spectra of raw counts
      - Signatures: signatures used in fitting
      - Exposures: fraction of calls attributed to each signature
      - Reconstructions: sigfit reconstruction of spectra using component signatures
14. mutation_signatures.RDS: R data file with sigfit analysis results

## Citation
If you use HiDEF-seq, please cite:

Liu MH, Costa B, Choi U, Bandler RC, Lassen E, Grońska-Pęski M, Schwing A, Murphy ZR, Rosenkjær D, Picciotto S, Bianchi V, Stengs L, Edwards M, Loh CA, Truong TK, Brand RE, Pastinen R, Wagner JR, Skytte AB, Tabori U, Shoag JE, Evrony GD. Single-strand mismatch and damage patterns revealed by single-molecule DNA sequencing. bioRxiv 2023.02.19.526140; doi: https://doi.org/10.1101/2023.02.19.526140 (2023).
