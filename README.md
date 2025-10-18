# HiDEF-seq

## Intro
HiDEF-seq is a single-molecule sequencing method with single-molecule fidelity. This repository contains a pipeline for analysis of HiDEF-seq data.

The latest HiDEF-seq v3 library preparation protocol (with random fragmentation and whole-genome coverage; link coming soon) and this HiDEF-seq v3 analysis pipeline are compatible with both PacBio Revio instruments (ccs consensus sequence data) and Sequel II instruments (subread sequence data). See the HiDEF-seq library preparation protocol for important settings required for sequencing on PacBio Revio instruments, without which the data cannot be analyzed.

HiDEF-seq analysis also requires standard germline sequencing data for filtering germline variants, which should also be sequenced on PacBio Revio instruments (i.e., same platform as HiDEF-seq data) to minimize false-positive calls due to missed germline variants. Once PacBio sequencing costs drop further, the HiDEF-seq data could feasibly also be used as the germline sequencing data.

The HiDEF-seq analysis pipeline is designed to be run by Nextflow and a pre-configured docker image. Below are instructions on how to run the pipeline.

## Outline
- [Computing environment](#computing-environment)
- [Reference genome](#reference-genome)
- [Germline sequencing data processing](#germline-sequencing-data-processing)
- [Run HiDEF-seq analysis](#run-hidef-seq-analysis)
- [Outputs](#outputs)
- [Citation](#citation)



## Computing environment

### Overall structure of the analysis pipeline

HiDEF-seq is orchestrated with [Nextflow](https://www.nextflow.io/), which can execute pipelines on local workstations, high-performance computing clusters, and the cloud.

Specifically, we configured this GitHub repository to be able to be directly run by Nextflow via the [main.nf](main.nf) script. This main pipeline script includes a set of workflows that cover each step of the analysis. Each workflow is comprised of a set of processes that execute either bash or [R](https://www.r-project.org/) scripts that Nextflow can be configured to run within our pre-configured HiDEF-seq docker image.

### Create a Nextflow configuration file

Create a [Nextflow configuration file](https://www.nextflow.io/docs/latest/config.html) tailored to your computing environment. For example, running on a SLURM-based cluster requires defining a [SLURM executor](https://www.nextflow.io/docs/latest/executor.html#slurm) in the configuration profile that sets options such as queue names, maximum memory, and CPU resources. When using Singularity, enable it explicitly with `singularity.enabled = true`.

Pass this configuration file to Nextflow with the `-config` flag when running the pipeline as described in [Run HiDEF-seq analysis](#run-hidef-seq-analysis).

Consult the <a href="https://www.nextflow.io/docs/latest/index.html" target="_blank" rel="noopener noreferrer">Nextflow documentation</a> for more details of Nextflow's capabilities and options.

### HiDEF-seq pipeline docker image

We provide a fully configured docker image for the pipeline at `docker://gevrony/hidef-seq:3.0`.

When using singularity, download the docker image into an .sif file with `singularity pull docker://gevrony/hidef-seq:3.0`, and set `hidefseq_container` parameter to that.

This docker image can be utilized by the Nextflow pipeline by setting the `hidefseq_container` parameter in the pipeline's [YAML parameters file](#yaml-parameters-file) to point either to the docker link or to the singularity. sif file.

## Reference genome
The pipeline requires reference genome files and multiple derivative files, which can be prepared per below.

### Script requirements
- <a href="http://www.htslib.org/" target="_blank" rel="noopener noreferrer">samtools</a>
- <a href="https://github.com/PacificBiosciences/pbmm2" target="_blank" rel="noopener noreferrer">pbmm2</a> — you may use pbmm2 that is already installed inside the HiDEF-seq docker image. To call it inside the container, first activate the bundled environment with `source /hidef/miniconda3/etc/profile.d/conda.sh` followed by `conda activate /hidef/bin/pbconda`.
- <a href="https://bedtools.readthedocs.io/" target="_blank" rel="noopener noreferrer">bedtools</a>
- <a href="http://hgdownload.soe.ucsc.edu/admin/exe/" target="_blank" rel="noopener noreferrer">bedGraphToBigWig</a>
- <a href="http://www.htslib.org/" target="_blank" rel="noopener noreferrer">bcftools</a>

### Preparing reference genome files
1. Download a single FASTA file containing sequences of all contigs for the reference genome of interest.

2. Create a FASTA index: `samtools faidx genome.fasta`

3. Create a pbmm2 index: `pbmm2 index genome.fasta genome.mmi --preset CCS` for Revio data. Use `--preset SUBREAD` for Sequel II subread data.

4. Extract chromosome contig sizes in a tab-delimited file (columns: `chrom\tchrom_size`): `cut -f 1,2 genome.fa.fai > genome.chromsizes.tsv` 

5. Prepare bigWig format genomic filter tracks used by the pipeline:
   - Create BED files for every desired filter you will use for `read_filters` and `genome_filters` (see the [YAML configuration documentation](config_templates/README.md#region-filter-configuration) for details). For example, public resources such as the <a href="https://genome.ucsc.edu/" target="_blank" rel="noopener noreferrer">UCSC Genome Browser</a> provide centromere, telomere, and segmental duplication annotations.
   - Convert each BED file to a sorted, merged bedGraph and finally to bigWig format:
     ```
     bedtools sort -i filter.bed -g genome.fa.fai | \
       bedtools merge -i stdin | \
       awk '{print $0 "\t1"}' | \
       sort -k1,1 -k2,2n > filter.bedgraph
     
     bedGraphToBigWig filter.bedgraph genome.chrsizes.tsv filter.bw
     ```

6. Prepare gnomAD data (required for human analyses) from the <a href="https://gnomad.broadinstitute.org/downloads" target="_blank" rel="noopener noreferrer">gnomAD download portal</a>:

   - Download gnomAD genomes sites vcf files.

   - For each chromosome's vcf file, remove extraneous tags and normalize alleles:

      ```
      bcftools annotate -x ^INFO/AF,FORMAT gnomAD.[chrom].vcf.gz | \
        bcftools norm -Oz -m- > gnomAD.[chrom].norm.vcf.gz
      ```
      For mitochondrial data, use `-x ^INFO/AF_hom,INFO/AF_het,FORMAT`.

   - Concatenate all chromosomes' vcf files, filter by allele frequency, and index:

     ```
     bcftools concat -Oz `ls gnomAD.chrom*.norm.vcf.gz | sort -V | xargs` > gnomAD.norm.vcf.gz
     bcftools view -Oz -i '(INFO/AF[*] >= 0.001 | INFO/AF_hom[*] >= 0.001 | INFO/AF_het[*] >= 0.001) & FILTER="PASS"' \
       gnomAD.norm.vcf.gz > gnomAD.norm.AFfiltered.vcf.gz
     bcftools index gnomAD.norm.AFfiltered.vcf.gz
     ```

   - Split single-nucleotide variants (SNVs) and indels and convert each set to BED and bigWig tracks:
     ```
     bcftools view gnomAD.norm.AFfiltered.vcf.gz | \
       grep -v '^#' | \
       awk -v OFS='\t' '$5!="*" {if(length($4)==length($5)){print $1,$2-1,$2}}' \
       > gnomAD.norm.AFfiltered.vcf.gz.snvs.bed
     
     bcftools view gnomAD.norm.AFfiltered.vcf.gz | \
       grep -v '^#' | \
       awk -v OFS='\t' '$5!="*" {if(length($4) > length($5)){print $1,$2,$2+length($4)-1} else if(length($5)>length($4)){print $1,$2-1,$2+1}}' \
       > gnomAD.norm.AFfiltered.vcf.gz.indels.bed
     ```
     Then convert the resulting BED files to bigWig as described above.

## Germline sequencing data processing
Germline variants need to be called from germline sequencing data before starting HiDEF-seq analysis. These variant calls are used to filter germline variants during HiDEF-seq analysis. We recommend calling variants with more than one caller, to reduce the probability of missing germline variants that then would be mis-called as somatic events.

The repository provides an example workflow ([`scripts/Process_PacBio_GermlineWGS_for_HiDEF-seq_v3.sh`](scripts/Process_PacBio_GermlineWGS_for_HiDEF-seq_v3.sh)) that aligns PacBio HiFi reads with pbmm2 and runs both DeepVariant and Clair3 variant calling inside Singularity containers.

### Script requirements
- <a href="https://www.docker.com/" target="_blank" rel="noopener noreferrer">Docker</a> or <a href="https://sylabs.io/singularity/" target="_blank" rel="noopener noreferrer">Singularity</a> (to execute DeepVariant and Clair3)
- <a href="https://github.com/PacificBiosciences/pbmm2" target="_blank" rel="noopener noreferrer">pbmm2</a>
- <a href="https://github.com/google/deepvariant" target="_blank" rel="noopener noreferrer">DeepVariant</a>
- <a href="https://github.com/HKU-BAL/Clair3" target="_blank" rel="noopener noreferrer">Clair3</a>

### Processing outline
Run the script as follows:

```
Process_PacBio_GermlineWGS_for_HiDEF-seq_v3.sh [input_bam] [output_basename] [reference.fasta] [reference.mmi] [hidef-seq .sif path] [clair3 .sif path] [deepvariant .sif path]
```



## Run HiDEF-seq analysis

### Requirements
- HiDEF-seq docker image (either accessible via docker or downloaded as a singularity .sif image file)
- <a href="https://www.nextflow.io/" target="_blank" rel="noopener noreferrer">Nextflow</a> v25.04.3 or newer
- Nextflow configuration file ([described above](#create-a-nextflow-configuration-file))
- YAML parameters file (see below)

### YAML parameters file
All configuration parameters for the HiDEF-seq pipeline reside in a YAML-format file that enumerates samples, reference resources, filters, and per-workflow options.

You will need to prepare a YAML parameters file for each run of the pipeline. Template YAML parameters files and detailed documentation of parameters are available in [`config_templates`](config_templates).

### Run pipeline
The pipeline is run using Nextflow, which pulls from this GitHub repository all of the pipeline scripts.

Launch the pipeline as follows:

```
HIDEFSEQ_GITREPO=evronylab/hidef-seq
HIDEFSEQ_GITTAG=3.0
NEXTFLOW_CONFIG=/path/to/nextflow.config
YAML=/path/to/analysis.yaml
WORK_DIR=/path/to/nextflow_work

nextflow -config $NEXTFLOW_CONFIG \
  run $HIDEFSEQ_GITREPO \
  -r $HIDEFSEQ_GITTAG \
  -latest \
  -params-file $YAML \
  -resume \ #optional if resuming prior runs
  -work-dir "$WORK_DIR"
```

Refer to the <a href="https://www.nextflow.io/docs/latest/cli.html#run" target="_blank" rel="noopener noreferrer">Nextflow CLI documentation</a> for additional Nextflow runtime options.

Each pipeline invocation executes the full analysis sequence end-to-end.

HiDEF-seq records lightweight configuration digests for the R-driven stages so that `-resume` only re-runs the tasks affected by YAML edits. Changes to sections used by `installBSgenome`, `processGermlineVCFs`, `extractCallsChunk`, `filterCallsChunkChromgroupFiltergroup`, `calculateBurdensChromgroupFiltergroup`, or `outputResultsSample` mutate their command scripts (via embedded signatures), prompting Nextflow to invalidate the cache for just those processes while leaving unrelated stages untouched.

### Outputs
Pipeline outputs are organized per sample within `analysis_output_dir` configured in the YAML parameters file. A comprehensive description of every final product generated by the pipeline (i.e., outputs of `processReads`, `calculateBurdens`, and `outputResults` workflows) is available in [docs/outputs.md](docs/outputs.md).

## Citation
If you use HiDEF-seq, please cite:

> Liu MH *, Costa B *, Bianchini EC, Choi U, Bandler RC, Lassen E, Grońska-Pęski M, Schwing A, Murphy ZR, Rosenkjær D, Picciotto S, Bianchi V, Stengs L, Edwards M, Nunes NM, Loh CA, Truong TK, Brand RE, Pastinen T, Wagner JR, Skytte AB, Tabori U, Shoag JE, Evrony GD. Single-strand mismatch and damage patterns revealed by single-molecule DNA sequencing. *Nature* 630, 752–761 (2024). <a href="https://www.nature.com/articles/s41586-024-07532-8" target="_blank" rel="noopener noreferrer">https://www.nature.com/articles/s41586-024-07532-8</a>.
