# HiDEF-seq

## Intro
HiDEF-seq is a single-molecule sequencing method with single-molecule accuracy for single base substitutions, when present in either one or both strands of DNA. This repository contains a pipeline for analysis of HiDEF-seq data.

The latest HiDEF-seq library preparation protocol (v3: with random fragmentation and whole-genome coverage; link coming soon) and the this HiDEF-seq analysis pipeline (v3) are compatible with PacBio Revio. Note: the HiDEF-seq library preparation protocol contains important required settings for sequencing on PacBio Revio instruments, without which the data cannot be analyzed. HiDEF-seq analysis also requires standard germline sequencing data for filtering germline variants, which can also be sequenced on PacBio Revio.

This repository is configured so it can be directly run by [Nextflow](https://www.nextflow.io/), which dispatches and manages jobs flexibly in diverse computing environments. We also provide a pre-configured docker image for the pipeline.

## Outline
- [Computing environment](#computing-environment)
- [Reference genome](#reference-genome)
- [Germline sequencing data processing](#germline-sequencing-data-processing)
- [Run HiDEF-seq pipeline](#run-hidef-seq-pipeline)
- [Outputs](#outputs)
- [Citation](#citation)

## Computing environment
HiDEF-seq is orchestrated with [Nextflow](https://www.nextflow.io/), which can execute pipelines on local workstations, high-performance computing clusters, and the cloud. Consult the <a href="https://www.nextflow.io/docs/latest/index.html" target="_blank" rel="noopener noreferrer">Nextflow documentation</a> for an overview of supported executors and configuration patterns.

### Create a Nextflow configuration file

Create a [Nextflow configuration file](https://www.nextflow.io/docs/latest/config.html) tailored to your computing environment. For example, running on a SLURM-based cluster requires defining a [SLURM executor](https://www.nextflow.io/docs/latest/executor.html#slurm) in the Nextflow configuration profile that sets options such as queue names, maximum memory, and CPU resources. When using Singularity, enable it explicitly with `singularity.enabled = true`.

Pass this Nextflow configuration file to Nextflow with the `-config` flag when running the pipeline as described in [Run HiDEF-seq pipeline](#run-hidef-seq-pipeline).

### HiDEF-seq pipeline docker image

We provide a fully configured docker image for the pipeline that can be found on [docker hub](https://hub.docker.com/r/gevrony/hidef-seq) (tag:3.0).

This image is automatically pull If using singularity

## Reference genome
Preparing the reference genome requires installing several command-line tools and generating multiple derivative files used throughout the workflow.

### Script requirements
- <a href="http://www.htslib.org/" target="_blank" rel="noopener noreferrer">samtools</a>
- <a href="https://github.com/PacificBiosciences/pbmm2" target="_blank" rel="noopener noreferrer">pbmm2</a> — you may use pbmm2 that is already installed inside the HiDEF-seq docker image. To call it inside the container, first activate the bundled environment with `source /hidef/miniconda3/etc/profile.d/conda.sh` followed by `conda activate /hidef/bin/pbconda`.
- <a href="https://bedtools.readthedocs.io/" target="_blank" rel="noopener noreferrer">bedtools</a>
- <a href="http://hgdownload.soe.ucsc.edu/admin/exe/" target="_blank" rel="noopener noreferrer">bedGraphToBigWig</a>
- <a href="http://www.htslib.org/" target="_blank" rel="noopener noreferrer">bcftools</a>

### Preparing reference genome files
1. Download the FASTA for the reference genome of interest.
2. Create the FASTA index by running `samtools faidx chm13.draft_v1.0.fasta`.
3. Build the minimap2 index used by pbmm2 (use `--preset SUBREAD` for Sequel II subread data) by running `pbmm2 index ref.fasta ref.mmi --preset CCS`.
4. Record chromosome sizes in a tab-delimited file named `genome.chrsizes.tsv` (columns: `chrom\tchrom_size`).
5. Prepare genomic filter tracks used by the pipeline:
   - Identify BED files for `read_filters` and `genome_filters` (see the [YAML configuration documentation](config_templates/README.md#region-filter-configuration) for details). Public resources such as the <a href="https://genome.ucsc.edu/" target="_blank" rel="noopener noreferrer">UCSC Genome Browser</a> provide centromere, telomere, and segmental duplication annotations.
   - Convert each BED file to a sorted, merged bedGraph and finally to bigWig format:
     ```
     bedtools sort -i filter.bed -g genome.fa.fai | \
       bedtools merge -i stdin | \
       awk '{print $0 "\t1"}' | \
       sort -k1,1 -k2,2n > filter.bedgraph
     
     bedGraphToBigWig filter.bedgraph genome.chrsizes.tsv filter.bw
     ```
6. Prepare gnomAD data (required for human analyses) from the <a href="https://gnomad.broadinstitute.org/downloads" target="_blank" rel="noopener noreferrer">gnomAD download portal</a>:
   ```
   bcftools annotate -x ^INFO/AF,FORMAT gnomAD.chrom.vcf.gz | \
     bcftools norm -Oz -m- > gnomAD.chrom.norm.vcf.gz
   ```
   For mitochondrial data, retain `AF_hom` and `AF_het` instead of `AF`. Concatenate all chromosomes, filter by allele frequency, and index:
   ```
   bcftools concat -Oz gnomAD.chrom.norm.vcf.gz > gnomAD.norm.vcf.gz
   bcftools view -Oz -i '(INFO/AF[*] >= 0.001 | INFO/AF_hom[*] >= 0.001 | INFO/AF_het[*] >= 0.001) & FILTER="PASS"' \
     gnomAD.norm.vcf.gz > gnomAD.norm.AFfiltered.vcf.gz
   bcftools index gnomAD.norm.AFfiltered.vcf.gz
   ```
   Split single-nucleotide variants (SNVs) and indels and convert each set to BED and bigWig tracks:
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
   Convert the resulting BED files to bigWig as in step 5.

## Germline sequencing data processing
Germline variant calling should be completed before starting HiDEF-seq analysis. The repository provides an example workflow ([`scripts/Process_PacBio_GermlineWGS_for_HiDEF-seq_v3.sh`](scripts/Process_PacBio_GermlineWGS_for_HiDEF-seq_v3.sh)) that aligns PacBio HiFi reads with <a href="https://github.com/PacificBiosciences/pbmm2" target="_blank" rel="noopener noreferrer">pbmm2</a> and runs both DeepVariant and Clair3 inside Singularity containers.

### Script requirements
- <a href="https://www.docker.com/" target="_blank" rel="noopener noreferrer">Docker</a> or <a href="https://sylabs.io/singularity/" target="_blank" rel="noopener noreferrer">Singularity</a> (to execute DeepVariant and Clair3)
- <a href="https://github.com/PacificBiosciences/pbmm2" target="_blank" rel="noopener noreferrer">pbmm2</a>
- <a href="https://github.com/google/deepvariant" target="_blank" rel="noopener noreferrer">DeepVariant</a>
- <a href="https://github.com/HKU-BAL/Clair3" target="_blank" rel="noopener noreferrer">Clair3</a>

### Processing outline
1. Align germline reads to the reference genome with pbmm2.
2. Call germline variants with DeepVariant and Clair3 using the [provided script](scripts/Process_PacBio_GermlineWGS_for_HiDEF-seq_v3.sh) as a template for batching jobs on your scheduler.

## Run HiDEF-seq pipeline

### Requirements
- HiDEF-seq container image. For Singularity run `singularity pull docker://gevrony/hidef-seq:3.0`.
- <a href="https://www.nextflow.io/" target="_blank" rel="noopener noreferrer">Nextflow</a> v25.04.3 or newer.

### YAML parameters file
All run-time configuration resides in a YAML file that enumerates samples, reference resources, filters, and per-workflow options. Template files and detailed documentation are available in [`config_templates/`](config_templates) and elaborated in [YAML parameters](config_templates/README.md).

### Run pipeline
The pipeline is run using nextflow, whichi pulls from this GitHub repository all of the pipeline scripts (main.nf nextflow pipeline script and R scripts in the bin folder) from GitHub. Set environment variables that describe the repository revision, Nextflow environment configuration, and run-specific inputs, then launch the workflow:

```
HIDEFSEQ_GITREPO=evronylab/HiDEF-seq
HIDEFSEQ_GITTAG=3.0  # Update to the desired release tag
NEXTFLOW_CONFIG=/path/to/nextflow.config
YAML=/path/to/analysis.yaml
WORK_DIR=/path/to/nextflow_work

nextflow -config "$NEXTFLOW_CONFIG" \
  run "$HIDEFSEQ_GITREPO" \
  -r "$HIDEFSEQ_GITTAG" \
  -latest \
  -params-file "$YAML" \
  --workflow all \
  -resume \
  -work-dir "$WORK_DIR"
```

Refer to the <a href="https://www.nextflow.io/docs/latest/cli.html#run" target="_blank" rel="noopener noreferrer">Nextflow CLI documentation</a> for additional runtime options.

The `--workflow` parameter controls which segments of the pipeline execute. Supported values are:
- `all` — runs the entire pipeline.
- Ordered, contiguous subsets drawn from `processReads`, `splitBAMs`, `prepareFilters`, `extractCalls`, `filterCalls`, `calculateBurdens`, `outputResults`, `removeIntermediateFiles`.

Rules enforced by `main.nf`:
- `all` cannot be combined with other entries.
- Individual workflow names must be listed in the canonical order above.
- When running a subset, the steps must be contiguous (for example `processReads,splitBAMs,prepareFilters` is valid; `processReads,extractCalls` is rejected because `splitBAMs` is skipped).
- `removeIntermediateFiles` executes only if `remove_intermediate_files: true` in the YAML configuration.

### Outputs
Pipeline outputs are organized per sample within `analysis_output_dir`. A comprehensive description of every final product generated by the `processReads` and `outputResults` workflows is available in [docs/outputs.md](docs/outputs.md).

## Citation
If you use HiDEF-seq, please cite:

> Liu MH *, Costa B *, Bianchini EC, Choi U, Bandler RC, Lassen E, Grońska-Pęski M, Schwing A, Murphy ZR, Rosenkjær D, Picciotto S, Bianchi V, Stengs L, Edwards M, Nunes NM, Loh CA, Truong TK, Brand RE, Pastinen R, Wagner JR, Skytte AB, Tabori U, Shoag JE, Evrony GD. Single-strand mismatch and damage patterns revealed by single-molecule DNA sequencing. *Nature* (2024). <a href="https://www.nature.com/articles/s41586-024-07532-8" target="_blank" rel="noopener noreferrer">https://www.nature.com/articles/s41586-024-07532-8</a>.
