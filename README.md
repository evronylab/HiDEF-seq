# HiDEF-seq

## Intro
HiDEF-seq is a single-molecule sequencing method with single-molecule fidelity. This repository contains a pipeline for analysis of HiDEF-seq data.

The latest HiDEF-seq v3 library preparation [protocol](https://dx.doi.org/10.17504/protocols.io.kxygxy9mwl8j/v3) (with random fragmentation and whole-genome coverage) and this HiDEF-seq v3 analysis pipeline are compatible with both PacBio Revio instruments (ccs consensus sequence data) and Sequel II instruments (subread sequence data). See the HiDEF-seq library preparation protocol for important settings required for sequencing on PacBio Revio instruments, without which the data cannot be analyzed.

HiDEF-seq analysis also requires standard germline sequencing data for filtering germline variants, which should also be sequenced on PacBio Revio instruments (i.e., same platform as HiDEF-seq data) to minimize false-positive calls due to missed germline variants. Once PacBio sequencing costs drop further, the HiDEF-seq data could feasibly also be used as the germline sequencing data. Note also that HiDEF-seq somatic insertion and deletion (indel) calls have not been validated to achieve single-molecule fidelity.

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

Specifically, we configured this GitHub repository to be able to be directly run by Nextflow via the [main.nf](main.nf) script. This main pipeline script includes a set of workflows that cover each step of the analysis. Each workflow is comprised of processes that execute shell commands, [R](https://www.r-project.org/) scripts, or the Python BAM splitter within the HiDEF-seq container.

### Create a Nextflow configuration file

Create a [Nextflow configuration file](https://www.nextflow.io/docs/latest/config.html) tailored to your computing environment. For example, running on a SLURM-based cluster requires defining a [SLURM executor](https://www.nextflow.io/docs/latest/executor.html#slurm) in the configuration profile that sets options such as queue names, maximum memory, and CPU resources. When using Singularity, enable it explicitly with `singularity.enabled = true`.

Pass this configuration file to Nextflow with the `-config` flag when running the pipeline as described in [Run HiDEF-seq analysis](#run-hidef-seq-analysis).

Consult the <a href="https://www.nextflow.io/docs/latest/index.html" target="_blank" rel="noopener noreferrer">Nextflow documentation</a> for more details of Nextflow's capabilities and options.

### HiDEF-seq pipeline docker image

We provide a fully configured docker image for the pipeline at `docker://gevrony/hidef-seq:3.0`.

When using singularity, download the docker image into an .sif file with `singularity pull docker://gevrony/hidef-seq:3.0`, and set `hidefseq_container` parameter to that.

This docker image can be utilized by the Nextflow pipeline by setting the `hidefseq_container` parameter in the pipeline's [YAML parameters file](#yaml-parameters-file) to point either to the docker link or to the singularity. sif file.

The `optimization` branch additionally requires `pysam` with libdeflate support
for its Python BAM splitter. The tested 3.0 container has system Python 3.12 at
`/usr/bin/python3`, but no system pip. Run the following as root while building
the writable container:

```bash
apt-get update
apt-get install --no-install-recommends python3-pip

env -u HTSLIB_LIBRARY_DIR -u HTSLIB_INCLUDE_DIR -u HTSLIB_LIBRARY_MODE \
  PYTHONNOUSERSITE=1 \
  CPPFLAGS=-I/hidef/miniconda3/include/python3.12 \
  HTSLIB_CONFIGURE_OPTIONS=--with-libdeflate \
  /usr/bin/python3 -m pip install \
    --no-cache-dir --no-deps --no-binary=pysam \
    --target /usr/local/lib/python3.12/dist-packages pysam==0.24.1

/usr/bin/python3 -c 'import pysam, pysam.config; print(pysam.__version__); assert pysam.config.HAVE_LIBDEFLATE == 1'
```

The existing Conda installation supplies only Python 3.12 header files for
compilation; no Conda packages are installed or changed. The source build uses
pysam's [bundled HTSlib](https://pysam.readthedocs.io/en/latest/installation.html)
and existing system compiler/compression libraries. Both tested upstream
0.23.3 and 0.24.1 Linux wheels lacked libdeflate, so the command builds pysam
from source. `--no-deps` avoids installing other runtime Python packages;
`--target` selects the existing system Python's local package directory.
There is no new environment or wrapper.

The 2026-10-07 apt simulation with current Ubuntu metadata showed that
`--no-install-recommends` adds pip, setuptools and wheel and updates only
`python3-pkg-resources`. Omitting that option also updates system Python and
native libraries used by existing pipeline tools. Package plans can change;
the observed minimal plan does not modify either Conda installation or the
pipeline's native libraries.

The native source build passed all splitter fixtures and dependency checks in
job `19321297`, using a project-local test target. All 13 extension libraries
resolved only package-local or existing system libraries under both normal
and PacBio-activated environments, without Conda runtime libraries. The test
used temporary pip 26.2.1; apt currently provides pip 24.0. The apt installation
was simulated, and the rebuilt SIF still requires inventory and full-workflow
validation. The earlier full LIB1 benchmark used pysam 0.23.3 in an isolated
test Conda installation; neither version is enforced by the splitter itself.

Save the rebuilt image under a new filename and set `hidefseq_container` to it.
The workflow uses `python3` from the container's PATH; no Python path entry is
needed in the run YAML. The pipeline does not install Python packages during
tasks.
Prepared-cache operations use R with `jsonlite`, `digest` and `openssl`, plus
Linux `flock`, `stat` and `sync`, already available in the tested 3.0 environment.

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
Process_PacBio_GermlineWGS_for_HiDEF-seq_v3.sh [input_bam] [output_basename] [reference.fasta] [reference.mmi] [hidef-seq .sif path] [clair3 .sif path] [deepvariant .sif path] [male/female] [PAR.bed]
```



## Run HiDEF-seq analysis

### Requirements
- HiDEF-seq docker image (either accessible via docker or downloaded as a singularity .sif image file)
- <a href="https://www.nextflow.io/" target="_blank" rel="noopener noreferrer">Nextflow</a> v26.04.0 or newer
- Nextflow configuration file ([described above](#create-a-nextflow-configuration-file))
- YAML parameters file (see below)

### YAML parameters file
All configuration parameters for the HiDEF-seq pipeline reside in a YAML-format file that enumerates samples, reference resources, filters, and per-workflow options.

You will need to prepare a YAML parameters file for each run of the pipeline. Template YAML parameters files and detailed documentation of parameters are available in [`config_templates`](config_templates).

### Run pipeline
The pipeline is run using Nextflow, which pulls from this GitHub repository all of the pipeline scripts.

Launch the pipeline as follows:

```bash
HIDEFSEQ_GITREPO=evronylab/HiDEF-seq
HIDEFSEQ_GITTAG=optimization
NEXTFLOW_CONFIG=/path/to/nextflow.config
YAML=/path/to/analysis.yaml
WORK_DIR=/path/to/nextflow_work

nextflow -config "$NEXTFLOW_CONFIG" \
  run "$HIDEFSEQ_GITREPO" \
  -r "$HIDEFSEQ_GITTAG" \
  -latest \
  -params-file "$YAML" \
  -resume \
  -work-dir "$WORK_DIR" \
  -with-report
```

`-r optimization -latest` fetches the current optimization branch for each test.
Omit `-resume` for a new execution; keep the same launch directory and work
directory when resuming. If several runs share a launch directory, use
`-resume <session-UUID>` to select the intended run. `-with-report` is optional.
See [Torch launch and measured resource settings](docs/workflow-optimization.md#torch-launch-and-resource-settings)
for the SLURM command, tested source revision, and workload-specific allocations.

Refer to the <a href="https://www.nextflow.io/docs/latest/cli.html#run" target="_blank" rel="noopener noreferrer">Nextflow CLI documentation</a> for additional Nextflow runtime options.

Each pipeline invocation executes the full analysis sequence end-to-end.

Prepared artifacts have separate cache identities based on their relevant inputs, settings, and code, allowing unaffected reference and germline preparation to be reused. Downstream R tasks receive an immutable effective YAML whose filename hashes its full contents. A change to any retained configuration field can therefore rerun downstream R tasks under `-resume`, even when their individual configuration signatures are unchanged. Downstream invalidation is currently conservative; it is not limited to the scientifically affected tasks. See [prepared caches and resume behavior](docs/workflow-optimization.md) for details.

These optimizations preserve scientific output schemas, formats, and filenames, including one final QS2 per sample. They introduce no call batching or output sharding. In the completed two-sample Torch workload, common downstream actual CPU fell from 447.2 to 230.2 hours (48.5%), and the highest burden-task Slurm RSS fell from 288.0 to 127.7 GiB (55.6%). All 711 required scientific output paths and 22 prepared-cache products passed validation under the agreed metadata, cross-chromgroup ordering and BAM coordinate-tie rules. The [validation ledger](docs/optimization-validation.md) records the full comparison scope, preparation costs and measurement limits; a [machine-readable result](docs/benchmarks/torch-2026-10-06.json) preserves evidence checksums.

### Outputs
A comprehensive description of every final product generated by the pipeline is available in the [outputs documentation](docs/outputs.md).

## Citation
If you use HiDEF-seq, please cite:

> Liu MH *, Costa B *, Bianchini EC, Choi U, Bandler RC, Lassen E, Grońska-Pęski M, Schwing A, Murphy ZR, Rosenkjær D, Picciotto S, Bianchi V, Stengs L, Edwards M, Nunes NM, Loh CA, Truong TK, Brand RE, Pastinen T, Wagner JR, Skytte AB, Tabori U, Shoag JE, Evrony GD. Single-strand mismatch and damage patterns revealed by single-molecule DNA sequencing. *Nature* 630, 752–761 (2024). <a href="https://www.nature.com/articles/s41586-024-07532-8" target="_blank" rel="noopener noreferrer">https://www.nature.com/articles/s41586-024-07532-8</a>.
