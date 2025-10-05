# YAML parameters

This document explains every configuration option consumed by the HiDEF-seq Nextflow pipeline (`main.nf`) and the downstream R scripts. Template files are available in [`config_templates/`](../config_templates) and mirror the structure described below.

## Table of contents
- [Global analysis identifiers](#global-analysis-identifiers)
- [Sequencing runs and samples](#sequencing-runs-and-samples)
- [Individuals and germline resources](#individuals-and-germline-resources)
- [Output locations](#output-locations)
- [Container and tool paths](#container-and-tool-paths)
- [Pipeline runtime parameters](#pipeline-runtime-parameters)
- [Reference genome resources](#reference-genome-resources)
- [Call extraction settings](#call-extraction-settings)
- [Germline VCF filter definitions](#germline-vcf-filter-definitions)
- [Filter group definitions](#filter-group-definitions)
- [Region filter configuration](#region-filter-configuration)
- [Sensitivity estimation](#sensitivity-estimation)

## Global analysis identifiers
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `analysis_id` | string | yes | Unique identifier used to prefix every output file and directory (`main.nf`, multiple R scripts). |
| `reads_type` | string (`subreads` or `ccs`) | yes | Declares the input data type for `processReads`. `subreads` triggers CCS generation (`ccsChunk` and `mergeCCSchunks` processes), whereas `ccs` skips CCS and starts from existing HiFi reads. All `runs[].reads_file` entries must match this type. |

## Sequencing runs and samples
The `runs` list describes each PacBio instrument run contributing HiDEF-seq reads.

| Key | Type | Required | Notes |
| --- | --- | --- | --- |
| `runs[].run_id` | string | yes | Run identifier used throughout logs and output names. |
| `runs[].reads_file` | path | yes | Absolute or relative path to the CCS or subread BAM file (`processReads`). A companion `.pbi` is expected at the same path. |
| `runs[].samples[]` | list | yes | Barcode definitions for each sample present in the run. |
| `runs[].samples[].sample_id` | string | yes | Sample identifier. Must match an entry in the top-level `samples` list. |
| `runs[].samples[].barcode_id` | string | yes | Lima barcode label used to name demultiplexed BAMs. |
| `runs[].samples[].barcode` | DNA sequence | yes | Forward barcode sequence written to a per-run FASTA (`makeBarcodesFasta`). |

Top-level sample metadata is defined in the `samples` list.

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `samples[].sample_id` | string | yes | Sample identifier referenced by `runs`, `individuals`, and downstream outputs. |
| `samples[].individual_id` | string | yes | Links each sample to a germline individual record. |
| `samples[].tissue` | string | optional | Free-form annotation stored with run metadata (`outputResults`). |

## Individuals and germline resources
Each entry ties a germline individual to BAM/VCF inputs used for filtering and sensitivity calculations.

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `individuals[].individual_id` | string | yes | Identifier referenced by `samples[].individual_id`. |
| `individuals[].sex` | string (`male` or `female`) | yes | Used when summarising burdens and sensitivity (sex chromosomes handled separately). |
| `individuals[].germline_bam_file` | path | yes | BAM aligned to the same reference; processed by `processGermlineBAMs` to produce coverage and VCF tracks. |
| `individuals[].germline_bam_type` | string (`Illumina` or `PacBio`) | yes | Chooses samtools/bcftools mpileup parameters in `processGermlineBAMs`. |
| `individuals[].germline_vcf_files[]` | list | optional | One or more VCFs describing germline variants. |
| `individuals[].germline_vcf_files[].germline_vcf_file` | path | yes (if list present) | VCF or BCF processed by `processGermlineVCFs`. Must be bgzipped and indexed. |
| `individuals[].germline_vcf_files[].germline_vcf_type` | string | yes (if list present) | Identifier matching an entry in `germline_vcf_types`, used to attach filter thresholds. |

## Output locations
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `analysis_output_dir` | path | yes | Root directory for per-sample results (all workflows). |
| `cache_dir` | path | yes | Shared cache for reference-dependent intermediates (BSgenome archives, region filters, germline BAM/VCF derivatives). Must be accessible inside the container and outside `/home` for Singularity bind mounts. |

## Container and tool paths
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `hidefseq_container` | string | yes | Docker image reference (`docker://...`) or Singularity `.sif` path used for every process. |
| `conda_base_script` | path | defaults to `/hidef/miniconda3/etc/profile.d/conda.sh` | Script sourced prior to activating the bundled conda environment in several processes. |
| `conda_pbbioconda_env` | path | defaults to `/hidef/bin/pbconda` | Conda environment that provides PacBio tools. |
| `samtools_bin`, `bcftools_bin`, `bedGraphToBigWig_bin`, `wiggletools_bin`, `wigToBigWig_bin`, `seqkit_bin`, `bedtools_bin`, `bgzip_bin`, `tabix_bin` | string | yes | Command names or absolute paths invoked inside the container. Override when custom installations are required. |
| `ccs_ld_preload` | path | optional | Library injected via `LD_PRELOAD` when running CCS (workaround for certain SLURM environments). Leave empty to disable. |

## Pipeline runtime parameters
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `ccs_chunks` | integer | required when `reads_type: subreads` | Number of CCS chunks processed per run (`ccsChunk`). |
| `lima_min_score` | integer | optional | Minimum score passed to Lima for barcode demultiplexing (`limaDemux`). |
| `analysis_chunks` | integer | yes | Number of BAM chunks produced per sample for downstream processing (`splitBAMs`, `extractCalls`, `filterCalls`). |
| `mem_extractCallsChunk`, `time_extractCallsChunk`, `maxRetries_extractCallsChunk` | string/integer | optional | Baseline memory, wall-clock time, and retry count for `extractCallsChunk`. Each retry increases memory/time as coded in `main.nf`. |
| `mem_filterCallsChunk`, `time_filterCallsChunk`, `maxRetries_filterCallsChunk` | string/integer | optional | Analogous settings for `filterCallsChunk`. |
| `mem_calculateBurdensChromgroupFiltergroup`, `time_calculateBurdensChromgroupFiltergroup`, `maxRetries_calculateBurdensChromgroupFiltergroup` | string/integer | optional | Resource controls for `calculateBurdensChromgroupFiltergroup`. |
| `remove_intermediate_files` | boolean | optional | When true and `--workflow` includes `removeIntermediateFiles`, purges intermediate per-sample directories after a successful run. |

Additional Nextflow parameters:

- `workflow` (CLI `--workflow`) selects pipeline stages (see [README](../README.md#run-hidef-seq-pipeline)).
- `paramsFileName` is captured automatically by `main.nf` and should not be set manually.

## Reference genome resources
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `genome_fasta` | path | yes | Reference FASTA used by pbmm2 and downstream scripts. |
| `genome_fai` | path | yes | FASTA index (output of `samtools faidx`). |
| `genome_mmi` | path | yes | pbmm2 minimap2 index (`pbmm2 index`). |
| `BSgenome.BSgenome_name` | string | yes | BSgenome package name (e.g., `BSgenome.Hsapiens.UCSC.hg38`) loaded by the R scripts. |
| `BSgenome.BSgenome_file` | path | optional | Tarball containing a custom BSgenome build; if supplied, `installBSgenome.R` installs it into `cache_dir`. |
| `sex_chromosomes` | comma-separated string | yes | Chromosome names treated as sex chromosomes (excluded from sensitivity calculations). |
| `mitochondrial_chromosome` | string | yes | Chromosome identifier for mitochondrial DNA. |
| `chromgroups[]` | list | yes | Groups of chromosomes analysed together. |
| `chromgroups[].chromgroup` | string | yes | Group name (used in output directories and reporting). |
| `chromgroups[].chroms` | comma-separated string | yes | Chromosome names included in the group. |

## Call extraction settings
These parameters control `extractCalls.R`, `filterCalls.R`, and `calculateBurdens.R`.

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `min_strand_overlap` | float (0–1) | yes | Minimum fraction of read overlap between forward and reverse strands required before extracting calls. |
| `call_types[]` | list | yes | Defines the call classes to process. Each entry contains: |
| `call_types[].call_type` | string | yes | Logical name for the call (e.g., `SBS`, `insertion`, custom MDB labels). |
| `call_types[].call_class` | string (`SBS`, `indel`, `MDB`) | yes | Governs which R code path handles the call. |
| `call_types[].analyzein_chromgroups` | string | yes | Either `all` or a comma-separated subset of `chromgroups[].chromgroup`. |
| `call_types[].SBSindel_call_types[]` | list | yes | One or more call subtypes per `call_type`. |
| `call_types[].SBSindel_call_types[].SBSindel_call_type` | string | yes | Subclass analysed (`mutation`, `mismatch-ss`, `mismatch-ds`, `mismatch-os`, `match`). |
| `call_types[].SBSindel_call_types[].filtergroup` | string | yes | Name of a `filtergroups[].filtergroup` used for this subtype. |
| MDB-specific keys | | | |
| `call_types[].MDB_bamtag` | string | required for MDB entries | BAM auxiliary tag containing MDB scores. |
| `call_types[].MDB_min_score` | numeric | required for MDB entries | Minimum per-read score to keep. |
| `call_types[].MDB_sensitivity` | float | optional | Sensitivity used to correct burdens for MDB call types. |
| `call_types[].MDB_base_opposite_strand` | string or `null` | optional | If set, requires the specified base on the opposite strand when evaluating MDBs. |

## Germline VCF filter definitions
Each entry describes thresholds applied to germline VCFs. Names must match `individuals[].germline_vcf_files[].germline_vcf_type`.

| Key | Type | Description |
| --- | --- | --- |
| `germline_vcf_type` | string | Identifier referenced by individuals. |
| `SBS_FILTERS[]` | list of strings | Filter names that, if present in the VCF `FILTER` column, cause the variant to be removed. Use quotes to preserve `.`. |
| `SBS_min_Depth`, `SBS_min_VAF`, `SBS_min_GQ`, `SBS_min_QUAL` | numeric | Minimum depth, variant allele fraction, genotype quality, and QUAL for SBS calls. |
| `indel_FILTERS[]` | list of strings | Filters to exclude for indels. |
| `indel_min_Depth`, `indel_min_VAF`, `indel_min_GQ`, `indel_min_QUAL` | numeric | Minimum thresholds for indels. |
| `indel_inspad`, `indel_delpad` | string | Indel padding settings passed to `bcftools` (e.g., `m2b15`). Set to `NA` to disable. |

## Filter group definitions
Filter groups are referenced by call types and specify per-molecule quality and coverage filters. Each group is a map containing:

| Key | Description |
| --- | --- |
| `filtergroup` | Name used by `call_types[].SBSindel_call_types[].filtergroup`. |
| `min_rq_eachstrand`, `min_rq_avgstrands` | Minimum read quality (rq tag) per strand and on average. |
| `min_ec_eachstrand`, `min_ec_avgstrands` | Minimum number of polymerase passes (ec tag). |
| `min_mapq_eachstrand`, `min_mapq_avgstrands` | Minimum mapping quality thresholds. |
| `max_num_SBScalls_eachstrand`, `max_num_SBScalls_stranddiff` | Maximum number of SBS calls allowed per strand and maximum imbalance between strands. |
| `max_num_SBSmutations` | Maximum SBS mutations per molecule before filtering. |
| `max_num_indelcalls_eachstrand`, `max_num_indelcalls_stranddiff`, `max_num_indelmutations` | Analogous controls for indels. |
| `max_num_softclipbases_eachstrand`, `max_num_softclipbases_avgstrands` | Maximum soft-clipped bases per strand and averaged. |
| `max_num_SBScalls_postVCF_eachstrand`, `max_num_SBSmutations_postVCF` | Post-germline filtering caps. |
| `max_num_indelcalls_postVCF_eachstrand`, `max_num_indelmutations_postVCF` | Post-germline caps for indels. |
| `min_qual`, `min_qual_method` | Minimum HiFi quality score (`qual` tag) and aggregation method (`mean`, `all`, or `any`) for multi-base events. |
| `read_trim_bp` | Number of bases trimmed from each read end before evaluating filters. |
| `ccsindel_inspad`, `ccsindel_delpad` | Padding strings used when intersecting indels with CCS subreads. |
| `min_BAMTotalReads`, `max_BAMVariantReads`, `max_BAMVAF` | Depth-based filters computed from germline BAMs (applied to SBS or indel calls as coded). |
| `min_frac_subreads_cvg`, `min_num_subreads_match`, `min_frac_subreads_match` | Subread coverage thresholds derived from `sa/sm/sx` tags. |
| `min_subreads_cvgmatch_method` | Aggregation method for subread filters (`mean`, `all`, `any`). |
| `max_finalcalls_eachstrand` | Maximum number of final calls per strand after all filters; if exceeded the molecule is dropped. |

## Region filter configuration
Region filters remove reads or genomic loci prior to burden calculation.

### Read filters (`region_filters[].read_filters[]`)
| Key | Description |
| --- | --- |
| `region_filter_file` | bigWig containing per-base scores (must correspond to prepared resources). |
| `binsize` | Integer bin size used when rescaling the bigWig. |
| `threshold` | Comparison and value (e.g., `gte0.1`, `lt0.5`). |
| `padding` | Number of bases expanded around filtered intervals. |
| `read_threshold` | Threshold applied to per-read annotations when intersecting reads (same syntax as `threshold`). |
| `applyto_chromgroups` | `all` or comma-separated list of chromgroups where the filter applies. |
| `applyto_filtergroups` | `all` or comma-separated list of filter groups that inherit the filter. |
| `is_germline_filter` | Boolean indicating if the filter targets germline sites. Germline filters are excluded from sensitivity calculations. |

### Genome filters (`region_filters[].genome_filters[]`)
| Key | Description |
| --- | --- |
| `region_filter_file` | bigWig of genome regions to mask. |
| `binsize`, `threshold`, `padding` | As above. |
| `applyto_chromgroups`, `applyto_filtergroups`, `is_germline_filter` | As above (without `read_threshold`). |

## Sensitivity estimation
Sensitivity settings influence burden correction in `calculateBurdens.R` and `outputResults.R`.

| Key | Type | Description |
| --- | --- | --- |
| `sensitivity_parameters.use_chromgroup` | string | Chromgroup from which sensitivity is derived for all groups. Leave blank to skip sensitivity adjustments (burdens default to 1.0 correction). |
| `sensitivity_parameters.sensitivity_vcf` | path | High-confidence germline VCF used to measure detection sensitivity. |
| `sensitivity_parameters.genotype` | string (`heterozygous` or `homozygous`) | Determines which genotypes are considered when tallying detection rates. |
| `sensitivity_parameters.SBS_min_Depth_quantile`, `SBS_min_VAF`, `SBS_max_VAF`, `SBS_min_GQ_quantile`, `SBS_min_QUAL_quantile`, `SBS_min_variant_detections` | numeric | Thresholds used to select SBS variants for sensitivity estimation. Quantile parameters are applied to the distribution of observed values. |
| `sensitivity_parameters.indel_min_Depth_quantile`, `indel_min_VAF`, `indel_max_VAF`, `indel_min_GQ_quantile`, `indel_min_QUAL_quantile`, `indel_min_variant_detections` | numeric | Corresponding thresholds for indel variants. |

MDB-specific sensitivity corrections are configured per call type using `call_types[].MDB_sensitivity`.
