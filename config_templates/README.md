# Configuration templates

Here, we describe the HiDEF-seq pipeline's parameters, which are configured in a [YAML-format](https://en.wikipedia.org/wiki/YAML) parameters file.

We provide example template YAML files in this directory: a generic template (`analysis.template.yaml`) and an hg38-focused template (`analysis.hg38.template.yaml`) that capture commonly used settings for the hg38 reference genome.

Below, the `parameter[].subparameter` notation indicates that `parameter` is a list that can contain multiple `subparameter` items.

## Table of contents
- [Global analysis identifier](#global-analysis-identifier)
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

## Global analysis identifier

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `analysis_id` | string | yes | Unique identifier prefixed onto every output path emitted by `main.nf` and the R summary scripts. |

## Sequencing runs and samples

### Run-level entries
| Key | Type | Required | Notes |
| --- | --- | --- | --- |
| `reads_type` | string (`subreads` or `ccs`) | yes | Type of reads: `subreads` data are split into `ccs_chunks` for CCS consensus calling, while `ccs` data skip CCS consensus calling. All `runs[].reads_file` entries must have the same `reads_type`. |
| `runs[].run_id` | string | yes | Run identifier propagated to data outputs. |
| `runs[].reads_file` | path | yes | Absolute path to the CCS or subread BAM. A companion BAM index file with suffix `.pbi` is expected. |
| `runs[].samples[]` | list | yes | Barcode definitions for each sample present in the run. |
| `runs[].samples[].sample_id` | string | yes | Sample identifier. Must match an entry in samples[].sample_id`. |
| `runs[].samples[].barcode_id` | string | yes | Barcode label used during BAM demultiplexing. |
| `runs[].samples[].barcode` | DNA sequence | yes | Barcode sequence. |

### Sample metadata
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `samples[].sample_id` | string | yes | Sample identifier. Referenced by `runs[].samples[].sample_id`  and propagated to data outputs. |
| `samples[].individual_id` | string | yes | Individual identifier. Must match an entry in `individuals[].individual_id`. |
| `samples[].tissue` | string | optional | Tissue type annotation propagated to run metadata table. |

## Individuals and germline resources

### Individual-level inputs
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `individuals[].individual_id` | string | yes | Individual identifier (i.e., identifier of the person/animal). Referenced by`samples[].individual_id` and propagated to data outputs. |
| `individuals[].sex` | string (`male` or `female`) | yes | Sex of the individual. Not currently used in analysis. |
| `individuals[].germline_bam_file` | path | yes | Germline sequencing aligned BAM. |
| `individuals[].germline_bam_type` | string (`Illumina` or `PacBio`) | yes | Type of germline sequencing. Affects coverage analysis in the `processGermlineBAMs` process. |
| `individuals[].germline_vcf_files[]` | list | optional | Germline VCF variant files for each individual. |
| `individuals[].germline_vcf_files[].germline_vcf_file` | path | yes (if list present) | Bgzipped VCF file. Must have a matching index. |
| `individuals[].germline_vcf_files[].germline_vcf_type` | string | yes (if list present) | Name of tool used to call the germline VCF file's variants. Must match an entry in  `germline_vcf_types[].germline_vcf_type`. |

## Output locations

### Paths
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `analysis_output_dir` | path | yes | Root directory for final output folders. |
| `cache_dir` | path | yes | Shared cache for reference-dependent intermediates (BSgenome installation, region filters, processed germline variant resources). |

## Container and tool paths

### Container and environment hooks
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `hidefseq_container` | string | yes | Docker image reference (`docker://…`) or Singularity `.sif` path used for all processes. |
| `conda_base_script` | path | defaults to `/hidef/miniconda3/etc/profile.d/conda.sh` | Script sourced prior to activating conda environment in the docker image. |
| `conda_pbbioconda_env` | path | defaults to `/hidef/bin/pbconda` | Conda environment housing PacBio command-line tools in the docker image. |
| `samtools_bin`, `bcftools_bin`, `bedGraphToBigWig_bin`, `wiggletools_bin`, `wigToBigWig_bin`, `seqkit_bin`, `bedtools_bin`, `bgzip_bin`, `tabix_bin` | string | yes | Absolute paths to binaries of tools inside the docker image. |
| `ccs_ld_preload` | path | optional | Shared library path exported via `LD_PRELOAD` before invoking the PacBio `ccs` binary, mitigating thread-affinity issues described in <a href="https://github.com/microsoft/onnxruntime/issues/10736" target="_blank" rel="noopener noreferrer">onnxruntime issue #10736</a>. Leave blank to disable. |

## Pipeline runtime parameters

### Nextflow chunking and retries
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `ccs_chunks` | integer | required when `reads_type: subreads` | Number of chunks in which CCS consensus sequence calling is performed for subreads data. Resulting CCS chunks are then merged. |
| `lima_min_score` | integer | yes | Minimum Lima barcode score enforced by `limaDemux`. See [lima documentation](https://lima.how/faq/filter-input.html#--min-score) for details. |
| `analysis_chunks` | integer | yes | Number of chunks to split processed BAM files into in the `splitBAMs` workflow for subsequent `extractCalls` and `filterCalls` workflows. Higher values increase parallelism. |
| `mem_extractCallsChunk`, `time_extractCallsChunk`, `maxRetries_extractCallsChunk` | string/integer | optional | Baseline memory, time limit, and retry count for `extractCallsChunk` processes. Each retry will increase the memory and time limit by one half of the baseline. |
| `mem_filterCallsChunk`, `time_filterCallsChunk`, `maxRetries_filterCallsChunk` | string/integer | optional | Analogous settings for `filterCallsChunk` processes. |
| `mem_calculateBurdensChromgroupFiltergroup`, `time_calculateBurdensChromgroupFiltergroup`, `maxRetries_calculateBurdensChromgroupFiltergroup` | string/integer | optional | Analogous settings for `calculateBurdensChromgroupFiltergroup` processes. |
| `remove_intermediate_files` | boolean, true or false | optional | When true, enables the `removeIntermediateFiles` workflow segment to delete intermediate per-sample directories once outputs are finalised. |

## Reference genome resources

### Required files
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `genome_fasta` | path | yes | Reference FASTA. |
| `genome_fai` | path | yes | FASTA index produced by `samtools faidx`. |
| `genome_mmi` | path | yes | pbmm2 index produced by `pbmm2 index`. |
| `BSgenome.BSgenome_name` | string | yes | BSgenome package name (for example `BSgenome.Hsapiens.UCSC.hg38`) used by R scripts. Run `BSGenome::available.genomes()` in R to see all publicly available genomes. If the desired reference genome is not available, [create a custom BSgenome package](https://bioconductor.org/packages/devel/bioc/manuals/BSgenomeForge/man/BSgenomeForge.pdf) and set `BSgenome.BSgenome_name` to the custom package's name. |
| `BSgenome.BSgenome_file` | path | optional | Tarball (.tar.gz) file containing a custom BSgenome build. When supplied, `installBSgenome.R` installs it into `cache_dir`. |
| `sex_chromosomes` | comma-separated string | yes | Names of sex chromosomes. Used to exclude them from sensitivity analysis. |
| `mitochondrial_chromosome` | string | yes | Name of the mitochondrial chromosome. Used to exclude it from sensitivity analysis. |

### Chromosome grouping
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `chromgroups[]` | list | yes | Declares groups of chromosomes to be analyzed together. All pipeline outputs are calculated per chromosome group. For example, you may wish to analyze chr1-22,chrX,chrY separately from chrM to obtain separate burdens and spectra for each. |
| `chromgroups[].chromgroup` | string | yes | Chromosome group name. |
| `chromgroups[].chroms` | comma-separated string | yes | Chromosomes to include in the group. Do not include a chromosome in more than one group. |

## Call extraction settings

### Overlap and call type structure
| Key | Type | Required | Description |
| --- | :-- | :-- | --- |
| `min_strand_overlap` | float (0–1) | yes | Minimum fraction of forward/reverse strand read overlap required for a molecule to be included in call extraction. Fraction of forward strand overlapping reverse read, and fraction of reverse read overlapping forward read must pass this filter. |
| `call_types[]` | list | yes | Settings for analysis of each type of call. |
| `call_types[].call_type` | string | yes | Name of the call type. Must be either `SBS` (single base substitution) when `call_class` = `SBS`, `insertion` or `deletion` when `call_class` = `indel`, or the name of an MDB (modified base; e.g., 5-methylcytosine). |
| `call_types[].call_class` | string (`SBS`, `indel`, `MDB`) | yes | Class of the call type. |
| `call_types[].analyzein_chromgroups` | string | yes | Either `all` or a comma-separated subset of `chromgroups[].chromgroup` values indicating in which chromgroups to analyze the call type. |
| `call_types[].SBSindel_call_types[]` | list | yes | Call subtypes that specify different options for strand match and mismatches. |
| `call_types[].SBSindel_call_types[].SBSindel_call_type` | string | yes | Call subtype.<br /><br />For `call_class` = `SBS` or `indel`, accepted values are `mutation` (a reverse complement call was identified in the opposite strand), `mismatch-ss` ('single-strand mismatch': the call was identified in only a single strand and there is no call on the opposite strand, and `mismatch-ds` ('double-strand mismatch': the call was identified in only a single strand and there is a non-reverse complement call on the opposite strand).<br /><br />For `call_class`= `MDB`, accepted values are `mutation` (the MDB is at the position of an SBS or indel mutation), `mismatch-ss` (the MDB is at a position and strand with an SBS or indel mismatch-ss call), `mismatch-ds` (the MDB is at the position of an SBS or indel mismatch-ds call), `mismatch-os` (the MDB is at a position and strand whose opposite strand has an SBS or indel mismatch-ss call), and `match` (the MDB is at a position with no SBS or indel call on either strand). |
| `call_types[].SBSindel_call_types[].filtergroup` | string | yes | Filtergroup whose filter settings will be used for this call type. Must match an entry in `filtergroups[].filtergroup`. |

### MDB-specific options
| Key | Type | Required | Description |
| --- | :-- | --- | --- |
| `call_types[].MDB_bamtag` | string | set only for MDB call_class entries | BAM auxiliary tag storing MDB scores. |
| `call_types[].MDB_min_score` | numeric | set only for MDB call_class entries | Minimum score filter for calls of this MDB type. |
| `call_types[].MDB_sensitivity` | float (0-1) | set only for MDB call_class entries | Sensitivity for calls of this MDB type. Used to calculate sensitivity-corrected burdens. |
| `call_types[].MDB_base_opposite_strand` | string or `null` | optional | When set, requires the specified base on the opposite strand for MDB calls; `null` disables this filter. |

## Germline VCF filter definitions

### Threshold fields
| Key | Type | Required | Description |
| --- | :-- | :-- | --- |
| `germline_vcf_types[]` | list | yes | Settings for filters that determine which germline variants from each germline VCF type are used for germline variant filtering. |
| `germline_vcf_types[].germline_vcf_type` | string | yes | Name of tool used to call the germline VCF's variants.  Referenced by`individuals[].germline_vcf_files[].germline_vcf_type`. |
| `germline_vcf_types[].SBS_FILTERS[]` | list of strings | yes | List of VCF `FILTER` column values, at least one of which must be present in order to include an SBS variant. Important: surround each entry with quotes. |
| `germline_vcf_types[]`.`SBS_min_Depth`, `SBS_min_VAF`, `SBS_min_GQ`, `SBS_min_QUAL` | numeric | yes | Minimum total read depth (for all alleles at the site), allele fraction of the variant, genotype quality of the variant, and VCF QUAL column value of the variant in order to include an SBS variant. |
| `germline_vcf_types[].indel_FILTERS[]` | list of strings | yes | List of VCF `FILTER` column values, at least one of which must be present in order to include an indel variant. Important: surround each entry with quotes. |
| `germline_vcf_types[]`.`indel_min_Depth`, `indel_min_VAF`, `indel_min_GQ`, `indel_min_QUAL` | numeric | yes | Minimum thresholds in order to include an indel variant. |
| `germline_vcf_types[]`.`indel_inspad`, `indel_delpad` | string | yes; set to `NA` or `m0b0` to disable | Specification of padding to add around germline insertion (`indel_inspad`) and deletion (`indel_delpad`) variants for filtering HiDEF-seq calls. Specified as `m<multiplier>b<offset>` (for example `m2b15`). `m` multiplies the insertion or deletion length, `b` adds a fixed base count, and the pipeline adds flanking bases on both sides of each variant's span where each flank size is the larger of these two numbers. |

## Filter group definitions

### Basic molecule- and call-level filters

Filter groups enable creation of different sets of filter thresholds, each of which can be applied to multiple call types.

| Key | Type | Required | Description |
| --- | :-- | --- | --- |
| `filtergroups[]` | list | yes | Settings for filters that are applied to molecules and calls for each call type assigned to the filtergroup. |
| `filtergroups[].filtergroup` | string | yes | Name referenced by `call_types[].SBSindel_call_types[].filtergroup`. |
| `filtergroups[]`.`min_rq_eachstrand`, `min_rq_avgstrands`<br />(molecule-level filter) | float (0-1), float (0-1) | yes | Minimum read-quality (calculated by `ccs` as the average of consensus base qualities; obtained from the BAM file's`rq` tag) thresholds that must pass for both strands (`_eachstrand`) and for the average across both strands (`_avgstrands`) to keep the molecule in the analysis. |
| `filtergroups[]`.`min_ec_eachstrand`, `min_ec_avgstrands`<br />(molecule-level filter) | numeric, numeric | yes | Minimum effective coverage  (calculated by `ccs` as the coverage of the consensus sequence by subreads; obtained from the BAM file's `ec` tag) thresholds that must pass for both strands and for the average across both strands. |
| `filtergroups[]`.`min_mapq_eachstrand`, `min_mapq_avgstrands`<br />(molecule-level filter) | numeric, numeric | yes | Minimum mapping quality thresholds that must pass for both strands and for the average across both strands. |
| `filtergroups[]`.`max_num_SBScalls_eachstrand`, `max_num_SBScalls_stranddiff`, `max_num_SBSmutations`<br />(molecule-level filter) | integer, integer, integer | yes | Maximum number of SBS calls per strand, allowable difference in the number of SBS calls between strands, and number of SBS mutations per molecule. |
| `filtergroups[]`.`max_num_indelcalls_eachstrand`, `max_num_indelcalls_stranddiff`, `max_num_indelmutations`<br />(molecule-level filter) | integer, integer, integer | yes | Analogous filters for indels. |
| `filtergroups[]`.`max_num_softclipbases_eachstrand`, `max_num_softclipbases_avgstrands`<br />(molecule-level filter) | integer, numeric | yes | Maximum number of soft-clipped bases thresholds that must pass for both strands and for the average across both strands. |
| `filtergroups[]`.`max_num_SBScalls_postVCF_eachstrand`, `max_num_SBSmutations_postVCF`<br />(molecule-level filter) | integer, integer | yes | Maximum number of SBS calls per strand and number of SBS mutations per molecule after germline VCF filtering is applied. |
| `filtergroups[]`.`max_num_indelcalls_postVCF_eachstrand`, `max_num_indelmutations_postVCF`<br />(molecule-level filter) | integer, integer | yes | Analogous filters for indels. |
| `filtergroups[]`.`min_qual`, `min_qual_method`<br />(call-level filters) | numeric, string (`mean`, `all`, `any`) | yes | Minimum base quality score filter for a call (`min_qual`) that must pass in both strands, even for single-strand calls. The thresholding method (`min_qual_method`) used for indels with length > 1 can be `mean`, `all`, or `any` , indicating that either the mean of base qualities, all base qualities, or any base qualities of a strand are compared to `min_qual`. If any base quality in the opposite strand is NA, the call fails this filter. |
| `filtergroups[].read_trim_bp`<br />(molecule-region filter) | integer | yes | Number of bases trimmed from each read end. |
| `filtergroups[]`.`ccsindel_inspad`, `ccsindel_delpad`<br />(molecule-region filter) | string, string | yes; set to `NA` or `m0b0` to disable | Padding strings using the same `m<multiplier>b<offset>` syntax described above, to filter regions near HiDEF-seq (i.e. non-germline) indel calls in each strand. This filter is not applied to indel calls. |
| `filtergroups[].min_BAMTotalReads`<br />(genomic-region filter) | integer | yes | Minimum number of reads required in the germline sequencing data at the site of the call (to avoid false-positives due to false-negative detection of a germline variant). For insertions, left and right flanking bases, and for deletions, the deleted bases, in genome reference space must pass the filter. |
| `filtergroups[]`.`max_BAMVariantReads`, `max_BAMVAF`<br />(call-level filters) | integer, float (0-1) | yes | Maximum number of germline sequencing variant reads and maximum VAF of germline sequencing variant reads that match the HiDEF-seq call as obtained by `bcftools mpileup` of the germline sequencing BAM file. Applied only to SBS calls. |
| `filtergroups[]`.`min_frac_subreads_cvg`, `min_num_subreads_match`, `min_frac_subreads_match`<br />(call-level filter) | float (0-1), integer, float (0-1) | yes | Subread coverage thresholds derived from `sa`, `sm`, and `sx` tags. |
| `filtergroups[].min_subreads_cvgmatch_method`<br />(call-level filter) | string (`mean`, `all`,  `any`) | yes | Aggregation method (`mean`, `all`, or `any`) applied when summarising subread coverage metrics across multi-base events. |
| `filtergroups[].max_finalcalls_eachstrand`<br />(molecule-level filter) | Integer | yes | Maximum allowed final calls per strand; molecules exceeding the cap are discarded. |

## Region filter configuration

### Read filters (`region_filters[].read_filters[]`)
| Key | Description |
| --- | --- |
| `region_filter_file` | bigWig track providing per-base scores. |
| `binsize` | Integer bin size used when rescaling the bigWig to genomic intervals. |
| `threshold` | Comparison encoded as `<operator><value>` using `lt` (<), `lte` (≤), `gt` (>), or `gte` (≥); for example `gte0.1`. |
| `padding` | Number of bases expanded around filtered intervals. |
| `read_threshold` | Optional per-read comparison using the same operator syntax as `threshold`. |
| `applyto_chromgroups` | `all` or a comma-separated list of chromgroups receiving the filter. |
| `applyto_filtergroups` | `all` or a comma-separated list of filter groups inheriting the filter. |
| `is_germline_filter` | Boolean indicating whether the filter targets germline contexts. Germline filters are excluded from sensitivity calculations. |

### Genome filters (`region_filters[].genome_filters[]`)
| Key | Description |
| --- | --- |
| `region_filter_file` | bigWig describing genome regions to mask during burden calculations. |
| `binsize`, `threshold`, `padding` | Same semantics as the read filters. |
| `applyto_chromgroups`, `applyto_filtergroups`, `is_germline_filter` | Behave as in the read filter configuration (without a `read_threshold`). |

## Sensitivity estimation

### Sensitivity parameters
| Key | Type | Description |
| --- | --- | --- |
| `sensitivity_parameters.use_chromgroup` | string | Chromgroup whose sensitivity estimates seed other groups. Leave blank to skip sensitivity correction (burdens default to 1.0). |
| `sensitivity_parameters.sensitivity_vcf` | path | External population reference germline VCF used to identify high-confidence heterozygous or homozygous variants for sensitivity analysis. |
| `sensitivity_parameters.genotype` | string (`heterozygous` or `homozygous`) | Determines which genotype states are tallied when measuring detection sensitivity. |
| `sensitivity_parameters.SBS_min_Depth_quantile`, `SBS_min_VAF`, `SBS_max_VAF`, `SBS_min_GQ_quantile`, `SBS_min_QUAL_quantile`, `SBS_min_variant_detections` | numeric | Thresholds applied to SBS variants; quantile parameters evaluate empirical distributions, whereas absolute thresholds apply directly. |
| `sensitivity_parameters.indel_min_Depth_quantile`, `indel_min_VAF`, `indel_max_VAF`, `indel_min_GQ_quantile`, `indel_min_QUAL_quantile`, `indel_min_variant_detections` | numeric | Analogous thresholds for indel sensitivity selection. |
| `sensitivity_parameters.default_sensitivity` | numeric | Optional fallback sensitivity applied when no qualifying variants are detected. |

MDB-specific sensitivity overrides can also be provided per call type through `call_types[].MDB_sensitivity`.
