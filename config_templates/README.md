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
- [Germline VCF filter settings](#germline-vcf-filter-settings)
- [Filter group settings](#filter-group-settings)
- [Region filter configuration](#region-filter-configuration)
- [Sensitivity estimation settings](#sensitivity-estimation-settings)

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

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `analysis_output_dir` | path | yes | Root directory for final output folders. |
| `cache_dir` | path | yes | Shared cache for reference-dependent intermediates (BSgenome installation, region filters, processed germline variant resources). Populate a shared cache with a single Nextflow run at a time (or point concurrent runs to distinct caches) so that parallel executions do not overwrite the same files while they are being written. |

## Container and tool paths

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `hidefseq_container` | string | yes | Docker image reference (`docker://…`) or Singularity `.sif` path used for all processes. |
| `conda_base_script` | path | defaults to `/hidef/miniconda3/etc/profile.d/conda.sh` | Script sourced prior to activating conda environment in the docker image. |
| `conda_pbbioconda_env` | path | defaults to `/hidef/bin/pbconda` | Conda environment housing PacBio command-line tools in the docker image. |
| `samtools_bin`, `bcftools_bin`, `bedGraphToBigWig_bin`, `wiggletools_bin`, `wigToBigWig_bin`, `seqkit_bin`, `bedtools_bin`, `bgzip_bin`, `tabix_bin` | string | yes | Absolute paths to binaries of tools inside the docker image. |
| `ccs_ld_preload` | path | optional | Shared library path exported via `LD_PRELOAD` before invoking the PacBio `ccs` binary, mitigating thread-affinity issues described in <a href="https://github.com/microsoft/onnxruntime/issues/10736" target="_blank" rel="noopener noreferrer">onnxruntime issue #10736</a>. Leave blank to disable. |

## Pipeline runtime parameters

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `ccs_chunks` | integer | required when `reads_type: subreads` | Number of chunks in which CCS consensus sequence calling is performed for subreads data. Resulting CCS chunks are then merged. |
| `lima_min_score` | integer | yes | Minimum Lima barcode score enforced by `limaDemux`. See [lima documentation](https://lima.how/faq/filter-input.html#--min-score) for details. |
| `analysis_chunks` | integer | yes | Number of chunks to split processed BAM files into in the `splitBAMs` workflow for subsequent `extractCalls` and `filterCalls` workflows. Higher values increase parallelism. |
| `mem_extractCallsChunk`, `time_extractCallsChunk`, `maxRetries_extractCallsChunk` | string/integer | optional | Baseline memory, time limit, and retry count for `extractCallsChunk` processes. Each retry will increase the memory and time limit by one half of the baseline. |
| `mem_filterCallsChunkChromgroupFiltergroup`, `time_filterCallsChunkChromgroupFiltergroup`, `maxRetries_filterCallsChunkChromgroupFiltergroup` | string/integer | optional | Analogous settings for `filterCallsChunkChromgroupFiltergroup` processes. |
| `mem_calculateBurdensChromgroupFiltergroup`, `time_calculateBurdensChromgroupFiltergroup`, `maxRetries_calculateBurdensChromgroupFiltergroup` | string/integer | optional | Analogous settings for `calculateBurdensChromgroupFiltergroup` processes. |
| `mem_outputResultsSample`, `time_outputResultsSample`, `maxRetries_outputResultsSample` | string/integer | optional | Analogous settings for `outputResultsSample` processes. |
| `output_intermediate_files` | boolean, `true` or `false` | optional | When true, publishes intermediate per-sample outputs (for example `processedReads`, `splitBAMs`, `extractCalls`, and `calculateBurdens`). When false, those directories are omitted from `analysis_output_dir`. |

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

### Read overlap filter and call type structure
| Key | Type | Required | Description |
| --- | :-- | :-- | --- |
| `min_strand_overlap` | float (0–1) | yes | Minimum fraction of forward/reverse strand read overlap required for a molecule to be included in call extraction. Fraction of forward strand read overlapping reverse strand read, and fraction of reverse strand read overlapping forward strand read must both pass this filter. |
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

## Germline VCF filter settings

| Key | Type | Required | Description |
| --- | :-- | :-- | --- |
| `germline_vcf_types[]` | list | yes | Settings for filters that determine which germline variants from each germline VCF type are used for germline variant filtering. |
| `germline_vcf_types[].germline_vcf_type` | string | yes | Name of tool used to call the germline VCF's variants.  Referenced by`individuals[].germline_vcf_files[].germline_vcf_type`. |
| `germline_vcf_types[].SBS_FILTERS[]` | list of strings | yes | List of VCF `FILTER` column values, at least one of which must be present in order to include an SBS variant. Important: surround each entry with quotes. |
| `germline_vcf_types[]`.`SBS_min_Depth`, `SBS_min_VAF`, `SBS_min_GQ`, `SBS_min_QUAL` | numeric | yes | Minimum total read depth (for all alleles at the site), allele fraction, genotype quality, and VCF QUAL column value in order to include an SBS variant. |
| `germline_vcf_types[].indel_FILTERS[]` | list of strings | yes | List of VCF `FILTER` column values, at least one of which must be present in order to include an indel variant. Important: surround each entry with quotes. |
| `germline_vcf_types[]`.`indel_min_Depth`, `indel_min_VAF`, `indel_min_GQ`, `indel_min_QUAL` | numeric | yes | Minimum thresholds in order to include an indel variant. |
| `germline_vcf_types[]`.`indel_ins_pad`, `indel_del_pad` | string | optional | Specification of padding to add around germline insertion (`indel_ins_pad`) and deletion (`indel_del_pad`) variants for filtering HiDEF-seq calls. Specified as `m<multiplier>b<offset>` (for example `m2b15`). `m` multiplies the indel length, `b` is a fixed number of bases for all indels, and the pipeline filters flanking bases on each side of each variant where the size of the flank on each side is the larger of these two numbers calculated for each variant. Leave `indel_ins_pad` and `indel_del_pad` blank to disable this filter (including disabling of filtering at the site of the variant itself without padding). |

## Filter group settings

Filter groups enable creation of different sets of thresholds for basic molecule- and call-level filters. Each filter group can be applied to multiple call types per the `call_types[].SBSindel_call_types[].filtergroup` setting.

| Key | Type | Required | Description |
| --- | :-- | --- | --- |
| `filtergroups[]` | list | yes | Settings for filters that are applied to molecules and calls for each call type assigned to the filtergroup. |
| `filtergroups[].filtergroup` | string | yes | Name referenced by `call_types[].SBSindel_call_types[].filtergroup`. |
| `filtergroups[]`.`min_rq_eachstrand`, `min_rq_avgstrands`<br />(molecule-level filter) | float (0-1), float (0-1) | yes | Minimum read-quality (calculated by `ccs` as the average of consensus base qualities; obtained from the BAM file's`rq` tag) thresholds that must pass for both strands (`_eachstrand`) and for the average across both strands (`_avgstrands`) to keep the molecule in the analysis. |
| `filtergroups[]`.`min_ec_eachstrand`, `min_ec_avgstrands`<br />(molecule-level filter) | numeric, numeric | yes | Minimum effective coverage  (calculated by `ccs` as the coverage of the consensus sequence by subreads; obtained from the BAM file's `ec` tag) thresholds that must pass for both strands and for the average across both strands. |
| `filtergroups[]`.`min_mapq_eachstrand`, `min_mapq_avgstrands`<br />(molecule-level filter) | numeric, numeric | yes | Minimum mapping quality thresholds that must pass for both strands and for the average across both strands. |
| `filtergroups[]`.`max_num_SBScalls_eachstrand`, `max_num_SBScalls_stranddiff`, `max_num_SBSmutations`<br />(molecule-level filter) | integer, integer, integer | yes | Maximum number of SBS calls per strand, maximum allowable difference in the number of SBS calls between strands, and maximum number of SBS mutations per molecule. SBS calls for these filters are counted irrespective of any other filters. |
| `filtergroups[]`.`max_num_indelcalls_eachstrand`, `max_num_indelcalls_stranddiff`, `max_num_indelmutations`<br />(molecule-level filter) | integer, integer, integer | yes | Analogous filters for indels. |
| `filtergroups[]`.`max_num_softclipbases_eachstrand`, `max_num_softclipbases_avgstrands`<br />(molecule-level filter) | integer, numeric | yes | Maximum number of soft-clipped bases thresholds that must pass for both strands and for the average across both strands. |
| `filtergroups[]`.`max_num_SBScalls_postVCF_eachstrand`, `max_num_SBSmutations_postVCF`<br />(molecule-level filter) | integer, integer | yes | Maximum number of SBS calls per strand and maximum number of SBS mutations per molecule when only germline VCF filtering is applied. |
| `filtergroups[]`.`max_num_indelcalls_postVCF_eachstrand`, `max_num_indelmutations_postVCF`<br />(molecule-level filter) | integer, integer | yes | Analogous filters for indels. |
| `filtergroups[]`.`min_qual`, `min_qual_method`<br />(call-level filters) | numeric, string (`mean`, `all`, `any`) | yes | Minimum base quality score filter for a call (`min_qual`) that must pass in both strands, even for single-strand calls. The thresholding method (`min_qual_method`) used for indels with length > 1 can be `mean`, `all`, or `any` , indicating that either the mean of base qualities, all base qualities, or any base qualities of a strand are compared to `min_qual`. If any base quality in the opposite strand is NA, the call fails this filter. |
| `filtergroups[].read_trim_bp`<br />(molecule-region filter) | integer | yes | Number of bases trimmed from each read end. |
| `filtergroups[]`.`ccs_sbs_flank`, `ccs_ins_pad`, `ccs_del_pad`<br />(molecule-region filter) | integer, string, string | optional | Specification of regions around each HiDEF-seq (i.e. non-germline) SBS (`ccs_sbs_flank`), insertion (`ccs_ins_pad`), and deletion (`ccs_del_pad`) call for filtering nearby HiDEF-seq calls. A region surrounding a call in one strand is filtered from both strands. `ccs_sbs_flank` is the number of bases in each immediately adjacent flank to filter (i.e. not including the site of the SBS itself). To disable `ccs_sbs_flank` filtering, leave it blank or set it to 0. `ccs_ins_pad` and `ccs_del_pad` include the site of the indel itself, but they are not applied to indel calls, and they are specified using the same `m<multiplier>b<offset>` padding syntax described above for `germline_vcf_types`. Note that `ccs_ins_pad` and `ccs_del_pad` filter out all SBS calls opposite from an indel within the indel span + padding range. To disable one or more of these filters, leave them blank, and note that setting them to `m0b0` still applies the filter for the span of the indel (i.e., in cases where there is an SBS call on the opposite strand of a deletion). |
| `filtergroups[].min_germlineBAM_TotalReads`<br />(genomic-region filter) | integer | yes | Minimum number of reads required in the germline sequencing data at the site of the call (to avoid false-positives due to false-negative detection of a germline variant) as obtained by `samtools mpileup`. For insertions, left and right flanking bases, and for deletions, the deleted bases, in genome reference space must pass the filter. |
| `filtergroups[]`.`max_germlineBAM_VariantReads`, `max_germlineBAM_VAF`<br />(call-level filters) | integer, float (0-1) | yes | Maximum allowed number of germline sequencing variant reads and maximum allowed VAF of germline sequencing variant reads that match the HiDEF-seq call as obtained by `bcftools mpileup` of the germline sequencing BAM file. Applied only to SBS calls. |
| `filtergroups[]`.`min_frac_subreads_cvg`, `min_num_subreads_match`, `min_frac_subreads_match`<br />(call-level filter) | float (0-1), integer, float (0-1) | yes | Minimum required subread coverage for the strand at the site of the call, subreads whose sequence matches the call, and fraction of subreads whose sequence matches the call (obtained from `sa`, `sm`, and `sx` BAM tags and calculated as `sa/max(sa)` for the strand, `sm`, and `sm/(sm+sx)`, respectively). Filters must pass in both strands for all call types. Note: `sm/(sm+sx)` may be slightly deflated in very rare molecules with > 255 subreads per strand, due to sm capping at 255 but sx continuing to increase. This would cause filtering out of these calls if sm/(sm+sx) falls below the threshold. |
| `filtergroups[].min_subreads_cvgmatch_method`<br />(call-level filter) | string (`mean`, `all`,  `any`) | yes | The thresholding method used for `min_frac_subreads_cvg`, `min_num_subreads_match`, and `min_frac_subreads_match` when applied to indels with length > 1. Can be `mean`, `all`, or `any` , indicating that either the mean of the values for each base, all the values of bases, or any values of bases of a strand are compared to the respective threshold. If any `sa`, `sm`, or `sx` value in the opposite strand is NA, the call fails the respective filter. |
| `filtergroups[].max_finalcalls_eachstrand`<br />(molecule-level filter) | Integer | yes | Maximum allowed number of final calls per strand, counting the number of calls separately for each `call_type` x `SBSindel_call_type` combination. If either strand fails the filter, the entire molecule is discarded from analysis of that specific `call_type`.`SBSindel_call_type`. |

## Region filter configuration

Processed region filter bigWig files are prepared by the `prepareFilters` workflow for analysis based on the below settings and saved in `cache_dir` for future analyses. 

### Region-based molecule filters

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `region_filters[].read_filters[]` | list | optional | Region-based molecule filters that are applied as follows: each molecule is assessed against the set of genomic regions defined by the region filter and the below `binsize`, `threshold`, and `padding` settings. If the fraction of the molecule covered by the region filter (average of each strand) passes the `read_threshold` setting, then the entire molecule is filtered out. Otherwise, the filter does not remove the molecule or any part of it. |
| `region_filters[].read_filters[].region_filter_file` | path | yes | bigWig file of the region filter providing per-base scores. |
| `region_filters[].read_filters[].binsize` | integer | yes | Size of bins (in bases) within each of which bigWig per-base scores are averaged before applying `threshold`. `binsize: 1` retains the original per-base scores. |
| `region_filters[].read_filters[].threshold` | string | yes | Threshold applied to the bigWig score to select regions that will comprise the set of filter regions to consider in the subsequent `read_threshold` setting. Encoded as `<operator><value>` using `lt` (<), `lte` (≤), `gt` (>), or `gte` (≥); for example `gte0.1`. |
| `region_filters[].read_filters[].padding` | integer | yes | Number of bases to add to each flank of each region of the region filter. |
| `region_filters[].read_filters[].read_threshold` | string | yes | Fraction of the molecule that must be covered by the region filter (average of both strands) in order to be filtered. Encoded as for `threshold`. Encoded as `<operator><value>` as in `threshold`, but with `<value>` limited to a range of 0-1. |
| `region_filters[].read_filters[].applyto_chromgroups` | string | yes | `all` or a comma-separated list of chromgroups to which the filter is applied. |
| `region_filters[].read_filters[].applyto_filtergroups` | string | yes | `all` or a comma-separated list of filter groups to which the filter is applied. |
| `region_filters[].read_filters[].is_germline_filter` | boolean, `true` or `false` | yes | Boolean indicating whether the filter targets germline variant contexts (for example, a filter based on gnomAD variants). If `true`, the filter is excluded from sensitivity calculations that rely on germline variant calls that pass non-germline filters. |

### Genome region filters 
| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `region_filters[].genome_filters[]` | list | optional | Filters based on genomic regions that filter calls. Deletions are filtered out even if they partially overlap filtered regions. |
| `region_filters[].genome_filters[].region_filter_file` | path | yes | bigWig file of the region filter providing per-base scores. |
| `region_filters[].genome_filters[]`.`binsize`, `threshold`, `padding` | integer, string, integer | yes | Same formats as in `region_filters[].read_filters[]`. |
| `region_filters[].genome_filters[]`.`applyto_chromgroups`, `applyto_filtergroups`, `is_germline_filter` | string, string, boolean | yes | Same formats as in `region_filters[].read_filters[]`. Note, `genome_filters[]` do not have a `read_threshold` setting as they filter calls rather than molecules. |

## Sensitivity estimation settings

Settings for estimating sensitivity for SBS and indels. Sensitivity is estimated separated for SBS and indels first for mutations (double-strand events) by counting the number of a set of high-confidence germline variants that are detected per the number of opportunities to detect them in the final interrogated bases (while including interrogated bases that were only filtered by germline-related filters), and correcting for zygosity of the germline variant.

After calculating sensitivity for mutations, sensitivity estimates are assigned to SBSindel_call_types `mismatch-ss` and `mismatch-ds` (single-strand events) as the square-root of the respective call type's estimated sensitivity for mutations. MDB sensitivities for all SBSindel_call_types, including `mismatch-os` and `match` calls, are not calculated and instead set manually in the `call_types` section.

Below also are filter settings to select high-confidence germline variants used for sensitivity analysis. These filters are calculated and applied separately for each `germline_vcf_file` and all filters must pass for the variant in all germline VCF files in order for the variant to be used for sensitivity estimation.

Sex chromosomes and mitochondrial chromosome are excluded from sensitivity analysis, due to the complication of pseudo-autosomal regions.

| Key | Type | Required | Description |
| --- | --- | --- | --- |
| `sensitivity_parameters.use_chromgroup` | string | optional | Chromgroup in which sensitivity is estimated, and these estimates are then used for all other chromgroups. Leave blank to skip sensitivity estimation (sets sensitivities for all SBS and indel call types to 1.0). |
| `sensitivity_parameters.sensitivity_vcf` | path | yes | External population reference germline VCF intersected with the analyzed individual's germline variants in order to retain only high-confidence variants for sensitivity analysis. |
| `sensitivity_parameters.genotype` | string (`heterozygous` or `homozygous`) | yes | Determines which germline variants are used for sensitivity estimation. Heterozygous germline variants are more numerous, but sensitivity estimates per site are noisier due to allele-sampling noise, while when using homozygous variants there is more noise due to fewer sites. |
| `sensitivity_parameters.SBS_min_Depth_quantile`, `SBS_min_VAF`, `SBS_max_VAF`, `SBS_min_GQ_quantile`, `SBS_min_QUAL_quantile` | numeric | yes | Thresholds applied to germline SBS variants to include them in sensitivity analysis. Quantile parameters are calculated separately for each `germline_vcf_file`. |
| `sensitivity_parameters.SBS_min_variant_detections` | numeric | yes | Minimum number of times high-confidence germline SBS variants were detected in the HiDEF-seq sample (counting 1 detection per molecule) in order to calculate sensitivity, and otherwise sensitivity is set to the default 1.0 |
| `sensitivity_parameters.indel_min_Depth_quantile`, `indel_min_VAF`, `indel_max_VAF`, `indel_min_GQ_quantile`, `indel_min_QUAL_quantile`, `indel_min_variant_detections` | numeric | yes | Analogous settings for indel sensitivity estimation. |
