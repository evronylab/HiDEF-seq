# HiDEF-seq pipeline outputs

This guide documents every final file produced by the `processReads` and `outputResults` workflows described in `main.nf`, including directory layouts, serialized objects, and column-level metadata.

## Outline
- [Directory layout](#directory-layout)
- [processReads outputs](#processreads-outputs)
- [outputResults outputs](#outputresults-outputs)

## Directory layout

### Sample root
For each sample (`sample_id`) belonging to an individual (`individual_id`), HiDEF-seq writes results beneath
`[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/`.
Workflow-specific subdirectories (`processReads`, `outputResults`, …) mirror the executed stages.

### Shared logs
Global logs and run metadata accumulate under `[analysis_output_dir]/[analysis_id].sharedLogs/`. Key contents include:
- `runParams.*.yaml` and `runInfo.*.txt` snapshots written at pipeline launch.
- `.command.log` files for workflow-wide tasks such as `makeBarcodesFasta`, `countZMWs`, `combineCCSChunks`, `prepareFilters`, `calculateBurdens`, and `outputResults` (see `publishDir "${sharedLogsDir}"` statements in `main.nf`).
- Molecule-count summaries from `countZMWs` (`*.zmwcount.txt`).
- CCS QC files copied from PacBio outputs (for example `statistics/*.ccs_report.json`, `statistics/*.summary.json`) and Lima demultiplexing summaries (`*.lima.summary`, `*.lima.counts`).
- Any additional helper tables emitted globally, including cached configuration manifests and chunk-level tracking TSVs.

### Per-sample logs
Each sample directory contains a `logs/` subfolder: `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/logs/`. This folder stores `.command.log` transcripts for sample-scoped processes (`pbmm2Align`, `mergeAlignedSampleBAMs`, `splitBAMs`, `filterCallsChunk`, `calculateBurdensChromgroupFiltergroup`, `outputResultsSample`, etc.), renamed with the full `analysis_id.individual_id.sample_id.processName.command.log` pattern for clarity.

## processReads outputs

### CCS BAM files
Location: `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/processReads/`

| File | Description |
| --- | --- |
| `${analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam` | HiFi CCS reads that passed adapter trimming, demultiplexing, and pbmm2 alignment; sorted by coordinate. |
| `${analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.bai` | BAM index created with `samtools index`. |
| `${analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.pbi` | PacBio BAM index created with `pbindex`. |

### BAM tag reference
The aligned CCS BAM retains PacBio-rich annotations used by downstream filters. The table below summarises notable tags and columns:

| Tag | Category | Description |
| --- | --- | --- |
| `ec` | General | Effective coverage: average subread coverage across windows, counting all polished subreads (full- and partial-length) that pass filters. Typically approximates `np + 1`. |
| `np` | General | Number of full-length subreads used during CCS polishing. |
| `rq` | General | Mean per-base quality across the CCS read. |
| `sn` | General | Signal-to-noise ratios for nucleotides A, C, G, T across the HQRegion. |
| `we` | General | Start of the last base (`qe - 1`) in approximate raw frame counts from movie start. |
| `ws` | General | Start of the first base in approximate raw frame counts from movie start. |
| `zm` | General | ZMW hole number. |
| `RG` | General | Read-group identifier propagated from alignment headers. |
| `ac` | Barcode/adapter | Array `[left_detected, left_missing, right_detected, right_missing]` counting adapters. |
| `bx` | Barcode/adapter | Pair of clipped barcode sequence lengths. |
| `ff` | Barcode/adapter | Failure flag (`4` indicates non-failed single-stranded CCS reads). |
| `ma` | Barcode/adapter | Adapter detection status (`0` both sides detected, `1` missing left, `2` missing right; adapter called when >25% of subreads detect it). |
| `qs` | Barcode/adapter | Query start position after barcode trimming. |
| `qe` | Barcode/adapter | Query end position before barcode trimming. |
| `bc` | Barcode/adapter | Integer barcode pair indices (0-based offsets within barcode FASTA). |
| `bq` | Barcode/adapter | Quality of the barcode call. |
| `cx` | Barcode/adapter | Sum of subread local context flags (adapter/barcode presence, pass orientation, etc.); CCS BAMs downstream of pbmm2 set `cx=12`. |
| `bl` | Barcode/adapter | Barcode sequence clipped from the leading end. |
| `bt` | Barcode/adapter | Barcode sequence clipped from the trailing end. |
| `ql` | Barcode/adapter | Qualities of barcode bases clipped from the leading end. |
| `qt` | Barcode/adapter | Qualities of barcode bases clipped from the trailing end. |
| `ls` | Barcode/adapter | Binary blob storing clipped adapter/barcode data. |
| `sa` | Rich HiFi | Run-length encoded counts of subread alignments spanning each consensus position. |
| `sm` | Rich HiFi | Count of aligned subread bases matching the consensus per position (values >255 capped). |
| `sx` | Rich HiFi | Count of aligned subread bases mismatching the consensus per position (values >255 capped). |
| `MM` | Methylation | Base modification calls per BAM specification. |
| `ML` | Methylation | Modification probabilities corresponding to `MM`. |
| `ip` | Kinetics | Inter-pulse width measurements. |
| `pw` | Kinetics | Pulse width measurements. |
| `mg` | Alignment | Gap-compressed sequence identity percentage (`100 * matches / (matches + mismatches + gaps)`). |
| `NM` | Alignment | Total mismatches plus indel events for the alignment. |
| `pos` | Alignment summary | 1-based genomic start position excluding clipping. |
| `qwidth` | Alignment summary | Length of the query sequence including soft clipping. |
| `isize` | Alignment summary | Insert size / template length (`TLEN`), sum of `M/D/=/X` CIGAR operations. |

Downstream filters often rely on `sa` counts and the consensus agreement ratio `sm / (sm + sx)`; the latter can slightly underestimate agreement when subreads per strand exceed 255 because `sm` is capped.
Consult the <a href="https://pacbiofileformats.readthedocs.io/en/13.1/BAM.html" target="_blank" rel="noopener noreferrer">PacBio BAM format reference</a> for exhaustive tag definitions.

### QC contributions to shared logs
During `processReads`, the pipeline adds to the shared logs directory:
- Stage-by-stage ZMW counts produced by `countZMWs`.
- Per-chunk CCS diagnostics (`statistics/*.ccs_report.*`, `statistics/*.summary.json`) and their corresponding `.command.log` files from `ccsChunk`, `mergeCCSchunks`, and `filterAdapter`.
- Lima demultiplexing summaries (`*.lima.summary`, `*.lima.counts`) plus associated `.command.log` transcripts.
- Any additional `.command.log` outputs from global helper steps executed within `processReads`.

## outputResults outputs

### Directory structure
`outputResults` organises every downstream result by chromgroup and filter group. The script defines helper suffixes such as
`filterStats_dir <- "/filterStats/"` and `finalCalls_dir <- "/finalCalls/"`; when concatenated with a chromgroup name the leading slash becomes a path separator. Consequently, each chromgroup has its own directory under `outputResults`:

```
outputResults/
  └─ [chromgroup]/
       ├─ filterStats/
       ├─ finalCalls/
       ├─ germlineVariantCalls/
       ├─ finalCalls.spectra/
       ├─ interrogatedBases.spectra/
       ├─ genome.spectra/
       ├─ sensitivity/
       ├─ finalCalls.burdens/
       └─ estimatedSBSMutationErrorProbability/
```

Files inside these folders are further keyed by `filtergroup`, `call_class`, `call_type`, and `SBSindel_call_type`.
All burden-related tables originate from the `calculateBurdensChromgroupFiltergroup` tasks and are published here by `outputResults`.

### Top-level files
Location: `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/outputResults/`

| File | Description |
| --- | --- |
| `${analysis_id}.${individual_id}.${sample_id}.yaml_config.tsv` | Snapshot of the YAML parameter file used for the run. |
| `${analysis_id}.${individual_id}.${sample_id}.run_metadata.tsv` | Read-group metadata joined with sample annotations. Columns reflect BAM header fields (`rg_id`, `movie_id`, `SM`, `PM`, `PL`, `DS`, `PU`, etc.) and may vary with instrument metadata. |
| `${analysis_id}.${individual_id}.${sample_id}.qs2` | Serialized object described in [Serialized object (.qs2)](#serialized-object-qs2). |

### Filter statistics (filterStats)
Each combination of chromgroup and filter group yields three TSV tables named
`[chromgroup]/filterStats/[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].{table}.tsv`.

| Table | Columns |
| --- | --- |
| `molecule_stats.by_run_id.tsv` | `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `run_id`, `stat`, `value`. `stat` encodes sequential filter labels (for example `num_molecules.min_rq_eachstrand.passfilter`). |
| `molecule_stats.by_analysis_id.tsv` | Same as above without `run_id`, aggregating across runs. |
| `region_genome_filter_stats.tsv` | `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `filter`, `binsize`, `threshold`, `padding`, `region_filter_threshold_file`, `num_genomebases_individually_filtered`, `num_genomebases_remaining`. |

### Final calls
Files live under `[chromgroup]/finalCalls/` and are named
`[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].[call_class].[call_type].[SBSindel_call_type].*`.
The pipeline produces four files per subtype:

- **`.finalCalls.tsv`** — One row per strand-specific call that passed all filters.
  - Metadata columns: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`, `analysis_chunk`, `run_id`, `zm`.
  - Coordinates and context: `seqnames`, `start_refspace`, `end_refspace`, `ref_plus_strand`, `alt_plus_strand`, `reftnc_plus_strand`, `alttnc_plus_strand`, `reftnc_pyr`, `alttnc_pyr`, `indel_width`, optional MDB-specific scores.
  - Strand-resolved attributes appear twice with suffixes `_refstrand_plus_read` and `_refstrand_minus_read` (for example `start_queryspace`, `end_queryspace`, `qual`, `qual.opposite_strand`, `sa`, `sm`, `sx`, `call_type.opposite_strand`, `deletion.bothstrands.startendmatch`). Multi-value fields are comma-delimited.
- **`.finalCalls_unique.tsv`** — One row per unique mutation (strand columns removed, retains metadata and context fields).
- **`.finalCalls.vcf.bgz`** — Bgzipped VCF containing the same calls as `.finalCalls.tsv`; INFO fields mirror TSV columns except those mapped to CHROM, POS, REF, and ALT. Indexed with `.tbi` courtesy of `writeVcf(index = TRUE)`.
- **`.finalCalls_unique.vcf.bgz`** — Bgzipped, indexed VCF of unique mutations.

### Germline variant calls
Stored under `[chromgroup]/germlineVariantCalls/` with filenames ending `.germlineVariantCalls.tsv` or `.germlineVariantCalls.vcf.bgz`.

- Metadata columns match the final-call tables (`analysis_id` through `indel_width`).
- Germline filter flags without strand suffixes: `germline_vcf.passfilter`, any `germline_vcf_indel_region_filter_*`, `max_BAMVariantReads.passfilter`, `max_BAMVAF.passfilter`, and `region_read_filter_*` / `region_genome_filter_*` columns marked `is_germline_filter: true` in the YAML.
- Strand-specific columns reuse the `_refstrand_plus_read` / `_refstrand_minus_read` suffix pattern for query positions, qualities, subread metrics, and match flags.
- VCFs are bgzipped/indexed and share the same INFO mappings as the TSVs.

### Final-call spectra
Located under `[chromgroup]/finalCalls.spectra/` with filenames suffixed by spectrum type.
Each TSV contains metadata columns (`analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`) and spectrum-specific fields:

- `.finalCalls.reftnc_pyr.tsv` / `.finalCalls_unique.reftnc_pyr.tsv`: `reftnc_pyr`, `count`, `fraction`.
- `.finalCalls.reftnc_template_strand.tsv`: `reftnc_template_strand`, `count`, `fraction`.
- `.finalCalls.reftnc_pyr_spectrum.tsv` / `.finalCalls_unique.reftnc_pyr_spectrum.tsv`: `channel`, `count`, `fraction`.
- `.finalCalls.reftnc_template_strand_spectrum.tsv`: `channel`, `count`, `fraction`.
- `.finalCalls.refindel_spectrum.sigfit.tsv` / `.finalCalls_unique.refindel_spectrum.sigfit.tsv`: `channel`, `count`, `fraction` for indel contexts.

Each TSV is paired with a `.pdf` spectrum plot generated via `plot_spectrum`.

### Interrogated-base spectra
Under `[chromgroup]/interrogatedBases.spectra/`:

- `.bam.gr.filtertrack.reftnc_pyr.tsv`: Metadata plus `reftnc_pyr`, `count`, `fraction` summarising duplex coverage.
- `.bam.gr.filtertrack.reftnc_both_strands.tsv`: Metadata plus `reftnc`, `count`, `fraction` across template orientations.

### Genome spectra
`[chromgroup]/genome.spectra/` stores reference-context summaries:

- `genome.reftnc_pyr.tsv`: `analysis_id`, `individual_id`, `sample_id`, `reftnc_pyr`, `count`, `fraction`.
- `genome.reftnc_both_strands.tsv`: Same metadata with `reftnc`, `count`, `fraction`.
- `[chromgroup].genome_chromgroup.reftnc_pyr.tsv`: Chromgroup-restricted `reftnc_pyr`, `count`, `fraction`.
- `[chromgroup].genome_chromgroup.reftnc_both_strands.tsv`: Chromgroup-restricted `reftnc`, `count`, `fraction`.

### Sensitivity summaries
`[chromgroup]/sensitivity/[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].sensitivity.tsv`
contains:

- Metadata: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`.
- Estimates: `sensitivity`, `sensitivity_source`, `high_confidence_germline_vcf_variants_sum_zm_detected`, `high_confidence_germline_vcf_variants_sum_duplex_coverage`.
- Nested column `high_confidence_germline_vcf_variants` (list-column describing contributing variants) empty when sensitivity falls back to defaults.

### Final-call burdens
`[chromgroup]/finalCalls.burdens/[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].finalCalls.burdens.tsv`
originates from `calculateBurdensChromgroupFiltergroup` and reports burden statistics before/after sensitivity correction:

- Metadata: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`.
- Counts: `num_calls`, `interrogated_bases_or_bp`, `sensitivity`, `sensitivity_source`, `sensitivity_corrected` (boolean).
- Rates and intervals: `burden_calls`, `burden_calls_lci`, `burden_calls_uci`, `num_calls_lci`, `num_calls_uci`.

### Estimated SBS mutation error probability
`[chromgroup]/estimatedSBSMutationErrorProbability/[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].estimatedSBSMutationErrorProbability.tsv`
contains:

- Metadata columns mirroring the burden tables.
- Nested tibble `by_channel_pyr` expanded during TSV export with `channel_pyr` and `error_prob` for each trinucleotide channel.
- Scalar column `total` representing the summed SBS mutation error probability across channels.

### Serialized object (.qs2)
The `.qs2` file captures the final data structures required to regenerate tables and plots without re-running the Nextflow processes. It is saved with the <a href="https://github.com/qsbase/qs2" target="_blank" rel="noopener noreferrer">qs2 R package</a> format and can be loaded directly with that library. Each component is listed below with its schema.

#### yaml.config
Top-level configuration loaded from the analysis YAML. Keys mirror the templates in [`config_templates/`](../config_templates).

| Key | Description |
| --- | --- |
| `analysis_id`, `analysis_output_dir`, `cache_dir` | Identifiers and directories governing output locations and shared caches. |
| `runs` | List of run entries (`run_id`, `reads_file`, barcode definitions). |
| `samples` | Sample metadata linking `sample_id` to `individual_id` and annotations. |
| `individuals` | Germline resources (sex, BAM/VCF paths, technology tags). |
| `chromgroups`, `call_types`, `filtergroups` | Workflow configuration describing chromosome groupings, call classes, and molecule-level filters. |
| `region_read_filters_config`, `region_genome_filters_config` | Region filter definitions expanded from the YAML `region_filters` blocks. |
| `germline_vcf_types_config` | Germline filter thresholds, including padding strings (`m<multiplier>b<offset>`). |
| `BSgenome` | Reference genome package name and optional tarball. |
| `sensitivity_parameters` | Default chromgroup and thresholds used to measure detection sensitivity. |

#### run_metadata
`run_metadata` matches the TSV emitted by `outputResults`. Columns originate from BAM read-group (`@RG`) headers and the YAML configuration.

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id` | Identifiers propagated into every output. |
| `run_id` | Run identifier resolved from the YAML configuration. |
| `rg_id` | Read-group identifier (`@RG:ID`). |
| `movie_id` | PacBio movie identifier (`@RG:PU`). |
| `SM`, `PL`, `PM` | Sample name, platform, and platform model tags copied from the BAM header. |
| `DS.*` | One column per key extracted from the `@RG:DS` field (for example `DS_SEQUENCING_KIT`, `DS_BINDING_KIT`, `DS_CCS_VERSION`); contents depend on instrument metadata. |
| Additional `@RG` tags | Any other read-group tags (e.g. `LB`, `CN`, `DT`) are surfaced verbatim when present. |

#### individual and sample identifiers

| Field | Description |
| --- | --- |
| `individual_id` | Scalar string naming the individual analysed in the fragment. |
| `sample_id` | Scalar string naming the sample. |

#### chromgroup and filtergroup

| Field | Description |
| --- | --- |
| `chromgroup` | Scalar string naming the chromgroup processed when `calculateBurdens` produced the source QS file. |
| `filtergroup` | Scalar string naming the molecule-level filter group associated with the chromgroup slice. |

#### call_types

| Column | Description |
| --- | --- |
| `call_type`, `call_class`, `SBSindel_call_type` | Call taxonomy as defined in the YAML configuration. |
| `filtergroup` | Molecule-level filter group applied to the subtype. |
| `analyzein_chromgroups` | Chromgroup scope (`all` or comma-separated subset). |
| MDB-specific fields | Optional columns (`MDB_bamtag`, `MDB_min_score`, `MDB_sensitivity`, `MDB_base_opposite_strand`) populated for MDB call types. |

#### molecule_stats.by_run_id

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup` | Identifiers for the molecule statistics slice. |
| `run_id` | Run contributing to the statistics row. |
| `stat` | Filter stage label (e.g. `num_molecules.passalignmentfilter`). |
| `value` | Integer count recorded at the corresponding filter stage. |

#### molecule_stats.by_analysis_id

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup` | Identifiers. |
| `stat` | Filter stage label aggregated across runs. |
| `value` | Count aggregated across runs. |

#### region_genome_filter_stats

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup` | Identifiers. |
| `filter` | Filter label derived from the threshold file basename. |
| `binsize`, `threshold`, `padding` | Region filter parameters (threshold strings use `lt`, `lte`, `gt`, `gte`; padding is measured in bases). |
| `region_filter_threshold_file` | Source bigWig used to derive the mask. |
| `num_genomebases_individually_filtered` | Bases removed by the filter prior to merging. |
| `num_genomebases_remaining` | Bases retained after filtering. |

#### finalCalls

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type` | Call metadata. |
| `analysis_chunk`, `run_id`, `zm` | Chunk identifier, sequencing run, and ZMW hole number. |
| `seqnames`, `start_refspace`, `end_refspace` | Reference contig and 1-based coordinates. |
| `ref_plus_strand`, `alt_plus_strand` | Forward-strand reference and alternate alleles. |
| `reftnc_plus_strand`, `alttnc_plus_strand`, `reftnc_pyr`, `alttnc_pyr` | Context nucleotides in forward and pyrimidine orientation. |
| `indel_width` | Event width in reference bases (0 for SBS). |
| `refstrand` | Indicates which strand supplied the call (`refstrand_plus_read`/`refstrand_minus_read`). |
| `start_queryspace`, `end_queryspace` | Comma-delimited query-space coordinates for the reporting strand. |
| `qual`, `qual.opposite_strand` | Base qualities on the reporting and opposite strand (comma-delimited Phred values). |
| `sa`, `sa.opposite_strand` | Run-length encoded subread coverage counts per position. |
| `sm`, `sm.opposite_strand`; `sx`, `sx.opposite_strand` | Run-length encoded match and mismatch counts per position. |
| `call_class.opposite_strand`, `call_type.opposite_strand`, `alt_plus_strand.opposite_strand` | Allele annotations observed on the opposite strand. |
| `deletion.bothstrands.startendmatch` | Flag indicating whether deletion endpoints agree across strands. |

#### finalCalls_for_tsv

| Column | Description |
| --- | --- |
| Metadata columns | Same metadata and context fields listed for `finalCalls`. |
| `*_refstrand_plus_read`, `*_refstrand_minus_read` | For each strand-resolved metric in `finalCalls` (`start_queryspace`, `end_queryspace`, `qual`, `sa`, `sm`, `sx`, `call_type`, `alt_plus_strand`, etc.), the data are pivoted into explicit columns suffixed with `_refstrand_plus_read` and `_refstrand_minus_read`. Values remain comma-delimited where appropriate. |

#### finalCalls_unique_for_tsv

| Column | Description |
| --- | --- |
| Metadata columns | Identical to `finalCalls_for_tsv`. |
| Strand-specific fields | Removed to ensure one row per unique mutation; coordinates and context fields remain. |

#### finalCalls.bytype

| Column | Description |
| --- | --- |
| Metadata columns | Same identifiers and context fields as `finalCalls`. |
| `finalCalls_for_tsv`, `finalCalls_unique_for_tsv` | List-columns containing the tables described above prior to unnesting. |
| `finalCalls_for_vcf`, `finalCalls_unique_for_vcf` | List-columns of VCF-ready tables (after left-normalisation of indels). |

#### germlineVariantCalls

| Column | Description |
| --- | --- |
| Metadata columns | Same as `finalCalls_for_tsv`. |
| Germline flags | `germline_vcf.passfilter`, `germline_vcf_indel_region_filter_*`, `max_BAMVariantReads.passfilter`, `max_BAMVAF.passfilter`, and any `region_*` filters marked as germline in the YAML. |
| Strand-specific fields | `_refstrand_plus_read` and `_refstrand_minus_read` columns mirroring the structure described for `finalCalls_for_tsv`. |

#### finalCalls.reftnc_spectra

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type` | Identifiers. |
| `finalCalls.reftnc_pyr`, `finalCalls_unique.reftnc_pyr` | List-columns of tibbles (`reftnc_pyr`, `count`, `fraction`). |
| `finalCalls.reftnc_template_strand` | List-column with template-strand contexts (`reftnc_template_strand`, `count`, `fraction`). |
| `finalCalls.reftnc_pyr_spectrum`, `finalCalls_unique.reftnc_pyr_spectrum` | List-columns containing spectrum tables (`channel`, `count`, `fraction`). |
| `finalCalls.reftnc_template_strand_spectrum` | List-column for strand-specific spectra (`channel`, `count`, `fraction`). |
| `finalCalls.refindel_spectrum`, `finalCalls_unique.refindel_spectrum` | List-columns of indel spectra matrices. |
| `finalCalls.*.sigfit` columns | Sigfit-formatted matrices generated from the corresponding spectra. |

#### bam.gr.filtertrack.bytype.coverage_tnc

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type` | Identifiers for the duplex coverage slice. |
| `bam.gr.filtertrack.coverage` | List-column of GRanges coverage objects for interrogated bases. |
| `bam.gr.filtertrack.reftnc_pyr`, `bam.gr.filtertrack.reftnc_both_strands` | List-columns of trinucleotide counts and fractions (with genome and chromgroup ratios attached). |

#### genome.reftnc and genome_chromgroup.reftnc

| Column | Description |
| --- | --- |
| `reftnc_pyr`/`reftnc` | Reference trinucleotide label. |
| `count`, `fraction` | Genome-wide or chromgroup-restricted counts and fractions. |

#### sensitivity

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type` | Identifiers. |
| `sensitivity`, `sensitivity_source` | Estimated detection rate and provenance (`use_chromgroup`, default fallback, etc.). |
| `high_confidence_germline_vcf_variants_sum_zm_detected`, `high_confidence_germline_vcf_variants_sum_duplex_coverage` | Counts of germline variants supporting the estimate. |
| `high_confidence_germline_vcf_variants` | List-column summarising the contributing variants (empty when defaults are used). |

#### finalCalls.burdens

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type` | Identifiers. |
| `num_calls`, `interrogated_bases_or_bp` | Raw counts of calls and interrogated sequence. |
| `sensitivity`, `sensitivity_source`, `sensitivity_corrected` | Sensitivity estimate applied to the burden. |
| `burden_calls`, `burden_calls_lci`, `burden_calls_uci`, `num_calls_lci`, `num_calls_uci` | Rates with confidence intervals. |
| `unique_calls`, `reftnc_corrected`, `reftnc_corrected_chromgroup` | Flags indicating whether corrections were applied. |

#### estimatedSBSMutationErrorProbability

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type` | Identifiers. |
| `by_channel_pyr` | List-column with `channel_pyr` and `error_prob` per trinucleotide. |
| `total` | Aggregate SBS mutation error probability. |
