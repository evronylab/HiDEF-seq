# HiDEF-seq pipeline outputs

This guide documents every final artifact produced by the `processReads` and `outputResults` workflows described in `main.nf`, including directory layouts, serialized objects, and column-level metadata.

## Outline
- [Directory layout](#directory-layout)
  - [A. Sample root](#a-sample-root)
  - [B. Shared logs](#b-shared-logs)
- [processReads outputs](#processreads-outputs)
  - [A. CCS BAM artifacts](#a-ccs-bam-artifacts)
  - [B. QC contributions to shared logs](#b-qc-contributions-to-shared-logs)
- [outputResults outputs](#outputresults-outputs)
  - [A. Top-level files](#a-top-level-files)
  - [B. Chromgroup directory structure](#b-chromgroup-directory-structure)
  - [C. Filter statistics (filterStats)](#c-filter-statistics-filterstats)
  - [D. Final calls](#d-final-calls)
  - [E. Germline variant calls](#e-germline-variant-calls)
  - [F. Final-call spectra](#f-final-call-spectra)
  - [G. Interrogated-base spectra](#g-interrogated-base-spectra)
  - [H. Genome spectra](#h-genome-spectra)
  - [I. Sensitivity summaries](#i-sensitivity-summaries)
  - [J. Final-call burdens](#j-final-call-burdens)
  - [K. Estimated SBS mutation error probability](#k-estimated-sbs-mutation-error-probability)
  - [L. Serialized object (.qs2)](#l-serialized-object-qs2)
  - [M. Shared logs from outputResults](#m-shared-logs-from-outputresults)

## Directory layout

### A. Sample root
For each sample (`sample_id`) belonging to an individual (`individual_id`), HiDEF-seq writes results beneath
`[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/`.
Workflow-specific subdirectories (`processReads`, `outputResults`, …) mirror the executed stages.

### B. Shared logs
All workflows publish command transcripts and stage-specific QC summaries to
`[analysis_output_dir]/[analysis_id].sharedLogs/`. Expect entries such as:
- `${analysis_id}.${individual_id}.${sample_id}.processReads.command.log` and matching logs for other processes (`splitBAMs`, `prepareFilters`, `extractCalls`, `filterCalls`, `calculateBurdens`, `outputResults`, `removeIntermediateFiles`).
- `*.zmwcount.txt` files capturing molecule counts after major read-processing checkpoints.
- CCS-specific metrics (for example `*.ccs_report.json`, `*.summary.json`), Lima barcode summaries (`*.lima.summary`, `*.lima.counts`), and intermediate QC tables emitted via the numerous `publishDir` statements in `main.nf`.

## processReads outputs

### A. CCS BAM artifacts
Location: `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/processReads/`

| File | Description |
| --- | --- |
| `${analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam` | HiFi CCS reads that passed adapter trimming, demultiplexing, and pbmm2 alignment; sorted by coordinate. |
| `${analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.bai` | BAM index created with `samtools index`. |
| `${analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.pbi` | PacBio BAM index created with `pbindex`. |

### B. QC contributions to shared logs
During `processReads`, the pipeline adds to the shared logs directory:
- Stage-by-stage ZMW counts produced by the `countZMWs` helper.
- `.command.log` transcripts for each Nextflow process touching the sample.
- If CCS generation is enabled, PacBio CCS statistics and Lima barcode summaries copied from the run workspaces.

## outputResults outputs

### A. Top-level files
Location: `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/outputResults/`

| File | Description |
| --- | --- |
| `${analysis_id}.${individual_id}.${sample_id}.yaml_config.tsv` | Snapshot of the YAML parameter file used for the run. |
| `${analysis_id}.${individual_id}.${sample_id}.run_metadata.tsv` | Read-group metadata joined with sample annotations. Columns reflect BAM header fields (`rg_id`, `movie_id`, `SM`, `PM`, `PL`, `DS`, `PU`, etc.) and may vary with instrument metadata. |
| `${analysis_id}.${individual_id}.${sample_id}.qs2` | Serialized object described in [L. Serialized object (.qs2)](#l-serialized-object-qs2). |

### B. Chromgroup directory structure
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

### C. Filter statistics (filterStats)
Each combination of chromgroup and filter group yields three TSV tables named
`[chromgroup]/filterStats/[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].{table}.tsv`.

| Table | Columns |
| --- | --- |
| `molecule_stats.by_run_id.tsv` | `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `run_id`, `stat`, `value`. `stat` encodes sequential filter labels (for example `num_molecules.min_rq_eachstrand.passfilter`). |
| `molecule_stats.by_analysis_id.tsv` | Same as above without `run_id`, aggregating across runs. |
| `region_genome_filter_stats.tsv` | `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `filter`, `binsize`, `threshold`, `padding`, `region_filter_threshold_file`, `num_genomebases_individually_filtered`, `num_genomebases_remaining`. |

### D. Final calls
Files live under `[chromgroup]/finalCalls/` and are named
`[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].[call_class].[call_type].[SBSindel_call_type].*`.
The pipeline produces four artifacts per subtype:

- **`.finalCalls.tsv`** — One row per strand-specific call that passed all filters.
  - Metadata columns: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`, `analysis_chunk`, `run_id`, `zm`.
  - Coordinates and context: `seqnames`, `start_refspace`, `end_refspace`, `ref_plus_strand`, `alt_plus_strand`, `reftnc_plus_strand`, `alttnc_plus_strand`, `reftnc_pyr`, `alttnc_pyr`, `indel_width`, optional MDB-specific scores.
  - Strand-resolved attributes appear twice with suffixes `_refstrand_plus_read` and `_refstrand_minus_read` (for example `start_queryspace`, `end_queryspace`, `qual`, `qual.opposite_strand`, `sa`, `sm`, `sx`, `call_type.opposite_strand`, `deletion.bothstrands.startendmatch`). Multi-value fields are comma-delimited.
- **`.finalCalls_unique.tsv`** — One row per unique mutation (strand columns removed, retains metadata and context fields).
- **`.finalCalls.vcf.bgz`** — Bgzipped VCF containing the same calls as `.finalCalls.tsv`; INFO fields mirror TSV columns except those mapped to CHROM, POS, REF, and ALT. Indexed with `.tbi` courtesy of `writeVcf(index = TRUE)`.
- **`.finalCalls_unique.vcf.bgz`** — Bgzipped, indexed VCF of unique mutations.

### E. Germline variant calls
Stored under `[chromgroup]/germlineVariantCalls/` with filenames ending `.germlineVariantCalls.tsv` or `.germlineVariantCalls.vcf.bgz`.

- Metadata columns match the final-call tables (`analysis_id` through `indel_width`).
- Germline filter flags without strand suffixes: `germline_vcf.passfilter`, any `germline_vcf_indel_region_filter_*`, `max_BAMVariantReads.passfilter`, `max_BAMVAF.passfilter`, and `region_read_filter_*` / `region_genome_filter_*` columns marked `is_germline_filter: true` in the YAML.
- Strand-specific columns reuse the `_refstrand_plus_read` / `_refstrand_minus_read` suffix pattern for query positions, qualities, subread metrics, and match flags.
- VCFs are bgzipped/indexed and share the same INFO mappings as the TSVs.

### F. Final-call spectra
Located under `[chromgroup]/finalCalls.spectra/` with filenames suffixed by spectrum type.
Each TSV contains metadata columns (`analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`) and spectrum-specific fields:

- `.finalCalls.reftnc_pyr.tsv` / `.finalCalls_unique.reftnc_pyr.tsv`: `reftnc_pyr`, `count`, `fraction`.
- `.finalCalls.reftnc_template_strand.tsv`: `reftnc_template_strand`, `count`, `fraction`.
- `.finalCalls.reftnc_pyr_spectrum.tsv` / `.finalCalls_unique.reftnc_pyr_spectrum.tsv`: `channel`, `count`, `fraction`.
- `.finalCalls.reftnc_template_strand_spectrum.tsv`: `channel`, `count`, `fraction`.
- `.finalCalls.refindel_spectrum.sigfit.tsv` / `.finalCalls_unique.refindel_spectrum.sigfit.tsv`: `channel`, `count`, `fraction` for indel contexts.

Each TSV is paired with a `.pdf` spectrum plot generated via `plot_spectrum`.

### G. Interrogated-base spectra
Under `[chromgroup]/interrogatedBases.spectra/`:

- `.bam.gr.filtertrack.reftnc_pyr.tsv`: Metadata plus `reftnc_pyr`, `count`, `fraction` summarising duplex coverage.
- `.bam.gr.filtertrack.reftnc_both_strands.tsv`: Metadata plus `reftnc`, `count`, `fraction` across template orientations.

### H. Genome spectra
`[chromgroup]/genome.spectra/` stores reference-context summaries:

- `genome.reftnc_pyr.tsv`: `analysis_id`, `individual_id`, `sample_id`, `reftnc_pyr`, `count`, `fraction`.
- `genome.reftnc_both_strands.tsv`: Same metadata with `reftnc`, `count`, `fraction`.
- `[chromgroup].genome_chromgroup.reftnc_pyr.tsv`: Chromgroup-restricted `reftnc_pyr`, `count`, `fraction`.
- `[chromgroup].genome_chromgroup.reftnc_both_strands.tsv`: Chromgroup-restricted `reftnc`, `count`, `fraction`.

### I. Sensitivity summaries
`[chromgroup]/sensitivity/[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].sensitivity.tsv`
contains:

- Metadata: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`.
- Estimates: `sensitivity`, `sensitivity_source`, `high_confidence_germline_vcf_variants_sum_zm_detected`, `high_confidence_germline_vcf_variants_sum_duplex_coverage`.
- Nested column `high_confidence_germline_vcf_variants` (list-column describing contributing variants) empty when sensitivity falls back to defaults.

### J. Final-call burdens
`[chromgroup]/finalCalls.burdens/[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].finalCalls.burdens.tsv`
originates from `calculateBurdensChromgroupFiltergroup` and reports burden statistics before/after sensitivity correction:

- Metadata: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`.
- Counts: `num_calls`, `interrogated_bases_or_bp`, `sensitivity`, `sensitivity_source`, `sensitivity_corrected` (boolean).
- Rates and intervals: `burden_calls`, `burden_calls_lci`, `burden_calls_uci`, `num_calls_lci`, `num_calls_uci`.

### K. Estimated SBS mutation error probability
`[chromgroup]/estimatedSBSMutationErrorProbability/[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].estimatedSBSMutationErrorProbability.tsv`
contains:

- Metadata columns mirroring the burden tables.
- Nested tibble `by_channel_pyr` expanded during TSV export with `channel_pyr` and `error_prob` for each trinucleotide channel.
- Scalar column `total` representing the summed SBS mutation error probability across channels.

### L. Serialized object (.qs2)
The `.qs2` file stores a named list of R objects for downstream programmatic consumption. Components and schemas are summarised below.

| Component | Type | Column summary |
| --- | --- | --- |
| `yaml.config` | named list | Parsed representation of the YAML parameters. |
| `run_metadata` | tibble | Same schema as `${analysis_id}.${individual_id}.${sample_id}.run_metadata.tsv`. |
| `individual_id` | character scalar | Individual identifier. |
| `sample_id` | character scalar | Sample identifier. |
| `chromgroup` | character scalar | Chromgroup represented by the current `.qs2` fragment. |
| `filtergroup` | character scalar | Filter group identifier tied to the current fragment. |
| `call_types` | tibble | Per-call-type configuration subset drawn from the YAML (`call_type`, `call_class`, `SBSindel_call_type`, `filtergroup`, `analyzein_chromgroups`, MDB-specific fields). |
| `molecule_stats.by_run_id` | tibble | Columns described in [C. Filter statistics (filterStats)](#c-filter-statistics-filterstats) row 1. |
| `molecule_stats.by_analysis_id` | tibble | Columns described in [C. Filter statistics (filterStats)](#c-filter-statistics-filterstats) row 2. |
| `region_genome_filter_stats` | tibble | Columns described in [C. Filter statistics (filterStats)](#c-filter-statistics-filterstats) row 3. |
| `bam.gr.filtertrack.bytype.coverage_tnc` | tibble | `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`, `reftnc_pyr`, `count`, `fraction`. |
| `finalCalls` | tibble | Strand-resolved calls with the column layout described in [D. Final calls](#d-final-calls). |
| `finalCalls.bytype` | nested tibble | List-columns containing per-subtype tibbles used to render TSV/VCF outputs; each inner tibble shares the schemas outlined for final calls and unique calls. |
| `germlineVariantCalls` | tibble | Strand-resolved germline calls with the schema described in [E. Germline variant calls](#e-germline-variant-calls). |
| `finalCalls.reftnc_spectra` | tibble | Spectrum tables described in [F. Final-call spectra](#f-final-call-spectra). |
| `finalCalls.burdens` | tibble | Burden metrics described in [J. Final-call burdens](#j-final-call-burdens) prior to TSV export. |
| `genome.reftnc` | tibble | Genome-wide spectra described in [H. Genome spectra](#h-genome-spectra). |
| `genome_chromgroup.reftnc` | tibble | Chromgroup-restricted spectra described in [H. Genome spectra](#h-genome-spectra). |
| `sensitivity` | tibble | Sensitivity summaries described in [I. Sensitivity summaries](#i-sensitivity-summaries). |
| `estimatedSBSMutationErrorProbability` | tibble | Error-probability tables described in [K. Estimated SBS mutation error probability](#k-estimated-sbs-mutation-error-probability). |

### M. Shared logs from outputResults
Within the chromgroup directories, `outputResults` records a `.command.log` per invocation (e.g., `${analysis_id}.${individual_id}.${sample_id}.outputResultsSample.command.log`) and copies spectrum-generation logs. These complement the global shared logs and provide traceability for each chromgroup/filtergroup combination processed during result packaging.
