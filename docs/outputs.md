# HiDEF-seq pipeline outputs

This guide documents every final file produced by the `processReads` and `outputResults` workflows described in `main.nf`, including directory layouts, serialized objects, and column-level metadata.

## Outline
- [Main analysis directory layout](#main-analysis-directory-layout)
  - [Sample root directory](#sample-root-directory)
  - [Shared logs](#shared-logs)
  - [Per-sample logs](#per-sample-logs)
- [processReads outputs](#processreads-outputs)
  - [CCS BAM files](#ccs-bam-files)
  - [BAM tag reference](#bam-tag-reference)
  - [QC contributions to shared logs](#qc-contributions-to-shared-logs)
- [outputResults outputs](#outputresults-outputs)
  - [Directory structure](#directory-structure)
  - [Top-level files](#top-level-files)
  - [Coverage and reference trinucleotides](#coverage-and-reference-trinucleotides)
  - [Filter statistics](#filter-statistics)
  - [Final calls](#final-calls)
  - [Germline variant calls](#germline-variant-calls)
  - [Final-call spectra](#final-call-spectra)
  - [Interrogated-base spectra](#interrogated-base-spectra)
  - [Genome spectra](#genome-spectra)
  - [Sensitivity summaries](#sensitivity-summaries)
  - [Final-call burdens](#final-call-burdens)
  - [Estimated SBS mutation error probability](#estimated-sbs-mutation-error-probability)
  - [Serialized object (.qs2)](#serialized-object-qs2)

## Main analysis directory layout

### Sample root directory
For each sample (`sample_id`) belonging to an individual (`individual_id`), HiDEF-seq writes results to the folder `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/`. Within each of these folders are folders for each of the Nextflow pipeline's workflow outputs (`processReads`, `outputResults`, …).

### Shared logs
Global logs and run metadata are saved in `[analysis_output_dir]/[analysis_id].sharedLogs/`. Key contents include:
- `runParams.*.yaml` and `runInfo.*.txt` snapshots written at pipeline launch.
- `.command.log` files for pipeline-wide tasks such as `makeBarcodesFasta`, `countZMWs`, `combineCCSChunks`, `prepareFilters`, `calculateBurdens`, and `outputResults` (see `publishDir "${sharedLogsDir}"` statements in `main.nf`).
- Molecule-count summaries from `countZMWs` (`*.zmwcount.txt`).
- CCS QC files copied from PacBio outputs (for example `statistics/*.ccs_report.json`, `statistics/*.summary.json`) and Lima demultiplexing summaries (`*.lima.summary`, `*.lima.counts`).
- Any additional helper tables emitted globally, including cached configuration manifests and chunk-level tracking TSVs.

### Per-sample logs
Each sample root directory contains a `logs/` subfolder: `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/logs/`. This folder stores `.command.log` transcripts for sample-scoped processes (`pbmm2Align`, `mergeAlignedSampleBAMs`, `splitBAMs`, `filterCallsChunk`, `calculateBurdensChromgroupFiltergroup`, `outputResultsSample`, etc.), renamed with the full `analysis_id.individual_id.sample_id.processName.command.log` pattern for clarity.

## processReads outputs

### CCS BAM files
Location: `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/processReads/`

| File | Description |
| --- | --- |
| `${analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam` | HiFi CCS reads that passed adapter trimming, demultiplexing, and pbmm2 alignment; sorted by coordinate. |
| `${analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.bai` | BAM index created with `samtools index`. |
| `${analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.pbi` | PacBio BAM index created with `pbindex`. |

### BAM tag reference
The table below summarises tags present in the aligned CCS BAM. See references for [PacBio BAM format](https://pacbiofileformats.readthedocs.io/en/13.1/BAM.html), [lima (demultiplexing)](https://lima.how/), [CCS](https://ccs.how/), and [pbmm2](https://github.com/PacificBiosciences/pbmm2) for more information.

| Tag | Category | Description |
| --- | --- | --- |
| `ec` | General | Effective coverage: average subread coverage of the consensus sequence across all [windows](https://ccs.how/how-does-ccs-work.html#4-windowing). |
| `np` | General | Number of full-length subreads used during CCS polishing. |
| `rq` | General | Mean per-base quality across the CCS read. |
| `sn` | General | Signal-to-noise ratios for nucleotides A, C, G, T across the HQRegion. |
| `we` | General | Start of the last base (`qe - 1`) in approximate raw frame counts from movie start. |
| `ws` | General | Start of the first base in approximate raw frame counts from movie start. |
| `zm` | General | ZMW hole number. |
| `RG` | General | Read-group identifier. |
| `ac` | Barcode/adapter | Array of adapter detection: `[left_detected, left_missing, right_detected, right_missing]`. |
| `bx` | Barcode/adapter | Pair of clipped barcode sequence lengths. |
| `ff` | Barcode/adapter | Failure flag. |
| `ma` | Barcode/adapter | Adapter detection status (`0` both sides detected, `1` missing left, `2` missing right). |
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
| `sa` | Rich HiFi | Run-length encoded counts of subread alignments spanning each consensus position. Ignores local details of alignments (i.e., indels) and only considers the start and end of each subread alignment. |
| `sm` | Rich HiFi | Count of aligned subread bases matching the consensus sequence per position (values >255 capped). |
| `sx` | Rich HiFi | Count of aligned subread bases mismatching the consensus sequence per position (values >255 capped). |
| `MM` | Methylation | Base modification calls per BAM specification. |
| `ML` | Methylation | Modification probabilities corresponding to `MM`. |
| `ip` | Kinetics | Inter-pulse width measurements. |
| `pw` | Kinetics | Pulse width measurements. |
| `mg` | Alignment | Gap-compressed sequence identity percentage (`100 * matches / (matches + mismatches + gaps)`). |
| `NM` | Alignment | Total number of mismatches and gaps in the alignment. |

### QC contributions to shared logs
During `processReads`, the pipeline adds to the sharedLogs directory:
- Stage-by-stage ZMW counts produced by `countZMWs`.
- Per-chunk CCS diagnostics (`statistics/*.ccs_report.*`, `statistics/*.summary.json`) and their corresponding `.command.log` files from `ccsChunk`, `mergeCCSchunks`, and `filterAdapter`.
- Lima demultiplexing summaries (`*.lima.summary`, `*.lima.counts`) plus associated `.command.log` transcripts.
- Any additional `.command.log` outputs from global helper steps executed within `processReads`.

## outputResults outputs

### Directory structure
The outputs of `outputResults` are saved in `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/outputResults`, with separate sub-folders for the outputs of each chromgroup as follows:

```
outputResults/
  └─ [chromgroup]/
       ├─ coverage_reftnc/
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

Files inside these folders are further keyed by subsets of `analysis_id`, `individual_id`, `sample_id`, `filtergroup`, `call_class`, `call_type`, and `SBSindel_call_type` fields relevant to the file.

### Top-level files
Location: `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/outputResults/`

| File | Description |
| --- | --- |
| `${analysis_id}.${individual_id}.${sample_id}.yaml_config.tsv` | Snapshot of the YAML parameter file used for the run. |
| `${analysis_id}.${individual_id}.${sample_id}.run_metadata.tsv` | Read-group metadata joined with sample annotations. Columns reflect BAM header fields (`rg_id`, `movie_id`, `SM`, `PM`, `PL`, `DS`, `PU`, etc.). |
| `${analysis_id}.${individual_id}.${sample_id}.qs2` | Serialized object enabling analysis of all outputs in R. Contents described in [Serialized object (.qs2)](#serialized-object-qs2). |

### Coverage and reference trinucleotides

Location: `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/outputResults/[chromgroup]/coverage_reftnc/`

For every combination of `call_class`, `call_type`, and `SBSindel_call_type`, the `calculateBurdensChromgroupFiltergroup` process that is run for each `sample_id`, `chromgroup`, and `filtergroup` writes a bgzipped BED file named:
 `[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].[call_class].[call_type].[SBSindel_call_type].coverage_reftnc.bed.gz`, with a companion Tabix index (`.tbi`) that contains genome-wide per-base HiDEF-seq duplex coverage (final interrogated bases) and reference trincucleotide sequences for all non-zero coverage positions:

| Column               | Description                                                  |
| -------------------- | ------------------------------------------------------------ |
| `seqnames`           | Reference sequence name.                                     |
| `start`              | Zero-based inclusive start position.                         |
| `end`                | Zero-based exclusive end position.                           |
| `duplex_coverage`    | Number of final interrogated duplex base pairs at the position (i.e., coverage by each duplex molecule is counted as 1). |
| `reftnc_plus_strand` | Reference trinucleotide context of the plus strand at the position. |

These files originate from `bin/calculateBurdens.R` and are moved into the `outputResults` folder hierarchy by the `calculateBurdensChromgroupFiltergroup` process so they accompany other per-chromgroup summaries.

### Filter statistics

Each combination of chromgroup and filter group yields three TSV tables named:
`[chromgroup]/filterStats/[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].{table}.tsv`.

| Table | Columns |
| --- | --- |
| `molecule_stats.by_run_id.tsv` | Number of molecules and number of sequenced bases (in reference space) remaining after each filter, for each run_id. Columns: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `run_id`, `stat`, `value`. `stat` contains the filter label and type of statistic, and `value` contains the value of the statistic. |
| `molecule_stats.by_analysis_id.tsv` | Same as above without aggregating across `run_id`. |
| `region_genome_filter_stats.tsv` | Number of genome bases remaining after each genome region filter. Columns: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `filter`, `binsize`, `threshold`, `padding`, `region_filter_threshold_file`, `num_genomebases_individually_filtered` (number of genome bases removed when applying only this filter), `num_genomebases_remaining` (number of genome bases remaining after applying this and all prior filters). |

### Final calls
Final calls files are output to `[chromgroup]/finalCalls/` and are named
`[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].[call_class].[call_type].[SBSindel_call_type].*`.
The pipeline produces four files per subtype:

- **`.finalCalls.tsv`** — Tables of the final calls that passed all filters.
  - Shared identifiers: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`, `analysis_chunk`, `run_id`, `zm`.
  - Coordinates in genome reference space and alleles: `seqnames`, `start_refspace`, `end_refspace`, `ref_plus_strand`, `alt_plus_strand`.
  - Strand-resolved coordinates in query space (i.e., read coordinates) and other measurements: `start_queryspace`, `end_queryspace`, `qual`, `qual.opposite_strand`, `sa`, `sm`, `sx`). Measurements for multi-base calls, such as `qual` are stored as comma-delimited values for every base. SBS and indel mutation calls are pivoted to one row per event with strand-specific columns suffixed `_refstrand_plus_read` and `_refstrand_minus_read`. Non-mutation call types remain one row per strand with a `refstrand` column annotating to which reference strand the call's read aligned.
  - Columns specific to each `call_class`:
    * **SBS** and **MDB** entries include:
      * `reftnc_plus_strand` and `alttnc_plus_strand` — reference and call trinucleotide sequence at the call position on the reference genome plus strand.
      * `reftnc_pyr` and `alttnc_pyr`— reference and call trinucleotide sequences collapsed to central pyrimidine trinucleotide sequences.
      * `reftnc_synthesized_strand` and `alttnc_synthesized_strand` — reference and call trinucleotide sequence at the call position on the strand synthesized by the sequencer polymerase.
      * `reftnc_template_strand` and `alttnc_template_strand` — reference and call trinucleotide sequence at the call position on the strand replicated by the sequencer polymerase (i.e., template strand).
    * **Indel** entries include: `indel_width`.
    * **MDB** entries include `MDB_score`.
- **`.finalCalls_unique.tsv`** — Tables of with only one row per unique SBS or indel mutation, i.e., retaining only one row for calls detected in > 1 molecule. Not created for single-strand call types.
  - Shared identifiers: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`.
  - Coordinates in genome reference space and alleles: `seqnames`, `start_refspace`, `end_refspace`, `ref_plus_strand`, `alt_plus_strand`.
  - **SBS** entries include: `reftnc_plus_strand`, `alttnc_plus_strand`, `reftnc_pyr` and `alttnc_pyr`.
  - **Indel** entries include: `indel_width`.

- **`.finalCalls.vcf.bgz`** — Bgzipped VCF containing the same calls as `.finalCalls.tsv`; INFO fields mirror TSV columns except those mapped to CHROM, POS, REF, and ALT. Indexed with `.tbi`.
- **`.finalCalls_unique.vcf.bgz`** — Bgzipped, indexed VCF containing the same calls as `.finalCalls_unique.tsv`. INFO fields mirror TSV columns except those mapped to CHROM, POS, REF, and ALT. Indexed with `.tbi`.

### Germline variant calls
Stored under `[chromgroup]/germlineVariantCalls/` with filenames ending `.germlineVariantCalls.tsv` or `.germlineVariantCalls.vcf.bgz`.

- Metadata columns mirror the final-call tables (`analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `SBSindel_call_type`, `analysis_chunk`, `run_id`, `zm`).
- Germline filter flags without strand suffixes: `germline_vcf.passfilter`, any `germline_vcf_indel_region_filter_*`, `max_BAMVariantReads.passfilter`, `max_BAMVAF.passfilter`, and `region_read_filter_*` / `region_genome_filter_*` columns marked `is_germline_filter: true` in the YAML.
- Strand-specific columns reuse the `_refstrand_plus_read` / `_refstrand_minus_read` suffix pattern for `start_queryspace`, `end_queryspace`, `qual`, `qual.opposite_strand`, `sa`, `sm`, and `sx`, with each field storing comma-delimited per-base values. Columns such as `deletion.bothstrands.startendmatch` and `MDB_score` are removed before export so each TSV focuses on germline confirmation metrics.
- VCFs are bgzipped/indexed (`*.germlineVariantCalls.vcf.bgz` plus `.tbi`) and share the same INFO mappings as the TSVs.

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
- When sensitivity is borrowed from another chromgroup, the table records `sensitivity_source = "calculated_other_chromgroup"`. If no calculated sensitivity exists for that chromgroup/filtergroup pair, the value defaults to `1` with `sensitivity_source = "default_other_chromgroup"`.

### Final-call burdens
`[chromgroup]/finalCalls.burdens/[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].finalCalls.burdens.tsv`
originates from `calculateBurdensChromgroupFiltergroup` and reports burden statistics before/after sensitivity correction:

- Metadata: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`.
- Counts: `num_calls`, `interrogated_bases_or_bp`, `sensitivity`, `sensitivity_source`, `sensitivity_corrected` (boolean). A default sensitivity of `1` with `sensitivity_source = "default_other_chromgroup"` signals that no chromgroup-specific estimate was available.
- Rates and intervals: `burden_calls`, `burden_calls_lci`, `burden_calls_uci`, `num_calls_lci`, `num_calls_uci`.

### Estimated SBS mutation error probability
`[chromgroup]/estimatedSBSMutationErrorProbability/` now holds two TSV files per chromgroup/filtergroup combination:

- `[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].estimatedSBSMutationErrorProbability.by_channel_pyr.tsv` — Metadata columns plus `channel_pyr` and `error_prob` for each trinucleotide context.
- `[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].estimatedSBSMutationErrorProbability.total.tsv` — The same metadata paired with aggregate fields summarising the total SBS mutation error probability.

### Serialized object (.qs2)
The `.qs2` file captures the final data structures required to regenerate tables and plots without re-running the Nextflow processes. It is saved with the <a href="https://github.com/qsbase/qs2" target="_blank" rel="noopener noreferrer">qs2 R package</a> format and can be loaded directly with that library. Each component is listed below with its schema.

#### yaml.config
Top-level configuration loaded from the analysis YAML. Keys mirror the templates in [`config_templates/`](../config_templates),
so their meanings follow the descriptions in that documentation (for example `runs`, `samples`, `chromgroups`, `region_filters`,
and `sensitivity_parameters`). The list-column stores the original nested structure and string padding codes exactly as
specified in the YAML.

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
| `analyzein_chromgroups` | Present when the YAML specifies chromgroup restrictions for the call type. |
| MDB-specific fields | Columns such as `MDB_bamtag`, `MDB_min_score`, `MDB_sensitivity`, and `MDB_base_opposite_strand` appear when MDB call types are configured. |

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
| `num_genomebases_remaining` | Bases remaining after filtering. |

#### finalCalls

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type` | Call metadata. |
| `analysis_chunk`, `run_id`, `zm` | Chunk identifier, sequencing run, and ZMW hole number. |
| `seqnames`, `start_refspace`, `end_refspace` | Reference contig and 1-based coordinates. |
| `ref_plus_strand`, `alt_plus_strand` | Forward-strand reference and alternate alleles. |
| `refstrand` | Indicates which strand supplied the call (`+` or `-`). |
| `start_queryspace`, `end_queryspace` | Comma-delimited query-space coordinates for the reporting strand. |
| `qual`, `qual.opposite_strand` | Base qualities on the reporting and opposite strand (comma-delimited Phred values). |
| `sa`, `sa.opposite_strand` | Run-length encoded subread coverage counts per position. |
| `sm`, `sm.opposite_strand`; `sx`, `sx.opposite_strand` | Run-length encoded match and mismatch counts per position. |
| `call_class.opposite_strand`, `call_type.opposite_strand`, `alt_plus_strand.opposite_strand` | Allele annotations observed on the opposite strand. |
| `reftnc_plus_strand`, `alttnc_plus_strand`, `reftnc_pyr`, `alttnc_pyr` | Present only for SBS entries. |
| `indel_width` | Present only for indel entries. |
| `MDB_score` and related mismatch-damage metrics | Present only for MDB entries. |

#### finalCalls.bytype

| Column | Description |
| --- | --- |
| Metadata columns | Same identifiers and context fields as `finalCalls`. |
| `finalCalls_for_tsv`, `finalCalls_unique_for_tsv` | List-columns containing the strand-resolved and unique call tables used for TSV export. For SBS and indel `SBSindel_call_type = "mutation"`, the nested tables include strand-split columns with `_refstrand_plus_read` and `_refstrand_minus_read` suffixes for `start_queryspace`, `end_queryspace`, `qual`, `qual.opposite_strand`, `sa`, `sm`, and `sx`. |

#### germlineVariantCalls

| Column | Description |
| --- | --- |
| Metadata columns | Same as the strand-resolved `finalCalls` table described above. |
| Germline flags | `germline_vcf.passfilter`, `germline_vcf_indel_region_filter_*`, `max_BAMVariantReads.passfilter`, `max_BAMVAF.passfilter`, and any `region_*` filters marked as germline in the YAML. |
| Strand-specific fields | `_refstrand_plus_read` and `_refstrand_minus_read` columns mirroring the structure described for `finalCalls`. |

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
| `sensitivity`, `sensitivity_source` | Estimated detection rate and provenance (`calculated`, `calculated_other_chromgroup`, `default_other_chromgroup`, or YAML overrides such as `use_chromgroup`). When the provenance is `default_other_chromgroup`, the sensitivity value is `1`. |
| `high_confidence_germline_vcf_variants_sum_zm_detected`, `high_confidence_germline_vcf_variants_sum_duplex_coverage` | Counts of germline variants supporting the estimate. |

#### finalCalls.burdens

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type` | Identifiers. |
| `num_calls`, `interrogated_bases_or_bp` | Raw counts of calls and interrogated sequence. |
| `sensitivity`, `sensitivity_source`, `sensitivity_corrected` | Sensitivity estimate applied to the burden, including fallback designations such as `default_other_chromgroup`; when that fallback is used the sensitivity equals `1`. |
| `burden_calls`, `burden_calls_lci`, `burden_calls_uci`, `num_calls_lci`, `num_calls_uci` | Rates with confidence intervals. |
| `unique_calls`, `reftnc_corrected`, `reftnc_corrected_chromgroup` | Flags indicating whether corrections were applied. |

#### estimatedSBSMutationErrorProbability

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type` | Identifiers. |
| `by_channel_pyr` | List-column with `channel_pyr` and `error_prob` per trinucleotide. |
| `total` | Aggregate SBS mutation error probability. |
