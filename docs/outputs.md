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
| `${analysis_id}.${individual_id}.${sample_id}.run_metadata.tsv` | `analysis_id`, `individual_id`, `sample_id` identifier columns and sequencing run and read group metadata derived from BAM header fields (`rg_id`, `movie_id`, `SM`, `PM`, `PL`, `DS`, `PU`, etc.). |
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
| `region_genome_filter_stats.tsv` | Number of genome bases remaining after each genome region filter. Columns: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `filter` (filter label derived from the threshold file basename), `binsize` (per YAML configuration), `threshold` (per YAML configuration), `padding` (per YAML configuration), `region_filter_threshold_file` (per YAML configuration), `num_genomebases_individually_filtered` (number of genome bases removed when applying only this filter), `num_genomebases_remaining` (number of genome bases remaining after applying this and all prior filters). |

### Final calls
Final calls files are output to `[chromgroup]/finalCalls/` and are named:
`[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].[call_class].[call_type].[SBSindel_call_type].*`.

The pipeline produces four files per subtype:

- **`.finalCalls.tsv`** — Tables of the final calls that passed all filters.
  - Shared identifiers: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`, `analysis_chunk`, `run_id`, `zm`.
  - Coordinates in genome reference space and alleles: `seqnames`, `start_refspace`, `end_refspace`, `ref_plus_strand`, `alt_plus_strand`.
  - Strand-specific coordinates in query space (i.e., read coordinates) and other measurements: `start_queryspace`, `end_queryspace`, `qual`, `qual.opposite_strand`, `sa`, `sm`, `sx`. Measurements for multi-base calls, such as `qual` are stored as comma-delimited values for every base. SBS and indel mutation call tables are pivoted to one row per mutation with these strand-specific columns suffixed with `_refstrand_plus_read` and `_refstrand_minus_read`. Non-mutation call types remain one row per strand with a `refstrand` column annotating to which reference strand the call's read aligned.
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
Germline variant call files are output to `[chromgroup]/germlineVariantCalls/` and are named:
`[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].*`.

The pipeline produces two files:

- `.germlineVariantCalls.tsv` — Tables of the final germline variant calls that passed all filters, with one row per mutation and strand-specific columns suffixed with `_refstrand_plus_read` and `_refstrand_minus_read`.
  - Shared identifiers: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `analysis_chunk`, `run_id`, `zm`.
  - Coordinates in genome reference space and alleles: `seqnames`, `start_refspace`, `end_refspace`, `ref_plus_strand`, `alt_plus_strand`.
  - Strand-specific coordinates in query space (i.e., read coordinates) and other strand-specific measurements: `start_queryspace`, `end_queryspace`, `qual`, `qual.opposite_strand`, `sa`, `sm`, `sx`, suffixed with `_refstrand_plus_read` and `_refstrand_minus_read`.
  - Additional call information: `reftnc_plus_strand`,  `alttnc_plus_strand`, `reftnc_pyr`,  `alttnc_pyr`, and `indel_width`.
  - Germline filter information: `region_read_filter_*.passfilter` and `region_genome_filter_*.passfilter` columns (TRUE/FALSE) for each region filter configured with `is_germline_filter: true` in the YAML; `germline_vcf_types_detected` and `germline_vcf_files_detected` columns annotating in which germline VCF types and files the calls were detected.
- `.germlineVariantCalls.vcf.bgz`—  Bgzipped VCF containing the same calls as `.germlineVariantCalls.tsv`; INFO fields mirror TSV columns except those mapped to CHROM, POS, REF, and ALT. Indexed with `.tbi`.

### Final-call spectra
Trinucleotide distributions and call spectra of final calls are output to `[chromgroup]/finalCalls.spectra/` and are named:
`[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].[call_class].[call_type].[SBSindel_call_type].*`.

The pipeline produces several files per per below, and each table contains metadata columns (`analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`) and additional fields:

- **SBS** `call_class` files:

  - `.finalCalls.reftnc_pyr.tsv` / `.finalCalls_unique.reftnc_pyr.tsv` (trinucelotide distributions of reference sequences at positions of all and unique calls, respectively, collapsed to central pyrimidine contexts): `reftnc_pyr`, `count`, `fraction`.

  - `.finalCalls.reftnc_template_strand.tsv` (trinucelotide distributions of reference sequences in the strand replicated by the sequencer polymerase, at positions of all calls, for single-strand call types): `reftnc_template_strand`, `count`, `fraction`.

  - `.finalCalls.reftnc_pyr_spectrum.tsv` / `.finalCalls_unique.reftnc_pyr_spectrum.tsv` (reference>call trinucleotide spectra, collapsed to central pyrmidine contexts, of all and unique calls, respectively): `channel`, `count`, `fraction`.

  - `.finalCalls.reftnc_template_strand_spectrum.tsv` (reference>call trinucleotide spectra of the strand replicated by the sequencer polymerase, of all calls, for single-strand call types): `channel`, `count`, `fraction`.


- **indel** `call_class` files:
  - `.finalCalls.refindel_spectrum.sigfit.tsv` / `.finalCalls_unique.refindel_spectrum.sigfit.tsv` (indel spectra of all and unique calls, respectively, per sigfit-style channel labels): `channel`, `count`, `fraction`.


Each `_spectrum` file is paired with a `.pdf` spectrum plot generated via `plot_spectrum`.

### Interrogated-base spectra
Trinucleotide distributions of interrogated bases are output to `[chromgroup]/interrogatedBases.spectra/` and are named:
`[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].[call_class].[call_type].[SBSindel_call_type].*`.

The pipeline produces several files per per below, and each table contains metadata columns (`analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`) and additional fields:

- `.bam.gr.filtertrack.reftnc_pyr.tsv`: `reftnc_pyr`, `count`, `fraction`, `fraction_ratio_to_genome`, `fraction_ratio_to_genome_chromgroup` summarizing counts of interrogated base pairs (i.e., duplex coverage) together with ratios of the interrogated-base trinucleotide fractions to the whole-genome and chromgroup-restricted fractions.
- `.bam.gr.filtertrack.reftnc_both_strands.tsv`: `reftnc`, `count`, `fraction`, `fraction_ratio_to_genome`, `fraction_ratio_to_genome_chromgroup` summarizing counts of interrogated bases across both strands and the corresponding genome/chromgroup fraction ratios.

### Genome spectra
Trinucleotide distributions of the genome are output to `[chromgroup]/genome.spectra/` and are named:

- `genome.reftnc_pyr.tsv`: `analysis_id`, `individual_id`, `sample_id`, `reftnc_pyr`, `count`, `fraction`.
- `genome.reftnc_both_strands.tsv`: Same metadata with `reftnc`, `count`, `fraction`.
- `[chromgroup].genome_chromgroup.reftnc_pyr.tsv`: Same as `genome.reftnc_pyr.tsv`, restricting to `[chromgroup]`.
- `[chromgroup].genome_chromgroup.reftnc_both_strands.tsv`: Same as `genome.reftnc_both_strands.tsv`, restricting to `[chromgroup]`.

### Sensitivity summaries
Summaries of sensitivity analyses are output to `[chromgroup]/sensitivity/` and are named:
`[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].sensitivity.tsv` with columns:

- Shared identifiers: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`.
- Number of molecules detecting high-confidence germline variants: `high_confidence_germline_vcf_variants_sum_zm_detected`
- Duplex coverage of interrogated base pairs at sites of high-confidence germline variants: `high_confidence_germline_vcf_variants_sum_duplex_coverage`.
- `sensitivity`, : Sensitivity estimate, calculated as `high_confidence_germline_vcf_variants_sum_zm_detected` /  `high_confidence_germline_vcf_variants_sum_duplex_coverage` for homozygous variants, and 2 * this value for heterozygous variants.
- `sensitivity_source`: Provenance of the sensitivity value. Possible values are:
  - `calculated` — sensitivity derived from high-confidence germline variants detected in the current chromgroup/filtergroup.
  - `default` — sensitivity left at the default value of 1 because thresholds were not met or sensitivity analysis was not configured for that call type.
  - `yaml.config` — sensitivity sourced directly from YAML configuration (MDB call types).
  - `calculated_other_chromgroup` — sensitivity borrowed from another chromgroup where it was successfully calculated for the same filtergroup.
  - `default_other_chromgroup` — no calculated sensitivity was available to borrow, so the default value of 1 was used.

### Final-call burdens
Final calls are output to `[chromgroup]/finalCalls.burdens/` and are named:
`[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].finalCalls.burdens.tsv`, with one row per type of burden calculation (see below), with columns:

- Shared identifiers: `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`.
- Annotations for type of burden calculation:
  - Unique calls (boolean TRUE/FALSE for mutations, NA otherwise): `unique_calls`
  - Corrected for the ratio of the trinucleotide distribution of interrogated bases or base pairs to the trinucleotide distribution of the whole genome (boolean TRUE/FALSE for `call_class` SBS and MDB when `unique_calls` is FALSE and for `call_class` SBS when `unique_calls` is TRUE, NA otherwise): `reftnc_corrected`
  - Baseline used for trinucleotide correction when `reftnc_corrected` is TRUE: `reftnc_corrected_chromgroup` (`genome` for whole-genome corrections, `genome_chromgroup` for chromgroup-restricted corrections, NA otherwise).
  - Corrected for sensitivity: `sensitivity_corrected` (boolean TRUE/FALSE).

- Number of calls corrected for above selected correction methods, and Poisson 95% lower (lci) and upper (uci) confidence intervals: `num_calls`,  `num_calls_lci`, `num_calls_uci`.
- Pre-correction call counts retained alongside corrected values (only populated for trinucleotide- and sensitivity-corrected rows): `num_calls_noncorrected`.
- Number of interrogated bases (single-strand calls) or base pairs (double-strand calls): `interrogated_bases_or_bp`.
- Sensitivity estimate and source (from the sensitivity summary table, annotated when `sensitivity_corrected` is TRUE, NA otherwise): `sensitivity`, `sensitivity_source`.
- Burden (`num_calls` / `interrogated_bases_or_bp`) and Poisson 95% lower (lci) and upper (uci) confidence intervals: `burden_calls`, `burden_calls_lci`, `burden_calls_uci`.

### Estimated SBS mutation error probability
Estimated SBS mutation error probabilities are output to `[chromgroup]/estimatedSBSMutationErrorProbability/` and are named:
`[analysis_id].[individual_id].[sample_id].[chromgroup].[filtergroup].estimatedSBSMutationErrorProbability.*`. 

The pipeline produces several files per per below, and each table contains metadata columns (`analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`) and additional fields:

- `by_channel_pyr.tsv` — `channel_pyr` (reference>call central pyrimidine trinucleotide channel) and `error_prob` (error probability).
- `total.tsv` — Total SBS mutation error probability across all trinucleotide channels.

### Serialized object (.qs2)
For each sample, a `.qs2` file is stored in `[analysis_output_dir]/[analysis_id].[individual_id].[sample_id]/outputResults/${analysis_id}.${individual_id}.${sample_id}.qs2`. It contains all the final data structures required for downstream analyses in R in a more convenient bundle than the many separate files described above. It is saved with the <a href="https://github.com/qsbase/qs2" target="_blank" rel="noopener noreferrer">qs2 R package</a> format and can be loaded with `qs2::qs_read`. Each component of the object is listed below with its schema.

#### yaml.config
Top-level configuration loaded from the analysis YAML, stored as nested lists per the structure and parameters documented in [`config_templates/`](../config_templates).

#### run_metadata
Columns as per `run_metadata.tsv` described above. 

#### individual and sample identifiers

| Field | Description |
| --- | --- |
| `individual_id` | String naming the individual analysed. |
| `sample_id` | Scalar string naming the sample. |

#### molecule_stats.by_run_id,  molecule_stats.by_analysis_id

Columns as per `molecule_stats.by_run_id.tsv` and `molecule_stats.by_analysis_id.tsv` described above. 

#### region_genome_filter_stats

Columns as per `region_genome_filter_stats.tsv` described above.

#### finalCalls

Unified table of all call types across all chromgroups and filtergroups, with one row per call at the strand level (i.e. each mutation has two rows, one for each strand)

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type` | Call metadata. |
| `analysis_chunk`, `run_id`, `zm` | Chunk identifier, sequencing run ID, and ZMW hole number. |
| `seqnames`, `start_refspace`, `end_refspace` | Reference contig and reference space 1-based coordinates. |
| `ref_plus_strand`, `alt_plus_strand` | Reference genome forward strand reference and alternate allele sequences. |
| `refstrand` | Reference genome strand that the read aligned to which supplied the call (`+` or `-`). |
| `start_queryspace`, `end_queryspace` | Query-space 1-based coordinates. |
| `qual`, `qual.opposite_strand` | Base qualities of the call and on the opposite strand per reference space alignment coordinates (Phred values; list of values for each call, as some calls span multiple bases). |
| `sa`, `sa.opposite_strand` | Subread coverage of the call and on the opposite strand per reference space alignment coordinates (list of values for each call, as some calls span multiple bases). |
| `sm`, `sm.opposite_strand`; `sx`, `sx.opposite_strand` | Subread match and mismatch counts of the call and on the opposite strand per reference space alignment coordinates (list of values for each call, as some calls span multiple bases). |
| `call_class.opposite_strand`, `call_type.opposite_strand`, `alt_plus_strand.opposite_strand` | Call annotations observed on the opposite strand, if any. |
| `reftnc_plus_strand`, `alttnc_plus_strand`, `reftnc_pyr`, `alttnc_pyr` | Present only for SBS and MDB entries. |
| `indel_width` | Present only for indel entries. |
| `MDB_score` | Present only for MDB entries. |

#### finalCalls.bytype

Contains final calls split into one row for each combination of `call_class`, `call_type`, and `SBSindel_call_type`. For each row, there are metadata identifiers (`analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`) and two tables:

| Table | Description |
| --- | --- |
| `finalCalls_for_tsv` | Same format as `finalCalls.tsv` described above. |
| `finalCalls_unique_for_tsv` | Same format as `finalCalls_unique.tsv` described above. |

#### germlineVariantCalls

Germline variant calls with the same schema as the unified finalCalls table.

#### finalCalls.reftnc_spectra

Table containing trinucleotide distributions and spectra of final calls, split into one row for each combination of `call_class`, `call_type`, and `SBSindel_call_type`. For each row, there are metadata identifiers (`analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`) and the following tables:

| Column | Description |
| --- | --- |
| `finalCalls.reftnc_pyr`, `finalCalls_unique.reftnc_pyr` | Same format as `finalCalls.reftnc_pyr.tsv` and `finalCalls_unique.reftnc_pyr.tsv` described above. |
| `finalCalls.reftnc_template_strand` | Same format as `finalCalls.reftnc_template_strand.tsv` described above. |
| `finalCalls.reftnc_pyr_spectrum`, `finalCalls_unique.reftnc_pyr_spectrum` | Same format as `finalCalls.reftnc_pyr_spectrum.tsv` and `finalCalls_unique.reftnc_pyr_spectrum.tsv` described above. |
| `finalCalls.reftnc_template_strand_spectrum` | Same format as `finalCalls.reftnc_template_strand_spectrum.tsv` described above. |
| `finalCalls.refindel_spectrum`, `finalCalls_unique.refindel_spectrum` | Indel spectra, with channel labels per the [indelwald package](https://github.com/MaximilianStammnitz/Indelwald) (created with the indel.spectrum() function in the HiDEF-seq bin/sharedFunctions.R script). |
| `finalCalls.*.sigfit` columns | Sigfit-formatted matrices generated from the corresponding spectra. |

#### bam.gr.filtertrack.bytype.coverage_tnc

Table containing genome coverage by interrogated base pairs (i.e., duplex coverage) and trinucleotide distributions of interrogated bases, split into one row for each combination of `call_class`, `call_type`, and `SBSindel_call_type`. For each row, there are metadata identifiers (`analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`) and the following tables:

| Column | Description |
| --- | --- |
| `bam.gr.filtertrack.coverage` | List-column of GRanges coverage objects for interrogated base pairs (i.e., duplex coverage counting once coverage by both strands of a molecule). |
| `bam.gr.filtertrack.reftnc_pyr`, `bam.gr.filtertrack.reftnc_both_strands` | Same format as `.bam.gr.filtertrack.reftnc_pyr.tsv` and `bam.gr.filtertrack.reftnc_both_strands.tsv` described above. |

#### genome.reftnc and genome_chromgroup.reftnc

`genome.reftnc` and `genome_chromgroup.reftnc` contain a metadata identifier `chromgroup` and the following tables:

- `reftnc_pyr`: Same format as `genome.reftnc_pyr.tsv` and `[chromgroup].genome_chromgroup.reftnc_pyr.tsv`.
- `reftnc_both_strands`: Same format as `genome.reftnc_both_strands.tsv` and `[chromgroup].genome_chromgroup.reftnc_both_strands.tsv`.

#### sensitivity

| Column | Description |
| --- | --- |
| `analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type` | Identifiers. |
| `sensitivity`, `sensitivity_source` | Estimated detection rate and provenance (`calculated`, `calculated_other_chromgroup`, `default_other_chromgroup`, or YAML overrides such as `use_chromgroup`). When the provenance is `default_other_chromgroup`, the sensitivity value is `1`. |
| `high_confidence_germline_vcf_variants_sum_zm_detected`, `high_confidence_germline_vcf_variants_sum_duplex_coverage` | Counts of germline variants supporting the estimate. |

#### finalCalls.burdens

Columns as per `.finalCalls.burdens.tsv` described above, in one table for all call types.

#### estimatedSBSMutationErrorProbability

Table containing estimated SBS mutation error probabilities, split into one row for each combination of `call_class`, `call_type`, and `SBSindel_call_type`. For each row, there are metadata identifiers (`analysis_id`, `individual_id`, `sample_id`, `chromgroup`, `filtergroup`, `call_class`, `call_type`, `SBSindel_call_type`) and the following tables:

| Column           | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| `by_channel_pyr` | Same format as `estimatedSBSMutationErrorProbability.*.channel_pyr.tsv` described above. |
| `total`          | Same format as `.estimatedSBSMutationErrorProbability.*.total.tsv` described above. |
