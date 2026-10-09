# Optimization measurement and scientific validation

Run large QS2 reads, writes, and comparisons **inside an allocated compute job**.
The host need not have R or GNU time; the pipeline container may provide both.
The measurement wrapper falls back to Linux Python child-process accounting when
GNU time is unavailable.
All examples write below the working `codex` directory. These tools do not submit
jobs or change pipeline output formats. Scratch component files used by validation
are disposable and are never published pipeline outputs.

Benchmark and validation drivers use R, Python and shell separately from
production code. The [implementation table](../../docs/workflow-optimization.md#implementation-languages)
distinguishes custom production helpers, orchestration, existing CLI tools and
compiled package internals; historical C++ reference arms are identified below.

## CPU accounting

The primary performance measure is total **actual user + system CPU time**, including
waited-for child processes. Wall time, allocated CPU hours, and peak RSS are separate
measures. Allocated CPU hours are not actual CPU consumed. There is no per-step 5%
acceptance threshold: evaluate total representative-workload CPU and inspect stage
regressions in context. Keep identical input/configuration/thread/resource settings,
record commit/container, and alternate reference/candidate runs when practical.

Wrap the actual foreground task, not `sbatch`, an asynchronous launcher, or the
Nextflow driver when its tasks execute elsewhere:

```sh
python3 scripts/benchmark/measure.py --label reference-outputResults \
  --out /projects/work/evrong01/HiDEF-seq/codex/bench/reference.outputResults \
  --allocated-cpus 4 -- Rscript --vanilla bin/outputResults.R [pipeline arguments]
```

The wrapper invokes GNU `/usr/bin/time` (override with `--time`) and writes JSON,
TSV, and raw timing output. When absent, it measures the delta of Python
`resource.getrusage(RUSAGE_CHILDREN)` around the foreground command and records this
method in the report. Linux child resource usage propagates from descendants that
their parents wait for; the same detached/remote-child limitation applies. The
command's output remains attached to the job's
stdout/stderr. Failed commands retain measurements and propagate their exit code.
Both methods count children only when the workload waits for them; detached or remote
jobs are not covered. Maximum RSS is GNU time's process/descendant high-water metric,
**not** a sum of simultaneous memory across a process tree. The profiler below
measures R's own phase high-water RSS.

After a Slurm job finishes, retain allocation and step accounting separately:

```sh
python3 scripts/benchmark/collect_sacct.py --jobs JOBID \
  --out /projects/work/evrong01/HiDEF-seq/codex/bench/reference.sacct
```

Use the allocation record for allocated CPU hours and appropriate completed task
or `.batch` records for actual CPU checks. **Do not sum an allocation and its steps.**
`TotalCPU` is converted to seconds; `CPUTimeRAW / 3600` is reported separately.
Accounting can lag completion and RSS fields may be blank on allocation records.
Preserve records rather than filling missing values with zero. Sum task-level GNU
time CPU over the same executed workflow scope, including retries; do not compare
a cached/resumed run against a complete run.

For precise standalone-job comparisons, also retain `sacct --json -j JOBID` and
read the normal allocation's CPU fields as integer microseconds:
`seconds * 1000000 + microseconds`. Slurm may return an unnormalized microsecond
field greater than one million; that is valid and must not be rejected or
counted twice. Verify that exact user plus system equals total, and keep step
memory records separate. The later storage and output-replay controls use this
check. The collectors above and frozen full-workflow reports retain their
original displayed-field precision; this documentation does not change them.

`summarize_workflow.py` joins one or more Nextflow traces to saved accounting.
It deduplicates original job IDs across resumed traces and uses allocation CPU
totals without adding their `.batch` or other steps a second time. Failed attempts
remain in total expenditure and are separated from successful-workload CPU.
Missing accounting stays unresolved rather than becoming zero. For example:

```sh
job_ids=$(python3 scripts/benchmark/summarize_workflow.py --trace TRACE.tsv --job-ids-only)
python3 scripts/benchmark/collect_sacct.py --jobs "$job_ids" --out ACCOUNTING_PREFIX
python3 scripts/benchmark/summarize_workflow.py --trace TRACE.tsv \
  --accounting ACCOUNTING_PREFIX.json --preparation-scope cold --report WORKFLOW_METRICS.json
```

Repeat `--trace` to include earlier attempts and `--accounting` to supply later
accounting snapshots. By default, cached task records retain their original
measured workload cost. `--exclude-cached` requires exactly one supplied trace so
that it measures recorded new work in a single invocation; an executed record in
any supplied trace otherwise counts once, even if a later invocation cached it.

Classify preparation explicitly with `--preparation-scope cold|verified-cache|historical-unverified|mixed|unknown`
(default `unknown`). Preparation reused by `storeDir` without trace records and
historical cache construction have no measured cost here. The report separates
preparation from common downstream CPU, and successful from unsuccessful attempts
within each scope. Total observed CPU includes failed retries and incomplete work;
success-only CPU is an additional diagnostic, not total expenditure. A baseline
importing historical caches and a cold candidate do not establish a comparable
end-to-end cost: report preparation separately, compare matched downstream work,
and account for preparation independently. No scope label alone establishes
comparability; reports explicitly leave that conclusion unasserted.

Reports remain partial unless `--workflow-complete` is supplied after verifying
termination and the recorded accounting is resolved. That flag does not establish
scientific correctness, cache equivalence, or completeness of unrecorded work.
Signal-interrupted jobs/steps remain CPU-incomplete even when numeric `TotalCPU`
exists: [Slurm documents](https://slurm.schedmd.com/sacct.html) that interrupted
steps may omit child-process CPU. Their observed cost is retained as a lower bound.
An allocation reporting less CPU than the sum of its deduplicated steps, or an
allocation with still-live steps, also remains unresolved. The consistency check
allows Slurm's loss of subsecond precision when displaying CPU durations of an
hour or longer; it never adds step CPU to the allocation total. Later accounting
snapshots can resolve
ordinary lag; evidence of signal interruption is retained across snapshots.
Peak RSS is reported per task and stage, not summed across tasks. Wall time
requires the workflow/controller records.

## Profile one existing final QS2

```sh
Rscript --vanilla scripts/benchmark/profile_qs2.R FINAL.outputResults.qs2 \
  /projects/work/evrong01/HiDEF-seq/codex/bench/LIB1.profile \
  /projects/work/evrong01/HiDEF-seq/codex/bench/LIB1.roundtrip.qs2
```

The final argument is optional. The script loads exactly one object, reports logical
sizes of top-level components and data-frame columns, and optionally times a save
of that object to one new QS2 file. CPU and wall timings are separate for load,
component inspection, and save. RSS is recorded before/after each phase and at its
peak. Linux `clear_refs=5` resets the RSS high-water counter before each phase when
permitted; `peak_is_phase_specific=FALSE` explicitly marks a cumulative peak when
resetting fails. Save-phase memory includes the already loaded object. Logical
`object.size()` can double-count shared storage and excludes some external storage;
component sizes are neither compressed sizes nor process RSS. The complete prior
`outputResults` task's peak memory is **not** the size of the final object.

## Compare scientific R objects without loading both full objects

```sh
Rscript --vanilla scripts/benchmark/compare_qs2.R stage REFERENCE.outputResults.qs2 \
  /projects/work/evrong01/HiDEF-seq/codex/bench/reference-components
Rscript --vanilla scripts/benchmark/compare_qs2.R compare \
  /projects/work/evrong01/HiDEF-seq/codex/bench/reference-components \
  CANDIDATE.outputResults.qs2 \
  /projects/work/evrong01/HiDEF-seq/codex/bench/objects.tsv --chromgroup-blocks
```

Reference staging reads one complete QS2, then saves/releases one top-level
component at a time. Comparing holds the candidate plus one reference component,
releasing both compared components after each step. Budget RAM for the candidate
plus its largest reference component and serialization/attribute overhead, and
scratch space for staged data. This is not a streaming QS2 reader.

Comparison preserves storage types, classes, all attributes (including factor
levels, dimensions, dimnames, and names), list order, and S4 slots such as Rle and
GRanges/Seqinfo. Root schema must match. The only default metadata whitelist is
`/run_metadata`. Add exact paths, for example `--ignore=/yaml.config/work_dir`, only
after verifying that the field is run-specific. Paths use slash-separated names
with JSON-pointer escaping; positional list entries use `[[N]]`. Every whitelist
hit is reported. No whole-configuration whitelist is supplied.

For final saved configurations, `--config-policy=/absolute/path/policy.json`
binds the already loaded `yaml.config` component to reviewed effective YAML on
each side before applying metadata exclusions; it does not load either QS again.
The component must exist. `config_helper_bindings.ready` must be true, and its
`reference` and `candidate` entries each specify `effective_config`, the file's
`sha256`, and the exact `helper_keys` set. The shared
`compare_config_bindings.R` helper checks every helper value (including nested
thresholds and artifact paths) plus all scientific configuration fields.
`config_added_metadata` names the permitted helper fields; listing a field does
not exempt its value from the binding check. Missing helpers, changed targets,
stale YAML hashes, or scientific drift fail even if a corresponding `--ignore`
path is present. `tests/test_config_bindings.R BENCHMARK_DIRECTORY OUTPUT_DIRECTORY`
exercises these checks and the actual disk-staged QS interface on tiny fixtures.

`--token-policy=/absolute/path/rules.json` permits a value-checked source-path
relocation only in explicitly named top-level table columns named
`germline_vcf_files_detected`. For example:

```json
{"columns":{"/germlineVariantCalls/germline_vcf_files_detected":{
  "reference":["/old/a.vcf.gz","/old/b.vcf.gz"],
  "candidate":["/new/a.vcf.gz","/new/b.vcf.gz"]}}}
```

The paired arrays define complete comma-delimited tokens. Every nonempty token
must occur in its side's array; unknown paths fail. Mapping preserves token order,
duplicate multiplicity, empty strings, NA, character-vector attributes and all
other table data. It does not replace substrings or assembly labels. Add separate
exact paths for other tables only after inspecting their actual provenance.
The policy also supports the one reviewed derived-cache path column
`/region_genome_filter_stats/region_filter_threshold_file`, with
`"mode":"exact-path"` and paired `reference`/`candidate` arrays. These entries
map complete field values, without comma splitting or substring replacement.
Build the mapping from independently verified region-filter products and inspect
the observed values first. Unknown paths, wrong products, changed NA positions,
attributes, row order, and all other table values still fail. No other cache-path
column is authorized by this mode.

The policy changes validation in memory; source QS files remain untouched.
Run its small regression fixture from the repository root:

```sh
Rscript --vanilla tests/test_provenance_tokens.R
Rscript --vanilla tests/test_provenance_paths.R .
```

`--chromgroup-blocks` permits the known cross-chromgroup combined-table ordering
variation by stably grouping data-frame rows on `chromgroup`. It preserves order
within each chromgroup, including filtergroup order, and does not sort scientific
rows within groups. Automatic row names are regenerated; explicit row names are
compared. Omit the flag for entirely strict order checks.

Exit codes: **0** scientific match (including explicitly allowed metadata changes),
**1** discrepancy, **2** floating-point investigation required. Changed doubles at
symmetric relative error `abs(a-b)/max(abs(a),abs(b)) <= 1e-12` are review findings,
never automatic passes. Zero receives no absolute-tolerance exception. NA, NaN,
and signed infinity must match. Published numeric precision must still match
exactly through the text comparator. Reports cap detail at 1000 rows per component
but count all findings. The CLI is designed for pipeline top-level lists; the
sourceable `scientific_compare()` also supports focused intermediate-object tests.

## Compare published files

```sh
python3 scripts/benchmark/compare_outputs.py REFERENCE_DIR CANDIDATE_DIR \
  --report /projects/work/evrong01/HiDEF-seq/codex/bench/published.json
```

The comparator inventories all relative paths and checks byte identity. For
TSV/CSV/text/BED/VCF/YAML (including gzip), an exact decompressed-byte fast path
avoids iterating through billions of matching BED lines. It checks gzip CRCs and
strict UTF-8 across block boundaries. Byte differences fall back to line-by-line
comparison without numerical tolerance or sorting. Compression and line-ending
differences are allowed. It excludes only
`*.run_metadata.tsv` and VCF `##fileDate=` by default; additional explicit metadata
globs/VCF keys are recorded in the report. Configuration files remain strict unless
specific run metadata is separately reviewed. Differing QS2, PDF, index, and other
binary files are reported as **review required**, not passed: resolve QS2 findings
with the R comparator, validate VCF indexing separately, and inspect/raster-compare
PDFs as needed. The report is a validation inventory, not proof that unexamined
binary differences are harmless. Its exit codes are 0 match, 1 discrepancy, 2 review.

Small harness checks (R/Bioconductor checks use the pipeline environment):

```sh
python3 -m unittest discover -s tests -p test_validation_harness.py -v
Rscript --vanilla tests/test_scientific_compare.R
```

### PDF content and binary indices

`compare_pdf.py` provides a conservative format-specific check for changed PDFs.
The pinned production container had no qpdf, Poppler, Ghostscript, Python PDF
parser, or R qpdf/pdftools package. Validation therefore uses **pypdf 6.19.0** and
**typing_extensions 4.15.0** in a separate workspace environment; these are not
production dependencies. Hash-pinned wheels are listed in
`requirements-pdf-validation.txt`. The current workspace installation and verified
download manifest are under `runtime/pdf-validation/`. The parser's strict reader
API is documented by [pypdf](https://pypdf.readthedocs.io/en/latest/modules/PdfReader.html).

```sh
PYTHONPATH=/projects/work/evrong01/HiDEF-seq/codex/runtime/pdf-validation/site-packages \
  python3 scripts/benchmark/compare_pdf.py REFERENCE.pdf CANDIDATE.pdf --report pdf.json
PYTHONPATH=/projects/work/evrong01/HiDEF-seq/codex/runtime/pdf-validation/site-packages \
  python3 -m unittest discover -s tests -p test_pdf_comparison.py -v
```

The check compares the rooted parsed PDF object graph, preserving indirect-object
sharing, arrays, dictionary values, page geometry and font attributes. It hashes
every decoded reachable stream exactly, including drawing commands, text and
embedded font programs. Object numbers, dictionary-key order and supported stream
compression encoding may differ. The pinned reader retains exact original decimal
values before conversion to binary floats, so small coordinate differences cannot
silently disappear. It is not a raster approximation or a general PDF conformance
validator, and does not compare obsolete/unreachable objects outside the rooted
document graph.

Only `/Info/CreationDate` and `/Info/ModDate` are ignored by default, and every
ignored key is reported. Additional proven provenance-only Info fields require
an explicit `--ignore-info-key EXACT_KEY`; a document identifier requires
`--ignore-document-id`. These options never replace paths, dates or text inside
page content, fonts or XMP streams. Non-whitelisted differences fail. Unsupported
filters (anything beyond unfiltered/Flate streams), encryption, detected
incremental revisions, active/interactive features and parser warnings require
review; no broad byte-regex normalization is used. Exit codes are 0 pass, 1
difference, 2 review. The general output inventory still reports differing PDFs
for review; retain the separate PDF reports as evidence resolving those entries.
Five adversarial fixtures cover metadata, compression, drawing/text/font changes,
indirect sharing and decimal precision. Compute job `19194972` additionally parsed
three historical PDFs with compressed object/XRef streams: date-only rewrites
preserving the PDF header passed, and all three changed-drawing versions failed.
Those validation artifacts are in workspace `runs/pdf-real-validation-v2/`.

`compare_bam.py` compares raw decompressed BAM record blocks in order, including
every tag byte and float bit. It does not reserialize through SAM or a BAM writer,
so it cannot hide differences through text precision or long-CIGAR transport
normalization. BGZF compression may differ. Header fields/order and the binary
reference dictionary remain exact, except the explicit `CL` field of `@PG`
records. Other `@PG` provenance/count/order differences require manual review;
they are never silently removed. Other header or raw-record differences fail.
The framing follows the [SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
`samtools quickcheck` checks header/EOF before streaming, gzip validates CRCs,
and malformed/truncated frames fail. Frames over 64 MiB or oversized reference
dictionaries require review rather than unbounded allocation.

```sh
python3 scripts/benchmark/compare_bam.py REFERENCE.bam CANDIDATE.bam \
  --output NEW_SCRATCH_DIRECTORY --samtools samtools --pbindex pbindex --rebuild-indexes
python3 -m unittest discover -s tests -p test_bam_comparison.py -v
python3 tests/test_bam_comparison_tools.py --output NEW_TINY_FIXTURE_DIRECTORY \
  --samtools samtools --pbindex pbindex --bgzip bgzip
```

The optional index check rebuilds BAI and PBI against **each BAM's own scratch
symlink**, then compares that BAM's published index with its rebuilt index (PBI
decompressed, BAI bytes). It never compares offsets across differently compressed
BAMs. Input indices must remain untouched; file metadata is checked before and
after. Mismatches require review because they can reflect an invalid index or a
different valid index encoding. Missing indices or rebuild failures fail.
Without `--rebuild-indexes`, BAM contents can pass but the overall report remains
review-required for unchecked indices. Exit codes are 0 pass, 1 failure, 2 review.
Scratch rebuilding adds full BAM reads and must run on compute **outside measured
pipeline work**. No full BAM validation job is launched by the utility itself.
Pinned-tool fixture job `19196140` passed: both BAI/PBI scratch rebuilds matched,
source BAMs/indices remained byte-identical, reordered/changed records and a
truncated BAM failed, and a stale PBI required review. A one-ULP float-tag mutation
rendered identically in SAM but was correctly rejected by raw-record comparison.
Artifacts and frozen source hashes are in workspace `runs/bam-comparison-fixtures-v2/`.

The original ordered comparator remains unchanged for reproducibility. The
completed Torch comparison also has an explicitly approved coordinate-tie rule:
all logical records, multiplicities, coordinates and coordinate ordering must
match, but records at identical coordinates may have a different relative order.
Full raw-record diagnostics checked that rule for both processed BAMs and rebuilt
all eight own-file BAI/PBI indexes. The original failed ordered report is retained.

`accept_bam_coordinate_diagnoses.py` validates the pinned supplemental evidence
for that completed two-sample workload. It accepts only its two exact BAM paths
and four owner-linked index paths; all 705 other required scientific paths must
already pass. This is an evidence gate, not a replacement BAM scanner. It requires
the explicit approval, pinned original reports and diagnostic code, exact input
bindings, zero record/coordinate differences, complete successful index checks,
input timestamps predating diagnosis, and candidate identity continuity with the
frozen producer audit. It does not relax any record field or numeric comparison.

```sh
python3 scripts/benchmark/accept_bam_coordinate_diagnoses.py PINNED_PLAN.json \
  --output NEW_ACCEPTANCE_REPORT.json
python3 -m unittest discover -s tests -p test_bam_coordinate_acceptance.py -v
```

The retained plan/report live under workspace
`runs/bam-coordinate-acceptance-3911047/`; the full diagnostic code and reports are
under `runs/bam-coordinate-diagnosis-lib1/` and `runs/bam-coordinate-diagnosis-lib2/`.
Twenty-one synthetic acceptance and rejection tests cover altered records,
missing/failed indexes, wrong files/owners, incomplete evidence, changed inputs,
and extra scientific failures. See the
[validation ledger](../../docs/optimization-validation.md) for final acceptance
and the unchanged original failure counts.

The BED harness already exercises Tabix boundary/empty-region queries. Matching
compressed scientific contents alone does not prove a changed index is valid;
retain these format-specific reports alongside the published-file inventory.

## Effective YAML verification

`compare_config.R ORIGINAL.yaml EFFECTIVE.yaml REPORT.tsv` parses both files using
the pipeline's `configr::read.config` and compares every original key and nested
value. It permits only the four explicit added artifact metadata fields used by
the workflow; unknown additions fail. Original run-metadata differences require
an explicit `--metadata-key=analysis_output_dir` (or another exact key).
Additional permitted new keys require `--allow-extra=EXACT_KEY`. Scientific
parameters are never covered by a blanket configuration exclusion.

## Experimental combined filtering benchmark

The prototype uses custom R and existing package internals; production adoption
would require Nextflow/shell changes. Python controls the benchmark only.

This experiment is **not connected to the production workflow**. It compares the
current filtering script in separate R processes with the same script evaluated
in isolated per-group environments inside one R session. The latter retains one
extraction QS object and intercepts only reads of that exact input path. Helpers
are sourced into each group's environment so their lexical scope matches the
normal standalone script. Each group has its own working directory and complete
QS output. Package initialization and the extraction read are shared, so measure
the full chain; do not attribute the entire difference to file loading alone.

Both arms require already prepared reference-summary and relevant germline
coverage-cache files, avoiding confounding by repeated preparation work. The
driver freezes R source and config, alternates arm order across paired repeats,
and records complete-chain actual CPU (including children), wall time, allocated
CPU hours and peak RSS through `measure.py`. Shared extraction residency may
increase memory. Complete scientific QS comparisons run outside timed regions;
any difference or floating-point review stops the experiment. Scratch reference
components and all outputs are retained for review. No output sharding or
production filtering behavior is changed.

```sh
python3 scripts/benchmark/benchmark_filter_chain.py \
  --config PREPARED_CONFIG.yaml --extract extractCalls.chunk1.qs2 \
  --sample 70007.1-CordBlood-LIB1 --chromgroup 1-22X \
  --filtergroups lenient strict --allocated-cpus 2 --repeats 3 \
  --output /projects/work/evrong01/HiDEF-seq/codex/runs/filter-chain-experiment
Rscript --vanilla tests/test_filter_chain_harness.R
```

The experiment retains both arms' intermediates for diagnosing differences;
only complete filtering QS objects are scientific outputs of these tasks.
Whole-process RSS is a process high-water metric, not summed simultaneous
memory across all child processes. Allocation/accounting for the complete job,
including validation overhead, can additionally be recorded with `collect_sacct.py`.

## Paired complete extraction or filtering

The measured extraction/filtering changes are custom R using existing packages
and CLI tools. The Python replay driver adds no production scientific code.

`benchmark_extraction.py` runs alternating pairs of independent R processes in
one allocation. Supply `--bam INPUT.bam` for extraction, or
`--extract-qs INPUT.qs2 --chromgroup GROUP --filtergroup GROUP` for filtering.
Both modes require `--baseline-repo`, `--candidate-repo`, `--config`, `--sample`,
`--output` (a new directory), and the appropriate `--allocated-cpus` value.
The default is two pairs. Source/configuration snapshots, input identities and
per-process metrics are retained; full scientific QS comparisons run afterward,
outside the timed operations.

Filtering gives both arms the same prepared configuration and the exact relative
input name `extractCalls.chunk1.qs2`, since the script records that argument in its
output. The original script ignores the new preparation fields and performs its
usual inline work. The comparison retains every configuration field; only the
normal `/run_metadata` exception applies. Run on compute with enough memory for
the original operation as well as the subsequent comparison.

## Germline VCF annotation block benchmark

The matching and summary optimization is custom R with existing package
internals and the loader's existing shell/`bcftools` pipeline. The block harness
is R; the optional Python measurement wrapper is benchmark infrastructure.

`benchmark_germline_vcf.R` extracts the actual germline VCF block from original
and candidate `filterCalls.R` scripts. It executes the candidate's preceding
filters once on a real extraction QS, then gives both arms the same resulting
calls and configuration. Each timed arm includes the complete VCF read,
conversion, quality filtering, deduplication, summary, and annotation join. The
candidate's additional matching scan is inside this timed region. The original
loader is redirected by basename to the same prepared artifact, without changing
its scientific expressions. No VCF read is cached between arms.

```sh
Rscript --vanilla tests/test_germline_vcf_annotation.R
Rscript --vanilla scripts/benchmark/benchmark_germline_vcf.R \
  BASELINE/bin/filterCalls.R bin/filterCalls.R PREPARED_CONFIG.yaml \
  extractCalls.chunk1.qs2 70007.1-CordBlood-LIB1 1-22X lenient NEW_OUTPUT_DIR 3
```

Run this on a compute node. It alternates arm order, reports actual CPU including
waited-for children, wall time and R phase RSS, and records whether the RSS
high-water reset succeeded. RSS includes the retained common preceding-filter
context. Arms use separate evaluation environments; their output objects are
released and garbage-collected before the next arm. They still share one R
process: allocator-resident pages and filesystem caches can persist. Alternating
order mitigates this effect but does not make the RSS values independent
fresh-process peaks. Use full-filter process replays for that claim.
Each pair compares both annotated calls and the complete filtered
germline table with `identical()` outside timing; the latter must remain intact
for downstream whole-genome statistics. Source ASTs, source/input fingerprints,
outputs, metrics and session information are retained. Wrap the command with
`measure.py` for the entire experiment's resource cost, keeping setup/validation
cost separate from the per-block results. No full-pipeline gain is inferred from
this focused benchmark.

## Germline VCF export allocation profile

Run `profile_germline_vcf.py` inside the pinned container on compute (the real
LIB1 profile requests 96 GiB, two CPUs and two hours):

```sh
python3 scripts/benchmark/profile_germline_vcf.py \
  --repo . --input ORIGINAL.outputResults.qs2 \
  --reference-vcf ORIGINAL.1-22X.lenient.germlineVariantCalls.vcf.bgz \
  --chromgroup 1-22X --filtergroup lenient --output NEW_PROFILE_DIRECTORY
```

The driver freezes the actual writer/normalizer source and harness with SHA256
hashes, and records the input path and size, installed writer methods, package version and
session details. A lightweight compute preflight precedes the large QS load.
The worker retains **all** original final QS components while exporting the
selected formatted table, including raw germline and coverage. This deliberately
also retains later-group formatted/final results that may not yet exist at the
first real export; the payload approximation is recorded explicitly.

Expression-entry and function-exit markers report cumulative CPU, elapsed time,
current RSS and process HWM. They do not retain result snapshots, force garbage
collection, reset HWM or change the native writer buffer. CPU differences between
nested events are inclusive and must not be summed across overlapping phases.
These instrumented diagnostic timings are separate from accepted performance
benchmarks. Whole-worker metrics include input loading; `export:begin/end` delimit
the complete normalization and export operation.

After the worker exits, the driver compares every decompressed VCF scientific
line with the published reference, ignoring only `fileDate`. It also checks index
chromosome lists and first/last-record-position tabix queries per chromosome.
It records those checks separately from byte identity of compressed files or
indices. No new pipeline batches or scientific output shards are introduced.

For the focused native writer-buffer experiment, add `--pairs 2` to the same
command and use a new output directory. This runs fresh workers in alternating
`default → nchunk=100000L`, then reversed order. Both workers retain the same
payload, and their actual writer AST differs only by that single argument.
Expression instrumentation is disabled; load/export boundary events and complete
worker CPU/RSS are retained. The original QS remains read-only. Pinned preflight
fixtures cover mixed floating-point magnitudes, missing values, flags, factors,
anchored indels, duplicate loci, empty tables and indexed queries. Each pair
also rebuilds every index against its own BGZF encoding with pinned Rsamtools and
requires identical decompressed index bytes. Compressed offsets are never
compared between different BGZF encodings. Comparisons run outside worker timing.

The production change is one R argument to the existing VariantAnnotation
writer, with existing compiled package/Rsamtools operations. Python drives the
benchmark and validation; no new non-R production helper is introduced.

Index rebuilds stage physical scratch BGZF copies: Rsamtools normalizes paths,
so symlinks could resolve to original files. Before real-data work, compute
preflight rebuilds four actual fixture indices and verifies the original
BGZF/index SHA256, device, inode, mode, size, mtime and ctime stay unchanged.
Read access time is excluded because reads may update it. The same invariant is
covered by `python3 tests/test_vcf_export_harness.py`.

The paired harness explicitly removes only an absent or known `nchunk=100000L`
argument for its default arm and sets `100000L` for the native arm. It rejects
unexpected configured values. Thus it remains a legacy-versus-native benchmark
after the production writer gains that argument. The manifest and each effective
writer AST record this normalization; fixtures cover configured production calls
and argument removal/restoration. Diagnostic mode retains the actual source
writer setting and records its effective value.
Add `--preflight-only --pairs 2` to run just the configured-source and index
fixtures on compute, without loading the real QS or running timed workers.

## Coverage annotation validation and replay

`calculateBurdens.R` now calls the R coverage functions in `sharedFunctions.R`.
The default `--coverage-annotation r` uses bounded reference windows and native
`data.table::fwrite` buffers piped to bgzip. `--coverage-annotation legacy` selects
the original annotation chain; compressed FASTA and references with ambiguous
contig names also use that fallback. The original R BED writer and
`chunk_runs=1e7` remain unchanged, as do the final QS and published BED schemas.

Production annotation is custom R using existing compiled package internals,
with shell invocations of existing `bgzip`/`tabix` and legacy CLI tools. Python
provides fixtures and replay control; the historical C++ arm is separate.

Run these small fixtures on an allocated compute node in the pinned container,
with `/hidef/bin` on PATH. Use fresh fixture output directories:

```sh
python3 tests/test_annotate_coverage.py \
  --output-dir test-results/coverage-helper
Rscript --vanilla tests/test_coverage_annotation_dispatch.R \
  bin/calculateBurdens.R bin/sharedFunctions.R \
  test-results/coverage-helper test-results/coverage-dispatch
```

The Python fixture produces independent expected outputs using the original
reference preprocessing and complete annotation chain. The R fixture extracts
the actual production functions, CLI options and dispatch expression. It checks
reference boundaries, context/count formatting, fractional and large depths,
small windows, wrapped FASTA, fallback, empty/counts-only output and propagation
of input, compression and indexing failures.

For historical C++ comparisons, obtain the removed helper from the revision that
was actually benchmarked. Compilation is only for this reference arm; the
production workflow no longer compiles it:

```sh
mkdir -p build
git show 604ce8c:bin/annotateCoverage.cpp > build/annotateCoverage.cpp
g++ -O3 -std=c++17 -Wall -Wextra build/annotateCoverage.cpp \
  -o build/annotateCoverage -lhts -Wl,-rpath,/usr/local/lib

Rscript --vanilla scripts/benchmark/benchmark_coverage_annotation.R prepare \
  /path/to/original.outputResults.qs2 /path/to/new-packet.qs2 1-22X strict all

python3 scripts/benchmark/benchmark_r_coverage_annotation.py \
  --packet /path/to/new-packet.qs2 \
  --baseline-script /path/to/baseline/bin/calculateBurdens.R \
  --candidate-script bin/calculateBurdens.R \
  --shared-functions bin/sharedFunctions.R \
  --helper build/annotateCoverage \
  --consumer-test tests/test_coverage_annotation_consumer.R \
  --chromosomes all --output /path/to/new-coverage-replay \
  --repeats 1 --allocated-cpus 4
```

The baseline script above must be the original `6b6598b` version; its unchanged
BED writer is checked against the candidate. The replay freezes and measures
complete C++ and R workers, compares decompressed BED bytes, context-count maps
and indexed queries, and rebuilds each output's own Tabix index. The consumer
test checks ordered downstream count tables and their QS2 round trip; it does
not replace a complete burden-task scientific comparison. Preparation and
validation costs are reported separately from worker CPU. Historical original
outputs may additionally be supplied with `--historical-directory`. Preserve
all historical reports and distinguish them from the new R-port measurements.

## R artifact-cache validation and replay

The production protocol is custom R, integrated through Nextflow/shell and
existing `flock`, `stat`, `sync` and `cp`. Python supplies tests, replay control
and the frozen historical comparison arm, not the current cache implementation.

Run the protocol tests inside the HiDEF container on an allocated node. The
optional frozen Python helper enables mixed Python/R locking and manifest
interoperability checks:

```sh
mkdir -p build
git show 604ce8c:bin/artifactCache.py > build/artifactCache.py
ARTIFACT_CACHE_PYTHON="$PWD/build/artifactCache.py" \
  python3 -m unittest discover -s tests -p test_artifact_cache.py -v

python3 scripts/benchmark/benchmark_r_artifact_cache.py \
  --python-helper build/artifactCache.py \
  --r-helper bin/artifactCache.R \
  --product reference-summary=/absolute/path/to/existing/summary.qs2 \
  --product library=/absolute/path/to/existing/library-directory \
  --out /absolute/path/to/new-cache-replay --repeats 3
```

The replay alternates language order and runs each operation in a fresh worker.
Cold measurements include the same `cp --reflink=never` closed-product builder,
verification, durable publication and startup. Warm measurements verify and
restore the same products with hard links; they fail if the builder executes.
The driver checks product bytes, manifests, hard links and unchanged source
metadata. It removes only the temporary product/cache copies it created.
Outputs must use a new directory. This measures cache overhead, not scientific
preparation; filesystem page-cache state is uncontrolled. CPU includes waited-for
children, and RSS is the maximum individual process/child value rather than a
sum of concurrent processes. Actual workflow behavior is separately exercised
by `tests/test_prepared_cache_workflow.py`.

## Follow-on regression tests

The ordering regressions run on the host with Nextflow available, inside a
one-CPU, 4 GiB, ten-minute compute allocation. Use a fresh directory for each:

```sh
module load nextflow/26.04.0
python3 tests/test_output_assembly_order.py --directory /path/to/new-output-fixture
python3 tests/test_merge_input_order.py --directory /path/to/new-merge-fixture
```

The output fixture extracts the actual channel fragment from `main.nf` and checks all 24 arrival
permutations across interleaved samples, using deliberately nonlexical configured
groups and reversed filenames. It checks exact file order and tuple contents;
it launches no scientific workers. The merge fixture checks 50 sample groups
and 14 demultiplex groups, including multiple barcodes per sample, round2
identifiers, paired BAM/PBI order and singleton lists. It exercises merge-list
grouping and sorting, not upstream demultiplex routing. Production changes are
Nextflow/Groovy; the regression drivers are Python. These correctness tests add
no CPU-saving claim.

The separate disabled layout-policy helper tests are recorded in the
[output-ordering ledger](../../docs/optimization-validation.md#historical-pythonr-run-ordering-failure-and-subsequent-repair).
They exercise R QS-column and Python TSV-row validation helpers, not this
Nextflow arrival-order regression, and do not activate an acceptance policy.

The subsequent [ordered output-stage replay](../../docs/optimization-validation.md#ordered-output-stage-repeatability)
passed two fresh workers per sample and comparisons with the actual pipeline
outputs. It checks strict QS/TSV data and layout without migration normalization.
This is separate from both the synthetic fixtures and original-versus-new
publication acceptance; it does not represent two complete pipeline reruns.

The integrated source passed job `19374216`: all four changed R expression trees were identical to the measured candidate, the 16 parser/168 aggregation/32 VCF cases passed, and all 26 pins and the container identity remained unchanged. These scientific regressions add no performance measurement. See the [compact report](../../docs/benchmarks/torch-2026-10-07-followon.json) for measured scopes and provenance.

From the repository root inside the pipeline environment, run:

```sh
Rscript --vanilla tests/test_generated_pair_parser.R .
Rscript --vanilla tests/test_region_read_aggregation.R .
python3 tests/test_vcf_region_lookup.py .
```

Use a compute allocation with one CPU, 4 GiB RAM and five minutes. The Python VCF driver invokes its paired R companion, so these four source files have three entry points.

| Source file | Invocation and coverage |
|---|---|
| `test_generated_pair_parser.R` | R entry point taking the candidate root; 16 cases covering values, types, missing tokens, expected rejection, grouping, row names and data-frame attributes |
| `test_region_read_aggregation.R` | R entry point taking the candidate root; 168 exact cases executing the actual grouping, threshold and ordered-join expression with controlled fractions |
| `test_vcf_region_lookup.py` | Python entry point taking the candidate root; creates small FASTA/VCFs, runs the real indexing tools, and invokes the paired R test |
| `test_vcf_region_lookup.R` | Companion invoked by the Python driver with candidate root and fixture directory; 32 direct/coalesced target comparisons with GT/no-GT, duplicate records, long deletions, multiallelic/MNV records, ordering and boundaries |

Run the following inside an existing compute allocation with at least one CPU, 4 GiB RAM and five minutes. The example uses the same container and three entry points as the completed suite. The VCF R companion requires the generated fixture directory, so the Python driver is its normal entry point.

```sh
apptainer exec --cleanenv -B /projects \
  /projects/rps/evrong01/evronylab/bin/HiDEF-seq/hidef-seq_3.0.sif \
  /bin/bash -s <<'SH'
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
FOLLOWON_TESTS=/projects/work/evrong01/HiDEF-seq/codex/runtime/followon-tests-v2
FOLLOWON_SOURCE=/projects/work/evrong01/HiDEF-seq/codex/runtime/followon-candidate-v3
Rscript --vanilla "$FOLLOWON_TESTS/test_generated_pair_parser.R" "$FOLLOWON_SOURCE"
Rscript --vanilla "$FOLLOWON_TESTS/test_region_read_aggregation.R" "$FOLLOWON_SOURCE"
/usr/bin/python3 "$FOLLOWON_TESTS/test_vcf_region_lookup.py" "$FOLLOWON_SOURCE"
SH
```

The container supplies R/Bioconductor, tidyverse, dtplyr, data.table, GenomicRanges, vcfR and jsonlite, plus Python, samtools, bcftools, bgzip and tabix. The aggregation test fixes data.table to one thread. The Python driver uses a temporary directory and removes its generated FASTA/VCFs after completion; each entry point prints its PASS result and exits nonzero on a failed assertion or tool call.

The frozen [runner plan](/projects/work/evrong01/HiDEF-seq/codex/runtime/followon-tests-v2/plan.json) binds all 15 candidate/test source hashes and the SIF stat identity; its [reviewed result](/projects/work/evrong01/HiDEF-seq/codex/runtime/followon-tests-v2/reviewed-result.json) records the successful run. The direct commands above rerun the fixtures. To repeat the hash-bound runner, prepare a new plan with a fresh output directory and its matching plan hash; its completed output directory is intentionally not reusable.

The original candidate-source suite passed job `19355705`. Its source hashes remain distinct from the integrated hashes: only an explanatory comment and changed-line whitespace were adjusted in `outputResults.R` and `sharedFunctions.R`. Integration job `19374216` verified parsed-expression identity and reran the relocated suites; its [reviewed receipt](/projects/work/evrong01/HiDEF-seq/codex/runtime/followon-integration-v1/reviewed-result.json) records `COMPLETED`, exit 0, 23 seconds and unchanged source/container guards.

The 432 aggregation comparisons in the exploratory benchmark cover three alternative mappings; the reusable source test covers the selected implementation with 168 production-expression cases. These counts describe different suites and must not be combined into an independent sample count. Output object lifetime and R-annotation compression retain separate full-output and full-row validation evidence; this small suite does not exercise them.

## Incremental burden coverage tests

The selected accumulator uses custom R with existing compiled Bioconductor
operations; it adds no custom non-R scientific helper. Its regression drivers
are R. The allocated integration runner uses Python and shell for source checks,
measurement and job control only.

From the repository root inside the pipeline container, with one CPU, 4 GiB and
ten minutes allocated, run:

```sh
Rscript --vanilla tests/test_burden_native_coverage.R .
Rscript --vanilla tests/test_burden_coverage_optimization.R . /path/to/installed/BSgenome.package
Rscript --vanilla tests/test_burden_sensitivity_optimization.R .
```

The native suite compares 70 cases against independent preserved functions from
`ea8d06f`, including coverage sharing and guarded sensitivity paths. The other
two suites retain their earlier assertions; their loaders also import the new
helpers. Supply the installed BSgenome package directory to the coverage suite
to exercise its actual chrM/chrY checks. These tests do not measure the complete
pipeline or validate the final combined run.

## Allele-key regression

Observed allele-key construction is custom R using existing compiled
vctrs/Bioconductor operations. It introduces no custom non-R production code,
new dependency, CLI or Nextflow/shell change. The regression is R; its allocated
runner uses Python/shell only for source guards, measurement and job control.

From the repository root in the pipeline container, within a compute allocation
with one CPU, 4 GiB RAM and five minutes, run:

```sh
Rscript --vanilla tests/test_allele_interaction.R . NEW_OUTPUT_DIRECTORY
```

The output directory must not exist. The test loads the actual installed encoder
and overlap function without running the filtering CLI. Its preserved 7e
baseline fixture and 60 cases compare exact factor codes/levels and overlap
results, including row/factor order, unused and NA levels, separator collisions,
attributes, fallback types, both strand/adjacency modes and zero-width ranges.
Warnings and errors retain their classes/messages under `warn=2`; inputs must
remain unchanged. The test also requires exactly the two intended helper-call
substitutions. It writes comparisons, session information and a JSON result.

The original 60-case job `19448420` passed in 15 seconds with a 1.17 GiB Slurm
batch peak. The relocated test passed job `19471735`: all 60 exact rows matched
the original fixture, source/container guards passed, and installed R bytes
matched the measured candidate. The integration allocation completed in 16
seconds using 12.693 CPU seconds and a 0.921 GiB Slurm batch peak.
The four nuclear contexts, fresh two-chunk ABBA experiment and
candidate-only mitochondrial checks are recorded in the
[allele-key report](../../docs/benchmarks/torch-2026-10-09-allele-interaction.json).
This fixture adds neither a performance measurement nor full-workflow acceptance.
