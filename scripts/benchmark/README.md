# Optimization measurement and scientific validation

Run large QS2 reads, writes, and comparisons **inside an allocated compute job**.
The host need not have R or GNU time; the pipeline container may provide both.
The measurement wrapper falls back to Linux Python child-process accounting when
GNU time is unavailable.
All examples write below the working `codex` directory. These tools do not submit
jobs or change pipeline output formats. Scratch component files used by validation
are disposable and are never published pipeline outputs.

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
An allocation reporting less CPU than one of its steps, or an allocation with
still-live steps, also remains unresolved. Later accounting snapshots can resolve
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

The comparator inventories all relative paths, checks byte identity, and compares
TSV/CSV/text/BED/VCF/YAML (including gzip) line by line without numerical tolerance
or sorting. Compression and line-ending differences are allowed. It excludes only
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
