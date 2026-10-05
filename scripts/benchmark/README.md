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
