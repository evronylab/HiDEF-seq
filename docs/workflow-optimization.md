# Computational optimizations

The optimization branch changes repeated work and intermediate allocations while
retaining the scientific output schemas and one final QS2 per sample. TSV, VCF,
PDF and indexed coverage BED outputs are retained. The
[validation ledger](optimization-validation.md) records completed component
benchmarks, complete-workload measurements and acceptance of all 711 required
scientific output paths under the agreed comparison rules for `358e4bf`.
Historical results for `3911047`, `ea8d06f` and `7e37ada` remain separately scoped.
Helper generalization and production cleanup at `358e4bf` have completed fresh
full-workflow execution, accounting and scientific acceptance. The output-only
R handoff installed later in `3dabe21` has separate both-sample and lifecycle
proofs; the full-workflow measurements remain those of `358e4bf`.

## Production implementation and validation status

Production now uses one R coverage annotator, one shared reference-summary
builder and one shared germline-coverage-filter builder. Standalone consumers
use those same builders when no prepared product is supplied. Sparse coverage,
allele grouping, generated-field parsing and compressed-tag lookup now support
the general input cases without dispatching to an older scientific
implementation. No selectable legacy scientific algorithm or experimental
storage path remains in production. Historical implementations remain only as
test or benchmark oracles. These changes add no custom non-R scientific code.

The frozen `358e4bf` source passed nine focused kernel suites, annotation and
preparation checks, cache protocol checks and a Nextflow configuration preview.
Its fresh two-sample workflow completed all **649 tasks**, with zero failures,
retries or cached tasks, and the controller completed successfully. Accounting
and scientific acceptance are independently reviewed. All 711 required
publications, 20 prepared products and two cache startup checks passed.
Component correctness does not establish whole-workflow performance: an
eight-input burden-kernel replay used 16.11% more operation CPU and 7.75% more
observed process high-water memory after generalization. Extraction and allele
kernel CPU fell 19.76% and 22.11% in their respective replays. Those measurements
exclude input loading, explicit GC and serialization from operation CPU.
High-water memory was sampled through the operations, before any later
serialization peak; these are not complete-task measurements. The completed
workflow accounting and task-memory observations below now supply that broader
measurement scope; the [finalization report](benchmarks/torch-python-r-finalization.json)
records the separate completed scientific acceptance.

The 631 common downstream tasks used **151.542205 actual CPU-hours**, versus
447.194170 for the original 749-task workload (**66.112661% less**). The 18
cold-preparation tasks used a separate **7.034264 CPU-hours**. All 649 workers
plus the controller used **158.609414 CPU-hours**, with preparation counted once;
validation and exploratory jobs are excluded. CPU here is each normal
allocation's reported UserCPU + SystemCPU, summed once at retained display
precision. Reserved cores times elapsed time totalled **218.801667 core-hours**
for workers plus controller; that is a different measure from actual CPU.
The original has no matched cold-preparation measurement or end-to-end wall-time
comparison. The [current accounting ledger](optimization-validation.md#current-accepted-workflow-358e4bf)
retains the stage breakdown and workspace evidence.

For the eight burden tasks, Nextflow's maximum sampled process-tree RSS was
**44.944 GiB** and its maximum sampled sum of process high-water marks was
**44.970 GiB**, versus 230.348/235.352 GiB originally. The two final-output
tasks measured **47.249/47.492 GiB**, versus 76.446/76.446 GiB originally.
The [reviewed process-memory supplement](benchmarks/torch-finalization-process-memory.json)
pins all ten task mappings and counters. The separate maximum task/step Slurm
RSS was **136.687 GiB** for burdens and **51.600 GiB** for outputs. These are
different instruments: sampled process-tree sums can miss short-lived children
and count shared pages repeatedly; process high-water marks need not occur
simultaneously. None is an R-heap measure, simultaneous pipeline peak or minimum
safe allocation. No resource request was reduced from these observations.

Nextflow publication again follows the original per-process rules: hard links
for BAMs and optional intermediate QS files, moves for coverage and final sample
outputs, and copies for logs, barcode FASTA and VerifyBAMID results. No global
publication-mode override remains, and default input staging is unchanged.
Prepared-cache restoration is a separate operation: it hard-links a verified
cache product into the consuming task on the same filesystem and copies it only
when the filesystem rejects that link with `EXDEV`. This does not change a
process's Nextflow publication policy.

The unused whole-genome reference trinucleotide BED producer is removed. Its
two private cache products disappear; all published coverage BEDs and their
indexes remain required. The completed scientific gates accept all 711
publications and the remaining 20 prepared products, with strict QS/TSV
ordering and the existing BAM coordinate-tie exception. No historical QS/TSV
layout allowance applies to the new comparison. Earlier accepted measurements
below describe their frozen revisions. The separately validated single-public-QS
output helper was installed in `3dabe21` and is excluded from the `358e4bf`
full-workflow measurements.

The reopened Arrow/DuckDB evaluation is complete. Exact storage and export
proofs passed, but the current SQL coverage integration is rejected for
production: complete five- and fifteen-input workers used 44.5% and 176.2% more
CPU without improving process peak memory. No disk-backed call representation
has been adopted. This experimental scientific logic uses SQL through existing
DuckDB/DBI packages; it adds no custom C++/Rcpp.

## Historical accepted full-workflow comparison: `7e37ada`

The `7e37ada` workflow completed all 650 tasks and its accounting. The comparable
downstream workload used 631 candidate tasks versus 749 original tasks, consuming
157.81978 CPU hours versus 447.19417 (64.71% less). The 19 cold-preparation tasks
used another 8.99136 CPU hours; the original run has no matching cold-preparation
measurement. These are complete-workload accounting totals, separate from the
component timings below.

All 22 prepared-cache comparisons and two startup checks passed. Both samples
also passed the separate deterministic output-stage replays described below.
The composed original-versus-new scientific acceptance gate passed all **711
required paths**: 687 direct strict passes, 18 approved layout proofs, and two
BAMs with four owner-linked indices. Original failed reports remain intact.
The [complete-workload report](benchmarks/torch-python-r-followon.json) binds
execution, accounting and acceptance; the [validation ledger](optimization-validation.md)
also records the separate strict output-stage repeatability proof.

This run includes the selected native burden accumulator and the five follow-on
R changes. It excludes the later allele-key helper installed in `00f6e33` and the
experimental temporary coverage storage described below. Production includes
custom Python BAM dispatch and Nextflow/Groovy orchestration alongside R and
existing compiled packages and CLI tools; benchmark-only Python/shell controls
are identified separately.

## Implementation languages

The labels below distinguish implemented, rejected and experimental approaches,
including orchestration.
An R implementation can use existing packages with compiled internals without
adding a custom C++ or Rcpp helper. Existing external tools are listed separately.
Benchmark and regression infrastructure also uses Python and shell; the VCF
regression test, for example, has a Python fixture driver and an R comparison.

| Optimization | Custom implementation | Non-R involvement |
|---|---|---|
| Prepared reference summaries and germline coverage caches | R preparers; Nextflow and shell integration | Nextflow/shell orchestration and existing reference/coverage tools; no new compiled helper |
| Artifact-cache protocol | R, including shared functions | Existing `flock`, `stat`, `sync` and `cp`; Nextflow/shell invokes the R entry point |
| Output publication, original per-process policy | Nextflow | Existing Nextflow hard-link, move and copy modes; no custom filesystem selector |
| Deterministic final-output assembly (output-stage reproducibility passed) | Nextflow/Groovy | **Yes: custom non-R orchestration**; orders the small input-file list without sorting scientific records |
| BAM dispatch | Python using `pysam`; Nextflow/shell integration | **Yes: custom Python**; existing BAM indexing tools and pysam's compiled library |
| Extraction quality lookup and generated-field parsing | R | Existing R packages; no new custom non-R implementation |
| Chromosome restriction, germline summaries, threshold and coordinate helpers, region-read aggregation | R | Existing R/Bioconductor packages, including dtplyr/data.table; no new custom non-R implementation |
| Observed allele-key construction, selected after context and two-chunk comparisons | R | Existing compiled vctrs/Bioconductor operations; no custom non-R code or new CLI |
| Coalesced VCF lookup | R constructs the query commands | Changes shell command construction using existing Bash/`bcftools`; no new standalone executable |
| Coverage accumulation and sensitivity calculations | R | Existing R/Bioconductor packages; no new custom non-R implementation |
| Coverage annotation and lower compression level | R | Existing `bgzip`/`tabix` for compression and indexing; no custom non-R annotation kernel |
| Final table formatting, VCF writer buffer and earlier object release | R | Existing R/Bioconductor packages; no new custom non-R implementation |
| Coordinate-sort reuse, rejected | Nextflow/shell command changes | Existing `pbmerge`/`samtools`; would involve non-R orchestration |
| Shared-session filtering, rejected | R worker plus workflow changes | Would involve Nextflow/shell orchestration |
| Flattened interval positions and direct logical masks, rejected | R | Existing Bioconductor packages; no new custom non-R implementation |
| Endpoint accumulation, tested and not selected | R | Existing Bioconductor packages; no new custom Python/C++/Rcpp implementation |
| Endpoint capacity fallback, tested and not selected | R | Existing Bioconductor coverage and Rle addition; no new custom non-R implementation |
| Native Rle addition and immediate sharing, selected after whole-worker proof | R | Uses existing Bioconductor compiled operations; no new custom Python/C++/Rcpp implementation |
| Temporary-QS coverage storage, experimental and not adopted | R | Existing compiled qs2/Bioconductor packages and annotation/index CLIs; Python/shell benchmark controls only, no custom non-R scientific kernel |
| DuckDB calls and coverage storage, current SQL integration rejected | R with embedded SQL | **Yes: SQL query/reduction logic** through the existing DuckDB package; no custom C++/Rcpp; Python/shell benchmark controls only |
| Single-public-QS output task, both samples validated; installed in `3dabe21` | R | Existing packages and export CLIs; no new package or custom non-R scientific logic; Python/shell measurement and validation controls only |

## Torch launch and resource settings

For each test, launch from its own run directory with its own work directory.
Keep both directories and the run's `.nextflow` state for resume. This example
uses the current branch each time; set the two run paths before submitting:

```bash
RUN_DIR=/path/to/test-run
RUN_PARAMETERS=/path/to/test-run/parameters.yaml
NEXTFLOW_TORCH_CONFIG=/projects/rps/evrong01/evronylab/bin/HiDEF-seq/nextflow.torch.config
HIDEFSEQ_GITTAG=optimization

cd "$RUN_DIR"
sbatch --time=48:00:00 --mem=4G --cpus-per-task=1 \
  --account=torch_pr_285_general \
  --wrap "module purge && module load nextflow/26.04.0 && nextflow -config \"$NEXTFLOW_TORCH_CONFIG\" run evronylab/HiDEF-seq -r \"$HIDEFSEQ_GITTAG\" -latest -params-file \"$RUN_PARAMETERS\" -resume -work-dir \"$RUN_DIR/work\" -with-report -with-trace"
```

The outer allocation runs the Nextflow controller; each process receives its own
SLURM allocation. Do not start a second controller for the same run while the
first controller or its workers are still active. An explicit
`-resume <session-UUID>` selects the intended session when several runs share a
launch directory. Use distinct report/trace filenames for repeated tests.

To reproduce the latest fully benchmarked and accepted workflow exactly, set
`HIDEFSEQ_GITTAG=358e4bfcc89e409456c6ec97d8fdd142475c07f6` instead of the branch
name. That accepted run includes the allele-key helper, generalized R
implementations, preparation consolidation and publication cleanup. The branch
additionally contains the output-only R handoff installed in `3dabe21`, with
separate both-sample and 42-check lifecycle validation; no full workflow was
measured at that later revision. The earlier accepted candidate was `3911047`; subsequent
reference-summary consolidation into `sharedFunctions.R` preserved its function
definitions and passed focused container checks in job `19289528`.
Final validation passed under the approved BAM rule: logical records, multiplicities,
coordinates and coordinate ordering remain exact; only relative ordering at
identical coordinates may differ. See the validation ledger for the complete
comparison scope and preserved original strict-order report.

The two-sample, 60-chunk-per-sample `7e37ada` Torch comparison used the rebuilt
`hidef-seq_3.0.sif`; the original image was retained for baseline validation.
These allocations matched between baseline and candidate:

| YAML process suffix | Memory | Time | CPUs |
| --- | ---: | ---: | ---: |
| `extractCallsChunk` | `32GB` | `1h` | 1 |
| `filterCallsChunkChromgroupFiltergroup` | `48GB` | `2h` | 1 |
| `calculateBurdensChromgroupFiltergroup` | `288GB` | `30h` | 2 |
| `outputResultsSample` | `96GB` | `4h` | 1 |

Set the corresponding `mem_<suffix>` and `time_<suffix>` keys in the run YAML;
CPU counts above come from the workflow process definitions. The template's
64 GB burden/output requests are not validated for this workload. The accepted
`7e37ada` comparison's maximum task Slurm memory was 135.261 GiB for burdens and 51.733 GiB
for final output. Requests were kept unchanged for the comparison. Size future
requests from representative full tasks with headroom;
these peaks are workload-dependent and do not establish a universal lower limit.
Slurm task RSS, Nextflow's sampled process RSS, and simultaneous pipeline memory
are different measurements.

Memory is not bounded independently of sequencing yield. With a fixed
`analysis_chunks`, larger inputs create larger extraction and filtering tasks.
The Python dispatcher streams BAM records but retains a ZMW-to-chunk lookup
whose size grows with the ID enumeration. Burden calculation loads each input
QS2 in full and retains calls, germline records and aggregate coverage; final
output loads each burden QS2 and constructs the combined sample object. More
chunks would be a separate configuration change requiring its own comparison,
and would not remove those sample-wide allocations. Windowed R annotation
limits its additional working buffers, not the memory of all retained data.
The 5/15/60-input storage experiments test accumulation across existing inputs;
they do not establish scaling for larger individual input QS2 files.

## Prepared caches and publication

Prepared reference libraries and summaries, per-individual VCF annotations,
germline BAM coverage/calls, and thresholded region tracks use
separate directories under `cache_dir/prepared/v1/<kind>/<identity>/`.
Each identity includes the relevant settings, process code, helper code, and
container/tool identity. Reference preparation does not depend on individual
VCF settings; a VCF change for one individual does not rebuild other prepared
artifacts. Existing downstream Nextflow configuration signatures remain
conservative and may rerun unrelated downstream analysis tasks.

Low-coverage germline filters additionally have their own scope for each
individual and distinct `min_germlineBAM_TotalReads` threshold. Preparation uses
the same whole-genome wiggletools comparison, fill-in, BigWig conversion and
import as the former per-chunk path, then stores the intervals as QS2. The
consumer performs its original reference Seqinfo assignment after loading.
Thresholds shared by multiple filtergroups reuse one preparation; genome-wide
filtered-base statistics are preserved. Standalone YAML without this mapping
uses the same shared builder as the preparation task.
The existing empty-stream failure is preserved: when a threshold comparison
produces no intervals, `wigToBigWig` rejects the empty input in both paths.

A reference-summary bundle is prepared after BSgenome installation and reused
across samples, chunks and chromosome groups. It contains full-reference N
intervals with their Seqinfo and integer trinucleotide counts for each chromosome.
Filtering retains the original whole-genome N-base statistic; burden tables sum
the requested chromosomes and keep the existing channel/fraction transformations.
Circular chromosomes retain the original linear counting boundaries. Standalone
YAML without `reference_summary_file` uses the same shared builders, requesting
only the components the consumer needs.
`tests/test_reference_summary.R` checks empty/N-free/all-N references, chromosome
order, repeated selections, types, and equivalence of the earlier chromosome
restriction in call loading.

These reference-summary functions live in the single common
`bin/sharedFunctions.R` file. The separate helper file has been removed, and its
production callers, benchmark and cache dependency list use the common file.
Every prepared-artifact identity includes the complete `sharedFunctions.R` hash.
Any change to that file, including the consolidation or the later allele-key
helper addition, therefore changes all prepared-cache identities on the next
launch, even when a preparer's own calculation is unchanged. Downstream R task
signatures also include this hash. Existing cache bundles are retained; they are
not silently reused under changed source hashes.

Large external files use an explicitly labeled identity comprising canonical
path, size, modification time, and Unix change time and inode where supported.
Node-local device numbers are excluded: the same file on Torch's shared filesystem
has different device numbers on different nodes. These are **metadata identities,
not input content hashes**. No full BAM,
reference, or container checksum is computed on the orchestration node. Small
scripts use SHA256. Existing external executable paths have metadata identities;
commands supplied by the container use its identity plus configured command.
Use immutable input files and pinned container images. Metadata-preserving
content changes or mutable container tags can evade this identity scheme.

`artifactCache.R` runs inside preparation tasks and uses the reusable cache
functions in `sharedFunctions.R`, with `jsonlite`, `digest` and `openssl`.
Linux `flock` holds the same lock protocol as the former Python helper; GNU
`stat` checks file types and `sync` provides file/directory fsync. No R process
retains complete prepared products in memory during copying or hashing.
It serializes simultaneous builds of the same identity, accepts only explicitly
declared closed products after a successful command, checksums those products,
and copies them into a temporary sibling directory. It writes
`manifest.complete.json` there and atomically renames the complete directory
into place. Failed builds/copies cannot create a complete entry. Existing corrupt
entries fail verification and are never silently overwritten. Inspect and move
a corrupt entry aside before rebuilding it.

The R port passed 16 protocol tests and an actual 50-task Nextflow fixture
(jobs `19299712` and `19299729`). Cold, warm, missing-entry and concurrent
workflows produced the expected build counts and identical scientific objects.
Three paired measurements on each of four existing product bundles measured
about 0.4–1.4 additional CPU seconds per cache operation versus Python. For the
1.63 GB bundle, median cold CPU was 2.96 seconds in Python and 4.35 seconds in R;
warm CPU was 1.12 and 1.92 seconds. R used roughly 73–129 MiB maximum individual
process/child RSS versus 16–19 MiB for Python. This port consolidates the
implementation in R; it does not itself reduce resource usage.

Those measurements include startup, checksumming, publication or restoration
and waited-for child CPU. Cold measurements use the same closed-product copy
builder in both arms, so they do not measure the unchanged scientific
preparation algorithms. Filesystem page-cache state was uncontrolled. Evidence
and source hashes are retained in `runs/r-cache-v1/result.json`. Accounting for
the completed `ea8d06f` run is complete; that run remains scientifically
unaccepted because LIB1 failed the strict comparison. The later combined
`7e37ada` run has completed accounting, cache checks and composed scientific
acceptance, as summarized [above](#historical-accepted-full-workflow-comparison-7e37ada).

Preparation tasks verify bundle contents on every workflow launch, including
resume. Cache hits restore products to the Nextflow work directory with hard
links, falling back to copies across filesystems. Published cache directories
are never modified by the helper; consumers must treat restored links as
read-only. Legacy cache files are left intact and are not trusted automatically.
The first run therefore prepares the new cache namespace.

R scripts receive a full effective YAML named by the SHA256 of its exact text.
It retains original input configuration keys with intentional runtime overrides
and adds only `reference_cache_dir`, `reference_summary_file`, `cache_artifacts`
and `germline_coverage_filters`. Run snapshots and effective YAML contain ordinary
scalars/collections, excluding injected runtime objects and unrelated defaults.
The source parameters file is read only. The shared resolver maps the original
product basenames to those paths. Standalone R invocation without these added
fields resolves inputs directly in its configured `cache_dir`. Conflicting prepared product basenames
fail early instead of sharing an ambiguous cached file. Scientific product
schemas and basenames are unchanged; absolute cache paths are metadata.
Because the immutable effective configuration path is included in R task scripts,
changes to any retained configuration field can invalidate downstream R tasks.
The existing configuration signatures therefore do not yet provide fully scoped
downstream resume. This conservative behavior preserves concurrent-run safety;
separating execution parameters from full provenance is a future design change.

Output publication follows the original pipeline's per-process policy:

- Merged BAMs and optional split BAM/extraction/filter/burden intermediates use
  `link` (hard links).
- Coverage BED.gz files and their indexes, and final per-sample results, use
  `move`.
- Logs, barcode FASTA files, and VerifyBamID outputs use `copy`.

The automatic filesystem-based publication override has been removed. A hard
link gives the work and results paths the same stored file; editing its contents
through either path affects both names. Removing one name leaves the other
usable. Moved outputs follow the original resume behavior and no longer remain
at their work paths. These publication rules are distinct from Nextflow input
staging and from preparation-cache restoration. No input-staging link mode has
been overridden. The existing preparation barrier remains, with the unused
whole-genome trinucleotide BED producer removed.

## BAM dispatch

`countAnalysisZMWs` retains one exact `zmwfilter --show-all` enumeration alongside
its existing count text. Empty samples are skipped as before; effective chunk
count remains the smaller of requested chunks and enumerated IDs. Nextflow
stages `splitBamByZmw.py` and calls it with `python3` from the container's PATH.
No Python path entry is needed in the run YAML. The interpreter requires pysam with libdeflate
support; initial benchmarks used pysam 0.23.3 without an exact-version check in
the splitter. The script SHA256 is an explicit task input, so changing
source invalidates dispatch even if its size and timestamp are preserved.
There is no compiler process or per-task dependency installation.

Each sample dispatch assigns the same contiguous quotient/remainder ID
partitions and streams original BAM records to those chunks. Numeric ZMW IDs
remain the grouping key across runs. If an ID occurs in multiple legacy
partitions, its records are emitted to each matching chunk, preserving the
legacy include behavior. Missing `zm` tags follow the pinned PacBio index's
observed zero-ID behavior. Headers, record order and tags are preserved.

Production dispatch inputs always pass through `mergeAlignedSampleBAMs`, which
unconditionally runs pinned `pbmerge` (including one input), then coordinate sort
and indexing. This boundary matters for noncanonical synthetic PacBio headers:
legacy `zmwfilter` adds missing `@HD pb:5.0.0` and `@RG PM:SEQUEL` defaults, while
the dispatcher preserves its input header. In the long-CIGAR boundary fixture,
one-input `pbmerge` added those defaults before dispatch; all non-`@PG` header
lines, the binary reference dictionary and all three ordered raw BAM records
then matched legacy exactly. Two malformed inputs lacking the PacBio version
were rejected by `pbmerge` before dispatch. The direct unmerged synthetic case
therefore remains a documented standalone-helper difference, with no observed
production mismatch. Program-header differences remain explicit review findings;
this is not a blanket header whitelist. Boundary evidence is in
`runs/dispatcher-header-boundary/results.json` (job 19201853); that focused probe
did not rebuild indices, which are covered by the separate dispatch fixtures.

Up to 128 output writers stream binary BAM records through pysam. Input
decompression uses the configured thread count; output compression is serial
so the writer count does not multiply the thread budget. Larger chunk counts
use additional sequential input passes without changing chunk assignments.
The ZMW-to-chunk mapping still scales with the number of enumerated IDs.
Each output receives
its PacBio and samtools indices sequentially. BAM/PBI/BAI filenames and the
seven-field downstream tuples are unchanged; tuple expansion explicitly matches
basenames, including the single-chunk case and chunk numbers above nine. One
sample-level dispatch log replaces separate chunk logs. `pbmerge`, coordinate
sorts and the global preparation barrier remain unchanged.

The adversarial helper harness covers duplicate IDs across groups, missing tags,
QNAME disagreement, one/many writers, original headers/tags/order, indexing and
truncated input. `tests/test_dispatch_workflow.py` extracts the actual process
and channel code for a small multi/single/empty-sample Nextflow replay; it checks
all resulting tuples and compares each chunk with legacy `zmwfilter --include`.
The Python Nextflow replay passed execution, publication and exact binary BAM
record comparisons for 12-chunk and one-chunk samples and skipped the empty
sample (job `19299386`). All five tasks were cached on unchanged resume. A
source-only edit preserving file size and timestamp reran exactly the two
dispatch tasks while retaining all three cached enumeration tasks. Evidence is
in `runs/python-splitter-v1/workflow/`. The earlier C++ workflow evidence remains
unchanged in `runs/dispatch-workflow-v3/`.

The historical C++ LIB1 replay produced all 60 chunks in 2,161.59 actual CPU seconds
(1,092.32 wall seconds; 526,472 KiB maximum process RSS). Quickcheck and PBI/BAI
indexing passed for all chunks and used another 1,499.20 actual CPU seconds.
Chunk 1's complete ordered SAM record stream matched the legacy chunk SHA256
`b09dc1de7a0fdc777326febf672a3dbfe1163190a119e1dc9f7e9d15948b8490`.
These measurements and hashes are retained in `runs/dispatch-benchmark/`. This
is a historical component measurement.

The initial Python full LIB1 paired replay passed all 60 chunks and 10,382,760 ordered raw
BAM records, reference dictionaries and headers under the existing program
command-line metadata exception (job `19299191`). Both implementations completed
all BAI/PBI indexing. The paired costs were:

| Operation | C++ CPU seconds | Python CPU seconds | C++ elapsed seconds | Python elapsed seconds |
| --- | ---: | ---: | ---: | ---: |
| Split | 1,856.50 | 1,891.36 | 938.50 | 1,630.67 |
| Index | 1,149.14 | 1,149.01 | 588.31 | 589.00 |
| Combined | 3,005.64 | 3,040.36 | 1,526.81 | 2,219.67 |

Python used 1.16% more combined CPU and 45.38% more elapsed time. Maximum stage
process/descendant RSS increased from 515 to 875 MiB. Its serial output
compression preserves the thread budget but increases elapsed time. This is
one full-sample pair with sequential arms, excluding compilation and equivalence
validation. It used isolated Conda Python 3.12.15/pysam 0.23.3 with libdeflate
inside the original SIF; it does not validate a rebuilt SIF or the complete
follow-up pipeline. Evidence is in
`runs/python-splitter-v1/full-conda/paired-comparison-reviewed.json`.

The final rebuilt-container pair (job `19322722`) uses system Python 3.12.3,
pysam 0.24.1 and HTSlib 1.24 with libdeflate, invoked through `python3` on PATH.
It compares the frozen C++ splitter under the old image with the production
Python splitter under the rebuilt image on the same node and input. All 60
chunks and 10,382,760 ordered raw records, complete headers and reference
dictionaries matched exactly. All 240 saved BAI/PBI comparisons against fresh
own-BAM rebuilds passed.

| Operation | C++ CPU seconds | Python CPU seconds | C++ elapsed seconds | Python elapsed seconds |
| --- | ---: | ---: | ---: | ---: |
| Split | 2,254.70 | 2,230.42 | 1,140.20 | 1,975.50 |
| Index | 1,405.79 | 1,458.13 | 734.84 | 767.62 |
| Combined | 3,660.56 | 3,688.62 | 1,875.25 | 2,743.29 |

Combined actual CPU increased 0.77%; elapsed time increased 46.29%. Individual
process peak RSS increased from 516 to 874 MiB. This port removes the custom
C++ implementation with nearly unchanged CPU cost; it does not save memory or
elapsed time relative to C++. This is one sequential pair, C++ then Python,
with no randomized-order or cold-cache claim. Validation costs are excluded.
The [component report](benchmarks/torch-2026-10-07-splitter.json) records metrics
and evidence checksums. Full-pipeline acceptance and performance are separate.

Against the original per-chunk `zmwfilter` implementation, completed two-library
accounting for `ea8d06f` measures the entire enumeration, splitting and indexing
stage at **66.3835 CPU-hours originally versus 2.20819 CPU-hours with Python**:
**64.1753 CPU-hours saved (96.67%)**. Python retains the one-pass dispatch saving.
These are allocation UserCPU + SystemCPU totals, excluding validation. This is
not a 30-fold elapsed-time claim: the original chunk jobs ran in parallel.
The source is workspace
`runs/candidate-python-r-20261006/accounting/comparison.json`,
`common_downstream_stages.bam_dispatch`. The run's unrelated output-assembly
failure remains recorded below; these stage measurements do not accept the
whole workflow. The user explicitly accepted the Python implementation.

## Processing inside the R scripts

### Extraction

`extractCalls.R` decodes the run-length encoded `sa` tag directly into an Rle and
reverses its runs when strand orientation requires it. It no longer expands all
read-level `sa` tags into dense arrays for quality lookup. Queried positions are
mapped to run values using cumulative endpoints. Compact index normalization
preserves ordinary R subscript semantics without expanding the tag values into
a dense vector. The published extraction object still contains its
original `sa` Rle representation, and `sm` and `sx` retain their original vectors.

For indel quality lookup, read names are matched to tag-list positions once for
each strand. The six existing assignments then use integer positions, avoiding
repeated linear searches through long named lists. All three tag lists originate
from the same BAM rows and share the same name order. The temporary index is
removed before aggregation; first-match, absent-name and output-order behavior
is unchanged. Run `tests/test_extract_sa_rle_optimization.R` and
`tests/test_indel_tag_indices.R` in the pipeline container for focused regression
checks against the actual production assignments.

### Filtering

`filterCalls.R` restricts calls to the requested chromosome group before expensive
processing. Whole-genome statistics and germline variant filters retain their
original scope. In the germline annotation step, only variants matching the call
keys need the expensive per-variant summary; the complete filtered germline table
is still retained where downstream calculations require it. The reusable reference
summary and germline coverage masks above remove repeated reference scans and
whole-genome coverage conversion from individual filter tasks. Filtergroups remain
separate jobs: the tested shared-load alternative did not justify its additional
memory and complexity.

The Python/R follow-up also removes two kinds of repeated R overhead inside
filtering. Strand-threshold checks keep the same `mean`, `all`, `any` and
short-circuit expressions while replacing generic missing-value replacement
with scalar checks. Query-to-reference conversion keeps `mapFromAlignments`
and constructs its GRanges inputs/results directly, using the returned read
indices to attach metadata instead of repeated data-frame conversions and a
join. Adversarial fixtures and all four real coordinate queries passed exact
values, order, types and Seqinfo checks. Two alternating pairs of complete strict
nuclear filtering workers also passed exact comparisons of all seven QS2
components, with zero ignored fields, failures or numerical reviews (job
`19300725`). CPU reductions were 6.95% and 3.36%; mean complete-worker CPU fell
from 1,407.89 to 1,336.09 seconds, a 5.10% reduction. Mean elapsed time fell
6.96%. Maximum RSS differed by only about 0.4%, which does not establish a
meaningful memory improvement. These measurements cover one real chunk, not a
whole-pipeline extrapolation; evidence is in `runs/filter-r-v1/summary.json`.

Early isolated coordinate-helper timings are not accepted performance evidence:
R's lazy argument evaluation included input construction in only the original
arm. The diagnostic harness now forces inputs before either timing. Scientific
comparisons were unaffected, and the separate complete-worker benchmark does
not have that timing problem. The changes are retained on the strength of those complete-worker comparisons and their simpler function bodies. Subsequent profiling selected a targeted dtplyr conversion for region-read aggregation: it preserves `base::mean`, threshold semantics and ordered joins while avoiding repeated group work. Its separate component and combined-chain evidence is recorded in [the follow-on results](#follow-on-r-optimization-results). This does not introduce a broad dtplyr rewrite or a custom Rcpp helper.

The subsequently selected allele-key helper constructs labels only for observed
combinations of the four allele columns, avoiding unobserved Cartesian labels.
It retains factor preparation and ordering, with a common encoder for missing
values, separator-containing labels, factor levels and recycled inputs. It no
longer delegates unusual inputs to `interaction` or allocates all unobserved
Cartesian combinations. This is custom R
using existing compiled vctrs/Bioconductor operations; it adds no custom non-R
production code, dependency, CLI or Nextflow/shell change.

Four nuclear sample/filtergroup ABBA comparisons showed mean complete-worker CPU
reductions of 4.39–15.80%, with context-dependent process-peak changes, including
increases. A fresh extraction of distinct LIB1 chunks 1–2 gave a 7.69% CPU
reduction and a 16.02% lower mean individual-process peak for strict nuclear
filtering. Ordered scientific comparisons passed in all four nuclear contexts,
four candidate-only mitochondrial contexts and the two-chunk comparison.
These scoped results do not establish whole-pipeline speed or a memory bound
for larger inputs. The [allele-key report](benchmarks/torch-2026-10-09-allele-interaction.json)
retains observations, source identities and limits. The earlier helper was
installed after integration job `19471735` passed its relocated 60-case regression
and source/container guards. This check adds no performance measurement.
The full-workflow measurements above remain attributed to `7e37ada`, which
predates this helper. Its completed acceptance does not extend the helper's
scoped performance evidence into a full-workflow measurement.
The separate scalar-bound extension was not integrated.

### Burden and sensitivity calculations

`calculateBurdens.R` accumulates coverage one category at a time and replaces that
category's prior value. Joins operate on category metadata instead of carrying
old and new genome-wide coverage list-columns through the join. The final coverage
tables retain their original structure.

Barcode-orientation tracks are constructed only for requested non-mutation call
types with asymmetric final-round barcodes, where the downstream summaries consume
them. The observed orientation incorporates both demultiplexing rounds; asymmetric
round-2 outputs remain supported. Aggregate duplex coverage and the existing
plus/minus and even-depth consistency checks remain in place.
Six barcode scenarios, including both rounds, passed regression fixtures. The
measured real-data workload had no second-round Lima; a real second-round run
was not benchmarked.

Sensitivity retains coverage only at selected high-confidence germline variant
positions instead of accumulating extra genome-wide coverage tracks. The original
per-VCF quality quantiles are calculated over the whole genome before chromosome
selection. Each indel flank is summed across all chunks before the minimum of its
flanks is taken. These details preserve the denominator, including variants whose
two flanks are covered by different chunks. Indel reference annotation also loads
only chromosomes containing calls.

The coverage BED writer retains its `chunk_runs = 1e7` setting. After it writes
the original coverage runs, `annotate_coverage_row()` in `sharedFunctions.R`
reads indexed FASTA windows and emits the same per-base BED records and context
counts. It uses `data.table::fread` for run tables and `fwrite` buffers piped
directly to bgzip, avoiding per-base R text construction and temporary expanded
BEDs. Reference windows and input buffers bound annotation's additional memory;
they do not partition calls, add chromosome jobs or divide the final QS2.

Coverage annotation has one R implementation and no method selector. The
annotation reader accepts ordinary indexed FASTA, gzip/BGZF, bzip2 and xz files.
Compressed references are expanded once per task into a temporary FASTA in the
task working directory using a 1 MiB buffer, then use the same window reader;
the temporary file is removed on success or failure. Other reference setup
steps retain their own format requirements. Empty contigs are permitted, and
chromosome names are read literally from the FAI. This fixes the old shell
parser's truncation of names containing colons or hyphens. Invalid indexes and
annotation, compression or indexing failures stop the task.

The obsolete whole-genome trinucleotide BED preparation and its cache entry
have been removed. Their unused `seqkit_bin` and `bedtools_bin` entries have
also been removed from configuration templates; existing YAML files may retain
those keys. Existing cache files are left untouched. Fixtures cover
ordinary and compressed references, literal names, empty contigs, exact context
arithmetic, formatting, ordering, reference edges and failure cleanup. Earlier
full-workflow and C++ comparisons below retain their original revision scopes;
they do not constitute validation of subsequent consolidation changes.

The paired chr22 replay passed all eight exact BED/context-count comparisons
and both implementations' own-file index checks (job `19299543`). C++ used
130.34 CPU seconds and 716 MiB maximum individual process/child RSS; R used
189.60 CPU seconds and 859 MiB. The actual downstream ordered tables and QS2
round-trip also passed (job `19299897`). These results establish a measured
45.5% CPU cost versus C++ for this component. The original pre-optimization chr22 operation
used substantially more CPU in earlier measurements, but those runs were not
this same-node paired comparison.

The full LIB1 strict nuclear replay also passed (job `19299718`): all eight
ordered decompressed BEDs and context-count files match both the paired C++
worker and the preserved original outputs. All 16 saved indexes match their
own BED's rebuilt index payload, and the actual downstream ordered context
tables, factors, fractions and QS2 round trips pass. The original writer's
`chunk_runs = 1e7` remains unchanged.

| Nuclear annotation worker | Actual CPU hours | Elapsed hours | Maximum individual process/child RSS |
| --- | ---: | ---: | ---: |
| Historical original | 7.1969 | 4.4383 | 16.730 GiB |
| Paired C++ | 3.4005 | 2.8181 | 16.729 GiB |
| R | 4.5582 | 4.3535 | 19.626 GiB |

R costs 34.05% more CPU than the paired C++ implementation and uses 17.31%
more peak individual-process memory. Compared with the earlier original run,
R uses 36.66% less CPU; that historical comparison was not a randomized paired
measurement. It does not show an individual-process memory saving. The sampled
sum of process RSS fell from 31.768 GiB in the original to 19.642 GiB in R,
but that separate metric can double-count shared pages and miss brief peaks.
These are annotation-worker measurements, including compression and indexing,
not measurements of the whole burden task. Evidence is retained in
`runs/r-coverage-v1/component-summary-after-nuclear.json`. Full-workflow accounting
and scientific acceptance are complete as summarized [above](#historical-accepted-full-workflow-comparison-7e37ada).

### Final output

The approved Nextflow change orders each sample's burden input files by the configured
chromosome group and filter-group order before passing them to R. Previously, task
completion order could reverse filter-group blocks in the combined QS2. This
correctness fix sorts only the small list of filenames and retains every record
within each group. The sample still starts output as soon as all its own burden
tasks finish. This involves custom Nextflow/Groovy code; it is not a claimed
performance optimization. Sample BAM merges also use YAML run and sample-entry
order, including multiple barcode assignments to the same sample within a run;
each BAM remains paired with its PBI. Demultiplex merge inputs use basename
order. These changes order filenames before the existing tools run.
The follow-on's completed [migration layout proof](optimization-validation.md#completed-migration-layout-proof)
confirmed both complete sample QS2 files and all 16 reviewed statistics TSVs,
including exact per-row name attributes after the existing chromosome-group
alignment. The original failed reports and a corrected validation-only check
are preserved. BAM and composed migration acceptance also passed.
The user approved this layout on 2026-10-08: comparison
with the original run may align these reviewed ordering differences, while
repeated optimized runs must preserve deterministic QS/TSV values, row order and
column order. Migration exceptions must not be used for new-versus-new checks.
The existing BAM coordinate-tie exception remains applicable to repeated runs;
it does not permit changes to records, multiplicities or coordinate ordering.
The original failed comparison reports remain intact; the
[validation ledger](optimization-validation.md#historical-pythonr-run-ordering-failure-and-subsequent-repair)
distinguishes the failed strict replay from the passing diagnostic.

Both samples have now passed two fresh output-stage executions, compared with
each other and with their actual `7e37ada` outputs. Across the four comparisons,
all 72 ordered QS component checks, 1,244 scientific-file comparisons and 224
own-index checks passed. QS/TSV values, types, attributes and row/column order
were exact, with no ignored QS fields or ordering normalization. Existing VCF
date/encoding and PDF date exceptions remain applicable. These are output-stage
replays, not two additional full-workflow runs. Scientific writers remain R;
the replay and validation controls use R, Python and shell.

`outputResults.R` skips expensive list-column conversions that the following
germline-table pivot discards. Its native `VariantAnnotation::writeVcf()` call uses
`nchunk = 100000L` to reduce temporary VCF formatting allocations. This is the
existing writer's export buffer; it does not divide analysis calls into separate
jobs or split the final QS2. The measured final object fits as a single object,
so no public split-QS2 format or new reader API is introduced.

The R-only output handoff installed in `3dabe21` uses a small R parent to run
the calculation/export worker, wait for it to exit, then start a fresh R
assembler. Private temporary QS parts transfer the existing ordered payload
between these sequential processes. Successful final rename publishes one
ordinary standalone QS2, then the private handoff is cleaned up. The existing
CLI, input-file ordering, scientific record order and Nextflow publication
modes are retained. Final assembly still materializes the complete ordinary
object; this does not provide constant memory or alter upstream burden memory.
The [complete-task evidence](#complete-r-only-single-public-qs-output-task-lib1-and-lib2)
records both samples, 42 lifecycle checks and the measured CPU/memory tradeoff.

## Follow-on R optimization results

Five integrated changes were assembled as candidate v3 against scientific baseline `ea8d06f3c0f29a7145e9e56dee1b0bb2424f4b77`: generated-field parsing, VCF target lookup, region-read aggregation, output object lifetime, and compression in R coverage annotation. The first three passed a combined extraction and filtering comparison on one real chunk from each library. Output lifetime and compression have separate component evidence. The complete `7e37ada` accounting and scientific acceptance are summarized [above](#historical-accepted-full-workflow-comparison-7e37ada).

### Selected changes and component measurements

| Change and measured scope | Implementation / non-R involvement | Actual CPU seconds, before → after | Process peak RSS, before → after | Scientific comparison |
|---|---|---:|---:|---|
| Generated-field parser; two alternating complete extraction pairs, LIB1 chunk 1 | R; no new custom non-R code | 414.84 → 279.75 (−32.56%) | 8.655 → 8.073 GiB (−6.73%) | 16 exact QS component checks, including baselines against the current pipeline; 16 parser fixtures |
| Coalesced VCF seeks; two alternating complete strict nuclear filter pairs, LIB1 chunk 1 | R changes shell command construction; existing Bash/`bcftools`, no new standalone executable | 1412.06 → 812.74 (−42.44%) | 12.315 → 12.359 GiB (+0.36%) | 28 exact QS component checks; 32 VCF fixtures |
| dtplyr region aggregation; grouping and ordered joins only | R; existing dtplyr/data.table compiled internals, no new custom non-R helper | 66.908 → 0.761 | No whole-filter memory claim | 432 exploratory fixture comparisons; all three alternative mappings exact on 124,090 reads and six region filters |
| Output object lifetime; one complete LIB1 output-stage pair | R; no new custom non-R code | 1930.72 → 1819.27 (−5.77%) | 59.529 → 43.284 GiB (−27.29%) | 311 scientific files passed, including all 18 QS components and run metadata |
| BGZF level 1 in R annotation; one complete nuclear coverage-row pair | R changes shell `bgzip` arguments; existing `bgzip`/`tabix`, no new standalone executable | 2004.91 → 1086.09 (−45.83%) | 15.859 → 15.858 GiB (effectively unchanged) | 1,509,632,319 ordered rows and 45,175,589,162 decompressed bytes exact; each output's own index validated |

The parser uses `tidyr::separate_wider_delim()` at 16 generated underscore-field sites. The measured historical version retained an older parser for colliding names and unusual data-frame attributes. Current production handles column placement, row names and grouping through one reconstruction path; it does not dispatch to that older parser. VCF lookup coalesces nearby seek intervals, then applies the original fine targets with matching overlap semantics before the unchanged normalization steps. Duplicate records, ordering, deletion overlap and target boundaries remain covered by exact fixtures.

Region aggregation uses `dtplyr::lazy_dt(immutable = TRUE)`, computes the original `base::mean` once per group, restores dplyr ordering, and retains the threshold expression and ordered join. This is a targeted conversion of the region-read grouping step. The selected dtplyr implementation took 0.761 seconds versus direct data.table's 0.711 seconds in the same operation-level comparison. Installed versions were R 4.4.2, tidyr 1.3.1, dtplyr 1.3.1 and data.table 1.17.6; no package installation was needed.

The output change saves the same single final QS object earlier and releases raw calls and coverage before VCF/export construction. Per-file traversal and formatted data remain unchanged. Failed tasks must still prevent publication, because the intermediate QS file now appears earlier. The output pilot used four completed historical candidate-optimization burdens inputs and their effective configuration. Its comparison allowed only VCF `fileDate` and PDF `CreationDate`/`ModDate` differences; it used no numeric tolerance or ignored QS metadata. Passing 311 scientific files does not imply identical compressed or PDF bytes.

Level-1 compression applies to BED output from the single R `annotate_coverage_row()` implementation. The tested R-annotation BED grew from 8,038,272,657 to 10,198,489,389 compressed bytes (+26.87%). The full-row timer includes packet loading, the original writer, annotation, compression and indexing. It is not a complete burdens-stage measurement. Compressor-only timings are attribution within that total and must not be added to it.

In the completed full workflow, the 36 published coverage BEDs grew from
151.27 to 192.15 GB (+27.0%). These are summed compressed file lengths in decimal
GB, excluding indexes and temporary files, rather than physical filesystem usage
or an isolated compression-only comparison. Final sample QS2 sizes were nearly
unchanged. This storage cost accompanies the CPU results above.

### Combined extraction and filtering

Jobs `19352943` and `19352944` completed successfully with frozen candidate v3. Each arm ran a fresh extraction worker and four separate filters: nuclear strict/lenient and mitochondrial strict/lenient. Both arms matched the current pipeline reference in all 128 QS component comparisons, with zero ignored fields and zero numeric tolerance. Source, executed-copy, input and producer bindings passed independent review.

| Library and scope | Actual CPU seconds, baseline → candidate | Highest individual worker RSS, baseline → candidate |
|---|---:|---:|
| LIB1 chunk 1; extraction plus four filters | 5173.631 → 3293.236 (−36.35%) | 14,826,940 → 14,299,480 KiB (−3.56%) |
| LIB2 chunk 1; extraction plus four filters | 4349.449 → 2790.702 (−35.84%) | 15,771,884 → 15,810,648 KiB (+0.25%) |

Across these two measured chunks, scientific-worker CPU totaled 9523.079 seconds before and 6083.938 seconds after, a 36.11% reduction. This sum includes 20 measured worker executions across the two arms and libraries. Exact comparisons and launch preflight are excluded. Arm order was baseline then candidate for LIB1, reversed for LIB2; filesystem cache state was uncontrolled.

Memory changes varied by operation. Strict nuclear filter peaks increased 2.63% for LIB1 and 2.58% for LIB2; LIB2 mitochondrial filter peaks also increased. These chains do not establish a general memory reduction. They exercise the parser, VCF and aggregation changes, and do not exercise output lifetime or coverage compression. They are not complete-workflow measurements. [Independent chain review](optimization-validation.md#follow-on-component-and-chain-validation).

### Rejected interval approaches

Candidate v3 retains the original subread interval construction. Flattened positions with native compressed reduction passed all 28 complete-filter component checks and reduced worker CPU by 9.22%, but increased process peak RSS by 8.94%. That whole-worker pair measured the earlier unguarded source; a separate integer-overflow guard later passed boundary and real-constructor tests. Dense failure masks remain a temporary-memory concern.

Both interval alternatives use R with existing Bioconductor operations; neither
would add a custom non-R production helper.

Direct logical masks passed 192 fixtures and four real constructor comparisons, but increased constructor CPU from 148.64 to 369.29 seconds for only 4.46% lower diagnostic-process peak RSS. Those CPU timings exclude loading and serialization, whereas the diagnostic RSS includes them. Neither interval approach is included in candidate v3.

### Interpretation and regression checks

CPU means actual user plus system time, including waited child processes for worker measurements. Component RSS is the individual process peak, averaged across arms where the table reports repeated pairs. Chain RSS is the maximum individual worker peak, never a sum or concurrent workflow memory estimate. Scheduler peaks may include separate comparison processes and are not substituted for worker peaks. The component gains are not additive and must not be extrapolated to all chunks or a complete pipeline.

Reusable fixtures passed in job `19355705`: 16 parser cases, 168 exact production region-aggregation/ordered-join cases, and 32 VCF target comparisons. The four test source files have three entry points, because the Python VCF driver prepares fixtures and invokes its paired R test. [Test invocation instructions](../scripts/benchmark/README.md#follow-on-regression-tests) describe the environment and repository commands. After integration, job `19374216` completed successfully in 23 seconds: all four integrated R expression trees matched the measured source, the relocated 16/168/32 fixture suites passed, and all 26 source pins and the SIF identity remained unchanged. Historical measured-source hashes and post-format integration hashes are retained separately; the identity check introduces no new performance measurement.

The [compact benchmark](benchmarks/torch-2026-10-07-followon.json) records the measured scopes and comparison rules. The [evidence ledger](optimization-validation.md#follow-on-component-and-chain-validation) retains source, configuration, receipt and report hashes. Full-workflow validation and accounting remain separate acceptance gates.

### Scaling limits

The current 60 chunks per sample divide work but do not cap chunk size as input grows. Extraction and filtering can therefore need more memory on larger inputs. Burdens processing still retains raw calls and germline data, and its final Rle memory depends on run complexity. Earlier object release reduces simultaneous live objects in output generation, but the final result remains one complete QS object. None of these measurements establishes constant memory with increasing input. The selected native accumulator is included in `7e37ada`; later temporary-storage experiments remain separate.

### Burden accumulator comparison and selection

The following comparisons describe the earlier accumulator selection. They do
not describe selectable alternatives in current production. The rejected
endpoint prototype and its capacity/legacy paths remain experimental evidence
only; `358e4bf` uses the generalized native R implementation.

Two R implementations have completed matched 60-chunk accumulation benchmarks.
Both use existing compiled Bioconductor operations; neither adds custom C++, Rcpp
or Python scientific code. Python and shell control the benchmarks only. These
measurements include chunk loading, ordinary coverage, sensitivity and result
serialization, but exclude raw-call/germline accumulation and coverage annotation.
Their baselines ran separately, so compare each candidate with its own baseline.

| Candidate and scope | Non-R involvement | Actual CPU seconds, before → after | Process peak RSS, before → after |
|---|---|---:|---:|
| Retain interval endpoints and finalize coverage once; 60-chunk accumulation | Existing compiled Bioconductor packages only | 10,417.00 → 3,773.84 (−63.77%) | 36.216 → 26.220 GiB (−27.60%) |
| Incremental native Rle addition with immediate sharing; 60-chunk accumulation | Existing compiled Bioconductor packages only | 11,741.48 → 9,473.69 (−19.31%) | 36.690 → 29.105 GiB (−20.67%) |

Both comparisons passed all ten exact component checks. These component gains
are not whole-stage or whole-pipeline results. The original 8-chunk pilots showed
different gains and are not extrapolated to the full workload.

The complete incremental Rle worker has now also passed: all 20 QS components,
including configuration and run metadata, are identical; all eight compressed
coverage BEDs are byte-identical, and all eight own-file index proofs and indexed
queries pass. This uses custom R code with existing compiled Bioconductor
operations. Complete-worker actual CPU was 36,149.475 → 29,084.372 seconds
(10.042 → 8.079 hours, −19.54%). Observed Slurm batch peak memory was
120.258 → 112.563 GiB (−6.40%). Validation used another 1,906.119 CPU seconds,
excluded from worker costs.

This complete-task comparison is unpaired across nodes and uses a 288 GiB
baseline allocation versus a 160 GiB replay allocation, both with two CPUs.
Torch's cgroup memory accounting can include filesystem cache, so the observed
peak difference is not an isolated estimate of R allocation savings. The final
combined pipeline comparison retains matched original resource requests. The
reviewed complete-task proof is workspace
`runtime/exploration-burdens-v2/whole-worker-native-postrun-v1/review.json`, SHA-256
`c25289659decc7e1cc342cc380704b098797112947894c8342b07205c66b40eb`;
its resource-context supplement is
`runtime/exploration-burdens-v2/whole-worker-resource-context-v1/review.json`.
The baseline here is `ea8d06f`, not the original `6b6598b` pipeline. The endpoint
worker and its strict validation have also completed. All 20 QS components,
including configuration and run metadata, match exactly; all eight compressed
BED files, decoded TBI payloads, own-index checks and indexed-query proofs pass.
Its observed actual CPU was 30,148.661 seconds (8.375 hours, 16.60% below that same
baseline), with a Slurm batch peak of 108.373 GiB. The separate process/descendant
maximum was 35.066 GiB for endpoints and 41.486 GiB for native; those are not the
Slurm metric or concurrent process totals. Endpoints used 3.66% more CPU than
native in these unpaired runs, with a lower observed process peak. The comparisons
do not isolate causal timing or memory differences between alternatives. Both
implementations use custom R and existing compiled packages/CLI tools. The
endpoint review is workspace
`runtime/exploration-burdens-v2/whole-worker-endpoints-postrun-v1/review.json`,
SHA-256 `ab2fd7839b9fa0ec64ca1abf06f465f2c801469d9e9640713379d16b13156554`.
Its separate validation cost was 2,128.765 CPU seconds, excluded from worker costs.

Endpoint storage grows with retained intervals. It also needs a guard against
concatenating more than 2,147,483,647 ranges for one membership group and
chromosome: the installed IRanges version uses an integer range count at that
boundary. The present 60-chunk dataset stays below 1.46% of that limit even when
summing all groups for its largest chromosome. The proposed R-only fallback uses
coverage of the already retained individual pieces and adds the resulting Rles
above the library limit, preserving the ordinary path otherwise. It introduces
no call batches, output shards, chromosome jobs or new tuning parameter.

The capacity guard passed 32 boundary checks using a tiny forced limit, all 126
previous endpoint cases and the existing coverage/sensitivity suites. These
fixtures do not allocate billions of intervals. The matched eight-chunk
normal-path benchmark also passed all ten exact component comparisons. Worker
CPU was 760.33 seconds without the guard and 701.39 seconds with it; peak RSS was
6,778,328 and 6,766,412 KiB, respectively. This single pair showed no resource
regression, but most timing variation occurred before the added helper ran; it
does not establish a speed improvement from the guard. Validation costs are
excluded from those worker measurements.

The guard avoids an aggregate library limit; it does not bound the retained
endpoints or make total memory constant. Neither the guard nor the alternative
accumulator has been selected for production.

The native alternative retains aggregate Rles and clears its incoming-coverage
and new-Rle sharing caches between chunks. Endpoints instead retain interval
history until finalization and add membership grouping, a capacity path and a
permanent legacy fallback. Any observed asymmetric barcode orientation on either
strand sends ordinary endpoint coverage to that legacy path, even when the
category does not request orientation summaries. A late transition materializes
the retained history while the triggering chunk is live. Native keeps its
ordinary coverage optimization for those configurations; both candidates apply
independent guards to sensitivity calculations. The observed orientation reflects
the relevant final barcode round, so eligibility is not simply a round1/round2
setting. These are implementation properties, not measured performance claims
for untested barcode configurations.

Native accumulation was selected after these complete comparisons. Its simpler
chunk-local lifecycle, lack of retained interval history and broader ordinary
barcode-path applicability outweigh the endpoint candidate's lower observed
peak for this dataset. Native also used less observed whole-worker CPU, although
the small unpaired difference does not establish a causal speed advantage.
Neither design makes whole-pipeline memory constant: raw calls, germline records,
Rle run complexity and the final single QS2 still grow with the input.
Native has been installed after the current-shared integration preflight passed
in job `19421102`: all seven parsed-source comparisons, 70 native cases and the
coverage/sensitivity suites, with actual chrM/chrY exercised in the coverage
suite. Three explanatory comment
blocks were clarified after installation; all other source lines remain identical
to the tested source. The combined `7e37ada` workflow completed execution,
accounting and scientific acceptance under the original resource requests.
The selection is recorded
in workspace `runtime/exploration-burdens-v2/accumulator-selection-v1/selection.json`.
Both implementations are R with existing compiled packages and unchanged external
annotation tools; neither adds custom non-R scientific code.
The [compact comparison report](benchmarks/torch-2026-10-08-burden-accumulator.json)
keeps measured-source, installed-source and integration evidence distinct.

### Experimental temporary coverage storage

A separate R prototype stores coverage groups temporarily as QS files, reduces
them, then assembles the ordinary self-contained final QS2 in a fresh R process.
Calls, statistics and sensitivity retain their existing calculations. It uses
existing compiled qs2/Bioconductor packages and annotation/index CLIs; custom
Python/shell code controls the benchmark only. This approach is **not adopted**,
and the production single-QS2 output remains unchanged.

For one complete LIB1 strict nuclear burden task using the same 60 inputs,
allocation CPU was 6.82530 hours versus 7.54434 for the native task, including
temporary I/O and fresh assembly. All 20 ordered QS components, eight compressed
coverage BEDs and their index proofs passed exact comparison. These descriptive
runs occurred separately; they are neither a paired causal estimate nor a
whole-pipeline result. Earlier 5/15/60-input coverage-only pilots omitted other
burden work, and the five-input temporary-storage pilot was more costly; their
results should not be substituted for this complete-task comparison. The full
prototype uses the native route for eight or fewer inputs; the five-input pilot
deliberately forced storage to measure its overhead.

Observed Slurm batch peaks were 98.723 GiB for the prototype and 128.668 GiB for
native, with different nodes/times and 160 versus 288 GiB memory requests.
Torch's cgroup accounting includes file cache. Native has no isolated R-process
peak comparable to the prototype's process measurements, so these figures do
not establish a causal R-memory reduction or a safe lower allocation. Stat-only
sampling observed a maximum of about 38.0 GiB across eight temporary raw BED files; that is
a sampled lower bound for those files, not total peak scratch usage.

The ordinary output consumer completed from a standalone copy of the final QS2
with the private parts hidden, producing 311 outputs. Its R worker used
2,098.585498 CPU seconds and peaked at 45,114,820 KiB (43.025 GiB). The separate
consumer comparison passed all 18 ordered QS components, all 311 exports and
56 own-file index checks. QS configuration, metadata and TSV order remained exact;
existing VCF date/encoding and PDF Info-date exceptions applied. Native output
replays had essentially the same command memory peak. This proves compatibility
for the tested consumer, with a substantial downstream memory requirement still
present.

Groups of eight limit input count, not bytes. One input QS2, the retained calls,
complete final assembly and later sample-level output still impose growing
memory requirements. The prototype retained 3.55 GB in 1,760 private files and
adds 455 lines of reduction/storage logic before production lifecycle integration.
The evaluation therefore weighs the further 9.53% observed task CPU saving
against retry, cleanup and maintenance costs. The recommendation is to retain
the current native default and single QS2. No production adoption or larger-input
memory guarantee follows from these experiments.

The first disk-backed call-table feasibility attempt stopped after Arrow failed
while applying table metadata and DuckDB failed startup configuration. That
did not establish that either backend was unsuitable. The reopened evaluation
is now complete: removing DuckDB's invalid configuration option permits startup, and
representing three row-associated name attributes as native columns permits
ordinary Arrow filtering. Both complete actual call tables now reconstruct
exactly, including their list columns and R attributes. Across both backends,
all 84 tested operation outputs passed on retained call tables grouped by the
selected 5, 15 and 60 analysis chunks. These were queries of retained tables,
not fresh burden calculations.
The tested Arrow joins join scalar row IDs, then collect and reattach list
columns in R. Native joins carrying those list columns remain unsupported in
the installed Arrow version; DuckDB supports the tested native list joins.
DuckDB used less CPU and peak process memory than Arrow in the tested
read/query/materialization suites. Both read prebuilt storage; the R baseline
loaded QS2 during its suite, and staging costs were measured separately. These
results favor DuckDB over Arrow for the tested disk read/query path, without
establishing faster resident-R operations or complete-burden savings. No further
full Arrow burden rewrite is warranted by this evaluation. The complete DuckDB
burden/output prototype uses R with custom SQL through the existing compiled
DuckDB/DBI packages; it adds no custom C++ or Rcpp. Python and shell control the
experiments outside the production repository.

A complete one-input burden replay passed all 20 ordered components and all
eight coverage BED/index/query checks. Its downstream output comparisons passed
all 18 final components exactly and 140 exports under the existing date/encoding
exceptions, both through disk-backed readers and after reconstructing an
ordinary QS2 for the unchanged consumer. That complete
path still peaked at about 20.32 GiB in both implementations. It therefore does
not demonstrate a whole-task memory reduction. The five-input burden and
its ordinary-QS assembly also pass all 20 components. Complete worker CPU
increased 44.5%, with essentially unchanged process peak memory; assembly adds
a separate 29.65 CPU seconds.

The fifteen-input pair and ordinary burden-QS assembly are complete. All 20
components match exactly both through the backed reader and after assembly.
The ordinary worker used 5,879.13 CPU seconds (1.63 hours), versus 16,238.18
(4.51 hours) for DuckDB: 176.2% more CPU. Both reached about 20.32 GiB process
peak, with greater sampled file-cache pressure in the database arm. Ordinary
assembly adds 111.02 CPU seconds and peaks at 7.55 GiB; its full burden path
therefore costs 16,349.20 CPU seconds. This single native-first pair used one
CPU per worker on the same node, with uncontrolled cache. It does not include
an external BED/index comparison or final consumer at this input size.
The separate 8 GiB-pool probes are complete and exact at five inputs
(211,430,455 coverage rows) and fifteen inputs (**673,395,946 ordered coverage
rows**). The fifteen-input source contained 707,349,842 events. Its isolated
query used 4,413.557 CPU seconds and 4,714.940 wall seconds; the complete fresh
process used 4,415.574 CPU seconds and peaked at 14.93 GiB RSS. Sampled temporary
file occupancy reached 46.735 GB, while cumulative cgroup peak was 60.446 GB.
Validation costs are separate. The effective 8.0 GiB setting limits the engine
pool, not the whole process or job. This isolated query and the earlier broader
reduction phase have different timing boundaries and cache/process states;
they do not establish a causal tuning speedup or a measured complete 8 GiB
worker. The 60-input successor remains disabled by the completed decision to
reject this SQL integration.

The separate four-source output test is complete. It starts from four accepted
ordinary burden results migrated to native storage, rather than four upstream
backed calculations. All 18 final components match through both the backed
reader and assembled ordinary QS2; all 310 exports and their own-file index
checks pass under the existing date/encoding exceptions.

| Four-source output path | Actual CPU seconds | Maximum stage process RSS, GiB |
|---|---:|---:|
| Ordinary consumer | 1,444.75 | 43.02 |
| Backed consumer, without ordinary final-QS assembly | 1,938.51 | 6.59 |
| Backed consumer plus ordinary final-QS assembly | 2,630.58 | 40.37 |

The backed consumer alone uses 34.18% more CPU with 84.68% lower process peak.
Reconstructing the self-contained QS2 adds 692.08 CPU seconds and largely
restores the observed process peak. First-use migration of all four
ordinary inputs adds another 2,180.19 CPU seconds separately; it would not be
charged to already-backed inputs. The final backed result requires a 12.81 GB
closure of database and small QS files, versus 2.57 GB for the ordinary final
QS2. Both also publish 3.44 GB across the same 310 external files. This compares
the final logical-data storage only, excluding retained ordinary burden
intermediates and peak scratch; it is not a whole-pipeline storage ratio.
These are single paired consumer measurements with native-first ordering and
uncontrolled filesystem cache, not whole-workflow measurements. The actual
60-input SQL-backed burden, ordinary burden assembly and downstream consumer
chain have not run and will not be activated for this finalization. The
complete five- and fifteen-input worker CPU increases, without a process-peak
improvement, do not justify the current integration's CPU, cache and I/O costs.
The observed 20.32 GiB process peaks are neither a proven memory floor nor an
attribution to input loading.

The experimental database memory setting is not a cap on the entire R task.
Filesystem cache and temporary query files also contribute to job memory and
I/O; process RSS and sampled cgroup accounting are recorded separately.
Input deserialization, R-side calculations, collected tables and ordinary QS2
assembly all require additional memory. Keeping a self-contained QS2 preserves
an output-sized allocation during assembly and later ordinary reads. A small
QS2 containing database references would change the storage contract and require
retaining the database and using its readers; that is a separate evaluated
option, not the production output format.

The prototype's core burden/export helpers contain 509 physical lines of R/SQL;
the accepted-output migration utility adds 37. Generated integration changes
and an unchanged scientific fragment reused by the exporter are separate.
The attribute adapter accounts for only a small part of that total. Production
adoption would also require complete handling of supported barcode configurations,
backing-file publication and restart rules, format versions and cleanup. Neither
exact native-column storage nor a successful dataset replay removes that work.

R-native coverage combined with native-column DuckDB calls remains technically
feasible, but that combined path is unvalidated and deferred. It would retain
custom SQL ordering/query logic for calls, without adding custom C++, and would
need combined readers/manifests, track mapping, lifecycle handling, ordinary
assembly and all-configuration validation. Earlier integer-only reduction
proofs do not permit reassociating fractional or double sums. Neither this
hybrid nor a few-QS handoff removes whole input and final-QS allocations. Current
finalization retains clean R calculation and a self-contained QS2. The separate
R-only single-public-QS output task has passed its LIB1 completion and LIB2
integration/repeat proofs and was installed in `3dabe21` after acceptance of
the `358e4bf` workflow. Its separate output-stage proofs do not establish a full
workflow at the later source revision. The final storage decision
and reviewed boundaries are in the
[completed storage evaluation](benchmarks/torch-storage-evaluation.json).
The readable companion and underlying reviews remain in workspace
`runtime/reopened-storage-evaluation-20261009/`:
`storage-evaluation-final-v1.txt`, `storage-evaluation-final-independent-review-v1.json`
and `coverage-pool-probe-terminal-root-review-v2.json`. No new storage rewrite
or allocation is proposed.

The [first evaluation](benchmarks/torch-2026-10-09-storage-evaluation.txt)
and [compact evidence](benchmarks/torch-2026-10-09-storage-evaluation.json) retain
the earlier R-coverage storage prototype's 5/15/60-input results, phase
measurements, exactness proofs and limitations.

### R-only private components with ordinary final-QS reconstruction

A separate downstream experiment used the same four ordinary burden inputs and
saved each annotated coverage portion before releasing it, plus one file holding
the other 17 final components. Both the component reader and reconstructed
ordinary QS2 reproduce all 18 components exactly; all 310 exports and 56
own-file index rebuilds pass. This tests ordinary R objects using existing
qs2/dplyr/Bioconductor packages, without SQL, a custom codec or new compiled code.
Python/shell are experimental measurement and validation controls only.

| Downstream boundary | Actual CPU seconds | Maximum stage process RSS, GiB |
|---|---:|---:|
| Ordinary final QS2 and exports | 2,012.88 | 43.02 |
| Five private payload files and exports | 1,891.31 | 26.58 |
| Private files plus fresh ordinary-QS reconstruction | 2,080.26 | 32.02 |

The first two rows are one native-first pair on the same node, with uncontrolled
cache/order effects. Reconstruction ran in another allocation on another node;
the last row adds its measured 188.96 CPU seconds and reports the larger process
peak. The resulting +3.35% CPU and -25.57% process peak are descriptive sums,
not a matched complete-task result or a required-allocation measurement.
The required five-file bundle plus manifest occupies 2,573,368,484 bytes, only
338,067 bytes more than the ordinary QS2. Retaining both representations doubles
that final-data storage; the prototype does not clean up the parts.

This earlier few-QS variant motivated the distinct complete-task experiment
below. Its component bundle and separately measured reconstruction remain
separate evidence; the descriptive sum above is not the completed task's
measurement. No public component format is adopted. Ordinary input loads and
complete final assembly still grow with data size; selecting a non-coverage
component in this earlier experiment loads all 17 remainder components.
Source-bound results and independent review are in workspace
`runtime/native-few-qs-output-replay-v1/postrun-review-v1/`.

### Complete R-only single-public-QS output task: LIB1 and LIB2

The LIB1 completion experiment passed retained candidate A versus native,
fresh candidate B versus matched native B, and B versus A. Each comparison
passed all 18 QS components with exact values, types, attributes and order,
all 310 exports under the existing date/encoding policy, and 56 own-index
rebuild checks. All 42 repository-relative lifecycle tests also passed in an
isolated checkout. Source, input, reference, SIF and product-identity guards
passed; the failed earlier attempt remains preserved.

| LIB1 complete paired output command | Actual CPU seconds | Largest-process high-water memory, GiB |
|---|---:|---:|
| Ordinary native B | 1,488.797765 | 43.022 |
| Single-public-QS candidate B | 1,599.939114 | 32.020 |

The candidate used **7.465% more CPU** with **25.572% lower largest-process
high-water memory**. Its entire lightweight R parent, calculation worker,
fresh assembly, required I/O and cleanup are included; there is no omitted
assembly cost to add. This is one same-node candidate-first pair with
uncontrolled cache/order effects. Candidate concurrent-memory sampling was
incomplete, so the individual-process high-water values do not establish a
complete concurrent-memory bound or a minimum allocation. Synchronous execution
and separate role logs establish that the calculation child exited before
assembly began. Allocation and validation overhead remain separate.

The LIB1 candidate and native QS2 files each occupy 2,573,030,417 bytes (2.573 GB), with
about 3.436 GB of common exports retained separately. Private parts are absent
after successful completion. Ordinary final assembly still materializes the
complete object; this is neither an upstream/whole-workflow memory improvement
nor a constant-memory claim.

LIB2 integration is also complete: candidate A against the actual `358e4bf`
native output and candidate B against A each passed all 18 ordered components,
310 exports and 56 own-index rebuild checks. All source, input, reference, SIF,
producer and product guards passed. Both candidate concurrent-memory samplers
were incomplete. Their individual-process peaks cannot be compared with native
Slurm batch memory, and the separate native task is not a matched timing arm.
This proves output-stage integration and repeatability for the current LIB2
inputs; it does not establish a later-source whole-workflow result.

The exact R-only candidate and 42 lifecycle checks were installed in `3dabe21`
after both-sample validation and acceptance of the `358e4bf` workflow. It uses
existing R packages and unchanged
export CLIs, with no new package or custom non-R scientific logic; Python/shell
remain measurement and validation controls. The
[combined single-QS report](benchmarks/torch-single-qs-output-handoff.json)
records both sample proofs, measured costs and source hashes. Underlying reviews
remain in workspace `runtime/native-single-qs-completion-v1/postrun-review-v1/`
and `runtime/native-single-qs-lib2-replay-v1/postrun-review-v1/`. LIB1 used
historical burden inputs; LIB2 used the four current `358e4bf` burden outputs.
Neither constitutes a full workflow executed with the candidate.

The 42-check lifecycle fixture was verified in an isolated checkout. In a
checkout containing the candidate and fixtures, run from its root on a compute
node in the pipeline R/Bioconductor environment, using a fresh output directory:

```sh
Rscript --vanilla tests/test_output_results_handoff.R /absolute/fresh-output
```

## Focused local validation

Coordinate-sort reuse was tested and rejected. The single-input fixture matched,
but every tested multi-input order (`AB`, `BA`, `ABC`, `CBA`) changed SAM record
order at coordinate ties. Record multisets and headers excluding program records
matched, and indices validated. That experiment preceded approval of the more
permissive BAM coordinate-tie comparison rule; it was rejected under the
then-current strict record-order rule. The agreed implementation retains the
existing coordinate sorts and `pbmerge`, with
no single-input exception. The retained fixture report is workspace
`runs/coordinate-merge-fixtures/results.json`.

Run `python3 -m unittest discover -s tests -p test_artifact_cache.py -v` from the
repository. Tests exercise identity changes, exact serialized identity handling,
concurrent publication, failure cleanup, cache hits, corruption, and symlink
rejection. Workflow/container smoke tests and scientific equivalence benchmarks
are separate requirements before production use.


`tests/test_prepared_cache_workflow.py` is a compute-only integration fixture for
actual region and germline-coverage preparation processes. It extracts their
production definitions, identity construction, `cachedBuild`, and completion
barrier, then uses tiny real BigWig tracks. The protocol checks cold builds, warm
`-resume` restores, a rebuild after renaming one entry inside its own test cache,
and concurrent launches sharing a fresh cache with separate work directories.
It counts real tool executions, verifies complete manifests and unchanged
unrelated entries, and compares imported BigWig/GRanges scientific objects.
Two individuals share one coverage input, with a duplicate threshold and a
fractional threshold to exercise scoping and identity serialization.

Fixture job 19194716 passed all phases: eight cold builds, zero warm-resume builds, one rebuild for the renamed entry, and eight total builds across two
concurrent workflows with 16 preparation tasks. All scientific objects and
manifest checks passed. Evidence is in
`runs/prepared-cache-workflow-v1/results.json`, the per-phase traces, and
`scientific.log` files. This exercises actual cache publication and restoration
without replacing the separate full-reference workflow validation.

## Python/R follow-up status

The follow-up was authorized on 2026-10-06. The rebuilt-container Python/R changes and the five additional R optimizations have component evidence. The parser, VCF lookup and targeted dtplyr changes also passed combined extraction/filter chains on one real chunk from each library, with 128 exact QS component comparisons and 36.11% lower summed scientific-worker CPU. Memory results were mixed. Output lifetime and R-annotation compression retain their separate component measurements.

The `ea8d06f` workflow completed all 650 tasks, its full comparison and CPU/memory
accounting, but its LIB1 final QS2 failed the original strict comparison; that run
is not scientifically accepted. The combined run is pinned to immutable revision
`7e37ada`; execution, accounting, cache checks and both deterministic output-stage
replays are complete, and all 711 required publication paths are accepted. It includes
the five later changes and selected native burden accumulator, but predates the
allele-key helper installed in `00f6e33`. The earlier `3911047` acceptance remains
a separate historical result. Retained benchmark JSON files are
historical snapshots: their earlier pending or under-evaluation status text does
not supersede this current status or the native-accumulator selection above.

1. **Python BAM splitter:** implemented, with exact full LIB1 comparisons and
   measured splitting/indexing costs, including the final system-Python/pysam
   installation in the rebuilt container. This adds custom Python and uses
   pysam's existing compiled library and existing index tools.
2. **R artifact cache:** implemented with common functions in
   `sharedFunctions.R`. Protocol and actual Nextflow cold/warm/concurrent tests
   passed. Paired closed-product cache operations measured its additional CPU
   and memory cost; unchanged scientific preparers were outside that benchmark.
   It uses existing Linux commands and Nextflow/shell invocation, with no custom
   non-R cache helper.
3. **R coverage annotation:** implemented in `sharedFunctions.R`; chr22 and full
   nuclear comparisons passed, including downstream context consumers and
   indexes. CPU and memory costs versus C++ are recorded above. Existing compiled
   R packages and `bgzip`/`tabix` remain; no custom C++ annotation helper remains.
4. **data.table/dtplyr evaluation:** `fread`/`fwrite` handle coverage I/O. The earlier two pure-R filtering changes retained their measured 5.10% mean complete-worker CPU reduction. Subsequent profiling selected an immutable dtplyr grouping step for region-read filters, preserving `base::mean`, threshold expressions and ordered joins. Combined extraction/filter results are recorded in the [follow-on section](#follow-on-r-optimization-results). A later real-interval ABBA comparison also evaluated `data.table::foverlaps` and non-equi joins. On the largest captured subject, `foverlaps` increased mean operation CPU from 10.38 to 58.31 seconds and mean whole-process peak RSS from 3.36 to 11.71 GiB; the non-equi alternative also regressed. Smaller cases were mixed, so neither became a broad interval-engine replacement. These are captured-operation results, not whole-filter performance claims. All implementations use R with existing compiled package internals, without a new custom non-R scientific helper; Python/shell controls remain outside production.
5. **Rcpp evaluation:** profiling did not identify a further candidate with
   demonstrated dramatic gains after removing the selected R overhead. No new
   custom compiled helper was added for incremental gains. An Rcpp implementation
   would involve custom C++ called from R; none is selected here.
6. **Incremental burden accumulation:** implemented in R using existing compiled
   Bioconductor operations, without a custom non-R kernel. A complete strict
   nuclear task and current-shared integration tests passed. The retained-endpoint
   R alternative also passed, but was not selected because its extra state and
   growth with interval history did not justify its observed memory advantage.

For all these items, retain the agreed scientific comparison rules and record
performance regressions as well as improvements. Judge aggregate performance by
actual user-plus-system CPU-hours across the pipeline, with significantly lower
peak memory in the currently expensive steps. An individual step may become
more than 5% slower; there is no per-step 5% acceptance limit. Preserve the single
final QS2 and the existing TSV, VCF, PDF and indexed BED outputs.
