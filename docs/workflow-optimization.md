# Computational optimizations

The optimization branch changes repeated work and intermediate allocations while
retaining the scientific output schemas and one final QS2 per sample. TSV, VCF,
PDF and indexed coverage BED outputs are retained. The
[validation ledger](optimization-validation.md) records completed component
benchmarks, complete-workload measurements and acceptance of all 711 required
scientific outputs under the agreed comparison rules for the earlier
`3911047` production revision. The Python/R follow-up described below is being
validated separately; those historical results do not establish acceptance of
the current follow-up code.

## Implementation languages

The labels below describe production implementation, including orchestration.
An R implementation can use existing packages with compiled internals without
adding a custom C++ or Rcpp helper. Existing external tools are listed separately.
Benchmark and regression infrastructure also uses Python and shell; the VCF
regression test, for example, has a Python fixture driver and an R comparison.

| Optimization | Custom implementation | Non-R involvement |
|---|---|---|
| Prepared reference summaries and germline coverage caches | R preparers; Nextflow and shell integration | Nextflow/shell orchestration and existing reference/coverage tools; no new compiled helper |
| Artifact-cache protocol | R, including shared functions | Existing `flock`, `stat`, `sync` and `cp`; Nextflow/shell invokes the R entry point |
| Output publication with hard links or copies | Nextflow/Groovy | **Yes: custom non-R orchestration** using existing JVM/Nextflow filesystem APIs; no new scientific kernel |
| BAM dispatch | Python using `pysam`; Nextflow/shell integration | **Yes: custom Python**; existing BAM indexing tools and pysam's compiled library |
| Extraction quality lookup and generated-field parsing | R | Existing R packages; no new custom non-R implementation |
| Chromosome restriction, germline summaries, threshold and coordinate helpers, region-read aggregation | R | Existing R/Bioconductor packages, including dtplyr/data.table; no new custom non-R implementation |
| Coalesced VCF lookup | R constructs the query commands | Existing shell pipeline and `bcftools`; no new Python/C++ helper |
| Coverage accumulation and sensitivity calculations | R | Existing R/Bioconductor packages; no new custom non-R implementation |
| Coverage annotation and lower compression level | R | Existing `bgzip` and `tabix`; legacy annotation also uses existing `bedtools` |
| Final table formatting, VCF writer buffer and earlier object release | R | Existing R/Bioconductor packages; no new custom non-R implementation |
| Coordinate-sort reuse, rejected | Nextflow/shell command changes | Existing `pbmerge`/`samtools`; would involve non-R orchestration |
| Shared-session filtering, rejected | R worker plus workflow changes | Would involve Nextflow/shell orchestration |
| Flattened interval positions and direct logical masks, rejected | R | Existing Bioconductor packages; no new custom non-R implementation |
| Endpoint accumulation, under evaluation | R | Existing Bioconductor packages; no new custom Python/C++/Rcpp implementation |
| Native Rle addition and immediate sharing, under evaluation | R | Uses existing Bioconductor compiled operations; no new custom Python/C++/Rcpp implementation |

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

The completed candidate ran production revision
`391104794588819601f89bf6acfe588fc6c579bc`. The branch subsequently consolidated
the three reference-summary functions into `sharedFunctions.R` without changing
their definitions; focused container checks passed in job `19289528`. That source
organization change was not another full-workload benchmark. To reproduce the
fully benchmarked source exactly, set `HIDEFSEQ_GITTAG` to the full revision
instead of the branch name. Final
validation passed under the approved BAM rule: logical records, multiplicities,
coordinates and coordinate ordering remain exact; only relative ordering at
identical coordinates may differ. See the validation ledger for the complete
comparison scope and preserved original strict-order report.

The two-sample, 60-chunk-per-sample Torch comparison used the existing
`hidef-seq_3.0.sif` and matched these allocations between baseline and candidate:

| YAML process suffix | Memory | Time | CPUs |
| --- | ---: | ---: | ---: |
| `extractCallsChunk` | `32GB` | `1h` | 1 |
| `filterCallsChunkChromgroupFiltergroup` | `48GB` | `2h` | 1 |
| `calculateBurdensChromgroupFiltergroup` | `288GB` | `30h` | 2 |
| `outputResultsSample` | `96GB` | `4h` | 1 |

Set the corresponding `mem_<suffix>` and `time_<suffix>` keys in the run YAML;
CPU counts above come from the workflow process definitions. The template's
64 GB burden/output requests are not validated for this workload. Candidate
strict nuclear burden tasks reached about 120–128 GiB maximum Slurm RSS, and
final output tasks reached about 69–76 GiB. Requests were kept unchanged for the
comparison. Size future requests from representative full tasks with headroom;
these peaks are workload-dependent and do not establish a universal lower limit.
Slurm task RSS, Nextflow's sampled process RSS, and simultaneous pipeline memory
are different measurements.

## Prepared caches and publication

Prepared reference libraries, reference trinucleotide BEDs, per-individual VCF
annotations, germline BAM coverage/calls, and thresholded region tracks now use
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
filtered-base statistics are preserved. Legacy YAML without this mapping still
uses the original inline calculation.
The existing empty-stream failure is preserved: when a threshold comparison
produces no intervals, `wigToBigWig` rejects the empty input in both paths.

A reference-summary bundle is prepared after BSgenome installation and reused
across samples, chunks and chromosome groups. It contains full-reference N
intervals with their Seqinfo and integer trinucleotide counts for each chromosome.
Filtering retains the original whole-genome N-base statistic; burden tables sum
the requested chromosomes and keep the existing channel/fraction transformations.
Circular chromosomes retain the original linear counting boundaries. Standalone
YAML without `reference_summary_file` retains the original reference scans.
`tests/test_reference_summary.R` checks empty/N-free/all-N references, chromosome
order, repeated selections, types, and equivalence of the earlier chromosome
restriction in call loading.

These reference-summary functions live in the single common
`bin/sharedFunctions.R` file. The separate helper file has been removed, and its
production callers, benchmark and cache dependency list use the common file.
The consolidation changes that file's hash, so affected prepared artifacts and
downstream tasks receive new cache identities on the next run. Existing cache
bundles are retained; they are not silently reused under changed source hashes.

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
and source hashes are retained in `runs/r-cache-v1/result.json`; complete
pipeline accounting for the follow-up remains pending.

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
fields retains the legacy cache layout. Conflicting prepared product basenames
fail early instead of sharing an ambiguous cached file. Scientific product
schemas and basenames are unchanged; absolute cache paths are metadata.
Because the immutable effective configuration path is included in R task scripts,
changes to any retained configuration field can invalidate downstream R tasks.
The existing configuration signatures therefore do not yet provide fully scoped
downstream resume. This conservative behavior preserves concurrent-run safety;
separating execution parameters from full provenance is a future design change.

Output publication retains files in Nextflow work directories. It uses hard
links when work/output directories are on the same filesystem and copies
otherwise. This replaces final-output `move` publication, which removed files
required by resume. The existing global preparation barrier remains in place.

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

## Processing inside the R scripts

### Extraction

`extractCalls.R` decodes the run-length encoded `sa` tag directly into an Rle and
reverses its runs when strand orientation requires it. It no longer expands all
read-level `sa` tags into dense arrays for quality lookup. Queried positions are
mapped to run values using cumulative endpoints; unusual indices use the original
dense-vector semantics. The published extraction object still contains its
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

The original annotation path remains available through `--coverage-annotation
legacy`; compressed FASTA and reference contig names that the legacy parser
treats ambiguously also use it. Ordinary workflow and standalone calls default
to R annotation. Annotation, compression and indexing failures stop the task
before results can be saved. Fixtures cover exact context arithmetic across
buffer boundaries, output formatting and ordering, empty output, reference
edges and failure propagation. The full nuclear component comparison has passed;
the full follow-up workflow comparison remains pending. The earlier 51.14%
nuclear CPU reduction in the validation ledger applies to the removed C++
implementation, not the R port.

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
`runs/r-coverage-v1/component-summary-after-nuclear.json`; full-pipeline CPU,
peak-memory and scientific acceptance still require the follow-up run.

### Final output

`outputResults.R` skips expensive list-column conversions that the following
germline-table pivot discards. Its native `VariantAnnotation::writeVcf()` call uses
`nchunk = 100000L` to reduce temporary VCF formatting allocations. This is the
existing writer's export buffer; it does not divide analysis calls into separate
jobs or split the final QS2. The measured final object fits as a single object,
so no separate QS2 components or new reader API are introduced.

## Follow-on R optimization results

Five integrated changes were assembled as candidate v3 against scientific baseline `ea8d06f3c0f29a7145e9e56dee1b0bb2424f4b77`: generated-field parsing, VCF target lookup, region-read aggregation, output object lifetime, and compression in R coverage annotation. The first three passed a combined extraction and filtering comparison on one real chunk from each library. Output lifetime and compression have separate component evidence. Full-workflow acceptance and whole-pipeline CPU and memory results remain pending.

### Selected changes and component measurements

| Change and measured scope | Implementation / non-R involvement | Actual CPU seconds, before → after | Process peak RSS, before → after | Scientific comparison |
|---|---|---:|---:|---|
| Generated-field parser; two alternating complete extraction pairs, LIB1 chunk 1 | R; no new custom non-R code | 414.84 → 279.75 (−32.56%) | 8.655 → 8.073 GiB (−6.73%) | 16 exact QS component checks, including baselines against the current pipeline; 16 parser fixtures |
| Coalesced VCF seeks; two alternating complete strict nuclear filter pairs, LIB1 chunk 1 | R; existing shell/`bcftools` commands | 1412.06 → 812.74 (−42.44%) | 12.315 → 12.359 GiB (+0.36%) | 28 exact QS component checks; 32 VCF fixtures |
| dtplyr region aggregation; grouping and ordered joins only | R; existing dtplyr/data.table packages | 66.908 → 0.761 | No whole-filter memory claim | 432 exploratory fixture comparisons; all three alternative mappings exact on 124,090 reads and six region filters |
| Output object lifetime; one complete LIB1 output-stage pair | R; no new custom non-R code | 1930.72 → 1819.27 (−5.77%) | 59.529 → 43.284 GiB (−27.29%) | 311 scientific files passed, including all 18 QS components and run metadata |
| BGZF level 1 in R annotation; one complete nuclear coverage-row pair | R; existing `bgzip`/`tabix` commands | 2004.91 → 1086.09 (−45.83%) | 15.859 → 15.858 GiB (effectively unchanged) | 1,509,632,319 ordered rows and 45,175,589,162 decompressed bytes exact; each output's own index validated |

The parser uses `tidyr::separate_wider_delim()` at 16 generated underscore-field sites. Colliding names and unusual data-frame attributes retain the original parser. VCF lookup coalesces nearby seek intervals, then applies the original fine targets with matching overlap semantics before the unchanged normalization steps. Duplicate records, ordering, deletion overlap and target boundaries remain covered by exact fixtures.

Region aggregation uses `dtplyr::lazy_dt(immutable = TRUE)`, computes the original `base::mean` once per group, restores dplyr ordering, and retains the threshold expression and ordered join. This is a targeted conversion of the region-read grouping step. The selected dtplyr implementation took 0.761 seconds versus direct data.table's 0.711 seconds in the same operation-level comparison. Installed versions were R 4.4.2, tidyr 1.3.1, dtplyr 1.3.1 and data.table 1.17.6; no package installation was needed.

The output change saves the same single final QS object earlier and releases raw calls and coverage before VCF/export construction. Per-file traversal and formatted data remain unchanged. Failed tasks must still prevent publication, because the intermediate QS file now appears earlier. The output pilot used four completed historical candidate-optimization burdens inputs and their effective configuration. Its comparison allowed only VCF `fileDate` and PDF `CreationDate`/`ModDate` differences; it used no numeric tolerance or ignored QS metadata. Passing 311 scientific files does not imply identical compressed or PDF bytes.

Level-1 compression applies only to BED output from the R `annotate_coverage_row()` path. Explicit `--coverage-annotation legacy` and automatic fallback retain their existing `bgzip -c` behavior. The tested R-annotation BED grew from 8,038,272,657 to 10,198,489,389 compressed bytes (+26.87%). The full-row timer includes packet loading, the original writer, annotation, compression and indexing. It is not a complete burdens-stage measurement. Compressor-only timings are attribution within that total and must not be added to it.

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

The current 60 chunks per sample divide work but do not cap chunk size as input grows. Extraction and filtering can therefore need more memory on larger inputs. Burdens processing still retains raw calls and germline data, and its final Rle memory depends on run complexity. Earlier object release reduces simultaneous live objects in output generation, but the final result remains one complete QS object. None of these measurements establishes constant memory with increasing input. Native burden-accumulation experiments are not included in these changes or claims.

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

The follow-up was authorized on 2026-10-06. The rebuilt-container Python/R changes and the five additional R optimizations have component evidence. The parser, VCF lookup and targeted dtplyr changes also passed combined extraction/filter chains on one real chunk from each library, with 128 exact QS component comparisons and 36.11% lower summed scientific-worker CPU. Memory results were mixed. Output lifetime and R-annotation compression retain their separate component measurements. Full-workflow scientific comparison and CPU/memory accounting remain pending; the earlier complete-workload acceptance applies to production revision `3911047`.

1. **Python BAM splitter:** implemented, with exact full LIB1 comparisons and
   measured splitting/indexing costs, including the final system-Python/pysam
   installation in the rebuilt container.
2. **R artifact cache:** implemented with common functions in
   `sharedFunctions.R`. Protocol and actual Nextflow cold/warm/concurrent tests
   passed. Paired closed-product cache operations measured its additional CPU
   and memory cost; unchanged scientific preparers were outside that benchmark.
3. **R coverage annotation:** implemented in `sharedFunctions.R`; chr22 and full
   nuclear comparisons passed, including downstream context consumers and
   indexes. CPU and memory costs versus C++ are recorded above.
4. **data.table/dtplyr evaluation:** `fread`/`fwrite` handle coverage I/O. The earlier two pure-R filtering changes retained their measured 5.10% mean complete-worker CPU reduction. Subsequent profiling selected an immutable dtplyr grouping step for region-read filters, preserving `base::mean`, threshold expressions and ordered joins. Combined extraction/filter results are recorded in the [follow-on section](#follow-on-r-optimization-results); no broad dtplyr conversion was made.
5. **Rcpp evaluation:** profiling did not identify a further candidate with
   demonstrated dramatic gains after removing the selected R overhead. No new
   custom compiled helper was added for incremental gains. An Rcpp implementation
   would involve custom C++ called from R; none is selected here.

For all five items, retain the agreed scientific comparison rules and record
performance regressions as well as improvements. Judge aggregate performance by
actual user-plus-system CPU-hours across the pipeline, with significantly lower
peak memory in the currently expensive steps. An individual step may become
more than 5% slower; there is no per-step 5% acceptance limit. Preserve the single
final QS2 and the existing TSV, VCF, PDF and indexed BED outputs.
