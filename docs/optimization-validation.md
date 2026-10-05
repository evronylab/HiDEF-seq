# Optimization validation ledger

This living ledger records the optimization measurements available on 2026-10-05.
Completed component results do **not** establish a whole-pipeline improvement.
Full scientific equivalence, total pipeline CPU, and high-memory-stage RSS
comparisons remain pending.

## Measurement and acceptance

The primary performance measure is actual **user + system CPU time**, summed over
the same executed workload and including children that their parents wait for.
Report wall time and allocated CPU-hours separately: allocated CPU-hours measure
reserved CPUs multiplied by allocation duration, not CPU actually consumed.
Measure tasks on compute nodes; timing a Nextflow driver does not include its
remote jobs. Keep inputs, configuration, thread counts, container and resource
settings comparable, record source revisions, and alternate paired run order.

GNU time and the Python child-accounting fallback report a process/descendant
maximum RSS, not the simultaneous sum of a process tree. This matters for shell
pipelines with concurrent children: that value can understate their combined
memory requirement. Preserve Slurm allocation and step records separately without
double-counting them. Distinguish intrinsic successful-workload performance from
total expenditure including retries, failed experiments, and validation itself.

Acceptance considers lower total actual CPU and lower memory in high-memory
steps. Individual stages may slow down when the full workload improves; there is
no per-step 5% threshold. Maintainability matters alongside measured gains.

Scientific data, schemas, records, formats, filenames, and within-group order must
be preserved. Existing nondeterminism in combined cross-chromgroup row order may
be normalized explicitly. Run metadata and compressed encoding may differ.
Internal floating-point differences at relative error up to `1e-12` require
investigation rather than automatic acceptance; published numeric precision must
match exactly. Configuration comparison checks every original scientific key
and uses an explicit metadata whitelist, never a whole-YAML exclusion.

Reproducible commands and comparator details are in the
[benchmark harness documentation](../scripts/benchmark/README.md). Prepared-cache
identity and current conservative downstream invalidation are documented in
[workflow optimization](workflow-optimization.md).

## Completed: full-reference summary preparation

Job `19190017` compared the original reference scans with the reusable reference
summary in three independent pairs on the full hg38 reference.

| Pair | Baseline actual CPU (s) | Candidate actual CPU (s) | Scientific comparison |
| --- | ---: | ---: | --- |
| 1 | 172.23 | 80.62 | Identical |
| 2 | 175.55 | 80.58 | Identical |
| 3 | 185.78 | 86.81 | Identical |

Reported peak RSS was **16.045 GiB baseline** and **3.190 GiB candidate**.
Every pair matched the complete N-interval GRanges and collapsed 64-element
trinucleotide count vector using R `identical()`. These measurements cover
reference preparation only. Wall time and allocated CPU-hours are not transcribed
in this ledger; retain the raw measurement/accounting records for those metrics.

## Completed: one real LIB1 extraction chunk

The original and candidate extraction scripts were replayed once each on LIB1
chunk 1 with the same configuration and a two-CPU allocation. Complete scientific
QS contents matched exactly; only the explicit `/run_metadata` whitelist applied.

| Metric | Baseline | Candidate |
| --- | ---: | ---: |
| Actual CPU (s) | 971.42 | 702.81 |
| Wall time (s) | 985.35 | 733.63 |
| Allocated CPU-hours during measured command | 0.54741 | 0.40757 |
| Peak RSS (KiB) | 12,205,592 | 9,047,372 |

This single real-chunk comparison used **27.65% less actual CPU** and **25.88%
less peak RSS**, but the arms ran on different nodes. The apparent CPU gain is
**not confirmed** by the completed same-allocation replay below; it must not be
used as the accepted speed estimate. This is not a whole-pipeline result.
Measurements are in workspace `runs/chunk-replay/extract.json` and
`runs/extract-candidate/extract.json`; the exact scientific comparison report is
`runs/extract-candidate/comparison.tsv`. Allocated CPU-hours above cover the timed
command, not its containing job's setup or validation overhead.

## Completed: matched extraction with compressed Rle slicing

Job `19193173` froze both source trees and alternated two pairs of independent
extraction processes within one allocation on `cl014`. Both complete scientific
QS comparisons passed exactly, with only `/run_metadata` ignored.

| Pair | Baseline CPU (s) | Candidate CPU (s) | Baseline peak RSS (KiB) | Candidate peak RSS (KiB) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 658.724 | 733.030 | 12,214,684 | 9,047,936 |
| 2 | 626.403 | 720.260 | 12,216,424 | 9,052,344 |

Mean actual CPU increased **13.09%**, from **642.563 to 726.645 s**. Mean wall
time increased **12.82%**, from **651.427 to 734.960 s**. Mean peak RSS decreased
**25.91%**, from **12,215,554 to 9,050,140 KiB**. These matched-node results
supersede the earlier unpaired CPU estimate. They measure the compressed-tag
implementation that constructs an Rle slice at each query, before the interval
lookup change below. Sources, per-arm metrics and both comparison reports are
retained in workspace `runs/extract-paired/`.

## Completed: real sa query lookup and interval extraction validation

Job `19195014` replayed **182,949 real read/category query groups**, containing
**1,533,808 own-strand positions**, from the original LIB1 chunk-1 extraction.
Reconstructed SBS/MDB, insertion and deletion queries first matched saved
`calls$sa` values exactly; every timed result then matched those same values.

| Repeat | Rle slicing CPU (s) | Interval lookup CPU (s) | Selected-tag dense lookup CPU (s) |
| --- | ---: | ---: | ---: |
| 1 | 44.532 | 2.677 | 2.133 |
| 2 | 43.285 | 2.991 | 2.469 |
| 3 | 43.047 | 2.764 | 5.544 |

Median interval lookup CPU was **2.764 s versus 43.285 s**, a **93.61%** reduction
for this lookup operation. The dense comparator expands only tags with selected
queries; it is a lower-bound lookup comparator, not the full original extraction.
This replay does not measure opposite-strand coordinate mapping or claim a
whole-extraction gain. Full fixtures also preserved types, reversed runs, unusual
indices and fallback errors, including cumulative run endpoints beyond `INT_MAX`.
Results, source/function snapshots, SHA256 hashes and the prior `fad3028`
extraction-source hash are in workspace `runs/sa-lookup/`.

The validated lookup now finds run values from cumulative run endpoints while
retaining direct Rle decoding and the original unusual-index dense fallback.
Job `19195412` passed both fixture suites against the **actual production helper**
and completed two alternating original-versus-interval full extraction pairs
within one allocation on `cl017`. Both complete scientific QS comparisons passed
exactly, with only `/run_metadata` ignored.

| Pair | Baseline CPU (s) | Interval CPU (s) | Baseline peak RSS (KiB) | Interval peak RSS (KiB) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 937.473 | 1,106.634 | 12,210,968 | 9,538,772 |
| 2 | 1,278.315 | 1,265.165 | 12,214,704 | 9,535,276 |

Mean actual CPU increased **7.04%**, from **1,107.894 to 1,185.899 s**; mean wall
time increased **7.00%**, from **1,120.690 to 1,199.123 s**. Mean peak RSS decreased
**21.91%**, from **12,212,836 to 9,537,024 KiB**. The first baseline preceded a
colocated filter benchmark, and baseline CPU differed substantially between the
two pairs; retain both results rather than treating either pair as a stable speed
estimate. The interval lookup microbenchmark does not establish a complete
extraction CPU improvement. Sources, metrics and exact comparison reports are in
workspace `runs/extract-paired-interval/`. The indexed lookup below addresses
the remaining measured extraction cost.

## Completed: indexed indel tag lookup, two whole-extraction pairs

Diagnostic profiling localized about 48% of sampled original extraction time to
six indel tag assignments. Each read group repeatedly searched named tag lists.
The candidate resolves read names once with `match()` and uses integer indices
for the existing sa/sm/sx assignments on both strands. It removes the temporary
index before aggregation. The three tag lists originate from the same BAM rows
with identical names and ordering. First duplicate-name matches and missing,
empty or NA name behavior remain the same.

Job `19200377` alternated two original/candidate pairs in separate R processes
within one allocation on `cs612`, replaying real LIB1 chunk 1. The candidate also
includes the direct Rle decoding and interval lookup described above. Both
complete scientific QS comparisons passed exactly, with only `/run_metadata`
ignored; all workers finished before comparison began.

| Pair | Original CPU (s) | Candidate CPU (s) | Original peak RSS (KiB) | Candidate peak RSS (KiB) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 731.129 | 391.982 | 12,214,900 | 9,529,492 |
| 2 | 707.680 | 394.897 | 12,213,496 | 9,535,212 |

Mean actual CPU decreased **45.31%**, from **719.405 to 393.440 s**; mean wall
time decreased **45.74%**, from **737.400 to 400.116 s**. Mean peak process RSS
decreased **21.96%**, from **12,214,198 to 9,532,352 KiB**. These measurements
include startup, extraction and QS serialization; they do not establish an
all-chunk or whole-pipeline improvement. Source snapshots, input/configuration
identities, per-arm metrics, exact comparisons and aggregate metrics are retained
in workspace `runs/extract-paired-index/`.

The actual assignment fixture passed successful homogeneous integer/double tags,
duplicate and absent names, empty queries, out-of-bounds positions, keys and row
order. Mixed tag types that trigger a legacy truncation error retain that error.
The additional integrated SA/extraction fixture remains queued at this entry;
production promotion awaits that check.

## Completed: germline VCF annotation block, three pairs

The full VCF load/conversion, quality filter, deduplication, summary, and call
annotation block was replayed on LIB1 chunk 1, nuclear chromosome group `1-22X`,
filtergroup `lenient`. Both arms used the same actual calls from the candidate
script's preceding filters and loaded the same germline artifact afresh. The
candidate restricts only annotation summaries to matching call keys; its extra
matching scan is included in these timings.

| Pair | Baseline CPU (s) | Candidate CPU (s) | Baseline wall (s) | Candidate wall (s) | Exact comparison |
| --- | ---: | ---: | ---: | ---: | --- |
| 1 | 275.089 | 13.674 | 278.905 | 13.755 | Both tables identical |
| 2 | 264.022 | 13.947 | 265.061 | 14.009 | Both tables identical |
| 3 | 255.846 | 15.116 | 256.769 | 15.179 | Both tables identical |

Each pair matched all **1,023,018 annotated calls** and the entire
**10,772,319-row filtered germline table** using `identical()`. The latter remains
unchanged for downstream per-file filters and whole-genome statistics. Mean
actual CPU was **264.986 s baseline** and **14.246 s candidate**, a **94.62%**
reduction for this block only. Full-filter and whole-pipeline gains remain pending.

Observed phase peak RSS was 11,746,588 / 11,394,356 / 10,184,904 KiB for baseline
and 8,903,732 / 8,945,940 / 8,967,980 KiB for candidate. The Linux high-water reset
succeeded for each arm. However, arms share one R process and retain a common
preceding-filter context; allocator-resident pages can persist even after each
arm's outputs are released and garbage-collected. These are not independent
fresh-process memory estimates. Alternating arm order limits systematic warm-up
effects; full-filter process replay is required for a standalone memory claim.
Serialization and exact comparison were outside block timing.

Actual ASTs, input/source fingerprints, session details, complete output tables,
`metrics.tsv`, and `comparisons.tsv` are retained in workspace
`runs/germline-vcf-block/`. The sourceable AST locator and benchmark commands are
documented in the [benchmark harness](../scripts/benchmark/README.md).

## Completed: one real nuclear lenient burden chunk

Job `19194158` replayed complete `calculateBurdens.R` baseline and candidate
processes within one allocation using the same original LIB1 chunk-1 nuclear
lenient filter output. Both arms included the unchanged full coverage BED writer.

| Metric | Baseline | Candidate |
| --- | ---: | ---: |
| Actual CPU (s) | 1,870.558 | 1,924.817 |
| Wall time (s) | 1,730.558 | 1,761.718 |
| Peak RSS (KiB) | 23,739,512 | 21,305,220 |

The complete scientific QS comparison passed exactly, with only the explicit
`/run_metadata` whitelist. The coverage BED gzip and its tabix index were both
byte-identical. Peak RSS decreased **10.25%**, while actual CPU increased
**2.90%**; this one-chunk replay does not establish a burden CPU improvement.
It also does not measure multi-chunk accumulator scaling. Frozen sources,
configuration, complete per-process metrics and comparison reports are retained
in workspace `runs/burdens-one-chunk/`.

## Completed: germline formatting, three fresh-process pairs

Job `19194048` alternated three pairs of independent R workers on `cl012`, each
loading the same benchmark packet containing only raw germline data and config.
The actual formatter ASTs came from original/candidate `outputResults.R`.
All three formatted **5,374,705-row × 56-column** tables matched with R
`identical()`, including values, classes, attributes and row/column order.

| Pair | Baseline formatting CPU (s) | Candidate formatting CPU (s) | Baseline process peak RSS (KiB) | Candidate process peak RSS (KiB) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 1,974.829 | 773.302 | 22,654,384 | 22,729,472 |
| 2 | 1,387.397 | 711.664 | 22,653,780 | 22,728,184 |
| 3 | 1,392.369 | 714.989 | 22,654,248 | 22,730,140 |

Mean formatting CPU decreased **53.73%**, from **1,584.865 to 733.318 s**.
Complete-worker CPU, including packet loading and formatted-table saving,
decreased **51.42%**, from **1,645.095 to 799.124 s**. Mean worker wall time was
**1,652.912 versus 802.874 s**. The first baseline was slower than the other two;
retain every pair rather than reporting only the largest gain.

Mean process peak RSS was **22,654,137 versus 22,729,265 KiB**, a **0.33%
increase**; this optimization has no demonstrated formatter memory benefit.
The profile resets HWM immediately before formatting and reconstructs complete
worker peak from pre-reset and remaining peaks, including serialization. These
fresh processes avoid the allocator-retention confound of the earlier shared-R
experiment, but they do not retain coverage and are not full-output-stage memory
measurements. Sources, phase/process metrics and all exact comparisons are in
workspace `runs/germline-formatting-fresh/`.

The earlier shared-process job `19190488` was intentionally canceled after one
exact pair when this independent-process replacement was adopted. Its partial
results are preserved; it was not a scientific failure and its memory values
are not used as accepted independent-worker estimates.

## Completed: existing final LIB1 QS2 profile

Job `19188581` profiled one existing final QS2. Its compressed size was **2.397
GiB** and its logical R object size was **31.751 GiB**. Principal logical components
were coverage (**20.049 GiB**), raw germline data (**9.456 GiB**), and formatted
germline data (**1.798 GiB**). Logical `object.size()` values can double-count
shared storage and are distinct from compressed size and process RSS.

| Phase | Actual CPU (s) | Wall time (s) | Peak RSS (GiB) |
| --- | ---: | ---: | ---: |
| Load | 137.56 | 138.81 | 31.8 |
| Save | 86.00 | 87.08 | approximately 32.2 |

Saving this loaded object did not reproduce the much larger historical
`outputResults` task peak. The measured final object does not explain that entire
task's memory use; assembly and export temporaries remain profiling targets.
The selected design retains **one final QS2 per sample**, along with existing
TSV, VCF, PDF, BED.gz, and Tabix products. No call batching or output sharding is
introduced, and no scientific schema change is authorized by these measurements.

## Completed: BAM dispatch and Nextflow integration

Job `19191730` dispatched the full LIB1 analysis BAM into 60 chunks. All chunks
passed BAM readability checks and received PacBio and samtools indices. Chunk 1
matched legacy ordered SAM records and tags exactly; the remaining 59 chunks
have not each been compared by independently rerunning the legacy splitter.

| Phase | Actual CPU (s) | Wall (s) | Process/descendant peak RSS (KiB) |
| --- | ---: | ---: | ---: |
| All 60 dispatch outputs | 2,161.586 | 1,092.321 | 526,472 |
| All 60 PBI/BAI index sets | 1,499.199 | 764.531 | 63,596 |

Dispatch plus indexing consumed **1.017 actual CPU-hours**. The original measured
one-chunk split consumed 1,666.483 CPU seconds before its separate indexing.
That single-chunk measurement is not an independently measured 60-chunk total.
Slurm's sampled memory accounting was substantially larger than the dispatcher's
private/process RSS; retain both rather than interpreting file-cache accounting
as the helper's live heap. The shared compression pool has two workers, but
HTSlib also creates background stream threads (64 OS threads were observed).

The fresh full baseline subsequently completed all **120 original split tasks**
for both samples. Their Slurm allocation records sum to **238,950.663 actual CPU
seconds (66.375 CPU-hours)**, including indexing within each task. The two
`countAnalysisZMWs` tasks add 30.031 CPU seconds. Native job IDs are deduplicated
across trace files, and allocation/step CPU is not double-counted. The frozen
trace and accounting snapshot is workspace
`runtime/accounting-checkpoints/20261005T051136Z/baseline-6b6598b/`.
This provides a measured two-sample baseline cost; comparison with the full
two-sample candidate and all-chunk scientific validation remains pending.

Follow-up job `19198407` strengthened the real chunk-1 check: all **173,355 ordered
decompressed BAM record blocks** and the binary reference dictionary matched
exactly. Both BAI/PBI payloads matched independent rebuilds against their own BAM.
Its raw report requested review solely for program-header provenance: the legacy
`zmwfilter` record is absent, and fields within one `samtools` program record are
reordered. A separate value-checked audit accepts these allowed metadata changes;
all remaining ordered program records have identical field values except command
provenance. The unmodified report, exact headers and audit are retained in
`runs/dispatch-bam-raw/`. This does not relax checks on non-program headers.

Nextflow fixture job `19194083` passed 12-chunk, single-chunk and empty-sample
cases, exact legacy comparisons for all 13 emitted BAMs, index and tuple checks,
and publication checks. Published files were hard links to retained work files.
All six tasks were reused on `-resume`. Preview job `19194084` also verified all
57 original configuration keys exactly. The integrated checkpoint is
`fad3028563f70e37dcb958d7d27e87ca4e78e3a6`; full remote candidate job `19194442`
uses `-r optimization -latest` and records its resolved revision separately.

## Completed: chromosome 22 coverage annotation experiment

Job `19192378` alternated two complete legacy/candidate operation pairs for all
eight strict coverage rows on actual LIB1 chromosome 22. Both used the unchanged
original R BED writer and `chunk_runs=1e7`. The candidate substitutes direct FASTA
annotation for the large reference-BED intersection; production still uses the
legacy annotation pending the full nuclear benchmark.

| Pair | Baseline process CPU (s) | Candidate process CPU (s) | Baseline wall (s) | Candidate wall (s) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 782.175 | 134.502 | 688.412 | 112.771 |
| 2 | 895.969 | 134.966 | 786.281 | 113.348 |

Every decompressed BED byte and numeric context count matched in both pairs;
Tabix contig lists and boundary queries also matched. Mean process CPU was
83.94% lower. Process RSS was approximately 0.70 GiB in both arms, so this
experiment establishes no annotation-memory reduction. Separate whole-final-QS
preparation peaked at 32.03 GiB and is excluded from these operation timings.
The process metrics include packet loading, writing, annotation, compression and
indexing; operation-only phase metrics and sampled concurrent process-group RSS
are also retained under workspace `runs/coverage-annotation-chr22/`. Full nuclear
job `19194115` will assess scaling before production adoption.

## Rejected: coordinate-sort reuse

The completed tiny fixture matched ordered SAM records for one input but failed
for every tested multi-input order (`AB`, `BA`, `ABC`, `CBA`) because coordinate
ties changed record order. Record multisets, headers excluding program records,
and index checks passed. These do not establish the required ordered equivalence.
Keep the current coordinate sorts and `pbmerge`; no single-input exception or
sort-reuse implementation is adopted. Results are retained in workspace
`runs/coordinate-merge-fixtures/results.json`.

## Rejected: mitochondrial shared-session filtering

Job `19191880` completed one paired mitochondrial replay of lenient and strict
filtergroups. Both complete QS outputs matched exactly between separate R
processes and the shared-session prototype. Direct inspection confirms this
mitochondrial source snapshot already included the germline VCF annotation
`semi_join` optimization; the earlier nuclear experiment used a different snapshot.

| Complete two-group chain | Actual CPU (s) | Wall (s) | Allocated CPU-hours during command | Peak RSS (KiB) |
| --- | ---: | ---: | ---: | ---: |
| Separate processes | 479.830 | 484.334 | 0.26907 | 7,063,536 |
| Shared session | 448.831 | 452.033 | 0.25113 | 9,918,568 |

The single pair saved 6.46% actual CPU but increased peak RSS by 40.42%.
Mitochondrial fusion is rejected; no production implementation is adopted.
Records and strict comparison reports are in workspace `runs/filter-chain-mito/`.
The separate nuclear experiment is described below.

## Rejected: nuclear shared-session filtering

Job `19191877` completed one paired replay of the full nuclear lenient→strict
chain on the original LIB1 extraction chunk. Both full QS outputs matched
exactly, with **no metadata or floating-point exceptions** in these comparisons.
The shared session read the extraction once and recorded one cache hit for the
second group.

| Complete two-group chain | Actual CPU (s) | Wall (s) | Peak RSS (KiB) |
| --- | ---: | ---: | ---: |
| Separate processes | 4,052.435 | 4,106.104 | 16,177,948 |
| Shared session | 4,005.754 | 4,002.416 | 16,464,072 |

Actual CPU decreased only **1.15%**, while peak RSS increased **1.77%**. This
single pair does not justify the additional shared-session state and workflow
complexity; nuclear fusion is rejected and production remains unfused.
The frozen `filterCalls.R` hash
`5dd7b3c00ddbb6f9261f3409b0cd5580cd0194e1bdf47b53d0f1f70a4618e162`
**predates the germline annotation `semi_join` optimization**. These results
therefore describe that frozen prototype, not an estimate of fusion on the
current filtering source. Sources, complete measurements, cache-use counts and
both exact comparisons are retained in workspace `runs/filter-chain-1-22X/`.

## Completed: all four complete chunk-filtering comparisons

Job `19194137` replayed all four chromosome/filtergroup combinations for the
original LIB1 chunk-1 extraction with 24 GiB and two CPUs. Both it and original
replay job `19188367` completed successfully. All four complete scientific QS
comparisons passed exactly.
Workspace `runs/filter-candidate-current/` contains the frozen `code/bin` and
benchmark scripts, per-file SHA256 manifest, prepared config, logs and metrics.
The extraction argument remains the exact relative `extractCalls.chunk1.qs2`
used by the original replay; its symlink resolves to the original extraction.

Preflight required every original parsed configuration value to match and exactly
two added preparation fields. Complete original outputs were compared with only
`/config/yaml.config/reference_summary_file` and
`/config/yaml.config/germline_coverage_filters` whitelisted. The full config,
inherited run metadata, scientific schemas and within-group row order otherwise
remained strict.

| Chromgroup / filtergroup | Baseline CPU (s) | Candidate CPU (s) | Baseline peak RSS (KiB) | Candidate peak RSS (KiB) |
| --- | ---: | ---: | ---: | ---: |
| Nuclear / lenient | 3,215.248 | 1,897.161 | 25,199,392 | 14,156,524 |
| Nuclear / strict | 2,597.530 | 1,594.728 | 25,121,932 | 13,037,608 |
| Mitochondrial / lenient | 843.808 | 171.237 | 21,726,896 | 7,061,936 |
| Mitochondrial / strict | 796.762 | 172.205 | 21,725,876 | 7,064,580 |

These replays used **different nodes**, so the CPU effects remain preliminary.
Nuclear CPU fell 40.99% / 38.61%, with peak RSS 43.82% / 48.10% lower;
mitochondrial CPU fell 79.71% / 78.39%, with peak RSS about 67.5% lower.
Same-allocation job `19196053` is running two alternating nuclear-lenient pairs
to assess the complete-filter CPU effect under better matched conditions.
The measured operations include script startup, filtering and QS serialization;
new shared-cache preparation and scientific comparison are outside these timings.
This is not a whole-pipeline or preparation-inclusive result. Wall times and
per-group results are retained in `comparisons/four-group-summary.json`.

## Pending before overall conclusions

- Complete comparable baseline and candidate runs for both samples; compare
  total actual CPU, wall time, allocated CPU-hours, and high-memory-stage RSS.
- Compare complete nested scientific QS objects and every published scientific
  output, resolving all discrepancies and binary-file review findings.
- Complete paired measurements and exact-output checks for the remaining
  implementations. Full nuclear BED annotation remains an experiment until
  separately accepted. BAM dispatch passed its component and workflow gates;
  complete pipeline validation remains pending. Coordinate-sort reuse and both
  mitochondrial and nuclear shared-session filtering were rejected above.
- Validate cache/resume behavior, effective YAML parsing and scientific keys,
  including missing prepared artifacts and concurrent launches.

The baseline source revision is
`6b6598b236d1f4a96e9e34597f7c7fa3c2cb2a3a`. Candidate revisions must be recorded with
each subsequent measurement. This initial ledger transcribes the completed
reference and QS2 measurements from the workspace `WORK_LOG.md`; isolated runs
and accounting artifacts are retained under the workspace `runs/` directory.

## Diagnostic: unchanged germline VCF export with resident final payload

Job `19197272` profiled the actual normalizer and VCF writer from the pinned
sources on **2,992,350 nuclear lenient germline records**, retaining the complete
original final QS object, including raw germline and coverage. This also retains
later-group formatted/final results absent at the first real export, so it is a
conservative payload approximation rather than a replay of the entire stage.

The installed VariantAnnotation **1.52.0** writer defaults to two half-table
chunks above one million rows. Expression markers localized large temporary
allocations to nonlogical INFO strings, the INFO matrix row collapse, and the
fixed output-line strings. RSS was **32.67 GiB after loading**, **33.78 GiB at
export entry**, and reached **59.09 GiB process peak**. This reproduces substantial
export growth without preceding formatter allocator residue. The same job's
Slurm batch MaxRSS was **69,231,192 KiB (66.02 GiB)**, which is a different
whole-job accounting measure from the R process maximum. It must not be
substituted into the paired per-worker RSS comparisons or compared directly
with historical full-stage peaks from another metric source. The residual
historical peak difference is not assigned wholly to formatter retention.
Comparable complete-stage process and Slurm memory measurements remain required
before lowering production memory requests.

The diagnostic normalization/export region took **248.469 actual CPU seconds**
and **252.735 wall seconds**; the complete worker, including loading, took
**417.719 CPU seconds** and **425.063 wall seconds**. These instrumented timings
are diagnostic, not an accepted optimized-versus-original performance comparison.
There was no explicit GC, HWM reset, changed writer buffer, or production edit.

Every decompressed scientific VCF line matched the published original, ignoring
only `fileDate`. Tabix chromosome names and first/last-record-position queries
matched on all **23 chromosomes**. Frozen source hashes, installed writer methods,
event metrics, payload disclosure and comparisons are in workspace
`runs/germline-vcf-profile/`.

## Completed: native VCF writer buffer, two fresh-process pairs

Job `19198553` alternated independent default→native and native→default workers
on `cs648`, retaining the same complete final QS payload. Only the existing
`writeVcf(..., nchunk=100000L)` argument changed; timed workers had no expression
instrumentation. No QS rewriting, scientific output shards, or pipeline call
batching was introduced.

| Pair | Default whole-worker CPU (s) | Native whole-worker CPU (s) | Default process peak RSS (KiB) | Native process peak RSS (KiB) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 405.871 | 374.196 | 61,973,644 | 51,508,468 |
| 2 | 405.012 | 374.630 | 61,973,256 | 51,509,896 |

Mean complete-worker CPU decreased **7.65%**, from **405.442 to 374.413 s**.
Mean normalization/export CPU decreased **13.17%**, from **253.595 to 220.192 s**;
input loading was outside that phase. Mean whole-worker wall time decreased
**9.04%**, from **414.415 to 376.969 s**. Mean process peak RSS decreased
**16.89%**, from **61,973,450 to 51,509,182 KiB**. These are focused export-worker
measurements, not measured whole-output-stage or whole-pipeline improvements.
Production memory requests remain unchanged pending comparable full-stage data.

Both pairs matched every decompressed scientific VCF line against the published
original and each other, ignoring only `fileDate`. All **2,992,350 records**,
index chromosome lists and first/last-record-position queries on all **23
chromosomes** matched. Every original and generated index also matched a pinned
Rsamtools rebuild against its own BGZF encoding. Pinned fixtures covered mixed
numeric magnitudes, missing values, flags, factors, anchored indels, duplicate
loci, empty tables and indexed queries.

The one-argument change is included in the production writer. The harness now
explicitly removes only the known configured argument for its default arm and
sets it for its native arm; it rejects unexpected settings and records effective
writer ASTs. Job `19201862` passed those guards and fixtures against the
actual configured production source (28 seconds, exit 0). Its saved default-arm
AST removes the configured argument, and all eight original BGZF/index fixture
hash/stat fingerprints remained unchanged. The configured-source preflight is in
`runs/germline-vcf-production-preflight/`; paired sources, measurements,
comparisons and original-fixture immutability proof are in
`runs/germline-vcf-native-v2/`.

The earlier attempt `19198210` was intentionally canceled during its first worker,
before comparisons, after review found that Rsamtools path normalization made
symlink staging unsafe for index rebuilding. Its partial artifacts remain in
`runs/germline-vcf-native/`. The corrected run stages physical scratch copies;
published inputs were not modified. Whole-job Slurm/cgroup accounting remains
separate from the per-worker process RSS figures above.
