# Optimization validation ledger

This living ledger records the optimization measurements available on 2026-10-06.
Completed component results do **not** establish a whole-pipeline improvement.
Full scientific equivalence, total pipeline CPU, and the remaining high-memory
stage comparisons remain pending. Both complete LIB1 nuclear burden groups have
passed strict scientific comparison and show substantial CPU/RSS reductions.
The full candidate run has completed; the original run and final publication
comparison are still in progress.

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

## Verified: cache identities across shared-filesystem mounts

The same immutable filtered input was recorded with device number 59 in
`runs/burdens-one-chunk/manifest.json` and 57 in
`runs/burdens-annotation-integrated/manifest.json`; inode, path, size and nanosecond
timestamps matched. Java's Unix file key includes this node-local device number,
which could unnecessarily invalidate prepared caches after moving the controller.
The identity now records the Unix inode separately and excludes the device.
It remains a metadata identity, with the limitations documented above.

Job `19207641` exercised the actual production Nextflow function: repeated reads
and symlink aliases matched; mtime changes, size changes, and a same-path,
same-size, same-mtime file replacement invalidated the identity. Three real input
probes also matched their signed Unix inode and size. Full-workflow preview
`19207702` completed and emitted 19 expected prepared-artifact identities.
Evidence is in `runs/prepared-input-identity/` and `runs/cache-identity-preview/`.
Migration of existing development-run cache entries requires separate checksum
verification and exact comparison of every other identity field before reuse.
The development-run migration passed 18 focused fixtures, then job `19208174`
verified and atomically published 17 completed bundles under the new identities.
All product checksums and the original manifests remained exact. The only allowed
identity changes were removal of the device number and explicit dependency-key
updates to other proven bundles. This is a one-time development migration, not a
production compatibility path. The two unfinished germline bundles were excluded
from that first migration. After they completed, job `19213931` verified and
published those two bundles (36,287,915,378 product bytes) in 132.495 actual CPU
seconds and 221.768 wall seconds. All 19 old manifests and the 17 prior migrated
targets remained unchanged; all 19 new prepared entries are now verified.
Evidence is in `runs/prepared-device-migration-v2/` and
`runs/prepared-device-migration-additional/`.

The old full candidate was checkpointed only after all valuable frontend tasks
and raw germline preparation completed. Its remaining derived-coverage worker
finished successfully after controller shutdown; its 362.168 CPU seconds are
retained in a separate accounting-only supplemental trace. No stale R analysis
was launched. Frozen source, logs, native-job accounting and completion proofs
are in `runtime/checkpoint-evidence/final-drained-20261005T0951Z/`. Historical
effective-YAML audit `19214158` confirmed all 57 original fields, four exact
helper fields and 21 artifact paths against the 19 retained old bundles.
Actual-run audit `19217206` then verified all 19 current product inventories and
bound the resumed YAML (`8319091851aa5c61ed331999973893d8096c1ff0bbe149830ecc519ebc689f1c`)
to all 57 unchanged source values and four exact helper fields/21 product paths.
The actual task identities matched the reviewed inventory, and all 22 completed
read-processing tasks were cache hits. Evidence is in
`runs/resumed-effective-config-audit/`.

All 22 prepared-cache comparison products are now scientifically validated.
The 16 earlier comparisons are joined to the six remaining products by their
verified immutable manifests and identical product checksums in
`runs/combined-cache-validation-evidence.json`. The raw germline VCF comparison
checked all **391,979,965,111 decompressed body bytes** exactly after its narrowly
audited header provenance differences. The reference BED also matched every
decompressed byte; the full germline BigWig and derived coverage QS matched.
The initial reference-index comparison incorrectly inferred `tabix -p bed`;
both actual pipeline versions use `-s1 -b2 -e3`. Corrected rebuilds in `19218720`
matched each original index exactly against its own compressed file. Existing
pipeline indexing semantics were preserved. The initial review remains recorded
in `19217283`; its other five products passed. This establishes preparation
equivalence, not final pipeline-publication equivalence.

## Completed stages: both-sample extraction and filtering

All 120 candidate extractions and 480 filters completed successfully before
the singleton-reference correction. That correction changes only burden code,
so the completed upstream work is reusable.

| Stage | Original actual CPU-hours | Candidate actual CPU-hours | CPU decrease | Original maximum Slurm RSS (KiB) | Candidate maximum Slurm RSS (KiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Extraction, 120 tasks | 29.314 | 17.394 | 40.66% | 15,020,132 | 13,359,224 |
| Filtering, 480 tasks | 257.643 | 161.463 | 37.33% | 38,828,464 | 29,138,900 |

These complete-stage totals have no missing or unresolved accounting. They span
different nodes and are distinct from the controlled paired process-RSS trials.
Allocated CPU-hours were 30.228 versus 18.568 for extraction and 263.865 versus
164.018 for filtering. Frozen accounting is in
`runtime/accounting-checkpoints/20261005T115110Z/` and
`runtime/accounting-checkpoints/20261005T135308Z/`, with the original completed
stage snapshot in `20261005T082229Z/`.

The candidate controller finished with failure after the discovered burden bug;
all healthy filtering tasks were allowed to finish and remain cached. Failed and
deliberately stopped development burden attempts are retained separately, with
signal-interrupted CPU marked as a lower bound. Full successful-workload CPU,
the remaining nuclear burden peaks, and final publication equivalence remain pending.

## Completed: first full 60-chunk nuclear burden comparison

The LIB1 nuclear lenient group completed successfully in baseline job `19210736`
and candidate job `19223739` (production revision `3911047`). Both use all 60
ordered filter outputs and the same scientific settings and two-CPU allocation.

| Metric | Baseline | Candidate | Decrease |
| --- | ---: | ---: | ---: |
| Actual user + system CPU-hours | 9.4044 | 3.9775 | 57.71% |
| Slurm-reported maximum RSS (KiB) | 178,074,992 | 64,809,888 | 63.61% |
| Slurm-reported maximum RSS (GiB) | 169.826 | 61.808 | 63.61% |
| Allocation elapsed time | 8:46:47 | 3:49:59 | 56.34% |

These are complete task measurements on different nodes (`cl011` and `cl015`),
not a repeated matched-node trial or a whole-pipeline result. Slurm RSS is distinct
from R heap size and Nextflow's sampled peak: the latter reported 132,119,540 KiB
for this baseline task. Frozen allocation and step accounting is in
`runtime/accounting-checkpoints/lenient-nuclear-first-pair-20261005/`.

Job `19234598` compared both complete burden QS objects, including real nuclear
sensitivity results and all accumulated coverage. Scientific values, types,
attributes and row order matched with **zero failures and zero numeric reviews**.
No cross-chromgroup row normalization was enabled. The nine explicit
configuration/run-metadata exclusions were accompanied by saved-YAML binding.
Table provenance used only the reviewed exact region-cache and germline-source
path mappings; no other table values were normalized.
Input file identities remained unchanged throughout comparison. Frozen source,
commands and reports are in `runs/full-nuclear-burden-lib1-lenient/`.

This establishes the complete lenient nuclear burden object for one sample.
Coverage BED/index publications, combined final QS/TSV/VCF/PDF outputs and
whole-pipeline totals still require their pending comparisons.

## Completed: full LIB1 strict nuclear burden comparison

Baseline job `19210451` and candidate `19223765` both completed the full 60-chunk
strict nuclear calculation successfully. Comparison job `19261320` then matched
the complete burden QS objects with **zero failures and zero numeric reviews**,
strict row order, the same nine reviewed metadata exclusions and saved-YAML/path
bindings used in the lenient test. Source files remained unchanged. Evidence is
in `runs/full-nuclear-burden-lib1-strict/`.

| Metric | Baseline | Candidate | Decrease |
| --- | ---: | ---: | ---: |
| Actual user + system CPU-hours | 20.9361 | 9.3028 | 55.57% |
| Slurm-reported maximum RSS (KiB) | 301,981,536 | 125,981,960 | 58.28% |
| Slurm-reported maximum RSS (GiB) | 287.992 | 120.146 | 58.28% |
| Nextflow sampled peak RSS (KiB) | 241,732,676 | 51,527,288 | 78.68% |
| Allocation elapsed time | 18:46:11 | 8:47:42 | 53.14% |

Both tasks requested two CPUs and 288 GiB, on `cl017` and `cl015` respectively.
The baseline's Slurm memory measurement reached its allocation limit; it finished
without an OOM error. Slurm accounting and Nextflow's sampled process metrics
measure different aspects of memory and must not be interpreted as R heap size
or used interchangeably. This is one completed pair on different nodes.
Allocation/step evidence is in
`runtime/accounting-checkpoints/lib1-strict-nuclear-20261005/`.

## Full candidate completed; final comparison pending

Candidate controller `19223438` completed production revision `3911047` with
29 newly completed and 623 cached tasks, zero failed tasks and zero retries.
All 652 original or current native task IDs have resolved successful accounting:
**230.2022 actual CPU-hours**, including **51.780 CPU seconds** for the final 19
cache verification/restoration tasks. Cached tasks retain their original CPU
cost; this is not the cost of only the final resumed invocation. Allocated CPU
time was 261.1728 hours. The final controller itself used another 65.831 CPU
seconds. Evidence is in
`runtime/accounting-checkpoints/candidate-complete-3911047/`.

The candidate's first successful construction of all 19 prepared bundles cost
**8.9972 additional CPU-hours**. The separately reviewed preparation plan includes
the derived-coverage worker that finished after its original controller stopped:
`runtime/accounting-preparation-review/cold-preparation-plan.json`. Thus cold
construction plus the final successful worker workload is 239.1994 CPU-hours.
This still excludes failed/canceled development, earlier restoration, migration,
other controller invocations and validation expenditure. The original cold cache
construction was not measured; final reporting must compare matched downstream
work separately. Charging the candidate's cold construction while charging no
original cold preparation can establish a conservative worker-CPU bound, not a
matched cold-versus-cold speedup.

Separate accounting reconciles all 703 distinct worker submissions across the
three candidate invocations, including the orphan preparation task. Their
observed cost is **239.8723 CPU-hours**, of which **0.6729 hours** is additional
work beyond the final successful workload plus first cold preparation. That
additional work comprises earlier restores, two completed burden jobs absent
from the final trace, two failed burden attempts and nine deliberately canceled
attempts. All three controllers add 172.523 CPU seconds; the two successful
cache-identity migration allocations add 211.053 seconds. Canceled workers and
one controller have incomplete child accounting, so totals containing them are
lower bounds. Validation, benchmarks and other development activity are excluded;
these figures are not a complete development-cost total. Evidence and disjoint
job categories are in `runtime/accounting-development-final/`.

The candidate has exactly the 711 expected scientific publication paths. A
provenance audit binds every candidate publication to its actual completed or
cached producer. Generic genome-spectrum filenames are shared across samples,
so the audit now also binds each output task's sample argument and output prefix
to its publication directory; inode or copied-byte checks remain unchanged.
The correction passed 38 focused fixtures and the complete candidate audit.
This proves inventory and producer identity, not scientific contents. The full
scientific publication/index comparison remains queued behind the original run.

## Corrected and verified: singleton indel reference and four mitochondrial burdens

The full candidate exposed a singleton-reference bug in strict mitochondrial
burden job `19219429`: BSgenome `getSeq()` returns a `DNAString` when exactly one
chromosome is requested, while indel context lookup requires a named
`DNAStringSet`. The targeted-reference optimization now explicitly constructs
that set and restores the requested chromosome names. This also covers nuclear
groups whose indel calls occur on only one chromosome; multiple-chromosome and
empty-reference behavior remain unchanged.

Fixture job `19221726` passed the actual installed BSgenome API checks for chrM,
chrY, reordered multiple chromosomes and empty input. It also passed the six
barcode scenarios, coverage accumulator and invariant checks, and exact
full/subset indel spectra with terminal mitochondrial and nuclear contexts. The
earlier miniature mock incorrectly always returned a set; it now reproduces
BSgenome's singleton simplification. Failed fixture records are retained:
`19221187` lacked a frozen fixture dependency, and `19221461` compared internal
sequence backing pools directly. The corrected actual-reference assertion
compares exact ordered bases, names, widths and class. Neither failure ran a
scientific replay.

Job `19221727` replayed both complete strict mitochondrial burdens using all 60
originally ordered candidate filter outputs, the unchanged effective YAML and
the compiled annotation helper. Both workers completed successfully. Their QS
objects and the two previously completed lenient mitochondrial objects were
compared with the corresponding original pipeline outputs. Initial comparisons
identified only the relocated `region_filter_threshold_file` provenance column;
the driver therefore returned failure under its original policy. The input
identity checks passed, and all initial reports remain preserved.

Independent audit `19222187` checked all four observed columns: each has 13 rows,
six unchanged missing values and seven relocated paths, with unchanged column
attributes and order. Complete paths were checked against the 13 scientifically
verified region-filter product pairs and the actual effective configuration.
Mutations covering wrong products, foreign paths, missing values, order, types
and attributes were rejected. Job `19222615` then passed all four comparisons
through the actual shared comparator and its current provenance fixtures:
**zero failures and zero numeric reviews** per object. Each comparison reports
nine reviewed configuration/run-metadata exclusions; saved configurations are
independently bound to their exact reviewed YAML. The cache-path rule checks
complete values in that one named column, and the existing germline provenance
rules remain unchanged.

Evidence is in `runs/mito-singleton-fix-v3/`, `runs/mito-cache-path-audit/` and
`runs/mito-reviewed-path-comparison/`; earlier failure snapshots remain in the
preceding directories. This is real 60-chunk mitochondrial downstream proof.
These groups skip the configured nuclear sensitivity calculation. The separate
complete LIB1 lenient nuclear check above covers that calculation; the other
nuclear groups, combined final publication and whole-pipeline performance
comparisons remain pending.

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
The additional integrated SA/extraction and actual indel-assignment fixtures
passed in job `19206115` against the same candidate source. The indexed lookup is
now included in the production script. The original queued fixture `19202242`
was cancelled before execution after its identical replacement was accepted on
Torch's short-job partition; no timed scientific worker was rerun for that change.

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
reduction for this block only. The complete filter comparisons below measure its
combined effect with the other filtering changes; whole-pipeline gains remain
pending.

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
The two-sample candidate has now completed both dispatches and their indices,
both enumerations and the shared compiler. That complete dispatch scope used
**6,959.049 actual CPU seconds (1.933 CPU-hours)** versus **238,980.694 seconds
(66.384 CPU-hours)** for the original splits plus enumerations: **97.09% less
CPU**. The compiler is included in the candidate total. Allocated CPU-hours were
1.993 versus 71.294; they are separate from actual CPU use. These are completed
stage totals from the full runs, not paired node-controlled trials or complete
pipeline totals. The updated frozen traces and accounting are in
`runtime/accounting-checkpoints/20261005T082229Z/`; neither dispatch scope has
unresolved accounting. Full downstream scientific validation remains pending.

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
annotation for the large reference-BED intersection; the legacy arm retained
the original annotation. The full nuclear gate and subsequent integration are
documented below.

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
are also retained under workspace `runs/coverage-annotation-chr22/`. This
chromosome-specific percentage is not projected to the full nuclear operation.

## Completed and integrated: full nuclear coverage annotation

Job `19194115` completed one legacy/candidate pair for all eight strict nuclear
coverage rows from the original LIB1 final QS. Every decompressed BED byte and
all eight numeric context-count maps matched exactly. Tabix contig lists and
first/last record boundary queries matched for every BED. Per-BED decompressed
SHA256 hashes, frozen sources/helper, input identity and measurements are in
workspace `runs/coverage-annotation-nuclear/`. The job completed with exit 0.
Both arms retained the original R BED writer and `chunk_runs=1e7`.

| Metric | Legacy | Candidate |
| --- | ---: | ---: |
| Whole-worker actual CPU (s) | 25,908.885 | 12,660.006 |
| Whole-worker wall (s) | 15,977.813 | 11,120.429 |
| Peak individual/descendant RSS (KiB) | 17,542,956 | 17,542,284 |
| Original writer CPU (s) | 251.892 | 263.565 |
| Annotation/compression/index CPU (s) | 25,624.548 | 12,363.233 |
| Complete operation CPU (s) | 25,876.440 | 12,626.798 |
| Complete operation wall (s) | 15,945.205 | 11,086.803 |

Whole-worker actual CPU fell **51.14%** (7.197 to 3.517 CPU-hours); wall time
fell 30.40%. Annotation/compression/index CPU fell 51.75%. **Process peak RSS
was effectively unchanged**, so this experiment establishes no process-memory
reduction. Sampled summed process-group RSS was 33,311,640 versus 17,559,052 KiB,
but that sum can double-count shared/forked pages and miss short peaks; it is not
a heap-memory estimate. Whole-job Slurm batch `MaxRSS` was separately
67,103,700 KiB, including preparation and validation, and must not be confused
with either worker's process peak. Shared final-QS preparation, sampling and
scientific comparison are outside the worker CPU measurements.

Legacy ran first and candidate second in the same allocation on `cl015`.
Contention was not identical: baseline burden jobs `19210990` and `19211105`
started on that node at 08:23:45 UTC during the candidate arm. Earlier and other
colocated work was not fully inventoried; `colocation-note.json` records this
limitation. This is one observed full-operation pair, not a projection to all
groups or total pipeline CPU. Complete pipeline scientific and high-memory-stage
validation remain required; production resource requests remain unchanged.

The accepted helper is integrated through an optional `--coverage-annotator`
argument, preserving the saved YAML schema. Without that argument, or for a
reference containing any legacy-ambiguous contig name, R uses the original
annotation path. The helper's only source difference from the benchmarked C++
is its first-line description comment; `runtime/coverage-annotation-promotion.json`
records both source hashes. The portable dispatch fixture was copied
byte-for-byte from the passed runtime source. The promoted R script differs
from its passed runtime source only by removing four inherited trailing-space
instances in the newly indented legacy fallback; the provenance records that
whitespace-only delta.

Pinned-container jobs `19206485` and `19207026` passed fixtures extracting the
actual candidate R functions and dispatch expression. They cover aggregate and
counts-only rows, empty indexed output, all 125 normalized contexts, terminal
bases, small contigs, fractional/large depths and quoted paths. An unsafe name
in an **uncovered** reference contig selects the unchanged legacy block before
the helper is invoked. A late malformed BED row, compressor failure and indexer
failure each stop R before a downstream save sentinel; source BEDs remain for
diagnosis. The final fixture also proves a valid FAI without a trailing newline
works under `options(warn=2)`. Frozen source, tests and results are in workspace
`runs/coverage-annotation-dispatch-v3/`. Initial job `19206382` failed before R
fixtures because the clean container PATH omitted `/hidef/bin/seqkit`; the
replacement explicitly includes `/hidef/bin`. These fixtures establish dispatch
and failure behavior; the full nuclear comparison above separately establishes
record/count/index equivalence. The portable test is retained as
`tests/test_coverage_annotation_dispatch.R`.

Runtime integration job `19207217` then completed the full nuclear-lenient
one-chunk burden operation with the fixture-passed R dispatch and compiled
helper. It used the unchanged original YAML and filtered chunk, reusing the
original baseline output from `runs/burdens-one-chunk/` without rerunning it.
Every scientific QS component matched exactly, including YAML, coverage,
spectra, burdens and sensitivity; only `/run_metadata` was ignored. The single
published coverage BED gzip and its tabix index were both **byte-identical**.
Original baseline file identities remained unchanged. Frozen source/helper,
measurements and comparison reports are in `runs/burdens-annotation-integrated/`.
The integrated worker used 737.208 CPU seconds, 694.923 wall seconds and
21,304,252 KiB peak process RSS; Slurm batch `MaxRSS` was separately 21,789,872 KiB.
It ran on `cs602`, whereas the reused baseline ran on `cl011`, so this is an
integration correctness result, **not a matched performance estimate**. The
complete pipeline/QS checks remain required before overall conclusions.

Combined workflow check `19216785` passed after integration with the inode-only
cache fix. Full-workflow preview emitted the same 19 prepared identities; a cold
compiler/dependency fixture completed four tasks, all four were cached on resume,
and changing only the helper source rebuilt the compiler and two burden-consumer
stubs while preserving an unrelated cached frontend. The actual preparation and
global-barrier workflow block remains byte-identical to the prior full candidate.
The fixture is a dependency/staging check, separate from the real R scientific
replays above. Evidence: `runs/coverage-annotation-integrated-workflow/`.

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
The completed same-allocation nuclear-lenient comparison is reported below.
The measured operations include script startup, filtering and QS serialization;
new shared-cache preparation and scientific comparison are outside these timings.
This is not a whole-pipeline or preparation-inclusive result. Wall times and
per-group results are retained in `comparisons/four-group-summary.json`.

## Completed: matched complete nuclear-lenient filtering, two pairs

Job `19196053` completed two alternating pairs of independent R workers on
`cl017`, using the same full chunk-1 extraction and identical prepared YAML for
both source versions. Pair 1 ran baseline then candidate; pair 2 reversed the
order. Sources, input identity, config hash, logs, measurements and comparisons
are frozen in workspace `runs/filter-paired/`.

| Pair | Baseline CPU (s) | Candidate CPU (s) | Baseline peak RSS (KiB) | Candidate peak RSS (KiB) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 4,097.287 | 1,971.739 | 26,488,152 | 14,499,304 |
| 2 | 2,923.138 | 2,223.242 | 26,487,048 | 14,047,936 |

Both complete scientific QS comparisons passed exactly, including the config,
calls, ordinary and germline-excluded coverage trackers, genome tracker and
region/molecule statistics. Each comparison reports zero failures, zero numeric
reviews and zero ignored differences. The comparator was configured with the
usual `/run_metadata` whitelist; it was not broadened for these comparisons.

Mean actual worker CPU fell from **3,510.213 to 2,097.491 seconds (40.25%)**;
mean wall time fell from 3,384.005 to 2,096.194 seconds (38.06%). Mean process
peak RSS fell from **26,487,600 to 14,273,620 KiB (46.11%)**. Both pairs are shown
because baseline runtime varied substantially even within the same allocation.
These measurements include R startup, filtering and QS serialization; shared
preparation and scientific comparison are outside worker timings. They do not
measure total pipeline CPU or preparation-inclusive speedup. Whole-job Slurm
batch `MaxRSS` was 33,548,272 KiB and is separate from the worker/descendant RSS
above; neither is a measurement of the complete pipeline's peak memory.

## Pending before overall conclusions

- Complete comparable baseline and candidate runs for both samples; compare
  total actual CPU, wall time, allocated CPU-hours, and high-memory-stage RSS.
- Compare complete nested scientific QS objects and every published scientific
  output, resolving all discrepancies and binary-file review findings.
- Complete paired measurements and exact-output checks for the remaining
  implementations. Full nuclear BED annotation passed its full-operation and
  integration checks. BAM dispatch passed its component and workflow gates;
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
