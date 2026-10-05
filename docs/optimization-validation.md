# Optimization validation ledger

This living ledger records the optimization measurements available on 2026-10-04.
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
less peak RSS**. It is not a repeated benchmark or a whole-pipeline result.
Measurements are in workspace `runs/chunk-replay/extract.json` and
`runs/extract-candidate/extract.json`; the exact scientific comparison report is
`runs/extract-candidate/comparison.tsv`. Allocated CPU-hours above cover the timed
command, not its containing job's setup or validation overhead.

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
Nuclear results remain a separate pending experiment. Records and strict
comparison reports are in workspace `runs/filter-chain-mito/`.

## Running: complete current-candidate chunk filtering

Job `19194137` replays all four chromosome/filtergroup combinations for the
original LIB1 chunk-1 extraction with 24 GiB, two CPUs, and a four-hour limit.
Workspace `runs/filter-candidate-current/` contains the frozen `code/bin` and
benchmark scripts, per-file SHA256 manifest, prepared config, logs and metrics.
The extraction argument remains the exact relative `extractCalls.chunk1.qs2`
used by the original replay; its symlink resolves to the original extraction.

Preflight requires every original parsed configuration value to match and exactly
two added preparation fields. After each candidate group, available complete
original outputs are compared with only
`/config/yaml.config/reference_summary_file` and
`/config/yaml.config/germline_coverage_filters` whitelisted. The full config,
inherited run metadata, scientific schemas and within-group row order otherwise
remain strict. `comparisons/status.json` distinguishes pending from passed groups;
the saved comparator can be rerun on compute as remaining originals arrive.
No complete-filter equivalence or performance result is claimed yet.

## Pending before overall conclusions

- Complete comparable baseline and candidate runs for both samples; compare
  total actual CPU, wall time, allocated CPU-hours, and high-memory-stage RSS.
- Compare complete nested scientific QS objects and every published scientific
  output, resolving all discrepancies and binary-file review findings.
- Complete paired measurements and exact-output checks for the remaining
  implementations. Combined filtering and BED annotation alternatives remain
  experiments until separately accepted; BAM dispatch integration requires its
  real-data and workflow checks. Coordinate-sort reuse was rejected above.
- Validate cache/resume behavior, effective YAML parsing and scientific keys,
  including missing prepared artifacts and concurrent launches.

The baseline source revision is
`6b6598b236d1f4a96e9e34597f7c7fa3c2cb2a3a`. Candidate revisions must be recorded with
each subsequent measurement. This initial ledger transcribes the completed
reference and QS2 measurements from the workspace `WORK_LOG.md`; isolated runs
and accounting artifacts are retained under the workspace `runs/` directory.
