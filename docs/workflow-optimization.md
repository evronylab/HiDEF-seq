# Computational optimizations

The optimization branch changes repeated work and intermediate allocations while
retaining the scientific output schemas and one final QS2 per sample. TSV, VCF,
PDF and indexed coverage BED outputs are retained. The
[validation ledger](optimization-validation.md) records completed component
benchmarks, complete-workload measurements and acceptance of all 711 required
scientific outputs under the agreed comparison rules.

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

`artifactCache.py` runs inside preparation tasks using Python's standard library.
It serializes simultaneous builds of the same identity, accepts only explicitly
declared closed products after a successful command, checksums those products,
and copies them into a temporary sibling directory. It writes
`manifest.complete.json` there and atomically renames the complete directory
into place. Failed builds/copies cannot create a complete entry. Existing corrupt
entries fail verification and are never silently overwritten. Inspect and move
a corrupt entry aside before rebuilding it.

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
count remains the smaller of requested chunks and enumerated IDs. A small
Nextflow process compiles `splitBamByZmw.cpp` once inside the pinned container,
with the source supplied as a content-hashed task input. All sample dispatches
reuse that executable.

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

Up to 128 output writers share an HTSlib compression thread pool. Larger chunk
counts use additional sequential input passes over bounded writer groups,
without changing chunk assignments. The pool limits compression workers;
HTSlib may also create background I/O threads per stream. Each output receives
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
The pinned Nextflow 26.04.0 replay passed actual execution and publication for
12-chunk and one-chunk samples, skipped the empty sample, and cached all six tasks
on `-resume` (job 19194083). All 13 published BAMs retained hardlinks to their
work outputs. Full workflow preview and the 57-field scientific configuration
comparison also passed (job 19194084). Evidence is in
`runs/dispatch-workflow-v3/{validation.json,first.trace.tsv,resume.trace.tsv}`.

The full LIB1 replay produced all 60 chunks in 2,161.59 actual CPU seconds
(1,092.32 wall seconds; 526,472 KiB maximum process RSS). Quickcheck and PBI/BAI
indexing passed for all chunks and used another 1,499.20 actual CPU seconds.
Chunk 1's complete ordered SAM record stream matched the legacy chunk SHA256
`b09dc1de7a0fdc777326febf672a3dbfe1163190a119e1dc9f7e9d15948b8490`.
These measurements and hashes are retained in `runs/dispatch-benchmark/`. This
is an exact real-data record comparison for chunk 1, not a claim that every
real-data chunk was compared or that a complete optimized pipeline was validated.

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
the original coverage runs, a helper reads the reference FASTA directly to emit
the same per-base BED records and context counts. Nextflow compiles the helper
once using the compiler and HTSlib already in the pinned container, then stages
the executable as an explicit burden-task dependency. No additional software
dependency or output format is introduced.

The original annotation path remains available for reference contig names that
the legacy parser treats ambiguously, and for standalone `calculateBurdens.R`
invocations without `--coverage-annotator`. Helper, compression and index failures
stop the task before results can be saved. Full nuclear validation matched all
eight BEDs, count tables and indexed queries, while the integrated one-chunk
replay matched the complete scientific QS and byte-identical BED/index files.
The full nuclear operation used 51.14% less worker CPU with effectively unchanged
process peak RSS; see the validation ledger for measurement scope and node
colocation. Complete pipeline validation remains required.

### Final output

`outputResults.R` skips expensive list-column conversions that the following
germline-table pivot discards. Its native `VariantAnnotation::writeVcf()` call uses
`nchunk = 100000L` to reduce temporary VCF formatting allocations. This is the
existing writer's export buffer; it does not divide analysis calls into separate
jobs or split the final QS2. The measured final object fits as a single object,
so no separate QS2 components or new reader API are introduced.

## Focused local validation

Coordinate-sort reuse was tested and rejected. The single-input fixture matched,
but every tested multi-input order (`AB`, `BA`, `ABC`, `CBA`) changed SAM record
order at coordinate ties. Record multisets and headers excluding program records
matched, and indices validated; those checks are insufficient because ordered
records must also match. Keep the existing coordinate sorts and `pbmerge`, with
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

## Pending work — awaiting user instruction

Recorded on 2026-10-06. These items are proposals only. Do not implement or
benchmark them until the user instructs us to proceed.

1. **Port `artifactCache.py` to R.** Put reusable functions in
   `sharedFunctions.R` where appropriate. Preserve cache identity, locking,
   verification, atomic publication, failure handling and restore behavior.
   Measure actual CPU time, elapsed time and peak memory against the Python
   implementation, including cold preparation and warm restoration.
2. **Move `annotateCoverage.cpp` functionality into R.** Integrate it into the
   relevant R scripts or `sharedFunctions.R`, as appropriate. Preserve coverage
   records, context counts, ordering and indexed output behavior. Measure CPU
   time, elapsed time and peak memory against the current C++ implementation,
   including the complete annotation/compression/index operation.
3. **Evaluate `data.table` / `dtplyr` in expensive R operations.** Explore
   material reductions in CPU-hours and peak memory from fewer copies and more
   efficient joins, grouping and aggregation. Preserve scientific values,
   output types and required ordering; avoid broad rewrites without measured
   benefit.
4. **Evaluate targeted Rcpp functions for heavy steps.** Consider only candidates
   with the potential for dramatic CPU-hour and memory reductions. Do not add
   compiled helpers for incremental gains. Compare complete affected operations,
   including conversion and allocation overhead, with the current implementation.

For all four items, retain the agreed scientific comparison rules and record
performance regressions as well as improvements. Judge aggregate performance by
actual user-plus-system CPU-hours across the pipeline, with significantly lower
peak memory in the currently expensive steps. An individual step may become
more than 5% slower; there is no per-step 5% acceptance limit. Preserve the single
final QS2 and the existing TSV, VCF, PDF and indexed BED outputs.
