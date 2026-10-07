# R coverage annotation replay

Run this benchmark inside the pinned pipeline container on a SLURM compute node.
It accepts the immutable real-coverage packet produced by the earlier complete
annotation benchmark, so packet preparation does not reload an entire final
result object. The benchmark freezes its sources and refuses an existing output
directory.

```bash
python3 scripts/benchmark/benchmark_r_coverage_annotation.py \
  --packet /absolute/path/to/coverage-annotation-nuclear/coverage.qs2 \
  --baseline-script /absolute/path/to/frozen/calculateBurdens.original.R \
  --candidate-script /absolute/path/to/current/bin/calculateBurdens.R \
  --shared-functions /absolute/path/to/current/bin/sharedFunctions.R \
  --helper /absolute/path/to/preserved/annotateCoverage \
  --output /absolute/path/to/new/replay-directory \
  --chromosomes all --repeats 1 --allocated-cpus 4 \
  --historical-directory /absolute/path/to/original/pair01/legacy
```

Start with `--chromosomes chr22`; omit `--historical-directory` unless that
directory covers exactly the same chromosomes and rows. Both measured arms use
the original unchanged BED writer, enforced by AST identity and `chunk_runs =
1e7`. The C++ arm is a preserved baseline binary; the R arm evaluates the actual
three shared coverage definitions from the frozen candidate. All reference
inputs and original coverage values are shared. No OS-cache reset is performed.

The R annotator reads at most 65,536 coverage runs at a time and expands at most
one million genomic bases per window. `data.table::fwrite` formats output into an
R sink connected to bgzip. No whole-chromosome sequence, per-base R string table,
or uncompressed annotation BED is retained. Sequential double accumulation is
preserved by seeding each `rowsum` with the previous context totals. Contig-edge
contexts, ambiguous-base normalization and original depth tokens are preserved.
References with legacy-delimiter contig names, empty indexed contigs or gzip
compression select the existing legacy reference-BED path.

Results report complete-worker user plus system CPU (including waited-for
children), wall time, maximum individual-process RSS, sampled summed process
RSS, and operation phases. Packet selection/loading and the unchanged writer
can dominate RSS; the maximum individual-process measurement is not simultaneous
whole-job memory. Keep SLURM accounting separately. Compare timing across the
paired workers before interpreting historical timings from different nodes.

Validation compares exact ordered decompressed BED bytes, observed context
counts, tabix chromosome lists and boundary queries. Every newly generated
index must also exactly match its own independent rebuild after decompression.
An optional historical directory adds a direct original-output comparison.
Source hashes and the input packet's size/mtime must remain unchanged.

For standalone correctness tests, add `/hidef/bin` to `PATH` inside the pinned
container so seqkit is available. `tests/test_annotate_coverage.py --output-dir
NEW_DIRECTORY` creates independent legacy goldens without a compiled helper.
Run `tests/test_coverage_annotation_dispatch.R CANDIDATE.R SHARED.R GOLDENS
NEW_RESULT_DIRECTORY` to exercise the actual production dispatch, fallback,
bounded windows, precision and failures. An optional `--helper` still tests a
preserved historical C++ binary.

`tests/test_coverage_annotation_consumer.R CANDIDATE.R SHARED.R CPP_DIRECTORY
R_DIRECTORY NEW_RESULT_DIRECTORY` evaluates the actual downstream burden
context-consumer expressions on the paired real counts. It checks exact ordered
tables, factors, fractions and QS2 round trips for aggregate rows and synthetic
asymmetric-strand combinations. This is a producer/consumer boundary test, not
an entire burden-task or full-workflow result comparison.
