# Complete coverage annotation replay

This experimental benchmark does not change the production workflow. Run it in
the pinned container on an allocated compute node. `annotateCoverage` must
already have been compiled against that container's HTSlib.

```bash
python3 scripts/benchmark/benchmark_coverage_annotation.py \
  --input /absolute/path/to/sample.outputResults.qs2 \
  --baseline-script /absolute/path/to/baseline/bin/calculateBurdens.R \
  --helper /absolute/path/to/annotateCoverage \
  --chromgroup 1-22X --filtergroup strict --chromosomes chr22 \
  --output /absolute/path/to/new/replay-directory \
  --repeats 2 --allocated-cpus 4
```

Start with chr22; `--chromosomes all` retains the selected chromgroup's original
coverage. Loading the full LIB1 final QS2 previously peaked above 40 GiB, so use
96 GiB for the initial run. The reference trinucleotide BED is reused from the
saved configuration, with no regeneration in either measured arm.

The preparation phase selects real coverage from
`bam.gr.filtertrack.bytype.coverage_tnc`, reports every available metadata row,
and saves one common packet. Final artifacts retain combined-orientation
coverage; selecting a strand-specific row without its original strand RleLists
fails explicitly. Preparation CPU and memory are reported separately.

Both arms execute the original BED-writing expression extracted from the frozen
baseline script. The extractor requires `chunk_runs = 1e7`. The legacy arm then
executes the original union, per-base expansion, reference intersection, row
intersection, awk counting, bgzip and tabix. The candidate replaces annotation
and counting with the experimental FASTA helper and retains bgzip and tabix.
Paired repeats alternate their execution order. The harness neither drops OS
caches nor claims a cold-cache comparison.

`results.json` records operation-only CPU/wall phases and complete-process
metrics, including packet loading. CPU means user plus system time for the
process and waited-for descendants. RSS includes both the maximum individual
process measurement and a separate 250 ms sampled sum across the process group;
the latter can double-count shared pages and miss short peaks. Retain scheduler
accounting independently when assessing the allocated job.

Validation runs outside the measured arms. It compares exact decompressed BED
bytes in order, numeric context-count maps including N contexts and dot, the
expected count-file set, tabix contig lists, and indexed first/last-base queries.
It permits gzip/index encoding differences. Tiny adversarial reference fixtures
are tested separately by `tests/test_annotate_coverage.py`.

A passing chr22 replay is a component result. Adoption still requires a full
nuclear operation replay and integration/output validation; no production
performance claim follows from the helper fixtures alone.
