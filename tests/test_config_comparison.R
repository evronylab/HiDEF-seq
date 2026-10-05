#!/usr/bin/env Rscript
# Exercise generated-YAML parsing and exact scientific-key comparison.
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if(length(args)) args[[1]] else ".", mustWork = TRUE)
directory <- tempfile(pattern = "config-compare-")
dir.create(directory)
original <- file.path(directory, "original.yaml")
lines <- c("analysis_id: fixture", "cache_dir: /original/cache", "genome_fasta: /reference.fa",
           "samples:", "  - sample_id: sample1", "    individual_id: individual1",
           "filtergroups:", "  - filtergroup: strict", "    threshold: 10")
writeLines(lines, original)
compare <- function(candidate_lines, expected_status, flags = character()) {
  candidate <- tempfile(tmpdir = directory, fileext = ".yaml")
  report <- tempfile(tmpdir = directory, fileext = ".tsv")
  writeLines(candidate_lines, candidate)
  command <- c("--vanilla", file.path(repo, "scripts", "benchmark", "compare_config.R"),
               original, candidate, report, flags)
  status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(command)))
  stopifnot(identical(as.integer(status), as.integer(expected_status)), file.exists(report))
  stopifnot(identical(readLines(original), lines))
}
extras <- c("reference_cache_dir: /prepared/library", "reference_summary_file: /prepared/reference.qs2",
            "cache_artifacts: {}", "germline_coverage_filters: []")
compare(c(lines, extras), 0L)
compare(c(sub("threshold: 10", "threshold: 11", lines, fixed = TRUE), extras), 1L)
compare(c(lines, "unexpected_torch_runtime: true"), 1L)
compare(c(sub("/original/cache", "/different/cache", lines, fixed = TRUE), extras), 1L)
compare(c(sub("/original/cache", "/different/cache", lines, fixed = TRUE), extras), 0L,
        "--metadata-key=cache_dir")
unlink(directory, recursive = TRUE)
cat("PASS: configr parses effective YAML, explicit metadata additions only, exact scientific keys, input unchanged\n")
