#!/usr/bin/env Rscript
# Run from repository root with the pipeline's R environment.
source("scripts/benchmark/compare_qs2.R")
check <- function(a, b, expected, ...) {
  result <- scientific_compare(a, b, ...)
  stopifnot(identical(unname(result$counts[c("failure", "review")]), as.integer(expected)))
}
check(list(x = 1:3), list(x = 1:3), c(0, 0))
check(list(x = 1:3), list(x = c(1L, 3L, 2L)), c(1, 0))
attribute_reference <- structure(1:3, alpha = "same", beta = "before")
attribute_candidate <- structure(1:3, beta = "after", alpha = "same")
check(attribute_reference, attribute_candidate, c(1, 0))
check(factor("a", levels = c("a", "b")), factor("a", levels = c("b", "a")), c(2, 0))
check(matrix(1:4, 2, dimnames = list(c("a", "b"), NULL)),
      matrix(1:4, 2, dimnames = list(c("b", "a"), NULL)), c(1, 0))
check(c(1, NA_real_, NaN, Inf), c(1, NaN, NA_real_, Inf), c(1, 0))
check(1, 1 + 1e-13, c(0, 1))
check(0, 1e-20, c(1, 0))
check(1, 1 + 1e-10, c(1, 0))
check(list(run_metadata = "old", result = 1), list(run_metadata = "new", result = 1), c(0, 0))
check(list(yaml.config = list(threshold = 1)),
      list(yaml.config = list(threshold = 1, reference_summary_file = "/prepared/ref.qs2")),
      c(0, 0), ignore = "/yaml.config/reference_summary_file")
check(list(yaml.config = list(threshold = 1)),
      list(yaml.config = list(threshold = 2, reference_summary_file = "/prepared/ref.qs2")),
      c(1, 0), ignore = "/yaml.config/reference_summary_file")
reference <- data.frame(chromgroup = c("b", "b", "a", "a"), value = c(1, 2, 3, 4))
candidate <- reference[c(3, 4, 1, 2), ]; rownames(candidate) <- NULL
check(reference, candidate, c(2, 0))
check(reference, candidate, c(0, 0), chromgroup_blocks = TRUE)
candidate <- reference[c(4, 3, 1, 2), ]; rownames(candidate) <- NULL
check(reference, candidate, c(1, 0), chromgroup_blocks = TRUE)
if (requireNamespace("GenomicRanges", quietly = TRUE)) {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 4), width = 2))
  check(ranges, ranges, c(0, 0))
  other <- ranges; BiocGenerics::start(other) <- c(2, 4)
  stopifnot(scientific_compare(ranges, other)$counts[["failure"]] > 0L)
  check(S4Vectors::Rle(c(1L, 1L, 2L)), S4Vectors::Rle(c(1L, 2L, 2L)), c(1, 0))
}
if (requireNamespace("qs2", quietly = TRUE)) {
  reference_path <- tempfile(fileext = ".qs2")
  candidate_path <- tempfile(fileext = ".qs2")
  stage_path <- tempfile(pattern = "reference-stage-")
  report_path <- tempfile(fileext = ".tsv")
  qs2::qs_save(list(run_metadata = "old", result = reference), reference_path)
  reordered <- reference[c(3, 4, 1, 2), ]; rownames(reordered) <- NULL
  qs2::qs_save(list(run_metadata = "new", result = reordered), candidate_path)
  stopifnot(comparison_main(c("stage", reference_path, stage_path)) == 0L)
  stopifnot(comparison_main(c("compare", stage_path, candidate_path, report_path,
                             "--chromgroup-blocks")) == 0L)
  stopifnot(file.exists(report_path))
  unlink(c(reference_path, candidate_path, stage_path, report_path), recursive = TRUE)
}
cat("Scientific comparator tests passed\n")
