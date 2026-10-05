#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
repo <- if(length(args)) args[[1]] else "."
suppressPackageStartupMessages(library(S4Vectors))
source(file.path(repo, "scripts/benchmark/sa_interval_lookup.R"))
#Exercise the actual production helper too; the benchmark mirror alone cannot
#prove that a later integration preserves its guards and exact return types.
expressions <- parse(file.path(repo, "bin/extractCalls.R"))
hits <- vapply(expressions, function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                 identical(x[[2]], as.name("subset_tag_positions")), logical(1))
stopifnot(sum(hits) == 1L)
eval(expressions[[which(hits)]])
production_lookup <- subset_tag_positions
benchmark_lookup <- subset_tag_positions_interval
subset_tag_positions_interval <- function(tag, positions) {
  actual <- production_lookup(tag, positions)
  stopifnot(identical(actual, benchmark_lookup(tag, positions)))
  actual
}
capture <- function(expr) tryCatch(list(value = force(expr)), error = function(e) list(error = conditionMessage(e)))
tags <- list(Rle(integer()), Rle(numeric()), Rle(c(1L, 2L, 3L), c(2L, 0L, 3L)),
             Rle(c(1.5, NA_real_, 8), c(2L, 1L, 3L)), Rle(c(TRUE, FALSE), c(3L, 2L)),
             Rle(factor(c("A", "B")), c(2L, 3L)), Rle(c("A", "C"), c(2L, 3L)))
for(tag in tags) for(tag in list(tag, rev(tag))) {
  dense <- as.vector(tag)
  indices <- list(integer(), numeric(), seq_along(dense), rev(seq_along(dense)), c(2L, 1L, 2L),
                  c(NA_integer_, 1L), c(0L, 1L), c(-1L, -2L), c(-1L, 1L),
                  c(1L, length(tag) + 1L), c(1.8, 2.1), c(TRUE, FALSE), c(Inf, 1),
                  NULL, setNames(c(1L, 2L), c("a", "b")))
  for(index in indices) stopifnot(identical(capture(subset_tag_positions_interval(tag, index)), capture(dense[index])))
}
# Large compressed vectors remain valid without integer cumulative overflow;
# validate sparse endpoints without ever constructing their dense expansion.
long <- Rle(c(7L, 8L), c(.Machine$integer.max, 2))
stopifnot(identical(subset_tag_positions_interval(long, c(1, .Machine$integer.max,
                                                       .Machine$integer.max + 1, .Machine$integer.max + 2)),
                    c(7L, 7L, 8L, 8L)))
for(tag in list(1:6, as.numeric(1:6), setNames(1:6, letters[1:6]), NULL)) {
  for(index in list(1:3, integer(), NA_integer_, -1L))
    stopifnot(identical(capture(subset_tag_positions_interval(tag, index)), capture(tag[index])))
}
cat("PASS: interval lookup equals dense semantics, including boundaries, reversal, classes, unusual indices and long run endpoints\n")
