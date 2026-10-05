#!/usr/bin/env Rscript
#Run on a compute node in the pipeline container:
#Rscript tests/test_extract_sa_rle_optimization.R [repository root]
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(tidyverse)
  library(data.table)
})
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
repo <- if(length(args)) args[[1]] else "."
script <- file.path(repo, "bin", "extractCalls.R")
expressions <- parse(script)
find_assignments <- function(node, name) {
  if(is.call(node) && identical(node[[1]], as.name("<-"))) {
    return(if(identical(node[[2]], as.name(name))) list(node) else list())
  }
  if(is.expression(node) || (is.call(node) && is.symbol(node[[1]]) &&
                           as.character(node[[1]]) %in% c("{", "if", "for"))) {
    return(unlist(lapply(as.list(node), find_assignments, name = name), recursive = FALSE))
  }
  list()
}
assignment <- function(node, name) {
  hits <- find_assignments(node, name)
  stopifnot(length(hits) == 1L)
  hits[[1]]
}
for(name in c("decode_sa_rle", "subset_tag_positions")) eval(assignment(expressions, name))
legacy_decode <- function(tag) Rle(inverse.rle(list(lengths = tag[c(TRUE, FALSE)], values = tag[c(FALSE, TRUE)])))
capture <- function(expr) {
  tryCatch(list(value = force(expr)), error = function(e) list(error = conditionMessage(e)))
}
valid_tags <- list(
  empty_integer = integer(), empty_double = numeric(),
  single_run = c(5L, 3L), adjacent_equal = c(2L, 3L, 4L, 3L, 2L, 7L),
  zero_runs = c(0L, 9L, 2L, 3L, 0L, 8L, 3L, 3L, 0L, 1L),
  all_zero_runs = c(0L, 5L, 0L, 8L),
  doubles = c(2, 1.5, 3, 7.25, 1, 0), fractional_lengths = c(2.9, 3, 1.1, 4),
  na_values = c(2, NA_real_, 1, NA_real_, 2, 4),
  named_tag = setNames(c(2L, 3L, 2L, 8L), c("l1", "v1", "l2", "v2"))
)
for(case_name in names(valid_tags)) {
  tag <- valid_tags[[case_name]]
  baseline_outcome <- capture(legacy_decode(tag))
  candidate_outcome <- capture(decode_sa_rle(tag))
  if(!identical(candidate_outcome, baseline_outcome)) stop(paste("sa decode differs for", case_name))
  #Some unusual length encodings (including fractional counts) are rejected by
  #rep.int in the legacy decoder. Compare that failure instead of assuming a
  #dense reference vector can always be constructed.
  if(!is.null(baseline_outcome$error)) next
  reference <- baseline_outcome$value
  candidate <- candidate_outcome$value
  stopifnot(identical(candidate, reference), identical(rev(candidate), Rle(rev(as.vector(reference)))))
  for(reverse_strand in c(FALSE, TRUE)) {
    compressed <- if(reverse_strand) rev(candidate) else candidate
    dense <- if(reverse_strand) rev(as.vector(reference)) else as.vector(reference)
    positions <- list(integer(), seq_along(dense), as.numeric(seq_along(dense)),
                      c(1L, 1L), c(0L, 1L), c(NA_integer_, 1L), c(1L, length(dense) + 1L),
                      c(-1L, -2L), c(-1L, 1L), c(1.8, 2.1), c(TRUE, FALSE),
                      c(Inf, 1), c(NaN, 1), NULL, setNames(c(1L, 2L), c("a", "b")))
    for(index in positions) stopifnot(identical(capture(subset_tag_positions(compressed, index)), capture(dense[index])))
  }
}
for(tag in list(NULL, c(1L, 5L, 2L), c(-1L, 5L), c(NA_integer_, 5L), c(Inf, 5))) {
  baseline <- capture(legacy_decode(tag))
  candidate <- capture(decode_sa_rle(tag))
  stopifnot(!is.null(baseline$error), identical(candidate, baseline))
}
#Plain sm/sx vectors retain their original classes, names and indexing behavior.
for(tag in list(1:6, as.numeric(1:6), setNames(1:6, letters[1:6]), numeric(), NULL)) {
  for(index in list(c(1L, 3L, 3L), c(NA_integer_, 7L), -1L, integer())) {
    stopifnot(identical(subset_tag_positions(tag, index), tag[index]))
  }
}
cat("PASS: direct/legacy sa decoding, integer/double types, adjacent/zero/empty runs, strand reversal, unusual indices and errors\n")

#Load the actual nested SBS/MDB extractor without executing BAM/CIGAR work.
extract_calls_definition <- assignment(expressions, "extract_calls")[[3]]
extract_body <- extract_calls_definition[[3]]
strand_levels <- c("+", "-")
eval(assignment(extract_body, "extract_sasmsx"))
tags <- list("101_+" = c(2L, 3L, 4L, 6L), "101_-" = c(3L, 8L, 3L, 2L),
             "102_+" = c(1, 2.5, 5, 7.5))
compressed <- map(tags, decode_sa_rle)
compressed[["101_-"]] <- rev(compressed[["101_-"]])
dense <- map(tags, ~ as.vector(legacy_decode(.x)))
dense[["101_-"]] <- rev(dense[["101_-"]])
queries <- list("101_+" = c(1L, 3L, 6L), "101_-" = c(2L, 4L), "102_+" = c(1L, 5L))
for(tag_name in c("sa", "sm", "sx")) {
  candidate_input <- if(tag_name == "sa") compressed else dense
  reference <- extract_sasmsx(tag_name, dense, queries)
  candidate <- extract_sasmsx(tag_name, candidate_input, queries)
  stopifnot(identical(candidate, reference), is.list(candidate[[tag_name]]))
}
#The opposite-strand lookup swaps keys but does not reverse already-oriented
#arrays again; coordinate conversion and resulting row/column order stay exact.
opposite_dense <- setNames(dense, chartr("+-", "-+", names(dense)))
opposite_compressed <- setNames(compressed, chartr("+-", "-+", names(compressed)))
opposite_queries <- list("101_-" = c(1L, 4L), "101_+" = c(2L, 6L))
conversion <- tibble(zm = c(101L, 101L, 101L, 101L),
                     strand = factor(c("-", "-", "+", "+"), levels = strand_levels),
                     start = c(1L, 4L, 2L, 6L), end = c(1L, 4L, 2L, 6L),
                     start_queryspace = c(5L, 2L, 4L, 1L), end_queryspace = c(5L, 2L, 4L, 1L))
stopifnot(identical(extract_sasmsx("sa.opposite_strand", opposite_compressed, opposite_queries, conversion),
                    extract_sasmsx("sa.opposite_strand", opposite_dense, opposite_queries, conversion)))
#Even legacy errors/empty outputs for no queries must remain the same.
stopifnot(identical(capture(extract_sasmsx("sa", compressed, list())),
                    capture(extract_sasmsx("sa", dense, list()))))
cat("PASS: actual SBS/MDB sa/sm/sx extraction and opposite-strand coordinate-join outputs\n")

#Run the actual six data.table indel tag assignments on ordinary vs compressed
#sa inputs, covering insertion-base ranges, deletion flanks, duplicates and OOB.
lines <- readLines(script)
indel_lines <- lines[grepl("indels_queryspace_pos[", lines, fixed = TRUE) |
                       grepl("indels_queryspace_pos.opposite_strand[", lines, fixed = TRUE)]
indel_lines <- indel_lines[grepl(":=", indel_lines, fixed = TRUE)]
stopifnot(length(indel_lines) == 12L, sum(grepl("_val :=", indel_lines, fixed = TRUE)) == 6L)
indel_code <- parse(text = indel_lines)
indel_queries <- data.table(zm_strand = c(rep("101_+", 4), rep("101_-", 2), rep("102_+", 2)),
                            start_end_queryspace = c(rep("2_4", 3), "1_6", "1_6", "1_6", "1_2", "1_2"),
                            pos_queryspace = c(2L, 3L, 4L, 6L, 1L, 6L, 1L, 2L))
run_indel <- function(sa, query, vector_tags) {
  env <- new.env(parent = globalenv())
  env$indels_queryspace_pos <- copy(query)
  env$indels_queryspace_pos.opposite_strand <- copy(query)[, zm_strand := chartr("+-", "-+", zm_strand)]
  env$sa.input <- sa
  env$sm.input <- vector_tags
  env$sx.input <- map(vector_tags, ~ .x * 0)
  env$sa.input.opposite_strand <- setNames(sa, chartr("+-", "-+", names(sa)))
  env$sm.input.opposite_strand <- setNames(env$sm.input, chartr("+-", "-+", names(dense)))
  env$sx.input.opposite_strand <- setNames(env$sx.input, chartr("+-", "-+", names(dense)))
  eval(indel_code, env)
  map(list(env$indels_queryspace_pos, env$indels_queryspace_pos.opposite_strand), function(x) {
    x[, .(sa = list(sa_val), sm = list(sm_val), sx = list(sx_val)), by = .(zm_strand, start_end_queryspace)] %>% as_tibble()
  })
}
for(vector_tags in list(map(dense, as.integer), map(dense, as.numeric))) {
  for(query in list(indel_queries, indel_queries[0],
                    copy(indel_queries)[, pos_queryspace := replace(pos_queryspace, c(1, 2), c(NA_integer_, 9L))])) {
    stopifnot(identical(capture(run_indel(map(vector_tags, Rle), query, vector_tags)),
                        capture(run_indel(vector_tags, query, vector_tags))))
  }
}
cat("PASS: actual insertion/deletion tag assignments and aggregated list-of-vector outputs, including empty/OOB queries\n")
