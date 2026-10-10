#!/usr/bin/env Rscript
# Run on a compute node in the pipeline container.
args <- commandArgs(TRUE)
repo <- if(length(args)) args[[1L]] else '.'
suppressPackageStartupMessages({library(GenomicRanges); library(plyranges); library(tidyverse)})
options(warn = 2)
load_functions <- function(path, environment) {
 for(expr in parse(path)) if(is.call(expr) && identical(expr[[1L]], as.name('<-')) &&
    is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name('function'))) eval(expr, environment)
}
candidate <- new.env(parent = globalenv())
load_functions(file.path(repo, 'bin/sharedFunctions.R'), candidate)
load_functions(file.path(repo, 'bin/calculateBurdens.R'), candidate)
load_functions(file.path(repo, 'bin/extractCalls.R'), candidate)
reference <- new.env(parent = globalenv())
sys.source(file.path(repo, 'tests/fixtures/burden_native_reference_ea8d06f.R'), reference)
capture <- function(expr) tryCatch(list(value = force(expr)), error = function(e) list(error = conditionMessage(e)))
#Recycling is part of interaction's contract even though data-frame columns
#normally have equal lengths; a lone factor preserves an explicit NA level.
for(columns in list(list(factor(c('a','b')),factor(c('w','x','y','z'))),
                    list(factor(c('x',NA),exclude=NULL)),
                    list(c('a','b'), c('x','y','z')))) {
  stopifnot(identical(capture(interaction(columns,drop=TRUE)),
                      capture(candidate$observed_allele_interaction(columns))))
}
set.seed(159)
labels <- c('', '.', 'A', 'A.B', 'B.C', 'C', 'NA', NA_character_)
for(i in seq_len(500L)) {
 columns <- replicate(sample(1:5, 1), factor(sample(labels, 20L, TRUE),
     levels = sample(labels), exclude = if(i %% 3L) NA else NULL), simplify = FALSE)
 expected <- capture(interaction(columns, drop = TRUE))
 actual <- capture(candidate$observed_allele_interaction(columns))
 if(is.null(expected$error)) {
   if(!identical(expected, actual)) {dput(columns); print(expected); print(actual); stop('allele factor mismatch')}
 } else {
   # base::interaction can fail during collision remapping when rows contain NA.
   # The sparse encoder still returns the observed labels with NA propagation.
   expected_labels <- do.call(paste, c(lapply(columns, as.character), sep = '.'))
   missing <- Reduce(`|`, lapply(columns, is.na))
   expected_labels[missing] <- NA_character_
   stopifnot(is.null(actual$error), identical(as.character(actual$value), expected_labels))
 }
}
# Cardinalities whose Cartesian product exceeds the integer limit are still
#valid observed factors. No Cartesian allocation or code multiplication occurs.
x <- factor(seq_len(50000L), levels = seq_len(50000L))
large <- candidate$observed_allele_interaction(list(x, x))
stopifnot(identical(as.integer(large), seq_len(50000L)), nlevels(large) == 50000L)
cat('PASS: 500 randomized separator/NA/factor combinations and Cartesian-overflow cardinality\n')
for(i in seq_len(150L)) {
 extent <- sample(c(1L, 5L, 20L), 1L)
 circular <- sample(c(TRUE, FALSE, NA), 1L)
 known <- sample(c(TRUE, FALSE), 1L)
 starts <- sample(-40:40, 6L)
 widths <- sample(0:80, 6L)
 gr <- suppressWarnings(GRanges(rep('a', 12L), IRanges(rep(starts, 2L), width = rep(widths, 2L)),
   strand = rep(c('+', '-'), each = 6L), seqinfo = Seqinfo('a', if(known) extent else NA_integer_, isCircular = circular),
   run_id = 'r', zm = rep(1:6, 2L), bc_orientation = 'x-x'))
 if(i %% 3L == 0L) {
   extra <- gr[1L]; strand(extra) <- '*'; gr <- suppressWarnings(c(gr, extra, extra))
 }
 cov <- reference$calc_duplex_coverage(gr)
 stopifnot(identical(candidate$calc_duplex_coverage(gr), cov))
 query <- GRanges(rep('a', length(cov[[1L]])), IRanges(seq_len(length(cov[[1L]])), width = 1L))
 expected <- list(start = reference$gr_1bp_cov(query, cov), end = reference$gr_1bp_cov(rev(query), cov))
 actual <- candidate$sensitivity_site_counts(gr, query, rev(query))
 if(!identical(expected, actual)) {print(gr);print(expected);print(actual);stop('circular/clipping mismatch')}
}
#Index normalization and assignment are independent for the two query sides.
track <- GRanges(rep('a',2L),IRanges(2L,5L),strand=c('+','-'),seqinfo=Seqinfo('a',10L),
                 bc_orientation='x-x',run_id='r',zm=1L)
cov <- reference$calc_duplex_coverage(track)
query <- function(p) suppressWarnings(GRanges(rep('a',length(p)),IRanges(p,width=1L)))
for(a in list(0L,c(0L,2L),c(0L,2L,3L),-2L,c(-2L,-3L),c(0L,-2L))) {
 for(b in list(1L,0L,c(1L,3L),-3L)) {
   expected <- capture(list(start=reference$gr_1bp_cov(query(a),cov),end=reference$gr_1bp_cov(query(b),cov)))
   actual <- capture(candidate$sensitivity_site_counts(track,query(a),query(b)))
   stopifnot(identical(expected,actual))
 }
}
cat('PASS: 150 randomized clipped/circular/multiturn/unknown-length cases and independent query-side indexing\n')
for(i in seq_len(100L)) {
 values <- sample(c(NA_real_, 1:8), 10L, TRUE)
 lengths <- sample(0:20, 10L, TRUE)
 tag <- Rle(values, lengths)
 positions <- list(sample(-250:250, 20L), sample(c(NA,TRUE,FALSE), 300L, TRUE),
                   sample(c(NA_real_, NaN, Inf, -Inf, runif(30L,-250,250)), 25L), c('a', NA_character_),
                   matrix(c(1L,2L,NA_integer_,400L),2L))
 for(index in positions) stopifnot(identical(capture(candidate$subset_tag_positions(tag,index)), capture(as.vector(tag)[index])))
}
cat('PASS: 500 randomized compressed-vector indexing cases\n')
