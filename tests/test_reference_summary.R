#!/usr/bin/env Rscript
# Tiny reference and consumer equivalence fixtures, run in the pipeline R image.
suppressPackageStartupMessages({
  library(BSgenome)
  library(GenomicRanges)
  library(Biostrings)
  library(tidyverse)
})
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
repo <- if(length(args)) args[[1]] else "."
source(file.path(repo, "bin", "sharedFunctions.R"))
source(file.path(repo, "bin", "referenceSummaryFunctions.R"))

sequences <- DNAStringSet(c(nfree = "ACGTACGT", allN = "NNNNNN",
                           linear = "NNACGNNTN", circular = "NACGTN"))
reference_info <- Seqinfo(names(sequences), width(sequences),
                         isCircular = c(FALSE, FALSE, FALSE, TRUE), genome = "fixture")
# A small list with the BSgenome accessors used by the production helper avoids
# forging/installing a package and still tests full Seqinfo and circular flags.
setOldClass("reference_summary_fixture")
setMethod("seqnames", "reference_summary_fixture", function(x) names(x))
setMethod("seqinfo", "reference_summary_fixture", function(x) attr(x, "seqinfo"))
genome <- structure(as.list(sequences), class = "reference_summary_fixture", seqinfo = reference_info)
expected_n <- GRanges(c("allN", "linear", "linear", "linear", "circular", "circular"),
                     IRanges(c(1L, 1L, 6L, 9L, 1L, 6L), c(6L, 2L, 7L, 9L, 1L, 6L)),
                     strand = "*", seqinfo = reference_info)
stopifnot(identical(reference_n_ranges(genome), expected_n))
nfree <- structure(as.list(sequences["nfree"]), class = "reference_summary_fixture",
                   seqinfo = reference_info["nfree"])
stopifnot(identical(reference_n_ranges(nfree), GRanges(seqinfo = reference_info["nfree"])))
empty <- structure(list(), class = "reference_summary_fixture", seqinfo = Seqinfo())
stopifnot(identical(reference_n_ranges(empty), GRanges(seqinfo = Seqinfo())))
counts <- reference_trinucleotide_counts(genome)
chromosome_sets <- list(names(sequences), rev(names(sequences)), "circular", "allN", "nfree",
                        c("linear", "linear", "circular"), character())
for(chromosomes in chromosome_sets) {
  expected <- trinucleotideFrequency(sequences[chromosomes], simplify.as = "collapsed")
  stopifnot(identical(reference_counts_for_chromosomes(counts, chromosomes), expected))
}
# Circular flags must not introduce extra end-to-start trinucleotides.
stopifnot(sum(counts$circular) == 2L)

assignment <- function(file, name) {
  expressions <- parse(file)
  matches <- vapply(expressions, function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                       identical(x[[2]], as.name(name)), logical(1))
  stopifnot(any(matches))
  expressions[[which(matches)[[1]]]]
}
env <- new.env(parent = globalenv())
env$test_reference <- sequences
env$getSeq <- function(x, names) x[names]
eval(assignment(file.path(repo, "bin", "calculateBurdens.R"), "get_genome_reftnc"), env)
for(chromosomes in chromosome_sets) {
  env$reference_summary <- NULL
  original <- env$get_genome_reftnc("test_reference", chromosomes)
  env$reference_summary <- list(trinucleotide_counts = counts)
  cached <- env$get_genome_reftnc("test_reference", chromosomes)
  stopifnot(identical(original, cached))
}
cat("PASS: reference N ranges/Seqinfo, empty/N-free/all-N references, linear/circular boundaries, exact count types/order and burden tables\n")

# Compare the production calls-loading expression against the previous join-
# then-filter order, including factor levels, duplicate rows, NA chromosomes,
# additional SBS/indel categories, and an empty chromosome group.
env$extractedCalls <- list(calls = tibble(
  seqnames = factor(c("chr1", "chrM", "chrM", "chr1", "chrM", NA, "chrM"), levels = c("chr1", "chrM", "unused")),
  call_type = factor(c("SBS", "SBS", "SBS", "insertion", "MDB2", "SBS", "deletion")),
  call_class = factor(c("SBS", "SBS", "SBS", "indel", "MDB", "SBS", "indel")),
  SBSindel_call_type = factor(c("mutation", "mutation", "mutation", "mutation", "match", "mutation", "mismatch-ss")),
  value = c(1L, 2L, 2L, 3L, 4L, 5L, 6L)))
env$call_types_toanalyze <- env$extractedCalls$calls %>%
  filter(call_class == "SBS") %>% distinct(call_type, call_class, SBSindel_call_type)
for(chromosomes in list("chrM", "chr1", character(), c("chr1", "chrM"))) {
  env$chroms_toanalyze <- chromosomes
  original <- env$extractedCalls$calls %>%
    left_join(env$call_types_toanalyze %>% mutate(call_toanalyze = TRUE),
              by = join_by(call_type, call_class, SBSindel_call_type)) %>%
    mutate(call_toanalyze = replace_na(call_toanalyze, FALSE)) %>%
    filter(seqnames %in% chromosomes, call_toanalyze | call_class %in% c("SBS", "indel"))
  eval(assignment(file.path(repo, "bin", "filterCalls.R"), "calls"), env)
  stopifnot(identical(env$calls, original))
}
cat("PASS: early chromosome filtering preserves calls, row order, duplicates, factors, NA and empty-group behavior\n")
