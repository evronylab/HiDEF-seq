#!/usr/bin/env Rscript
# Differential benchmark of reference preparation, run on a compute node.
# Usage: reference_summary.R CACHE_LIBRARY PACKAGE MODE OUTPUT.qs2 REPO
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L, args[[3]] %in% c("baseline", "candidate"))
suppressPackageStartupMessages({
  library(BSgenome)
  library(GenomicRanges)
  library(Biostrings)
  library(qs2)
  library(args[[2]], character.only = TRUE, lib.loc = args[[1]])
})
source(file.path(args[[5]], "bin", "referenceSummaryFunctions.R"))
genome <- get(args[[2]])
if (args[[3]] == "baseline") {
  n_ranges <- GenomicRanges::reduce(vmatchPattern("N", genome), ignore.strand = TRUE)
  counts <- trinucleotideFrequency(DNAStringSet(getSeq(genome)), simplify.as = "collapsed")
} else {
  n_ranges <- reference_n_ranges(genome)
  counts <- Reduce(`+`, reference_trinucleotide_counts(genome))
}
cat("N intervals:", length(n_ranges), "N bases:", sum(width(n_ranges)), "\n")
cat("Trinucleotide counts type:", typeof(counts), "\n")
qs_save(list(n_ranges = n_ranges, counts = counts), args[[4]])
