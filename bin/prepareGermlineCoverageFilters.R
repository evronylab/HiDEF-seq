#!/usr/bin/env -S Rscript --vanilla

# Prepare the same whole-genome low-coverage intervals previously rebuilt in
# every filterCalls task. Reference Seqinfo is assigned by the consumer, exactly
# where the reference is available; no chromosome restriction is applied.
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(qs2))

source(Sys.which("sharedFunctions.R"))

options(warn=2)
option_list <- list(
  make_option(c("-b", "--bigwig"), type="character"),
  make_option(c("-f", "--fai"), type="character"),
  make_option(c("-t", "--threshold"), type="double"),
  make_option(c("-w", "--wiggletools"), type="character"),
  make_option(c("-u", "--wig_to_bigwig"), type="character"),
  make_option(c("-o", "--output"), type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))
required <- c("bigwig", "fai", "threshold", "wiggletools", "wig_to_bigwig", "output")
if(any(vapply(required, function(name) is.null(opt[[name]]), logical(1)))) {
  stop("All input, tool, threshold and output options are required", call.=FALSE)
}
intervals <- prepare_germline_coverage_filter(
  opt$bigwig, opt$fai, opt$threshold, opt$wiggletools, opt$wig_to_bigwig)
qs2::qs_save(intervals, opt$output)
