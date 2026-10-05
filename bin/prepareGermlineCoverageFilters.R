#!/usr/bin/env -S Rscript --vanilla

# Prepare the same whole-genome low-coverage intervals previously rebuilt in
# every filterCalls task. Reference Seqinfo is assigned by the consumer, exactly
# where it was assigned in the legacy path; no chromosome restriction is applied.
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(qs2))

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
if(length(opt$threshold) != 1L || !is.finite(opt$threshold)) {
  stop("threshold must be one finite number", call.=FALSE)
}

tmpchromsizes <- tempfile(tmpdir=getwd(), pattern=".germline-coverage-", fileext=".bed")
tmpbw <- tempfile(tmpdir=getwd(), pattern=".germline-coverage-", fileext=".bw")

commands <- c(
  paste("awk '{print $1 \"\t0\t\" $2}'", shQuote(opt$fai),
        "| sort -k1,1 -k2,2n >", shQuote(tmpchromsizes)),
  paste(shQuote(opt$wiggletools), "lt", opt$threshold,
        "trim", shQuote(tmpchromsizes), "fillIn", shQuote(tmpchromsizes),
        shQuote(opt$bigwig), "|", shQuote(opt$wig_to_bigwig),
        "stdin <(cut -f 1,2", shQuote(opt$fai), ")", shQuote(tmpbw))
)
status <- system2("/bin/bash", args="-s", input=c("set -euo pipefail", commands))
if(status != 0L) stop("Failed to prepare germline coverage filter", call.=FALSE)

intervals <- rtracklayer::import(tmpbw, format="bigWig")
intervals <- plyranges::select(intervals, -score)
qs2::qs_save(intervals, opt$output)
invisible(file.remove(tmpchromsizes, tmpbw))
