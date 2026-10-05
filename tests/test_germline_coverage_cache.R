#!/usr/bin/env -S Rscript --vanilla
# Run inside the pinned container with wiggletools and wigToBigWig on PATH.
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(qs2))

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2L) stop("Usage: test_germline_coverage_cache.R HELPER OUTPUT_DIR")
helper <- normalizePath(args[[1]], mustWork=TRUE)
dir.create(args[[2]], recursive=TRUE, showWarnings=FALSE)
setwd(args[[2]])
wiggletools <- Sys.which("wiggletools")
wigToBigWig <- Sys.which("wigToBigWig")
stopifnot(nzchar(wiggletools), nzchar(wigToBigWig))

reference <- Seqinfo(c("chr1", "chr2", "chrEmpty"), c(30L, 20L, 8L),
                     isCircular=c(FALSE, TRUE, FALSE), genome="fixture")
coverage <- GRanges(c("chr1", "chr1", "chr1", "chr2", "chr2"),
                    IRanges(c(1L, 8L, 20L, 2L, 14L), c(4L, 13L, 28L, 7L, 20L)),
                    score=c(0, 15, 16, 1, 25), seqinfo=reference)
export(coverage, "coverage.bw", format="BigWig")
writeLines(c("chr1\t30\t0\t30\t31", "chr2\t20\t0\t20\t21", "chrEmpty\t8\t0\t8\t9"), "reference.fa.fai")
query <- GRanges(c("chr1", "chr1", "chr2", "chrEmpty"),
                 IRanges(c(4L, 13L, 7L, 1L), c(5L, 14L, 8L, 8L)), seqinfo=reference)

for(threshold in c(0, 1, 15, 16, 26)) {
  # Literal legacy shell/import/Seqinfo/select sequence on a tiny sparse track.
  commands <- c("set -euo pipefail",
    "awk '{print $1 \"\t0\t\" $2}' reference.fa.fai | sort -k1,1 -k2,2n > legacy.chromsizes.bed",
    paste(shQuote(wiggletools), "lt", threshold,
          "trim legacy.chromsizes.bed fillIn legacy.chromsizes.bed coverage.bw |",
          shQuote(wigToBigWig), "stdin <(cut -f 1,2 reference.fa.fai) legacy.bw"))
  legacy_status <- system2("/bin/bash", args="-s", input=commands)
  output <- paste0("threshold", threshold, ".qs2")
  status <- system2(file.path(R.home("bin"), "Rscript"), args=c("--vanilla", shQuote(helper),
    "--bigwig coverage.bw --fai reference.fa.fai --threshold", threshold,
    "--wiggletools", shQuote(wiggletools), "--wig_to_bigwig", shQuote(wigToBigWig),
    "--output", shQuote(output)))
  if(threshold == 0) {
    # Existing wigToBigWig rejects the empty comparison stream. Preserve that
    # failure rather than silently changing the pipeline's empty-filter policy.
    stopifnot(legacy_status != 0L, status != 0L, !file.exists(output))
    next
  }
  stopifnot(legacy_status == 0L, status == 0L)
  expected <- import("legacy.bw", format="bigWig")
  seqlevels(expected) <- seqlevels(reference)
  seqinfo(expected) <- reference
  expected <- plyranges::select(expected, -score)
  actual <- qs_read(output)
  seqlevels(actual) <- seqlevels(reference)
  seqinfo(actual) <- reference
  stopifnot(identical(actual, expected), identical(sum(width(actual)), sum(width(expected))),
            identical(overlapsAny(query, actual), overlapsAny(query, expected)))
}
cat("Germline coverage cache matches legacy GRanges, Seqinfo, widths and overlaps at four thresholds; both paths reject the legacy empty-stream edge.\n")
