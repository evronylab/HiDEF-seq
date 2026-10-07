#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({
  library(GenomicRanges); library(vcfR); library(tidyverse); library(jsonlite)
})
options(warn = 2)
baseline <- new.env(parent = globalenv())
candidate <- new.env(parent = globalenv())
stopifnot(length(args) == 2L)
sys.source(file.path(args[[1]], 'bin/sharedFunctions.R'), envir = baseline)
sys.source(file.path(args[[1]], 'bin/sharedFunctions.R'), envir = candidate)
# Control: query the unchanged exact target intervals directly, disabling only
# coalescing. Normalization and parsing remain the same for both arms.
baseline$coalesce_vcf_regions <- function(regions, ...) regions
folder <- normalizePath(args[[2]])
setwd(folder)
fake_genome <- GRanges(c("chrA", "chrB"), IRanges(1L, 50000L),
                      seqinfo = Seqinfo(c("chrA", "chrB"), c(50000L, 50000L)))
rows <- function(chr, start, end = start) {
  tibble(seqnames = chr, start_refspace = start, end_refspace = end)
}
cases <- list(
  empty = rows(character(), integer()),
  nearby_sbs = rows(rep("chrA", 4), c(100, 101, 900, 901)),
  deletion_interior = rows(rep("chrA", 4), c(500, 502, 504, 510)),
  insertion = rows(rep("chrA", 3), c(699, 700, 701)),
  multiallelic_and_mnv = rows(rep("chrA", 7), c(900, 901, 902, 999, 1000, 1001, 1002)),
  overlap_duplicate_targets = rows(rep("chrA", 5), c(100, 100, 495, 501, 699), c(100, 102, 505, 510, 702)),
  reverse_order = rows(c("chrB", "chrA", "chrA"), c(100, 1000, 100)),
  factor_contigs = rows(factor(c("chrA", "chrB", "chrA"), levels=c("chrB", "chrA")), c(100, 100, 900)),
  no_records = rows("chrA", 40000),
  all_records = rows(c("chrA", "chrB"), c(1, 1), c(50000, 50000)),
  exact_gap_boundary = rows(rep("chrA", 4), c(100, 1100, 2101, 3101)),
  nearby_ranges = rows(rep("chrA", 3), c(95, 498, 899), c(105, 510, 1003)),
  long_deletion_separate_groups = rows(rep("chrA", 2), c(5100, 7200)),
  long_deletion_and_sbs = rows(rep("chrA", 4), c(5000, 5100, 7200, 7300)),
  same_coordinate_metadata = rows(rep("chrA", 4), c(100, 101, 900, 1100)),
  deeply_nested_ranges = rows(rep("chrA", 5), c(95, 96, 100, 499, 502), c(1000, 900, 100, 7100, 505))
)
capture <- function(fun, vcffile, regions) {
  tryCatch(list(ok = TRUE, value = fun(vcffile, regions,
    genome_fasta=file.path(folder, "reference.fa"), BSgenome_name="fake_genome",
    bcftools_bin=Sys.which("bcftools"))),
    error = function(e) list(ok = FALSE, message=conditionMessage(e)))
}
result <- list()
for(vcf in c("with_gt.vcf.gz", "without_gt.vcf.gz")) {
  for(name in names(cases)) {
    old <- capture(baseline$load_vcf, file.path(folder, vcf), cases[[name]])
    new <- capture(candidate$load_vcf, file.path(folder, vcf), cases[[name]])
    if(name %in% c('nearby_sbs', 'deletion_interior')) {
      stopifnot(old$ok, new$ok, length(old$value) > 0L, length(new$value) > 0L)
      if(name == 'nearby_sbs') stopifnot(sum(start(old$value) == 100L) == 2L)
      if(name == 'deletion_interior') stopifnot(any(old$value$call_type == 'deletion'))
    }
    exact <- identical(old$ok, new$ok) && (!old$ok || identical(old$value, new$value))
    result[[length(result)+1L]] <- list(vcf=vcf, case=name, exact=exact,
       baseline_success=old$ok, candidate_success=new$ok,
       baseline_error=old$message, candidate_error=new$message)
    if(!exact) saveRDS(list(baseline=old,candidate=new), paste0(vcf,".",name,".mismatch.rds"))
  }
}
# A gap of exactly 1000 joins; nested intervals must not shorten the group's
# reach; chromosome boundaries remain independent.
input <- rows(c('chrA','chrA','chrA','chrA','chrB'),
              c(1L,2L,1100L,2200L,1L), c(100L,3L,1100L,2200L,1L))
expected <- rows(c('chrA','chrA','chrB'), c(1L,2200L,1L), c(1100L,2200L,1L))
stopifnot(identical(candidate$coalesce_vcf_regions(input), expected))
stopifnot(all(map_lgl(result,'exact')),
          all(map_lgl(keep(result, ~ .x$case != 'no_records'), 'baseline_success')))
cat('PASS: 32 coalesced/direct VCF target comparisons, including GT/no-GT, long deletions, duplicate records, ordering and interval boundaries\n')
