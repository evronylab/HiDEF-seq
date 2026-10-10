#!/usr/bin/env Rscript
#Run in the pipeline container on a compute node, with optional repository root.
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(plyranges)
  library(tidyverse)
})
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
repo <- if(length(args)) args[[1]] else "."
reference <- new.env(parent = globalenv())
sys.source(file.path(repo, "tests", "fixtures", "burden_native_reference_ea8d06f.R"), reference)
gr_1bp_cov <- reference$gr_1bp_cov
expressions <- parse(file.path(repo, "bin", "calculateBurdens.R"))
assignment <- function(name) {
  hits <- vapply(expressions, function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                   identical(x[[2]], as.name(name)), logical(1))
  stopifnot(any(hits))
  expressions[[which(hits)[[1]]]]
}
for(name in c("validate_bam.gr.filtertrack", "calc_duplex_coverage", "sum_RleList",
              "make_sensitivity_coverage_queries", "duplex_coverage_ranges", "sensitivity_site_counts",
              "sum_filtertrack_sensitivity_coverage")) {
  eval(assignment(name))
}
si <- Seqinfo(c("chr1", "chr2"), c(100L, 100L))
variants <- tribble(
  ~seqnames, ~start, ~end, ~call_class, ~call_type,
  "chr1", 10L, 12L, "indel", "deletion",
  "chr1", 11L, 10L, "indel", "insertion",
  "chr1", 5L, 5L, "SBS", "SBS",
  "chr2", 20L, 20L, "SBS", "SBS"
)
original_variants <- variants
queries <- make_sensitivity_coverage_queries(variants, si)
stopifnot(identical(variants, original_variants))
keys <- tribble(
  ~call_type, ~call_class, ~SBSindel_call_type, ~filtergroup,
  "deletion", "indel", "mutation", "test",
  "deletion", "indel", "mismatch-ss", "test",
  "insertion", "indel", "mutation", "test",
  "SBS", "SBS", "mutation", "test",
  "SBS", "SBS", "mismatch-ss", "test",
  "MDB2", "MDB", "mismatch-os", "test"
) %>% mutate(across(everything(), factor))
paired_ranges <- function(starts, ends = starts) {
  gr <- GRanges(rep("chr1", 2 * length(starts)), IRanges(rep(starts, 2), rep(ends, 2)),
                strand = rep(c("+", "-"), each = length(starts)), seqinfo = si)
  gr$run_id <- rep("run", length(gr))
  gr$zm <- rep(seq_along(starts), 2)
  gr$bc_orientation <- factor(rep("A-A", length(gr)))
  #A queried chromosome absent from this coverage object must get numeric zero.
  keepSeqlevels(gr, "chr1")
}
left <- paired_ranges(c(9L, 9L, 5L, 10L))
right <- paired_ranges(c(13L, 13L, 5L, 11L))
chunks <- list(
  mutate(keys, bam.gr.filtertrack = list(left, left[0], left, left, left[0], left)),
  mutate(keys, bam.gr.filtertrack = list(right, right, right, right, right, right))
)
legacy <- map2(chunks[[1]]$bam.gr.filtertrack, chunks[[2]]$bam.gr.filtertrack,
               ~ sum_RleList(coverage(.x) %/% 2L, coverage(.y) %/% 2L))
sparse <- sum_filtertrack_sensitivity_coverage(chunks[[1]], queries)
sparse <- sum_filtertrack_sensitivity_coverage(chunks[[2]], queries, previous = sparse)
stopifnot(identical(select(sparse, all_of(names(keys))), keys))
for(i in seq_len(nrow(keys))) {
  idx <- match(keys$call_type[i], queries$call_type)
  expected_start <- if(is.na(idx)) numeric() else gr_1bp_cov(queries$query_start[[idx]], legacy[[i]])
  expected_end <- if(is.na(idx)) numeric() else gr_1bp_cov(queries$query_end[[idx]], legacy[[i]])
  stopifnot(identical(sparse$coverage_start[[i]], expected_start),
            identical(sparse$coverage_end[[i]], expected_end),
            identical(pmin(sparse$coverage_start[[i]], sparse$coverage_end[[i]]), pmin(expected_start, expected_end)))
}
#Left=2/right=0 then left=0/right=2 must produce two, not zero.
stopifnot(identical(pmin(sparse$coverage_start[[1]], sparse$coverage_end[[1]]), 2))
stopifnot(identical(pmin(sparse$coverage_start[[2]], sparse$coverage_end[[2]]), 0))
stopifnot(identical(pmin(sparse$coverage_start[[3]], sparse$coverage_end[[3]]), 1))
stopifnot(identical(pmin(sparse$coverage_start[[4]], sparse$coverage_end[[4]]), c(2, 0)))

#Check the final nested annotations and scalar count column classes, not just
#values, against the former aggregate-coverage query method.
annotated_variants <- variants %>%
  mutate(start = if_else(call_class == "indel", start - 1, start),
         end = if_else(call_class == "indel", end + 1, end), num_zm_detected = 1L)
sparse_annotated <- sparse %>% nest_join(annotated_variants, by = "call_type", name = "data") %>%
  mutate(data = pmap(list(data, coverage_start, coverage_end), ~ mutate(..1, duplex_coverage = pmin(..2, ..3)))) %>%
  select(-coverage_start, -coverage_end)
legacy_annotated <- keys %>% mutate(cov = legacy) %>% nest_join(annotated_variants, by = "call_type", name = "data") %>%
  mutate(data = map2(data, cov, function(x, cov) {
    gr <- makeGRangesFromDataFrame(x, seqinfo = si)
    mutate(x, duplex_coverage = pmin(gr_1bp_cov(resize(gr, 1, fix = "start"), cov),
                                    gr_1bp_cov(resize(gr, 1, fix = "end"), cov)))
  })) %>% select(-cov)
stopifnot(identical(sparse_annotated, legacy_annotated))
counts <- function(x) mutate(x, num_zm_detected = map_int(data, ~ sum(.x$num_zm_detected)),
                             duplex_coverage = map_int(data, ~ sum(.x$duplex_coverage))) %>% select(-data)
stopifnot(identical(counts(sparse_annotated), counts(legacy_annotated)))
for(empty_queries in list(NULL, make_sensitivity_coverage_queries(variants[0,], si))) {
  empty <- sum_filtertrack_sensitivity_coverage(chunks[[1]], empty_queries)
  empty <- sum_filtertrack_sensitivity_coverage(chunks[[2]], empty_queries, empty)
  stopifnot(all(lengths(empty$coverage_start) == 0L), all(lengths(empty$coverage_end) == 0L))
}
expect_error <- function(expr, pattern) {
  message <- tryCatch({force(expr); NA_character_}, error = function(e) conditionMessage(e))
  stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE))
}
bad <- paired_ranges(80L)
start(bad)[2] <- 81L
for(q in list(NULL, queries)) {
  expect_error(sum_filtertrack_sensitivity_coverage(mutate(keys[1,], bam.gr.filtertrack = list(bad)), q),
               "Mismatched plus and minus strand ranges")
}
bad <- paired_ranges(80L)
seqlevels(bad) <- c("chr1", "chr2")
seqlengths(bad) <- c(100L, 100L)
seqnames(bad) <- c("chr1", "chr2")
for(q in list(NULL, queries)) {
  expect_error(sum_filtertrack_sensitivity_coverage(mutate(keys[1,], bam.gr.filtertrack = list(bad)), q),
               "Non-even strand coverage")
}
cat("PASS: sparse/aggregate sensitivity, cross-chunk flank minimum, missing chromosomes, empty/disabled sensitivity, invariant checks, types\n")

#Terminal indels retain the original out-of-bounds failure rather than silently
#clipping a queried flank. Suppress the same GRanges construction warnings in
#both paths here so the actual coordinate-0/length+1 coverage lookup is tested.
capture_lookup <- function(expr) {
  tryCatch(list(value = force(expr)), error = function(e) list(error = conditionMessage(e)))
}
for(pos in c(1L, 100L)) {
  terminal <- tibble(seqnames = "chr1", start = pos, end = pos,
                     call_class = "indel", call_type = "deletion")
  terminal_queries <- suppressWarnings(make_sensitivity_coverage_queries(terminal, si))
  legacy_range <- suppressWarnings(terminal %>%
    mutate(start = start - 1, end = end + 1) %>%
    makeGRangesFromDataFrame(seqinfo = si))
  side <- if(pos == 1L) "start" else "end"
  candidate_query <- terminal_queries[[paste0("query_", side)]][[1]]
  legacy_query <- suppressWarnings(resize(legacy_range, 1, fix = side))
  stopifnot(identical(candidate_query, legacy_query),
            identical(start(candidate_query), if(pos == 1L) 0L else 101L))
  baseline_failure <- capture_lookup(gr_1bp_cov(legacy_query, legacy[[1]]))
  candidate_failure <- capture_lookup(gr_1bp_cov(candidate_query, legacy[[1]]))
  stopifnot(!is.null(baseline_failure$error), identical(candidate_failure, baseline_failure))
}
cat("PASS: unchanged coordinate-0 and seqlength+1 sensitivity flank failures\n")


#Exercise the moved preparation function with in-memory I/O stand-ins. Its
#global quantile must include the excluded chromosomes and both call classes.
fixture <- tibble(
  seqnames = c(rep("chr1", 8), "chrX", "chrM", "chr2", "chr2"),
  start = seq_len(12), end = seq_len(12), ref_plus_strand = "A", alt_plus_strand = "T",
  call_class = "SBS", call_type = "SBS", SBSindel_call_type = "mutation",
  Depth = c(10, 20, 30, 40, 160, 180, 200, 220, 300, 400, 500, 600),
  VAF = 0.5, GT = "1/1"
) %>% mutate(
  GQ = Depth, QUAL = Depth,
  call_class = if_else(start %in% c(3L, 5L), "indel", call_class),
  call_type = if_else(call_class == "indel", "deletion", call_type),
  ref_plus_strand = if_else(call_class == "indel", "AA", ref_plus_strand),
  alt_plus_strand = if_else(call_class == "indel", "A", alt_plus_strand),
  end = if_else(call_class == "indel", start + 1L, end)
)
variant_keys <- c("seqnames", "start", "end", "ref_plus_strand", "alt_plus_strand",
                  "call_class", "call_type", "SBSindel_call_type")
config <- list(cache_dir = "/unused", genome_fasta = "/unused/ref.fa", bcftools_bin = "/unused/bcftools",
               individuals = list(list(individual_id = "individual", germline_bam_file = "germline.bam",
                                       germline_vcf_files = list(list(file = "a"), list(file = "b")))))
thresholds <- list(sensitivity_vcf = "/unused/sensitivity.vcf", genotype = "homozygous",
                   SBS_min_Depth_quantile = 0.6, indel_min_Depth_quantile = 0.3,
                   SBS_min_GQ_quantile = 0.6, indel_min_GQ_quantile = 0.3,
                   SBS_min_QUAL_quantile = 0.6, indel_min_QUAL_quantile = 0.3,
                   SBS_min_VAF = 0, indel_min_VAF = 0, SBS_max_VAF = 1, indel_max_VAF = 1)
for(genotype in c("homozygous", "heterozygous")) {
  thresholds$genotype <- genotype
  fixture$GT <- if(genotype == "homozygous") "1/1" else "0/1"
  env <- new.env(parent = globalenv())
  env$germline <- bind_rows(mutate(fixture, germline_vcf_file = "a"), mutate(fixture, germline_vcf_file = "b"))
  env$allowed <- select(fixture, all_of(variant_keys))
  env$qs_read <- function(...) env$germline
  env$load_vcf <- function(...) env$allowed
  env$cache_file <- function(file, config) file
  eval(assignment("prepare_high_confidence_germline_vcf_variants"), env)
  prepare <- function() env$prepare_high_confidence_germline_vcf_variants(
    config, "individual", thresholds, "unused", c("chr1", "chrX", "chrM"), "chrX", "chrM")
  stopifnot(identical(prepare(), select(fixture[c(5, 8),], all_of(variant_keys))))
  env$germline <- filter(env$germline, !(germline_vcf_file == "b" & start %in% c(5L, 8L)))
  stopifnot(nrow(prepare()) == 0L) #Must be in every VCF, despite per-VCF quantiles.
  env$germline <- bind_rows(mutate(fixture, germline_vcf_file = "a"), mutate(fixture, germline_vcf_file = "b"))
  env$allowed <- filter(env$allowed, !start %in% c(5L, 8L))
  stopifnot(nrow(prepare()) == 0L) #External sensitivity VCF restriction.
}
cat("PASS: genome-wide per-VCF quantiles, both genotypes, all-VCF intersection, sex/mito and external-VCF filtering\n")
