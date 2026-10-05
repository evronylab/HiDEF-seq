#!/usr/bin/env Rscript
#Run baseline and candidate as separate compute-node processes with identical
#parameters, then compare their single output QS2 files with identical(qs_read()).
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3L || length(args) > 6L) {
  stop("Usage: benchmark_coverage_accumulator.R REPO baseline|candidate NEW_PREFIX [INTERVALS=1000000] [CATEGORIES=8] [CHUNKS=6]")
}
suppressPackageStartupMessages({library(GenomicRanges); library(plyranges); library(tidyverse); library(qs2)})
options(warn = 2)
repo <- args[[1]]
mode <- match.arg(args[[2]], c("baseline", "candidate"))
prefix <- args[[3]]
intervals <- if(length(args) >= 4L) as.integer(args[[4]]) else 1000000L
categories <- if(length(args) >= 5L) as.integer(args[[5]]) else 8L
chunks <- if(length(args) >= 6L) as.integer(args[[6]]) else 6L
stopifnot(intervals > 0L, categories > 0L, chunks > 0L)
paths <- paste0(prefix, c(".qs2", ".metrics.tsv", ".chunks.tsv", ".session.txt"))
if(any(file.exists(paths))) stop("Refusing to overwrite benchmark outputs")
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
expressions <- parse(file.path(repo, "bin", "calculateBurdens.R"))
for(name in c("sum_RleList", "bc_orientation_is_asymmetric", "validate_bam.gr.filtertrack",
              "calc_duplex_coverage", "accumulate_bam.gr.filtertracks", "filtertrack_coverage_result")) {
  hits <- vapply(expressions, function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                   identical(x[[2]], as.name(name)), logical(1))
  stopifnot(sum(hits) == 1L)
  eval(expressions[[which(hits)]])
}
sys.source(file.path(repo, "tests", "fixtures", "sum_filtertrack_coverage_before_accumulator.R"), envir = globalenv())
stride <- chunks * 12 + 32
chrom_length <- as.double(intervals) * stride + 32
stopifnot(chrom_length < .Machine$integer.max, as.double(intervals) * chunks < .Machine$integer.max)
si <- Seqinfo("synthetic", as.integer(chrom_length))
keys <- tibble(call_type = factor(paste0("category", seq_len(categories))), call_class = factor("SBS"),
               SBSindel_call_type = factor("mutation"), filtergroup = factor("synthetic"))
rss <- function(field) {
  line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
  if(!length(line)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line))
}
reset_hwm <- function() {
  tryCatch({con <- file("/proc/self/clear_refs", "w"); on.exit(close(con)); writeLines("5", con); TRUE},
           warning = function(w) FALSE, error = function(e) FALSE)
}
state <- if(mode == "candidate") new.env(parent = emptyenv()) else NULL
invisible(gc())
reset <- reset_hwm()
before_rss <- rss("VmRSS")
before <- proc.time()
chunk_metrics <- list()
for(chunk in seq_len(chunks)) {
  #Interleaved, non-overlapping chunks create steadily growing Rle run counts.
  #Each category independently computes coverage from the same input intervals,
  #as real call categories often cover most of the same molecules.
  starts <- as.integer(seq.int(from = 1 + (chunk - 1) * 12, by = stride, length.out = intervals))
  gr <- GRanges(rep("synthetic", 2 * intervals), IRanges(rep(starts, 2), width = 5L),
                strand = rep(c("+", "-"), each = intervals), seqinfo = si)
  gr$run_id <- factor(rep("synthetic", length(gr)))
  gr$zm <- rep((chunk - 1L) * intervals + seq_len(intervals), 2)
  gr$bc_orientation <- factor(rep("A-A", length(gr)))
  incoming <- mutate(keys, bam.gr.filtertrack = rep(list(gr), nrow(keys)))
  if(mode == "candidate") {
    accumulate_bam.gr.filtertracks(state, incoming)
  } else {
    state <- if(chunk == 1L) sum_bam.gr.filtertracks(incoming) else sum_bam.gr.filtertracks(state, incoming)
  }
  rm(gr, incoming, starts)
  #Match the existing production cleanup after each input chunk.
  invisible(gc())
  timing <- proc.time() - before
  chunk_metrics[[chunk]] <- data.frame(chunk = chunk, actual_cpu_seconds = unname(sum(timing[c(1, 2, 4, 5)])),
                                       wall_seconds = unname(timing[[3]]), rss_kib = rss("VmRSS"), peak_rss_kib = rss("VmHWM"))
  write.table(bind_rows(chunk_metrics), paths[[3]], sep = "\t", row.names = FALSE, quote = FALSE)
  cat(mode, "chunk", chunk, "complete\n")
}
result <- if(mode == "candidate") filtertrack_coverage_result(state) else state
rm(state)
timing <- proc.time() - before
metrics <- data.frame(mode = mode, intervals = intervals, categories = categories, chunks = chunks,
                       user_seconds = unname(timing[[1]] + timing[[4]]), system_seconds = unname(timing[[2]] + timing[[5]]),
                       actual_cpu_seconds = unname(sum(timing[c(1, 2, 4, 5)])), wall_seconds = unname(timing[[3]]),
                       rss_before_kib = before_rss, rss_after_kib = rss("VmRSS"), peak_rss_kib = rss("VmHWM"),
                       peak_is_phase_specific = reset)
write.table(metrics, paths[[2]], sep = "\t", row.names = FALSE, quote = FALSE)
#Write one complete final accumulator outside the measured accumulation phase.
qs_save(result, paths[[1]])
capture.output(sessionInfo(), file = paths[[4]])
