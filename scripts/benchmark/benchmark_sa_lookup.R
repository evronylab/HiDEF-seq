#!/usr/bin/env Rscript
# Real own-strand queries from an existing extraction QS. This isolates lookup
# cost; it does not claim to measure complete extraction or opposite-strand CIGAR mapping.
args <- commandArgs(TRUE)
if(length(args) < 3L || length(args) > 5L) stop("Usage: benchmark_sa_lookup.R REPO EXTRACT.qs2 NEW_PREFIX [REPEATS=3] [REFERENCE_EXTRACT.R]")
repo <- args[[1]]
prefix <- args[[3]]
repeats <- if(length(args) >= 4L) as.integer(args[[4]]) else 3L
reference_source <- if(length(args) == 5L) args[[5]] else file.path(repo, "bin/extractCalls.R")
stopifnot(repeats > 0L)
if(any(file.exists(paste0(prefix, c(".metrics.tsv", ".queries.tsv", ".session.txt", ".code"))))) stop("Output exists")
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
code_dir <- paste0(prefix, ".code")
dir.create(code_dir)
file.copy(reference_source, file.path(code_dir, "reference_extractCalls.R"))
file.copy(file.path(repo, "scripts/benchmark/sa_interval_lookup.R"), code_dir)
file.copy(file.path(repo, "scripts/benchmark/benchmark_sa_lookup.R"), code_dir)
hashes <- vapply(list.files(code_dir, full.names=TRUE), function(path) {
  system2("sha256sum", shQuote(path), stdout=TRUE)
}, character(1))
writeLines(hashes, file.path(code_dir, "source.sha256"))
suppressPackageStartupMessages({library(qs2); library(GenomicRanges); library(data.table)})
options(warn = 2)
source(file.path(code_dir, "sa_interval_lookup.R"))
expressions <- parse(file.path(code_dir, "reference_extractCalls.R"))
hits <- vapply(expressions, function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                 identical(x[[2]], as.name("subset_tag_positions")), logical(1))
stopifnot(sum(hits) == 1L)
eval(expressions[[which(hits)]])
object <- qs_read(args[[2]])
bam <- as.data.table(object$bam)[, .(run_id, zm, strand, sa)]
calls <- as.data.table(object$calls)[, .(run_id, zm, strand, call_class, call_type,
                                      start_queryspace, end_queryspace, sa)]
rm(object)
invisible(gc())
stopifnot(!anyDuplicated(bam[, .(run_id, zm, strand)]))
bam[, tag_row := .I]
calls <- merge(calls, bam[, .(run_id, zm, strand, tag_row)], by=c("run_id", "zm", "strand"), sort=FALSE)
stopifnot(!anyNA(calls$tag_row))
# SBS/MDB rows are already one base each; insertion ranges include the inserted
# bases. Deletions retain zero-width start>end labels, and extraction queries
# the two flanks by swapping start/end exactly as its forqualdata block does.
calls[, positions := Map(function(a, b) seq.int(min(a, b), max(a, b)), start_queryspace, end_queryspace)]
queries <- calls[, .(positions=list(unlist(positions, use.names=FALSE)),
                     expected=list(unlist(sa, use.names=FALSE))), by=.(tag_row, call_class, call_type)]
tags <- bam$sa[queries$tag_row]
positions <- queries$positions
expected <- queries$expected
rm(bam, calls)
invisible(gc())
for(i in seq_along(tags)) {
  baseline <- subset_tag_positions(tags[[i]], positions[[i]])
  stopifnot(identical(baseline, expected[[i]]),
            identical(subset_tag_positions_interval(tags[[i]], positions[[i]]), expected[[i]]))
}
queries[, `:=`(tag_row=NULL, positions=NULL, expected=NULL)]
queries[, `:=`(tag_length=vapply(tags, length, numeric(1)),
               run_count=vapply(tags, function(x) length(runValue(x)), integer(1)),
               query_count=lengths(positions))]
write.table(queries, paste0(prefix, ".queries.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
lookup <- function(mode) {
  if(mode == "dense_selected") {
    dense <- lapply(tags, as.vector)
    return(Map(function(tag, index) tag[index], dense, positions))
  }
  Map(if(mode == "source_lookup") subset_tag_positions else subset_tag_positions_interval, tags, positions)
}
metrics <- list()
for(iteration in seq_len(repeats)) {
  modes <- c("source_lookup", "interval", "dense_selected")
  if(iteration %% 2L == 0L) modes <- rev(modes)
  for(mode in modes) {
    invisible(gc())
    before <- proc.time()
    result <- lookup(mode)
    elapsed <- proc.time() - before
    stopifnot(identical(result, expected))
    metrics[[length(metrics) + 1L]] <- data.frame(iteration=iteration, mode=mode,
      groups=length(tags), positions=sum(lengths(positions)),
      actual_cpu_seconds=unname(sum(elapsed[c(1, 2, 4, 5)])), wall_seconds=unname(elapsed[[3]]))
    write.table(do.call(rbind, metrics), paste0(prefix, ".metrics.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
    rm(result)
  }
}
capture.output(sessionInfo(), file=paste0(prefix, ".session.txt"))
cat("PASS: real query reconstruction and every measured lookup matched saved own-strand sa values exactly\n")
