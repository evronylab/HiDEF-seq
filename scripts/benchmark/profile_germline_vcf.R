#!/usr/bin/env Rscript
# Profile actual export expressions; keep the original final object resident.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(5L, 6L))
mode <- if(length(args) == 6L) args[[6]] else "profile"
stopifnot(mode %in% c("profile", "default", "100000"))
source_dir <- args[[1]]
input <- args[[2]]
chromgroup <- args[[3]]
filtergroup <- args[[4]]
prefix <- args[[5]]
source(file.path(source_dir, "germline_formatting_helpers.R"))
parsed <- parse(file.path(source_dir, "outputResults.R"))
writer <- assignment(parsed, "write_vcf_from_calls")
shared <- parse(file.path(source_dir, "sharedFunctions.R"))
for(name in c("normalize_indels_for_vcf", "get_bsgenome_name", "reference_cache_dir")) eval(assignment(shared, name))
suppressPackageStartupMessages({
  library(GenomicAlignments); library(GenomicRanges); library(BSgenome)
  library(plyr); library(plyranges); library(configr); library(rtracklayer)
  library(VariantAnnotation); library(qs2); library(tidyverse)
})
options(warn = 2)
eval(writer)
# Absolute event times allow inclusive phase CPU differences without retaining
# results or forcing collection. VmHWM is never reset during this process.
events <- paste0(prefix, ".events.tsv")
if(file.exists(events)) stop("Profile output already exists")
writeLines("event\tuser_seconds\tsystem_seconds\twall_seconds\trss_kib\thwm_kib", events)
profile_event <- function(label) {
  timing <- proc.time()
  cat(label, sum(timing[c(1, 4)]), sum(timing[c(2, 5)]), timing[[3]],
      rss("VmRSS"), rss("VmHWM"), sep = "\t", file = events, append = TRUE)
  cat("\n", file = events, append = TRUE)
}
# Insert BEFORE markers, preserving the original last expression and return.
# on.exit supplies an exit marker without holding the function result.
instrument <- function(fn, label) {
  original <- as.list(body(fn))[-1L]
  entry <- substitute(on.exit(.GlobalEnv$profile_event(LABEL), add = TRUE),
                      list(LABEL = paste0(label, ":exit")))
  parts <- list(as.name("{"), entry)
  for(i in seq_along(original)) {
    description <- substr(paste(deparse(original[[i]]), collapse = " "), 1L, 100L)
    marker <- substitute(.GlobalEnv$profile_event(LABEL),
      list(LABEL = paste0(label, ":", i, ":", gsub("[\t\r\n]", " ", description))))
    parts <- c(parts, list(marker, original[[i]]))
  }
  body(fn) <- as.call(parts)
  fn
}
namespace <- asNamespace("VariantAnnotation")
internal_names <- c(".chunkIndex", ".makeVcfInfo", ".makeVcfMatrix", ".makeVcfGeno")
# Save installed code BEFORE instrumentation; this records native chunk rules.
capture.output({
  print(packageVersion("VariantAnnotation"))
  print(methods::getMethod("writeVcf", c("VCF", "character")))
  print(methods::getMethod("writeVcf", c("VCF", "connection")))
  for(name in internal_names) {cat("\n", name, "\n"); print(get(name, namespace))}
}, file = paste0(prefix, ".installed-writer.R.txt"))
# Benchmark arms differ by this single existing API argument only.
with_native_buffer <- function(fn, rows) {
  expressions <- as.list(body(fn))
  last <- expressions[[length(expressions)]]
  stopifnot(is.call(last), identical(last[[1]], quote(VariantAnnotation::writeVcf)))
  arguments <- as.list(last)
  stopifnot(sum(names(arguments) == "nchunk") <= 1L,
            is.null(arguments$nchunk) || identical(arguments$nchunk, 100000L))
  if(!is.null(rows)) stopifnot(is.integer(rows), length(rows) == 1L, !is.na(rows), rows > 0L)
  # Explicit legacy normalization: remove only the known configured argument.
  # Assigning NULL to a list removes the element rather than passing NULL to R.
  arguments$nchunk <- rows
  expressions[[length(expressions)]] <- as.call(arguments)
  body(fn) <- as.call(expressions)
  fn
}
if(mode == "default") write_vcf_from_calls <- with_native_buffer(write_vcf_from_calls, NULL)
if(mode == "100000") write_vcf_from_calls <- with_native_buffer(write_vcf_from_calls, 100000L)
effective_nchunk <- as.list(body(write_vcf_from_calls))[[length(body(write_vcf_from_calls))]]$nchunk
capture.output(print(write_vcf_from_calls), file = paste0(prefix, ".effective-writer.R.txt"))
if(mode == "profile") {
  for(name in internal_names[-1L]) {
    fn <- get(name, namespace)
    assignInNamespace(name, instrument(fn, paste0("VariantAnnotation", name)), "VariantAnnotation")
  }
  write_vcf_from_calls <- instrument(write_vcf_from_calls, "write_vcf_from_calls")
  normalize_indels_for_vcf <- instrument(normalize_indels_for_vcf, "normalize_indels_for_vcf")
}
if(input == "--preflight") {
  fn <- function(x) { y <- x + 1L; if(x < 0L) return(y); y * 2L }
  changed <- instrument(fn, "fixture")
  stopifnot(identical(fn(3L), changed(3L)), identical(fn(-2L), changed(-2L)))
  if(mode != "profile") {
    source(file.path(source_dir, "test_vcf_native_buffer.R"))
  }
  cat("AST, installed writer methods, instrumentation returns, and requested fixtures passed\n")
  quit(status = 0L)
}
profile_event("load:begin")
resident <- qs_read(input)
profile_event("load:end")
yaml.config <- resident$yaml.config
BSgenome_name <- get_bsgenome_name(yaml.config)
suppressPackageStartupMessages(library(BSgenome_name, character.only = TRUE,
                                       lib.loc = reference_cache_dir(yaml.config)))
formatted <- resident$germlineVariantCalls_for_tsv %>%
  filter(.data$chromgroup == .env$chromgroup, .data$filtergroup == .env$filtergroup)
stopifnot(nrow(formatted) > 0L)
writeLines(c(paste("selected_rows", nrow(formatted)),
             paste("selected_columns", ncol(formatted)),
             "All original final QS components remain resident, including the complete formatted table.",
             "This conservatively retains later-group formatted/final results absent at the first real export.",
             "No explicit GC or HWM reset is performed.",
             paste("writer_mode", mode),
             paste("native_nchunk", if(is.null(effective_nchunk)) "default" else effective_nchunk),
             paste("expression_instrumentation", mode == "profile")),
           paste0(prefix, ".payload.txt"))
profile_event("export:begin")
formatted %>% rename(start = start_refspace, end = end_refspace) %>%
  normalize_indels_for_vcf(BSgenome_name = BSgenome_name) %>%
  rename(start_refspace = start) %>%
  write_vcf_from_calls(BSgenome_name = BSgenome_name, out_vcf = paste0(prefix, ".vcf"))
profile_event("export:end")
stopifnot(file.exists(paste0(prefix, ".vcf.bgz")), file.exists(paste0(prefix, ".vcf.bgz.tbi")))
capture.output(sessionInfo(), file = paste0(prefix, ".session.txt"))
