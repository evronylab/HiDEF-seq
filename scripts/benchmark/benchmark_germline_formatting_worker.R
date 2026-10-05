#!/usr/bin/env Rscript
# Each command runs in a fresh process; prepare runs only once per benchmark.
args <- commandArgs(TRUE)
if(length(args) < 1L) stop("Expected prepare, format, or compare")
options(warn = 2)
script_file <- sub("^--file=", "", commandArgs(FALSE)[startsWith(commandArgs(FALSE), "--file=")][[1]])
source(file.path(dirname(normalizePath(script_file)), "germline_formatting_helpers.R"))
command <- args[[1]]
if(command == "prepare") {
  stopifnot(length(args) == 5L)
  # Resolve all formatter nodes before loading the large final QS object.
  for(path in args[4:5]) for(name in required_names) assignment(parse(path), name)
  assignment(parse(file.path(dirname(args[[5]]), "sharedFunctions.R")), "cache_file")
  if(file.exists(args[[3]])) stop("Packet already exists")
  suppressPackageStartupMessages(library(qs2))
  object <- qs_read(args[[2]])
  packet <- object[c("yaml.config", "germlineVariantCalls")]
  stopifnot(identical(names(packet), c("yaml.config", "germlineVariantCalls")))
  rm(object)
  invisible(gc())
  qs_save(packet, args[[3]])
  cat("Prepared raw germline packet:", nrow(packet$germlineVariantCalls), "rows\n")
} else if(command == "format") {
  stopifnot(length(args) == 6L)
  packet_path <- args[[2]]
  formatter_path <- args[[3]]
  candidate_shared_path <- args[[4]]
  prefix <- args[[5]]
  mode <- args[[6]]
  outputs <- paste0(prefix, c(".qs2", ".phase.tsv", ".session.txt", ".process_peak.tsv"))
  if(any(file.exists(outputs))) stop("Worker outputs already exist")
  parsed <- parse(formatter_path)
  resolved <- setNames(lapply(required_names, function(name) assignment(parsed, name)), required_names)
  eval(assignment(parse(candidate_shared_path), "cache_file"))
  suppressPackageStartupMessages({library(qs2); library(tidyverse)})
  packet <- qs_read(packet_path)
  yaml.config <- packet$yaml.config
  germlineVariantCalls <- packet$germlineVariantCalls
  rm(packet)
  invisible(gc())
  preload_peak <- rss("VmHWM")
  reset <- reset_hwm()
  before_rss <- rss("VmRSS")
  before <- proc.time()
  formatted <- format_all(resolved)
  timing <- proc.time() - before
  metrics <- data.frame(mode = mode, rows = nrow(formatted), columns = ncol(formatted),
    user_seconds = unname(timing[[1]] + timing[[4]]),
    system_seconds = unname(timing[[2]] + timing[[5]]),
    actual_cpu_seconds = unname(sum(timing[c(1, 2, 4, 5)])), wall_seconds = unname(timing[[3]]),
    rss_before_kib = before_rss, rss_after_kib = rss("VmRSS"), peak_rss_kib = rss("VmHWM"),
    peak_is_phase_specific = reset)
  write.table(metrics, outputs[[2]], sep = "\t", row.names = FALSE, quote = FALSE)
  qs_save(formatted, outputs[[1]])
  capture.output(sessionInfo(), file = outputs[[3]])
  # Reconstruct the full-process high-water mark across the phase reset; GNU
  # time/getrusage alone may lose the earlier load peak after clear_refs=5.
  write.table(data.frame(preformat_peak_rss_kib = preload_peak,
                         remaining_peak_rss_kib = rss("VmHWM")), outputs[[4]],
              sep = "\t", row.names = FALSE, quote = FALSE)
} else if(command == "compare") {
  stopifnot(length(args) == 4L, !file.exists(args[[4]]))
  suppressPackageStartupMessages(library(qs2))
  baseline <- qs_read(args[[2]])
  candidate <- qs_read(args[[3]])
  same <- identical(baseline, candidate)
  write.table(data.frame(identical = same, rows = nrow(baseline), columns = ncol(baseline)),
              args[[4]], sep = "\t", row.names = FALSE, quote = FALSE)
  if(!same) stop("Formatted tables differ (including values, types, attributes or order)")
} else stop("Unknown command")
