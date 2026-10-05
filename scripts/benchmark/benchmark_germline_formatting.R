#!/usr/bin/env Rscript
#Compute-node benchmark; measures the full existing germline formatting workload,
#without retaining the final QS2 object's large coverage component in memory.
args <- commandArgs(trailingOnly = TRUE)
preflight <- length(args) == 3L && identical(args[[1]], "--preflight")
if(!preflight && (length(args) < 4L || length(args) > 5L)) {
  stop(paste(
    "Usage: benchmark_germline_formatting.R INPUT.qs2 BASELINE_outputResults.R CANDIDATE_outputResults.R NEW_OUTPUT_PREFIX [PAIRED_REPEATS=3]",
    "or: benchmark_germline_formatting.R --preflight BASELINE_outputResults.R CANDIDATE_outputResults.R", sep = "\n"))
}
options(warn = 2)
script_paths <- setNames(args[2:3], c("baseline", "candidate"))

#Read assignments from the actual old/new scripts, without running the pipeline.
#Recurse only into control-flow blocks; do not interpret unrelated expressions.
find_assignments <- function(node, name) {
  if(is.call(node) && identical(node[[1]], as.name("<-"))) {
    return(if(identical(node[[2]], as.name(name))) list(node) else list())
  }
  if(is.expression(node) || (is.call(node) && is.symbol(node[[1]]) && as.character(node[[1]]) %in% c("{", "for", "if"))) {
    return(unlist(lapply(as.list(node), find_assignments, name = name), recursive = FALSE))
  }
  list()
}
assignment <- function(expressions, name) {
  hits <- find_assignments(expressions, name)
  if(length(hits) != 1L) stop(paste("Expected one assignment for", name))
  hits[[1]]
}
#Resolve every required node before loading the large QS2 object. Preflight
#uses base R only, reads no scientific data and creates no output files.
parsed <- lapply(script_paths, parse)
initial_names <- c("chromgroups", "filtergroups", "strand_identical_cols_keep",
                   "strand_identical_cols_discard", "strand_redundant_cols_discard")
config_names <- c("region_read_filters_cols_keep", "region_genome_filters_cols_keep", "germline_filter_cols_keep")
required_names <- c(initial_names, config_names, "germlineVariantCalls.out")
resolved <- lapply(parsed, function(expressions) {
  setNames(lapply(required_names, function(name) assignment(expressions, name)), required_names)
})
cache_expressions <- parse(file.path(dirname(script_paths[["candidate"]]), "sharedFunctions.R"))
eval(assignment(cache_expressions, "cache_file"))
if(preflight) {
  cat("AST preflight passed:", length(required_names), "assignments in each script plus cache_file\n")
  print(script_paths)
  quit(save = "no", status = 0L)
}
suppressPackageStartupMessages({library(qs2); library(tidyverse)})
input <- args[[1]]
prefix <- args[[4]]
repeats <- if(length(args) == 5L) as.integer(args[[5]]) else 3L
stopifnot(!is.na(repeats), repeats > 0L)
metrics_path <- paste0(prefix, ".metrics.tsv")
comparison_path <- paste0(prefix, ".comparisons.tsv")
output_files <- unlist(lapply(seq_len(repeats), function(i) paste0(prefix, ".", i, ".", names(script_paths), ".qs2")))
if(any(file.exists(c(metrics_path, comparison_path, output_files)))) stop("Refusing to overwrite benchmark outputs")
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
cat("Loading final QS2; retaining only config and raw germline calls...\n")
object <- qs_read(input)
yaml.config <- object$yaml.config
germlineVariantCalls <- object$germlineVariantCalls
rm(object)
invisible(gc())
cat("Raw germline rows:", nrow(germlineVariantCalls), "\n")

format_all <- function(code_by_name) {
  env <- new.env(parent = globalenv())
  for(name in initial_names) eval(code_by_name[[name]], env)
  config_code <- code_by_name[config_names]
  format_code <- code_by_name[["germlineVariantCalls.out"]]
  formatted <- list()
  for(chromgroup in env$chromgroups) for(filtergroup in env$filtergroups) {
    env$i <- chromgroup
    env$j <- filtergroup
    for(code in config_code) eval(code, env)
    eval(format_code, env)
    formatted[[length(formatted) + 1L]] <- env$germlineVariantCalls.out
    rm("germlineVariantCalls.out", envir = env)
  }
  bind_rows(formatted)
}
rss <- function(field) {
  if(!file.exists("/proc/self/status")) return(NA_real_)
  line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
  if(!length(line)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line))
}
reset_hwm <- function() {
  tryCatch({
    con <- file("/proc/self/clear_refs", open = "w")
    on.exit(close(con))
    writeLines("5", con)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
}
metrics <- list()
comparisons <- list()
for(iteration in seq_len(repeats)) {
  #Alternate ordering to reduce systematic warm-up/order effects.
  modes <- if(iteration %% 2L) names(script_paths) else rev(names(script_paths))
  for(mode in modes) {
    invisible(gc())
    reset <- reset_hwm()
    before_rss <- rss("VmRSS")
    before <- proc.time()
    cat("Pair", iteration, mode, "formatting started\n")
    formatted <- format_all(resolved[[mode]])
    timing <- proc.time() - before
    metrics[[length(metrics) + 1L]] <- data.frame(
      iteration = iteration, mode = mode, rows = nrow(formatted), columns = ncol(formatted),
      user_seconds = unname(timing[[1]] + timing[[4]]),
      system_seconds = unname(timing[[2]] + timing[[5]]),
      actual_cpu_seconds = unname(sum(timing[c(1, 2, 4, 5)])), wall_seconds = unname(timing[[3]]),
      rss_before_kib = before_rss, rss_after_kib = rss("VmRSS"), peak_rss_kib = rss("VmHWM"),
      peak_is_phase_specific = reset)
    write.table(bind_rows(metrics), metrics_path, sep = "\t", row.names = FALSE, quote = FALSE)
    #Serialization and exact comparisons occur outside formatting measurements.
    qs_save(formatted, paste0(prefix, ".", iteration, ".", mode, ".qs2"))
    rm(formatted)
    invisible(gc())
  }
  baseline <- qs_read(paste0(prefix, ".", iteration, ".baseline.qs2"))
  candidate <- qs_read(paste0(prefix, ".", iteration, ".candidate.qs2"))
  same <- identical(baseline, candidate)
  comparisons[[length(comparisons) + 1L]] <- data.frame(iteration = iteration, identical = same,
                                                     rows = nrow(baseline), columns = ncol(baseline))
  write.table(bind_rows(comparisons), comparison_path, sep = "\t", row.names = FALSE, quote = FALSE)
  if(!same) stop(paste("Baseline/candidate differ in pair", iteration))
  rm(baseline, candidate)
  invisible(gc())
  cat("Pair", iteration, "exact QS2-loaded table comparison passed\n")
}
capture.output(sessionInfo(), file = paste0(prefix, ".session.txt"))
