#!/usr/bin/env Rscript
# Experimental shared-session runner; never invoked by the production workflow.
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 8L) stop(paste(
  "Usage: benchmark_filter_chain.R SCRIPT_BIN CONFIG EXTRACT.qs2 SAMPLE CHROMGROUP OUTPUT_DIR FILTERGROUP1 FILTERGROUP2 [FILTERGROUP...]"))
script_bin <- normalizePath(args[[1]], mustWork = TRUE)
config <- normalizePath(args[[2]], mustWork = TRUE)
extract <- normalizePath(args[[3]], mustWork = TRUE)
sample_id <- args[[4]]
chromgroup <- args[[5]]
output_dir <- normalizePath(args[[6]], mustWork = TRUE)
filtergroups <- args[-seq_len(6L)]
if(anyDuplicated(filtergroups)) stop("Filtergroups must be distinct")
script <- parse(file.path(script_bin, "filterCalls.R"), keep.source = FALSE)
helper_path <- normalizePath(file.path(script_bin, "sharedFunctions.R"), mustWork = TRUE)
Sys.setenv(PATH = paste(script_bin, Sys.getenv("PATH"), sep = .Platform$path.sep))
cache <- new.env(parent = emptyenv())
cache$loaded <- FALSE
cache$reads <- 0L
cache$hits <- 0L
group_stats <- list()

run_group <- function(group, index) {
  directory <- file.path(output_dir, sprintf("group%03d", index))
  if(dir.exists(directory)) stop("Group output directory already exists: ", directory)
  dir.create(directory)
  prior_directory <- getwd()
  prior_options <- options()
  on.exit({setwd(prior_directory); options(prior_options)}, add = TRUE)
  setwd(directory)
  log <- file("filter.log", open = "wt")
  sink(log)
  sink(log, type = "message")
  on.exit({sink(type = "message"); sink(); close(log)}, add = TRUE)
  group_env <- new.env(parent = globalenv())
  filter_args <- c("-c", config, "-f", extract, "-s", sample_id,
                   "-g", chromgroup, "-v", group, "-o", "output.qs2")
  group_env$parse_args <- function(object, ...) {
    optparse::parse_args(object, args = filter_args, ...)
  }
  group_env$source <- function(file, local = FALSE, ...) {
    # In a normal Rscript run these definitions and pipeline variables share
    # .GlobalEnv. Bind them to this group's environment to preserve that scope
    # while preventing state leaking into the next group.
    if(!identical(normalizePath(file, mustWork = TRUE), helper_path) || !identical(local, FALSE)) {
      stop("Experimental runner supports only the filter script's sharedFunctions source call")
    }
    base::source(file, local = group_env, ...)
  }
  group_env$qs_read <- function(file, ...) {
    if(length(file) == 1L && identical(normalizePath(file, mustWork = TRUE), extract)) {
      if(!cache$loaded) {
        cache$object <- qs2::qs_read(file, ...)
        cache$loaded <- TRUE
        cache$reads <- cache$reads + 1L
      } else {
        cache$hits <- cache$hits + 1L
      }
      return(cache$object)
    }
    qs2::qs_read(file, ...)
  }
  before <- proc.time()
  eval(script, envir = group_env)
  timing <- proc.time() - before
  if(!file.exists("output.qs2")) stop("Filter script did not produce its complete QS2 output")
  result <- data.frame(index = index, filtergroup = group,
                       user_seconds = unname(timing[[1]] + timing[[4]]),
                       system_seconds = unname(timing[[2]] + timing[[5]]),
                       wall_seconds = unname(timing[[3]]))
  # Explicitly release all per-group objects, preserving only the extraction
  # cache. R copy-on-modify semantics apply; final QS comparisons remain the
  # authority for detecting any accidental shared-object mutation.
  rm(list = ls(group_env, all.names = TRUE), envir = group_env)
  invisible(gc())
  result
}

for(index in seq_along(filtergroups)) {
  group_stats[[index]] <- run_group(filtergroups[[index]], index)
}
if(cache$reads != 1L || cache$hits != length(filtergroups) - 1L) {
  stop("Unexpected extraction-read pattern; cache interception needs review")
}
write.table(do.call(rbind, group_stats), file.path(output_dir, "group-timings.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(extraction_reads = cache$reads, extraction_cache_hits = cache$hits),
            file.path(output_dir, "cache-usage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
