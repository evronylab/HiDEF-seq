#!/usr/bin/env Rscript
# Tiny test of isolation/cached-read interception; no production data needed.
suppressPackageStartupMessages({library(qs2); library(optparse)})
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if(length(args)) args[[1]] else ".", mustWork = TRUE)
temporary <- tempfile(pattern = "filter-chain-test-")
dir.create(temporary)
on.exit_cleanup <- function() unlink(temporary, recursive = TRUE)
script_bin <- file.path(temporary, "bin")
dir.create(script_bin)
helper <- file.path(script_bin, "sharedFunctions.R")
writeLines("helper <- function() list(group = opt$filtergroup, marker = group_marker)", helper)
Sys.chmod(helper, "0755")
script <- file.path(script_bin, "filterCalls.R")
writeLines(c(
  "suppressPackageStartupMessages({library(qs2); library(optparse)})",
  "source(Sys.which('sharedFunctions.R'))",
  "option_list <- list(make_option(c('-c','--config')), make_option(c('-f','--file')),",
  "  make_option(c('-s','--sample')), make_option(c('-g','--chromgroup')),",
  "  make_option(c('-v','--filtergroup')), make_option(c('-o','--output')))",
  "opt <- parse_args(OptionParser(option_list = option_list))",
  "options(warn = 2)",
  "stopifnot(!file.exists('group-marker.txt'))",
  "writeLines(opt$filtergroup, 'group-marker.txt')",
  "group_marker <- paste0('private-', opt$filtergroup)",
  "data <- qs_read(opt$file)",
  "stopifnot(identical(data$values, 1:3))",
  "data$values[1] <- if(opt$filtergroup == 'first') 100L else 200L",
  "qs_save(list(data = data, context = helper()), opt$output)"
), script)
config <- file.path(temporary, "config.yaml")
writeLines("unused: true", config)
input <- file.path(temporary, "input.qs2")
qs_save(list(values = 1:3), input)
output <- file.path(temporary, "shared")
dir.create(output)
rscript <- file.path(R.home("bin"), "Rscript")
status <- system2(rscript, shQuote(c("--vanilla", file.path(repo, "scripts", "benchmark", "benchmark_filter_chain.R"),
                                   script_bin, config, input, "sample", "chr1", output, "first", "second")))
stopifnot(status == 0L)
first <- qs_read(file.path(output, "group001", "output.qs2"))
second <- qs_read(file.path(output, "group002", "output.qs2"))
stopifnot(identical(first, list(data = list(values = c(100L, 2L, 3L)),
                              context = list(group = "first", marker = "private-first"))))
stopifnot(identical(second, list(data = list(values = c(200L, 2L, 3L)),
                               context = list(group = "second", marker = "private-second"))))
usage <- read.delim(file.path(output, "cache-usage.tsv"))
stopifnot(usage$extraction_reads == 1L, usage$extraction_cache_hits == 1L)
on.exit_cleanup()
cat("PASS: shared filter runner isolates helper globals/options/work directories and caches one immutable extraction read\n")
