#!/usr/bin/env Rscript
# Rscript tests/test_optional_region_filters.R [repository root]
suppressPackageStartupMessages(library(tidyverse))
options(warn = 2)
args <- commandArgs(TRUE)
repo <- if (length(args)) args[[1L]] else "."
expressions <- parse(file.path(repo, "bin", "filterCalls.R"))
assignment <- function(name) {
  found <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
    identical(x[[2L]], as.name(name)), expressions)
  stopifnot(length(found) >= 1L)
  found[[1L]]
}
individual_stats <- Filter(function(x)
  grepl('individual_filter_pattern = "region_read_filter_.*passfilter$"',
        paste(deparse(x, width.cutoff = 500L), collapse = " "), fixed = TRUE), expressions)
stopifnot(length(individual_stats) == 1L)
individual_stats <- individual_stats[[1L]]
env <- new.env(parent = globalenv())
eval(assignment("calculate_molecule_stats_frombam"), env)
for (expr in parse(file.path(repo, "bin", "sharedFunctions.R"))) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% c("region_filter_config", "cache_file")) {
    eval(expr, env)
  }
}
env$chromgroup_toanalyze <- "nuclear"
env$filtergroup_toanalyze <- "strict"
env$extractCallsFile <- "fixture.extractCalls.qs2"
env$bam <- tibble(run_id = factor(rep("run1", 4L)), zm = c(1L, 1L, 2L, 2L),
  start = c(1L, 1L, 11L, 11L), end = c(10L, 10L, 20L, 20L))
initial <- tibble(run_id = factor("run1"), chromgroup = "nuclear", filtergroup = "strict",
  extractCallsFile = env$extractCallsFile, stat = "previous", value = 4L)

# No applicable filters: the existing stats must survive exactly, and no
# individual-filter statistic is invented when there is no filter to report.
env$region_read_filters_config <- tibble(region_filter_file = character())
env$molecule_stats <- initial
eval(individual_stats, env)
stopifnot(identical(env$molecule_stats, initial))

# Nonempty path: preserve the exact original helper output and row order.
env$region_read_filters_config <- tibble(region_filter_file = "fixture.bw")
env$bam$region_read_filter_fixture.passfilter <- c(TRUE, TRUE, FALSE, FALSE)
expected <- bind_rows(initial, env$calculate_molecule_stats_frombam(env$bam,
  individual_filter_pattern = "region_read_filter_.*passfilter$", all_filter_suffix = NULL,
  chromgroup_toanalyze.input = env$chromgroup_toanalyze,
  filtergroup_toanalyze.input = env$filtergroup_toanalyze,
  extractCallsFile.input = env$extractCallsFile))
env$molecule_stats <- initial
eval(individual_stats, env)
stopifnot(identical(env$molecule_stats, expected),
  identical(tail(env$molecule_stats$value, 2L), c(1L, 20L)))
cat("PASS: no applicable read-filter stats and unchanged nonempty stats, values, schema and order\n")

# Preserve populated configuration inference exactly; optional empty lists get
# a usable zero-row schema for both filtering and germline export selection.
read_one <- list(region_filter_file = "/input/read1.bw", binsize = 1L, threshold = "gte0.1",
  padding = 2L, read_threshold = "gte0.5", applyto_chromgroups = "nuclear",
  applyto_filtergroups = "strict", is_germline_filter = TRUE)
read_two <- modifyList(read_one, list(region_filter_file = "/input/read2.bw", binsize = 5,
  applyto_chromgroups = "mitochondrial", is_germline_filter = FALSE))
genome_one <- read_one[setdiff(names(read_one), "read_threshold")]
genome_one$region_filter_file <- "/input/genome1.bw"
configured <- list(cache_dir = "/cache", region_filters = list(
  list(read_filters = list(read_one, read_two)), list(genome_filters = list(genome_one))))
legacy <- function(config, type) config$region_filters %>% map(type) %>% flatten %>%
  enframe(name = NULL) %>% unnest_wider(value)
for (type in c("read_filters", "genome_filters")) {
  stopifnot(identical(env$region_filter_config(configured, type), legacy(configured, type)))
}
env$cache_dir <- "/cache/"
env$i <- "nuclear"
env$j <- "strict"
output <- parse(file.path(repo, "bin", "outputResults.R"))
find_assignment <- function(nodes, name) {
  found <- list()
  visit <- function(node) {
    if (!is.call(node)) return(invisible(NULL))
    if (identical(node[[1L]], as.name("<-")) && identical(node[[2L]], as.name(name))) {
      found[[length(found) + 1L]] <<- node
    } else if (is.symbol(node[[1L]]) && as.character(node[[1L]]) %in% c("for", "if", "{")) {
      for (part in as.list(node)[-1L]) visit(part)
    }
    invisible(NULL)
  }
  for (node in nodes) visit(node)
  stopifnot(length(found) == 1L)
  found[[1L]]
}
configurations <- list(
  configured = configured,
  zero_applicable = modifyList(configured, list()),
  omitted = list(cache_dir = "/cache"),
  empty_parent = list(cache_dir = "/cache", region_filters = list()),
  empty_children = list(cache_dir = "/cache", region_filters = list(
    list(read_filters = list(), genome_filters = list()))),
  null_children = list(cache_dir = "/cache", region_filters = list(
    list(read_filters = NULL, genome_filters = NULL)))
)
for (case in names(configurations)) {
  env$yaml.config <- configurations[[case]]
  env$chromgroup_toanalyze <- env$i <- if (case == "zero_applicable") "unselected" else "nuclear"
  for (type in c("read", "genome")) {
    target <- paste0("region_", type, "_filters_config")
    eval(assignment(target), env)
    expected_rows <- if (case == "configured") 1L else 0L
    stopifnot(nrow(env[[target]]) == expected_rows)
    export_target <- paste0("region_", type, "_filters_cols_keep")
    eval(find_assignment(output, export_target), env)
    expected_export <- if (case == "configured") paste0("region_", type, "_filter_",
      if (type == "read") "read1" else "genome1", ".bw.bin1.gte0.1.bw.passfilter") else character()
    stopifnot(identical(env[[export_target]], expected_export))
  }
}
cat("PASS: populated/empty/omitted/zero-applicable region configuration and output column selection\n")
