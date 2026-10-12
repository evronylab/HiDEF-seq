#!/usr/bin/env Rscript
# Run in the pipeline R image: Rscript tests/test_burden_optional_sensitivity.R [repo or calculateBurdens.R]
# Evaluate the production configuration and dispatch expressions without loading
# references, calls or coverage. Blank/omitted use_chromgroup is a valid setting.
suppressPackageStartupMessages(library(tidyverse))
options(warn = 2)
args <- commandArgs(TRUE)
source_path <- if (length(args)) args[[1L]] else "."
if (dir.exists(source_path)) source_path <- file.path(source_path, "bin", "calculateBurdens.R")
expressions <- parse(source_path)
assignment <- function(name) {
  found <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
    identical(x[[2L]], as.name(name)), expressions)
  stopifnot(length(found) >= 1L)
  found[[1L]]
}
dispatch <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("if")) &&
  grepl("All non-MDB sensitivities set to 1", paste(deparse(x), collapse = " "), fixed = TRUE), expressions)
stopifnot(length(dispatch) == 1L)
dispatch <- dispatch[[1L]]
# Enabled sensitivity science is covered by the coverage/sensitivity suites.
# Keep its production condition and both real disabled/other-group branches.
dispatch[[3L]] <- quote(selected_branch <- "enabled")
call_types <- tibble(call_type = c("SBS", "SBS", "insertion", "deletion"),
  call_class = c("SBS", "SBS", "indel", "indel"),
  SBSindel_call_type = c("mutation", "mismatch-ss", "mutation", "mutation"),
  filtergroup = factor("strict"), analyzein_chromgroups = "all")
settings <- list(
  blank_yaml = list(use_chromgroup = NULL, genotype = "heterozygous"),
  omitted_key = list(genotype = "heterozygous"),
  current_group = list(use_chromgroup = "nuclear", genotype = "heterozygous"),
  other_group = list(use_chromgroup = "nuclear", genotype = "heterozygous")
)
for (case in names(settings)) {
  env <- new.env(parent = globalenv())
  env$yaml.config <- list(sensitivity_parameters = settings[[case]])
  env$chromgroup_toanalyze <- if (case == "other_group") "mitochondrial" else "nuclear"
  env$call_types_toanalyze <- call_types
  env$selected_branch <- "disabled_or_other"
  eval(assignment("sensitivity_parameters"), env)
  eval(assignment("sensitivity_enabled_for_chromgroup"), env)
  eval(assignment("sensitivity"), env)
  eval(dispatch, env)
  stopifnot(identical(env$sensitivity_enabled_for_chromgroup, case == "current_group"),
    identical(env$selected_branch, if (case == "current_group") "enabled" else "disabled_or_other"))
  expected <- call_types[, setdiff(names(call_types), "analyzein_chromgroups")]
  expected$sensitivity <- rep(if (case == "other_group") NA_real_ else 1, nrow(expected))
  expected$sensitivity_source <- rep(if (case == "other_group") "other_chromgroup" else "default", nrow(expected))
  stopifnot(identical(env$sensitivity, expected))
  cat("PASS:", case, "sensitivity dispatch, values, source, types and schema\n")
}
