#!/usr/bin/env Rscript
# Rscript tests/test_optional_germline_vcf_paths.R [repository root]
# Exercise the real filter-path assignments without loading any files.
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
env <- new.env(parent = globalenv())
for (expr in parse(file.path(repo, "bin", "sharedFunctions.R"))) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], as.name("cache_file"))) eval(expr, env)
}
env$cache_dir <- "/cache/"
env$individual_id_toanalyze <- "sample1"
env$qs_read <- function(path) {
  env$germline_variants_path <- path
  tibble()
}
base_individual <- list(individual_id = "sample1", germline_bam_file = "/data/person.bam", sex = "female")
vcf <- list(germline_vcf_file = "/data/person.vcf.gz", germline_vcf_type = "germline")
cases <- list(
  omitted = base_individual,
  null = c(base_individual, list(germline_vcf_files = NULL)),
  empty = c(base_individual, list(germline_vcf_files = list())),
  one = c(base_individual, list(germline_vcf_files = list(vcf))),
  multiple = c(base_individual, list(germline_vcf_files = list(vcf, vcf)))
)
expected <- c("/cache/sample1.person.bam.germline_vcf_variants.qs2",
  "/cache/person.bam.bw", "/cache/person.bam.vcf.gz")
for (case in names(cases)) {
  original_config <- list(cache_dir = "/cache", individuals = list(
    list(individual_id = "other", germline_bam_file = "/data/other.bam", germline_vcf_files = list()),
    cases[[case]]))
  for (prepared in c(FALSE, TRUE)) {
    env$yaml.config <- original_config
    expected_paths <- expected
    if (prepared) {
      expected_paths <- file.path("/prepared", basename(expected))
      env$yaml.config$cache_artifacts <- setNames(as.list(expected_paths), basename(expected))
    }
    before <- env$yaml.config
    eval(assignment("germline_vcf_variants"), env)
    eval(assignment("germline_bam_samtools_mpileup_file"), env)
    eval(assignment("germline_bam_bcftools_mpileup_file"), env)
    observed <- c(env$germline_variants_path, env$germline_bam_samtools_mpileup_file,
      env$germline_bam_bcftools_mpileup_file)
    stopifnot(identical(observed, expected_paths), identical(env$yaml.config, before))
  }
  cat("PASS:", case, "germline VCF list with standalone and prepared cache paths\n")
}
