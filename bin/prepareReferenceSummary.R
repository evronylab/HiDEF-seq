#!/usr/bin/env -S Rscript --vanilla

# Reference-wide data are deliberately retained here: filter statistics count
# all N bases, and burden normalization also needs the complete reference.
suppressPackageStartupMessages({
  library(optparse)
  library(BSgenome)
  library(GenomicRanges)
  library(Biostrings)
  library(configr)
  library(qs2)
  library(tidyverse)
})
source(Sys.which("sharedFunctions.R"))
source(Sys.which("referenceSummaryFunctions.R"))
options(warn = 2)
option_list <- list(
  make_option(c("-c", "--config"), type = "character"),
  make_option(c("-o", "--output"), type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$config) || is.null(opt$output)) stop("--config and --output are required")
yaml.config <- suppressWarnings(read.config(opt$config))
BSgenome_name <- get_bsgenome_name(yaml.config)
suppressPackageStartupMessages(library(BSgenome_name, character.only = TRUE,
                                       lib.loc = reference_cache_dir(yaml.config)))
genome <- get(BSgenome_name)
summary <- list(n_ranges = reference_n_ranges(genome),
                trinucleotide_counts = reference_trinucleotide_counts(genome))
qs2::qs_save(summary, opt$output)
