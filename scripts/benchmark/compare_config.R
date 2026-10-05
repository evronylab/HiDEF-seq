#!/usr/bin/env Rscript
# Prove scientific configuration equality independently of run metadata.
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3L) stop(paste(
  "Usage: compare_config.R ORIGINAL.yaml EFFECTIVE.yaml REPORT.tsv",
  "[--metadata-key=EXACT_KEY] [--allow-extra=EXACT_KEY]"))
flags <- if(length(args) > 3L) args[4:length(args)] else character()
if(any(!startsWith(flags, "--metadata-key=") & !startsWith(flags, "--allow-extra="))) stop("Unknown option")
if(file.exists(args[[3]])) stop("Report already exists")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(script_file)), "compare_qs2.R"))
suppressPackageStartupMessages(library(configr))
original <- suppressWarnings(read.config(args[[1]]))
effective <- suppressWarnings(read.config(args[[2]]))
if(!is.list(original) || !is.list(effective) || is.null(names(original)) || is.null(names(effective))) {
  stop("Both configurations must parse as named mappings through configr::read.config")
}
metadata <- sub("^--metadata-key=", "", flags[startsWith(flags, "--metadata-key=")])
allowed_extra <- c("germline_coverage_filters", "reference_cache_dir", "reference_summary_file", "cache_artifacts",
                   sub("^--allow-extra=", "", flags[startsWith(flags, "--allow-extra=")]))
if(any(!metadata %in% names(original))) stop("Each metadata-key must name an actual original configuration key")
missing <- setdiff(names(original), names(effective))
extras <- setdiff(names(effective), names(original))
unknown <- setdiff(extras, allowed_extra)
reports <- list()
add <- function(key, status, reason) {
  reports[[length(reports) + 1L]] <<- data.frame(path = paste0("/yaml.config/", key), status = status,
                                               reason = reason, max_relative = NA_real_)
}
for(key in missing) add(key, "failure", "original configuration key missing")
for(key in unknown) add(key, "failure", "extra key outside explicit metadata whitelist")
for(key in intersect(extras, allowed_extra)) add(key, "ignored", "explicit added-metadata whitelist")
for(key in metadata) add(key, "ignored", "explicit original run-metadata whitelist")
# Top-level YAML mapping order is not scientific data. Preserve the original
# order when looking up the same keys, then compare all their nested values,
# classes and attributes exactly; nested scientific list order remains strict.
keys <- setdiff(intersect(names(original), names(effective)), metadata)
result <- scientific_compare(original[keys], effective[keys], path = "/yaml.config", ignore = character())
reports[[length(reports) + 1L]] <- result$reports
report <- do.call(rbind, reports)
write.table(report, args[[3]], sep = "\t", quote = TRUE, row.names = FALSE)
failure <- length(missing) + length(unknown) + result$counts[["failure"]]
review <- result$counts[["review"]]
cat("Original keys compared:", length(keys), "failures:", failure, "numeric reviews:", review, "\n")
quit(status = if(failure) 1L else if(review) 2L else 0L)
