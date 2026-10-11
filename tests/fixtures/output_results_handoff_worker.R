# Tiny data only; exercise the real launcher, save, assembly and cleanup helpers.
# Test-only loader: evaluate the installed script's private function definitions.
# The scientific body is replaced by tiny explicit data, never a legacy runtime.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
stopifnot(length(script) == 1L)
root <- normalizePath(file.path(dirname(normalizePath(script)), "../.."), mustWork = TRUE)
expressions <- parse(file.path(root, "bin", "outputResults.R"), keep.source = FALSE)
entry <- which(vapply(expressions, function(x) identical(x,
  quote(output_results_entry(commandArgs(trailingOnly = TRUE)))), logical(1)))
stopifnot(length(entry) == 1L, entry > 1L)
definitions <- expressions[seq_len(entry - 1L)]
stopifnot(all(vapply(definitions, function(x) is.call(x) &&
  identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) &&
  startsWith(as.character(x[[2L]]), "output_results_") &&
  is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("function")), logical(1))))
for (x in definitions) eval(x, .GlobalEnv)
output_results_entry(commandArgs(trailingOnly = TRUE))
suppressPackageStartupMessages({library(GenomicRanges); library(tidyverse); library(qs2)})
opt <- output_results_parse_cli(commandArgs(trailingOnly = TRUE))
spec <- readRDS(opt$config)
inputs <- stringr::str_split_1(opt$files, ",") |> trimws()
stopifnot(length(inputs) == length(spec$rows))
launch <- output_results_launch()
saveRDS(list(worker_pid = Sys.getpid(), launch = launch,
             stage = Sys.getenv("HIDEF_OUTPUT_RESULTS_STAGE")), spec$probe)
if (identical(spec$fail, "worker")) stop("INTENTIONAL_TINY_WORKER_FAILURE")
parts <- lapply(seq_along(inputs), function(i) {
  n <- spec$rows[[i]]
  tibble(scope = factor(rep(inputs[[i]], n), levels = unique(c("unused", inputs))),
         cell = rep(list(structure(1:2, names = c("", "a"))), n),
         bam.gr.filtertrack.coverage = rep(list(IRanges::RleList(
           chr1 = Rle(c(0L, 1L, 1L, 0L)), chrM = Rle(integer()), compress = FALSE)), n))
})
thin <- lapply(seq_along(inputs), function(i) output_results_save_coverage(parts[[i]], i, inputs))
component <- "bam.gr.filtertrack.bytype.coverage_tnc"
# Public component names/order are an explicit oracle, independent of the helper.
component_names <- c(
  "yaml.config", "run_metadata", "individual_id", "sample_id",
  "molecule_stats.by_run_id", "molecule_stats.by_analysis_id",
  "region_genome_filter_stats", "finalCalls", "finalCalls.bytype",
  "germlineVariantCalls", "germlineVariantCalls_for_tsv",
  "finalCalls.reftnc_spectra", "bam.gr.filtertrack.bytype.coverage_tnc",
  "genome.reftnc", "genome_chromgroup.reftnc", "sensitivity",
  "finalCalls.burdens", "estimatedSBSMutationErrorProbability"
)
expected <- setNames(rep(list(structure(1:2, names = c(NA_character_, ""))), 18L),
                     component_names)
expected[[component]] <- bind_rows(parts)
saveRDS(expected, spec$expected)
small <- expected
small[[component]] <- bind_rows(thin)
output_results_save_remainder(small, paste0(opt$output_basename, ".outputResults.qs2"))
if (identical(spec$fail, "exports")) stop("INTENTIONAL_TINY_EXPORT_FAILURE")
writeLines("unchanged tiny export", paste0(opt$output_basename, ".yaml_config.tsv"))
output_results_worker_done()
if (identical(spec$fail, "assembly")) {
  writeLines("INTENTIONAL_INVALID_QS", file.path(launch$directory, "remainder.qs2"))
}
