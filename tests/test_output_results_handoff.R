#!/usr/bin/env Rscript
# Run in the pipeline R/Bioconductor environment on a compute node:
# Rscript --vanilla tests/test_output_results_handoff.R /absolute/fresh-output
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
stopifnot(length(script) == 1L)
root <- normalizePath(file.path(dirname(normalizePath(script)), ".."), mustWork = TRUE)
suppressPackageStartupMessages({library(GenomicRanges); library(tidyverse); library(qs2)})
options(warn = 2)
out <- args[[1L]]
stopifnot(!file.exists(out), dir.create(out))
out <- normalizePath(out)
checks <- list()
check <- function(name, condition) {
  stopifnot(!name %in% names(checks), isTRUE(condition))
  checks[[name]] <<- TRUE
}
rscript <- file.path(R.home("bin"), "Rscript")
worker <- file.path(root, "tests", "fixtures", "output_results_handoff_worker.R")
quote_arg <- function(x) shQuote(as.character(x))
launch <- function(argv, log, env = character()) {
  system2(rscript, c("--vanilla", quote_arg(worker), quote_arg(argv)),
          stdout = log, stderr = log, env = env, wait = TRUE)
}
run_case <- function(name, inputs, rows, fail = "", existing_final = FALSE) {
  directory <- file.path(out, name)
  stopifnot(dir.create(directory))
  basename <- file.path(directory, "result with spaces ' $x; literal")
  final <- paste0(basename, ".outputResults.qs2")
  config <- file.path(directory, "config with spaces.rds")
  spec <- list(rows = rows, fail = fail, expected = file.path(directory, "expected.rds"),
               probe = file.path(directory, "worker-probe.rds"))
  saveRDS(spec, config)
  if (existing_final) qs_save(list(old = TRUE), final)
  status <- launch(c("-c", config, "-s", "sample ' $x; literal", "-f", paste(inputs, collapse = ","),
                     "-o", basename), file.path(directory, "worker.log"))
  private <- list.files(directory, pattern = "^\\.hidef-output-parts-", all.files = TRUE, full.names = TRUE)
  if (!nzchar(fail)) {
    check(paste0(name, "_exact"), status == 0L && identical(qs_read(final), readRDS(spec$expected)))
    check(paste0(name, "_clean_private"), length(private) == 0L)
    public_qs <- list.files(directory, pattern = "\\.outputResults\\.qs2$", full.names = TRUE)
    check(paste0(name, "_one_public_QS"), identical(public_qs, final))
    probe <- readRDS(spec$probe)
    log <- readLines(file.path(directory, "worker.log"))
    roles <- c("launcher", "worker", "assembly")
    role_lines <- lapply(roles, function(role) grep(paste0("^## outputResults ", role, " PID: "), log))
    stopifnot(all(lengths(role_lines) == 1L))
    pids <- as.integer(vapply(role_lines, function(i) sub(".* PID: ", "", log[[i]]), character(1)))
    exited <- grep("^## outputResults calculation/export child exited successfully$", log)
    check(paste0(name, "_sequential_distinct_R_processes"), identical(probe$stage, "worker") &&
            !anyDuplicated(pids) && !Sys.getpid() %in% pids && probe$worker_pid == pids[[2L]] &&
            length(exited) == 1L && role_lines[[1L]] < role_lines[[2L]] &&
            role_lines[[2L]] < exited && exited < role_lines[[3L]] &&
            identical(probe$launch$args, c("-c", config, "-s", "sample ' $x; literal", "-f",
                                         paste(inputs, collapse = ","), "-o", basename)))
  } else {
    check(paste0(name, "_nonzero_no_final"), status != 0L && !file.exists(final))
    check(paste0(name, "_evidence_retained"), length(private) == 1L && dir.exists(private))
  }
}
run_case("one", "a", 2L)
run_case("two_with_empty_table", c("input one ' $x; literal", "b"), c(0L, 2L))
run_case("four", letters[1:4], c(2L, 1L, 0L, 2L))
run_case("seven", letters[1:7], c(0L, 1L, 2L, 0L, 1L, 0L, 1L))
run_case("duplicate_inputs", c("a", "a", "b", "a"), c(2L, 2L, 1L, 2L))
run_case("all_empty_tables", c("a", "b"), c(0L, 0L))
run_case("replace_existing_final", "a", 1L, existing_final = TRUE)
run_case("worker_failure", "a", 1L, fail = "worker")
run_case("export_failure", "a", 1L, fail = "exports")
run_case("assembly_failure", "a", 1L, fail = "assembly")
empty <- file.path(out, "empty-input")
stopifnot(dir.create(empty))
status <- launch(c("-c", "unused", "-s", "sample", "-f", "", "-o", file.path(empty, "result")),
                 file.path(empty, "worker.log"))
check("empty_input_fails_before_parts", status != 0L &&
        !length(list.files(empty, pattern = "^\\.hidef-output-parts-", all.files = TRUE)))
missing <- file.path(out, "missing-argument")
stopifnot(dir.create(missing))
status <- launch(c("-c", "unused", "-s", "sample", "-o", file.path(missing, "result")),
                 file.path(missing, "worker.log"))
check("missing_argument_fails_before_parts", status != 0L &&
        !length(list.files(missing, pattern = "^\\.hidef-output-parts-", all.files = TRUE)))
check("unknown_stage_rejected", launch(c("-h"), file.path(out, "stage.log"),
      "HIDEF_OUTPUT_RESULTS_STAGE=unexpected") != 0L)
check("stage_without_handoff_rejected", launch(c("-h"), file.path(out, "handoff.log"),
      "HIDEF_OUTPUT_RESULTS_STAGE=worker") != 0L)
check("inherited_handoff_rejected", launch(c("-h"), file.path(out, "inherited.log"),
      "HIDEF_OUTPUT_RESULTS_PARTS=/nonexistent") != 0L)
check("public_help_preserved", launch("--help", file.path(out, "help.log")) == 0L &&
      any(grepl("--output_basename", readLines(file.path(out, "help.log")), fixed = TRUE)))
check("parent_environment_unchanged", identical(Sys.getenv("HIDEF_OUTPUT_RESULTS_STAGE"), "") &&
      identical(Sys.getenv("HIDEF_OUTPUT_RESULTS_PARTS"), ""))
# Private artifacts use only rds/qs2, so none match the existing scientific export glob.
failure_parts <- unlist(lapply(c("worker_failure", "export_failure", "assembly_failure"), function(n) {
  list.files(file.path(out, n), pattern = "^\\.hidef-output-parts-", all.files = TRUE, full.names = TRUE)
}), use.names = FALSE)
check("private_files_outside_public_export_glob", !any(grepl("\\.(tsv|vcf\\.bgz|vcf\\.bgz\\.tbi|pdf)$",
      unlist(lapply(failure_parts, list.files, recursive = TRUE, all.files = TRUE), use.names = FALSE))))
stopifnot(length(checks) == 42L)
jsonlite::write_json(list(exact = TRUE, count = length(checks), checks = checks),
                    file.path(out, "completion.json"), pretty = TRUE, auto_unbox = TRUE)
cat("COMPLETE42 CLI/lifecycle/ordered-attribute fixtures\n")
