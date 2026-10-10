#!/usr/bin/env Rscript
# Regression of actual installed helpers against the preserved 7e function.
# Run in a compute allocation with the pipeline R/Bioconductor environment.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
stopifnot(!file.exists(output), dir.create(output))
suppressPackageStartupMessages({library(GenomicRanges); library(plyranges); library(tidyverse)})
options(warn = 2)
read_definition <- function(path, name) {
  nodes <- as.list(parse(path, keep.source = FALSE))
  found <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
                    identical(x[[2L]], as.name(name)), nodes)
  stopifnot(length(found) == 1L)
  found[[1L]]
}
read_helper <- function(path, name) {
  env <- new.env(parent = globalenv())
  eval(read_definition(path, name), env)
  env[[name]]
}
baseline_path <- file.path(root, "tests", "fixtures", "overlapsAny_bymcols_7e37ada.R")
filter_path <- file.path(root, "bin", "filterCalls.R")
shared_path <- file.path(root, "bin", "sharedFunctions.R")
allele_key_names <- c("call_class", "call_type", "ref_plus_strand", "alt_plus_strand")
for (name in "observed_allele_interaction") {
  eval(read_definition(shared_path, name), globalenv())
}
# The overlap function may differ only by its two calls to the tested encoder.
replacements <- 0L
rewrite <- function(x) {
  if (identical(x, quote(interaction(drop = TRUE)))) {
    replacements <<- replacements + 1L
    return(quote(observed_allele_interaction()))
  }
  if (is.call(x)) return(as.call(lapply(as.list(x), rewrite)))
  x
}
expected <- rewrite(read_definition(baseline_path, "overlapsAny_bymcols"))
stopifnot(replacements == 2L,
          identical(expected, read_definition(filter_path, "overlapsAny_bymcols")))
oracle <- read_helper(baseline_path, "overlapsAny_bymcols")
candidate <- read_helper(filter_path, "overlapsAny_bymcols")
capture <- function(expr) {
  warnings <- list()
  value <- tryCatch(withCallingHandlers(list(value = force(expr)),
    warning = function(w) warnings[[length(warnings) + 1L]] <<-
      list(class = class(w), message = conditionMessage(w))),
    error = function(e) list(error = list(class = class(e), message = conditionMessage(e))))
  list(value = value, warnings = warnings)
}
cases <- list()
add <- function(name, columns, fast) cases[[name]] <<- list(columns = columns, fast = fast)
a <- list(call_class = factor(c("SBS", "indel", "SBS", "SBS", "indel"), levels = c("unused", "indel", "SBS")),
          call_type = factor(c("SBS", "insertion", "SBS", "SBS", "deletion"), levels = c("deletion", "unused", "insertion", "SBS")),
          ref_plus_strand = c("A", "", "T", "A", "GG"),
          alt_plus_strand = c("C", "AA", "G", "C", ""))
add("ordinary", a, TRUE)
b <- lapply(a, rev); add("reversed_rows", b, TRUE)
b <- a; b[1:2] <- lapply(b[1:2], function(x) factor(x, levels = rev(levels(x))))
add("reversed_factor_levels", b, TRUE)
b <- lapply(a, as.character); add("all_character", b, TRUE)
b <- lapply(a, function(x) x[c(1L, 1L, 5L, 2L, 1L)]); add("duplicate_rows", b, TRUE)
b <- a; b$ref_plus_strand[1] <- NA_character_; add("missing_value", b, FALSE)
b <- a; b$call_class <- addNA(b$call_class); add("explicit_unused_NA_level", b, FALSE)
b$call_class[1] <- NA; b$call_class <- addNA(b$call_class); add("explicit_used_NA_level", b, FALSE)
b <- a; b$ref_plus_strand[1] <- "NA"; add("literal_NA_string", b, TRUE)
b <- a; b$ref_plus_strand[1] <- "A.C"; add("separator_value", b, FALSE)
b <- a; b$ref_plus_strand <- c("A.B", "A", "T", "A", "GG")
b$alt_plus_strand <- c("C", "B.C", "G", "C", ""); add("separator_collision", b, FALSE)
b <- a; b$call_class <- ordered(b$call_class); add("ordered_factor", b, FALSE)
b <- a; attr(b$ref_plus_strand, "note") <- "preserve"; add("custom_attribute", b, FALSE)
b <- a; names(b$ref_plus_strand) <- letters[1:5]; add("named_vector", b, FALSE)
b <- a; b$ref_plus_strand <- Rle(b$ref_plus_strand); add("Rle_column", b, FALSE)
b <- a; b$ref_plus_strand <- as.list(b$ref_plus_strand); add("list_column", b, FALSE)
b <- lapply(a, function(x) x[0L]); add("empty_columns", b, FALSE)
b <- a[4:1]; add("reversed_key_order", b, FALSE)
b <- a; names(b)[1] <- "other"; add("other_key_names", b, FALSE)
b <- a; b$ref_plus_strand <- 1:5; add("numeric_column", b, FALSE)

rows <- list()
check <- function(label, before, after, require_success = TRUE) {
  expected <- capture(before); actual <- capture(after)
  same <- identical(expected, actual) && (!require_success || is.null(expected$value$error))
  rows[[length(rows) + 1L]] <<- data.frame(label = label, exact = same,
    baseline_success = is.null(expected$value$error))
  if (!same) saveRDS(list(expected = expected, actual = actual),
    file.path(output, paste0("mismatch-", length(rows), ".rds")))
}
for (name in names(cases)) {
  f <- cases[[name]]
  before <- serialize(f$columns, NULL)
  check(paste0("factor:", name), base::interaction(f$columns, drop = TRUE),
        observed_allele_interaction(f$columns), require_success = f$fast)
  stopifnot(identical(before, serialize(f$columns, NULL)))
}

make_gr <- function(columns, query = TRUE) {
  x <- GRanges("chr1", IRanges(c(10L, 20L, 30L, 10L, 50L),
    width = if (query) c(1L, 0L, 1L, 1L, 2L) else c(1L, 1L, 1L, 1L, 2L)),
    strand = if (query) c("+", "-", "*", "+", "-") else "*")
  mcols(x) <- DataFrame(columns)
  x$run_id <- factor(rep("runA", length(x)), levels = c("unused", "runA"))
  x$zm <- seq_along(x)
  names(x) <- paste0("row", seq_along(x)); x
}
# Full helper comparisons cover independently ordered factor levels, reversed
# query/subject rows, missing values and separator labels.
for (name in c("ordinary", "reversed_rows", "reversed_factor_levels", "all_character",
               "separator_value", "missing_value", "explicit_unused_NA_level", "explicit_used_NA_level")) {
  query <- make_gr(cases[[name]]$columns)
  subject <- make_gr(a, query = FALSE)
  subject$call_class <- factor(subject$call_class, levels = rev(levels(subject$call_class)))
  for (ignore in c(FALSE, TRUE)) for (adjacent in c(FALSE, TRUE)) {
    frozen <- serialize(list(query, subject), NULL)
    check(paste("helper", name, ignore, adjacent, sep = ":"),
      oracle(query, subject, allele_key_names, ignore, adjacent),
      candidate(query, subject, allele_key_names, ignore, adjacent),
      require_success = !grepl("NA_level", name, fixed = TRUE))
    stopifnot(identical(frozen, serialize(list(query, subject), NULL)))
  }
}
query <- make_gr(a); subject <- make_gr(a, query = FALSE)
for (keys in list(character(), c("run_id", "zm"))) {
  for (ignore in c(FALSE, TRUE)) for (adjacent in c(FALSE, TRUE)) {
    check(paste("other-keys", paste(keys, collapse = ","), ignore, adjacent, sep = ":"),
      oracle(query, subject, keys, ignore, adjacent),
      candidate(query, subject, keys, ignore, adjacent))
  }
}
results <- do.call(rbind, rows)
stopifnot(nrow(results) == 60L, !anyDuplicated(results$label))
write.table(results, file.path(output, "comparisons.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
completion <- list(status = if (all(results$exact)) "PASS" else "FAIL", checks = nrow(results),
  overlap_call_replacements = replacements, performance_claim = FALSE,
  implementation = "Custom R with existing compiled vctrs/Bioconductor internals; no custom non-R scientific code")
jsonlite::write_json(completion, file.path(output, "completion.json"), auto_unbox = TRUE, pretty = TRUE)
cat(completion$status, ": ", nrow(results), " exact allele-factor/overlap cases\n", sep = "")
quit(status = if (all(results$exact)) 0L else 1L)
