#!/usr/bin/env Rscript
# Run in the pipeline container: Rscript tests/test_indel_tag_indices.R [repository root]
args <- commandArgs(TRUE)
suppressPackageStartupMessages({library(S4Vectors); library(data.table)})
options(warn=2)
repo <- if(length(args)) args[[1]] else "."
source_file <- file.path(repo, "bin/extractCalls.R")
expressions <- parse(source_file)
for (expr in expressions) {
  if(is.call(expr) && identical(expr[[1]], as.name("<-")) &&
     identical(expr[[2]], as.name("subset_tag_positions"))) eval(expr)
}
lines <- readLines(source_file)
actual <- lines[grepl("indels_queryspace_pos[", lines, fixed=TRUE) |
                grepl("indels_queryspace_pos.opposite_strand[", lines, fixed=TRUE)]
actual <- actual[grepl(":=", actual, fixed=TRUE)]
stopifnot(length(actual)==12L)
legacy <- actual[grepl("_val :=", actual, fixed=TRUE)]
legacy <- gsub("tag_index[1]", "zm_strand[1]", legacy, fixed=TRUE)
stopifnot(length(legacy)==6L)
run <- function(code, tags, query) {
  env <- new.env(parent=globalenv())
  env$sa.input <- lapply(tags, Rle)
  env$sm.input <- tags
  env$sx.input <- lapply(tags, function(x) x * 0)
  for (tag in c("sa", "sm", "sx")) {
    value <- get(paste0(tag, ".input"), env)
    assign(paste0(tag, ".input.opposite_strand"),
           setNames(value, chartr("+-", "-+", names(value))), env)
  }
  env$indels_queryspace_pos <- copy(query)
  env$indels_queryspace_pos.opposite_strand <- copy(query)[, zm_strand := chartr("+-", "-+", zm_strand)]
  for (table in c("indels_queryspace_pos", "indels_queryspace_pos.opposite_strand")) {
    setkeyv(get(table, env), c("zm_strand", "start_end_queryspace"))
  }
  eval(parse(text=code), env)
  lapply(c("indels_queryspace_pos", "indels_queryspace_pos.opposite_strand"), function(table) {
    value <- get(table, env)
    list(columns=as.list(value), key=key(value))
  })
}
capture <- function(expr) tryCatch(list(value=force(expr)), error=function(e) list(error=conditionMessage(e)))
tags <- setNames(list(1:6, 6:1, c(1, 2.5, NA_real_, 4, 5, 6), rep(99L, 6), rep(88L, 6), rep(77L, 6)),
                 c("101_+", "101_-", "102_+", "101_+", "", NA_character_))
query <- data.table(zm_strand=c("102_+", "101_+", "101_+", "101_-", "102_+", "101_-"),
                    start_end_queryspace=c("1_2", "2_4", "2_4", "1_6", "1_2", "1_6"),
                    pos_queryspace=c(1L, 2L, 4L, 1L, 2L, 6L))
for (valid_tags in list(lapply(tags, as.numeric), lapply(tags, as.integer))) {
  stopifnot(identical(run(legacy, valid_tags, query), run(actual, valid_tags, query)))
}
for (vector_tags in list(lapply(tags, as.numeric), lapply(tags, as.integer), tags)) {
 for (q in list(query, query[0], copy(query)[, pos_queryspace := c(NA_integer_, 9L, 2L, 3L, 0L, 6L)],
               copy(query)[, zm_strand := c("", NA_character_, "missing", "101_+", "101_+", "101_+")])) {
  old <- capture(run(legacy, vector_tags, q))
  new <- capture(run(actual, vector_tags, q))
  if(!identical(old, new)) {str(old); str(new); stop("indel tag index differs")}
 }
}
# Confirm the matching rule also retains base-list first-match and missing-key semantics.
keys <- c(names(tags), "missing")
indices <- match(keys, names(tags))
indices[is.na(keys) | keys == ""] <- NA_integer_
stopifnot(identical(lapply(keys, function(key) tags[[key]]), lapply(indices, function(index) tags[[index]])))
cat("PASS: actual six indel assignments preserve values, types, keys, ordering, duplicate/absent names and empty queries\n")
