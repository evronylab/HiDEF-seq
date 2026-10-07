#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
options(warn = 2)
data.table::setDTthreads(1L)
args <- commandArgs(trailingOnly = TRUE)
repo <- if(length(args)) args[[1L]] else "."
nodes <- as.list(parse(file.path(repo, "bin/filterCalls.R")))
loops <- Filter(function(node) is.call(node) && identical(node[[1L]], as.name("for")) &&
  identical(node[[3L]], quote(seq_len(nrow(region_read_filters_config)))), nodes)
stopifnot(length(loops) == 1L)
assignments <- Filter(function(node) is.call(node) && identical(node[[1L]], as.name("<-")) &&
  identical(node[[2L]], as.name("bam")), as.list(loops[[1L]][[4L]])[-1L])
stopifnot(length(assignments) == 1L)

# Feed controlled fractions at the boundary of the unchanged binnedAverage
# operation. Execute the actual remaining production grouping and ordered join.
replacements <- 0L
replace_fractions <- function(node) {
  if(!is.call(node)) return(node)
  if(identical(node[[1L]], as.name("%>%")) && is.call(node[[3L]]) &&
     identical(node[[3L]][[1L]], as.name("binnedAverage"))) {
    replacements <<- replacements + 1L
    return(as.name("fractions"))
  }
  for(i in seq_along(node)[-1L]) node[[i]] <- replace_fractions(node[[i]])
  node
}
annotation <- replace_fractions(assignments[[1L]])
stopifnot(replacements == 1L)
legacy <- function(fractions, bam, type, value) {
  flags <- fractions %>% group_by(run_id, zm) %>%
    summarize(passfilter = !case_when(
      !!type == "gt" ~ mean(frac) > value,
      !!type == "gte" ~ mean(frac) >= value,
      !!type == "lt" ~ mean(frac) < value,
      !!type == "lte" ~ mean(frac) <= value), .groups = "drop")
  left_join(bam, flags, by = join_by(run_id, zm))
}
cases <- list(
  ordinary = tibble(run_id = factor(c("b","a","b","a"), levels=c("b","a","unused")),
    zm=c(2L,1L,2L,1L), frac=c(.1,.5,.3,.5)),
  missing_keys = tibble(run_id=factor(c("a",NA,"a",NA)), zm=c(1L,2L,1L,2L), frac=c(0,1,.5,NA_real_)),
  nonfinite = tibble(run_id=factor(rep("a",8)), zm=rep(1:4,each=2), frac=c(NA,NaN,Inf,Inf,-Inf,Inf,-0,0)),
  empty = tibble(run_id=factor(levels=c("b","a")), zm=integer(), frac=numeric()),
  numeric_edges = tibble(run_id=factor(rep("a",6)), zm=rep(1:3,each=2), frac=c(.5,.5,.5-2^-53,.5,.5+2^-52,.5)),
  integer_fraction = tibble(run_id=factor(c("a","a","b")), zm=c(1L,1L,2L), frac=c(0L,1L,1L)))
cases$character_keys <- cases$ordinary %>% mutate(run_id=as.character(run_id))
checks <- 0L
for(fractions in cases) for(type in c("gt","gte","lt","lte")) for(value in c(-Inf,0,.5,1,Inf,NA_real_)) {
  bam <- fractions %>% select(-frac) %>% mutate(order=seq_len(n()), payload=map(order, ~c(.x,.x)))
  expected <- legacy(fractions, bam, type, value)
  env <- new.env(parent=globalenv())
  env$fractions <- fractions
  env$bam <- bam
  env$read_threshold_type <- type
  env$read_threshold_value <- value
  env$passfilter_label <- "passfilter"
  eval(annotation, env)
  stopifnot(identical(env$bam, expected), identical(env$fractions, fractions))
  checks <- checks + 1L
}
cat("PASS:", checks, "exact production region aggregation and ordered join cases\n")
