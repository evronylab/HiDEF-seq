#!/usr/bin/env Rscript
#Validation-only rules must never mask scientific changes or arbitrary paths.
args <- commandArgs(trailingOnly=TRUE)
repo <- if(length(args)) args[[1]] else "."
source(file.path(repo,"scripts","benchmark","compare_qs2.R"))
rule <- list(mode="exact-path",reference=c("/old/a.bw","/old/b.bw"),
             candidate=c("/new/A/a.bw","/new/B/b.bw"))
path <- "/region_genome_filter_stats/region_filter_threshold_file"
rules <- setNames(list(rule),path)
a <- data.frame(filter=c("none","a","b","empty"),
                region_filter_threshold_file=c(NA,rule$reference,""),count=c(0L,10L,20L,0L))
attr(a$region_filter_threshold_file,"provenance_attr") <- "retained"
b <- a
b$region_filter_threshold_file[2:3] <- rule$candidate
passes <- function(candidate,policy=rules,component="region_genome_filter_stats") {
  tryCatch({
    left <- normalize_provenance_column(a,component,policy,"reference")
    right <- normalize_provenance_column(candidate,component,policy,"candidate")
    comparison <- scientific_compare(left,right,ignore=character())
    !comparison$counts[["failure"]] && !comparison$counts[["review"]]
  },error=function(e) FALSE)
}
stopifnot(passes(b))
mutations <- list()
x <- b; x$region_filter_threshold_file[[2]] <- rule$candidate[[2]]; mutations$wrong_product <- x
x <- b; x$region_filter_threshold_file[[2]] <- "/new/foreign/a.bw"; mutations$foreign <- x
x <- b; x$count[[2]] <- 11L; mutations$count <- x
x <- b; x$region_filter_threshold_file[[2]] <- NA; mutations$missing <- x
x <- b; x$region_filter_threshold_file[[1]] <- rule$candidate[[1]]; mutations$added <- x
x <- b; attr(x$region_filter_threshold_file,"provenance_attr") <- "changed"; mutations$attribute <- x
x <- b; x$region_filter_threshold_file <- factor(x$region_filter_threshold_file); mutations$factor <- x
x <- b; x$region_filter_threshold_file[[2]] <- paste(rule$candidate,collapse=","); mutations$joined <- x
mutations$reordered <- b[c(1,3,2,4),]
mutations$deleted <- b[-2,]
stopifnot(!any(vapply(mutations,passes,logical(1))))
bad_mode <- rules; bad_mode[[path]]$mode <- "tokens"; stopifnot(!passes(b,bad_mode))
bad_scope <- setNames(list(rule),"/finalCalls/region_filter_threshold_file")
stopifnot(!passes(b,bad_scope))
stopifnot(identical(normalize_cache_paths(c(NA,NA),rule,"candidate"),c(NA,NA)))
#The existing comma-separated germline provenance format remains unchanged.
tokens <- list(reference=c("/old/DV","/old/Clair"),candidate=c("/new/DV","/new/Clair"))
stopifnot(identical(normalize_provenance_tokens(c(NA,"","/new/Clair,/new/Clair,/new/DV"),tokens,"candidate"),
                    c(NA,"","/old/Clair,/old/Clair,/old/DV")))
cat("PASS: exact region paths; wrong products, scientific values, attributes, NA/order and scope rejected; germline tokens unchanged\n")
