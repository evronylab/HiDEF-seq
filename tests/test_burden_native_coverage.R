#!/usr/bin/env Rscript
# Run in the pipeline R environment on a compute node:
# Rscript tests/test_burden_native_coverage.R [repository root] [candidate calculateBurdens.R]
# The 70 checks preserve the exact case bodies that passed isolated job19365179.
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(plyranges)
  library(tidyverse)
})
options(warn = 2)
args <- commandArgs(TRUE)
repo <- if(length(args)) args[[1]] else "."
candidate_source <- if(length(args) >= 2L) args[[2]] else file.path(repo, "bin", "calculateBurdens.R")
baseline <- new.env(parent = globalenv())
sys.source(file.path(repo, "tests", "fixtures", "burden_native_reference_ea8d06f.R"), envir = baseline)
load_candidate <- function(path) {
  env <- new.env(parent = globalenv())
  wanted <- c("sum_RleList", "bc_orientation_is_asymmetric", "validate_bam.gr.filtertrack",
              "calc_duplex_coverage", "accumulate_bam.gr.filtertracks", "filtertrack_coverage_result",
              "gr_1bp_cov", "sum_filtertrack_sensitivity_coverage", "duplex_plus_ranges",
              "make_chunk_coverage_kernels", "sensitivity_site_counts")
  for(expr in parse(path)) {
    if(is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
       as.character(expr[[2]]) %in% wanted) eval(expr, env)
  }
  stopifnot(all(wanted %in% ls(env)))
  env
}
candidate <- load_candidate(candidate_source)
report_directory <- tempfile("burden-native-regressions-")
dir.create(report_directory)

# Ordinary coverage, original row joins, factors and orientation sums (14).
local({
  args <- c("", file.path(report_directory, "ordinary.tsv"))
  report_file <- args[[2]]
si<-Seqinfo(c("a","b"),seqlengths=c(20L,20L))
gr<-GRanges(c("a","a"),IRanges(c(2L,2L),c(6L,6L)),strand=c("+","-"),seqinfo=si,run_id=c("r","r"),zm=c(1L,1L),bc_orientation=factor(c("x-x","x-x")))
table<-function(x,keys=c("one","two"))tibble(call_type=keys,call_class="SBS",SBSindel_call_type="mutation",filtergroup=factor("strict"),bam.gr.filtertrack=rep(list(x),length(keys)))
capture<-function(fun)tryCatch(fun(),error=function(e)conditionMessage(e))
run<-function(chunks,intern=FALSE) {
  state<-new.env(parent=emptyenv())
  for(x in chunks) {
    orientation<-if(any(x$SBSindel_call_type=="mismatch-ss")) x %>% select(-bam.gr.filtertrack) else NULL
    if(intern) {
      candidate$accumulate_bam.gr.filtertracks(state,x,orientation)
    } else baseline$accumulate_bam.gr.filtertracks(state,x,orientation)
  }
  baseline$filtertrack_coverage_result(state)
}
cases<-list(ordinary=list(table(gr),table(gr)),empty=list(table(gr[0]),table(gr[0])),
  repeated_keys=list(table(gr,c("one","one")),table(gr,c("one","one"))),
  reordered_keys=list(table(gr),table(gr,c("two","one"))),
  unjoined_category=list(table(gr,"one"),table(gr,c("one","new"))))
other<-gr;seqnames(other)[2]<-"b"
cases$seqnames_odd<-list(table(other));cases$seqnames_even<-list(table(c(other,other)))
star<-gr[1];strand(star)<-"*"
cases$star_odd<-list(table(c(gr,star)));cases$star_even<-list(table(c(gr,star,star)))
cases$late_fallback<-list(table(gr),table(c(gr,star,star)))
unknown<-gr;seqlengths(unknown)<-rep(NA_integer_,length(seqlevels(unknown)));cases$unknown_lengths<-list(table(unknown),table(unknown))
circular<-gr;isCircular(circular)[1]<-TRUE;cases$circular<-list(table(circular))
asymmetric<-gr;asymmetric$bc_orientation<-factor(c("x-y","y-x"));cases$asymmetric<-list(table(asymmetric))
orientation<-table(asymmetric);orientation$SBSindel_call_type<-"mismatch-ss"
cases$asymmetric_orientation_sums<-list(orientation,orientation)
results<-list()
for(label in names(cases)) {
  original<-suppressWarnings(capture(function()run(cases[[label]])))
  if(label=="ordinary") stopifnot(is.data.frame(original),nrow(original)==2L)
  actual<-suppressWarnings(capture(function()run(cases[[label]],TRUE)))
  ok<-identical(original,actual)
  results[[length(results)+1L]]<-data.frame(case=label,identical=ok,outcome=if(is.character(original))original else "exact result")
  write.table(do.call(rbind,results),args[2],sep="\t",row.names=FALSE,quote=FALSE)
  if(!ok){print(original);print(actual);stop(label)}
}

})

# Original sensitivity valid/invalid inputs and warning behavior (24).
local({
  args <- c("", file.path(report_directory, "sensitivity.tsv"))
  report_file <- args[[2]]
si <- Seqinfo(c("a","b"),seqlengths=c(20L,20L))
gr <- GRanges(c("a","a"),IRanges(c(2L,2L),c(6L,6L)),strand=c("+","-"),
  seqinfo=si,run_id=c("r","r"),zm=c(1L,1L),bc_orientation=factor(c("x-x","x-x")))
q <- GRanges(c("a","a","b"),IRanges(c(1L,4L,3L),width=1L),seqinfo=si)
cases <- list()
add <- function(name,track=gr,start=q,end=rev(q)) {
  cases[[name]] <<- list(track=track,start=start,end=end)
}
add("ordinary")
for(s in c("+","-","*")) {x<-q;strand(x)<-s;add(paste0("query_strand_",s),start=x,end=rev(x))}
add("empty_queries",start=q[0],end=q[0]);add("empty_track",track=gr[0])
outside <- GRanges("unknown",IRanges(2L,width=1L),seqinfo=Seqinfo("unknown",20L))
add("absent_sequence",start=outside,end=outside)
x<-q;seqlengths(x)<-c(21L,20L);add("query_length_conflict",start=x,end=x)
x<-q;genome(x)<-"other";z<-gr;genome(z)<-"reference";add("query_genome_conflict",track=z,start=x,end=x)
add("query_sides_differ",start=q,end=x)
x<-q;width(x)<-2L;add("wide_query",start=x)
x<-q;width(x)<-0L;add("zero_width_query",start=x)
x<-q;suppressWarnings(ranges(x)[1]<-IRanges(0L,width=1L));add("below_reference",start=x)
x<-q;suppressWarnings(ranges(x)[1]<-IRanges(21L,width=1L));add("above_reference",start=x)
x<-gr;seqlengths(x)<-rep(NA_integer_,2L);add("unknown_lengths",track=x)
x<-gr;isCircular(x)[1]<-TRUE;add("circular",track=x)
x<-gr;seqnames(x)[2]<-"b";add("seqnames_odd",track=x);add("seqnames_even",track=c(x,x))
x<-gr[1];strand(x)<-"*";add("star_odd",track=c(gr,x));add("star_even",track=c(gr,x,x))
x<-gr;width(x)<-0L;add("zero_width_track",track=x)
x<-gr;end(x)[2]<-7L;add("mismatched_pair",track=x)
capture <- function(fn) {
  warnings <- character()
  result <- tryCatch(withCallingHandlers(list(value=fn()),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  }),error=function(e)list(error=conditionMessage(e)))
  list(result=result,warnings=warnings)
}
rows <- list()
for(name in names(cases)) {
  x <- cases[[name]]
  expected <- capture(function() {
    baseline$validate_bam.gr.filtertrack(x$track)
    cov <- baseline$calc_duplex_coverage(x$track)
    list(start=baseline$gr_1bp_cov(x$start,cov),end=baseline$gr_1bp_cov(x$end,cov))
  })
  actual <- capture(function()candidate$sensitivity_site_counts(x$track,x$start,x$end))
  ok <- identical(expected,actual)
  rows[[length(rows)+1L]] <- tibble(case=name,identical=ok,
    outcome=if(is.null(expected$result$error))"exact result and warnings" else expected$result$error)
  write_tsv(bind_rows(rows),args[2])
  if(!ok){print(expected);print(actual);stop(name)}
}
bytype <- tibble(call_type="SBS",call_class="SBS",SBSindel_call_type="mutation",
  filtergroup=factor("strict"),bam.gr.filtertrack=list(gr))
for(label in c("disabled","unmatched")) {
  queries <- if(label=="disabled")NULL else tibble(call_type="insertion",query_start=list(q),query_end=list(q))
  expected <- baseline$sum_filtertrack_sensitivity_coverage(bytype,queries)
  actual <- candidate$sum_filtertrack_sensitivity_coverage(bytype,queries)
  stopifnot(identical(expected,actual))
  rows[[length(rows)+1L]] <- tibble(case=label,identical=TRUE,outcome="exact result")
}
write_tsv(bind_rows(rows),args[2])

})

# Circularity and unused metadata at both helper/wrapper levels, warn0/warn2 (32).
local({
  args <- c("", file.path(report_directory, "review-regressions.tsv"))
  report_file <- args[[2]]
capture <- function(fn, warn_as_error = FALSE) {
  old <- options(warn = if(warn_as_error) 2L else 0L)
  on.exit(options(old))
  warnings <- character()
  result <- tryCatch(withCallingHandlers(list(value = fn()), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    if(!warn_as_error) invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e)))
  list(result = result, warnings = warnings)
}

track <- function(start = 18L, end = 23L, circular = NA) {
  suppressWarnings(GRanges(c("a", "a"), IRanges(rep(start, 2), rep(end, 2)),
    strand = c("+", "-"), seqinfo = Seqinfo("a", 20L, isCircular = circular),
    run_id = c("r", "r"), zm = c(1L, 1L), bc_orientation = factor(c("x-x", "x-x"))))
}
query <- function(circular = TRUE) GRanges(rep("a", 5L),
  IRanges(c(1L, 2L, 3L, 18L, 20L), width = 1L),
  seqinfo = Seqinfo("a", 20L, isCircular = circular))
cases <- list(
  query_true_track_na_outside = list(gr = track(), qs = query(), qe = query()),
  query_true_track_na_inside = list(gr = track(end = 20L), qs = query(), qe = query()),
  query_na_track_na_outside = list(gr = track(), qs = query(NA), qe = query(NA)),
  query_false_track_na_outside = list(gr = track(), qs = query(FALSE), qe = query(FALSE)),
  query_true_track_false_outside = list(gr = track(circular = FALSE), qs = query(), qe = query()),
  query_true_track_true_outside = list(gr = track(circular = TRUE), qs = query(), qe = query()),
  query_true_track_na_below = list(gr = track(start = -2L, end = 3L), qs = query(), qe = query())
)
qs_meta <- query(NA)
qe_meta <- query(NA)
qs_meta$unused <- IRanges(seq_along(qs_meta), width = 1L)
qe_meta$unused <- letters[seq_along(qe_meta)]
cases$query_mcols_incompatible <- list(gr = track(end = 20L), qs = qs_meta, qe = qe_meta)
results <- list()
for(label in names(cases)) for(warn_as_error in c(FALSE, TRUE)) {
  x <- cases[[label]]
  bytype <- tibble(call_type="SBS",call_class="SBS",SBSindel_call_type="mutation",
    filtergroup=factor("strict"),bam.gr.filtertrack=list(x$gr))
  queries <- tibble(call_type="SBS",query_start=list(x$qs),query_end=list(x$qe))
  for(scope in c("site_counts", "whole_sensitivity_wrapper")) {
    expected <- capture(function() {
      if(scope == "whole_sensitivity_wrapper") return(baseline$sum_filtertrack_sensitivity_coverage(bytype,queries))
      baseline$validate_bam.gr.filtertrack(x$gr)
      cov <- baseline$calc_duplex_coverage(x$gr)
      list(start=baseline$gr_1bp_cov(x$qs,cov),end=baseline$gr_1bp_cov(x$qe,cov))
    }, warn_as_error)
    actual <- capture(function() {
      if(scope == "whole_sensitivity_wrapper") return(candidate$sum_filtertrack_sensitivity_coverage(bytype,queries))
      candidate$sensitivity_site_counts(x$gr,x$qs,x$qe)
    }, warn_as_error)
    stopifnot(is.null(expected$result$error), !is.null(expected$result$value))
    ok <- identical(expected, actual)
    results[[length(results)+1L]] <- data.frame(case=label,warn_as_error=warn_as_error,scope=scope,identical=ok)
    write.table(do.call(rbind,results),report_file,sep="\t",row.names=FALSE,quote=FALSE)
    if(!ok) {print(expected);print(actual);stop(paste(label,scope,warn_as_error))}
  }
}
cat("PASS: 16 reviewed guard cases, each helper and whole wrapper identical\n")

})

counts <- vapply(list.files(report_directory, full.names = TRUE), function(path) nrow(read.delim(path)), integer(1))
stopifnot(sum(counts) == 70L)
cat("PASS: 70 native coverage/sensitivity regression checks, exact original reference\n")
unlink(report_directory, recursive = TRUE)
