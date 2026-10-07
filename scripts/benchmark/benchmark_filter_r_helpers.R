#!/usr/bin/env Rscript
# Replay the actual four query-to-reference helper calls on real filter inputs.
# All timing is diagnostic; complete fresh-process filter pairs are authoritative.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(5L,6L))
if(length(args)==6L) stopifnot(args[[6L]] == "coordinates-only")
baseline_bin <- normalizePath(args[[1]],mustWork=TRUE)
candidate <- normalizePath(args[[2]],mustWork=TRUE)
config <- normalizePath(args[[3]],mustWork=TRUE)
extract <- normalizePath(args[[4]],mustWork=TRUE)
dir.create(args[[5]],recursive=TRUE)
setwd(normalizePath(args[[5]],mustWork=TRUE))
Sys.setenv(PATH=paste(baseline_bin,Sys.getenv("PATH"),sep=.Platform$path.sep))
original_nodes <- as.list(parse(file.path(baseline_bin,"filterCalls.R")))
candidate_nodes <- as.list(parse(candidate))
assignment <- function(node,name) is.call(node) && identical(node[[1]],as.name("<-")) &&
  identical(node[[2]],as.name(name))
find_assignment <- function(nodes,name){
  out <- Filter(function(x) assignment(x,name),nodes)
  stopifnot(length(out)==1L);out[[1]]
}
environment <- new.env(parent=globalenv())
environment$parse_args <- function(object,...) optparse::parse_args(object,
  args=c("-c",config,"-f",extract,"-s","70007.1-CordBlood-LIB1","-g","1-22X",
         "-v","strict","-o","unused.qs2"),...)
environment$source <- function(file,local=FALSE,...) base::source(file,local=environment,...)
timings <- list()
check_coordinate <- function(irangeslist.input,bam.input,BSgenome_name.input){
  index <- length(timings)+1L
  # Force the caller's query construction before timing either function. R's
  # lazy promises otherwise charge that shared work only to the first arm.
  force(irangeslist.input); force(bam.input); force(BSgenome_name.input)
  before <- proc.time()
  expected <- original_coordinate(irangeslist.input,bam.input,BSgenome_name.input)
  middle <- proc.time()
  actual <- candidate_coordinate(irangeslist.input,bam.input,BSgenome_name.input)
  after <- proc.time()
  if(!identical(expected,actual)){
    saveRDS(list(expected=expected,actual=actual),paste0("coordinate-failure-",index,".rds"))
    stop("Actual query-to-reference helper mismatch: ",index)
  }
  timings[[index]] <<- data.frame(helper="coordinates",index=index,
    input_reads=nrow(bam.input),query_ranges=sum(elementNROWS(irangeslist.input)),
    output_ranges=length(actual),original_cpu=sum((middle-before)[1:2]),
    candidate_cpu=sum((after-middle)[1:2]))
  write.table(do.call(rbind,timings),"helper-timings.tsv",sep="\t",row.names=FALSE,quote=FALSE)
  cat("PASS actual coordinate helper",index,"reads",nrow(bam.input),"ranges",length(actual),"\n")
  actual
}
for(node in original_nodes){
  eval(node,environment)
  if(assignment(node,"convert_query_to_refspace")){
    original_coordinate <- environment$convert_query_to_refspace
    eval(find_assignment(candidate_nodes,"convert_query_to_refspace"),environment)
    candidate_coordinate <- environment$convert_query_to_refspace
    environment$convert_query_to_refspace <- check_coordinate
  }
  if(assignment(node,"min_qual.fail.gr")) break
}
# bam is unchanged between base-quality and subread coordinate conversion in
# production. Reuse its actual post-molecule-filter rows for all three queries.
for(name in c("frac_subreads_cvg.fail.gr","num_subreads_match.fail.gr","frac_subreads_match.fail.gr")){
  eval(find_assignment(original_nodes,name),environment)
  rm(list=name,envir=environment)
  invisible(gc())
}
stopifnot(length(timings)==4L)
if(length(args)==6L){
  writeLines("PASS all four actual coordinate helper calls; argument construction excluded from helper timings", "PASS.txt")
  quit(status=0L)
}
original_threshold <- environment$min_threshold_eachstrand
eval(find_assignment(candidate_nodes,"min_threshold_eachstrand"),environment)
candidate_threshold <- environment$min_threshold_eachstrand
calls <- environment$calls
settings <- environment$filtergroup_toanalyze_config
# Own/opposite base qualities and integer subread counts use the original
# materialized list columns; compare every actual call in all three modes.
threshold_timings <- list()
for(mode in c("mean","all","any")) for(column in c("qual","sa","sm","sx")){
  threshold <- if(column=="qual") settings$min_qual else settings$min_num_subreads_match
  before <- proc.time()
  expected <- mapply(function(x,y) original_threshold(x,y,threshold,mode),
                     calls[[column]],calls[[paste0(column,".opposite_strand")]],USE.NAMES=FALSE)
  middle <- proc.time()
  actual <- mapply(function(x,y) candidate_threshold(x,y,threshold,mode),
                   calls[[column]],calls[[paste0(column,".opposite_strand")]],USE.NAMES=FALSE)
  after <- proc.time()
  stopifnot(identical(expected,actual))
  threshold_timings[[length(threshold_timings)+1L]] <- data.frame(column=column,mode=mode,
    calls=length(actual),original_cpu=sum((middle-before)[1:2]),candidate_cpu=sum((after-middle)[1:2]))
  write.table(do.call(rbind,threshold_timings),"threshold-timings.tsv",sep="\t",row.names=FALSE,quote=FALSE)
  cat("PASS actual threshold helper",column,mode,"calls",length(actual),"\n")
}
writeLines("PASS all four actual coordinate helper calls and all-call threshold checks", "PASS.txt")
