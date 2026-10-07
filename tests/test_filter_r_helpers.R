#!/usr/bin/env Rscript
# Compare the actual two function definitions with a frozen pre-change script.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
suppressPackageStartupMessages({
  library(GenomicAlignments); library(GenomicRanges); library(plyranges)
  library(tidyverse)
})
options(warn=2)
dir.create(args[[3]], recursive=TRUE, showWarnings=FALSE)
read_helpers <- function(path){
  env <- new.env(parent=globalenv())
  nodes <- as.list(parse(path))
  for(name in c("min_threshold_eachstrand", "convert_query_to_refspace")){
    found <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                      identical(x[[2]], as.name(name)), nodes)
    stopifnot(length(found) == 1L)
    eval(found[[1]], env)
  }
  env
}
original <- read_helpers(args[[1]])
candidate <- read_helpers(args[[2]])
capture <- function(expression) tryCatch(list(value=force(expression)),
  error=function(error) list(error=conditionMessage(error)))
same <- function(before, after, label, allow_errors=FALSE){
  if(allow_errors && !is.null(before$error) && !is.null(after$error)) return(invisible(NULL))
  if(!allow_errors && (!is.null(before$error) || !is.null(after$error))){
    print(before); print(after); stop("Expected successful helper results: ",label)
  }
  if(!identical(before,after)){
    saveRDS(list(before=before,after=after),file.path(args[[3]],paste0(label,".rds")))
    print(before); print(after); print(all.equal(before,after))
    stop("Helper mismatch: ",label)
  }
}
vectors <- list(integer(),numeric(),0L,1L,-1L,NA_integer_,NA_real_,NaN,Inf,-Inf,
  c(0L,1L),c(1,NA_real_),c(NA_real_,1),c(0,NA_real_),c(Inf,-Inf),
  c(9007199254740992,1,-9007199254740992),c(1e-300,1e300),
  setNames(c(1,2),c("a","b")))
cases <- 0L
for(mode in c("mean","all","any")) for(threshold in c(-Inf,-1,0,.5,1,Inf,NA_real_,NaN)){
  for(x in vectors) for(y in vectors){
    cases <- cases+1L
    same(capture(original$min_threshold_eachstrand(x,y,threshold,mode)),
         capture(candidate$min_threshold_eachstrand(x,y,threshold,mode)),
         paste0("threshold-",cases))
  }
}
# A failed first strand must not force/evaluate the opposite-strand expression.
for(mode in c("mean","all","any")){
  same(capture(original$min_threshold_eachstrand(0,stop("forced"),1,mode)),
       capture(candidate$min_threshold_eachstrand(0,stop("forced"),1,mode)),
       paste0("short-circuit-",mode))
}
reference <- GRanges(seqinfo=Seqinfo(c("chr2","chr1","unused"),c(1000L,1000L,50L),
                                    isCircular=c(FALSE,TRUE,FALSE),genome="fixture"))
bam <- tibble(seqnames=factor(c("chr1","chr2","chr1","chr2"),
                             levels=c("chr2","chr1","unused")),
  start=c(10L,100L,210L,400L),end=c(19L,111L,221L,409L),
  cigar=c("10M","2S4M2I3M2D3M","3M2D3M3N1M","10M"),
  strand=factor(c("+","-","*","-"),levels=c("+","-","*")),
  run_id=factor(c("runA","runB","runA","runB"),levels=c("runB","runA","unused")),
  zm=c(10L,10L,2L,3L))
sets <- list(
  ordinary=IRangesList(IRanges(c(1,4),c(2,7)),IRanges(c(2,5,8),c(4,8,12)),
                      IRanges(c(1,3,5),c(2,5,7)),IRanges()),
  empty=IRangesList(rep(list(IRanges()),4L)),
  duplicated=IRangesList(IRanges(c(1,1,2),c(3,3,4)),IRanges(),IRanges(),IRanges(3,5)),
  zero_width=IRangesList(IRanges(start=c(0,1,5,11),width=0),IRanges(6,5),IRanges(),IRanges()),
  outside=IRangesList(IRanges(start=c(-2,15),width=c(4,3)),IRanges(),IRanges(),IRanges()),
  named=IRangesList(IRanges(c(1,4),c(2,7),names=c("a","b")),IRanges(),IRanges(),IRanges()))
sets$metadata <- sets$named
for(i in seq_along(sets$metadata)) mcols(sets$metadata[[i]]) <-
  DataFrame(extra=seq_len(length(sets$metadata[[i]])))
for(label in names(sets)){
  same(capture(original$convert_query_to_refspace(sets[[label]],bam,"reference")),
       capture(candidate$convert_query_to_refspace(sets[[label]],bam,"reference")),label,
       allow_errors=label %in% c("zero_width","outside"))
}
for(label in c("empty_bam","reversed","numeric_zm","character_fields")){
  b <- bam; ranges <- sets$ordinary
  if(label == "empty_bam"){b <- bam[0,];ranges <- IRangesList()}
  if(label == "reversed"){b <- bam[4:1,];ranges <- ranges[4:1]}
  if(label == "numeric_zm") b$zm <- as.double(b$zm)
  if(label == "character_fields"){
    b$seqnames <- as.character(b$seqnames); b$strand <- as.character(b$strand)
    b$run_id <- as.character(b$run_id)
  }
  same(capture(original$convert_query_to_refspace(ranges,b,"reference")),
       capture(candidate$convert_query_to_refspace(ranges,b,"reference")),label)
}
writeLines(paste("PASS",cases,"threshold cases, short-circuit and coordinate fixtures"),
           file.path(args[[3]],"PASS.txt"))
cat(readLines(file.path(args[[3]],"PASS.txt")),"\n")
