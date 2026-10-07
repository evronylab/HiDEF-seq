#!/usr/bin/env Rscript
# Run the actual burden context-consumer expressions on paired real count files.
# This checks the changed producer boundary, not an entire burden-task QS2.
args <- commandArgs(TRUE)
stopifnot(length(args) == 5L)
candidate <- normalizePath(args[[1]], mustWork=TRUE)
shared <- normalizePath(args[[2]], mustWork=TRUE)
directories <- vapply(args[3:4], normalizePath, character(1), mustWork=TRUE)
dir.create(args[[5]], recursive=TRUE)
output <- normalizePath(args[[5]], mustWork=TRUE)
suppressPackageStartupMessages({library(Biostrings);library(qs2);library(tidyverse)})
options(warn=2)
source(shared)
nodes <- as.list(parse(candidate))
named <- function(node,name) is.call(node) && identical(node[[1]],as.name("<-")) &&
  identical(node[[2]],as.name(name))
first <- which(vapply(nodes,named,logical(1),"reftnc_plus_strand_by_annotation_row"))
last <- which(vapply(nodes,named,logical(1),"reference_summary"))
stopifnot(length(first)==1L,length(last)==1L,last>first)
consumer <- Filter(function(node) is.call(node) && identical(node[[1]],as.name("<-")),
                   nodes[seq.int(first,last-1L)])
stopifnot(length(consumer)==5L)
paths <- lapply(directories,function(path) sort(list.files(path,pattern="^[0-9]+\\.reftnc_plus_strand\\.tsv$")))
stopifnot(length(paths[[1]])>0L,identical(paths[[1]],paths[[2]]))
for(file in paths[[1]]) stopifnot(identical(readLines(file.path(directories[[1]],file)),
                                          readLines(file.path(directories[[2]],file))))
original_wd <- getwd()
for(orientations in c(FALSE,TRUE)){
  results <- list()
  for(arm in 1:2){
    directory <- file.path(output,paste(arm,orientations,sep="-"));dir.create(directory)
    setwd(directory)
    n <- length(paths[[arm]])
    rows <- tibble(annotation_row_id=seq_len(n),parent_row_id=seq_len(n),
                   bc_orientation="all_bc_orientations",strand=NA_character_,source_row=seq_len(n))
    if(orientations) rows <- bind_rows(rows,tibble(annotation_row_id=n+seq_len(2L*n),
      parent_row_id=n+rep(seq_len(n),each=2L),bc_orientation="A-B",strand=rep(c("+","-"),n),
      source_row=rep(seq_len(n),each=2L)))
    for(i in seq_len(nrow(rows))){
      lines <- readLines(file.path(directories[[arm]],paths[[arm]][rows$source_row[[i]]]))
      lines <- sub("^[^\t]+",as.character(rows$annotation_row_id[[i]]),lines)
      writeLines(lines,paste0(rows$annotation_row_id[[i]],".reftnc_plus_strand.tsv"))
    }
    env <- new.env(parent=globalenv())
    env$coverage_annotation_rows <- select(rows,-source_row)
    env$coverage_rows <- distinct(rows,row_id=parent_row_id,bc_orientation) %>%
      mutate(bam.gr.filtertrack.by_bc_orientation_strand.coverage=rep(list(NULL),n()))
    for(node in consumer) eval(node,env)
    result <- list(coverage_rows=env$coverage_rows,
      annotation_counts=env$reftnc_plus_strand_by_annotation_row,
      plus_counts=env$reftnc_plus_strand_by_row,
      template_counts=env$reftnc_template_strand_by_row)
    qs_save(result,"consumer-summary.qs2")
    results[[arm]] <- qs_read("consumer-summary.qs2")
    stopifnot(identical(result,results[[arm]]))
  }
  stopifnot(identical(results[[1]],results[[2]]))
  setwd(original_wd)
}
writeLines("PASS ordered actual burden context-consumer tables, factors, fractions and QS2 round trips; aggregate and asymmetric-strand cases",
           file.path(output,"PASS.txt"))
cat(readLines(file.path(output,"PASS.txt")),"\n")
