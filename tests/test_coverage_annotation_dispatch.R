#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4L)
candidate <- normalizePath(args[[1]], mustWork=TRUE)
shared <- normalizePath(args[[2]], mustWork=TRUE)
fixture <- normalizePath(args[[3]], mustWork=TRUE)
dir.create(args[[4]], recursive=TRUE, showWarnings=FALSE)
output <- normalizePath(args[[4]], mustWork=TRUE)
suppressPackageStartupMessages({library(tidyverse); library(optparse)})
options(warn=2)
nodes <- as.list(parse(candidate))
assignment <- function(name){
  found <- Filter(function(n) is.call(n) && as.character(n[[1]]) %in% c("<-", "=") &&
    identical(n[[2]], as.name(name)), Filter(function(n) is.call(n) && is.symbol(n[[1]]), nodes))
  stopifnot(length(found) == 1L)
  found[[1]]
}
shared_nodes <- as.list(parse(shared))
for(node in shared_nodes){
  if(is.call(node) && identical(node[[1]], as.name("<-")) &&
     as.character(node[[2]]) %in% c("coverage_annotation_index", "use_r_coverage_annotation", "annotate_coverage_row")) eval(node)
}
eval(assignment("option_list"))
parsed <- parse_args(OptionParser(option_list=option_list), args=c("--coverage-annotation", "legacy"))
stopifnot(identical(parsed$coverage_annotation, "legacy"))
stopifnot(identical(parse_args(OptionParser(option_list=option_list), args=character())$coverage_annotation, "r"))
dispatch <- Filter(function(n) is.call(n) && identical(n[[1]], as.name("if")) &&
  is.call(n[[2]]) && identical(n[[2]][[1]], as.name("use_r_coverage_annotation")), nodes)
stopifnot(length(dispatch) == 1L)
dispatch <- dispatch[[1]]
fasta <- file.path(fixture, "reference.fa")
fai <- paste0(fasta, ".fai")
stopifnot(use_r_coverage_annotation(fasta, fai), !use_r_coverage_annotation(fasta, fai, "legacy"))
expect_error <- function(expression){
  err <- tryCatch({force(expression); NULL}, error=identity)
  stopifnot(inherits(err, "error"))
  conditionMessage(err)
}
expect_error(use_r_coverage_annotation(file.path(output,"missing"), fai))
expect_error(use_r_coverage_annotation(fasta, fai, "unknown"))
read_counts <- function(path){
  x <- read.delim(path, header=FALSE, col.names=c("row", "context", "count"),
    colClasses=c("character", "character", "numeric"), na.strings=character())
  x <- x[order(x$row, x$context), ]; rownames(x) <- NULL; x
}
read_bed <- function(path){con <- gzfile(path); on.exit(close(con)); readLines(con)}
for(row in 1:4){
  # The Python legacy fixture covers all contexts, original depth tokens,
  # noninteger/large depths, small contigs, chrM edges and empty input.
  for(aggregate in c(TRUE,FALSE)){
    directory <- file.path(output,paste(row,aggregate,sep="-")); dir.create(directory)
    input <- file.path(directory,"input.bed")
    stopifnot(file.copy(file.path(fixture,paste0(row,".bed")),input))
    counts <- file.path(directory,"counts.tsv")
    bed <- if(aggregate) file.path(directory,"result.bed.gz") else NULL
    annotate_coverage_row(input,fasta,fai,row,counts,bed,"bgzip","tabix", window_bases=7L,input_rows=2L)
    stopifnot(!file.exists(input), identical(read_counts(counts),
      read_counts(file.path(fixture,paste0(row,".legacy.counts.tsv")))))
    if(aggregate) stopifnot(file.exists(paste0(bed,".tbi")),
      identical(read_bed(bed),read_bed(file.path(fixture,paste0(row,".legacy.bed.gz")))))
    else stopifnot(length(list.files(directory,pattern="bed.gz")) == 0L)
  }
}
# A valid FAI may omit its final newline; options(warn=2) must not reject it.
no_lf <- file.path(output,"reference-no-final-lf.fai")
writeBin(charToRaw(paste(readLines(fai),collapse="\n")),no_lf)
stopifnot(use_r_coverage_annotation(fasta,no_lf))
no_lf_bed <- file.path(output,"no-final-lf.bed")
stopifnot(file.copy(file.path(fixture,"1.bed"),no_lf_bed))
no_lf_counts <- file.path(output,"no-final-lf.counts.tsv")
annotate_coverage_row(no_lf_bed,fasta,no_lf,1L,no_lf_counts)
stopifnot(identical(read_counts(no_lf_counts),read_counts(file.path(fixture,"1.legacy.counts.tsv"))))
# Exercise the actual conditional and untouched legacy block, including an
# unsafe UNCOVERED contig; its length 1 yields no legacy trinucleotide rows.
fallback <- file.path(output,"fallback"); dir.create(fallback)
oldwd <- setwd(fallback)
writeLines(c(">chr1","ACGTACGT",">chr-uncovered","A"), "reference.fa")
stopifnot(system2("samtools", c("faidx", "reference.fa")) == 0L)
writeLines(c("chr1\t1\t2\tACG","chr1\t2\t3\tCGT","chr1\t3\t4\tGTA",
  "chr1\t4\t5\tTAC","chr1\t5\t6\tACG","chr1\t6\t7\tCGT"), "reference.trinuc.bed")
stopifnot(system2("bgzip",c("-c","reference.trinuc.bed"),stdout="reference.trinuc.bed.gz") == 0L)
yaml.config <- list(genome_fasta="reference.fa",genome_fai="reference.fa.fai",
  bedtools_bin="bedtools",bgzip_bin="bgzip",tabix_bin="tabix",cache_dir=".")
cache_file <- function(...) "reference.trinuc.bed.gz"
coverage_annotation_rows <- tibble(annotation_row_id=1L,bc_orientation="all_bc_orientations",
  call_class="SBS",call_type="SBS",SBSindel_call_type="mutation")
get_coverage_reftnc_output_file <- function(...) "result.bed.gz"
writeLines(c("#!/bin/bash", "touch invoked-unexpectedly", "exit 19"),"forbidden-helper")
Sys.chmod("forbidden-helper", "0755")
opt <- list(coverage_annotation="r")
stopifnot(!use_r_coverage_annotation(yaml.config$genome_fasta,yaml.config$genome_fai,opt$coverage_annotation))
writeLines("chr1\t0\t8\t2","1.bed")
eval(dispatch)
stopifnot(!file.exists("invoked-unexpectedly"),file.exists("result.bed.gz.tbi"))
expected <- c("chr1\t0\t1\t2\t.",paste0("chr1\t",1:6,"\t",2:7,"\t2\t",c("ACG","CGT","GTA","TAC","ACG","CGT")),"chr1\t7\t8\t2\t.")
stopifnot(identical(read_bed("result.bed.gz"),expected))
fallback_counts <- read_counts("1.reftnc_plus_strand.tsv")
opt$coverage_annotation <- "legacy"
file.remove("result.bed.gz","result.bed.gz.tbi","1.reftnc_plus_strand.tsv")
writeLines("chr1\t0\t8\t2","1.bed")
eval(dispatch)
stopifnot(identical(read_bed("result.bed.gz"),expected),identical(read_counts("1.reftnc_plus_strand.tsv"),fallback_counts))
setwd(oldwd)

# Actual production conditional must propagate late helper, compressor and
# indexer failures before a downstream save sentinel. Source BED is retained.
for(kind in c("late-aggregate","late-orientation","compressor","indexer")){
  directory <- file.path(output,kind);dir.create(directory);setwd(directory)
  input <- if(startsWith(kind,"late")) "chr1\t1\t3\t2\nchr1\t2\t4\t1" else "chr1\t1\t3\t2"
  writeLines(input,"1.bed")
  yaml.config <- list(genome_fasta=fasta,genome_fai=fai,bgzip_bin="bgzip",tabix_bin="tabix")
  if(kind %in% c("compressor","indexer")){
    writeLines(c("#!/bin/bash","cat >/dev/null","exit 17"),"fail-tool")
    Sys.chmod("fail-tool","0755")
    yaml.config[[if(kind == "compressor") "bgzip_bin" else "tabix_bin"]] <- "./fail-tool"
  }
  coverage_annotation_rows$bc_orientation <- if(kind == "late-orientation") "bc1" else "all_bc_orientations"
  opt <- list(coverage_annotation="r")
  error <- expect_error({eval(dispatch); writeLines("invalid success","saved-qs-sentinel")})
  stopifnot(grepl("Coverage",error,fixed=TRUE),
    file.exists("1.bed"),!file.exists("saved-qs-sentinel"),!file.exists("result.bed.gz.tbi"))
  if(startsWith(kind,"late")) stopifnot(!file.exists("1.reftnc_plus_strand.tsv"))
  setwd(oldwd)
}
# Helper paths and all input/output path arguments may contain spaces.
space <- file.path(output,"quoted paths");dir.create(space)
for(name in c("reference.fa","reference.fa.fai","1.bed")) stopifnot(file.copy(file.path(fixture,name),file.path(space,name)))
annotate_coverage_row(file.path(space,"1.bed"),
  file.path(space,"reference.fa"),file.path(space,"reference.fa.fai"),1L,
  file.path(space,"counts table.tsv"),file.path(space,"coverage output.bed.gz"),"bgzip","tabix")
stopifnot(identical(read_bed(file.path(space,"coverage output.bed.gz")),read_bed(file.path(fixture,"1.legacy.bed.gz"))))

# CRLF-wrapped reference and missing final LF exercise byte offsets and raw
# normalization independently of Rsamtools/BSgenome. Vary both block boundaries.
wrapped <- file.path(output, "wrapped"); dir.create(wrapped)
sequence <- paste(rep(c("aCGTrYsW", "NNACgtXA"), 200L), collapse="")
reference <- file.path(wrapped, "reference.fa")
pieces <- substring(sequence, seq(1L, nchar(sequence), 13L),
                    pmin(seq(1L, nchar(sequence), 13L) + 12L, nchar(sequence)))
writeBin(charToRaw(paste0(">chr1\r\n", paste(pieces, collapse="\r\n"))), reference)
stopifnot(system2("samtools", c("faidx", shQuote(reference))) == 0L)
normalized <- strsplit(gsub("[^ACGT]", "N", toupper(sequence)), "", fixed=TRUE)[[1L]]
expected_context <- c(".", vapply(2:(length(normalized)-1L), function(i)
  paste(normalized[(i-1L):(i+1L)], collapse=""), character(1)), ".")
for(window in c(1L, 11L, 1000L)){
  input <- file.path(wrapped, paste0(window, ".bed"))
  writeLines(paste("chr1", 0, nchar(sequence), "02", sep="\t"), input)
  bed <- file.path(wrapped, paste0(window, ".bed.gz"))
  annotate_coverage_row(input,reference,paste0(reference,".fai"),1L,
    paste0(input,".counts"),bed,"bgzip","tabix",window_bases=window)
  stopifnot(identical(read_bed(bed),paste("chr1",0:(nchar(sequence)-1L),
    1:nchar(sequence),"02",expected_context,sep="\t")))
}

# Addition must remain left-to-right within each context across both windows
# and input batches; multiplying run width by depth changes these answers.
precision <- file.path(output,"precision");dir.create(precision)
ref <- file.path(precision,"reference.fa")
writeLines(c(">chr1",paste(rep("A",12003L),collapse="")),ref)
stopifnot(system2("samtools", c("faidx",shQuote(ref))) == 0L)
input_lines <- c("chr1\t1\t2\t9007199254740992", "chr1\t2\t2001\t1",
                 "chr1\t2001\t4001\t0.1", "chr1\t4001\t6001\t3",
                 "chr1\t6001\t12002\t1.5")
awk_input <- file.path(precision,"per-base.tsv")
writeLines(rep(c("9007199254740992","1","0.1","3","1.5"),c(1L,1999L,2000L,2000L,6001L)),awk_input)
expected_file <- file.path(precision,"expected.tsv")
stopifnot(system2("awk",c(shQuote('{s+=$1} END{print s}'),shQuote(awk_input)),stdout=expected_file) == 0L)
expected_sum <- as.numeric(readLines(expected_file))
for(window in c(7L, 2000L, 1000000L)){
  input <- file.path(precision,paste0(window,".bed"));writeLines(input_lines,input)
  counts <- paste0(input,".counts")
  annotate_coverage_row(input,ref,paste0(ref,".fai"),1L,counts,
                        window_bases=window,input_rows=2L)
  actual <- read_counts(counts)
  stopifnot(nrow(actual)==1L,actual$context=="AAA",actual$count==expected_sum)
}

# Compressed FASTA deliberately chooses the existing reference BED fallback.
compressed <- file.path(output,"reference.fa.gz")
con <- gzfile(compressed,"wb");writeBin(readBin(fasta,"raw",file.info(fasta)$size),con);close(con)
stopifnot(!use_r_coverage_annotation(compressed,fai))
empty_contig_index <- file.path(output,"reference-empty-contig.fai")
writeLines(c(readLines(fai), "chrEmpty\t0\t0\t0\t0"),empty_contig_index)
stopifnot(!use_r_coverage_annotation(fasta,empty_contig_index))

writeLines("PASS R coverage annotation, actual production dispatch, legacy fallback and failure propagation",file.path(output,"PASS.txt"))
cat(readLines(file.path(output,"PASS.txt")),"\n")
