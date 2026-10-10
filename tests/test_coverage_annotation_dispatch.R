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
     as.character(node[[2]]) %in% c("coverage_annotation_index", "with_coverage_annotation_reference", "annotate_coverage_row")) eval(node)
}
eval(assignment("option_list"))
stopifnot(!"coverage_annotation" %in% names(parse_args(OptionParser(option_list=option_list), args=character())))
dispatch <- Filter(function(n) is.call(n) &&
  identical(n[[1]], as.name("with_coverage_annotation_reference")), nodes)
stopifnot(length(dispatch) == 1L)
dispatch <- dispatch[[1]]
fasta <- file.path(fixture, "reference.fa")
fai <- paste0(fasta, ".fai")
expect_error <- function(expression){
  err <- tryCatch({force(expression); NULL}, error=identity)
  stopifnot(inherits(err, "error"))
  conditionMessage(err)
}
expect_error(coverage_annotation_index(file.path(output,"missing"), fai))
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
stopifnot(identical(coverage_annotation_index(fasta,no_lf),coverage_annotation_index(fasta,fai)))
no_lf_bed <- file.path(output,"no-final-lf.bed")
stopifnot(file.copy(file.path(fixture,"1.bed"),no_lf_bed))
no_lf_counts <- file.path(output,"no-final-lf.counts.tsv")
annotate_coverage_row(no_lf_bed,fasta,no_lf,1L,no_lf_counts)
stopifnot(identical(read_counts(no_lf_counts),read_counts(file.path(fixture,"1.legacy.counts.tsv"))))
# Production entry point is now a single annotator, independent of a reference
# trinucleotide BED or method selector, for aggregate and orientation rows.
oldwd <- getwd()
coverage_annotation_rows <- tibble(annotation_row_id=1L,bc_orientation="all_bc_orientations",
  call_class="SBS",call_type="SBS",SBSindel_call_type="mutation")
get_coverage_reftnc_output_file <- function(...) "result.bed.gz"
# Actual production entry point must propagate late annotation, compressor and
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

# gzip and BGZF references (with or without a .gzi sidecar) use exactly the
# same ordinary indexed reader after one bounded-memory expansion per pass.
for(format in c("gzip", "bgzf", "bgzf-indexed", "bzip2", "xz")){
  compressed <- file.path(output,paste0("reference-",format,".fa.gz"))
  if(format %in% c("gzip", "bzip2", "xz")){
    writer <- switch(format,gzip=gzfile,bzip2=bzfile,xz=xzfile)
    con <- writer(compressed,"wb");writeBin(readBin(fasta,"raw",file.info(fasta)$size),con);close(con)
  }else{
    stopifnot(system2("bgzip",c("-c",shQuote(fasta)),stdout=compressed) == 0L)
    if(format == "bgzf-indexed") stopifnot(system2("bgzip",c("-r",shQuote(compressed))) == 0L)
  }
  directory <- file.path(output,format);dir.create(directory);setwd(directory)
  yaml.config <- list(genome_fasta=compressed,genome_fai=fai,bgzip_bin="bgzip",tabix_bin="tabix")
  coverage_annotation_rows <- tibble(annotation_row_id=1:4,
    bc_orientation=c("all_bc_orientations","bc1","bc1","bc1"),
    call_class="SBS",call_type="SBS",SBSindel_call_type="mutation")
  for(row in 1:4) stopifnot(file.copy(file.path(fixture,paste0(row,".bed")),paste0(row,".bed")))
  before <- list.files(getwd(),pattern="^coverage-reference-",full.names=TRUE)
  eval(dispatch)
  stopifnot(identical(before,list.files(getwd(),pattern="^coverage-reference-",full.names=TRUE)),
    identical(read_bed("result.bed.gz"),read_bed(file.path(fixture,"1.legacy.bed.gz"))))
  for(row in 1:4) stopifnot(identical(read_counts(paste0(row,".reftnc_plus_strand.tsv")),
    read_counts(file.path(fixture,paste0(row,".legacy.counts.tsv")))))
  setwd(oldwd)
  # Cleanup also applies when a later row fails.
  before <- list.files(getwd(),pattern="^coverage-reference-",full.names=TRUE)
  bad <- file.path(directory,"invalid.bed");writeLines("unknown\t0\t1\t2",bad)
  expect_error(annotate_coverage_row(bad,compressed,fai,1L,paste0(bad,".counts")))
  stopifnot(file.exists(bad),identical(before,list.files(getwd(),pattern="^coverage-reference-",full.names=TRUE)))
}

# Correct the historical seqkit/awk delimiter bug: contig names are literal
# FAI/BED fields, never split at colon, hyphen or space. Include an uncovered
# empty contig, wrapped lines, ambiguity/lowercase, and a multi-block BGZF.
named <- file.path(output,"literal-contigs");dir.create(named)
sequences <- setNames(c("", "aCGTRYacgtN", "n", "aR", paste(rep("AcGTN",15000L),collapse="")),
                     c("chrEmpty", "chr-1", "chr:2", "chr space", "chr-long:3"))
ref <- file.path(named,"reference.fa")
con <- file(ref,"wb");index_lines <- character();expected <- character();input_lines <- character()
for(chromosome in names(sequences)){
  value <- sequences[[chromosome]];len <- nchar(value)
  writeBin(charToRaw(paste0(">",chromosome,"\n")),con)
  offset <- seek(con)
  if(len){
    chunks <- substring(value,seq.int(1L,len,17L),pmin(seq.int(1L,len,17L)+16L,len))
    writeBin(charToRaw(paste0(paste(chunks,collapse="\n"),"\n")),con)
    normalized <- gsub("[^ACGT]","N",toupper(value))
    contexts <- rep(".",len)
    if(len > 2L) contexts[2:(len-1L)] <- substring(normalized,1:(len-2L),3:len)
    expected <- c(expected,paste(chromosome,0:(len-1L),1:len,"02",contexts,sep="\t"))
    input_lines <- c(input_lines,paste(chromosome,0,len,"02",sep="\t"))
  }
  index_lines <- c(index_lines,paste(chromosome,len,offset,min(len,17L),if(len) min(len,17L)+1L else 0L,sep="\t"))
}
close(con);writeLines(index_lines,paste0(ref,".fai"))
for(format in c("plain","gzip","bgzf")){
  reference <- if(format == "plain") ref else paste0(ref,".",format)
  if(format == "gzip"){
    con <- gzfile(reference,"wb");writeBin(readBin(ref,"raw",file.info(ref)$size),con);close(con)
  }else if(format == "bgzf") stopifnot(system2("bgzip",c("-c",shQuote(ref)),stdout=reference) == 0L)
  bed <- file.path(named,paste0(format,".bed"));writeLines(input_lines,bed)
  annotated <- paste0(bed,".gz")
  annotate_coverage_row(bed,reference,paste0(ref,".fai"),1L,paste0(bed,".counts"),
    annotated,"bgzip","tabix",window_bases=30001L,input_rows=2L)
  stopifnot(identical(read_bed(annotated),expected))
  con <- textConnection(expected);expected_table <- read.delim(con,header=FALSE,colClasses="character");close(con)
  expected_counts <- aggregate(as.numeric(expected_table$V4),list(context=expected_table$V5),sum)
  actual <- read_counts(paste0(bed,".counts"))
  stopifnot(identical(actual$context,expected_counts$context),identical(actual$count,expected_counts$x))
  indexed_contigs <- system2("tabix",c("-l",shQuote(annotated)),stdout=TRUE)
  stopifnot(identical(indexed_contigs,names(sequences)[-1L]))
}
# Only zero-length records and empty coverage are valid, including gzip.
empty_ref <- file.path(named,"empty.fa");writeLines(">empty",empty_ref)
empty_fai <- paste0(empty_ref,".fai");writeLines("empty\t0\t7\t0\t0",empty_fai)
empty_bed <- file.path(named,"empty.bed");file.create(empty_bed)
annotate_coverage_row(empty_bed,empty_ref,empty_fai,1L,paste0(empty_bed,".counts"),
  paste0(empty_bed,".gz"),"bgzip","tabix")
stopifnot(identical(readLines(paste0(empty_bed,".counts")),"1\tNA\t0"),
  identical(read_bed(paste0(empty_bed,".gz")),character()))

# Reject invalid indices and input intervals explicitly, without an alternate
# annotation path. Fractional/nonfinite offsets must not reach file seeking.
for(line in c("chr1\t8\t0.5\t8\t9","chr1\tInf\t0\t8\t9",
              "chr1\t8\t0\t0\t9","chr1\t8\t0\t8\t7",
              "chr1\t8\t0\t8\tNA","chr1\t-1\t0\t8\t9")){
  invalid <- file.path(named,"invalid.fai");writeLines(line,invalid)
  stopifnot(grepl("Invalid FASTA index",expect_error(coverage_annotation_index(ref,invalid)),fixed=TRUE))
}
for(line in c("chrEmpty\t0\t1\t1", "chr-1\t0\t99\t2", "chr-1\t-1\t1\t2",
              "chr-1\t0\t1\tInf", "chr-1\t0\t1\t0", "chr-1\t0\t1\t2\textra",
              "chr:2\t0\t1\t2\nchr-1\t0\t1\t2")){
  invalid <- file.path(named,"invalid.bed");writeLines(line,invalid)
  expect_error(annotate_coverage_row(invalid,ref,paste0(ref,".fai"),1L,paste0(invalid,".counts")))
  stopifnot(file.exists(invalid),!file.exists(paste0(invalid,".counts")))
}
writeLines("PASS single R coverage annotation: legacy oracle, production call, compression, literal contigs, empty contigs and failures",file.path(output,"PASS.txt"))
cat(readLines(file.path(output,"PASS.txt")),"\n")
