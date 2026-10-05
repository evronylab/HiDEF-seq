#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
repo <- if(length(args)) args[[1]] else "."
source(file.path(repo, "scripts/benchmark/germline_vcf_block.R"))
annotation <- germline_vcf_block(file.path(repo, "bin/filterCalls.R"))$annotation

# Exact original annotation expression, independently retained as the oracle.
legacy <- function(calls, germline_vcf_variants) {
  calls %>% left_join(
    germline_vcf_variants %>%
      group_by(seqnames,start,end,ref_plus_strand,alt_plus_strand) %>%
      summarize(germline_vcf_types_detected = str_c(germline_vcf_type,collapse=","),
                germline_vcf_files_detected = str_c(germline_vcf_file,collapse=","),
                .groups="drop") %>% mutate(germline_vcf.passfilter = FALSE),
    by = join_by(seqnames,start,end,ref_plus_strand,alt_plus_strand)
  ) %>% mutate(germline_vcf.passfilter = replace_na(germline_vcf.passfilter, TRUE),
               germline_vcf.passfilter = if_else(call_class == "MDB", TRUE, germline_vcf.passfilter))
}
variants <- tibble(
  seqnames = factor(c("chr1","chr2","chr1","chr1",NA,"chr1"), levels=c("chr2","chr1","unused")),
  start = c(10L,20L,10L,10L,NA_integer_,30L),
  end = c(10L,20L,10L,10L,NA_integer_,31L),
  ref_plus_strand = c("A","G","A","A",NA,"AT"),
  alt_plus_strand = c("C","T","C","C",NA,"A"),
  germline_vcf_type = c("second","other","first","second",NA,"indel"),
  germline_vcf_file = factor(c("z.vcf","other.vcf","a.vcf","z.vcf",NA,"indel.vcf")),
  call_class = factor(c("SBS","SBS","SBS","SBS","SBS","indel")),
  call_type = factor(c("A>C","G>T","A>C","A>C",NA,"deletion")))
calls <- variants[c(1L,1L,5L,6L),] %>%
  select(seqnames,start,end,ref_plus_strand,alt_plus_strand) %>%
  mutate(call_class = factor(c("SBS","MDB","SBS","indel"), levels=c("MDB","indel","SBS","unused")),
         row_id = 1:4)
calls <- bind_rows(calls, calls[2L,] %>% mutate(start=99L,end=99L,row_id=5L))
check <- function(input_calls, input_variants) {
  env <- new.env(parent=globalenv())
  env$calls <- input_calls
  env$germline_vcf_variants <- input_variants
  eval(annotation, env)
  expected <- legacy(input_calls, input_variants)
  stopifnot(identical(env$calls, expected), identical(env$germline_vcf_variants, input_variants))
  env$calls
}
result <- check(calls, variants)
stopifnot(identical(result$germline_vcf_types_detected[1:2], rep("second,first,second",2L)),
          identical(result$germline_vcf_files_detected[1:2], rep("z.vcf,a.vcf,z.vcf",2L)),
          identical(result$germline_vcf.passfilter, c(FALSE,TRUE,FALSE,FALSE,TRUE)))
check(calls[0,], variants)
check(calls, variants[0,])
check(calls[0,], variants[0,])
check(calls %>% filter(call_class == "MDB"), variants)
check(calls %>% filter(start == 99L), variants)
check(calls %>% mutate(seqnames=as.character(seqnames)), variants)
check(calls, variants %>% mutate(seqnames=as.character(seqnames)))
cat("PASS: germline annotation preserves exact types, factor levels, NA matching, empty tables, duplicates, source order, MDB annotations and full variant table\n")
