#!/usr/bin/env Rscript
# Experimental complete-operation replay. Production code is never sourced.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2L) stop("Use preflight SOURCE.R; prepare INPUT.qs2 PACKET.qs2 CHROMGROUP FILTERGROUP CHROMS; or run MODE SOURCE.R PACKET.qs2 HELPER OUTPUT_DIR")
mode <- args[[1]]

resolve_code <- function(path) {
  expressions <- as.list(parse(path))
  text <- vapply(expressions, function(node) paste(deparse(node, width.cutoff=500L), collapse="\n"), character(1))
  assignment <- function(name) {
    matches <- Filter(function(node) is.call(node) && identical(node[[1]], as.name("<-")) &&
                        identical(node[[2]], as.name(name)), expressions)
    if(length(matches) != 1L) stop(paste("Expected one source assignment:", name))
    matches[[1]]
  }
  containing <- function(token) {
    hits <- which(grepl(token, text, fixed=TRUE))
    if(length(hits) != 1L) stop(paste("Expected one source expression containing:", token))
    expressions[[hits]]
  }
  walks <- which(vapply(expressions, function(node) {
    is.call(node) && identical(node[[1]], as.name("%>%")) &&
      identical(node[[2]], as.name("coverage_annotation_rows")) &&
      is.call(node[[3]]) && identical(node[[3]][[1]], as.name("pwalk"))
  }, logical(1)))
  if(length(walks) != 2L) stop("Expected the original writer and annotation pwalk expressions")
  chunks <- assignment("chunk_runs")
  if(!identical(eval(chunks[[3]]), 1e7)) stop("The original chunk_runs must remain exactly 1e7")
  list(chunk_runs=chunks, writer=expressions[[walks[[1]]]],
       output_name=assignment("get_coverage_reftnc_output_file"),
       union=containing("sort -m -k1,1n -k3,3n -k4,4n"),
       reference_path=assignment("genome_trinuc_file"),
       reference_intersection=containing("makewindows -w 1 -b all.bed"),
       annotation=expressions[[walks[[2]]]])
}

if(mode == "preflight") {
  code <- resolve_code(args[[2]])
  cat("Resolved original operation AST:", paste(names(code), collapse=", "), "\n")
  quit(save="no")
}
suppressPackageStartupMessages({
  library(qs2)
  library(tidyverse)
  library(GenomicRanges)
  library(data.table)
})
options(warn=2)

if(mode == "prepare") {
  if(length(args) != 6L) stop("prepare INPUT.qs2 PACKET.qs2 CHROMGROUP FILTERGROUP CHROMS(comma separated or all)")
  if(file.exists(args[[3]])) stop("Refusing to overwrite prepared benchmark packet")
  object <- qs_read(args[[2]])
  component <- "bam.gr.filtertrack.bytype.coverage_tnc"
  if(is.null(object[[component]])) stop(paste("Final QS2 lacks", component))
  available <- object[[component]] %>%
    select(any_of(c("analysis_id", "individual_id", "sample_id", "chromgroup", "filtergroup",
                    "bc_orientation", "call_class", "call_type", "SBSindel_call_type")))
  write_tsv(available, paste0(args[[3]], ".available-rows.tsv"))
  selected <- object[[component]] %>%
    filter(as.character(chromgroup) == args[[4]], as.character(filtergroup) == args[[5]]) %>%
    select(all_of(c(names(available), "bam.gr.filtertrack.coverage")))
  if(nrow(selected) == 0L) stop("No coverage rows matched requested chromgroup/filtergroup")
  # Final artifacts do not retain orientation-specific strand RleLists. Never
  # substitute combined orientation coverage for missing strand-specific data.
  if(any(selected$bc_orientation != "all_bc_orientations"))
    stop("Selected final QS2 contains orientation rows without original strand RleLists; use a richer intermediate for that benchmark")
  configuration <- object$yaml.config
  rm(object)
  invisible(gc())
  chromosomes <- if(args[[6]] == "all") NULL else strsplit(args[[6]], ",", fixed=TRUE)[[1]]
  if(!is.null(chromosomes)) {
    selected$bam.gr.filtertrack.coverage <- lapply(selected$bam.gr.filtertrack.coverage, function(coverage) {
      if(!all(chromosomes %in% names(coverage))) stop("Requested chromosome absent from a selected coverage RleList")
      coverage[names(coverage) %in% chromosomes]
    })
  }
  selected <- selected %>% mutate(annotation_row_id=row_number())
  metadata <- list(source=normalizePath(args[[2]]), component=component,
                   chromgroup=args[[4]], filtergroup=args[[5]], chromosomes=chromosomes,
                   original_writer_chunk_runs=1e7, original_rows_retained=nrow(selected))
  qs_save(list(yaml.config=configuration, coverage_annotation_rows=selected, metadata=metadata), args[[3]])
  write_tsv(selected %>% select(-bam.gr.filtertrack.coverage), paste0(args[[3]], ".selected-rows.tsv"))
  cat("Prepared", nrow(selected), "real coverage rows; chromosomes:", args[[6]], "\n")
  quit(save="no")
}

if(mode != "run" || length(args) != 6L) stop("run MODE SOURCE.R PACKET.qs2 HELPER OUTPUT_DIR")
arm <- args[[2]]
if(!arm %in% c("legacy", "candidate")) stop("Mode must be legacy or candidate")
code <- resolve_code(args[[3]])
packet <- qs_read(args[[4]])
helper <- normalizePath(args[[5]], mustWork=TRUE)
output <- normalizePath(args[[6]], mustWork=TRUE)
setwd(output)
yaml.config <- packet$yaml.config
coverage_annotation_rows <- packet$coverage_annotation_rows
individual_id <- unique(as.character(coverage_annotation_rows$individual_id))
sample_id_toanalyze <- unique(as.character(coverage_annotation_rows$sample_id))
chromgroup_toanalyze <- unique(as.character(coverage_annotation_rows$chromgroup))
filtergroup_toanalyze <- unique(as.character(coverage_annotation_rows$filtergroup))
stopifnot(length(individual_id)==1L, length(sample_id_toanalyze)==1L,
          length(chromgroup_toanalyze)==1L, length(filtergroup_toanalyze)==1L)
rm(packet)
invisible(gc())
# Supports either an original cache configuration or a resolved effective YAML.
cache_file <- function(path, yaml.config) {
  if(is.null(yaml.config$cache_artifacts)) return(path)
  as.character(yaml.config$cache_artifacts[[basename(path)]])
}
eval(code$output_name)
eval(code$reference_path)
stopifnot(file.exists(genome_trinuc_file), file.exists(yaml.config$genome_fasta), file.exists(yaml.config$genome_fai))
output_paths <- coverage_annotation_rows %>%
  transmute(annotation_row_id, file=pmap_chr(list(call_class, call_type, SBSindel_call_type), get_coverage_reftnc_output_file))
stopifnot(!anyDuplicated(output_paths$file))
write_tsv(output_paths, "output-paths.tsv")
write_tsv(tibble(tabix=as.character(yaml.config$tabix_bin)), "tool-paths.tsv")

# Time only the requested operation. The outer Python wrapper separately records
# startup/packet-load-inclusive CPU/RSS. proc.time includes waited-for children.
before <- proc.time()
eval(code$chunk_runs)
eval(code$writer)
writer_done <- proc.time()
if(arm == "legacy") {
  eval(code$union)
  eval(code$reference_intersection)
  invisible(file.remove("all.bed"))
  eval(code$annotation)
  invisible(file.remove("all.trinuc.bed"))
} else {
  for(index in seq_len(nrow(coverage_annotation_rows))) {
    row <- coverage_annotation_rows[index, ]
    input <- paste0(row$annotation_row_id, ".bed")
    counts <- paste0(row$annotation_row_id, ".reftnc_plus_strand.tsv")
    command <- paste(shQuote(helper), "--bed", shQuote(input),
      "--fasta", shQuote(yaml.config$genome_fasta), "--fai", shQuote(yaml.config$genome_fai),
      "--row-id", row$annotation_row_id, "--counts", shQuote(counts))
    if(row$bc_orientation == "all_bc_orientations") {
      destination <- output_paths$file[[index]]
      command <- paste(command, "--bed-output - |", shQuote(yaml.config$bgzip_bin), "-c >", shQuote(destination),
        "&&", shQuote(yaml.config$tabix_bin), "-@2 -s1 -b2 -e3", shQuote(destination))
    }
    status <- system2("/bin/bash", args="-s", input=c("set -euo pipefail", command))
    if(status != 0L) stop("Candidate annotation pipeline failed")
    invisible(file.remove(input))
  }
}
finished <- proc.time()
metric <- function(start, end, phase) {
  delta <- end-start
  tibble(phase=phase, user_seconds=unname(delta[["user.self"]]+delta[["user.child"]]),
    system_seconds=unname(delta[["sys.self"]]+delta[["sys.child"]]),
    actual_cpu_seconds=user_seconds+system_seconds, wall_seconds=unname(delta[["elapsed"]]))
}
write_tsv(bind_rows(metric(before, writer_done, "original_bed_writer"),
                   metric(writer_done, finished, "annotation_compression_index"),
                   metric(before, finished, "complete_operation")), "operation-metrics.tsv")
cat("Complete", arm, "operation finished for", nrow(coverage_annotation_rows), "rows\n")
