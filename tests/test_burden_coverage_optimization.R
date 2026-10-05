#!/usr/bin/env Rscript
#Run inside the pipeline R environment on a compute node:
#Rscript tests/test_burden_coverage_optimization.R [repository root] [installed BSgenome package directory]
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(plyranges)
  library(Biostrings)
  library(tidyverse)
})
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
repo <- if(length(args)) args[[1]] else "."
expressions <- parse(file.path(repo, "bin", "calculateBurdens.R"))
assignment <- function(expressions, name, occurrence = 1L) {
  matches <- vapply(expressions, function(x) {
    is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(name))
  }, logical(1))
  stopifnot(sum(matches) >= occurrence)
  expressions[[which(matches)[[occurrence]]]]
}
for(name in c("sum_RleList", "bc_orientation_is_asymmetric", "validate_bam.gr.filtertrack",
              "calc_duplex_coverage", "accumulate_bam.gr.filtertracks", "filtertrack_coverage_result",
              "make_per_bc_orientation_coverage")) {
  eval(assignment(expressions, name))
}

#Adapter keeps the existing numerical/consumer checks concise; production keeps
#the environment alive across chunks and materializes a table only at the end.
sum_bam.gr.filtertracks <- function(first, second = NULL, orientation_call_types = NULL) {
  state <- new.env(parent = emptyenv())
  if(is.null(second)) {
    accumulate_bam.gr.filtertracks(state, first, orientation_call_types)
  } else {
    state$metadata <- select(first, -bam.gr.filtertrack.coverage, -bam.gr.filtertrack.by_bc_orientation_strand.coverage)
    state$coverage <- first$bam.gr.filtertrack.coverage
    state$orientation <- first$bam.gr.filtertrack.by_bc_orientation_strand.coverage
    accumulate_bam.gr.filtertracks(state, second, orientation_call_types)
  }
  filtertrack_coverage_result(state)
}

keys <- tribble(
  ~call_type, ~call_class, ~SBSindel_call_type, ~filtergroup,
  "SBS", "SBS", "mutation", "test",
  "SBS", "SBS", "mismatch-ss", "test",
  "indel", "indel", "mutation", "test",
  "indel", "indel", "mismatch-ss", "test",
  "MDB2", "MDB", "mismatch-ss", "test",
  "SBS", "SBS", "mismatch-os", "other"
)
si <- Seqinfo(c("chr1", "chrM"), c(100L, 60L), isCircular = c(FALSE, TRUE))

#The paired orientation is exactly the annotation emitted by extractCalls for
#the final round. First-round orientation must not influence a round-2 sample.
make_molecule <- function(round1, round2 = NULL, id = 1L, reverse = FALSE) {
  barcode <- if(is.null(round2)) round1 else round2
  pair <- strsplit(barcode, "-", fixed = TRUE)[[1]]
  orientation <- c(paste(pair, collapse = "-"), paste(rev(pair), collapse = "-"))
  if(reverse) orientation <- rev(orientation)
  gr <- GRanges(rep(if(id %% 2L) "chr1" else "chrM", 2),
                IRanges(rep(3L + id, 2), width = 8L), strand = c("+", "-"), seqinfo = si)
  gr$run_id <- rep(paste0("run", id %% 2L), 2)
  gr$zm <- rep(id, 2)
  gr$bc_orientation <- factor(orientation)
  gr
}
tracks <- function(gr) mutate(keys, bam.gr.filtertrack = rep(list(gr), nrow(keys)))

#Independent baseline: compute every observed orientation, as before optimization.
legacy_tracks <- function(gr) {
  orientations <- as_tibble(gr) %>% distinct(bc_orientation, strand)
  orientations$bam.gr.filtertrack.coverage <- map2(
    orientations$bc_orientation, orientations$strand,
    ~ coverage(gr[gr$bc_orientation == .x & as.character(strand(gr)) == as.character(.y)]))
  keys %>% mutate(bam.gr.filtertrack.coverage = rep(list(coverage(gr) %/% 2L), n()),
                  bam.gr.filtertrack.by_bc_orientation_strand.coverage = rep(list(orientations), n()))
}
consumer_rows <- function(cov, requested) {
  env <- new.env(parent = globalenv())
  env$bam.gr.filtertrack.bytype <- cov
  env$call_types_toanalyze <- requested
  eval(assignment(expressions, "coverage_rows"), env)
  eval(assignment(expressions, "coverage_annotation_rows"), env)
  #The aggregate rows' orientation list is never consumed or serialized; remove
  #that internal column, comparing every actual aggregate/template annotation.
  list(rows = select(env$coverage_rows, -bam.gr.filtertrack.by_bc_orientation_strand.coverage),
       annotation = select(env$coverage_annotation_rows, -bam.gr.filtertrack.by_bc_orientation_strand.coverage))
}

scenarios <- list(
  round1_symmetric = list("A-A", NULL),
  round1_asymmetric = list("A-B", NULL),
  symmetric_then_asymmetric = list("A-A", "C-D"),
  asymmetric_then_symmetric = list("A-B", "C-C"),
  both_asymmetric = list("A-B", "C-D"),
  both_symmetric = list("A-A", "C-C")
)
requested_sets <- list(
  all = filter(keys, filtergroup == "test"),
  mutation_only = filter(keys, SBSindel_call_type == "mutation"),
  implicit_sbs_mismatch = filter(keys, call_type != "SBS", filtergroup == "test")
)
for(scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  chunks <- list(make_molecule(scenario[[1]], scenario[[2]], 1L),
                 make_molecule(scenario[[1]], scenario[[2]], 2L, reverse = TRUE))
  #Empty chunks and symmetric/asymmetric runs in the same sample exercise typed
  #empty joins, newly encountered orientations and preserved row ordering.
  variants <- list(pure = chunks,
                   mixed = c(list(chunks[[1]][0]), chunks,
                             list(make_molecule("X-X", NULL, 3L), make_molecule("E-F", NULL, 4L))))
  for(variant in variants) {
    combined <- do.call(c, variant)
    baseline <- legacy_tracks(combined)
    for(requested in requested_sets) {
      optimized <- sum_bam.gr.filtertracks(tracks(variant[[1]]), orientation_call_types = requested)
      for(gr in variant[-1]) {
        optimized <- sum_bam.gr.filtertracks(optimized, tracks(gr), orientation_call_types = requested)
      }
      stopifnot(identical(optimized$bam.gr.filtertrack.coverage, baseline$bam.gr.filtertrack.coverage))
      stopifnot(identical(consumer_rows(optimized, requested), consumer_rows(baseline, requested)))
      eligible <- keys %>% mutate(row = row_number()) %>%
        filter(SBSindel_call_type != "mutation") %>% semi_join(requested, by = names(keys)) %>% pull(row)
      for(i in seq_len(nrow(keys))) {
        observed <- optimized$bam.gr.filtertrack.by_bc_orientation_strand.coverage[[i]]
        if(!i %in% eligible) stopifnot(nrow(observed) == 0L)
        stopifnot(all(bc_orientation_is_asymmetric(observed$bc_orientation)))
      }
    }
    sensitivity <- sum_bam.gr.filtertracks(tracks(variant[[1]]))
    for(gr in variant[-1]) sensitivity <- sum_bam.gr.filtertracks(sensitivity, tracks(gr))
    stopifnot(identical(sensitivity$bam.gr.filtertrack.coverage, baseline$bam.gr.filtertrack.coverage),
              all(map_int(sensitivity$bam.gr.filtertrack.by_bc_orientation_strand.coverage, nrow) == 0L))
  }
  cat("PASS:", scenario_name, "(pure/mixed, empty chunks, requested categories, sensitivity)\n")
}


#Compare the real persistent accumulator with the previous whole-table helper.
#No table alias is materialized from the candidate until all chunks finish.
legacy_env <- new.env(parent = globalenv())
sys.source(file.path(repo, "tests", "fixtures", "sum_filtertrack_coverage_before_accumulator.R"), envir = legacy_env)
chunks <- list(tracks(make_molecule("A-A", id = 1L)),
               tracks(make_molecule("C-D", id = 2L, reverse = TRUE)),
               tracks(make_molecule("A-A", id = 3L)[0]),
               tracks(make_molecule("E-F", id = 4L)))
#Reordered rows and factor-level growth must have the former join semantics.
chunks[[1]]$filtergroup <- factor(chunks[[1]]$filtergroup, levels = c("test", "other"))
for(i in 2:length(chunks)) {
  chunks[[i]] <- chunks[[i]][rev(seq_len(nrow(keys))),]
  chunks[[i]]$filtergroup <- factor(chunks[[i]]$filtergroup, levels = c("other", "test", "unused"))
}
state <- new.env(parent = emptyenv())
old <- NULL
for(chunk in chunks) {
  accumulate_bam.gr.filtertracks(state, chunk, requested_sets$all)
  old <- if(is.null(old)) legacy_env$sum_bam.gr.filtertracks(chunk, orientation_call_types = requested_sets$all) else
    legacy_env$sum_bam.gr.filtertracks(old, chunk, orientation_call_types = requested_sets$all)
}
stopifnot(identical(filtertrack_coverage_result(state), old))
#A one-to-many metadata match preserves legacy row expansion, too.
state <- new.env(parent = emptyenv())
one <- chunks[[1]][1,,drop = FALSE]
two <- chunks[[2]][c(nrow(chunks[[2]]), nrow(chunks[[2]])),,drop = FALSE]
accumulate_bam.gr.filtertracks(state, one)
accumulate_bam.gr.filtertracks(state, two)
old <- legacy_env$sum_bam.gr.filtertracks(legacy_env$sum_bam.gr.filtertracks(one), two)
stopifnot(identical(filtertrack_coverage_result(state), old))
cat("PASS: persistent category accumulator, metadata order/factors, orientation growth and legacy row expansion\n")

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE))
}
bad <- make_molecule("A-A")
start(bad)[2] <- start(bad)[2] + 1L
expect_error(sum_bam.gr.filtertracks(tracks(bad)), "Mismatched plus and minus strand ranges")
#Matching ranges/metadata but different seqnames reaches the independent even
#coverage assertion. This must still fail when orientation tracks are disabled.
bad <- make_molecule("A-A")
seqnames(bad) <- c("chr1", "chrM")
expect_error(sum_bam.gr.filtertracks(tracks(bad)), "Non-even strand coverage")
cat("PASS: plus/minus and even-coverage invariants\n")

#Compare spectra against an independent full miniature reference. Boundary
#contexts on chrM keep the pre-existing clipping/N-padding behavior (no new
#circular wrapping). The unused chromosome must never be requested.
shared <- parse(file.path(repo, "bin", "sharedFunctions.R"))
eval(assignment(shared, "indel.spectrum"))
reference <- DNAStringSet(c(chr1 = paste(rep("ACGT", 50), collapse = ""),
                           chrM = paste(rep("CATG", 30), collapse = ""),
                           unused = paste(rep("G", 1000), collapse = "")))
indels <- tribble(
  ~CHROM, ~POS, ~REF, ~ALT, ~TEMPLATE_STRAND,
  "chr1", 1L, "AC", "A", "+",
  "chr1", 199L, "GT", "G", "-",
  "chr1", 40L, "T", "TACGTAC", "+",
  "chr1", 80L, "TACGTAC", "T", "-",
  "chrM", 1L, "C", "CA", "-",
  "chrM", 119L, "TG", "T", "+",
  "chrM", 30L, "ATGC", "A", "-",
  "chrM", 70L, "A", "ATGC", "+"
) %>% as.data.frame()
for(calls in list(indels, filter(indels, CHROM == "chrM"), indels[0,])) {
  env <- new.env(parent = globalenv())
  env$toy_reference <- reference
  env$BSgenome_name <- "toy_reference"
  env$requested_chroms <- NULL
  env$getSeq <- function(x, names) {
    env$requested_chroms <- names
    #Match the actual BSgenome method: a singleton character query simplifies
    #to DNAString and loses its chromosome name. The old mock hid this case.
    if(length(names) == 1L) x[[names]] else x[names]
  }
  env$finalCalls.reftnc_spectra <- tibble(
    call_class = c("SBS", "indel", "indel"),
    finalCalls_for_vcf = list(tibble(seqnames = "unused"),
                             tibble(seqnames = factor(calls$CHROM, levels = names(reference))),
                             tibble(seqnames = calls$CHROM)))
  eval(assignment(expressions, "indel_spectrum_chroms"), env)
  eval(assignment(expressions, "BSgenome_for_indel.spectrum"), env)
  stopifnot(is(env$BSgenome_for_indel.spectrum, "DNAStringSet"))
  stopifnot(identical(as.character(names(env$BSgenome_for_indel.spectrum)), unique(calls$CHROM)))
  if(nrow(calls) == 0L) stopifnot(is.null(env$requested_chroms))
  for(spectrum_type in c("pyr", "template")) {
    stopifnot(identical(indel.spectrum(calls, reference, spectrum_type = spectrum_type),
                        indel.spectrum(calls, env$BSgenome_for_indel.spectrum, spectrum_type = spectrum_type)))
  }
}
cat("PASS: full/subset indel spectra, absent indels, nuclear/mitochondrial terminal contexts\n")

#When an installed reference is supplied, verify the pinned BSgenome API itself
#rather than only its faithful miniature mock. Only chrM and chrY are loaded.
if(length(args) >= 2L) {
  suppressPackageStartupMessages(library(BSgenome))
  package_dir <- normalizePath(args[[2]], mustWork = TRUE)
  package_name <- basename(package_dir)
  suppressPackageStartupMessages(library(package_name, character.only = TRUE,
                                         lib.loc = dirname(package_dir)))
  genome <- get(package_name)
  scalar <- getSeq(genome, names = "chrM")
  stopifnot(is(scalar, "DNAString"), !is(scalar, "DNAStringSet"))
  legacy_reference <- getSeq(genome, names = c("chrM", "chrY"))
  stopifnot(is(legacy_reference, "DNAStringSet"),
            identical(names(legacy_reference), c("chrM", "chrY")))
  for(chroms in list("chrM", "chrY", c("chrY", "chrM"), character())) {
    env <- new.env(parent = globalenv())
    env$BSgenome_name <- package_name
    env$indel_spectrum_chroms <- chroms
    eval(assignment(expressions, "BSgenome_for_indel.spectrum"), env)
    expected <- if(length(chroms)) legacy_reference[chroms] else DNAStringSet()
    observed <- env$BSgenome_for_indel.spectrum
    #XStringSet subsets can retain the unused backing sequence pool. Compare
    #the exact ordered bases, names, widths and class, not pool allocation.
    stopifnot(identical(class(observed), class(expected)),
              identical(names(observed), names(expected)),
              identical(width(observed), width(expected)),
              identical(as.character(observed), as.character(expected)))
  }
  cat("PASS: actual BSgenome singleton simplification; named singleton/multiple/empty references exact\n")
}
