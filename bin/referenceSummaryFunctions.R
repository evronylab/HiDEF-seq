# Small reusable reference summaries. Callers load BSgenome, GenomicRanges and
# Biostrings. Reference sequence boundaries and literal N matching are unchanged.
reference_n_ranges <- function(genome) {
  reference_info <- seqinfo(genome)
  chromosomes <- seqnames(genome)
  pieces <- lapply(chromosomes, function(chromosome) {
    positions <- IRanges::reduce(ranges(matchPattern("N", genome[[chromosome]])))
    GRanges(seqnames = rep(chromosome, length(positions)), ranges = positions,
            strand = rep("*", length(positions)), seqinfo = reference_info)
  })
  if (!length(pieces)) return(GRanges(seqinfo = reference_info))
  GenomicRanges::reduce(do.call(c, unname(pieces)), ignore.strand = TRUE)
}

reference_trinucleotide_counts <- function(genome) {
  chromosomes <- seqnames(genome)
  setNames(lapply(chromosomes, function(chromosome) {
    trinucleotideFrequency(genome[[chromosome]])
  }), chromosomes)
}

reference_counts_for_chromosomes <- function(counts, chromosomes) {
  chromosomes <- as.character(chromosomes)
  if (!length(chromosomes)) {
    return(trinucleotideFrequency(DNAStringSet(), simplify.as = "collapsed"))
  }
  if (anyNA(chromosomes) || any(!chromosomes %in% names(counts))) {
    stop("Requested chromosome is absent from the prepared reference summary", call. = FALSE)
  }
  # Preserve requested order and repeated chromosomes, as getSeq(chromosomes)
  # did. Integer counts and their channel names retain their original types.
  Reduce(`+`, counts[chromosomes])
}
