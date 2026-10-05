# Locate the unchanged section boundaries in the actual pipeline AST. Sourceable
# by tiny fixtures and the compute-node benchmark; does not run pipeline code.
germline_vcf_block <- function(path) {
  expressions <- parse(path, keep.source = FALSE)
  marker <- function(text) {
    matches <- which(vapply(expressions, function(node) {
      is.call(node) && identical(node[[1]], as.name("cat")) &&
        length(node) >= 2L && identical(node[[2]], text)
    }, logical(1)))
    if(length(matches) != 1L) stop("Expected one section marker: ", text)
    matches[[1]]
  }
  first <- marker("## Applying germline VCF variant filters...")
  next_section <- marker("## Applying post-germline VCF variant filtering molecule filters...")
  if(next_section <= first) stop("Unexpected germline VCF section order")
  block <- expressions[seq.int(first + 1L, next_section - 1L)]
  annotation <- which(vapply(block, function(node) {
    is.call(node) && identical(node[[1]], as.name("<-")) &&
      identical(node[[2]], as.name("calls"))
  }, logical(1)))
  if(length(annotation) != 1L) stop("Expected one germline annotation assignment")
  list(prefix = expressions[seq_len(first - 1L)], block = block,
       annotation = block[[annotation]])
}
