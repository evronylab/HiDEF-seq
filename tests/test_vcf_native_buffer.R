# Sourced in the pinned compute preflight, with the actual writer and buffer
# modifier already resolved. Fake seqinfo suffices: this tests export, not getSeq.
original_writer <- eval(writer)
small_writer <- with_native_buffer(original_writer, 1L)
medium_writer <- with_native_buffer(original_writer, 3L)
stopifnot(identical(as.list(body(original_writer))[-length(body(original_writer))],
                    as.list(body(small_writer))[-length(body(small_writer))]))
fixture_genome <- GRanges(seqinfo = Seqinfo("chr1", 1000L))
calls_fixture <- tibble(
  seqnames = factor(rep("chr1", 9L)),
  start_refspace = c(6L, 1L, 1L, 2L, 3L, 4L, 4L, 5L, 7L),
  ref_plus_strand = c("A", "AT", "A", "C", "T", "G", "G", "T", "C"),
  alt_plus_strand = c("C", "A", "AGG", "G", "A", "T", "A", "C", "G"),
  numeric_mixed = c(1, 1e-15, 1e20, 1234567890.123456, NA_real_, 0, -1e-15,
                    .Machine$double.xmin, .Machine$double.xmax),
  numeric_zero = rep(0, 9L),
  integer_values = c(NA_integer_, 1:8),
  logical_values = c(TRUE, FALSE, NA, TRUE, FALSE, TRUE, FALSE, NA, TRUE),
  factor_values = factor(c("x", "y", NA, "x", "z", "y", "x", "y", "z")),
  string_values = c("1,2", "-3,4", NA, "", "test", "5e-10", ".", "x", "y"))
read_scientific <- function(path) {
  con <- gzfile(path)
  on.exit(close(con))
  lines <- readLines(con)
  lines[!startsWith(lines, "##fileDate=")]
}
for(empty in c(FALSE, TRUE)) {
  input_fixture <- if(empty) calls_fixture[FALSE, ] else calls_fixture
  paths <- paste0(prefix, ".fixture.", empty, ".", c("default", "one", "three"), ".vcf")
  original_writer(input_fixture, "fixture_genome", paths[[1]])
  small_writer(input_fixture, "fixture_genome", paths[[2]])
  medium_writer(input_fixture, "fixture_genome", paths[[3]])
  for(path in paths[-1L]) {
    stopifnot(identical(read_scientific(paste0(paths[[1]], ".bgz")),
                        read_scientific(paste0(path, ".bgz"))))
    if(!empty) for(pos in c(1L, 4L, 7L)) {
      region <- GRanges("chr1", IRanges(pos, width = 1L))
      stopifnot(identical(scanTabix(paste0(paths[[1]], ".bgz"), param = region),
                          scanTabix(paste0(path, ".bgz"), param = region)))
    }
  }
}
cat("PASS: actual writer native buffer preserves floats, NA/flags, factors, empty tables, anchored indels, duplicate loci and indexed queries\n")
