#!/usr/bin/env Rscript
# Run on a compute node. Save phase writes one optional disposable QS2 file.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: profile_qs2.R INPUT.qs2 NEW_OUTPUT_PREFIX [ROUNDTRIP.qs2]")
}
suppressPackageStartupMessages(library(qs2))
input <- args[[1L]]
prefix <- args[[2L]]
outputs <- paste0(prefix, c(".phases.tsv", ".components.tsv", ".session.txt"))
if (any(file.exists(outputs)) || (length(args) == 3L && file.exists(args[[3L]]))) {
  stop("Refusing to replace an existing report or roundtrip file")
}
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
rss <- function(field) {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
  if (!length(line)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line))
}
reset_hwm <- function() {
  # Linux supports resetting VmHWM to current RSS with clear_refs=5.
  tryCatch({
    con <- file("/proc/self/clear_refs", open = "w")
    on.exit(close(con))
    writeLines("5", con)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
}
phases <- list()
measure <- function(label, expr) {
  invisible(gc())
  reset <- reset_hwm()
  before_rss <- rss("VmRSS")
  before <- proc.time()
  value <- force(expr)
  timing <- proc.time() - before
  phases[[length(phases) + 1L]] <<- data.frame(
    phase = label, user_seconds = unname(timing[[1]] + timing[[4]]),
    system_seconds = unname(timing[[2]] + timing[[5]]),
    actual_cpu_seconds = unname(sum(timing[c(1, 2, 4, 5)])),
    wall_seconds = unname(timing[[3]]), rss_before_kib = before_rss,
    rss_after_kib = rss("VmRSS"), peak_rss_kib = rss("VmHWM"),
    peak_is_phase_specific = reset)
  write.table(do.call(rbind, phases), outputs[[1]], sep = "\t", row.names = FALSE, quote = FALSE)
  value
}
object <- measure("load", qs2::qs_read(input))
# object.size describes logical in-memory size, not peak process RSS, and can
# double count shared objects. It does not report compressed component sizes.
component <- function(value, path) data.frame(
  path = path, class = paste(class(value), collapse = ";"),
  length = length(value), rows = if (is.null(dim(value))) NA_integer_ else dim(value)[[1]],
  columns = if (length(dim(value)) < 2L) NA_integer_ else dim(value)[[2]],
  logical_bytes = as.numeric(object.size(value)))
sizes <- measure("component_sizes", {
  answer <- list(component(object, "/"))
  if (is.list(object)) for (i in seq_along(object)) {
    nm <- if (is.null(names(object))) as.character(i) else names(object)[[i]]
    answer[[length(answer) + 1L]] <- component(object[[i]], paste0("/", nm))
    if (is.data.frame(object[[i]])) for (column in names(object[[i]])) {
      answer[[length(answer) + 1L]] <- component(object[[i]][[column]], paste0("/", nm, "/", column))
    }
  }
  do.call(rbind, answer)
})
write.table(sizes, outputs[[2]], sep = "\t", row.names = FALSE, quote = TRUE)
if (length(args) == 3L) invisible(measure("save", qs2::qs_save(object, args[[3L]])))
capture.output(list(input = normalizePath(input), input_bytes = file.info(input)$size,
                    qs2 = as.character(packageVersion("qs2")), session = sessionInfo()),
               file = outputs[[3]])
