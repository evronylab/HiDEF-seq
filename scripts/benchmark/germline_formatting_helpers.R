# Shared by the same-process and fresh-process measurement drivers.
#Read assignments from the actual old/new scripts, without running the pipeline.
#Recurse only into control-flow blocks; do not interpret unrelated expressions.
find_assignments <- function(node, name) {
  if(is.call(node) && identical(node[[1]], as.name("<-"))) {
    return(if(identical(node[[2]], as.name(name))) list(node) else list())
  }
  if(is.expression(node) || (is.call(node) && is.symbol(node[[1]]) && as.character(node[[1]]) %in% c("{", "for", "if"))) {
    return(unlist(lapply(as.list(node), find_assignments, name = name), recursive = FALSE))
  }
  list()
}
assignment <- function(expressions, name) {
  hits <- find_assignments(expressions, name)
  if(length(hits) != 1L) stop(paste("Expected one assignment for", name))
  hits[[1]]
}
format_all <- function(code_by_name) {
  env <- new.env(parent = globalenv())
  for(name in initial_names) eval(code_by_name[[name]], env)
  config_code <- code_by_name[config_names]
  format_code <- code_by_name[["germlineVariantCalls.out"]]
  formatted <- list()
  for(chromgroup in env$chromgroups) for(filtergroup in env$filtergroups) {
    env$i <- chromgroup
    env$j <- filtergroup
    for(code in config_code) eval(code, env)
    eval(format_code, env)
    formatted[[length(formatted) + 1L]] <- env$germlineVariantCalls.out
    rm("germlineVariantCalls.out", envir = env)
  }
  bind_rows(formatted)
}
rss <- function(field) {
  if(!file.exists("/proc/self/status")) return(NA_real_)
  line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
  if(!length(line)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line))
}
reset_hwm <- function() {
  tryCatch({
    con <- file("/proc/self/clear_refs", open = "w")
    on.exit(close(con))
    writeLines("5", con)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
}

initial_names <- c("chromgroups", "filtergroups", "strand_identical_cols_keep",
                   "strand_identical_cols_discard", "strand_redundant_cols_discard")
config_names <- c("region_read_filters_cols_keep", "region_genome_filters_cols_keep", "germline_filter_cols_keep")
required_names <- c(initial_names, config_names, "germlineVariantCalls.out")
