#!/usr/bin/env Rscript
# Content identities and immutable prepared-artifact publication. The reusable
# implementation lives in sharedFunctions.R; the CLI needs only jsonlite/digest/openssl
# plus GNU stat/sync/cp, util-linux flock and Linux /proc/self/fd.
arguments <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", arguments[startsWith(arguments, "--file=")][[1L]])
source(file.path(dirname(normalizePath(script, mustWork = TRUE)), "sharedFunctions.R"))

artifact_cache_main <- function(arguments) {
  usage <- "Usage: artifactCache.R identify --spec FILE --output FILE | {verify,publish,run} --identity FILE --root DIR [--source DIR | --product PATH ... -- COMMAND ...]"
  if (!length(arguments) || arguments[[1L]] %in% c("-h", "--help")) {
    cat(usage, "\n")
    return(invisible(NULL))
  }
  operation <- arguments[[1L]]
  allowed <- switch(operation, identify = c("spec", "output"), verify = c("identity", "root"),
                    publish = c("identity", "root", "source"), run = c("identity", "root", "product"),
                    stop(usage, call. = FALSE))
  options <- list()
  command <- character()
  i <- 2L
  while (i <= length(arguments)) {
    option <- arguments[[i]]
    if (option == "--") {
      if (i < length(arguments)) command <- arguments[seq.int(i + 1L, length(arguments))]
      break
    }
    if (!startsWith(option, "--") && operation == "run") {
      command <- arguments[seq.int(i, length(arguments))]
      break
    }
    name <- substring(option, 3L)
    if (!name %in% allowed || i == length(arguments)) stop(usage, call. = FALSE)
    if (name == "product") options[[name]] <- c(options[[name]], arguments[[i + 1L]]) else options[[name]] <- arguments[[i + 1L]]
    i <- i + 2L
  }
  if (length(setdiff(allowed, names(options)))) stop(usage, call. = FALSE)
  if (operation == "identify") {
    identity <- artifact_cache_identify(artifact_cache_read_json(options$spec))
    connection <- file(options$output, "wb")
    tryCatch(writeBin(charToRaw(paste0(artifact_cache_canonical(identity), "\n")), connection),
             finally = close(connection))
    cat(artifact_cache_key(identity), "\n", sep = "")
  } else {
    identity <- artifact_cache_read_json(options$identity)
    destination <- switch(operation,
      verify = artifact_cache_verify(options$root, identity),
      publish = artifact_cache_publish(options$root, identity, options$source),
      run = {
        if (!length(command)) stop("run requires a build command after --", call. = FALSE)
        artifact_cache_run(options$root, identity, options$product, command)
      })
    cat(destination, "\n", sep = "")
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  tryCatch(artifact_cache_main(commandArgs(trailingOnly = TRUE)), error = function(error) {
    message(conditionMessage(error))
    quit(status = 1L)
  })
}
