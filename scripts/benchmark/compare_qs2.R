#!/usr/bin/env Rscript
# Sourceable comparator plus disk-staged CLI. Scratch component files are
# validation artifacts only; pipeline output remains one final QS2 per sample.
scientific_compare <- function(reference, candidate, path = "", ignore = "/run_metadata",
                               chromgroup_blocks = FALSE, tolerance = 1e-12,
                               max_reports = 1000L) {
  reports <- list()
  counts <- c(failure = 0L, review = 0L, ignored = 0L)
  report <- function(path, status, reason, max_relative = NA_real_) {
    counts[[status]] <<- counts[[status]] + 1L
    if (length(reports) < max_reports) reports[[length(reports) + 1L]] <<-
      data.frame(path = path, status = status, reason = reason, max_relative = max_relative)
  }
  child_path <- function(path, name) paste0(path, "/", gsub("/", "~1", gsub("~", "~0", name, fixed = TRUE), fixed = TRUE))
  sorted_attributes <- function(value) {
    # R attributes form a named set (identical's default attrib.as.set=TRUE),
    # unlike scientific named-list elements, whose order remains significant.
    attrs <- attributes(value)
    if (is.null(attrs)) NULL else attrs[order(names(attrs), method = "radix")]
  }
  canonical <- function(x) {
    if (chromgroup_blocks && is.data.frame(x) && "chromgroup" %in% names(x)) {
      # Stable order among rows in the same chromgroup; never sort scientific
      # calls or other keys within a group. Automatic row names carry no data.
      automatic <- .row_names_info(x, type = 1L) < 0L
      x <- x[order(as.character(x$chromgroup), na.last = TRUE, method = "radix"), , drop = FALSE]
      if (automatic) attr(x, "row.names") <- .set_row_names(nrow(x))
    }
    x
  }
  compare <- function(a, b, where) {
    if (where %in% ignore) {
      report(where, "ignored", "explicit metadata whitelist")
      return(invisible(NULL))
    }
    if (identical(a, b)) return(invisible(NULL))
    if (!identical(typeof(a), typeof(b)) || !identical(class(a), class(b))) {
      report(where, "failure", "storage type or class differs")
      return(invisible(NULL))
    }
    # Exact whitelisted metadata keys may be added to named configuration
    # mappings. Remove only those named entries, never an entire configuration
    # or a scientific data-frame column. All remaining schema/attributes/order
    # and values are still compared normally.
    if (is.list(a) && !is.data.frame(a) && is.null(dim(a)) && is.null(dim(b)) &&
        !is.null(names(a)) && !is.null(names(b)) &&
        !anyDuplicated(names(a)) && !anyDuplicated(names(b))) {
      candidates <- union(names(a), names(b))
      ignored_names <- candidates[vapply(candidates, function(nm) child_path(where, nm) %in% ignore, logical(1))]
      if (length(ignored_names)) {
        for (nm in ignored_names) report(child_path(where, nm), "ignored", "explicit metadata whitelist")
        remove_metadata <- function(value) {
          attrs <- attributes(value)
          value <- value[!names(value) %in% ignored_names]
          attrs$names <- names(value)
          attributes(value) <- attrs
          value
        }
        a <- remove_metadata(a)
        b <- remove_metadata(b)
      }
    }
    a <- canonical(a)
    b <- canonical(b)
    if (identical(a, b)) return(invisible(NULL))
    if (isS4(a)) {
      if (!identical(slotNames(a), slotNames(b))) {
        report(where, "failure", "S4 slot names differ")
        return(invisible(NULL))
      }
      for (slot_name in slotNames(a)) compare(slot(a, slot_name), slot(b, slot_name), paste0(where, "/@", slot_name))
      aa <- sorted_attributes(a); bb <- sorted_attributes(b)
      aa[slotNames(a)] <- NULL; bb[slotNames(b)] <- NULL
      compare(aa, bb, paste0(where, "/@attributes"))
      return(invisible(NULL))
    }
    if (!identical(length(a), length(b))) {
      report(where, "failure", "length differs")
      return(invisible(NULL))
    }
    compare(sorted_attributes(a), sorted_attributes(b), paste0(where, "/@attributes"))
    if (is.list(a) || is.pairlist(a)) {
      for (i in seq_along(a)) {
        nm <- if (is.null(names(a)) || is.na(names(a)[[i]]) || !nzchar(names(a)[[i]])) paste0("[[", i, "]]") else names(a)[[i]]
        compare(a[[i]], b[[i]], child_path(where, nm))
      }
    } else if (is.double(a)) {
      # Numeric differences are review-required even within tolerance. Relative
      # error is symmetric, with no absolute tolerance around zero.
      for (offset in seq.int(1, max(1, length(a)), by = 1000000)) {
        if (!length(a)) break
        idx <- seq.int(offset, min(length(a), offset + 999999))
        av <- a[idx]; bv <- b[idx]
        special_a <- !is.finite(av); special_b <- !is.finite(bv)
        if (!identical(is.na(av), is.na(bv)) || !identical(is.nan(av), is.nan(bv)) ||
            !identical(special_a, special_b) || !identical(av[special_a], bv[special_b])) {
          report(where, "failure", paste("NA, NaN, or infinity differs in chunk starting", offset))
          next
        }
        changed <- which(!special_a & av != bv)
        if (length(changed)) {
          relative <- abs(av[changed] - bv[changed]) / pmax(abs(av[changed]), abs(bv[changed]))
          worst <- max(relative)
          report(where, if (is.finite(worst) && worst <= tolerance) "review" else "failure",
                 paste(length(changed), "numeric values differ; first index", idx[changed[[1]]]), worst)
        }
      }
    } else {
      # Attributes have already been compared independently.
      attributes(a) <- NULL; attributes(b) <- NULL
      if (!identical(a, b)) report(where, "failure", "values or order differ")
    }
    invisible(NULL)
  }
  compare(reference, candidate, path)
  list(counts = counts, reports = if (length(reports)) do.call(rbind, reports) else
    data.frame(path = character(), status = character(), reason = character(), max_relative = double()))
}

# Validate and normalize only complete, ordered tokens in explicitly named
# provenance columns. Scientific values and column/vector attributes stay intact.
normalize_provenance_tokens <- function(x, rule, side) {
  stopifnot(is.character(x), !is.object(x), side %in% c("reference","candidate"))
  old <- unlist(rule$reference,use.names=FALSE)
  current <- unlist(rule$candidate,use.names=FALSE)
  stopifnot(length(old)>0L, length(old)==length(current), !anyDuplicated(old), !anyDuplicated(current))
  allowed <- if(side=="reference") old else current
  values <- unique(x)
  normalized <- vapply(values,function(value) {
    if(is.na(value) || identical(value,"")) return(value)
    tokens <- strsplit(value,",",fixed=TRUE)[[1]]
    stopifnot(identical(paste(tokens,collapse=","),value), all(tokens %in% allowed))
    paste(old[match(tokens,allowed)],collapse=",")
  },character(1),USE.NAMES=FALSE)
  x[] <- normalized[match(x,values)]
  x
}

normalize_provenance_column <- function(x, component, rules, side) {
  for(path in names(rules)) {
    pieces <- strsplit(sub("^/","",path),"/",fixed=TRUE)[[1]]
    stopifnot(length(pieces)==2L, identical(pieces[[2]],"germline_vcf_files_detected"))
    if(identical(pieces[[1]],component)) {
      stopifnot(is.data.frame(x), sum(names(x)==pieces[[2]])==1L)
      x[[pieces[[2]]]] <- normalize_provenance_tokens(x[[pieces[[2]]]],rules[[path]],side)
    }
  }
  x
}

comparison_main <- function(args) {
  if (length(args) < 3L) stop(paste(
    "Usage: compare_qs2.R stage REFERENCE.qs2 NEW_SCRATCH_DIR",
    "or: compare_qs2.R compare SCRATCH_DIR CANDIDATE.qs2 REPORT.tsv",
    "[--chromgroup-blocks] [--ignore=/exact/path] [--token-policy=rules.json]", sep = "\n"))
  suppressPackageStartupMessages(library(qs2))
  if (args[[1]] == "stage") {
    if (dir.exists(args[[3]])) stop("Scratch directory already exists")
    dir.create(args[[3]], recursive = TRUE)
    reference <- qs2::qs_read(args[[2]])
    if (!is.list(reference) || is.data.frame(reference)) stop("Expected top-level list")
    saveRDS(list(attributes = attributes(reference), length = length(reference),
                 type = typeof(reference), class = class(reference)), file.path(args[[3]], "manifest.rds"))
    for (i in seq_along(reference)) {
      qs2::qs_save(reference[[i]], file.path(args[[3]], sprintf("%04d.qs2", i)))
      reference[i] <- list(NULL)
      invisible(gc())
    }
    writeLines("complete", file.path(args[[3]], "COMPLETE"))
    return(0L)
  }
  if (args[[1]] != "compare" || length(args) < 4L) stop("Invalid command")
  if (!file.exists(file.path(args[[2]], "COMPLETE"))) stop("Reference staging incomplete")
  if (file.exists(args[[4]])) stop("Report already exists")
  flags <- if (length(args) > 4L) args[5:length(args)] else character()
  if (any(!flags %in% "--chromgroup-blocks" & !startsWith(flags, "--ignore=") & !startsWith(flags, "--token-policy="))) stop("Unknown option")
  token_flags <- flags[startsWith(flags,"--token-policy=")]
  stopifnot(length(token_flags)<=1L)
  token_rules <- if(length(token_flags)) jsonlite::fromJSON(sub("^--token-policy=","",token_flags),simplifyVector=FALSE)$columns else list()
  ignore <- c("/run_metadata", sub("^--ignore=", "", flags[startsWith(flags, "--ignore=")]))
  manifest <- readRDS(file.path(args[[2]], "manifest.rds"))
  candidate <- qs2::qs_read(args[[3]])
  candidate_manifest <- list(attributes = attributes(candidate), length = length(candidate),
                             type = typeof(candidate), class = class(candidate))
  if (!identical(manifest, candidate_manifest)) stop("Root list schema/attributes differ")
  counts <- c(failure = 0L, review = 0L, ignored = 0L)
  reports <- list()
  for (i in seq_along(candidate)) {
    reference_part <- qs2::qs_read(file.path(args[[2]], sprintf("%04d.qs2", i)))
    nm <- if (is.null(names(candidate))) paste0("[[", i, "]]") else names(candidate)[[i]]
    reference_part <- normalize_provenance_column(reference_part,nm,token_rules,"reference")
    candidate[[i]] <- normalize_provenance_column(candidate[[i]],nm,token_rules,"candidate")
    result <- scientific_compare(reference_part, candidate[[i]], path = paste0("/", nm),
                                 ignore = ignore, chromgroup_blocks = "--chromgroup-blocks" %in% flags)
    counts <- counts + result$counts
    reports[[i]] <- result$reports
    cat(nm, paste(names(result$counts), result$counts, collapse = "; "), "\n")
    rm(reference_part)
    candidate[i] <- list(NULL)
    invisible(gc())
  }
  write.table(do.call(rbind, reports), args[[4]], sep = "\t", quote = TRUE, row.names = FALSE)
  print(counts)
  if (counts[["failure"]]) 1L else if (counts[["review"]]) 2L else 0L
}

if (sys.nframe() == 0L) quit(status = comparison_main(commandArgs(trailingOnly = TRUE)))
