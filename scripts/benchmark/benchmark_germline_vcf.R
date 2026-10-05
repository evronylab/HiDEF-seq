#!/usr/bin/env Rscript
# Full VCF load/filter/summary/join benchmark using actual original/current ASTs.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 8L || length(args) > 9L) stop(paste(
  "Usage: benchmark_germline_vcf.R BASELINE_filterCalls.R CANDIDATE_filterCalls.R",
  "PREPARED_CONFIG EXTRACT.qs2 SAMPLE CHROMGROUP FILTERGROUP NEW_OUTPUT_DIR [REPEATS=3]"))
paths <- normalizePath(args[1:4], mustWork=TRUE)
output <- args[[8]]
if(dir.exists(output) || file.exists(output)) stop("Refusing existing output: ", output)
repeats <- if(length(args)==9L) as.integer(args[[9]]) else 3L
stopifnot(length(repeats)==1L, !is.na(repeats), repeats>0L)
dir.create(output, recursive=TRUE)
output <- normalizePath(output, mustWork=TRUE)
script_argument <- grep("^--file=", commandArgs(), value=TRUE)
source(file.path(dirname(sub("^--file=", "", script_argument[[1]])), "germline_vcf_block.R"))
code <- lapply(paths[1:2], germline_vcf_block)
names(code) <- c("baseline","candidate")
writeLines(capture.output(dput(list(arguments=args, source_md5=tools::md5sum(paths)))), file.path(output,"inputs.txt"))
writeLines(vapply(code$baseline$block, function(x) paste(deparse(x),collapse="\n"), character(1)), file.path(output,"baseline-block.R"))
writeLines(vapply(code$candidate$block, function(x) paste(deparse(x),collapse="\n"), character(1)), file.path(output,"candidate-block.R"))

# Execute the current real preceding filters exactly once, with helpers bound to
# the same lexical environment as their pipeline globals. Both arms inherit the
# same resulting calls/configuration. Preparation is outside block timing.
script_bin <- dirname(paths[[2]])
Sys.setenv(PATH=paste(script_bin,Sys.getenv("PATH"),sep=.Platform$path.sep))
setwd(output)
context <- new.env(parent=globalenv())
filter_args <- c("-c",paths[[3]],"-f",paths[[4]],"-s",args[[5]],"-g",args[[6]],
                 "-v",args[[7]],"-o",file.path(output,"unused-full-filter.qs2"))
context$parse_args <- function(object, ...) optparse::parse_args(object,args=filter_args,...)
context$source <- function(file, local=FALSE, ...) base::source(file,local=context,...)
eval(code$candidate$prefix, context)
cat("Prepared actual pre-VCF calls:", nrow(context$calls), "\n")

rss <- function(field) {
  if(!file.exists("/proc/self/status")) return(NA_real_)
  line <- grep(paste0("^",field,":"),readLines("/proc/self/status"),value=TRUE)
  if(!length(line)) return(NA_real_)
  as.numeric(gsub("[^0-9]","",line))
}
reset_hwm <- function() tryCatch({
  con <- file("/proc/self/clear_refs",open="w")
  on.exit(close(con))
  writeLines("5",con)
  TRUE
},error=function(e) FALSE,warning=function(w) FALSE)
metrics <- list()
comparisons <- list()
for(pair in seq_len(repeats)) {
  modes <- if(pair %% 2L) names(code) else rev(names(code))
  for(mode in modes) {
    env <- new.env(parent=context)
    env$calls <- context$calls
    # Original AST predates immutable-cache resolution. Resolve its unchanged
    # requested basename to the same prepared file used by the candidate.
    env$qs_read <- function(file, ...) qs2::qs_read(context$cache_file(file,context$yaml.config),...)
    invisible(gc())
    reset <- reset_hwm()
    before_rss <- rss("VmRSS")
    before <- proc.time()
    eval(code[[mode]]$block,env)
    elapsed <- proc.time()-before
    metrics[[length(metrics)+1L]] <- data.frame(pair=pair,mode=mode,
      user_seconds=unname(elapsed[[1]]+elapsed[[4]]),
      system_seconds=unname(elapsed[[2]]+elapsed[[5]]),
      actual_cpu_seconds=unname(sum(elapsed[c(1,2,4,5)])),wall_seconds=unname(elapsed[[3]]),
      rss_before_kib=before_rss,rss_after_kib=rss("VmRSS"),peak_rss_kib=rss("VmHWM"),
      peak_is_phase_specific=reset,calls=nrow(env$calls),filtered_variants=nrow(env$germline_vcf_variants))
    write.table(dplyr::bind_rows(metrics),"metrics.tsv",sep="\t",quote=FALSE,row.names=FALSE)
    # Persist BOTH outputs; the full variant table feeds later global stats.
    qs2::qs_save(list(calls=env$calls,germline_vcf_variants=env$germline_vcf_variants),
                  sprintf("pair%02d.%s.qs2",pair,mode))
    rm(env)
    invisible(gc())
  }
  baseline <- qs2::qs_read(sprintf("pair%02d.baseline.qs2",pair))
  candidate <- qs2::qs_read(sprintf("pair%02d.candidate.qs2",pair))
  comparisons[[length(comparisons)+1L]] <- data.frame(pair=pair,
    calls_identical=identical(baseline$calls,candidate$calls),
    full_variants_identical=identical(baseline$germline_vcf_variants,candidate$germline_vcf_variants))
  write.table(dplyr::bind_rows(comparisons),"comparisons.tsv",sep="\t",quote=FALSE,row.names=FALSE)
  if(!identical(baseline,candidate)) stop("Exact full-block differential comparison failed in pair ",pair)
  rm(baseline,candidate)
  invisible(gc())
}
capture.output(sessionInfo(),file="session.txt")
