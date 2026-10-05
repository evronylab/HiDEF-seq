#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
bench <- normalizePath(args[[1]])
out <- normalizePath(args[[2]])
source(file.path(bench,'compare_qs2.R'))
source(file.path(bench,'compare_config_bindings.R'))
reference <- list(threshold=0.5,scientific=list(order=c('a','b')))
helpers <- list(reference_cache_dir='/cache/reference',reference_summary_file='/cache/summary.qs2',
                germline_coverage_filters=list(list(individual_id='one',threshold=15L,file='/cache/mask.qs2')),
                cache_artifacts=list('referenceSummary.qs2'='/cache/summary.qs2'))
candidate <- c(reference,helpers)
ref_path <- file.path(out,'reference.yaml'); cand_path <- file.path(out,'candidate.yaml')
yaml::write_yaml(reference,ref_path); yaml::write_yaml(candidate,cand_path)
# Use the same parser as the real pipeline for both expected and actual values.
reference <- read_bound_yaml(ref_path); candidate <- read_bound_yaml(cand_path)
anchor <- function(path,keys) list(effective_config=path,sha256=digest::digest(file=path,algo='sha256'),helper_keys=as.list(keys))
policy <- list(reference_config=ref_path,candidate_config=ref_path,config_replacements=list(),
               config_added_metadata=as.list(names(helpers)),qs_ignore_paths=as.list(paste0('/yaml.config/',names(helpers))),
               config_helper_bindings=list(ready=TRUE,reference=anchor(ref_path,character()),
                 candidate=anchor(cand_path,names(helpers)),historical_effective_configs=list()))
reject <- function(expr) stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),'error'))
validate_config_binding(reference,'reference',policy)
validate_config_binding(candidate,'candidate',policy)
validate_config_binding(reference,'candidate',policy,source_only=TRUE)
for(key in names(helpers)) {
  bad <- candidate
  if(key=='germline_coverage_filters') bad[[key]][[1]]$threshold <- 16L
  else if(key=='cache_artifacts') bad[[key]][[1]] <- '/wrong/product.qs2'
  else bad[[key]] <- '/wrong/path'
  reject(validate_config_binding(bad,'candidate',policy))
}
bad <- candidate; bad$reference_summary_file <- NULL
reject(validate_config_binding(bad,'candidate',policy))
bad <- candidate; bad$threshold <- 0.6
reject(validate_config_binding(bad,'candidate',policy))
bad <- reference; bad$threshold <- 0.6
reject(validate_config_binding(bad,'reference',policy))
bad_policy <- policy; bad_policy$config_helper_bindings$candidate$sha256 <- paste(rep('0',64),collapse='')
reject(validate_config_binding(candidate,'candidate',bad_policy))
bad_policy <- policy; bad_policy$config_helper_bindings$ready <- FALSE
reject(validate_config_binding(candidate,'candidate',bad_policy))
old <- candidate; old$reference_summary_file <- '/old/summary.qs2'
old_path <- file.path(out,'historical-effective.yaml');yaml::write_yaml(old,old_path)
old <- read_bound_yaml(old_path)
policy$config_helper_bindings$historical_effective_configs <- list(c(anchor(old_path,names(helpers)),
  list(path=old_path,side='candidate')))
validate_config_binding(old,'candidate',policy,artifact_path=old_path)
reject(validate_config_binding(old,'candidate',policy))
policy_path <- file.path(out,'policy.json')
jsonlite::write_json(policy,policy_path,auto_unbox=TRUE,pretty=TRUE)
qs2::qs_save(list(yaml.config=reference,values=c(1,2)),file.path(out,'reference.qs2'))
qs2::qs_save(list(yaml.config=candidate,values=c(1,2)),file.path(out,'candidate.qs2'))
run <- function(arguments,log) system2(file.path(R.home('bin'),'Rscript'),
  vapply(c('--vanilla',file.path(bench,'compare_qs2.R'),arguments),shQuote,character(1)),
  stdout=file.path(out,log),stderr=file.path(out,log))
stage <- file.path(out,'components')
stopifnot(run(c('stage',file.path(out,'reference.qs2'),stage),'stage.log')==0L)
flags <- c(paste0('--config-policy=',policy_path),paste0('--ignore=/yaml.config/',names(helpers)))
stopifnot(run(c('compare',stage,file.path(out,'candidate.qs2'),file.path(out,'exact.tsv'),flags),'exact.log')==0L)
for(key in names(helpers)) {
  bad <- candidate
  if(key=='germline_coverage_filters') bad[[key]][[1]]$threshold <- 16L
  else if(key=='cache_artifacts') bad[[key]][[1]] <- '/wrong/product.qs2'
  else bad[[key]] <- '/wrong/path'
  path <- file.path(out,paste0(key,'.qs2'))
  qs2::qs_save(list(yaml.config=bad,values=c(1,2)),path)
  stopifnot(run(c('compare',stage,path,paste0(path,'.tsv'),flags),paste0(key,'.log'))!=0L)
}
qs2::qs_save(list(values=c(1,2)),file.path(out,'missing.qs2'))
stopifnot(run(c('stage',file.path(out,'missing.qs2'),file.path(out,'missing-components')),'missing-stage.log')==0L)
stopifnot(run(c('compare',file.path(out,'missing-components'),file.path(out,'missing.qs2'),
                file.path(out,'missing.tsv'),flags),'missing.log')!=0L)
cat('PASS: exact source/helper bindings; all four helper mutations/drop rejected; scientific drift and stale hashes rejected; historical exact binding; disk-staged QS binding and missing-component gate\n')
