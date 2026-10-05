# Sourceable validation-only binding of saved configuration to reviewed input.
# The policy pins complete effective YAML bytes and the exact helper-key set.
read_bound_yaml <- function(path) {
  # Publication suffixes can be .yaml_config.tsv; force YAML format explicitly.
  suppressWarnings(configr::read.config(path, file.type='yaml'))
}

validate_config_binding <- function(actual, side, policy, source_only=FALSE, artifact_path=NULL) {
  stopifnot(side %in% c('reference','candidate'), is.list(actual), !anyDuplicated(names(actual)))
  bindings <- policy$config_helper_bindings
  stopifnot(isTRUE(bindings$ready))
  expected <- bindings[[side]]
  if(!is.null(artifact_path)) {
    selected <- Filter(function(entry) identical(entry$side,side) &&
      identical(normalizePath(entry$path,mustWork=TRUE),normalizePath(artifact_path,mustWork=TRUE)),
      bindings$historical_effective_configs)
    stopifnot(length(selected)<=1L)
    if(length(selected)) expected <- selected[[1]]
  }
  stopifnot(is.character(expected$effective_config),length(expected$effective_config)==1L,
            is.character(expected$sha256),length(expected$sha256)==1L,
            grepl('^[0-9a-f]{64}$',expected$sha256),
            identical(digest::digest(file=expected$effective_config,algo='sha256'),expected$sha256))
  reviewed <- read_bound_yaml(expected$effective_config)
  metadata_keys <- unlist(policy$config_added_metadata,use.names=FALSE)
  keys <- unlist(expected$helper_keys,use.names=FALSE)
  if(is.null(keys)) keys <- character()
  stopifnot(!anyDuplicated(keys),all(keys %in% metadata_keys),
            setequal(intersect(names(reviewed),metadata_keys),keys))
  if(source_only) {
    stopifnot(!any(names(actual) %in% metadata_keys))
  } else {
    stopifnot(setequal(intersect(names(actual),metadata_keys),keys),
              identical(actual[keys],reviewed[keys]))
  }
  scientific <- setdiff(names(reviewed),metadata_keys)
  stopifnot(setequal(setdiff(names(actual),metadata_keys),scientific))
  bound <- scientific_compare(reviewed[scientific],actual[scientific],ignore=character())
  stopifnot(!bound$counts[['failure']],!bound$counts[['review']])
  invisible(TRUE)
}
