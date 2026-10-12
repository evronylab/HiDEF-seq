#!/usr/bin/env Rscript
# Run in the pipeline R environment on a compute node:
# Rscript tests/test_burden_call_store.R [repository root]
# Exact native-R contracts from the qualified DuckDB call-store fixtures.
# These tests exercise successful schema, ordering, counting, and QS2 handoff.
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(DBI)
  library(duckdb)
  library(qs2)
})
options(warn = 2)
args <- commandArgs(TRUE)
repo <- normalizePath(if(length(args)) args[[1L]] else ".")
for(expr in parse(file.path(repo, "bin", "sharedFunctions.R"))) {
  if(is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
     identical(expr[[2L]], as.name("new_burden_call_store"))) eval(expr)
}
stopifnot(exists("new_burden_call_store", mode = "function"))

local({
previous_directory <- getwd()
test_directory <- tempfile("burden-call-store-tests-", tmpdir = previous_directory)
stopifnot(dir.create(test_directory))
setwd(test_directory)
on.exit(setwd(previous_directory))
checks <- 0L
record <- function(label, actual, expected) {
  if(!identical(actual, expected)) {
    print(list(case = label, actual = actual, expected = expected))
    stop("Burden call-store contract differs: ", label, call. = FALSE)
  }
  checks <<- checks + 1L
  invisible(NULL)
}
check_case <- function(label, fun) fun()
roundtrip <- function(label, chunks) check_case(label, function() {
  store <- new_burden_call_store(buffer_rows = 2L)
  on.exit({store$close(); store$cleanup()})
  for(x in chunks) store$append_table("finalCalls", x)
  expected <- bind_rows(chunks)
  actual <- store$read_table("finalCalls")
  record(label, actual, expected)
  record(paste0(label, "_repeat"), store$read_table("finalCalls"), actual)
})

make <- function(f, named = TRUE, scalar = seq_along(f)) {
  x <- tibble(group = f, scalar = scalar, alt_plus_strand = rep('T', length(f)),
    qual = lapply(seq_along(f), function(i) if (i == 1L) integer() else c(i, NA_integer_)))
  if (named) names(x$alt_plus_strand) <- paste0('row', seq_along(f))
  x
}
a <- make(factor('a', levels = c('unused', 'a')))
b <- make(factor(c('b', 'a', NA), levels = c('b', 'a', 'next')))
cases <- list(
  list(a[0, ], a, b), list(a, b), list(a, make(c('z', 'a'))), list(a[0, ]),
  list(a[0, ], make(c('c', 'd'), FALSE)), list(make(c('c', 'd'), FALSE), b),
  list(a[0, ], b[0, ]), list(a[0, ], make(character(), FALSE)),
  list(make(character(), FALSE), a[0, ]),
  list(make(c('a', 'b'), scalar = c(TRUE, NA)), make('c', scalar = 3L), make('d', scalar = 4.5)))
for (i in seq_along(cases)) {
  store <- new_burden_call_store(buffer_rows = 1L)
  for (x in cases[[i]]) store$append_table('finalCalls', x)
  record(paste0('raw_', i), store$read_table('finalCalls'), bind_rows(cases[[i]]))
  store$close()
  store$cleanup()
}

# Exact original raw-chunk bind/analysis_chunk/projection/filter boundaries.
call_chunk <- function(index, empty = FALSE) {
  x <- tibble(run_id = factor(c('r2', 'r1'), levels = c('unused', if (index == 1L) c('r2', 'r1') else c('r1', 'r3', 'r2'))),
    zm = c(2L, 1L), call_class = factor(c('SBS', 'SBS'), levels = c('indel', 'SBS')),
    call_type = factor(c('SBS', 'SBS')), SBSindel_call_type = factor(c('mutation', 'mismatch-ss')),
    seqnames = factor(c('chr2', 'chr1'), levels = c('chr1', 'chr2', 'chrX')),
    start = c(22L, 11L), end = c(22L, 11L), ref_plus_strand = c('A', NA_character_),
    alt_plus_strand = c('T', 'C'), call_toanalyze = c(TRUE, FALSE),
    germline_vcf_types_detected = c('', ''), germline_vcf_files_detected = c('', ''),
    germline_vcf.passfilter = c(TRUE, TRUE), unusual.PASSFILTER = c(TRUE, TRUE),
    qual = list(integer(), c(NA_integer_, 19L, 7L)))
  names(x$alt_plus_strand) <- c('original-2', 'original-1')
  if (empty) x[0, ] else x
}
for (empty in c(FALSE, TRUE)) {
  chunks <- list(call_chunk(1L, empty), call_chunk(2L, empty))
  store <- new_burden_call_store(buffer_rows = 1L)
  for (i in seq_along(chunks)) store$append_calls(chunks[[i]], chunks[[i]], i)
  combined <- bind_rows(chunks, .id = 'analysis_chunk') %>% mutate(analysis_chunk = as.integer(analysis_chunk))
  projected <- combined %>% select(-c('germline_vcf_types_detected', 'germline_vcf_files_detected',
    'call_toanalyze', contains('passfilter')))
  record(paste0('projection_', empty), store$project_final(), projected)
  record(paste0('public_final_', empty), store$read_table('finalCalls',
    columns = setdiff(names(combined), 'call_toanalyze'), where = '"call_toanalyze" = TRUE'),
    combined %>% filter(call_toanalyze == TRUE) %>% select(-call_toanalyze))
  record(paste0('public_germline_', empty), store$read_table('germlineVariantCalls'), combined %>% select(-call_toanalyze))
  store$close()
  store$cleanup()
}

# SQL grouping must match dplyr distinct/count including NA keys and duplicates
# across chunks. Then exercise the original left_join's exact result contract.
for (empty in c(FALSE, TRUE)) {
  x <- tibble(run_id = factor(c('r1', 'r1', 'r2', 'r1', NA, NA)), zm = c(1L, 1L, 1L, 2L, 3L, 3L),
    call_type = factor(c('SBS', 'SBS', 'SBS', 'indel', 'SBS', 'SBS'), levels = c('unused', 'indel', 'SBS')),
    seqnames = factor(c('chr2', 'chr2', 'chr2', 'chr1', NA, NA), levels = c('chr1', 'chr2', 'chrX')),
    start = c(20L, 20L, 20L, 10L, NA_integer_, NA_integer_),
    ref_plus_strand = c('A', 'A', 'A', 'G', NA_character_, NA_character_))
  names(x$ref_plus_strand) <- paste0('allele-', seq_len(nrow(x)))
  if (empty) x <- x[0, ]
  keys <- c('call_type', 'seqnames', 'start', 'ref_plus_strand')
  store <- new_burden_call_store(buffer_rows = 1L)
  halves <- if (empty) list(x, x) else list(x[c(1L, 3L, 5L), ], x[c(2L, 4L, 6L), ])
  halves <- map2(halves, seq_along(halves), function(value, chunk)
    add_column(value, analysis_chunk = rep.int(as.integer(chunk), nrow(value)), .before = 1L))
  for (part in halves) store$append_table('germlineVariantCalls', part)
  combined <- bind_rows(halves)
  expected <- combined %>% distinct(across(all_of(c('run_id', 'zm', keys)))) %>%
    count(across(all_of(keys)), name = 'num_zm_detected')
  actual <- store$sensitivity_counts(keys)
  record(paste0('counts_', empty), actual, expected)
  variants <- distinct(combined, across(all_of(keys)))
  record(paste0('sensitivity_left_join_', empty),
    variants %>% left_join(actual, by = keys) %>% mutate(num_zm_detected = replace_na(num_zm_detected, 0)),
    variants %>% left_join(expected, by = keys) %>% mutate(num_zm_detected = replace_na(num_zm_detected, 0)))
  store$close()
  store$cleanup()
}

# Empty producer schemas, scalar promotion, names, and duplicate multiplicity.
named <- function(x,n=rep('',length(x))) setNames(x,n)
empty0 <- tibble(group=factor(), value=integer())
some <- tibble(group=factor(c('b',NA,'a'),levels=c('unused','b','a')),value=1:3)
roundtrip('producer_zero_level_factor_all_empty',list(empty0,empty0))
roundtrip('producer_zero_level_factor_then_nonempty',list(empty0,some))
roundtrip('producer_nonempty_then_zero_level_factor',list(some,empty0))
a <- tibble(group=factor(c('b',NA,'a'),levels=c('unused','b','a')),value=named(c('T',NA,'C'),c('b','','a')),
  quality=list(integer(),c(NA_integer_,2L,2L),c(-1L,.Machine$integer.max)))
b <- tibble(group=factor(c('c','b'),levels=c('new-unused','c','b')),value=c('A','T'),quality=list(NA_integer_,3:1))
cases <- list(empty_before=list(a[0,],a,b),empty_after=list(a,b,a[0,]),
  reverse_union=list(b,a),all_empty_unions=list(a[0,],b[0,]),
  named_empty_unnamed_data=list(a[0,],b),unnamed_data_named_empty=list(b,a[0,]),
  unnamed_empty_named_data=list(b[0,],a),named_data_unnamed_empty=list(a,b[0,]),
  empty_named_then_empty_unnamed=list(a[0,],b[0,]),empty_unnamed_then_empty_named=list(b[0,],a[0,]),
  factor_to_character=list(tibble(v=factor(c('a',NA),levels=c('unused','a'))),tibble(v=c('b',NA))),
  character_to_factor=list(tibble(v=c('b',NA)),tibble(v=factor(c('a',NA),levels=c('unused','a')))),
  logical_integer_double=list(tibble(v=c(TRUE,NA)),tibble(v=c(2L,NA_integer_)),tibble(v=c(2.5,NA_real_))),
  double_integer_logical=list(tibble(v=c(2.5,NA_real_)),tibble(v=c(2L,NA_integer_)),tibble(v=c(TRUE,NA))),
  named_every_scalar=list(tibble(l=named(c(TRUE,NA),c('yes','missing')),i=named(c(1L,NA_integer_),c('one','na')),
    d=named(c(1.5,NA_real_),c('x','y')),f=named(factor(c('b',NA),levels=c('unused','b')),c('fb','fn')),
    q=named(list(integer(),c(1L,NA_integer_)),c('q0','q1')))),
  duplicate_multiplicity=list(a,a,b,b))
for(label in names(cases)) roundtrip(label,cases[[label]])

# Order differs deliberately from sorted key order, with NA identities and keys,
# duplicated molecules across chunks, unused factor levels, and named alleles.
check_case('grouping_projection',function() {
  x <- tibble(run_id=factor(c('r2','r1','r2',NA,NA,'r3','r1','r1'),levels=c('unused','r3','r2','r1')),
    zm=c(2L,1L,2L,NA_integer_,NA_integer_,1L,1L,2L),
    call_type=factor(c('SBS','indel','SBS','SBS','SBS','indel','indel','indel'),levels=c('unused','SBS','indel')),
    seqnames=factor(c('chr2','chr1','chr2',NA,NA,'chr1','chr1','chr1'),levels=c('unused','chr2','chr1')),
    start=c(20L,10L,20L,NA_integer_,NA_integer_,10L,10L,10L),
    ref_plus_strand=named(c('A','G','A',NA,NA,'G','G','G'),paste0('n',1:8)),
    call_toanalyze=c(TRUE,FALSE,TRUE,NA,TRUE,FALSE,FALSE,TRUE),
    germline_vcf_types_detected=rep('',8),germline_vcf_files_detected=rep('',8),
    xxx.PASSFILTER=rep(TRUE,8),qual=rep(list(c(1L,NA_integer_)),8))
  chunks <- list(x[0,],x[c(3,5,8,2),],x[c(1,4,6,7),],x[0,])
  store <- new_burden_call_store(buffer_rows=1L)
  on.exit({store$close(); store$cleanup()})
  for(i in seq_along(chunks))store$append_calls(chunks[[i]],chunks[[i]],i)
  combined <- bind_rows(lapply(seq_along(chunks),function(i)add_column(chunks[[i]],analysis_chunk=rep.int(as.integer(i),nrow(chunks[[i]])),.before=1L)))
  record('projection_exact',store$project_final(),select(combined,-c(germline_vcf_types_detected,germline_vcf_files_detected,call_toanalyze,contains('passfilter'))))
  record('public_filter_with_na',store$read_table('finalCalls',setdiff(names(combined),'call_toanalyze'),'"call_toanalyze" = TRUE'),select(filter(combined,call_toanalyze==TRUE),-call_toanalyze))
  record('raw_duplicate_multiplicity',store$read_table('finalCalls'),combined)
  for(keys in list(c('call_type','seqnames','start','ref_plus_strand'),c('ref_plus_strand','start'),c('call_type'))) {
    tag <- paste(keys,collapse='_')
    expected <- combined %>% distinct(across(all_of(c('run_id','zm',keys)))) %>% count(across(all_of(keys)),name='num_zm_detected')
    actual <- store$sensitivity_counts(keys)
    record(paste0('count_',tag),actual,expected)
    record(paste0('count_repeat_',tag),store$sensitivity_counts(keys),actual)
  }
})

# Factor names and levels must match native bind_rows, including empty chunks.
f<-function(value,levels,n=NULL) {x<-factor(value,levels=levels);if(!is.null(n))names(x)<-n;x}
named_factor_cases<-list(
 single_named=list(f('a','a','first')),
 zero_to_one=list(f(NA,character(),'missing'),f('a','a','present')),
 one_to_zero=list(f('a','a','present'),f(NA,character(),'missing')),
 same_levels=list(f('a',c('a','b'),'first'),f('b',c('a','b'),'second')),
 different_levels=list(f('a','a','first'),f('b','b','second')),
 grow_levels=list(f('a','a','first'),f('b',c('a','b'),'second')),
 shrink_levels=list(f('a',c('a','b'),'first'),f('b','b','second')),
 different_unused=list(f('a',c('unused1','a'),'first'),f('a',c('unused2','a'),'second')),
 reversed_level_order=list(f('a',c('b','a'),'first'),f('b',c('a','b'),'second')),
 empty_grow=list(f(character(),'a',character()),f('b',c('a','b'),'second')),
 named_then_unnamed=list(f('a','a','first'),f('b','b')),
 unnamed_then_named=list(f('a','a'),f('b','b','second')),
 factor_to_character=list(f('a','a','first'),setNames('b','second')),
 character_to_factor=list(setNames('b','first'),f('a','a','second')),
 three_factors=list(f('a','a','first'),f('b','b','second'),f('c',c('c','unused'),'third')),
 factor_factor_character=list(f('a','a','first'),f('b','b','second'),setNames('c','third')),
 character_factor_factor=list(setNames('c','first'),f('a','a','second'),f('b','b','third')),
 factor_character_factor=list(f('a','a','first'),setNames('c','second'),f('b','b','third')),
 named_na_empty_named=list(f(NA,character(),'na'),f(character(),'unused',character()),f('a','a','third')),
 all_empty_named_factors=list(f(character(),'a',character()),f(character(),'b',character()),f(character(),character(),character())),
 factor_empty_character=list(f('a','a','first'),setNames(character(),character())),
 character_empty_factor=list(setNames('a','first'),f(character(),'unused',character())),
 factor_factor_empty_character=list(f('a','a','first'),f('b','b','second'),setNames(character(),character())))

named_factor_cases$single_empty_named_factor<-list(f(character(),'a',character()))
named_factor_cases$single_empty_zero_named_factor<-list(f(character(),character(),character()))
named_factor_cases$empty_factor_named_character<-list(f(character(),'unused',character()),setNames('a','last'))
named_factor_cases$two_empty_zero_named_factors<-list(f(character(),character(),character()),f(character(),character(),character()))

shared_factor<-f(c('a','b'),c('unused','a','b'),c('first','second'))
named_factor_cases$shared_factor_twice<-list(shared_factor,shared_factor)
named_factor_cases$identical_separate_factors<-list(f(c('a','b'),c('unused','a','b'),c('first','second')),f(c('a','b'),c('unused','a','b'),c('first','second')))
named_factor_cases$same_names_different_values<-list(f(c('a','b'),c('unused','a','b'),c('first','second')),f(c('b','a'),c('unused','a','b'),c('first','second')))
named_factor_cases$same_names_same_levels_different_order<-list(f(c('a','b'),c('unused','a','b'),c('first','second')),f(c('a','b'),c('a','unused','b'),c('first','second')))
named_factor_cases$shared_factor_three<-list(shared_factor,shared_factor,shared_factor)
named_factor_cases$shared_factor_then_char<-list(shared_factor,shared_factor,setNames('c','third'))
named_factor_cases$named_missing_column_names<-list(f(c('a','b'),c('a','b'),c(NA_character_,'second')),setNames('c','third'))

unnamed_seed<-factor(c('a','b'),levels=c('unused','a','b'))
named_factor_cases$unnamed_then_separate_then_named_shared_levels<-list(unnamed_seed,factor(c('b','a'),levels=c('unused','a','b')),setNames(unnamed_seed,c('late1','late2')))
named_factor_cases$empty_then_separate_then_named_shared_levels<-list(unnamed_seed[0],factor(character(),levels=c('unused','a','b')),setNames(unnamed_seed,c('late1','late2')))

for(label in names(named_factor_cases)) roundtrip(paste0('factor_names_',label),lapply(named_factor_cases[[label]],function(x)tibble(v=x)))
for(parts in 1:2) check_case(paste0('named_factor_calls_',parts),function(){
  x<-tibble(run_id=factor(c('r2','r1','r2')),zm=c(2L,1L,2L),
    call_type=setNames(factor(c('SBS','indel','SBS'),levels=c('unused','SBS','indel')),c('s1','i1','s2')),
    seqnames=factor(c('chr2','chr1','chr2')),start=c(20L,10L,20L),
    call_toanalyze=c(TRUE,FALSE,TRUE),germline_vcf_types_detected=rep('',3),
    germline_vcf_files_detected=rep('',3),foo.passfilter=rep(TRUE,3))
  chunks<-rep(list(x),parts)
  store<-new_burden_call_store(buffer_rows=1L)
  on.exit({store$close(); store$cleanup()})
  for(i in seq_along(chunks))store$append_calls(chunks[[i]],chunks[[i]],i)
  expected<-bind_rows(chunks,.id='analysis_chunk') %>% mutate(analysis_chunk=as.integer(analysis_chunk))
  record(paste0('named_factor_append_calls_',parts),store$read_table('finalCalls'),expected)
  record(paste0('named_factor_projection_',parts),store$project_final(),select(expected,-c(call_toanalyze,germline_vcf_types_detected,germline_vcf_files_detected,contains('passfilter'))))
  keys<-c('call_type','seqnames','start')
  counts<-expected %>% distinct(across(all_of(c('run_id','zm',keys)))) %>% count(across(all_of(keys)),name='num_zm_detected')
  record(paste0('named_factor_sensitivity_',parts),store$sensitivity_counts(keys),counts)
})

roundtrip('genuine_na_character_names',list(tibble(v=setNames(c('a','b'),c(NA_character_,'second'))),tibble(v='c'),tibble(v=setNames('d','fourth'))))
roundtrip('genuine_na_character_names_factor_middle',list(tibble(v=setNames(c('a','b'),c(NA_character_,'second'))),tibble(v=factor('c')),tibble(v=setNames('d','fourth'))))

# The backing database and all DuckDB directories are private to this task.
# The saved QS2 remains independently readable after successful cleanup.
local({
  store <- new_burden_call_store(buffer_rows = 1L)
  other <- new_burden_call_store(buffer_rows = 1L)
  path <- environment(store$append_table)$path
  other_path <- environment(other$append_table)$path
  record("private_directory_in_task", dirname(path), normalizePath(getwd()))
  record("independent_private_directories", path != other_path, TRUE)
  con <- environment(store$append_table)$con
  for(setting in c("extension_directory", "secret_directory", "temp_directory")) {
    actual <- DBI::dbGetQuery(con, paste0("SELECT current_setting('", setting, "')"))[[1L]]
    suffix <- c(extension_directory = "extensions", secret_directory = "secrets", temp_directory = "temp")[[setting]]
    record(paste0("private_", setting), actual, file.path(path, suffix))
  }
  input <- tibble(value = factor(c("a", NA), levels = c("unused", "a")),
                  quality = list(integer(), c(1L, NA_integer_)))
  store$append_table("finalCalls", input)
  result <- list(finalCalls = store$read_table("finalCalls"))
  store$close()
  qs2::qs_save(result, "result.qs2")
  store$cleanup()
  record("closed_database_removed", dir.exists(path), FALSE)
  record("other_store_preserved", file.exists(file.path(other_path, "calls.duckdb")), TRUE)
  record("standalone_qs2_after_cleanup", qs2::qs_read("result.qs2"), list(finalCalls = bind_rows(input)))
  other$close()
  other$cleanup()
  record("only_public_output_remains", list.files(all.files = TRUE, no.. = TRUE), "result.qs2")
})
cat("PASS", checks, "exact burden call-store checks\n")
setwd(previous_directory)
stopifnot(unlink(test_directory, recursive = TRUE) == 0L)
})
