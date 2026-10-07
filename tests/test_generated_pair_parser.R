#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
options(warn=2)
args <- commandArgs(trailingOnly=TRUE)
repo <- if(length(args)) args[[1L]] else '.'
old <- function(x) separate(x, token, sep='_', into=c('left','right'))
nodes <- as.list(parse(file.path(repo, 'bin/extractCalls.R')))
helper <- Filter(function(node) is.call(node) && identical(node[[1L]], as.name('<-')) &&
  identical(node[[2L]], as.name('separate_generated_pair')), nodes)
stopifnot(length(helper) == 1L)
eval(helper[[1L]])
new <- function(x) separate_generated_pair(x, token, into=c('left','right'))
capture <- function(fun,x) tryCatch(list(success=TRUE,value=fun(x)),error=function(e) list(success=FALSE,message=conditionMessage(e)))
cases <- list(
  ordinary=tibble(before=1:4,token=c('123_+','9_-','-2_17',NA_character_),after=factor(c('a','b','a','b'))),
  empty=tibble(before=integer(),token=character(),after=factor(levels=c('a','b'))),
  leading_empty=tibble(token=c('_+','1_')),
  no_separator=tibble(token='12'),
  extra_separator=tibble(token='12_+_bad'),
  empty_string=tibble(token=''),
  all_missing=tibble(token=c(NA_character_,NA_character_)),
  factor_token=tibble(token=factor(c('1_+','2_-',NA))),
  collision=tibble(left=c('old','old'),token=c('1_+','2_-'),right=1:2),
  data_frame=data.frame(before=1:2,token=c('1_+','2_-'),after=3:4),
  grouped=tibble(group=factor(c('b','a','b'),levels=c('b','a','c')),token=c('1_+','2_-',NA_character_)) %>% group_by(group,.drop=FALSE),
  grouped_token=tibble(token=c('1_+','2_-'),value=1:2) %>% group_by(token),
  rowwise=tibble(token=c('1_+','2_-'),value=1:2) %>% rowwise()
)
cases$named_rows <- data.frame(token=c('1_+','2_-'),row.names=c('a','b'))
cases$custom_attribute <- structure(data.frame(token=c('1_+','2_-')), note='preserve me')
cases$custom_subclass <- structure(data.frame(token=c('1_+','2_-')), class=c('custom','data.frame'))
results <- lapply(names(cases),function(name) {
  a <- capture(old,cases[[name]]); b <- capture(new,cases[[name]])
  list(name=name,old_success=a$success,new_success=b$success,
       exact=identical(a$success,b$success) && (!a$success || identical(a$value,b$value)),
       old_class=if(a$success) class(a$value) else a$message,
       new_class=if(b$success) class(b$value) else b$message,
       details=if(a$success && b$success && !identical(a$value,b$value)) capture.output(all.equal(a$value,b$value)) else character())
})
stopifnot(all(vapply(results[!names(cases) %in% c('no_separator','extra_separator','empty_string')],
  function(x) x$old_success && x$new_success, logical(1))))
stopifnot(all(vapply(results,`[[`,logical(1),'exact')))
cat('PASS: generated-field parser preserves values, types, missing data, grouping and data-frame attributes\n')
