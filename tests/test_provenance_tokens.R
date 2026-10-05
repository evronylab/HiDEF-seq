#!/usr/bin/env Rscript
# Tiny validation-only regression; no pipeline data or external tools required.
args <- commandArgs(trailingOnly=TRUE)
source(if(length(args)) args[[1]] else 'scripts/benchmark/compare_qs2.R')
rule <- list(reference=c('/old/a.vcf.gz','/old/b.vcf.gz'),
             candidate=c('/new/a.vcf.gz','/new/b.vcf.gz'))
path <- '/germlineVariantCalls/germline_vcf_files_detected'
rules <- setNames(list(rule),path)
a <- data.frame(position=1:4,germline_vcf_files_detected=c(
  '/old/a.vcf.gz,/old/b.vcf.gz,/old/a.vcf.gz','/old/b.vcf.gz','',NA_character_))
b <- a
b$germline_vcf_files_detected <- c(
  '/new/a.vcf.gz,/new/b.vcf.gz,/new/a.vcf.gz','/new/b.vcf.gz','',NA_character_)
check <- function(x,y) {
  xn <- normalize_provenance_column(x,'germlineVariantCalls',rules,'reference')
  yn <- normalize_provenance_column(y,'germlineVariantCalls',rules,'candidate')
  stopifnot(all(scientific_compare(xn,yn,ignore=character())$counts==0L))
  TRUE
}
fails <- function(x,y) stopifnot(inherits(try(check(x,y),silent=TRUE),'try-error'))
stopifnot(check(a,b))
for(value in c('/third/a.vcf.gz','/new/b.vcf.gz,/new/a.vcf.gz,/new/a.vcf.gz',
               '/new/a.vcf.gz,/new/b.vcf.gz','/new/a.vcf.gz,/new/b.vcf.gz,/new/a.vcf.gz,')) {
  changed <- b; changed$germline_vcf_files_detected[[1]] <- value; fails(a,changed)
}
changed <- a; changed$germline_vcf_files_detected[[1]] <- '/third/a.vcf.gz'; fails(changed,b)
changed <- b; changed$position[[1]] <- 100L; fails(a,changed)
changed <- b; changed$germline_vcf_files_detected[[4]] <- ''; fails(a,changed)
changed <- b; attr(changed,'assembly') <- 'different'; fails(a,changed)
changed <- b; changed$germline_vcf_files_detected <- factor(changed$germline_vcf_files_detected); fails(a,changed)
changed <- b; changed$germline_vcf_files_detected <- NULL; fails(a,changed)
# Identical substrings elsewhere receive no normalization.
a$other <- '/old/a.vcf.gz'; b$other <- '/new/a.vcf.gz'; fails(a,b)
stopifnot(identical(normalize_provenance_column(a,'unlistedComponent',rules,'reference'),a))
# Vector names and table attributes are retained and compared.
named <- setNames(c('/new/a.vcf.gz',NA_character_),c('first','second'))
stopifnot(identical(normalize_provenance_tokens(named,rule,'candidate'),
                    setNames(c('/old/a.vcf.gz',NA_character_),c('first','second'))))
cat('PASS: exact provenance column tokens; order, duplicate multiplicity, NA, types, attributes and all other science preserved\n')
