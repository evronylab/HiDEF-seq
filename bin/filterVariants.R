#!/usr/bin/env -S Rscript --vanilla

#filterVariants.R:
# Filters variants

cat("#### Running filterVariants ####\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(tidyverse))

######################
### Load configuration
######################
cat("## Loading configuration...")

#General options
options(datatable.showProgress = FALSE)
options(warn=2) #Stop script for any warnings

#Command line arguments
option_list = list(
	make_option(c("-c", "--config"), type = "character", default=NULL,
							help="path to YAML configuration file"),
	make_option(c("-s", "--sample_id"), type = "character", default=NULL,
	            help="sample_id to analyze"),
	make_option(c("-f", "--file"), type = "character", default=NULL,
							help="path to extractVariants qs2 file"),
	make_option(c("-o", "--output"), type = "character", default=NULL,
							help="output qs2 file")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$config) | is.na(opt$file) | is.na(opt$output) ){
	stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))
sample_id <- opt$sample_id
extractVariantsFile <- opt$file
outputFile <- opt$output

#Load the BSgenome reference
suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name,character.only=TRUE,lib.loc=yaml.config$cache_dir))

#Load miscellaneous configuration parameters
 #individual_id of this sample_id
individual_id <- yaml.config$samples %>%
  bind_rows %>%
  filter(sample_id == !!sample_id) %>%
  pull(individual_id)

 #cache_dir
cache_dir <- yaml.config$cache_dir

 #germline vcf filter parameters
germline_vcf_filters_config <- yaml.config$germline_vcf_filtergroups %>%
  enframe(name=NULL) %>%
  unnest_wider(value,transform=list(vcf_SNV_FILTERS=list, vcf_INDEL_FILTERS=list))

 #variant types
variant_types <- yaml.config$variant_types %>%
  enframe(name=NULL) %>%
  unnest_wider(value)
  
 #variant filter group parameters
variant_filtergroups_config <- yaml.config$variant_filtergroups %>%
  enframe(name=NULL) %>%
  unnest_wider(value)
  
 #region filter parameters
region_filters_config <- list()
region_filters_config[["read_filters"]] <- yaml.config$region_filters %>%
  map("read_filters") %>%
  flatten %>%
  enframe(name=NULL) %>%
  unnest_wider(value) %>%
  mutate(region_filter_threshold_file = str_c(cache_dir,"/",basename(region_filter_file),".bin",binsize,".",threshold,".bw"))

region_filters_config[["bin_filters"]] <- yaml.config$region_filters %>%
  map("bin_filters") %>%
  flatten %>%
  enframe(name=NULL) %>%
  unnest_wider(value) %>%
  mutate(region_filter_threshold_file = str_c(cache_dir,"/",basename(region_filter_file),".bin",binsize,".",threshold,".bw"))

#Output data lists
read_region_filters <- list()
genome_region_filters <- list()
stats <- list() 
cat("DONE\n")

######################
### Load extracted variants
######################
cat("## Loading extracted variants:",sample_id,"->",extractVariantsFile,"...")

extractedVariants <- qs_read(extractVariantsFile)

cat("DONE\n")

######################
### Initial whole-molecule filters
######################
cat("## Applying initial whole-molecule filters...")


cat("DONE\n")

######################
### Germline VCF variant filters
######################
cat("## Applying germline VCF variant filters...")

#Load germline VCF filter data
germline_vcf_variants <- qs_read(str_c(cache_dir,"/",individual_id,".germline_vcf_variants.qs2"))

cat("DONE\n")

######################
### Post-germline VCF variant filtering whole-molecule filters
######################
cat("## Applying post-germline VCF variant filtering whole-molecule filters...")


cat("DONE\n")

######################
### Region filters
######################
cat("## Applying region filters...")

#Also noN filter

cat("DONE\n")

######################
### Base-quality filters
######################
cat("## Applying base-quality filters...")

cat("DONE\n")

######################
### Molecule position filters
######################
cat("## Applying molecule position filters...")

cat("DONE\n")

######################
### Germline VCF indel region filters
######################
cat("## Applying germline VCF indel region filters...")

cat("DONE\n")

######################
### Read indel region filters
######################
cat("## Applying read indel region filters...")

cat("DONE\n")

######################
### Germline BAM filter
######################
cat("## Applying germline BAM filter...")

cat("DONE\n")

######################
### Subread filters
######################
cat("## Applying subread filters...")

cat("DONE\n")