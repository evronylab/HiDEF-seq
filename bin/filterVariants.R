#!/usr/bin/env -S Rscript --vanilla

#extractVariants.R:
# Loads and formats aligned ccs bamFile in HiDEF-seq format RDS file that includes all required alignment and variant information for analysis.  
# Usage: extractVariants.R -c [configuration.yaml] -b [bamFile]

cat("#### Running extractVariants ####\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs2))

######################
### Load configuration
######################
cat("## Loading configuration...")

#General options
options(datatable.showProgress = FALSE)

#Command line arguments
option_list = list(
	make_option(c("-c", "--config"), type = "character", default=NULL,
							help="path to YAML configuration file"),
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
extractVariantsFile <- opt$file
outputFile <- opt$output

#Load the BSgenome reference
suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name,character.only=TRUE,lib.loc=yaml.config$cache_dir))

#Comma-separated chromosomes to analyze across all chromgroups
chroms_to_analyze <- yaml.config$chromgroups %>%
	map_vec("chroms") %>%
	map(str_split_1,",") %>%
	unlist

cat("DONE\n")