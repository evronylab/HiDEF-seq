#!/usr/bin/env -S Rscript --vanilla

#processGermlineVCF.R:
# Processes germline VCF.
# Usage: processGermlineVCF.R -c [configuration.yaml]

cat("#### Running processGermlineVCF ####\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs2))

######################
### Load custom shared functions
######################
source(Sys.which("sharedFunctions.R"))

######################
### Load configuration
######################
cat("## Loading configuration...")

#General options
options(warn=2) #Stop script for any warnings

#Command line arguments
option_list = list(
  make_option(c("-c", "--config"), type = "character", default=NULL,
              help="path to YAML configuration file"),
  make_option(c("-i", "--individual_id"), type = "character", default=NULL,
              help="individual_id to process"),
  make_option(c("-o", "--output"), type = "character", default=NULL,
  						help="output qs2 file")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$config) | is.null(opt$individual_id) | is.null(opt$output)){
  stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))

individual_id_toprocess <- opt$individual_id
outputFile <- opt$output

cache_dir <- yaml.config$cache_dir

#Load the BSgenome reference
suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name,character.only=TRUE,lib.loc=yaml.config$cache_dir))

cat("DONE\n")

######################
### Load VCF variants
######################
#Load all VCF variants, not just analyzed chromosomes, since this will be used for all chunks

cat("#### Loading VCF variants...\n")

#Load list of germline VCFs for all configured individuals
vcf_files <- yaml.config$individuals %>%
	modify_tree(leaf = as.character) %>%
	bind_rows %>%
	unnest_wider(germline_vcf_files)

#Load VCF files
cat("> Processing individual:",individual_id_toprocess,"\n")

germline_vcf_variants <- list()

vcf_files_individual <- vcf_files %>% filter(individual_id == individual_id_toprocess)

for(i in 1:nrow(vcf_files_individual)){
	
	germline_vcf_file <- vcf_files_individual %>% pluck("germline_vcf_file",i)
	germline_vcf_type <- vcf_files_individual %>% pluck("germline_vcf_type",i)

  cat(">> Processing VCF file:",germline_vcf_file,"...")
  
  germline_vcf_variants[[i]] <- load_vcf(
	  	vcf_file = germline_vcf_file,
	  	genome_fasta = yaml.config$genome_fasta,
	  	BSgenome_name = yaml.config$BSgenome$BSgenome_name,
	  	bcftools_bin = yaml.config$bcftools_bin
  	) %>%
  	mutate(
  		germline_vcf_file = factor(!!germline_vcf_file),
  		germline_vcf_type = factor(!!germline_vcf_type)
  	)

  cat("DONE\n")
}

#Save variants to file
qs_save(
	germline_vcf_variants %>% GRangesList %>% unlist,
	outputFile
  )

cat("DONE\n")