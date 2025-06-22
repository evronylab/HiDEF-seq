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
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs2))

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
              help="individual_id to process")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$config)){
  stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))

individual_id_toprocess <- opt$individual_id
cache_dir <- yaml.config$cache_dir

cat("DONE\n")

######################
### Load VCF variants
######################
#Load all VCF variants, not just analyzed chromosomes, since this will be used for all chunks

cat("#### Loading VCF variants...\n")

#Load list of germline VCFs for all configured individuals
vcf_files <- yaml.config$individuals %>%
  map(~ .x %>% keep(names(.x) %in% c("individual_id","germline_vcf_files"))) %>%
  enframe(name=NULL) %>%
  unnest_wider(value) %>%
  unnest_longer(germline_vcf_files) %>%
  unnest_wider(germline_vcf_files)

#Load VCF files
cat("> Processing individual:",individual_id_toprocess,"\n")
    
germline_vcf_variants <- list()

vcf_files_individual <- vcf_files %>% filter(individual_id == individual_id_toprocess)

for(i in 1:nrow(vcf_files_individual)){
	
	germline_vcf_file <- vcf_files_individual %>% pluck("germline_vcf_file",i)
	germline_vcf_filtergroup <- vcf_files_individual %>% pluck("germline_vcf_filtergroup",i)

  cat(">> Processing VCF file:",germline_vcf_file,"...")

	#Split multi-allelic sites (bcftools norm -m -both -f [fastaref]), and keep only variants in chromosomes present in this bamfile chunk.
	tmpvcf <- tempfile(tmpdir=getwd(),pattern=".")
	system(paste("/bin/bash -c",shQuote(paste(yaml.config$bcftools_bin,"norm -m -both -f",yaml.config$genome_fasta,germline_vcf_file," 2>/dev/null >",tmpvcf))))

  #Load vcf file
	 #Remove ALT == "*" alleles that indicate overlapping deletions (not needed because every deletion already has a separate vcf entry)
	 #Annotate with Depth (sum of AD) and VAF (AD of variant / sum of Depth).
	 #Annotate SNVs vs indels
	 #For deletions, change ranges to reflect position of deletion. For insertions: width=-1, and the insertion is immediately to the right of the first REF base.
	 #Keep only columns CHROM, POS, REF, ALT, QUAL, FILTER, GT, GQ.
	
  vcf <- read.vcfR(tmpvcf,convertNA=FALSE,verbose=FALSE)
  file.remove(tmpvcf) %>% invisible

  germline_vcf_variants[[i]] <- cbind(
  	vcf@fix,
  	data.frame(
	  	GT=as.character(extract.gt(vcf,element="GT",IDtoRowNames=FALSE)),
	  	GQ=as.numeric(extract.gt(vcf,element="GQ",as.numeric=TRUE,IDtoRowNames=FALSE)),
	  	AD=as.character(extract.gt(vcf,element="AD",IDtoRowNames=FALSE))
  	)) %>%
  	as_tibble %>%
  	filter(ALT!="*") %>%
  	separate(AD,c("AD1","AD2"),",",convert=TRUE) %>%
  	mutate(Depth=AD1+AD2, VAF=AD2/Depth) %>%
  	mutate(
  		variant_type = case_when(
  			nchar(REF)==1 & nchar(ALT)==1 ~ "SBS_mutation",
  			nchar(REF) != nchar(ALT) ~ "indel_mutation"
  			) %>%
  			factor,
  		indel_type = case_when(
  			nchar(REF) - nchar(ALT) > 0 ~ "deletion",
  			nchar(REF) - nchar(ALT) < 0 ~ "insertion"
  		) %>%
  			factor,
  		CHROM = factor(CHROM),
  		POS = as.numeric(POS),
  		end_refspace = if_else(indel_type == "deletion", POS + nchar(REF) - nchar(ALT), POS),
  		start_refspace = if_else(variant_type == "indel_mutation", POS + 1, POS),
  		QUAL = as.numeric(QUAL),
  		germline_vcf_file = factor(!!germline_vcf_file),
  		germline_vcf_filtergroup = factor(!!germline_vcf_filtergroup)
  	) %>%
  	select(-c(POS,ID,INFO,AD1,AD2)) %>%
  	rename(
  		seqnames = CHROM,
  		ref = REF,
  		alt = ALT
  	)
  
  rm(vcf)
  cat("DONE\n")
}

#Save variants to file
qs_save(
	germline_vcf_variants %>% bind_rows,
  str_c(individual_id_toprocess,".germline_vcf_variants.qs2")
  )

cat("DONE\n")