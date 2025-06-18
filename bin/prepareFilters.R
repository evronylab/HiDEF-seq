#!/usr/bin/env -S Rscript --vanilla

#prepareFilters.R:
# Prepares all required reference genome and filter files for downstream analysis steps.
# Usage: prepareFilters.R -c [configuration.yaml]

cat("#### Running prepareFilters ####\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs2))

######################
### Load configuration
######################
cat("## Loading configuration...")

#Command line arguments
option_list = list(
  make_option(c("-c", "--config"), type = "character", default=NULL,
              help="path to YAML configuration file")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$config)){
  stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))

cache_dir <- yaml.config$prepareFilters_cache_dir

cat("DONE\n")

######################
### Configure BSgenome reference
######################
cat("## Configuring BSgenome reference...")

#Configure cache dir as the path for the BSgenome installation
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths(cache_dir)

#Check if BSgenome_name is already installed or available to be installed in Bioconductor, and if not, install the genome from BSgenome_file
if(!yaml.config$BSgenome$BSgenome_name %in% installed.genomes()){
	if(yaml.config$BSgenome$BSgenome_name %in% available.genomes()){
		BiocManager::install(yaml.config$BSgenome$BSgenome_name, lib=cache_dir)
	}else if(!is.null(yaml.config$BSgenome$BSgenome_file)){
		install.packages(yaml.config$BSgenome$BSgenome_file, repos = NULL, lib = cache_dir)
	}else{
		stop("ERROR: Must specify either BSgenome_name that is in available.genomes() or a BSgenome_file!", call.=FALSE)
	}
}

######################
### Load VCF variants
######################
#Load all variants, not just analyzed chromosomes, since this will be used for all chunks
# 
# cat("#### Loading VCF variants...")
# 
# #Load list of germline VCFs for all configured individuals
# vcf_files <- yaml.config$individuals %>%
#   map(~ .x %>% keep(names(.x) %in% c("individual_id","germline_vcf_files"))) %>%
#   enframe(name=NULL) %>%
#   unnest_wider(value) %>%
#   unnest_longer(germline_vcf_files) %>%
#   unnest_wider(germline_vcf_files)
# 
# #Load germline VCF filter parameters
# vcf_config <- yaml.config$germline_vcf_filtergroups %>%
#   bind_rows %>%
#   group_by(across(-c(vcf_SNV_FILTERS,vcf_INDEL_FILTERS))) %>%
#   summarize(vcf_SNV_FILTERS=list(vcf_SNV_FILTERS),vcf_INDEL_FILTERS=list(vcf_INDEL_FILTERS), .groups="drop")
# 
# #Annotate with each VCF file (make a list of tables, each table with the same number of rows as sbses.df)
# vcf_indels <- list()
# vcf_snvs <- list()
# 
# for(i in germline_vcf_files){
# 
#   cat("     Loading VCF File:",i$germline_vcf_file,"...")
# 
# 	#Split multi-allelic sites (bcftools norm -m -both -f [fastaref]), and keep only variants in chromosomes present in this bamfile chunk.
# 	tmpvcf <- tempfile(tmpdir=getwd(),pattern=".")
# 	chunkchroms <- paste(unique(fwd.ccs.df$rname),collapse=",")
# 	system(paste("/bin/bash -c",shQuote(paste(yaml.config$bcftools_bin,"norm -m -both -f",fastaref,"-r",chunkchroms,i," 2>/dev/null >",tmpvcf))))
# 
#   #Load vcf file, filter to keep only chromosomes being analyzed, and keep only columns CHROM, POS, REF, ALT, QUAL, FILTER, GT, GQ
#   vcffile <- read.vcfR(tmpvcf,convertNA=FALSE,verbose=FALSE)
#   file.remove(tmpvcf)
# 
#   vcftable <- cbind(vcffile@fix,
#   	data.frame(
# 	  	GT=as.character(extract.gt(vcffile,element="GT",convertNA=FALSE,IDtoRowNames=FALSE)),
# 	  	GQ=as.numeric(extract.gt(vcffile,element="GQ",as.numeric=TRUE,convertNA=FALSE,IDtoRowNames=FALSE)),
# 	  	AD=as.character(extract.gt(vcffile,element="AD",convertNA=FALSE,IDtoRowNames=FALSE))
#   	)) %>%
#   	filter(ALT!="*") %>%
#     filter(CHROM %in% chroms) %>%
#   	as.data.frame() %>%
#   	mutate_at(c("POS","QUAL"),as.numeric) %>%
#   	separate(AD,c("AD1","AD2"),",") %>%
#   	mutate_at(c("AD1","AD2"),as.numeric) %>%
#   	mutate(Depth=AD1+AD2, VAF=AD2/Depth) %>%
#   	rename(Ref = REF, Alt = ALT) %>%
#   	select(-c(ID,INFO,AD1,AD2))
# 
#   rm(vcffile)
# 
#   #Separate indels from SNVs and remove '*' Alt alleles, which indicates overlapping deletions, which are not needed because every deletion already has a separate vcf entry. Then annotate with Depth (sum of AD) and VAF (AD of variant / sum of Depth).
#   vcfIndels[[i]] <- vcftable %>%
#   	filter(nchar(Ref)!=nchar(Alt))
# 
#   vcfSNVs[[i]] <- vcftable %>%
#   	filter(nchar(Ref)==1 & nchar(Alt)==1)
# 
#   rm(vcftable)
# 
#   #Convert Indels to Granges object. Also, change ranges for deletion variants to accurately reflect position of deletion. Insertion variants are annotated with width=-1 at the insertion site (insertion is immediately to the right of the first Ref base).
#   vcfIndels[[i]] <- makeGRangesFromDataFrame(vcfIndels[[i]],keep.extra.columns=TRUE,ignore.strand=TRUE,seqnames.field="CHROM",start.field="POS",end.field="POS")
#   lengthdiff <- nchar(vcfIndels[[i]]$Ref) - nchar(vcfIndels[[i]]$Alt)
#   end(vcfIndels[[i]])[lengthdiff>0] <- start(vcfIndels[[i]])[lengthdiff>0] + lengthdiff[lengthdiff>0]
#   start(vcfIndels[[i]]) <- start(vcfIndels[[i]]) + 1
# 
#   cat("DONE\n")
# }
# 
# 
# ######################
# ### Prepare region filters
# ######################
# cat("## Preparing region filters...")
# 
# cat("DONE\n")