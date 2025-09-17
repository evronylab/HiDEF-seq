#!/usr/bin/env -S Rscript --vanilla

#installBSgenome.R:
# Installs required BSgenome.

cat("#### Running installBSgenome ####\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(tidyverse))

######################
### Load configuration
######################
cat("## Loading configuration...")

#General options
options(warn=2) #Stop script for any warnings

#Command line arguments
option_list = list(
  make_option(c("-c", "--config"), type = "character", default=NULL,
              help="path to YAML configuration file")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$config)){
  stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))

cache_dir <- yaml.config$cache_dir

cat("DONE\n")

######################
### Install BSgenome reference
######################
cat("## Installing BSgenome reference...")

#Configure cache dir as the path for the BSgenome installation
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths(cache_dir)

#Check if BSgenome_name is already installed or available to be installed in Bioconductor, and if not, install the genome from BSgenome_file
if(!yaml.config$BSgenome$BSgenome_name %in% installed.genomes()){
	if(yaml.config$BSgenome$BSgenome_name %in% available.genomes() %>% suppressMessages){
		BiocManager::install(yaml.config$BSgenome$BSgenome_name, lib = cache_dir)
	}else if(!is.null(yaml.config$BSgenome$BSgenome_file)){
		install.packages(yaml.config$BSgenome$BSgenome_file, repos = NULL, lib = cache_dir)
	}else{
		stop("ERROR: Must specify either BSgenome_name that is in available.genomes() or a BSgenome_file!", call.=FALSE)
	}
	
	cat("DONE\n")
}else{
  cat("Previously installed\n")
}

######################
### Output duckdb of genome trinucleotide sequences
######################
cat("## Extracting trinucleotide sequences...")

genome_trinuc_file <- str_c(yaml.config$BSgenome$BSgenome_name,".bed.gz")

if(!file.exists(str_c(cache_dir,"/",genome_trinuc_file))){
	
	#Extract sequences for all bases (except contig edges) from genome with seqkit and write to BED format (bgzipped, tabix indexed)
	tmpseqs <- tempfile(tmpdir=getwd(),pattern=".")
	
	invisible(
		system(
			paste(
				yaml.config$seqkit_bin,"sliding -S '' -s1 -W3",yaml.config$genome_fasta,"|",
				yaml.config$seqkit_bin,"seq -u |", #convert to upper case
				yaml.config$seqkit_bin,"fx2tab -Q |",
				"awk -F '[:\\-\\t]' 'BEGIN {OFS=\"\t\"}{print $1, $2, $2+1, $4}' |", #$2 is first position of 1-based trinucleotide position
				yaml.config$bgzip_bin,"-c >",genome_trinuc_file,"&&",
				yaml.config$tabix_bin,"-s 1 -b 2 -e 3",genome_trinuc_file
			),
		intern = FALSE)
	)
	
	#Copy to cache_dir
	if(!file.exists(str_c(cache_dir,"/",genome_trinuc_file))){
		file_move(genome_trinuc_file, cache_dir)
	}
	
	cat("DONE\n")
}else{
	cat("Previously extracted\n")
}