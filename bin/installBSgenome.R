#!/usr/bin/env -S Rscript --vanilla

#installBSgenome.R:
# Installs required BSgenome.

cat("#### Running installBSgenome ####\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(BSgenome))
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

cache_dir <- yaml.config$cache_dir

cat("DONE\n")

######################
### Install BSgenome reference
######################
cat("## Installing BSgenome reference...\n")

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
}else{
  cat("> Already installed\n")
}

cat("DONE\n")