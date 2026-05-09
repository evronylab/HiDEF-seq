#!/usr/bin/env -S Rscript --vanilla

#installBSgenome.R:
# Installs required BSgenome.

cat("#### Running installBSgenome ####\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(BSgenomeForge))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(tidyverse))

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
              help="path to YAML configuration file")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$config)){
  stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))

cache_dir <- yaml.config$cache_dir
BSgenome_name <- get_bsgenome_name(yaml.config)
writeLines(BSgenome_name, "BSgenome_name.txt")

cat("DONE\n")

######################
### Install BSgenome reference
######################
cat("## Installing BSgenome reference...")

#Configure cache dir as the path for the BSgenome installation
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths(cache_dir)

#Check only the configured cache. Downstream scripts load the forged package from cache_dir.
if(length(find.package(BSgenome_name, lib.loc=cache_dir, quiet=TRUE)) == 0){
	cat("\n")
	cat("> BSgenome package not found in cache:", BSgenome_name, "\n")
	cat("> Forging BSgenome package from:", yaml.config$genome_fasta, "\n")

	genome_fasta_basename <- basename(yaml.config$genome_fasta)
	genome_part <- genome_fasta_basename %>%
		str_replace_all(pattern = "[^0-9a-zA-Z.]", replacement = "")

	circular_chromosomes <- yaml.config$circular_chromosomes
	if(is.null(circular_chromosomes)){
		circular_chromosomes <- character(0)
	}else{
		circular_chromosomes <- circular_chromosomes %>%
			str_split_1(",") %>%
			str_trim
	}

	twobit_path <- tempfile(tmpdir=getwd(), pattern=str_c(genome_part, "."), fileext=".2bit")
	BSgenomeForge::fastaTo2bit(
		origfile = yaml.config$genome_fasta,
		destfile = twobit_path
	)

	pkg_dir <- BSgenomeForge::forgeBSgenomeDataPkgFromTwobitFile(
		filepath = twobit_path,
		organism = yaml.config$genome_organism,
		provider = "user",
		genome = genome_fasta_basename,
		pkg_maintainer = "HiDEF-seq <hidef-seq@example.invalid>",
		pkg_author = "HiDEF-seq",
		circ_seqs = circular_chromosomes
	)

	if(basename(pkg_dir) != BSgenome_name){
		stop(str_c(
			"ERROR: Forged BSgenome package name did not match derived package name. Expected ",
			BSgenome_name,
			" but got ",
			basename(pkg_dir),
			"."
		), call.=FALSE)
	}

	pkg_tarball <- devtools::build(pkg_dir, path=getwd(), vignettes=FALSE, manual=FALSE, quiet=FALSE)
	devtools::check_built(pkg_tarball, cran=FALSE, force_suggests=FALSE, manual=FALSE, error_on="error", quiet=FALSE)

	devtools::install_local(pkg_tarball, dependencies=FALSE, upgrade="never", build=FALSE, quiet=FALSE)

	if(length(find.package(BSgenome_name, lib.loc=cache_dir, quiet=TRUE)) == 0){
		stop(str_c("ERROR: ", BSgenome_name, " was not installed into cache_dir: ", cache_dir), call.=FALSE)
	}

	cat("DONE\n")
}else{
  cat("Previously installed\n")
}
