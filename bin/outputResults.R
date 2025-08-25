#!/usr/bin/env -S Rscript --vanilla

#outputResults.R:
# Outputs results

cat("#### Running outputResults ####\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(sigfit))
suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(tidyverse))

######################
### Load configuration
######################
cat("## Loading configuration...\n")

#General options
options(warn=2) #Stop script for any warnings

#Command line arguments
option_list = list(
	make_option(c("-c", "--config"), type = "character", default=NULL,
							help="path to YAML configuration file"),
	make_option(c("-s", "--sample_id_toanalyze"), type = "character", default=NULL,
							help="sample_id to analyze"),
	make_option(c("-g", "--chromgroup_toanalyze"), type = "character", default=NULL,
	            help="chromgroup to analyze"),
	make_option(c("-v", "--filtergroup_toanalyze"), type = "character", default=NULL,
	            help="filtergroup to analyze"),
	make_option(c("-f", "--files"), type = "character", default=NULL,
							help="comma-separated filterCalls qs2 files"),
		make_option(c("-o", "--output_basename"), type = "character", default=NULL,
							help="output basename")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$config) | is.na(opt$sample_id_toanalyze) | is.na(opt$chromgroup_toanalyze) | is.na(opt$filtergroup_toanalyze) | is.na(opt$files) | is.na(opt$output_basename) ){
	stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))
sample_id_toanalyze <- opt$sample_id_toanalyze
chromgroup_toanalyze <- opt$chromgroup_toanalyze
filtergroup_toanalyze <- opt$filtergroup_toanalyze
filterCallsFiles <- opt$files %>% str_split_1(",") %>% str_trim
output_basename <- opt$output_basename

#Load the BSgenome reference
suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name,character.only=TRUE,lib.loc=yaml.config$cache_dir))

#Load miscellaneous configuration parameters
 #chromosomes to analyze
chroms_toanalyze <- yaml.config$chromgroups %>%
	enframe %>%
	unnest_wider(value) %>%
	filter(chromgroup == !!chromgroup_toanalyze) %>%
	pull(chroms) %>%
	str_split_1(",") %>%
	str_trim

 #sensitivity thresholds
sensitivity_thresholds <- yaml.config$chromgroups %>%
	enframe %>%
	unnest_wider(value) %>%
	filter(chromgroup == !!chromgroup_toanalyze) %>%
	unnest_wider(sensitivity_thresholds) %>%
	select(-name,-chroms)

 #individual_id of this sample_id
individual_id <- yaml.config$samples %>%
  bind_rows %>%
  filter(sample_id == sample_id_toanalyze) %>%
  pull(individual_id)

#Display basic configuration parameters
cat("> Processing:\n")
cat("    individual_id:",individual_id,"\n")
cat("    sample_id:",sample_id_toanalyze,"\n")
cat("    chromgroup:",chromgroup_toanalyze,"\n")
cat("    filtergroup:",filtergroup_toanalyze,"\n")

cat("DONE\n")

######################
### Define custom functions
######################
#Function to sum two RleLists, including handling of chromosomes only present in one of the RleLists.
sum_RleList <- function(a, b) {
	seqs_union <- union(names(a), names(b))
	
	seqs_union %>%
		map(function(nm){
			seq_in_a <- nm %in% names(a)
			seq_in_b <- nm %in% names(b)
			if(seq_in_a && seq_in_b){
				ra <- a[[nm]]
				rb <- b[[nm]]
				if (length(ra) != length(rb)) {
					stop(sprintf("Length mismatch for '%s': %d vs %d", nm, length(ra), length(rb)))
				}
				ra + rb
			}else if(seq_in_a){
				a[[nm]]
			}else{
				b[[nm]]
			}
		}) %>%
			set_names(seqs_union) %>%
			RleList(compress=TRUE)
}

#Order of trinucleotide context labels
trint_subs_labels <- c(
	"ACA>AAA","ACC>AAC","ACG>AAG","ACT>AAT","CCA>CAA","CCC>CAC","CCG>CAG","CCT>CAT",
	"GCA>GAA","GCC>GAC","GCG>GAG","GCT>GAT","TCA>TAA","TCC>TAC","TCG>TAG","TCT>TAT",
	"ACA>AGA","ACC>AGC","ACG>AGG","ACT>AGT","CCA>CGA","CCC>CGC","CCG>CGG","CCT>CGT",
	"GCA>GGA","GCC>GGC","GCG>GGG","GCT>GGT","TCA>TGA","TCC>TGC","TCG>TGG","TCT>TGT",
	"ACA>ATA","ACC>ATC","ACG>ATG","ACT>ATT","CCA>CTA","CCC>CTC","CCG>CTG","CCT>CTT",
	"GCA>GTA","GCC>GTC","GCG>GTG","GCT>GTT","TCA>TTA","TCC>TTC","TCG>TTG","TCT>TTT",
	"ATA>AAA","ATC>AAC","ATG>AAG","ATT>AAT","CTA>CAA","CTC>CAC","CTG>CAG","CTT>CAT",
	"GTA>GAA","GTC>GAC","GTG>GAG","GTT>GAT","TTA>TAA","TTC>TAC","TTG>TAG","TTT>TAT",
	"ATA>ACA","ATC>ACC","ATG>ACG","ATT>ACT","CTA>CCA","CTC>CCC","CTG>CCG","CTT>CCT",
	"GTA>GCA","GTC>GCC","GTG>GCG","GTT>GCT","TTA>TCA","TTC>TCC","TTG>TCG","TTT>TCT",
	"ATA>AGA","ATC>AGC","ATG>AGG","ATT>AGT","CTA>CGA","CTC>CGC","CTG>CGG","CTT>CGT",
	"GTA>GGA","GTC>GGC","GTG>GGG","GTT>GGT","TTA>TGA","TTC>TGC","TTG>TGG","TTT>TGT"
	)

genome_freqs_labels <- str_sub(trint_subs_labels,1,3)

#Order of indel context labels
indel_labels <- c(
	"1:Del:C:0","1:Del:C:1","1:Del:C:2","1:Del:C:3","1:Del:C:4","1:Del:C:5",
	"1:Del:T:0","1:Del:T:1","1:Del:T:2","1:Del:T:3","1:Del:T:4","1:Del:T:5",
	"1:Ins:C:0","1:Ins:C:1","1:Ins:C:2","1:Ins:C:3","1:Ins:C:4","1:Ins:C:5",
	"1:Ins:T:0","1:Ins:T:1","1:Ins:T:2","1:Ins:T:3","1:Ins:T:4","1:Ins:T:5",
	"2:Del:R:0","2:Del:R:1","2:Del:R:2","2:Del:R:3","2:Del:R:4","2:Del:R:5",
	"3:Del:R:0","3:Del:R:1","3:Del:R:2","3:Del:R:3","3:Del:R:4","3:Del:R:5",
	"4:Del:R:0","4:Del:R:1","4:Del:R:2","4:Del:R:3","4:Del:R:4","4:Del:R:5",
	"5:Del:R:0","5:Del:R:1","5:Del:R:2","5:Del:R:3","5:Del:R:4","5:Del:R:5",
	"2:Ins:R:0","2:Ins:R:1","2:Ins:R:2","2:Ins:R:3","2:Ins:R:4","2:Ins:R:5",
	"3:Ins:R:0","3:Ins:R:1","3:Ins:R:2","3:Ins:R:3","3:Ins:R:4","3:Ins:R:5",
	"4:Ins:R:0","4:Ins:R:1","4:Ins:R:2","4:Ins:R:3","4:Ins:R:4","4:Ins:R:5",
	"5:Ins:R:0","5:Ins:R:1","5:Ins:R:2","5:Ins:R:3","5:Ins:R:4","5:Ins:R:5",
	"2:Del:M:1","3:Del:M:1","3:Del:M:2","4:Del:M:1","4:Del:M:2","4:Del:M:3",
	"5:Del:M:1","5:Del:M:2","5:Del:M:3","5:Del:M:4","5:Del:M:5")

#All possible trinucleotides
trinucleotides_64 <- expand.grid(
	c("A","C","G","T"),
	c("A","C","G","T"),
	c("A","C","G","T")
) %>%
	unite("tri",everything(),sep="") %>%
	pull(tri)

trinucleotides_32_pyr <- expand.grid(
	c("A","C","G","T"),
	c("C","T"),
	c("A","C","G","T")
) %>%
	unite("tri",everything(),sep="") %>%
	pull(tri)

trinucleotides_32_pur <- trinucleotides_64 %>% setdiff(trinucleotides_32_pyr)

#Function to reduce 64 to 32 trinucleotide frequency with central pyrimidine. Input is integer array with named elements that results from the trinucleotideFrequency function of Biostrings.
trinucleotide64to32 <- function(x){
	left_join(
		tibble(
			tri = trinucleotides_32_pyr,
			count_pyr = x[trinucleotides_32_pyr]
		),
		tibble(
			tri = trinucleotides_32_pur %>%
				DNAStringSet %>%
				reverseComplement %>%
				as.character,
			count_pur = x[trinucleotides_32_pur]
		),
		by = "tri"
	) %>%
		mutate(count = count_pyr + count_pur) %>%
		select(tri,count) %>%
		deframe
}

#Function to convert indelwald spectrum (produced by indel.spectrum) to sigfit format
indelwald.to.sigfit <- function(indelwald.spectrum){
	indelwald.spectrum %>%
		map(c) %>%
		unlist %>%
		as_tibble %>%
		na.omit %>% ###CHECK
		set_names("count") %>%
		bind_cols(label=indel_labels,.) %>%
		pivot_wider(names_from=label,values_from=count)
}

######################
### Create set of high-quality germline variants for sensitivity analysis
######################
cat("## Creating set of high-quality germline variants for sensitivity analysis...\n")

#Load germline VCF variants and filter to retain high-confidence variants
germline_vcf_variants <- qs_read(str_c(yaml.config$cache_dir,"/",individual_id,".germline_vcf_variants.qs2")) %>%
	as_tibble %>%
	
	#filter per sensitivity_threshold settings
	filter(
		(Depth >= quantile(Depth, sensitivity_thresholds$SBS_min_Depth_quantile) & call_class == "SBS") |
			(Depth >= quantile(Depth, sensitivity_thresholds$indel_min_Depth_quantile) & call_class == "indel"),
		(VAF >= sensitivity_thresholds$SBS_min_VAF & call_class == "SBS") |
			(VAF >= sensitivity_thresholds$indel_min_VAF & call_class == "indel"),
		(VAF <= sensitivity_thresholds$SBS_max_VAF & call_class == "SBS") |
			VAF <= sensitivity_thresholds$indel_max_VAF & call_class == "indel",
		(GQ >= sensitivity_thresholds$SBS_min_GQ_quantile & call_class == "SBS") | (GQ >= sensitivity_thresholds$indel_min_GQ_quantile & call_class == "indel"),
		(QUAL >= sensitivity_thresholds$SBS_min_QUAL_quantile & call_class == "SBS") |
			(QUAL >= sensitivity_thresholds$indel_min_QUAL_quantile & call_class == "indel")
	)

sensitivity_thresholds$sensitivity_genotype

#variants in all germline vcfs

#Remove chrX and chrY
yaml.config$sex_chromosomes

if is.null(gnomad_sensitivity_vcf) -> set sensitivity to 1.0

Variants in gnomad_sensitivity_vcf
 -> use load_vcf yaml.config$gnomad_sensitivity_vcf

cat("DONE\n")

######################
### Load data from filterCalls files
######################
cat("## Loading data from filterCalls files...\n")

#Load basic configuration only from first chunk, since identical in all chunks.
yaml.config <- filterCallsFiles[1] %>% pluck("config","yaml.config")
run_metadata <- filterCallsFiles[1] %>% pluck("config","run_metadata")
genome_chromgroup.gr <- filterCallsFiles[1] %>% pluck("config","genome_chromgroup.gr")
call_types_toanalyze <- filterCallsFiles[1] %>% pluck("config","call_types")

#Create list for finalCalls
finalCalls <- list()

#Loop over all filterCallsFiles
for(i in seq_along(filterCallsFiles)){
	cat("> Chunk:",i,"\n")
	
	filterCallsFile <- qs_read(filterCallsFiles[i])
	
	#Load final calls that pass all filters
	finalCalls[[i]] <- filterCallsFile %>%
		pluck("calls") %>%
		filter(
			call_toanalyze == TRUE,
			if_all(contains("passfilter"), ~ .x == TRUE)
		)
	
	#Filtered read coverage of the genome
	if(i == 1){
		bam.gr.filtertrack.bytype <- filterCallsFile %>%
			pluck("bam.gr.filtertrack.bytype") %>%
			mutate(
				bam.gr.filtertrack.coverage = bam.gr.filtertrack %>% map(coverage)
			) %>%
			select(-bam.gr.filtertrack)
	}else{
		bam.gr.filtertrack.bytype <- bam.gr.filtertrack.bytype %>%
			left_join(
				filterCallsFile %>%
					pluck("bam.gr.filtertrack.bytype") %>%
					mutate(
						bam.gr.filtertrack.coverage_add = bam.gr.filtertrack %>% map(coverage)
					) %>%
					select(-bam.gr.filtertrack),
				by = names(.) %>% setdiff("bam.gr.filtertrack.coverage")
			) %>%
			mutate(
				bam.gr.filtertrack.coverage = map2(
					bam.gr.filtertrack.coverage,
					bam.gr.filtertrack.coverage_add,
					function(x,y){sum_RleList(x,y)}
				)
			) %>%
			select(-bam.gr.filtertrack.coverage_add)
	}
	
	#Count how many opportunities there was to detect each sensitivity variant set in this chunk's non-gnomad filter tracker
	
	#Load germline calls that pass all filters and intersect with sensitivity variant set to calculate running sensitivity

	
	#molecule_stats
	
	#region_genome_filter_stats
	
}

rm(filterCallsFile)
invisible(gc())

cat("DONE\n")

######################
### Output configuration parameters
######################
cat("> Outputting:\n")

cat("    Configuration parameters...")
#yaml.config, run_metadata

cat("DONE\n")

######################
### Output analysis statistics
######################
cat("    Analysis statistics...:")

cat("DONE\n")

######################
### Output coverage per genome site
######################
cat("    Coverage per genome site...")

#Coverage per genome site

cat("DONE\n")

######################
### Output trinucleotide background counts
######################
cat("    Trinucleotide background counts...")

#Number of trinucleotide counts of interrogated bases and base pairs

cat("DONE\n")

######################
### Output filtered calls
######################
cat("    Filtered calls...")

#Germline and somatic
cat("DONE\n")

######################
### Calculate and output sensitivity
######################
cat("    Sensitivity...")

if sensitivity_threshold$min_sensitivity_variants or set sensitivity to 1

Correct final sensitivity if used het sensitivity_threshold$sensitivity_genotype, up to max of 1.0

cat("DONE\n")

######################
### Output call frequencies
######################
cat("    Call frequencies...")
#Number of interrogated bases and base pairs
#all and unique, observed vs genome corrected vs sensitivity corrected vs both corrected
#upper and lower poisson conf int
#Interrogated bases, Interrogated base pairs

cat("DONE\n")

######################
### Output call spectra
######################
cat("    Call spectra...")
#Plots and tables of observed and corrected, for all and for unique counts

cat("DONE\n")

######################
### Output estimated mutation error rate
######################
cat("    Estimated mutation error rate...")
#For ssDNA mismatches only, for each channel and total

cat("DONE\n")