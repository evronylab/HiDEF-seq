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
suppressPackageStartupMessages(library(vcfR))
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

 #sex chromosomes
sex_chromosomes <- yaml.config$sex_chromosomes %>%
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
			RleList(compress=FALSE)
}

#Function to calculate genome coverage for bam.gr.filtertrack_bytype. If two bam.gr.filtertrack_bytypes are provided, it calculates genome coverage only for the second and adds it to the first. Also removes the bam.gr.filtertrack column that is no longer necessary.
sum_bam.gr.filtertracks <- function(bam.gr.filtertrack1, bam.gr.filtertrack2=NULL){
	
	if(is.null(bam.gr.filtertrack2)){
		bam.gr.filtertrack1 %>%
			mutate(
				bam.gr.filtertrack.coverage = bam.gr.filtertrack %>% map(coverage)
			) %>%
			select(-bam.gr.filtertrack)
		
	}else{
		bam.gr.filtertrack1 %>%
			left_join(
				bam.gr.filtertrack2 %>%
					mutate(
						bam.gr.filtertrack.coverage = bam.gr.filtertrack %>% map(coverage)
					) %>%
					select(-bam.gr.filtertrack),
				by = names(.) %>% setdiff("bam.gr.filtertrack.coverage"),
				suffix = c("",".2")
			) %>%
			mutate(
				bam.gr.filtertrack.coverage = map2(
					bam.gr.filtertrack.coverage,
					bam.gr.filtertrack.coverage.2,
					function(x,y){sum_RleList(x,y)}
				)
			) %>%
			select(-bam.gr.filtertrack.coverage.2)
	}
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
### Load data from filterCalls files
######################
cat("## Loading data from filterCalls files...\n> chunk:")

#Create list for finalCalls and molecule_stats
finalCalls <- list()
molecule_stats <- list()

#Loop over all filterCallsFiles
for(i in seq_along(filterCallsFiles)){
	
	#Extract chunk number for logging in case filterCallsFiles are not supplied in order
	chunk_num <- filterCallsFiles[i] %>% str_extract("(?<=chunk)\\d+")
		
	cat(" ", chunk_num, sep="")
	
	#Load filterCallsFile
	filterCallsFile <- qs_read(filterCallsFiles[i])
	
	#Load basic configuration only from first chunk, since identical in all chunks.
	if(i == 1){
		#Basic configuration parameters
		run_metadata <- filterCallsFile %>% pluck("config","run_metadata")
		genome_chromgroup.gr <- filterCallsFile %>% pluck("config","genome_chromgroup.gr")
		call_types_toanalyze <- filterCallsFile %>% pluck("config","call_types")
		region_read_filters_config <- filterCallsFile %>% pluck("config","region_read_filters_config")
		region_genome_filters_config <- filterCallsFile %>% pluck("config","region_genome_filters_config")
		
		#gnomad-related filters that are excluded when calculating sensitivity
		gnomad_filters <- c(
			"germline_vcf.passfilter",
			str_subset(filterCallsFile %>% pluck("calls") %>% colnames, "germline_vcf_indel_region_filter"),
			"max_BAMVariantReads.passfilter",
			"max_BAMVAF.passfilter",
			region_read_filters_config %>%
				filter(is_gnomad_filter==TRUE) %>%
				pull(region_filter_threshold_file) %>%
				basename %>%
				str_c("region_read_filter_",.,".passfilter"),
			region_genome_filters_config %>%
				filter(is_gnomad_filter==TRUE) %>%
				pull(region_filter_threshold_file) %>%
				basename %>%
				str_c("region_genome_filter_",.,".passfilter")
		)
		
		#region_genome_filter_stats
		region_genome_filter_stats <- filterCallsFile %>% pluck("region_genome_filter_stats")
	}

	#Load final calls that pass all filters, ignoring the gnomAD-related filters, since we will later also need the calls filtered only by gnomAD-related filters to calculate sensitivity
	finalCalls[[i]] <- filterCallsFile %>%
		pluck("calls") %>%
		filter(
			call_toanalyze == TRUE,
			if_all(contains("passfilter") & !all_of(gnomad_filters), ~ .x == TRUE)
		)
	
	#Filtered read coverage of the genome, with and without gnomad_filters
	if(i == 1){
		bam.gr.filtertrack.bytype <- sum_bam.gr.filtertracks(
			filterCallsFile %>% pluck("bam.gr.filtertrack.bytype")
		)
		
		bam.gr.filtertrack.except_gnomad_filters.bytype <- sum_bam.gr.filtertracks(
			filterCallsFile %>% pluck("bam.gr.filtertrack.except_gnomad_filters.bytype")
		)
		
	}else{
		bam.gr.filtertrack.bytype <- sum_bam.gr.filtertracks(
			bam.gr.filtertrack.bytype,
			filterCallsFile %>% pluck("bam.gr.filtertrack.bytype")
		)
		
		bam.gr.filtertrack.except_gnomad_filters.bytype <- sum_bam.gr.filtertracks(
			bam.gr.filtertrack.except_gnomad_filters.bytype,
			filterCallsFile %>% pluck("bam.gr.filtertrack.except_gnomad_filters.bytype")
		)
	}

	#molecule_stats
	molecule_stats[[i]] <- filterCallsFile %>% pluck("molecule_stats")

}

#Remove temp objects
rm(filterCallsFile)
invisible(gc())

#Combine finalCalls to one tibble
finalCalls <- finalCalls %>%
	bind_rows(.id = "chunk") %>%
	mutate(chunk = chunk %>% as.integer)

#Combine molecule_stats to one tibble and sum stats. Do this by left_join to molecule_stats[[1]] without value column to preserve the same stats row order in the final tibble
molecule_stats <- molecule_stats[[1]] %>%
	select(run_id,chromgroup,filtergroup,stat) %>%
	left_join(
		molecule_stats %>%
			bind_rows %>%
			group_by(run_id,chromgroup,filtergroup,stat) %>%
			summarize(value = sum(value), .groups = "drop"),
		by = join_by(run_id,chromgroup,filtergroup,stat)
	)

cat(" DONE\n")

######################
### Output configuration parameters
######################
cat("> Analyzing and outputting:\n")

cat("    Configuration parameters...")
#yaml.config, run_metadata

cat("DONE\n")

######################
### Output analysis statistics
######################
cat("    Statistics...:")

molecule_stats

region_genome_filter_stats

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
cat("    Sensitivity...\n")

cat("      > Setting sensitivity to 1 since gnomad_sensitivity_vcf is NULL.\n")

#Create set of high-quality germline variants for sensitivity analysis
if(!is.null(yaml.config$gnomad_sensitivity_vcf)){
	
	#Count number of germline VCFs for the analyzed individual
	num_germline_vcf_files <- yaml.config$individuals %>%
		map(~ .x %>% keep(names(.x) %in% c("individual_id","germline_vcf_files"))) %>%
		enframe(name=NULL) %>%
		unnest_wider(value) %>%
		unnest_longer(germline_vcf_files) %>%
		unnest_wider(germline_vcf_files) %>%
		filter(individual_id == !!individual_id) %>%
		nrow
	
	#Load gnomad_sensitivity_vcf
	gnomad_sensitivity_vcf <- load_vcf(
		vcf_file = yaml.config$gnomad_sensitivity_vcf,
		genome_fasta = yaml.config$genome_fasta,
		BSgenome_name = yaml.config$BSgenome$BSgenome_name,
		bcftools_bin = yaml.config$bcftools_bin
	) %>%
		as_tibble
	
	#Load germline VCF variants and filter to retain high-confidence variants
	germline_vcf_variants <- qs_read(str_c(yaml.config$cache_dir,"/",individual_id,".germline_vcf_variants.qs2")) %>%
		as_tibble %>%
		
		#separate genotypes of each allele
		separate_wider_delim(GT, regex("[^[:digit:].]"), names=c("GT1","GT2"), cols_remove = FALSE) %>%
		
		#group by germline_vcf_file so that filters are calculated separately for each germline VCF.
		group_by(germline_vcf_file) %>% 
		
		#filter per sensitivity_threshold settings
		filter(
			(Depth >= quantile(Depth, sensitivity_thresholds$SBS_min_Depth_quantile) & call_class == "SBS") |
				(Depth >= quantile(Depth, sensitivity_thresholds$indel_min_Depth_quantile) & call_class == "indel"),
			(VAF >= sensitivity_thresholds$SBS_min_VAF & call_class == "SBS") |
				(VAF >= sensitivity_thresholds$indel_min_VAF & call_class == "indel"),
			(VAF <= sensitivity_thresholds$SBS_max_VAF & call_class == "SBS") |
				(VAF <= sensitivity_thresholds$indel_max_VAF & call_class == "indel"),
			(GQ >= quantile(GQ, sensitivity_thresholds$SBS_min_GQ_quantile) & call_class == "SBS") |
				(GQ >= quantile(GQ, sensitivity_thresholds$indel_min_GQ_quantile) & call_class == "indel"),
			(QUAL >= quantile(QUAL, sensitivity_thresholds$SBS_min_QUAL_quantile) & call_class == "SBS") |
				(QUAL >= quantile(QUAL, sensitivity_thresholds$indel_min_QUAL_quantile) & call_class == "indel"),
			(sensitivity_thresholds$sensitivity_genotype == "het" & ((GT1 == "1" & GT2 != "1") | (GT1 != "1" & GT2 == "1"))) |
				(sensitivity_thresholds$sensitivity_genotype == "hom" & (GT1 == "1" & GT2 == "1"))
		) %>%
		ungroup %>%
		
		#keep variants detected in all germline vcfs
		count(seqnames, start, end, ref_plus_strand, alt_plus_strand, call_class, call_type, SBSindel_call_type) %>%
		filter(n == num_germline_vcf_files) %>%
		select(-n) %>%
		
		#keep variants in gnomad_sensitivity_vcf
		semi_join(
			gnomad_sensitivity_vcf,
			by = names(.)
		)
	
	rm(num_germline_vcf_files, gnomad_sensitivity_vcf)
	
	if sensitivity_threshold$min_sensitivity_variants or set sensitivity to 1
	
	#Load germline calls that pass all filters and intersect with sensitivity variant set to calculate sensitivity
	
	#Count how many opportunities there was to detect each sensitivity variant set in this chunk's non-gnomad filter tracker
	
	Take into account sex_chromosomes
	
	Correct final sensitivity if used het sensitivity_threshold$sensitivity_genotype, up to max of 1.0
	
	cat("DONE\n")
	
}else{
	cat("      > Setting sensitivity to 1 since gnomad_sensitivity_vcf is NULL.\n")
	
	sensitivity_sbs <- 1.0
	sensitivity_indel <- 1.0
}

cat("    DONE\n")

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