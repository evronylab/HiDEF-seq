#!/usr/bin/env -S Rscript --vanilla

#outputResults.R:
# Output results.
	#Also calculates sensitivity-corrected burdens.

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
suppressPackageStartupMessages(library(VariantAnnotation))
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
	make_option(c("-f", "--files"), type = "character", default=NULL,
							help="comma-separated calculateBurdens qs2 files"),
	make_option(c("-o", "--output_basename"), type = "character", default=NULL,
							help="output basename")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$config) | is.na(opt$sample_id_toanalyze) | is.na(opt$files) | is.na(opt$output_basename) ){
	stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))
sample_id_toanalyze <- opt$sample_id_toanalyze
calculateBurdensFiles <- opt$files %>% str_split_1(",") %>% str_trim
output_basename <- opt$output_basename

#Load the BSgenome reference
suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name,character.only=TRUE,lib.loc=yaml.config$cache_dir))

#Load miscellaneous configuration parameters

 #individual_id of this sample_id
individual_id <- yaml.config$samples %>%
  bind_rows %>%
  filter(sample_id == sample_id_toanalyze) %>%
  pull(individual_id)

#Display basic configuration parameters
cat("> Processing:\n")
cat("    individual_id:",individual_id,"\n")
cat("    sample_id:",sample_id_toanalyze,"\n")

cat("DONE\n")

######################
### Define custom functions
######################
#Function to write vcf from finalCalls/germlineCalls
write_vcf_from_calls <- function(calls, BSgenome_name, out_vcf){
	
	#Fixed fields
	CHROM <- calls$seqnames
	POS <- calls$start
	REF <- calls$ref_plus_strand %>% DNAStringSet
	ALT <- calls$alt_plus_strand %>% DNAStringSet
	QUAL <- rep(NA_real_, nrow(calls))
	FILTER <- rep("PASS", nrow(calls))
	
	# Row ranges
	if(nrow(calls)>0){
		calls.gr <- GRanges(
			seqnames = CHROM,
			ranges = IRanges(start = POS, width = REF %>% as.character %>% nchar),
			strand = "*",
			seqinfo = BSgenome_name %>% get %>% seqinfo
		)
	}else{
		calls.gr <- GRanges(seqinfo = BSgenome_name %>% get %>% seqinfo)
	}
	
	fixed <- DataFrame(REF = REF, ALT = ALT, QUAL = QUAL, FILTER = FILTER)
	
	# INFO: everything except the basic VCF columns and fields you want to drop
	info_cols <- names(calls) %>%
		setdiff(c("seqnames","start","ref_plus_strand","alt_plus_strand","qual","strand")) %>%
		str_subset("passfilter$", negate=TRUE)
	
	#Convert factor columns to character
	info_df <- calls %>%
		select(all_of(info_cols)) %>%
		mutate(across(where(is.factor), as.character)) %>%
		as.data.frame %>%
		DataFrame
	
	#Sort by seqnames, start, end
	calls.order <- order(
		calls.gr %>% seqnames %>% as.integer,
		calls.gr %>% start,
		calls.gr %>% end
	)
	calls.gr <- calls.gr[calls.order]
	fixed <- fixed[calls.order, , drop=FALSE]
	info_df <- info_df[calls.order, , drop=FALSE]
	
	#INFO header
	number_map <- function(x){
		if(is.logical(x)){0L}else{1L}
	}
	
	type_map <- function(x){
		if(is.logical(x)){return("Flag")}
		if(is.integer(x)){return("Integer")}
		if(is.double(x)){return("Float")}
		"String"
	}
	
	info_header <- DataFrame(
		Number = vapply(info_df, number_map, integer(1)),
		Type = vapply(info_df, type_map, character(1)),
		Description = info_cols,
		row.names = info_cols
	)
	
	# Build header (fileformat/source/reference/contigs via Seqinfo)
	hdr <- VariantAnnotation::VCFHeader()
	info(hdr) <- info_header
	meta(hdr)$fileformat <- DataFrame(Value = "VCFv4.3", row.names = "fileformat")
	meta(hdr)$source <- DataFrame(Value = "HiDEF-seq", row.names = "source")
	meta(hdr)$reference <- DataFrame(Value = BSgenome_name, row.names = "reference")
	
	vcf <- VariantAnnotation::VCF(
		rowRanges = calls.gr,
		fixed = fixed,
		info = info_df,
		collapsed = FALSE
	)
	header(vcf) <- hdr
	
	VariantAnnotation::writeVcf(vcf, filename = out_vcf, index = TRUE)
}

######################
### Load data from calculateBurdens files
######################
cat("## Loading data from calculateBurdens files...\n> chromgroup/filtergroup:")

#Create lists for data loading

#Loop over all calculateBurdens
for(i in seq_along(calculateBurdensFiles)){
	
	cat(" ", **, sep="")
	
	#Load calculateBurdensFile
	calculateBurdensFile <- qs_read(calculateBurdensFiles[i])
	
	#Load basic configuration only from first chunk, since identical in all chunks.
	if(i == 1){
		#Basic configuration parameters
		run_metadata -> same for all filtergroup/chromgroup
		
	}

}

#Remove temp objects
rm(calculcateBurdensFile)
invisible(gc())

######################
### Output configuration parameters
######################
cat("## Outputting configuration parameters...")

opt$config %>% file.copy(str_c(output_basename,".yaml.config.tsv"))

run_metadata %>% write_tsv(str_c(output_basename,".run_metadata.tsv"))

cat("DONE\n")

######################
### Output filtering statistics
######################
cat("## Outputting filtering statistics...:")

molecule_stats_by_run_id %>% write_tsv(str_c(output_basename,".molecule_stats_by_run_id.tsv"))

molecule_stats_by_analysis_id %>% write_tsv(str_c(output_basename,".molecule_stats_by_analysis_id.tsv"))

region_genome_filter_stats %>% write_tsv(str_c(output_basename,".region_genome_filter_stats.tsv"))

cat("DONE\n")


### output coverage and reference sequences of interrogated genome bases
#Output, with a separate folder for each call_class, call_type, SBSindel_call_type combination
for(i in seq_len(nrow(bam.gr.filtertrack.bytype))){
	
	prefix <- 
		str_c(
			bam.gr.filtertrack.bytype$call_class[i],
			bam.gr.filtertrack.bytype$call_type[i],
			bam.gr.filtertrack.bytype$SBSindel_call_type[i],
			sep="."
		)
	
	dir.create(prefix, showWarnings = FALSE)
	
	bam.gr.filtertrack.bytype %>%
		pluck("bam.gr.filtertrack.coverage",i) %>%
		mutate(name = reftnc_pyr, score = coverage) %>% #rename to to allow output by 'export'
		export(
			con = str_c(
				str_c(prefix,"/",output_basename),
				prefix,
				"bed",
				sep="."
			),
			format = "bed",
			index = TRUE
		)
}

cat("DONE\n")


#Output trinucleotide counts
#bam.gr.filtertracks
for(i in c("bam.gr.filtertrack.reftnc_duplex_pyr","bam.gr.filtertrack.reftnc_both_strands","bam.gr.filtertrack.reftnc_both_strands_pyr")){
	for(j in seq_len(nrow(bam.gr.filtertrack.bytype))){
		
		prefix <- str_c(
			bam.gr.filtertrack.bytype$call_class[j],
			bam.gr.filtertrack.bytype$call_type[j],
			bam.gr.filtertrack.bytype$SBSindel_call_type[j],
			sep = "."
		)
		
		bam.gr.filtertrack.bytype %>%
			pluck(i,j) %>%
			write_tsv(
				file = str_c(
					str_c(prefix,"/",output_basename),
					prefix,
					i,
					"tsv",
					sep="."
				)
			)
	}
}

#genome
for(i in c("genome.reftnc_both_strands","genome.reftnc_duplex_pyr","genome_chromgroup.reftnc_both_strands","genome_chromgroup.reftnc_duplex_pyr")){
	i %>%
		get %>%
		write_tsv(str_c(output_basename,i,"tsv",sep="."))
}

######################
### Output final calls and germline variant calls
######################

cat("## Outputting final calls and detected germline variants...")

#Output final calls to tsv and vcf, separately for each combination of call_class, call_type, SBSindel_call_type
finalCalls.bytype %>%
	pwalk(
		function(...){
			x <- list(...)
			
			prefix <- str_c(
				x$call_class,
				x$call_type,
				x$SBSindel_call_type,
				sep="."
			)
			
			#Format list columns to comma-delimited
			x$finalCalls <- x$finalCalls %>%
				mutate(
					across(
						where(is.list),
						function(x){map_chr(x, function(v){str_c(v, collapse = ",")})}
					)
				)
			
			#tsv
			x$finalCalls %>%
				write_tsv(
					str_c(
						str_c(prefix,"/",output_basename),
						prefix,
						"finalCalls",
						"tsv",
						sep="."
					)
				)
			
			#vcf
			x$finalCalls_for_vcf %>%
				write_vcf_from_calls(
					BSgenome_name = yaml.config$BSgenome$BSgenome_name,
					out_vcf = str_c(
						str_c(prefix,"/",output_basename),
						prefix,
						"finalCalls",
						"vcf",
						sep="."
					)
				)
			
		}
	)

##TODO - OUTPUT UNIQUE CALLS to TSV and VCF

#Output germline variant calls
 #Format list columns to comma-delimited
germlineVariantCalls.out <- germlineVariantCalls %>%
	mutate(
		across(
			where(is.list),
			function(x){map_chr(x, function(v){str_c(v, collapse = ",")})}
		)
	)

 #tsv
germlineVariantCalls.out %>%
	write_tsv(
		str_c(output_basename,"germlineVariantCalls","tsv",sep=".")
	)

 #vcf
germlineVariantCalls.out %>% 
	normalize_indels_for_vcf(
		BSgenome_name = yaml.config$BSgenome$BSgenome_name
	) %>%
	write_vcf_from_calls(
		BSgenome_name = yaml.config$BSgenome$BSgenome_name,
		out_vcf = str_c(output_basename,"germlineVariantCalls","vcf",sep=".")
	)

rm(germlineVariantCalls.out)

cat("DONE\n")

####
Output Sensitivity


#Call burdens
 -> Use sensitivity from matching filtergroup and from use_chromgroup, or from any analysis if use_chromgroup is null (since all identically set SBS and indel senstivity to 1)
 - calculate for uncorrected and tnc_corrected burdens
	
#Call spectra
 #plots

#Estimated mutation error rate