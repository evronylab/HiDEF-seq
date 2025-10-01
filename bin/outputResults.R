#!/usr/bin/env -S Rscript --vanilla

#outputResults.R:
# Output results. Also calculates sensitivity-corrected burdens and estimated SBS mutation error rate.

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

if(is.null(opt$config) | is.null(opt$sample_id_toanalyze) | is.null(opt$files) | is.null(opt$output_basename) ){
	stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))
sample_id_toanalyze <- opt$sample_id_toanalyze
calculateBurdensFiles <- opt$files %>% str_split_1(",") %>% str_trim
output_basename <- opt$output_basename

#Load the BSgenome reference
suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name,character.only=TRUE,lib.loc=yaml.config$cache_dir))

#Load miscellaneous configuration parameters
 #analysis_id
analysis_id <- yaml.config$analysis_id

 #individual_id of this sample_id
individual_id <- yaml.config$samples %>%
  bind_rows %>%
  filter(sample_id == sample_id_toanalyze) %>%
  pull(individual_id)

 #call types
call_types <- yaml.config$call_types %>%
	enframe(name=NULL) %>%
	unnest_wider(value) %>%
	unnest_longer(SBSindel_call_types) %>%
	unnest_wider(SBSindel_call_types) %>%
	select(-starts_with("MDB"))

chromgroups <-  yaml.config$chromgroups %>%
	enframe(name=NULL) %>%
	unnest_wider(value) %>%
	pull(chromgroup)

filtergroups <-  yaml.config$filtergroups %>%
	enframe(name=NULL) %>%
	unnest_wider(value) %>%
	pull(filtergroup)

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
		setdiff(c("seqnames","start","ref_plus_strand","alt_plus_strand","qual"))
	
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
cat("## Loading data from calculateBurdens files...\n")

#Create lists for data loading
run_metadata <- list()
molecule_stats.by_run_id <- list()
molecule_stats.by_analysis_id <- list()
region_genome_filter_stats <- list()
finalCalls <- list()
finalCalls.bytype = list()
germlineVariantCalls <- list()
finalCalls.reftnc_spectra <- list()
finalCalls.burdens <- list()
bam.gr.filtertrack.bytype.coverage_tnc <- list()
genome.reftnc <- list()
genome_chromgroup.reftnc <- list()
sensitivity <- list()
estimatedSBSMutationErrorProbability <- list()

#Loop over all calculateBurdens
for(i in seq_along(calculateBurdensFiles)){
	
	#Load calculateBurdensFile
	calculateBurdensFile <- qs_read(calculateBurdensFiles[i])
	
	#Annotations to add to tables
	sample_annotations <- list(
		analysis_id = analysis_id %>% factor,
		individual_id = individual_id %>% factor,
		sample_id = sample_id_toanalyze %>% factor
	)
	
	chromgroup_filtergroup_annotations <- list(
		chromgroup = calculateBurdensFile$chromgroup %>%
			factor(levels = chromgroups),
		filtergroup = calculateBurdensFile$filtergroup %>%
			factor(levels = filtergroups)
	)
	
	cat(" > chromgroup:", calculateBurdensFile$chromgroup, "/", "filtergroup:", calculateBurdensFile$filtergroup, "\n")
	
	#Run metadata
	run_metadata[[i]] <- calculateBurdensFile %>%
		pluck("run_metadata") %>%
		mutate(!!!sample_annotations, .before = 1)
	
	#Filtering stats
	molecule_stats.by_run_id[[i]] <- calculateBurdensFile %>%
		pluck("molecule_stats.by_run_id") %>%
		mutate(
			!!!sample_annotations,
			chromgroup = chromgroup %>% factor(levels = c("all_chroms","all_chromgroups",chromgroups)),
			filtergroup = filtergroup %>% factor(levels = filtergroups),
			.before = 1
		) %>%
		relocate(run_id, .after = filtergroup)
	
	molecule_stats.by_analysis_id[[i]] <- calculateBurdensFile %>%
		pluck("molecule_stats.by_analysis_id") %>%
		mutate(
			!!!sample_annotations,
			chromgroup = chromgroup %>% factor(levels = c("all_chroms","all_chromgroups",chromgroups)),
			filtergroup = filtergroup %>% factor(levels = filtergroups),
			.before = 1
		)
	
	region_genome_filter_stats[[i]] <- calculateBurdensFile %>%
		pluck("region_genome_filter_stats") %>%
		mutate(
			!!!sample_annotations,
			!!!chromgroup_filtergroup_annotations,
			.before = 1
		)
	
	#Final calls
	finalCalls[[i]] <- calculateBurdensFile %>%
		pluck("finalCalls") %>%
		mutate(
			!!!sample_annotations,
			!!!chromgroup_filtergroup_annotations,
			.before = 1
		)
	
	#finalCalls for tsv and vcf output
	finalCalls.bytype[[i]] <- calculateBurdensFile %>%
		pluck("finalCalls.bytype") %>%
		mutate(
			!!!sample_annotations,
			chromgroup = !!calculateBurdensFile$chromgroup %>% factor(levels = c(chromgroups)),
			filtergroup = filtergroup %>% factor(levels = filtergroups),
			.before = 1
		) %>%
		relocate(filtergroup, call_class, .after = chromgroup) %>%
		select(-analyzein_chromgroups)
	
	#Germline variant calls
	germlineVariantCalls[[i]] <- calculateBurdensFile %>%
		pluck("germlineVariantCalls") %>%
		mutate(
			!!!sample_annotations,
			!!!chromgroup_filtergroup_annotations,
			.before = 1
		)
	
	#Trinucleotide context counts and fractions of finalCalls
	finalCalls.reftnc_spectra[[i]] <- calculateBurdensFile %>%
		pluck("finalCalls.reftnc_spectra") %>%
		mutate(
			!!!sample_annotations,
			chromgroup = !!calculateBurdensFile$chromgroup %>% factor(levels = c(chromgroups)),
			filtergroup = filtergroup %>% factor(levels = filtergroups),
			.before = 1
		) %>%
		relocate(filtergroup, call_class, .after = chromgroup) %>%
		select(-analyzein_chromgroups)
	
	#Burdens
	finalCalls.burdens[[i]] <- calculateBurdensFile %>%
		pluck("finalCalls.burdens") %>%
		mutate(
			!!!sample_annotations,
			chromgroup = !!calculateBurdensFile$chromgroup %>% factor(levels = c(chromgroups)),
			filtergroup = filtergroup %>% factor(levels = filtergroups),
			.before = 1
		) %>%
		relocate(filtergroup, call_class, .after = chromgroup) %>%
		select(-analyzein_chromgroups)

	#Genome coverage and trinucleotide counts, fractions, and ratio to genome
	bam.gr.filtertrack.bytype.coverage_tnc[[i]] <- calculateBurdensFile %>%
		pluck("bam.gr.filtertrack.bytype.coverage_tnc") %>%
		mutate(
			!!!sample_annotations,
			chromgroup = !!calculateBurdensFile$chromgroup %>% factor(levels = c(chromgroups)),
			filtergroup = filtergroup %>% factor(levels = filtergroups),
			.before = 1
		) %>%
		relocate(filtergroup, call_class, .after = chromgroup) %>%
		select(-analyzein_chromgroups)
	
	#Whole genome trinucleotide counts and fractions
	genome.reftnc[[i]] <- calculateBurdensFile %>%
		pluck("genome.reftnc") %>%
		enframe %>%
		pivot_wider(names_from = name, values_from = value)
	
	#Genome chromgroup trinucleotide counts and fractions
		genome_chromgroup.reftnc[[i]] <- calculateBurdensFile %>%
			pluck("genome_chromgroup.reftnc") %>%
			enframe %>%
			pivot_wider(names_from = name, values_from = value) %>%
				mutate(
					chromgroup = !!calculateBurdensFile$chromgroup %>% factor(levels = chromgroups),
					.before = 1
				)
	
	#Sensitivity
	if(! calculateBurdensFile %>% pluck("sensitivity") %>% is.null){
		sensitivity[[i]] <- calculateBurdensFile %>%
			pluck("sensitivity") %>%
			mutate(
				!!!sample_annotations,
				chromgroup = !!calculateBurdensFile$chromgroup %>% factor(levels = c(chromgroups)),
				filtergroup = filtergroup %>% factor(levels = filtergroups),
				.before = 1
			) %>%
			relocate(filtergroup, call_class, .after = chromgroup) %>%
			select(-analyzein_chromgroups)
	}
	
	#estimated SBS mutation error probability
	if(! calculateBurdensFile %>% pluck("estimatedSBSMutationErrorProbability") %>% is.null){
		estimatedSBSMutationErrorProbability[[i]] <- calculateBurdensFile %>%
			pluck("estimatedSBSMutationErrorProbability") %>%
			enframe %>%
			pivot_wider(names_from = name, values_from = value) %>%
			mutate(
				!!!sample_annotations,
				!!!chromgroup_filtergroup_annotations,
				.before = 1
			)
	}
	
	#Remove temporary objects
	rm(calculateBurdensFile)
	invisible(gc())
}

#Combine data
run_metadata <- run_metadata %>%
	bind_rows %>%
	distinct

molecule_stats.by_run_id <- molecule_stats.by_run_id %>%
	bind_rows %>%
	distinct

molecule_stats.by_analysis_id <- molecule_stats.by_analysis_id %>%
	bind_rows %>%
	distinct

region_genome_filter_stats <- region_genome_filter_stats %>% bind_rows

finalCalls <- finalCalls %>% bind_rows

finalCalls.bytype <- finalCalls.bytype %>% bind_rows
	
germlineVariantCalls <- germlineVariantCalls %>% bind_rows
	
finalCalls.reftnc_spectra <- finalCalls.reftnc_spectra %>% bind_rows
	
finalCalls.burdens <- finalCalls.burdens %>% bind_rows

bam.gr.filtertrack.bytype.coverage_tnc <- bam.gr.filtertrack.bytype.coverage_tnc %>% bind_rows

genome.reftnc <- genome.reftnc %>%
	bind_rows %>%
	distinct

genome_chromgroup.reftnc <- genome_chromgroup.reftnc %>%
	bind_rows %>%
	distinct

sensitivity <- sensitivity %>% bind_rows
	
estimatedSBSMutationErrorProbability <- estimatedSBSMutationErrorProbability %>% bind_rows

invisible(gc())

cat(" DONE\n")

######################
### Output configuration parameters
######################
cat("## Outputting configuration parameters...")

opt$config %>% file.copy(str_c(output_basename,".yaml_config.tsv"))

run_metadata %>% write_tsv(str_c(output_basename,".run_metadata.tsv"))

cat("DONE\n")

######################
### Output filtering statistics
######################
cat("## Outputting filtering statistics...")

for(i in chromgroups){
	#Create output folders once
	dir_create(i)
	
	for(j in filtergroups){
		
		#molecule_stats.by_run_id
		 #Outupt all_chroms and all_chromgroups stats to each output file
		molecule_stats.by_run_id %>%
			filter(
				chromgroup %in% c("all_chroms", "all_chromgroups") |
					(chromgroup == i & filtergroup == j)
				) %>%
			write_tsv(str_c(i,"/",output_basename,i,j,"molecule_stats.by_run_id.tsv",sep="."))
		
		#molecule_stats.by_analysis_id
		molecule_stats.by_analysis_id %>%
			filter(
				chromgroup %in% c("all_chroms", "all_chromgroups") |
					(chromgroup == i & filtergroup == j)
			) %>%
			write_tsv(str_c(i,"/",output_basename,i,j,"molecule_stats.by_analysis_id.tsv",sep="."))
		
		#region_genome_filter_stats
		region_genome_filter_stats %>%
			filter(chromgroup == i & filtergroup == j) %>%
			write_tsv(str_c(i,"/",output_basename,i,j,"region_genome_filter_stats.tsv",sep="."))

	}
}

cat("DONE\n")

######################
### Output final calls and germline variant calls
######################

cat("## Outputting final calls and germline variant calls...")

#Output final calls to tsv and vcf, separately for each combination of call_class, call_type, SBSindel_call_type

finalCalls.bytype %>%
	##**ADD here to rename strand in all the finalCalls tibbles (for tsv and for vcf, all and unique) to aligned_synthesized_strand, to help users understand more easily what the 'strand' column is
	pwalk(
		function(...){
			x <- list(...)
			
			prefix <- str_c(
				x$call_class,
				x$call_type,
				x$SBSindel_call_type,
				sep="."
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
#Format list columns to comma-delimited *** UPDATE FROM PRIOR calcburdens CODE used for finalCalls.
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


#Output trinucleotide counts
genome.reftnc_pyr
genome.reftnc_both_strands

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



##Call spectra
#tables

#plots

##burdens

####
Output Sensitivity
- Assign the use_chromgroup sensitivity to all chromgroups, or if use_chromgroup is null, every chromgroup will already have the default sensitivity tibble assigned by the calculateBurdens script

#Call burdens
-> Use sensitivity from matching filtergroup and from use_chromgroup, or from any analysis if use_chromgroup is null (since all identically set SBS and indel senstivity to 1)
- calculate for uncorrected and tnc_corrected burdens

#Estimated SBS mutation error rate

#Output RDS with configuration parameters and stats
qs_save(
	lst(
		yaml.config,
		run_metadata,
		individual_id,
		sample_id = sample_id_toanalyze,
		molecule_stats.by_run_id,
		molecule_stats.by_analysis_id,
		region_genome_filter_stats,
		finalCalls,
		finalCalls.bytype,
		germlineVariantCalls,
		finalCalls.reftnc_spectra,
		finalCalls.burdens,
		bam.gr.filtertrack.bytype.coverage_tnc,
		genome.reftnc,
		genome_chromgroup.reftnc,
		sensitivity,
		estimatedSBSMutationErrorProbability
	)
	str_c(output_basename,".qs2")
)