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
suppressPackageStartupMessages(library(survival))
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

#Create output directories
filterStats_dir <- "/filterStats/"
finalCalls_dir <- "/finalCalls/"
germlineVariantCalls_dir <- "/germlineVariantCalls/"
finalCalls.spectra_dir <- "/finalCalls.spectra/"
interrogatedBases.spectra_dir <- "/interrogatedBases.spectra/"
genome.spectra_dir <- "/genome.spectra/"
sensitivity_dir <- "/sensitivity/"
finalCalls.burdens_dir <- "/finalCalls.burdens/"
estimatedSBSMutationErrorProbability_dir <- "/estimatedSBSMutationErrorProbability/"

for(i in chromgroups){
	for(j in c(filterStats_dir, finalCalls_dir, germlineVariantCalls_dir, finalCalls.spectra_dir, interrogatedBases.spectra_dir, genome.spectra_dir, sensitivity_dir, finalCalls.burdens_dir, estimatedSBSMutationErrorProbability_dir))
	fs::dir_create(str_c(i,"/",j))
}

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
	
	#Repair names in 'calls' so they are compatible with DataFrame objects required by VariantAnnotation (i.e. change '-' to '.' in column names)
	names(calls) <- names(calls) %>% vctrs::vec_as_names(repair = "universal", quiet = TRUE)
	
	#Fixed fields
	CHROM <- calls$seqnames
	POS <- calls$start_refspace
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
	
	# INFO: everything except the basic VCF columns
	info_cols <- names(calls) %>%
		setdiff(c("seqnames","start_refspace","ref_plus_strand","alt_plus_strand"))
	
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
bam.gr.filtertrack.bytype.coverage_tnc <- list()
genome.reftnc <- list()
genome_chromgroup.reftnc <- list()
sensitivity <- list()
finalCalls.burdens <- list()
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
			chromgroup = !!calculateBurdensFile$chromgroup %>% factor(levels = chromgroups),
			filtergroup = filtergroup %>% factor(levels = filtergroups),
			.before = 1
		) %>%
		relocate(filtergroup, call_class, .after = chromgroup)
	
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
			chromgroup = !!calculateBurdensFile$chromgroup %>% factor(levels = chromgroups),
			filtergroup = filtergroup %>% factor(levels = filtergroups),
			.before = 1
		) %>%
		relocate(filtergroup, call_class, .after = chromgroup)
	
	#Genome coverage and trinucleotide counts, fractions, and ratio to genome
	bam.gr.filtertrack.bytype.coverage_tnc[[i]] <- calculateBurdensFile %>%
		pluck("bam.gr.filtertrack.bytype.coverage_tnc") %>%
		mutate(
			!!!sample_annotations,
			chromgroup = !!calculateBurdensFile$chromgroup %>% factor(levels = c(chromgroups)),
			filtergroup = filtergroup %>% factor(levels = filtergroups),
			.before = 1
		) %>%
		relocate(filtergroup, call_class, .after = chromgroup)
	
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
				chromgroup = !!calculateBurdensFile$chromgroup %>% factor(levels = chromgroups),
				filtergroup = filtergroup %>% factor(levels = filtergroups),
				.before = 1
			) %>%
			relocate(filtergroup, call_class, .after = chromgroup)
	}

	#Burdens
	finalCalls.burdens[[i]] <- calculateBurdensFile %>%
		pluck("finalCalls.burdens") %>%
		mutate(
			!!!sample_annotations,
			chromgroup = !!calculateBurdensFile$chromgroup %>% factor(levels = chromgroups),
			filtergroup = filtergroup %>% factor(levels = filtergroups),
			.before = 1
		) %>%
		relocate(filtergroup, call_class, .after = chromgroup)
		
	#estimated SBS mutation error probability
	if(! calculateBurdensFile %>% pluck("estimatedSBSMutationErrorProbability") %>% is.null){
		estimatedSBSMutationErrorProbability[[i]] <- calculateBurdensFile %>%
			pluck("estimatedSBSMutationErrorProbability") %>%
			enframe %>%
			pivot_wider(names_from = name, values_from = value) %>%
			unnest(total) %>%
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

finalCalls <- finalCalls %>%
	bind_rows %>%
	#Rename columns for greater clarity in final output
	rename(
		refstrand = strand,
		start_refspace = start,
		end_refspace = end
	)

finalCalls.bytype <- finalCalls.bytype %>%
	bind_rows %>%
	#Rename columns for greater clarity in final output
	mutate(
		across(
			c(finalCalls_for_tsv, finalCalls_unique_for_tsv, finalCalls_for_vcf, finalCalls_unique_for_vcf),
			function(x){
				x %>%
					map(function(y){
						if(!is.null(y)){
							y %>%
								rename(
									refstrand = any_of("strand"),
									start_refspace = start,
									end_refspace = any_of("end")
								)
						}else{
							y
						}
					})
			}
		)
	)

germlineVariantCalls <- germlineVariantCalls %>%
	bind_rows %>%
	#Rename columns for greater clarity in final output
	rename(
		refstrand = strand,
		start_refspace = start,
		end_refspace = end
	)
	
finalCalls.reftnc_spectra <- finalCalls.reftnc_spectra %>% bind_rows
	
bam.gr.filtertrack.bytype.coverage_tnc <- bam.gr.filtertrack.bytype.coverage_tnc %>% bind_rows

genome.reftnc <- genome.reftnc %>%
	bind_rows %>%
	distinct

genome_chromgroup.reftnc <- genome_chromgroup.reftnc %>%
	bind_rows %>%
	distinct

sensitivity <- sensitivity %>% bind_rows

finalCalls.burdens <- finalCalls.burdens %>% bind_rows

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
	for(j in filtergroups){
		
		output_basename_full <- str_c(str_c(i,filterStats_dir,output_basename), i, j, sep=".")
		
		#molecule_stats.by_run_id
		 #Outupt all_chroms and all_chromgroups stats to each output file
		molecule_stats.by_run_id %>%
			filter(
				chromgroup %in% c("all_chroms", "all_chromgroups") |
					(chromgroup == i & filtergroup == j)
				) %>%
			write_tsv(str_c(output_basename_full,".molecule_stats.by_run_id.tsv"))
		
		#molecule_stats.by_analysis_id
		molecule_stats.by_analysis_id %>%
			filter(
				chromgroup %in% c("all_chroms", "all_chromgroups") |
					(chromgroup == i & filtergroup == j)
			) %>%
			write_tsv(str_c(output_basename_full,".molecule_stats.by_analysis_id.tsv"))
		
		#region_genome_filter_stats
		region_genome_filter_stats %>%
			filter(chromgroup == i & filtergroup == j) %>%
			write_tsv(str_c(output_basename_full,".region_genome_filter_stats.tsv"))

	}
}

rm(output_basename_full)

cat("DONE\n")

######################
### Output final calls
######################

cat("## Outputting final calls...")

#Output final calls to tsv and vcf, separately for each combination of chromgroup, filtergroup, call_class, call_type, SBSindel_call_type
finalCalls.bytype %>%
	pwalk(
		function(...){
			x <- list(...)
			
			output_basename_full <- str_c(
				str_c(x$chromgroup,finalCalls_dir,output_basename),
				x$chromgroup,
				x$filtergroup,
				x$call_class,
				x$call_type,
				x$SBSindel_call_type,
				sep="."
			)
			
			metadata <- x %>%
				keep(
					names(.) %in%
						c(
							"analysis_id", "individual_id", "sample_id",
							"chromgroup", "filtergroup",
							"call_class", "call_type", "SBSindel_call_type"
						)
				) %>%
				as_tibble
			
			#all calls tsv
			metadata %>%
				mutate(finalCalls_for_tsv = list(x$finalCalls_for_tsv)) %>%
				unnest(finalCalls_for_tsv) %>%
				write_tsv(str_c(output_basename_full,".finalCalls.tsv"))
			
			#unique calls tsv
			if(!is.null(x$finalCalls_unique_for_tsv)){
				metadata %>%
					mutate(finalCalls_unique_for_tsv = list(x$finalCalls_unique_for_tsv)) %>%
					unnest(finalCalls_unique_for_tsv) %>%
					write_tsv(str_c(output_basename_full,".finalCalls_unique.tsv"))
			}
			
			#all calls vcf
			metadata %>%
				mutate(finalCalls_for_vcf = list(x$finalCalls_for_vcf)) %>%
				unnest(finalCalls_for_vcf) %>%
				write_vcf_from_calls(
					BSgenome_name = yaml.config$BSgenome$BSgenome_name,
					out_vcf = str_c(output_basename_full,".finalCalls.vcf")
				)
			
			#unique calls vcf
			if(!is.null(x$finalCalls_unique_for_vcf)){
				metadata %>%
					mutate(finalCalls_unique_for_vcf = list(x$finalCalls_unique_for_vcf)) %>%
					unnest(finalCalls_unique_for_vcf) %>%
					write_vcf_from_calls(
						BSgenome_name = yaml.config$BSgenome$BSgenome_name,
						out_vcf = str_c(output_basename_full,".finalCalls_unique.vcf")
					)
			}
		}
	)

#Remove 'for_vcf' tables since not needed anymore.
finalCalls.bytype <- finalCalls.bytype %>%
	select(-c(finalCalls_for_vcf, finalCalls_unique_for_vcf))

invisible(gc())

cat("DONE\n")

######################
### Output germline variant calls
######################

cat("## Outputting germline variant calls...")

#Define columns that are identical between strands to either '_keep' or '_discard' in the subsequent pivot_wider.
strand_identical_cols_keep <- c(
	"analysis_id", "individual_id", "sample_id", "chromgroup", "filtergroup",
	"analysis_chunk", "run_id", "zm",
	"call_class", "call_type",
	"seqnames", "start_refspace", "end_refspace", "ref_plus_strand", "alt_plus_strand",
	"reftnc_plus_strand", "alttnc_plus_strand", "reftnc_pyr", "alttnc_pyr",
	"indel_width"
)

strand_identical_cols_discard <- c(
	"refstrand",
	"call_class.opposite_strand", "call_type.opposite_strand",
	"alt_plus_strand.opposite_strand",
	"deletion.bothstrands.startendmatch", "MDB_score"
)

#Output germline variant calls to tsv and vcf, separately for each combination of chromgroup x filtergroup
for(i in chromgroups){
	for(j in filtergroups){

		#Output path
		output_basename_full <- str_c(str_c(i,germlineVariantCalls_dir,output_basename),i,j,sep=".")
		
		#Define germline filter columns to keep in the output
		region_read_filters_cols_keep <- yaml.config$region_filters %>%
			map("read_filters") %>%
			flatten %>%
			enframe(name=NULL) %>%
			unnest_wider(value) %>%
			mutate(region_filter_threshold_file = str_c(yaml.config$cache_dir,"/",basename(region_filter_file),".bin",binsize,".",threshold,".bw")) %>%
			filter(
				applyto_chromgroups == "all" | (applyto_chromgroups %>% str_split(",") %>% map(str_trim) %>% map_lgl(~ !!i %in% .x)),
				applyto_filtergroups == "all" | (applyto_filtergroups %>% str_split(",") %>% map(str_trim) %>% map_lgl(~ !!j %in% .x))
			) %>%
			filter(is_germline_filter==TRUE) %>%
			pull(region_filter_threshold_file) %>%
			basename %>%
			str_c("region_read_filter_",.,".passfilter")
		
		region_genome_filters_cols_keep <- yaml.config$region_filters %>%
			map("genome_filters") %>%
			flatten %>%
			enframe(name=NULL) %>%
			unnest_wider(value) %>%
			mutate(region_filter_threshold_file = str_c(yaml.config$cache_dir,"/",basename(region_filter_file),".bin",binsize,".",threshold,".bw")) %>%
			filter(
				applyto_chromgroups == "all" | (applyto_chromgroups %>% str_split(",") %>% map(str_trim) %>% map_lgl(~ !!i %in% .x)),
				applyto_filtergroups == "all" | (applyto_filtergroups %>% str_split(",") %>% map(str_trim) %>% map_lgl(~ !!j %in% .x))
			) %>%
			filter(is_germline_filter==TRUE) %>%
			pull(region_filter_threshold_file) %>%
			basename %>%
			str_c("region_genome_filter_",.,".passfilter")
		
		germline_filter_cols_keep <- c( #not including germline_vcf.passfilter, since FALSE for all calls, and not including germline_vcf_indel_region_filter, max_germlineBAM_VariantReads.passfilter, max_germlineBAM_VAF.passfilter, since these are not relevant annotation to output for germline calls themselves (only relevant for somatic calls that are filtered based on germline calls).
			region_read_filters_cols_keep,
			region_genome_filters_cols_keep,
			"germline_vcf_types_detected", "germline_vcf_files_detected"
		)
		
		#Create and format output tibble
		germlineVariantCalls.out <- germlineVariantCalls %>%
			#Remove all passfilter columns (since all true) except germline filter columns
			select(-(contains("passfilter") & !all_of(germline_filter_cols_keep))) %>%
			
			#Select current chromgroup/filtergroup
			filter(chromgroup == i, filtergroup == j) %>%
			
			#Reformat list columns to be comma-delimited
			mutate(
				across(
					where(is.list),
					function(x){x %>% map_chr(function(v){str_c(v, collapse = ",")})}
				)
			) %>%
			
			#Change strand levels so that new column names after pivot_wider have suffixes ref_strand_(plus/minus)_read instead of +/-.
			mutate(
				refstrand = refstrand %>%
					fct_recode(
						refstrand_plus_read  = "+",
						refstrand_minus_read = "-"
					)
			) %>%
			
			#Collapse to one row per mutation
			pivot_wider(
				id_cols = all_of(c(strand_identical_cols_keep, germline_filter_cols_keep)),
				names_from = refstrand,
				values_from = -all_of(c(strand_identical_cols_keep, germline_filter_cols_keep, strand_identical_cols_discard)),
				names_glue = "{.value}_{refstrand}",
				names_expand = TRUE
			)
		
		#tsv
		germlineVariantCalls.out %>%
			write_tsv(str_c(output_basename_full,".germlineVariantCalls.tsv"))
		
		#vcf
		germlineVariantCalls.out %>%
			#Rename back columns for normalize_indels_for_vcf
			rename(
				start = start_refspace,
				end = end_refspace
			) %>%
			normalize_indels_for_vcf(
				BSgenome_name = yaml.config$BSgenome$BSgenome_name
			) %>%
			#Rename back columns for greater clarity in final output
			rename(start_refspace = start) %>%
			write_vcf_from_calls(
				BSgenome_name = yaml.config$BSgenome$BSgenome_name,
				out_vcf = str_c(output_basename_full,".germlineVariantCalls.vcf")
			)
		
		rm(germlineVariantCalls.out)
		invisible(gc())
		
	}
}

cat("DONE\n")

######################
### Output spectra of calls, interrogated bases, and the genome
######################
cat("## Outputting spectra of calls, interrogated bases, and the genome...")

#Helper function: write TSV for a given tibble col_name
write_col <- function(df.input, col_name.input, metadata.input, output_basename_full.input) {
	metadata.input %>%
		mutate(!!sym(col_name.input) := list(df.input)) %>%
		unnest(all_of(col_name.input)) %>%
		write_tsv(str_c(output_basename_full.input, ".", col_name.input, ".tsv"))
}

#Helper function: plots spectrum for a given sigfit col_name
plot_col <- function(df.input, col_name.input, output_basename_full.input){
	#Output name
	output_name <- str_c(output_basename_full.input, ".", col_name.input)
	
	#Count number of calls
	sum_counts <- sum(df.input)
	
	#Normalize to number of calls
	if(sum_counts > 0){
		df.input <- df.input / sum(df.input)
	}
	
	df.input %>%
		plot_spectrum(
			pdf_path = str_c(output_name, ".pdf"),
			name = str_c(output_name %>% basename, " (", sum_counts, "calls)")
			)
}

#Output spectra of finalCalls for each combination of chromgroup, filtergroup, call_class, call_type, SBSindel_call_type
#Note: this outputs SBS/mismatch-ss trinucleotide distributions that were retained in finalCalls.reftnc_spectra from upstream calculateBurdens for every chromgroup/filtergroup, even if not configured to call SBS/mismatch-ss in a chromgroup/filtergroup. This aids assessment of mismatch patterns for every chromgroup/filtergroup for the purpose of evaluating SBS mutation error probability (which is also later calculated directly).
finalCalls.reftnc_spectra %>%
	pwalk(
		function(...){
			x <- list(...)
			
			output_basename_full <- str_c(
				str_c(x$chromgroup,finalCalls.spectra_dir,output_basename),
				x$chromgroup,
				x$filtergroup,
				x$call_class,
				x$call_type,
				x$SBSindel_call_type,
				sep="."
			)
			
			metadata <- x %>%
				keep(
					names(.) %in%
						c(
							"analysis_id", "individual_id", "sample_id",
							"chromgroup", "filtergroup",
							"call_class", "call_type", "SBSindel_call_type"
						)
				) %>%
				as_tibble
			
			sbs_tables_to_output <- c(
				"finalCalls.reftnc_pyr",
				"finalCalls_unique.reftnc_pyr",
				"finalCalls.reftnc_template_strand",
				"finalCalls.reftnc_pyr_spectrum",
				"finalCalls_unique.reftnc_pyr_spectrum",
				"finalCalls.reftnc_template_strand_spectrum"
			)
			
			indel_tables_to_output <- c(
				"finalCalls.refindel_spectrum.sigfit",
				"finalCalls_unique.refindel_spectrum.sigfit"
			)
			
			plots_to_output <- c(
				"finalCalls.reftnc_pyr_spectrum.sigfit",
				"finalCalls_unique.reftnc_pyr_spectrum.sigfit",
				"finalCalls.reftnc_template_strand_spectrum.sigfit",
				"finalCalls.refindel_spectrum.sigfit",
				"finalCalls_unique.refindel_spectrum.sigfit"
			)
			
			sbs_tables_to_output %>%
				walk(function(y){
					if(!is.null(x[[y]])){
						write_col(
							df.input = x[[y]],
							col_name.input = y, metadata.input = metadata, output_basename_full.input = output_basename_full
						)
					}
				})
			
			indel_tables_to_output %>%
				walk(function(y){
					if(!is.null(x[[y]])){
						write_col(
							df.input = x[[y]] %>%
								as_tibble %>%
								pivot_longer(cols=everything(), names_to = "channel", values_to = "count") %>%
								mutate(fraction = count / sum(count)),
							col_name.input = y, metadata.input = metadata, output_basename_full.input = output_basename_full
						)
					}
				})
			
			plots_to_output %>%
				walk(function(y){
					if(!is.null(x[[y]])){
						plot_col(
							df.input = x[[y]],
							col_name.input = y, output_basename_full.input = output_basename_full
						)
					}
				})
			
		}
	)

#Output spectra of interrogated bases for each combination of chromgroup, filtergroup, call_class, call_type, SBSindel_call_type
bam.gr.filtertrack.bytype.coverage_tnc %>%
	pwalk(
		function(...){
			x <- list(...)
			
			output_basename_full <- str_c(
				str_c(x$chromgroup,interrogatedBases.spectra_dir,output_basename),
				x$chromgroup,
				x$filtergroup,
				x$call_class,
				x$call_type,
				x$SBSindel_call_type,
				sep="."
			)
			
			metadata <- x %>%
				keep(
					names(.) %in%
						c(
							"analysis_id", "individual_id", "sample_id",
							"chromgroup", "filtergroup",
							"call_class", "call_type", "SBSindel_call_type"
						)
				) %>%
				as_tibble
			
			sbs_tables_to_output <- c(
				"bam.gr.filtertrack.reftnc_pyr",
				"bam.gr.filtertrack.reftnc_both_strands"
			)
			
			sbs_tables_to_output %>%
				walk(function(y){
					write_col(
						df.input = x[[y]],
						col_name.input = y, metadata.input = metadata, output_basename_full.input = output_basename_full
					)
				})
			
		}
	)

#Output spectra of genome bases
for(i in chromgroups){
	genome.reftnc %>%
		select(reftnc_pyr) %>%
		unnest(reftnc_pyr) %>%
		write_tsv(str_c(i,genome.spectra_dir,"genome.reftnc_pyr.tsv"))
	
	genome_chromgroup.reftnc %>%
		filter(chromgroup == i) %>%
		select(chromgroup, reftnc_pyr) %>%
		unnest(reftnc_pyr) %>%
		write_tsv(str_c(i,genome.spectra_dir,i,".genome_chromgroup.reftnc_pyr.tsv"))
		
	genome.reftnc %>%
		select(reftnc_both_strands) %>%
		unnest(reftnc_both_strands) %>%
		write_tsv(str_c(i,genome.spectra_dir,"genome.reftnc_both_strands.tsv"))
	
	genome_chromgroup.reftnc %>%
		filter(chromgroup == i) %>%
		select(chromgroup, reftnc_both_strands) %>%
		unnest(reftnc_both_strands) %>%
		write_tsv(str_c(i,genome.spectra_dir,i,".genome_chromgroup.reftnc_both_strands.tsv"))
}

cat("DONE\n")

######################
### Output sensitivity
######################
cat("## Outputting sensitivity...")

#For each call_class, call_type, SBSindel_call_type, for rows with sensitivity_source == "other_chromgroup", copy sensitivity column value from row with the same filtergroup that has sensitivity_source = "calculated", and change the row's sensitivity source to "calculated_other_chromgroup". If there is no such "calculated" row, then set the row's sensitivity to 1 and the sensitivity source to "default_other_chromgroup".
sensitivity <- sensitivity %>%
	group_by(call_class, call_type, SBSindel_call_type, filtergroup) %>%
	mutate(
		sensitivity.new = if_else(
			sensitivity_source == "other_chromgroup",
			first(sensitivity[sensitivity_source == "calculated"], default = NA),
			sensitivity
		),
		
		sensitivity_source = case_when(
			sensitivity_source == "other_chromgroup" & is.na(sensitivity.new) ~ "default_other_chromgroup",
			sensitivity_source == "other_chromgroup" & !is.na(sensitivity.new) ~ "calculated_other_chromgroup",
			.default = sensitivity_source
		),
		
		sensitivity = sensitivity.new %>% replace_na(1)
	) %>%
	ungroup %>% 
	select(-sensitivity.new)

#Output as tsv
sensitivity %>%
	nest_by(chromgroup, filtergroup, .keep = TRUE) %>%
	pwalk(
		function(...){
			x <- list(...)
			
			output_basename_full <- str_c(
				str_c(x$chromgroup, sensitivity_dir, output_basename),
				x$chromgroup,
				x$filtergroup,
				x$call_class,
				x$call_type,
				x$SBSindel_call_type,
				sep="."
			)
			
			x$data %>% write_tsv(str_c(output_basename_full, ".sensitivity.tsv"))
			
		}
	)

cat("DONE\n")

######################
### Output call burdens
######################
cat("## Outputting call burdens...")

#Add burdens corrected for sensitivity, and retain sensitivity and sensitivity_source columns in finalCalls.burdens tibble
finalCalls.burdens <- finalCalls.burdens %>%
	bind_rows(
		finalCalls.burdens %>%
			left_join(
				sensitivity %>%
					select(
						analysis_id, individual_id, sample_id, chromgroup, filtergroup, call_class, call_type, SBSindel_call_type,
						sensitivity, sensitivity_source
						),
				by = join_by(analysis_id, individual_id, sample_id, chromgroup, filtergroup, call_class, call_type, SBSindel_call_type)
			) %>%
			
			#Calculate corrected counts, burdens, and Poisson 95% confidence intervals for number of calls and burdens
			mutate(
				num_calls = num_calls / sensitivity,
				sensitivity_corrected = TRUE,
				
				ci = num_calls %>% map( function(x){cipoisson(x)} ),
				
				num_calls_lci = map_dbl(ci,1),
				num_calls_uci = map_dbl(ci,2),
				
				burden_calls = num_calls / interrogated_bases_or_bp,
				burden_calls_lci = num_calls_lci / interrogated_bases_or_bp,
				burden_calls_uci = num_calls_uci / interrogated_bases_or_bp
			) %>%
			select(-ci)
	) %>%
	
	#Set sensitivity and sensitivity_source to NA when sensitivity_corrected = FALSE
	mutate(
		sensitivity = if_else(sensitivity_corrected == TRUE, sensitivity, NA),
		sensitivity_source = if_else(sensitivity_corrected == TRUE, sensitivity_source, NA)
	)

#Output as tsv
finalCalls.burdens %>%
	nest_by(chromgroup, filtergroup, .keep = TRUE) %>%
	pwalk(
		function(...){
			x <- list(...)
			
			output_basename_full <- str_c(
				str_c(x$chromgroup, finalCalls.burdens_dir, output_basename),
				x$chromgroup,
				x$filtergroup,
				x$call_class,
				x$call_type,
				x$SBSindel_call_type,
				sep="."
			)
			
			x$data %>% write_tsv(str_c(output_basename_full, ".finalCalls.burdens.tsv"))
			
		}
	)

cat("DONE\n")

######################
### Output estimated SBS mutation error probability
######################
cat("## Outputting estimated SBS mutation error probability...")

#Output as tsv
estimatedSBSMutationErrorProbability %>%
	pwalk(
		function(...){
			x <- list(...)
			
			output_basename_full <- str_c(
				str_c(x$chromgroup, estimatedSBSMutationErrorProbability_dir, output_basename),
				x$chromgroup,
				x$filtergroup,
				x$call_class,
				x$call_type,
				x$SBSindel_call_type,
				sep="."
			)
			
			metadata <- x %>%
				keep(
					names(.) %in%
						c(
							"analysis_id", "individual_id", "sample_id",
							"chromgroup", "filtergroup"
						)
				) %>%
				as_tibble
			
			#by_channel_pyr error probability tsv
			metadata %>%
				mutate(by_channel_pyr = list(x$by_channel_pyr)) %>%
				unnest(by_channel_pyr) %>%
				write_tsv(str_c(output_basename_full, ".estimatedSBSMutationErrorProbability.by_channel_pyr.tsv"))
			
			#total error probability tsv
			metadata %>%
				mutate(total = list(x$total)) %>%
				unnest(total) %>%
				write_tsv(str_c(output_basename_full, ".estimatedSBSMutationErrorProbability.total.tsv"))
		}
	)

cat("DONE\n")

######################
### Output qs2 data object
######################
cat("## Outputting qs2 data object...")

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
		bam.gr.filtertrack.bytype.coverage_tnc,
		genome.reftnc,
		genome_chromgroup.reftnc,
		sensitivity,
		finalCalls.burdens,
		estimatedSBSMutationErrorProbability
	)
	str_c(output_basename,".qs2")
)

cat("DONE\n")