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

#Create output directories
filterStats_dir <- "/filterStats/"
finalCalls_dir <- "/finalCalls/"
germlineVariantCalls_dir <- "/germlineVariantCalls/"
finalCalls.spectra_dir <- "/finalCalls.spectra/"
interrogatedBases.spectra_dir <- "/interrogatedBases.spectra/"
genome.spectra_dir <- "/genome.spectra/"

for(i in chromgroups){
	for(j in c(filterStats_dir, finalCalls_dir, germlineVariantCalls_dir, finalCalls.spectra_dir, interrogatedBases.spectra_dir, genome.spectra_dir))
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
	
	# INFO: everything except the basic VCF columns and fields you want to drop
	info_cols <- names(calls) %>%
		setdiff(c("seqnames","start_refspace","ref_plus_strand","alt_plus_strand","qual"))
	
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
		finalCalls_for_tsv = finalCalls_for_tsv %>%
			map(function(x){
				x %>%
					rename(
						refstrand.refstrand_plus_minus_read = strand.refstrand_plus_minus_read,
						start_refspace = start,
						end_refspace = end
					)
			}),
		
		finalCalls_unique_for_tsv = finalCalls_unique_for_tsv %>%
			map(function(x){
				if(!is.null(x)){
					x %>%
						rename(
							start_refspace = start,
							end_refspace = end
						)
				}else{
					x
				}
			}),
		
		finalCalls_for_vcf = finalCalls_for_vcf %>%
			map(function(x){
				x %>%
					rename(
						refstrand.refstrand_plus_minus_read = strand.refstrand_plus_minus_read,
						start_refspace = start
					)
			}),
		
		finalCalls_unique_for_vcf = finalCalls_unique_for_vcf %>%
			map(function(x){
				if(!is.null(x)){
					x %>%
						rename(
							start_refspace = start
						)
				}else{
					x
				}
			})
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
### Output final calls and germline variant calls
######################

cat("## Outputting final calls and germline variant calls...")

#Output final calls to tsv and vcf, separately for each combination of call_class, call_type, SBSindel_call_type
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

#Create new tibbles without 'for_vcf' tables since not needed anymore, and unnest 'for_tsv' tables, since now that tsv files have been output per call type, we don't need them separated by type.
finalCalls_for_tsv <- finalCalls.bytype %>%
	select(-c(finalCalls_unique_for_tsv, finalCalls_for_vcf, finalCalls_unique_for_vcf)) %>%
	unnest(finalCalls_for_tsv)

finalCalls_unique_for_tsv <- finalCalls.bytype %>%
	select(-c(finalCalls_for_tsv, finalCalls_for_vcf, finalCalls_unique_for_vcf)) %>%
	unnest(finalCalls_unique_for_tsv)

rm(finalCalls.bytype)
invisible(gc())

#Output germline variant calls
 #Reformat list columns to be comma-delimited bounded by square brackets.
for(i in chromgroups){
	for(j in filtergroups){

		#Output path
		output_basename_full <- str_c(str_c(i,germlineVariantCalls_dir,output_basename),i,j,sep=".")
		
		#Create and format output tibble
		germlineVariantCalls.out <- germlineVariantCalls %>%
			filter(chromgroup == i, filtergroup == j) %>%
			mutate(
				across(
					where(is.list),
					function(x){x %>% map_chr(function(v){str_c("[",str_c(v, collapse = ","),"]")})}
				)
			) %>%
			group_by( #Fields that are identical between strands. Not including call_class,call_type,SBSindel_call_type as these are identical for all rows within each finalCalls after nest_join. Not including deletion.bothstrands.startendmatch, since not a field of interest.
				analysis_id,individual_id,sample_id,chromgroup,filtergroup,
				analysis_chunk,run_id,zm,
				seqnames,start_refspace,end_refspace,ref_plus_strand,alt_plus_strand,
				reftnc_plus_strand,alttnc_plus_strand,reftnc_pyr,alttnc_pyr,
				indel_width
			) %>%
			arrange(refstrand) %>% #Sort by reference genome aligned strand
			summarize(
				across( #Collapse to one row fields that differ between strands
					c(refstrand,start_queryspace,end_queryspace,
						qual,qual.opposite_strand,sa,sa.opposite_strand,
						sm,sm.opposite_strand,sx,sx.opposite_strand,
						ref_synthesized_strand,alt_synthesized_strand,
						ref_template_strand,alt_template_strand,
						reftnc_synthesized_strand,alttnc_synthesized_strand,
						reftnc_template_strand,alttnc_template_strand),
					function(y){
						y %>% as.character %>% replace_na("NA") %>% str_c(collapse=",")
					},
					.names="{.col}.refstrand_plus_minus_read"
				),
				.groups = "drop"
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
				out_vcf = str_c(output_basename_full,"germlineVariantCalls.vcf")
			)
		
		rm(germlineVariantCalls.out)
		invisible(gc())
		
	}
}

cat("DONE\n")

######################
### Output trinucleotide distributions
######################
#Output trinucleotide counts

cat("## Outputting trinucleotide distributions...")

#Helper function: write TSV for a given tibble col_name
write_col <- function(df.input, col_name.input, metadata.input, output_basename_full.input) {
	metadata.input %>%
		mutate(!!sym(col_name.input) := list(df.input)) %>%
		unnest(all_of(col_name.input)) %>%
		write_tsv(str_c(output_basename_full.input, ".", col_name.input, ".tsv"))
}

#Helper function: plots spectrum for a given sigfit col_name
plot_col <- function(df.input, col_name.input, output_basename_full.input){
	df.input %>%
		plot_spectrum(str_c(output_basename_full.input, ".", col_name.input, ".pdf"))
}

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

- Assign the use_chromgroup sensitivity to all chromgroups, or if use_chromgroup is null, every chromgroup will already have the default sensitivity tibble assigned by the calculateBurdens script

cat("DONE\n")

######################
### Output call burdens
######################
cat("## Outputting call burdens...")

-> Use sensitivity from matching filtergroup and from use_chromgroup, or from any analysis if use_chromgroup is null (since all identically set SBS and indel senstivity to 1)
- calculate for uncorrected and tnc_corrected burdens

cat("DONE\n")

######################
### Output estimated SBS mutation error probability
######################
cat("## Outputting estimated SBS mutation error probability...")

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
		finalCalls_for_tsv,
		finalCalls_unique_for_tsv,
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