#!/usr/bin/env -S Rscript --vanilla

#filterCalls.R:
# Filters calls
# Note: for all newly created 'passfilter' columns, TRUE = read/call successfully passed the filter and should be analyzed.

cat("#### Running filterCalls ####\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(rtracklayer))
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
	make_option(c("-f", "--file"), type = "character", default=NULL,
							help="path to extractCalls qs2 file"),
	make_option(c("-s", "--sample_id_toanalyze"), type = "character", default=NULL,
							help="sample_id to analyze"),
	make_option(c("-g", "--chromgroup_toanalyze"), type = "character", default=NULL,
	            help="chromgroup to analyze"),
	make_option(c("-v", "--filtergroup_toanalyze"), type = "character", default=NULL,
	            help="filtergroup to analyze"),
	make_option(c("-o", "--output"), type = "character", default=NULL,
							help="output qs2 file")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$config) | is.null(opt$file) | is.null(opt$sample_id_toanalyze) | is.null(opt$chromgroup_toanalyze) | is.null(opt$filtergroup_toanalyze) | is.null(opt$output) ){
	stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))
extractCallsFile <- opt$file
sample_id_toanalyze <- opt$sample_id_toanalyze
chromgroup_toanalyze <- opt$chromgroup_toanalyze
filtergroup_toanalyze <- opt$filtergroup_toanalyze
outputFile <- opt$output

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

 #GRanges of chroms to analyze, required as input for some analyses
genome_chromgroup.gr <- yaml.config$BSgenome$BSgenome_name %>%
	get %>%
	seqinfo %>%
	GRanges %>%
	unname %>%
	filter(seqnames %in% !!chroms_toanalyze)

 #individual_id of this sample_id
individual_id_toanalyze <- yaml.config$samples %>%
  bind_rows %>%
  filter(sample_id == sample_id_toanalyze) %>%
  pull(individual_id)

 #cache_dir
cache_dir <- str_c(yaml.config$cache_dir,"/")

 #germline vcf filter parameters
germline_vcf_types_config <- yaml.config$germline_vcf_types %>%
  enframe(name=NULL) %>%
  unnest_wider(value)

 #call types (restrict to selected chromgroup_toanalyze and filtergroup_toanalyze
call_types_toanalyze <- yaml.config$call_types %>%
  enframe(name=NULL) %>%
  unnest_wider(value) %>%
  unnest_longer(SBSindel_call_types) %>%
  unnest_wider(SBSindel_call_types) %>%
	filter(
		analyzein_chromgroups == "all" | (analyzein_chromgroups %>% str_split(",") %>% map(str_trim) %>% map_lgl(~ !!chromgroup_toanalyze %in% .x)),
		filtergroup == filtergroup_toanalyze
		) %>%
	mutate(filtergroup = filtergroup %>% factor)
  
 #filter group parameters (restrict to selected filtergroup_toanalyze)
filtergroup_toanalyze_config <- yaml.config$filtergroups %>%
  enframe(name=NULL) %>%
  unnest_wider(value) %>%
	filter(filtergroup == filtergroup_toanalyze)
  
 #region filter parameters (restrict to selected chromgroup_toanalyze and filtergroup_toanalyze)
region_read_filters_config <- yaml.config$region_filters %>%
  map("read_filters") %>%
  flatten %>%
  enframe(name=NULL) %>%
  unnest_wider(value) %>%
  mutate(region_filter_threshold_file = str_c(cache_dir,"/",basename(region_filter_file),".bin",binsize,".",threshold,".bw")) %>%
	filter(
	  applyto_chromgroups == "all" | (applyto_chromgroups %>% str_split(",") %>% map(str_trim) %>% map_lgl(~ !!chromgroup_toanalyze %in% .x)),
	  applyto_filtergroups == "all" | (applyto_filtergroups %>% str_split(",") %>% map(str_trim) %>% map_lgl(~ !!filtergroup_toanalyze %in% .x))
	)

region_genome_filters_config <- yaml.config$region_filters %>%
  map("genome_filters") %>%
  flatten %>%
  enframe(name=NULL) %>%
  unnest_wider(value) %>%
  mutate(region_filter_threshold_file = str_c(cache_dir,"/",basename(region_filter_file),".bin",binsize,".",threshold,".bw")) %>%
	filter(
	  applyto_chromgroups == "all" | (applyto_chromgroups %>% str_split(",") %>% map(str_trim) %>% map_lgl(~ !!chromgroup_toanalyze %in% .x)),
	  applyto_filtergroups == "all" | (applyto_filtergroups %>% str_split(",") %>% map(str_trim) %>% map_lgl(~ !!filtergroup_toanalyze %in% .x))
	)

#Display basic configuration parameters
cat("> Processing:\n")
cat("    extractCalls File:",extractCallsFile,"\n")
cat("    individual_id:",individual_id_toanalyze,"\n")
cat("    sample_id:",sample_id_toanalyze,"\n")
cat("    chromgroup:",chromgroup_toanalyze,"\n")
cat("    filtergroup:",filtergroup_toanalyze,"\n")

cat("DONE\n")

######################
### Load custom shared functions
######################
source(Sys.which("sharedFunctions.R"))

######################
### Define custom functions
######################

#Function to calculate number of molecules and reference space bases passing each individual filter matching a pattern (starting from the upstream filtered bam of extractCalls, i.e. starting after minstrandoverlapfilter), and passing all filters so far (if all_filter_suffix is not NULL), from bam input.
calculate_molecule_stats_frombam <- function(bam.input, individual_filter_pattern, all_filter_suffix=NULL, chromgroup_toanalyze.input, filtergroup_toanalyze.input, extractCallsFile.input){
	
	#Number of molecules and reference space bases passing each basic molecule filter, each applied separately
	filter_stats <- bam.input %>%
		group_by(run_id) %>%
		summarize(
			across(
				matches(individual_filter_pattern),
				list(
					num_molecules = ~ n_distinct(zm[.x]),
					num_refspacebases = ~ sum(end[.x] - start[.x])
				),
				.names = "{.fn}_individualfilter.{.col}"
			)
		)
	
	#Number of molecules and reference space bases passing all the basic molecule filters
	if(!is.null(all_filter_suffix)){
		filter_stats <- filter_stats %>%
			left_join(
				bam.input %>%
					filter(if_all(matches("passfilter$"), ~ .x == TRUE)) %>%
					group_by(run_id) %>%
					summarize(
						!!str_c("num_molecules_remaining.",all_filter_suffix) := n_distinct(zm),
						!!str_c("num_refspacebases_remaining.",all_filter_suffix) := sum(end-start)
					),
				by = "run_id"
			)
	}
	
	filter_stats %>%
		complete(run_id) %>%
		mutate(across(-run_id, ~ .x %>% replace_na(0))) %>%
		
		#Pivot longer
		pivot_longer(
			-run_id,
			names_to = "stat",
			values_to = "value"
		) %>%
		
		#Annotate extractCallsFile, chromgroup, filtergroup of this process
		mutate(
			chromgroup = chromgroup_toanalyze.input,
			filtergroup = filtergroup_toanalyze.input,
			extractCallsFile = extractCallsFile.input %>% basename
		) %>%
		
		#Reorder columns
		relocate(run_id,chromgroup,filtergroup,extractCallsFile)
}

#Function to calculate number of molecules and refspace bases remaining in bam filtertracker.
calculate_molecule_stats_frombamfiltertrack <- function(bam.gr.filtertrack.input, stat_label.suffix){
	
	bam.gr.filtertrack.input %>%
		as_tibble %>%
		group_by(run_id) %>%
		summarize(
			!!str_c("num_molecules_remaining.",stat_label.suffix) := n_distinct(zm),
			!!str_c("num_refspacebases_remaining.",stat_label.suffix) := sum(end-start)
		) %>%
		complete(run_id) %>%
		mutate(across(-run_id, ~ .x %>% replace_na(0))) %>%
		pivot_longer(
			-run_id,
			names_to = "stat",
			values_to = "value"
		) %>%
		mutate(
			chromgroup = chromgroup_toanalyze,
			filtergroup = filtergroup_toanalyze,
			extractCallsFile = extractCallsFile %>% basename
		)
}

#Function returning a TRUE/FALSE vector signifying for each element of a query GRanges if it overlaps any element of a subject GRanges, joining on optional columns. If overlap_adjacent_query_insertion = TRUE, then insertions (zero width) in query will be considered an overlap to an immediately adjacent non-insertion subject range.
overlapsAny_bymcols <- function(query, subject, join_mcols = character(), ignore.strand = TRUE, overlap_adjacent_query_insertion = FALSE) {
	
	stopifnot(inherits(query, "GenomicRanges"), inherits(subject, "GenomicRanges"))
	
	nq <- length(query)
	
	if(nq == 0){return(logical())} #Empty query
	if(length(subject) == 0){return(rep(FALSE, nq))} #Empty subject
	
	hits <- findOverlaps(query, subject, ignore.strand = ignore.strand)
	
	if(overlap_adjacent_query_insertion == TRUE){
		hits_adjacent_query_insertion <- findOverlaps(query, subject, ignore.strand = ignore.strand, maxgap = 0)
		
		q_hits_adjacent_query_insertion <- queryHits(hits_adjacent_query_insertion)
		s_hits_adjacent_query_insertion <- subjectHits(hits_adjacent_query_insertion)
		
		hits_adjacent_query_insertion <- hits_adjacent_query_insertion[
			query[q_hits_adjacent_query_insertion] %>% width == 0 &
				subject[s_hits_adjacent_query_insertion] %>% width != 0
			]
		
		hits <- GenomicRanges::union(hits, hits_adjacent_query_insertion)
		rm(hits_adjacent_query_insertion)
	}
	
	if (length(hits) == 0) return(rep(FALSE, nq)) #No overlap
	
	qh <- queryHits(hits)
	sh <- subjectHits(hits)
	
	#Check join_mcols equality
	if (length(join_mcols) == 0) {
	  return(tabulate(qh, nbins = nq) > 0L) #No keys, so every hit valid
	}
	
	q_mcols <- mcols(query)
	s_mcols <- mcols(subject)
	
	if (!all(join_mcols %in% colnames(q_mcols)) || !all(join_mcols %in% colnames(s_mcols))){
		stop("All `join_mcols` must exist as metadata columns in both query and subject.")
	}
	
	#Find shared join_mcols; note: any NA in any join_mcols column ⇒ non-match
	same_join_mcols <- vctrs::vec_equal(
	  as.data.frame(q_mcols[qh, join_mcols, drop = FALSE]),
	  as.data.frame(s_mcols[sh, join_mcols, drop = FALSE]),
	  na_equal = FALSE
	  )
	same_join_mcols[is.na(same_join_mcols)] <- FALSE

	#Output result
	tabulate(qh[same_join_mcols], nbins = nq) > 0L
}

#Convert positions of bases from query space to reference space, while retaining run_id and zm info
#iranges.input: IRangesList with query space positions, with the IRangesList length the same as the number of rows in bam.input
convert_query_to_refspace <- function(irangeslist.input, bam.input, BSgenome_name.input){
	
	#Format irangeslist.input to GRanges of positions failing subreads filters in query space, with annotated bam read id
	irangeslist.input <- irangeslist.input %>%
		setNames(bam.input %>% pull(seqnames)) %>%
		as.data.frame %>%
		as_tibble %>%
		makeGRangesFromDataFrame(
			seqnames = "group_name",
			keep.extra.columns=TRUE
		) %>%
		setNames(.$group) %>%
		select(-group)
	
	if(length(irangeslist.input) == 0){
		return(
			GRanges(run_id=factor(), zm=integer())
		)
	}
	
	mapFromAlignments( #Main function
		
		irangeslist.input,
		
		#GAlignments of reads in reference space
		bam.input %>%
			select(seqnames, start, end, cigar) %>% 
			makeGRangesFromDataFrame(
				keep.extra.columns=TRUE,
				seqinfo=BSgenome_name.input %>% get %>% seqinfo
			) %>%
			as("GAlignments") %>%
			setNames(1:nrow(bam.input))
	) %>%
		
		#Annotate result with run_id, zm, and strand
		as_tibble %>%
		left_join(
			bam.input %>%
				select(run_id, zm, strand) %>%
				mutate(row_id = row_number()),
			by = join_by(alignmentsHits == row_id)
		) %>%
		select(-c(width, xHits, alignmentsHits, strand.x)) %>%
		rename(strand = strand.y) %>%
		makeGRangesFromDataFrame(
			keep.extra.columns=TRUE,
			seqinfo=BSgenome_name.input %>% get %>% seqinfo
		)
}

#Function to check a parameter on both strands (x and y), per a specified minimum threshold and method (mean, all, any)
min_threshold_eachstrand <- function(x, y, threshold, mode = c("mean","all","any")){
	chosenmode <- switch(
		mode,
		mean = function(v) (mean(v) %>% replace_na(0)) >= threshold,
		all = function(v) all(v >= threshold) %>% replace_na(0),
		any = function(v) any(v >= threshold) %>% replace_na(0)
	)

	chosenmode(x) && chosenmode(y) && !any(is.na(y))
}

######################
### Load read and extracted call information
######################
cat("## Loading reads and extracted calls...")

#Load extractedCalls RDS file
extractedCalls <- qs_read(extractCallsFile)

#Load run metadata
run_metadata <- extractedCalls %>% pluck("run_metadata")

#Load bam reads, keeping only reads in selected chroms_toanalyze.
bam <- extractedCalls %>%
	pluck("bam") %>%
	filter(seqnames %in% chroms_toanalyze)

#Load calls, keeping:
# a. calls in selected chroms_toanalyze with call_type/call_class/SBSindel_call_type in call_types_toanalyze (i.e. call_types in selected chromgroup_toanalyze and filtergroup_toanalyze), and mark these with call_toanalyze = TRUE
# b. call_class = SBS or indel (needed for later max calls/mutations postVCF, read indel region filters, and downstream sensitivity and estimated SBS mutation error calculations).
calls <- extractedCalls %>%
	pluck("calls") %>%
	left_join(
		call_types_toanalyze %>%
			distinct(call_type,call_class,SBSindel_call_type) %>%
			mutate(call_toanalyze = TRUE),
		by = join_by(call_type,call_class,SBSindel_call_type)
	) %>%
  mutate(call_toanalyze = call_toanalyze %>% replace_na(FALSE)) %>%
	filter(
	  seqnames %in% chroms_toanalyze,
	  call_toanalyze == TRUE | call_class %in% c("SBS","indel")
		)

#Load prior molecule stats for this chromgroup, for all_chromgroups, and for all_chroms
molecule_stats <- extractedCalls %>%
	pluck("molecule_stats") %>%
	filter(chromgroup %in% c(chromgroup_toanalyze,"all_chromgroups","all_chroms"))

#Create tibble to track genome region filter stats and initialize with stats for whole genome and for chromgroup.
genome_numbases <- yaml.config$BSgenome$BSgenome_name %>%
	get %>%
	seqinfo %>%
	GRanges %>%
	width %>%
	sum

chromgroup_toanalyze_numbases <- genome_chromgroup.gr %>%
	width %>%
	sum

region_genome_filter_stats <- tibble(
	filter = "None (whole genome)",
	binsize = NA,
	threshold = NA,
	padding = NA,
	region_filter_threshold_file = NA,
	num_genomebases_individually_filtered = 0L,
	num_genomebases_remaining = genome_numbases
) %>%
	bind_rows(
		tibble(
			filter = "chromgroup",
			binsize = NA,
			threshold = NA,
			padding = NA,
			region_filter_threshold_file = NA,
			num_genomebases_individually_filtered = genome_numbases - chromgroup_toanalyze_numbases,
			num_genomebases_remaining = chromgroup_toanalyze_numbases
		)
	)

rm(extractedCalls)
invisible(gc())

cat("DONE\n")

######################
### Basic molecule filters
######################
cat("## Applying basic molecule filters...")

#Annotate reads for filters: rq, ec, mapq, num_softclipbases
bam <- bam %>%
	mutate(num_softclipbases = cigar %>% cigarOpTable %>% as_tibble %>% pull(S)) %>%
	group_by(run_id,zm) %>%
	mutate(
		
		#rq
		min_rq_eachstrand.passfilter = all(rq >= filtergroup_toanalyze_config$min_rq_eachstrand),
		min_rq_avgstrands.passfilter = mean(rq) >= filtergroup_toanalyze_config$min_rq_avgstrands,
		
		#ec
		min_ec_eachstrand.passfilter = all(ec >= filtergroup_toanalyze_config$min_ec_eachstrand),
		min_ec_avgstrands.passfilter = mean(ec) >= filtergroup_toanalyze_config$min_ec_avgstrands,
		
		#mapq
		min_mapq_eachstrand.passfilter = all(mapq >= filtergroup_toanalyze_config$min_mapq_eachstrand),
		min_mapq_avgstrands.passfilter = mean(mapq) >= filtergroup_toanalyze_config$min_mapq_avgstrands,
		
		#num_softclipbases
		max_num_softclipbases_eachstrand.passfilter = all(num_softclipbases <= filtergroup_toanalyze_config$max_num_softclipbases_eachstrand),
		max_num_softclipbases_avgstrands.passfilter = mean(num_softclipbases) <= filtergroup_toanalyze_config$max_num_softclipbases_avgstrands
		
	) %>%
	select(-num_softclipbases) %>%
	ungroup

#Annotate reads for filters: max_num_SBScalls_eachstrand, max_num_SBScalls_stranddiff, max_num_indelcalls_eachstrand, max_num_indelcalls_stranddiff
bam <- bam %>%
	left_join(
			
			#Filter for SBS and indel calls
			calls %>%
				filter(
				  call_class %in% c("SBS","indel")
				  ) %>%
			  
				#Change call_class to factor so that count results are listed for both SBS and indel below.
				mutate(call_class = call_class %>% factor(levels=c("SBS","indel"))) %>% 
				
				#Count number of SBS and indel calls per strand while completing missing strand and SBS/indel values for each molecule so that num_{call_class}calls are calculated correctly
				count(run_id,zm,strand,call_class) %>%
			  complete(run_id,zm,strand,call_class,fill = list(n=0)) %>%
				pivot_wider(
					names_from = call_class,
					values_from = n,
					names_glue = "num_{call_class}calls",
					names_expand = TRUE #Necessary if there are no calls
					) %>%
				
				#Calculate filters for each molecule
				group_by(run_id,zm) %>%
				summarize(
					#max_num_SBScalls_eachstrand.passfilter
					max_num_SBScalls_eachstrand.passfilter = all(num_SBScalls <= filtergroup_toanalyze_config$max_num_SBScalls_eachstrand),
					
					#max_num_indelcalls_eachstrand.passfilter
					max_num_indelcalls_eachstrand.passfilter = all(num_indelcalls <= filtergroup_toanalyze_config$max_num_indelcalls_eachstrand),
					
					#max_num_SBScalls_stranddiff.passfilter
					max_num_SBScalls_stranddiff.passfilter = (num_SBScalls * if_else(strand == "+",  1L, -1L)) %>%
						sum %>%
						abs %>%
						#Curly braces ensure the '.' refers to the immediately preceding result
						{. <= filtergroup_toanalyze_config$max_num_SBScalls_stranddiff},
					
					#max_num_indelcalls_stranddiff.passfilter
					max_num_indelcalls_stranddiff.passfilter = (num_indelcalls * if_else(strand == "+",  1L, -1L)) %>%
						sum %>%
						abs %>%
					  {. <= filtergroup_toanalyze_config$max_num_indelcalls_stranddiff},
	
					.groups = "drop"
				),
		by = join_by(run_id,zm)
	) %>%
	mutate(across(
			c(
				max_num_SBScalls_eachstrand.passfilter,
				max_num_indelcalls_eachstrand.passfilter,
				max_num_SBScalls_stranddiff.passfilter,
				max_num_indelcalls_stranddiff.passfilter
			),
			~ . %>% replace_na(TRUE)
	))

#Annotate reads for filters: max_num_SBSmutations, max_num_indelmutations
bam <- bam %>%
	left_join(
			
			#Filter for SBS and indel mutations. Not filtering here on only '+' or only '-' strand since the below code won't double count mutations, because it counts mutations for each strand and those numbers will be identical for each strand of a molecule.
		calls %>%
			filter(
				call_class %in% c("SBS","indel"),
				SBSindel_call_type == "mutation"
				) %>%
		  
		  #Count number of SBS and indel mutations per molecule while completing missing SBS/indel values for each molecule so that num_{call_class}mutations are calculated correctly.
		  mutate(call_class = call_class %>% factor(levels=c("SBS","indel"))) %>% 
		  count(run_id,zm,strand,call_class) %>%
		  complete(run_id,zm,strand,call_class,fill = list(n=0)) %>%
			pivot_wider(
				names_from=call_class,
				values_from=n,
				names_glue = "num_{call_class}mutations",
				names_expand = TRUE #Necessary if there are no calls
			) %>%
			
			#Calculate filters for each molecule
			group_by(run_id,zm) %>%
			summarize(
				#max_num_SBSmutations
				max_num_SBSmutations.passfilter = all(num_SBSmutations <= filtergroup_toanalyze_config$max_num_SBSmutations),
				
				#max_num_indelmutations
				max_num_indelmutations.passfilter = all(num_indelmutations <= filtergroup_toanalyze_config$max_num_indelmutations),
				
				.groups = "drop"
			),
		by = join_by(run_id,zm)
	) %>%
	mutate(across(
		c(
			max_num_SBSmutations.passfilter,
			max_num_indelmutations.passfilter
		),
		~ . %>% replace_na(TRUE)
	))

#Annotate calls with read filters
calls <- calls %>%
  left_join(
    bam %>%
      distinct(
        run_id,
        zm,
        pick(contains("passfilter"))
      ),
    by = join_by(run_id,zm)
  )

#Calculate basic molecule filter stats
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombam(
			bam.input = bam,
			individual_filter_pattern = "passfilter$",
			all_filter_suffix = "passallbasicmoleculefilters",
			chromgroup_toanalyze.input = chromgroup_toanalyze,
			filtergroup_toanalyze.input = filtergroup_toanalyze,
			extractCallsFile.input = extractCallsFile
		)
	) %>%
	relocate(run_id,chromgroup,filtergroup,extractCallsFile)

cat("DONE\n")

######################
### Germline VCF variant filters
######################
cat("## Applying germline VCF variant filters...")

#Load germline VCF filter data
germline_vcf_variants <- qs_read(str_c(cache_dir,"/",individual_id_toanalyze,".germline_vcf_variants.qs2")) %>%
  as_tibble
  
#Filter to keep germline VCF variants that pass configured germline VCF filters and keep only columns necessary for downstream filtering
germline_vcf_variants <- germline_vcf_variants  %>%
  group_by(germline_vcf_type) %>% #split by germline_vcf_type
  nest %>%
  deframe %>%
  imap(
    function(x,idx){
      filters <- germline_vcf_types_config %>% filter(germline_vcf_type == idx)
      
      x %>%
        filter(
          (
            call_class == "SBS" &
              FILTER %in% (filters$SBS_FILTERS %>% unlist) &
              Depth >= filters$SBS_min_Depth &
              GQ >= filters$SBS_min_GQ &
              VAF >= filters$SBS_min_VAF &
              QUAL >= filters$SBS_min_QUAL
          ) |
          (
            call_class == "indel" &
              FILTER %in% (filters$indel_FILTERS %>% unlist) &
              Depth >= filters$indel_min_Depth &
              GQ >= filters$indel_min_GQ &
              VAF >= filters$indel_min_VAF &
              QUAL >= filters$indel_min_QUAL
          )
        ) %>%
        select(seqnames,start,end,ref_plus_strand,alt_plus_strand,call_class,call_type,germline_vcf_file) %>%
        distinct #in case atomization/multi-allelic records created duplicate records
    }
  ) %>%
  bind_rows(.id="germline_vcf_type") #Combine back to one tibble

#Annotate calls matching germline VCF variants, and for MDB calls set germline_vcf.passfilter = TRUE regardless if there is a matching germline VCF variant, since we do not want to filter out MDBs just because they are in the location of a germline variant.
calls <- calls %>%
  left_join(
    germline_vcf_variants %>%
      group_by(seqnames,start,end,ref_plus_strand,alt_plus_strand) %>%
      summarize(
        germline_vcf_types_detected = str_c(germline_vcf_type,collapse=","),
        germline_vcf_files_detected = str_c(germline_vcf_file,collapse=","),
        .groups="drop"
      ) %>%
      mutate(germline_vcf.passfilter = FALSE),
    by = join_by(seqnames,start,end,ref_plus_strand,alt_plus_strand)
  ) %>%
  mutate(
    germline_vcf.passfilter = germline_vcf.passfilter %>% replace_na(TRUE),
    germline_vcf.passfilter = if_else(call_class == "MDB", TRUE, germline_vcf.passfilter)
    )

cat("DONE\n")

######################
### Post-germline VCF variant filtering molecule filters
######################
cat("## Applying post-germline VCF variant filtering molecule filters...")

#Annotate reads for filters: max_num_SBScalls_postVCF_eachstrand, max_num_indelcalls_postVCF_eachstrand
bam <- bam %>%
  left_join(
    
    #Filter to keep SBS and indel calls that pass all filters applied so far
    calls %>%
      filter(
        call_class %in% c("SBS","indel"),
        if_all(contains("passfilter"), ~ .x == TRUE)
      ) %>%
      
      #Count number of SBS and indel calls per strand while completing missing strand and SBS/indel values for each molecule so that num_{call_class}calls are calculated correctly. Fix strand levels to '+/-' to avoid expanding strand = '*'.
      mutate(
        call_class = call_class %>% factor(levels=c("SBS","indel")),
        strand = strand %>% factor(levels=c("+","-"))
      ) %>% 
      count(run_id,zm,strand,call_class) %>%
      complete(run_id, zm, strand, call_class, fill = list(n=0)) %>%
      pivot_wider(
        names_from=call_class,
        values_from=n,
        names_glue = "num_{call_class}calls",
        names_expand = TRUE #Necessary if there are no calls
      ) %>%
      
      #Calculate filters for each molecule
      group_by(run_id,zm) %>%
      summarize(
        #max_num_SBScalls_postVCF_eachstrand
        max_num_SBScalls_postVCF_eachstrand.passfilter = all(num_SBScalls <= filtergroup_toanalyze_config$max_num_SBScalls_postVCF_eachstrand),
        
        #max_num_indelcalls_postVCF_eachstrand
        max_num_indelcalls_postVCF_eachstrand.passfilter = all(num_indelcalls <= filtergroup_toanalyze_config$max_num_indelcalls_postVCF_eachstrand),
        ,.groups = "drop"
      )
    ,by = join_by(run_id,zm)
  ) %>%
	mutate(across(
		c(
			max_num_SBScalls_postVCF_eachstrand.passfilter,
			max_num_indelcalls_postVCF_eachstrand.passfilter
		),
		~ . %>% replace_na(TRUE)
	))

#Annotate reads for filters: max_num_SBSmutations_postVCF, max_num_indelmutations_postVCF
bam <- bam %>%
  left_join(
    
    #Filter to keep SBS and indel calls that pass all filters applied so far
    calls %>%
      filter(
        call_class %in% c("SBS","indel"),
        SBSindel_call_type == "mutation",
        if_all(contains("passfilter"), ~ .x == TRUE)
      ) %>%
      
      #Count number of SBS and indel mutations per molecule while completing missing SBS/indel values for each molecule so that num_{call_class}mutations are calculated correctly
      mutate(call_class = call_class %>% factor(levels=c("SBS","indel"))) %>% 
      count(run_id,zm,strand,call_class) %>%
      complete(run_id,zm,strand,call_class,fill = list(n=0)) %>%
      pivot_wider(
        names_from=call_class,
        values_from=n,
        names_glue = "num_{call_class}mutations",
        names_expand = TRUE #Necessary if there are no calls
      ) %>%
      
      #Calculate filters for each molecule
      group_by(run_id,zm) %>%
      summarize(
        #max_num_SBSmutations_postVCF
        max_num_SBSmutations_postVCF.passfilter = all(num_SBSmutations <= filtergroup_toanalyze_config$max_num_SBSmutations_postVCF),
        
        #max_num_indelmutations_postVCF
        max_num_indelmutations_postVCF.passfilter = all(num_indelmutations <= filtergroup_toanalyze_config$max_num_indelmutations_postVCF),
        
        .groups = "drop"
      ),
    by = join_by(run_id,zm)
  )  %>%
	mutate(across(
		c(
			max_num_SBSmutations_postVCF.passfilter,
			max_num_indelmutations_postVCF.passfilter
		),
		~ . %>% replace_na(TRUE)
	))

#Update molecule_stats
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombam(
			bam.input = bam,
			individual_filter_pattern = "postVCF.*passfilter$",
			all_filter_suffix = "passallbasicmoleculeandpostgermlineVCFfilters",
			chromgroup_toanalyze.input = chromgroup_toanalyze,
			filtergroup_toanalyze.input = filtergroup_toanalyze,
			extractCallsFile.input = extractCallsFile
		)
	)

#Annotate calls with new read filters
calls <- calls %>%
  left_join(
    bam %>%
      distinct(
        run_id,
        zm,
        max_num_SBScalls_postVCF_eachstrand.passfilter,
        max_num_indelcalls_postVCF_eachstrand.passfilter,
        max_num_SBSmutations_postVCF.passfilter,
        max_num_indelmutations_postVCF.passfilter
      )
    ,by = join_by(run_id,zm)
  )

cat("DONE\n")

######################
### Create filter trackers
######################

cat("## Creating filter trackers...")

#Create bam read GRanges to track genome region filtering, containing only molecules that have passed all filters so far
bam.gr.filtertrack <- bam %>%
	filter(if_all(matches("passfilter$"), ~ .x == TRUE)) %>%
	select(run_id, zm, seqnames, start, end, strand) %>%
	makeGRangesFromDataFrame(
		keep.extra.columns = TRUE,
		seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
	)

#Create genome GRanges of selected chroms_toanalyze to track genome region filtering
genome_chromgroup.gr.filtertrack <- genome_chromgroup.gr

#Convert calls tibble back to GRanges, keeping only metadata columns necessary for its subsequent use
calls.gr <- calls %>%
	makeGRangesFromDataFrame(
		keep.extra.columns=TRUE,
		seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
	) %>%
	select(run_id, zm, call_class, call_type, SBSindel_call_type, ref_plus_strand, alt_plus_strand)

cat("DONE\n")

######################
### Genome 'N' base sequences
######################
cat("## Applying genome 'N' base sequence filter...")

passfilter_label <- "region_genome_filter_Nbases.passfilter"

#Extract 'N' base sequence ranges
region_genome_filter <- yaml.config$BSgenome$BSgenome_name %>%
	get %>%
	vmatchPattern("N",.) %>%
	GenomicRanges::reduce(ignore.strand=TRUE)

#Subtract regions from filter trackers
bam.gr.filtertrack <- bam.gr.filtertrack %>%
	GRanges_subtract(region_genome_filter, ignore.strand=TRUE)

genome_chromgroup.gr.filtertrack <- genome_chromgroup.gr.filtertrack %>%
	GRanges_subtract(region_genome_filter)

#Annotate filtered calls. Annotates even if partial overlap (relevant for deletions).
calls[[passfilter_label]] <- ! overlapsAny_bymcols(
	calls.gr, region_genome_filter,
	ignore.strand = TRUE
)

#Record number of filtered genome bases to stats
region_genome_filter_stats <- region_genome_filter_stats %>%
  bind_rows(
  	tibble(
  		filter = "Nbases",
  		binsize = NA,
  		threshold = NA,
  		padding = NA,
  		num_genomebases_individually_filtered = region_genome_filter %>%
  			width %>%
  			sum,
  		num_genomebases_remaining = genome_chromgroup.gr.filtertrack %>%
  			width %>%
  			sum
  	)
	)

#Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack,
			stat_label.suffix = "passNbasesfilter"
		)
	)

rm(region_genome_filter)
invisible(gc())

cat("DONE\n")

######################
### Base-quality filters
######################
cat("## Applying base-quality filters...")

#Convert positions of bases failing min_qual filter from query space to reference space, while retaining run_id and zm info
min_qual.fail.gr <- bam %>%
	pull(qual) %>%
	PhredQuality %>%
	as("IntegerList") %>%
	as("RleList") %>%
	#Handle empty RleList due to empty bam
	{
		if(length(.) > 0){
			IRanges::slice(., upper=filtergroup_toanalyze_config$min_qual, includeUpper=FALSE, rangesOnly=TRUE)
		}else{
			IRangesList()
		}
	} %>%
	convert_query_to_refspace(bam.input = bam, BSgenome_name.input = yaml.config$BSgenome$BSgenome_name)

#Subtract bases below min_qual threshold from bam reads filter tracker, joining on run_id and zm. Set ignore.strand = TRUE so that if a position is filtered in one strand, also filter the same reference space position on the opposite strand, since that is how calls are filtered.
bam.gr.filtertrack <- bam.gr.filtertrack %>%
  GRanges_subtract_bymcols(min_qual.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

#Annotate calls below min_qual threshold. If either strand fails the filter, both strands are set to fail. If any qual value in the opposite strand contains NA, consider that failure to pass the filter.
if(! filtergroup_toanalyze_config$min_qual_method %in% c("mean","all","any")){
  stop("Incorrect min_qual_method setting!")
}

calls <- calls %>%
  mutate(
  	min_qual.passfilter = map2_lgl(
	  	qual,
	  	qual.opposite_strand,
	  	~ min_threshold_eachstrand(
	  		.x, .y,
	  		filtergroup_toanalyze_config$min_qual,
	  		filtergroup_toanalyze_config$min_qual_method
	  	)
	  )
  )

#Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack,
			stat_label.suffix = "passminqualfilter"
		)
	)

rm(min_qual.fail.gr)
invisible(gc())

cat("DONE\n")

######################
### Molecule position filters
######################
cat("## Applying molecule position filters...")

#Create GRanges of read ends
bam.gr.onlyranges <- bam %>%
  select(seqnames, start, end, strand, run_id, zm) %>%
  makeGRangesFromDataFrame(
    keep.extra.columns=TRUE,
    seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
  )

bam.gr.trim <- c(
    bam.gr.onlyranges %>%
	    mutate(
	      temp_end = start + filtergroup_toanalyze_config$read_trim_bp - 1,
	      end = if_else(temp_end < end, temp_end, end)  #prevents exceeding genome boundaries
	      ) %>%
	    select(-temp_end),
    
    bam.gr.onlyranges %>%
      mutate(
        temp_start = end - filtergroup_toanalyze_config$read_trim_bp + 1,
        start = if_else(temp_start > start, temp_start, start)
      ) %>%
      select(-temp_start)
  )

rm(bam.gr.onlyranges)
invisible(gc())

#Subtract bases filtered by read_trim_bp from bam reads filter tracker, joining on run_id and zm. If a position is filtered in one strand, also filter the same reference space position on the opposite strand, since that is how calls are filtered.
bam.gr.filtertrack <- bam.gr.filtertrack %>%
  GRanges_subtract_bymcols(bam.gr.trim, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

#Update molecule stats
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack,
			stat_label.suffix = "passreadtrimbpfilter"
		)
	)

#Annotate calls filtered by read_trim_bp without taking strand into account.
calls[["read_trim_bp.passfilter"]] <- ! overlapsAny_bymcols(
	calls.gr, bam.gr.trim,
	join_mcols = c("run_id","zm"),
	ignore.strand = TRUE
)

rm(bam.gr.trim)
invisible(gc())

cat("DONE\n")

######################
### Read SBS region filter
######################
cat("## Applying read SBS region filter...")

#Extract padding configuration
sbs_flank <- filtergroup_toanalyze_config %>% pull(ccs_sbs_flank)

if(!is.null(sbs_flank) && sbs_flank == 0){sbs_flank <- NULL}

#Create GRanges of sbs flanks with configured padding
if(!is.null(sbs_flank)){
	read_sbs_region_filter <- c(
		#Left flank
		calls.gr %>%
			filter(call_class=="SBS") %>%
			flank(width=sbs_flank, start = TRUE, ignore.strand = TRUE) %>%
			suppressWarnings %>% #remove warnings of out of bounds regions due to resize
			trim,
		
		#Right flank
		calls.gr %>%
			filter(call_class=="SBS") %>%
			flank(width=sbs_flank, start = FALSE, ignore.strand = TRUE) %>%
			suppressWarnings %>% #remove warnings of out of bounds regions due to resize
			trim
	) %>%
		select(run_id, zm)
}else{
	read_sbs_region_filter <- calls.gr %>%
		select(run_id, zm) %>%
		slice(0)
}
	
#Annotate calls filtered by read_sbs_region_filter without taking strand into account, so that if a call on either strand fails the filter, calls on both strands fail the filter. Applied to both indels and SBS calls. Insertions immediately adjacent to an SBS will be filtered using the overlap_adjacent_query_insertion = TRUE option.
calls[["read_sbs_region_filter.passfilter"]] <- ! overlapsAny_bymcols(
	calls.gr, read_sbs_region_filter,
	join_mcols = c("run_id","zm"),
	ignore.strand = TRUE,
	overlap_adjacent_query_insertion = TRUE
)

#Subtract read_indel_region_filter bases from bam reads filter tracker, joining on run_id and zm. Set ignore.strand = TRUE so that if a position is filtered in one strand, also filter the same reference space position on the opposite strand, since that is how calls are filtered.
bam.gr.filtertrack <- bam.gr.filtertrack %>%
	GRanges_subtract_bymcols(read_sbs_region_filter, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

#Record number of molecule reference space bases remaining after filters
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack,
			stat_label.suffix = "passreadsbsregionfilter"
		)
	)

rm(read_sbs_region_filter)
invisible(gc())

cat("DONE\n")

######################
### Read indel region filters
######################
cat("## Applying read indel region filters (not applied to indels)...")

#Extract padding configuration
ccs_ins_pad <- filtergroup_toanalyze_config %>% pull(ccs_ins_pad)
ccs_del_pad <- filtergroup_toanalyze_config %>% pull(ccs_del_pad)

#Create GRanges of insertions with configured padding
if(!is.null(ccs_ins_pad)){
	ccs_ins_pad <- ccs_ins_pad %>%
	  tibble(pad=.) %>%
	  extract(pad, into = c("m", "b"), regex = "m(\\d+)b(\\d+)", convert = TRUE)
	
	read_ins_region_filter <- calls.gr %>%
		filter(call_class=="indel", call_type == "insertion") %>%
		mutate(
			padding_m = ccs_ins_pad$m * nchar(alt_plus_strand) %>% round %>% as.integer,
			padding_b = rep(ccs_ins_pad$b, length(.)),
			start = start - pmax(padding_m,padding_b),
			end = end + pmax(padding_m,padding_b)
		) %>%
		suppressWarnings %>% #remove warnings of out of bounds regions due to resize
		trim %>%
		select(run_id, zm)
	
}else{
	read_ins_region_filter <- calls.gr %>%
		select(run_id, zm) %>%
		slice(0)
}

#Create GRanges of deletions with configured padding
if(!is.null(ccs_del_pad)){
	ccs_del_pad <- ccs_del_pad %>%
	  tibble(pad=.) %>%
	  extract(pad, into = c("m", "b"), regex = "m(\\d+)b(\\d+)", convert = TRUE)
	
	read_del_region_filter <- calls.gr %>%
		filter(call_class=="indel", call_type == "deletion") %>%
		mutate(
			padding_m = ccs_del_pad$m * nchar(ref_plus_strand) %>% round %>% as.integer,
			padding_b = rep(ccs_del_pad$b, length(.)),
			start = start - pmax(padding_m,padding_b),
			end = end + pmax(padding_m,padding_b)
		) %>%
		suppressWarnings %>% #remove warnings of out of bounds regions due to resize
		trim %>%
		select(run_id, zm)
	
}else{
	read_del_region_filter <- calls.gr %>%
		select(run_id, zm) %>%
		slice(0)
}

#Combine insertion and deletion GRanges with configured padding
read_indel_region_filter <- c(read_ins_region_filter, read_del_region_filter)
 
#Annotate calls filtered by read_indel_region_filter without taking strand into account, so that if a call on either strand fails the filter, calls on both strands fail the filter. Not applied to indels.
calls[["read_indel_region_filter.passfilter"]] <- ! overlapsAny_bymcols(
	calls.gr, read_indel_region_filter,
	join_mcols = c("run_id","zm"),
	ignore.strand = TRUE
)

calls <- calls %>%
  mutate(
  	read_indel_region_filter.passfilter = if_else(call_class == "indel", TRUE, read_indel_region_filter.passfilter)
  	)

#Subtract read_indel_region_filter bases from bam reads filter tracker, joining on run_id and zm. Set ignore.strand = TRUE so that if a position is filtered in one strand, also filter the same reference space position on the opposite strand, since that is how calls are filtered. Create new bam.gr.filtertrack.nonindels and .indels, since these regions are not used to filter indels.
bam.gr.filtertrack.indelanalysis <- bam.gr.filtertrack

bam.gr.filtertrack.nonindelanalysis <- bam.gr.filtertrack %>%
  GRanges_subtract_bymcols(read_indel_region_filter, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

rm(bam.gr.filtertrack)
invisible(gc())

#Record number of molecule reference space bases remaining after filters
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.indelanalysis,
			stat_label.suffix = "indelanalysis.passreadindelregionfilter"
		),
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.nonindelanalysis,
			stat_label.suffix = "nonindelanalysis.passreadindelregionfilter"
		)
	)

rm(read_ins_region_filter, read_del_region_filter, read_indel_region_filter)
invisible(gc())

cat("DONE\n")

######################
### Germline BAM total reads (coverage) filter
######################
cat("## Applying germline BAM total reads (coverage) filter...")

passfilter_label <- "min_germlineBAM_TotalReads.passfilter"

#Load germline BAM samtools mpileup bigwig; only bases failing filter due to memory constraints.
germline_bam_samtools_mpileup_file <- cache_dir %>%
	str_c(
		yaml.config$individuals %>%
			bind_rows %>%
			filter(individual_id == individual_id_toanalyze) %>%
			pull(germline_bam_file) %>%
			unique %>%
			basename,
		".bw"
	)

tmpchromsizes <- tempfile(tmpdir=getwd(),pattern=".",fileext=".bed")
system(paste("/bin/bash -c",shQuote(paste(
	"awk '{print $1 \"\t0\t\" $2}'", yaml.config$genome_fai,
	"| sort -k1,1 -k2,2n >",
	tmpchromsizes
)
)))

tmpbw <- tempfile(tmpdir=getwd(),pattern=".")

system(paste("/bin/bash -c",shQuote(paste(
	yaml.config$wiggletools_bin, "lt",
	filtergroup_toanalyze_config$min_germlineBAM_TotalReads,
	"trim", tmpchromsizes, "fillIn", tmpchromsizes,
	germline_bam_samtools_mpileup_file, "|",
	yaml.config$wigToBigWig_bin, "stdin <(cut -f 1,2",
	yaml.config$genome_fai,")",
	tmpbw
	)
)))

germline_bam_samtools_mpileup_filter <- tmpbw %>%
	import(format = "bigWig") %>%
	{
		seqlevels(.) <- seqlevels(genome_chromgroup.gr)
		seqinfo(.) <- seqinfo(genome_chromgroup.gr)
		.
	} %>%
	select(-score)

invisible(file.remove(tmpchromsizes, tmpbw))

#Subtract filtered regions from filter trackers
bam.gr.filtertrack.indelanalysis <- bam.gr.filtertrack.indelanalysis %>%
	GRanges_subtract(germline_bam_samtools_mpileup_filter,ignore.strand = TRUE)

bam.gr.filtertrack.nonindelanalysis <- bam.gr.filtertrack.nonindelanalysis %>%
	GRanges_subtract(germline_bam_samtools_mpileup_filter,ignore.strand = TRUE)

genome_chromgroup.gr.filtertrack <- genome_chromgroup.gr.filtertrack %>%
	GRanges_subtract(germline_bam_samtools_mpileup_filter,ignore.strand = TRUE)

#Annotate calls with germline BAM samtools mpileup bigwig filter. For insertions, left and right flanking bases, and for deletions, the deleted bases, in genome reference space must pass the filter.
calls[[passfilter_label]] <- ! overlapsAny_bymcols(
	calls.gr %>%
		mutate(
			startnew = if_else(call_type == "insertion", end, start),
			endnew = if_else(call_type == "insertion", start, end),
			start = startnew,
			end = endnew
		),
	germline_bam_samtools_mpileup_filter,
	ignore.strand = TRUE
)

#Record number of filtered genome bases to stats
region_genome_filter_stats <- region_genome_filter_stats %>%
	bind_rows(
		tibble(
			filter = "min_germlineBAM_TotalReads",
			binsize = NA,
			threshold = NA,
			padding = NA,
			num_genomebases_individually_filtered = germline_bam_samtools_mpileup_filter %>%
				width %>%
				sum,
			num_genomebases_remaining = genome_chromgroup.gr.filtertrack %>%
				width %>%
				sum
		)
	)

#Record number of molecule reference space bases remaining after filters
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.indelanalysis,
			stat_label.suffix = "indelanalysis.passgermlinebamtotalreadsfilter"
		),
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.nonindelanalysis,
			stat_label.suffix = "nonindelanalysis.passgermlinebamtotalreadsfilter"
		)
	)

rm(germline_bam_samtools_mpileup_filter)
invisible(gc())

cat("DONE\n")

######################
### Germline BAM variant read filters
######################
#Apply only to SBS calls, since germline BAM mpileup calling of indels is too noisy

cat("## Applying germline BAM variant read filters (only applied to SBS calls)...")

#Format variant regions to VCF POS coordinates for extracting variants from germline BAM bcftools mpileup. Retrieve only SBS calls since only annotating SBS calls. Filter is applied both to SBS mutation and SBS non-mutation calls.
variant_regions_bcftools_mpileup <- calls %>%
	filter(call_class == "SBS") %>%
	rename(
		start_refspace = start,
		end_refspace = end
	) %>%
	select(seqnames,start_refspace,end_refspace)

#Load germline BAM bcftools mpileup VCF
germline_bam_bcftools_mpileup_file <- cache_dir %>%
	str_c(
		yaml.config$individuals %>%
			bind_rows %>%
			filter(individual_id == individual_id_toanalyze) %>%
			pull(germline_bam_file) %>%
			unique %>%
			basename,
		".vcf.gz"
	)

germline_bam_bcftools_mpileup_filter <- load_vcf(
	vcf_file = germline_bam_bcftools_mpileup_file,
	regions = variant_regions_bcftools_mpileup,
	genome_fasta = yaml.config$genome_fasta,
	BSgenome_name = yaml.config$BSgenome$BSgenome_name,
	bcftools_bin = yaml.config$bcftools_bin
)

#Extract variants failing filters
germline_bam_bcftools_mpileup_BAMVariantReads_filter <- germline_bam_bcftools_mpileup_filter %>%
	filter(AD2 > filtergroup_toanalyze_config$max_germlineBAM_VariantReads)

germline_bam_bcftools_mpileup_BAMVAF_filter <- germline_bam_bcftools_mpileup_filter %>%
	filter(VAF > filtergroup_toanalyze_config$max_germlineBAM_VAF)

#Annotate SBS calls with germline BAM bcftools mpileup VCF filters.
calls[["max_germlineBAM_VariantReads.passfilter"]] <- ! overlapsAny_bymcols(
	calls.gr, germline_bam_bcftools_mpileup_BAMVariantReads_filter,
	join_mcols = c("call_class","call_type","ref_plus_strand","alt_plus_strand"), #Not joining by SBSindel_call_type so filter is also applied to mismatch-ss and mismatch-ds calls 
	ignore.strand = TRUE
)

calls[["max_germlineBAM_VAF.passfilter"]] <- ! overlapsAny_bymcols(
	calls.gr, germline_bam_bcftools_mpileup_BAMVAF_filter,
	join_mcols = c("call_class","call_type","ref_plus_strand","alt_plus_strand"), #Not joining by SBSindel_call_type so filter is also applied to mismatch-ss and mismatch-ds calls
	ignore.strand = TRUE
)

rm(variant_regions_bcftools_mpileup, germline_bam_bcftools_mpileup_filter)
invisible(gc())

cat("DONE\n")

######################
### Subread filters
######################
cat("## Applying subread filters...")

#Convert positions of bases failing subread filters from query space to reference space, while retaining run_id and zm info
 #min_frac_subreads_cvg (sa/max(sa)) 
frac_subreads_cvg.fail.gr <- bam %>%
		pull(sa) %>%
		map(
			function(x){
				vals <- x %>% as.vector
				(vals/max(vals)) %>%
					{. < filtergroup_toanalyze_config$min_frac_subreads_cvg} %>%
					which %>%
					IRanges
				}
		) %>%
		IRangesList %>%
		convert_query_to_refspace(bam.input = bam, BSgenome_name.input = yaml.config$BSgenome$BSgenome_name)

 #min_num_subreads_match (sm)
num_subreads_match.fail.gr <- bam %>%
	pull(sm) %>%
	map(
		function(x){
			(x < filtergroup_toanalyze_config$min_num_subreads_match) %>%
				which %>%
				IRanges
		}
	) %>%
	IRangesList %>%
	convert_query_to_refspace(bam.input = bam, BSgenome_name.input = yaml.config$BSgenome$BSgenome_name)

 #min_frac_subreads_match (sm / [sm+sx])
frac_subreads_match.fail.gr <- map2(
		.x = bam$sm,
		.y = bam$sx,
		function(x,y){
			(x/(x+y)) %>%
				replace_na(0) %>%
				{. < filtergroup_toanalyze_config$min_frac_subreads_match} %>%
				which %>%
				IRanges
			}
	) %>%
	IRangesList %>%
	convert_query_to_refspace(bam.input = bam, BSgenome_name.input = yaml.config$BSgenome$BSgenome_name)

#Subtract bases below filter thresholds from bam reads filter tracker, joining on run_id and zm. Set ignore.strand = TRUE so that if a position is filtered in one strand, also filter the same reference space position on the opposite strand, since that is how calls are filtered.
bam.gr.filtertrack.indelanalysis <- bam.gr.filtertrack.indelanalysis %>%
	GRanges_subtract_bymcols(frac_subreads_cvg.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

bam.gr.filtertrack.nonindelanalysis <- bam.gr.filtertrack.nonindelanalysis %>%
	GRanges_subtract_bymcols(frac_subreads_cvg.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

 #Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.indelanalysis,
			stat_label.suffix = "indelanalysis.passminfracsubreadscvgfilter"
		),
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.nonindelanalysis,
			stat_label.suffix = "nonindelanalysis.passminfracsubreadscvgfilter"
		)
	)

bam.gr.filtertrack.indelanalysis <- bam.gr.filtertrack.indelanalysis %>%
	GRanges_subtract_bymcols(num_subreads_match.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

bam.gr.filtertrack.nonindelanalysis <- bam.gr.filtertrack.nonindelanalysis %>%
	GRanges_subtract_bymcols(num_subreads_match.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

 #Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.indelanalysis,
			stat_label.suffix = "indelanalysis.passminnumsubreadsmatchfilter"
		),
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.nonindelanalysis,
			stat_label.suffix = "nonindelanalysis.passminnumsubreadsmatchfilter"
		)
	)

bam.gr.filtertrack.indelanalysis <- bam.gr.filtertrack.indelanalysis %>%
	GRanges_subtract_bymcols(frac_subreads_match.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

bam.gr.filtertrack.nonindelanalysis <- bam.gr.filtertrack.nonindelanalysis %>%
	GRanges_subtract_bymcols(frac_subreads_match.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

 #Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.indelanalysis,
			stat_label.suffix = "indelanalysis.passminfracsubreadsmatchfilter"
		),
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.nonindelanalysis,
			stat_label.suffix = "nonindelanalysis.passminfracsubreadsmatchfilter"
		)
	)

#Annotate calls below filter thresholds. If either strand fails the filter, both strands are set to fail. If sa, sm, or sx in the opposite strand contain NA, consider that failure to pass the filter.
if(! filtergroup_toanalyze_config$min_subreads_cvgmatch_method %in% c("mean","all","any")){
	stop("Incorrect min_qual_method setting!")
}

#min_frac_subreads_cvg
calls <- calls %>%
	#Extract max sa of strand and opposite strand, required for frac_subreads_cvg calculations
	left_join(
		bam %>%
			mutate(
				max_sa = map_int(sa, function(x){x %>% as.vector %>% max})
			) %>%
			select(run_id,zm,strand,max_sa),
		by = join_by(run_id,zm,strand)
	) %>%
	left_join(
		bam %>%
			mutate(
				strand = if_else(strand == "+","-","+") %>% factor,
				max_sa.opposite_strand = map_int(sa, function(x){x %>% as.vector %>% max})
				) %>%
			select(run_id,zm,strand,max_sa.opposite_strand),
		by = join_by(run_id,zm,strand)
	) %>%
	
	#Calculate frac_subreads_cvg for both strands
	mutate(
		frac_subreads_cvg = map2(
			.x = sa,
			.y = max_sa,
			function(x,y){x/y}
		),
		frac_subreads_cvg.opposite_strand = map2(
			.x = sa.opposite_strand,
			.y = max_sa.opposite_strand,
			function(x,y){x/y}
		)
	) %>%
	
	#Perform filtering. NA in any base on the opposite strand is assigned passfilter = FALSE.
	mutate(
		min_frac_subreads_cvg.passfilter = map2_lgl(
			frac_subreads_cvg,
			frac_subreads_cvg.opposite_strand,
			~ min_threshold_eachstrand(
				.x, .y,
				filtergroup_toanalyze_config$min_frac_subreads_cvg,
				filtergroup_toanalyze_config$min_subreads_cvgmatch_method
			)
		)
	) %>%
	select(-c(max_sa,max_sa.opposite_strand,frac_subreads_cvg,frac_subreads_cvg.opposite_strand))

#min_num_subreads_match. NA in any base on the opposite strand is assigned passfilter = FALSE.
calls <- calls %>%
	mutate(
		min_num_subreads_match.passfilter = map2_lgl(
			sm,
			sm.opposite_strand,
			~ min_threshold_eachstrand(
				.x, .y,
				filtergroup_toanalyze_config$min_num_subreads_match,
				filtergroup_toanalyze_config$min_subreads_cvgmatch_method
			)
		)
	)
	
#min_frac_subreads_match
calls <- calls %>%
	#Calculate frac_subreads_match for both strands
	mutate(
		frac_subreads_match = map2(
			.x = sm,
			.y = sx,
			function(x,y){x/(x+y)}
		),
		frac_subreads_match.opposite_strand = map2(
			.x = sm.opposite_strand,
			.y = sx.opposite_strand,
			function(x,y){x/(x+y)}
		)
	) %>%
	
	#Perform filtering. If any base on the call strand or opposite strand has (sm+sx) = 0, evaluate frac_subreads_match as 0. If any opposite strand base has sm or sx = NA, then assign passfilter = FALSE.
	mutate(
		min_frac_subreads_match.passfilter = map2_lgl(
			frac_subreads_match,
			frac_subreads_match.opposite_strand,
			~ min_threshold_eachstrand(
				.x, .y,
				filtergroup_toanalyze_config$min_frac_subreads_match,
				filtergroup_toanalyze_config$min_subreads_cvgmatch_method
			)
		)
	) %>%
	select(-c(frac_subreads_match,frac_subreads_match.opposite_strand))
	
rm(frac_subreads_cvg.fail.gr, num_subreads_match.fail.gr, frac_subreads_match.fail.gr)
invisible(gc())

cat("DONE\n")

######################
### Duplex coverage filters
######################
cat("## Applying duplex coverage filters...")

#Extract read regions remaining after all filtering that do not have duplex (i.e. both strand) coverage

#Regions only in plus or minus strand reads
#indel analysis
duplex_coverage.fail.indelanalysis.gr <- c(
	#Plus strand only
	bam.gr.filtertrack.indelanalysis %>%
		filter(strand=="+") %>%
		GRanges_subtract_bymcols(bam.gr.filtertrack.indelanalysis %>% filter(strand=="-"), join_mcols = c("run_id","zm"), ignore.strand = TRUE),
	
	#Minus strand only
	bam.gr.filtertrack.indelanalysis %>%
		filter(strand=="-") %>%
		GRanges_subtract_bymcols(bam.gr.filtertrack.indelanalysis %>% filter(strand=="+"), join_mcols = c("run_id","zm"), ignore.strand = TRUE)
)

#non indel analysis
duplex_coverage.fail.nonindelanalysis.gr <- c(
	#Plus strand only
	bam.gr.filtertrack.nonindelanalysis %>%
		filter(strand=="+") %>%
		GRanges_subtract_bymcols(bam.gr.filtertrack.nonindelanalysis %>% filter(strand=="-"), join_mcols = c("run_id","zm"), ignore.strand = TRUE),
	
	#Minus strand only
	bam.gr.filtertrack.nonindelanalysis %>%
		filter(strand=="-") %>%
		GRanges_subtract_bymcols(bam.gr.filtertrack.nonindelanalysis %>% filter(strand=="+"), join_mcols = c("run_id","zm"), ignore.strand = TRUE)
)

#Subtract these regions from the bam.gr filter tracker
bam.gr.filtertrack.indelanalysis <- bam.gr.filtertrack.indelanalysis %>%
	GRanges_subtract_bymcols(duplex_coverage.fail.indelanalysis.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

bam.gr.filtertrack.nonindelanalysis <- bam.gr.filtertrack.nonindelanalysis %>%
	GRanges_subtract_bymcols(duplex_coverage.fail.nonindelanalysis.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

#Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.indelanalysis,
			stat_label.suffix = "indelanalysis.passduplexcoveragefilter"
		),
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.nonindelanalysis,
			stat_label.suffix = "nonindelanalysis.passduplexcoveragefilter"
		)
	)

##Annotate which calls have duplex coverage after all filtering by comparing to final filtered bam.gr filter tracker without taking strand into account, so that if a call on either strand fails the filter, calls on both strands fail the filter. Use appropriate indel or nonindel analysis filter tracker.
calls[["duplex_coverage.indelanalysis.passfilter"]] <- ! overlapsAny_bymcols(
	calls.gr, duplex_coverage.fail.indelanalysis.gr,
	join_mcols = c("run_id","zm"),
	ignore.strand = TRUE
)

calls[["duplex_coverage.nonindelanalysis.passfilter"]] <- ! overlapsAny_bymcols(
	calls.gr, duplex_coverage.fail.nonindelanalysis.gr,
	join_mcols = c("run_id","zm"),
	ignore.strand = TRUE
)

calls <- calls %>%
	mutate(
		duplex_coverage.passfilter = if_else(
			call_class == "indel",
			duplex_coverage.indelanalysis.passfilter,
			duplex_coverage.nonindelanalysis.passfilter,
		)
	) %>%
	select(-c(duplex_coverage.indelanalysis.passfilter,duplex_coverage.nonindelanalysis.passfilter))

rm(duplex_coverage.fail.indelanalysis.gr,duplex_coverage.fail.nonindelanalysis.gr)
invisible(gc())

cat("DONE\n")

######################
### Germline VCF indel region filters
######################
cat("## Applying germline VCF indel region filters...")

#Create a copy of the bam filter trackers from which this and subsequent germline-related filters will be excluded, for use in downstream sensitivity analysis
bam.gr.filtertrack.indelanalysis.except_germline_filters <- bam.gr.filtertrack.indelanalysis
bam.gr.filtertrack.nonindelanalysis.except_germline_filters <- bam.gr.filtertrack.nonindelanalysis

#Split germline_vcf_variants by germline_vcf_file
germline_vcf_variants <- germline_vcf_variants  %>%
	group_by(germline_vcf_file) %>% 
	nest %>%
	deframe

for(i in names(germline_vcf_variants)){
	
	#Construct label for the newly created filter column
	passfilter_label <- i %>%
		basename %>%
		str_c("germline_vcf_indel_region_filter_",.,".passfilter")
	
	#Extract germline VCF type and padding configuration
	vcf_type <- germline_vcf_variants %>% pluck(i,"germline_vcf_type",1)
	indel_ins_pad <- germline_vcf_types_config %>% filter(germline_vcf_type == vcf_type) %>% pull(indel_ins_pad)
	indel_del_pad <- germline_vcf_types_config %>% filter(germline_vcf_type == vcf_type) %>% pull(indel_del_pad)
	
	#Create GRanges of insertions with configured padding
	if(!is.null(indel_ins_pad)){
		indel_ins_pad <- indel_ins_pad %>%
			tibble(pad=.) %>%
			extract(pad, into = c("m", "b"), regex = "m(\\d+)b(\\d+)", convert = TRUE)
		
		germline_vcf_ins_region_filter <- germline_vcf_variants %>%
			pluck(i) %>%
			filter(call_class=="indel", call_type == "insertion") %>%
			mutate(
				padding_m = indel_ins_pad$m * nchar(alt_plus_strand) %>% round %>% as.integer,
				padding_b = indel_ins_pad$b,
				start = start - pmax(padding_m,padding_b),
				end = end + pmax(padding_m,padding_b)
			)
	}else{
		germline_vcf_ins_region_filter <- germline_vcf_variants %>%
			pluck(i) %>%
			slice(0)
	}
	
	#Create GRanges of deletions with configured padding
	if(!is.null(indel_del_pad)){
		indel_del_pad <- indel_del_pad %>%
			tibble(pad=.) %>%
			extract(pad, into = c("m", "b"), regex = "m(\\d+)b(\\d+)", convert = TRUE)
		
		germline_vcf_del_region_filter <- germline_vcf_variants %>%
			pluck(i) %>%
			filter(call_class=="indel", call_type == "deletion") %>%
			mutate(
				padding_m = indel_del_pad$m * nchar(ref_plus_strand) %>% round %>% as.integer,
				padding_b = indel_del_pad$b,
				start = start - pmax(padding_m,padding_b),
				end = end + pmax(padding_m,padding_b)
			)
	}else{
		germline_vcf_del_region_filter <- germline_vcf_variants %>%
			pluck(i) %>%
			slice(0)
	}
	
	#Combine insertion and deletion GRanges with configured padding
	germline_vcf_indel_region_filter <- c(germline_vcf_ins_region_filter, germline_vcf_del_region_filter) %>%
		makeGRangesFromDataFrame(
			seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
		) %>%
		suppressWarnings %>% #remove warnings of out of bounds regions due to resize
		trim %>%
		GenomicRanges::reduce(ignore.strand=TRUE)
	
	#Subtract filtered regions from filter trackers
	bam.gr.filtertrack.indelanalysis <- bam.gr.filtertrack.indelanalysis %>%
		GRanges_subtract(germline_vcf_indel_region_filter, ignore.strand = TRUE)
	
	bam.gr.filtertrack.nonindelanalysis <- bam.gr.filtertrack.nonindelanalysis %>%
		GRanges_subtract(germline_vcf_indel_region_filter, ignore.strand = TRUE)
	
	genome_chromgroup.gr.filtertrack <- genome_chromgroup.gr.filtertrack %>%
		GRanges_subtract(germline_vcf_indel_region_filter)
	
	#Annotate filtered calls. Annotates even if partial overlap (relevant for deletions).
	calls[[passfilter_label]] <- ! overlapsAny_bymcols(
		calls.gr, germline_vcf_indel_region_filter,
		ignore.strand = TRUE
	)
	
	#Record number of filtered genome bases to stats
	region_genome_filter_stats <- region_genome_filter_stats %>%
		bind_rows(
			tibble(
				filter = str_c(i %>% basename,":indelregion"),
				binsize = NA,
				threshold = NA,
				padding = NA,
				num_genomebases_individually_filtered = germline_vcf_indel_region_filter %>%
					width %>%
					sum,
				num_genomebases_remaining = genome_chromgroup.gr.filtertrack %>%
					width %>%
					sum
			)
		)
}

#Update molecule stats
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.indelanalysis,
			stat_label.suffix = "indelanalysis.passgermlinevcfindelregionfilters"
		),
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.nonindelanalysis,
			stat_label.suffix = "nonindelanalysis.passgermlinevcfindelregionfilters"
		)
	)

rm(germline_vcf_ins_region_filter, germline_vcf_del_region_filter, germline_vcf_indel_region_filter)
invisible(gc())

cat("DONE\n")

######################
### Region-based molecule filters
######################
cat("## Applying region-based molecule filters...")

for(i in seq_len(nrow(region_read_filters_config))){
	
	#Construct label for the newly created filter column
	passfilter_label <- region_read_filters_config %>%
		pluck("region_filter_threshold_file",i) %>%
		basename %>%
		str_c("region_read_filter_",.,".passfilter")
	
	#Extract read threshold
	read_threshold_type <- region_read_filters_config %>%
		pluck("read_threshold",i) %>%
		str_extract("^(gt|gte|lt|lte)")
	
	read_threshold_value <- region_read_filters_config %>%
		pluck("read_threshold",i) %>%
		str_extract("[0-9.]+$")
	
	#Load region read filter bigwig
	region_read_filter <- region_read_filters_config %>%
		pluck("region_filter_threshold_file",i) %>%
		import(format = "bigWig") %>%
		{
			seqlevels(.) <- seqlevels(genome_chromgroup.gr)
			seqinfo(.) <- seqinfo(genome_chromgroup.gr)
			.
		} %>%
		
		#Add padding
		resize(
			width = width(.) + 2 * (region_read_filters_config %>% pluck("padding",i)),
			fix="center"
		) %>%
		suppressWarnings %>% #remove warnings of out of bounds regions due to resize
		trim %>%
		GenomicRanges::reduce(ignore.strand=TRUE) %>%
		
		#Create a logical TRUE/FALSE RLE of genomic bases covered so that the below binnedAverage function calculates, for each read, the fraction covered by the region_read_filter
		coverage %>%
		{. > 0}
	
	#Annotate which molecules pass the filter
	bam <- bam %>%
		left_join(
			#Create ranges-only copy of bam read GRanges on which to check filter
			bam %>%
				select(seqnames, start, end, strand, run_id, zm) %>%
				makeGRangesFromDataFrame(
					keep.extra.columns=TRUE,
					seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
				) %>%
				
				#Calculate fraction of coverage by the region read filter
				binnedAverage(region_read_filter,varname="frac") %>%
				
				#Annotate which molecules have a fraction of their length filtered by the above bigwig (averaged across plus and minus strands) that is either greater than (gt), greater than or equal (gte), less than (lt), or less than or equal (lte) to a read_threshold_value
				as_tibble %>%
				group_by(run_id, zm) %>%
				summarize(
					!!passfilter_label := ! case_when(
						!!read_threshold_type == "gt" ~ mean(frac) > read_threshold_value,
						!!read_threshold_type == "gte" ~ mean(frac) >= read_threshold_value,
						!!read_threshold_type == "lt" ~ mean(frac) < read_threshold_value,
						!!read_threshold_type == "lte" ~ mean(frac) <= read_threshold_value
					),
					.groups = "drop"
				),
			by = join_by(run_id,zm)
		)
	
	#Remove filtered molecules from bam.gr.filtertrack
	bam.gr.filtertrack.indelanalysis <- bam.gr.filtertrack.indelanalysis[
		!vctrs::vec_in(
			bam.gr.filtertrack.indelanalysis %>%
				mcols %>%
				as_tibble %>%
				select(run_id,zm),
			bam %>%
				filter(!!sym(passfilter_label) == FALSE) %>%
				distinct(run_id,zm)
			)
		]
		
		bam.gr.filtertrack.nonindelanalysis <- bam.gr.filtertrack.nonindelanalysis[
			!vctrs::vec_in(
				bam.gr.filtertrack.nonindelanalysis %>%
					mcols %>%
					as_tibble %>%
					select(run_id,zm),
				bam %>%
					filter(!!sym(passfilter_label) == FALSE) %>%
					distinct(run_id,zm)
			)
		]
		
		#Remove filtered molecules from bam.gr.filtertrack.except_germline_filters only if it is a non-germline filter
		if(region_read_filters_config %>% pluck("is_germline_filter",i) == FALSE){
			bam.gr.filtertrack.indelanalysis.except_germline_filters <- bam.gr.filtertrack.indelanalysis.except_germline_filters[
				!vctrs::vec_in(
					bam.gr.filtertrack.indelanalysis.except_germline_filters %>%
						mcols %>%
						as_tibble %>%
						select(run_id,zm),
					bam %>%
						filter(!!sym(passfilter_label) == FALSE) %>%
						distinct(run_id,zm)
				)
			]
			
			bam.gr.filtertrack.nonindelanalysis.except_germline_filters <- bam.gr.filtertrack.nonindelanalysis.except_germline_filters[
				!vctrs::vec_in(
					bam.gr.filtertrack.nonindelanalysis.except_germline_filters %>%
						mcols %>%
						as_tibble %>%
						select(run_id,zm),
					bam %>%
						filter(!!sym(passfilter_label) == FALSE) %>%
						distinct(run_id,zm)
				)
			]
		}
		
		#Update molecule stats of remaining molecules and bases after each filter is added
		molecule_stats <- molecule_stats %>%
			bind_rows(
				calculate_molecule_stats_frombamfiltertrack(
					bam.gr.filtertrack.input = bam.gr.filtertrack.indelanalysis,
					stat_label.suffix = str_c("indelanalysis.",passfilter_label)
				),
				calculate_molecule_stats_frombamfiltertrack(
					bam.gr.filtertrack.input = bam.gr.filtertrack.nonindelanalysis,
					stat_label.suffix = str_c("nonindelanalysis.",passfilter_label)
				)
			)

	rm(region_read_filter)
	invisible(gc())		
}

#Update molecule stats of each filter applied individually regardless of any prior filters. This is possible in this case since these are molecule-level filters. Note, all_filter_suffix = NULL, because can't calculate all_filter stats from bam at this point since non-molecule filters have already been applied upstream.
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombam(
			bam.input = bam,
			individual_filter_pattern = "region_read_filter_.*passfilter$",
			all_filter_suffix = NULL,
			chromgroup_toanalyze.input = chromgroup_toanalyze,
			filtergroup_toanalyze.input = filtergroup_toanalyze,
			extractCallsFile.input = extractCallsFile
		)
	)

#Update molecule stats of remaining molecules and bases after all region read filters are applied
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.indelanalysis,
			stat_label.suffix = "indelanalysis.passregionreadfilters"
		),
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.nonindelanalysis,
			stat_label.suffix = "nonindelanalysis.passregionreadfilters"
		)
	)

#Annotate calls with new read filters
calls <- calls %>%
	left_join(
		bam %>%
			distinct(
				run_id,
				zm,
				pick(matches("region_read_filter.*passfilter"))
			)
		,by = join_by(run_id,zm)
	)

cat("DONE\n")

######################
### Genome region filters
######################
cat("## Applying genome region filters...")

#Perform filtering
for(i in seq_len(nrow(region_genome_filters_config))){
	
	#Construct label for the newly created filter column
	passfilter_label <- region_genome_filters_config %>%
		pluck("region_filter_threshold_file",i) %>%
		basename %>%
		str_c("region_genome_filter_",.,".passfilter")
	
	#Load region read filter bigwig
	region_genome_filter <- region_genome_filters_config %>%
		pluck("region_filter_threshold_file",i) %>%
		import(format = "bigWig") %>%
		{
			seqlevels(.) <- seqlevels(genome_chromgroup.gr)
			seqinfo(.) <- seqinfo(genome_chromgroup.gr)
			.
		} %>%
		
		#Add padding
		resize(
			width = width(.) + 2 * (region_genome_filters_config %>% pluck("padding",i)),
			fix="center"
		) %>%
		suppressWarnings %>% #remove warnings of out of bounds regions due to resize
		trim %>%
		GenomicRanges::reduce(ignore.strand=TRUE)
	
	#Subtract genome region filters from filter trackers
	bam.gr.filtertrack.indelanalysis <- bam.gr.filtertrack.indelanalysis %>%
		GRanges_subtract(region_genome_filter, ignore.strand=TRUE) #GRanges_subtract to retain run_id and zm columns
	
	bam.gr.filtertrack.nonindelanalysis <- bam.gr.filtertrack.nonindelanalysis %>%
		GRanges_subtract(region_genome_filter, ignore.strand=TRUE)
	
	#Remove filtered regions from bam.gr.filtertrack.except_germline_filters only if it is a non-germline filter
	if(region_genome_filters_config %>% pluck("is_germline_filter",i) == FALSE){
		bam.gr.filtertrack.indelanalysis.except_germline_filters <- bam.gr.filtertrack.indelanalysis.except_germline_filters %>%
			GRanges_subtract(region_genome_filter, ignore.strand=TRUE)
		
		bam.gr.filtertrack.nonindelanalysis.except_germline_filters <- bam.gr.filtertrack.nonindelanalysis.except_germline_filters %>%
			GRanges_subtract(region_genome_filter, ignore.strand=TRUE)
	}
	
	genome_chromgroup.gr.filtertrack <- genome_chromgroup.gr.filtertrack %>%
		GRanges_subtract(region_genome_filter)
	
	#Annotate filtered calls. Annotates even if partial overlap (relevant for deletions).
	calls[[passfilter_label]] <- ! overlapsAny_bymcols(
		calls.gr, region_genome_filter,
		ignore.strand = TRUE
	)
	
	#Update molecule stats
	molecule_stats <- molecule_stats %>%
		bind_rows(
			calculate_molecule_stats_frombamfiltertrack(
				bam.gr.filtertrack.input = bam.gr.filtertrack.indelanalysis,
				stat_label.suffix = str_c("indelanalysis.",passfilter_label)
			),
			calculate_molecule_stats_frombamfiltertrack(
				bam.gr.filtertrack.input = bam.gr.filtertrack.nonindelanalysis,
				stat_label.suffix = str_c("nonindelanalysis.",passfilter_label)
			)
		)
	
	#Record number of filtered genome bases to stats
	region_genome_filter_stats <- region_genome_filter_stats %>%
		bind_rows(
			region_genome_filters_config %>%
				slice(i) %>%
				rename(filter = region_filter_file) %>%
				mutate(
					filter = filter %>% basename,
					num_genomebases_individually_filtered = region_genome_filter %>%
						width %>%
						sum,
					num_genomebases_remaining = genome_chromgroup.gr.filtertrack %>%
						width %>%
						sum
				) %>%
				select(-c(applyto_chromgroups,applyto_filtergroups))
		)
	
	rm(region_genome_filter)
	invisible(gc())
}

#Record number of molecule reference space bases remaining after all genome region filters
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.indelanalysis,
			stat_label.suffix = "indelanalysis.passgenomeregionfilters"
		),
		calculate_molecule_stats_frombamfiltertrack(
			bam.gr.filtertrack.input = bam.gr.filtertrack.nonindelanalysis,
			stat_label.suffix = "nonindelanalysis.passgenomeregionfilters"
		)
	)

cat("DONE\n")

######################
### Max final calls per molecule filter
######################
cat("## Applying max final calls per molecule filters...")

#Extract which molecules fail the filter for each call_type / SBSindel_call_type combination
max_finalcalls_eachstrand.filter <- calls %>%
	mutate(passallfilters = if_all(contains("passfilter"), ~ . == TRUE)) %>%
	group_by(run_id,zm,strand,call_type,SBSindel_call_type) %>%
	summarize(num_calls = sum(passallfilters), .groups="drop") %>%
	group_by(run_id,zm,call_type,SBSindel_call_type) %>%
	summarize(
		max_finalcalls_eachstrand.passfilter = all(num_calls <= filtergroup_toanalyze_config$max_finalcalls_eachstrand),
		.groups="drop"
	)

#Subtract filtered molecules from bam.gr.filtertrack, creating a new bam.gr.filtertrack.bytype for each call_type x SBSindel_call_type combination being analyzed. For only bam.gr.filtertrack.bytype, also add a row for SBS/mismatch-ss (if not already there) for later calculation of SBS mutation error probability. Perform for both standard and 'except_germline_filters' filter trackers.
bam.gr.filtertrack.bytype <- call_types_toanalyze %>%
	
	#Remove analyzein_chromgroups and MDB configuration parameters from call_types_toanalyze that are not needed here
	select(-analyzein_chromgroups, -starts_with("MDB")) %>%
	
	#Add row for SBS/mismatch-ss
	bind_rows(
		tibble(
			call_type = "SBS",
			call_class = "SBS",
			SBSindel_call_type = "mismatch-ss",
			filtergroup = filtergroup_toanalyze
		)
	) %>%
	distinct %>%
	
	mutate(
		bam.gr.filtertrack = if_else(
			call_class == "indel",
			bam.gr.filtertrack.indelanalysis %>% list,
			bam.gr.filtertrack.nonindelanalysis %>% list
		),
		bam.gr.filtertrack = pmap(
			list(bam.gr.filtertrack, call_type, SBSindel_call_type),
			function(x,y,z){
				molecules_to_filter <- !!max_finalcalls_eachstrand.filter %>%
					filter(
						call_type == y,
						SBSindel_call_type == z,
						max_finalcalls_eachstrand.passfilter == FALSE
						) %>%
					distinct(run_id,zm)
					
				x[
					!vctrs::vec_in(
						x %>% mcols %>% as_tibble %>% select(run_id,zm),
						molecules_to_filter
					)
				]
				
			}
		)
	)

bam.gr.filtertrack.except_germline_filters.bytype <- call_types_toanalyze %>%
	select(-analyzein_chromgroups, -starts_with("MDB")) %>% #remove columns that are not needed here
	mutate(
		bam.gr.filtertrack = if_else(
			call_class == "indel",
			bam.gr.filtertrack.indelanalysis.except_germline_filters %>% list,
			bam.gr.filtertrack.nonindelanalysis.except_germline_filters %>% list
		),
		bam.gr.filtertrack = pmap(
			list(bam.gr.filtertrack, call_type, SBSindel_call_type),
			function(x,y,z){
				molecules_to_filter <- !!max_finalcalls_eachstrand.filter %>%
					filter(
						call_type == y,
						SBSindel_call_type == z,
						max_finalcalls_eachstrand.passfilter == FALSE
					) %>%
					distinct(run_id,zm)
				
				x[
					!vctrs::vec_in(
						x %>% mcols %>% as_tibble %>% select(run_id,zm),
						molecules_to_filter
					)
				]
				
			}
		)
	)

#Update molecule stats for each call_type x SBSindel_call_type combination being analyzed.
molecule_stats <- molecule_stats %>%
	bind_rows(
		bam.gr.filtertrack.bytype %>%
			#Exclude SBS/mismatch-ss row if it was only added to bam.gr.filtertrack.bytype for downstream calculation of SBS mutation error probability. 
			semi_join(
				call_types_toanalyze,
				by = join_by(call_type, call_class, SBSindel_call_type, filtergroup)
			) %>%
			mutate(
				numbases = pmap(
						list(bam.gr.filtertrack, call_type, SBSindel_call_type),
						function(x,y,z){
							calculate_molecule_stats_frombamfiltertrack(
								bam.gr.filtertrack.input = x,
								stat_label.suffix = str_c(y,".",z,".passmaxfinalcallseachstrandfilter")
							)
						}
				)
			) %>%
			select(numbases) %>%
			deframe %>%
			bind_rows %>%
			relocate(run_id,chromgroup,filtergroup,extractCallsFile)
		)

#Annotate calls with filter
calls <- calls %>%
	left_join(
		max_finalcalls_eachstrand.filter,
		by = join_by(run_id,zm,call_type,SBSindel_call_type)
	) %>%
	mutate(max_finalcalls_eachstrand.passfilter = max_finalcalls_eachstrand.passfilter %>% replace_na(TRUE))

rm(max_finalcalls_eachstrand.filter, bam.gr.filtertrack.indelanalysis, bam.gr.filtertrack.nonindelanalysis, bam.gr.filtertrack.indelanalysis.except_germline_filters, bam.gr.filtertrack.nonindelanalysis.except_germline_filters)
invisible(gc())

cat("DONE\n")

######################
### Save output data
######################
cat("## Saving output data...")

#Remove unecessary indelanalysis or nonindelanalysis stats from molecule_stats
molecule_stats <- molecule_stats %>%
	filter(
		if(!("indel" %in% call_types_toanalyze$call_class)){
			stat %>% str_detect("\\.indelanalysis\\.",negate=TRUE)
		}else{
			TRUE
		}
	) %>%
	filter(
		if(!any(c("SBS","MDB") %in% call_types_toanalyze$call_class)){
			stat %>% str_detect("\\.nonindelanalysis\\.",negate=TRUE)
		}else{
			TRUE
		}
	)

#Save output
qs_save(
	list(
		config = list(
			yaml.config = yaml.config,
			run_metadata = run_metadata,
			extractCallsFile = extractCallsFile,
			individual = individual_id_toanalyze,
			sample_id = sample_id_toanalyze,
			chromgroup = chromgroup_toanalyze,
			genome_chromgroup.gr = genome_chromgroup.gr,
			filtergroup = filtergroup_toanalyze,
			call_types = call_types_toanalyze,
			region_read_filters_config = region_read_filters_config,
			region_genome_filters_config = region_genome_filters_config
		),
		calls = calls,
		bam.gr.filtertrack.bytype = bam.gr.filtertrack.bytype,
		bam.gr.filtertrack.except_germline_filters.bytype = bam.gr.filtertrack.except_germline_filters.bytype,
		genome_chromgroup.gr.filtertrack = genome_chromgroup.gr.filtertrack,
		region_genome_filter_stats = region_genome_filter_stats,
		molecule_stats = molecule_stats
	),
	outputFile
)

cat("DONE\n")