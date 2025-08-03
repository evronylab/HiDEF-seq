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
options(datatable.showProgress = FALSE)
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

if(is.na(opt$config) | is.na(opt$file) | is.na(opt$sample_id_toanalyze) | is.na(opt$chromgroup_toanalyze) | is.na(opt$filtergroup_toanalyze) | is.na(opt$output) ){
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

#General parameters
strand_levels <- c("+","-")

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

 #call types (restrict to selected chromgroup_toanalyze and filtergroup_toanalyze)
call_types_toanalyze <- yaml.config$call_types %>%
  enframe(name=NULL) %>%
  unnest_wider(value) %>%
  unnest_longer(SBSindel_call_types) %>%
  unnest_wider(SBSindel_call_types) %>%
	filter(
		analyzein_chromgroups == "all" | (analyzein_chromgroups %>% str_split(",") %>% map(str_trim) %>% map_lgl(~ !!chromgroup_toanalyze %in% .x)),
		filtergroup == filtergroup_toanalyze
		)
  
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

######################
### Load read and extracted call information
######################
cat("## Loading reads and extracted calls...")

#Load extractedCalls RDS file
extractedCalls <- qs_read(extractCallsFile)

#Load bam reads as tibble, keeping only reads in selected chroms_toanalyze. Note, this transforms qual from PhredQuality to list.
bam <- extractedCalls %>%
	pluck("bam.gr") %>%
	filter(
		seqnames %in% chroms_toanalyze
	) %>%
	as_tibble %>%
	mutate(strand = strand %>% factor(levels=strand_levels))

#Load calls as tibble, keeping only calls in selected chroms_toanalyze with call_type in call_types_toanalyze (i.e. call_types in selected chromgroup_toanalyze and filtergroup_toanalyze) or call_class = SBS or indel (needed for later max calls/mutations postVCF filters and downstream sensitivity calculations). Mark calls with call_type in call_types_toanalyze with a new column call_toanalyze = TRUE. Also fix strand_levels back to '+/-'.
calls <- extractedCalls %>%
	pluck("calls.gr") %>%
  mutate(
  	call_toanalyze = call_type %in% (!!call_types_toanalyze %>% pull(call_type) %>% unique),
  	strand = strand %>% factor(strand_levels)
  	) %>%
	filter(
	  seqnames %in% chroms_toanalyze,
	  call_toanalyze == TRUE | call_class %in% c("SBS","indel")
		) %>%
  as_tibble

#Load prior molecule stats
molecule_stats <- extractedCalls %>%
	pluck("molecule_stats")

#Create tibble to track genome region filter stats
region_genome_filters_stats <- region_genome_filters_config %>%
	mutate(num_bases_filtered = NA_integer_)

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
		  
			#Change call_class to factor so that count results are listed for both SBS and indel below. Fix strand levels to '+/-' to avoid expanding strand = '*'.
			mutate(
				call_class = call_class %>% factor(levels=c("SBS","indel")),
				strand = strand %>% factor(strand_levels)
			) %>% 
			
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
		  
		  #Count number of SBS and indel mutations per molecule while completing missing SBS/indel values for each molecule so that num_{call_class}mutations are calculated correctly. Fix strand levels to '+/-' to avoid expanding strand = '*'.
		  mutate(
		    call_class = call_class %>% factor(levels=c("SBS","indel")),
		    strand = strand %>% factor(levels=strand_levels)
		  ) %>% 
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
			)
		,by = join_by(run_id,zm)
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
      )
    ,by = join_by(run_id,zm)
  )

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
        .groups="drop") %>%
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
      mutate(
        call_class = call_class %>% factor(levels=c("SBS","indel")),
        strand = strand %>% factor(levels=strand_levels)
      ) %>% 
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
### Region-based molecule filters
######################
cat("## Applying region-based molecule filters...")

#Ranges-only copy of bam read GRanges on which to check filters to reduce memory requirements
bam.gr.onlyranges <- bam %>%
	select(seqnames, start, end, strand, run_id, zm) %>%
	makeGRangesFromDataFrame(
		keep.extra.columns=TRUE,
		seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
	)

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
  	subsetByOverlaps(genome_chromgroup.gr) %>%
    
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
      bam.gr.onlyranges %>%
        
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
    
}

rm(bam.gr.onlyranges, region_read_filter)

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
### Calculate molecule filter stats
######################

cat("## Calculating molecule filter stats...")
molecule_stats <- molecule_stats %>%
	bind_rows(
		#Number of molecules passing each individual filter, each applied separately
		bam %>%
			group_by(run_id) %>%
			summarize(across(
				matches("passfilter$"),
				~ n_distinct(zm[.x]),
				.names = "num_molecules_{.col}"
			)) %>%
		
		#Number of molecules passing all the molecule filters
		left_join(
			bam %>%
				filter(if_all(matches("passfilter$"), ~ .x == TRUE)) %>%
				group_by(run_id) %>%
				summarize(num_molecules.passallmoleculefilters = n_distinct(zm)),
			by = "run_id"
		) %>%
	
		#Number of query space bases passing each individual filter each applied separately
		left_join(
			bam %>%
				group_by(run_id) %>%
				summarize(across(
					matches("passfilter$"),
					~ sum(qwidth[.x]),
					.names = "num_queryspacebases_{.col}"
				)),
			by = "run_id"
		) %>%
			
		#Number of reference space bases passing each individual filter each applied separately
		left_join(
			bam %>%
				group_by(run_id) %>%
				summarize(across(
					matches("passfilter$"),
					~ sum(end[.x] - start[.x]),
					.names = "num_refspacebases_{.col}"
				)),
			by = "run_id"
		) %>%
		
		#Number of query space bases passing all the molecule filters
		left_join(
			bam %>%
				filter(if_all(matches("passfilter$"), ~ .x == TRUE)) %>%
				group_by(run_id) %>%
				summarize(num_queryspacebases.passallmoleculefilters = sum(qwidth)),
			by = "run_id"
		) %>%
			
		#Number of reference space bases passing all the molecule filters
		left_join(
			bam %>%
				filter(if_all(matches("passfilter$"), ~ .x == TRUE)) %>%
				group_by(run_id) %>%
				summarize(num_refspacebases.passallmoleculefilters = sum(end-start)),
			by = "run_id"
		) %>%
		
		#Pivot longer
		pivot_longer(
			-run_id,
			names_to = "stat",
			values_to = "value"
		) %>%
			
		#Annotate extractCallsFile, chromgroup, filtergroup of this process
		mutate(
			chromgroup = chromgroup_toanalyze,
			filtergroup = filtergroup_toanalyze,
			extractCallsFile = extractCallsFile %>% basename
		)
	) %>%
	
	#Reorder columns
	relocate(
		c(run_id,filtergroup,extractCallsFile),
		.after = chromgroup
	)

cat("DONE\n")

######################
### Format data for subsequent filters
######################

cat("## Formatting data for subsequent filters...")

#Convert calls tibble back to GRanges
calls.gr <- calls %>%
	makeGRangesFromDataFrame(
		keep.extra.columns = TRUE,
		seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
	)

#Get list of columns by which to join calls to calls.gr (all original columns before filter annotations)
calls.joincols <- extractedCalls %>%
  pluck("calls.gr") %>%
  as_tibble %>%
  colnames

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

rm(extractedCalls)

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
		subsetByOverlaps(genome_chromgroup.gr) %>%
		
		#Add padding
		resize(
			width = width(.) + 2 * (region_genome_filters_config %>% pluck("padding",i)),
			fix="center"
		) %>%
		suppressWarnings %>% #remove warnings of out of bounds regions due to resize
		trim %>%
		GenomicRanges::reduce(ignore.strand=TRUE) %>%
		mutate(passfilter = rep(FALSE,length(.))) #convoluted assignment so it works with empty GRanges
	
	#Subtract genome region filters from filter trackers
	bam.gr.filtertrack <- bam.gr.filtertrack %>%
		GRanges_subtract(region_genome_filter, ignore.strand=TRUE) #GRanges_subtract to retain run_id and zm columns
	
	genome_chromgroup.gr.filtertrack <- genome_chromgroup.gr.filtertrack %>%
		GRanges_subtract(region_genome_filter)
	
	#Annotate filtered calls. Annotates even if partial overlap (relevant for deletions). Annotates final result with a tibble join instead of with the GRanges join due to some calls that overlap > 1 region genome filter range.
	calls <- calls %>%
		left_join(
			calls.gr %>%
				select(run_id,zm) %>%
				join_overlap_inner(region_genome_filter,suffix=c("",".filter")) %>%
				as_tibble %>%
				select(-width) %>%
				rename(!!passfilter_label := passfilter) %>%
				distinct,
			by = join_by(run_id,zm,strand,seqnames,start,end) #ok to join by strand here since calls on both strands were annotated in the join_overlap_inner step
		) %>%
		mutate(across(all_of(passfilter_label), ~ . %>% replace_na(TRUE)))
	
	#Record number of filtered genome bases to stats
	region_genome_filters_stats[i,"num_bases_filtered"] <- region_genome_filter %>%
		width %>%
		sum
}

#Record number of molecule reference space bases remaining after filters
molecule_stats <- molecule_stats %>%
  left_join(
    bam.gr.filtertrack %>%
      as_tibble %>%
      group_by(run_id) %>%
      summarize(num_refspacebases_passgenomeregionfilters = sum(end-start)),
    by = "run_id"
  )

rm(region_genome_filter)

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
	GenomicRanges::reduce(ignore.strand=TRUE) %>%
	mutate(!!passfilter_label := FALSE)

#Subtract genome region filters from filter trackers
bam.gr.filtertrack <- bam.gr.filtertrack %>%
	GRanges_subtract(region_genome_filter, ignore.strand=TRUE)

genome_chromgroup.gr.filtertrack <- genome_chromgroup.gr.filtertrack %>%
	GRanges_subtract(region_genome_filter)

#Annotate filtered calls. Annotates even if partial overlap (relevant for deletions). Performed with a tibble join instead of with GRanges join due to some calls that overlap > 1 region genome filter range.
calls <- calls %>%
	left_join(
		calls.gr %>%
			select(run_id,zm) %>%
			join_overlap_inner(region_genome_filter,suffix=c("",".filter")) %>%
			as_tibble %>%
		  select(-width) %>%
			distinct,
		by = join_by(run_id,zm,strand,seqnames,start,end) #ok to join by strand here since calls on both strands were annotated in the join_overlap_inner step
	) %>%
	mutate(across(all_of(passfilter_label), ~ . %>% replace_na(TRUE)))

#Record number of filtered genome bases to stats
region_genome_filters_stats <- region_genome_filters_stats %>%
  bind_rows(
  	tibble(
  		region_filter_file = "Nbases",
  		binsize = NA,
  		threshold = NA,
  		padding = NA,
  		applyto_chromgroups = "all",
  		applyto_filtergroups = "all",
  		num_bases_filtered = region_genome_filter %>%
  			width %>%
  			sum
  	)
)

#Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
  left_join(
    bam.gr.filtertrack %>%
      as_tibble %>%
      group_by(run_id) %>%
      summarize(num_refspacebases_passNbasesfilter = sum(end-start)),
    by = "run_id"
  )

rm(region_genome_filter)

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
  mutate(min_qual.passfilter =
           if(filtergroup_toanalyze_config$min_qual_method=="mean"){
           	map2_lgl(qual, qual.opposite_strand, ~
           					 	(mean(.x) >= filtergroup_toanalyze_config$min_qual) &
           					 	(mean(.y) >= filtergroup_toanalyze_config$min_qual) &
           					 	!(any(is.na(.y)))
           	)
           	
           }else if(filtergroup_toanalyze_config$min_qual_method=="all"){
           	map2_lgl(qual, qual.opposite_strand, ~
           					 	(all(.x >= filtergroup_toanalyze_config$min_qual)) &
           					 	(all(.y >= filtergroup_toanalyze_config$min_qual)) &
           					 	!(any(is.na(.y)))
           	)
           }else if(filtergroup_toanalyze_config$min_qual_method=="any"){
           	map2_lgl(qual, qual.opposite_strand, ~
           					 	(all(.x >= filtergroup_toanalyze_config$min_qual)) &
           					 	(all(.y >= filtergroup_toanalyze_config$min_qual)) &
           					 	!(any(is.na(.y)))
           	)
           }
  )

#Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
  left_join(
    bam.gr.filtertrack %>%
      as_tibble %>%
      group_by(run_id) %>%
      summarize(num_refspacebases_passminqualfilter = sum(end-start)),
    by = "run_id"
  )

rm(min_qual.fail.gr)

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
  ) %>%
  mutate(read_trim_bp.passfilter = rep(FALSE,length(.)))

rm(bam.gr.onlyranges)

#Subtract bases filtered by read_trim_bp from bam reads filter tracker, joining on run_id and zm. If a position is filtered in one strand, also filter the same reference space position on the opposite strand, since that is how calls are filtered.
bam.gr.filtertrack <- bam.gr.filtertrack %>%
  GRanges_subtract_bymcols(bam.gr.trim, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

#Annotate calls filtered by read_trim_bp without taking strand into account, so that if a call on either strand fails the filter, calls on both strands fail the filter.
calls <- calls %>%
  left_join(
    calls.gr %>%
      select(run_id,zm) %>%
      join_overlap_inner(bam.gr.trim, suffix = c("",".trim")) %>% #strand ignored
      filter(
        run_id == run_id.trim,
        zm == zm.trim
      ) %>%
      as_tibble %>%
      select(-run_id.trim,-zm.trim,-width) %>%
      distinct,
    by = join_by(run_id,zm,strand,seqnames,start,end) #ok to join by strand here since calls on both strands were annotated in the join_overlap_inner step
  ) %>%
  mutate(read_trim_bp.passfilter = read_trim_bp.passfilter %>% replace_na(TRUE))

#Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
  left_join(
    bam.gr.filtertrack %>%
      as_tibble %>%
      group_by(run_id) %>%
      summarize(num_refspacebases_passreadtrimbpfilter = sum(end-start)),
    by = "run_id"
  )

rm(bam.gr.trim)

cat("DONE\n")

######################
### Germline VCF indel region filters
######################
cat("## Applying germline VCF indel region filters...")

#Split by germline_vcf_variants by germline_vcf_file
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
  indel_inspad <- germline_vcf_types_config %>% filter(germline_vcf_type == vcf_type) %>% pull(indel_inspad)
  indel_delpad <- germline_vcf_types_config %>% filter(germline_vcf_type == vcf_type) %>% pull(indel_delpad)
  
  if(!is.na(indel_inspad)){
    indel_inspad <- indel_inspad %>%
      tibble(pad=.) %>%
      extract(pad, into = c("m", "b"), regex = "m(\\d+)b(\\d+)", convert = TRUE)
  }else{
    indel_inspad <- tibble(m = 0, b = 0)
  }
  
  if(!is.na(indel_delpad)){
    indel_delpad <- indel_delpad %>%
      tibble(pad=.) %>%
      extract(pad, into = c("m", "b"), regex = "m(\\d+)b(\\d+)", convert = TRUE)
  }else{
    indel_delpad <- tibble(m = 0, b = 0)
  }
  
  #Create GRanges of variants to filter with configured padding
  germline_vcf_indel_region_filter <- germline_vcf_variants %>%
    pluck(i) %>%
    filter(call_class == "indel") %>%
    mutate(
      padding_m = case_when(
        call_type == "insertion" ~ indel_inspad$m * nchar(alt_plus_strand) %>% round %>% as.integer,
        call_type == "deletion" ~ indel_delpad$m * nchar(ref_plus_strand) %>% round %>% as.integer
      ),
      padding_b = case_when(
        call_type == "insertion" ~ indel_inspad$b,
        call_type == "deletion" ~ indel_delpad$b
      ),
      start = start - pmax(padding_m,padding_b),
      end = end + pmax(padding_m,padding_b)
    ) %>%
    select(-padding_m,-padding_b) %>%
    makeGRangesFromDataFrame(
      seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
    ) %>%
    trim %>%
    GenomicRanges::reduce(ignore.strand=TRUE) %>%
  	mutate(passfilter = FALSE)
  
  #Subtract filtered regions from filter trackers
  bam.gr.filtertrack <- bam.gr.filtertrack %>%
    GRanges_subtract(germline_vcf_indel_region_filter, ignore.strand = TRUE)
  
  genome_chromgroup.gr.filtertrack <- genome_chromgroup.gr.filtertrack %>%
    GRanges_subtract(germline_vcf_indel_region_filter)
  
  #Annotate filtered calls. Annotates even if partial overlap (relevant for deletions). Annotates final result with a tibble join instead of with the GRanges join due to some calls that overlap > 1 region genome filter range.
  calls <- calls %>%
    left_join(
      calls.gr %>%
      	select(run_id,zm) %>%
        join_overlap_inner(germline_vcf_indel_region_filter,suffix=c("",".filter")) %>%
        as_tibble %>%
        select(-width) %>%
      	rename(!!passfilter_label := passfilter) %>%
      	distinct,
      by = join_by(run_id,zm,strand,seqnames,start,end) #ok to join by strand here since calls on both strands were annotated in the join_overlap_inner step
    ) %>%
  	mutate(across(all_of(passfilter_label), ~ . %>% replace_na(TRUE)))
  
  #Record number of filtered genome bases to stats
  region_genome_filters_stats <- region_genome_filters_stats %>%
    bind_rows(
      tibble(
        region_filter_file = i,
        binsize = NA,
        threshold = NA,
        padding = NA,
        applyto_chromgroups = "all",
        applyto_filtergroups = "all",
        num_bases_filtered = germline_vcf_indel_region_filter %>%
          width %>%
          sum
      )
  )
}

#Record number of molecule reference space bases remaining after filters
molecule_stats <- molecule_stats %>%
  left_join(
    bam.gr.filtertrack %>%
      as_tibble %>%
      group_by(run_id) %>%
      summarize(num_refspacebases_passgermlinevcfindelregionfilters = sum(end-start)),
    by = "run_id"
  )

rm(germline_vcf_indel_region_filter)

cat("DONE\n")

######################
### Read indel region filters
######################
cat("## Applying read indel region filters (not applied to indels)...")

#Extract padding configuration
indel_inspad <- filtergroup_toanalyze_config %>% pull(ccsindel_inspad)
indel_delpad <- filtergroup_toanalyze_config %>% pull(ccsindel_delpad)

indel_inspad <- indel_inspad %>%
  tibble(pad=.) %>%
  extract(pad, into = c("m", "b"), regex = "m(\\d+)b(\\d+)", convert = TRUE)

indel_delpad <- indel_delpad %>%
  tibble(pad=.) %>%
  extract(pad, into = c("m", "b"), regex = "m(\\d+)b(\\d+)", convert = TRUE)

#Create GRanges of indels with configured padding
read_indel_region_filter <- calls %>%
  filter(call_class=="indel") %>%
  mutate(
    padding_m = case_when(
      call_type == "insertion" ~ indel_inspad$m * nchar(alt_plus_strand) %>% round %>% as.integer,
      call_type == "deletion" ~ indel_delpad$m * nchar(ref_plus_strand) %>% round %>% as.integer
    ),
    padding_b = case_when(
      call_type == "insertion" ~ indel_inspad$b,
      call_type == "deletion" ~ indel_delpad$b
    ),
    start = start - pmax(padding_m,padding_b),
    end = end + pmax(padding_m,padding_b)
  ) %>%
  select(run_id, zm, seqnames, start, end, strand,-padding_m,-padding_b) %>%
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
  ) %>%
  mutate(read_indel_region_filter.passfilter = rep(FALSE,length(.)))
 
#Annotate calls filtered by read_indel_region_filter without taking strand into account, so that if a call on either strand fails the filter, calls on both strands fail the filter. Not applied to indels
calls <- calls %>%
  left_join(
    calls.gr %>%
      select(run_id,zm) %>%
      join_overlap_inner(read_indel_region_filter, suffix = c("",".filter")) %>% #strand ignored
      filter(
        run_id == run_id.filter,
        zm == zm.filter
      ) %>%
      as_tibble %>%
      select(-run_id.filter,-zm.filter,-width) %>%
      distinct,
    by = join_by(run_id,zm,strand,seqnames,start,end) #ok to join by strand here since calls on both strands were annotated in the join_overlap_inner step
  ) %>%
  mutate(
  	read_indel_region_filter.passfilter = read_indel_region_filter.passfilter %>% replace_na(TRUE),
  	read_indel_region_filter.passfilter = if_else(call_class == "indel", TRUE, read_indel_region_filter.passfilter)
  	)

#Subtract bases below min_qual threshold from bam reads filter tracker, joining on run_id and zm. Set ignore.strand = TRUE so that if a position is filtered in one strand, also filter the same reference space position on the opposite strand, since that is how calls are filtered.
bam.gr.filtertrack <- bam.gr.filtertrack %>%
  GRanges_subtract_bymcols(read_indel_region_filter, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

#Record number of molecule reference space bases remaining after filters
molecule_stats <- molecule_stats %>%
  left_join(
    bam.gr.filtertrack %>%
      as_tibble %>%
      group_by(run_id) %>%
      summarize(num_refspacebases_passreadindelregionfilter = sum(end-start)),
    by = "run_id"
  )

rm(read_indel_region_filter)

cat("DONE\n")

######################
### Germline BAM total reads (coverage) filter
######################
cat("## Applying germline BAM total reads (coverage) filter...")

passfilter_label <- "min_BAMTotalReads.passfilter"

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

tmpbw <- tempfile(tmpdir=getwd(),pattern=".")

system(paste("/bin/bash -c",shQuote(paste(
	yaml.config$wiggletools_bin,"lt",
	filtergroup_toanalyze_config$min_BAMTotalReads,
	germline_bam_samtools_mpileup_file,"|",
	yaml.config$wigToBigWig_bin,"stdin <(cut -f 1,2",
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
	subsetByOverlaps(genome_chromgroup.gr) %>%
	select(-score) %>%
	mutate(!!passfilter_label := rep(FALSE,length(.)))

file.remove(tmpbw) %>% invisible

#Subtract filtered regions from filter trackers
bam.gr.filtertrack <- bam.gr.filtertrack %>%
	GRanges_subtract(germline_bam_samtools_mpileup_filter,ignore.strand = TRUE)

genome_chromgroup.gr.filtertrack <- genome_chromgroup.gr.filtertrack %>%
	GRanges_subtract(germline_bam_samtools_mpileup_filter,ignore.strand = TRUE)

#Annotate calls with germline BAM samtools mpileup bigwig filter. For insertions, left and right flanking bases, and for deletions, the deleted bases, in genome reference space must pass the filter.
calls <- calls %>%
	left_join(
		calls.gr %>%
			mutate(
				startnew = if_else(call_type == "insertion", end, start),
				endnew = if_else(call_type == "insertion", start, end),
				start = startnew,
				end = endnew
			) %>%
			join_overlap_inner(germline_bam_samtools_mpileup_filter) %>%
			as_tibble %>%
			mutate(
				startnew = if_else(call_type == "insertion", end, start),
				endnew = if_else(call_type == "insertion", start, end),
				start = startnew,
				end = endnew,
				width = if_else(call_type == "insertion", 0, width)
			) %>%
			select(all_of(c(calls.joincols,passfilter_label))) %>%
			distinct,
		by = calls.joincols
	) %>%
	mutate(across(all_of(passfilter_label), ~ . %>% replace_na(TRUE)))

#Record number of filtered genome bases to stats
region_genome_filters_stats <- region_genome_filters_stats %>%
	bind_rows(
		tibble(
			region_filter_file = "min_BAMTotalReads",
			binsize = NA,
			threshold = NA,
			padding = NA,
			applyto_chromgroups = "all",
			applyto_filtergroups = "all",
			num_bases_filtered = germline_bam_samtools_mpileup_filter %>%
				width %>%
				sum
		)
	)

#Record number of molecule reference space bases remaining after filters
molecule_stats <- molecule_stats %>%
	left_join(
		bam.gr.filtertrack %>%
			as_tibble %>%
			group_by(run_id) %>%
			summarize(num_refspacebases_passgermlinebamtotalreadsfilter = sum(end-start)),
		by = "run_id"
	)

rm(germline_bam_samtools_mpileup_filter)

cat("DONE\n")

######################
### Germline BAM variant read filters
######################
#Apply only to SBS calls, since germline BAM mpileup calling of indels is too noisy

cat("## Applying germline BAM variant read filters (only applied to SBS calls)...")

#Format variant regions to VCF POS coordinates for extracting variants from germline BAM bcftools mpileup. Retrieve only SBS calls since only annotating SBS calls.
variant_regions_bcftools_mpileup <- calls.gr %>%
	as_tibble %>%
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

#Annotate mpileup with variants passing filters
germline_bam_bcftools_mpileup_filter <- germline_bam_bcftools_mpileup_filter %>%
	mutate(
		max_BAMVariantReads.passfilter = AD2 <= filtergroup_toanalyze_config$max_BAMVariantReads,
		max_BAMVAF.passfilter = VAF <= filtergroup_toanalyze_config$max_BAMVAF
	)

#Annotate SBS calls with germline BAM bcftools mpileup VCF filters. Then set non-SBS calls to pass these filters.
calls <- calls %>%
	left_join(
		calls.gr %>%
			join_overlap_inner(germline_bam_bcftools_mpileup_filter, suffix = c("",".mpileup")) %>%
			as_tibble %>%
			filter(
				call_class %>% as.character == call_class.mpileup %>% as.character,
				call_type %>% as.character == call_type.mpileup %>% as.character,
				SBSindel_call_type %>% as.character == SBSindel_call_type.mpileup %>% as.character,
				ref_plus_strand == ref_plus_strand.mpileup,
				alt_plus_strand == alt_plus_strand.mpileup
			) %>%
			select(all_of(c(calls.joincols,"max_BAMVariantReads.passfilter","max_BAMVAF.passfilter"))) %>%
			distinct,
		by = calls.joincols
	) %>%
	mutate(
		max_BAMVariantReads.passfilter = max_BAMVariantReads.passfilter %>% replace_na(TRUE),
		max_BAMVAF.passfilter = max_BAMVAF.passfilter %>% replace_na(TRUE)
	)

rm(variant_regions_bcftools_mpileup, germline_bam_bcftools_mpileup_filter)

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
bam.gr.filtertrack <- bam.gr.filtertrack %>%
	GRanges_subtract_bymcols(frac_subreads_cvg.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

 #Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
	left_join(
		bam.gr.filtertrack %>%
			as_tibble %>%
			group_by(run_id) %>%
			summarize(num_refspacebases_passminfracsubreadscvgfilter = sum(end-start)),
		by = "run_id"
	)

bam.gr.filtertrack <- bam.gr.filtertrack %>%
	GRanges_subtract_bymcols(num_subreads_match.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

 #Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
	left_join(
		bam.gr.filtertrack %>%
			as_tibble %>%
			group_by(run_id) %>%
			summarize(num_refspacebases_passminnumsubreadsmatchfilter = sum(end-start)),
		by = "run_id"
	)

bam.gr.filtertrack <- bam.gr.filtertrack %>%
	GRanges_subtract_bymcols(frac_subreads_match.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

 #Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
	left_join(
		bam.gr.filtertrack %>%
			as_tibble %>%
			group_by(run_id) %>%
			summarize(num_refspacebases_passminfracsubreadsmatchfilter = sum(end-start)),
		by = "run_id"
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
	mutate(min_frac_subreads_cvg.passfilter =
				 	if(filtergroup_toanalyze_config$min_subreads_cvgmatch_method=="mean"){
				 		map2_lgl(frac_subreads_cvg, frac_subreads_cvg.opposite_strand, ~
				 						 	(mean(.x) >= filtergroup_toanalyze_config$min_frac_subreads_cvg) &
				 						 	(mean(.y) >= filtergroup_toanalyze_config$min_frac_subreads_cvg) &
				 						 	!(any(is.na(.y)))
				 		)
				 	}else if(filtergroup_toanalyze_config$min_subreads_cvgmatch_method=="all"){
				 		map2_lgl(frac_subreads_cvg, frac_subreads_cvg.opposite_strand, ~
				 						 	(all(.x >= filtergroup_toanalyze_config$min_frac_subreads_cvg)) &
				 						 	(all(.y >= filtergroup_toanalyze_config$min_frac_subreads_cvg)) &
				 						 	!(any(is.na(.y)))
				 		)
				 	}else if(filtergroup_toanalyze_config$min_subreads_cvgmatch_method=="any"){
				 		map2_lgl(frac_subreads_cvg, frac_subreads_cvg.opposite_strand, ~
				 						 	(any(.x >= filtergroup_toanalyze_config$min_frac_subreads_cvg)) &
				 						 	(any(.y >= filtergroup_toanalyze_config$min_frac_subreads_cvg)) &
				 						 	!(any(is.na(.y)))
				 		)
				 	}
	) %>%
	select(-c(max_sa,max_sa.opposite_strand,frac_subreads_cvg,frac_subreads_cvg.opposite_strand))


#min_num_subreads_match. NA in any base on the opposite strand is assigned passfilter = FALSE.
calls <- calls %>%
	mutate(min_num_subreads_match.passfilter =
				 	if(filtergroup_toanalyze_config$min_subreads_cvgmatch_method=="mean"){
				 		map2_lgl(sm, sm.opposite_strand, ~
				 						 	(mean(.x) >= filtergroup_toanalyze_config$min_num_subreads_match) &
				 						 	(mean(.y) >= filtergroup_toanalyze_config$min_num_subreads_match) &
				 						 	!(any(is.na(.y)))
				 		)
				 	}else if(filtergroup_toanalyze_config$min_subreads_cvgmatch_method=="all"){
				 		map2_lgl(sm, sm.opposite_strand, ~
				 						 	(all(.x >= filtergroup_toanalyze_config$min_num_subreads_match)) &
				 						 	(all(.y >= filtergroup_toanalyze_config$min_num_subreads_match)) &
				 						 	!(any(is.na(.y)))
				 		)
				 	}else if(filtergroup_toanalyze_config$min_subreads_cvgmatch_method=="any"){
				 		map2_lgl(sm, sm.opposite_strand, ~
				 						 	(any(.x >= filtergroup_toanalyze_config$min_num_subreads_match)) &
				 						 	(any(.y >= filtergroup_toanalyze_config$min_num_subreads_match)) &
				 						 	!(any(is.na(.y)))
				 		)
				 	}
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
	mutate(min_frac_subreads_match.passfilter =
				 	if(filtergroup_toanalyze_config$min_subreads_cvgmatch_method=="mean"){
				 		map2_lgl(frac_subreads_match, frac_subreads_match.opposite_strand, ~
				 						 	((mean(.x) %>% replace_na(0)) >= filtergroup_toanalyze_config$min_frac_subreads_match) &
				 						 	((mean(.y) %>% replace_na(0)) >= filtergroup_toanalyze_config$min_frac_subreads_match) &
				 						 	!(any(is.na(.y)))
				 		)
				 	}else if(filtergroup_toanalyze_config$min_subreads_cvgmatch_method=="all"){
				 		map2_lgl(frac_subreads_match, frac_subreads_match.opposite_strand, ~
				 						 	(all(.x >= filtergroup_toanalyze_config$min_frac_subreads_match) %>% replace_na(0)) &
				 						 	(all(.y >= filtergroup_toanalyze_config$min_frac_subreads_match) %>% replace_na(0)) &
				 						 	!(any(is.na(.y)))
				 		)
				 	}else if(filtergroup_toanalyze_config$min_subreads_cvgmatch_method=="any"){
				 		map2_lgl(frac_subreads_match, frac_subreads_match.opposite_strand, ~
				 						 	(any(.x >= filtergroup_toanalyze_config$min_frac_subreads_match) %>% replace_na(0)) &
				 						 	(any(.y >= filtergroup_toanalyze_config$min_frac_subreads_match) %>% replace_na(0)) &
				 						 	!(any(is.na(.y)))
				 		)
				 	}
	) %>%
	select(-c(frac_subreads_match,frac_subreads_match.opposite_strand))
	
rm(frac_subreads_cvg.fail.gr, num_subreads_match.fail.gr, frac_subreads_match.fail.gr)

cat("DONE\n")


######################
### Duplex coverage filters
######################
cat("## Applying duplex coverage filters...")

#Extract read regions remaining after all filtering that do not have duplex (i.e. both strand) coverage

 #Regions only in plus or minus strand reads
duplex_coverage.fail.gr <- c(
	#Plus strand only
	bam.gr.filtertrack %>%
		filter(strand=="+") %>%
		GRanges_subtract_bymcols(bam.gr.filtertrack %>% filter(strand=="-"), join_mcols = c("run_id","zm"), ignore.strand = TRUE),
	
	#Minus strand only
	bam.gr.filtertrack %>%
		filter(strand=="-") %>%
		GRanges_subtract_bymcols(bam.gr.filtertrack %>% filter(strand=="+"), join_mcols = c("run_id","zm"), ignore.strand = TRUE)
) %>%
	mutate(duplex_coverage.passfilter = rep(FALSE,length(.)))

 #Subtract these regions from the bam.gr filter tracker
bam.gr.filtertrack <- bam.gr.filtertrack %>%
	GRanges_subtract_bymcols(duplex_coverage.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

#Record number of molecule reference space bases remaining after filter
molecule_stats <- molecule_stats %>%
	left_join(
		bam.gr.filtertrack %>%
			as_tibble %>%
			group_by(run_id) %>%
			summarize(num_refspacebases_passduplexcoveragefilter = sum(end-start)),
		by = "run_id"
	)

##Annotate which calls have duplex coverage after all filtering by comparing to final filtered bam.gr filter tracker without taking strand into account, so that if a call on either strand fails the filter, calls on both strands fail the filter.
calls <- calls %>%
	left_join(
		calls.gr %>%
			select(run_id,zm) %>%
			join_overlap_inner(duplex_coverage.fail.gr,suffix = c("",".filter")) %>% #strand ignored
			filter(
				run_id == run_id.filter,
				zm == zm.filter
			) %>%
			as_tibble %>%
			select(-run_id.filter,-zm.filter,-width) %>%
			distinct,
		by = join_by(run_id,zm,strand,seqnames,start,end) #ok to join by strand here since calls on both strands were annotated in the join_overlap_inner step
	) %>%
	mutate(duplex_coverage.passfilter = duplex_coverage.passfilter %>% replace_na(TRUE))

rm(duplex_coverage.fail.gr)

cat("DONE\n")

######################
### Max final calls per molecule filter
######################
cat("## Applying max final calls per molecule filters...")

#Extract which molecules fail the filter for each call_type / SBSindel_call_type being analyzed
max_finalcalls_eachstrand.filter <- calls %>%
	mutate(passallfilters = if_all(contains("passfilter"), ~ . == TRUE)) %>%
	group_by(run_id,zm,strand,call_type,SBSindel_call_type) %>%
	summarize(num_calls = sum(passallfilters), .groups = "drop") %>%
	group_by(run_id,zm,call_type,SBSindel_call_type) %>% #not grouping by strand so entire molecule is filtered if either strand fails
	summarize(
		max_finalcalls_eachstrand.passfilter =
			if(n() == 0){
				TRUE
			}else{
				all(num_calls <= filtergroup_toanalyze_config$max_finalcalls_eachstrand)
			},
		.groups="drop"
		)

#Subtract filtered molecules from bam.gr.filtertrack, creating a new copy of bam.gr.filtertrack for each call_type x SBSindel_call_type combination.
bam.gr.filtertrack <- call_types_toanalyze %>%
	mutate(
		bam.gr.filtertrack.call_type = list(!!bam.gr.filtertrack),
		bam.gr.filtertrack.call_type = pmap(
			list(bam.gr.filtertrack.call_type, call_type, SBSindel_call_type),
			function(x,y,z){
				molecules_to_filter <- !!max_finalcalls_eachstrand.filter %>%
					filter(
						max_finalcalls_eachstrand.passfilter == FALSE,
						call_type == y,
						SBSindel_call_type == z
						)
				
				x %>%
					filter(
						!(run_id %in% molecules_to_filter$run_id &
								zm %in% molecules_to_filter$zm)
						) %>%
					return
			}
		)
	)

#Record number of molecule reference space bases remaining after filter, for each call_type x SBSindel_call_type combination.
molecule_stats <- list(
	pre_max_finalcalls_eachstrand_filter = molecule_stats,
	post_max_finalcalls_eachstrand_filter = bam.gr.filtertrack %>%
		mutate(
			numbases = 
				map(
					bam.gr.filtertrack.call_type,
					function(x){
						x %>%
							as_tibble %>%
							group_by(run_id) %>%
							summarize(num_refspacebases_passmaxfinalcallseachstrandfilter = sum(end-start))
					}
				)
		) %>%
		select(call_type,SBSindel_call_type,numbases) %>%
		unnest_wider(numbases)
)

#Annotate calls with filter
calls <- calls %>%
	left_join(
		max_finalcalls_eachstrand.filter,
		by = join_by(run_id,zm,call_type,SBSindel_call_type)
	) %>%
	mutate(max_finalcalls_eachstrand.passfilter = max_finalcalls_eachstrand.passfilter %>% replace_na(TRUE))

rm(max_finalcalls_eachstrand.filter)

cat("DONE\n")

######################
### Save output data
######################
cat("## Saving output data...")
#Also save total number of bases in chromgroups to analyze

qs_save(
	list(
		config = list(
			yaml.config = yaml.config,
			run_metadata = extractedCalls %>% pluck("run_metadata"),
			extractCallsFile = extractCallsFile,
			individual = individual_id_toanalyze,
			sample_id = sample_id_toanalyze,
			chromgroup = chromgroup_toanalyze,
			genome_chromgroup.gr = genome_chromgroup.gr,
			filtergroup = filtergroup_toanalyze,
			call_types = call_types_toanalyze
		),
		calls = calls,
		bam.gr.filtertrack = bam.gr.filtertrack,
		genome_chromgroup.gr.filtertrack = genome_chromgroup.gr.filtertrack,
		region_genome_filters_stats = region_genome_filters_stats,
		molecule_stats = molecule_stats
	),
	outputFile
)

cat("DONE\n")