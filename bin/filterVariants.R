#!/usr/bin/env -S Rscript --vanilla

#filterVariants.R:
# Filters variants
# Note: for all newly created 'passfilter' columns, TRUE = read/variant successfully passed the filter and should be analyzed.

cat("#### Running filterVariants ####\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
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
							help="path to extractVariants qs2 file"),
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
extractVariantsFile <- opt$file
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
cache_dir <- yaml.config$cache_dir

 #germline vcf filter parameters
germline_vcf_types_config <- yaml.config$germline_vcf_types %>%
  enframe(name=NULL) %>%
  unnest_wider(value)

 #variant types (restrict to selected chromgroup_toanalyze and filtergroup_toanalyze)
variant_types_toanalyze <- yaml.config$variant_types %>%
  enframe(name=NULL) %>%
  unnest_wider(value) %>%
  unnest_longer(SBSindel_call_types) %>%
  unnest_wider(SBSindel_call_types) %>%
	filter(
		analyzein_chromgroups == "all" | (analyzein_chromgroups %>% str_split(",") %>% map(str_trim) %>% map_lgl(~ !!chromgroup_toanalyze %in% .x)),
		filtergroup == filtergroup_toanalyze
		) %>%
	pull(variant_type) %>%
  unique
  
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

#Output data lists

#Display basic configuration parameters
cat("> Processing:\n")
cat("    extractVariants File:",extractVariantsFile,"\n")
cat("    individual_id:",individual_id_toanalyze,"\n")
cat("    sample_id:",sample_id_toanalyze,"\n")
cat("    chromgroup:",chromgroup_toanalyze,"\n")
cat("    filtergroup:",filtergroup_toanalyze,"\n")

cat("DONE\n")


######################
### Define custom functions
######################
# Function to subtract two granges (x and y) from each other, without reducing overlaps in the final output. Modified from code written by Herve Pages.
GRanges_subtract <- function(x, y, ignore.strand=FALSE){
	y <- GenomicRanges::reduce(y, ignore.strand=ignore.strand)
	hits <- findOverlaps(x, y, ignore.strand=ignore.strand)
	ans <- psetdiff(x, extractList(y, as(hits, "IntegerList")))
	unlisted_ans <- unlist(ans, use.names=FALSE)
	mcols(unlisted_ans) <- extractROWS(mcols(x),Rle(seq_along(ans), lengths(ans)))
	unlist(setNames(relist(unlisted_ans, ans), names(x)))
}
	
# Function to subtract two granges (x and y) from each other, if they match on join_mcols (character vector of all columns to match), without reducing overlaps in the final output.
GRanges_subtract_bymcols <- function(x, y, join_mcols, ignore.strand = FALSE) {
  #Encode join keys as a factor
  key_x <- x %>%
    as_tibble %>%
    select(all_of(join_mcols)) %>%
    as.list %>%
    interaction(drop = TRUE)
  
  key_y <- y %>%
    as_tibble %>%
    select(all_of(join_mcols)) %>%
    as.list %>%
    interaction(drop = TRUE)
  
  #Find all overlaps, then pick only those where the metadata‐keys match
  hits <- findOverlaps(x, y, ignore.strand = ignore.strand)
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  sameKey <- as.character(key_x[qh]) == as.character(key_y[sh])
  if (!any(sameKey)) return(x)
  
  #Obtain intersecting segments for matched hits
  qh <- qh[sameKey]
  sh <- sh[sameKey]
  segs_to_remove <- pintersect(x[qh], y[sh], ignore.strand = ignore.strand)
  
  #Split x and segs_to_remove by the index of x,
  x_list  <- split(x, factor(seq_along(x), levels = seq_along(x)), drop = FALSE)
  rem_list <- split(segs_to_remove, factor(qh, levels = seq_along(x)), drop = FALSE)
  
  #Subtract for each range in x (x_list) the ranges in y (rem_list) that intersect it
  out_list <- GenomicRanges::setdiff(x_list, rem_list, ignore.strand = ignore.strand)
  
  #Flatten result and restore all metadata from the original x
   #How many fragments came from each original x[i]
  nfrags <- elementNROWS(out_list)
   #Unlist into a single GRanges
  out <- unlist(out_list, use.names = FALSE)
   #Build an index mapping each fragment back to its source x[i]
  src_idx  <- rep(seq_along(nfrags), times = nfrags)
   #Re‐attach starnd and metadata columns from x[src_idx, ]
  strand(out) <- strand(x)[src_idx]
  mcols(out) <- mcols(x)[src_idx, , drop = FALSE]
  
  return(out)
}

######################
### Load read and extracted variant information
######################
cat("## Loading reads and extracted variants...")

#Load extractedVariants RDS file
extractedVariants <- qs_read(extractVariantsFile)

#Load bam reads as tibble, keeping only reads in selected chroms_toanalyze. Note, this transforms qual from PhredQuality to list.
bam <- extractedVariants %>%
	pluck("bam.gr") %>%
	filter(
		seqnames %in% chroms_toanalyze
	) %>%
	as_tibble

#Load variants as tibble, keeping only variants in selected chroms_toanalyze with variant_type in variant_types_toanalyze (i.e. variant_types in selected chromgroup_toanalyze and filtergroup_toanalyze) or variant_class = SBS or indel (needed for later max calls/mutations postVCF filters and for sensitivity estimation using germline SBS and indel mutations). Mark variants with variant_type in variant_types_toanalyze with a new column variant_type_toanalyze = TRUE.
variants <- extractedVariants %>%
	pluck("variants.gr") %>%
  mutate(variant_type_toanalyze = if_else(variant_type %in% !!variant_types_toanalyze, TRUE, FALSE)) %>%
	filter(
	  seqnames %in% chroms_toanalyze,
	  variant_type_toanalyze == TRUE | variant_class %in% c("SBS","indel")
		) %>%
  as_tibble

#Load prior molecule stats
molecule_stats <- extractedVariants %>%
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
		min_rq_eachstrand.passfilter = if_else(
			all(rq >= filtergroup_toanalyze_config$min_rq_eachstrand),
			TRUE,
			FALSE
		),
		min_rq_avgstrands.passfilter = if_else(
			mean(rq) >= filtergroup_toanalyze_config$min_rq_avgstrands,
			TRUE,
			FALSE
		),
		
		
		#ec
		min_ec_eachstrand.passfilter = if_else(
			all(ec >= filtergroup_toanalyze_config$min_ec_eachstrand),
			TRUE,
			FALSE
		),
		min_ec_avgstrands.passfilter = if_else(
			mean(ec) >= filtergroup_toanalyze_config$min_ec_avgstrands,
			TRUE,
			FALSE
		),
		
		#mapq
		min_mapq_eachstrand.passfilter = if_else(
			all(mapq >= filtergroup_toanalyze_config$min_mapq_eachstrand),
			TRUE,
			FALSE
		),
		min_mapq_avgstrands.passfilter = if_else(
			mean(mapq) >= filtergroup_toanalyze_config$min_mapq_avgstrands,
			TRUE,
			FALSE
		),
		
		#num_softclipbases
		max_num_softclipbases_eachstrand.passfilter = if_else(
			all(num_softclipbases <= filtergroup_toanalyze_config$max_num_softclipbases_eachstrand),
			TRUE,
			FALSE
		),
		max_num_softclipbases_avgstrands.passfilter = if_else(
			mean(num_softclipbases) <= filtergroup_toanalyze_config$max_num_softclipbases_avgstrands,
			TRUE,
			FALSE
		)
		
	) %>%
	select(-num_softclipbases) %>%
	ungroup

#Annotate reads for filters: max_num_SBScalls_eachstrand, max_num_SBScalls_stranddiff, max_num_indelcalls_eachstrand, max_num_indelcalls_stranddiff
bam <- bam %>%
	left_join(
		
		#Input all SBS and indel variants, not just variant_types being analyzed in this run, as the filters being calculated are used for general molecule filters that require info on both SBS and indels
		extractedVariants %>%
			pluck("variants.gr") %>%
			as_tibble %>%
			
			#Count number of SBS and indel calls per strand
			filter(
			  seqnames %in% chroms_toanalyze,
			  variant_class %in% c("SBS","indel")
			  ) %>%
		  
		  #Count number of SBS and indel calls per strand while completing missing strand and SBS/indel values for each molecule so that num_{variant_class}calls are calculated correctly
			mutate(
			  strand = strand %>% factor(levels=c("+","-")),
			  variant_class = variant_class %>% factor(levels=c("SBS","indel"))
			  ) %>% 
			count(run_id,zm,strand,variant_class) %>%
		  complete(nesting(run_id,zm), strand, variant_class, fill = list(n=0)) %>%
			pivot_wider(
				names_from=variant_class,
				values_from=n,
				names_glue = "num_{variant_class}calls"
				) %>%
			
			#Calculate filters for each molecule
			group_by(run_id,zm) %>%
			summarize(
				#max_num_SBScalls_eachstrand.passfilter
				max_num_SBScalls_eachstrand.passfilter = if_else(
					all(num_SBScalls <= filtergroup_toanalyze_config$max_num_SBScalls_eachstrand),
					TRUE,
					FALSE
				),
				
				#max_num_indelcalls_eachstrand.passfilter
				max_num_indelcalls_eachstrand.passfilter = if_else(
					all(num_indelcalls <= filtergroup_toanalyze_config$max_num_indelcalls_eachstrand),
					TRUE,
					FALSE
				),
				
				#max_num_SBScalls_stranddiff.passfilter
				max_num_SBScalls_stranddiff.passfilter = (num_SBScalls * if_else(strand == "+",  1L, -1L)) %>%
					sum %>%
					abs %>%
					#Curly braces ensure the '.' refers to the immediately preceding result
					{if_else(. <= filtergroup_toanalyze_config$max_num_SBScalls_stranddiff, TRUE, FALSE)},
				
				#max_num_indelcalls_stranddiff.passfilter
				max_num_indelcalls_stranddiff.passfilter = (num_indelcalls * if_else(strand == "+",  1L, -1L)) %>%
					sum %>%
					abs %>%
				  {if_else(. <= filtergroup_toanalyze_config$max_num_indelcalls_stranddiff, TRUE, FALSE)}
				
				,.groups = "drop"
			)
		,by = join_by(run_id,zm)
	)
				

#Annotate reads for filters: max_num_SBSmutations, max_num_indelmutations
bam <- bam %>%
	left_join(
		
		#Input all variants, not just variant_types being analyzed, as this is used for general molecule filters that require info on SBS and indels
		extractedVariants %>%
			pluck("variants.gr") %>%
			as_tibble %>%
			
			#Count number of SBS and indel mutations per molecule
			filter(
			  seqnames %in% chroms_toanalyze,
				variant_class %in% c("SBS","indel"),
				SBSindel_call_type == "mutation"
				) %>%
		  
		  #Count number of SBS and indel mutations per molecule while completing missing SBS/indel values for each molecule so that num_{variant_class}mutations are calculated correctly
		  mutate(
		    variant_class = variant_class %>% factor(levels=c("SBS","indel")),
		    SBSindel_call_type = SBSindel_call_type %>% factor(levels="mutation")
		  ) %>% 
		  count(run_id,zm,variant_class,SBSindel_call_type) %>%
		  complete(nesting(run_id,zm), variant_class, SBSindel_call_type, fill = list(n=0)) %>%
			pivot_wider(
				names_from=variant_class,
				values_from=n,
				names_glue = "num_{variant_class}mutations"
			) %>%
			
			#Calculate filters for each molecule
			group_by(run_id,zm) %>%
			summarize(
				#max_num_SBSmutations
				max_num_SBSmutations.passfilter = if_else(
					num_SBSmutations <= filtergroup_toanalyze_config$max_num_SBSmutations,
					TRUE,
					FALSE
				),
				
				#max_num_indelmutations
				max_num_indelmutations.passfilter = if_else(
				  num_indelmutations <= filtergroup_toanalyze_config$max_num_indelmutations,
				  TRUE,
				  FALSE
				)
				
				,.groups = "drop"
			)
		,by = join_by(run_id,zm)
	)

#Replace NA with TRUE for read filters, since some molecules had no SBS or indel variants called.
bam <- bam %>%
  mutate(
    across(
      contains("passfilter"),~ replace_na(.x, TRUE)
    )
  )

#Annotate variants with read filters
variants <- variants %>%
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

#Load germline VCF filter data, and split by germline_vcf_type
germline_vcf_variants <- qs_read(str_c(cache_dir,"/",individual_id_toanalyze,".germline_vcf_variants.qs2")) %>%
  as_tibble %>%
  group_by(germline_vcf_type) %>%
  nest %>%
  deframe
  
#Filter to keep germline VCF variants that pass configured germline VCF filters and keep only columns necessary for downstream filtering
germline_vcf_variants <- germline_vcf_variants  %>%
  imap(
    function(x,idx){
      filters <- germline_vcf_types_config %>% filter(germline_vcf_type == idx)
      
      x %>%
        filter(
          (
            variant_class == "SBS" &
              FILTER %in% (filters$SBS_FILTERS %>% unlist) &
              Depth >= filters$SBS_min_Depth &
              GQ >= filters$SBS_min_GQ &
              VAF >= filters$SBS_min_VAF &
              QUAL >= filters$SBS_min_QUAL
          ) |
          (
            variant_class == "indel" &
              FILTER %in% (filters$indel_FILTERS %>% unlist) &
              Depth >= filters$indel_min_Depth &
              GQ >= filters$indel_min_GQ &
              VAF >= filters$indel_min_VAF &
              QUAL >= filters$indel_min_QUAL
          )
        ) %>%
        select(seqnames,start,end,ref_plus_strand,alt_plus_strand,germline_vcf_file) %>%
        distinct #in case atomization/multi-allelic records created duplicate records
    }
  )

#Annotate calls matching germline VCF variants, and for MDB variants set germline_vcf.passfilter = TRUE regardless if there is a matching germline VCF variant, since we do not want to filter out MDBs just because they are in the location of a germline variant.
variants <- variants %>%
  left_join(
    germline_vcf_variants %>%
      bind_rows(.id="germline_vcf_type") %>%
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
    germline_vcf.passfilter = if_else(variant_class == "MDB", TRUE, germline_vcf.passfilter)
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
    variants %>%
      filter(
        variant_class %in% c("SBS","indel"),
        if_all(contains("passfilter"), ~ .x == TRUE)
      ) %>%
      
      #Count number of SBS and indel calls per strand while completing missing strand and SBS/indel values for each molecule so that num_{variant_class}calls are calculated correctly
      mutate(
        strand = strand %>% factor(levels=c("+","-")),
        variant_class = variant_class %>% factor(levels=c("SBS","indel"))
      ) %>% 
      count(run_id,zm,strand,variant_class) %>%
      complete(nesting(run_id,zm), strand, variant_class, fill = list(n=0)) %>%
      pivot_wider(
        names_from=variant_class,
        values_from=n,
        names_glue = "num_{variant_class}calls"
      ) %>%
      
      #Calculate filters for each molecule
      group_by(run_id,zm) %>%
      summarize(
        #max_num_SBScalls_postVCF_eachstrand
        max_num_SBScalls_postVCF_eachstrand.passfilter = if_else(
          all(num_SBScalls <= filtergroup_toanalyze_config$max_num_SBScalls_postVCF_eachstrand),
          TRUE,
          FALSE
        ),
        
        #max_num_indelcalls_postVCF_eachstrand
        max_num_indelcalls_postVCF_eachstrand.passfilter = if_else(
          all(num_indelcalls <= filtergroup_toanalyze_config$max_num_indelcalls_postVCF_eachstrand),
          TRUE,
          FALSE
        )
        ,.groups = "drop"
      )
    ,by = join_by(run_id,zm)
  )


#Annotate reads for filters: max_num_SBSmutations_postVCF, max_num_indelmutations_postVCF
bam <- bam %>%
  left_join(
    
    #Filter to keep SBS and indel calls that pass all filters applied so far
    variants %>%
      filter(
        variant_class %in% c("SBS","indel"),
        SBSindel_call_type == "mutation",
        if_all(contains("passfilter"), ~ .x == TRUE)
      ) %>%
      
      #Count number of SBS and indel mutations per molecule while completing missing SBS/indel values for each molecule so that num_{variant_class}mutations are calculated correctly
      mutate(
        variant_class = variant_class %>% factor(levels=c("SBS","indel")),
        SBSindel_call_type = SBSindel_call_type %>% factor(levels="mutation")
      ) %>% 
      count(run_id,zm,variant_class,SBSindel_call_type) %>%
      complete(nesting(run_id,zm), variant_class, SBSindel_call_type, fill = list(n=0)) %>%
      pivot_wider(
        names_from=variant_class,
        values_from=n,
        names_glue = "num_{variant_class}mutations"
      ) %>%
      
      #Calculate filters for each molecule
      group_by(run_id,zm) %>%
      summarize(
        #max_num_SBSmutations_postVCF
        max_num_SBSmutations_postVCF.passfilter = if_else(
          num_SBSmutations <= filtergroup_toanalyze_config$max_num_SBSmutations_postVCF,
          TRUE,
          FALSE
        ),
        
        #max_num_indelmutations_postVCF
        max_num_indelmutations_postVCF.passfilter = if_else(
          num_indelmutations <= filtergroup_toanalyze_config$max_num_indelmutations_postVCF,
          TRUE,
          FALSE
        )
        
        ,.groups = "drop"
      )
    ,by = join_by(run_id,zm)
  )

#Replace NA with TRUE for read filters, since some molecules had no SBS or indel variants called.
bam <- bam %>%
  mutate(
    across(
      matches("postVCF.*passfilter"),~ replace_na(.x, TRUE)
    )
  )

#Annotate variants with new read filters
variants <- variants %>%
  left_join(
    bam %>%
      distinct(
        run_id,
        zm,
        pick(matches("postVCF.*passfilter"))
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
  filter_label <- region_read_filters_config %>%
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
    import(format = "bigWig", which = genome_chromgroup.gr) %>%
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
      bam.gr.onlyranges %>%
        
        #Calculate fraction of coverage by the region read filter
        binnedAverage(region_read_filter,varname="frac") %>%
        
        #Annotate which molecules have a fraction of their length filtered by the above bigwig (averaged across plus and minus strands) that is either greater than (gt), greater than or equal (gte), less than (lt), or less than or equal (lte) to a read_threshold_value
        as_tibble %>%
        group_by(run_id, zm) %>%
        summarize(
          !!filter_label := ! case_when(
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

#Annotate variants with new read filters
variants <- variants %>%
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
	
	#Annotate number of molecules passing each individual filter each applied separately
	left_join(
		bam %>%
			group_by(run_id) %>%
			summarize(across(
				matches("passfilter$"),
				~ n_distinct(zm[.x]),
				.names = "num_molecules_{.col}"
			)),
		by = "run_id"
	) %>%
	
	#Annotate number of molecules passing all the molecule filters
	left_join(
		bam %>%
			filter(if_all(matches("passfilter$"), ~ .x == TRUE)) %>%
			group_by(run_id) %>%
			summarize("num_molecules.passallfilters" = n_distinct(zm)),
		by = "run_id"
	)

cat("DONE\n")

######################
### Format data for subsequent filters
######################

cat("## Formatting data for subsequent filters...")

#Convert variants tibble back to GRanges
variants.gr <- variants %>%
	makeGRangesFromDataFrame(
		keep.extra.columns = TRUE,
		seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
	)

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

cat("DONE\n")

######################
### Genome region filters
######################
cat("## Applying genome region filters...")

#Perform filtering
for(i in seq_len(nrow(region_genome_filters_config))){
	
	#Construct label for the newly created filter column
	filter_label <- region_genome_filters_config %>%
		pluck("region_filter_threshold_file",i) %>%
		basename %>%
		str_c("region_genome_filter_",.,".passfilter")
	
	#Load region read filter bigwig
	region_genome_filter <- region_genome_filters_config %>%
		pluck("region_filter_threshold_file",i) %>%
		import(format = "bigWig", which = genome_chromgroup.gr) %>%
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
		GenomicRanges::reduce(ignore.strand=TRUE) %>%
		mutate(!!filter_label := FALSE)
	
	#Subtract genome region filters from filter trackers
	bam.gr.filtertrack <- bam.gr.filtertrack %>%
		GRanges_subtract(region_genome_filter) #GRanges_subtract to retain run_id and zm columns
	
	genome_chromgroup.gr.filtertrack <- genome_chromgroup.gr.filtertrack %>%
		GRanges_subtract(region_genome_filter)
	
	#Annotate filtered variants. Annotates even if partial overlap (relevant for deletions). Annotates final result with a tibble join instead of with the GRanges join due to some variants that overlap > 1 region genome filter range.
	variants <- variants %>%
		left_join(
			variants.gr %>%
				join_overlap_left(region_genome_filter) %>%
				as_tibble %>%
				mutate(across(
					any_of(filter_label),
					~ replace_na(.x,TRUE)
					)) %>%
				distinct,
			by = colnames(.) %>% str_subset("^region_genome_filter", negate=TRUE)
		)
	
	#Record number of filtered bases to stats
	region_genome_filters_stats[i,"num_bases_filtered"] <- region_genome_filter %>%
		width %>%
		sum
}

rm(region_genome_filter)

cat("DONE\n")

######################
### Genome 'N' base sequences
######################
cat("## Applying genome 'N' base sequence filter...")

filter_label <- "region_genome_filter_Nbases.passfilter"

#Extract 'N' base sequence ranges
region_genome_filter <- yaml.config$BSgenome$BSgenome_name %>%
	get %>%
	vmatchPattern("N",.) %>%
	GenomicRanges::reduce(ignore.strand=TRUE) %>%
	mutate(!!filter_label := FALSE)

#Subtract genome region filters from filter trackers
bam.gr.filtertrack <- bam.gr.filtertrack %>%
	GRanges_subtract(region_genome_filter)

genome_chromgroup.gr.filtertrack <- genome_chromgroup.gr.filtertrack %>%
	GRanges_subtract(region_genome_filter)

#Annotate filtered variants. Annotates even if partial overlap (relevant for deletions). Performed with a tibble join instead of with GRanges join due to some variants that overlap > 1 region genome filter range.
variants <- variants %>%
	left_join(
		variants.gr %>%
			join_overlap_left(region_genome_filter) %>%
			as_tibble %>%
			mutate(across(
				any_of(filter_label),
				~ replace_na(.x,TRUE)
			)) %>%
			distinct,
		by = colnames(.) %>% str_subset("^region_genome_filter", negate=TRUE)
	)

#Record number of filtered bases to stats
region_genome_filters_stats <- region_genome_filters_stats %>% bind_rows(
	tibble(
		region_filter_file = "Nbases",
		binsize = 1,
		threshold = "gte0.1",
		padding = 0,
		applyto_chromgroups = "all",
		applyto_filtergroups = "all",
		num_bases_filtered = region_genome_filter %>%
			width %>%
			sum
	)
)

rm(region_genome_filter)

cat("DONE\n")

######################
### Base-quality filters
######################
cat("## Applying base-quality filters...")

#Convert positions of bases failing min_qual filter from query space to reference space, while retaining run_id and zm info
min_qual.fail.gr <- mapFromAlignments( #Main function
  #GRanges of positions failing min_qual filter in query space
  bam %>%
    pull(qual) %>%
    PhredQuality %>%
    as("IntegerList") %>%
    as("RleList") %>%
    IRanges::slice(upper=filtergroup_toanalyze_config$min_qual, includeUpper=FALSE, rangesOnly=TRUE) %>%
    setNames(bam %>% pull(seqnames)) %>%
    as.data.frame %>%
    as_tibble %>%
    makeGRangesFromDataFrame(
      seqnames = "group_name",
      keep.extra.columns=TRUE
    ) %>%
    setNames(.$group) %>%
    select(-group),
  
  #GAlignments of reads in reference space
  bam %>%
    select(seqnames, start, end, cigar) %>% 
    makeGRangesFromDataFrame(
      keep.extra.columns=TRUE,
      seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
    ) %>%
    as("GAlignments") %>%
    setNames(1:nrow(bam))
  ) %>%
  
  #Annotate result with run_id, zm, and strand
  as_tibble %>%
  left_join(
    bam %>%
      select(run_id, zm, strand) %>%
      mutate(row_id = row_number()),
    by = join_by(alignmentsHits == row_id)
  ) %>%
  select(-c(width, xHits, alignmentsHits, strand.x)) %>%
  rename(strand = strand.y) %>%
  makeGRangesFromDataFrame(
    keep.extra.columns=TRUE,
    seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
  )

#Subtract bases below min_qual threshold from bam reads filter tracker, joining on run_id and zm. Set ignore.strand = TRUE so that if a position is filtered in one strand, also filter the same reference space position on the opposite strand, since that is how variants are filtered.
bam.gr.filtertrack <- bam.gr.filtertrack %>%
  GRanges_subtract_bymcols(min_qual.fail.gr, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

#Annotate variants below min_qual threshold. If either strand fails the filter, both strands are set to fail.
if(! filtergroup_toanalyze_config$min_qual_method %in% c("mean","all","any")){
  stop("Incorrect min_qual_method setting!")
}

###FIX THIS CODE TO ALSO ANALYZE QUAL ON OPPOSITE STRAND and take into account possibility of NA in qual of opposite strand (situations where there wasn't a reference-space-matched base on the opposite strand)
variants <- variants %>%
  mutate(min_qual.passfilter =
           if(filtergroup_toanalyze_config$min_qual_method=="mean"){
             qual %>% sapply(mean) >= filtergroup_toanalyze_config$min_qual
           }else if(filtergroup_toanalyze_config$min_qual_method=="all"){
             qual %>% sapply(function(x){all(x >= filtergroup_toanalyze_config$min_qual)})
           }else if(filtergroup_toanalyze_config$min_qual_method=="any"){
             qual %>% sapply(function(x){any(x >= filtergroup_toanalyze_config$min_qual)})
           }
  )

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
  mutate(read_trim_bp.passfilter = FALSE)

rm(bam.gr.onlyranges)

#Subtract bases filtered by read_trim_bp from bam reads filter tracker, joining on run_id and zm. If a position is filtered in one strand, also filter the same reference space position on the opposite strand, since that is how variants are filtered.
bam.gr.filtertrack <- bam.gr.filtertrack %>%
  GRanges_subtract_bymcols(bam.gr.trim, join_mcols = c("run_id","zm"), ignore.strand = TRUE)

#Annotate variants filtered by read_trim_bp. If either strand fails the filter, both strands are set to fail.
variants <- variants %>%
  left_join(
    variants.gr %>%
      join_overlap_left(bam.gr.trim) %>%
      as_tibble %>%
      mutate(read_trim_bp.passfilter = read_trim_bp.passfilter ~ replace_na(.x,TRUE)) %>%
      distinct,
    by = colnames(.) %>% str_subset("^region_genome_filter", negate=TRUE)
  ) %>%
  group_by(seqnames,start,end,run_id) %>%
  mutate(read_trim_bp.passfilter = if_else(all(read_trim_bp.passfilter), TRUE, FALSE)) %>%
  ungroup

cat("DONE\n")

######################
### Germline VCF indel region filters
######################
cat("## Applying germline VCF indel region filters...")

cat("DONE\n")

######################
### Read indel region filters
######################
cat("## Applying read indel region filters...")

cat("DONE\n")

######################
### Germline BAM filter
######################
cat("## Applying germline BAM filter...")

cat("DONE\n")

######################
### Subread filters
######################
cat("## Applying subread filters...")

#Take into account possibility of values of 0 in the denominator, which will lead to value of Inf or 0/0 = NaN, etc.

cat("DONE\n")

######################
### Duplex coverage filters
######################
cat("## Applying duplex coverage filters...")

##Annotate which variants have duplex coverage after all filtering by comparing to final filtered bam.gr filter tracker.

cat("DONE\n")

######################
### Save output data
######################

#Also save total number of bases in chromgroups to analyze