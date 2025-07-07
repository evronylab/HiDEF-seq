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

region_bin_filters_config <- yaml.config$region_filters %>%
  map("bin_filters") %>%
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
### Load read and extracted variant information
######################
cat("## Loading reads and extracted variants...")

#Load data
extractedVariants <- qs_read(extractVariantsFile)

#Keep only reads in selected chroms_toanalyze
bam.gr.filtered <- extractedVariants %>%
	pluck("bam.gr") %>%
	filter(
		seqnames %in% chroms_toanalyze
	)

#Keep only variants in selected chroms_toanalyze with variant_type in variant_types_toanalyze (i.e. variant_types in selected chromgroup_toanalyze and filtergroup_toanalyze) or variant_class = SBS or indel (needed for later max calls/mutations postVCF filters and for sensitivity estimation using germline SBS and indel mutations). Mark variants with variant_type in variant_types_toanalyze with a new column variant_type_toanalyze = TRUE.
variants.gr.filtered <- extractedVariants %>%
	pluck("variants.gr") %>%
  mutate(variant_type_toanalyze = if_else(variant_type %in% !!variant_types_toanalyze, TRUE, FALSE)) %>%
	filter(
	  seqnames %in% chroms_toanalyze,
	  variant_type_toanalyze == TRUE | variant_class %in% c("SBS","indel")
		) 

#Convert reads and variants to tibbles since much faster than plyranges analysis for initial filters
bam.gr.filtered <- bam.gr.filtered %>%
  as_tibble

variants.gr.filtered <- variants.gr.filtered %>%
  as_tibble

#Create GRanges of selected chroms_toanalyze in the genome to track genome filtering
genome_chromgroup.gr.filtered <- yaml.config$BSgenome$BSgenome_name %>%
	get %>%
	seqinfo %>%
	GRanges %>%
	filter(seqnames %in% !!chroms_toanalyze)

#Molecule stats
molecule_stats <- extractedVariants %>%
	pluck("molecule_stats")

cat("DONE\n")

######################
### Basic whole-molecule filters
######################
cat("## Applying basic whole-molecule filters...")

#Annotate reads for filters: rq, ec, mapq, num_softclipbases
bam.gr.filtered <- bam.gr.filtered %>%
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
bam.gr.filtered <- bam.gr.filtered %>%
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
bam.gr.filtered <- bam.gr.filtered %>%
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
bam.gr.filtered <- bam.gr.filtered %>%
  mutate(
    across(
      contains("passfilter"),~ replace_na(.x, TRUE)
    )
  )

#Annotate variants with read filters
variants.gr.filtered <- variants.gr.filtered %>%
  left_join(
    bam.gr.filtered %>%
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
variants.gr.filtered <- variants.gr.filtered %>%
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
### Post-germline VCF variant filtering whole-molecule filters
######################
cat("## Applying post-germline VCF variant filtering whole-molecule filters...")

#Annotate reads for filters: max_num_SBScalls_postVCF_eachstrand, max_num_indelcalls_postVCF_eachstrand
bam.gr.filtered <- bam.gr.filtered %>%
  left_join(
    
    #Filter to keep SBS and indel calls that pass all filters applied so far
    variants.gr.filtered %>%
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
bam.gr.filtered <- bam.gr.filtered %>%
  left_join(
    
    #Filter to keep SBS and indel calls that pass all filters applied so far
    variants.gr.filtered %>%
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
bam.gr.filtered <- bam.gr.filtered %>%
  mutate(
    across(
      matches("postVCF.*passfilter"),~ replace_na(.x, TRUE)
    )
  )

#Annotate variants with new read filters
variants.gr.filtered <- variants.gr.filtered %>%
  left_join(
    bam.gr.filtered %>%
      distinct(
        run_id,
        zm,
        pick(matches("postVCF.*passfilter"))
      )
    ,by = join_by(run_id,zm)
  )

cat("DONE\n")

######################
### Region filters
######################
cat("## Applying region filters...")

#Also noN filter

cat("DONE\n")

######################
### Base-quality filters
######################
cat("## Applying base-quality filters...")

cat("DONE\n")

######################
### Molecule position filters
######################
cat("## Applying molecule position filters...")

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

cat("DONE\n")

#Convert back to GRanges
%>%
  makeGRangesFromDataFrame(
    keep.extra.columns=TRUE,
    seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
  )