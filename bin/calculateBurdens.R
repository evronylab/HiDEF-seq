#!/usr/bin/env -S Rscript --vanilla

#calculateBurdens.R:
# Calculate burdens.
	#Calculates non-sensitivity-corrected burdens and the sensitivity values (if sensitivity analysis is configured). The calculated sensitivity values are then used to calculate sensitivity-corrected burdens in the outputResults script rather than in this script, since this script operates per chromgroup x filtergroup whereas the sensitivity correction uses for each filtergroup the sensitivity value of one chromgroup (specified in sensitivity_parameters$use_chromgroup) for sensitivity-correction of all chromgroups.

cat("#### Running calculateBurdens ####\n")

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
	make_option(c("-g", "--chromgroup_toanalyze"), type = "character", default=NULL,
	            help="chromgroup to analyze"),
	make_option(c("-v", "--filtergroup_toanalyze"), type = "character", default=NULL,
	            help="filtergroup to analyze"),
	make_option(c("-f", "--files"), type = "character", default=NULL,
							help="comma-separated filterCalls qs2 files"),
		make_option(c("-o", "--output"), type = "character", default=NULL,
							help="output qs2 file")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$config) | is.na(opt$sample_id_toanalyze) | is.na(opt$chromgroup_toanalyze) | is.na(opt$filtergroup_toanalyze) | is.na(opt$files) | is.na(opt$output) ){
	stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))
sample_id_toanalyze <- opt$sample_id_toanalyze
chromgroup_toanalyze <- opt$chromgroup_toanalyze
filtergroup_toanalyze <- opt$filtergroup_toanalyze
filterCallsFiles <- opt$files %>% str_split_1(",") %>% str_trim
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

 #sex chromosomes
sex_chromosomes <- yaml.config$sex_chromosomes %>%
	str_split_1(",") %>%
	str_trim

#mitochondrial chromosome
mitochondrial_chromosome <- yaml.config$mitochondrial_chromosome

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

 #sensitivity thresholds
sensitivity_parameters <- yaml.config$sensitivity_parameters %>%
	as_tibble

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
### Load custom shared functions
######################
source(Sys.which("sharedFunctions.R"))

######################
### Define custom functions
######################
#Order of sbs trinucleotide context labels
sbs96_labels <- c(
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

sbs192_labels <- c(
	sbs96_labels,
	str_c(
		str_sub(sbs96_labels,1,3) %>% DNAStringSet %>% reverseComplement %>% as.character,
		">",
		str_sub(sbs96_labels,5,7) %>% DNAStringSet %>% reverseComplement %>% as.character
	)
)

sbs192_labels.sigfit <- c(
	str_c("T:",sbs96_labels),
	str_c("U:",sbs96_labels)
)

#Order of indel context labels
indel_labels.sigfit <- c(
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
	"5:Del:M:1","5:Del:M:2","5:Del:M:3","5:Del:M:4","5:Del:M:5"
)

#Function to sum two RleLists, including handling of chromosomes only present in one of the RleLists.
sum_RleList <- function(a, b) {
	seqs_union <- union(names(a), names(b))
	
	seqs_union %>%
		map(function(nm){
			seq_in_a <- nm %in% names(a)
			seq_in_b <- nm %in% names(b)
			if(seq_in_a && seq_in_b){
				a[[nm]] + b[[nm]]
			}else if(seq_in_a){
				a[[nm]]
			}else{
				b[[nm]]
			}
		}) %>%
			set_names(seqs_union) %>%
			RleList(compress=FALSE)
}

#Function to calculate duplex genome coverage for bam.gr.filtertrack_bytype. If two bam.gr.filtertrack_bytypes are provided, it calculates duplex genome coverage only for the second and adds it to the first. Also removes the bam.gr.filtertrack column that is no longer necessary.
sum_bam.gr.filtertracks <- function(bam.gr.filtertrack1, bam.gr.filtertrack2=NULL){

	if(is.null(bam.gr.filtertrack2)){
		bam.gr.filtertrack1 %>%
			mutate(
				bam.gr.filtertrack.coverage = bam.gr.filtertrack %>%
					map(function(x){GenomicRanges::reduce(x,ignore.strand=TRUE) %>% coverage}) #Collapse separate + and - reads
			) %>%
			select(-bam.gr.filtertrack)
		
	}else{
		bam.gr.filtertrack1 %>%
			left_join(
				bam.gr.filtertrack2 %>%
					mutate(
						bam.gr.filtertrack.coverage = bam.gr.filtertrack %>%
							map(function(x){GenomicRanges::reduce(x,ignore.strand=TRUE) %>% coverage}) #Collapse separate + and - reads
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

#Function to extract coverage for a GRanges object with only 1 bp ranges from a SimpleRleList coverage object. Coverage = 0 for seqnames in the GRanges that are not in the coverage object.
gr_1bp_cov <- function(gr, cov){
	
	stopifnot(all(width(gr) == 1))
	
	n <- length(gr)
	result <- rep.int(0, n) #pre-fill zeros
	chr <- gr %>% seqnames %>% as.character
	pos <- gr %>% start
	
	# indices per chromosome (keeps order in 'gr')
	idx_by_chr <- split(seq_len(n), chr)
	
	for (ch in names(idx_by_chr)) {
		idx <- idx_by_chr[[ch]]
		if (ch %in% names(cov)){
			result[idx] <- cov[[ch]][pos[idx]]
		}
	}
	
	return(result)
}

#Function to re-anchor indels from HiDEF-seq calls to the style of VCF format so that ref_plus_strand and alt_plus_strand contain the base preceding the indel
normalize_indels_for_vcf <- function(df, BSgenome_name) {
	
	#Identify deletions and insertions
	# insertions: ref_plus_strand == ""; alt_plus_strand = inserted bases
	# deletions: ref_plus_strand = deleted bases; alt_plus_strand == ""
	ins <- !nzchar(df$ref_plus_strand) & nzchar(df$alt_plus_strand)
	del <- nzchar(df$ref_plus_strand) & !nzchar(df$alt_plus_strand)
	ins_or_del <- ins | del
	
	#Normalize start position of insertions and deletions to start = start - 1, and extract sequence of the new start position
	if(any(ins_or_del)){
		df$start[ins_or_del] <- df$start[ins_or_del] - 1
		
		gr <- GRanges(
			df$seqnames[ins_or_del],
			IRanges(df$start[ins_or_del],	width = 1),
			seqinfo = BSgenome_name %>% get %>% seqinfo
		)
		
		start_base_ref <- df %>% nrow %>% character
		start_base_ref[ins_or_del] <- getSeq(BSgenome_name %>% get, gr) %>% as.character
	}
	
	#Normalize insertion sequences: ref_plus_strand = new start position base; alt_plus_strand = new start position base + inserted bases
	if(any(ins)){
		df$ref_plus_strand[ins] <- start_base_ref[ins]
		df$alt_plus_strand[ins] <- str_c(start_base_ref[ins], df$alt_plus_strand[ins])
	}
	
	#Normalize deletion sequences: ref_plus_strand = new start position base + deleted bases; alt_plus_strand = new start position base
	if(any(del)){
		df$ref_plus_strand[del] <- str_c(start_base_ref[del], df$ref_plus_strand[del])
		df$alt_plus_strand[del] <- start_base_ref[del]
	}
	
	#Remove 'end' column that is not needed for vcf, to avoid confusion
	df$end <- NULL
	
	return(df)
}

######################
### Load data from filterCalls files
######################
cat("## Loading data from filterCalls files...\n> analysis chunk:")

#Create lists for data loading
finalCalls <- list()
germlineVariantCalls <- list()
molecule_stats <- list()

#Loop over all filterCallsFiles
for(i in seq_along(filterCallsFiles)){
	
	#Extract analysis chunk number for logging in case filterCallsFiles are not supplied in order
	analysis_chunk_num <- filterCallsFiles[i] %>% str_extract("(?<=chunk)\\d+")
		
	cat(" ", analysis_chunk_num, sep="")
	
	#Load filterCallsFile
	filterCallsFile <- qs_read(filterCallsFiles[i])
	
	#Load basic configuration only from first chunk, since identical in all chunks.
	if(i == 1){
		#Basic configuration parameters
		run_metadata <- filterCallsFile %>% pluck("config","run_metadata")
		genome_chromgroup.gr <- filterCallsFile %>% pluck("config","genome_chromgroup.gr")
		region_read_filters_config <- filterCallsFile %>% pluck("config","region_read_filters_config")
		region_genome_filters_config <- filterCallsFile %>% pluck("config","region_genome_filters_config")
		
		#germline-related filters that are excluded when calculating sensitivity
		germline_filters <- c(
			"germline_vcf.passfilter",
			str_subset(filterCallsFile %>% pluck("calls") %>% colnames, "germline_vcf_indel_region_filter"),
			"max_BAMVariantReads.passfilter",
			"max_BAMVAF.passfilter",
			region_read_filters_config %>%
				filter(is_germline_filter==TRUE) %>%
				pull(region_filter_threshold_file) %>%
				basename %>%
				str_c("region_read_filter_",.,".passfilter"),
			region_genome_filters_config %>%
				filter(is_germline_filter==TRUE) %>%
				pull(region_filter_threshold_file) %>%
				basename %>%
				str_c("region_genome_filter_",.,".passfilter")
		)
		
		#region_genome_filter_stats
		region_genome_filter_stats <- filterCallsFile %>% pluck("region_genome_filter_stats")
	}

	#Load final calls that pass all filters, ignoring the germline filters, since we will later also need the calls filtered by germline filters to calculate sensitivity
	finalCalls[[i]] <- filterCallsFile %>%
		pluck("calls") %>%
		filter(
			call_toanalyze == TRUE,
			if_all(contains("passfilter"), ~ .x == TRUE)
		)
	
	germlineVariantCalls[[i]] <- filterCallsFile %>%
		pluck("calls") %>%
		filter(
			call_class %in% c("SBS","indel"),
			SBSindel_call_type == "mutation",
			if_all(contains("passfilter") & !all_of(germline_filters), ~ .x == TRUE), #All non-germline filters == TRUE
			germline_vcf.passfilter == FALSE #Detected in at least one germline VCF
		)
	
	#Filtered duplex read coverage of the genome, with and without germline_filters
	if(i == 1){
		bam.gr.filtertrack.bytype <- sum_bam.gr.filtertracks(
			filterCallsFile %>% pluck("bam.gr.filtertrack.bytype")
		)
		
		bam.gr.filtertrack.except_germline_filters.bytype <- sum_bam.gr.filtertracks(
			filterCallsFile %>% pluck("bam.gr.filtertrack.except_germline_filters.bytype")
		)
		
	}else{
		bam.gr.filtertrack.bytype <- sum_bam.gr.filtertracks(
			bam.gr.filtertrack.bytype,
			filterCallsFile %>% pluck("bam.gr.filtertrack.bytype")
		)
		
		bam.gr.filtertrack.except_germline_filters.bytype <- sum_bam.gr.filtertracks(
			bam.gr.filtertrack.except_germline_filters.bytype,
			filterCallsFile %>% pluck("bam.gr.filtertrack.except_germline_filters.bytype")
		)
	}

	#molecule_stats
	molecule_stats[[i]] <- filterCallsFile %>% pluck("molecule_stats")
}

#Remove temp objects
rm(filterCallsFile)
invisible(gc())

#Combine finalCalls and germlineVariantCalls each to one tibble
finalCalls <- finalCalls %>%
	bind_rows(.id = "analysis_chunk") %>%
	mutate(analysis_chunk = analysis_chunk %>% as.integer) %>%
	select(-call_toanalyze)

germlineVariantCalls <- germlineVariantCalls %>%
	bind_rows(.id = "analysis_chunk") %>%
	mutate(analysis_chunk = analysis_chunk %>% as.integer) %>%
	select(-call_toanalyze)

#Combine molecule_stats to one tibble and sum stats. Also do a left_join to molecule_stats[[1]] without the value column to preserve the original row order of stats. Also calculate stats summed across all run_ids.
molecule_stats.by_run_id <- molecule_stats[[1]] %>%
	select(run_id,chromgroup,filtergroup,stat) %>%
	left_join(
		molecule_stats %>%
			bind_rows %>%
			group_by(run_id,chromgroup,filtergroup,stat) %>%
			summarize(value = sum(value), .groups = "drop"),
		by = join_by(run_id,chromgroup,filtergroup,stat)
	)

molecule_stats.by_analysis_id <- molecule_stats[[1]] %>%
	select(chromgroup,filtergroup,stat) %>%
	left_join(
		molecule_stats %>%
			bind_rows %>%
			group_by(chromgroup,filtergroup,stat) %>%
			summarize(value = sum(value), .groups = "drop"),
		by = join_by(chromgroup,filtergroup,stat)
	)

rm(molecule_stats)
invisible(gc())

cat(" DONE\n")

######################
### Calculate coverage and extract reference sequences of interrogated genome bases
######################
cat("## Calculating coverage and extracting reference sequences of interrogated genome bases...")

#Convert duplex coverage RleList of final interrogated genome bases to GRanges for every base with a 'coverage' column, excluding zero coverage bases
bam.gr.filtertrack.bytype <- bam.gr.filtertrack.bytype %>%
	mutate(
		bam.gr.filtertrack.coverage = bam.gr.filtertrack.coverage %>%
			map(
				function(x){
					if(length(x) == 0){
						return(GRanges(
							duplex_coverage=integer(),
							seqinfo = yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
							))
						}
					
					x %>%
						imap(
							function(r,chr){
								pos <- which(r != 0L) 
								
								if(length(pos) == 0L){
									return(GRanges(
										duplex_coverage=integer(),
										seqinfo = yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
										))
									}
								
								rl <- runLength(r)
								rv <- runValue(r)
								nz <- rv != 0L
								
								GRanges(
									seqnames = chr,
									ranges = IRanges(start = pos, width = 1L),
									duplex_coverage = rep.int(rv[nz], rl[nz]),
									seqinfo = yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
								)
							}
						) %>%
						GRangesList %>%
						unlist(use.names = FALSE)
				}
			)
	)

#Annotate trinucleotide reference sequences for final interrogated genome bases

 #Extract all regions for which to obtain sequence
regions_to_getseq <- bam.gr.filtertrack.bytype %>%
	pull(bam.gr.filtertrack.coverage) %>%
	GRangesList %>%
	unlist(use.names = FALSE) %>%
	select(-duplex_coverage) %>%
	sort %>%
	unique %>%
	resize(width = 3, fix="center") %>%
	trim %>% #Remove trinucleotide contexts that extend past a non-circular chromosome edge (after trim, width < 3)
	filter(width == 3)

 #Extract reference sequences with seqkit (faster than R getSeq)
tmpregions <- tempfile(tmpdir=getwd(),pattern=".")
tmpseqs <- tempfile(tmpdir=getwd(),pattern=".")

str_c(
	regions_to_getseq %>% seqnames %>% as.character,
	":", regions_to_getseq %>% start,
	"-", regions_to_getseq %>% end
) %>%
	write_lines(tmpregions)

invisible(system(paste(
	yaml.config$seqkit_bin,"faidx --quiet -l",tmpregions,yaml.config$genome_fasta,"|",
	yaml.config$seqkit_bin,"seq -u |", #convert to upper case
	yaml.config$seqkit_bin,"fx2tab -Q |",
	"sed -E 's/:([0-9]+)-([0-9]+)\\t/\\t\\1\\t\\2\\t/' >", #change : and - to tab
	tmpseqs
),intern=FALSE))

seqkit_seqs <- tmpseqs %>%
	read_tsv(
		col_names=c("seqnames","start","end","reftnc_plus_strand"),
		col_types = cols(
			seqnames = col_factor(),
			start = col_integer(),
			end = col_integer(),
			reftnc_plus_strand = col_factor(levels = trinucleotides_64)
		),
	) %>%
	makeGRangesFromDataFrame(
		keep.extra.columns = TRUE,
		seqinfo = yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
	)

invisible(file.remove(tmpregions,tmpseqs))

hits <- findOverlaps(regions_to_getseq, seqkit_seqs, type = "equal")
regions_to_getseq$reftnc_plus_strand <- factor(NA_character_, levels = trinucleotides_64)
regions_to_getseq$reftnc_plus_strand[queryHits(hits)] <- seqkit_seqs$reftnc_plus_strand[subjectHits(hits)]

rm(seqkit_seqs, hits)
invisible(gc())

regions_to_getseq <- regions_to_getseq %>%
	resize(width = 1, fix = "center")

 #Join extracted sequences back to each bam.gr.filtertrack.coverage, and use that to also annotate reftnc_minus_strand and reftnc_pyr. 'for' loop uses less memory than 'map'
for(i in seq_len(nrow(bam.gr.filtertrack.bytype))){
	h <- findOverlaps(
		bam.gr.filtertrack.bytype$bam.gr.filtertrack.coverage[[i]],
		regions_to_getseq,
		type = "equal"
	)
	
	mcols(bam.gr.filtertrack.bytype$bam.gr.filtertrack.coverage[[i]])$reftnc_plus_strand <- factor(NA_character_, levels = trinucleotides_64)
	
	mcols(bam.gr.filtertrack.bytype$bam.gr.filtertrack.coverage[[i]])$reftnc_plus_strand[queryHits(h)] <- regions_to_getseq$reftnc_plus_strand[subjectHits(h)]
	
	bam.gr.filtertrack.bytype$bam.gr.filtertrack.coverage[[i]] <- bam.gr.filtertrack.bytype$bam.gr.filtertrack.coverage[[i]] %>%
		mutate(
			reftnc_minus_strand = reftnc_plus_strand %>%
				fct_relabel(
					function(x){
						x %>% DNAStringSet %>% reverseComplement %>% as.character
					}
				) %>%
				fct_relevel(trinucleotides_64),
			
			reftnc_pyr = reftnc_plus_strand %>%
				fct_collapse(!!!trinucleotides_64_32_pyr_list) %>%
				factor(levels = trinucleotides_32_pyr)
		)
}

rm(regions_to_getseq, h)
invisible(gc())

######################
### Calculate trinucleotide distributions of interrogated bases and the genome
######################
cat("## Calculating trinucleotide distributions of interrogated bases and the genome...")

#Calculate trinucleotide distribution of final interrogated genome bases
bam.gr.filtertrack.bytype <- bam.gr.filtertrack.bytype %>%
	mutate(
		bam.gr.filtertrack.reftnc_pyr = bam.gr.filtertrack.coverage %>%
			map(
				function(x){
					x %>%
						as_tibble %>%
						count(reftnc_pyr, wt = duplex_coverage, name = "count") %>%
						complete(reftnc_pyr, fill = list(count = 0)) %>%
						arrange(reftnc_pyr) %>%
						filter(!is.na(reftnc_pyr)) %>%
						mutate(fraction = count / sum(count)) %>%
						select(-count)
				}
			),
		
		bam.gr.filtertrack.reftnc_both_strands = bam.gr.filtertrack.coverage %>%
			map(
				function(x){
					bind_rows(
						x %>%
							as_tibble %>%
							count(reftnc_plus_strand, wt = duplex_coverage, name = "count") %>%
							rename(reftnc = reftnc_plus_strand),
						x %>%
							as_tibble %>%
							count(reftnc_minus_strand, wt = duplex_coverage, name = "count") %>%
							rename(reftnc = reftnc_minus_strand)
					) %>%
						group_by(reftnc) %>%
						summarize(count = sum(count), .groups = "drop") %>%
						complete(reftnc, fill = list(count = 0)) %>%
						arrange(reftnc) %>%
						filter(!is.na(reftnc)) %>%
						mutate(fraction = count / sum(count)) %>%
						select(-count)
				}
			),
		
		#Remove reftnc_plus_strand and reftnc_minus_strand that are no longer needed
		bam.gr.filtertrack.coverage = bam.gr.filtertrack.coverage %>%
			map(
				function(x){
					x %>% select(-reftnc_plus_strand,-reftnc_minus_strand)
				}
			)
	)

#Calculate trinucleotide distributions of the whole genome and of the genome in the analyzed chromgroup
 #Function to extract trinucleotide distribution for selected chromosomes
get_genome_reftnc <- function(BSgenome_name, chroms){
	
	reftnc_plus_strand <- BSgenome_name %>%
		get %>%
		getSeq(chroms) %>%
		trinucleotideFrequency(simplify.as = "collapsed") %>%
		enframe(name = "reftnc", value = "count") %>%
		mutate(reftnc = reftnc %>% factor(levels = trinucleotides_64))
	
	reftnc_minus_strand <- reftnc_plus_strand %>%
		mutate(
			reftnc = reftnc %>%
				DNAStringSet %>%
				reverseComplement %>%
				as.character %>%
				factor(levels = trinucleotides_64)
		)

	reftnc_pyr <- reftnc_plus_strand %>%
		trinucleotides_64to32(tri_column = "reftnc", count_column = "count") %>%
		rename(reftnc_pyr = reftnc) %>%
		filter(!is.na(reftnc_pyr)) %>%
		mutate(fraction = count / sum(count)) %>%
		select(-count)
	
	reftnc_both_strands <- bind_rows(
		reftnc_plus_strand,
		reftnc_minus_strand
	) %>%
		group_by(reftnc) %>%
		summarize(count = sum(count), .groups = "drop") %>%
		arrange(reftnc) %>%
		filter(!is.na(reftnc)) %>%
		mutate(fraction = count / sum(count)) %>%
		select(-count)
	
	return(
		list(
			reftnc_pyr = reftnc_pyr,
			reftnc_both_strands = reftnc_both_strands
		)
	)
}

	#Whole genome
genome.reftnc <- get_genome_reftnc(
	BSgenome_name = yaml.config$BSgenome$BSgenome_name,
	chroms = yaml.config$BSgenome$BSgenome_name %>% get %>% seqnames
	)

 #genome in the analyzed chromgroup
genome_chromgroup.reftnc <- get_genome_reftnc(
	BSgenome_name = yaml.config$BSgenome$BSgenome_name,
	chroms = chroms_toanalyze
)

cat("DONE\n")

#Calculate ratio of trinucleotide fractions of final interrogated genome bases to the trinucleotide fractions of the whole genome and the genome in the analyzed chromgroup
bam.gr.filtertrack.bytype <- bam.gr.filtertrack.bytype %>%
	mutate(
		bam.gr.filtertrack.reftnc_pyr = bam.gr.filtertrack.reftnc_pyr %>%
			map(
				function(x){
					x %>%
						left_join(
							genome.reftnc$reftnc_pyr,
							by = "reftnc_pyr",
							suffix = c("",".genome")
						) %>%
						left_join(
							genome_chromgroup.reftnc$reftnc_pyr,
							by = "reftnc_pyr",
							suffix = c("",".genome_chromgroup")
						)
				}
			),
		
		bam.gr.filtertrack.reftnc_both_strands = bam.gr.filtertrack.reftnc_both_strands %>%
			map(
				function(x){
					x %>%
						left_join(
							genome.reftnc$reftnc_both_strands,
							by = "reftnc",
							suffix = c("",".genome")
						) %>%
						left_join(
							genome_chromgroup.reftnc$reftnc_both_strands,
							by = "reftnc",
							suffix = c("",".genome_chromgroup")
						)
				}
			),
		
		across(
			starts_with("bam.gr.filtertrack.reftnc"),
			function(x){
				x %>%
					map(
						function(y){
							y %>% 
								mutate(
									fraction_ratio_to_genome = fraction / fraction.genome,
									fraction_ratio_to_genome_chromgroup = fraction / fraction.genome_chromgroup
								) %>%
								select(-contains(".genome"))
						}
					)
			}
		)
	)

######################
### Calculate trinucleotide distributions (SBS and MDB calls) and spectra (SBS and indel calls)
######################
cat("## Calculating trinucleotide disributions (SBS and MDB calls) and spectra (SBS and indel calls)...")
	
#Nest_join finalCalls for each call_class x call_type x SBSindel_call_type combination and collapse to distinct calls ignoring strand for call_class = "SBS" and "indel" with SBSindel_call_type = "mutation" so that each of these events is only counted once, while all other SBS and indel SBSindel_call_types and all MDB SBSindel_call_types will be counted once for each strand on which they occur. Also extract unique calls for SBSindel_call_type = "mutation". Also create a table formatted for vcf output. Used for burden, spectra, and vcf output.
finalCalls.bytype <- call_types_toanalyze %>%
	nest_join(
		finalCalls,
		by = join_by(call_class, call_type, SBSindel_call_type),
		name = "finalCalls"
	) %>%
	mutate(
		finalCalls = if_else(
			call_class %in% c("SBS","indel") & SBSindel_call_type == "mutation",
			finalCalls %>% map(
				function(x){
					x %>%
						distinct( #Not including call_class,call_type,SBSindel_call_type as these are identical for all rows within each finalCalls after nest_join
							run_id,zm,
							seqnames,start,end,ref_plus_strand,alt_plus_strand,reftnc_pyr,alttnc_pyr,reftnc_template_strand,alttnc_template_strand
						)
				}
			),
			finalCalls %>% map(
				function(x){
					x %>%
						distinct(
							run_id,zm,
							seqnames,strand,start,end,ref_plus_strand,alt_plus_strand,reftnc_pyr,alttnc_pyr,reftnc_template_strand,alttnc_template_strand
						)
				}
			)
		),
		
		finalCalls_unique = if_else(
			SBSindel_call_type == "mutation",
			finalCalls %>% map(
				function(x){
					x %>%
						distinct(
							seqnames,start,end,ref_plus_strand,alt_plus_strand,reftnc_pyr,alttnc_pyr,reftnc_template_strand,alttnc_template_strand
						)
				}
			),
			NA
		),
		
		finalCalls_for_vcf = finalCalls %>% map(
			function(x){
				x %>% 
					normalize_indels_for_vcf(
						BSgenome_name = yaml.config$BSgenome$BSgenome_name
					)
			}
		),
		
		finalCalls_unique_for_vcf = if_else(
			SBSindel_call_type == "mutation",
			finalCalls_unique %>%
				map(
					function(x){
						x %>% 
							normalize_indels_for_vcf(
								BSgenome_name = yaml.config$BSgenome$BSgenome_name
							)
					}
				),
			NA
		)
	)

#Calculate trinucleotide counts and fractions for call_class = "SBS" and "MDB". For all call types, calculate reftnc_pyr. For call_class "SBS" with SBSindel_call_type = "mutation", also calculate reftnc_pyr for unique calls, and for all other call types, calculate reftnc_template_strand.
finalCalls.reftnc_spectra <- finalCalls.bytype %>%
	mutate(
		
		finalCalls.reftnc_pyr = if_else(
			call_class %in% c("SBS","MDB"),
			finalCalls %>%
				map(
					function(x){
						x %>%
							count(reftnc_pyr, name = "count") %>%
							complete(reftnc_pyr, fill = list(count = 0)) %>%
							arrange(reftnc_pyr) %>%
							filter(!is.na(reftnc_pyr)) %>%
							mutate(fraction = count / sum(count))
					}
				),
			NA
		),
		
		finalCalls_unique.reftnc_pyr = if_else(
			call_class == "SBS" & SBSindel_call_type == "mutation",
			finalCalls_unique %>%
				map(
					function(x){
						if(is.null(x)){return(NA)}
						x %>%
							count(reftnc_pyr, name = "count") %>%
							complete(reftnc_pyr, fill = list(count = 0)) %>%
							arrange(reftnc_pyr) %>%
							filter(!is.na(reftnc_pyr)) %>%
							mutate(fraction = count / sum(count))
					}
				),
			NA
		),
		
		finalCalls.reftnc_template_strand = if_else(
			(call_class == "SBS" & SBSindel_call_type != "mutation") | call_class == "MDB",
			finalCalls %>%
				map(
					function(x){
						x %>%
							count(reftnc_template_strand, name = "count") %>%
							complete(reftnc_template_strand, fill = list(count = 0)) %>%
							arrange(reftnc_template_strand) %>%
							filter(!is.na(reftnc_template_strand)) %>%
							mutate(fraction = count / sum(count))
					}
				),
			NA
		)
	)

#Extract SBS spectra
finalCalls.reftnc_spectra <- finalCalls.reftnc_spectra %>%
	mutate(
		finalCalls.reftnc_pyr_spectrum = if_else(
			call_class == "SBS",
			finalCalls %>%
				map(
					function(x){
						x %>%
							mutate(
								channel = str_c(reftnc_pyr,">",alttnc_pyr) %>%
									factor(levels = sbs96_labels)
							) %>%
							count(channel, name = "count") %>%
							complete(channel, fill = list(count = 0)) %>%
							arrange(channel) %>%
							filter(!is.na(channel))
					}
				),
			NA
		),
		
		finalCalls_unique.reftnc_pyr_spectrum = if_else(
			call_class == "SBS" & SBSindel_call_type == "mutation",
			finalCalls_unique %>%
				map(
					function(x){
						if(is.null(x)){return(NA)}
						x %>%
							mutate(
								channel = str_c(reftnc_pyr,">",alttnc_pyr) %>%
									factor(levels = sbs96_labels)
							) %>%
							count(channel, name = "count") %>%
							complete(channel, fill = list(count = 0)) %>%
							arrange(channel) %>%
							filter(!is.na(channel))
					}
				),
			NA
		),
		
		finalCalls.reftnc_template_strand_spectrum = if_else(
			call_class == "SBS" & SBSindel_call_type != "mutation",
			finalCalls %>%
				map(
					function(x){
						x %>%
							mutate(
								channel = str_c(reftnc_template_strand,">",alttnc_template_strand) %>%
									factor(levels = sbs192_labels)
							) %>%
							count(channel, name = "count") %>%
							complete(channel, fill = list(count = 0)) %>%
							arrange(channel) %>%
							filter(!is.na(channel))
					}
				),
			NA
		)
	)

#Extract indel spectra (only for mutations)
BSgenome_for_indel.spectrum <- yaml.config$BSgenome$BSgenome_name %>%
	get %>%
	getSeq

finalCalls.reftnc_spectra <- finalCalls.reftnc_spectra %>%
	mutate(
		finalCalls.refindel_spectrum = pmap(
			list(call_class, SBSindel_call_type, finalCalls_for_vcf),
			function(x,y,z){
				if(x == "indel" & y == "mutation"){
					z %>%
						rename(
							CHROM = seqnames, POS = start,
							REF = ref_plus_strand, ALT = alt_plus_strand
						) %>%
						select(CHROM, POS, REF, ALT) %>%
						as.data.frame %>%
						indel.spectrum(BSgenome_for_indel.spectrum)
				}else{
					NULL
				}
			}
		),
		
		finalCalls_unique.refindel_spectrum = pmap(
			list(call_class, SBSindel_call_type, finalCalls_unique_for_vcf),
			function(x,y,z){
				if(x == "indel" & y == "mutation"){
					z %>%
						rename(
							CHROM = seqnames, POS = start,
							REF = ref_plus_strand, ALT = alt_plus_strand
						) %>%
						select(CHROM, POS, REF, ALT) %>%
						as.data.frame %>%
						indel.spectrum(BSgenome_for_indel.spectrum)
				}else{
					NULL
				}
			}
		)
	)

rm(BSgenome_for_indel.spectrum)
invisible(gc())

#Create sigfit-format spectra tables
finalCalls.reftnc_spectra <- finalCalls.reftnc_spectra %>%
	mutate(
		rowname_col = str_c(
			yaml.config$analysis_id, sample_id_toanalyze,
			chromgroup_toanalyze, filtergroup_toanalyze,
			call_class, call_type, SBSindel_call_type,
			sep="."
		),
		
		across(
			c(finalCalls.reftnc_pyr_spectrum, finalCalls_unique.reftnc_pyr_spectrum),
			function(x){
				map2(
					x, rowname_col,
					function(y, rn){
						if(is.null(y)){return(NULL)}
						y %>%
							pivot_wider(names_from = channel, values_from = count) %>%
							mutate(rowname = rn) %>%
							column_to_rownames("rowname") %>%
							as.data.frame
					}
				)
			},
			.names = "{.col}.sigfit"
		),
		
		across(
			c(finalCalls.reftnc_template_strand_spectrum),
			function(x){
				map2(
					x, rowname_col,
					function(y, rn){
						if(is.null(y)){return(NULL)}
						y %>%
							mutate(
								channel = if_else(
									str_sub(channel,2,2) %in% c("C","T"),
									str_c("T:",channel),
									str_c(
										"U:",
										str_sub(channel,1,3) %>% DNAStringSet %>% reverseComplement %>% as.character,
										">",
										str_sub(channel,5,7) %>% DNAStringSet %>% reverseComplement %>% as.character
									)
								) %>%
									factor(levels = sbs192_labels.sigfit)
							) %>%
							pivot_wider(names_from = channel, values_from = count) %>%
							mutate(rowname = rn) %>%
							column_to_rownames("rowname") %>%
							as.data.frame
					}
				)
			},
			.names = "{.col}.sigfit"
		),
		
		across(
			c(finalCalls.refindel_spectrum, finalCalls_unique.refindel_spectrum),
			function(x){
				map2(
					x, rowname_col,
					function(y, rn){
						if(is.null(y)){return(NULL)}
						y %>%
							indelspectrum.to.sigfit %>%
							mutate(rowname = rn) %>%
							column_to_rownames("rowname") %>%
							as.data.frame
					}
				)
			}
		)
	) %>%
	select(-rowname_col)

#Remove unnecessary columns
finalCalls.reftnc_spectra <- finalCalls.reftnc_spectra %>%
	select(-c(finalCalls, finalCalls_unique, finalCalls_for_vcf, finalCalls_unique_for_vcf))

cat("DONE\n")

######################
### Calculate call burdens
######################
cat("## Calculating call burdens...")
	
#Calculate number of interrogated base pairs (SBSindel_call_type = mutation) or bases (mismatch-ss, mismatch-ds, mismatch-os, match), for each call_type to analyze. Note, for mismatch-ds, each mismatch-ds is counted as 2 mismatches.
finalCalls.burdens <- call_types_toanalyze %>%
	select(-starts_with("MDB")) %>%
	left_join(
		bam.gr.filtertrack.bytype %>%
			mutate(
				interrogated_bases_or_bp = if_else(
					SBSindel_call_type == "mutation",
					bam.gr.filtertrack.coverage %>%
						map_dbl(
							function(x){
								x %>% mcols %>% with(duplex_coverage) %>% sum
							}
						),
					bam.gr.filtertrack.coverage %>%
						map_dbl(
							function(x){
								2 * (x %>% mcols %>% with(duplex_coverage) %>% sum)
							}
						)
				)
			) %>%
			select(names(call_types_toanalyze),interrogated_bases_or_bp),
		by = names(call_types_toanalyze)
	)

#Calculate number of calls for each call_type to analyze.
finalCalls.burdens <- finalCalls.burdens  %>%
	left_join(
		finalCalls.bytype,
		by = join_by(call_type,call_class,analyzein_chromgroups,SBSindel_call_type,filtergroup)
	) %>%
	mutate(
		num_calls = finalCalls %>%
			map_dbl(nrow),
		num_calls_unique = finalCalls_unique %>%
			map_dbl(
				function(x){
					if(!is.null(x)){x %>% nrow} else{NA_real_}
				}
			)
	)

#Calculate Poisson 95% confidence intervals for number of calls and number of unique calls
finalCalls.burdens <- finalCalls.burdens %>%
	mutate(
		ci = num_calls %>% map( function(x){poisson.test(x)$conf.int} ),
		ci_unique = map2(
			SBSindel_call_type,
			num_calls_unique,
			function(x,y){
				if(x == "mutation"){
					poisson.test(y)$conf.int
				}else{
					c(NA_real_,NA_real_)
				}
			}
		),
		
		num_calls_lci = map_dbl(ci,1),
		num_calls_uci = map_dbl(ci,2),
		
		num_calls_unique_lci = map_dbl(ci_unique,1),
		num_calls_unique_uci = map_dbl(ci_unique,2)
	) %>%
	select(-ci, -ci_unique)

#Calculate burdens for all calls and unique calls, with Poisson 95% confidence intervals
finalCalls.burdens <- finalCalls.burdens %>%
	mutate(
		burden_calls = num_calls / interrogated_bases_or_bp,
		burden_calls_lci = num_calls_lci / interrogated_bases_or_bp,
		burden_calls_uci = num_calls_uci / interrogated_bases_or_bp,
		
		burden_calls_unique = num_calls_unique / interrogated_bases_or_bp,
		burden_calls_unique_lci = num_calls_unique_lci / interrogated_bases_or_bp,
		burden_calls_unique_uci = num_calls_unique_uci / interrogated_bases_or_bp
	)

##TODO: Calculate burdens corrected for trinucleotide background

cat("DONE\n")

######################
### Calculate SBS and indel sensitivity
######################
#Create sensitivity tibble, set default sensitivity to 1 and source to 'default' for all call_types, except for call_class == 'MDB' whose source is set to 'yaml.config' and sensitivity is set per the yaml.config.
sensitivity <- call_types_toanalyze %>%
	mutate(
		sensitivity = if("MDB_sensitivity" %in% names(.)){
			if_else(call_class == "MDB", MDB_sensitivity, 1, ptype = numeric())
		}else{
			1
		},
		sensitivity_source = if_else(call_class == "MDB", "yaml.config", "default")
	) %>%
	select(-starts_with("MDB"))

#Create set of high-quality germline variants for sensitivity analysis if use_chromgroup is defined and is the current chromgroup
if(!is.null(sensitivity_parameters$use_chromgroup) & sensitivity_parameters$use_chromgroup == chromgroup_toanalyze){
	
	cat("## Calculating SBS and indel sensitivity...")
	
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
		vcf_file = sensitivity_parameters$gnomad_sensitivity_vcf,
		genome_fasta = yaml.config$genome_fasta,
		BSgenome_name = yaml.config$BSgenome$BSgenome_name,
		bcftools_bin = yaml.config$bcftools_bin
	) %>%
		as_tibble
	
	#Load germline VCF variants and filter to retain high-confidence variants
	high_confidence_germline_vcf_variants <- qs_read(str_c(yaml.config$cache_dir,"/",individual_id,".germline_vcf_variants.qs2")) %>%
		as_tibble %>%
		
		#Separate genotypes of each allele
		separate_wider_delim(GT, regex("[^[:digit:].]"), names=c("GT1","GT2"), cols_remove = FALSE) %>%
		
		#Group by germline_vcf_file so that filters are calculated separately for each germline VCF.
		group_by(germline_vcf_file) %>% 
		
		#Filter per sensitivity_threshold settings
		filter(
			(Depth >= quantile(Depth, sensitivity_parameters$SBS_min_Depth_quantile) & call_class == "SBS") |
				(Depth >= quantile(Depth, sensitivity_parameters$indel_min_Depth_quantile) & call_class == "indel"),
			(VAF >= sensitivity_parameters$SBS_min_VAF & call_class == "SBS") |
				(VAF >= sensitivity_parameters$indel_min_VAF & call_class == "indel"),
			(VAF <= sensitivity_parameters$SBS_max_VAF & call_class == "SBS") |
				(VAF <= sensitivity_parameters$indel_max_VAF & call_class == "indel"),
			(GQ >= quantile(GQ, sensitivity_parameters$SBS_min_GQ_quantile) & call_class == "SBS") |
				(GQ >= quantile(GQ, sensitivity_parameters$indel_min_GQ_quantile) & call_class == "indel"),
			(QUAL >= quantile(QUAL, sensitivity_parameters$SBS_min_QUAL_quantile) & call_class == "SBS") |
				(QUAL >= quantile(QUAL, sensitivity_parameters$indel_min_QUAL_quantile) & call_class == "indel"),
			(sensitivity_parameters$genotype == "heterozygous" & ((GT1 == "1" & GT2 != "1") | (GT1 != "1" & GT2 == "1"))) |
				(sensitivity_parameters$genotype == "homozygous" & (GT1 == "1" & GT2 == "1"))
		) %>%
		ungroup %>%
		
		#Keep variants detected in all germline vcfs
		count(seqnames, start, end, ref_plus_strand, alt_plus_strand, call_class, call_type, SBSindel_call_type) %>%
		filter(n == !!num_germline_vcf_files) %>%
		select(-n) %>%
		
		#Keep variants in chromgroup, and exclude variants in sex chromosomes and mitochondrial chromosome
		filter(
			seqnames %in% chroms_toanalyze,
			! seqnames %in% sex_chromosomes,
			! seqnames %in% mitochondrial_chromosome
		) %>%
		
		#Keep variants in gnomad_sensitivity_vcf
		semi_join(
			gnomad_sensitivity_vcf,
			by = names(.)
		)
	
	rm(num_germline_vcf_files, gnomad_sensitivity_vcf)
	invisible(gc())
	
	#Annotate for each variant in high_confidence_germline_vcf_variants how many times it was detected in germlineVariantCalls, counting 1 for each zm in which it was detected. 
	high_confidence_germline_vcf_variants <- high_confidence_germline_vcf_variants %>%
		left_join(
			germlineVariantCalls %>%
				distinct(
					across(all_of(c("run_id","zm",names(high_confidence_germline_vcf_variants))))
				) %>%
				count(
					across(all_of(names(high_confidence_germline_vcf_variants))),
					name = "num_zm_detected"
				),
			by = names(high_confidence_germline_vcf_variants)
		) %>%
		mutate(num_zm_detected = num_zm_detected %>% replace_na(0))
	
	#Annotate duplex coverage of high_confidence_germline_vcf_variants for each non-germline filter tracker. For insertions and deletions, annotate the minimum coverage at the left and right bases immediately flanking the insertion site/deleted bases.
	#For high_confidence_germline_vcf_variants, change insertion and deletion ranges to span immediately flanking bases
	high_confidence_germline_vcf_variants <- high_confidence_germline_vcf_variants %>%
		mutate(
			start = if_else(call_class == "indel", start - 1, start),
			end = if_else(call_class == "indel", end + 1, end)
		)
	
	#Get coverage for start and end positions so that coverage is calculated for those bases specifically, and then set coverage to the minimum between these.
	bam.gr.filtertrack.except_germline_filters.bytype <- bam.gr.filtertrack.except_germline_filters.bytype %>%
		nest_join(high_confidence_germline_vcf_variants, by = "call_type") %>%
		mutate(
			high_confidence_germline_vcf_variants = map2(
				high_confidence_germline_vcf_variants,
				bam.gr.filtertrack.coverage,
				function(x,y){
					
					gr <- x %>%
						makeGRangesFromDataFrame(
							seqinfo = yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
						)
					
					x %>%
						mutate(
							duplex_coverage = pmin(
								gr %>%
									resize(width = 1, fix = "start") %>%
									gr_1bp_cov(y),
								gr %>%
									resize(width = 1, fix = "end") %>%
									gr_1bp_cov(y)
							)
						)
				}
			)
		)
	
	#Join high_confidence_germline_vcf_variants.annotated to sensitivity tibble 
	sensitivity <- sensitivity %>%
		left_join(
			bam.gr.filtertrack.except_germline_filters.bytype %>%
				select(-bam.gr.filtertrack.coverage),
			by = join_by(call_type,call_class,analyzein_chromgroups,SBSindel_call_type,filtergroup)
		)
	
	rm(high_confidence_germline_vcf_variants, bam.gr.filtertrack.except_germline_filters.bytype)
	invisible(gc())
	
	#Sum number of high confidence germline VCF variant detections and coverage
	sensitivity <- sensitivity %>%
		mutate(
			high_confidence_germline_vcf_variants_sum_zm_detected = high_confidence_germline_vcf_variants %>%
				map_int( function(x){x$num_zm_detected %>% sum} ),
			high_confidence_germline_vcf_variants_sum_duplex_coverage = high_confidence_germline_vcf_variants %>%
				map_int( function(x){x$duplex_coverage %>% sum} )
		)
	
	#Calculate sensitivity.
	#Calculate sensitivity only if number of high-confidence germline variant detections is above the configured minimum required, otherwise keep as default of 1. 
	#If sensitivity_parameters$genotype = 'homozygous' (i.e. used only homozygous variants for sensitivity calculation): calculate sensitivity as sum(num_zm_detected) / sum(duplex_coverage)
	#If sensitivity_parameters$genotype = 'heterozygous', calculate sensitivity as 2 * sum(num_zm_detected) / sum(duplex_coverage).
	#Cap sensitivity to a max of 1, since it is possible to exceed 1 due to the above 'heterozygous' correction and edge cases. 
	sensitivity <- sensitivity %>% 
		mutate(
			
			#Check based on thresholds if sensitivity should be calculated
			calculate_sensitivity = (call_class == "SBS" &
															 	high_confidence_germline_vcf_variants_sum_zm_detected >= sensitivity_parameters$SBS_min_variant_detections &
															 	high_confidence_germline_vcf_variants_sum_duplex_coverage > 0) |
				
				(calculate_indel_sensitivity = call_class == "indel" &
				 	high_confidence_germline_vcf_variants_sum_zm_detected >= sensitivity_parameters$indel_min_variant_detections &
				 	high_confidence_germline_vcf_variants_sum_duplex_coverage > 0),
			
			sensitivity = case_when(
				calculate_sensitivity & sensitivity_parameters$genotype == "homozygous" ~
					high_confidence_germline_vcf_variants_sum_zm_detected / high_confidence_germline_vcf_variants_sum_duplex_coverage,
				
				calculate_sensitivity & sensitivity_parameters$genotype == "heterozygous" ~
					2 * (high_confidence_germline_vcf_variants_sum_zm_detected / high_confidence_germline_vcf_variants_sum_duplex_coverage),
				
				.default = sensitivity
			),
			
			sensitivity = pmin(1, sensitivity),
			
			sensitivity_source = if_else(calculate_sensitivity, "calculated", sensitivity_source)
			
		) %>%
		select(-calculate_sensitivity)
	
	#Change sensitivity to sqrt(sensitivity) for SBSindel_call_type = 'mismatch-ss' and 'mismatch-ds'. Do not do this transformation for SBSindel_call_type = 'mutation' since this involves detection in both strands. Although 'mismatch-ds' is detected in both strands, these are counted for burdens as two events for simplicity in output of mismatch-ds calls. SBSindel_call_type = 'mismatch-os' and 'match' are also not changed since these are for type MDB that is set in the yaml.config.
	sensitivity <- sensitivity %>%
		mutate(
			sensitivity = if_else(
				SBSindel_call_type %in% c("mismatch-ss","mismatch-ds"),
				sqrt(sensitivity),
				sensitivity
			)
		)
	
	cat("DONE\n")
	
}else if(!is.null(sensitivity_parameters$use_chromgroup) & sensitivity_parameters$use_chromgroup != chromgroup_toanalyze){{
	cat("Skipping sensitivity calculation since use_chromgroup not in currently analyzed chromgroup.\n")
	sensitivity <- NULL #Remove sensitivity tibble so it is not used.
}else if(is.null(sensitivity_parameters$use_chromgroup)){{
	cat("All sensitivities set to 1 since use_chromgroup not defined.\n")
}

######################
### Calculate estimated mutation error rate
######################
cat("## Calculating estimated mutation error rate...")

estimatedMutationErrorRate

#For ssDNA mismatches only, for each channel and total

- review Nanoseq code for calculating error rate again and update HiDEF-seq accordingly. See https://github.com/cancerit/NanoSeq/issues/92 and https://github.com/cancerit/NanoSeq/blob/4136d3ca943b96b2cf9013ac14183f098b8234be/R/nanoseq_results_plotter.R#L730. Specifically, do I need to divide the interrogated ssDNA bases by 2 in the denominator, since the two false positive calls must happen in opposite strands. But think about it carefully. Not sure about this.

cat("DONE\n")

######################
### Save results
######################
cat("## Saving results...")
qs_save(
	list(
		yaml.config = yaml.config,
		run_metadata = run_metadata,
		individual_id = individual_id,
		sample_id = sample_id_toanalyze,
		chromgroup = chromgroup_toanalyze,
		filtergroup = filtergroup_toanalyze,
		call_types = call_types_toanalyze
		molecule_stats.by_run_id = molecule_stats.by_run_id,
		molecule_stats.by_analysis_id = molecule_stats.by_analysis_id
		region_genome_filter_stats = region_genome_filter_stats,
		bam.gr.filtertrack.bytype.coverage_tnc = bam.gr.filtertrack.bytype,
		genome.reftnc = genome.reftnc,
		genome_chromgroup.reftnc = genome_chromgroup.reftnc,
		finalCalls = finalCalls,
		germlineVariantCalls = germlineVariantCalls,
		finalCalls.bytype = finalCalls.bytype,
		finalCalls.reftnc_spectra = finalCalls.reftnc_spectra,
		finalCalls.burdens = finalCalls.burdens,
		sensitivity = sensitivity,
		estimatedMutationErrorRate = estimatedMutationErrorRate
	),
	outputFile
)

cat("DONE\n")