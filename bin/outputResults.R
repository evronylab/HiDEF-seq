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
#Order of sbs trinucleotide context labels
sbs_labels <- c(
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
	"5:Del:M:1","5:Del:M:2","5:Del:M:3","5:Del:M:4","5:Del:M:5"
)

#All possible trinucleotides
trinucleotides_64 <- mkAllStrings(c("A","C","G","T"), 3)

#central pyrimidine trinucleotides
trinucleotides_32_pyr <- trinucleotides_64 %>% str_subset(".[CT].")

#List named with trinucleotides_32_pyr sequences, each containing the corresponding trinucleotides_64 sequences
trinucleotides_64_32_pyr_list <- split(
	trinucleotides_64,
	ifelse(
		substr(trinucleotides_64, 2, 2) %in% c("C","T"),
		trinucleotides_64,
		trinucleotides_64 %>% DNAStringSet %>% reverseComplement %>% as.character
	) %>%
		factor(levels = trinucleotides_32_pyr)
)

#Function to reduce 64 to 32 trinucleotide frequency with central pyrimidine. Input is a 2-column tibble.
trinucleotides_64to32 <- function(x, tri_column, count_column){
	x %>%
		mutate(
			!!sym(tri_column) := !!sym(tri_column) %>%
				factor(levels = trinucleotides_64) %>%
				fct_collapse(!!!trinucleotides_64_32_pyr_list) %>%
				factor(levels = trinucleotides_32_pyr)
		)%>%
		group_by(!!sym(tri_column)) %>%
		summarize(!!sym(count_column) := sum(!!sym(count_column)), .groups = "drop") %>%
		arrange(!!sym(tri_column))
}

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

#Function to write vcf from finalCalls/germlineCalls
write_vcf_from_calls <- function(calls, BSgenome_name, out_vcf){
	
	#Fixed fields
	CHROM <- calls$seqnames
	POS <- calls$start
	REF <- calls$ref_plus_strand %>% DNAStringSet
	ALT <- calls$alt_plus_strand %>% DNAStringSet
	QUAL <- rep(NA_real_, nrow(calls))
	FILTER <- rep("PASS", nrow(calls))
	
	# Row ranges; width from REF (>=1)
	calls.gr <- GRanges(
		seqnames = CHROM,
		ranges = IRanges(start = POS, width = REF %>% as.character %>% nchar),
		strand = "*",
		seqinfo = BSgenome_name %>% get %>% seqinfo
		)
	
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
		calls.gr %>% start
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
		call_types_toanalyze <- filterCallsFile %>% pluck("config","call_types")
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
			call_toanalyze == TRUE &
				if_all(contains("passfilter"), ~ .x == TRUE)
		)
	
	germlineVariantCalls[[i]] <- filterCallsFile %>%
		pluck("calls") %>%
		filter(
			call_class %in% c("SBS","indel") &
				SBSindel_call_type == "mutation" &
				if_all(contains("passfilter") & !all_of(germline_filters), ~ .x == TRUE) & #All non-germline filters == TRUE
				germline_vcf.passfilter == FALSE #Detected in at least one germline VCF
		)
	
	#Filtered read coverage of the genome, with and without germline_filters
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

#Combine molecule_stats to one tibble and sum stats. Also do a left_join to molecule_stats[[1]] without the value column to preserve the original row order of stats
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
cat("## Outputting configuration parameters...")

opt$config %>% file.copy(str_c(output_basename,".yaml.config.tsv"))

run_metadata %>% write_tsv(str_c(output_basename,".run_metadata.tsv"))

cat("DONE\n")

######################
### Output filtering statistics
######################
cat("## Outputting filtering statistics...:")

molecule_stats %>% write_tsv(str_c(output_basename,".molecule_stats.tsv"))

region_genome_filter_stats %>% write_tsv(str_c(output_basename,".region_genome_filter_stats.tsv"))

cat("DONE\n")

######################
### Output coverage and reference sequences of interrogated genome bases
######################
cat("## Outputting coverage and reference sequences of interrogated genome bases...")

#Convert duplex coverage RleList of final interrogated genome bases to GRanges for every base with a 'coverage' column, excluding zero coverage bases
bam.gr.filtertrack.bytype <- bam.gr.filtertrack.bytype %>%
	mutate(
		bam.gr.filtertrack.coverage = bam.gr.filtertrack.coverage %>%
			map(
				function(x){
					if(length(x) == 0){
						return(GRanges(
							coverage=integer(),
							seqinfo = yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
							))
						}
					
					x %>%
						imap(
							function(r,chr){
								pos <- which(r != 0L) 
								
								if(length(pos) == 0L){
									return(GRanges(
										coverage=integer(),
										seqinfo = yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
										))
									}
								
								rl <- runLength(r)
								rv <- runValue(r)
								nz <- rv != 0L
								
								GRanges(
									seqnames = chr,
									ranges = IRanges(start = pos, width = 1L),
									coverage = rep.int(rv[nz], rl[nz]),
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
	select(-coverage) %>%
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
	yaml.config$seqkit_bin,"faidx --quiet -l",tmpregions,genome_fasta,"|",
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

file.remove(tmpregions,tmpseqs)

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

#Output
for(i in seq_len(nrow(bam.gr.filtertrack.bytype))){
	bam.gr.filtertrack.bytype %>%
		pluck("bam.gr.filtertrack.coverage",i) %>%
		mutate(name = reftnc_pyr, score = coverage) %>% #rename to to allow output by 'export'
		export(
			con = str_c(
				output_basename,
				bam.gr.filtertrack.bytype$call_class[i],
				bam.gr.filtertrack.bytype$call_type[i],
				bam.gr.filtertrack.bytype$SBSindel_call_type[i],
				"bed",
				sep="."
			),
			format = "bed",
			index = TRUE
		)
}

cat("DONE\n")

######################
### Output trinucleotide background counts
######################
cat("## Outputting trinucleotide background counts...")

#Calculate trinucleotide counts for final interrogated genome bases
bam.gr.filtertrack.bytype <- bam.gr.filtertrack.bytype %>%
	mutate(
		bam.gr.filtertrack.reftnc_duplex_pyr = bam.gr.filtertrack.coverage %>%
			map(
				function(x){
					x %>%
						as_tibble %>%
						count(reftnc_pyr, wt = coverage, name = "count") %>%
						complete(reftnc_pyr, fill = list(count = 0)) %>%
						arrange(reftnc_pyr) %>%
						filter(!is.na(reftnc_pyr))
				}
			),
		
		bam.gr.filtertrack.reftnc_both_strands = bam.gr.filtertrack.coverage %>%
			map(
				function(x){
					bind_rows(
						x %>%
							as_tibble %>%
							count(reftnc_plus_strand, wt = coverage, name = "count") %>%
							rename(reftnc = reftnc_plus_strand),
						x %>%
							as_tibble %>%
							count(reftnc_minus_strand, wt = coverage, name = "count") %>%
							rename(reftnc = reftnc_minus_strand)
					) %>%
						group_by(reftnc) %>%
						summarize(count = sum(count), .groups = "drop") %>%
						complete(reftnc, fill = list(count = 0)) %>%
						arrange(reftnc) %>%
						filter(!is.na(reftnc))
				}
			),
		
		bam.gr.filtertrack.reftnc_both_strands_pyr = bam.gr.filtertrack.reftnc_both_strands %>%
			map(
				function(x){
					x %>%
						trinucleotides_64to32(
							tri_column = "reftnc",
							count_column = "count"
						) %>%
						rename(reftnc_pyr = reftnc)
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

#Calculate trinucleotide counts for the whole genome
genome.reftnc_plus_strand <- yaml.config$BSgenome$BSgenome_name %>%
	get %>%
	getSeq %>%
	trinucleotideFrequency(simplify.as = "collapsed") %>%
	enframe(name = "reftnc", value = "count") %>%
	mutate(reftnc = reftnc %>% factor(levels = trinucleotides_64))

genome.reftnc_minus_strand <- genome.reftnc_plus_strand %>%
	mutate(
		reftnc = reftnc %>%
			DNAStringSet %>%
			reverseComplement %>%
			as.character %>%
			factor(levels = trinucleotides_64)
	)

genome.reftnc_both_strands <- bind_rows(
		genome.reftnc_plus_strand,
		genome.reftnc_minus_strand
	) %>%
	group_by(reftnc) %>%
	summarize(count = sum(count), .groups = "drop") %>%
	arrange(reftnc) %>%
	filter(!is.na(reftnc))
	
genome.reftnc_duplex_pyr <- genome.reftnc_plus_strand %>%
	trinucleotides_64to32(tri_column = "reftnc", count_column = "count") %>%
	rename(reftnc_pyr = reftnc) %>%
	filter(!is.na(reftnc_pyr))

rm(genome.reftnc_plus_strand,genome.reftnc_minus_strand)

#Calculate trinucleotide counts for the analyzed chromgroup
genome_chromgroup.reftnc_plus_strand <- yaml.config$BSgenome$BSgenome_name %>%
	get %>%
	getSeq(chroms_toanalyze) %>%
	trinucleotideFrequency(simplify.as = "collapsed") %>%
	enframe(name = "reftnc", value = "count") %>%
	mutate(reftnc = reftnc %>% factor(levels = trinucleotides_64))

genome_chromgroup.reftnc_minus_strand <- genome_chromgroup.reftnc_plus_strand %>%
	mutate(
		reftnc = reftnc %>%
			DNAStringSet %>%
			reverseComplement %>%
			as.character %>%
			factor(levels = trinucleotides_64)
	)

genome_chromgroup.reftnc_both_strands <- bind_rows(
	genome_chromgroup.reftnc_plus_strand,
	genome_chromgroup.reftnc_minus_strand
) %>%
	group_by(reftnc) %>%
	summarize(count = sum(count), .groups = "drop") %>%
	arrange(reftnc) %>%
	filter(!is.na(reftnc))

genome_chromgroup.reftnc_duplex_pyr <- genome_chromgroup.reftnc_plus_strand %>%
	trinucleotides_64to32(tri_column = "reftnc", count_column = "count") %>%
	rename(reftnc_pyr = reftnc) %>%
	filter(!is.na(reftnc_pyr))

rm(genome_chromgroup.reftnc_plus_strand,genome_chromgroup.reftnc_minus_strand)

#Output trinucleotide counts
 #bam.gr.filtertracks
for(j in c("bam.gr.filtertrack.reftnc_duplex_pyr","bam.gr.filtertrack.reftnc_both_strands","bam.gr.filtertrack.reftnc_both_strands_pyr")){
	for(i in seq_len(nrow(bam.gr.filtertrack.bytype))){
		bam.gr.filtertrack.bytype %>%
			pluck(j,i) %>%
			write_tsv(
				file = str_c(
					output_basename,
					bam.gr.filtertrack.bytype$call_class[i],
					bam.gr.filtertrack.bytype$call_type[i],
					bam.gr.filtertrack.bytype$SBSindel_call_type[i],
					j,
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

cat("DONE\n")

######################
### Output final calls and germline variant calls
######################

cat("## Outputting final calls and detected germline variants...")

#Final calls
 #Extract final calls, split by type, and format list columns to comma-delimited
finalCalls.out <- call_types_toanalyze %>%
	distinct(call_type, call_class, SBSindel_call_type) %>%
	nest_join(
		finalCalls %>%
			mutate(
				across(
					where(is.list),
					function(x){map_chr(x, function(v){str_c(v, collapse = ",")})}
				)
			),
		by = join_by(call_type, call_class, SBSindel_call_type),
		name = "finalCalls"
	)

 #Output final calls to tsv, separately for each combination of call_class, call_type, SBSindel_call_type
finalCalls.out %>%
	pwalk(
		function(...){
			x <- list(...)
			
			x$finalCalls %>%
				write_tsv(
					str_c(output_basename,x$call_class,x$call_type,x$SBSindel_call_type,"finalCalls","tsv",sep=".")
					)
		}
	)

 #Output final calls to vcf, separately for each combination of call_class, call_type, SBSindel_call_type
finalCalls.out %>%
	pwalk(
		function(...){
			x <- list(...)
			
			x$finalCalls %>%
				normalize_indels_for_vcf(
					BSgenome_name = yaml.config$BSgenome$BSgenome_name
				) %>%
				write_vcf_from_calls(
					BSgenome_name = yaml.config$BSgenome$BSgenome_name,
					out_vcf = str_c(output_basename,x$call_class,x$call_type,x$SBSindel_call_type,"finalCalls","vcf",sep=".")
				)
		}
	)

#Germline variant calls
 #Extract germline variant calls and format list columns to comma-delimited
germlineVariantCalls.out <- germlineVariantCalls %>%
	mutate(
		across(
			where(is.list),
			function(x){map_chr(x, function(v){str_c(v, collapse = ",")})}
		)
	)

 #Output germline variant calls to tsv
germlineVariantCalls.out %>%
	write_tsv(
		str_c(output_basename,"germlineVariantCalls","tsv",sep=".")
	)

 #Output germline variant calls to vcf
germlineVariantCalls.out %>% 
	normalize_indels_for_vcf(
		BSgenome_name = yaml.config$BSgenome$BSgenome_name
	) %>%
	write_vcf_from_calls(
		BSgenome_name = yaml.config$BSgenome$BSgenome_name,
		out_vcf = str_c(output_basename,"germlineVariantCalls","vcf",sep=".")
	)

cat("DONE\n")

######################
### Calculate sensitivity
######################
cat("## Calculating sensitivity...")

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
	
	if sensitivity_threshold$SBS_min_variants or set sensitivity to 1
	if sensitivity_threshold$indel_min_variants or set sensitivity to 1
	
	#Load germline calls that pass all filters and intersect with sensitivity variant set to calculate sensitivity
	
	#Count how many opportunities there was to detect each sensitivity variant set in this chunk's non-germline filter tracker
	
	Take into account sex_chromosomes
	
	Correct final sensitivity if used het sensitivity_threshold$sensitivity_genotype, up to max of 1.0
	
	cat("DONE\n")
	
}else{
	cat("Setting sensitivity to 1 since gnomad_sensitivity_vcf is NULL.\n")
	
	sensitivity_sbs <- 1.0
	sensitivity_indel <- 1.0
}

######################
### Output call frequencies
######################
cat("## Outputting call frequencies...")
#Number of interrogated bases and base pairs -> calculate bases interrogated  = 2 x base pairs interrogated!
#all and unique, observed vs genome corrected vs sensitivity corrected vs both corrected
#upper and lower poisson conf int
#Interrogated bases, Interrogated base pairs

cat("DONE\n")

######################
### Output call spectra
######################
cat("## Outputting call spectra...")
#Plots and tables of observed and corrected, only for unique counts

cat("DONE\n")

######################
### Output estimated mutation error rate
######################
cat("## Outputting estimated mutation error rate...")
#For ssDNA mismatches only, for each channel and total

- review Nanoseq code for calculating error rate again and update HiDEF-seq accordingly. See https://github.com/cancerit/NanoSeq/issues/92 and https://github.com/cancerit/NanoSeq/blob/4136d3ca943b96b2cf9013ac14183f098b8234be/R/nanoseq_results_plotter.R#L730. Specifically, do I need to divide the interrogated ssDNA bases by 2 in the denominator, since the two false positive calls must happen in opposite strands. But think about it carefully. Not sure about this.

cat("DONE\n")

######################
### Output results in RDS file
######################
cat("## Outputting data in RDS file...")
qs_save(
	list(
		yaml.config = yaml.config,
		run_metadata = run_metadata,
		individual_id = individual_id,
		sample_id = sample_id_toanalyze,
		chromgroup = chromgroup_toanalyze,
		filtergroup = filtergroup_toanalyze,
		call_types = call_types_toanalyze,
		molecule_stats = molecule_stats,
		region_genome_filter_stats = region_genome_filter_stats,
		sensitivity_sbs = sensitivity_sbs,
		sensitivity_indel = sensitivity_indel,
		bam.gr.filtertrack.bytype = bam.gr.filtertrack.bytype,
		genome.reftnc_both_strands = genome.reftnc_both_strands,
		genome.reftnc_duplex_pyr = genome.reftnc_duplex_pyr,
		genome_chromgroup.reftnc_both_strands = genome_chromgroup.reftnc_both_strands,
		genome_chromgroup.reftnc_duplex_pyr = genome_chromgroup.reftnc_duplex_pyr,
		finalCalls = finalCalls,
		germlineVariantCalls = germlineVariantCalls,
		
	)
	str_c(output_basename,".output.RDS")
)

cat("DONE\n")