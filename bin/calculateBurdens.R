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

if(is.null(opt$config) | is.null(opt$sample_id_toanalyze) | is.null(opt$chromgroup_toanalyze) | is.null(opt$filtergroup_toanalyze) | is.null(opt$files) | is.null(opt$output) ){
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
	) %>%
	mutate(filtergroup = filtergroup %>% factor)

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
sbs96_labels.sigfit <- c(
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
	sbs96_labels.sigfit,
	str_c(
		str_sub(sbs96_labels.sigfit,1,3) %>% DNAStringSet %>% reverseComplement %>% as.character,
		">",
		str_sub(sbs96_labels.sigfit,5,7) %>% DNAStringSet %>% reverseComplement %>% as.character
	)
)

sbs192_labels.sigfit <- c(
	str_c("T:",sbs96_labels.sigfit),
	str_c("U:",sbs96_labels.sigfit)
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

######################
### Load data from filterCalls files
######################
cat("## Loading data from filterCalls files...\n > analysis chunk:")

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
			filterCallsFile %>% pluck("calls") %>% colnames %>% str_subset("germline_vcf_indel_region_filter"),
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

	#Load final calls that pass all filters, and germline variant calls that pass all filters ignoring the germline filters, since we will later also need the calls filtered by germline filters to calculate sensitivity
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
	
	#Remove temporary objects
	rm(filterCallsFile)
	invisible(gc())
}

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
### Calculate trinucleotide distributions of interrogated bases and the genome
######################
cat("## Calculating trinucleotide distributions of interrogated bases and the genome...")
#Convert coverage Rle to coverage runs and write to disk
 
chunk_runs <- 1e7 #Size of blocks to write to disk. Lower number takes longer, but decreases peak RAM.

bam.gr.filtertrack.bytype %>%
	mutate(row_id = row_number()) %>%
	pwalk(
		function(...){
			x <- list(...)
			
			output_file <- str_c(x$row_id,".bed")
			if(file.exists(output_file)){file.remove(output_file)}
			file.create(output_file)
			
			x$bam.gr.filtertrack.coverage %>%
				iwalk(
					function(rle, seqname){
						lens <- rle %>% runLength
						vals <- rle %>% runValue
						n <- lens %>% length
						
						offset <- 0L
						
						for(lo in seq.int(1L, n, by = chunk_runs)){
							hi <- min(lo + chunk_runs - 1L, n)
							len_slice <- lens[lo:hi]
							val_slice <- vals[lo:hi]
							
							nz <- which(val_slice != 0L)
							if(length(nz) == 0){
								offset <- offset + sum(len_slice)
								next
							}
							
							ends_chunk <- offset + cumsum(len_slice)
							
							starts_chunk_sub <- ends_chunk[nz]
							gt1 <- nz > 1L
							if(any(gt1)){starts_chunk_sub[gt1] <- ends_chunk[nz[gt1] - 1L]}
							if(nz[1L] == 1L){starts_chunk_sub[1L] <- offset}

							data.frame(
								seqname = seqname,
								start = starts_chunk_sub,
								end = ends_chunk[nz],
								value = val_slice[nz]
							) %>%
								data.table::fwrite(output_file, sep="\t", col.names = FALSE, append = TRUE, showProgress = FALSE)
							
							offset <- ends_chunk[length(ends_chunk)]
							
							invisible(gc())
						}
					}
				)
		}
	)

#Extract unique merged coverage runs
paste(
	"sort -m -k1,1n -k3,3n -k4,4n",
	paste(
		"<(awk -v OFS='\t' 'FNR==NR{ord[$1]=NR-1;next}{print ord[$1],$1,$2,$3}'",
		yaml.config$genome_fai,str_c(bam.gr.filtertrack.bytype %>% nrow %>% seq_len,".bed"),
		")",
		collapse=" "
	),
	"| cut -f2- |",
	yaml.config$bedtools_bin,"merge -i stdin > all.bed"
) %>%
	system2(command="/bin/bash", args="-s", input=.) %>% 
	invisible

#Expand to per-base coverage runs and annotate with genome trinucleotide sequences
genome_trinuc_file <- str_c(yaml.config$cache_dir,"/",yaml.config$BSgenome$BSgenome_name, ".bed.gz")

paste(
	yaml.config$bedtools_bin, "makewindows -w 1 -b all.bed |",
	yaml.config$bedtools_bin, "intersect -sorted -loj -wa -wb -a stdin -b",genome_trinuc_file, "-g", yaml.config$genome_fai,
	"| cut -f 1-3,7 > all.trinuc.bed"
) %>%
	system2(command="/bin/bash", args="-s", input=.) %>% 
	invisible

file.remove("all.bed") %>% invisible

#Intersect individual coverage run BED files with per-base trinucleotide-annotated genome, and output: 1) final per-base BED output (bgzipped, tabix indexed): seqnames, start, end, duplex_coverage, reftnc_plus_strand; 2) counts for every reftnc_plus_strand trinucleotide context
bam.gr.filtertrack.bytype %>%
	mutate(row_id = row_number()) %>%
	pwalk(
		function(...){
			x <- list(...)
			
			cov_output_file <- str_c(
				yaml.config$analysis_id,
				individual_id,
				sample_id_toanalyze,
				chromgroup_toanalyze,
				filtergroup_toanalyze,
				x$call_class,
				x$call_type,
				x$SBSindel_call_type,
				"coverage_reftnc",
				"bed.gz",
				sep="."
			)
			
			paste(
				yaml.config$bedtools_bin, "intersect -sorted -wa -wb",
				"-a", str_c(x$row_id,".bed"), "-b all.trinuc.bed",
				"-g", yaml.config$genome_fai, "|",
				"awk -v OFS='\t'", str_c("-v row_id=",x$row_id),
				"'{print $5,$6,$7,$4,$8; sum[$8]+=$4}
					END{
						if(length(sum)==0){
							print row_id, \"NA\", 0 > (row_id \".reftnc_plus_strand.tsv\")
						}else{
							for(k in sum){print row_id, k, sum[k] > (row_id \".reftnc_plus_strand.tsv\")}
						}
					}' |",
				yaml.config$bgzip_bin, "-c >",
				cov_output_file,
				"&&",
				yaml.config$tabix_bin, "-@2 -s1 -b2 -e3", cov_output_file
			) %>%
				system2(command="/bin/bash", args="-s", input=.) %>% 
				invisible
			
			file.remove(str_c(x$row_id,".bed")) %>% invisible
		}
	)

file.remove("all.trinuc.bed") %>% invisible

#Remove coverage_reftnc.bed.gz and coverage_reftnc.bed.gz.tbi files for SBS/mismatch-ss if it is not in call_types_toanalyze, since the data is only being retained for calculation of SBS mutation error probability.
if(
	call_types_toanalyze %>%
		filter(call_class=="SBS", SBSindel_call_type=="mismatch-ss") %>%
		nrow == 0
	){
	cov_output_file <- str_c(
		yaml.config$analysis_id,
		individual_id,
		sample_id_toanalyze,
		chromgroup_toanalyze,
		filtergroup_toanalyze,
		"SBS",
		"SBS",
		"mismatch-ss",
		"coverage_reftnc",
		"bed.gz",
		sep="."
	)
	
	file.remove(cov_output_file, str_c(cov_output_file,".tbi")) %>% invisible
}

#Add reftnc_plus_strand results to bam.gr.filtertrack.bytype and annotate reftnc_minus_strand
bam.gr.filtertrack.bytype <- bam.gr.filtertrack.bytype %>%
	rowid_to_column("row_id") %>%
	nest_join(
		list.files(pattern="*.reftnc_plus_strand.tsv", full.names=TRUE) %>%
			read_tsv(col_names = c("row_id", "reftnc", "count"), col_types="ici") %>%
			mutate(reftnc = reftnc %>% factor(levels = trinucleotides_64)) %>%
			group_by(row_id) %>%
			complete(reftnc, fill = list(count = 0)) %>%
			arrange(reftnc) %>%
			filter(!is.na(reftnc)) %>% #This filters any reftnc's (such as '.' from contig edges) that are not in trinucleotides_64, as setting levels = trinucleotides_64 makes anything not in trinucleotides_64 become NA
			mutate(
				count = count %>% as.integer,
				fraction = count / sum(count)
			) %>%
			ungroup,
		by = "row_id",
		name = "bam.gr.filtertrack.reftnc_plus_strand"
	) %>%
	select(-row_id) %>%
	
	mutate(
		bam.gr.filtertrack.reftnc_minus_strand = bam.gr.filtertrack.reftnc_plus_strand %>%
			map(function(x){
				x %>%
					mutate(
						reftnc = reftnc %>%
							DNAStringSet %>%
							reverseComplement %>%
							as.character %>%
							factor(levels = trinucleotides_64)
					)
			})
	)

file.remove(
	str_c(bam.gr.filtertrack.bytype %>% nrow %>% seq_len, ".reftnc_plus_strand.tsv")
) %>%
	invisible

#Annotate reftnc_pyr and reftnc_both_strands
bam.gr.filtertrack.bytype <- bam.gr.filtertrack.bytype %>%
	mutate(
		bam.gr.filtertrack.reftnc_pyr = bam.gr.filtertrack.reftnc_plus_strand %>%
			map(function(x){
				x %>%
					trinucleotides_64to32(tri_column = "reftnc", count_column = "count") %>%
					rename(reftnc_pyr = reftnc) %>%
					mutate(fraction = count / sum(count))
			}),
		
		bam.gr.filtertrack.reftnc_both_strands = map2(
			bam.gr.filtertrack.reftnc_plus_strand,
			bam.gr.filtertrack.reftnc_minus_strand,
			function(x,y){
				bind_rows(x,y) %>%
					group_by(reftnc) %>%
					summarize(count = sum(count), .groups = "drop") %>%
					arrange(reftnc) %>%
					mutate(fraction = count / sum(count))
			}
		)
	) %>%
	select(-bam.gr.filtertrack.reftnc_plus_strand, -bam.gr.filtertrack.reftnc_minus_strand)

#Calculate trinucleotide distributions of the whole genome and of the genome in the analyzed chromgroup
 #Function to extract trinucleotide distribution for selected chromosomes
get_genome_reftnc <- function(BSgenome_name, chroms){
	
	reftnc_plus_strand <- BSgenome_name %>%
		get %>%
		getSeq(chroms) %>%
		DNAStringSet %>%
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
		mutate(fraction = count / sum(count))
	
	reftnc_both_strands <- bind_rows(
		reftnc_plus_strand,
		reftnc_minus_strand
	) %>%
		group_by(reftnc) %>%
		summarize(count = sum(count), .groups = "drop") %>%
		arrange(reftnc) %>%
		filter(!is.na(reftnc)) %>%
		mutate(fraction = count / sum(count))
	
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

#Define columns to discard from each finalCalls table after the nest join, depending on the call_class, as these are not relevant for the final tsv and vcf outputs
alltypes_cols_discard <- c("deletion.bothstrands.startendmatch")

sbs_indel_cols_discard <- c("MDB_score")

sbs_mdb_cols_discard <- c("indel_width")

indel_cols_discard <- c(
	"reftnc_plus_strand", "alttnc_plus_strand",
	"reftnc_pyr", "alttnc_pyr",
	"reftnc_synthesized_strand", "reftnc_template_strand",
	"alttnc_synthesized_strand", "alttnc_template_strand"
)

#Define columns that are identical between zm and strands to either '_keep' or '_discard' in the subsequent pivot_wider for SBSindel_call_type == "mutation" call tables. Not including call_class, call_type, SBSindel_call_type as these are removed from finalCalls after the below nest_join.
zm_identical_cols_keep <- c(
	"analysis_chunk", "run_id", "zm"
)

strand_identical_cols_keep <- c(
	"seqnames", "start", "end", "ref_plus_strand", "alt_plus_strand",
	"reftnc_plus_strand", "alttnc_plus_strand", "reftnc_pyr", "alttnc_pyr"
)

strand_identical_cols_discard <- c(
	"strand",
	"call_class.opposite_strand", "call_type.opposite_strand",
	"alt_plus_strand.opposite_strand"
)
	
#Nest_join finalCalls for each call_class x call_type x SBSindel_call_type combination, and add SBS/mismatch-ss if not in call_types_toanalyze for downstream SBS mutation error rate calculation, and collapse to distinct calls ignoring strand for call_class = "SBS" and "indel" with SBSindel_call_type = "mutation" so that each of these events is only counted once, while all other SBS and indel SBSindel_call_types and all MDB SBSindel_call_types will be counted once for each strand on which they occur. Also extract unique calls for SBSindel_call_type = "mutation". Then convert to a table formatted for tsv output. Before the nest join, remove "germline_vcf_types_detected", "germline_vcf_files_detected", and "passfilter" columns as these are empty, empty, and TRUE, respectively for all finalCalls. Used for burden, spectra, and vcf output.
finalCalls.bytype <- call_types_toanalyze %>%
	
	#remove columns that are not needed here
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
	
	#Nest join finalCalls
	nest_join(
		finalCalls %>%
			#Remove filter annotation columns since all finalCalls pass filters so this information is not needed for final tsv and vcf outputs
			select(-c("germline_vcf_types_detected", "germline_vcf_files_detected", contains("passfilter"))) %>%
			#Reformat list columns to be comma-delimited.
			mutate(
				across(
					where(is.list),
					function(x){x %>% map_chr(function(v){str_c(v, collapse = ",")})}
				)
			),
		by = join_by(call_class, call_type, SBSindel_call_type),
		name = "finalCalls"
	) %>%
	
	#Discard columns that are not needed for each call_class
	mutate(
		finalCalls = map2(
			call_class, finalCalls,
			function(x,y){
				if(x == "SBS"){
					y %>% select(-all_of(c(alltypes_cols_discard, sbs_indel_cols_discard, sbs_mdb_cols_discard)))
				}else if(x == "indel"){
					y %>% select(-all_of(c(alltypes_cols_discard, sbs_indel_cols_discard, indel_cols_discard)))
				}else if(x == "MDB"){
					y %>% select(-all_of(c(alltypes_cols_discard, sbs_mdb_cols_discard)))
				}
			}
		)
	) %>%
	
	#Make tables for tsv and vcf output
	mutate(
		
		#Calls for tsv output
		finalCalls_for_tsv = if_else(
			call_class %in% c("SBS","indel") & SBSindel_call_type == "mutation",
			finalCalls %>% map(
				function(x){
					x %>%
						
						#Change strand levels so that new column names after pivot_wider have suffixes ref_strand_(plus/minus)_read instead of +/-.
						mutate(
							strand = strand %>%
								fct_recode(
									refstrand_plus_read  = "+",
									refstrand_minus_read = "-"
								)
						) %>%
						
						#Collapse to one row per mutation
						pivot_wider(
							id_cols = all_of(c(zm_identical_cols_keep, strand_identical_cols_keep)),
							names_from = strand,
							values_from = -all_of(c(zm_identical_cols_keep, strand_identical_cols_keep, strand_identical_cols_discard)),
							names_glue = "{.value}_{strand}"
						)
				}
			),
			finalCalls
		),
		
		#Unique calls for tsv output
		finalCalls_unique_for_tsv = if_else(
			SBSindel_call_type == "mutation",
			finalCalls_for_tsv %>% map(
				function(x){
					x %>% distinct(across(all_of(strand_identical_cols_keep)))
				}
			),
			NA
		),
		
		#Calls for vcf output
		finalCalls_for_vcf = finalCalls_for_tsv %>% map(
			function(x){
				x %>% 
					normalize_indels_for_vcf(
						BSgenome_name = yaml.config$BSgenome$BSgenome_name
					)
			}
		),
		
		#Unique calls for vcf output
		finalCalls_unique_for_vcf = if_else(
			SBSindel_call_type == "mutation",
			finalCalls_unique_for_tsv %>%
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
	) %>%
	select(-finalCalls)

#Calculate trinucleotide counts and fractions for call_class = "SBS" and "MDB". For all call types (including SBS/mismatch-ss even if not in call_types_toanalyze, for downstream SBS mutation error calculation), calculate reftnc_pyr. For call_class "SBS" with SBSindel_call_type = "mutation", also calculate reftnc_pyr for unique calls, and for all other call types, calculate reftnc_template_strand.
finalCalls.reftnc_spectra <- finalCalls.bytype %>%
	mutate(
		
		finalCalls.reftnc_pyr = pmap(
			list(call_class, finalCalls_for_tsv),
			function(x,z){
				if(x %in% c("SBS","MDB")){
					z %>%
						count(reftnc_pyr, name = "count") %>%
						complete(reftnc_pyr, fill = list(count = 0)) %>%
						arrange(reftnc_pyr) %>%
						filter(!is.na(reftnc_pyr)) %>%
						mutate(fraction = count / sum(count))
				}else{
					NULL
				}
			}
		),
		
		finalCalls_unique.reftnc_pyr = pmap(
			list(call_class, SBSindel_call_type, finalCalls_unique_for_tsv),
			function(x,y,z){
				if(x == "SBS" & y == "mutation"){
					z %>%
						count(reftnc_pyr, name = "count") %>%
						complete(reftnc_pyr, fill = list(count = 0)) %>%
						arrange(reftnc_pyr) %>%
						filter(!is.na(reftnc_pyr)) %>%
						mutate(fraction = count / sum(count))
				}else{
					NULL
				}
			}
		),
		
		finalCalls.reftnc_template_strand = pmap(
			list(call_class, SBSindel_call_type, finalCalls_for_tsv),
			function(x,y,z){
				if((x == "SBS" & y != "mutation") | x == "MDB"){
					z %>%
						count(reftnc_template_strand, name = "count") %>%
						complete(reftnc_template_strand, fill = list(count = 0)) %>%
						arrange(reftnc_template_strand) %>%
						filter(!is.na(reftnc_template_strand)) %>%
						mutate(fraction = count / sum(count))
				}else{
					NULL
				}
			}
		)
	)

#Remove SBS/mismatch-ss from finalCalls.bytype if it is not part of call_types_toanalyze, since it was only retained until here in order to calculate SBS mutation error probability.
finalCalls.bytype <- finalCalls.bytype %>%
	semi_join(
		call_types_toanalyze,
		by = join_by(call_type, call_class, SBSindel_call_type, filtergroup)
	)

#Extract SBS spectra
finalCalls.reftnc_spectra <- finalCalls.reftnc_spectra %>%
	mutate(
		finalCalls.reftnc_pyr_spectrum = pmap(
			list(call_class, finalCalls_for_tsv),
			function(x,z){
				if(x == "SBS"){
					z %>%
						mutate(
							channel = str_c(reftnc_pyr,">",alttnc_pyr) %>%
								factor(levels = sbs96_labels.sigfit)
						) %>%
						count(channel, name = "count") %>%
						complete(channel, fill = list(count = 0)) %>%
						arrange(channel) %>%
						filter(!is.na(channel)) %>%
						mutate(fraction = count / sum(count))
				}else{
					NULL
				}
			}
		),
		
		finalCalls_unique.reftnc_pyr_spectrum = pmap(
			list(call_class, SBSindel_call_type, finalCalls_unique_for_tsv),
			function(x,y,z){
				if(x == "SBS" & y == "mutation"){
					z %>%
						mutate(
							channel = str_c(reftnc_pyr,">",alttnc_pyr) %>%
								factor(levels = sbs96_labels.sigfit)
						) %>%
						count(channel, name = "count") %>%
						complete(channel, fill = list(count = 0)) %>%
						arrange(channel) %>%
						filter(!is.na(channel)) %>%
						mutate(fraction = count / sum(count))
				}else{
					NULL
				}
			}
		),
		
		finalCalls.reftnc_template_strand_spectrum = pmap(
			list(call_class, SBSindel_call_type, finalCalls_for_tsv),
			function(x,y,z){
				if(x == "SBS" & y != "mutation"){
					z %>%
						mutate(
							channel = str_c(reftnc_template_strand,">",alttnc_template_strand) %>%
								factor(levels = sbs192_labels)
						) %>%
						count(channel, name = "count") %>%
						complete(channel, fill = list(count = 0)) %>%
						arrange(channel) %>%
						filter(!is.na(channel)) %>%
						mutate(fraction = count / sum(count))
				}else{
					NULL
				}
			}
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
							select(channel, count) %>%
							pivot_wider(names_from = channel, values_from = count) %>%
							mutate(rowname = rn) %>%
							column_to_rownames("rowname")
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
							select(channel, count) %>%
							pivot_wider(names_from = channel, values_from = count) %>%
							mutate(rowname = rn) %>%
							column_to_rownames("rowname")
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
							column_to_rownames("rowname")
					}
				)
			},
			.names = "{.col}.sigfit"
		)
	) %>%
	select(-rowname_col)

#Remove unnecessary columns
finalCalls.reftnc_spectra <- finalCalls.reftnc_spectra %>%
	select(-c(finalCalls_for_tsv, finalCalls_unique_for_tsv, finalCalls_for_vcf, finalCalls_unique_for_vcf))

cat("DONE\n")

######################
### Calculate call burdens
######################
cat("## Calculating call burdens...")
#Note: Sensitivity-corrected burdens are not calculated in this calculateBurdens script. They are calculated later in the outputResults script, since sensitivity values are calculated for one chromgroup and then shared across chromgroups, while this script only analyzes one chromgroup.

#A. Call burdens not corrected for interrogated vs genome trinucleotide distributions

#Calculate number of interrogated base pairs (SBSindel_call_type = mutation) or bases (mismatch-ss, mismatch-ds, mismatch-os, match), for each call_type to analyze. Note, for mismatch-ds, each mismatch-ds is counted as 2 mismatches so burdens use interrogated bases rather than base pairs.
finalCalls.burdens <- bam.gr.filtertrack.bytype %>%
	
	#Exclude SBS/mismatch-ss row if it was only added to bam.gr.filtertrack.bytype for downstream calculation of SBS mutation error probability. 
	semi_join(
		call_types_toanalyze,
		by = join_by(call_type, call_class, SBSindel_call_type, filtergroup)
	) %>%
	
	mutate(
		cov_sum = bam.gr.filtertrack.coverage %>%
			map_dbl(function(x){
					x %>% map_dbl(sum) %>% sum
			}),
		
		interrogated_bases_or_bp = if_else(
			SBSindel_call_type == "mutation",
			cov_sum,
			cov_sum * 2
		)
	) %>%
	select(-cov_sum, -starts_with("bam.gr.filtertrack"))

#Calculate number of calls for each call_type to analyze.
finalCalls.burdens <- finalCalls.burdens %>%
	left_join(
		bind_rows(
			#Not unique calls
			finalCalls.bytype %>%
				mutate(
					num_calls = finalCalls_for_tsv %>% map_dbl(nrow),
					num_calls_noncorrected = NA_real_,
					unique_calls = if_else(SBSindel_call_type == "mutation", FALSE, NA),
					reftnc_corrected = if_else(call_class %in% c("SBS","MDB"), FALSE, NA),
					reftnc_corrected_chromgroup = NA,
					sensitivity_corrected = FALSE
				),
		
			#Unique calls (only mutations)
			finalCalls.bytype %>%
				filter(!map_lgl(finalCalls_unique_for_tsv,is.null)) %>%
				mutate(
					num_calls = finalCalls_unique_for_tsv %>% map_dbl(nrow),
					num_calls_noncorrected = NA_real_,
					unique_calls = TRUE,
					reftnc_corrected = if_else(call_class =="SBS", FALSE, NA),
					reftnc_corrected_chromgroup = NA,
					sensitivity_corrected = FALSE
				)
		),
		by = join_by(call_type, call_class, SBSindel_call_type, filtergroup)
	) %>%
	select(-starts_with("finalCalls"))

##B. Call burdens corrected for interrogated vs genome trinucleotide distributions
#Performed for SBS and MDB.
	#Note: interrogated bases/base pairs are calculated here as the sum of counts of each trinucleotide in bam.gr.filtertrack.bytype rather than directly from the bam.gr.filtertrack.bytype coverage as before, since there may be rare edge cases of interrogated bases being at contig boundaries so they won't contribute to trinucleotide count, and likewise the analyzed calls also could have this issue. So it is more consistent for trinucleotide distribution-corrected analysis to calculate everything using trinucleotide counts.

 #Join trinucleotide distributions of calls and genome
finalCalls.reftnc_spectra.genome_correction <- finalCalls.reftnc_spectra %>%
	#Exclude SBS/mismatch-ss row if it was only added to finalCalls.reftnc_spectra for downstream calculation of SBS mutation error probability. 
	semi_join(
		call_types_toanalyze,
		by = join_by(call_type, call_class, SBSindel_call_type, filtergroup)
	) %>%
	
	left_join(
		bam.gr.filtertrack.bytype %>% select(-bam.gr.filtertrack.coverage),
		by = join_by(call_type, call_class, SBSindel_call_type, filtergroup)
	) %>%
	filter(call_class %in% c("SBS","MDB"))

 #Helper function to get the call counts corrected for the distribution of trinucleotides in interrogated bases/base pairs vs the genome.
#		x and y: calls reftnc and bam.gr.filtertrack reftnc, respectively
#		reftnc_cols: reftnc columns in x and y to join by (passed as expr(col1 == col2))
#		ratio_col: column by which to divide the counts column of x
#		unique_col_annotation: Label to give the unique_calls column in the resulting tibble
get_burden_data <- function(x, y, reftnc_cols, ratio_col, unique_col_annotation){
	df <- left_join(x, y, by = join_by(!!reftnc_cols)) %>%
		mutate(count_corrected = count.x / !!sym(ratio_col))
	
	return(
		tibble(
			interrogated_bases_or_bp = df %>% pull(count.y) %>% sum,
			num_calls = df %>% pull(count_corrected) %>% sum,
			num_calls_noncorrected = df %>% pull(count.x) %>% sum,
			unique_calls = unique_col_annotation,
			reftnc_corrected = TRUE,
			reftnc_corrected_chromgroup = ratio_col %>% str_remove("fraction_ratio_to_"),
			sensitivity_corrected = FALSE
		)
	)
}

 #SBS mutations: Not unique calls and unique calls, reftnc_pyr_corrected for whole genome and for chromgroup
finalCalls.reftnc_spectra.genome_correction.SBSmutations <- finalCalls.reftnc_spectra.genome_correction %>%
	filter(call_class == "SBS" & SBSindel_call_type == "mutation")

finalCalls.burdens <- finalCalls.burdens %>%
	bind_rows(
		bind_rows(
			
			#Not unique calls, whole-genome corrected
			finalCalls.reftnc_spectra.genome_correction.SBSmutations %>%
				mutate(
					burden_data = map2(
						finalCalls.reftnc_pyr,
						bam.gr.filtertrack.reftnc_pyr,
						function(x,y){get_burden_data(x,y,expr(reftnc_pyr==reftnc_pyr),"fraction_ratio_to_genome",FALSE)}
					)
				),
				
			#Unique calls, whole-genome corrected
			finalCalls.reftnc_spectra.genome_correction.SBSmutations %>%
				mutate(
					burden_data = map2(
						finalCalls_unique.reftnc_pyr,
						bam.gr.filtertrack.reftnc_pyr,
						function(x,y){get_burden_data(x,y,expr(reftnc_pyr==reftnc_pyr),"fraction_ratio_to_genome",TRUE)}
					)
				),
				
			#Not unique calls, genome chromgroup corrected
			finalCalls.reftnc_spectra.genome_correction.SBSmutations %>%
				mutate(
					burden_data = map2(
							finalCalls.reftnc_pyr,
							bam.gr.filtertrack.reftnc_pyr,
							function(x,y){get_burden_data(x,y,expr(reftnc_pyr==reftnc_pyr),"fraction_ratio_to_genome_chromgroup",FALSE)}
					)
				),
				
			#Unique calls, genome chromgroup corrected
			finalCalls.reftnc_spectra.genome_correction.SBSmutations %>%
				mutate(
					burden_data = map2(
						finalCalls_unique.reftnc_pyr,
						bam.gr.filtertrack.reftnc_pyr,
						function(x,y){get_burden_data(x,y,expr(reftnc_pyr==reftnc_pyr),"fraction_ratio_to_genome_chromgroup",TRUE)}
					)
				)
		) %>%
		select(call_type, call_class, SBSindel_call_type, filtergroup, burden_data) %>%
		unnest_wider(burden_data)
	)

 #SBS non-mutations and MDB: Not unique calls, reftnc for whole genome and for chromgroup
finalCalls.reftnc_spectra.genome_correction.SBSnonmutations <- finalCalls.reftnc_spectra.genome_correction %>%
	filter((call_class == "SBS" & SBSindel_call_type != "mutation") | call_class == "MDB")

finalCalls.burdens <- finalCalls.burdens %>%
	bind_rows(
		bind_rows(
			
			#Not unique calls, whole-genome corrected
			finalCalls.reftnc_spectra.genome_correction.SBSnonmutations %>%
				mutate(
					burden_data = map2(
						finalCalls.reftnc_template_strand,
						bam.gr.filtertrack.reftnc_both_strands,
						function(x,y){get_burden_data(x,y,expr(reftnc_template_strand==reftnc),"fraction_ratio_to_genome",FALSE)}
					)
				),
			
			#Not unique calls, genome chromgroup corrected
			finalCalls.reftnc_spectra.genome_correction.SBSnonmutations %>%
				mutate(
					burden_data = map2(
						finalCalls.reftnc_template_strand,
						bam.gr.filtertrack.reftnc_both_strands,
						function(x,y){get_burden_data(x,y,expr(reftnc_template_strand==reftnc),"fraction_ratio_to_genome_chromgroup",FALSE)}
					)
				)
		) %>%
		select(call_type, call_class, SBSindel_call_type, filtergroup, burden_data) %>%
		unnest_wider(burden_data)
	)

rm(finalCalls.reftnc_spectra.genome_correction, finalCalls.reftnc_spectra.genome_correction.SBSmutations, finalCalls.reftnc_spectra.genome_correction.SBSnonmutations)
invisible(gc())

#Calculate burdens and Poisson 95% confidence intervals for number of calls and burdens
finalCalls.burdens <- finalCalls.burdens %>%
	mutate(
		ci = num_calls %>% map( function(x){cipoisson(x)} ),
		
		num_calls_lci = map_dbl(ci,1),
		num_calls_uci = map_dbl(ci,2),
		
		burden_calls = num_calls / interrogated_bases_or_bp,
		burden_calls_lci = num_calls_lci / interrogated_bases_or_bp,
		burden_calls_uci = num_calls_uci / interrogated_bases_or_bp
	) %>%
	select(-ci)

cat("DONE\n")

######################
### Calculate SBS and indel sensitivity
######################
#Create sensitivity tibble, set default sensitivity to 1 and source to 'default' for all call_types, except for call_class == 'MDB' whose source is set to 'yaml.config' and sensitivity is set per the yaml.config.
sensitivity <- call_types_toanalyze %>%
	select(-analyzein_chromgroups, -starts_with("MDB")) %>% #remove columns that are not needed here
	mutate(
		sensitivity = if("MDB_sensitivity" %in% names(.)){
			if_else(call_class == "MDB", MDB_sensitivity, 1, ptype = numeric())
		}else{
			1
		},
		sensitivity_source = if_else(call_class == "MDB", "yaml.config", "default")
	)

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
	
	#Load sensitivity_vcf
	sensitivity_vcf <- load_vcf(
		vcf_file = sensitivity_parameters$sensitivity_vcf,
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
		
		#Keep variants in sensitivity_vcf
		semi_join(
			sensitivity_vcf,
			by = names(.)
		)
	
	rm(num_germline_vcf_files, sensitivity_vcf)
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
			by = join_by(call_type, call_class, SBSindel_call_type, filtergroup)
		)
	
	rm(high_confidence_germline_vcf_variants, bam.gr.filtertrack.except_germline_filters.bytype)
	invisible(gc())
	
	#Sum number of high confidence germline VCF variant detections and coverage, and remove high_confidence_germline_vcf_variants that is no longer needed
	sensitivity <- sensitivity %>%
		mutate(
			high_confidence_germline_vcf_variants_sum_zm_detected = high_confidence_germline_vcf_variants %>%
				map_int( function(x){x$num_zm_detected %>% sum} ),
			high_confidence_germline_vcf_variants_sum_duplex_coverage = high_confidence_germline_vcf_variants %>%
				map_int( function(x){x$duplex_coverage %>% sum} )
		) %>% 
		select(-high_confidence_germline_vcf_variants)
	
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
				
				(call_class == "indel" &
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
	
}else if(!is.null(sensitivity_parameters$use_chromgroup) & sensitivity_parameters$use_chromgroup != chromgroup_toanalyze){
	cat("## All non-MDB sensitivities set to NA since the currently analyzed chromgroup is not equal to sensitivity_parameters$use_chromgroup.\n")
	sensitivity <- sensitivity %>%
		mutate(
			sensitivity = if_else(call_class == "MDB", sensitivity, NA_real_),
			sensitivity_source = if_else(call_class == "MDB", sensitivity_source, "other_chromgroup")
		)
}else if(is.null(sensitivity_parameters$use_chromgroup)){
	cat("## All non-MDB sensitivities set to 1 since use_chromgroup not defined.\n") #Already set to default of 1 upon above initialization
}

######################
### Calculate estimated SBS mutation error rates per trinucleotide call channel and total error rate
######################
cat("## Calculating estimated SBS mutation error probability per trinucleotide call channel and total error probability...")
#The SBS mutation error probability is estimated for each of the 192 trinucleotide call channels using SBS mismatch-ss calls as [burden(tri-call_i) * burden(rev complement tri-call_i)], where burden(tri-call_i) = [# tri_i calls] / [# interrogated bases with tri_i's trinucleotide context], and likewise for burden(rev complement tri-call_i). Then these error probabilities are summed across the 96 central pyrimidine trinucelotide call channels. The total error probability is then calculated as the sum of the 96 trinucleotide call channel error probabilities. Note: assumes that trinucleotide contexts that are not in interrogated bases have 0 error probability.
estimatedSBSMutationErrorProbability <- list()

#Calculate burden of each trinucleotide call channel
estimatedSBSMutationErrorProbability_by_channel <- left_join(
	#SBS mismatch-ss calls template_strand spectrum
	finalCalls.reftnc_spectra %>%
		filter(call_class=="SBS", SBSindel_call_type=="mismatch-ss") %>%
		pluck("finalCalls.reftnc_template_strand_spectrum",1) %>%
		mutate(reftnc = channel %>% str_sub(1,3)),

	#Interrogated bases spectrum
	bam.gr.filtertrack.bytype %>%
		filter(call_class=="SBS", SBSindel_call_type=="mismatch-ss") %>%
		pluck("bam.gr.filtertrack.reftnc_both_strands",1) %>%
		select(reftnc, count),
	
	by = "reftnc",
	suffix = c(".calls",".interrogated_bases")
) %>%
	mutate(
		burden = count.calls / count.interrogated_bases
	) %>%
	select(channel, burden)

#Remove SBS/mismatch-ss from bam.gr.filtertrack.bytype if it is not part of call_types_toanalyze, since it was only retained until here in order to calculate SBS mutation error probability.
bam.gr.filtertrack.bytype <- bam.gr.filtertrack.bytype %>%
	semi_join(
		call_types_toanalyze,
		by = join_by(call_type, call_class, SBSindel_call_type, filtergroup)
	)

#Multiply each trinucleotide call channel burden by its reverse complement's burden to obtain the error probability for each channel
estimatedSBSMutationErrorProbability$by_channel_pyr <- estimatedSBSMutationErrorProbability_by_channel %>%
	mutate(
		channel_rc = str_c(
			channel %>% str_sub(1,3) %>% DNAStringSet %>% reverseComplement %>% as.character,
			">",
			channel %>% str_sub(5,7) %>% DNAStringSet %>% reverseComplement %>% as.character
		) %>%
			factor(levels = sbs192_labels)
	) %>%
	left_join(
		x = select(., channel, channel_rc, burden),
		y = select(., channel_rc, burden),
		by = join_by(channel == channel_rc),
		suffix = c("",".rc")
	) %>%
	mutate(
		error_prob = burden * burden.rc,
		channel_pyr = if_else(channel %>% str_sub(2,2) %in% c("C","T"), channel, channel_rc) %>%
			factor(levels = sbs96_labels.sigfit)
	) %>%
	group_by(channel_pyr) %>%
	summarize(error_prob = sum(error_prob, na.rm = TRUE)) %>%
	arrange(channel_pyr)

estimatedSBSMutationErrorProbability$total <- estimatedSBSMutationErrorProbability$by_channel_pyr %>%
	pull(error_prob) %>% 
	sum

rm(estimatedSBSMutationErrorProbability_by_channel)

cat("DONE\n")

######################
### Save results
######################
cat("## Saving results...")
qs_save(
	lst(
		yaml.config,
		run_metadata,
		individual_id,
		sample_id = sample_id_toanalyze,
		chromgroup = chromgroup_toanalyze,
		filtergroup = filtergroup_toanalyze,
		call_types = call_types_toanalyze,
		molecule_stats.by_run_id,
		molecule_stats.by_analysis_id,
		region_genome_filter_stats,
		bam.gr.filtertrack.bytype.coverage_tnc = bam.gr.filtertrack.bytype,
		finalCalls,
		finalCalls.bytype,
		germlineVariantCalls,
		finalCalls.reftnc_spectra,
		finalCalls.burdens,
		genome.reftnc,
		genome_chromgroup.reftnc,
		sensitivity,
		estimatedSBSMutationErrorProbability
	),
	outputFile
)

cat("DONE\n")