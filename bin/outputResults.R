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
	make_option(c("-s", "--sample_id"), type = "character", default=NULL,
							help="sample_id to analyze"),
	make_option(c("-f", "--files"), type = "character", default=NULL,
							help="comma-separated calculateBurdens qs2 files"),
	make_option(c("-d", "--duckdbs"), type = "character", default=NULL,
							help="comma-separated coverage_tnc.duckdb files"),
	make_option(c("-o", "--output_basename"), type = "character", default=NULL,
							help="output basename")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$config) | is.null(opt$sample_id) | is.null(opt$files) | is.null(opt$duckdbs) | is.null(opt$output_basename) ){
	stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))
sample_id <- opt$sample_id
calculateBurdensFiles <- opt$files %>% str_split_1(",") %>% str_trim
calculateBurdensDuckDbs <- opt$duckdbs %>% str_split_1(",") %>% str_trim
output_basename <- opt$output_basename

#Load the BSgenome reference
suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name,character.only=TRUE,lib.loc=yaml.config$cache_dir))

#Load miscellaneous configuration parameters
 #analysis_id
analysis_id <- yaml.config$analysis_id

 #individual_id of this sample_id
individual_id <- yaml.config$samples %>%
  bind_rows %>%
  filter(sample_id == sample_id) %>%
  pull(individual_id)

 #call types
call_type <- yaml.config$call_types %>%
	enframe(name=NULL) %>%
	unnest_wider(value) %>%
	unnest_longer(SBSindel_call_types) %>%
	unnest_wider(SBSindel_call_types) %>%
	select(-starts_with("MDB"))

#Display basic configuration parameters
cat("> Processing:\n")
cat("    individual_id:",individual_id,"\n")
cat("    sample_id:",sample_id,"\n")

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
molecule_stats.by_run_id <- list()
molecule_stats.by_analysis_id <- list()
region_genome_filter_stats <- list()
bam.gr.filtertrack.bytype.coverage_tnc <- list()
finalCalls <- list()
germlineVariantCalls <- list()
finalCalls.bytype <- list()
finalCalls.reftnc_spectra <- list()
finalCalls.burdens <- list()
sensitivity <- list()

#Loop over all calculateBurdens
for(i in seq_along(calculateBurdensFiles)){
	
	#Load calculateBurdensFile
	calculateBurdensFile <- qs_read(calculateBurdensFiles[i])
	chromgroup <- calculateBurdensFile$chromgroup
	filtergroup <- calculateBurdensFile$filtergroup
	
	cat(" > chromgroup/filtergroup:", chromgroup, "/", filtergroup, "...")
	
	#Load data that is identical in all calculateBurdensFiles
	if(i == 1){
		#Run metadata
		run_metadata <- calculateBurdensFile %>% pluck("run_metadata")
		
		#Whole genome trinucleotide counts and fractions
		genome.reftnc_pyr <- calculateBurdensFile %>% pluck("genome.reftnc","reftnc_pyr")
		genome.reftnc_both_strands <- calculateBurdensFile %>% pluck("genome.reftnc","reftnc_both_strands")
	}
	
	#Load filtering stats
	molecule_stats.by_run_id[[i]] <- calculateBurdensFile$molecule_stats.by_run_id %>%
		mutate(
			analysis_id = !!analysis_id %>% factor,
			individual_id = !!individual_id %>% factor,
			sample_id = !!sample_id %>% factor,
			.before = 1
		) %>%
		mutate(
			chromgroup = chromgroup %>% factor,
			filtergroup = filtergroup %>% factor
		)
	
	molecule_stats.by_analysis_id[[i]] <- calculateBurdensFile$molecule_stats.by_analysis_id %>%
		mutate(
			analysis_id = !!analysis_id %>% factor,
			individual_id = !!individual_id %>% factor,
			sample_id = !!sample_id %>% factor,
			.before = 1
		) %>%
		mutate(
			chromgroup = chromgroup %>% factor,
			filtergroup = filtergroup %>% factor
		)
	
	region_genome_filter_stats[[i]] <- calculateBurdensFile$region_genome_filter_stats %>%
		mutate(
			analysis_id = !!analysis_id %>% factor,
			individual_id = !!individual_id %>% factor,
			sample_id = !!sample_id %>% factor,
			chromgroup = !!chromgroup %>% factor,
			filtergroup = !!filtergroup %>% factor,
			.before = 1
		)
	
	#Load HiDEF-seq bam genome coverage and trinucleotide counts, fractions, and ratio to genome
	bam.gr.filtertrack.bytype.coverage_tnc[[i]] <- calculateBurdensFile$bam.gr.filtertrack.bytype.coverage_tnc %>%
		mutate(
			analysis_id = !!analysis_id %>% factor,
			individual_id = !!individual_id %>% factor,
			sample_id = !!sample_id %>% factor,
			chromgroup = !!chromgroup %>% factor,
			.before = 1
		) %>%
		relocate(filtergroup, .after = chromgroup) %>%
		select(-analyzein_chromgroups)
	
	#Final calls for single table output
	finalCalls[[i]] <- calculateBurdensFile$finalCalls %>%
		mutate(
			analysis_id = !!analysis_id %>% factor,
			individual_id = !!individual_id %>% factor,
			sample_id = !!sample_id %>% factor,
			chromgroup = !!chromgroup %>% factor,
			filtergroup = !!filtergroup %>% factor
			.before = 1
		)
	
	#Germline variant calls
	germlineVariantCalls[[i]] <- calculateBurdensFile$germlineVariantCalls %>%
		mutate(
			analysis_id = !!analysis_id %>% factor,
			individual_id = !!individual_id %>% factor,
			sample_id = !!sample_id %>% factor,
			chromgroup = !!chromgroup %>% factor,
			filtergroup = !!filtergroup %>% factor
			.before = 1
		)
	
	#finalCalls for tsv and VCF output
	finalCalls.bytype[[i]] <- calculateBurdensFile$finalCalls.bytype %>%
		mutate(
			analysis_id = !!analysis_id %>% factor,
			individual_id = !!individual_id %>% factor,
			sample_id = !!sample_id %>% factor,
			chromgroup = !!chromgroup %>% factor,
			.before = 1
		) %>%
		relocate(filtergroup, .after = chromgroup) %>%
		select(-analyzein_chromgroups)
	
	#Trinucleotide context counts and fractions of finalCalls
	finalCalls.reftnc_spectra[[i]] <- calculateBurdensFile$finalCalls.reftnc_spectra %>%
		mutate(
			analysis_id = !!analysis_id %>% factor,
			individual_id = !!individual_id %>% factor,
			sample_id = !!sample_id %>% factor,
			chromgroup = !!chromgroup %>% factor,
			.before = 1
		) %>%
		relocate(filtergroup, .after = chromgroup) %>%
		select(-analyzein_chromgroups)
	
	#Burdens
	finalCalls.burdens[[i]] <- calculateBurdensFile$finalCalls.burdens %>%
		mutate(
			analysis_id = !!analysis_id %>% factor,
			individual_id = !!individual_id %>% factor,
			sample_id = !!sample_id %>% factor,
			chromgroup = !!chromgroup %>% factor,
			.before = 1
		) %>%
		relocate(filtergroup, .after = chromgroup) %>%
		select(-analyzein_chromgroups)
	
	#Sensitivity
	if(!is.null(calculateBurdensFile$sensitivity)){
		sensitivity[[i]] <- calculateBurdensFile$sensitivity %>%
			mutate(
				analysis_id = !!analysis_id %>% factor,
				individual_id = !!individual_id %>% factor,
				sample_id = !!sample_id %>% factor,
				chromgroup = !!chromgroup %>% factor,
				.before = 1
			) %>%
			relocate(filtergroup, .after = chromgroup) %>%
			select(-analyzein_chromgroups)
	}
	
	#estiamted SBS mutation error probability
	if(!is.null(calculateBurdensFile$estimatedSBSMutationErrorProbability)){
		estimatedSBSMutationErrorProbability[[i]] <- calculateBurdensFile$estimatedSBSMutationErrorProbability
	}
	
	#Remove temporary objects
	rm(calculcateBurdensFile)
	invisible(gc())
	
	cat("DONE\n")

}

#Combine data across chromgroups/filtergroups
molecule_stats.by_run_id <-
molecule_stats.by_analysis_id <-
region_genome_filter_stats <-
bam.gr.filtertrack.bytype.coverage_tnc <-
finalCalls <-
germlineVariantCalls <-
finalCalls.bytype <-
finalCalls.reftnc_spectra <-
finalCalls.burdens <-
sensitivity <-

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

molecule_stats_by_run_id %>% write_tsv(str_c(output_basename,".molecule_stats_by_run_id.tsv"))

molecule_stats_by_analysis_id %>% write_tsv(str_c(output_basename,".molecule_stats_by_analysis_id.tsv"))

region_genome_filter_stats %>% write_tsv(str_c(output_basename,".region_genome_filter_stats.tsv"))

cat("DONE\n")


### output coverage and reference sequences of interrogated genome bases
#Output, with a separate folder for each call_class, call_type, SBSindel_call_type combination

**con %>% dbExecute("CREATE INDEX idx_cov_perbase_rowid ON coverage_perbase(row_id);") %>% invisible

calculateBurdensDuckDbs

Order by rowid:
	con %>% dbExecute("
  CREATE OR REPLACE TABLE coverage_runs_sorted AS
  SELECT * FROM coverage_runs
  ORDER BY row_id, seqnames, start;
  DROP TABLE coverage_runs;
  ALTER TABLE coverage_runs_sorted RENAME TO coverage_runs;
  CHECKPOINT;
") %>% invisible

Then materialize a temp table for each rowid, then do the jointo genome_trinuc and stream to BED

#Attach previously created genome trinucleotide duckdb
genome_trinuc_duckdb_file <- str_c(yaml.config$cache_dir,"/",yaml.config$BSgenome$BSgenome_name,".trinuc.duckdb")

con %>%
	dbExecute(
		sprintf(
			"ATTACH '%s' AS genome_trinuc (READ_ONLY);",
			genome_trinuc_duckdb_file
		)
	) %>%
	invisible

DBI::dbExecute(
	con,
	sprintf("
    COPY (
      SELECT
        t.seqnames               AS chrom,           -- col1
        t.pos - 1                AS start0,          -- col2 (BED start)
        t.pos                    AS end1,            -- col3 (BED end)
        t.reftnc                 AS reftnc,          -- col4
        '.'                      AS dot,             -- col5
        r.duplex_coverage        AS duplex_coverage  -- col6
      FROM coverage_runs r
      JOIN genome_trinuc.reftnc_plus_strand t
        ON t.seqnames = r.seqnames
       AND t.pos BETWEEN r.start AND r.end_pos
      ORDER BY t.seqnames, t.pos
    )
    TO '%s' (DELIMITER '\t', HEADER FALSE);
  ", out_tsv)
)

**stream output to bgzip, then do tabix - check that this actually streams.
con <- dbConnect(duckdb::duckdb(), "mydb.duckdb")

cmd <- pipe("gzip > out.tsv.gz", "w")
dbExecute(con, "COPY (SELECT * FROM my_table) TO '/dev/stdout' (DELIMITER '\t', HEADER FALSE);", immediate=TRUE, resultFormat="raw") %>%
	writeBin(cmd)
close(cmd)

for(i in seq_len(nrow(bam.gr.filtertrack.bytype))){
	
	prefix <- 
		str_c(
			bam.gr.filtertrack.bytype$call_class[i],
			bam.gr.filtertrack.bytype$call_type[i],
			bam.gr.filtertrack.bytype$SBSindel_call_type[i],
			sep="."
		)
	
	dir.create(prefix, showWarnings = FALSE)
	
	bam.gr.filtertrack.bytype %>%
		pluck("bam.gr.filtertrack.coverage",i) %>%
		** Fix to output each chromosome sequentially, since it is now a List
		mutate(name = reftnc_pyr, score = coverage) %>% #rename to to allow output by 'export'
		export(
			con = str_c(
				str_c(prefix,"/",output_basename),
				prefix,
				"bed",
				sep="."
			),
			format = "bed",
			index = TRUE
		)
}

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