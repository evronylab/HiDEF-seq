######################
### Custom shared functions
######################

#Function to derive the BSgenome package name created from the configured FASTA
get_bsgenome_name <- function(yaml.config){
	genome_suffix <- yaml.config$genome_fasta %>%
		basename %>%
		str_replace_all(pattern = "[^0-9a-zA-Z.]", replacement = "")

	str_c("BSgenome", yaml.config$genome_organism, "user", genome_suffix, sep=".")
}

#Function to load and format VCF
# regions: optional tibble with columns: seqnames start_refspace end_refspace
load_vcf <- function(vcf_file, regions = NULL, genome_fasta, BSgenome_name, bcftools_bin){

  #Check if AD, GT, and GQ tags exist
	AD_exists <- system(paste("/bin/bash -c",shQuote(paste(
		bcftools_bin,"view -h",
		vcf_file,"|",
		"grep '^##FORMAT=<ID=AD' | wc -l"
	)
	)), intern=TRUE) != "0"

	GT_exists <- system(paste("/bin/bash -c",shQuote(paste(
		bcftools_bin,"view -h",
		vcf_file,"|",
		"grep '^##FORMAT=<ID=GT' | wc -l"
		)
	 )), intern=TRUE) != "0"
	
	GQ_exists <- system(paste("/bin/bash -c",shQuote(paste(
		bcftools_bin,"view -h",
		vcf_file,"|",
		"grep '^##FORMAT=<ID=GQ' | wc -l"
	)
	)), intern=TRUE) != "0"
	
	if(!is.null(regions) && nrow(regions) == 0){
		return(
			GRanges(
				ref_plus_strand = character(),
				alt_plus_strand = character(),
				QUAL = numeric(),
				FILTER = character(),
				GT = character(),
				GQ = numeric(),
				AD1 = integer(),
				AD2 = integer(),
				Depth = integer(),
				VAF = numeric(),
				call_class = factor(),
				call_type = factor(),
				SBSindel_call_type = factor()
			)
		)
	}
	
	#Create regions files if specified
	if(!is.null(regions)){
		tmpregions <- tempfile(tmpdir=getwd(),pattern=".")
		regions %>%
			arrange(seqnames,start_refspace,end_refspace) %>%
			distinct %>%
			write_tsv(tmpregions, col_names=FALSE)
	}
	
	#Atomize and split multi-allelic sites (bcftools norm -a -f [fastaref] | bcftools norm -m -both -f [fastaref]), and filter for records containing ALT alleles.
  tmpvcf <- tempfile(tmpdir=getwd(),pattern=".")
  
  system(paste("/bin/bash -c",shQuote(paste(
    bcftools_bin,"view",
    if(!is.null(regions)){paste("-R",tmpregions)},
    vcf_file,"|",
    bcftools_bin,"norm -a -f",genome_fasta,"2>/dev/null |",
    bcftools_bin,"norm -m -both -f",genome_fasta,"2>/dev/null",
    if(GT_exists){paste("|",bcftools_bin,"view -i 'GT=\"alt\"'")},
    ">",
    tmpvcf
    )
   )))
  
  if(!is.null(regions)){
  	file.remove(tmpregions) %>% invisible
  }

  #Load vcf file
   #Remove ALT == "*" alleles that indicate overlapping deletions (not needed because every deletion already has a separate vcf entry; seen in germline VCF files) and ALT = "<*>" alleles that are non-variant gVCF blocks (seen in bcftools mpileup files).
   #If the AD tag exists, annotate with allele depths (AD1, AD2), Depth (AD1+AD2) and VAF (AD2 / Depth).
   #Annotate SBS vs insertions vs deletions
   #For deletions, change ranges to reflect position of deletion. For insertions: change start position to the base on the right of the insertion position and end to the base on the left of the insertion position.
   #For deletions, change to: REF = deleted bases and ALT = "". For insertions, change to: REF = "" and ALT = inserted bases
   #Keep only columns CHROM, POS, REF, ALT, QUAL, FILTER, GT (if exists), GQ (if exists).
   #Convert to GRanges
  
  vcf <- read.vcfR(tmpvcf,convertNA=FALSE,verbose=FALSE)
  file.remove(tmpvcf) %>% invisible
  
  vcf <- vcf@fix %>%
  	as_tibble %>%
  	mutate(
  		GT = if(GT_exists){
  			as.character(extract.gt(vcf,element="GT",IDtoRowNames=FALSE))
  		}else{NULL},
  		GQ = if(GQ_exists){
  			as.numeric(extract.gt(vcf,element="GQ",as.numeric=TRUE,IDtoRowNames=FALSE))
  		}else{NULL},
  		QUAL = if_else(QUAL == ".",NA,QUAL)
  	) %>%
  	{
  		if(AD_exists){
	  			mutate(
	  				.,
	  				AD = as.character(extract.gt(vcf,element="AD",IDtoRowNames=FALSE))
					)
  		}else{
  			.
  		}
  	} %>%
  	filter(! ALT %in% c("*","<*>")) %>%
  	{
  		if(AD_exists){
  				separate(.,AD,c("AD1","AD2"),",",convert=TRUE) %>%
  				mutate(Depth=AD1+AD2, VAF=AD2/Depth)
  		}else{
  			.
  		}
  	} %>%
    mutate(
      call_class = if_else(nchar(REF)==1 & nchar(ALT)==1,"SBS","indel") %>% factor,
      call_type = case_when(
        call_class == "SBS" ~ "SBS",
        call_class == "indel" & (nchar(REF) - nchar(ALT) < 0) ~ "insertion",
        call_class == "indel" & (nchar(REF) - nchar(ALT) > 0) ~ "deletion"
        ) %>%
        factor,
      SBSindel_call_type = "mutation" %>% factor,
      CHROM = factor(CHROM),
      POS = as.numeric(POS),
      start_refspace = if_else(call_type %in% c("deletion","insertion"), POS + 1, POS),
      end_refspace = if_else(call_type == "deletion", POS + nchar(REF) - nchar(ALT), POS),
      REF = case_when(
        call_type == "insertion" ~ "",
        call_type == "deletion" ~ str_sub(REF,2),
        .default = REF
        ),
      ALT = case_when(
        call_type == "insertion" ~ str_sub(ALT,2),
        call_type == "deletion" ~ "",
        .default = ALT
      ),
      QUAL = as.numeric(QUAL)
    ) %>%
    select(-c(POS,ID,INFO)) %>%
    rename(
      seqnames = CHROM,
      ref_plus_strand = REF,
      alt_plus_strand = ALT
    ) %>%
  	distinct %>% # select distinct rows since -R regions can lead to duplicated entries
    makeGRangesFromDataFrame(
      seqnames.field="seqnames",
      start.field="start_refspace",
      end.field="end_refspace",
      keep.extra.columns=TRUE,
      seqinfo=BSgenome_name %>% get %>% seqinfo
    ) %>%
  	sort

  return(vcf)
}

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
GRanges_subtract_bymcols <- function(x, y, join_mcols, ignore.strand = FALSE){
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
  
  keys_all <- factor(c(key_x,key_y))
  id_x <- as.integer(keys_all)[seq_along(x)]
  id_y <- as.integer(keys_all)[length(x) + seq_along(y)]
  
  #Change seqnames to contain key
  x2 <- GRanges(
  	seqnames = str_c(as.character(seqnames(x)), "|", id_x),
  	ranges = ranges(x),
  	strand = strand(x)
  )
  
  y2 <- GRanges(
  	seqnames = str_c(as.character(seqnames(y)), "|", id_y),
  	ranges = ranges(y),
  	strand = strand(y)
  ) %>%
  	GNCList
  
  #Find all overlaps
  hits <- findOverlaps(x2, y2, ignore.strand = ignore.strand) %>%
  	suppressWarnings #Hide warning about unshared seqlevels between x2 and y2
  
  if(length(hits) == 0){return(x)}
  
  #Obtain intersecting segments for matched hits
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  
  segs_to_remove <- pintersect(x[qh], y[sh], ignore.strand = ignore.strand)
  
  #Split segs_to_remove by the index of x
  rem_list <- split(segs_to_remove, factor(qh, levels = seq_along(x)), drop = FALSE)
  
  #Subtract for each range in x (x_list) the ranges in y (rem_list) that intersect it
  out_list <- psetdiff(x, rem_list, ignore.strand = ignore.strand)
  
  #Flatten result and restore all metadata from the original x
   #How many fragments came from each original x[i]
  nfrags <- elementNROWS(out_list)
   #Unlist into a single GRanges
  out <- unlist(out_list, use.names = FALSE)
   #Build an index mapping each fragment back to its source x[i]
  src_idx  <- rep(seq_along(nfrags), times = nfrags)
   #Re‐attach metadata columns
  mcols(out) <- mcols(x)[src_idx, , drop = FALSE]
  
  return(out)
}

#All possible trinucleotides
trinucleotides_64 <- mkAllStrings(c("A","C","G","T"), 3)

#central pyrimidine trinucleotides
trinucleotides_32_pyr <- trinucleotides_64 %>% str_subset(".[CT].")

#List named with trinucleotides_32_pyr sequences, each containing the corresponding trinucleotides_64 sequences
trinucleotides_64_32_pyr_list <- split(
	trinucleotides_64,
	ifelse(
		str_sub(trinucleotides_64, 2, 2) %in% c("C","T"),
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
		) %>%
		group_by(!!sym(tri_column)) %>%
		summarize(!!sym(count_column) := sum(!!sym(count_column)), .groups = "drop") %>%
		arrange(!!sym(tri_column))
}

#Function to re-anchor indels from HiDEF-seq calls to the style of VCF format so that ref_plus_strand and alt_plus_strand contain the base preceding the indel
normalize_indels_for_vcf <- function(df, BSgenome_name){
	
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

#Function to import indels for spectrum analysis, modified from INDELWALD package:
## Max Stammnitz; maxrupsta@gmail.com; University of Cambridge  ##
## Citation: The evolution of two transmissible cancers in Tasmanian devils (Stammnitz et al. 2023, Science 380:6642)
indel.spectrum <- function(x, reference, context_bp = 1000, spectrum_type = c("pyr", "template")){

	spectrum_type <- match.arg(spectrum_type)
	one_bp_bases <- if(spectrum_type == "pyr") c("C", "T") else c("C", "T", "G", "A")
	
	# context_bp: total number of bp extracted around ≥5 bp events (flank on each side ~ context_bp/2)
	.flank <- as.integer(context_bp / 2L)  # symmetric flank size used in ≥5 bp sections
	
	## 1. split VCF indels into:
	# (i) 1-bp deletions
	dels.1bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-1,,drop=F]
	
	# (ii) 1-bp insertions
	ins.1bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-1,,drop=F]
	
	# (iii) >1 bp deletions
	
	## 2
	dels.2bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-2,,drop=F]
	
	## 3
	dels.3bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-3,,drop=F]
	
	## 4
	dels.4bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-4,,drop=F]
	
	## 5+
	dels.5bp <- x[nchar(as.character(x[,'ALT'])) <= nchar(as.character(x[,'REF']))-5,,drop=F]
	
	# (iv) >1 bp insertions
	
	## 2
	ins.2bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-2,,drop=F]
	
	## 3
	ins.3bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-3,,drop=F]
	
	## 4
	ins.4bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-4,,drop=F]
	
	## 5+
	ins.5bp <- x[nchar(as.character(x[,'REF'])) <= nchar(as.character(x[,'ALT']))-5,,drop=F]
	
	## 2. classify 1 bp events into:
	
	# (i) 1 bp deletions at homopolymers (length 1 == "no neighbouring homopolymers")
	if(nrow(dels.1bp) > 0){
		
		## extract 10 bp upstream/downstream sequence context from reference
		dels.1bp.context <- as.character(subseq(x = reference[as.character(dels.1bp[,'CHROM'])], 
																						start = as.numeric(dels.1bp[,'POS']) - 9, 
																						end = as.numeric(dels.1bp[,'POS']) + 11))
		dels.1bp.context.middle <- paste0('[', str_split_fixed(dels.1bp.context, '', 21)[,11,drop=F], ']')
		dels.1bp.context.start <- str_split_fixed(dels.1bp.context, '', 11)[,1:10,drop=F]
		dels.1bp.context.start <- paste(dels.1bp.context.start[,1], dels.1bp.context.start[,2], dels.1bp.context.start[,3],
																		dels.1bp.context.start[,4], dels.1bp.context.start[,5], dels.1bp.context.start[,6],
																		dels.1bp.context.start[,7], dels.1bp.context.start[,8], dels.1bp.context.start[,9],
																		dels.1bp.context.start[,10], sep = '')
		dels.1bp.context.end <- str_split_fixed(dels.1bp.context, '', 12)[,12,drop=F]
		dels.1bp[,'CONTEXT FW'] <- paste(dels.1bp.context.start, dels.1bp.context.middle, dels.1bp.context.end, sep = '')
		
		## For pyr spectra, pyrimidine-center 1-bp indels; for template spectra, orient by TEMPLATE_STRAND.
		dels.1bp.context.rc <- as.character(reverseComplement(subseq(x = reference[as.character(dels.1bp[,'CHROM'])], 
																																 start = as.numeric(dels.1bp[,'POS']) - 9, 
																																 end = as.numeric(dels.1bp[,'POS']) + 11)))
		dels.1bp.context.rc.middle <- paste0('[', str_split_fixed(dels.1bp.context.rc, '', 21)[,11,drop=F], ']')
		dels.1bp.context.rc.start <- str_split_fixed(dels.1bp.context.rc, '', 11)[,1:10,drop=F]
		dels.1bp.context.rc.start <- paste(dels.1bp.context.rc.start[,1], dels.1bp.context.rc.start[,2], dels.1bp.context.rc.start[,3],
																			 dels.1bp.context.rc.start[,4], dels.1bp.context.rc.start[,5], dels.1bp.context.rc.start[,6],
																			 dels.1bp.context.rc.start[,7], dels.1bp.context.rc.start[,8], dels.1bp.context.rc.start[,9],
																			 dels.1bp.context.rc.start[,10], sep = '')
		dels.1bp.context.rc.end <- str_split_fixed(dels.1bp.context.rc, '', 12)[,12,drop=F]
		dels.1bp[,'CONTEXT RC'] <- paste(dels.1bp.context.rc.start, dels.1bp.context.rc.middle, dels.1bp.context.rc.end, sep = '')
		
		## summarise 1 bp deletions in matrix format
		dels.1bp.summary <- matrix(0, ncol = length(one_bp_bases), nrow = 6)
		colnames(dels.1bp.summary) <- one_bp_bases
		rownames(dels.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
		fw <- str_split_fixed(str_split_fixed(dels.1bp[,'CONTEXT FW'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		rc <- str_split_fixed(str_split_fixed(dels.1bp[,'CONTEXT RC'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		fw.if <- if(spectrum_type == "pyr"){
			fw == 'C' | fw == 'T'
		}else{
			dels.1bp[,'TEMPLATE_STRAND'] == "+"
		}
		for (i in 1:nrow(dels.1bp)){
			
			if(fw.if[i] == T){
				
				## look at forward context
				upstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT FW']), '\\[', 2)[,1], '', 10)
				upstream.tmp <- which(upstream.tmp == fw[i])
				
				### check all 6 categories for upstream bases
				if(any(upstream.tmp %in% 10)){
					
					if(any(upstream.tmp %in% 9)){
						
						if(any(upstream.tmp %in% 8)){
							
							if(any(upstream.tmp %in% 7)){
								
								if(any(upstream.tmp %in% 6)){
									
									upstream.tmp <- 6
									
								}else{
									upstream.tmp <- 5
								}
								
							}else{
								upstream.tmp <- 4
							}
							
						}else{
							upstream.tmp <- 3
						}
						
					}else{
						upstream.tmp <- 2
					}
					
				}else{
					upstream.tmp <- 1
				}
				
				downstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT FW']), '\\]', 2)[,2], '', 10)
				downstream.tmp <- which(downstream.tmp == fw[i])
				
				### check all 6 categories for downstream bases
				if(any(downstream.tmp %in% 1)){
					
					if(any(downstream.tmp %in% 2)){
						
						if(any(downstream.tmp %in% 3)){
							
							if(any(downstream.tmp %in% 4)){
								
								if(any(downstream.tmp %in% 5)){
									
									downstream.tmp <- 6
									
								}else{
									downstream.tmp <- 5
								}
								
							}else{
								downstream.tmp <- 4
							}
							
						}else{
							downstream.tmp <- 3
						}
						
					}else{
						downstream.tmp <- 2
					}
					
				}else{
					downstream.tmp <- 1
				}
				
				## summarise, which homopolymer (upstream vs. downstream context) is longer
				dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] <- dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] + 1
				
			} else {
				
				## look at reverse complement context
				upstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT RC']), '\\[', 2)[,1], '', 10)
				upstream.tmp <- which(upstream.tmp == rc[i])
				
				### check all 6 categories for upstream bases
				if(any(upstream.tmp %in% 10)){
					
					if(any(upstream.tmp %in% 9)){
						
						if(any(upstream.tmp %in% 8)){
							
							if(any(upstream.tmp %in% 7)){
								
								if(any(upstream.tmp %in% 6)){
									
									upstream.tmp <- 6
									
								}else{
									upstream.tmp <- 5
								}
								
							}else{
								upstream.tmp <- 4
							}
							
						}else{
							upstream.tmp <- 3
						}
						
					}else{
						upstream.tmp <- 2
					}
					
				}else{
					upstream.tmp <- 1
				}
				
				downstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT RC']), '\\]', 2)[,2], '', 10)
				downstream.tmp <- which(downstream.tmp == rc[i])
				
				### check all 6 categories for downstream bases
				if(any(downstream.tmp %in% 1)){
					
					if(any(downstream.tmp %in% 2)){
						
						if(any(downstream.tmp %in% 3)){
							
							if(any(downstream.tmp %in% 4)){
								
								if(any(downstream.tmp %in% 5)){
									
									downstream.tmp <- 6
									
								}else{
									downstream.tmp <- 5
								}
								
							}else{
								downstream.tmp <- 4
							}
							
						}else{
							downstream.tmp <- 3
						}
						
					}else{
						downstream.tmp <- 2
					}
					
				}else{
					downstream.tmp <- 1
				}
				
				## summarise, which homopolymer (upstream vs. downstream context) is longer
				dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] <- dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] + 1
			}
			
		}
		
	}else{
		
		## summarise 1 bp deletions in matrix format
		dels.1bp.summary <- matrix(0, ncol = length(one_bp_bases), nrow = 6)
		colnames(dels.1bp.summary) <- one_bp_bases
		rownames(dels.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
		
	}
	
	# (ii) 1 bp insertions at homopolymers (length 0 == "no neighbouring homopolymers")
	if(nrow(ins.1bp) > 0){
		
		## extract 10 bp upstream/downstream sequence context from reference
		ins.1bp.context <- as.character(subseq(x = reference[as.character(ins.1bp[,'CHROM'])], 
																					 start = as.numeric(ins.1bp[,'POS']) - 9, 
																					 end = as.numeric(ins.1bp[,'POS']) + 10))
		## insertion after base pos. 10
		ins.1bp.context.middle <- paste0('[', str_split_fixed(as.character(ins.1bp[,'ALT']), '', 2)[,2,drop=F], ']')
		ins.1bp.context.start <- str_split_fixed(ins.1bp.context, '', 11)[,1:10,drop=F]
		ins.1bp.context.start <- paste(ins.1bp.context.start[,1], ins.1bp.context.start[,2], ins.1bp.context.start[,3],
																	 ins.1bp.context.start[,4], ins.1bp.context.start[,5], ins.1bp.context.start[,6],
																	 ins.1bp.context.start[,7], ins.1bp.context.start[,8], ins.1bp.context.start[,9],
																	 ins.1bp.context.start[,10], sep = '')
		ins.1bp.context.end <- str_split_fixed(ins.1bp.context, '', 11)[,11,drop=F]
		ins.1bp[,'CONTEXT FW'] <- paste(ins.1bp.context.start, ins.1bp.context.middle, ins.1bp.context.end, sep = '')
		
		## For pyr spectra, pyrimidine-center 1-bp indels; for template spectra, orient by TEMPLATE_STRAND.
		ins.1bp.context.rc <- as.character(reverseComplement(subseq(x = reference[as.character(ins.1bp[,'CHROM'])], 
																																start = as.numeric(ins.1bp[,'POS']) - 9, 
																																end = as.numeric(ins.1bp[,'POS']) + 10)))
		ins.1bp.context.rc.middle <- paste0('[', as.character(reverseComplement(DNAStringSet(str_split_fixed(as.character(ins.1bp[,'ALT']), '', 2)[,2,drop=F]))), ']')
		ins.1bp.context.rc.start <- str_split_fixed(ins.1bp.context.rc, '', 11)[,1:10,drop=F]
		ins.1bp.context.rc.start <- paste(ins.1bp.context.rc.start[,1], ins.1bp.context.rc.start[,2], ins.1bp.context.rc.start[,3],
																			ins.1bp.context.rc.start[,4], ins.1bp.context.rc.start[,5], ins.1bp.context.rc.start[,6],
																			ins.1bp.context.rc.start[,7], ins.1bp.context.rc.start[,8], ins.1bp.context.rc.start[,9],
																			ins.1bp.context.rc.start[,10], sep = '')
		ins.1bp.context.rc.end <- str_split_fixed(ins.1bp.context.rc, '', 11)[,11,drop=F]
		ins.1bp[,'CONTEXT RC'] <- paste(ins.1bp.context.rc.start, ins.1bp.context.rc.middle, ins.1bp.context.rc.end, sep = '')
		
		## summarise 1 bp insertions in matrix format
		ins.1bp.summary <- matrix(0, ncol = length(one_bp_bases), nrow = 6)
		colnames(ins.1bp.summary) <- one_bp_bases
		rownames(ins.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
		fw <- str_split_fixed(str_split_fixed(ins.1bp[,'CONTEXT FW'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		rc <- str_split_fixed(str_split_fixed(ins.1bp[,'CONTEXT RC'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		fw.if <- if(spectrum_type == "pyr"){
			fw == 'C' | fw == 'T'
		}else{
			ins.1bp[,'TEMPLATE_STRAND'] == "+"
		}
		for (i in 1:nrow(ins.1bp)){
			
			if(fw.if[i] == T){
				
				## look at forward context
				upstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT FW']), '\\[', 2)[,1], '', 10)
				upstream.tmp <- which(upstream.tmp == fw[i])
				
				### check all 6 categories for upstream bases
				if(any(upstream.tmp %in% 10)){
					
					if(any(upstream.tmp %in% 9)){
						
						if(any(upstream.tmp %in% 8)){
							
							if(any(upstream.tmp %in% 7)){
								
								if(any(upstream.tmp %in% 6)){
									
									upstream.tmp <- 6
									
								}else{
									upstream.tmp <- 5
								}
								
							}else{
								upstream.tmp <- 4
							}
							
						}else{
							upstream.tmp <- 3
						}
						
					}else{
						upstream.tmp <- 2
					}
					
				}else{
					upstream.tmp <- 1
				}
				
				downstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT FW']), '\\]', 2)[,2], '', 10)
				downstream.tmp <- which(downstream.tmp == fw[i])
				
				### check all 6 categories for downstream bases
				if(any(downstream.tmp %in% 1)){
					
					if(any(downstream.tmp %in% 2)){
						
						if(any(downstream.tmp %in% 3)){
							
							if(any(downstream.tmp %in% 4)){
								
								if(any(downstream.tmp %in% 5)){
									
									downstream.tmp <- 6
									
								}else{
									downstream.tmp <- 5
								}
								
							}else{
								downstream.tmp <- 4
							}
							
						}else{
							downstream.tmp <- 3
						}
						
					}else{
						downstream.tmp <- 2
					}
					
				}else{
					downstream.tmp <- 1
				}
				
				## summarise, which homopolymer is longer
				ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] <- ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] + 1
				
			} else {
				
				## look at reverse complement context
				upstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT RC']), '\\[', 2)[,1], '', 10)
				upstream.tmp <- which(upstream.tmp == rc[i])
				
				### check all 6 categories for upstream bases
				if(any(upstream.tmp %in% 10)){
					
					if(any(upstream.tmp %in% 9)){
						
						if(any(upstream.tmp %in% 8)){
							
							if(any(upstream.tmp %in% 7)){
								
								if(any(upstream.tmp %in% 6)){
									
									upstream.tmp <- 6
									
								}else{
									upstream.tmp <- 5
								}
								
							}else{
								upstream.tmp <- 4
							}
							
						}else{
							upstream.tmp <- 3
						}
						
					}else{
						upstream.tmp <- 2
					}
					
				}else{
					upstream.tmp <- 1
				}
				
				downstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT RC']), '\\]', 2)[,2], '', 10)
				downstream.tmp <- which(downstream.tmp == rc[i])
				
				### check all 6 categories for downstream bases
				if(any(downstream.tmp %in% 1)){
					
					if(any(downstream.tmp %in% 2)){
						
						if(any(downstream.tmp %in% 3)){
							
							if(any(downstream.tmp %in% 4)){
								
								if(any(downstream.tmp %in% 5)){
									
									downstream.tmp <- 6
									
								}else{
									downstream.tmp <- 5
								}
								
							}else{
								downstream.tmp <- 4
							}
							
						}else{
							downstream.tmp <- 3
						}
						
					}else{
						downstream.tmp <- 2
					}
					
				}else{
					downstream.tmp <- 1
				}
				
				## summarise, which homopolymer is longer
				ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] <- ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] + 1
			}
			
		} 
		
	}else{
		
		## summarise 1 bp insertions in matrix format
		ins.1bp.summary <- matrix(0, ncol = length(one_bp_bases), nrow = 6)
		colnames(ins.1bp.summary) <- one_bp_bases
		rownames(ins.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
		
	}
	
	## 3. classify >=2 bp deletions into:
	
	# (i) 2 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
	if(nrow(dels.2bp) > 0){
		
		## extract 5 x 2 bp upstream/downstream sequence context from reference
		dels.2bp.context <- as.character(subseq(x = reference[as.character(dels.2bp[,'CHROM'])], 
																						start = as.numeric(dels.2bp[,'POS']) - 9, 
																						end = as.numeric(dels.2bp[,'POS']) + 2 + 10))
		dels.2bp.context.middle <- str_split_fixed(dels.2bp.context, '', 22)[,11:12,drop=F]
		dels.2bp.context.middle <- paste0('[', dels.2bp.context.middle[,1], dels.2bp.context.middle[,2], ']')
		dels.2bp.context.start <- str_split_fixed(dels.2bp.context, '', 11)[,1:10,drop=F]
		dels.2bp.context.start <- paste(dels.2bp.context.start[,1], dels.2bp.context.start[,2], dels.2bp.context.start[,3],
																		dels.2bp.context.start[,4], dels.2bp.context.start[,5], dels.2bp.context.start[,6],
																		dels.2bp.context.start[,7], dels.2bp.context.start[,8], dels.2bp.context.start[,9],
																		dels.2bp.context.start[,10], sep = '')
		dels.2bp.context.end <- str_split_fixed(dels.2bp.context, '', 13)[,13,drop=F]
		dels.2bp[,'CONTEXT'] <- paste(dels.2bp.context.start, dels.2bp.context.middle, dels.2bp.context.end, sep = '')
		
	}
	
	# (ii) 3 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
	if(nrow(dels.3bp) > 0){
		
		## extract 5 x 3 bp upstream/downstream sequence context from reference
		dels.3bp.context <- as.character(subseq(x = reference[as.character(dels.3bp[,'CHROM'])], 
																						start = as.numeric(dels.3bp[,'POS']) - 14, 
																						end = as.numeric(dels.3bp[,'POS']) + 3 + 15))
		dels.3bp.context.middle <- str_split_fixed(dels.3bp.context, '', 33)[,16:18,drop=F]
		dels.3bp.context.middle <- paste0('[', dels.3bp.context.middle[,1], dels.3bp.context.middle[,2], dels.3bp.context.middle[,3], ']')
		dels.3bp.context.start <- str_split_fixed(dels.3bp.context, '', 33)[,1:15,drop=F]
		dels.3bp.context.start <- paste(dels.3bp.context.start[,1], dels.3bp.context.start[,2], dels.3bp.context.start[,3],
																		dels.3bp.context.start[,4], dels.3bp.context.start[,5], dels.3bp.context.start[,6],
																		dels.3bp.context.start[,7], dels.3bp.context.start[,8], dels.3bp.context.start[,9],
																		dels.3bp.context.start[,10], dels.3bp.context.start[,11], dels.3bp.context.start[,12], 
																		dels.3bp.context.start[,13], dels.3bp.context.start[,14], dels.3bp.context.start[,15], sep = '')
		dels.3bp.context.end <- str_split_fixed(dels.3bp.context, '', 19)[,19,drop=F]
		dels.3bp[,'CONTEXT'] <- paste(dels.3bp.context.start, dels.3bp.context.middle, dels.3bp.context.end, sep = '')
		
	}
	
	# (iii) 4 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
	if(nrow(dels.4bp) > 0){
		
		## extract 5 x 4 bp upstream/downstream sequence context from reference
		dels.4bp.context <- as.character(subseq(x = reference[as.character(dels.4bp[,'CHROM'])], 
																						start = as.numeric(dels.4bp[,'POS']) - 19, 
																						end = as.numeric(dels.4bp[,'POS']) + 4 + 20))
		dels.4bp.context.middle <- str_split_fixed(dels.4bp.context, '', 41)[,21:24,drop=F]
		dels.4bp.context.middle <- paste0('[', dels.4bp.context.middle[,1], 
																			dels.4bp.context.middle[,2], 
																			dels.4bp.context.middle[,3], 
																			dels.4bp.context.middle[,4], ']')
		dels.4bp.context.start <- str_split_fixed(dels.4bp.context, '', 21)[,1:20,drop=F]
		dels.4bp.context.start <- paste(dels.4bp.context.start[,1], dels.4bp.context.start[,2], dels.4bp.context.start[,3],
																		dels.4bp.context.start[,4], dels.4bp.context.start[,5], dels.4bp.context.start[,6],
																		dels.4bp.context.start[,7], dels.4bp.context.start[,8], dels.4bp.context.start[,9],
																		dels.4bp.context.start[,10], dels.4bp.context.start[,11], dels.4bp.context.start[,12], 
																		dels.4bp.context.start[,13], dels.4bp.context.start[,14], dels.4bp.context.start[,15], 
																		dels.4bp.context.start[,16], dels.4bp.context.start[,17], dels.4bp.context.start[,18], 
																		dels.4bp.context.start[,19], dels.4bp.context.start[,20], sep = '')
		dels.4bp.context.end <- str_split_fixed(dels.4bp.context, '', 25)[,25,drop=F]
		dels.4bp[,'CONTEXT'] <- paste(dels.4bp.context.start, dels.4bp.context.middle, dels.4bp.context.end, sep = '')
		
	}
	
	# (iv) 5+ bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
	if(nrow(dels.5bp) > 0){
		
		# lengths of the deleted segments
		dels.5bp.context.middle.lengths <- nchar(as.character(dels.5bp[,'REF'])) - 1
		
		# extract symmetric upstream/downstream flanks sized by context_bp and extend end by deletion length
		dels.5bp.context <- as.character(
			subseq(
				x = reference[as.character(dels.5bp[,'CHROM'])],
				start = as.numeric(dels.5bp[,'POS']) - (.flank - 1),
				end = as.numeric(dels.5bp[,'POS']) + dels.5bp.context.middle.lengths + .flank
				)
			)
		
		# build [deleted] middle
		dels.5bp.context.middle <- as.character(dels.5bp[,'REF'])
		dels.5bp.context.middle <- str_split_fixed(dels.5bp.context.middle, '', 2)[,2,drop=F]
		dels.5bp.context.middle <- paste0('[', dels.5bp.context.middle, ']')
		
		#upstream
		dels.5bp.context.start <- rep(NA, nrow(dels.5bp))
		for (i in 1:length(dels.5bp.context.start)){
			l <- dels.5bp.context.middle.lengths[i]
			nC <- 2*.flank + l # total length of context string for this deletion
			s <- str_split_fixed(dels.5bp.context[i], '', nC)
			i1 <- max(1, .flank - (5 * l) + 1) # start index (clamped)
			i2 <- .flank # end index (upstream flank end)
			dels.5bp.context.start[i] <- paste(s[, i1:i2, drop=FALSE], collapse='')
		}
		
		#downstream
		dels.5bp.context.end <- rep(NA, nrow(dels.5bp))
		for (i in 1:length(dels.5bp.context.end)){
			l <- dels.5bp.context.middle.lengths[i]
			nC <- 2*.flank + l
			s <- str_split_fixed(dels.5bp.context[i], '', nC)
			i1 <- .flank + 1 + l # first base after the deleted block
			i2 <- min(.flank + (6*l), nC) # clamp to available sequence
			dels.5bp.context.end[i] <- paste(s[, i1:i2, drop=FALSE], collapse='')
		}
		
		dels.5bp[,'CONTEXT'] <- paste(dels.5bp.context.start, dels.5bp.context.middle, dels.5bp.context.end, sep = '')
		
	}
	
	## summarise >=2 bp deletions at simple repeats in matrix format
	## also create sub-tables for repeat length = 1 hits, to later assess these for microhomologies
	dels.greater.2bp.summary <- matrix(0, ncol = 4, nrow = 6)
	colnames(dels.greater.2bp.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## deletion size
	rownames(dels.greater.2bp.summary) <- c('1 or MH', '2', '3', '4', '5', '6+') ## number of repeats
	
	## 2 bp
	if(nrow(dels.2bp) > 0){
		
		repeat.nts <- str_split_fixed(str_split_fixed(dels.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		downstream.context <- str_split_fixed(str_split_fixed(dels.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
		dels.2bp.pot.MH.ind <- c()
		for (i in 1:nrow(dels.2bp)){
			
			## look at repeat
			tmp.repeat.nts <- repeat.nts[i]
			tmp.repeat.length <- nchar(tmp.repeat.nts)
			
			## how often does it match consecutively in the immediate downstream context?
			tmp.downstream.context <- downstream.context[i]
			tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
			
			## group
			tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
																	paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
			
			### check all 6 categories for upstream bases
			if(tmp.downstream.context[1] %in% tmp.repeat.nts){
				
				if(tmp.downstream.context[2] %in% tmp.repeat.nts){
					
					if(tmp.downstream.context[3] %in% tmp.repeat.nts){
						
						if(tmp.downstream.context[4] %in% tmp.repeat.nts){
							
							if(tmp.downstream.context[5] %in% tmp.repeat.nts){
								
								dels.greater.2bp.summary['6+', '2 bp'] <- dels.greater.2bp.summary['6+', '2 bp'] + 1
								
							} else{
								
								dels.greater.2bp.summary['5', '2 bp'] <- dels.greater.2bp.summary['5', '2 bp'] + 1
								
							}
							
						}else{
							
							dels.greater.2bp.summary['4', '2 bp'] <- dels.greater.2bp.summary['4', '2 bp'] + 1
							
						}
						
					}else{
						
						dels.greater.2bp.summary['3', '2 bp'] <- dels.greater.2bp.summary['3', '2 bp'] + 1
						
					}
					
				}else{
					
					dels.greater.2bp.summary['2', '2 bp'] <- dels.greater.2bp.summary['2', '2 bp'] + 1
					
				}
				
			}else{
				
				dels.greater.2bp.summary['1 or MH', '2 bp'] <- dels.greater.2bp.summary['1 or MH', '2 bp'] + 1
				dels.2bp.pot.MH.ind <- c(dels.2bp.pot.MH.ind, i)
				
			}
			
		}
		dels.2bp.pot.MH <- dels.2bp[dels.2bp.pot.MH.ind,,drop=F] 
		
	}
	
	## 3 bp
	if(nrow(dels.3bp) > 0){
		
		repeat.nts <- str_split_fixed(str_split_fixed(dels.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		downstream.context <- str_split_fixed(str_split_fixed(dels.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
		dels.3bp.pot.MH.ind <- c()
		for (i in 1:nrow(dels.3bp)){
			
			## look at repeat
			tmp.repeat.nts <- repeat.nts[i]
			tmp.repeat.length <- nchar(tmp.repeat.nts)
			
			## how often does it match consecutively in the immediate downstream context?
			tmp.downstream.context <- downstream.context[i]
			tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
			
			## group
			tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
																	paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
			
			### check all 6 categories for upstream bases
			if(tmp.downstream.context[1] %in% tmp.repeat.nts){
				
				if(tmp.downstream.context[2] %in% tmp.repeat.nts){
					
					if(tmp.downstream.context[3] %in% tmp.repeat.nts){
						
						if(tmp.downstream.context[4] %in% tmp.repeat.nts){
							
							if(tmp.downstream.context[5] %in% tmp.repeat.nts){
								
								dels.greater.2bp.summary['6+', '3 bp'] <- dels.greater.2bp.summary['6+', '3 bp'] + 1
								
							} else{
								
								dels.greater.2bp.summary['5', '3 bp'] <- dels.greater.2bp.summary['5', '3 bp'] + 1
								
							}
							
						}else{
							
							dels.greater.2bp.summary['4', '3 bp'] <- dels.greater.2bp.summary['4', '3 bp'] + 1
							
						}
						
					}else{
						
						dels.greater.2bp.summary['3', '3 bp'] <- dels.greater.2bp.summary['3', '3 bp'] + 1
						
					}
					
				}else{
					
					dels.greater.2bp.summary['2', '3 bp'] <- dels.greater.2bp.summary['2', '3 bp'] + 1
					
				}
				
			}else{
				
				dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] + 1
				dels.3bp.pot.MH.ind <- c(dels.3bp.pot.MH.ind, i)
				
			}
			
		}
		dels.3bp.pot.MH <- dels.3bp[dels.3bp.pot.MH.ind,,drop=F] 
		
	}
	
	## 4 bp
	if(nrow(dels.4bp) > 0){
		
		repeat.nts <- str_split_fixed(str_split_fixed(dels.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		downstream.context <- str_split_fixed(str_split_fixed(dels.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
		dels.4bp.pot.MH.ind <- c()
		for (i in 1:nrow(dels.4bp)){
			
			## look at repeat
			tmp.repeat.nts <- repeat.nts[i]
			tmp.repeat.length <- nchar(tmp.repeat.nts)
			
			## how often does it match consecutively in the immediate downstream context?
			tmp.downstream.context <- downstream.context[i]
			tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
			
			## group
			tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
																	paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
			
			### check all 6 categories for upstream bases
			if(tmp.downstream.context[1] %in% tmp.repeat.nts){
				
				if(tmp.downstream.context[2] %in% tmp.repeat.nts){
					
					if(tmp.downstream.context[3] %in% tmp.repeat.nts){
						
						if(tmp.downstream.context[4] %in% tmp.repeat.nts){
							
							if(tmp.downstream.context[5] %in% tmp.repeat.nts){
								
								dels.greater.2bp.summary['6+', '4 bp'] <- dels.greater.2bp.summary['6+', '4 bp'] + 1
								
							} else{
								
								dels.greater.2bp.summary['5', '4 bp'] <- dels.greater.2bp.summary['5', '4 bp'] + 1
								
							}
							
						}else{
							
							dels.greater.2bp.summary['4', '4 bp'] <- dels.greater.2bp.summary['4', '4 bp'] + 1
							
						}
						
					}else{
						
						dels.greater.2bp.summary['3', '4 bp'] <- dels.greater.2bp.summary['3', '4 bp'] + 1
						
					}
					
				}else{
					
					dels.greater.2bp.summary['2', '4 bp'] <- dels.greater.2bp.summary['2', '4 bp'] + 1
					
				}
				
			}else{
				
				dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] + 1
				dels.4bp.pot.MH.ind <- c(dels.4bp.pot.MH.ind, i)
				
			}
			
		}
		dels.4bp.pot.MH <- dels.4bp[dels.4bp.pot.MH.ind,,drop=F]
		
	}
	
	## 5+ bp
	if(nrow(dels.5bp) > 0){
		
		repeat.nts <- str_split_fixed(str_split_fixed(dels.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		downstream.context <- str_split_fixed(str_split_fixed(dels.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
		dels.5bp.pot.MH.ind <- c()
		for (i in 1:nrow(dels.5bp)){
			
			## look at repeat
			tmp.repeat.nts <- repeat.nts[i]
			tmp.repeat.length <- nchar(tmp.repeat.nts)
			
			## how often does it match consecutively in the immediate downstream context?
			tmp.downstream.context <- downstream.context[i]
			tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
			
			## group
			tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
																	paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
			
			### check all 6 categories for upstream bases
			if(tmp.downstream.context[1] %in% tmp.repeat.nts){
				
				if(tmp.downstream.context[2] %in% tmp.repeat.nts){
					
					if(tmp.downstream.context[3] %in% tmp.repeat.nts){
						
						if(tmp.downstream.context[4] %in% tmp.repeat.nts){
							
							if(tmp.downstream.context[5] %in% tmp.repeat.nts){
								
								dels.greater.2bp.summary['6+', '5+ bp'] <- dels.greater.2bp.summary['6+', '5+ bp'] + 1
								
							} else{
								
								dels.greater.2bp.summary['5', '5+ bp'] <- dels.greater.2bp.summary['5', '5+ bp'] + 1
								
							}
							
						}else{
							
							dels.greater.2bp.summary['4', '5+ bp'] <- dels.greater.2bp.summary['4', '5+ bp'] + 1
							
						}
						
					}else{
						
						dels.greater.2bp.summary['3', '5+ bp'] <- dels.greater.2bp.summary['3', '5+ bp'] + 1
						
					}
					
				}else{
					
					dels.greater.2bp.summary['2', '5+ bp'] <- dels.greater.2bp.summary['2', '5+ bp'] + 1
					
				}
				
			}else{
				
				dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] + 1
				dels.5bp.pot.MH.ind <- c(dels.5bp.pot.MH.ind, i)
				
			}
			
		}
		dels.5bp.pot.MH <- dels.5bp[dels.5bp.pot.MH.ind,,drop=F] 
		
	}
	
	# (v) >=2 bp deletions at microhomologies
	## go directly into in matrix format, thereby also re-classify repeat length = 1 deletions
	dels.greater.2bp.MH.summary <- matrix(0, ncol = 4, nrow = 5)
	colnames(dels.greater.2bp.MH.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## deletion size
	rownames(dels.greater.2bp.MH.summary) <- c('1 bp MH', '2 bp MH', '3 bp MH', '4 bp MH', '5+ bp MH') ## microhomology length
	dels.greater.2bp.MH.summary[2:5, '2 bp'] <- NA
	dels.greater.2bp.MH.summary[3:5, '3 bp'] <- NA
	dels.greater.2bp.MH.summary[4:5, '4 bp'] <- NA
	
	## 2 bp
	if(nrow(dels.2bp) > 0){
		
		if(nrow(dels.2bp.pot.MH) > 0){
			
			repeat.nts <- str_split_fixed(str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
			upstream.context <- str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
			downstream.context <- str_split_fixed(str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
			for (i in 1:nrow(dels.2bp.pot.MH)){
				
				## look at repeat
				tmp.repeat.nts <- repeat.nts[i]
				
				## isolate first nucleotide in the immediate upstream context
				tmp.upstream.context <- upstream.context[i]
				tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][length(strsplit(tmp.upstream.context, '')[[1]])]
				
				## isolate first nucleotide in the immediate downstream context
				tmp.downstream.context <- downstream.context[i]
				tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1]
				
				if(strsplit(tmp.repeat.nts, '')[[1]][2] == tmp.upstream.context | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context){
					
					dels.greater.2bp.MH.summary['1 bp MH', '2 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '2 bp'] + 1
					dels.greater.2bp.summary['1 or MH', '2 bp'] <- dels.greater.2bp.summary['1 or MH', '2 bp'] - 1
					
				}
			} 
			
		}
		
	}
	
	## 3 bp
	if(nrow(dels.3bp) > 0){
		
		if(nrow(dels.3bp.pot.MH) > 0){
			
			repeat.nts <- str_split_fixed(str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
			upstream.context <- str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
			downstream.context <- str_split_fixed(str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
			for (i in 1:nrow(dels.3bp.pot.MH)){
				
				## look at repeat
				tmp.repeat.nts <- repeat.nts[i]
				
				## isolate first two nucleotides in the immediate upstream context
				tmp.upstream.context <- upstream.context[i]
				tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 1):length(strsplit(tmp.upstream.context, '')[[1]])]
				
				## isolate first two nucleotides in the immediate downstream context
				tmp.downstream.context <- downstream.context[i]
				tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:2]
				
				### MH length = 2 (important to start with the highest possible MH length)
				
				if(all(c(strsplit(tmp.repeat.nts, '')[[1]][2:3] == tmp.upstream.context) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context) == T)){
					
					dels.greater.2bp.MH.summary['2 bp MH', '3 bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '3 bp'] + 1
					dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] - 1
					
				}else{
					
					### MH length = 1
					
					if(strsplit(tmp.repeat.nts, '')[[1]][3] == tmp.upstream.context[2] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
						
						dels.greater.2bp.MH.summary['1 bp MH', '3 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '3 bp'] + 1
						dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] - 1
						
					}
					
				}
				
			} 
			
		}
		
	}
	
	## 4 bp
	if(nrow(dels.4bp) > 0){
		
		if(nrow(dels.4bp.pot.MH) > 0){
			
			repeat.nts <- str_split_fixed(str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
			upstream.context <- str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
			downstream.context <- str_split_fixed(str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
			for (i in 1:nrow(dels.4bp.pot.MH)){
				
				## look at repeat
				tmp.repeat.nts <- repeat.nts[i]
				
				## isolate first three nucleotides in the immediate upstream context
				tmp.upstream.context <- upstream.context[i]
				tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 2):length(strsplit(tmp.upstream.context, '')[[1]])]
				
				## isolate first three nucleotides in the immediate downstream context
				tmp.downstream.context <- downstream.context[i]
				tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:3]
				
				### MH length = 3 (important to start with the highest possible MH length)
				
				if(all(c(strsplit(tmp.repeat.nts, '')[[1]][2:4] == tmp.upstream.context) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:3] == tmp.downstream.context) == T)){
					
					dels.greater.2bp.MH.summary['3 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['3 bp MH', '4 bp'] + 1
					dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
					
				}else{
					
					### MH length = 2
					
					if(all(c(strsplit(tmp.repeat.nts, '')[[1]][3:4] == tmp.upstream.context[2:3]) ==T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context[1:2]) == T)){
						
						dels.greater.2bp.MH.summary['2 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '4 bp'] + 1
						dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
						
					}else{
						
						### MH length = 1
						
						if(strsplit(tmp.repeat.nts, '')[[1]][4] == tmp.upstream.context[3] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
							
							dels.greater.2bp.MH.summary['1 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '4 bp'] + 1
							dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
							
						}
						
					}
					
				}
				
			} 
			
		}
		
	}
	
	## 5+ bp
	if(nrow(dels.5bp) > 0){
		
		if(nrow(dels.5bp.pot.MH) > 0){
			
			repeat.nts <- str_split_fixed(str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
			upstream.context <- str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
			downstream.context <- str_split_fixed(str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
			for (i in 1:nrow(dels.5bp.pot.MH)){
				
				## look at repeat
				tmp.repeat.nts <- repeat.nts[i]
				tmp.repeat.length <- nchar(tmp.repeat.nts)
				
				## isolate first X (= repeat length - 1) nucleotides in the immediate upstream context
				tmp.upstream.context <- upstream.context[i]
				tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 4):length(strsplit(tmp.upstream.context, '')[[1]])]
				
				## isolate first X (= repeat length - 1) nucleotides in the immediate downstream context
				tmp.downstream.context <- downstream.context[i]
				tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:5]
				
				### MH length = 5+ (important to start with the highest possible MH length); only need to look at first 5 up/downstream NTs
				
				if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-4):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 4):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:5] == tmp.downstream.context[1:5]) == T)){
					
					dels.greater.2bp.MH.summary['5+ bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['5+ bp MH', '5+ bp'] + 1
					dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
					
				}else{
					
					### MH length = 4
					
					if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-3):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 3):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:4] == tmp.downstream.context[1:4]) == T)){
						
						dels.greater.2bp.MH.summary['4 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['4 bp MH', '5+ bp'] + 1
						dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
						
					}else{
						
						### MH length = 3
						
						if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-2):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 2):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:3] == tmp.downstream.context[1:3]) == T)){
							
							dels.greater.2bp.MH.summary['3 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['3 bp MH', '5+ bp'] + 1
							dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
							
						}else{
							
							### MH length = 2
							
							if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-1):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 1):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context[1:2]) == T)){
								
								dels.greater.2bp.MH.summary['2 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '5+ bp'] + 1
								dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
								
							}else{
								
								if(strsplit(tmp.repeat.nts, '')[[1]][tmp.repeat.length] == tmp.upstream.context[length(tmp.upstream.context)] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
									
									dels.greater.2bp.MH.summary['1 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '5+ bp'] + 1
									dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
									
								}
								
							}
							
						}
						
					}
					
				}
				
			} 
			
		}
		
	}
	
	## after MH cases have been taken out, rename row of deletions at 1-unit repeats in original table
	rownames(dels.greater.2bp.summary)[1] <- '1'
	
	## 4. classify >=2 bp insertions into:
	
	# (i) 2 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
	if(nrow(ins.2bp) > 0){
		
		ins.2bp.context <- as.character(subseq(x = reference[as.character(ins.2bp[,'CHROM'])], 
																					 start = as.numeric(ins.2bp[,'POS']) - 9, 
																					 end = as.numeric(ins.2bp[,'POS']) + 10))
		ins.2bp.context.middle <- as.character(ins.2bp[,'ALT'])
		ins.2bp.context.middle <- paste0('[', str_split_fixed(ins.2bp.context.middle, '', 2)[,2,drop=F], ']')
		ins.2bp.context.start <- str_split_fixed(ins.2bp.context, '', 11)[,1:10,drop=F]
		ins.2bp.context.start <- paste(ins.2bp.context.start[,1], ins.2bp.context.start[,2], ins.2bp.context.start[,3],
																	 ins.2bp.context.start[,4], ins.2bp.context.start[,5], ins.2bp.context.start[,6],
																	 ins.2bp.context.start[,7], ins.2bp.context.start[,8], ins.2bp.context.start[,9],
																	 ins.2bp.context.start[,10], sep = '')
		ins.2bp.context.end <- str_split_fixed(ins.2bp.context, '', 11)[,11,drop=F]
		ins.2bp[,'CONTEXT'] <- paste(ins.2bp.context.start, ins.2bp.context.middle, ins.2bp.context.end, sep = '')
		
	}
	
	# (ii) 3 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
	if(nrow(ins.3bp) > 0){
		
		ins.3bp.context <- as.character(subseq(x = reference[as.character(ins.3bp[,'CHROM'])], 
																					 start = as.numeric(ins.3bp[,'POS']) - 14, 
																					 end = as.numeric(ins.3bp[,'POS']) + 15))
		ins.3bp.context.middle <- as.character(ins.3bp[,'ALT'])
		ins.3bp.context.middle <- paste0('[', str_split_fixed(ins.3bp.context.middle, '', 2)[,2,drop=F], ']')
		ins.3bp.context.start <- str_split_fixed(ins.3bp.context, '', 16)[,1:15,drop=F]
		ins.3bp.context.start <- paste(ins.3bp.context.start[,1], ins.3bp.context.start[,2], ins.3bp.context.start[,3],
																	 ins.3bp.context.start[,4], ins.3bp.context.start[,5], ins.3bp.context.start[,6],
																	 ins.3bp.context.start[,7], ins.3bp.context.start[,8], ins.3bp.context.start[,9],
																	 ins.3bp.context.start[,10], ins.3bp.context.start[,11], ins.3bp.context.start[,12], 
																	 ins.3bp.context.start[,13], ins.3bp.context.start[,14], ins.3bp.context.start[,15], sep = '')
		ins.3bp.context.end <- str_split_fixed(ins.3bp.context, '', 16)[,16,drop=F]
		ins.3bp[,'CONTEXT'] <- paste(ins.3bp.context.start, ins.3bp.context.middle, ins.3bp.context.end, sep = '')
		
	}
	
	# (iii) 4 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
	if(nrow(ins.4bp) > 0){
		
		ins.4bp.context <- as.character(subseq(x = reference[as.character(ins.4bp[,'CHROM'])], 
																					 start = as.numeric(ins.4bp[,'POS']) - 19, 
																					 end = as.numeric(ins.4bp[,'POS']) + 20))
		ins.4bp.context.middle <- as.character(ins.4bp[,'ALT'])
		ins.4bp.context.middle <- paste0('[', str_split_fixed(ins.4bp.context.middle, '', 2)[,2,drop=F], ']')
		ins.4bp.context.start <- str_split_fixed(ins.4bp.context, '', 21)[,1:20,drop=F]
		ins.4bp.context.start <- paste(ins.4bp.context.start[,1], ins.4bp.context.start[,2], ins.4bp.context.start[,3],
																	 ins.4bp.context.start[,4], ins.4bp.context.start[,5], ins.4bp.context.start[,6],
																	 ins.4bp.context.start[,7], ins.4bp.context.start[,8], ins.4bp.context.start[,9],
																	 ins.4bp.context.start[,10], ins.4bp.context.start[,11], ins.4bp.context.start[,12], 
																	 ins.4bp.context.start[,13], ins.4bp.context.start[,14], ins.4bp.context.start[,15], 
																	 ins.4bp.context.start[,16], ins.4bp.context.start[,17], ins.4bp.context.start[,18], 
																	 ins.4bp.context.start[,19], ins.4bp.context.start[,20], sep = '')
		ins.4bp.context.end <- str_split_fixed(ins.4bp.context, '', 21)[,21,drop=F]
		ins.4bp[,'CONTEXT'] <- paste(ins.4bp.context.start, ins.4bp.context.middle, ins.4bp.context.end, sep = '')
		
	}
	
	# (iv) 5+ bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
	if(nrow(ins.5bp) > 0){
		
		ins.5bp.context <- as.character(
			subseq(
				x = reference[as.character(ins.5bp[,'CHROM'])],
				start = as.numeric(ins.5bp[,'POS']) - (.flank - 1L),
				end = as.numeric(ins.5bp[,'POS']) + .flank
			)
		)
		
		# inserted segment and its length
		ins.5bp.context.middle <- as.character(ins.5bp[,'ALT'])
		ins.5bp.context.middle <- str_split_fixed(ins.5bp.context.middle, '', 2)[,2,drop=F]
		ins.5bp.context.middle <- paste0('[', ins.5bp.context.middle, ']')
		ins.5bp.context.middle.lengths <- nchar(as.character(ins.5bp[,'ALT'])) - 1

		# upstream: take up to 5×repeat-length
		ins.5bp.context.start <- rep(NA, nrow(ins.5bp))
		for (i in 1:length(ins.5bp.context.start)){
			l <- ins.5bp.context.middle.lengths[i]
			s <- str_split_fixed(ins.5bp.context[i], '', 2*.flank)
			i1 <- max(1, .flank - (5 * l) + 1)
			i2 <- .flank
			ins.5bp.context.start[i] <- paste(s[, i1:i2, drop=FALSE], collapse='')
		}
		
		# downstream: take up to 6×repeat-length starting immediately after POS
		ins.5bp.context.end <- rep(NA, nrow(ins.5bp))
		for (i in 1:length(ins.5bp.context.end)){
			l <- ins.5bp.context.middle.lengths[i]
			s <- str_split_fixed(ins.5bp.context[i], '', 2*.flank)
			i1 <- .flank + 1
			i2 <- min(.flank + (6 * l), 2*.flank) # clamp to available sequence
			ins.5bp.context.end[i] <- paste(s[, i1:i2, drop=FALSE], collapse='')
		}
		
		ins.5bp[,'CONTEXT'] <- paste(ins.5bp.context.start, ins.5bp.context.middle, ins.5bp.context.end, sep = '')
		
	}
	
	## summarise >=2 bp insertions at simple repeats in matrix format
	ins.greater.2bp.summary <- matrix(0, ncol = 4, nrow = 6)
	colnames(ins.greater.2bp.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## insertion size
	rownames(ins.greater.2bp.summary) <- c('0', '1', '2', '3', '4', '5+') ## number of repeats
	
	## 2 bp
	if(nrow(ins.2bp) > 0){
		
		repeat.nts <- str_split_fixed(str_split_fixed(ins.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		downstream.context <- str_split_fixed(str_split_fixed(ins.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
		for (i in 1:nrow(ins.2bp)){
			
			## look at repeat
			tmp.repeat.nts <- repeat.nts[i]
			tmp.repeat.length <- nchar(tmp.repeat.nts)
			
			## how often does it match consecutively in the immediate downstream context?
			tmp.downstream.context <- downstream.context[i]
			tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
			
			## group
			tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
																	paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
			
			### check all 6 categories for upstream bases
			if(tmp.downstream.context[1] %in% tmp.repeat.nts){
				
				if(tmp.downstream.context[2] %in% tmp.repeat.nts){
					
					if(tmp.downstream.context[3] %in% tmp.repeat.nts){
						
						if(tmp.downstream.context[4] %in% tmp.repeat.nts){
							
							if(tmp.downstream.context[5] %in% tmp.repeat.nts){
								
								ins.greater.2bp.summary['5+', '2 bp'] <- ins.greater.2bp.summary['5+', '2 bp'] + 1
								
							} else{
								
								ins.greater.2bp.summary['4', '2 bp'] <- ins.greater.2bp.summary['4', '2 bp'] + 1
								
							}
							
						}else{
							
							ins.greater.2bp.summary['3', '2 bp'] <- ins.greater.2bp.summary['3', '2 bp'] + 1
							
						}
						
					}else{
						
						ins.greater.2bp.summary['2', '2 bp'] <- ins.greater.2bp.summary['2', '2 bp'] + 1
						
					}
					
				}else{
					
					ins.greater.2bp.summary['1', '2 bp'] <- ins.greater.2bp.summary['1', '2 bp'] + 1
					
				}
				
			}else{
				
				ins.greater.2bp.summary['0', '2 bp'] <- ins.greater.2bp.summary['0', '2 bp'] + 1
				
			}
			
		}
		
	}
	
	## 3 bp
	if(nrow(ins.3bp) > 0){
		
		repeat.nts <- str_split_fixed(str_split_fixed(ins.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		downstream.context <- str_split_fixed(str_split_fixed(ins.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
		for (i in 1:nrow(ins.3bp)){
			
			## look at repeat
			tmp.repeat.nts <- repeat.nts[i]
			tmp.repeat.length <- nchar(tmp.repeat.nts)
			
			## how often does it match consecutively in the immediate downstream context?
			tmp.downstream.context <- downstream.context[i]
			tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
			
			## group
			tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
																	paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
			
			### check all 6 categories for upstream bases
			if(tmp.downstream.context[1] %in% tmp.repeat.nts){
				
				if(tmp.downstream.context[2] %in% tmp.repeat.nts){
					
					if(tmp.downstream.context[3] %in% tmp.repeat.nts){
						
						if(tmp.downstream.context[4] %in% tmp.repeat.nts){
							
							if(tmp.downstream.context[5] %in% tmp.repeat.nts){
								
								ins.greater.2bp.summary['5+', '3 bp'] <- ins.greater.2bp.summary['5+', '3 bp'] + 1
								
							} else{
								
								ins.greater.2bp.summary['4', '3 bp'] <- ins.greater.2bp.summary['4', '3 bp'] + 1
								
							}
							
						}else{
							
							ins.greater.2bp.summary['3', '3 bp'] <- ins.greater.2bp.summary['3', '3 bp'] + 1
							
						}
						
					}else{
						
						ins.greater.2bp.summary['2', '3 bp'] <- ins.greater.2bp.summary['2', '3 bp'] + 1
						
					}
					
				}else{
					
					ins.greater.2bp.summary['1', '3 bp'] <- ins.greater.2bp.summary['1', '3 bp'] + 1
					
				}
				
			}else{
				
				ins.greater.2bp.summary['0', '3 bp'] <- ins.greater.2bp.summary['0', '3 bp'] + 1
				
			}
			
		} 
		
	}
	
	## 4 bp
	if(nrow(ins.4bp) > 0){
		
		repeat.nts <- str_split_fixed(str_split_fixed(ins.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		downstream.context <- str_split_fixed(str_split_fixed(ins.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
		for (i in 1:nrow(ins.4bp)){
			
			## look at repeat
			tmp.repeat.nts <- repeat.nts[i]
			tmp.repeat.length <- nchar(tmp.repeat.nts)
			
			## how often does it match consecutively in the immediate downstream context?
			tmp.downstream.context <- downstream.context[i]
			tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
			
			## group
			tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
																	paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
			
			### check all 6 categories for upstream bases
			if(tmp.downstream.context[1] %in% tmp.repeat.nts){
				
				if(tmp.downstream.context[2] %in% tmp.repeat.nts){
					
					if(tmp.downstream.context[3] %in% tmp.repeat.nts){
						
						if(tmp.downstream.context[4] %in% tmp.repeat.nts){
							
							if(tmp.downstream.context[5] %in% tmp.repeat.nts){
								
								ins.greater.2bp.summary['5+', '4 bp'] <- ins.greater.2bp.summary['5+', '4 bp'] + 1
								
							} else{
								
								ins.greater.2bp.summary['4', '4 bp'] <- ins.greater.2bp.summary['4', '4 bp'] + 1
								
							}
							
						}else{
							
							ins.greater.2bp.summary['3', '4 bp'] <- ins.greater.2bp.summary['3', '4 bp'] + 1
							
						}
						
					}else{
						
						ins.greater.2bp.summary['2', '4 bp'] <- ins.greater.2bp.summary['2', '4 bp'] + 1
						
					}
					
				}else{
					
					ins.greater.2bp.summary['1', '4 bp'] <- ins.greater.2bp.summary['1', '4 bp'] + 1
					
				}
				
			}else{
				
				ins.greater.2bp.summary['0', '4 bp'] <- ins.greater.2bp.summary['0', '4 bp'] + 1
				
			}
			
		} 
		
	}
	
	## 5+ bp
	if(nrow(ins.5bp) > 0){
		
		repeat.nts <- str_split_fixed(str_split_fixed(ins.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
		downstream.context <- str_split_fixed(str_split_fixed(ins.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
		for (i in 1:nrow(ins.5bp)){
			
			## look at repeat
			tmp.repeat.nts <- repeat.nts[i]
			tmp.repeat.length <- nchar(tmp.repeat.nts)
			
			## how often does it match consecutively in the immediate downstream context?
			tmp.downstream.context <- downstream.context[i]
			tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
			
			## group
			tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
																	paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
																	paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
			
			### check all 6 categories for upstream bases
			if(tmp.downstream.context[1] %in% tmp.repeat.nts){
				
				if(tmp.downstream.context[2] %in% tmp.repeat.nts){
					
					if(tmp.downstream.context[3] %in% tmp.repeat.nts){
						
						if(tmp.downstream.context[4] %in% tmp.repeat.nts){
							
							if(tmp.downstream.context[5] %in% tmp.repeat.nts){
								
								ins.greater.2bp.summary['5+', '5+ bp'] <- ins.greater.2bp.summary['5+', '5+ bp'] + 1
								
							} else{
								
								ins.greater.2bp.summary['4', '5+ bp'] <- ins.greater.2bp.summary['4', '5+ bp'] + 1
								
							}
							
						}else{
							
							ins.greater.2bp.summary['3', '5+ bp'] <- ins.greater.2bp.summary['3', '5+ bp'] + 1
							
						}
						
					}else{
						
						ins.greater.2bp.summary['2', '5+ bp'] <- ins.greater.2bp.summary['2', '5+ bp'] + 1
						
					}
					
				}else{
					
					ins.greater.2bp.summary['1', '5+ bp'] <- ins.greater.2bp.summary['1', '5+ bp'] + 1
					
				}
				
			}else{
				
				ins.greater.2bp.summary['0', '5+ bp'] <- ins.greater.2bp.summary['0', '5+ bp'] + 1
				
			}
			
		} 
		
	}
	
	## 5. summarise outputs as five lists, each featuring one table per major ID category
	out <- list('1bp del' = dels.1bp.summary,
							'1bp ins' = ins.1bp.summary,
							'>2bp del' = dels.greater.2bp.summary,
							'>2bp ins' = ins.greater.2bp.summary,
							'MH dels' = dels.greater.2bp.MH.summary)
	return(out)
	
}

#Function to convert indel spectrum produced by indel.spectrum() to sigfit format
indelspectrum.to.sigfit <- function(indelspectrum, spectrum_type = c("pyr", "template")){

	spectrum_type <- match.arg(spectrum_type)
	indel_levels <- if(spectrum_type == "pyr") indel_pyr_labels.sigfit else indel_template_strand_labels.sigfit
	
	#Helpers to convert from indelwald to sigfit labels
	map_1bp_rep  <- function(x){
		case_match(x, "0 bp" ~ "0", "1 bp" ~ "1", "2 bp" ~ "2", "3 bp" ~ "3", "4 bp" ~ "4", "5+ bp" ~ "5")
	}
	
	map_sz <- function(x){
		case_match(x, "2 bp" ~ "2", "3 bp" ~ "3", "4 bp" ~ "4", "5+ bp" ~ "5")
	}
	
	map_R_del <- function(x){
		case_match(x, "1" ~ "0", "2" ~ "1", "3" ~ "2", "4" ~ "3", "5" ~ "4", "6+" ~ "5")
	}
	
	map_R_ins <- function(x){
		case_match(x,"0" ~ "0", "1" ~ "1", "2" ~ "2", "3" ~ "3", "4" ~ "4", "5+" ~ "5")
	}
	
	map_MH_len <- function(x){
		case_match(x,"1 bp MH" ~ "1", "2 bp MH" ~ "2", "3 bp MH" ~ "3", "4 bp MH" ~ "4", "5+ bp MH" ~ "5")
	}
	
	indelspectrum %>%
		imap_dfr(function(mat, grp){
			mat %>%
				as.data.frame %>%
				rownames_to_column("rowname") %>%
				pivot_longer(-rowname, names_to="colname", values_to="count") %>%
				mutate(group = grp)
		}) %>%
		filter(!is.na(count)) %>%
		mutate(
			label = case_when(
				group == "1bp del" ~ str_c("1","Del", colname, map_1bp_rep(rowname), sep=":"),
				group == "1bp ins" ~ str_c("1","Ins", colname, map_1bp_rep(rowname), sep=":"),
				group == ">2bp del" ~ str_c(map_sz(colname), "Del", "R", map_R_del(rowname), sep=":"),
				group == ">2bp ins" ~ str_c(map_sz(colname), "Ins", "R", map_R_ins(rowname), sep=":"),
				group == "MH dels" ~ str_c(map_sz(colname), "Del", "M", map_MH_len(rowname), sep=":")
			) %>%
				factor(levels = indel_levels)
		)%>%
		select(label, count) %>%
		arrange(label) %>%
		pivot_wider(names_from = label, values_from = count)
}

# sigfit spectrum plotting helpers.
# Source: sigfit R package, https://github.com/kgori/sigfit, upstream default
# branch "master" at commit 46e89a1e861048e43ea1ce2367dcddfeff9d695e.
# Copied/adapted from R/sigfit_plotting.R and R/sigfit_utility.R.
# Authors: Kevin Gori and Adrian Baez-Ortega. License: GPL (>=3).
# Citation: Gori K, Baez-Ortega A. sigfit: flexible Bayesian inference of
# mutational signatures. bioRxiv 372896 (2018). doi: 10.1101/372896.
# HiDEF-seq uses sigfit-style plots for SBS (96), strand-wise SBS (192), and
# indel (83/107) spectra, so the DBS helper that depends on sigfit package data is
# intentionally not copied.

#' Initial coercion to matrix for signatures/exposures/counts
to_matrix <- function(x, int = FALSE) {
    # If x is coming from retrieve_pars, get mean
    if (is.list(x) & "mean" %in% names(x))
        x <- x$mean
    # If x is a vector, transform to 1-row matrix
    if (is.vector(x))
        x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
    # Otherwise, try coercing to matrix
    if (!is.matrix(x))
        x <- as.matrix(x)
    # For counts matrix: if real-valued, round
    if (int) {
        x <- round(x)
    }
    x
}


#' Generates character vector of 96 (or 192, if strand=TRUE) trinucleotide mutation types
mut_types <- function(strand = FALSE) {
    bases <- c("A", "C", "G", "T")
    muts <- paste0(rep(rep(bases, each = 4), 6),
                   rep(bases[c(2, 4)], each = 48),
                   rep(bases, 6 * 16 / 4),
                   ">",
                   rep(rep(bases, each = 4), 6),
                   c(rep(bases[-2], each = 16), rep(bases[-4], each = 16)),
                   rep(bases, 6 * 16 / 4))
    if (strand) {
        paste(c(rep("T", 96), rep("U", 96)), muts, sep = ":")
    }
    else {
        muts
    }
}


#' Plot mutational spectra
plot_spectrum <- function(spectra, pdf_path = NULL, pdf_width = 24, pdf_height = 8, name = NULL,
                          max_y = NULL, colors = NULL, boxes = TRUE, generic = FALSE) {
    # Fetch HPD interval values, if present
    if (is.list(spectra) & "mean" %in% names(spectra)) {
        spec <- to_matrix(spectra$mean)
        lwr <- to_matrix(spectra$lower)
        upr <- to_matrix(spectra$upper)
    }
    else {
        spec <- to_matrix(spectra)
        lwr <- NULL
        upr <- NULL
    }
    if (is.null(rownames(spec))) {
        rownames(spec) <- paste("Sample", 1:nrow(spec))
    }

    NTYPES <- c("SBS"=96, "TSW"=192, "ID"=83, "IDT"=107)  # number of categories per spectrum type
    NCAT <- ncol(spec)                          # number of categories in spectrum

    # Initialise PDF
    if (!is.null(pdf_path)) {
        cairo_pdf(pdf_path, width = pdf_width, height = pdf_height, onefile = TRUE)
        if (NCAT %in% NTYPES) {
            par(mar = c(5.5, 7, 7.5, 2))
        }
        else {
            par(mar = c(10, 7, 7.5, 2))
        }
    }

    # Generic spectrum (NCAT!={83,96,192})
    if (!(NCAT %in% NTYPES) | generic) {
        plot_spectrum_generic(spec, lwr, upr, name, max_y, colors, boxes)
    }
    else {
        # Standard SBS spectrum (NCAT=96)
        if (NCAT == NTYPES["SBS"]) {
            plot_spectrum_sbs(spec, lwr, upr, name, max_y, colors, boxes)
        }
        # Strand-wise SBS spectrum (NCAT=192)
        if (NCAT == NTYPES["TSW"]) {
            plot_spectrum_tsw(spec, lwr, upr, name, max_y, colors, boxes)
        }
        # Indel spectrum (NCAT=83)
        if (NCAT == NTYPES["ID"]) {
            plot_spectrum_id(spec, lwr, upr, name, max_y, colors, boxes)
        }
        # Template-strand indel spectrum (NCAT=107)
        if (NCAT == NTYPES["IDT"]) {
            plot_spectrum_id_template(spec, lwr, upr, name, max_y, colors, boxes)
        }
    }
    if (!is.null(pdf_path)) {
        invisible(dev.off())
    }
}


#' Plot SBS spectrum
plot_spectrum_sbs <- function(spec, lwr, upr, name, max_y, colors, boxes) {
    TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    LINECOL <- "gray60"
    XL <- c(0.2, 19.4, 38.6, 57.8, 77, 96.2)
    XR <- c(19.2, 38.4, 57.6, 76.8, 96, 115.2)
    FACTOR <- 1.095
    NCAT <- ncol(spec)   # number of categories
    NSAMP <- nrow(spec)  # number of samples
    for (i in 1:NSAMP) {
        if (is.null(max_y)) {
            samp_max_y <- max(0.05,
                              ifelse(is.null(upr), max(spec[i,]) * FACTOR, max(upr[i,]) * FACTOR))
        }
        else {
            samp_max_y <- max_y
        }
        if (boxes) {
            xlim <- c(-0.105, 115.5)
        }
        else {
            xlim <- c(-1, 116)
        }
        # Plot spectrum bars
        if (is.null(colors)) {
            colors = rep(COLORS, each = 16)
        }
        else if ((length(colors) > 1) & (length(colors) != NCAT)) {
            stop("'colors' must contain either a single value, or one value per mutation type.")
        }
        bars <- barplot(spec[i, ],
                        names.arg = substr(mut_types(), 1, 3), mgp = c(3, 0.8, 0),
                        col = colors, border = "white",
                        las = 2, ylim = c(0, samp_max_y), xlim = xlim,
                        yaxt = "n", cex.names = 1.6, xaxs = "i", family = "mono")
        # Highlight trinucleotide middle bases
        for (j in 1:length(COLORS)) {
            idx <- ((j-1) * 16 + 1):(j * 16)
            axis(side = 1, at = bars[idx], tick = FALSE, cex.axis = 1.6,
                 mgp = c(3, 0.8, 0), las = 2, family = "mono", font = 2,
                 col.axis = COLORS[j], labels = paste0(" ", substr(mut_types()[idx], 2, 2), " "))
        }
        # Plot axis
        if (any(spec > 1)) {
            axis(side = 2, cex.axis = 1.9, lwd = 2)
            label <- "Mutations"
            n_text <- paste0(" (", prettyNum(sum(spec[i,]), big.mark = ","), " mutations)")
        }
        else {
            axis(side = 2, at = seq(0, samp_max_y, ifelse(samp_max_y > 0.25, 0.1, 0.05)),
                 cex.axis = 1.9, lwd = 2)
            label <- "Mutation probability"
            n_text <- ""
        }
        if (is.null(name)) {
            nme <- rownames(spec)[i]
        }
        else {
            nme <- name
        }
        mtext(label, side = 2, cex = 2.4, line = 3.5)
        title(paste0(nme, n_text), line = 4, cex.main = 2.5)
        # Plot HPD intervals
        if (!is.null(lwr)) {
            arrows(bars, upr[i,],
                   bars, lwr[i,],
                   length = 0, lwd = 3, col = LINECOL)
        }
        # Plot mutation type labels
        text(x = (XL + XR) / 2, y = 1.055 * samp_max_y,
             labels = TYPES, cex = 2.4, xpd = TRUE)
        rect(xleft = XL, xright = XR, ybottom = 0.945 * samp_max_y,
             ytop = samp_max_y, col = COLORS, border = "white")
        # Plot box
        if (boxes) {
            box(lwd = 2)
        }
    }
}


#' Plot strandwise SBS spectrum
plot_spectrum_tsw <- function(spec, lwr, upr, name, max_y, colors, boxes) {
    TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    STRANDCOL <- c("deepskyblue3", "red3")
    BACKCOL <- c("#00BFFF33", "#00000033", "#EE2C2C33", "#C2C2C24D", "#A2CD5A4D", "#EEB4B44D")
    LINECOL <- "gray60"
    XL <- c(0.2, 19.4, 38.6, 57.8, 77, 96.2)
    XR <- c(19.2, 38.4, 57.6, 76.8, 96, 115.2)
    BACKLIM <- c(0, 46.8, 93.2, 139.55, 186, 232.35, 279.2)
    FACTOR <- 1.095
    NCAT <- ncol(spec)   # number of categories
    NSAMP <- nrow(spec)  # number of samples
    for (i in 1:NSAMP) {
        if (is.null(max_y)) {
            samp_max_y <- max(0.05,
                              ifelse(is.null(upr), max(spec[i,]) * FACTOR, max(upr[i,]) * FACTOR))
        }
        else {
            samp_max_y <- max_y
        }
        if (boxes) {
            xlim <- c(0, 279.2)
        }
        else {
            xlim <- c(-3, 280)
        }
        # Plot background panes and mutation type labels
        barplot(rbind(spec[i, 1:(NCAT/2)], spec[i, (NCAT/2+1):NCAT]), beside = TRUE,
                col = NA, border = NA, space = c(0.1, 0.8), xaxs = "i", yaxt = "n",
                xaxt = "n", ylim = c(0, samp_max_y), xlim = xlim)
        for (j in 1:length(COLORS)) {
            rect(xleft = BACKLIM[j], xright = BACKLIM[j+1], ybottom = 0,
                 ytop = samp_max_y, col = BACKCOL[j], border = "white")
            text(x = (BACKLIM[j] + BACKLIM[j+1]) / 2, y = 1.055 * samp_max_y,
                 labels = TYPES[j], cex = 2.4, xpd = TRUE)
            rect(xleft = BACKLIM[j], xright = BACKLIM[j+1], ybottom = 0.945 * samp_max_y,
                 ytop = samp_max_y, col = COLORS[j], border = "white")
        }
        # Plot legend
        legend("topright", bty = "n", inset = c(0.008, 0.03),
               legend = c("Central pyrimidine", "Central purine"),
               cex = 2.1, fill = NA, border = NA)
        legend("topright", bty = "n", inset = c(0.135, 0.03), pch=15, pt.cex=3.75,
               col=STRANDCOL, legend = c("", ""), cex = 2.1)
        # Plot spectrum bars
        bars <- barplot(rbind(spec[i, 1:(NCAT/2)],
                              spec[i, (NCAT/2+1):NCAT]),
                        names.arg = substr(mut_types(), 1, 3), beside = TRUE,
                        space = c(0.1, 0.8), mgp = c(3, 0.8, 0), las = 2,
                        col = STRANDCOL, border = "white", yaxt = "n",
                        cex.names = 1.6, xaxs = "i", family = "mono", add = TRUE)
        # Highlight trinucleotide middle bases
        for (j in 1:length(COLORS)) {
            idx <- ((j-1) * 16 + 1):(j * 16)
            axis(side = 1, at = colMeans(bars)[idx], tick = FALSE, cex.axis = 1.6,
                 mgp = c(3, 0.8, 0), las = 2, family = "mono", font = 2,
                 col.axis = COLORS[j], labels = paste0(" ", substr(mut_types()[idx], 2, 2), " "))
        }
        # Plot axis
        if (any(spec > 1)) {
            axis(side = 2, cex.axis = 1.9, lwd = 2)
            label <- "Mutations"
            n_text <- paste0(" (", prettyNum(sum(spec[i,]), big.mark = ","), " mutations)")
        }
        else {
            axis(side = 2, at = seq(0, samp_max_y, ifelse(samp_max_y > 0.25, 0.1, 0.05)),
                 cex.axis = 1.9, lwd = 2)
            label <- "Mutation probability"
            n_text <- ""
        }
        if (is.null(name)) {
            nme <- rownames(spec)[i]
        }
        else {
            nme <- name
        }
        mtext(label, side = 2, cex = 2.4, line = 3.5)
        title(paste0(nme, n_text), line = 4, cex.main = 2.5)
        # Plot HPD intervals
        if (!is.null(lwr)) {
            bars <- as.numeric(t(bars))
            arrows(bars, upr[i,],
                   bars, lwr[i,],
                   length = 0, lwd = 2.5, col = LINECOL)
        }
        # Plot box
        if (boxes) {
            box(lwd = 2)
        }
    }
}


#' Plot ID (indel) spectrum
plot_spectrum_id <- function(spec, lwr, upr, name, max_y, colors, boxes) {
    TYPES <- c("1-bp deletion", "1-bp insertion", ">1-bp deletion at repeat\n(deletion length)",
               ">1-bp insertion at repeat\n(insertion length)", "Microhomology\n(deletion length)")
    COLORS <- c(rep(c("#FCBD6F", "#FE8002", "#AFDC8A", "#36A02E", "#FCC9B4", "#FB896A", "#F04432",
                      "#BB191A", "#CFE0F1", "#93C3DE", "#4A97C8", "#1764AA"), each = 6), "#E1E1EE",
                rep("#B5B5D7", 2), rep("#8582BC", 3), rep("#62409A", 5))
    LINECOL <- "gray60"
    XL <- c(0.2, 7.4, 14.6, 21.8, 29, 36.2, 43.4, 50.6, 57.8, 65, 72.2, 79.4, 86.6, 87.8, 90.2, 93.8)
    XR <- c(7.2, 14.4, 21.6, 28.8, 36, 43.2, 50.4, 57.6, 64.8, 72, 79.2, 86.4, 87.6, 90, 93.6, 99.6)
    XM <- c(7.3, 21.7, 43.3, 72.1, 93.1)
    FACTOR <- 1.095
    NCAT <- ncol(spec)   # number of categories
    NSAMP <- nrow(spec)  # number of samples
    colnames(spec) <- c(rep(c(paste0(1:5, "  "), "6+"), 2), rep(c(paste0(0:4, "  "), "5+"), 2),
                        rep(c(paste0(1:5, "  "), "6+"), 4), rep(c(paste0(0:4, "  "), "5+"), 4),
                        paste0(c(1, 1:2, 1:3, 1:4), "  "), "5+")
    par(mar = c(5.5, 7, 11, 2))
    for (i in 1:NSAMP) {
        if (is.null(max_y)) {
            samp_max_y <- max(0.05,
                              ifelse(is.null(upr), max(spec[i,]) * FACTOR, max(upr[i,]) * FACTOR))
        }
        else {
            samp_max_y <- max_y
        }
        if (boxes) {
            xlim <- c(-0.105, 99.9)
        }
        else {
            xlim <- c(-1, 100.2)
        }
        # Plot spectrum bars
        if (is.null(colors)) {
            colors <- COLORS
        }
        else if ((length(colors) > 1) & (length(colors) != NCAT)) {
            stop("'colors' must contain either a single value, or one value per mutation type.")
        }
        bars <- barplot(spec[i, ], names.arg = colnames(spec),
                        mgp = c(3, 0.3, 0), col = colors, border = "white",
                        ylim = c(0, samp_max_y), xlim = xlim, yaxt = "n", xaxs = "i",
                        las = 2, cex.names = 1.5, adj = 0.5)
        # Plot axis
        if (any(spec > 1)) {
            axis(side = 2, cex.axis = 1.9, lwd = 2)
            label <- "Mutations"
            n_text <- paste0(" (", prettyNum(sum(spec[i,]), big.mark = ","), " mutations)")
        }
        else {
            axis(side = 2, at = seq(0, samp_max_y, ifelse(samp_max_y > 0.25, 0.1, 0.05)),
                 cex.axis = 1.9, lwd = 2)
            label <- "Mutation probability"
            n_text <- ""
        }
        if (is.null(name)) {
            nme <- rownames(spec)[i]
        }
        else {
            nme <- name
        }
        mtext(label, side = 2, cex = 2.4, line = 3.5)
        title(paste0(nme, n_text), line = 7.8, cex.main = 2.5)
        # Plot HPD intervals
        if (!is.null(lwr)) {
            arrows(bars, upr[i,],
                   bars, lwr[i,],
                   length = 0, lwd = 3, col = LINECOL)
        }
        # Plot mutation type labels
        rect(xleft = XL, xright = XR, ybottom = 0.945 * samp_max_y,
             ytop = samp_max_y, col = unique(COLORS), border = "white")
        mtext(TYPES, side = 3, at = XM, line = 2.4, cex = 2.2, xpd = TRUE)
        mtext(c(rep(c("C", "T"), 2), rep(c("2", "3", "4", "5+"), 3)),
              at = (XL + XR) / 2, side = 3, line = 0.35, cex = 2.1)
        mtext(c(rep(c("Homopolymer length", "Number of repeat units"), each = 2),
                "Microhomology length"),
              side = 1, at = XM + c(0, 0, 0, 0, 0.2), line = 2.8, cex = c(rep(1.95, 4), 1.85))
        # Plot box
        if (boxes) {
            box(lwd = 2)
        }
    }
}


#' Plot template-strand ID (indel) spectrum
plot_spectrum_id_template <- function(spec, lwr, upr, name, max_y, colors, boxes) {
    TYPES <- c("1-bp deletion", "1-bp insertion",
               ">1-bp deletion at repeat\n(deletion length)",
               ">1-bp insertion at repeat\n(insertion length)", "Microhomology\n(deletion length)")
    GROUP_WIDTHS <- c(rep(6, 16), 1, 2, 3, 5)
    TYPE_WIDTHS <- c(24, 24, 24, 24, 11)
    GROUP_LABELS <- c("C", "T", "G", "A", "C", "T", "G", "A",
                      rep(c("2", "3", "4", "5+"), 3))
    GROUP_COLORS <- c("#FCBD6F", "#FE8002", "#F9D6A8", "#C45D00",
                      "#AFDC8A", "#36A02E", "#D3EFB8", "#177D22",
                      "#FCC9B4", "#FB896A", "#F04432", "#BB191A",
                      "#CFE0F1", "#93C3DE", "#4A97C8", "#1764AA",
                      "#E1E1EE", "#B5B5D7", "#8582BC", "#62409A")
    COLORS <- c(rep(GROUP_COLORS[1:16], each = 6), GROUP_COLORS[17],
                rep(GROUP_COLORS[18], 2), rep(GROUP_COLORS[19], 3),
                rep(GROUP_COLORS[20], 5))
    LINECOL <- "gray60"
    FACTOR <- 1.095
    NCAT <- ncol(spec)   # number of categories
    NSAMP <- nrow(spec)  # number of samples
    colnames(spec) <- c(rep(c(paste0(1:5, "  "), "6+"), 4),
                        rep(c(paste0(0:4, "  "), "5+"), 4),
                        rep(c(paste0(1:5, "  "), "6+"), 4),
                        rep(c(paste0(0:4, "  "), "5+"), 4),
                        paste0(c(1, 1:2, 1:3, 1:4), "  "), "5+")
    par(mar = c(5.5, 7, 11, 2))
    for (i in 1:NSAMP) {
        if (is.null(max_y)) {
            samp_max_y <- max(0.05,
                              ifelse(is.null(upr), max(spec[i,]) * FACTOR, max(upr[i,]) * FACTOR))
        }
        else {
            samp_max_y <- max_y
        }
        bars_template <- barplot(spec[i, ], plot = FALSE)
        if (boxes) {
            xlim <- c(min(bars_template) - 0.805, max(bars_template) + 0.8)
        }
        else {
            xlim <- c(min(bars_template) - 1.7, max(bars_template) + 1.1)
        }
        # Plot spectrum bars
        if (is.null(colors)) {
            bar_colors <- COLORS
        }
        else if ((length(colors) > 1) & (length(colors) != NCAT)) {
            stop("'colors' must contain either a single value, or one value per mutation type.")
        }
        else {
            bar_colors <- colors
        }
        bars <- barplot(spec[i, ], names.arg = colnames(spec), axisnames = FALSE,
                        mgp = c(3, 0.3, 0), col = bar_colors, border = "white",
                        ylim = c(0, samp_max_y), xlim = xlim, yaxt = "n", xaxs = "i",
                        las = 2, cex.names = 1.35, adj = 0.5)
        axis(side = 1, at = bars, labels = colnames(spec), tick = FALSE,
             mgp = c(3, 0.3, 0), las = 2, cex.axis = 1.15, gap.axis = -1)
        group_end <- cumsum(GROUP_WIDTHS)
        group_start <- c(1, head(group_end, -1) + 1)
        group_xl <- bars[group_start] - 0.5
        group_xr <- bars[group_end] + 0.5
        group_xm <- (group_xl + group_xr) / 2
        type_end <- cumsum(TYPE_WIDTHS)
        type_start <- c(1, head(type_end, -1) + 1)
        type_xm <- ((bars[type_start] - 0.5) + (bars[type_end] + 0.5)) / 2
        # Plot axis
        if (any(spec > 1)) {
            axis(side = 2, cex.axis = 1.9, lwd = 2)
            label <- "Mutations"
            n_text <- paste0(" (", prettyNum(sum(spec[i,]), big.mark = ","), " mutations)")
        }
        else {
            axis(side = 2, at = seq(0, samp_max_y, ifelse(samp_max_y > 0.25, 0.1, 0.05)),
                 cex.axis = 1.9, lwd = 2)
            label <- "Mutation probability"
            n_text <- ""
        }
        if (is.null(name)) {
            nme <- rownames(spec)[i]
        }
        else {
            nme <- name
        }
        mtext(label, side = 2, cex = 2.4, line = 3.5)
        title(paste0(nme, n_text), line = 7.8, cex.main = 2.5)
        # Plot HPD intervals
        if (!is.null(lwr)) {
            arrows(bars, upr[i,],
                   bars, lwr[i,],
                   length = 0, lwd = 3, col = LINECOL)
        }
        # Plot mutation type labels
        rect(xleft = group_xl, xright = group_xr, ybottom = 0.945 * samp_max_y,
             ytop = samp_max_y, col = GROUP_COLORS, border = "white")
        mtext(TYPES, side = 3, at = type_xm, line = 2.4, cex = 1.75, xpd = TRUE)
        mtext(GROUP_LABELS, at = group_xm, side = 3, line = 0.35, cex = 1.75)
        mtext(c("Homopolymer length", "Homopolymer length",
                "Number of repeat units", "Number of repeat units",
                "Microhomology length"),
              side = 1, at = type_xm + c(0, 0, 0, 0, 0.2), line = 2.8,
              cex = c(rep(1.55, 4), 1.45))
        # Plot box
        if (boxes) {
            box(lwd = 2)
        }
    }
}


#' Plot generic spectrum
plot_spectrum_generic <- function(spec, lwr, upr, name, max_y, colors, boxes) {
    COLOR <- "orangered3"
    LINECOL <- "gray60"
    FACTOR <- 1.05
    NCAT <- ncol(spec)   # number of categories
    NSAMP <- nrow(spec)  # number of samples
    if (is.null(colnames(spec))) {
        types <- paste("Mut. type", 1:NCAT)
    }
    else {
        types <- colnames(spec)
    }

    for (i in 1:NSAMP) {
        if (is.null(max_y)) {
            samp_max_y <- ifelse(is.null(upr), max(spec[i,]) * FACTOR, max(upr[i,]) * FACTOR)
        }
        else {
            samp_max_y <- max_y
        }
        if (is.null(colors)) {
            colors = COLOR
        }
        else if ((length(colors) > 1) & (length(colors) != NCAT)) {
            stop("'colors' must contain either a single value, or one value per mutation type.")
        }
        # Plot spectrum bars
        bars <- barplot(spec[i,], names.arg = types, col = colors, mgp = c(3, 0.8, 0),
                        border = "white", las = 2, cex.names = 1,
                        ylim = c(0, samp_max_y), yaxt = "n", xaxs="i")
        # Plot axis
        if (any(spec > 1)) {
            axis(side = 2, cex.axis = 1.9, lwd = 2)
            label <- "Mutations"
            n_text <- paste0(" (", prettyNum(sum(spec[i,]), big.mark = ","), " mutations)")
        }
        else {
            axis(side = 2, cex.axis = 1.9, lwd = 2)
            label <- "Mutation probability"
            n_text <- ""
        }
        if (is.null(name)) {
            nme <- rownames(spec)[i]
        }
        else {
            nme <- name
        }
        mtext(label, side = 2, cex = 2.4, line = 3.5)
        title(paste0(nme, n_text), line = 4, cex.main = 2.5)
        # Plot HPD intervals
        if (!is.null(lwr)) {
            arrows(bars, upr[i,],
                   bars, lwr[i,],
                   length = 0, lwd = 3, col = LINECOL)
        }
        # Plot box
        if (boxes) {
            box(lwd = 2)
        }
    }
}
