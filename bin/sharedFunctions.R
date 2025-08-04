######################
### Custom shared functions
######################

#Function to load and format VCF
# regions: optional tibble with columns: seqnames start_refspace end_refspace
load_vcf <- function(vcf_file, regions = NULL, genome_fasta, BSgenome_name, bcftools_bin){

  #Check if GT and GQ tags exist
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
				AD1 = integer(),
				AD2 = integer(),
				GT = character(),
				GQ = numeric(),
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
	
	#Filter for records containing ALT alleles, and atomize and split multi-allelic sites (bcftools norm -a -f [fastaref] | bcftools norm -m -both -f [fastaref])
  tmpvcf <- tempfile(tmpdir=getwd(),pattern=".")
  
  system(paste("/bin/bash -c",shQuote(paste(
    bcftools_bin,"view",
    if(GT_exists){"-i 'GT=\"alt\"'"},
    if(!is.null(regions)){paste("-R",tmpregions)},
    vcf_file,"|",
    bcftools_bin,"norm -a -f",genome_fasta,"2>/dev/null |",
    bcftools_bin,"norm -m -both -f",genome_fasta,"2>/dev/null >",
    tmpvcf
    )
   )))
  
  if(!is.null(regions)){
  	file.remove(tmpregions) %>% invisible
  }

  #Load vcf file
   #Remove ALT == "*" alleles that indicate overlapping deletions (not needed because every deletion already has a separate vcf entry; seen in germline VCF files) and ALT = "<*>" alleles that are non-variant gVCF blocks (seen in bcftools mpileup files).
   #Annotate with allele depths (AD1, AD2), Depth (AD1+AD2) and VAF (AD2 / Depth).
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
  		AD=as.character(extract.gt(vcf,element="AD",IDtoRowNames=FALSE)),
  		GT = if(GT_exists){
  			as.character(extract.gt(vcf,element="GT",IDtoRowNames=FALSE))
  		}else{NULL},
  		GQ = if(GQ_exists){
  			as.numeric(extract.gt(vcf,element="GQ",as.numeric=TRUE,IDtoRowNames=FALSE))
  		}else{NULL}
  	) %>%
    filter(! ALT %in% c("*","<*>")) %>%
    separate(AD,c("AD1","AD2"),",",convert=TRUE) %>%
    mutate(Depth=AD1+AD2, VAF=AD2/Depth) %>%
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
    makeGRangesFromDataFrame(
      seqnames.field="seqnames",
      start.field="start_refspace",
      end.field="end_refspace",
      keep.extra.columns=TRUE,
      seqinfo=BSgenome_name %>% get %>% seqinfo
    ) %>%
  	sort %>% #Sort and select unique rows since -R regions can lead to out of order and duplicated entries
  	unique

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