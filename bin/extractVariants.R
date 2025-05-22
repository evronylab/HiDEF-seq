#!/usr/bin/env -S Rscript --vanilla

#extractVariants.R:
# Loads and formats aligned ccs bamFile in HiDEF-seq format RDS file that includes all required alignment and variant information for analysis.  
# Usage: extractVariants.R [bamFile] [configuration.yaml]

cat("#### Running extractVariants...\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs))

######################
### Load configuration
######################
cat("#### Loading configuration...")

#Command line arguments
args <- commandArgs(trailingOnly = TRUE)
bamFile <- args[1]
yaml.config <- suppressWarnings(read.config(args[2]))

#Load BSgenome reference, and install it first if needed
if(!yaml.config$BSgenome$BSgenome_name %in% installed.genomes()){
	if(yaml.config$BSgenome$BSgenome_name %in% available.genomes()){
		BiocManager::install(yaml.config$BSgenome$BSgenome_name)
	}else if(!is.null(yaml.config$BSgenome$BSgenome_file)){
		dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
		install.packages(yaml.config$BSgenome$BSgenome_file, repos = NULL,  lib = Sys.getenv("R_LIBS_USER"))
	}else{
		stop("ERROR: Must specify either BSgenome_name that is in available.genomes() or a BSgenome_file!", call.=FALSE)
	}
}
suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name,character.only=TRUE))

#Comma-separated chromosomes to analyze across all chromgroups
chroms_to_analyze <- yaml.config$chromgroups %>%
  map_vec("chroms") %>%
  map(str_split_1,",") %>%
  unlist

#Output stats
stats <- list()

cat("DONE\n")

######################
### Load BAM file
######################
cat("#### Loading BAM file:",bamFile,"...")

bam <- bamFile %>%
  scanBam(
    param=ScanBamParam(
      what=setdiff(scanBamWhat(),c("mrnm","mpos","groupid","mate_status")),
      tag=c("ec","np","rq","zm","sa","sm","sx")
    )
  ) %>%
  pluck(1)

#Create DataFrame of reads
bam.df <- cbind(
  DataFrame(bam[names(bam) != "tag"] %>% lapply(I)),
  DataFrame(bam$tag %>% lapply(I))
)
rm(bam)

#Extract ccs strand
bam.df$ccs_strand <- bam.df$qname %>% str_extract("(fwd|rev)$")

# Convert to GRanges
bam.df$end <- bam.df$pos + bam.df$isize - 1

bam.gr <- bam.df %>%
  makeGRangesFromDataFrame(
    seqnames.field="rname",
    start.field="pos",
    end.field="end",
    strand.field="strand",
    keep.extra.columns=TRUE,
    seqinfo=seqinfo(get(yaml.config$BSgenome$BSgenome_name))
  )
rm(bam.df)

#Count number of molecules
stats$num_molecules_initial <- bam.gr$zm %>% n_distinct

cat("DONE\n")

######################
### Initial molecule filtering
######################
cat("#### Initial molecule filtering...")

# Keep only molecules with 1 forward and 1 reverse primary alignment on the same chromosome (flags = 0 and 16 are plus and minus strand primary alignments, respectively; flags =  2048 and 2064 are plus and minus strand supplementary alignments, respectively).
zmwstokeep <- bam.gr %>%
  as_tibble %>%
  select(zm,ccs_strand,strand,flag,seqnames) %>%
  group_by(zm) %>%
  filter(
    n() == 2,
    sum(ccs_strand=="fwd") == 1,
    sum(ccs_strand=="rev") == 1,
    sum(strand=="+" & flag==0) == 1,
    sum(strand=="-" & flag==16) == 1,
    n_distinct(seqnames) == 1
    ) %>%
  pluck("zm")

bam.gr <- bam.gr[bam.gr$zm %in% zmwstokeep,]
rm(zmwstokeep)

# Count number of molecules
stats$num_molecules_postalignmentfilter <- bam.gr$zm %>% n_distinct

# Keep only molecules aligned to analyzed chromosomes
bam.gr <- bam.gr[seqnames(bam.gr) %in% chroms_to_analyze]

# Count number of molecules
stats$num_molecules_postanalyzedchroms <- bam.gr$zm %>% n_distinct

# Arrange bam.gr by ZMW id and strand
bam.gr <- bam.gr[order(bam.gr$zm, strand(bam.gr))]

# Keep only molecules with plus and minus strand alignment overlap >= min_strand_overlap (reciprocal or both plus and minus strand alignments)

 # Create simplified GRanges to reduce memory use
bam.gr.onlyranges <- GRanges(
  seqnames = seqnames(bam.gr),
  ranges = ranges(bam.gr),
  strand = strand(bam.gr),
  zm = bam.gr$zm
  )

 # Extract plus and minus strand reads separately. Due to prior sorting by ZMW id and strand, these are guaranteed to have the same ZMW id order.
bam.gr.onlyranges.plus <- bam.gr.onlyranges[strand(bam.gr.onlyranges) == "+"]
bam.gr.onlyranges.minus <- bam.gr.onlyranges[strand(bam.gr.onlyranges) == "-"]

 # Confirm ZMW id order is the same for plus and minus strand reads
if(! identical(bam.gr.onlyranges.plus$zm, bam.gr.onlyranges.minus$zm)){
  stop("ERROR: ZMW ids of plus and minus strand reads do not match!")
}

 # Compute reciprocal overlap fractions
bam.gr.onlyranges.overlap <- pintersect(bam.gr.onlyranges.plus, bam.gr.onlyranges.minus, ignore.strand=TRUE)
bam.gr.onlyranges.overlap.plus  <- width(bam.gr.onlyranges.overlap) / width(bam.gr.onlyranges.plus)
bam.gr.onlyranges.overlap.minus  <- width(bam.gr.onlyranges.overlap) / width(bam.gr.onlyranges.minus)

zmwstokeep <- bam.gr.onlyranges.overlap$zm[
  (bam.gr.onlyranges.overlap.plus >= yaml.config$min_strand_overlap) &
  (bam.gr.onlyranges.overlap.minus >= yaml.config$min_strand_overlap)
  ]

 # Keep only molecules passing reciprocal overlap filter
bam.gr <- bam.gr[bam.gr$zm %in% zmwstokeep,]

 # Remove intermediate objects
rm(bam.gr.onlyranges, bam.gr.onlyranges.plus, bam.gr.onlyranges.minus, bam.gr.onlyranges.overlap, bam.gr.onlyranges.overlap.plus, bam.gr.onlyranges.overlap.minus, zmwstokeep)

# Count number of molecules
stats$num_molecules_postminstrandoverlap <- bam.gr$zm %>% n_distinct

cat("DONE\n")

######################
### Extract read indels
######################
cat("#### Extracting read indels...")

#Extract indel positions in query space, and name them with the position to retain this in the final data frame.
cigar_query_indels <- cigarRangesAlongQuerySpace(bam.gr$cigar,ops=c("I","D"))
cigar_query_indels <- cigar_query_indels %>%
	unlist(use.names=FALSE) %>%
	setNames(str_c(start(.),end(.),sep="_")) %>%
	relist(.,cigar_query_indels)

#Check if no indels
if(cigar_query_indels %>% as.data.frame %>% nrow == 0){
	cat("  No indels in selected chromosomes!\n")
	indels.df <- tibble(zm=integer(),seqnames=factor(),strand=factor(),start_queryspace=integer(),end_queryspace=integer(),start_refspace=integer(),end_refspace=integer(),ref=character(),alt=character(),qual=integer(),sa=integer(),sm=integer(),sx=integer())
	
}else{

	#Indel positions in query space. Note, in query space, deletion width = 0.
	indels_queryspace <- cigar_query_indels %>%
		setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_")) %>%
		as.data.frame %>%
		separate(group_name,sep="_",into=c("zm","strand")) %>%
		mutate(
			zm=as.integer(zm),
			strand = as.factor(strand)
		) %>%
		rename(
			start_queryspace = start,
			end_queryspace = end,
			insertion_width = width
			) %>%
		select(-group,-names)
	
	#Indel positions in reference space, and also annotate with positions in query space for later joining. Note, in reference space, insertion width = 0.
	indels_refspace <- cigarRangesAlongReferenceSpace(bam.gr$cigar,ops=c("I","D"),pos=start(bam.gr))
	indels_refspace <- indels_refspace %>%
		unlist(use.names=FALSE) %>%
		setNames(str_c(cigar_query_indels %>% unlist(use.names=FALSE) %>% start,cigar_query_indels %>% unlist(use.names=FALSE) %>% end,sep="_")) %>%
		relist(.,indels_refspace) %>%
		setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_")) %>%
		as.data.frame %>%
		separate(group_name,sep="_",into=c("zm","strand")) %>%
		separate(names,sep="_",into=c("start_queryspace","end_queryspace")) %>%
		rename(
			start_refspace = start,
			end_refspace = end
		) %>%
		mutate(
			zm=as.integer(zm),
			strand = as.factor(strand),
			start_queryspace = as.integer(start_queryspace),
			end_queryspace = as.integer(end_queryspace),
			deletion_width = width
		) %>%
		left_join(
			data.frame(seqnames=seqnames(bam.gr), zm=bam.gr$zm) %>% distinct,
			by = join_by(zm)
		) %>%
		select(-group)
	
	#Insertion base sequences from query space, relative to reference plus strand
	indels_refspace$ref <- indels_refspace %>%
		mutate(strand="+") %>%
		makeGRangesFromDataFrame(
			seqnames.field="seqnames",
			start.field="start_refspace",
			end.field="start_refspace",
			strand.field="strand"
		) %>%
		getSeq(eval(parse(text=yaml.config$BSgenome$BSgenome_name)),.) %>%
		as.character
	
	#Single base substitution ALT base sequences, relative to reference plus strand
	sbs_alt <- extractAt(bam.gr$seq,cigar_query_sbs) %>%
		as("CharacterList") %>%
		setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_")) %>%
		as.list %>%
		enframe %>%
		unnest_longer(col=value,indices_to="start_queryspace",values_to="alt") %>%
		separate_longer_position(alt,width=1) %>%
		separate(name,sep="_",into=c("zm","strand")) %>%
		mutate(
			zm = as.integer(zm),
			strand = as.factor(strand),
			start_queryspace = as.integer(start_queryspace)
		) %>%
		group_by(zm,strand,start_queryspace) %>%
		mutate(start_queryspace = start_queryspace + row_number() - 1L) %>%
		ungroup
	
	#Single base substitution base qualities
	sbs_qual <- extractAt(bam.gr$qual,cigar_query_sbs) %>%
		as("CharacterList") %>%
		setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_")) %>%
		as.list %>%
		enframe %>%
		unnest_longer(col=value,indices_to="start_queryspace",values_to="qual") %>%
		separate_longer_position(qual,width=1) %>%
		separate(name,sep="_",into=c("zm","strand")) %>%
		mutate(
			zm = as.integer(zm),
			strand = as.factor(strand),
			start_queryspace = as.integer(start_queryspace),
			qual = qual %>% PhredQuality %>% as("IntegerList") %>% unlist
		) %>%
		group_by(zm,strand,start_queryspace) %>%
		mutate(start_queryspace = start_queryspace + row_number() - 1L) %>%
		ungroup
	
	rm(cigar_query_sbs)
	
	#sa/sm/sx: subread alignment tags
	
	#Convert sbs positions in query space to list for faster retrieval below of sa, sm, sx tags
	sbs_queryspace.list <- sbs_queryspace %>%
		mutate(zm_strand = str_c(zm,strand,sep="_")) %>%
		select(zm_strand,start_queryspace) %>%
		as.data.table %>%
		split(by="zm_strand",keep.by=F) %>%
		lapply(unlist,use.names=F)
	
	#sa: the number of subread alignments that span each CCS read position
	# sa is run-length encoded (rle) in the BAM file as length,value pairs. We inverse the rle encoding back to standard per-position values
	# sa needs to be reversed for reads aligned to genome minus strand
	
	#Extract all sa values and inverse rle.
	sa.lengths <- lapply(bam.gr$sa,function(x){x[c(TRUE, FALSE)]})
	sa.values <- lapply(bam.gr$sa,function(x){x[c(FALSE, TRUE)]})
	sa <- mapply(
		function(x,y){inverse.rle(list(lengths=x,values=y))},
		x=sa.lengths,
		y=sa.values,
		USE.NAMES=FALSE
	) %>%
		setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_"))
	rm(sa.lengths,sa.values)
	
	#Reverse order for reads aligned to genome minus strand
	sa[strand(bam.gr) %>% as.vector == "-"] <- lapply(sa[strand(bam.gr) %>% as.vector == "-"], rev)
	
	#Extract sa for all single base substitutions
	sbs_sa <- map2(sa[sbs_queryspace.list %>% names], sbs_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
		enframe %>%
		separate(name,sep="_",into=c("zm","strand")) %>%
		unnest_longer(col=value, values_to="sa", indices_to="start_queryspace") %>%
		mutate(
			zm = as.integer(zm),
			strand = as.factor(strand),
			start_queryspace = as.integer(start_queryspace)
		)
	
	rm(sa)
	
	#sm: the number of subreads that align as a match to each CCS read position
	# sm needs to be reversed for reads aligned to genome minus strand
	sm <- bam.gr$sm %>%
		setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_"))
	sm[strand(bam.gr) %>% as.vector == "-"] <- lapply(sm[strand(bam.gr) %>% as.vector == "-"], rev)
	
	sbs_sm <- map2(sm[sbs_queryspace.list %>% names], sbs_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
		enframe %>%
		separate(name,sep="_",into=c("zm","strand")) %>%
		unnest_longer(col=value, values_to="sm", indices_to="start_queryspace") %>%
		mutate(
			zm = as.integer(zm),
			strand = as.factor(strand),
			start_queryspace = as.integer(start_queryspace)
		)
	
	rm(sm)
	
	#sx: the number of subreads that align as a sbs to each CCS read position
	# sx needs to be reversed for reads aligned to genome minus strand
	sx <- bam.gr$sx %>%
		setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_"))
	sx[strand(bam.gr) %>% as.vector == "-"] <- lapply(sx[strand(bam.gr) %>% as.vector == "-"], rev)
	
	sbs_sx <- map2(sx[sbs_queryspace.list %>% names], sbs_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
		enframe %>%
		separate(name,sep="_",into=c("zm","strand")) %>%
		unnest_longer(col=value, values_to="sx", indices_to="start_queryspace") %>%
		mutate(
			zm = as.integer(zm),
			strand = as.factor(strand),
			start_queryspace = as.integer(start_queryspace)
		)
	
	rm(sx)
	
	#Combine sbs info
	sbs.df <- sbs_queryspace %>%
		left_join(sbs_alt, by=join_by(zm,strand,start_queryspace)) %>%
		left_join(sbs_qual, by=join_by(zm,strand,start_queryspace)) %>%
		left_join(sbs_refspace, by=join_by(zm,strand,start_queryspace)) %>%
		left_join(sbs_sa, by=join_by(zm,strand,start_queryspace)) %>%
		left_join(sbs_sm, by=join_by(zm,strand,start_queryspace)) %>%
		left_join(sbs_sx, by=join_by(zm,strand,start_queryspace)) %>%
		select(zm,seqnames,strand,start_queryspace,start_refspace,ref,alt,qual,sa,sm,sx)
	
	rm(sbs_queryspace,sbs_queryspace.list,sbs_alt,sbs_qual,sbs_refspace,sbs_sa,sbs_sm,sbs_sx)
	
	#Annotate each SBS as a change from the reference in reference space in only one strand ("SBS_mismatch"), non-complementary changes from the reference in both strands ("SBS_mismatch_ds"), or complementary changes from the reference in both strands.
	sbs.df <- sbs.df %>%
		group_by(zm,start_refspace) %>%
		mutate(
			count = n(),
			count_distinct_alt = n_distinct(alt),
			variant_type = case_when(
				count == 1 ~ "SBS_mismatch",
				count == 2 & count_distinct_alt == 1 ~ "SBS_mutation",
				count == 2 & count_distinct_alt == 2 ~ "SBS_mismatch_ds"
			) %>%
				as.factor
		) %>%
		select(-count,-count_distinct_alt) %>%
		ungroup
}

cat("DONE\n")

######################
### Extract read single base substitutions (sbs)
######################
cat("#### Extracting read single base substitutions...")

#Extract substitution positions in query space, and name them with the position to retain this in the final data frame
cigar_query_sbs <- cigarRangesAlongQuerySpace(bam.gr$cigar,ops="X")
cigar_query_sbs <- cigar_query_sbs %>%
  unlist(use.names=FALSE) %>%
  setNames(start(.)) %>%
  relist(.,cigar_query_sbs)

#Check if no single base substitutions
if(cigar_query_sbs %>% as.data.frame %>% nrow == 0){
  cat("  No single base substitutions in selected chromosomes!\n")
  sbs.df <- tibble(zm=integer(),seqnames=factor(),strand=factor(),start_queryspace=integer(),start_refspace=integer(),ref=character(),alt=character(),qual=integer(),sa=integer(),sm=integer(),sx=integer())
  
}else{
  #Single base substitution positions in query space
  sbs_queryspace <- cigar_query_sbs %>%
    setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_")) %>%
    as.data.frame %>%
    separate(group_name,sep="_",into=c("zm","strand")) %>%
    mutate(
      zm=as.integer(zm),
      strand = as.factor(strand)
    ) %>%
    rename(start_queryspace = start) %>%
    uncount(weights=width) %>%
    select(zm,strand,start_queryspace) %>%
    group_by(zm,strand,start_queryspace) %>%
    mutate(start_queryspace = start_queryspace + row_number() - 1L) %>%
    ungroup
  
  #Single base substitution positions in reference space, and also annotate with positions in query space for later joining
  sbs_refspace <- cigarRangesAlongReferenceSpace(bam.gr$cigar,ops="X",pos=start(bam.gr))
  sbs_refspace <- sbs_refspace %>%
    unlist(use.names=FALSE) %>%
    setNames(cigar_query_sbs %>% unlist(use.names=FALSE) %>% start) %>%
    relist(.,sbs_refspace) %>%
    setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_")) %>%
    as.data.frame %>%
    separate(group_name,sep="_",into=c("zm","strand")) %>%
    mutate(
      zm=as.integer(zm),
      strand = as.factor(strand),
      names = as.integer(names)
    ) %>%
    rename(
      start_refspace = start,
      start_queryspace = names
    ) %>%
    uncount(weights=width) %>%
    select(zm,strand,start_refspace,start_queryspace) %>%
    group_by(zm,strand,start_refspace,start_queryspace) %>%
    mutate(
      start_refspace = start_refspace + row_number() - 1L,
      start_queryspace = start_queryspace + row_number() - 1L
    ) %>%
    ungroup %>%
    left_join(
      data.frame(seqnames=seqnames(bam.gr), zm=bam.gr$zm) %>% distinct,
      by = join_by(zm)
    )
  
  #Single base substitution REF base sequences, relative to reference plus strand
  sbs_refspace$ref <- sbs_refspace %>%
    mutate(strand="+") %>%
    makeGRangesFromDataFrame(
      seqnames.field="seqnames",
      start.field="start_refspace",
      end.field="start_refspace",
      strand.field="strand"
    ) %>%
    getSeq(eval(parse(text=yaml.config$BSgenome$BSgenome_name)),.) %>%
    as.character

  #Single base substitution ALT base sequences, relative to reference plus strand
  sbs_alt <- extractAt(bam.gr$seq,cigar_query_sbs) %>%
    as("CharacterList") %>%
    setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_")) %>%
    as.list %>%
    enframe %>%
    unnest_longer(col=value,indices_to="start_queryspace",values_to="alt") %>%
    separate_longer_position(alt,width=1) %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    mutate(
      zm = as.integer(zm),
      strand = as.factor(strand),
      start_queryspace = as.integer(start_queryspace)
      ) %>%
    group_by(zm,strand,start_queryspace) %>%
    mutate(start_queryspace = start_queryspace + row_number() - 1L) %>%
    ungroup
  
  #Single base substitution base qualities
  sbs_qual <- extractAt(bam.gr$qual,cigar_query_sbs) %>%
    as("CharacterList") %>%
    setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_")) %>%
    as.list %>%
    enframe %>%
    unnest_longer(col=value,indices_to="start_queryspace",values_to="qual") %>%
    separate_longer_position(qual,width=1) %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    mutate(
      zm = as.integer(zm),
      strand = as.factor(strand),
      start_queryspace = as.integer(start_queryspace),
      qual = qual %>% PhredQuality %>% as("IntegerList") %>% unlist
    ) %>%
    group_by(zm,strand,start_queryspace) %>%
    mutate(start_queryspace = start_queryspace + row_number() - 1L) %>%
    ungroup
  
  rm(cigar_query_sbs)

  #sa/sm/sx: subread alignment tags

  #Convert sbs positions in query space to list for faster retrieval below of sa, sm, sx tags
  sbs_queryspace.list <- sbs_queryspace %>%
    mutate(zm_strand = str_c(zm,strand,sep="_")) %>%
    select(zm_strand,start_queryspace) %>%
    as.data.table %>%
    split(by="zm_strand",keep.by=F) %>%
    lapply(unlist,use.names=F)
  
  #sa: the number of subread alignments that span each CCS read position
  # sa is run-length encoded (rle) in the BAM file as length,value pairs. We inverse the rle encoding back to standard per-position values
  # sa needs to be reversed for reads aligned to genome minus strand

   #Extract all sa values and inverse rle.
  sa.lengths <- lapply(bam.gr$sa,function(x){x[c(TRUE, FALSE)]})
  sa.values <- lapply(bam.gr$sa,function(x){x[c(FALSE, TRUE)]})
  sa <- mapply(
      function(x,y){inverse.rle(list(lengths=x,values=y))},
      x=sa.lengths,
      y=sa.values,
      USE.NAMES=FALSE
    ) %>%
    setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_"))
  rm(sa.lengths,sa.values)
  
   #Reverse order for reads aligned to genome minus strand
  sa[strand(bam.gr) %>% as.vector == "-"] <- lapply(sa[strand(bam.gr) %>% as.vector == "-"], rev)
  
   #Extract sa for all single base substitutions
  sbs_sa <- map2(sa[sbs_queryspace.list %>% names], sbs_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
    enframe %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    unnest_longer(col=value, values_to="sa", indices_to="start_queryspace") %>%
    mutate(
      zm = as.integer(zm),
      strand = as.factor(strand),
      start_queryspace = as.integer(start_queryspace)
    )
  
  rm(sa)
  
  #sm: the number of subreads that align as a match to each CCS read position
  # sm needs to be reversed for reads aligned to genome minus strand
  sm <- bam.gr$sm %>%
    setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_"))
  sm[strand(bam.gr) %>% as.vector == "-"] <- lapply(sm[strand(bam.gr) %>% as.vector == "-"], rev)
  
  sbs_sm <- map2(sm[sbs_queryspace.list %>% names], sbs_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
    enframe %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    unnest_longer(col=value, values_to="sm", indices_to="start_queryspace") %>%
    mutate(
      zm = as.integer(zm),
      strand = as.factor(strand),
      start_queryspace = as.integer(start_queryspace)
    )
  
  rm(sm)
  
  #sx: the number of subreads that align as a sbs to each CCS read position
  # sx needs to be reversed for reads aligned to genome minus strand
  sx <- bam.gr$sx %>%
    setNames(str_c(bam.gr$zm,strand(bam.gr) %>% as.character,sep="_"))
  sx[strand(bam.gr) %>% as.vector == "-"] <- lapply(sx[strand(bam.gr) %>% as.vector == "-"], rev)
  
  sbs_sx <- map2(sx[sbs_queryspace.list %>% names], sbs_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
    enframe %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    unnest_longer(col=value, values_to="sx", indices_to="start_queryspace") %>%
    mutate(
      zm = as.integer(zm),
      strand = as.factor(strand),
      start_queryspace = as.integer(start_queryspace)
    )
  
  rm(sx)
  
  #Combine sbs info
  sbs.df <- sbs_queryspace %>%
    left_join(sbs_alt, by=join_by(zm,strand,start_queryspace)) %>%
    left_join(sbs_qual, by=join_by(zm,strand,start_queryspace)) %>%
    left_join(sbs_refspace, by=join_by(zm,strand,start_queryspace)) %>%
    left_join(sbs_sa, by=join_by(zm,strand,start_queryspace)) %>%
    left_join(sbs_sm, by=join_by(zm,strand,start_queryspace)) %>%
    left_join(sbs_sx, by=join_by(zm,strand,start_queryspace)) %>%
    select(zm,seqnames,strand,start_queryspace,start_refspace,ref,alt,qual,sa,sm,sx)
  
  rm(sbs_queryspace,sbs_queryspace.list,sbs_alt,sbs_qual,sbs_refspace,sbs_sa,sbs_sm,sbs_sx)
  
  #Annotate each SBS as a change from the reference in reference space in only one strand ("SBS_mismatch"), non-complementary changes from the reference in both strands ("SBS_mismatch_ds"), or complementary changes from the reference in both strands.
  sbs.df <- sbs.df %>%
    group_by(zm,start_refspace) %>%
    mutate(
      count = n(),
      count_distinct_alt = n_distinct(alt),
      variant_type = case_when(
        count == 1 ~ "SBS_mismatch",
        count == 2 & count_distinct_alt == 1 ~ "SBS_mutation",
        count == 2 & count_distinct_alt == 2 ~ "SBS_mismatch_ds"
      ) %>%
        as.factor
    ) %>%
    select(-count,-count_distinct_alt) %>%
    ungroup
}

cat("DONE\n")

######################
### Save output
######################
cat("#### Saving output...")

qsave(
  list(bam.gr = bam.gr, sbs.df = sbs.df, stats = stats),
  bamFile %>% basename %>% str_remove(".bam$") %>% str_c(".RDS")
  )

cat("DONE\n")