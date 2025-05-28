#!/usr/bin/env -S Rscript --vanilla

#extractVariants.R:
# Loads and formats aligned ccs bamFile in HiDEF-seq format RDS file that includes all required alignment and variant information for analysis.  
# Usage: extractVariants.R [bamFile] [configuration.yaml]

cat("#### Running extractVariants ####\n")

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
suppressPackageStartupMessages(library(qs2))

######################
### Load configuration
######################
cat("#### Loading configuration...")

#General options
options(datatable.showProgress = FALSE)
strand_levels <- c("+","-")

#Command line arguments
args <- commandArgs(trailingOnly = TRUE)
bamFile <- args[1]
yaml.config <- suppressWarnings(read.config(args[2]))

#Load the BSgenome reference
suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name,character.only=TRUE,lib.loc=yaml.config$prepareFilters_cache_dir))

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

#Get run metadata: movie, read group (rg), and run IDs from BAM file and config file
movie_ids <- bamFile %>%
	scanBamHeader %>%
	pluck(1,"text") %>%
	keep(names(.)=="@RG") %>%
	map(str_subset,"^PU:") %>%
	str_remove("^PU:")

rg_ids <- bamFile %>%
	scanBamHeader %>%
	pluck(1,"text") %>%
	keep(names(.)=="@RG") %>%
	map(str_subset,"^ID:") %>%
	str_remove("^ID:")

run_ids <- movie_ids %>%
	map_chr(
		function(x){
			map_lgl(yaml.config$runs, function(y){
				y$reads_file %>% str_detect(x)
			}) %>%
				keep(yaml.config$runs,.) %>%
				pluck(1,"run_id")
		}
	)

if(length(movie_ids) != length (run_ids)){
	stop("ERROR: Movie IDs in BAM file do not match yaml configuration file!")
}

run_metadata <- tibble(
	movie_id = factor(movie_ids),
	rg_id = factor(rg_ids),
	run_id = factor(run_ids)
)

rm(movie_ids,rg_ids,run_ids)

#Load BAM file reads
bam <- bamFile %>%
  scanBam(
    param=ScanBamParam(
      what=setdiff(scanBamWhat(),c("mrnm","mpos","groupid","mate_status")),
      tag=c("ec","np","rq","zm","sa","sm","sx","RG")
    )
  ) %>%
  pluck(1)

#Create DataFrame of reads
bam.df <- cbind(
  DataFrame(bam[names(bam) != "tag"] %>% lapply(I)),
  DataFrame(bam$tag %>% lapply(I))
)
rm(bam)

#Add movie_ids and run_ids using run metadata, while maintaining DataFrame format
bam.df <- bam.df %>%
	cbind(run_metadata[match(bam.df$RG, run_metadata$rg_id),c("movie_id","run_id")])

#Extract ccs strand
bam.df$ccs_strand <- bam.df$qname %>% str_extract("(fwd|rev)$")

#Convert some columns to factor
for(i in c("flag","RG","movie_id","run_id","ccs_strand")){
	bam.df[,i] <- factor(bam.df[,i])
}

#Convert to GRanges
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

# Arrange bam.gr by run id, ZMW id and strand
bam.gr <- bam.gr[order(bam.gr$run_id,bam.gr$zm, strand(bam.gr))]

#Count number of molecules per run
stats <- mcols(bam.gr)[,c("run_id","movie_id","zm")] %>%
	as_tibble %>%
	group_by(run_id,movie_id) %>%
	summarize(
		num_molecules_initial = n_distinct(zm),
		.groups="drop"
	)

cat("DONE\n")

######################
### Initial molecule filtering
######################
cat("#### Initial molecule filtering...")

# Keep only molecules with 1 forward and 1 reverse primary alignment on the same chromosome (flags = 0 and 16 are plus and minus strand primary alignments, respectively; flags =  2048 and 2064 are plus and minus strand supplementary alignments, respectively).
zmwstokeep <- bam.gr %>%
  as_tibble %>%
  select(run_id,zm,ccs_strand,strand,flag,seqnames) %>%
  group_by(run_id,zm) %>%
  filter(
    n() == 2,
    sum(ccs_strand=="fwd") == 1,
    sum(ccs_strand=="rev") == 1,
    sum(strand=="+" & flag==0) == 1,
    sum(strand=="-" & flag==16) == 1,
    n_distinct(seqnames) == 1
    ) %>%
  select(run_id,zm) %>%
	distinct

bam.gr <- bam.gr[bam.gr$run_id %in% zmwstokeep$run_id & bam.gr$zm %in% zmwstokeep$zm,]
rm(zmwstokeep)

# Count number of molecules per run
stats <- stats %>%
	left_join(
		mcols(bam.gr)[,c("run_id","movie_id","zm")] %>%
		as_tibble %>%
		group_by(run_id,movie_id) %>%
		summarize(
			num_molecules_postalignmentfilter = n_distinct(zm),
			.groups="drop"
		),
		by = join_by(run_id,movie_id)
)

# Keep only molecules aligned to analyzed chromosomes
bam.gr <- bam.gr[seqnames(bam.gr) %in% chroms_to_analyze]

# Count number of molecules per run
stats <- stats %>%
	left_join(
		mcols(bam.gr)[,c("run_id","movie_id","zm")] %>%
			as_tibble %>%
			group_by(run_id,movie_id) %>%
			summarize(
				num_molecules_postanalyzedchroms = n_distinct(zm),
				.groups="drop"
			),
		by = join_by(run_id,movie_id)
	)

# Keep only molecules with plus and minus strand alignment overlap >= min_strand_overlap (reciprocal or both plus and minus strand alignments)

 # Create simplified GRanges to reduce memory use
bam.gr.onlyranges <- GRanges(
  seqnames = seqnames(bam.gr),
  ranges = ranges(bam.gr),
  strand = strand(bam.gr),
  zm = bam.gr$zm,
  run_id=bam.gr$run_id
  )

 # Extract plus and minus strand reads separately. Due to prior sorting by ZMW id and strand, these are guaranteed to have the same ZMW id order.
bam.gr.onlyranges.plus <- bam.gr.onlyranges[strand(bam.gr.onlyranges) == "+"]
bam.gr.onlyranges.minus <- bam.gr.onlyranges[strand(bam.gr.onlyranges) == "-"]

 # Confirm run id and ZMW id order is the same for plus and minus strand reads
if(! identical(bam.gr.onlyranges.plus$zm, bam.gr.onlyranges.minus$zm) | ! identical(bam.gr.onlyranges.plus$run_id, bam.gr.onlyranges.minus$run_id) ){
  stop("ERROR: ZMW ids of plus and minus strand reads do not match!")
}

 # Compute reciprocal overlap fractions
bam.gr.onlyranges.overlap <- pintersect(bam.gr.onlyranges.plus, bam.gr.onlyranges.minus, ignore.strand=TRUE)
bam.gr.onlyranges.overlap$plus_overlap_frac <- width(bam.gr.onlyranges.overlap) / width(bam.gr.onlyranges.plus)
bam.gr.onlyranges.overlap$minus_overlap_frac <- width(bam.gr.onlyranges.overlap) / width(bam.gr.onlyranges.minus)

 # Keep only molecules passing reciprocal overlap filter
zmwstokeep <- bam.gr.onlyranges.overlap %>%
	as_tibble %>%
	select(run_id,zm,plus_overlap_frac,minus_overlap_frac) %>%
	filter(
		plus_overlap_frac >= yaml.config$min_strand_overlap,
		minus_overlap_frac >= plus_overlap_frac
	)

bam.gr <- bam.gr[bam.gr$run_id %in% zmwstokeep$run_id & bam.gr$zm %in% zmwstokeep$zm,]

 # Remove intermediate objects
rm(bam.gr.onlyranges, bam.gr.onlyranges.plus, bam.gr.onlyranges.minus, bam.gr.onlyranges.overlap, zmwstokeep)

# Count number of molecules per run
stats <- stats %>%
	left_join(
		mcols(bam.gr)[,c("run_id","movie_id","zm")] %>%
			as_tibble %>%
			group_by(run_id,movie_id) %>%
			summarize(
				num_molecules_postminstrandoverlap = n_distinct(zm),
				.groups="drop"
			),
		by = join_by(run_id,movie_id)
	)

cat("DONE\n")

#########################################
### Loop to extract variants from each run_id/movie_id separately
#########################################

#Split bam.gr by run_id
bam.gr <- bam.gr %>% split(bam.gr$run_id)

#Create output lists
sbs.df <- list()
indels.df <- list()

for(i in bam.gr %>% names){

cat("> Processing run:",i,"\n")

######################
### Extract subread alignment tags
######################
cat(" ## Extracting subread alignment tags...")

#sa: the number of subread alignments that span each CCS read position
# sa is run-length encoded (rle) in the BAM file as length,value pairs. We inverse the rle encoding back to standard per-position values
# sa needs to be reversed for reads aligned to genome minus strand

#Extract all sa values and inverse rle.
sa.lengths <- lapply(bam.gr[[i]]$sa,function(x){x[c(TRUE, FALSE)]})
sa.values <- lapply(bam.gr[[i]]$sa,function(x){x[c(FALSE, TRUE)]})
sa <- mapply(
	function(x,y){inverse.rle(list(lengths=x,values=y))},
	x=sa.lengths,
	y=sa.values,
	USE.NAMES=FALSE
) %>%
	setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_"))
rm(sa.lengths,sa.values)

#Reverse order for reads aligned to genome minus strand
sa[strand(bam.gr[[i]]) %>% as.vector == "-"] <- lapply(sa[strand(bam.gr[[i]]) %>% as.vector == "-"], rev)

#sm: the number of subreads that align as a match to each CCS read position
# sm needs to be reversed for reads aligned to genome minus strand
sm <- bam.gr[[i]]$sm %>%
	setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_"))
sm[strand(bam.gr[[i]]) %>% as.vector == "-"] <- lapply(sm[strand(bam.gr[[i]]) %>% as.vector == "-"], rev)

#sx: the number of subreads that align as a sbs to each CCS read position
# sx needs to be reversed for reads aligned to genome minus strand
sx <- bam.gr[[i]]$sx %>%
	setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_"))
sx[strand(bam.gr[[i]]) %>% as.vector == "-"] <- lapply(sx[strand(bam.gr[[i]]) %>% as.vector == "-"], rev)

cat("DONE\n")

######################
### Extract single base substitutions
######################
cat(" ## Extracting single base substitutions...")

#Extract substitution positions in query space, and name them with the position to retain this in the final data frame
cigar_query_sbs <- cigarRangesAlongQuerySpace(bam.gr[[i]]$cigar,ops="X")
cigar_query_sbs <- cigar_query_sbs %>%
  unlist(use.names=FALSE) %>%
  setNames(start(.)) %>%
  relist(.,cigar_query_sbs)

#Check if no single base substitutions
if(cigar_query_sbs %>% as.data.frame %>% nrow == 0){
  cat("  No single base substitutions in selected chromosomes!\n")
  sbs.df[[i]] <- tibble(zm=integer(),seqnames=factor(),strand=factor(),start_queryspace=integer(),start_refspace=integer(),ref=character(),alt=character(),qual=integer(),sa=integer(),sm=integer(),sx=integer(),variant_type=factor())
  
}else{
	
  #Single base substitution positions in query space
  sbs_queryspace <- cigar_query_sbs %>%
    setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_")) %>%
    as.data.frame %>%
    separate(group_name,sep="_",into=c("zm","strand")) %>%
    mutate(
      zm = as.integer(zm),
      strand = factor(strand,levels=strand_levels)
    ) %>%
    rename(start_queryspace = start) %>%
    uncount(weights=width) %>%
    select(zm,strand,start_queryspace) %>%
    group_by(zm,strand,start_queryspace) %>%
    mutate(start_queryspace = start_queryspace + row_number() - 1L) %>%
    ungroup
  
  #Single base substitution positions in reference space, annotated with positions in query space for later joining
  sbs_refspace <- cigarRangesAlongReferenceSpace(bam.gr[[i]]$cigar,ops="X",pos=start(bam.gr[[i]]))
  sbs_refspace <- sbs_refspace %>%
    unlist(use.names=FALSE) %>%
    setNames(cigar_query_sbs %>% unlist(use.names=FALSE) %>% start) %>%
    relist(.,sbs_refspace) %>%
    setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_")) %>%
    as.data.frame %>%
    separate(group_name,sep="_",into=c("zm","strand")) %>%
    mutate(
      zm=as.integer(zm),
      strand = factor(strand,levels=strand_levels),
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
      data.frame(seqnames=seqnames(bam.gr[[i]]), zm=bam.gr[[i]]$zm) %>% distinct,
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
  sbs_alt <- extractAt(bam.gr[[i]]$seq,cigar_query_sbs) %>%
    as("CharacterList") %>%
    setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_")) %>%
    as.list %>%
    enframe %>%
    unnest_longer(col=value,indices_to="start_queryspace",values_to="alt") %>%
    separate_longer_position(alt,width=1) %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    mutate(
      zm = as.integer(zm),
      strand = factor(strand,levels=strand_levels),
      start_queryspace = as.integer(start_queryspace)
      ) %>%
    group_by(zm,strand,start_queryspace) %>%
    mutate(start_queryspace = start_queryspace + row_number() - 1L) %>%
    ungroup
  
  #Single base substitution base qualities
  sbs_qual <- extractAt(bam.gr[[i]]$qual,cigar_query_sbs) %>%
    as("CharacterList") %>%
    setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_")) %>%
    as.list %>%
    enframe %>%
    unnest_longer(col=value,indices_to="start_queryspace",values_to="qual") %>%
    separate_longer_position(qual,width=1) %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    mutate(
      zm = as.integer(zm),
      strand = factor(strand,levels=strand_levels),
      start_queryspace = as.integer(start_queryspace),
      qual = qual %>% PhredQuality %>% as("IntegerList") %>% unlist
    ) %>%
    group_by(zm,strand,start_queryspace) %>%
    mutate(start_queryspace = start_queryspace + row_number() - 1L) %>%
    ungroup
  
  rm(cigar_query_sbs)

  #Extract sa/sm/sx tags

   #Convert sbs positions in query space to list for faster retrieval below of sa, sm, sx tags
  sbs_queryspace.list <- sbs_queryspace %>%
    mutate(zm_strand = str_c(zm,strand,sep="_")) %>%
    select(zm_strand,start_queryspace) %>%
    as.data.table %>%
    split(by="zm_strand",keep.by=F) %>%
    lapply(unlist,use.names=F)
  
   #sa
  sbs_sa <- map2(sa[sbs_queryspace.list %>% names], sbs_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
    enframe %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    unnest_longer(col=value, values_to="sa", indices_to="start_queryspace") %>%
    mutate(
      zm = as.integer(zm),
      strand = factor(strand,levels=strand_levels),
      start_queryspace = as.integer(start_queryspace)
    )
  
   #sm
  sbs_sm <- map2(sm[sbs_queryspace.list %>% names], sbs_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
    enframe %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    unnest_longer(col=value, values_to="sm", indices_to="start_queryspace") %>%
    mutate(
      zm = as.integer(zm),
      strand = factor(strand,levels=strand_levels),
      start_queryspace = as.integer(start_queryspace)
    )
  
   #sx
  sbs_sx <- map2(sx[sbs_queryspace.list %>% names], sbs_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
    enframe %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    unnest_longer(col=value, values_to="sx", indices_to="start_queryspace") %>%
    mutate(
      zm = as.integer(zm),
      strand = factor(strand,levels=strand_levels),
      start_queryspace = as.integer(start_queryspace)
    )
  
  #Combine sbs info
  sbs.df[[i]] <- sbs_queryspace %>%
    left_join(sbs_alt, by=join_by(zm,strand,start_queryspace)) %>%
    left_join(sbs_qual, by=join_by(zm,strand,start_queryspace)) %>%
    left_join(sbs_refspace, by=join_by(zm,strand,start_queryspace)) %>%
    left_join(sbs_sa, by=join_by(zm,strand,start_queryspace)) %>%
    left_join(sbs_sm, by=join_by(zm,strand,start_queryspace)) %>%
    left_join(sbs_sx, by=join_by(zm,strand,start_queryspace)) %>%
    select(zm,seqnames,strand,start_queryspace,start_refspace,ref,alt,qual,sa,sm,sx)
  
  rm(sbs_queryspace,sbs_queryspace.list,sbs_alt,sbs_qual,sbs_refspace,sbs_sa,sbs_sm,sbs_sx)
  
  #Annotate each SBS as a change from the reference in reference space in only one strand ("SBS_mismatch-ss"), non-complementary changes from the reference in both strands ("SBS_mismatch-ds"), or complementary changes from the reference in both strands ("SBS_mutation").
  sbs.df[[i]] <- sbs.df[[i]] %>%
    group_by(zm,start_refspace) %>%
    mutate(
      count = n(),
      count_distinct_alt = n_distinct(alt),
      variant_type = case_when(
        count == 1 ~ "SBS_mismatch-ss",
        count == 2 & count_distinct_alt == 1 ~ "SBS_mutation",
        count == 2 & count_distinct_alt == 2 ~ "SBS_mismatch-ds"
      ) %>%
        factor
    ) %>%
    select(-count,-count_distinct_alt) %>%
    ungroup %>%
    arrange(zm,start_refspace,strand)
}

cat("DONE\n")

######################
### Extract indels
######################
cat(" ## Extracting indels...")

#Extract indel positions in query space, and name them with the position to retain this in the final data frame.
cigar_query_indels <- cigarRangesAlongQuerySpace(bam.gr[[i]]$cigar,ops=c("I","D"))
cigar_query_indels <- cigar_query_indels %>%
	unlist(use.names=FALSE) %>%
	setNames(str_c(start(.),end(.),sep="_")) %>%
	relist(.,cigar_query_indels)

#Check if no indels
if(cigar_query_indels %>% as.data.frame %>% nrow == 0){
	cat("  No indels in selected chromosomes!\n")
	indels.df[[i]] <- tibble(zm=integer(),seqnames=factor(),strand=factor(),start_queryspace=integer(),end_queryspace=integer(),start_refspace=integer(),end_refspace=integer(),indel_type=factor(),ref=character(),alt=character(),indel_width=integer(),qual=list(),sa=list(),sm=list(),sx=list(),variant_type=factor())
	
}else{
	
	#Indel positions in query space. Note, in query space, deletion width = 0.
	indels_queryspace <- cigar_query_indels %>%
		setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_")) %>%
		as.data.frame %>%
		separate(group_name,sep="_",into=c("zm","strand")) %>%
		mutate(
			zm = as.integer(zm),
			strand = factor(strand,levels=strand_levels)
		) %>%
		rename(
			start_queryspace = start,
			end_queryspace = end,
			insertion_width = width
		) %>%
		mutate(insertion_width = insertion_width %>% na_if(0)) %>%
		select(-group,-names) %>%
	  as_tibble
	
	#Indel positions in reference space, annotated with positions in query space for later joining. Note, in reference space, insertion width = 0.
	indels_refspace <- cigarRangesAlongReferenceSpace(bam.gr[[i]]$cigar,ops=c("I","D"),pos=start(bam.gr[[i]]))
	indels_refspace <- indels_refspace %>%
		unlist(use.names=FALSE) %>%
		setNames(str_c(cigar_query_indels %>% unlist(use.names=FALSE) %>% start,cigar_query_indels %>% unlist(use.names=FALSE) %>% end,sep="_")) %>%
		relist(.,indels_refspace) %>%
		setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_")) %>%
		as.data.frame %>%
		separate(group_name,sep="_",into=c("zm","strand")) %>%
		separate(names,sep="_",into=c("start_queryspace","end_queryspace")) %>%
		rename(
			start_refspace = start,
			end_refspace = end
		) %>%
		mutate(
			zm=as.integer(zm),
			strand = factor(strand,levels=strand_levels),
			start_queryspace = as.integer(start_queryspace),
			end_queryspace = as.integer(end_queryspace),
			deletion_width = width
		) %>%
		mutate(deletion_width = deletion_width %>% na_if(0)) %>%
		left_join(
			data.frame(seqnames=seqnames(bam.gr[[i]]), zm=bam.gr[[i]]$zm) %>% distinct,
			by = join_by(zm)
		) %>%
		select(-group) %>%
	  as_tibble
	
	#Deletion base sequences from reference space, relative to reference plus strand
	indels_refspace$ref <- indels_refspace %>%
		mutate(strand="+") %>%
		makeGRangesFromDataFrame(
			seqnames.field="seqnames",
			start.field="start_refspace",
			end.field="end_refspace",
			strand.field="strand"
		) %>%
		getSeq(eval(parse(text=yaml.config$BSgenome$BSgenome_name)),.) %>%
		as.character
	
	#Insertion base sequences from query space, relative to reference plus strand
	insertions_queryspace_alt <- extractAt(bam.gr[[i]]$seq,cigar_query_indels) %>%
		as("CharacterList") %>%
		setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_")) %>%
		as.list %>%
		enframe %>%
		unnest_longer(col=value,indices_to="start_end_queryspace",values_to="alt") %>%
		separate(name,sep="_",into=c("zm","strand")) %>%
		separate(start_end_queryspace,sep="_",into=c("start_queryspace","end_queryspace")) %>%
		mutate(
			zm = as.integer(zm),
			strand = factor(strand,levels=strand_levels),
			start_queryspace = as.integer(start_queryspace),
			end_queryspace = as.integer(end_queryspace)
		)
	
	#Indel base qualities. For insertions, get base qualities of inserted bases in query space, relative to reference plus strand. For deletions, get base qualities of left and right flanking bases in query space, relative to reference plus strand.
	cigar_query_indels_forqualdata <- cigar_query_indels %>% unlist
	deletionidx <- start(cigar_query_indels_forqualdata) > end(cigar_query_indels_forqualdata)
	newstart <- end(cigar_query_indels_forqualdata)
	newend <- start(cigar_query_indels_forqualdata)
	start(cigar_query_indels_forqualdata[deletionidx]) <- newstart[deletionidx]
	end(cigar_query_indels_forqualdata[deletionidx]) <- newend[deletionidx]
	cigar_query_indels_forqualdata <- cigar_query_indels_forqualdata %>% relist(.,cigar_query_indels)
	
	rm(cigar_query_indels,deletionidx,newstart,newend)
	
	indels_qual <- extractAt(bam.gr[[i]]$qual,cigar_query_indels_forqualdata) %>%
		as("CharacterList") %>%
		setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_")) %>%
		as.list %>%
		enframe %>%
		unnest_longer(col=value,indices_to="start_end_queryspace",values_to="qual") %>%
		separate(name,sep="_",into=c("zm","strand")) %>%
		separate(start_end_queryspace,sep="_",into=c("start_queryspace","end_queryspace")) %>%
		mutate(
			zm = as.integer(zm),
			strand = factor(strand,levels=strand_levels),
			start_queryspace = as.integer(start_queryspace),
			end_queryspace = as.integer(end_queryspace),
			qual = qual %>% PhredQuality %>% as("IntegerList") %>% as.list %>% set_names(NULL)
		)
	
	#Extract sa/sm/sx tags
	
	 #Convert indel positions in query space to a format for faster retrieval below of sa, sm, sx tags.
	 # Flatten indel ranges, annotate with start_end_queryspace, then expand the ranges into the needed query_pos, preserving start_end_queryspace.
	indels_query_pos <- cigar_query_indels_forqualdata %>%
		setNames(str_c(bam.gr[[i]]$zm,strand(bam.gr[[i]]) %>% as.character,sep="_")) %>%
		as.data.frame %>%
		rename(
			zm_strand = group_name,
			start_queryspace = start,
			end_queryspace = end,
			start_end_queryspace = names
		) %>%
		select(zm_strand,start_queryspace,end_queryspace,start_end_queryspace) %>%
		group_by(zm_strand) %>%
		mutate(width=end_queryspace - start_queryspace + 1) %>%
		uncount(width) %>%
		group_by(zm_strand,start_end_queryspace) %>%
		mutate(pos_queryspace = start_queryspace + row_number() - 1) %>%
		select(zm_strand,start_end_queryspace,pos_queryspace) %>%
		as.data.table
	
	setkey(indels_query_pos, zm_strand, start_end_queryspace)
	
	 #Extract data
	indels_query_pos[, sa_val := sa[[zm_strand[1]]][pos_queryspace], by = zm_strand]
	indels_query_pos[, sm_val := sm[[zm_strand[1]]][pos_queryspace], by = zm_strand]
	indels_query_pos[, sx_val := sx[[zm_strand[1]]][pos_queryspace], by = zm_strand]
	
	rm(sa,sm,sx,cigar_query_indels_forqualdata)
	
	 #Aggregate back into a list-of-vectors per (name, start_end_queryspace), and reformat columns for later joining
	indels_query_pos <- indels_query_pos[, .(sa = list(sa_val), sm = list(sm_val), sx = list(sx_val)), by = .(zm_strand, start_end_queryspace)] %>%
	  as_tibble %>%
  	separate(zm_strand,sep="_",into=c("zm","strand")) %>%
	  separate(start_end_queryspace,sep="_",into=c("start_queryspace","end_queryspace")) %>%
  	  mutate(
  	    zm = as.integer(zm),
  	    strand = factor(strand,levels=strand_levels),
  	    start_queryspace = as.integer(start_queryspace),
  	    end_queryspace = as.integer(end_queryspace)
  	  )
	
	#Combine indel info
	indels.df[[i]] <- indels_queryspace %>%
		left_join(indels_refspace, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
		left_join(insertions_queryspace_alt, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
		left_join(indels_qual, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
		left_join(indels_query_pos, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
		mutate(
			indel_type = if_else(start_queryspace > end_queryspace, "deletion", "insertion") %>% factor,
			indel_width = if_else(is.na(insertion_width), deletion_width, insertion_width)
		) %>%
		select(zm,seqnames,strand,start_queryspace,end_queryspace,start_refspace,end_refspace,indel_type,ref,alt,indel_width,qual,sa,sm,sx)
	
	rm(indels_queryspace,indels_refspace,insertions_queryspace_alt,indels_qual,indels_query_pos)
	
	#Annotate each indel as a change from the reference in reference space in only one strand ("indel_mismatch-ss"), non-complementary changes from the reference in both strands ("indel_mismatch-ds"; only possible for insertions with different sequences), or complementary changes from the reference in both strands ("indel_mutation").
	indels.df[[i]] <- indels.df[[i]] %>%
		group_by(zm,start_refspace,end_refspace) %>%
		mutate(
			count = n(),
			count_distinct_alt = n_distinct(alt),
			variant_type = case_when(
				count == 1 ~ "indel_mismatch-ss",
				count == 2 & count_distinct_alt == 1 ~ "indel_mutation",
				count == 2 & count_distinct_alt == 2 ~ "indel_mismatch-ds"
			) %>%
				factor
		) %>%
		select(-count,-count_distinct_alt) %>%
		ungroup %>%
	  arrange(zm,start_refspace,end_refspace,strand)
}

cat("DONE\n")

} #End loop over each run

######################
### Save output
######################
cat("#### Saving output...")

#Collapse bam.gr back to a single GRanges object and remove seq and qual data that is no longer needed
bam.gr <- bam.gr %>% unlist
bam.gr$seq <- NULL
bam.gr$qual <- NULL

qs_save(
  list(
  	run_metadata = run_metadata,
  	bam.gr = bam.gr,
  	sbs.df = sbs.df %>% bind_rows(.id="run_id") %>% mutate(run_id=factor(run_id)),
  	indels.df = indels.df %>% bind_rows(.id="run_id") %>% mutate(run_id=factor(run_id)),
  	stats = stats
  	),
  bamFile %>% basename %>% str_replace(".bam$",".qs2")
  )

cat("DONE\n")