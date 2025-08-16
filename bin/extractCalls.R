#!/usr/bin/env -S Rscript --vanilla

#extractCalls.R:
# Perform initial molecule filtering and extract calls, including all call information required for analysis.  
# Note: extracts calls of all call_types across all chromosomes, not just those in configured chromgroups. filterCalls then removes calls outside configured chromgroups.

cat("#### Running extractCalls ####\n")

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(tidyverse))

######################
### Load configuration
######################
cat("## Loading configuration...")

#General options
options(datatable.showProgress = FALSE)
options(warn=2) #Stop script for any warnings

#Command line arguments
option_list = list(
  make_option(c("-c", "--config"), type = "character", default=NULL,
              help="path to YAML configuration file"),
  make_option(c("-b", "--bam"), type = "character", default=NULL,
              help="path to BAM file"),
  make_option(c("-o", "--output"), type = "character", default=NULL,
  						help="output qs2 file")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$config) | is.na(opt$bam) | is.na(opt$output)){
  stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))
bamFile <- opt$bam
outputFile <- opt$output

#Load the BSgenome reference
suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name,character.only=TRUE,lib.loc=yaml.config$cache_dir))

#General parameters
strand_levels <- c("+","-")

#Load call types, and filter to make a table only for MDB call types
call_types <- yaml.config$call_types %>%
  enframe(name=NULL) %>%
  unnest_wider(value) %>%
  unnest_longer(SBSindel_call_types) %>%
  unnest_wider(SBSindel_call_types)

call_types_MDB <- call_types %>%
  filter(call_class == "MDB")

#Load chromgroups (used for calculating molecule_stats per chromgroup) and add a chromgroup that combines all chromgroups
chroms_toanalyze <- yaml.config$chromgroups %>%
	enframe %>%
	unnest_wider(value) %>%
	mutate(
		chroms = chroms %>%
			map(function(x){
				x %>% str_split_1(",") %>% str_trim
			})
	) %>%
	select(-name, -sensitivity_thresholds) %>%
	bind_rows(
		tibble(
			chromgroup = "all_chroms",
			chroms = yaml.config$BSgenome$BSgenome_name %>% get %>% seqnames %>% list
		)
	)

cat("DONE\n")

######################
### Define custom functions
######################
#Function to count number of remaining molecules and molecule query space and reference space bases per run, and also per run x chromgroup
calculate_molecule_stats <- function(bam.gr.input, chroms_toanalyze.input, stat_label_suffix) {
	chroms_toanalyze.input %>%
		mutate(
			stats = map(chroms, ~{
				bam.gr.input %>%
					as_tibble %>%
					filter(seqnames %in% .x) %>%
						group_by(run_id) %>%
						summarize(
							num_molecules = n_distinct(zm),
							num_queryspacebases = sum(qwidth),
							num_refpacebases = sum(end - start),
							.groups = "drop"
						) %>%
					complete(run_id,fill = list(
						num_molecules=0, num_queryspacebases=0, num_refpacebases=0
						)
					)
				})
		) %>%
		unnest(stats) %>%
		pivot_longer(
			cols = starts_with("num_"),
			names_to = "stat",
			values_to = "value"
		) %>%
		mutate(
			stat = str_c(stat, ".", stat_label_suffix)
		) %>%
		select(-chroms) %>%
		arrange(stat) %>%
		relocate(run_id)
}

######################
### Load BAM file
######################
cat("## Loading BAM file:",bamFile,"...")

#Get run metadata from BAM file and join to run_id from config file
run_metadata <- bamFile %>%
  scanBamHeader %>%
  pluck(1,"text") %>%
  keep(names(.)=="@RG") %>%
  map(
    . %>%
      as_tibble %>%
      separate_wider_delim(value,delim=":",names=c("name","value"),too_many="merge") %>%
      pivot_wider
    ) %>%
  bind_rows %>%
  separate_rows(DS,sep=";") %>%
  separate_wider_delim(DS,delim="=",names=c("name","value"),too_many="merge") %>%
  pivot_wider %>%
  rename(
    movie_id = PU,
    rg_id = ID
    ) %>%
  mutate(
    run_id = movie_id %>%
      map_chr(
        function(x){
          yamlruns <- yaml.config$runs %>%
            map_df(~ tibble(run_id = .x$run_id, reads_file = .x$reads_file))
          
          yamlruns %>%
            pluck("run_id") %>%
            keep(yamlruns %>% pluck("reads_file") %>% str_detect(x))
        }
      )
  ) %>%
  type_convert %>%
  suppressMessages

#Load BAM file reads
bam <- bamFile %>%
  scanBam(
    param=ScanBamParam(
      what=setdiff(scanBamWhat(),c("mrnm","mpos","groupid","mate_status")),
      tag=c(
        "ec","np","rq","zm","sa","sm","sx","RG",
        call_types_MDB %>% pluck("call_bam_scoretag")
        )
    )
  ) %>%
  pluck(1)

#Create DataFrame of reads. DataFrame is required to store seq and qual at smaller size, since tibble cannot store S4Vector columns.
bam.df <- cbind(
  DataFrame(bam[names(bam) != "tag"] %>% lapply(I)),
  DataFrame(bam$tag %>% lapply(I))
)
rm(bam)
invisible(gc())

#Add movie_ids and run_ids using run metadata, while maintaining DataFrame format
bam.df <- bam.df %>%
	cbind(run_metadata[match(bam.df$RG, run_metadata$rg_id),c("movie_id","run_id")])

#Extract ccs strand
bam.df$ccs_strand <- bam.df$qname %>% str_extract("(fwd|rev)$")

#Convert some columns to factor
for(i in c("flag","RG","movie_id","run_id","ccs_strand")){
	bam.df[,i] <- factor(bam.df[,i])
}

#Recalculate isize, since pbmm2 has a bug that miscalculates it
bam.df$isize <- cigarWidthAlongReferenceSpace(bam.df$cigar)

#Reformat sa, sm, sx tags
 #sa: the number of subread alignments that span each CCS read position. sa is run-length encoded (rle) in the BAM file as length, value pairs. We inverse the rle encoding back to standard per-position values, then reverse the value orders for reads aligned to genome minus strand, then convert back to rle (as a list of rles instead of RleList as the former is much faster to convert to tibble).
sa.lengths <- lapply(bam.df$sa,function(x){x[c(TRUE, FALSE)]})
sa.values <- lapply(bam.df$sa,function(x){x[c(FALSE, TRUE)]})
bam.df$sa <- mapply(
		function(x,y){inverse.rle(list(lengths=x,values=y))},
		x=sa.lengths,
		y=sa.values,
		USE.NAMES=FALSE
	)

rm(sa.lengths,sa.values)
invisible(gc())

bam.df$sa[bam.df$strand %>% as.vector == "-"] <- lapply(bam.df$sa[bam.df$strand %>% as.vector == "-"], rev)

bam.df$sa <- bam.df$sa %>% lapply(Rle)

 #sm: the number of subreads that align as a match to each CCS read position. Reverse the value orders for reads aligned to genome minus strand
bam.df$sm[bam.df$strand %>% as.vector == "-"] <- lapply(bam.df$sm[bam.df$strand %>% as.vector == "-"], rev)

 #sx: the number of subreads that align as a mismatch to each CCS read position. Reverse the value orders for reads aligned to genome minus strand
bam.df$sx[bam.df$strand %>% as.vector == "-"] <- lapply(bam.df$sx[bam.df$strand %>% as.vector == "-"], rev)

#Convert to GRanges
bam.df$end <- bam.df$pos + bam.df$isize - 1

bam.gr <- bam.df %>%
  makeGRangesFromDataFrame(
    seqnames.field="rname",
    start.field="pos",
    end.field="end",
    strand.field="strand",
    keep.extra.columns=TRUE,
    seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
  )

rm(bam.df)
invisible(gc())

# Arrange bam.gr by run id, ZMW id and strand
bam.gr <- bam.gr[order(bam.gr$run_id,bam.gr$zm, strand(bam.gr))]

# Calculate molecule_stats
molecule_stats <- calculate_molecule_stats(
	bam.gr.input = bam.gr,
	chroms_toanalyze.input = chroms_toanalyze,
	stat_label_suffix = "prefilter"
)

cat("DONE\n")

######################
### Initial molecule filtering
######################
cat("## Initial molecule filtering...")

# Keep only molecules with 1 forward and 1 reverse ccs strand reads that are primary alignments to different strands on the same chromosome (flags = 0 and 16 are plus and minus strand primary alignments, respectively; flags =  2048 and 2064 are plus and minus strand supplementary alignments, respectively).
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
  distinct(run_id,zm)

bam.gr <- bam.gr[
	vctrs::vec_in(
		bam.gr %>% mcols %>% as_tibble %>% select(run_id,zm),
		zmwstokeep
		)
	]

# Calculate molecule_stats
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats(
			bam.gr.input = bam.gr,
			chroms_toanalyze.input = chroms_toanalyze,
			stat_label_suffix = "passalignmentfilter"
		)
)

# Keep only molecules with plus and minus strand alignment overlap >= min_strand_overlap (reciprocal or both plus and minus strand alignments)

 # Extract plus and minus strand reads separately. Due to prior sorting by ZMW id and strand, these are guaranteed to have the same ZMW id order.
bam.gr.onlyranges.plus <- bam.gr %>%
	filter(strand == "+") %>%
	select(run_id,zm)
bam.gr.onlyranges.minus <- bam.gr %>%
	filter(strand == "-") %>%
	select(run_id,zm)

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
	) %>%
	distinct(run_id,zm)

bam.gr <- bam.gr[
		vctrs::vec_in(
			bam.gr %>% mcols %>% as_tibble %>% select(run_id,zm),
			zmwstokeep
		)
	]

 # Remove intermediate objects
rm(bam.gr.onlyranges.plus, bam.gr.onlyranges.minus, bam.gr.onlyranges.overlap, zmwstokeep)
invisible(gc())

# Calculate molecule_stats
molecule_stats <- molecule_stats %>%
	bind_rows(
		calculate_molecule_stats(
			bam.gr.input = bam.gr,
			chroms_toanalyze.input = chroms_toanalyze,
			stat_label_suffix = "passminstrandoverlapfilter"
		)
	)

cat("DONE\n")

######################
### Define function to extract calls
######################
extract_calls <- function(bam.gr.input, call_class.input, call_type.input, cigar.ops.input = NULL, MDB_score_bamtag.input = NULL, MDB_score_min.input = NULL){
  #cigar.ops.input for the call_class: "X" for SBS; "I" for insertion; "D" for deletion; NULL for MDB
  #score_bamtag.input: tag in the BAM file containing the MDB score data (MDB only)
  #score_min.input: minimum score of MDBs to extract (MDB only)
  
  #Initialize empty calls tibble
  calls.out <- tibble(
    zm=integer(),
    seqnames=factor(),
    strand=factor(),
    start_queryspace=integer(),
    end_queryspace=integer(),
    start_refspace=integer(),
    end_refspace=integer(),
    ref_plus_strand=character(),
    alt_plus_strand=character(),
    qual=list(),
    qual.opposite_strand=list(),
    MDB_score=numeric(),
    sa=list(),
    sa.opposite_strand=list(),
    sm=list(),
    sm.opposite_strand=list(),
    sx=list(),
    sx.opposite_strand=list(),
    call_class=factor(),
    call_type=factor(),
    indel_width=integer()
    )
  
  #Extract call positions
  if(call_class.input %in% c("SBS","indel")){
    
    #Extract call positions in query space, and name calls with the position to retain this in the final data frame
    cigar.queryspace.var <- cigarRangesAlongQuerySpace(bam.gr.input$cigar,ops=cigar.ops.input)
    cigar.queryspace.var <- cigar.queryspace.var %>%
      unlist(use.names=FALSE) %>%
      setNames(str_c(start(.),end(.),sep="_")) %>%
      relist(.,cigar.queryspace.var)
    
    #Extract call positions in reference space
    cigar.refspace.var <- cigarRangesAlongReferenceSpace(bam.gr.input$cigar,ops=cigar.ops.input,pos=start(bam.gr.input))
    
  }else if(call_class.input == "MDB"){
    
  }
  
  #Check if no calls
  if(cigar.queryspace.var %>% as.data.frame %>% nrow == 0){
    cat("  No calls extracted of type",vartype_label,"!\n")
    return(calls.out)
  }
  
  #Call positions in query space
  # Note for indels: in query space, deletion width = 0.
  var_queryspace <- cigar.queryspace.var %>%
    setNames(str_c(bam.gr.input$zm,strand(bam.gr.input) %>% as.character,sep="_")) %>%
    as.data.frame %>%
    separate(group_name,sep="_",into=c("zm","strand")) %>%
    mutate(
      zm = as.integer(zm),
      strand = factor(strand,levels=strand_levels)
    ) %>%
    rename(
      start_queryspace = start,
      end_queryspace = end
      )
  
  if(call_class.input %in% c("SBS","MDB")){
    var_queryspace <- var_queryspace %>%
      uncount(weights=width) %>%
      select(zm,strand,start_queryspace) %>%
      group_by(zm,strand,start_queryspace) %>%
      mutate(
        start_queryspace = start_queryspace + row_number() - 1L,
        end_queryspace = start_queryspace
      ) %>%
      ungroup
  }else if(call_class.input == "indel"){
    var_queryspace <- var_queryspace %>%
      rename(insertion_width = width) %>%
      mutate(insertion_width = insertion_width %>% na_if(0)) %>%
      select(-group,-names) %>%
      as_tibble
  }
  
  #Call positions in reference space, annotated with positions in query space for later joining
  # Note for indels: in reference space, insertion width = 0.
  var_refspace <- cigar.refspace.var %>%
    unlist(use.names=FALSE) %>%
    setNames(
      str_c(
        cigar.queryspace.var %>% unlist(use.names=FALSE) %>% start,
        cigar.queryspace.var %>% unlist(use.names=FALSE) %>% end,
        sep="_"
        )
      ) %>%
    relist(.,cigar.refspace.var) %>%
    setNames(str_c(bam.gr.input$zm,strand(bam.gr.input) %>% as.character,sep="_")) %>%
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
      end_queryspace = as.integer(end_queryspace)
    ) %>%
    left_join(
      data.frame(seqnames=seqnames(bam.gr.input), zm=bam.gr.input$zm) %>% distinct,
      by = join_by(zm)
    )
  
  if(call_class.input %in% c("SBS","MDB")){
    var_refspace <- var_refspace %>%
      uncount(weights=width) %>%
      select(zm,strand,start_refspace,start_queryspace,seqnames) %>%
      group_by(zm,strand,start_refspace,start_queryspace,seqnames) %>%
      mutate(
        start_refspace = start_refspace + row_number() - 1L,
        end_refspace = start_refspace,
        start_queryspace = start_queryspace + row_number() - 1L,
        end_queryspace = start_queryspace
      ) %>%
      ungroup
  }else if(call_class.input == "indel"){
    var_refspace <- var_refspace %>%
      rename(deletion_width = width) %>%
      mutate(deletion_width = deletion_width %>% na_if(0)) %>%
      select(-group) %>%
      as_tibble
  }
  
  rm(cigar.refspace.var)
  invisible(gc())
  
  #Call REF base sequences, relative to reference plus strand
  var_refspace$ref_plus_strand <- var_refspace %>%
    mutate(strand="+") %>%
    makeGRangesFromDataFrame(
      seqnames.field="seqnames",
      start.field="start_refspace",
      end.field="end_refspace",
      strand.field="strand",
      seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
    ) %>%
    getSeq(eval(parse(text=yaml.config$BSgenome$BSgenome_name)),.) %>%
    as.character
  
  #Call ALT base sequences, relative to reference plus strand
  var_alt <- extractAt(bam.gr.input$seq,cigar.queryspace.var) %>%
    as("CharacterList") %>%
    setNames(str_c(bam.gr.input$zm,strand(bam.gr.input) %>% as.character,sep="_")) %>%
    as.list %>%
    enframe %>%
    unnest_longer(col=value,indices_to="start_end_queryspace",values_to="alt_plus_strand") %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    separate(start_end_queryspace,sep="_",into=c("start_queryspace","end_queryspace")) %>%
    mutate(
      zm = as.integer(zm),
      strand = factor(strand,levels=strand_levels),
      start_queryspace = as.integer(start_queryspace),
      end_queryspace = as.integer(end_queryspace)
    )
  
  if(call_class.input %in% c("SBS","MDB")){
    var_alt <- var_alt %>%
      separate_longer_position(alt_plus_strand,width=1) %>%
      group_by(zm,strand,start_queryspace) %>%
      mutate(
        start_queryspace = start_queryspace + row_number() - 1L,
        end_queryspace = start_queryspace
      ) %>%
      ungroup
  }
  
  #Call base qualities
  if(call_class.input %in% c("SBS","MDB")){
    cigar.queryspace.var.forqualdata <- cigar.queryspace.var
  }else if(call_class.input == "indel"){
    #Indel base qualities. For insertions, get base qualities of inserted bases in query space, relative to reference plus strand. For deletions, get base qualities of left and right flanking bases in query space, relative to reference plus strand.
    cigar.queryspace.var.forqualdata <- cigar.queryspace.var %>% unlist
    deletionidx <- start(cigar.queryspace.var.forqualdata) > end(cigar.queryspace.var.forqualdata)
    newstart <- end(cigar.queryspace.var.forqualdata)
    newend <- start(cigar.queryspace.var.forqualdata)
    start(cigar.queryspace.var.forqualdata[deletionidx]) <- newstart[deletionidx]
    end(cigar.queryspace.var.forqualdata[deletionidx]) <- newend[deletionidx]
    cigar.queryspace.var.forqualdata <- cigar.queryspace.var.forqualdata %>% relist(.,cigar.queryspace.var)
    
    rm(deletionidx,newstart,newend)
    invisible(gc())
  }
  
  var_qual <- extractAt(bam.gr.input$qual,cigar.queryspace.var.forqualdata) %>%
    as("CharacterList") %>%
    setNames(str_c(bam.gr.input$zm,strand(bam.gr.input) %>% as.character,sep="_")) %>%
    as.list %>%
    enframe %>%
    unnest_longer(col=value,indices_to="start_end_queryspace",values_to="qual") %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    separate(start_end_queryspace,sep="_",into=c("start_queryspace","end_queryspace")) %>%
    mutate(
      zm = as.integer(zm),
      strand = factor(strand,levels=strand_levels),
      start_queryspace = as.integer(start_queryspace),
      end_queryspace = as.integer(end_queryspace)
    )
  
  if(call_class.input %in% c("SBS","MDB")){
    var_qual <- var_qual %>%
      separate_longer_position(qual,width=1) %>%
      group_by(zm,strand,start_queryspace) %>%
      mutate(
        start_queryspace = start_queryspace + row_number() - 1L,
        end_queryspace = start_queryspace
      ) %>%
      ungroup
  }
  
  var_qual <- var_qual %>%
    mutate(
      qual = qual %>%
      	PhredQuality %>%
      	as("IntegerList") %>%
      	as.list %>%
      	set_names(NULL)
    )
  
  rm(cigar.queryspace.var)
  invisible(gc())
  
  #Call base qualities in opposite strand (matched via reference space).
  #For insertions, get the base qualities on the opposite strand from the left and right flanking bases in opposite strand query space that correspond to the call in reference space. For deletions, get base qualities from all bases in opposite strand query space that correspond to the call in reference space, or left and right flanking bases if opposite strand query space width = 0. For mutations, this temporarily assigns wrong values since for these the correct opposite strand bases to get quality from would be the inserted bases on the opposite strand for insertions and the bases flanking the deletion location for deletions. The issue is that we're going from call strand query space to reference space to opposite strand queryspace, rather than through direct alignment of the two strands to each other, which is too computationally complicated. For mutations, these imperfect values are later corrected to the precise opposite strand data: when we annotate which calls are mutations (i.e., on both strands), we reassign opposite strand data for each call from the call call on the opposite strand. For mismatches, this is not possible, so we use the above heuristics for filtering.
  #Set missing values to NA (for example, an SBS across from a deletion).
   #Make GRanges of call locations in reference space.
  var_refspace.gr <- var_refspace %>%
    select(-ref_plus_strand) %>%
    makeGRangesFromDataFrame(
      seqnames.field="seqnames",
      start.field="start_refspace",
      end.field="end_refspace",
      strand.field="strand",
      keep.extra.columns=TRUE,
      seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
    )
  
   #Make GRanges of bam reads parallel to call refspace GRanges, but matching to the opposite strand read from the call
  bam.gr.input.opposite_strand.parallel_to_var_refspace <- bam.gr.input %>%
    select(zm, cigar, qual) %>%
    mutate(
      opposite_strand = if_else(strand %>% as.character == "+", "-", "+"),
      strand = "*"
    ) %>%
    setNames(str_c(.$zm, .$opposite_strand, sep="_")) %>%
    select(-opposite_strand) %>%
    slice(
      match(
        str_c(var_refspace.gr$zm,strand(var_refspace.gr) %>% as.character,sep="_"),
        names(.)
      )
    )
  
   #Convert calls from reference space to query space of opposite strand
  vars_queryspace.opposite_strand <- var_refspace.gr %>%
    pmapToAlignments(bam.gr.input.opposite_strand.parallel_to_var_refspace %>% as("GAlignments")) %>%
    setNames(str_c(var_refspace.gr$start_queryspace, var_refspace.gr$end_queryspace,sep="_")) %>%
    (function(x) IRanges(start=start(x), end=end(x), names=names(x), seqnames = seqnames(x) %>% as.character)) %>%
    filter(seqnames != "UNMAPPED")
  
  #For indels, if opposite strand query space width = 0, swap start/end to get flanking bases.
  if(call_class.input == "indel"){
    vars_queryspace.opposite_strand <- vars_queryspace.opposite_strand %>%
      mutate(
        start_tmp = start,
        end_tmp = end,
        start = if_else(start_tmp > end_tmp, end, start),
        end = if_else(start_tmp > end_tmp, start, end)
      ) %>%
      select(-start_tmp,-end_tmp)
  }
  
  vars_queryspace.opposite_strand <- vars_queryspace.opposite_strand %>%
    split(mcols(.)$seqnames)
  
   #Extract qual from opposite strand. Set missing values to NA.
  var_qual.opposite_strand <- extractAt(
    bam.gr.input.opposite_strand.parallel_to_var_refspace %>%
        slice(
          match(
            names(vars_queryspace.opposite_strand),
            names(.)
            )
          ) %>%
        .$qual,
      vars_queryspace.opposite_strand %>% setNames(NULL)
    ) %>%
    as("CharacterList") %>%
    setNames(names(vars_queryspace.opposite_strand)) %>%
    as.list %>%
    enframe %>%
    unnest_longer(col=value,indices_to="start_end_queryspace",values_to="qual.opposite_strand") %>%
    separate(name,sep="_",into=c("zm","strand")) %>%
    separate(start_end_queryspace,sep="_",into=c("start_queryspace","end_queryspace")) %>%
    mutate(
      zm = as.integer(zm),
      strand = factor(strand,levels=strand_levels),
      start_queryspace = as.integer(start_queryspace),
      end_queryspace = as.integer(end_queryspace),
      qual.opposite_strand = qual.opposite_strand %>%
        PhredQuality %>%
        as("IntegerList") %>%
        as.list %>%
        set_names(NULL) %>%
        map_if(~ length(.x) == 0, ~NA_integer_)
    )
  
  rm(var_refspace.gr, bam.gr.input.opposite_strand.parallel_to_var_refspace)
  invisible(gc())
  
  
  #Extract sa/sm/sx tags
  # For insertions, get base qualities of inserted bases in query space, relative to reference plus strand. For deletions, get base qualities of left and right flanking bases in query space, relative to reference plus strand.
  #For insertions, get the base qualities on the opposite strand from the left and right flanking bases in opposite strand query space that correspond to the call in reference space. For deletions, get base qualities from all bases in opposite strand query space that correspond to the call in reference space. For mutations, this temporarily assigns wrong values since for these the correct opposite strand bases to get quality from would be the inserted bases on the opposite strand for insertions and the bases flanking the deletion location for deletions. The issue is that we're going from call strand query space to reference space to opposite strand queryspace, rather than through direct alignment of the two strands to each other, which is too computationally complicated. For mutations, these imperfect values are later corrected to the precise opposite strand data: when we annotate which calls are mutations (i.e., on both strands), we reassign opposite strand data for each call from the call on the opposite strand. For mismatches, this is not possible, so we use the above heuristics for filtering.
  
  #Format sa, sm, sx inputs
  sa.input <- bam.gr.input$sa %>%
  	lapply(as.vector) %>%
  	setNames(str_c(bam.gr.input$zm, strand(bam.gr.input) %>% as.character,sep="_"))
  
  sm.input <- bam.gr.input$sm %>%
  	setNames(str_c(bam.gr.input$zm, strand(bam.gr.input) %>% as.character,sep="_"))
  
  sx.input <- bam.gr.input$sx %>%
  	setNames(str_c(bam.gr.input$zm, strand(bam.gr.input) %>% as.character,sep="_"))
  
  if(call_class.input %in% c("SBS","MDB")){
    #Convert call positions in query space to list for faster retrieval below of sa, sm, sx tags
    var_queryspace.list <- var_queryspace %>%
      mutate(zm_strand = str_c(zm,strand,sep="_")) %>%
      select(zm_strand,start_queryspace) %>%
      as.data.table %>%
      split(by="zm_strand",keep.by=F) %>%
      lapply(unlist,use.names=F)
    
    var_queryspace.opposite_strand.list <- vars_queryspace.opposite_strand %>%
      as.data.table %>%
      filter(start==end) %>% #Remove sites whose opposite strand aligned to width = -1, meaning there is no reference-aligned base on the opposite strand. These will then become NA after the final join.
      select(seqnames, start) %>%
      split(by="seqnames",keep.by=F) %>%
      lapply(unlist,use.names=F)
    
    vars_queryspace.coordconversion <- vars_queryspace.opposite_strand %>%
      as.data.table %>%
      as_tibble %>%
      filter(start==end) %>%
      separate(seqnames,sep="_",into=c("zm","strand")) %>%
      separate(names,sep="_",into=c("start_queryspace","end_queryspace")) %>%
      mutate(
        zm = as.integer(zm),
        strand = factor(strand,levels=strand_levels),
        start_queryspace = as.integer(start_queryspace),
        end_queryspace = start_queryspace
      ) %>%
      select(-group,-group_name,-width)
    
    #Define extraction function
    extract_sasmsx <- function(tag.input, tagdata.input, var_queryspace.list.input, vars_queryspace.coordconversion.input=NULL){
      result <- map2(tagdata.input[var_queryspace.list.input %>% names], var_queryspace.list.input, ~ .x[.y] %>% set_names(.y)) %>%
        enframe %>%
        separate(name,sep="_",into=c("zm","strand")) %>%
        unnest_longer(col=value, values_to=tag.input, indices_to="start_queryspace") %>%
          mutate(
            zm = as.integer(zm),
            strand = factor(strand,levels=strand_levels),
            start_queryspace = as.integer(start_queryspace),
            end_queryspace = start_queryspace,
            !!tag.input := .data[[tag.input]] %>% as.list %>% set_names(NULL)
          )
      
      #Convert query space coordinate of opposite strand to call strand (if extracting opposite strand data)
      if(!is.null(vars_queryspace.coordconversion.input)){
        result <- result  %>%
          rename(
            start = start_queryspace,
            end = end_queryspace
          ) %>%
          left_join(
            vars_queryspace.coordconversion.input,
            by = join_by(zm,strand,start,end)
          ) %>%
          select(-start,-end)
      }
      
      return(result)
    }
    
    #Extract data
    var_sa <- extract_sasmsx(
      "sa",
      sa.input,
      var_queryspace.list
      )
    
    var_sm <- extract_sasmsx(
      "sm",
      sm.input,
      var_queryspace.list
    )
 
    var_sx <- extract_sasmsx(
      "sx",
      sx.input,
      var_queryspace.list
    )

    var_sa.opposite_strand <- extract_sasmsx(
      "sa.opposite_strand",
      sa.input %>% setNames(names(.) %>% chartr("+-","-+",.)),
      var_queryspace.opposite_strand.list,
      vars_queryspace.coordconversion
    )
    
    var_sm.opposite_strand <- extract_sasmsx(
      "sm.opposite_strand",
      sm.input %>% setNames(names(.) %>% chartr("+-","-+",.)),
      var_queryspace.opposite_strand.list,
      vars_queryspace.coordconversion
    )
    
    var_sx.opposite_strand <- extract_sasmsx(
      "sx.opposite_strand",
      sx.input %>% setNames(names(.) %>% chartr("+-","-+",.)),
      var_queryspace.opposite_strand.list,
      vars_queryspace.coordconversion
    )
    
    rm(var_queryspace.list, var_queryspace.opposite_strand.list, vars_queryspace.coordconversion)
    invisible(gc())
    
  }else if(call_class.input == "indel"){
    #Convert indel positions in query space to a format for faster retrieval below of sa, sm, sx tags.
    # Flatten indel ranges, annotate with start_end_queryspace, then expand the ranges into the needed query_pos, preserving start_end_queryspace.
    indels_queryspace_pos <- cigar.queryspace.var.forqualdata %>%
      setNames(str_c(bam.gr.input$zm,strand(bam.gr.input) %>% as.character,sep="_")) %>%
      as.data.frame %>%
      rename(
        zm_strand = group_name,
        start_queryspace = start,
        end_queryspace = end,
        start_end_queryspace = names
      ) %>%
      select(zm_strand,start_queryspace,end_queryspace,start_end_queryspace,width) %>%
      group_by(zm_strand) %>%
      uncount(width) %>%
      group_by(zm_strand,start_end_queryspace) %>%
      mutate(pos_queryspace = start_queryspace + row_number() - 1) %>%
      select(zm_strand,start_end_queryspace,pos_queryspace) %>%
      as.data.table
    
    indels_queryspace_pos.opposite_strand <- vars_queryspace.opposite_strand %>%
      as.data.frame %>%
      select(-seqnames) %>%
      rename(
        zm_strand = group_name,
        start_queryspace = start,
        end_queryspace = end,
        start_end_queryspace = names
      ) %>%
      select(zm_strand,start_queryspace,end_queryspace,start_end_queryspace,width) %>%
      group_by(zm_strand) %>%
      uncount(width) %>%
      group_by(zm_strand,start_end_queryspace) %>%
      mutate(pos_queryspace = start_queryspace + row_number() - 1) %>%
      select(zm_strand,start_end_queryspace,pos_queryspace) %>%
      as.data.table
    
    setkey(indels_queryspace_pos, zm_strand, start_end_queryspace)
    setkey(indels_queryspace_pos.opposite_strand, zm_strand, start_end_queryspace)
    
    #Extract data
    indels_queryspace_pos[, sa_val := sa.input[[zm_strand[1]]][pos_queryspace], by = zm_strand]
    indels_queryspace_pos[, sm_val := sm.input[[zm_strand[1]]][pos_queryspace], by = zm_strand]
    indels_queryspace_pos[, sx_val := sx.input[[zm_strand[1]]][pos_queryspace], by = zm_strand]
    
    sa.input.opposite_strand <- sa.input %>% setNames(names(.) %>% chartr("+-","-+",.))
    sm.input.opposite_strand <- sm.input %>% setNames(names(.) %>% chartr("+-","-+",.))
    sx.input.opposite_strand <- sx.input %>% setNames(names(.) %>% chartr("+-","-+",.))

    indels_queryspace_pos.opposite_strand[, sa_val := sa.input.opposite_strand[[zm_strand[1]]][pos_queryspace], by = zm_strand]
    indels_queryspace_pos.opposite_strand[, sm_val := sm.input.opposite_strand[[zm_strand[1]]][pos_queryspace], by = zm_strand]
    indels_queryspace_pos.opposite_strand[, sx_val := sx.input.opposite_strand[[zm_strand[1]]][pos_queryspace], by = zm_strand]
    
    rm(cigar.queryspace.var.forqualdata, vars_queryspace.opposite_strand)
    invisible(gc())
    
    #Aggregate back into a list-of-vectors per (name, start_end_queryspace), and reformat columns for later joining
    indels_queryspace_pos <- indels_queryspace_pos[, .(sa = list(sa_val), sm = list(sm_val), sx = list(sx_val)), by = .(zm_strand, start_end_queryspace)] %>%
      as_tibble %>%
      separate(zm_strand,sep="_",into=c("zm","strand")) %>%
      separate(start_end_queryspace,sep="_",into=c("start_queryspace","end_queryspace")) %>%
      mutate(
        zm = as.integer(zm),
        strand = factor(strand,levels=strand_levels),
        start_queryspace = as.integer(start_queryspace),
        end_queryspace = as.integer(end_queryspace)
      )
    
    indels_queryspace_pos.opposite_strand <- indels_queryspace_pos.opposite_strand[, .(sa.opposite_strand = list(sa_val), sm.opposite_strand = list(sm_val), sx.opposite_strand = list(sx_val)), by = .(zm_strand, start_end_queryspace)] %>%
      as_tibble %>%
      separate(zm_strand,sep="_",into=c("zm","strand")) %>%
      separate(start_end_queryspace,sep="_",into=c("start_queryspace","end_queryspace")) %>%
      mutate(
        zm = as.integer(zm),
        strand = factor(strand,levels=strand_levels),
        start_queryspace = as.integer(start_queryspace),
        end_queryspace = as.integer(end_queryspace)
      )
    
    rm(sa.input.opposite_strand, sm.input.opposite_strand, sx.input.opposite_strand)
    invisible(gc())
  }
  
  rm(sa.input, sm.input, sx.input)
  invisible(gc())
  
  #Combine call info with prior empty initialized tibble and return result
  if(call_class.input %in% c("SBS","MDB")){
    calls.out <- calls.out %>%
      bind_rows(
        var_queryspace %>%
          left_join(var_refspace, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_alt, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_qual, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_qual.opposite_strand, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_sa, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_sa.opposite_strand, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_sm, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_sm.opposite_strand, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_sx, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_sx.opposite_strand, by=join_by(zm,strand,start_queryspace,end_queryspace))
      )
  }else if(call_class.input == "indel"){
    calls.out <- calls.out %>%
      bind_rows(
        var_queryspace %>%
          left_join(var_refspace, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_alt, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_qual, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_qual.opposite_strand, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(indels_queryspace_pos, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(indels_queryspace_pos.opposite_strand, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          mutate(
            indel_width = if_else(is.na(insertion_width), deletion_width, insertion_width)
          ) %>%
          select(-deletion_width,-insertion_width)
      )
  }
  
  #Replace NA/NULL opposite strand values with list(NA)
  calls.out <- calls.out %>%
    mutate(
      across(
        c(qual.opposite_strand, sa.opposite_strand, sm.opposite_strand, sx.opposite_strand),
        ~ replace(.x, vapply(.x, function(el) is.null(el) || identical(el, NA), logical(1)), list(NA))
      )
    )
  
  #Set call_class and call_type
  calls.out <- calls.out %>%
    mutate(
      call_class = !!call_class.input %>% factor,
      call_type = !!call_type.input %>% factor
      )
  
  return(calls.out)
}


#########################################
### Loop to extract calls from each run_id separately
#########################################

#Split bam.gr by run_id
bam.gr <- bam.gr %>% split(bam.gr$run_id)

#Create output list
calls <- list()

for(i in bam.gr %>% names){

  cat("> Processing run:",i,"\n")
  
  ######################
  ### Extract calls
  ######################
  cat("    Extracting single base substitutions...")
  
  #SBS
  calls[[i]] <- extract_calls(
    bam.gr.input = bam.gr[[i]],
    call_class.input = "SBS",
    call_type.input = "SBS",
    cigar.ops.input = "X"
    )
  
  cat("DONE\n")
  
  #Insertions
  cat("    Extracting insertions...")
  
  calls[[i]] <- calls[[i]] %>%
    bind_rows(
      extract_calls(
        bam.gr.input = bam.gr[[i]],
        call_class.input = "indel",
        call_type.input = "insertion",
        cigar.ops.input = c("I")
      )
    )

  cat("DONE\n")
  
  #Deletions
  cat("    Extracting deletions...")
  
  calls[[i]] <- calls[[i]] %>%
    bind_rows(
      extract_calls(
        bam.gr.input = bam.gr[[i]],
        call_class.input = "indel",
        call_type.input = "deletion",
        cigar.ops.input = c("D")
      )
    )
  
  cat("DONE\n")
  
  
  #MDBs
  for(j in seq_len(nrow(call_types_MDB))){
    
    cat("    Extracting",call_types_MDB %>% pluck("call_type",j),"...")
    
    calls[[i]] <- calls[[i]] %>%
      bind_rows(
        extract_calls(
          bam.gr.input = bam.gr[[i]],
          call_class.input = "MDB",
          call_type.input = call_types_MDB %>% pluck("call_type",j),
          call_bam_scoretag.input = call_types_MDB %>% pluck("call_bam_scoretag",j),
          call_bam_scoretag_min.input = call_types_MDB %>% pluck("call_bam_scoretag_min",j)
        )
      )
    
    cat("DONE\n")
  }
  
}

######################
### Further annotate calls
######################
cat("## Further annotating calls...")

#Collapse bam.gr back to a single GRanges object, remove seq data that is no longer needed, and convert to a tibble. Note, this transforms qual from PhredQuality to list.
bam.gr <- bam.gr %>% unlist(use.names=F)
bam.gr$seq <- NULL
bam <- bam.gr %>%
	as_tibble %>%
	mutate(strand = strand %>% factor(levels=strand_levels))

rm(bam.gr)
invisible(gc())

#Combine calls of all runs
calls <- calls %>%
  bind_rows(.id="run_id") %>%
  mutate(
  	run_id=factor(run_id),
  	strand = strand %>% factor(levels=strand_levels)
  	)

#Annotate for each calls its trinucleotide context when call_class = SBS or MDB. For indels, this is set to NA.
calls <- calls %>%
	left_join(
		calls %>%
			filter(call_class %in% c("SBS","MDB")) %>%
			select(seqnames,start_refspace,end_refspace,call_class) %>%
			distinct %>%

			#Make GRanges of trinculeotide positions
			makeGRangesFromDataFrame(
				seqnames.field="seqnames",
				start.field="start_refspace",
				end.field="end_refspace",
				keep.extra.columns=TRUE,
				seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
			) %>%
			resize(width=3,fix="center") %>%
			
			#Remove trinucleotide contexts that extend past a non-circular chromosome edge (after trim, width < 3) to yield NA tnc value
			trim %>%
			filter(width == 3) %>%
			
			#Get trinucleotide context sequence for the reference plus and minus strands. The latter is calculated here to avoid error when calculating it for NA values, and it is used for later assignment to synthesized and template strand trinucleotide columns.
			mutate(
				reftnc_plus_strand = getSeq(eval(parse(text=yaml.config$BSgenome$BSgenome_name)),.) %>%
				  as.character,
				reftnc_minus_strand = reftnc_plus_strand %>% DNAStringSet %>% reverseComplement %>% as.character
			) %>%
			
			#Resize back to original start/end, and make tibble
			resize(width=1,fix="center") %>%
			as_tibble %>%
			rename(
				start_refspace = start,
				end_refspace = end
			) %>%
			select(-width,-strand)
		,
		by = join_by(seqnames,start_refspace,end_refspace,call_class)
	)

#Annotate for each calls its ref, alt, and trinucleotide context sequences relative to the reference plus strand, synthesized strand, and template strand.
calls <- calls %>%
	mutate(
		ref_minus_strand = ref_plus_strand %>% DNAStringSet %>% reverseComplement %>% as.character,
		ref_synthesized_strand = if_else(strand=="+",ref_plus_strand,ref_minus_strand),
		ref_template_strand = if_else(strand=="+",ref_minus_strand,ref_plus_strand),
		
		alt_minus_strand = alt_plus_strand %>% DNAStringSet %>% reverseComplement %>% as.character,
		alt_synthesized_strand = if_else(strand=="+",alt_plus_strand,alt_minus_strand),
		alt_template_strand = if_else(strand=="+",alt_minus_strand,alt_plus_strand),
		
		reftnc_synthesized_strand = if_else(strand=="+",reftnc_plus_strand,reftnc_minus_strand),
		reftnc_template_strand = if_else(strand=="+",reftnc_minus_strand,reftnc_plus_strand)
	) %>%
	select(-c(ref_minus_strand,alt_minus_strand,reftnc_minus_strand))

#Annotate for each call (for all call classes: SBS, indel, MDB) the base sequence on the opposite strand based on reference space coordinates, taking into account if there was an SBS or an indel called on the opposite strand, regardless if there was an MDB called on the opposite strand.
# Possible overlaps are SBS/MDB with SBS, SBS/MDB with deletion, insertion with insertion, insertion with deletion, deletion with SBS, deletion with insertion, deletion with deletion.
# If a deletion on the opposite strand only partially overlaps an analyzed deletion, annotate it as a deletion on the opposite strand.
# For deletions called on both strands, annotate if the deletions have the same coordinates on both strands (deletion.bothstrands.startendmatch).
# Note, for an analyzed deletion, it is possible for > 1 SBS, insertion or deletion to overlap it on the opposite strand.
# Whether and which MDBs were called on the same or opposite strand will be annotated later in the filterCalls step after further filtering MDBs, since filtering out an MDB won't trigger filtering of SBS and indel calls, whereas filtering an SBS or indel on any one strand will anyway cause filtering of calls on the opposite strand, so it doesn't matter if we annotate here opposite strand SBS and indel calls before filtering.

 #Make GRanges of calls. Swap start and end for insertions, since the subsequent join otherwise cannot join two insertions with the same coordinates (this won't affect later joining, since original start_refspace and end_refspace are used to filter the joins).
calls.gr <- calls %>%
  mutate(
    start = if_else(call_type=="insertion",end_refspace,start_refspace),
    end = if_else(call_type=="insertion",start_refspace,end_refspace),
    strand_copy = strand
    ) %>%
  makeGRangesFromDataFrame(
    seqnames.field="seqnames",
    start.field="start",
    end.field="end",
    strand.field="strand",
    keep.extra.columns=TRUE,
    seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
  )

 #Annotate calls GRanges with opposite strand calls
calls.gr <- calls.gr %>%
  join_overlap_left(
    calls.gr %>% filter(call_class %in% c("SBS","indel")),
    suffix = c("",".opposite_strand")
  ) %>%
  filter(
    run_id == run_id.opposite_strand,
    zm == zm.opposite_strand,
    strand_copy != strand_copy.opposite_strand,
    (
      call_class %in% c("SBS","MDB") & call_type.opposite_strand == "SBS" &
      start_refspace == start_refspace.opposite_strand &
      end_refspace == end_refspace.opposite_strand
    ) |
    (
      call_class %in% c("SBS","MDB") & call_type.opposite_strand == "deletion" &
        start_refspace >= start_refspace.opposite_strand &
        end_refspace <= end_refspace.opposite_strand
    ) |
    (
      call_type == "insertion" & call_type.opposite_strand == "insertion" &
        start_refspace == start_refspace.opposite_strand &
        end_refspace == end_refspace.opposite_strand
    ) |
    (
      call_type == "insertion" & call_type.opposite_strand == "deletion" &
        start_refspace <= end_refspace.opposite_strand &
        end_refspace >= start_refspace.opposite_strand
    ) |
    (
      call_type == "deletion" & call_type.opposite_strand == "SBS" &
      start_refspace <= start_refspace.opposite_strand &
      end_refspace >= end_refspace.opposite_strand
    ) |
    (
      call_type == "deletion" & call_type.opposite_strand == "insertion" &
        start_refspace <= end_refspace.opposite_strand &
        end_refspace >= start_refspace.opposite_strand
    ) |
    (
      call_type == "deletion" & call_type.opposite_strand == "deletion" &
        start_refspace <= end_refspace.opposite_strand &
        end_refspace >= start_refspace.opposite_strand
    )
  ) %>%
  mutate(
    deletion.bothstrands.startendmatch = case_when(
      call_type != "deletion" | call_type.opposite_strand != "deletion" ~ NA,
    	start_refspace == start_refspace.opposite_strand & end_refspace == end_refspace.opposite_strand ~ TRUE,
      start_refspace != start_refspace.opposite_strand | end_refspace != end_refspace.opposite_strand ~ FALSE
    )
  )

 #Join to calls with opposite strand information.
 # Note that for deletions with one partially overlapping deletion on the opposite strand, alt_plus_strand.opposite_strand will still equal "" even though part of the analyzed deletion's sequence is still present in the opposite strand. Likewise for deletions with > 1 overlapping SBS, insertion, or deletion on the opposite strand, alt_plus_strand.opposite_strand will contain the sequence of one random one of these. These two issues would be difficult to fix, and are not critical for downstream analysis.
calls <- calls %>%
	left_join(
		calls.gr %>%
			as_tibble %>%
			select(
				run_id,zm,start_refspace,end_refspace,strand,
				call_class.opposite_strand,
				call_type.opposite_strand,
				alt_plus_strand.opposite_strand,
				deletion.bothstrands.startendmatch
			),
		by=join_by(run_id,zm,start_refspace,end_refspace,strand),
		multiple="any" # Required due to possibility of > 1 call on the opposite strand overlapping an analyzed deletion
	)

rm(calls.gr)
invisible(gc())

#Annotate for each call (for all call classes: SBS, indel, MDB) its 'SBSindel_call_type'.
  #  Call strand | Opposite strand (SBS or indel) overlapping based on reference space coordinates
  #       N      |       N                       			=> "match" (possible only for MDB)
  #       N      |       Y                       			=> "mismatch-os" (possible only for MDB)
  #       Y      |       N                     				=> "mismatch-ss"
  #       Y      |       Y (non-complementary change)  => "mismatch-ds"
  #       Y      |       Y (complementary change)      => "mutation"

calls <- calls %>%
	mutate(SBSindel_call_type = case_when(
			ref_plus_strand == alt_plus_strand & is.na(alt_plus_strand.opposite_strand) ~ "match",
			ref_plus_strand == alt_plus_strand & !is.na(alt_plus_strand.opposite_strand) ~ "mismatch-os",
			ref_plus_strand != alt_plus_strand & is.na(alt_plus_strand.opposite_strand) ~ "mismatch-ss",
			ref_plus_strand != alt_plus_strand & !is.na(alt_plus_strand.opposite_strand) &
				(
				  alt_plus_strand != alt_plus_strand.opposite_strand |
				  (
				    alt_plus_strand == alt_plus_strand.opposite_strand &
				      !is.na(deletion.bothstrands.startendmatch) &
				      deletion.bothstrands.startendmatch == FALSE
				  )
				) ~ "mismatch-ds",
			ref_plus_strand != alt_plus_strand & !is.na(alt_plus_strand.opposite_strand) &
				(
				  alt_plus_strand == alt_plus_strand.opposite_strand &
				  (
				    is.na(deletion.bothstrands.startendmatch) |
				      deletion.bothstrands.startendmatch == TRUE
				  )
				) ~ "mutation"
		) %>%
		  factor
	) %>%
  arrange(run_id,zm,start_refspace,end_refspace,call_type,strand,SBSindel_call_type)

#For mutations, reannotate each strand's qual.opposite_strand, sa.opposite_strand, sm.opposite_strand, and sx.opposite_strand based on the opposite strand's qual, sa, sm, sx. This is more accurate than the previously calculated opposite strand data, because the prior data utilized a query space -> reference space -> opposite strand query space transformation, whereas annotated mutations reflect corresponding deletion coordinates and insertion sequences that are lost in the prior transformation. Though note that mismatches still have imperfect opposite strand qual, sa, sm, sx, and to obtain this would require aligning strands to each other which would be significantly more complex to incorporate.

calls <- calls %>% left_join(
  #Extract only mutations, and reverse strand orientation so that opposite strand info is annotated
  calls %>%
    filter(SBSindel_call_type == "mutation") %>%
    select(
      run_id,zm,seqnames,strand,start_refspace,end_refspace,ref_plus_strand,alt_plus_strand,
      qual,sa,sm,sx
    ) %>%
    rename(
      qual.tmp = qual,
      sa.tmp = sa,
      sm.tmp = sm,
      sx.tmp = sx,
    ) %>%
    mutate(strand = if_else(strand %>% as.character == "+", "-", "+")),
  by = join_by(run_id,zm,seqnames,strand,start_refspace,end_refspace,ref_plus_strand,alt_plus_strand)
) %>%
  #If there is a non-NULL value in the joined qual/sa/sm/sx columns, use those as the new corresponding opposite_strand values.
  mutate(
    qual.opposite_strand = if_else(map_lgl(qual.tmp, is.null), qual.opposite_strand, qual.tmp),
    sa.opposite_strand = if_else(map_lgl(sa.tmp, is.null), sa.opposite_strand, sa.tmp),
    sm.opposite_strand = if_else(map_lgl(sm.tmp, is.null), sm.opposite_strand, sm.tmp),
    sx.opposite_strand = if_else(map_lgl(sx.tmp, is.null), sx.opposite_strand, sx.tmp)
  ) %>%
  select(-c(qual.tmp,sa.tmp,sm.tmp,sx.tmp))

cat("DONE\n")

######################
### Save output
######################
cat("## Saving output...")

#In calls, change start/end_refspace to start/end for consistency with bam.
calls <- calls %>%
	rename(
		start = start_refspace,
		end = end_refspace
	)

#Save output
qs_save(
  list(
  	run_metadata = run_metadata,
  	molecule_stats = molecule_stats,
  	bam = bam,
  	calls = calls
  	),
  outputFile
  )

cat("DONE\n")