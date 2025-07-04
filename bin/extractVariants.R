#!/usr/bin/env -S Rscript --vanilla

#extractVariants.R:
# Perform initial molecule filtering and extract variants, including all variant information required for analysis.  
# Note: extracts variants of all variant_types across all chromosomes, not just those in configured chromgroups. filterVariants then removes variants outside configured chromgroups.

cat("#### Running extractVariants ####\n")

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
strand_levels <- c("+","-")

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

#Load variant types, and filter to make a table only for MDB variant types
variant_types <- yaml.config$variant_types %>%
  enframe(name=NULL) %>%
  unnest_wider(value) %>%
  unnest_longer(SBSindel_call_types) %>%
  unnest_wider(SBSindel_call_types)

variant_types_MDB <- variant_types %>%
  filter(variant_class == "MDB")

#Load the BSgenome reference
suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name,character.only=TRUE,lib.loc=yaml.config$cache_dir))

cat("DONE\n")

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
        variant_types_MDB %>% pluck("variant_bam_scoretag")
        )
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

#Recalculate isize, since pbmm2 has a bug that miscalculates it
bam.df$isize <- cigarWidthAlongReferenceSpace(bam.df$cigar)

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

# Arrange bam.gr by run id, ZMW id and strand
bam.gr <- bam.gr[order(bam.gr$run_id,bam.gr$zm, strand(bam.gr))]

#Count number of molecules per run
molecule_stats <- mcols(bam.gr)[,c("run_id","movie_id","zm")] %>%
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
  select(run_id,zm) %>%
	distinct

bam.gr <- bam.gr[bam.gr$run_id %in% zmwstokeep$run_id & bam.gr$zm %in% zmwstokeep$zm,]
rm(zmwstokeep)

# Count number of molecules per run
molecule_stats <- molecule_stats %>%
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

# Keep only molecules with plus and minus strand alignment overlap >= min_strand_overlap (reciprocal or both plus and minus strand alignments)

 # Extract plus and minus strand reads separately. Due to prior sorting by ZMW id and strand, these are guaranteed to have the same ZMW id order.
bam.gr.onlyranges.plus <- bam.gr[strand(bam.gr) == "+", c("zm","run_id")]
bam.gr.onlyranges.minus <- bam.gr[strand(bam.gr) == "-", c("zm","run_id")]

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
rm(bam.gr.onlyranges.plus, bam.gr.onlyranges.minus, bam.gr.onlyranges.overlap, zmwstokeep)

# Count number of molecules per run
molecule_stats <- molecule_stats %>%
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

######################
### Define function to extract variants
######################
extract_variants <- function(bam.gr.input, sa.input, sm.input, sx.input, variant_class.input, variant_type.input, cigar.ops.input = NULL, variant_bam_scoretag.input = NULL, variant_bam_scoretag_min.input = NULL){
  #cigar.ops.input for the variant_class: "X" for SBS; "I" for insertion; "D" for deletion; NULL for MDB
  #variant_bam_scoretag.input TAG in the MDB containing the MDB score data
  #variant_bam_scoretag_min.input: minimum score of MDBs to extract
  
  #Initialize empty variants tibble
  variants.out <- tibble(
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
    sa=list(),
    sm=list(),
    sx=list(),
    variant_class=factor(),
    variant_type=factor(),
    indel_width=integer()
    )
  
  #Extract variant positions
  if(variant_class.input %in% c("SBS","indel")){
    
    #Extract variant positions in query space, and name variants with the position to retain this in the final data frame
    cigar.queryspace.var <- cigarRangesAlongQuerySpace(bam.gr.input$cigar,ops=cigar.ops.input)
    cigar.queryspace.var <- cigar.queryspace.var %>%
      unlist(use.names=FALSE) %>%
      setNames(str_c(start(.),end(.),sep="_")) %>%
      relist(.,cigar.queryspace.var)
    
    #Extract variant positions in reference space
    cigar.refspace.var <- cigarRangesAlongReferenceSpace(bam.gr.input$cigar,ops=cigar.ops.input,pos=start(bam.gr.input))
    
  }else if(variant_class.input == "MDB"){
    
  }
  
  #Check if no variants
  if(cigar.queryspace.var %>% as.data.frame %>% nrow == 0){
    cat("  No variants extracted of type",vartype_label,"!\n")
    return(variants.out)
  }
  
  #Variant positions in query space
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
  
  if(variant_class.input %in% c("SBS","MDB")){
    var_queryspace <- var_queryspace %>%
      uncount(weights=width) %>%
      select(zm,strand,start_queryspace) %>%
      group_by(zm,strand,start_queryspace) %>%
      mutate(
        start_queryspace = start_queryspace + row_number() - 1L,
        end_queryspace = start_queryspace
      ) %>%
      ungroup
  }else if(variant_class.input == "indel"){
    var_queryspace <- var_queryspace %>%
      rename(insertion_width = width) %>%
      mutate(insertion_width = insertion_width %>% na_if(0)) %>%
      select(-group,-names) %>%
      as_tibble
  }
  
  #Variant positions in reference space, annotated with positions in query space for later joining
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
  
  if(variant_class.input %in% c("SBS","MDB")){
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
  }else if(variant_class.input == "indel"){
    var_refspace <- var_refspace %>%
      rename(deletion_width = width) %>%
      mutate(deletion_width = deletion_width %>% na_if(0)) %>%
      select(-group) %>%
      as_tibble
  }
  
  rm(cigar.refspace.var)
  
  #Variant REF base sequences, relative to reference plus strand
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
  
  #Variant ALT base sequences, relative to reference plus strand
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
  
  if(variant_class.input %in% c("SBS","MDB")){
    var_alt <- var_alt %>%
      separate_longer_position(alt_plus_strand,width=1) %>%
      group_by(zm,strand,start_queryspace) %>%
      mutate(
        start_queryspace = start_queryspace + row_number() - 1L,
        end_queryspace = start_queryspace
      ) %>%
      ungroup
  }
  
  #Variant base qualities
  if(variant_class.input %in% c("SBS","MDB")){
    cigar.queryspace.var.forqualdata <- cigar.queryspace.var
  }else if(variant_class.input == "indel"){
    #Indel base qualities. For insertions, get base qualities of inserted bases in query space, relative to reference plus strand. For deletions, get base qualities of left and right flanking bases in query space, relative to reference plus strand.
    cigar.queryspace.var.forqualdata <- cigar.queryspace.var %>% unlist
    deletionidx <- start(cigar.queryspace.var.forqualdata) > end(cigar.queryspace.var.forqualdata)
    newstart <- end(cigar.queryspace.var.forqualdata)
    newend <- start(cigar.queryspace.var.forqualdata)
    start(cigar.queryspace.var.forqualdata[deletionidx]) <- newstart[deletionidx]
    end(cigar.queryspace.var.forqualdata[deletionidx]) <- newend[deletionidx]
    cigar.queryspace.var.forqualdata <- cigar.queryspace.var.forqualdata %>% relist(.,cigar.queryspace.var)
    
    rm(deletionidx,newstart,newend)
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
  
  if(variant_class.input %in% c("SBS","MDB")){
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
      qual = qual %>% PhredQuality %>% as("IntegerList") %>% as.list %>% set_names(NULL)
    )
  
  rm(cigar.queryspace.var)
  
  #Extract sa/sm/sx tags
  # For insertions, get base qualities of inserted bases in query space, relative to reference plus strand. For deletions, get base qualities of left and right flanking bases in query space, relative to reference plus strand.
  
  if(variant_class.input %in% c("SBS","MDB")){
    #Convert variant positions in query space to list for faster retrieval below of sa, sm, sx tags
    var_queryspace.list <- var_queryspace %>%
      mutate(zm_strand = str_c(zm,strand,sep="_")) %>%
      select(zm_strand,start_queryspace) %>%
      as.data.table %>%
      split(by="zm_strand",keep.by=F) %>%
      lapply(unlist,use.names=F)
    
    #sa
    var_sa <- map2(sa.input[var_queryspace.list %>% names], var_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
      enframe %>%
      separate(name,sep="_",into=c("zm","strand")) %>%
      unnest_longer(col=value, values_to="sa", indices_to="start_queryspace") %>%
      mutate(
        zm = as.integer(zm),
        strand = factor(strand,levels=strand_levels),
        start_queryspace = as.integer(start_queryspace),
        end_queryspace = start_queryspace,
        sa = sa %>% as.list %>% set_names(NULL)
      )
    
    #sm
    var_sm <- map2(sm.input[var_queryspace.list %>% names], var_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
      enframe %>%
      separate(name,sep="_",into=c("zm","strand")) %>%
      unnest_longer(col=value, values_to="sm", indices_to="start_queryspace") %>%
      mutate(
        zm = as.integer(zm),
        strand = factor(strand,levels=strand_levels),
        start_queryspace = as.integer(start_queryspace),
        end_queryspace = start_queryspace,
        sm = sm %>% as.list %>% set_names(NULL)
      )
    
    #sx
    var_sx <- map2(sx.input[var_queryspace.list %>% names], var_queryspace.list, ~ .x[.y] %>% set_names(.y)) %>%
      enframe %>%
      separate(name,sep="_",into=c("zm","strand")) %>%
      unnest_longer(col=value, values_to="sx", indices_to="start_queryspace") %>%
      mutate(
        zm = as.integer(zm),
        strand = factor(strand,levels=strand_levels),
        start_queryspace = as.integer(start_queryspace),
        end_queryspace = start_queryspace,
        sx = sx %>% as.list %>% set_names(NULL)
      )
    
  }else if(variant_class.input == "indel"){
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
    
    setkey(indels_queryspace_pos, zm_strand, start_end_queryspace)
    
    #Extract data
    indels_queryspace_pos[, sa_val := sa.input[[zm_strand[1]]][pos_queryspace], by = zm_strand]
    indels_queryspace_pos[, sm_val := sm.input[[zm_strand[1]]][pos_queryspace], by = zm_strand]
    indels_queryspace_pos[, sx_val := sx.input[[zm_strand[1]]][pos_queryspace], by = zm_strand]
    
    rm(cigar.queryspace.var.forqualdata)
    
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
  }
  
  #Combine variant info with prior empty initialized tibble and return result
  if(variant_class.input %in% c("SBS","MDB")){
    variants.out <- variants.out %>%
      bind_rows(
        var_queryspace %>%
          left_join(var_refspace, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_alt, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_qual, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_sa, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_sm, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_sx, by=join_by(zm,strand,start_queryspace,end_queryspace))
      )
  }else if(variant_class.input == "indel"){
    variants.out <- variants.out %>%
      bind_rows(
        var_queryspace %>%
          left_join(var_refspace, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_alt, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(var_qual, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          left_join(indels_queryspace_pos, by=join_by(zm,strand,start_queryspace,end_queryspace)) %>%
          mutate(
            indel_width = if_else(is.na(insertion_width), deletion_width, insertion_width)
          ) %>%
          select(-deletion_width,-insertion_width)
      )
  }
  
  #Set variant_class and variant_type
  variants.out <- variants.out %>%
    mutate(
      variant_class = !!variant_class.input %>% factor,
      variant_type = !!variant_type.input %>% factor
      )
  
  return(variants.out)

}


#########################################
### Loop to extract variants from each run_id/movie_id separately
#########################################

#Split bam.gr by run_id
bam.gr <- bam.gr %>% split(bam.gr$run_id)

#Create output list
variants.df <- list()

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
  ### Extract variants
  ######################
  cat(" ## Extracting single base substitutions...")
  
  #SBS
  variants.df[[i]] <- extract_variants(
    bam.gr.input = bam.gr[[i]],
    sa.input = sa,
    sm.input = sm,
    sx.input = sx,
    variant_class.input = "SBS",
    variant_type.input = "SBS",
    cigar.ops.input = "X"
    )
  
  cat("DONE\n")
  
  #Insertions
  cat(" ## Extracting insertions...")
  
  variants.df[[i]] <- variants.df[[i]] %>%
    bind_rows(
      extract_variants(
        bam.gr.input = bam.gr[[i]],
        sa.input = sa,
        sm.input = sm,
        sx.input = sx,
        variant_class.input = "indel",
        variant_type.input = "insertion",
        cigar.ops.input = c("I")
      )
    )

  cat("DONE\n")
  
  #Deletions
  cat(" ## Extracting deletions...")
  
  variants.df[[i]] <- variants.df[[i]] %>%
    bind_rows(
      extract_variants(
        bam.gr.input = bam.gr[[i]],
        sa.input = sa,
        sm.input = sm,
        sx.input = sx,
        variant_class.input = "indel",
        variant_type.input = "deletion",
        cigar.ops.input = c("D")
      )
    )
  
  cat("DONE\n")
  
  
  #MDBs
  for(j in seq_len(nrow(variant_types_MDB))){
    
    cat(" ## Extracting",variant_types_MDB %>% pluck("variant_type",j),"...")
    
    variants.df[[i]] <- variants.df[[i]] %>%
      bind_rows(
        extract_variants(
          bam.gr.input = bam.gr[[i]],
          sa.input = sa,
          sm.input = sm,
          sx.input = sx,
          variant_class.input = "MDB",
          variant_type.input = variant_types_MDB %>% pluck("variant_type",j),
          variant_bam_scoretag.input = variant_types_MDB %>% pluck("variant_bam_scoretag",j),
          variant_bam_scoretag_min.input = variant_types_MDB %>% pluck("variant_bam_scoretag_min",j)
        )
      )
    
    cat("DONE\n")
  }
  
}

######################
### Format and save output
######################
cat("## Formatting and saving output...")

#Collapse bam.gr back to a single GRanges object and remove seq and qual data that is no longer needed
bam.gr <- bam.gr %>% unlist(use.names=F)
bam.gr$seq <- NULL
bam.gr$qual <- NULL

#Combine variant calls of all runs
variants.df <- variants.df %>%
  bind_rows(.id="run_id") %>%
  mutate(run_id=factor(run_id))

#Annotate for each variant position if it has duplex coverage
variants.df <- variants.df %>% left_join(
  variants.df %>%
  #Convert to GRanges
  select(run_id,zm,seqnames,start_refspace,end_refspace) %>%
  distinct %>% #In case there is an SBS and MDB at the same location, only want to join one row for the position to bam.gr
  makeGRangesFromDataFrame(
    seqnames.field="seqnames",
    start.field="start_refspace",
    end.field="end_refspace",
    keep.extra.columns=TRUE,
    seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
  ) %>%

  #Count for each variant position how many bam.gr reads from its zmw it overlaps (1 or 2)
  join_overlap_left_within(bam.gr[,c("run_id","zm")],suffix=c("",".bam.gr")) %>%
  filter(
    run_id == run_id.bam.gr,
    zm == zm.bam.gr
  ) %>%
  as_tibble %>%
  rename(
    start_refspace = start,
    end_refspace = end
    ) %>%
  group_by(run_id,zm,start_refspace,end_refspace) %>%
  summarize(
    duplex_coverage = if_else(n() == 2, TRUE, FALSE),
    .groups="drop"
    ),
  by = join_by(run_id,zm,start_refspace,end_refspace)
)

#Annotate for each variants its trinucleotide context when variant_class = SBS or MDB. For indels, this is set to NA.
variants.df <- variants.df %>%
	left_join(
		variants.df %>%
			filter(variant_class %in% c("SBS","MDB")) %>%
			select(seqnames,start_refspace,end_refspace,variant_class) %>%
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
		by = join_by(seqnames,start_refspace,end_refspace,variant_class)
	)

#Annotate for each variants its ref, alt, and trinucleotide context sequences relative to the reference plus strand, synthesized strand, and template strand.
variants.df <- variants.df %>%
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

#Annotate for each variant (for all variant classes: SBS, indel, MDB) the base sequence on the opposite strand based on reference space coordinates, taking into account if there was an SBS or an indel called on the opposite strand, regardless if there was an MDB called on the opposite strand.
# Possible overlaps are SBS/MDB with SBS, SBS/MDB with deletion, insertion with insertion, insertion with deletion, deletion with SBS, deletion with insertion, deletion with deletion.
# If a deletion on the opposite strand only partially overlaps an analyzed deletion, annotate it as a deletion on the opposite strand.
# For deletions called on both strands, annotate if the deletions have the same coordinates on both strands (deletion.bothstrands.startendmatch).
# Note, for an analyzed deletion, it is possible for > 1 SBS, insertion or deletion to overlap it on the opposite strand.
# Whether and which MDBs were called on the same or opposite strand will be annotated later in the filterVariants step after further filtering MDBs, since filtering out an MDB won't trigger filtering of SBS and indel calls, whereas filtering an SBS or indel on any one strand will anyway cause filtering of variants on the opposite strand, so it doesn't matter if we annotate here opposite strand SBS and indel variants before filtering.

 #Make GRanges of variants. Swap start and end for insertions, since the subsequent join otherwise cannot join two insertions with the same coordinates (this won't affect later joining, since original start_refspace and end_refspace are used to filter the joins).
variants.gr <- variants.df %>%
  mutate(
    start = if_else(variant_type=="insertion",end_refspace,start_refspace),
    end = if_else(variant_type=="insertion",start_refspace,end_refspace),
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

 #Annotate variants GRanges with opposite strand variants
variants.gr <- variants.gr %>%
  join_overlap_left(
    variants.gr %>% filter(variant_class %in% c("SBS","indel")),
    suffix = c("",".opposite_strand")
  ) %>%
  filter(
    run_id == run_id.opposite_strand,
    zm == zm.opposite_strand,
    strand_copy != strand_copy.opposite_strand,
    (
      variant_class %in% c("SBS","MDB") & variant_type.opposite_strand == "SBS" &
      start_refspace == start_refspace.opposite_strand &
      end_refspace == end_refspace.opposite_strand
    ) |
    (
      variant_class %in% c("SBS","MDB") & variant_type.opposite_strand == "deletion" &
        start_refspace >= start_refspace.opposite_strand &
        end_refspace <= end_refspace.opposite_strand
    ) |
    (
      variant_type == "insertion" & variant_type.opposite_strand == "insertion" &
        start_refspace == start_refspace.opposite_strand &
        end_refspace == end_refspace.opposite_strand
    ) |
    (
      variant_type == "insertion" & variant_type.opposite_strand == "deletion" &
        start_refspace <= end_refspace.opposite_strand &
        end_refspace >= start_refspace.opposite_strand
    ) |
    (
      variant_type == "deletion" & variant_type.opposite_strand == "SBS" &
      start_refspace <= start_refspace.opposite_strand &
      end_refspace >= end_refspace.opposite_strand
    ) |
    (
      variant_type == "deletion" & variant_type.opposite_strand == "insertion" &
        start_refspace <= end_refspace.opposite_strand &
        end_refspace >= start_refspace.opposite_strand
    ) |
    (
      variant_type == "deletion" & variant_type.opposite_strand == "deletion" &
        start_refspace <= end_refspace.opposite_strand &
        end_refspace >= start_refspace.opposite_strand
    )
  ) %>%
  mutate(
    deletion.bothstrands.startendmatch = case_when(
      variant_type != "deletion" | variant_type.opposite_strand != "deletion" ~ NA,
    	start_refspace == start_refspace.opposite_strand & end_refspace == end_refspace.opposite_strand ~ TRUE,
      start_refspace != start_refspace.opposite_strand | end_refspace != end_refspace.opposite_strand ~ FALSE
    )
  )

 #Join to variants.df with opposite strand information.
 # Note that for deletions with one partially overlapping deletion on the opposite strand, alt_plus_strand.opposite_strand will still equal "" even though part of the analyzed deletion's sequence is still present in the opposite strand. Likewise for deletions with > 1 overlapping SBS, insertion, or deletion on the opposite strand, alt_plus_strand.opposite_strand will contain the sequence of one random one of these. These two issues would be difficult to fix, and are not critical for downstream analysis.
variants.df <- variants.df %>%
	left_join(
		variants.gr %>%
			as_tibble %>%
			select(
				run_id,zm,start_refspace,end_refspace,strand,
				variant_class.opposite_strand,
				variant_type.opposite_strand,
				alt_plus_strand.opposite_strand,
				deletion.bothstrands.startendmatch
			),
		by=join_by(run_id,zm,start_refspace,end_refspace,strand),
		multiple="any" # Required due to possibility of > 1 variant on the opposite strand overlapping an analyzed deletion
	)

rm(variants.gr)
  
#Annotate for each variant (for all variant classes: SBS, indel, MDB) its 'SBSindel_call_type'.
 #No duplex coverage ("nonduplex")
 #Duplex coverage:
  #  Variant strand | Opposite strand (SBS or indel) overlapping based on reference space coordinates
  #       N         |       N                       			=> "match" (possible only for MDB)
  #       N         |       Y                       			=> "mismatch-os" (possible only for MDB)
  #       Y         |       N                     				=> "mismatch-ss"
  #       Y         |       Y (non-complementary change)  => "mismatch-ds"
  #       Y         |       Y (complementary change)      => "mutation"

variants.df <- variants.df %>%
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
  arrange(run_id,zm,start_refspace,end_refspace,variant_type,strand)

#Convert variants.df to GRanges
variants.gr <- variants.df %>%
	makeGRangesFromDataFrame(
		seqnames.field="seqnames",
		start.field="start_refspace",
		end.field="end_refspace",
		strand.field="strand",
		keep.extra.columns=TRUE,
		seqinfo=yaml.config$BSgenome$BSgenome_name %>% get %>% seqinfo
	)

rm(variants.df)

#Save output
qs_save(
  list(
  	run_metadata = run_metadata,
  	molecule_stats = molecule_stats,
  	bam.gr = bam.gr,
  	variants.gr = variants.gr
  	),
  outputFile
  )

cat("DONE\n")