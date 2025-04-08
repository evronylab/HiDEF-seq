#!/usr/bin/env Rscript

#extractVariants.R:
# Loads and formats aligned ccs bamFile in HiDEF-seq format RDS file that includes all required alignment and variant information for analysis.  
# Usage: extractVariants.R [bamFile] [sample_id] [configuration.yaml]

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
sample_id <- args[2]
yaml.config <- suppressWarnings(read.config(args[3]))

#Load BSgenome reference, and install it first if needed
if(! yaml.config$BSgenome$BSgenome_name %in% installed.genomes()){
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

#Output objects
output <- list()

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

#Extract ccs strand
bam.gr$ccs_strand <- bam.gr$qname %>% str_extract("(fwd|rev)$")

#Count number of molecules
output$stats$num_molecules_initial <- bam.gr$zm %>% n_distinct

cat("DONE\n")

######################
### Initial molecule filtering
######################
cat("#### Initial molecule filtering...")

# Keep only molecules with 1 forward and 1 reverse primary alignment on the same chromosome (flags = 0 and 16 are plus and minus strand primary alignments; flags =  2048 and 2064 are plus and minus strand supplementary alignments).
zmwstokeep <- bam.gr %>%
  as_tibble %>%
  select(zm,strand,flag,seqnames) %>%
  group_by(zm) %>%
  filter(sum(strand=="+") == 1 & sum(strand=="-") == 1) %>%
  ungroup %>%
  filter( flag==0 | flag==16 ) %>%
  group_by(zm) %>%
  filter(n_distinct(seqnames) == 1) %>%
  pluck("zm")

bam.gr <- bam.gr[bam.gr$zm %in% zmwstokeep,]
rm(zmwstokeep)

# Count number of molecules
output$stats$num_molecules_postalignmentfilter <- bam.gr$zm %>% n_distinct

# Keep only molecules aligned to analyzed chromosomes
bam.gr <- bam.gr[seqnames(bam.gr) %in% chroms_to_analyze]

# Count number of molecules
output$stats$num_molecules_postanalyzedchroms <- bam.gr$zm %>% n_distinct

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

 # Compute reciprocal overlap ractions
bam.gr.onlyranges.overlap <- pintersect(bam.gr.onlyranges.plus, bam.gr.onlyranges.minus, ignore.strand=TRUE)
bam.gr.onlyranges.overlap.plus  <- width(bam.gr.onlyranges.overlap) / width(bam.gr.onlyranges.plus)
bam.gr.onlyranges.overlap.minus  <- width(bam.gr.onlyranges.overlap) / width(bam.gr.onlyranges.minus)

zmwstokeep <- bam.gr.onlyranges.overlap$zm[
  (bam.gr.onlyranges.overlap.plus >= yaml.config$min_strand_overlap) & (bam.gr.onlyranges.overlap.minus >= yaml.config$min_strand_overlap)
  ]

 # Keep only molecules passing reciprocal overlap filter
bam.gr <- bam.gr[bam.gr$zm %in% zmwstokeep,]

 # Remove intermediate objects
rm(bam.gr.onlyranges, bam.gr.onlyranges.plus, bam.gr.onlyranges.minus, bam.gr.onlyranges.overlap, bam.gr.onlyranges.overlap.plus, bam.gr.onlyranges.overlap.minus, zmwstokeep)

# Count number of molecules
output$stats$num_molecules_postminstrandoverlap <- bam.gr$zm %>% n_distinct

cat("DONE\n")

######################
### Extract read indels
######################
#For each read, save an IRanges objects with the locations of insertions and deletions relative to reference space, and translated to reference coordinates using the read start position. Also save sizes of insertions and deletions relative to the reference space, since otherwise that information is lost for deletions.

cat("#### Extracting read indels...")
ccsIndels.fwd <- cigarRangesAlongReferenceSpace(fwd.ccs.df$cigar,ops=c("I","D"),with.ops=T)
ccsIndels.rev <- cigarRangesAlongReferenceSpace(rev.ccs.df$cigar,ops=c("I","D"),with.ops=T)
ccsIndels.fwd <- shift(ccsIndels.fwd,fwd.ccs.df$pos-1)
ccsIndels.rev <- shift(ccsIndels.rev,rev.ccs.df$pos-1)
ccsIndels.insertionlengths.fwd <- width(cigarRangesAlongQuerySpace(fwd.ccs.df$cigar,ops=c("I","D"),with.ops=T))
ccsIndels.insertionlengths.rev <- width(cigarRangesAlongQuerySpace(rev.ccs.df$cigar,ops=c("I","D"),with.ops=T))
cat("DONE\n")


####
#D. Extracting mismatch data from reads
####

cat("  Extracting mismatch data from reads...")

#Helper function to extract positions with specific CIGAR ops
extractpos_cigarops <- function(cigar,ops,query_or_reference,refpos){
	if(query_or_reference=="query"){
		x <- as.data.frame(cigarRangesAlongQuerySpace(cigar,ops="X"))
		if(nrow(x)==0){
		  return(data.table(ccsidx=numeric(),value=numeric()))
		}
	}else if(query_or_reference=="reference"){
		x <- as.data.frame(cigarRangesAlongReferenceSpace(cigar,ops="X",pos=refpos))
		if(nrow(x)==0){
		  return(data.table(ccsidx=numeric(),value=numeric()))
		}
	}
	
	x <- data.table(value=rep(x$start,times=x$width),group=rep(x$group,times=x$width))
	x[,onemore:= value + (0:(.N-1)), by=.(group, value)]
	x <- x %>% select(-value) %>% rename("value" = "onemore","ccsidx" = "group")
	
	return(x)
}

#qp: Query mismatch positions
qp <- extractpos_cigarops(bamfile$cigar,"X","query") %>% rename("qp" = "value")

#Check if no mutations remaining
if(nrow(qp)==0){
	cat("  Note: no mismatches in selected chromosomes!\n")
  
  mismatches.df <- data.table(chrom=factor(),qn=character(),qp=integer(),qq=integer(),rn=character(),rp=integer(),sa=integer(),sm=integer(),strand=factor(),sx=integer(),zm=integer(),zmwidx=integer())
  
  rm(bamfile,qp)

}else{
  #qn: Query mismatch base sequences, relative to reference plus strand
  qn <- extractAt(bamfile$seq,cigarRangesAlongQuerySpace(bamfile$cigar,ops="X")) %>%
  			as("CharacterList") %>%
  			unstrsplit("") %>%
  			str_c(collapse="") %>%
  			str_split(pattern="") %>%
  			unlist()
  
  #If keepseqdata=FALSE, remove seq data, which is not necessary for the subsequent analysis
  if(!keepseqdata){
    bamfile$seq <- NULL
  }
  
  #qq: Query mismatch base qualities
  qq <- extractAt(bamfile$qual,cigarRangesAlongQuerySpace(bamfile$cigar,ops="X")) %>%
  			as("CharacterList") %>%
  			unstrsplit("") %>%
  			str_c(collapse="") %>%
  			str_split(pattern="") %>%
  			unlist() %>%
  			PhredQuality() %>%
  			as("IntegerList") %>%
  			unlist()
  
  #rp: Reference mismatch positions
  rp <- extractpos_cigarops(bamfile$cigar,"X","reference",bamfile$pos) %>% rename("rp" = "value")
  
  #rn: Reference mismatch base sequences, relative to reference plus strand
  mismatch.counts <- rp %>%
  	group_by(ccsidx) %>% 
  	summarize(count=n()) %>%
  	left_join(data.frame(ccsidx=1:length(bamfile$rname)),.,by="ccsidx") %>%
  	replace_na(list(count=0)) %>%
  	select(count) %>%
  	pull()
  
  rn <- data.frame(seqnames=rep(bamfile$rname,times=mismatch.counts),start=rp$rp,end=rp$rp,strand="+") %>%
  	makeGRangesFromDataFrame() %>%
  	getSeq(eval(parse(text=yaml.config$BSgenomepackagename)),.) %>%
  	as.character()
  
  #sa: the number of subread alignments that span each CCS read position
  # sa is run-length encoded (rle) in the BAM file as length,value pairs. We inverse the rle encoding back to standard per-position values
  # sa needs to be reversed for reads aligned to genome minus strand
  
   #Extract all sa values and inverse rle.
  sa.lengths <- lapply(bamfile$tag$sa,function(x){x[c(TRUE, FALSE)]})
  sa.values <- lapply(bamfile$tag$sa,function(x){x[c(FALSE, TRUE)]})
  sa.full <- mapply(function(x,y){inverse.rle(list(lengths=x,values=y))},x=sa.lengths,y=sa.values,USE.NAMES=FALSE)
  
   #Reverse order for reads aligned to genome minus strand
  sa.full[bamfile$strand=="-"] <- lapply(sa.full[bamfile$strand=="-"],rev)
  
   #Extract sa for all mismatches
  sa <- map2(qp$ccsidx,qp$qp,~ sa.full[[.x]][.y]) %>% unlist()
  
   #Store final sa.full object back into bamfile object
  bamfile$tag$sa <- sa.full
  
  rm(sa.lengths,sa.values,sa.full)
  
  #sm: the number of subreads that align as a match to each CCS read position
  # sm needs to be reversed for reads aligned to genome minus strand
  bamfile$tag$sm[bamfile$strand=="-"] <- lapply(bamfile$tag$sm[bamfile$strand=="-"],rev)
  sm <- map2(qp$ccsidx,qp$qp,~ bamfile$tag$sm[[.x]][.y]) %>% unlist()
  
  #sx: the number of subreads that align as a mismatch to each CCS read position
  # sx needs to be reversed for reads aligned to genome minus strand
  bamfile$tag$sx[bamfile$strand=="-"] <- lapply(bamfile$tag$sx[bamfile$strand=="-"],rev)
  sx <- map2(qp$ccsidx,qp$qp,~ bamfile$tag$sx[[.x]][.y]) %>% unlist()
  
  #Combine mismatch info
  mismatches.df <- qp %>%
  	mutate(qn=qn,qq=qq,rp=rp$rp,rn=rn,sa=sa,sm=sm,sx=sx,
  	       strand=rep(bamfile$strand,times=mismatch.counts),
  	       zm=rep(bamfile$tag$zm,times=mismatch.counts)) %>%
    select(-ccsidx) %>%
  	left_join(data.frame(zm=bamfile$tag$zm,chrom=bamfile$rname) %>% distinct(), by="zm") %>%
    left_join(data.frame(zm=fwd.ccs.df$zm,zmwidx=fwd.ccs.df$zmwidx) %>% distinct(), by="zm") %>%
    arrange(zmwidx,rp) %>%
    select(order(colnames(.)))
  
  rm(bamfile,qp,qn,qq,rp,rn,sa,sm,sx)
}

cat("DONE\n")


####
#E. Load VCF files
####
cat("  Loading VCF files...\n")

#Annotate with each VCF file (make a list of tables, each table with the same number of rows as mismatches.df)
vcfIndels <- list()
vcfSNVs <- list()

for(i in vcffiles){
  cat("     Loading VCF File:",i,"...")
  
	#Split multi-allelic sites (bcftools norm -m -both -f [fastaref]), and keep only variants in chromosomes present in this bamfile chunk.
	tmpvcf <- tempfile(tmpdir=getwd(),pattern=".")
	chunkchroms <- paste(unique(fwd.ccs.df$rname),collapse=",")
	system(paste("/bin/bash -c",shQuote(paste(yaml.config$bcftools_bin,"norm -m -both -f",fastaref,"-r",chunkchroms,i," 2>/dev/null >",tmpvcf))))
	
  #Load vcf file, filter to keep only chromosomes being analyzed, and keep only columns CHROM, POS, REF, ALT, QUAL, FILTER, GT, GQ
  vcffile <- read.vcfR(tmpvcf,convertNA=FALSE,verbose=FALSE)
  file.remove(tmpvcf)
  
  vcftable <- cbind(vcffile@fix,
  	data.frame(
	  	GT=as.character(extract.gt(vcffile,element="GT",convertNA=FALSE,IDtoRowNames=FALSE)),
	  	GQ=as.numeric(extract.gt(vcffile,element="GQ",as.numeric=TRUE,convertNA=FALSE,IDtoRowNames=FALSE)),
	  	AD=as.character(extract.gt(vcffile,element="AD",convertNA=FALSE,IDtoRowNames=FALSE))
  	)) %>%
  	filter(ALT!="*") %>%
    filter(CHROM %in% chroms) %>%
  	as.data.frame() %>%
  	mutate_at(c("POS","QUAL"),as.numeric) %>%
  	separate(AD,c("AD1","AD2"),",") %>%
  	mutate_at(c("AD1","AD2"),as.numeric) %>%
  	mutate(Depth=AD1+AD2, VAF=AD2/Depth) %>%
  	rename(Ref = REF, Alt = ALT) %>%
  	select(-c(ID,INFO,AD1,AD2))
  
  rm(vcffile)

  #Separate indels from SNVs and remove '*' Alt alleles, which indicates overlapping deletions, which are not needed because every deletion already has a separate vcf entry. Then annotate with Depth (sum of AD) and VAF (AD of variant / sum of Depth).
  vcfIndels[[i]] <- vcftable %>%
  	filter(nchar(Ref)!=nchar(Alt))
  
  vcfSNVs[[i]] <- vcftable %>%
  	filter(nchar(Ref)==1 & nchar(Alt)==1)

  rm(vcftable)
  
  #Convert Indels to Granges object. Also, change ranges for deletion variants to accurately reflect position of deletion. Insertion variants are annotated with width=-1 at the insertion site (insertion is immediately to the right of the first Ref base).
  vcfIndels[[i]] <- makeGRangesFromDataFrame(vcfIndels[[i]],keep.extra.columns=TRUE,ignore.strand=TRUE,seqnames.field="CHROM",start.field="POS",end.field="POS")
  lengthdiff <- nchar(vcfIndels[[i]]$Ref) - nchar(vcfIndels[[i]]$Alt)
  end(vcfIndels[[i]])[lengthdiff>0] <- start(vcfIndels[[i]])[lengthdiff>0] + lengthdiff[lengthdiff>0]
  start(vcfIndels[[i]]) <- start(vcfIndels[[i]]) + 1
  
  cat("DONE\n")
}


####
#F. Collecting extracted data into bamfilezmw_all
####
cat("  Collecting extracted data into bamfilezmw_all...")

bamfilezmw <- list()

#Mismatches/mutations: Create fwd/rev/zmw data frames for mismatches found in one or both strands. Then make granges and add back rp as a column for simplicity even though it is already encoded in the start field. Note, strand ("+" or "-") corresponds to the strand of the CCS read in which the mismatch was found, i.e. the reference genome strand that the CCS read aligned to, which is also the synthesized strand. For zmw.mismatches.granges, strand is set to *. However, qn and rn always correspond to the reference genome plus strand.

fwd.mismatches.df <- mismatches.df %>% filter(strand=="+")
rev.mismatches.df <- mismatches.df %>% filter(strand=="-")
zmw.mismatches.df <- inner_join(fwd.mismatches.df,rev.mismatches.df,by=c("zm","zmwidx","chrom","rp","rn","qn"),suffix=c(".fwd",".rev")) %>% select(order(colnames(.)))

rm(mismatches.df)

	#Remove zmw mutations from fwd and rev data frames
fwd.mismatches.df <- anti_join(fwd.mismatches.df,zmw.mismatches.df,by=c("zm","zmwidx","chrom","rp","rn","qn"))
rev.mismatches.df <- anti_join(rev.mismatches.df,zmw.mismatches.df,by=c("zm","zmwidx","chrom","rp","rn","qn"))

fwd.mismatches.granges <- makeGRangesFromDataFrame(fwd.mismatches.df,seqnames.field="chrom",start.field="rp",end.field="rp",strand="strand",keep.extra.columns=TRUE)
rev.mismatches.granges <- makeGRangesFromDataFrame(rev.mismatches.df,seqnames.field="chrom",start.field="rp",end.field="rp",strand="strand",keep.extra.columns=TRUE)
zmw.mismatches.granges <- makeGRangesFromDataFrame(zmw.mismatches.df,seqnames.field="chrom",start.field="rp",end.field="rp",ignore.strand=TRUE,keep.extra.columns=TRUE)

fwd.mismatches.granges$rp <- start(fwd.mismatches.granges)
rev.mismatches.granges$rp <- start(rev.mismatches.granges)
zmw.mismatches.granges$rp <- start(zmw.mismatches.granges)

bamfilezmw$fwd.mismatches.granges <- fwd.mismatches.granges
bamfilezmw$rev.mismatches.granges <- rev.mismatches.granges
bamfilezmw$zmw.mismatches.granges <- zmw.mismatches.granges

rm(fwd.mismatches.granges,rev.mismatches.granges,zmw.mismatches.granges,zmw.mismatches.df.strand.fwd)

#vcfSNV.fwd/rev/zmw: mismatches annotated with vcfSNVs

annotate.mismatches.vcfSNV <- function(mismatches,vcfSNV,colstokeep){
  return(left_join(mismatches,vcfSNV,by=c("chrom"="CHROM","rp"="POS","rn"="Ref","qn"="Alt")) %>%
           select(all_of(colstokeep))
         )
}

vcfSNV.fwd <- list()
vcfSNV.rev <- list()
vcfSNV.zmw <- list()

for(i in names(vcfSNVs)){
  colstokeep <- c("QUAL","FILTER","GT","GQ","Depth","VAF")
	vcfSNV.fwd[[i]] <- annotate.mismatches.vcfSNV(fwd.mismatches.df,vcfSNVs[[i]],colstokeep)
	vcfSNV.rev[[i]] <- annotate.mismatches.vcfSNV(rev.mismatches.df,vcfSNVs[[i]],colstokeep)
	vcfSNV.zmw[[i]] <- annotate.mismatches.vcfSNV(zmw.mismatches.df,vcfSNVs[[i]],colstokeep)
}

bamfilezmw$vcfSNV.fwd <- vcfSNV.fwd
bamfilezmw$vcfSNV.rev <- vcfSNV.rev
bamfilezmw$vcfSNV.zmw <- vcfSNV.zmw

rm(vcfSNVs,vcfSNV.fwd,vcfSNV.rev,vcfSNV.zmw,fwd.mismatches.df,rev.mismatches.df,zmw.mismatches.df)

#vcfIndels
bamfilezmw$vcfIndels <- vcfIndels

rm(vcfIndels)

#CCS Indels
bamfilezmw$ccsIndels.fwd <- ccsIndels.fwd
bamfilezmw$ccsIndels.rev <- ccsIndels.rev
bamfilezmw$ccsIndels.insertionlengths.fwd <- ccsIndels.insertionlengths.fwd
bamfilezmw$ccsIndels.insertionlengths.rev <- ccsIndels.insertionlengths.rev

rm(ccsIndels.fwd,ccsIndels.rev,ccsIndels.insertionlengths.fwd,ccsIndels.insertionlengths.rev)

cat("DONE\n")

cat("#### Saving output RDS file...")
qsave(bamfilezmw_all,bamfilezmw_all_outputpath)
cat("DONE\n")