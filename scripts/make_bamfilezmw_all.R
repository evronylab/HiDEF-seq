#!/usr/bin/env Rscript
# Usage: make_bamfilezmw_all.R [configuration.yaml]
#  - Uses [configuration.yaml] file per documentation
#  - Outputs bamfilezmw_all for subsequent analyses
#
# Prerequisites:
#  - R libraries: Rsamtools, GenomicAlignments, GenomicRanges, vcfR, plyr, configr, qs
#
# Configuration YAML: Configuration YAML shared with the other analysis scripts

######################
### Load required libraries
######################
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs))

######################
### Define function to load bamfiles and annotate them with VCF files
######################
load_bamfile <- function(bamfile,vcffiles,chroms,keepseqdata=FALSE){
  # keepseqdata: whether to store the full sequence of each read (seq column of BAM file) in the output, which takes up almost half the output size and is not used in the analysis since the mutation sequences are already extracted in a separate tag in the prior steps of the pipeline.

####
#A. Load file, including all tags.
####

cat("  Loading BAM file:",bamfile,"...")
bamfile <- scanBam(bamfile,param=ScanBamParam(what=scanBamWhat(),tag=c("bx","ec","np","rq","sn","we","ws","zm","qs","qe","bc","bq","cx","bl","bt","ql","qt","RG","mg","qp","qn","qq","rp","rn")))[[1]]
cat("DONE\n")

#Keep only reads from selected chromosomes
keepreads <- bamfile$rname %in% chroms
for(i in setdiff(names(bamfile),"tag"))
{
  bamfile[[i]] <- bamfile[[i]][keepreads]
}

for(i in names(bamfile$tag)){
  bamfile$tag[[i]] <- bamfile$tag[[i]][keepreads]
}

if(length(bamfile$qname)==0 | length(unlist(bamfile$tag$qp))==0){
	cat("  No reads or mutations in selected chromosomes! Proceeding to next sample...\n")
	return("No reads or mutations in selected chromosomes")
}

#Remove full read sequence data, which is not necessary for the subsequent analysis, if keepseqdata=FALSE
if(!keepseqdata){
  bamfile$seq <- NULL
}

cat("  Formatting bamfile data...")
#Split qn and rn tags by comma delimiter
bamfile$tag$qn <- strsplit(bamfile$tag$qn,",")
bamfile$tag$rn <- strsplit(bamfile$tag$rn,",")

#Split qq tag by keeping odd numbered characters, since some quality values can equal ','.
#Convert base qualities to numeric.
bamfile$tag$qq <- lapply(as(PhredQuality(gsub("(.).","\\1",bamfile$tag$qq)),"IntegerList"),as.numeric)

#Add tag that equals TRUE for all positions that are not 'N' in the reference.
bamfile$tag$noN <- lapply(bamfile$tag$rn,function(x){x!="N"})

#For each CCS, save an IRanges objects with the locations of insertions and deletions relative to reference space, and translated to reference coordinates using the read start position. Also save sizes of insertions and deletions relative to the reference space, since otherwise that information is lost for deletions.
bamfile$ccsIndels <- cigarRangesAlongReferenceSpace(bamfile$cigar,ops=c("I","D"),with.ops=T)
bamfile$ccsIndels <- shift(bamfile$ccsIndels,bamfile$pos-1)
bamfile$ccsIndels.insertionlengths <- width(cigarRangesAlongQuerySpace(bamfile$cigar,ops=c("I","D"),with.ops=T))

cat("DONE\n")


####
#B. Load vcf files and annotate bamfile with variants.
####
cat("  Annotating bamfile with VCF data...\n")

#Format bamfile variants as chrom_pos_ref_alt in column 2, and an index in column 1 corresponding to the row number of the read in bamfile.
bamvariants <- mapply(paste0,bamfile$rname,"_",bamfile$tag$rp,"_",bamfile$tag$rn,"_",bamfile$tag$qn,SIMPLIFY=FALSE)
bamvariants <- lapply(bamvariants,function(x){if(grepl("___$",x[1])){character(0)}else{x}})

bamvariants <- data.frame(matrix(unlist(strsplit(unlist(lapply(seq_along(bamvariants),function(i){if(length(bamvariants[[i]])>0){paste(i,1:length(bamvariants[[i]]),bamvariants[[i]])}} ))," ")),ncol=3,byrow=TRUE))
bamvariants[,1] <- as.numeric(bamvariants[,1])
bamvariants[,2] <- as.numeric(bamvariants[,2])
colnames(bamvariants) <- c("readindex","varindex","id")

#Annotate each VCF file
vcfIndels <- list()
for(vcffilename in vcffiles){
  cat("     Annotating VCF File:",vcffilename,"...")
  
  #Load vcf file
  vcffile <- read.vcfR(vcffilename,convertNA=FALSE,verbose=FALSE)
  vcftable <- cbind(vcffile@fix,vcffile@gt)
  
  #Create table of vcf variants with a separate row for each variant (including splitting each allele of multi-allelic variants into a separate row), annotated with id = chrom_pos_ref_alt, FILTER, index of ALT allele (i.e. 1, 2, 3, needed for multi-allelic variants), GT (genotype), VAF (AD of variant / sum of AD), Depth (sum of AD), GQ. This step also removes indels by requiring the length of REF and ALT allele to be equal to 1.
  variantalt <- strsplit(vcftable[,5],",")
  variantaltnum <- lapply(variantalt,length)
  variantchrom <- mapply(function(x,y){rep(x,y)},vcftable[,"CHROM"],variantaltnum,USE.NAMES=FALSE)
  variantpos <- mapply(function(x,y){rep(x,y)},as.numeric(vcftable[,"POS"]),variantaltnum,USE.NAMES=FALSE)
  variantqual <- mapply(function(x,y){rep(x,y)},as.numeric(vcftable[,"QUAL"]),variantaltnum,USE.NAMES=FALSE)
  variantref <- mapply(function(x,y){rep(x,y)},vcftable[,"REF"],variantaltnum,USE.NAMES=FALSE)
  variantFILTER <- mapply(function(x,y){rep(x,y)},vcftable[,"FILTER"],variantaltnum,USE.NAMES=FALSE)
  variantGT <- mapply(function(x,y){rep(x,y)},extract.gt(vcffile,element="GT",convertNA=FALSE,IDtoRowNames=FALSE),variantaltnum,USE.NAMES=FALSE)
  variantGQ <- mapply(function(x,y){rep(x,y)},extract.gt(vcffile,element="GQ",as.numeric=TRUE,convertNA=FALSE,IDtoRowNames=FALSE),variantaltnum,USE.NAMES=FALSE)
  
  variantAD <- lapply(strsplit(extract.gt(vcffile,element="AD",convertNA=FALSE,IDtoRowNames=FALSE),","),function(x){as.numeric(x)})
  variantDepth <- lapply(variantAD,sum)
  variantVAF <- mapply(function(x,y){x[-1]/y},variantAD,variantDepth)
  variantDepth <- mapply(function(x,y){rep(x,y)},variantDepth,variantaltnum,USE.NAMES=FALSE)
  
  vcfvariants <- data.frame(Chrom=unlist(variantchrom),Pos=unlist(variantpos),QUAL=unlist(variantqual),Ref=unlist(variantref),Alt=unlist(variantalt),FILTER=unlist(variantFILTER),GTindex=unlist(lapply(variantaltnum,function(x){seq(x)})),GT=unlist(variantGT),VAF=unlist(variantVAF),Depth=unlist(variantDepth),GQ=unlist(variantGQ))
  
  vcfvariants$id <- paste0(vcfvariants$Chrom,"_",vcfvariants$Pos,"_",vcfvariants$Ref,"_",vcfvariants$Alt)
  
  rm(variantalt,variantaltnum,variantchrom,variantpos,variantqual,variantref,variantFILTER,variantGT,variantGQ,variantAD,variantDepth,variantVAF)
  
  #Match bamfile variants to vcf SNV variants using chrom_position_ref_alt.
  #Exclude '*' ALT allele and indels from vcf. '*' alleles indicate overlapping deletions, which are not needed.
  bamvcfvariants <- join(bamvariants,vcfvariants[vcfvariants$Alt!="*" & nchar(vcfvariants$Ref)==1 & nchar(vcfvariants$Alt)==1,],by="id")
  
  #Add the variant information to the main bamfile table, under the 'tag' section.
  for(i in c("QUAL","FILTER","GTindex","GT","VAF","Depth","GQ")){
    bamfile$vcfSNV[[vcffilename]][[i]] <- list()
    bamfile$vcfSNV[[vcffilename]][[i]][[length(bamfile$tag$rn)+1]] <- ""
    bamfile$vcfSNV[[vcffilename]][[i]][[length(bamfile$tag$rn)+1]] <- NULL
    
    apply(bamvcfvariants,1,function(x){
      if(i %in% c("QUAL","GTindex","VAF","Depth","GQ")){
        bamfile$vcfSNV[[vcffilename]][[i]][[as.numeric(x[1])]][as.numeric(x[2])] <<- as.numeric(x[i])
      }
      else{
        bamfile$vcfSNV[[vcffilename]][[i]][[as.numeric(x[1])]][as.numeric(x[2])] <<- x[i]
      }
    })
  }
  
    #Save indels as Granges object. Remove '*' ALT alleles, which indicates overlapping deletions, which are not needed because every deletion already has a separate vcf entry.
    #Also, change ranges for deletion variants to accurately reflect position of deletion. Insertion variants are annotated with width=-1 at the insertion site (insertion is immediately to the right of the first Ref base).
    #Load the result later into bamfilezmw.
    vcfIndels[[vcffilename]] <- subset(vcfvariants[vcfvariants$Alt!="*" & nchar(vcfvariants$Ref)!=nchar(vcfvariants$Alt),],select=-id)
    vcfIndels[[vcffilename]] <- makeGRangesFromDataFrame(vcfIndels[[vcffilename]],keep.extra.columns=TRUE,ignore.strand=TRUE,seqnames.field="Chrom",start.field="Pos",end.field="Pos")
    lengthdiff <- nchar(vcfIndels[[vcffilename]]$Ref) - nchar(vcfIndels[[vcffilename]]$Alt)
    end(vcfIndels[[vcffilename]])[lengthdiff>0] <- start(vcfIndels[[vcffilename]])[lengthdiff>0] + lengthdiff[lengthdiff>0]
    start(vcfIndels[[vcffilename]]) <- start(vcfIndels[[vcffilename]]) + 1
  
  cat("DONE\n")
}
cat("  Annotating bamfile with VCF data...DONE\n")


####
#C. Make a new bamfilezmw dataset that combines fwd and rev ccs reads of each zmw, with ".fwd" and ".rev" added to each column name. The naming of ".fwd" ".rev" for annotations corresponds to the reference genome strand the read aligned to, not the "/fwd" and "/rev" strand designation in the read name of the ccs consensus. For example, a "/fwd" ccs read may be aligned to the reverse strand of the genome (and vice versa), and in this case the read's information will be in the ".rev" annotations per the strand of the genome that the read aligned to.
####
#Extract zmwname and ccsstrand and place in new columns
cat("  Combining reads into ZMWs...")

bamfile$zmwname <- sapply(bamfile$qname,function(x){sub("(/fwd|/rev)$","",x)},USE.NAMES=FALSE)
bamfile$ccsstrand <- sapply(bamfile$qname,function(x){if(grepl("/fwd$",x)){"fwd"}else{"rev"}},USE.NAMES=FALSE)

#Create two datasets, bamfilefwd for strand="+" and bamfilerev for strand="-". Note, the 'strand' value is the direction the read aligned to relative to the genome reference (+ = fwd direction, - = rev direction), and this corresponds to the synthesized strand of the molecule.
bamfilefwd <- list()
bamfilerev <- list()
for(i in setdiff(names(bamfile),c("tag","vcfSNV"))){
  bamfilefwd[[i]] <- bamfile[[i]][bamfile$strand=="+"]
  bamfilerev[[i]] <- bamfile[[i]][bamfile$strand=="-"]
}
for(i in names(bamfile$tag)){
  bamfilefwd$tag[[i]] <- bamfile$tag[[i]][bamfile$strand=="+"]
  bamfilerev$tag[[i]] <- bamfile$tag[[i]][bamfile$strand=="-"]
}
for(i in names(bamfile$vcfSNV)){
  for(j in names(bamfile$vcfSNV[[i]])){
    bamfilefwd$vcfSNV[[i]][[j]] <- bamfile$vcfSNV[[i]][[j]][bamfile$strand=="+"]
    bamfilerev$vcfSNV[[i]][[j]] <- bamfile$vcfSNV[[i]][[j]][bamfile$strand=="-"]
  }
}

#Confirm that every zmw has both fwd and rev ccs reads (should return TRUE)
if(!all(apply(cbind(sort(bamfilefwd$zmwname),sort(bamfilerev$zmwname)),1,function(x){x[1]==x[2]}))){
  stop("ERROR: Not every ZMW has both fwd and rev ccs reads!")
}

#Confirm that there is only one alignment for each zmw (should return FALSE)
if(any(duplicated(bamfilefwd$zmwname)) | any(duplicated(bamfilerev$zmwname))){
  stop("ERROR: Some ZMWs have multiple alignments!")
}

#Create a table that joins the indices of the fwd and rev reads of each zmw
zmwjoin <- merge(data.frame(zmwname=bamfilefwd$zmwname,idx=1:length(bamfilefwd$zmwname)),data.frame(zmwname=bamfilerev$zmwname,idx=1:length(bamfilerev$zmwname)),by="zmwname")
row.names(zmwjoin) <- zmwjoin$zmwname

#Join bamfilefwd and bamfilerev horizontally in one dataset.
bamfilezmw <- list()
for(i in setdiff(names(bamfile),c("tag","vcfSNV"))){
  bamfilezmw[[paste0(i,".fwd")]] <- bamfilefwd[[i]]
  bamfilezmw[[paste0(i,".rev")]] <- bamfilerev[[i]][zmwjoin[bamfilefwd$zmwname,"idx.y"]]
}
for(i in names(bamfilefwd$tag)){
  bamfilezmw$tag.fwd[[i]] <- bamfilefwd$tag[[i]]
  bamfilezmw$tag.rev[[i]] <- bamfilerev$tag[[i]][zmwjoin[bamfilefwd$zmwname,"idx.y"]]
}
for(i in names(bamfilefwd$vcfSNV)){
  for(j in names(bamfilefwd$vcfSNV[[i]])){
    bamfilezmw$vcfSNV.fwd[[i]][[j]] <- bamfilefwd$vcfSNV[[i]][[j]]
    bamfilezmw$vcfSNV.rev[[i]][[j]] <- bamfilerev$vcfSNV[[i]][[j]][zmwjoin[bamfilefwd$zmwname,"idx.y"]]
  }
}

#Confirm zmwname.fwd and zmwname.rev match in bamfilezmw. Should return TRUE.
if(!all(bamfilezmw$zmwname.fwd==bamfilezmw$zmwname.rev)){
  stop("ERROR: Paired ZMW fwd and rev query names do not match!")
}

#Confirm that every zmw has one ccs read aligning to fwd and rev genome reference strands (should return TRUE)
if(!all(bamfilezmw$strand.fwd != bamfilezmw$strand.rev)){
  stop("ERROR: Not every ZMW has ccs reads aligning to both fwd and rev genome strands!")
}

cat("DONE\n")


####
#D. Add columns for each zmw with tag information for positions where fwd and rev ccs reads have the same reference mismatch position and the same ccs read base (i.e. the same mutation in both fwd and rev ccs reads):
#  - zmw.isize: the size in bp of the overlap (intersection) of fwd and reverse reference span.
#  - zmw.rp: Reference mismatch positions
#  - zmw.rn: Reference mismatch base sequences
#  - zmw.qn: Query mismatch base sequences
#  - zmw.fwdqp, zmw.revqp: Query mismatch base positions for shared mutations (from fwd and rev reads); useful for later query position filtering.
#  - zmw.fwdqq, zmw.revqq, zmw.avgqq: Query mismatch base qualities for shared mutations (from fwd read, rev read, and average)
#  - zmw.noN: TRUE for positions that are 'N' in the reference.
#  - zmw.vcfSNV[["vcffile"]]: For each VCF file-- QUAL, FILTER, GTindex, GT, VAF, Depth, GQ.
####

#Calculate zmw.isize (intersection between reference spans of fwd and rev reads)
bamfilezmw$zmw.isize <- apply(cbind(bamfilezmw$pos.fwd+bamfilezmw$isize.fwd,bamfilezmw$pos.rev+bamfilezmw$isize.rev),1,min)-apply(cbind(bamfilezmw$pos.fwd,bamfilezmw$pos.rev),1,max)

#Concatenate reference mismatch positions and query mismatch bases into one string to enable easy matching
rpandqn.fwd <- apply(cbind(bamfilezmw$tag.fwd$rp,bamfilezmw$tag.fwd$qn),1,function(x){  paste0(x[[1]],x[[2]])  },simplify=FALSE)
rpandqn.rev <- apply(cbind(bamfilezmw$tag.rev$rp,bamfilezmw$tag.rev$qn),1,function(x){  paste0(x[[1]],x[[2]])  },simplify=FALSE)

#Make list of all positions in fwd reads that match rev reads for both reference mismatch position and query mismatch base, and vice versa
rpandqn.fwd.match <- apply(cbind(rpandqn.fwd,rpandqn.rev),1,function(x){  x[[1]] %in% x[[2]]  },simplify=FALSE)
rpandqn.rev.match <- apply(cbind(rpandqn.fwd,rpandqn.rev),1,function(x){  x[[2]] %in% x[[1]]  },simplify=FALSE)

#Add columns with details of matching positions to bamfilezmw.
for(i in c("rp","rn","qn","noN")){
  bamfilezmw[[paste0("zmw.",i)]] <- apply(cbind(bamfilezmw$tag.fwd[[i]],rpandqn.fwd.match),1,function(x){  x[[1]][x[[2]]]  },simplify=FALSE)
}

#Annotate zmw.fwdqp and zmw.revqp.
bamfilezmw[["zmw.fwdqp"]] <- apply(cbind(bamfilezmw$tag.fwd$qp,rpandqn.fwd.match),1,function(x){x[[1]][x[[2]]]},simplify=FALSE)
bamfilezmw[["zmw.revqp"]] <- apply(cbind(bamfilezmw$tag.rev$qp,rpandqn.rev.match),1,function(x){x[[1]][x[[2]]]},simplify=FALSE)

#Annotate zmw.fwdqq and zmw.revqq, and zmw.avgqq.
bamfilezmw[["zmw.fwdqq"]] <- apply(cbind(bamfilezmw$tag.fwd$qq,rpandqn.fwd.match),1,function(x){x[[1]][x[[2]]]},simplify=FALSE)
bamfilezmw[["zmw.revqq"]] <- apply(cbind(bamfilezmw$tag.rev$qq,rpandqn.rev.match),1,function(x){x[[1]][x[[2]]]},simplify=FALSE)
bamfilezmw[["zmw.avgqq"]] <- apply(
  cbind(bamfilezmw[["zmw.fwdqq"]],bamfilezmw[["zmw.revqq"]]),1,
  function(x){rowMeans(cbind(x[1][[1]],x[2][[1]]))}
,simplify=FALSE)

#Annotate zmw.vcfSNV.
for(i in names(bamfilezmw$vcfSNV.fwd)){
  for(j in names(bamfilezmw$vcfSNV.fwd[[i]])){
    bamfilezmw[["zmw.vcfSNV"]][[i]][[j]] <- apply(cbind(bamfilezmw$vcfSNV.fwd[[i]][[j]],rpandqn.fwd.match),1,function(x){  x[[1]][x[[2]]]  },simplify=FALSE)
  }
}

#Add vcfIndels annotation
bamfilezmw$vcfIndels <- vcfIndels

cat("DONE\n")
return(bamfilezmw)
}


######################
### Load configuration
######################
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=1) {
  stop("Usage: Rscript make_bamfilezmw_all.R [configuration.yaml]\n", call.=FALSE)
}
cat("#### Loading configuration...")
yaml.config <- suppressWarnings(read.config(args[1]))

bamfilezmw_samplenames <- as.character(unlist(yaml.config$samplenames))

bamfilezmw_genomeversions <- list()
bamfile_filenames <- list()
vcffile_filenames <- list()

for(i in names(yaml.config$make_bamfilezmw_all_config)){
	sampleindex <- match(i,bamfilezmw_samplenames)
  bamfilezmw_genomeversions[[i]] <- yaml.config$genome
  bamfile_filenames[[i]] <- paste0(yaml.config$process_subreads_output_path,"/",yaml.config$ccs_BAM_prefix,".",bamfilezmw_samplenames[sampleindex],".ccs.demux.",sub(":.*","",yaml.config$barcodes[sampleindex]),".postccsfilter.aligned.final.bam")
  vcffile_filenames[[i]] <- yaml.config$make_bamfilezmw_all_config[[i]]$vcffile_filenames
}

chrs_to_analyze <- read.table(yaml.config$chrs_to_analyze)[,1]

bamfilezmw_all_outputpath <- yaml.config$bamfilezmw_all_filename

cat("DONE\n\n")

######################
### Run function to load and annotate bamfiles, and then save bamfilezmw_all into RDS file
######################
#Duplicate output to text file
sink(paste0(dirname(bamfilezmw_all_outputpath),"/make_bamfilezmw_all.output.",format(Sys.time(),"%Y-%m-%d_%H%M%OS2"),".txt"),split=TRUE)

bamfilezmw_all <- list()
for(i in bamfilezmw_samplenames){
  cat("### Processing",i,"-",bamfilezmw_genomeversions[[i]],"###\n")
  bamfilezmw_all[[i]][[bamfilezmw_genomeversions[[i]]]] <- load_bamfile(bamfile_filenames[[i]],vcffile_filenames[[i]],chrs_to_analyze,keepseqdata=FALSE)
}
cat("Saving bamfilezmw_all RDS file...")
qsave(bamfilezmw_all,bamfilezmw_all_outputpath)
cat("DONE\n")

#Close output stream
sink()
