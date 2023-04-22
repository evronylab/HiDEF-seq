#!/usr/bin/env Rscript
# Usage: mutation_frequencies.R [configuration.yaml]
#  - Create [configuration.yaml] per documentation
#
# Prerequisites
#   - R packages: GenomicAlignments,GenomicRanges, plyranges, MutationalPatterns, configr, magrittr, readr, digest, rtracklayer, qs, BSgenome, BSgenome package specified in the YAML configuration
#   - samtools
#   - kmc (https://github.com/refresh-bio/KMC)
#   - seqkit (https://bioinf.shenwei.me/seqkit/)
#		- print_mutations script must have been previously run
#
##Goal: Loop through every sample/genome/filterindex in configuration.yaml file to: 1. Calculate mutation frequencies; 2. Process mutations to save information (for ssDNA mutations, it saves template strand) necessary for mutational signature analysis formatted for MutationalPatterns, SigFit, signature.tools.lib, and general trinucleotide counts for any other signature analysis package.
#
# Saves:
#   1. mutation_frequencies.RDS in analysisoutput_path defined by configuration.yaml containing a list with:
# 		- tissues: List of tissues for samples, in same order as stored in lists
#			- mutations[[samplename]]: GRanges of mutations
#					- fwd
#					- rev
#					- dsDNA
#     - trinucleotide_contexts: list of sample names (sampleid) containing information for fwd/rev/dsDNA trinucleotide context analysis (for fwd/rev/fwdrev: trinucleotide contexts are without collapsing, i.e. all 64 possibilities; for dsDNA: trinucleotide contexts are collapsed to 32 with central pyrmidine)
#       - fwd
#       - rev
#				- fwdrev (combined fwd and rev strands)
#       - dsDNA
#       	=> For each:
#       	  genomefreq: trinucleotide frequencies of full genome
#       	  filteredgenomefreq: trinucleotide frequencies of genome excluding filtered regions
#       	  mutationreadsfreq: trinucleotide frequencies of interrogated portions of reads
#       	  mutationcounts: raw trinucleotide counts of mutations
#       	  mutationcounts.genomecorrected: trinucleotide mutation counts for ratio of trinucleotide context distributions of full genome to interrogated bases in the reads.
#       	  mutationcounts.filteredgenomecorrected: trinucleotide mutation counts for ratio of trinucleotide context distributions of interrogated bases in the genome to interrogated bases in the reads.
#			- signature.tools.lib: Contains data frames for all samples formatted for signature.tools.lib package. These ARE NOT yet collapsed into 96 mutation contexts. Each sample is in a column, with column name = sampleid.genome.filterid.
#				- ssDNA
#				- dsDNA (mutation contexts are based on the reference plus strand)
#			- signature.tools.libcollapsed: Contains data frames for all samples formatted for signature.tools.lib package. These ARE collapsed into 96 mutation contexts, so they are actually only useful for dsDNA analysis. Each sample is in a column, with column name = sampleid.genome.filterid.
#				- ssDNA
#				- dsDNA
# 		- SigFit: Contains data frames for all samples formatted for sigfit package. Note: These ARE NOT yet collapsed into 96 mutation contexts.
# 			- ssDNA
# 			- dsDNA (mutation contexts are based on the reference plus strand)
# 		- SigFitcollapsed: Contains data frames for all samples formatted for sigfit package. Note: These ARE collapsed into 96 mutation contexts.
# 			- ssDNA
# 			- dsDNA
# 		- MutationalPatterns[[samplename]]: Granges formatted for MutationalPatterns package. These are ARE NOT collapsed into 96 mutation contexts. Note that mutation names, REF, and ALT, are all relative to the reference PLUS strand and the synthesized strand, and the GS column contains strand information. Therefore, to get the template strand, need to reverse complement when GS = "+", and keep the same when GS = "-".
# 			- ssDNA
# 			- dsDNA
		
#   2. mutation_frequencies.output.[timestamp].txt log of script that includes all mutation frequencies and other output, but does not include conf. intervals (which are only output into the 'table' format output); in the analysisoutput_path directory.
#
#   3. mutation_frequencies.table.output.[timestamp].txt recording calculated mutation frequencies (including Poisson 95% confidence intervals for frequencies and haploid/diploid counts) and other output pipeline as a table; in the analysisoutput_path directory.

#Stop script for any warnings
options(warn=2)

######
#Load packages
######
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(MutationalPatterns))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(digest))

######
#Define custom functions
######

##Function to set the seed for random sampling (used during germline sensitivity if > max variants used) so that for each sampleid/genome/filterindex combination on any given machine (with the same .Machine$integer.max value), the analysis is reproducible.
set.seed.string <- function(x) {
	hexval <- paste0("0x",digest(x,"crc32"))
	intval <- type.convert(hexval,as.is=TRUE) %% .Machine$integer.max
	set.seed(intval)
}

#Function to merge 2 numeric lists, with correct formatting
merge2way <- function(x,y){
  xy <- merge(x,y,by=0)
  row.names(xy) <- xy[,1]
  return(xy[,c(2,3)])
}

#Function to merge 3 numeric lists, with correct formatting
merge3way <- function(x,y,z){
  xy <- merge(x,y,by=0)
  yz <- merge(y,z,by=0)
  xyz <- merge(xy,yz,by="Row.names")
  row.names(xyz) <- xyz[,1]
  return(xyz[,c(2,3,5)])
}

#Function to subtract two granges from each other, without reducing overlaps in the final output. Modified from code written by Herve Pages.
GRanges_subtract <- function(gr1, gr2, ignore.strand=FALSE){
    gr2 <- reduce(gr2, ignore.strand=ignore.strand)
    hits <- findOverlaps(gr1, gr2, ignore.strand=ignore.strand)
    ans <- psetdiff(gr1, extractList(gr2, as(hits, "IntegerList")))
    unlisted_ans <- unlist(ans, use.names=FALSE)
    mcols(unlisted_ans) <- extractROWS(mcols(gr1),Rle(seq_along(ans), lengths(ans)))
    unlist(setNames(relist(unlisted_ans, ans), names(gr1)))
}

#Function to subtract one granges from another (x and y), i.e. setdiff operation, using a different field for each (xfield and yfield) as the seqnames instead of the current seqnames. For example, use readnames, instead of chromosome names. Performs the operation after removing "/fwd" and "/rev" from xfield and yfield values, which will allow zmw (dsDNA) filtering using qq, bpends, and ccsindel filters derived from both /fwd and /rev reads (because if this filter is present in either fwd or rev reads, it should filter out that position from dsDNA interrogated bases). Removing "/fwd" and "/rev" for fwd and rev regionalfilter doesn't matter, because they are only applied to fwd and rev strands, respectively, so they will only be applied to the reads from which they derive.
setdiffgranges_by_otherfields <- function(x,y,xfield,yfield){
	if(length(y==0)){
		#filter granges is empty
		return(x)
	}else{
	  mcols(x)[xfield] <- sub("/fwd$|/rev$","",as.data.frame(mcols(x)[xfield])[,1])
	  mcols(y)[yfield] <- sub("/fwd$|/rev$","",as.data.frame(mcols(y)[yfield])[,1])
	  x <- split(x,mcols(x)[xfield])
	  y <- split(y,mcols(y)[yfield])
	  x <- x[sort(names(x))]
	  y <- y[sort(names(y))]
	  
	  readsinboth <- intersect(names(x),names(y))
	  readsonlyinx <- setdiff(names(x),names(y))
	  
	  xy.both <- setdiff(x[readsinboth],y[readsinboth])
	  rm(y)
	  xy.both <- c(xy.both,x[readsonlyinx])
	  rm(x)
	  
	  xy.both <- unlist(xy.both)
	  mcols(xy.both)[xfield] <- names(xy.both)
	  names(xy.both) <- NULL
	  
	  return(xy.both)
	}
}

#Function to rapidly calculate trinucleotide counts or distributions (all 64 possible trinucleotides). Takes into account strand; i.e. for reverse strand ranges, it returns the reverse complement using seqkit faidx that returns the reverse complement if start < end coordinate.
# Inputs are:
#   x: granges object to analyze
#   reffasta: path of genome reference fasta file
#   kmcdir: path of directory containing kmc binaries
#   seqkitbin: path of seqkit binary
#   as.prob: return absolute counts (as.prob=FALSE) or distribution fractions (as.prob=TRUE)
trinucleotideFrequency_kmc <- function(x,reffasta,kmcdir,seqkitbin,as.prob=FALSE){
  trinucleotides_64 <- apply(expand.grid(c("A","C","G","T"),c("A","C","G","T"),c("A","C","G","T")),1,paste,collapse="")
  
  tmpregions.tsv <- tempfile()
  tmpregions <- tempfile()
  tmpfasta <- tempfile()
  tmpkmcout <- tempfile()
  tmpkmcdumpout <- tempfile()
  
  x <- as.data.frame(x)
  minusstrandranges <- x[,"strand"]=="-"
  x[minusstrandranges,c("start","end")] <- x[minusstrandranges,c("end","start")]
  
  write_tsv(x[,c("seqnames","start","end")],tmpregions.tsv,col_names=F,progress=F)
  invisible(system(paste("/bin/bash -c",shQuote(paste("awk '{print $1 \":\" $2 \"-\" $3}'",tmpregions.tsv,">",tmpregions))),intern=TRUE))
  invisible(system(paste(seqkitbin,"faidx","-l",tmpregions,reffasta,">",tmpfasta),intern=TRUE,ignore.stderr=TRUE))
  
  invisible(system(paste(paste0(kmcdir,"/kmc"),"-k3 -b -m30 -fm -ci0 -cs9999999999999",tmpfasta,tmpkmcout,tempdir()),intern=TRUE,ignore.stdout=TRUE,ignore.stderr=TRUE))
  invisible(system(paste(paste0(kmcdir,"/kmc_dump"),tmpkmcout,tmpkmcdumpout),intern=TRUE))
  kmcoutput <- read.delim(tmpkmcdumpout,header=FALSE,sep="\t")
  
  file.remove(tmpregions.tsv,tmpregions,tmpfasta,list.files(tempdir(),paste0(tmpkmcout,"*")),tmpkmcdumpout)
  
  result <- rep(0,length(trinucleotides_64))
  names(result) <- trinucleotides_64
  result[kmcoutput$V1] <- kmcoutput$V2
  
  if(as.prob){
    return(result/sum(result))
  }else{
    return(result)
  }
}

#Function to reduce 64 to 32 trinucleotide frequency with central pyrimidine. Input is integer array with named elements that results from the trinucleotideFrequency function of Biostrings.
trinucleotide64to32 <- function(x){
  suppressPackageStartupMessages(library(Biostrings))
  trinucleotides_64 <- apply(expand.grid(c("A","C","G","T"),c("A","C","G","T"),c("A","C","G","T")),1,paste,collapse="")
  trinucleotides_32_pyr <- apply(expand.grid(c("A","C","G","T"),c("C","T"),c("A","C","G","T")),1,paste,collapse="")
  trinucleotides_32_pur <- setdiff(trinucleotides_64,trinucleotides_32_pyr)
  
  y <- x[trinucleotides_32_pur]
  x <- x[trinucleotides_32_pyr]
  names(y) <- reverseComplement(DNAStringSet(names(y)))
  result <- merge(x,y,by=0)
  row.names(result) <- result[,1]
  return(apply(result[,-1],1,sum))
}

######
#Begin program
######
ssDNAsamplenames <- c()
dsDNAsamplenames <- c()
ssDNAvcffiles <- c()
dsDNAvcffiles <- c()
trinucleotide_contexts <- list()
mutations <- list()
output.df <- c()

#Lists of all possible 64 trinucleotides.
trinucleotides_64 <- apply(expand.grid(c("A","C","G","T"),c("A","C","G","T"),c("A","C","G","T")),1,paste,collapse="")

######
#Load configuration
######
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=1) {
  stop("Usage: mutation_frequencies.R [configuration.yaml]\n", call.=FALSE)
}
yaml.config <- suppressWarnings(read.config(args[1]))

#Save log output to file
timestamp <- format(Sys.time(),"%Y-%m-%d_%H%M%OS2")
if(!dir.exists(yaml.config$analysisoutput_path)){dir.create(yaml.config$analysisoutput_path,recursive=TRUE)}
sink(paste0(yaml.config$analysisoutput_path,"/mutation_frequencies.output.",timestamp,".txt"),split=TRUE)

cat("####### Mutation Frequencies Analysis #######\n\n")

cat("#### Loading configuration...")
tissues <- yaml.config$tissues
genomes <- yaml.config$BSgenomepackagename
chrs_to_analyze <- read.table(yaml.config$chrs_to_analyze)[,1]

#Load BSgenome package
suppressPackageStartupMessages(library(yaml.config$BSgenomepackagename,character.only=TRUE))
cat("DONE\n\n")

cat("#### Loading data files...")

bamfilezmw_all <- qread(yaml.config$bamfilezmw_all_filename)
#To preserve memory, clear elements in bamfilezmw_all that are not needed for the analysis.
for(i in names(bamfilezmw_all)){
  for(j in names(bamfilezmw_all[[i]])){
  	if(class(bamfilezmw_all[[i]][[j]]) != "list"){
  		#Skip this sample
  		next
  	}
  	
    elementstoremove <- setdiff(names(bamfilezmw_all[[i]][[j]]),grep(paste0("^",c("rname","pos","isize","qname","zmwname","tag","zmw.rp","zmw.rn","zmw.qn","zmw.vcfSNV"),collapse="|"),names(bamfilezmw_all[[i]][[j]]),value=TRUE))
    for(k in elementstoremove){
      bamfilezmw_all[[i]][[j]][[k]] <- NULL
    }
  }
}

bamfilezmw_filtered <- qread(paste0(yaml.config$analysisoutput_path,"/",yaml.config$analysisoutput_basename,".RDS"))
#To preserve memory, clear elements in bamfilezmw_all that are not needed for the analysis.
for(i in names(bamfilezmw_filtered)){
  for(j in names(bamfilezmw_filtered[[i]])){
    for(k in 1:length(bamfilezmw_filtered[[i]][[j]])){
    	if(length(bamfilezmw_filtered[[i]][[j]][[k]])==0){
    		next
    	}
      elementstoremove <- setdiff(names(bamfilezmw_filtered[[i]][[j]][[k]]),grep(paste0("^",c("includezmws","includemutations.fwd","includemutations.rev","includemutations.zmw","readregionalfilters.fwd","readregionalfilters.rev","readregionalfilters.zmw","genomeregionalfilters.fwd","genomeregionalfilters.rev","genomeregionalfilters.zmw"),collapse="|"),names(bamfilezmw_filtered[[i]][[j]][[k]]),value=TRUE))
      for(l in elementstoremove){
        bamfilezmw_filtered[[i]][[j]][[k]][[l]] <- NULL
      }
    }
  }
}

cat("DONE\n\n")

#Calculate trinucleotide context counts and distributions of full genome fwd and rev strands, excluding N sequences, both 64 (all) and 32 (central pyrimidine) contexts. Note 32 (central pyrimidine) trinucleotide context is identical for both fwd and rev strands.
cat("#### Calculating genome trinucleotide context distribution...")
genome.granges <- as.data.frame(read.table(yaml.config$chrsizes))
genome.granges$V3 <- rep(1,nrow(genome.granges))
genome.granges <- makeGRangesFromDataFrame(genome.granges,seqnames.field="V1",start.field="V3",end.field="V2")
genome.granges <- genome.granges[seqnames(genome.granges) %in% chrs_to_analyze]

Nref.granges <- import(yaml.config$Nrefbed,format="BED")
Nref.granges <- Nref.granges[seqnames(Nref.granges) %in% chrs_to_analyze]
genome.noN.granges <- GenomicRanges::setdiff(genome.granges,Nref.granges)

genome.noN.trinucleotide_counts.fwd64 <- trinucleotideFrequency(getSeq(eval(parse(text=yaml.config$BSgenomepackagename)),genome.noN.granges),simplify.as="collapsed")[trinucleotides_64]

genome.noN.trinucleotide_counts.32 <- trinucleotide64to32(genome.noN.trinucleotide_counts.fwd64)

genome.noN.trinucleotide_counts.rev64 <- genome.noN.trinucleotide_counts.fwd64
names(genome.noN.trinucleotide_counts.rev64) <- reverseComplement(DNAStringSet(names(genome.noN.trinucleotide_counts.rev64)))
genome.noN.trinucleotide_counts.rev64 <- genome.noN.trinucleotide_counts.rev64[trinucleotides_64]

#Now swap fwd and rev genome trinucleotide counts, for consistence with interrogated base contexts and call contexts that reflect the template strand (i.e., - strand for fwd calls/reads and + strand for rev calls/reads).
genome.noN.trinucleotide_counts.originalrev64 <- genome.noN.trinucleotide_counts.rev64
genome.noN.trinucleotide_counts.rev64 <- genome.noN.trinucleotide_counts.fwd64
genome.noN.trinucleotide_counts.fwd64 <- genome.noN.trinucleotide_counts.originalrev64
rm(genome.noN.trinucleotide_counts.originalrev64)

#Calculate trinucleotide counts for both genome strands
genome.noN.trinucleotide_counts.fwdrev64 <- genome.noN.trinucleotide_counts.fwd64 + genome.noN.trinucleotide_counts.rev64

genome.noN.trinucleotide_freq.fwd64 <- genome.noN.trinucleotide_counts.fwd64/sum(genome.noN.trinucleotide_counts.fwd64)
genome.noN.trinucleotide_freq.rev64 <- genome.noN.trinucleotide_counts.rev64/sum(genome.noN.trinucleotide_counts.rev64)
genome.noN.trinucleotide_freq.fwdrev64 <- genome.noN.trinucleotide_counts.fwdrev64/sum(genome.noN.trinucleotide_counts.fwdrev64)
genome.noN.trinucleotide_freq.32 <- genome.noN.trinucleotide_counts.32/sum(genome.noN.trinucleotide_counts.32)
  
cat("DONE\n\n")
  
for(sampleid in names(bamfilezmw_all)){
  for(genome in names(bamfilezmw_all[[sampleid]])){
  	
    if(class(bamfilezmw_all[[sampleid]][[genome]]) != "list"){
      cat("     No reads or mutations in sample",paste0(sampleid,".",genome),"in selected chromosomes. Proceeding to next sample...\n")
    	#Skip this sample
    	next
  	}
    
    #fwd/rev/zmw.granges are Granges objects of all reads, which is used to calculate final number of interrogated bases for mutation frequency calculation. These are the same for all filter indexes.
    fwd.granges <- makeGRangesFromDataFrame(data.frame(chrom=bamfilezmw_all[[sampleid]][[genome]]$rname.fwd,start=bamfilezmw_all[[sampleid]][[genome]]$pos.fwd,end=bamfilezmw_all[[sampleid]][[genome]]$pos.fwd+bamfilezmw_all[[sampleid]][[genome]]$isize.fwd-1,qname=bamfilezmw_all[[sampleid]][[genome]]$qname.fwd),keep.extra.columns=TRUE)
    rev.granges <- makeGRangesFromDataFrame(data.frame(chrom=bamfilezmw_all[[sampleid]][[genome]]$rname.rev,start=bamfilezmw_all[[sampleid]][[genome]]$pos.rev,end=bamfilezmw_all[[sampleid]][[genome]]$pos.rev+bamfilezmw_all[[sampleid]][[genome]]$isize.rev-1,qname=bamfilezmw_all[[sampleid]][[genome]]$qname.rev),keep.extra.columns=TRUE)
    zmw.granges <- makeGRangesFromDataFrame(data.frame(chrom=seqnames(fwd.granges),start=apply(cbind(start(fwd.granges),start(rev.granges)),1,max),end=apply(cbind(end(fwd.granges),end(rev.granges)),1,min),qname=bamfilezmw_all[[sampleid]][[genome]]$zmwname.fwd),keep.extra.columns=TRUE)
    
    for(filterindex in 1:length(bamfilezmw_filtered[[sampleid]][[genome]])){
      samplename <- paste0(sampleid,".",genome,".",filterindex)
      
      #Set seed to make analysis reproducible on any given machine for any sampeid/genome/filterindex combination, due to random sampling steps for germline sensitivity calculation.
      set.seed.string(samplename)
      
      cat("#### Analyzing sample:",sampleid,"-",tissues[grep(sampleid,names(bamfilezmw_filtered))],"-",genome,"- filter index",filterindex,"...\n")

      if(length(bamfilezmw_filtered[[sampleid]][[genome]][[filterindex]])==0){
    		cat("\n    No ZMWs remain in filtered sample. Skipping to next sample/filter index!\n\n")
      	next
    	}
            
      ##Output filter settings
      cat("    => Filter settings\n")
      
      cat("\n")
      cat("       == Genome reference ==\n")
      for(yamlname in c("genome","chrs_to_analyze","gnomad_sensitivity_ref")){
      	cat("       ",paste0(yamlname,": ",yaml.config[[yamlname]]),"\n",sep="")
      }
      
      cat("\n")
      cat("       == Basic thresholds ==\n")
      print(unlist(yaml.config$thresholds[[filterindex]]))

      cat("\n")
      cat("       == VCF filters ==\n")
      if(yaml.config$vcffilters_sameforallbasicfiltersets){
        for(configindex in seq_along(yaml.config$vcffilters[[1]][[sampleid]])){
          cat(paste0("       ",names(yaml.config$vcffilters[[1]][[sampleid]][[configindex]]),": ",unlist(yaml.config$vcffilters[[1]][[sampleid]][[configindex]])),sep="\n")
          cat("\n")
        }
      }else{
        for(configindex in seq_along(yaml.config$vcffilters[[filterindex]][[sampleid]])){
          cat(paste0("       ",names(yaml.config$vcffilters[[filterindex]][[sampleid]][[configindex]]),": ",unlist(yaml.config$vcffilters[[filterindex]][[sampleid]][[configindex]])),sep="\n")
          cat("\n")
        }
      }
      
      cat("\n")
      cat("       == Bigwig filters ==\n")
      if(yaml.config$bigwigfilters_sameforallbasicfiltersets){
        for(configindex in seq_along(yaml.config$bigwigfilters[[1]])){
          cat(paste0("       ",names(yaml.config$bigwigfilters[[1]][[configindex]]),": ",unlist(yaml.config$bigwigfilters[[1]][[configindex]])),sep="\n")
          cat("\n")
        }
      }else{
        for(configindex in seq_along(yaml.config$bigwigfilters[[filterindex]])){
          cat(paste0("       ",names(yaml.config$bigwigfilters[[filterindex]][[configindex]]),": ",unlist(yaml.config$bigwigfilters[[filterindex]][[configindex]])),sep="\n")
          cat("\n")
        }
      }
      
      cat("\n")
      cat("       == Germline reference filters ==\n")
      if(yaml.config$bamfilters_sameforallbasicfiltersets){
        cat(paste0("       ",names(yaml.config$bamfilters[[1]][[sampleid]]),": ",unlist(yaml.config$bamfilters[[1]][[sampleid]])),sep="\n")
        cat("\n")
      }else{
        cat(paste0("       ",names(yaml.config$bamfilters[[filterindex]][[sampleid]]),": ",unlist(yaml.config$bamfilters[[filterindex]][[sampleid]])),sep="\n")
        cat("\n")
      }
      
      cat("\n")
      cat("       == Mutation rate calculation filters ==\n")
      cat("       Maximum mutations per single strand to include molecule in analysis:",yaml.config$mutratefilters$maxmutationsperssdna,"\n")
      cat("       Maximum mutations per ZMW to include molecule in analysis:",yaml.config$mutratefilters$maxmutationsperzmw,"\n")
      
      cat("\n\n    => Performing filtering for interrogated bases...")
      
      ##Get mutation filtering information
      includezmws <- bamfilezmw_filtered[[sampleid]][[genome]][[filterindex]]$includezmws
      includemutations.fwd <- bamfilezmw_filtered[[sampleid]][[genome]][[filterindex]]$includemutations.fwd
      includemutations.rev <- bamfilezmw_filtered[[sampleid]][[genome]][[filterindex]]$includemutations.rev
      includemutations.zmw <- bamfilezmw_filtered[[sampleid]][[genome]][[filterindex]]$includemutations.zmw
      readregionalfilters.fwd <- bamfilezmw_filtered[[sampleid]][[genome]][[filterindex]]$readregionalfilters.fwd
      readregionalfilters.rev <- bamfilezmw_filtered[[sampleid]][[genome]][[filterindex]]$readregionalfilters.rev
      readregionalfilters.zmw <- bamfilezmw_filtered[[sampleid]][[genome]][[filterindex]]$readregionalfilters.zmw
      genomeregionalfilters.fwd <- bamfilezmw_filtered[[sampleid]][[genome]][[filterindex]]$genomeregionalfilters.fwd
      genomeregionalfilters.rev <- bamfilezmw_filtered[[sampleid]][[genome]][[filterindex]]$genomeregionalfilters.rev
      genomeregionalfilters.zmw <- bamfilezmw_filtered[[sampleid]][[genome]][[filterindex]]$genomeregionalfilters.zmw
      
      ##Calculate raw mutation frequencies
      
			#Number of mutations per ssDNA and ZMW.
      rp.fwd <- unlist(mapply(function(x,y){length(x[y])},bamfilezmw_all[[sampleid]][[genome]]$tag.fwd$rp[includezmws],includemutations.fwd,SIMPLIFY=FALSE))
      rp.rev <- unlist(mapply(function(x,y){length(x[y])},bamfilezmw_all[[sampleid]][[genome]]$tag.rev$rp[includezmws],includemutations.rev,SIMPLIFY=FALSE))
      rp.zmw <- unlist(mapply(function(x,y){length(x[y])},bamfilezmw_all[[sampleid]][[genome]]$zmw.rp[includezmws],includemutations.zmw,SIMPLIFY=FALSE))
      
      #Total number of mutations in reads after max mutations per ssDNA/ZMW filtering
      nummutations.fwd <- sum(rp.fwd[rp.fwd<=yaml.config$mutratefilters$maxmutationsperssdna])
      nummutations.rev <- sum(rp.rev[rp.rev<=yaml.config$mutratefilters$maxmutationsperssdna])
      nummutations.zmw <- sum(rp.zmw[rp.zmw<=yaml.config$mutratefilters$maxmutationsperzmw])
      
      if(sum(nummutations.fwd,nummutations.rev,nummutations.zmw)==0){
        cat("\n    No mutations remain after filtering. Skipping to next sample/filter index!\n\n")
      	next
      }
      
      #Reads after max mutations per ssDNA/ZMW filtering
      mutationreads.fwd <- fwd.granges[includezmws][rp.fwd<=yaml.config$mutratefilters$maxmutationsperssdna]
      mutationreads.rev <- rev.granges[includezmws][rp.rev<=yaml.config$mutratefilters$maxmutationsperssdna]
      mutationreads.zmw <- zmw.granges[includezmws][rp.zmw<=yaml.config$mutratefilters$maxmutationsperzmw]
      
      #Number of bases before filtering
      rawmutationreadbases.fwd <- sum(width(mutationreads.fwd))
      rawmutationreadbases.rev <- sum(width(mutationreads.rev))
      rawmutationreadbases.zmw <- sum(width(mutationreads.zmw))
      
      #Apply regional filters to granges of reads and genome (excluding N regions). For any regional filter for reads that has the column 'appliestoread' defined (i.e., qq, bpends, and ccsindel filters), apply the filter only to the specified mutation reads, because these are filters for specific reads.
      #Perform filtering one chromosome at a time to avoid memory overload.
      
       #Function to filter readregionalfilters, one chromosome at a time. Also used for sensitivitymutationreads.zmw filtering, so it has an option to ignore the gnomad filter.
      filterreadregions <- function(mutationreads,readregionalfilters,chromstoanalyze,dontfiltergnomad=FALSE,gnomadsensitivityref=NULL,mingnomadsimilarity=0.5){
        mutationreads <- split(mutationreads,seqnames(mutationreads))
        
        for(filterchrom in chromstoanalyze){
          filter <- GRanges()
          for(i in 1:length(readregionalfilters)){
            if(!is.null(readregionalfilters[[i]]$appliestoread)){
              seqlevels(readregionalfilters[[i]]) <- seqlevels(mutationreads[[filterchrom]])
              filter <- c(filter,readregionalfilters[[i]] %>% filter(seqnames==filterchrom))
            }
          }
          if(length(filter)>0){
            mutationreads[[filterchrom]] <- setdiffgranges_by_otherfields(mutationreads[[filterchrom]],filter,"qname","appliestoread")
          }
          
          filter <- GRangesList(compress=FALSE)
          loopindex <- 1
          for(i in 1:length(readregionalfilters)){
            if(is.null(readregionalfilters[[i]]$appliestoread)){
              if(dontfiltergnomad){
                if(length(subsetByOverlaps(gnomadsensitivityref,readregionalfilters[[i]],type="equal"))/length(gnomadsensitivityref)>=mingnomadsimilarity){
                  next
                }
              }
              filter[[loopindex]] <- readregionalfilters[[i]] %>% filter(seqnames==filterchrom)
              seqlevels(filter[[loopindex]]) <- seqlevels(mutationreads[[filterchrom]])
              loopindex <- loopindex+1
            }
          }
          filter <- reduce(unlist(filter))
          mutationreads[[filterchrom]] <- GRanges_subtract(mutationreads[[filterchrom]],filter,ignore.strand=TRUE)
        }
        
        return(unlist(mutationreads,use.names=FALSE))
      }
      
       #Function to filter genomeregionalfilters, one chromosome at a time
      filtergenomeregions <- function(genomeranges,genomeregionalfilters,chromstoanalyze){
        genomeranges <- split(genomeranges,seqnames(genomeranges))
        
        for(filterchrom in chromstoanalyze){
          filter <- GRangesList(compress=FALSE)
          for(i in 1:length(genomeregionalfilters)){
            filter[[i]] <- genomeregionalfilters[[i]] %>% filter(seqnames==filterchrom)
            seqlevels(filter[[i]]) <- seqlevels(genomeranges)
          }
          filter <- reduce(unlist(filter))
          genomeranges[[filterchrom]] <- GenomicRanges::setdiff(genomeranges[[filterchrom]],filter)
        }

        return(unlist(genomeranges,use.names=FALSE))
      }
      
       # Filter readregionalfilteres
      mutationreads.fwd <- filterreadregions(mutationreads.fwd,readregionalfilters.fwd,chrs_to_analyze,dontfiltergnomad=FALSE)
      mutationreads.rev <- filterreadregions(mutationreads.rev,readregionalfilters.rev,chrs_to_analyze,dontfiltergnomad=FALSE)
      mutationreads.zmw <- filterreadregions(mutationreads.zmw,readregionalfilters.zmw,chrs_to_analyze,dontfiltergnomad=FALSE)
      
       # Filter genomeregionalfilteres
      genome.noN.granges.fwdfiltered <- filtergenomeregions(genome.noN.granges,genomeregionalfilters.fwd,chrs_to_analyze)
      genome.noN.granges.revfiltered <- filtergenomeregions(genome.noN.granges,genomeregionalfilters.rev,chrs_to_analyze)
      genome.noN.granges.zmwfiltered <- filtergenomeregions(genome.noN.granges,genomeregionalfilters.zmw,chrs_to_analyze)
      
      rm(readregionalfilters.fwd,readregionalfilters.rev,genomeregionalfilters.fwd,genomeregionalfilters.rev,genomeregionalfilters.zmw)
      
      #Configure .fwd/rev granges with "-" and "+" strands, respectively, because these are used for trinucleotide context correction, and we correct for the template strand trinucleotide context (the strand where the mutation actually occurs), not the synthesized strand, whereas the ".fwd"/".rev" designation refers to the synthesized strand.
      strand(mutationreads.fwd)[] <- "-"
      strand(mutationreads.rev)[] <- "+"
      strand(mutationreads.zmw)[] <- "*"
      strand(genome.noN.granges.fwdfiltered)[] <- "-"
      strand(genome.noN.granges.revfiltered)[] <- "+"
      strand(genome.noN.granges.zmwfiltered)[] <- "*"
      
      #Total number of bases in reads after filtering
      numbases.fwd <- sum(width(mutationreads.fwd))
      numbases.rev <- sum(width(mutationreads.rev))
      numbases.zmw <- sum(width(mutationreads.zmw))
      
      cat("DONE\n\n")
      
      ##Calculate sensitivity for detecting high quality gnomAD heterozygous germline variants among interrogated bases, specifically due to the minZMWsubreadsVariantReads, minZMWsubreadsVAF, and minsubreads_cvg_fraction filters, which are the only filters applied to final interrogated bases. Used for dsDNA mutation frequency calculations, and the square root of this sensitivity is used for ssDNA mutation frequency calculations since for dsDNA it is applied to each strand separately so the ssDNA sensitivity is the square root of the dsDNA sensitivity.
      
      #Get All high-quality heterozygous variants detected across all ZMWs, for each VCF. Or if only analyzing chrM, all high quality high-VAF variants.
      #Require a minimum number of variants for sensitivity calculation, otherwise set sensitivity to 1.0.
      minsensitivityvariantsrequired <- 100
      minQUALquantile <- 0.5
      minGQquantile <- 0.5
      minDepthquantile <- 0.5
      if(all(chrs_to_analyze=="chrM")){
        cat("    => Calculating sensitivity for known high VAF chrM variants...")
        minVAF <- 0.9
        maxVAF <- 1.0
        sensitivitygenotype <- "1.1"
      }else{
        cat("    => Calculating sensitivity for known heterozygous germline variants...")
        minVAF <- 0.3
        maxVAF <- 0.7
        sensitivitygenotype <- "0.1"
      }


      hetvariants.gr <- GRangesList()
      for(vcf in names(bamfilezmw_all[[sampleid]][[genome]]$zmw.vcfSNV)){
        hetvariants <- grepl(sensitivitygenotype,unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.vcfSNV[[vcf]]$GT))
        
        QUALquantile <- quantile(unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.vcfSNV[[vcf]]$QUAL),minQUALquantile,na.rm=TRUE)
        QUALvariants <- unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.vcfSNV[[vcf]]$QUAL)>=QUALquantile
        
        GQquantile <- quantile(unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.vcfSNV[[vcf]]$GQ),minGQquantile,na.rm=TRUE)
        GQvariants <- unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.vcfSNV[[vcf]]$GQ)>=GQquantile
        
        Depthquantile <- quantile(unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.vcfSNV[[vcf]]$Depth),minDepthquantile,na.rm=TRUE)
        Depthvariants <- unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.vcfSNV[[vcf]]$Depth)>=Depthquantile
        
        VAFvariants <- unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.vcfSNV[[vcf]]$VAF)>=minVAF & unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.vcfSNV[[vcf]]$VAF)<=maxVAF
        
        selectedvariants <- hetvariants & QUALvariants & GQvariants & Depthvariants & VAFvariants
        
        allchroms <- mapply(function(x,y){rep(x,length(y))},x=bamfilezmw_all[[sampleid]][[genome]]$rname.fwd,y=bamfilezmw_all[[sampleid]][[genome]]$zmw.rp)
        allzmwnames <- mapply(function(x,y){rep(x,length(y))},x=bamfilezmw_all[[sampleid]][[genome]]$zmwname.fwd,y=bamfilezmw_all[[sampleid]][[genome]]$zmw.rp)
        
        hetvariants.gr[[vcf]] <- makeGRangesFromDataFrame(data.frame(chrom=unlist(allchroms)[selectedvariants],start=unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.rp)[selectedvariants],end=unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.rp)[selectedvariants],strand="*",ref=unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.rn)[selectedvariants],alt=unlist(bamfilezmw_all[[sampleid]][[genome]]$zmw.qn)[selectedvariants],zmw=unlist(allzmwnames)[selectedvariants],row.names=NULL),keep.extra.columns=TRUE)
      }

      #Filter to keep variants identified in all VCF files
      hetvariants.gr.join <- hetvariants.gr[[1]]
      if(length(hetvariants.gr)>1){
	      for(vcf in 2:length(hetvariants.gr)){
	      	hetvariants.gr.join <- makeGRangesFromDataFrame(merge(data.frame(hetvariants.gr.join),data.frame(hetvariants.gr[[vcf]])),keep.extra.columns=TRUE)
	      }
      }
      hetvariants.gr <- sort(hetvariants.gr.join)
      rm(hetvariants.gr.join)
      
      #Filter to keep only autosomal variants, unless only analyzing chrM. This simplifies needing to handle males vs females differently, and it has no effect on sensitivity calculations, since the sensitivity is at a single molecule level and irrespective of zygosity.
      if(!all(chrs_to_analyze=="chrM")){
        hetvariants.gr <- hetvariants.gr[seqnames(hetvariants.gr) %in% setdiff(chrs_to_analyze,c("chrX","chrY","chrM"))]
      }
      
      #Filter to keep only variants detected in included ZMWs
      hetvariants.gr <- hetvariants.gr[hetvariants.gr$zmw %in% bamfilezmw_all[[sampleid]][[genome]]$zmwname.fwd[includezmws]]

      #Filter to keep only variants also detected by gnomAD
      gnomad_sensitivity_ref <- import(yaml.config$gnomad_sensitivity_ref)
      hetvariants.gr <- subsetByOverlaps(hetvariants.gr,gnomad_sensitivity_ref,type="equal")
      
      #Filter to keep only variants detected in final interrogated bases of dsDNA (zmw) analysis. Need to re-run regional filters but without any gnomad regional filter if one was used. Otherwise, this will filter out gnomad variants that we want for this sensitivity analysis. If a regional filter is found that removes > 50% of gnomad_sensitivity_ref, then we remove it from the sensitivity variant set filtering process, because it must be a gnomad-based filter.
      sensitivitymutationreads.zmw <- zmw.granges[includezmws][rp.zmw<=yaml.config$mutratefilters$maxmutationsperzmw]
      sensitivitymutationreads.zmw <- filterreadregions(sensitivitymutationreads.zmw,readregionalfilters.zmw,chrs_to_analyze,dontfiltergnomad=TRUE,gnomadsensitivityref=gnomad_sensitivity_ref,mingnomadsimilarity=0.5)
      
      hetvariants.gr <- subsetByOverlaps(hetvariants.gr,sensitivitymutationreads.zmw)
      
      #Filter to keep a maximum of 10,000 variants, since that will produce a reliable sensitivity estimate and it is computationally too intensive to process all variants.
      maxnum.hetvariants.gr <- 10000
      if(length(hetvariants.gr)>maxnum.hetvariants.gr){
        hetvariants.gr <- hetvariants.gr[sort(sample.int(length(hetvariants.gr),maxnum.hetvariants.gr))]
      }
      
      #If less than minsensitivityvariantsrequired number of variants, set sensitivity to 1.0, and user should know that sensitivity correction is applied.
      if(length(hetvariants.gr) < minsensitivityvariantsrequired){
        cat("       ",length(hetvariants.gr),"high-quality germline variants is less than required minimum of",minsensitivityvariantsrequired,"for reliable sensitivity correction. Setting sensitivity to 1.0, and therefore, sensitivity corrected-values output below are not usable.\n\n")
        numsensitivity_variants.output <- 0
      	sensitivity <- 1.0
        sensitivity.output <- sensitivity
      }
      else{
        #Change full zmw names to just hole numbers for later matching to bcftools mpileup output
        hetvariants.gr$zmw <- gsub(".*?/(.*)/.*","\\1",hetvariants.gr$zmw)
        
        #Make list of regions for mpileup
        hetvariants.gr.regionsfile <- tempfile()
        write.table(unique(cbind(as.character(seqnames(hetvariants.gr)),start(hetvariants.gr))),hetvariants.gr.regionsfile,quote=F,row.names=F,col.names=F,sep="\t")
        
        #Extract subreads of ZMW holes and align.
        zmwstempfile <- tempfile()
        write.table(unique(sort(as.numeric(hetvariants.gr$zmw))),zmwstempfile,quote=F,row.names=F,col.names=F)
        
        tmpsubreadsbam <- paste0(tempfile(),".bam")
        invisible(system(paste("/bin/bash -c",shQuote(paste("source",yaml.config$condabase_script,"; conda activate",yaml.config$conda_pbbioconda_env,"; zmwfilter --include",zmwstempfile,yaml.config$subreads_filename,tmpsubreadsbam,"; pbmm2 align --log-level INFO -j 4",yaml.config$genomemmiindex,tmpsubreadsbam,paste0(zmwstempfile,".subreads.aligned.bam"),"--preset SUBREAD --sort 2>/dev/null; rm",tmpsubreadsbam))),intern=T))
        
        #Reformat BAM file with each ZMW hole having a different read group and sample name, so that bcftools mpileup counts reads separately for each ZMW.
        invisible(system(paste("/bin/bash -c",shQuote(paste(yaml.config$samtools_bin,"view -H",paste0(zmwstempfile,".subreads.aligned.bam"),"| grep -v \"^@RG\" >",paste0(zmwstempfile,".subreads.aligned.bam.header"),"; for i in `cat",zmwstempfile,"`; do printf \"@RG\tID:$i\tSM:$i\n\"; done >>",paste0(zmwstempfile,".subreads.aligned.bam.header"),";",yaml.config$samtools_bin,"view",paste0(zmwstempfile,".subreads.aligned.bam"),"| awk '{split($1,zmw,\"/\"); sub(\"RG:Z:[[:alnum:]]+\t\",\"RG:Z:\" zmw[2] \"\t\"); print}' >",paste0(zmwstempfile,".subreads.aligned.bam.sam"),"; cat",paste0(zmwstempfile,".subreads.aligned.bam.header"),paste0(zmwstempfile,".subreads.aligned.bam.sam"),"|",yaml.config$samtools_bin,"view -b - >",paste0(zmwstempfile,".subreads.aligned.rg.bam"),";",yaml.config$samtools_bin,"index",paste0(zmwstempfile,".subreads.aligned.rg.bam")))),intern=T))
        
        #Do mpileup of fwd and rev subreads, each separately
        pileuptempfile <- tempfile()
        invisible(system(paste("/bin/bash -c",shQuote(paste(yaml.config$bcftools_bin,"mpileup --ff 2048 -I -A -B -Q 0 -d 9999999999 -a \"INFO/AD,FORMAT/DP,FORMAT/ADF,FORMAT/ADR\" -f",yaml.config$fastaref,"-R",hetvariants.gr.regionsfile,paste0(zmwstempfile,".subreads.aligned.rg.bam"),"2>/dev/null | grep -v ^\\#\\#",">",pileuptempfile))),intern=T))
        
        #Load VCF and annotate hetvariants.gr with relevant data
        pileupvcf <- makeGRangesFromDataFrame(read.delim(pileuptempfile),seqnames.field="X.CHROM",start.field="POS",end.field="POS",keep.extra.columns=TRUE)
        pileupvcf.qn <- strsplit(pileupvcf$ALT,",")
        
        FORMATsplit <- strsplit(pileupvcf$FORMAT,":")
        ADFindex <- unlist(lapply(FORMATsplit,function(x){grep("ADF",x)}))
        ADRindex <- unlist(lapply(FORMATsplit,function(x){grep("ADR",x)}))
        
        vcfhits <- findOverlaps(hetvariants.gr,pileupvcf,type="equal",select="first")
        
        #Rarely, no subreads will be aligned to site of variant called in CCS (subreads poor quality and don't align well to site of the CCS constructed by those subreads). However, bcftools mpileup will not output a record for these sites. Therefore, we identify these sites, and set their subread coverage values to 0 so they are filtered out, because their subread coverage is indeed 0.
        vcfhits.na <- which(is.na(vcfhits))
        vcfhits <- vcfhits[!is.na(vcfhits)]
          
        zmwcolmatch <- match(paste0("X",hetvariants.gr$zmw[setdiff(seq_along(hetvariants.gr$zmw),vcfhits.na)]),colnames(mcols(pileupvcf)))
        zmwcol <- as.data.frame(mcols(pileupvcf))[cbind(vcfhits,zmwcolmatch)]
        zmwcolsplit <- strsplit(zmwcol,":")
        
        ADFvalues <- lapply(strsplit(mapply(function(x,y){x[y]},x=zmwcolsplit,y=ADFindex[vcfhits]),","),as.numeric)
        ADRvalues <- lapply(strsplit(mapply(function(x,y){x[y]},x=zmwcolsplit,y=ADRindex[vcfhits]),","),as.numeric)
        
        altindex <- mapply(function(x,y){grep(x,y)+1},x=hetvariants.gr$alt[setdiff(seq_along(hetvariants.gr$alt),vcfhits.na)],y=pileupvcf.qn[vcfhits],USE.NAMES=F)
        
        #Function to add zero values at specific vector indices (postoadd), for sites without an mpileup record
        addzerovalues <- function(x,postoadd){
          postoadd <- postoadd-1
          x <- rep(x, (1:2)[(1:length(x) %in% postoadd) + 1])
          x[postoadd + 1:length(postoadd)] <- 0
          return(x)
        }
        
        hetvariants.gr$fwdsubreadsVariantReads <- addzerovalues(mapply(function(x,y){x[y]},x=ADFvalues,y=altindex),vcfhits.na)
        hetvariants.gr$revsubreadsVariantReads <- addzerovalues(mapply(function(x,y){x[y]},x=ADRvalues,y=altindex),vcfhits.na)
        hetvariants.gr$fwdsubreadsVAF <- addzerovalues(mapply(function(x,y){x[y]/sum(x)},x=ADFvalues,y=altindex),vcfhits.na)
        hetvariants.gr$revsubreadsVAF <- addzerovalues(mapply(function(x,y){x[y]/sum(x)},x=ADRvalues,y=altindex),vcfhits.na)
        
        #Set to 0 any variants not found in the pileup due to differences in positions of variants in gnomad versus subread alignments. This should not occur, but doing this check just in case
        for(annotatedcolumn in c("fwdsubreadsVariantReads","revsubreadsVariantReads","fwdsubreadsVAF","revsubreadsVAF")){
            varsnotfound <- unlist(lapply(mcols(hetvariants.gr)[[annotatedcolumn]],length))==0
            mcols(hetvariants.gr)[[annotatedcolumn]][varsnotfound] <- rep(list(0),length(which(varsnotfound==TRUE)))
            
            varsnotfound <- unlist(lapply(mcols(hetvariants.gr)[[annotatedcolumn]],is.na))
            mcols(hetvariants.gr)[[annotatedcolumn]][varsnotfound] <- rep(list(0),length(which(varsnotfound==TRUE)))
            mcols(hetvariants.gr)[[annotatedcolumn]] <- unlist(mcols(hetvariants.gr)[[annotatedcolumn]])
        }
        
        #Annotate hetvariants.gr with subreads coverage fraction for each strand separately, for use in minsubreads_cvg_fraction filtering
        
          #Load subreads BAM file
          subreadsbamfile <- scanBam(paste0(zmwstempfile,".subreads.aligned.bam"),param=ScanBamParam(what=c("rname","pos","isize","strand"),tag=c("zm")))[[1]]
          subreadsbamfile.gr <- makeGRangesFromDataFrame(data.frame(chrom=subreadsbamfile$rname,start=subreadsbamfile$pos,end=subreadsbamfile$pos+subreadsbamfile$isize-1,strand=subreadsbamfile$strand,zm=subreadsbamfile$tag$zm),keep.extra.columns=TRUE)
          fwdsubreadsbamfile.gr <- subreadsbamfile.gr[strand(subreadsbamfile.gr)=="+"]
          revsubreadsbamfile.gr <- subreadsbamfile.gr[strand(subreadsbamfile.gr)=="-"]
          
          #Annotate number of subreads covering each mutation
          annotatesubreadscvg <- function(gr,subreads.gr,annotationcolumn){
            subreadscvg.temp <- join_overlap_inner(gr,subreads.gr) %>% filter(zmw==zm) %>% group_by(seqnames,start,end,zmw) %>% summarise(n=n()) %>% as_granges()
            gr <- join_overlap_left(gr,subreadscvg.temp) %>% filter(zmw.x==zmw.y) %>% select(!zmw.y,zmw=zmw.x,!!annotationcolumn:=n)
            mcols(gr)[[annotationcolumn]][is.na(mcols(gr)[[annotationcolumn]])] <- 0
            return(gr)
          }
          hetvariants.gr <- annotatesubreadscvg(hetvariants.gr,fwdsubreadsbamfile.gr,"fwdsubreadscvg")
          hetvariants.gr <- annotatesubreadscvg(hetvariants.gr,revsubreadsbamfile.gr,"revsubreadscvg")
        
          #Annotate total number of subreads for each mutation's zmw
          fwdsubreadsperzmw <- table(fwdsubreadsbamfile.gr$zm)
          revsubreadsperzmw <- table(revsubreadsbamfile.gr$zm)
          hetvariants.gr$totalfwdsubreads <- as.numeric(fwdsubreadsperzmw[match(hetvariants.gr$zmw,names(fwdsubreadsperzmw))])
          hetvariants.gr$totalrevsubreads <- as.numeric(revsubreadsperzmw[match(hetvariants.gr$zmw,names(revsubreadsperzmw))])
          
          #Calculate subread coverage fraction
          hetvariants.gr$fwdsubreadscvgfraction <- hetvariants.gr$fwdsubreadscvg/hetvariants.gr$totalfwdsubreads
          hetvariants.gr$revsubreadscvgfraction <- hetvariants.gr$revsubreadscvg/hetvariants.gr$totalrevsubreads
          
          #Set to 0 for Na values (i.e. if total reads is 0, leads to division by zero, which leads to NA)
          hetvariants.gr$fwdsubreadscvgfraction[is.na(hetvariants.gr$fwdsubreadscvgfraction)] <- 0
          hetvariants.gr$revsubreadscvgfraction[is.na(hetvariants.gr$revsubreadscvgfraction)] <- 0
          
        #Remove temporary variables
        rm(pileupvcf,pileupvcf.qn,FORMATsplit,ADFindex,ADRindex,ADFvalues,ADRvalues,vcfhits,zmwcolmatch,zmwcol,zmwcolsplit,altindex,varsnotfound,subreadsbamfile,subreadsbamfile.gr,fwdsubreadsbamfile.gr,revsubreadsbamfile.gr,fwdsubreadsperzmw,revsubreadsperzmw)
        
        #Calculate sensitivity for variants detected above the minZMWsubreadsVariantReads, minZMWsubreadsVAF, and minsubreads_cvg_fraction thresholds. Note, both strands must separately pass these thresholds, which is the same way these filters are applied for somatic variants.
        sensitivity <- length(which(hetvariants.gr$fwdsubreadsVAF>=yaml.config$thresholds[[filterindex]]$minZMWsubreadsVAF & hetvariants.gr$revsubreadsVAF>=yaml.config$thresholds[[filterindex]]$minZMWsubreadsVAF & hetvariants.gr$fwdsubreadsVariantReads>=yaml.config$thresholds[[filterindex]]$minZMWsubreadsVariantReads & hetvariants.gr$revsubreadsVariantReads>=yaml.config$thresholds[[filterindex]]$minZMWsubreadsVariantReads & hetvariants.gr$fwdsubreadscvgfraction>=yaml.config$thresholds[[filterindex]]$minsubreads_cvg_fraction & hetvariants.gr$revsubreadscvgfraction>=yaml.config$thresholds[[filterindex]]$minsubreads_cvg_fraction))/length(hetvariants.gr)
        
        numsensitivity_variants.output <- length(hetvariants.gr)
        sensitivity.output <- sensitivity

        #Remove temp files
        invisible(file.remove(pileuptempfile,zmwstempfile,paste0(zmwstempfile,".subreads.aligned.bam"),paste0(zmwstempfile,".subreads.aligned.bam.bai"),paste0(zmwstempfile,".subreads.aligned.bam.header"),paste0(zmwstempfile,".subreads.aligned.bam.sam"),paste0(zmwstempfile,".subreads.aligned.rg.bam"),paste0(zmwstempfile,".subreads.aligned.rg.bam.bai")))
      }
      
      cat("DONE\n\n")
      
      ##Calculate statistics
      cat("    => Calculating output statistics...")
      
      #Output fraction of total ZMWs that were filtered by ZMW filters
      numZMW_filteredbyZMWfilters_beforemaxmutperZMW.output <- length(which(includezmws==FALSE))
	    numZMW.output <- length(zmw.granges)
	    fracZMW_filteredbyZMWfilters_beforemaxmutperZMW.output <- numZMW_filteredbyZMWfilters_beforemaxmutperZMW.output/numZMW.output
	      
      #Output fraction of total ZMWs|strands that were filtered due to number of mutations greater than threshold
      numFwdRevStrand_filtered_bymaxmutperStrand.output <- sum(length(fwd.granges[includezmws])-length(which(rp.fwd<=yaml.config$mutratefilters$maxmutationsperssdna)),length(rev.granges[includezmws])-length(which(rp.rev<=yaml.config$mutratefilters$maxmutationsperssdna)))
      numFwdRevStrand.output <- sum(length(fwd.granges),length(rev.granges))
      fracFwdRevStrand_filtered_bymaxmutperStrand.output <- numFwdRevStrand_filtered_bymaxmutperStrand.output/numFwdRevStrand.output
      numZMW_filtered_bymaxmutperZMW.output <- length(zmw.granges[includezmws])-length(which(rp.zmw<=yaml.config$mutratefilters$maxmutationsperzmw))
      fracZMW_filtered_bymaxmutperZMW.output <- numZMW_filtered_bymaxmutperZMW.output/numZMW.output

      #Output fraction of total sequenced bases that were filtered by the ZMW filters
			numFwdRevStrandBases_filteredbyZMWfilters_beforemaxmutperStrand.output <- sum(sum(width(fwd.granges[!includezmws])),sum(width(rev.granges[!includezmws])))
			numFwdRevStrandBases.output <- sum(sum(width(fwd.granges)),sum(width(rev.granges)))
			fracFwdRevStrandBases_filteredbyZMWfilters_beforemaxmutperStrand.output <- numFwdRevStrandBases_filteredbyZMWfilters_beforemaxmutperStrand.output/numFwdRevStrandBases.output
			numZMWBases_filteredbyZMWfilters_beforemaxmutperZMW.output <- sum(width(zmw.granges[!includezmws]))
			numZMWBases.output <- sum(width(zmw.granges))
			fracZMWBases_filteredbyZMWfilters_beforemaxmutperZMW.output <- numZMWBases_filteredbyZMWfilters_beforemaxmutperZMW.output/numZMWBases.output

			numFwdRevStrandBases_filtered_bymaxmutperStrand.output <- sum(sum(width(fwd.granges[includezmws]))-rawmutationreadbases.fwd,sum(width(rev.granges[includezmws]))-rawmutationreadbases.rev)
			fracFwdRevStrandBases_filtered_bymaxmutperStrand.output <- numFwdRevStrandBases_filtered_bymaxmutperStrand.output/numFwdRevStrandBases.output
			numZMWBases_filtered_bymaxmutperZMW.output <- sum(width(zmw.granges[includezmws]))-rawmutationreadbases.zmw
			fracZMWBases_filtered_bymaxmutperZMW.output <- numZMWBases_filtered_bymaxmutperZMW.output/numZMWBases.output

      #Output fraction of final ZMW bases (i.e. after basic ZMW filters and filters removing ZMWs with more than max number of mutations) that were filtered.
			numFwdRevStrandBases_filtered_afterZMWfiltersandmaxmutperStrand.output <- sum(rawmutationreadbases.fwd-numbases.fwd,rawmutationreadbases.rev-numbases.rev)
			numFwdRevStrandBases_afterZMWfilterandmaxmutperStrand.output <- sum(rawmutationreadbases.fwd,rawmutationreadbases.rev)
			fracFwdRevStrandBases_filtered_afterZMWfiltersandmaxmutperStrand.output <- numFwdRevStrandBases_filtered_afterZMWfiltersandmaxmutperStrand.output/numFwdRevStrandBases_afterZMWfilterandmaxmutperStrand.output
			
			numZMWBases_filtered_afterZMWfiltersandmaxmutperZMW.output <- rawmutationreadbases.zmw-numbases.zmw
			numZMWBases_afterZMWfilterandmaxmutperZMW.output <- rawmutationreadbases.zmw
			fracZMWBases_filtered_afterZMWfiltersandmaxmutperZMW.output <- numZMWBases_filtered_afterZMWfiltersandmaxmutperZMW.output/numZMWBases_afterZMWfilterandmaxmutperZMW.output
			
      #Output efficiency (i.e. fraction of total bases remaining after ZMW-level and regional filters)
			numFwdRevStrandBases_interrogated.output <- sum(numbases.fwd,numbases.rev)
			fracFwdRevStrandBases_interrogated_from_numFwdRevStrandBases.output <- numFwdRevStrandBases_interrogated.output/numFwdRevStrandBases.output
			
			numZMWBases_interrogated.output <- numbases.zmw
			fracZMWBases_interrogated_from_numZMWBases.output <- numZMWBases_interrogated.output/numZMWBases.output
			
			numSubreadBases.output <- sum(bamfilezmw_all[[sampleid]][[genome]]$tag.fwd$ec*bamfilezmw_all[[sampleid]][[genome]]$isize.fwd,bamfilezmw_all[[sampleid]][[genome]]$tag.rev$ec*bamfilezmw_all[[sampleid]][[genome]]$isize.rev)
			fracZMWBases_interrogated_from_numSubreadBases.output <- numZMWBases_interrogated.output/numSubreadBases.output

      #Output raw call burdens
			num_raw_FwdRevmutations.output <- sum(nummutations.fwd,nummutations.rev)
			num_raw_FwdRevmutations_lci.output <- poisson.test(num_raw_FwdRevmutations.output)$conf.int[1]
			num_raw_FwdRevmutations_uci.output <- poisson.test(num_raw_FwdRevmutations.output)$conf.int[2]
			
			#Calculate ratio of raw lci and raw uci FwdRev mutation counts to raw FwdRev mutation count, for calculating lci and uci of corrected FwdRev burdens
			raw_FwdRev_lci_to_calc_ratio <- num_raw_FwdRevmutations_lci.output/num_raw_FwdRevmutations.output
			raw_FwdRev_uci_to_calc_ratio <- num_raw_FwdRevmutations_uci.output/num_raw_FwdRevmutations.output
			
			freq_raw_FwdRevmutations.output <- num_raw_FwdRevmutations.output/numFwdRevStrandBases_interrogated.output
			freq_raw_FwdRevmutations_lci.output <- freq_raw_FwdRevmutations.output*raw_FwdRev_lci_to_calc_ratio
			freq_raw_FwdRevmutations_uci.output <- freq_raw_FwdRevmutations.output*raw_FwdRev_uci_to_calc_ratio
			
			num_sensCorrected_FwdRevmutations.output <- num_raw_FwdRevmutations.output/sqrt(sensitivity)
			freq_sensCorrected_FwdRevmutations.output <- num_sensCorrected_FwdRevmutations.output/numFwdRevStrandBases_interrogated.output
			freq_sensCorrected_FwdRevmutations_lci.output <- freq_sensCorrected_FwdRevmutations.output*raw_FwdRev_lci_to_calc_ratio
			freq_sensCorrected_FwdRevmutations_uci.output <- freq_sensCorrected_FwdRevmutations.output*raw_FwdRev_uci_to_calc_ratio
			
			num_raw_ZMWmutations.output <- nummutations.zmw
			num_raw_ZMWmutations_lci.output <- poisson.test(num_raw_ZMWmutations.output)$conf.int[1]
			num_raw_ZMWmutations_uci.output <- poisson.test(num_raw_ZMWmutations.output)$conf.int[2]
			
			#Calculate ratio of raw lci and raw uci ZMW mutation counts to raw ZMW mutation count, for calculating lci and uci of corrected ZMW burdens
			raw_ZMW_lci_to_calc_ratio <- num_raw_ZMWmutations_lci.output/num_raw_ZMWmutations.output
			raw_ZMW_uci_to_calc_ratio <- num_raw_ZMWmutations_uci.output/num_raw_ZMWmutations.output
			
			freq_raw_ZMWmutations.output <- num_raw_ZMWmutations.output/numZMWBases_interrogated.output
			freq_raw_ZMWmutations_lci.output <- freq_raw_ZMWmutations.output*raw_ZMW_lci_to_calc_ratio
			freq_raw_ZMWmutations_uci.output <- freq_raw_ZMWmutations.output*raw_ZMW_uci_to_calc_ratio
			
			num_sensCorrected_ZMWmutations.output <- num_raw_ZMWmutations.output/sensitivity
			freq_sensCorrected_ZMWmutations.output <- num_sensCorrected_ZMWmutations.output/numZMWBases_interrogated.output
			freq_sensCorrected_ZMWmutations_lci.output <- freq_sensCorrected_ZMWmutations.output*raw_ZMW_lci_to_calc_ratio
			freq_sensCorrected_ZMWmutations_uci.output <- freq_sensCorrected_ZMWmutations.output*raw_ZMW_uci_to_calc_ratio

      ## Create granges objects of mutations
      fwdzmwstoinclude <- rp.fwd>0 & rp.fwd<=yaml.config$mutratefilters$maxmutationsperssdna
      revzmwstoinclude <- rp.rev>0 & rp.rev<=yaml.config$mutratefilters$maxmutationsperssdna
      zmwzmwstoinclude <- rp.zmw>0 & rp.zmw<=yaml.config$mutratefilters$maxmutationsperzmw
      
      #Retrieve info for each mutation
      fwdchroms <- mapply(function(x,y,z){rep(x,length(y[z]))},as.character(bamfilezmw_all[[sampleid]][[genome]]$rname.fwd[includezmws][fwdzmwstoinclude]),bamfilezmw_all[[sampleid]][[genome]]$tag.fwd$rp[includezmws][fwdzmwstoinclude],includemutations.fwd[fwdzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)
      revchroms <- mapply(function(x,y,z){rep(x,length(y[z]))},as.character(bamfilezmw_all[[sampleid]][[genome]]$rname.rev[includezmws][revzmwstoinclude]),bamfilezmw_all[[sampleid]][[genome]]$tag.rev$rp[includezmws][revzmwstoinclude],includemutations.rev[revzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)
      zmwchroms <- mapply(function(x,y,z){rep(x,length(y[z]))},as.character(bamfilezmw_all[[sampleid]][[genome]]$rname.fwd[includezmws][zmwzmwstoinclude]),bamfilezmw_all[[sampleid]][[genome]]$zmw.rp[includezmws][zmwzmwstoinclude],includemutations.zmw[zmwzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)
      
      fwdpos <- mapply(function(x,y){x[y]},bamfilezmw_all[[sampleid]][[genome]]$tag.fwd$rp[includezmws][fwdzmwstoinclude],includemutations.fwd[fwdzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)
      revpos <- mapply(function(x,y){x[y]},bamfilezmw_all[[sampleid]][[genome]]$tag.rev$rp[includezmws][revzmwstoinclude],includemutations.rev[revzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)
      zmwpos <- mapply(function(x,y){x[y]},bamfilezmw_all[[sampleid]][[genome]]$zmw.rp[includezmws][zmwzmwstoinclude],includemutations.zmw[zmwzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)
      
      fwdrn <- mapply(function(x,y){x[y]},bamfilezmw_all[[sampleid]][[genome]]$tag.fwd$rn[includezmws][fwdzmwstoinclude],includemutations.fwd[fwdzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)
      revrn <- mapply(function(x,y){x[y]},bamfilezmw_all[[sampleid]][[genome]]$tag.rev$rn[includezmws][revzmwstoinclude],includemutations.rev[revzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)
      zmwrn <- mapply(function(x,y){x[y]},bamfilezmw_all[[sampleid]][[genome]]$zmw.rn[includezmws][zmwzmwstoinclude],includemutations.zmw[zmwzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)

      fwdqn <- mapply(function(x,y){x[y]},bamfilezmw_all[[sampleid]][[genome]]$tag.fwd$qn[includezmws][fwdzmwstoinclude],includemutations.fwd[fwdzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)
      revqn <- mapply(function(x,y){x[y]},bamfilezmw_all[[sampleid]][[genome]]$tag.rev$qn[includezmws][revzmwstoinclude],includemutations.rev[revzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)
      zmwqn <- mapply(function(x,y){x[y]},bamfilezmw_all[[sampleid]][[genome]]$zmw.qn[includezmws][zmwzmwstoinclude],includemutations.zmw[zmwzmwstoinclude],USE.NAMES=F,SIMPLIFY=FALSE)
      
      #Make granges. Configure the strand field as the strand of the template, not the synthesized strand, because that is the context where the mutation is actually happening, and we want trinucleotide contexts to be corrected relative to that. Therefore, fwdmutations and revmutations are assigned strands "-" and "+", respectively, since 'fwd'/'rev'mutations refers to the synthesized strand.
      fwdmutations <- data.frame(chrom=unlist(fwdchroms),start=unlist(fwdpos),end=unlist(fwdpos),rn=unlist(fwdrn),qn=unlist(fwdqn),strand=rep("-",length(unlist(fwdpos))))
      revmutations <- data.frame(chrom=unlist(revchroms),start=unlist(revpos),end=unlist(revpos),rn=unlist(revrn),qn=unlist(revqn),strand=rep("+",length(unlist(revpos))))
      dsDNAmutations <- data.frame(chrom=unlist(zmwchroms),start=unlist(zmwpos),end=unlist(zmwpos),rn=unlist(zmwrn),qn=unlist(zmwqn),strand=rep("*",length(unlist(zmwpos))))

      if(nrow(fwdmutations)>0){
        fwdmutations.granges <- makeGRangesFromDataFrame(fwdmutations,keep.extra.columns=TRUE,seqinfo=seqnames(get(yaml.config$BSgenomepackagename)))
      }else{
        fwdmutations.granges <- GRanges()
      }
      if(nrow(revmutations)>0){
        revmutations.granges <- makeGRangesFromDataFrame(revmutations,keep.extra.columns=TRUE,seqinfo=seqnames(get(yaml.config$BSgenomepackagename)))
      }else{
        revmutations.granges <- GRanges()
      }
      if(nrow(dsDNAmutations)>0){
        dsDNAmutations.granges <- makeGRangesFromDataFrame(dsDNAmutations,keep.extra.columns=TRUE,seqinfo=seqnames(get(yaml.config$BSgenomepackagename)))
      }else{
        dsDNAmutations.granges <- GRanges()
      }

 
      ## Calculate mutation frequencies corrected for trinucleotide context of analyzed portion of the genome (genomic regions remaining after regional filters).
      
      #A. Calculate trinucleotide context counts and distributions of mutations (64 trinucleotide contexts for fwd and rev ssDNA mutations, 64 trinucleotide contexts for combined fwd/rev mutations, and 32 trinucleotide contexts for dsDNA mutations)
      if(nrow(fwdmutations)>0){
        fwdmutations.trinucleotide_counts64 <- trinucleotideFrequency(getSeq(eval(parse(text=yaml.config$BSgenomepackagename)),resize(fwdmutations.granges,3,fix="center")),simplify.as="collapsed")[trinucleotides_64]
        fwdmutations.trinucleotide_freq64 <- fwdmutations.trinucleotide_counts64/sum(fwdmutations.trinucleotide_counts64)
      }else{
        fwdmutations.trinucleotide_counts64 <- rep(0,length(trinucleotides_64))
        names(fwdmutations.trinucleotide_counts64) <- trinucleotides_64
        fwdmutations.trinucleotide_freq64 <- 0
      }
      
      if(nrow(revmutations)>0){
        revmutations.trinucleotide_counts64 <- trinucleotideFrequency(getSeq(eval(parse(text=yaml.config$BSgenomepackagename)),resize(revmutations.granges,3,fix="center")),simplify.as="collapsed")[trinucleotides_64]
        revmutations.trinucleotide_freq64 <- revmutations.trinucleotide_counts64/sum(revmutations.trinucleotide_counts64)
      }else{
        revmutations.trinucleotide_counts64 <- rep(0,length(trinucleotides_64))
        names(revmutations.trinucleotide_counts64) <- trinucleotides_64
        revmutations.trinucleotide_freq64 <- 0
      }
      
      fwdrevmutations.trinucleotide_counts64 <- fwdmutations.trinucleotide_counts64 + revmutations.trinucleotide_counts64
      
      if(nrow(dsDNAmutations)>0){
        dsDNAmutations.trinucleotide_counts32 <- trinucleotide64to32(trinucleotideFrequency(getSeq(eval(parse(text=yaml.config$BSgenomepackagename)),resize(dsDNAmutations.granges,3,fix="center")),simplify.as="collapsed")[trinucleotides_64])
        dsDNAmutations.trinucleotide_freq32 <- dsDNAmutations.trinucleotide_counts32/sum(dsDNAmutations.trinucleotide_counts32)
      }else{
        dsDNAmutations.trinucleotide_counts32 <- rep(0,length(trinucleotides_64))
        names(dsDNAmutations.trinucleotide_counts32) <- trinucleotides_64
        dsDNAmutations.trinucleotide_counts32 <- trinucleotide64to32(dsDNAmutations.trinucleotide_counts32)
      }
      
      #B. Calculate trinucleotide context distributions of interrogated bases in genome. 64 trinucleotide contexts for fwd and rev strand regions, 64 trinucleotide contexts for combined fwd/rev strand regions, and 32 trinucleotide contexts for dsDNA regions. Note: this ignores regions that are < 3 base pairs.
      genome.noN.fwdfiltered.trinucleotide_counts.fwd64 <- trinucleotideFrequency_kmc(genome.noN.granges.fwdfiltered,yaml.config$fastaref,yaml.config$kmc_bindir,yaml.config$seqkit_bin,as.prob=FALSE)[trinucleotides_64]
      genome.noN.revfiltered.trinucleotide_counts.rev64 <- trinucleotideFrequency_kmc(genome.noN.granges.revfiltered,yaml.config$fastaref,yaml.config$kmc_bindir,yaml.config$seqkit_bin,as.prob=FALSE)[trinucleotides_64]
      genome.noN.fwdrevfiltered.trinucleotide_counts.fwdrev64 <- genome.noN.fwdfiltered.trinucleotide_counts.fwd64 + genome.noN.revfiltered.trinucleotide_counts.rev64
      genome.noN.zmwfiltered.trinucleotide_counts.32 <- trinucleotide64to32(trinucleotideFrequency_kmc(genome.noN.granges.zmwfiltered,yaml.config$fastaref,yaml.config$kmc_bindir,yaml.config$seqkit_bin,as.prob=FALSE)[trinucleotides_64])

      genome.noN.fwdfiltered.trinucleotide_freq.fwd64 <- genome.noN.fwdfiltered.trinucleotide_counts.fwd64/sum(genome.noN.fwdfiltered.trinucleotide_counts.fwd64)
      genome.noN.revfiltered.trinucleotide_freq.rev64 <- genome.noN.revfiltered.trinucleotide_counts.rev64/sum(genome.noN.revfiltered.trinucleotide_counts.rev64)
      genome.noN.fwdrevfiltered.trinucleotide_freq.fwdrev64 <- genome.noN.fwdrevfiltered.trinucleotide_counts.fwdrev64/sum(genome.noN.fwdrevfiltered.trinucleotide_counts.fwdrev64)
      genome.noN.zmwfiltered.trinucleotide_freq.32 <- genome.noN.zmwfiltered.trinucleotide_counts.32/sum(genome.noN.zmwfiltered.trinucleotide_counts.32)
      
      #C. Calculate trinucleotide context counts and distributions of interrogated bases in reads. 64 trinucleotide context for fwd and rev strand mutation reads, 64 trinucleotide context for combined fwd/rev strand mutation reads, and 32 context for dsDNA mutation reads. Note: this ignores regions that are < 3 base pairs.
      mutationreads.fwd.trinucleotide_counts.fwd64 <- trinucleotideFrequency_kmc(mutationreads.fwd,yaml.config$fastaref,yaml.config$kmc_bindir,yaml.config$seqkit_bin,as.prob=FALSE)[trinucleotides_64]
      mutationreads.rev.trinucleotide_counts.rev64 <- trinucleotideFrequency_kmc(mutationreads.rev,yaml.config$fastaref,yaml.config$kmc_bindir,yaml.config$seqkit_bin,as.prob=FALSE)[trinucleotides_64]
      mutationreads.fwdrev.trinucleotide_counts.fwdrev64 <- mutationreads.fwd.trinucleotide_counts.fwd64 + mutationreads.rev.trinucleotide_counts.rev64
      mutationreads.zmw.trinucleotide_counts.32 <- trinucleotide64to32(trinucleotideFrequency_kmc(mutationreads.zmw,yaml.config$fastaref,yaml.config$kmc_bindir,yaml.config$seqkit_bin,as.prob=FALSE)[trinucleotides_64])
      
      mutationreads.fwd.trinucleotide_freq.fwd64 <- mutationreads.fwd.trinucleotide_counts.fwd64/sum(mutationreads.fwd.trinucleotide_counts.fwd64)
      mutationreads.rev.trinucleotide_freq.rev64 <- mutationreads.rev.trinucleotide_counts.rev64/sum(mutationreads.rev.trinucleotide_counts.rev64)
      mutationreads.fwdrev.trinucleotide_freq.fwdrev64 <- mutationreads.fwdrev.trinucleotide_counts.fwdrev64/sum(mutationreads.fwdrev.trinucleotide_counts.fwdrev64)
      mutationreads.zmw.trinucleotide_freq.32 <- mutationreads.zmw.trinucleotide_counts.32/sum(mutationreads.zmw.trinucleotide_counts.32)

      #D. Correct fwd, rev, fwd/rev combined, and zmw trinucleotide mutation counts for ratio of trinucleotide context distributions of interrogated bases in the genome to interrogated bases in the reads. For any trinucleotide context for which the fraction of the trinucleotide in the interrogated bases in the reads is 0, the final corrected mutation count will result in NaN, because when dividing the interrogated bases in the genome by an interrogated bases in reads fraction of 0 will lead to 'Inf', which then multiplied by the raw mutation count will equal NaN. However, in these edge cases, the correct mutation count is set to 0, because it should not be possible to detect a mutation in a trinucleotide context that was not interrogated in the reads.
        fwdmutations.trinucleotide_counts64.filteredgenomecorrected <- apply(merge3way(fwdmutations.trinucleotide_counts64,genome.noN.fwdfiltered.trinucleotide_freq.fwd64,mutationreads.fwd.trinucleotide_freq.fwd64),1,function(x){x[1]*x[2]/x[3]})
        fwdmutations.trinucleotide_counts64.filteredgenomecorrected[is.na(fwdmutations.trinucleotide_counts64.filteredgenomecorrected)] <- 0
      
        revmutations.trinucleotide_counts64.filteredgenomecorrected <- apply(merge3way(revmutations.trinucleotide_counts64,genome.noN.revfiltered.trinucleotide_freq.rev64,mutationreads.rev.trinucleotide_freq.rev64),1,function(x){x[1]*x[2]/x[3]})
        revmutations.trinucleotide_counts64.filteredgenomecorrected[is.na(revmutations.trinucleotide_counts64.filteredgenomecorrected)] <- 0
        
        fwdrevmutations.trinucleotide_counts64.filteredgenomecorrected <- apply(merge3way(fwdrevmutations.trinucleotide_counts64,genome.noN.fwdrevfiltered.trinucleotide_freq.fwdrev64,mutationreads.fwdrev.trinucleotide_freq.fwdrev64),1,function(x){x[1]*x[2]/x[3]})
        fwdrevmutations.trinucleotide_counts64.filteredgenomecorrected[is.na(fwdrevmutations.trinucleotide_counts64.filteredgenomecorrected)] <- 0
      
        dsDNAmutations.trinucleotide_counts32.filteredgenomecorrected <- apply(merge3way(dsDNAmutations.trinucleotide_counts32,genome.noN.zmwfiltered.trinucleotide_freq.32,mutationreads.zmw.trinucleotide_freq.32),1,function(x){x[1]*x[2]/x[3]})
        dsDNAmutations.trinucleotide_counts32.filteredgenomecorrected[is.na(dsDNAmutations.trinucleotide_counts32.filteredgenomecorrected)] <- 0
      
      #E. Output mutation frequencies corrected for trinucleotide context of interrogated bases in the genome
      num_interrogatedgenomeCorrected_FwdRevmutations.output <- sum(fwdrevmutations.trinucleotide_counts64.filteredgenomecorrected)
			freq_interrogatedgenomeCorrected_FwdRevmutations.output <- num_interrogatedgenomeCorrected_FwdRevmutations.output/numFwdRevStrandBases_interrogated.output
			freq_interrogatedgenomeCorrected_FwdRevmutations_lci.output <- freq_interrogatedgenomeCorrected_FwdRevmutations.output*raw_FwdRev_lci_to_calc_ratio
			freq_interrogatedgenomeCorrected_FwdRevmutations_uci.output <- freq_interrogatedgenomeCorrected_FwdRevmutations.output*raw_FwdRev_uci_to_calc_ratio
			
      num_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations.output <- num_interrogatedgenomeCorrected_FwdRevmutations.output/sqrt(sensitivity)
			freq_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations.output <- num_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations.output/numFwdRevStrandBases_interrogated.output
			freq_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations_lci.output <- freq_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations.output*raw_FwdRev_lci_to_calc_ratio
			freq_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations_uci.output <- freq_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations.output*raw_FwdRev_uci_to_calc_ratio

      num_interrogatedgenomeCorrected_ZMWmutations.output <- sum(dsDNAmutations.trinucleotide_counts32.filteredgenomecorrected)
			freq_interrogatedgenomeCorrected_ZMWmutations.output <- num_interrogatedgenomeCorrected_ZMWmutations.output/numZMWBases_interrogated.output
			freq_interrogatedgenomeCorrected_ZMWmutations_lci.output <- freq_interrogatedgenomeCorrected_ZMWmutations.output*raw_ZMW_lci_to_calc_ratio
			freq_interrogatedgenomeCorrected_ZMWmutations_uci.output <- freq_interrogatedgenomeCorrected_ZMWmutations.output*raw_ZMW_uci_to_calc_ratio
			
      num_sensCorrected_interrogatedgenomeCorrected_ZMWmutations.output <- num_interrogatedgenomeCorrected_ZMWmutations.output/sensitivity
			freq_sensCorrected_interrogatedgenomeCorrected_ZMWmutations.output <- num_sensCorrected_interrogatedgenomeCorrected_ZMWmutations.output/numZMWBases_interrogated.output
			freq_sensCorrected_interrogatedgenomeCorrected_ZMWmutations_lci.output <- freq_sensCorrected_interrogatedgenomeCorrected_ZMWmutations.output*raw_ZMW_lci_to_calc_ratio
			freq_sensCorrected_interrogatedgenomeCorrected_ZMWmutations_uci.output <- freq_sensCorrected_interrogatedgenomeCorrected_ZMWmutations.output*raw_ZMW_uci_to_calc_ratio

      ## Calculate: 1) Mutation frequencies corrected for trinucleotide context of full genome; 2) Genome-wide mutation burdens per haploid and diploid cell, which assumes that the mutation rate in the regions of the genome that are not analyzed (i.e. filtered out) are the same as regions of the genome that were analyzed; and, 3) Haploid equivalents that were interrogated for each sample.

      #A. Correct fwd, rev, zmw trinucleotide mutation counts for ratio of trinucleotide context distributions of full genome to interrogated bases in the reads. For any trinucleotide context for which the fraction of the trinucleotide in the interrogated bases in the reads is 0, the final corrected mutation count will result in NaN, because when dividing the interrogated bases in the genome by an interrogated bases in reads fraction of 0 will lead to 'Inf', which then multiplied by the raw mutation count will equal NaN. However, in these cases, the correct mutation count is set to 0, because it should not be possible to detect a mutation in a trinucleotide context that was not interrogated in the reads. Nevertheless, this will skew the results, because technically in a full genome, these sites are present in the genome. Therefore, we will output in the results which trinucleotides contexts in the interrogated base distribution had a count of 0.
      fwdmutations.trinucleotide_counts64.genomecorrected <- apply(merge3way(fwdmutations.trinucleotide_counts64,genome.noN.trinucleotide_freq.fwd64,mutationreads.fwd.trinucleotide_freq.fwd64),1,function(x){x[1]*x[2]/x[3]})
      fwdmutations.trinucleotide_counts64.genomecorrected[is.na(fwdmutations.trinucleotide_counts64.genomecorrected)] <- 0
      mutationreads.fwd.trinucleotidesmissing <- names(mutationreads.fwd.trinucleotide_freq.fwd64)[mutationreads.fwd.trinucleotide_freq.fwd64==0]
      
      revmutations.trinucleotide_counts64.genomecorrected <- apply(merge3way(revmutations.trinucleotide_counts64,genome.noN.trinucleotide_freq.rev64,mutationreads.rev.trinucleotide_freq.rev64),1,function(x){x[1]*x[2]/x[3]})
      revmutations.trinucleotide_counts64.genomecorrected[is.na(revmutations.trinucleotide_counts64.genomecorrected)] <- 0
      mutationreads.rev.trinucleotidesmissing <- names(mutationreads.rev.trinucleotide_freq.rev64)[mutationreads.rev.trinucleotide_freq.rev64==0]

			fwdrevmutations.trinucleotide_counts64.genomecorrected <- apply(merge3way(fwdrevmutations.trinucleotide_counts64,genome.noN.trinucleotide_freq.fwdrev64,mutationreads.fwdrev.trinucleotide_freq.fwdrev64),1,function(x){x[1]*x[2]/x[3]})
      fwdrevmutations.trinucleotide_counts64.genomecorrected[is.na(fwdrevmutations.trinucleotide_counts64.genomecorrected)] <- 0
      mutationreads.fwdrev.trinucleotidesmissing <- names(mutationreads.fwdrev.trinucleotide_freq.fwdrev64)[mutationreads.fwdrev.trinucleotide_freq.fwdrev64==0]
      
      dsDNAmutations.trinucleotide_counts32.genomecorrected <- apply(merge3way(dsDNAmutations.trinucleotide_counts32,genome.noN.trinucleotide_freq.32,mutationreads.zmw.trinucleotide_freq.32),1,function(x){x[1]*x[2]/x[3]})
      dsDNAmutations.trinucleotide_counts32.genomecorrected[is.na(dsDNAmutations.trinucleotide_counts32.genomecorrected)] <- 0
      mutationreads.zmw.trinucleotidesmissing <- names(mutationreads.zmw.trinucleotide_freq.32)[mutationreads.zmw.trinucleotide_freq.32==0]
      
      #B. Output mutation frequencies corrected for trinucleotide context of full genome
			num_fullgenomeCorrected_FwdRevmutations.output <- sum(fwdrevmutations.trinucleotide_counts64.genomecorrected)
			freq_fullgenomeCorrected_FwdRevmutations.output <- num_fullgenomeCorrected_FwdRevmutations.output/numFwdRevStrandBases_interrogated.output
			freq_fullgenomeCorrected_FwdRevmutations_lci.output <- freq_fullgenomeCorrected_FwdRevmutations.output*raw_FwdRev_lci_to_calc_ratio
			freq_fullgenomeCorrected_FwdRevmutations_uci.output <- freq_fullgenomeCorrected_FwdRevmutations.output*raw_FwdRev_uci_to_calc_ratio
			
			num_sensCorrected_fullgenomeCorrected_FwdRevmutations.output <- num_fullgenomeCorrected_FwdRevmutations.output/sqrt(sensitivity)
			freq_sensCorrected_fullgenomeCorrected_FwdRevmutations.output <- num_sensCorrected_fullgenomeCorrected_FwdRevmutations.output/numFwdRevStrandBases_interrogated.output
			freq_sensCorrected_fullgenomeCorrected_FwdRevmutations_lci.output <- freq_sensCorrected_fullgenomeCorrected_FwdRevmutations.output*raw_FwdRev_lci_to_calc_ratio
			freq_sensCorrected_fullgenomeCorrected_FwdRevmutations_uci.output <- freq_sensCorrected_fullgenomeCorrected_FwdRevmutations.output*raw_FwdRev_uci_to_calc_ratio
			
			num_fullgenomeCorrected_ZMWmutations.output <- sum(dsDNAmutations.trinucleotide_counts32.genomecorrected)
			freq_fullgenomeCorrected_ZMWmutations.output <- num_fullgenomeCorrected_ZMWmutations.output/numZMWBases_interrogated.output
			freq_fullgenomeCorrected_ZMWmutations_lci.output <- freq_fullgenomeCorrected_ZMWmutations.output*raw_ZMW_lci_to_calc_ratio
			freq_fullgenomeCorrected_ZMWmutations_uci.output <- freq_fullgenomeCorrected_ZMWmutations.output*raw_ZMW_uci_to_calc_ratio
				
			num_sensCorrected_fullgenomeCorrected_ZMWmutations.output <- num_fullgenomeCorrected_ZMWmutations.output/sensitivity
			freq_sensCorrected_fullgenomeCorrected_ZMWmutations.output <- num_sensCorrected_fullgenomeCorrected_ZMWmutations.output/numZMWBases_interrogated.output
			freq_sensCorrected_fullgenomeCorrected_ZMWmutations_lci.output <- freq_sensCorrected_fullgenomeCorrected_ZMWmutations.output*raw_ZMW_lci_to_calc_ratio
			freq_sensCorrected_fullgenomeCorrected_ZMWmutations_uci.output <- freq_sensCorrected_fullgenomeCorrected_ZMWmutations.output*raw_ZMW_uci_to_calc_ratio
			
			FwdRev_interrogatedbases_trinucleotides_missing.output <- paste(mutationreads.fwdrev.trinucleotidesmissing,collapse=",")
			ZMW_interrogatedbases_trinucleotides_missing.output <- paste(mutationreads.zmw.trinucleotidesmissing,collapse=",")
			
      #C. Extrapolate fwd, rev, zmw mutation frequencies to haploid and diploid genome-corrected mutation burden counts. Note, this only extrapolates to chromosomes being analyzed per the configuration file (i.e. usually chr1-22 + chrX). For any trinucleotide context for which the number of counts of a trinucleotide context in the interrogated bases in the reads is 0, the final corrected mutation count will result in NaN, because dividing the interrogated bases in the genome by an interrogated bases in reads fraction of 0 will lead to 'Inf', which then multiplied by the raw mutation count will equal NaN. However, in these cases, the correct mutation count for that trinucleotide context is set to 0, because it should not be possible to detect a mutation in a trinucleotide context that was not interrogated in the reads. Nevertheless, this will skew the results, because technically in a full genome, these sites are present in the genome. Therefore, we will output in the results which trinucleotides contexts in the interrogated bases distribution had a count of 0.
      fwdmutations.genomecorrected.haploidcounts64 <- apply(merge3way(fwdmutations.trinucleotide_counts64,genome.noN.trinucleotide_counts.fwd64,mutationreads.fwd.trinucleotide_counts.fwd64),1,function(x){x[1]*x[2]/x[3]})
      fwdmutations.genomecorrected.haploidcounts64[is.na(fwdmutations.genomecorrected.haploidcounts64)] <- 0

      revmutations.genomecorrected.haploidcounts64 <- apply(merge3way(revmutations.trinucleotide_counts64,genome.noN.trinucleotide_counts.rev64,mutationreads.rev.trinucleotide_counts.rev64),1,function(x){x[1]*x[2]/x[3]})
      revmutations.genomecorrected.haploidcounts64[is.na(revmutations.genomecorrected.haploidcounts64)] <- 0

      fwdrevmutations.genomecorrected.haploidcounts64 <- apply(merge3way(fwdrevmutations.trinucleotide_counts64,genome.noN.trinucleotide_counts.fwdrev64,mutationreads.fwdrev.trinucleotide_counts.fwdrev64),1,function(x){x[1]*x[2]/x[3]})
      fwdrevmutations.genomecorrected.haploidcounts64[is.na(fwdrevmutations.genomecorrected.haploidcounts64)] <- 0
      
      dsDNAmutations.genomecorrected.haploidcounts32 <- apply(merge3way(dsDNAmutations.trinucleotide_counts32,genome.noN.trinucleotide_counts.32,mutationreads.zmw.trinucleotide_counts.32),1,function(x){x[1]*x[2]/x[3]})
      dsDNAmutations.genomecorrected.haploidcounts32[is.na(dsDNAmutations.genomecorrected.haploidcounts32)] <- 0
      
      #D. Output haploid and diploid genome-wide mutation burdens per cell, and haploid equivalents that were sequenced, as well as (95%) confidence intervals, corrected for full genome
			num_fullgenomeCorrected_FwdRevmutations.perHaploid.output <- sum(fwdrevmutations.genomecorrected.haploidcounts64)
			num_fullgenomeCorrected_FwdRevmutations.perHaploid_lci.output <- num_fullgenomeCorrected_FwdRevmutations.perHaploid.output*raw_FwdRev_lci_to_calc_ratio
			num_fullgenomeCorrected_FwdRevmutations.perHaploid_uci.output <- num_fullgenomeCorrected_FwdRevmutations.perHaploid.output*raw_FwdRev_uci_to_calc_ratio
			
			num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid.output <- num_fullgenomeCorrected_FwdRevmutations.perHaploid.output/sqrt(sensitivity)
			num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid_lci.output <- num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid.output*raw_FwdRev_lci_to_calc_ratio
			num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid_uci.output <- num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid.output*raw_FwdRev_uci_to_calc_ratio
			
			num_fullgenomeCorrected_FwdRevmutations.perDiploid.output <- num_fullgenomeCorrected_FwdRevmutations.perHaploid.output*2
			num_fullgenomeCorrected_FwdRevmutations.perDiploid_lci.output <- num_fullgenomeCorrected_FwdRevmutations.perHaploid_lci.output*2
			num_fullgenomeCorrected_FwdRevmutations.perDiploid_uci.output <- num_fullgenomeCorrected_FwdRevmutations.perHaploid_uci.output*2
			
			num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perDiploid.output <- num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid.output*2
			num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perDiploid_lci.output <- num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid_lci.output*2
			num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perDiploid_uci.output <- num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid_uci.output*2
			
			num_fullgenomeCorrected_ZMWmutations.perHaploid.output <- sum(dsDNAmutations.genomecorrected.haploidcounts32)
			num_fullgenomeCorrected_ZMWmutations.perHaploid_lci.output <- num_fullgenomeCorrected_ZMWmutations.perHaploid.output*raw_ZMW_lci_to_calc_ratio
			num_fullgenomeCorrected_ZMWmutations.perHaploid_uci.output <- num_fullgenomeCorrected_ZMWmutations.perHaploid.output*raw_ZMW_uci_to_calc_ratio
			
			num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid.output <- num_fullgenomeCorrected_ZMWmutations.perHaploid.output/sensitivity
			num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid_lci.output <- num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid.output*raw_ZMW_lci_to_calc_ratio
			num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid_uci.output <- num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid.output*raw_ZMW_uci_to_calc_ratio
				
			num_fullgenomeCorrected_ZMWmutations.perDiploid.output <- num_fullgenomeCorrected_ZMWmutations.perHaploid.output*2
			num_fullgenomeCorrected_ZMWmutations.perDiploid_lci.output <- num_fullgenomeCorrected_ZMWmutations.perHaploid_lci.output*2
			num_fullgenomeCorrected_ZMWmutations.perDiploid_uci.output <- num_fullgenomeCorrected_ZMWmutations.perHaploid_uci.output*2
			
			num_sensCorrected_fullgenomeCorrected_ZMWmutations.perDiploid.output <- num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid.output*2
			num_sensCorrected_fullgenomeCorrected_ZMWmutations.perDiploid_lci.output <- num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid_lci.output*2
			num_sensCorrected_fullgenomeCorrected_ZMWmutations.perDiploid_uci.output <- num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid_uci.output*2
			
			num_dsDNAhaploidequivalents_interrogated_in_FwdRevreads.output <- (numFwdRevStrandBases_interrogated.output/2)/sum(width(genome.noN.granges))
			num_dsDNAhaploidequivalents_interrogated_in_ZMWreads.output <- numZMWBases_interrogated.output/sum(width(genome.noN.granges))

      ##Save output calculations in data frame
      temp.df <- data.frame(sampleid=sampleid,genome=genome,filterid=filterindex,
				numsensitivity_variants=numsensitivity_variants.output,sensitivity=sensitivity.output,
				numZMW_filteredbyZMWfilters_beforemaxmutperZMW=numZMW_filteredbyZMWfilters_beforemaxmutperZMW.output,
				numZMW=numZMW.output,
				fracZMW_filteredbyZMWfilters_beforemaxmutperZMW=fracZMW_filteredbyZMWfilters_beforemaxmutperZMW.output,
				numFwdRevStrand_filtered_bymaxmutperStrand=numFwdRevStrand_filtered_bymaxmutperStrand.output,
				numFwdRevStrand=numFwdRevStrand.output,
				fracFwdRevStrand_filtered_bymaxmutperStrand=fracFwdRevStrand_filtered_bymaxmutperStrand.output,
				numZMW_filtered_bymaxmutperZMW=numZMW_filtered_bymaxmutperZMW.output,
				fracZMW_filtered_bymaxmutperZMW=fracZMW_filtered_bymaxmutperZMW.output,
				numFwdRevStrandBases_filteredbyZMWfilters_beforemaxmutperStrand=numFwdRevStrandBases_filteredbyZMWfilters_beforemaxmutperStrand.output,
				numFwdRevStrandBases=numFwdRevStrandBases.output,
				fracFwdRevStrandBases_filteredbyZMWfilters_beforemaxmutperStrand=fracFwdRevStrandBases_filteredbyZMWfilters_beforemaxmutperStrand.output,
				numZMWBases_filteredbyZMWfilters_beforemaxmutperZMW=numZMWBases_filteredbyZMWfilters_beforemaxmutperZMW.output,
				numZMWBases=numZMWBases.output,
				fracZMWBases_filteredbyZMWfilters_beforemaxmutperZMW=fracZMWBases_filteredbyZMWfilters_beforemaxmutperZMW.output,
				numFwdRevStrandBases_filtered_bymaxmutperStrand=numFwdRevStrandBases_filtered_bymaxmutperStrand.output,
				fracFwdRevStrandBases_filtered_bymaxmutperStrand=fracFwdRevStrandBases_filtered_bymaxmutperStrand.output,
				numZMWBases_filtered_bymaxmutperZMW=numZMWBases_filtered_bymaxmutperZMW.output,
				fracZMWBases_filtered_bymaxmutperZMW=fracZMWBases_filtered_bymaxmutperZMW.output,
				numFwdRevStrandBases_filtered_afterZMWfiltersandmaxmutperStrand=numFwdRevStrandBases_filtered_afterZMWfiltersandmaxmutperStrand.output,
				numFwdRevStrandBases_afterZMWfilterandmaxmutperStrand=numFwdRevStrandBases_afterZMWfilterandmaxmutperStrand.output,
				fracFwdRevStrandBases_filtered_afterZMWfiltersandmaxmutperStrand=fracFwdRevStrandBases_filtered_afterZMWfiltersandmaxmutperStrand.output,
				numZMWBases_filtered_afterZMWfiltersandmaxmutperZMW=numZMWBases_filtered_afterZMWfiltersandmaxmutperZMW.output,
				numZMWBases_afterZMWfilterandmaxmutperZMW=numZMWBases_afterZMWfilterandmaxmutperZMW.output,
				fracZMWBases_filtered_afterZMWfiltersandmaxmutperZMW=fracZMWBases_filtered_afterZMWfiltersandmaxmutperZMW.output,
				numFwdRevStrandBases_interrogated=numFwdRevStrandBases_interrogated.output,
				fracFwdRevStrandBases_interrogated_from_numFwdRevStrandBases=fracFwdRevStrandBases_interrogated_from_numFwdRevStrandBases.output,
				numZMWBases_interrogated=numZMWBases_interrogated.output,
				fracZMWBases_interrogated_from_numZMWBases=fracZMWBases_interrogated_from_numZMWBases.output,
				numSubreadBases=numSubreadBases.output,
				fracZMWBases_interrogated_from_numSubreadBases=fracZMWBases_interrogated_from_numSubreadBases.output,
				num_raw_FwdRevmutations=num_raw_FwdRevmutations.output,
				num_raw_FwdRevmutations_lci=num_raw_FwdRevmutations_lci.output,
				num_raw_FwdRevmutations_uci=num_raw_FwdRevmutations_uci.output,
				freq_raw_FwdRevmutations=freq_raw_FwdRevmutations.output,
				freq_raw_FwdRevmutations_lci=freq_raw_FwdRevmutations_lci.output,
				freq_raw_FwdRevmutations_uci=freq_raw_FwdRevmutations_uci.output,
				num_sensCorrected_FwdRevmutations=num_sensCorrected_FwdRevmutations.output,
				freq_sensCorrected_FwdRevmutations=freq_sensCorrected_FwdRevmutations.output,
				freq_sensCorrected_FwdRevmutations_lci=freq_sensCorrected_FwdRevmutations_lci.output,
				freq_sensCorrected_FwdRevmutations_uci=freq_sensCorrected_FwdRevmutations_uci.output,
				num_raw_ZMWmutations=num_raw_ZMWmutations.output,
				num_raw_ZMWmutations_lci=num_raw_ZMWmutations_lci.output,
				num_raw_ZMWmutations_uci=num_raw_ZMWmutations_uci.output,
				freq_raw_ZMWmutations=freq_raw_ZMWmutations.output,
				freq_raw_ZMWmutations_lci=freq_raw_ZMWmutations_lci.output,
				freq_raw_ZMWmutations_uci=freq_raw_ZMWmutations_uci.output,
				num_sensCorrected_ZMWmutations=num_sensCorrected_ZMWmutations.output,
				freq_sensCorrected_ZMWmutations=freq_sensCorrected_ZMWmutations.output,
				freq_sensCorrected_ZMWmutations_lci=freq_sensCorrected_ZMWmutations_lci.output,
				freq_sensCorrected_ZMWmutations_uci=freq_sensCorrected_ZMWmutations_uci.output,
				FwdRev_interrogatedbases_trinucleotides_missing=FwdRev_interrogatedbases_trinucleotides_missing.output,
				ZMW_interrogatedbases_trinucleotides_missing=ZMW_interrogatedbases_trinucleotides_missing.output,
				num_interrogatedgenomeCorrected_FwdRevmutations=num_interrogatedgenomeCorrected_FwdRevmutations.output,
				freq_interrogatedgenomeCorrected_FwdRevmutations=freq_interrogatedgenomeCorrected_FwdRevmutations.output,
				freq_interrogatedgenomeCorrected_FwdRevmutations_lci=freq_interrogatedgenomeCorrected_FwdRevmutations_lci.output,
				freq_interrogatedgenomeCorrected_FwdRevmutations_uci=freq_interrogatedgenomeCorrected_FwdRevmutations_uci.output,
				num_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations=num_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations.output,
				freq_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations=freq_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations.output,
				freq_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations_lci=freq_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations_lci.output,
				freq_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations_uci=freq_sensCorrected_interrogatedgenomeCorrected_FwdRevmutations_uci.output,
				num_interrogatedgenomeCorrected_ZMWmutations=num_interrogatedgenomeCorrected_ZMWmutations.output,
				freq_interrogatedgenomeCorrected_ZMWmutations=freq_interrogatedgenomeCorrected_ZMWmutations.output,
				freq_interrogatedgenomeCorrected_ZMWmutations_lci=freq_interrogatedgenomeCorrected_ZMWmutations_lci.output,
				freq_interrogatedgenomeCorrected_ZMWmutations_uci=freq_interrogatedgenomeCorrected_ZMWmutations_uci.output,
				num_sensCorrected_interrogatedgenomeCorrected_ZMWmutations=num_sensCorrected_interrogatedgenomeCorrected_ZMWmutations.output,
				freq_sensCorrected_interrogatedgenomeCorrected_ZMWmutations=freq_sensCorrected_interrogatedgenomeCorrected_ZMWmutations.output,
				freq_sensCorrected_interrogatedgenomeCorrected_ZMWmutations_lci=freq_sensCorrected_interrogatedgenomeCorrected_ZMWmutations_lci.output,
				freq_sensCorrected_interrogatedgenomeCorrected_ZMWmutations_uci=freq_sensCorrected_interrogatedgenomeCorrected_ZMWmutations_uci.output,
				num_fullgenomeCorrected_FwdRevmutations=num_fullgenomeCorrected_FwdRevmutations.output,
				freq_fullgenomeCorrected_FwdRevmutations=freq_fullgenomeCorrected_FwdRevmutations.output,
				freq_fullgenomeCorrected_FwdRevmutations_lci=freq_fullgenomeCorrected_FwdRevmutations_lci.output,
				freq_fullgenomeCorrected_FwdRevmutations_uci=freq_fullgenomeCorrected_FwdRevmutations_uci.output,
				num_sensCorrected_fullgenomeCorrected_FwdRevmutations=num_sensCorrected_fullgenomeCorrected_FwdRevmutations.output,
				freq_sensCorrected_fullgenomeCorrected_FwdRevmutations=freq_sensCorrected_fullgenomeCorrected_FwdRevmutations.output,
				freq_sensCorrected_fullgenomeCorrected_FwdRevmutations_lci=freq_sensCorrected_fullgenomeCorrected_FwdRevmutations_lci.output,
				freq_sensCorrected_fullgenomeCorrected_FwdRevmutations_uci=freq_sensCorrected_fullgenomeCorrected_FwdRevmutations_uci.output,
				num_fullgenomeCorrected_ZMWmutations=num_fullgenomeCorrected_ZMWmutations.output,
				freq_fullgenomeCorrected_ZMWmutations=freq_fullgenomeCorrected_ZMWmutations.output,
				freq_fullgenomeCorrected_ZMWmutations_lci=freq_fullgenomeCorrected_ZMWmutations_lci.output,
				freq_fullgenomeCorrected_ZMWmutations_uci=freq_fullgenomeCorrected_ZMWmutations_uci.output,
				num_sensCorrected_fullgenomeCorrected_ZMWmutations=num_sensCorrected_fullgenomeCorrected_ZMWmutations.output,
				freq_sensCorrected_fullgenomeCorrected_ZMWmutations=freq_sensCorrected_fullgenomeCorrected_ZMWmutations.output,
				freq_sensCorrected_fullgenomeCorrected_ZMWmutations_lci=freq_sensCorrected_fullgenomeCorrected_ZMWmutations_lci.output,
				freq_sensCorrected_fullgenomeCorrected_ZMWmutations_uci=freq_sensCorrected_fullgenomeCorrected_ZMWmutations_uci.output,
				num_fullgenomeCorrected_FwdRevmutations.perHaploid=num_fullgenomeCorrected_FwdRevmutations.perHaploid.output,
				num_fullgenomeCorrected_FwdRevmutations.perHaploid_lci=num_fullgenomeCorrected_FwdRevmutations.perHaploid_lci.output,
				num_fullgenomeCorrected_FwdRevmutations.perHaploid_uci=num_fullgenomeCorrected_FwdRevmutations.perHaploid_uci.output,
				num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid=num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid.output,
				num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid_lci=num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid_lci.output,
				num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid_uci=num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perHaploid_uci.output,
				num_fullgenomeCorrected_FwdRevmutations.perDiploid=num_fullgenomeCorrected_FwdRevmutations.perDiploid.output,
				num_fullgenomeCorrected_FwdRevmutations.perDiploid_lci=num_fullgenomeCorrected_FwdRevmutations.perDiploid_lci.output,
				num_fullgenomeCorrected_FwdRevmutations.perDiploid_uci=num_fullgenomeCorrected_FwdRevmutations.perDiploid_uci.output,
				num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perDiploid=num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perDiploid.output,
				num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perDiploid_lci=num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perDiploid_lci.output,
				num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perDiploid_uci=num_sensCorrected_fullgenomeCorrected_FwdRevmutations.perDiploid_uci.output,
				num_fullgenomeCorrected_ZMWmutations.perHaploid=num_fullgenomeCorrected_ZMWmutations.perHaploid.output,
				num_fullgenomeCorrected_ZMWmutations.perHaploid_lci=num_fullgenomeCorrected_ZMWmutations.perHaploid_lci.output,
				num_fullgenomeCorrected_ZMWmutations.perHaploid_uci=num_fullgenomeCorrected_ZMWmutations.perHaploid_uci.output,
				num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid=num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid.output,
				num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid_lci=num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid_lci.output,
				num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid_uci=num_sensCorrected_fullgenomeCorrected_ZMWmutations.perHaploid_uci.output,
				num_fullgenomeCorrected_ZMWmutations.perDiploid=num_fullgenomeCorrected_ZMWmutations.perDiploid.output,
				num_fullgenomeCorrected_ZMWmutations.perDiploid_lci=num_fullgenomeCorrected_ZMWmutations.perDiploid_lci.output,
				num_fullgenomeCorrected_ZMWmutations.perDiploid_uci=num_fullgenomeCorrected_ZMWmutations.perDiploid_uci.output,
				num_sensCorrected_fullgenomeCorrected_ZMWmutations.perDiploid=num_sensCorrected_fullgenomeCorrected_ZMWmutations.perDiploid.output,
				num_sensCorrected_fullgenomeCorrected_ZMWmutations.perDiploid_lci=num_sensCorrected_fullgenomeCorrected_ZMWmutations.perDiploid_lci.output,
				num_sensCorrected_fullgenomeCorrected_ZMWmutations.perDiploid_uci=num_sensCorrected_fullgenomeCorrected_ZMWmutations.perDiploid_uci.output,
				num_dsDNAhaploidequivalents_interrogated_in_FwdRevreads=num_dsDNAhaploidequivalents_interrogated_in_FwdRevreads.output,
				num_dsDNAhaploidequivalents_interrogated_in_ZMWreads=num_dsDNAhaploidequivalents_interrogated_in_ZMWreads.output,
				row.names=NULL)
      
      output.df <- rbind(output.df,temp.df)
      
      ##Save fwd, rev, combined fwd/rev, dsDNA mutations, trinucleotide context distributions of mutations, interrogated read bases, interrogated genome, and full genome. This will be used in subsequent mutational signature analyses.
      mutations[[samplename]]$fwd <- fwdmutations.granges
      mutations[[samplename]]$rev <- revmutations.granges
      mutations[[samplename]]$fwdrev <- c(fwdmutations.granges,revmutations.granges)
      mutations[[samplename]]$dsDNA <- dsDNAmutations.granges
      
      trinucleotide_contexts[[samplename]]$fwd$genomefreq <- genome.noN.trinucleotide_freq.fwd64
      trinucleotide_contexts[[samplename]]$rev$genomefreq <- genome.noN.trinucleotide_freq.rev64
      trinucleotide_contexts[[samplename]]$fwdrev$genomefreq <- genome.noN.trinucleotide_freq.fwdrev64
      trinucleotide_contexts[[samplename]]$dsDNA$genomefreq <- genome.noN.trinucleotide_freq.32
      
      trinucleotide_contexts[[samplename]]$fwd$filteredgenomefreq <- genome.noN.fwdfiltered.trinucleotide_freq.fwd64
      trinucleotide_contexts[[samplename]]$rev$filteredgenomefreq <- genome.noN.revfiltered.trinucleotide_freq.rev64
      trinucleotide_contexts[[samplename]]$fwdrev$filteredgenomefreq <- genome.noN.fwdrevfiltered.trinucleotide_freq.fwdrev64
      trinucleotide_contexts[[samplename]]$dsDNA$filteredgenomefreq <- genome.noN.zmwfiltered.trinucleotide_freq.32
      
      trinucleotide_contexts[[samplename]]$fwd$mutationreadsfreq <- mutationreads.fwd.trinucleotide_freq.fwd64
      trinucleotide_contexts[[samplename]]$rev$mutationreadsfreq <- mutationreads.rev.trinucleotide_freq.rev64
      trinucleotide_contexts[[samplename]]$fwdrev$mutationreadsfreq <- mutationreads.fwdrev.trinucleotide_freq.fwdrev64
      trinucleotide_contexts[[samplename]]$dsDNA$mutationreadsfreq <- mutationreads.zmw.trinucleotide_freq.32

      trinucleotide_contexts[[samplename]]$fwd$mutationcounts <- fwdmutations.trinucleotide_counts64
      trinucleotide_contexts[[samplename]]$rev$mutationcounts <- revmutations.trinucleotide_counts64
      trinucleotide_contexts[[samplename]]$fwdrev$mutationcounts <- fwdrevmutations.trinucleotide_counts64
      trinucleotide_contexts[[samplename]]$dsDNA$mutationcounts <- dsDNAmutations.trinucleotide_counts32
                  
      trinucleotide_contexts[[samplename]]$fwd$mutationcounts.genomecorrected <- fwdmutations.trinucleotide_counts64.genomecorrected
      trinucleotide_contexts[[samplename]]$rev$mutationcounts.genomecorrected <- revmutations.trinucleotide_counts64.genomecorrected
      trinucleotide_contexts[[samplename]]$fwdrev$mutationcounts.genomecorrected <- fwdrevmutations.trinucleotide_counts64.genomecorrected
      trinucleotide_contexts[[samplename]]$dsDNA$mutationcounts.genomecorrected <- dsDNAmutations.trinucleotide_counts32.genomecorrected

      trinucleotide_contexts[[samplename]]$fwd$mutationcounts.filteredgenomecorrected <- fwdmutations.trinucleotide_counts64.filteredgenomecorrected
      trinucleotide_contexts[[samplename]]$rev$mutationcounts.filteredgenomecorrected <- revmutations.trinucleotide_counts64.filteredgenomecorrected
      trinucleotide_contexts[[samplename]]$fwdrev$mutationcounts.filteredgenomecorrected <- fwdrevmutations.trinucleotide_counts64.filteredgenomecorrected
      trinucleotide_contexts[[samplename]]$dsDNA$mutationcounts.filteredgenomecorrected <- dsDNAmutations.trinucleotide_counts32.filteredgenomecorrected
      
      cat("DONE\n\n")
    
    }
  }
}

#Print table-formatted output
write.table(output.df,file=paste0(yaml.config$analysisoutput_path,"/mutation_frequencies.table.output.",timestamp,".txt"),quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)

#If no variants, quit and save mutation_frequencies.RDS with note
if(is.null(output.df)){
  cat("\n#### No variants in any sample/filter analysis... ending script!\n")
	qsave("No reads or mutations in selected chromosomes",paste0(yaml.config$analysisoutput_path,"/mutation_frequencies.RDS"))
	quit(save="no")
}

cat("#### Creating SNV catalogue for signature.tools.lib package from prior print_mutations output...")
#Note creating it from prior print_mutations output because that is simpler than needlessly recreating all the code for extracting trinucleotide contexts, synthesized vs template strand transformations, etc. ssDNA mutations are extracted from the PacBio template strand direction.
print_mutations_output <- qread(paste0(yaml.config$analysisoutput_path,"/print_mutations.output.RDS"))
samplenames <- paste0(print_mutations_output$sampleid,".",print_mutations_output$genome,".",print_mutations_output$filterid)

#Make list of all trinucleotide substitions, and also collapsed to central pyrimidine
trinucleotides_64_subs <- grep("A>A|C>C|G>G|T>T",apply(expand.grid(c("A","C","G","T"),"[",c("A","C","G","T"),">",c("A","C","G","T"),"]",c("A","C","G","T")),1,paste,collapse=""),value=TRUE,invert=TRUE)
trinucleotides_64_subs <- trinucleotides_64_subs[order(substr(trinucleotides_64_subs,3,5),substr(trinucleotides_64_subs,1,1),substr(trinucleotides_64_subs,7,7))]

trinucleotides_32_pyr_subs <- grep("A>A|C>C|G>G|T>T",apply(expand.grid(c("A","C","G","T"),"[",c("C","T"),">",c("A","C","G","T"),"]",c("A","C","G","T")),1,paste,collapse=""),value=TRUE,invert=TRUE)
trinucleotides_32_pyr_subs <- trinucleotides_32_pyr_subs[order(substr(trinucleotides_32_pyr_subs,3,5),substr(trinucleotides_32_pyr_subs,1,1),substr(trinucleotides_32_pyr_subs,7,7))]

ssDNAcatalogue <- data.frame(temp=rep(0,length(trinucleotides_64_subs)),row.names=trinucleotides_64_subs)
dsDNAcatalogue <- data.frame(temp=rep(0,length(trinucleotides_64_subs)),row.names=trinucleotides_64_subs)
ssDNAcataloguecollapsed <- data.frame(temp=rep(0,length(trinucleotides_32_pyr_subs)),row.names=trinucleotides_32_pyr_subs)
dsDNAcataloguecollapsed <- data.frame(temp=rep(0,length(trinucleotides_32_pyr_subs)),row.names=trinucleotides_32_pyr_subs)

for(i in unique(samplenames)){
  ssDNAcatalogue[,i] <- 0
  ssDNAcataloguecollapsed[,i] <- 0
  dsDNAcatalogue[,i] <- 0
  dsDNAcataloguecollapsed[,i] <- 0

	selectmutations <- samplenames==i & print_mutations_output$strandtype=="ssDNA"
  if(any(selectmutations)){
		
  	counts <- table(paste0(substr(print_mutations_output[selectmutations,"tnc_ssDNAtemplatestrand"],1,1),"[",print_mutations_output[selectmutations,"ref_ssDNAtemplatestrand"],">",print_mutations_output[selectmutations,"alt_ssDNAtemplatestrand"],"]",substr(print_mutations_output[selectmutations,"tnc_ssDNAtemplatestrand"],3,3)))
  	ssDNAcatalogue[names(counts),i] <- counts
	
	  mutations_to_collapse <- print_mutations_output$ref_ssDNAtemplatestrand %in% c("A","G")
	
	  counts <- table(c(
						paste0(substr(print_mutations_output[selectmutations & !mutations_to_collapse,"tnc_ssDNAtemplatestrand"],1,1),
						"[",
						print_mutations_output[selectmutations & !mutations_to_collapse,"ref_ssDNAtemplatestrand"],">",
						print_mutations_output[selectmutations & !mutations_to_collapse,"alt_ssDNAtemplatestrand"],"]",
						substr(print_mutations_output[selectmutations & !mutations_to_collapse,"tnc_ssDNAtemplatestrand"],3,3)),  
						
						paste0(as.character(reverseComplement(DNAStringSet(substr(print_mutations_output[selectmutations & mutations_to_collapse,"tnc_ssDNAtemplatestrand"],3,3)))),
						"[",
						as.character(reverseComplement(DNAStringSet(print_mutations_output[selectmutations & mutations_to_collapse,"ref_ssDNAtemplatestrand"]))),
						">",
						as.character(reverseComplement(DNAStringSet(print_mutations_output[selectmutations & mutations_to_collapse,"alt_ssDNAtemplatestrand"]))),
						"]",
						as.character(reverseComplement(DNAStringSet(substr(print_mutations_output[selectmutations & mutations_to_collapse,"tnc_ssDNAtemplatestrand"],1,1)))))
						))
	  	counts <- counts[names(counts)!="[>]"] #This removes empty mutation table counts when there is no mutations to reverse complement
	  ssDNAcataloguecollapsed[names(counts),i] <- counts
  }

	
	selectmutations <- samplenames==i & print_mutations_output$strandtype=="dsDNA"
	if(any(selectmutations)){
  	counts <- table(paste0(substr(print_mutations_output[selectmutations,"tnc_refplusstrand"],1,1),"[",print_mutations_output[selectmutations,"ref_refplusstrand"],">",print_mutations_output[selectmutations,"alt_refplusstrand"],"]",substr(print_mutations_output[selectmutations,"tnc_refplusstrand"],3,3)))
  	dsDNAcatalogue[names(counts),i] <- counts
  
  	counts <- table(paste0(substr(print_mutations_output[selectmutations,"tnc32_refplusstrand"],1,1),"[",print_mutations_output[selectmutations,"ref32_refplusstrand"],">",print_mutations_output[selectmutations,"alt32_refplusstrand"],"]",substr(print_mutations_output[selectmutations,"tnc32_refplusstrand"],3,3)))
  	dsDNAcataloguecollapsed[names(counts),i] <- counts
	}
}

ssDNAcatalogue$temp <- NULL
dsDNAcatalogue$temp <- NULL
ssDNAcataloguecollapsed$temp <- NULL
dsDNAcataloguecollapsed$temp <- NULL

cat("DONE\n")

##Converting signature.tools.lib catalogue to SigFit data frame
cat("#### Converting signature.tools.lib catalogue to SigFit data frame...")
ssDNAdf <- matrix(nrow=0,ncol=4)
dsDNAdf <- matrix(nrow=0,ncol=4)
ssDNAdfcollapsed <- matrix(nrow=0,ncol=4)
dsDNAdfcollapsed <- matrix(nrow=0,ncol=4)

ref_ssDNAcatalogue <- substr(rownames(ssDNAcatalogue),3,3)
alt_ssDNAcatalogue <- substr(rownames(ssDNAcatalogue),5,5)
tnc_ssDNAcatalogue <- paste0(substr(rownames(ssDNAcatalogue),1,1),ref_ssDNAcatalogue,substr(rownames(ssDNAcatalogue),7,7))
ref_ssDNAcataloguecollapsed <- substr(rownames(ssDNAcataloguecollapsed),3,3)
alt_ssDNAcataloguecollapsed <- substr(rownames(ssDNAcataloguecollapsed),5,5)
tnc_ssDNAcataloguecollapsed <- paste0(substr(rownames(ssDNAcataloguecollapsed),1,1),ref_ssDNAcataloguecollapsed,substr(rownames(ssDNAcataloguecollapsed),7,7))

ref_dsDNAcatalogue <- substr(rownames(dsDNAcatalogue),3,3)
alt_dsDNAcatalogue <- substr(rownames(dsDNAcatalogue),5,5)
tnc_dsDNAcatalogue <- paste0(substr(rownames(dsDNAcatalogue),1,1),ref_dsDNAcatalogue,substr(rownames(dsDNAcatalogue),7,7))
ref_dsDNAcataloguecollapsed <- substr(rownames(dsDNAcataloguecollapsed),3,3)
alt_dsDNAcataloguecollapsed <- substr(rownames(dsDNAcataloguecollapsed),5,5)
tnc_dsDNAcataloguecollapsed <- paste0(substr(rownames(dsDNAcataloguecollapsed),1,1),ref_dsDNAcataloguecollapsed,substr(rownames(dsDNAcataloguecollapsed),7,7))

for(i in colnames(ssDNAcatalogue)){
	ssDNAdf <- rbind(ssDNAdf,cbind(rep(i,sum(ssDNAcatalogue[,i])),rep(ref_ssDNAcatalogue,ssDNAcatalogue[,i]),rep(alt_ssDNAcatalogue,ssDNAcatalogue[,i]),rep(tnc_ssDNAcatalogue,ssDNAcatalogue[,i])))
	ssDNAdfcollapsed <- rbind(ssDNAdfcollapsed,cbind(rep(i,sum(ssDNAcataloguecollapsed[,i])),rep(ref_ssDNAcataloguecollapsed,ssDNAcataloguecollapsed[,i]),rep(alt_ssDNAcataloguecollapsed,ssDNAcataloguecollapsed[,i]),rep(tnc_ssDNAcataloguecollapsed,ssDNAcataloguecollapsed[,i])))
	dsDNAdf <- rbind(dsDNAdf,cbind(rep(i,sum(dsDNAcatalogue[,i])),rep(ref_dsDNAcatalogue,dsDNAcatalogue[,i]),rep(alt_dsDNAcatalogue,dsDNAcatalogue[,i]),rep(tnc_dsDNAcatalogue,dsDNAcatalogue[,i])))
	dsDNAdfcollapsed <- rbind(dsDNAdfcollapsed,cbind(rep(i,sum(dsDNAcataloguecollapsed[,i])),rep(ref_dsDNAcataloguecollapsed,dsDNAcataloguecollapsed[,i]),rep(alt_dsDNAcataloguecollapsed,dsDNAcataloguecollapsed[,i]),rep(tnc_dsDNAcataloguecollapsed,dsDNAcataloguecollapsed[,i])))
}

ssDNAdf <- as.data.frame(ssDNAdf) 
dsDNAdf <- as.data.frame(dsDNAdf)
ssDNAdfcollapsed <- as.data.frame(ssDNAdfcollapsed)
dsDNAdfcollapsed <- as.data.frame(dsDNAdfcollapsed)
colnames(ssDNAdf) <- c("Sample","Ref","Alt","Trinuc")
colnames(dsDNAdf) <- c("Sample","Ref","Alt","Trinuc")
colnames(ssDNAdfcollapsed) <- c("Sample","Ref","Alt","Trinuc")
colnames(dsDNAdfcollapsed) <- c("Sample","Ref","Alt","Trinuc")

cat("DONE\n")

##Creating MutationalPatterns GRanges
cat("#### Creating MutationalPatterns GRanges...")
ssDNAgranges <- GRangesList()
dsDNAgranges <- GRangesList()
for(i in unique(samplenames)){
	sstemp.granges <- c(mutations[[i]]$fwd,mutations[[i]]$rev)
	dstemp.granges <- mutations[[i]]$dsDNA
	
	if(length(sstemp.granges)>0){
  	names(sstemp.granges) <- paste0(seqnames(sstemp.granges),":",start(sstemp.granges),"_",sstemp.granges$rn,"/",sstemp.granges$qn)
  	sstemp.granges$paramRangeID <- factor(NA)
  	sstemp.granges$REF <- sstemp.granges$rn
  	sstemp.granges$REF <- DNAStringSet(sstemp.granges$REF)
  	sstemp.granges$ALT <- DNAStringSetList(as.list(sstemp.granges$qn))
  	sstemp.granges$rn <- NULL
  	sstemp.granges$qn <- NULL
  	sstemp.granges$QUAL <- 10
  	sstemp.granges$FILTER <- "PASS"
  	sstemp.granges$GS <- strand(sstemp.granges)
  	strand(sstemp.granges) <- as.character("*")
  	ssDNAgranges[[i]] <- sstemp.granges
	}
	
	if(length(dstemp.granges)>0){
	  names(dstemp.granges) <- paste0(seqnames(dstemp.granges),":",start(dstemp.granges),"_",dstemp.granges$rn,"/",dstemp.granges$qn)
  	dstemp.granges$paramRangeID <- factor(NA)
  	dstemp.granges$REF <- dstemp.granges$rn
  	dstemp.granges$REF <- DNAStringSet(dstemp.granges$REF)
  	dstemp.granges$ALT <- DNAStringSetList(as.list(dstemp.granges$qn))
  	dstemp.granges$rn <- NULL
  	dstemp.granges$qn <- NULL
  	dstemp.granges$QUAL <- 10
  	dstemp.granges$FILTER <- "PASS"
  	dstemp.granges$GS <- strand(dstemp.granges)
  	strand(dstemp.granges) <- as.character("*")
  	dsDNAgranges[[i]] <- dstemp.granges
	}
	
}
cat("DONE\n")

cat("#### Saving to RDS...")
output <- list()
output[["tissues"]] <- tissues
output[["mutations"]] <- mutations
output[["trinucleotide_contexts"]] <- trinucleotide_contexts
output[["signature.tools.lib"]][["ssDNA"]] <- ssDNAcatalogue
output[["signature.tools.lib"]][["dsDNA"]] <- dsDNAcatalogue
output[["signature.tools.libcollapsed"]][["ssDNA"]] <- ssDNAcataloguecollapsed
output[["signature.tools.libcollapsed"]][["dsDNA"]] <- dsDNAcataloguecollapsed
output[["SigFit"]][["ssDNA"]] <- ssDNAdf
output[["SigFit"]][["dsDNA"]] <- dsDNAdf
output[["SigFitcollapsed"]][["ssDNA"]] <- ssDNAdfcollapsed
output[["SigFitcollapsed"]][["dsDNA"]] <- dsDNAdfcollapsed
output[["MutationalPatterns"]][["ssDNA"]] <- ssDNAgranges
output[["MutationalPatterns"]][["dsDNA"]] <- dsDNAgranges

qsave(output,paste0(yaml.config$analysisoutput_path,"/mutation_frequencies.RDS"))
cat("DONE\n")
sink()
