#!/usr/bin/env Rscript
# Usage: mutation_filtering.R [configuration.yaml]
#  - Uses [configuration.yaml] file per documentation
#  - Outputs to analysisoutput_path directory defined in configuration.yaml:
#     1. bamfilezmw_filtered RDS file that saves the final list of mutations and additional data.
#     2. mutation_filtering.output.[timestamp].txt recording pipeline output.
#
# Prerequisites:
#  - R packages: GenomicAlignments, GenomicRanges, vcfR, Rsamtools, plyr, dplyr, plyranges, stringr, rtracklayer, configr, qs, BSgenome, BSgenome package for the reference genome specified in the configuration (e.g., BSgenome.Hsapiens.T2T.CHM13v1)
#  - Conda
#  - wiggletools (v1.2.11) in system path
#  - Conda environment installed with PacBio tools: pbmm2, zmwfilter
#  - bcftools
#  - samtools
#  - wigToBigWig (from Kent tools)
#  - SLURM for converting large wig to bigwig files in bw.granges.thresholding function
#

#Stop script for any warnings
options(warn=2)

######
#Load required libraries
######
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs))

######
#Load configuration
######
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=1) {
  stop("Usage: Rscript mutation_filtering.R [configuration.yaml]\n", call.=FALSE)
}
cat("####### Loading configuration...")
yaml.config <- suppressWarnings(read.config(args[1]))

# Configure environment
condabase_script <- yaml.config$condabase_script
conda_pbbioconda_env <- yaml.config$conda_pbbioconda_env
wigToBigWig_bin <- yaml.config$wigToBigWig_bin
bcftools_bin <- yaml.config$bcftools_bin
samtools_bin <- yaml.config$samtools_bin

#Load BSgenome package
suppressPackageStartupMessages(library(yaml.config$BSgenomepackagename,character.only=TRUE))

#Configure bigwig cache directory
bigwig_cachedir <- yaml.config$bigwig_cachedir

#Configure bamfilezmw_all and subreads bam input filenames. Note: subreads bam is used for extraction and alignment of subreads with mutations.
bamfilezmw_all_filename <- yaml.config$bamfilezmw_all_filename
subreads_filename <- yaml.config$subreads_filename

#Configure output path and output base name (i.e. without file extension).
analysisoutput_path <- yaml.config$analysisoutput_path
analysisoutput_basename <- yaml.config$analysisoutput_basename

#Configure sample names
bamfilezmw_samplenames <- as.character(unlist(yaml.config$samplenames))

#Configure genome reference information and files
genome <- yaml.config$genome
referenceinfo <- data.frame(chrsizes= yaml.config$chrsizes,Nrefbed=yaml.config$Nrefbed,fastaref=yaml.config$fastaref,genomemmiindex=yaml.config$genomemmiindex)
chroms_to_analyze <- read.table(yaml.config$chrs_to_analyze)[,1]

 ###################
 ## Basic filters ##
thresholds <- data.frame(rqfwd=numeric(),rqrev=numeric(),rqavg=numeric(),
                         ecfwd=numeric(),ecrev=numeric(),ecavg=numeric(),
                         mapqfwd=numeric(),mapqrev=numeric(),mapqavg=numeric(),
                         numsnvsfwd=numeric(),numsnvsrev=numeric(),numsnvszmw=numeric(),
                         numsnvsfwdrevdiff=numeric(),
                         numindelsfwd=numeric(),numindelsrev=numeric(),numindelszmw=numeric(),
                         numsnvsfwdpostVCF=numeric(),numsnvsrevpostVCF=numeric(),numsnvszmwpostVCF=numeric(),
                         numsoftclipfwd=numeric(),numsoftcliprev=numeric(),numsoftclipavg=numeric(),
                         minqq=numeric(),
                         bpends=numeric(),
                         ccsindelfilter=logical(),ccsindelinspad=character(),ccsindeldelpad=character(),
                         minsubreads_cvg_fraction=numeric(),
                         minZMWsubreadsVariantReads=numeric(),minZMWsubreadsVAF=numeric(),
                         minssDNAsubreadsVariantReads=numeric(),minssDNAsubreadsVAF=numeric())

for(i in names(yaml.config$thresholds)){
  for(j in names(yaml.config$thresholds[[i]])){
    thresholds[as.numeric(i),j] <- yaml.config$thresholds[[i]][[j]]
  }
}


###### *** NEW CODE for loading VCF filtering config ***###
#Load germline VCF filter parameters
vcf_config <- yaml.config$germline_vcf_filtergroups %>%
  bind_rows %>%
  group_by(across(-c(vcf_SNV_FILTERS,vcf_INDEL_FILTERS))) %>%
  summarize(vcf_SNV_FILTERS=list(vcf_SNV_FILTERS),vcf_INDEL_FILTERS=list(vcf_INDEL_FILTERS), .groups="drop")

 ################
 ## vcffilters ##
vcffilters <- list()
vcffilters.thresholds.emptydf <- data.frame(vcfSNVfilter=logical(),vcfSNVDepth=numeric(),vcfSNVGQ=numeric(),vcfSNVVAF=numeric(),vcfSNVQUAL=numeric(),vcfINDELfilter=logical(),vcfINDELDepth=numeric(),vcfINDELGQ=numeric(),vcfINDELVAF=numeric(),vcfINDELQUAL=numeric(),vcfINDELinspad=character(),vcfINDELdelpad=character())

for(i in names(yaml.config$vcffilters)){
  vcffilters[[as.numeric(i)]] <- list()
  
  for(j in names(yaml.config$vcffilters[[i]])){
    for(k in names(yaml.config$vcffilters[[i]][[j]])){
      vcffilters[[as.numeric(i)]][[j]][[as.numeric(k)]] <- yaml.config$vcffilters[[i]][[j]][[k]]
    }
  }

  if(yaml.config$vcffilters_sameforallbasicfiltersets==TRUE && nrow(thresholds)>1){
    for(q in 2:nrow(thresholds)){
      vcffilters[[q]] <- vcffilters[[1]]
    }
    break
  }
}

 ###################
 ## bigwigfilters ##
bigwigfilters <- list()
for(i in names(yaml.config$bigwigfilters)){
  bigwigfilters[[as.numeric(i)]] <- list()
  
  for(j in names(yaml.config$bigwigfilters[[i]])){
    bigwigfilters[[as.numeric(i)]][[as.numeric(j)]] <- yaml.config$bigwigfilters[[i]][[j]]
  }
  
  if(yaml.config$bigwigfilters_sameforallbasicfiltersets==TRUE && nrow(thresholds)>1){
    for(q in 2:nrow(thresholds)){
      bigwigfilters[[q]] <- bigwigfilters[[1]]
    }
    break
  }
}


 ################
 ## bamfilters ##
bamfilters <- list()
for(i in names(yaml.config$bamfilters)){
  bamfilters[[as.numeric(i)]] <- list()
  
  for(j in names(yaml.config$bamfilters[[i]])){
    bamfilters[[as.numeric(i)]][[j]] <- yaml.config$bamfilters[[i]][[j]]
  }
  
  if(yaml.config$bamfilters_sameforallbasicfiltersets==TRUE && nrow(thresholds)>1){
    for(q in 2:nrow(thresholds)){
      bamfilters[[q]] <- bamfilters[[1]]
    }
    break
  }
}

cat("DONE\n\n")



######
#Define analysis function
######
analyzebamfilezmw_all <- function(bamfilezmw_all_filename,thresholds,vcffilters,bigwigfilters,bamfilters,referenceinfo,outputpath,subreadsbam,bw_cachedir,chrs_to_analyze){

#Function that returns a list containing two Granges objects of regions above or equal, and regions below, a threshold value calculated as the mean of a wiggle file in a specified bin size (in bp) across the genome. Uses wiggletools. Saves result (bigwig) in bwcachedir for later viewing in IGV and to save time in future analyses.
bw.granges.thresholding <- function(bigwigfile,genomechrsizestsv,binsize,threshold,bwcachedir=bw_cachedir){
  
  #Generate temp filename and bwcachedir
  tmpoutput <- tempfile()
  if(!dir.exists(bwcachedir)){dir.create(bwcachedir,recursive=TRUE)}
  cachefile <- paste0(bwcachedir,"/",basename(bigwigfile),".bin",binsize)
  
  if(file.size(paste0(cachefile,".bw"))==0 | is.na(file.size(paste0(cachefile,".bw")))){
    cacheexists <- FALSE
    
    #Make genome BED file to use to fill in zero values for regions not in bigwig.
    system(paste("/bin/bash -c",shQuote(paste("awk '{print $1 \"\\t0\\t\" $2}'", genomechrsizestsv,"| sort -k1,1 -k2,2n >",paste0(tmpoutput,".genome.bed")))))
    
    #Run wiggletools
    if(binsize==1){
      system(paste("/bin/bash -c",shQuote(paste("wiggletools trim",paste0(tmpoutput,".genome.bed"),"fillIn",paste0(tmpoutput,".genome.bed"),bigwigfile,">",paste0(tmpoutput,".wig")))))
    }else{
      system(paste("/bin/bash -c",shQuote(paste("wiggletools trim",paste0(tmpoutput,".genome.bed"),"fillIn",paste0(tmpoutput,".genome.bed"), "scale",1/binsize,"bin",binsize,bigwigfile,">",paste0(tmpoutput,".wig")))))
    }
    
    #Convert wig to bigwig for later viewing in IGV. Save in bwcachedir. If wig file is large, run it as a SLURM job.
    if(file.size(paste0(tmpoutput,".wig"))>5000000000){
    	cat("Converting wig to bigwig as SLURM job...")
    	system(paste("mv",paste0(tmpoutput,".wig"),paste0(cachefile,".wig"))) #Move wig file to scratch space, since otherwise the SLURM job won't be able to access it.
    	cmd <- tempfile(pattern=".",tmpdir=getwd())
    	cat(paste(wigToBigWig_bin,paste0(cachefile,".wig"),genomechrsizestsv,paste0(cachefile,".bw")),file=cmd)
    	Sys.chmod(cmd,mode="700")
    	jobindex <- system(paste("sbatch",yaml.config$slurm_add_options,"-t 1920 --mem=64G --parsable -o /dev/null --wrap \"$HIDEF_SINGULARITY_WRAPPER",cmd,"; rm",cmd,"\""),intern=T)
    	
    	#Wait for job to finish
    	while(!(system(paste("sacct -o State -j",jobindex,"| tail -n 1 | awk '{print $1}'"),intern=T) %in% c("COMPLETED","FAILED"))){
        Sys.sleep(60) #Check SLURM job status every 1 minute
      }
      if(system(paste("sacct -o State -j",jobindex,"| tail -n 1 | awk '{print $1}'"),intern=T)=="FAILED"){
        cat("\n...SLURM job failed!")
      	system(paste("rm",paste0(cachefile,".wig")),ignore.stderr=TRUE)
        quit()
      }
    	system(paste("rm",paste0(cachefile,".wig")),ignore.stderr=TRUE)

    }else{
    	system(paste(wigToBigWig_bin,paste0(tmpoutput,".wig"),genomechrsizestsv,paste0(cachefile,".bw")))
    }
    
  }else{
    cacheexists <- TRUE
    cat("Loading from cache...")
  }

  system(paste("/bin/bash -c",shQuote(paste("wiggletools write_bg - compress gte",threshold,paste0(cachefile,".bw"),"| cut -f 1-3 >",paste0(tmpoutput,".gte")))))
  system(paste("/bin/bash -c",shQuote(paste("wiggletools write_bg - compress lt",threshold,paste0(cachefile,".bw"),"| cut -f 1-3 >",paste0(tmpoutput,".lt")))))
  
  #Make list of results
  listresult <- list(gte=import(paste0(tmpoutput,".gte"),format="BED"),lt=import(paste0(tmpoutput,".lt"),format="BED"))
  
  #Delete temp files
  system(paste("rm",paste0(tmpoutput,".gte"),"; rm",paste0(tmpoutput,".lt")))
  if(!cacheexists){system(paste("rm",paste0(tmpoutput,".genome.bed"),"; rm",paste0(tmpoutput,".wig")),ignore.stderr=TRUE)}
  
  return(listresult)
}

#Function to filter granges for specified chromosomes vector
filter.granges.by.chroms <- function(gr,chroms){
  return(gr[seqnames(gr) %in% chroms])
}

if(!dir.exists(outputpath)){dir.create(outputpath,recursive=TRUE)}

#Duplicate output to text file
sink(paste0(outputpath,"/mutation_filtering.output.",format(Sys.time(),"%Y-%m-%d_%H%M%OS2"),".txt"),split=TRUE)

cat("####### Mutation Filtering Pipeline #######\n\n")

#Load bamfilezmw_all
cat("#### Loading bamfilezmw_all file...")
bamfilezmw_all <- qread(bamfilezmw_all_filename)
cat("DONE\n\n")

bamfilezmw_filtered <- list()

for(i in names(bamfilezmw_all)){
  for(j in names(bamfilezmw_all[[i]])){
    
    cat("\n#### Analyzing sample:",i,j,"####\n")
    
  	if(class(bamfilezmw_all[[i]][[j]]) != "list"){
  		#Skip this sample
  		cat("     No reads or mutations in sample in selected chromosomes. Proceeding to next sample...\n")
  		next
  	}
  	
    #Variables for basic statistics calculations. genome.granges is used to calculate regional filtering stats, excluding 'N' sequence regions. fwd/rev/zmw.granges is used only for outputting statistics on bases filtered by base quality and bpends filters.
    totalnumzmws <- length(bamfilezmw_all[[i]][[j]]$tag.fwd$rq)
    genome.granges <- as.data.frame(read.table(referenceinfo$chrsizes))
    genome.granges$V3 <- rep(1,nrow(genome.granges))
    genome.granges <- makeGRangesFromDataFrame(genome.granges,seqnames.field="V1",start.field="V3",end.field="V2")
    genome.granges <- filter.granges.by.chroms(genome.granges,chrs_to_analyze)
    
    fwd.granges <- makeGRangesFromDataFrame(data.frame(chrom=bamfilezmw_all[[i]][[j]]$rname.fwd,start=bamfilezmw_all[[i]][[j]]$pos.fwd,end=bamfilezmw_all[[i]][[j]]$pos.fwd+bamfilezmw_all[[i]][[j]]$isize.fwd-1))
    rev.granges <- makeGRangesFromDataFrame(data.frame(chrom=bamfilezmw_all[[i]][[j]]$rname.rev,start=bamfilezmw_all[[i]][[j]]$pos.rev,end=bamfilezmw_all[[i]][[j]]$pos.rev+bamfilezmw_all[[i]][[j]]$isize.rev-1))
    zmw.granges <- makeGRangesFromDataFrame(data.frame(chrom=seqnames(fwd.granges),start=apply(cbind(start(fwd.granges),start(rev.granges)),1,max),end=apply(cbind(end(fwd.granges),end(rev.granges)),1,min)))

    for(k in 1:nrow(thresholds)){
      cat("\n ## Filter set",k,"\n")
      
      #Set up list to store filtering results
      bamfilezmw_filtered[[i]][[j]][[k]] <- list()
      
      cat(" -> Basic filter settings:\n")
      print(thresholds[k,],row.names=FALSE)
    
      #Lists to store read regional filters and genome regional filters for use in final mutation frequency calculation. Note: for base quality (qq), bpends, and ccsindel filters in read regional filters, add a column 'appliestoread' that species the full name (qname) of the consensus read that it applies to, since these filters are read-specific and we don't want to filter out these locations across all reads that overlap the filtered region.
      readregionalfilters.fwd <- list()
      readregionalfilters.rev <- list()
      readregionalfilters.zmw <- list()
      genomeregionalfilters.fwd <- list()
      genomeregionalfilters.rev <- list()
      genomeregionalfilters.zmw <- list()
      
      ##Apply ZMW filters
      cat(" Applying basic ZMW filters...")
      
      #Calculate average rq, ec, mapq, softclipped bases
      rq.avg <- apply(cbind(bamfilezmw_all[[i]][[j]]$tag.fwd$rq,bamfilezmw_all[[i]][[j]]$tag.rev$rq),1,mean)
      ec.avg <- apply(cbind(bamfilezmw_all[[i]][[j]]$tag.fwd$ec,bamfilezmw_all[[i]][[j]]$tag.rev$ec),1,mean)
      mapq.avg <- apply(cbind(bamfilezmw_all[[i]][[j]]$mapq.fwd,bamfilezmw_all[[i]][[j]]$mapq.rev),1,mean)
      softclip.fwd <- unlist(lapply(bamfilezmw_all[[i]][[j]]$cigar.fwd,function(x){sum(cigarOpTable(x)[,"S"])}))
      softclip.rev <- unlist(lapply(bamfilezmw_all[[i]][[j]]$cigar.rev,function(x){sum(cigarOpTable(x)[,"S"])}))
      softclip.avg <- apply(cbind(softclip.fwd,softclip.rev),1,mean)
      
      filteredzmws.rq <- bamfilezmw_all[[i]][[j]]$tag.fwd$rq < thresholds[k,"rqfwd"] | bamfilezmw_all[[i]][[j]]$tag.rev$rq < thresholds[k,"rqrev"] | rq.avg < thresholds[k,"rqavg"]
      filteredzmws.ec <- bamfilezmw_all[[i]][[j]]$tag.fwd$ec < thresholds[k,"ecfwd"] | bamfilezmw_all[[i]][[j]]$tag.rev$ec < thresholds[k,"ecrev"] | ec.avg < thresholds[k,"ecavg"]
      filteredzmws.mapq <- bamfilezmw_all[[i]][[j]]$mapq.fwd < thresholds[k,"mapqfwd"] | bamfilezmw_all[[i]][[j]]$mapq.rev < thresholds[k,"mapqrev"] | mapq.avg < thresholds[k,"mapqavg"]
      filteredzmws.numsnvs <- unlist(lapply(bamfilezmw_all[[i]][[j]]$tag.fwd$rp,length)) > thresholds[k,"numsnvsfwd"] | unlist(lapply(bamfilezmw_all[[i]][[j]]$tag.rev$rp,length)) > thresholds[k,"numsnvsrev"] | unlist(lapply(bamfilezmw_all[[i]][[j]]$zmw.rp,length)) > thresholds[k,"numsnvszmw"]
      
      numindels.fwd <- sapply(bamfilezmw_all[[i]][[j]]$cigar.fwd,function(x){str_count(x,"I|D")},USE.NAMES=F)
      numindels.rev <- sapply(bamfilezmw_all[[i]][[j]]$cigar.rev,function(x){str_count(x,"I|D")},USE.NAMES=F)
      numindels.zmw <- apply(cbind(numindels.fwd,numindels.rev),1,mean)
      filteredzmws.numindels <- numindels.fwd > thresholds[k,"numindelsfwd"] | numindels.rev > thresholds[k,"numindelsrev"] | numindels.zmw > thresholds[k,"numindelszmw"]
      
      filteredzmws.numsnvsfwdrevdiff <- abs(unlist(lapply(bamfilezmw_all[[i]][[j]]$tag.fwd$rp,length)) - unlist(lapply(bamfilezmw_all[[i]][[j]]$tag.rev$rp,length))) > thresholds[k,"numsnvsfwdrevdiff"]
      filteredzmws.numsoftclip <- softclip.fwd > thresholds[k,"numsoftclipfwd"] | softclip.rev > thresholds[k,"numsoftcliprev"] | softclip.avg > thresholds[k,"numsoftclipavg"]
      
      filteredzmws <- filteredzmws.rq | filteredzmws.ec | filteredzmws.mapq | filteredzmws.numsnvs | filteredzmws.numsnvsfwdrevdiff | filteredzmws.numindels | filteredzmws.numsoftclip
      includezmws <- !filteredzmws
      cat("DONE\n")
      
      if(!any(includezmws)){
      	bamfilezmw_filtered[[i]][[j]][[k]] <- list()
      	cat("\n     No ZMWs remain. Proceeding to next sample/filterset!\n\n")
      	next
      }
      
      cat(paste("        'rq' filtered:",length(which(filteredzmws.rq)),"/",totalnumzmws,"total ZMWs =",round(length(which(filteredzmws.rq))/totalnumzmws*100,1),"% of total ZMWs\n"))
      cat(paste("        'ec' filtered:",length(which(filteredzmws.ec)),"/",totalnumzmws,"total ZMWs =",round(length(which(filteredzmws.ec))/totalnumzmws*100,1),"% of total ZMWs\n"))
      cat(paste("        'mapq' filtered:",length(which(filteredzmws.mapq)),"/",totalnumzmws,"total ZMWs =",round(length(which(filteredzmws.mapq))/totalnumzmws*100,1),"% of total ZMWs\n"))
      cat(paste("        'numsnvs' filtered:",length(which(filteredzmws.numsnvs)),"/",totalnumzmws,"total ZMWs =",round(length(which(filteredzmws.numsnvs))/totalnumzmws*100,1),"% of total ZMWs\n"))
      cat(paste("        'numsnvsfwdrevdiff' filtered:",length(which(filteredzmws.numsnvsfwdrevdiff)),"/",totalnumzmws,"total ZMWs =",round(length(which(filteredzmws.numsnvsfwdrevdiff))/totalnumzmws*100,1),"% of total ZMWs\n"))
      cat(paste("        'numindels' filtered:",length(which(filteredzmws.numindels)),"/",totalnumzmws,"total ZMWs =",round(length(which(filteredzmws.numindels))/totalnumzmws*100,1),"% of total ZMWs\n"))
      cat(paste("        'numsoftclip' filtered:",length(which(filteredzmws.numsoftclip)),"/",totalnumzmws,"total ZMWs =",round(length(which(filteredzmws.numsoftclip))/totalnumzmws*100,1),"% of total ZMWs\n"))
      cat(paste("        Total filtered:",length(which(filteredzmws)),"/",totalnumzmws,"total ZMWs =",round(length(which(filteredzmws))/totalnumzmws*100,1),"% of total ZMWs\n"))
      cat(paste("        Remaining:",length(which(includezmws)),"ZMWs =",round(length(which(includezmws))/totalnumzmws*100,1),"% of total ZMWs\n"))
      
      #Clear temporary variables
      rm(rq.avg,ec.avg,mapq.avg,softclip.fwd,softclip.rev,softclip.avg,filteredzmws.rq,filteredzmws.ec,filteredzmws.mapq,filteredzmws.numsnvs,numindels.fwd,numindels.rev,numindels.zmw,filteredzmws.numsnvsfwdrevdiff,filteredzmws)
      
      ##Setup mutation filters to all true
      includemutations.fwd <- lapply(bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws],function(x){x | TRUE})
      includemutations.rev <- lapply(bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws],function(x){x | TRUE})
      includemutations.zmw <- lapply(bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws],function(x){x | TRUE})
      
      ##Apply noN filter
      cat(" Applying noN filters...")
      includemutations.fwd <- mapply(function(x,y){x & y},includemutations.fwd,bamfilezmw_all[[i]][[j]]$tag.fwd$noN[includezmws],SIMPLIFY=FALSE)
      includemutations.rev <- mapply(function(x,y){x & y},includemutations.rev,bamfilezmw_all[[i]][[j]]$tag.rev$noN[includezmws],SIMPLIFY=FALSE)
      includemutations.zmw <- mapply(function(x,y){x & y},includemutations.zmw,bamfilezmw_all[[i]][[j]]$zmw.noN[includezmws],SIMPLIFY=FALSE)
      cat("DONE\n")
      
      #Make grange of N regions and output filtering statistics. Also, subtract N regions from genome.granges.
      Nref.granges <- import(referenceinfo$Nrefbed,format="BED")
      Nref.granges <- filter.granges.by.chroms(Nref.granges,chrs_to_analyze)
      
      cat(paste("        'N' regions:",sum(width(Nref.granges)),"/",sum(width(genome.granges)),"bases =",round(sum(width(Nref.granges))/sum(width(genome.granges))*100,1),"% of the genome\n"))
      readregionalfilters.fwd <- c(readregionalfilters.fwd,Nref.granges)
      readregionalfilters.rev <- c(readregionalfilters.rev,Nref.granges)
      readregionalfilters.zmw <- c(readregionalfilters.zmw,Nref.granges)
      genomeregionalfilters.fwd <- c(genomeregionalfilters.fwd,Nref.granges)
      genomeregionalfilters.rev <- c(genomeregionalfilters.rev,Nref.granges)
      genomeregionalfilters.zmw <- c(genomeregionalfilters.zmw,Nref.granges)
      genome.noN.granges <- GenomicRanges::setdiff(genome.granges,Nref.granges)
      
      #Create genome.granges.filtered to track subsequent genome regions filtering for dsDNA mutations.
      genome.granges.filtered <- genome.noN.granges
      genome.granges.filtered.ssDNA <- genome.noN.granges
      
      ##Apply germline SNV filters
      cat(" Applying VCF SNV filters...\n")
      for(vcffileindex in 1:length(vcffilters[[k]][[i]])){
        if(vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVfilter){
          vcffile <- vcffilters[[k]][[i]][[vcffileindex]]$vcffilename
          cat("  VCF file:",basename(vcffile),"...")
            
          filterstoinclude <- paste0(paste0("^",vcffilters[[k]][[i]][[vcffileindex]]$SNVFILTERS,"$"),collapse="|")
          
          includemutations.fwd <- mapply(function(x,a,b,c,d,e){x & !(a>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVDepth & b>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVGQ & c>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVVAF & d>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVQUAL & grepl(filterstoinclude,e))},includemutations.fwd,bamfilezmw_all[[i]][[j]]$vcfSNV.fwd[[vcffile]]$Depth[includezmws],bamfilezmw_all[[i]][[j]]$vcfSNV.fwd[[vcffile]]$GQ[includezmws],bamfilezmw_all[[i]][[j]]$vcfSNV.fwd[[vcffile]]$VAF[includezmws],bamfilezmw_all[[i]][[j]]$vcfSNV.fwd[[vcffile]]$QUAL[includezmws],bamfilezmw_all[[i]][[j]]$vcfSNV.fwd[[vcffile]]$FILTER[includezmws],SIMPLIFY=FALSE)
      
          includemutations.rev <- mapply(function(x,a,b,c,d,e){x & !(a>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVDepth & b>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVGQ & c>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVVAF & d>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVQUAL & grepl(filterstoinclude,e))},includemutations.rev,bamfilezmw_all[[i]][[j]]$vcfSNV.rev[[vcffile]]$Depth[includezmws],bamfilezmw_all[[i]][[j]]$vcfSNV.rev[[vcffile]]$GQ[includezmws],bamfilezmw_all[[i]][[j]]$vcfSNV.rev[[vcffile]]$VAF[includezmws],bamfilezmw_all[[i]][[j]]$vcfSNV.rev[[vcffile]]$QUAL[includezmws],bamfilezmw_all[[i]][[j]]$vcfSNV.rev[[vcffile]]$FILTER[includezmws],SIMPLIFY=FALSE)
          
          includemutations.zmw <- mapply(function(x,a,b,c,d,e){x & !(a>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVDepth & b>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVGQ & c>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVVAF & d>=vcffilters[[k]][[i]][[vcffileindex]]$vcfSNVQUAL & grepl(filterstoinclude,e))},includemutations.zmw,bamfilezmw_all[[i]][[j]]$zmw.vcfSNV[[vcffile]]$Depth[includezmws],bamfilezmw_all[[i]][[j]]$zmw.vcfSNV[[vcffile]]$GQ[includezmws],bamfilezmw_all[[i]][[j]]$zmw.vcfSNV[[vcffile]]$VAF[includezmws],bamfilezmw_all[[i]][[j]]$zmw.vcfSNV[[vcffile]]$QUAL[includezmws],bamfilezmw_all[[i]][[j]]$zmw.vcfSNV[[vcffile]]$FILTER[includezmws],SIMPLIFY=FALSE)
          
          ##Set variants that equal NA after filtering to TRUE, since they did not meet germline filtering.
          includemutations.fwd <- lapply(includemutations.fwd,function(x){x==TRUE | is.na(x)})
          includemutations.rev <- lapply(includemutations.rev,function(x){x==TRUE | is.na(x)})
          includemutations.zmw <- lapply(includemutations.zmw,function(x){x==TRUE | is.na(x)})
      
          cat("DONE\n")
        }else{
          cat("NOT USED\n")
        }
      }
      cat("  VCF SNV filters...DONE\n")
      
      #Number of postVCF SNVs ZMW filter
      cat(" Applying postVCF SNV ZMW filters...")
      zmwstokeep <- unlist(lapply(includemutations.fwd,function(x){length(which(x))})) <= thresholds[k,"numsnvsfwdpostVCF"] & unlist(lapply(includemutations.rev,function(x){length(which(x))})) <= thresholds[k,"numsnvsrevpostVCF"] & unlist(lapply(includemutations.zmw,function(x){length(which(x))})) <= thresholds[k,"numsnvszmwpostVCF"]
      includezmws[includezmws==TRUE] <- zmwstokeep
      includemutations.fwd <- includemutations.fwd[zmwstokeep]
      includemutations.rev <- includemutations.rev[zmwstokeep]
      includemutations.zmw <- includemutations.zmw[zmwstokeep]

      cat("DONE\n")
      
      if(!any(includezmws)){
      	bamfilezmw_filtered[[i]][[j]][[k]] <- list()
      	cat("\n     No ZMWs remain. Proceeding to next sample/filterset!\n\n")
      	next
      }
      
      #Output statistics
      cat(paste("        Filtered:",length(which(!zmwstokeep)),"/",totalnumzmws,"additional ZMWs =",round(length(which(!zmwstokeep))/totalnumzmws*100,1),"% of ZMWs\n"))
      cat(paste("        Remaining:",length(which(includezmws)),"ZMWs =",round(length(which(includezmws))/totalnumzmws*100,1),"% of ZMWs\n"))

      #ZMW Bigwig filters
      if(any(unlist(lapply(bigwigfilters[[k]],function(x){x$bigwigfiltertype=="ZMW"})))){
        cat(" Applying ZMW Bigwig filters...\n")
        for(bigwigfileindex in which(unlist(lapply(bigwigfilters[[k]],function(x){x$bigwigfiltertype=="ZMW"})))){
          cat(paste0("  Bigwig file: ",bigwigfilters[[k]][[bigwigfileindex]]$bigwigfile,"..."))
          
          #Extract bigwig threshold
          if(grepl("^gte",bigwigfilters[[k]][[bigwigfileindex]]$threshold)){
            bigwigthresholdtype <- "gte"
            bigwigthreshold <- as.numeric(sub("^gte","",bigwigfilters[[k]][[bigwigfileindex]]$threshold))
          }else if(grepl("^lt",bigwigfilters[[k]][[bigwigfileindex]]$threshold)){
            bigwigthresholdtype <- "lt"
            bigwigthreshold <- as.numeric(sub("^lt","",bigwigfilters[[k]][[bigwigfileindex]]$threshold))
          }else{
            stop(paste("Wrong format of bigwigfilter",bigwigfilters[[k]][[bigwigfileindex]]$threshold))
          }
          
          #Extract ZMW threshold
          if(grepl("^gte",bigwigfilters[[k]][[bigwigfileindex]]$zmwthreshold)){
            zmwthresholdtype <- "gte"
            zmwthreshold <- as.numeric(sub("^gte","",bigwigfilters[[k]][[bigwigfileindex]]$zmwthreshold))
          }else if(grepl("^lt",bigwigfilters[[k]][[bigwigfileindex]]$zmwthreshold)){
            zmwthresholdtype <- "lt"
            zmwthreshold <- as.numeric(sub("^lt","",bigwigfilters[[k]][[bigwigfileindex]]$zmwthreshold))
          }else{
            stop(paste("Wrong format of ZMW bigwigfilter",bigwigfilters[[k]][[bigwigfileindex]]$zmwthreshold))
          }
          
          #Make bigwig and ccs granges
          bigwig.granges <- bw.granges.thresholding(bigwigfilters[[k]][[bigwigfileindex]]$bigwigfile,referenceinfo$chrsizes,bigwigfilters[[k]][[bigwigfileindex]]$binsize,bigwigthreshold)
          bigwig.granges <- bigwig.granges[[bigwigthresholdtype]]
          bigwig.granges <- filter.granges.by.chroms(bigwig.granges,chrs_to_analyze)
          
          fwdccs.granges <- makeGRangesFromDataFrame(data.frame(chrom=bamfilezmw_all[[i]][[j]]$rname.fwd[includezmws],start=bamfilezmw_all[[i]][[j]]$pos.fwd[includezmws],end=bamfilezmw_all[[i]][[j]]$pos.fwd[includezmws]+bamfilezmw_all[[i]][[j]]$isize.fwd[includezmws]-1))
          revccs.granges <- makeGRangesFromDataFrame(data.frame(chrom=bamfilezmw_all[[i]][[j]]$rname.rev[includezmws],start=bamfilezmw_all[[i]][[j]]$pos.rev[includezmws],end=bamfilezmw_all[[i]][[j]]$pos.rev[includezmws]+bamfilezmw_all[[i]][[j]]$isize.rev[includezmws]-1))
          
          #Find fraction overlap
          fwdccs.bw.overlaps <- findOverlaps(fwdccs.granges,bigwig.granges)
          revccs.bw.overlaps <- findOverlaps(revccs.granges,bigwig.granges)
            
          fwdccs.bw.intersect <- pintersect(fwdccs.granges[queryHits(fwdccs.bw.overlaps)], bigwig.granges[subjectHits(fwdccs.bw.overlaps)])
          revccs.bw.intersect <- pintersect(revccs.granges[queryHits(revccs.bw.overlaps)], bigwig.granges[subjectHits(revccs.bw.overlaps)])
            
          options(dplyr.summarise.inform=F)
          fwdccs.bw.intersect.fraction <- width(fwdccs.bw.intersect)/width(fwdccs.granges[queryHits(fwdccs.bw.overlaps)])
          fwdccs.bw.intersect.fraction.df <- data.frame(queryHits=queryHits(fwdccs.bw.overlaps),fraction=fwdccs.bw.intersect.fraction)
          fwdccs.bw.intersect.fraction.df.summary <- fwdccs.bw.intersect.fraction.df %>% group_by(queryHits) %>% summarise("sum_fraction"  = sum(fraction))
          
          revccs.bw.intersect.fraction <- width(revccs.bw.intersect)/width(revccs.granges[queryHits(revccs.bw.overlaps)])
          revccs.bw.intersect.fraction.df <- data.frame(queryHits=queryHits(revccs.bw.overlaps),fraction=revccs.bw.intersect.fraction)
          revccs.bw.intersect.fraction.df.summary <- revccs.bw.intersect.fraction.df %>% group_by(queryHits) %>% summarise("sum_fraction"  = sum(fraction))
            
          mcols(fwdccs.granges)$fractionfiltered <- rep(0,length(fwdccs.granges))
          mcols(fwdccs.granges)$fractionfiltered[fwdccs.bw.intersect.fraction.df.summary$queryHits] <- fwdccs.bw.intersect.fraction.df.summary$sum_fraction
          
          mcols(revccs.granges)$fractionfiltered <- rep(0,length(revccs.granges))
          mcols(revccs.granges)$fractionfiltered[revccs.bw.intersect.fraction.df.summary$queryHits] <- revccs.bw.intersect.fraction.df.summary$sum_fraction
          
          meanfractionfiltered <- apply(cbind(fwdccs.granges$fractionfiltered,revccs.granges$fractionfiltered),1,mean)
          
          #Perform filtering
          if(zmwthresholdtype == "gte"){
            zmwstokeep <- meanfractionfiltered<zmwthreshold
          }else if(zmwthresholdtype == "lt"){
            zmwstokeep <- meanfractionfiltered>=zmwthreshold
          }
          includezmws[includezmws==TRUE] <- zmwstokeep
          includemutations.fwd <- includemutations.fwd[zmwstokeep]
          includemutations.rev <- includemutations.rev[zmwstokeep]
          includemutations.zmw <- includemutations.zmw[zmwstokeep]
          cat("DONE\n")          
          
          cat(paste("        Filtered:",length(which(!zmwstokeep)),"/",totalnumzmws," additional ZMWs =",round(length(which(!zmwstokeep))/totalnumzmws*100,1),"% of ZMWs\n"))
          cat(paste("        Remaining:",length(which(includezmws)),"ZMWs =",round(length(which(includezmws))/totalnumzmws*100,1),"% of ZMWs\n"))
          
          rm(bigwig.granges,fwdccs.granges,revccs.granges,fwdccs.bw.overlaps,revccs.bw.overlaps,fwdccs.bw.intersect,revccs.bw.intersect)
          
          if(!any(includezmws)){
		      	bamfilezmw_filtered[[i]][[j]][[k]] <- list()
		      	cat("\n     No ZMWs remain. Skipping remaining ZMW Bigwig filters.\n\n")
		      	break
          }

        }
        cat("  ZMW Bigwig filters...DONE\n")
      }
      
      if(!any(includezmws)){
      	bamfilezmw_filtered[[i]][[j]][[k]] <- list()
      	cat("\n  No ZMWs remain. Proceeding to next sample/filterset!\n\n")
      	next
      }

      ##Make granges of SNV mutation positions.
      cat(" Preparing GRanges of SNV mutations...")
      fwdchroms <- mapply(function(x,y){rep(x,length(y))},as.character(bamfilezmw_all[[i]][[j]]$rname.fwd[includezmws]),bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws],USE.NAMES=F,SIMPLIFY=FALSE)
      revchroms <- mapply(function(x,y){rep(x,length(y))},as.character(bamfilezmw_all[[i]][[j]]$rname.rev[includezmws]),bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws],USE.NAMES=F,SIMPLIFY=FALSE)
      zmwchroms <- mapply(function(x,y){rep(x,length(y))},as.character(bamfilezmw_all[[i]][[j]]$rname.fwd[includezmws]),bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws],USE.NAMES=F,SIMPLIFY=FALSE)
      
      fwdindex <- mapply(function(x,y){rep(x,y)},seq_along(bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws]),lapply(bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws],length),SIMPLIFY=FALSE)
      revindex <- mapply(function(x,y){rep(x,y)},seq_along(bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws]),lapply(bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws],length),SIMPLIFY=FALSE)
      zmwindex <- mapply(function(x,y){rep(x,y)},seq_along(bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws]),lapply(bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws],length),SIMPLIFY=FALSE)
      
      fwdvariantindex <- lapply(bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws],function(x){if(length(x)>0){seq(length(x))}else{return(c())}})
      revvariantindex <- lapply(bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws],function(x){if(length(x)>0){seq(length(x))}else{return(c())}})
      zmwvariantindex <- lapply(bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws],function(x){if(length(x)>0){seq(length(x))}else{return(c())}})
      
      fwdholes <- mapply(function(x,y){rep(x,length(y))},bamfilezmw_all[[i]][[j]]$zmwname.fwd[includezmws],bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws],USE.NAMES=F,SIMPLIFY=FALSE)
      revholes <- mapply(function(x,y){rep(x,length(y))},bamfilezmw_all[[i]][[j]]$zmwname.rev[includezmws],bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws],USE.NAMES=F,SIMPLIFY=FALSE)
      zmwholes <- mapply(function(x,y){rep(x,length(y))},bamfilezmw_all[[i]][[j]]$zmwname.fwd[includezmws],bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws],USE.NAMES=F,SIMPLIFY=FALSE)
      
      fwdmutations <- data.frame(chrom=unlist(fwdchroms),start=unlist(bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws]),end=unlist(bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws]),zmwindex=unlist(fwdindex),variantindex=unlist(fwdvariantindex),qn=unlist(bamfilezmw_all[[i]][[j]]$tag.fwd$qn[includezmws]),zmw=unlist(fwdholes))
      revmutations <- data.frame(chrom=unlist(revchroms),start=unlist(bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws]),end=unlist(bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws]),zmwindex=unlist(revindex),variantindex=unlist(revvariantindex),qn=unlist(bamfilezmw_all[[i]][[j]]$tag.rev$qn[includezmws]),zmw=unlist(revholes))
      zmwmutations <- data.frame(chrom=unlist(zmwchroms),start=unlist(bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws]),end=unlist(bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws]),zmwindex=unlist(zmwindex),variantindex=unlist(zmwvariantindex),qn=unlist(bamfilezmw_all[[i]][[j]]$zmw.qn[includezmws]),zmw=unlist(zmwholes))
      
      #Convert full ZMW names to hole numbers for later subreads analysis
      fwdmutations$zmw <- gsub(".*?/(.*)/.*","\\1",fwdmutations$zmw)
      revmutations$zmw <- gsub(".*?/(.*)/.*","\\1",revmutations$zmw)
      zmwmutations$zmw <- gsub(".*?/(.*)/.*","\\1",zmwmutations$zmw)
  
      fwdmutations.granges <- makeGRangesFromDataFrame(fwdmutations,keep.extra.columns=TRUE,ignore.strand=TRUE)
      revmutations.granges <- makeGRangesFromDataFrame(revmutations,keep.extra.columns=TRUE,ignore.strand=TRUE)
      zmwmutations.granges <- makeGRangesFromDataFrame(zmwmutations,keep.extra.columns=TRUE,ignore.strand=TRUE)
      
      #Remove temporary objects
      rm(fwdchroms,revchroms,zmwchroms,fwdindex,revindex,zmwindex,fwdvariantindex,revvariantindex,zmwvariantindex,fwdholes,revholes,zmwholes,fwdmutations,revmutations,zmwmutations)
      cat("DONE\n")
      
      ##Change includemutations.fwd and includemutations.rev so that they only include ssDNA-only mutations, and not ZMW mutations. This will remove all ZMW mutations identified in the original CCS BAM file, not just the final ZMW mutations remaining after filtering, because we don't want to only keep ssDNA mutations that do not overlap the highest quality ZMW mutations, but rather any ZMW mutation. Up until this point ssDNA mutation lists included both ssDNA and ZMW mutations. The filter is simple: just subtract from Fwd and Rev mutation positions the positions of ZMW mutations. This filter will also retain as ssDNA mutations those positions where Fwd and Rev strand are mutated to different bases, because they will not be a ZMW mutation (because the original ZMW mutation calling required both Fwd and Rev strand to be mutated to the same base). Note: We do NOT need to add the positions of all the original ZMW mutations to the ssDNA region filters, because technically those bases are still being interrogated for ssDNA mutations, because even if there is a ZMW (dsDNA) mutation (including all germline sites), a ssDNA mutation of one of the strands to a different base (or back to a reference base) will still be detected.
      cat(" Removing ssDNA mutations that are dsDNA mutations...")
      includemutations.fwd <- mapply(function(x,y,q){y & !(x %in% intersect(x[y],q)) },bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws],includemutations.fwd,bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws],SIMPLIFY=FALSE)
      includemutations.rev <- mapply(function(x,y,q){y & !(x %in% intersect(x[y],q)) },bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws],includemutations.rev,bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws],SIMPLIFY=FALSE)
      cat("DONE\n")
      
      ##Apply SNV Bigwig filters
      if(any(unlist(lapply(bigwigfilters[[k]],function(x){x$bigwigfiltertype=="SNV"})))){
        cat(" Applying SNV Bigwig filters...\n")
        for(bigwigfileindex in which(unlist(lapply(bigwigfilters[[k]],function(x){x$bigwigfiltertype=="SNV"})))){
          cat("  Bigwig file:",basename(bigwigfilters[[k]][[bigwigfileindex]]$bigwigfile),"...")
          cat("apply to mutation types (ss/ds/both):",bigwigfilters[[k]][[bigwigfileindex]]$applyto,"...")
          
          #Extract threshold
          if(grepl("^gte",bigwigfilters[[k]][[bigwigfileindex]]$threshold)){
            bigwigthresholdtype <- "gte"
            bigwigthreshold <- as.numeric(sub("^gte","",bigwigfilters[[k]][[bigwigfileindex]]$threshold))
          }else if(grepl("^lt",bigwigfilters[[k]][[bigwigfileindex]]$threshold)){
            bigwigthresholdtype <- "lt"
            bigwigthreshold <- as.numeric(sub("^lt","",bigwigfilters[[k]][[bigwigfileindex]]$threshold))
          }else{
            stop(paste("Wrong format of bigwigfilter",bigwigfilters[[k]][[bigwigfileindex]]$threshold))
          }
          
          bigwig.granges <- bw.granges.thresholding(bigwigfilters[[k]][[bigwigfileindex]]$bigwigfile,referenceinfo$chrsizes,bigwigfilters[[k]][[bigwigfileindex]]$binsize,bigwigthreshold)
          bigwig.granges <- bigwig.granges[[bigwigthresholdtype]]
          bigwig.granges <- filter.granges.by.chroms(bigwig.granges,chrs_to_analyze)
          
          if(length(bigwig.granges)==0){
          	cat("Skipping: no ranges remain in bigwigfilter in selected chromosomes.\n")
          	next
          }

          #Apply padding: Expand bigwigthresholdtype (filtered regions), and reduce.
          bigwig.granges <- resize(bigwig.granges,width=width(bigwig.granges)+(bigwigfilters[[k]][[bigwigfileindex]]$padding*2),fix="center")
          bigwig.granges <- reduce(bigwig.granges)
          
          #Find overlaps - apply to ssDNA or dsDNA mutations, depending on 'applyto' settings.
          if(bigwigfilters[[k]][[bigwigfileindex]]$applyto %in% c("ss","both")){
          	
            fwdmutations.bigwig <- countOverlaps(fwdmutations.granges,bigwig.granges,ignore.strand=TRUE)
            revmutations.bigwig <- countOverlaps(revmutations.granges,bigwig.granges,ignore.strand=TRUE)
          }
          if(bigwigfilters[[k]][[bigwigfileindex]]$applyto %in% c("ds","both")){
            zmwmutations.bigwig <- countOverlaps(zmwmutations.granges,bigwig.granges,ignore.strand=TRUE)
          }

          #Filter out mutations that overlap bigwig ranges
          if(bigwigfilters[[k]][[bigwigfileindex]]$applyto %in% c("ss","both")){
            fwdmutations.bigwig.include <- lapply(includemutations.fwd,function(x){x|TRUE})
            if(any(fwdmutations.bigwig>0)){ #Perform filtering if any overlaps found with at least one mutation
            	invisible(apply(cbind(fwdmutations.granges$zmwindex[fwdmutations.bigwig>0],fwdmutations.granges$variantindex[fwdmutations.bigwig>0]),1,function(x){fwdmutations.bigwig.include[[x[1]]][x[2]] <<- FALSE}))
            }
            includemutations.fwd <- mapply(function(x,y){x&y},includemutations.fwd,fwdmutations.bigwig.include,SIMPLIFY=FALSE)
            
            revmutations.bigwig.include <- lapply(includemutations.rev,function(x){x|TRUE})
            if(any(revmutations.bigwig>0)){ #Perform filtering if any overlaps found with at least one mutation
            	invisible(apply(cbind(revmutations.granges$zmwindex[revmutations.bigwig>0],revmutations.granges$variantindex[revmutations.bigwig>0]),1,function(x){revmutations.bigwig.include[[x[1]]][x[2]] <<- FALSE}))
            }
            includemutations.rev <- mapply(function(x,y){x&y},includemutations.rev,revmutations.bigwig.include,SIMPLIFY=FALSE)
            rm(fwdmutations.bigwig,revmutations.bigwig)
          }
          if(bigwigfilters[[k]][[bigwigfileindex]]$applyto %in% c("ds","both")){
            sharedmutations.bigwig.include <- lapply(includemutations.zmw,function(x){x|TRUE})
            if(any(zmwmutations.bigwig>0)){ #Perform filtering if any overlaps found with at least one mutation
            	invisible(apply(cbind(zmwmutations.granges$zmwindex[zmwmutations.bigwig>0],zmwmutations.granges$variantindex[zmwmutations.bigwig>0]),1,function(x){sharedmutations.bigwig.include[[x[1]]][x[2]] <<- FALSE}))
            }
            includemutations.zmw <- mapply(function(x,y){x&y},includemutations.zmw,sharedmutations.bigwig.include,SIMPLIFY=FALSE)
            rm(zmwmutations.bigwig)
          }
          
          cat("DONE\n")
          
          #Save regional filter for later mutation frequency calculation
          if(bigwigfilters[[k]][[bigwigfileindex]]$applyto %in% c("ss","both")){
            readregionalfilters.fwd <- c(readregionalfilters.fwd,bigwig.granges)
            readregionalfilters.rev <- c(readregionalfilters.rev,bigwig.granges)
            genomeregionalfilters.fwd <- c(genomeregionalfilters.fwd,bigwig.granges)
            genomeregionalfilters.rev <- c(genomeregionalfilters.rev,bigwig.granges)
          }
          if(bigwigfilters[[k]][[bigwigfileindex]]$applyto %in% c("ds","both")){
            readregionalfilters.zmw <- c(readregionalfilters.zmw,bigwig.granges)
            genomeregionalfilters.zmw <- c(genomeregionalfilters.zmw,bigwig.granges)
          }
          
          #Output filtering statistics
          if(bigwigfilters[[k]][[bigwigfileindex]]$applyto == "ss"){
            genome.granges.filteredtemp.ssDNA <- GenomicRanges::setdiff(genome.granges.filtered.ssDNA,bigwig.granges)
            sumfiltergranges.ssDNA <- sum(width(bigwig.granges))
            sumadditionalfiltergranges.ssDNA <- sum(width(GenomicRanges::setdiff(genome.granges.filtered.ssDNA,genome.granges.filteredtemp.ssDNA)))
            
            genome.granges.filteredtemp <- genome.granges.filtered
            sumfiltergranges <- 0
            sumadditionalfiltergranges <- 0
          }else if(bigwigfilters[[k]][[bigwigfileindex]]$applyto == "ds"){
            genome.granges.filteredtemp <- GenomicRanges::setdiff(genome.granges.filtered,bigwig.granges)
            sumfiltergranges <- sum(width(bigwig.granges))
            sumadditionalfiltergranges <- sum(width(GenomicRanges::setdiff(genome.granges.filtered,genome.granges.filteredtemp)))
            
            genome.granges.filteredtemp.ssDNA <- genome.granges.filtered.ssDNA
            sumfiltergranges.ssDNA <- 0
            sumadditionalfiltergranges.ssDNA <- 0
          }else if(bigwigfilters[[k]][[bigwigfileindex]]$applyto == "both"){
            genome.granges.filteredtemp.ssDNA <- GenomicRanges::setdiff(genome.granges.filtered.ssDNA,bigwig.granges)
            sumfiltergranges.ssDNA <- sum(width(bigwig.granges))
            sumadditionalfiltergranges.ssDNA <- sum(width(GenomicRanges::setdiff(genome.granges.filtered.ssDNA,genome.granges.filteredtemp.ssDNA)))
            
            genome.granges.filteredtemp <- GenomicRanges::setdiff(genome.granges.filtered,bigwig.granges)
            sumfiltergranges <- sum(width(bigwig.granges))
            sumadditionalfiltergranges <- sum(width(GenomicRanges::setdiff(genome.granges.filtered,genome.granges.filteredtemp)))
          }

          cat(paste("        Filter covers:",sumfiltergranges.ssDNA,"/",sum(width(genome.granges)),"total genome bases  =",round(sumfiltergranges.ssDNA/sum(width(genome.granges))*100,1),"% of total genome for ssDNA mutations\n"))
          cat(paste("        Filter removes:",sumadditionalfiltergranges.ssDNA,"/",sum(width(genome.granges)),"total genome bases beyond prior filters =",round(sumadditionalfiltergranges.ssDNA/sum(width(genome.granges))*100,1),"% of total genome for ssDNA mutations\n"))
          
          cat(paste("        Filter covers:",sumfiltergranges,"/",sum(width(genome.granges)),"total genome bases  =",round(sumfiltergranges/sum(width(genome.granges))*100,1),"% of total genome for dsDNA mutations\n"))
          cat(paste("        Filter removes:",sumadditionalfiltergranges,"/",sum(width(genome.granges)),"total genome bases beyond prior filters =",round(sumadditionalfiltergranges/sum(width(genome.granges))*100,1),"% of total genome for dsDNA mutations\n"))
          
          genome.granges.filtered.ssDNA <- genome.granges.filteredtemp.ssDNA
          genome.granges.filtered <- genome.granges.filteredtemp
          
          sumremaininggranges <- sum(width(genome.granges.filtered.ssDNA))
          cat(paste("        Remaining:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for ssDNA mutations\n"))
          sumremaininggranges <- sum(width(genome.granges.filtered))
          cat(paste("        Remaining:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for dsDNA mutations\n"))
          
          rm(bigwig.granges)

        }
        cat("  SNV Bigwig filters...DONE\n")
      }

      ##Apply base quality filters (set to TRUE for any variant that passes filter). For dsDNA mutations, bases on both fwd and rev strand reads must pass filter.
      cat(" Applying base quality filters...")
      
      fwd.filter <- lapply(bamfilezmw_all[[i]][[j]]$tag.fwd$qq[includezmws],function(x){x < thresholds[k,"minqq"]})
      rev.filter <- lapply(bamfilezmw_all[[i]][[j]]$tag.rev$qq[includezmws],function(x){x < thresholds[k,"minqq"]})
      zmw.filter <- mapply(function(a,b){a < thresholds[k,"minqq"] | b < thresholds[k,"minqq"] },bamfilezmw_all[[i]][[j]]$zmw.fwdqq[includezmws],bamfilezmw_all[[i]][[j]]$zmw.revqq[includezmws],SIMPLIFY=FALSE)
      
      includemutations.fwd <- mapply(function(x,y){x & !y},includemutations.fwd,fwd.filter,SIMPLIFY=FALSE)
      includemutations.rev <- mapply(function(x,y){x & !y},includemutations.rev,rev.filter,SIMPLIFY=FALSE)
      includemutations.zmw <- mapply(function(x,y){x & !y},includemutations.zmw,zmw.filter,SIMPLIFY=FALSE)
      rm(fwd.filter,rev.filter,zmw.filter)
      
      #Save regional filters for later mutation frequency calculation.
      fwd.qqposlist <- sapply(as.character(sequenceLayer(bamfilezmw_all[[i]][[j]]$qual.fwd[includezmws],bamfilezmw_all[[i]][[j]]$cigar.fwd[includezmws],D.letter="!",N.letter="!")),function(x){which( (as.integer(charToRaw(x))-33) < thresholds[k,"minqq"])},USE.NAMES=FALSE,simplify=FALSE)
			fwd.qqposlist <- mapply(function(x,y){x+y-1},x=fwd.qqposlist,y=bamfilezmw_all[[i]][[j]]$pos.fwd[includezmws],SIMPLIFY=FALSE)
			fwd.qqposlist.notisempty <- lapply(fwd.qqposlist,length)>0
			fwd.qqposlist <- fwd.qqposlist[fwd.qqposlist.notisempty]
      
      rev.qqposlist <- sapply(as.character(sequenceLayer(bamfilezmw_all[[i]][[j]]$qual.rev[includezmws],bamfilezmw_all[[i]][[j]]$cigar.rev[includezmws],D.letter="!",N.letter="!")),function(x){which( (as.integer(charToRaw(x))-33) < thresholds[k,"minqq"])},USE.NAMES=FALSE,simplify=FALSE)
			rev.qqposlist <- mapply(function(x,y){x+y-1},x=rev.qqposlist,y=bamfilezmw_all[[i]][[j]]$pos.rev[includezmws],SIMPLIFY=FALSE)
			rev.qqposlist.notisempty <- lapply(rev.qqposlist,length)>0
			rev.qqposlist <- rev.qqposlist[rev.qqposlist.notisempty]
      
			if(length(fwd.qqposlist)>0){
	      names(fwd.qqposlist) <- bamfilezmw_all[[i]][[j]]$rname.fwd[includezmws][fwd.qqposlist.notisempty]
	      fwd.qqfilter.chrom <- as.character(stack(fwd.qqposlist)$ind)
	      names(fwd.qqposlist) <- bamfilezmw_all[[i]][[j]]$qname.fwd[includezmws][fwd.qqposlist.notisempty]
	      fwd.qqfilter.qname <- as.character(stack(fwd.qqposlist)$ind)
	      
	      fwd.qqposlist <- unlist(fwd.qqposlist,use.names=FALSE)
	      
	      fwd.qqfilter.granges <- makeGRangesFromDataFrame(data.frame(seqnames=fwd.qqfilter.chrom,start=fwd.qqposlist,end=fwd.qqposlist,appliestoread=fwd.qqfilter.qname),keep.extra.columns=TRUE)
	      
	      rm(fwd.qqfilter.chrom,fwd.qqfilter.qname)
			}else{
				#No mutations filtered - make empty granges
				fwd.qqfilter.granges <- makeGRangesFromDataFrame(data.frame(seqnames=1,start=1,end=1,appliestoread=1),keep.extra.columns=TRUE)
				fwd.qqfilter.granges[1] <- NULL
			}
      
			if(length(rev.qqposlist)>0){
	      names(rev.qqposlist) <- bamfilezmw_all[[i]][[j]]$rname.rev[includezmws][rev.qqposlist.notisempty]
	      rev.qqfilter.chrom <- as.character(stack(rev.qqposlist)$ind)
	      names(rev.qqposlist) <- bamfilezmw_all[[i]][[j]]$qname.rev[includezmws][rev.qqposlist.notisempty]
	      rev.qqfilter.qname <- as.character(stack(rev.qqposlist)$ind)
	      
	      rev.qqposlist <- unlist(rev.qqposlist,use.names=FALSE)
	      
	      rev.qqfilter.granges <- makeGRangesFromDataFrame(data.frame(seqnames=rev.qqfilter.chrom,start=rev.qqposlist,end=rev.qqposlist,appliestoread=rev.qqfilter.qname),keep.extra.columns=TRUE)
	      
	      rm(rev.qqfilter.chrom,rev.qqfilter.qname)
			}else{
				#No mutations filtered - make empty granges
				rev.qqfilter.granges <- makeGRangesFromDataFrame(data.frame(seqnames=1,start=1,end=1,appliestoread=1),keep.extra.columns=TRUE)
				rev.qqfilter.granges[1] <- NULL
			}
			
			rm(fwd.qqposlist,rev.qqposlist,fwd.qqposlist.notisempty,rev.qqposlist.notisempty)
      
      #qqfilter is only stored in readregionalfilters, because it is read-specific and should not be used as a genome region filter in mutation frequency calculations.
      readregionalfilters.fwd <- c(readregionalfilters.fwd,fwd.qqfilter.granges)
      readregionalfilters.rev <- c(readregionalfilters.rev,rev.qqfilter.granges)
      readregionalfilters.zmw <- c(readregionalfilters.zmw,fwd.qqfilter.granges,rev.qqfilter.granges)
      
      cat("DONE\n")      
      #Output filtering statistics
      cat(paste("        Filter covers:",sum(width(c(fwd.qqfilter.granges,rev.qqfilter.granges))),"fwd and rev strand bases below threshold out of",sum(width(zmw.granges))*2,"total (unfiltered) ZMW fwd and rev bases (ZMW lengths * 2)\n"))
      
      #Delete temporary objects
      rm(fwd.qqfilter.granges,rev.qqfilter.granges)

      
      ##Apply read location based filters
      cat(" Applying read location-based filters...")
      
      fwd.filter <- mapply(function(x,a,b){(x-a+1)<=thresholds[k,"bpends"] | (a+b-x)<=thresholds[k,"bpends"] },bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws],bamfilezmw_all[[i]][[j]]$pos.fwd[includezmws],bamfilezmw_all[[i]][[j]]$isize.fwd[includezmws],SIMPLIFY=FALSE)
      rev.filter <- mapply(function(x,a,b){(x-a+1)<=thresholds[k,"bpends"] | (a+b-x)<=thresholds[k,"bpends"] },bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws],bamfilezmw_all[[i]][[j]]$pos.rev[includezmws],bamfilezmw_all[[i]][[j]]$isize.rev[includezmws],SIMPLIFY=FALSE)
      zmw.filter <- mapply(function(x,a,b,c,d){(x-a+1)<=thresholds[k,"bpends"] | (x-c+1)<=thresholds[k,"bpends"] | (a+b-x)<=thresholds[k,"bpends"] | (c+d-x)<=thresholds[k,"bpends"]},bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws],bamfilezmw_all[[i]][[j]]$pos.fwd[includezmws],bamfilezmw_all[[i]][[j]]$isize.fwd[includezmws],bamfilezmw_all[[i]][[j]]$pos.rev[includezmws],bamfilezmw_all[[i]][[j]]$isize.rev[includezmws],SIMPLIFY=FALSE)
      
      includemutations.fwd <- mapply(function(x,y){x&!y},includemutations.fwd,fwd.filter,SIMPLIFY=FALSE)
      includemutations.rev <- mapply(function(x,y){x&!y},includemutations.rev,rev.filter,SIMPLIFY=FALSE)
      includemutations.zmw <- mapply(function(x,y){x&!y},includemutations.zmw,zmw.filter,SIMPLIFY=FALSE)
      
      #Save regional filters for later mutation frequency calculation
      fwd.bpends.left.start <- bamfilezmw_all[[i]][[j]]$pos.fwd[includezmws]
      fwd.bpends.left.end <- bamfilezmw_all[[i]][[j]]$pos.fwd[includezmws]+thresholds[k,"bpends"]-1
      fwd.bpends.right.start <- fwd.bpends.left.start+bamfilezmw_all[[i]][[j]]$isize.fwd[includezmws]-thresholds[k,"bpends"]
      fwd.bpends.right.end <- fwd.bpends.right.start+thresholds[k,"bpends"]-1
      
      rev.bpends.left.start <- bamfilezmw_all[[i]][[j]]$pos.rev[includezmws]
      rev.bpends.left.end <- bamfilezmw_all[[i]][[j]]$pos.rev[includezmws]+thresholds[k,"bpends"]-1
      rev.bpends.right.start <- rev.bpends.left.start+bamfilezmw_all[[i]][[j]]$isize.rev[includezmws]-thresholds[k,"bpends"]
      rev.bpends.right.end <- rev.bpends.right.start+thresholds[k,"bpends"]-1
      
      fwd.bpends.filter <- c(makeGRangesFromDataFrame(data.frame(seqnames=bamfilezmw_all[[i]][[j]]$rname.fwd[includezmws],start=fwd.bpends.left.start,end=fwd.bpends.left.end,appliestoread=bamfilezmw_all[[i]][[j]]$qname.fwd[includezmws]),keep.extra.columns=TRUE),makeGRangesFromDataFrame(data.frame(seqnames=bamfilezmw_all[[i]][[j]]$rname.fwd[includezmws],start=fwd.bpends.right.start,end=fwd.bpends.right.end,appliestoread=bamfilezmw_all[[i]][[j]]$qname.fwd[includezmws]),keep.extra.columns=TRUE))
      rev.bpends.filter <- c(makeGRangesFromDataFrame(data.frame(seqnames=bamfilezmw_all[[i]][[j]]$rname.rev[includezmws],start=rev.bpends.left.start,end=rev.bpends.left.end,appliestoread=bamfilezmw_all[[i]][[j]]$qname.rev[includezmws]),keep.extra.columns=TRUE),makeGRangesFromDataFrame(data.frame(seqnames=bamfilezmw_all[[i]][[j]]$rname.rev[includezmws],start=rev.bpends.right.start,end=rev.bpends.right.end,appliestoread=bamfilezmw_all[[i]][[j]]$qname.rev[includezmws]),keep.extra.columns=TRUE))
      
      #bpends is only stored in readregionalfilters, because it is read-specific and should not be used as a genome region filter in mutation frequency calculations.
      readregionalfilters.fwd <- c(readregionalfilters.fwd,fwd.bpends.filter)
      readregionalfilters.rev <- c(readregionalfilters.rev,rev.bpends.filter)
      readregionalfilters.zmw <- c(readregionalfilters.zmw,fwd.bpends.filter,rev.bpends.filter)

      cat("DONE\n")      
      
      #Output filtering statistics
      cat(paste("        Filter covers:",sum(width(GenomicRanges::union(fwd.bpends.filter,rev.bpends.filter))),"bases out of",sum(width(zmw.granges)),"total (unfiltered) ZMW bases\n"))
      
      #Delete temporary objects
      rm(fwd.bpends.filter,rev.bpends.filter,fwd.bpends.left.start,fwd.bpends.left.end,rev.bpends.left.start,rev.bpends.left.end)

      ##Apply VCF indel filters
      cat(" Applying VCF indel filters...\n")
      for(vcffileindex in 1:length(vcffilters[[k]][[i]])){
          if(vcffilters[[k]][[i]][[vcffileindex]]$vcfINDELfilter){
            vcffile <- vcffilters[[k]][[i]][[vcffileindex]]$vcffilename
            cat("  VCF file:",basename(vcffile),"...\n")
            
            filterstoinclude <- paste0(paste0("^",vcffilters[[k]][[i]][[vcffileindex]]$INDELFILTERS,"$"),collapse="|")
            
            #Get all indels that pass vcf indel filters
            includevcfindels <- bamfilezmw_all[[i]][[j]]$vcfIndels[[vcffile]]$Depth>=vcffilters[[k]][[i]][[vcffileindex]]$vcfINDELDepth & bamfilezmw_all[[i]][[j]]$vcfIndels[[vcffile]]$GQ>=vcffilters[[k]][[i]][[vcffileindex]]$vcfINDELGQ & bamfilezmw_all[[i]][[j]]$vcfIndels[[vcffile]]$VAF>=vcffilters[[k]][[i]][[vcffileindex]]$vcfINDELVAF & bamfilezmw_all[[i]][[j]]$vcfIndels[[vcffile]]$QUAL>=vcffilters[[k]][[i]][[vcffileindex]]$vcfINDELQUAL & grepl(filterstoinclude,bamfilezmw_all[[i]][[j]]$vcfIndels[[vcffile]]$FILTER)
            
            ##Adjust vcf indels ranges with specified padding
            vcfinsertions <- bamfilezmw_all[[i]][[j]]$vcfIndels[[vcffile]][includevcfindels & ((nchar(bamfilezmw_all[[i]][[j]]$vcfIndels[[vcffile]]$Ref)-nchar(bamfilezmw_all[[i]][[j]]$vcfIndels[[vcffile]]$Alt))<0)]
            vcfdeletions <- bamfilezmw_all[[i]][[j]]$vcfIndels[[vcffile]][includevcfindels & ((nchar(bamfilezmw_all[[i]][[j]]$vcfIndels[[vcffile]]$Ref)-nchar(bamfilezmw_all[[i]][[j]]$vcfIndels[[vcffile]]$Alt))>0)]
            
            #Keep only vcf indels in analyzed chromosomes
            vcfinsertions <- filter.granges.by.chroms(vcfinsertions,chrs_to_analyze)
            vcfdeletions <- filter.granges.by.chroms(vcfdeletions,chrs_to_analyze)
            
          if(length(c(vcfinsertions,vcfdeletions))==0){
          	cat("        Skipping: no indels in VCF in selected chromosomes.\n")
          	next
          }
            
            #Insertion padding
            vcfINDELinspad <- as.numeric(unlist(str_split(vcffilters[[k]][[i]][[vcffileindex]]$vcfINDELinspad,"m|b"))[2:3])
            multpadding <- vcfINDELinspad[1]
            bppadding <- vcfINDELinspad[2]
            padding <- apply(cbind(round(multpadding*(nchar(vcfinsertions$Alt)-nchar(vcfinsertions$Ref))),rep(bppadding,length(vcfinsertions))),1,max)
            start(vcfinsertions) <- start(vcfinsertions)-padding
            end(vcfinsertions) <- end(vcfinsertions)+padding
            
            #Deletion padding
            vcfINDELdelpad <- as.numeric(unlist(str_split(vcffilters[[k]][[i]][[vcffileindex]]$vcfINDELdelpad,"m|b"))[2:3])
            multpadding <- vcfINDELdelpad[1]
            bppadding <- vcfINDELdelpad[2]
            padding <- apply(cbind(round(multpadding*width(vcfdeletions)),rep(bppadding,length(vcfdeletions))),1,max)
            start(vcfdeletions) <- start(vcfdeletions) - padding
            end(vcfdeletions) <- end(vcfdeletions) + padding
            
            #Find overlaps
            fwdmutations.vcfindels <- countOverlaps(fwdmutations.granges,c(vcfinsertions,vcfdeletions),ignore.strand=TRUE)
            revmutations.vcfindels <- countOverlaps(revmutations.granges,c(vcfinsertions,vcfdeletions),ignore.strand=TRUE)
            zmwmutations.vcfindels <- countOverlaps(zmwmutations.granges,c(vcfinsertions,vcfdeletions),ignore.strand=TRUE)
            
            #Filter out mutations that overlap vcf indel ranges
            fwdmutations.vcfindels.include <- lapply(includemutations.fwd,function(x){x|TRUE})
            if(any(fwdmutations.vcfindels>0)){ #Perform filtering if any overlaps found with at least one mutation
            	invisible(apply(cbind(fwdmutations.granges$zmwindex[fwdmutations.vcfindels>0],fwdmutations.granges$variantindex[fwdmutations.vcfindels>0]),1,function(x){fwdmutations.vcfindels.include[[x[1]]][x[2]] <<- FALSE}))
            }
            includemutations.fwd <- mapply(function(x,y){x&y},includemutations.fwd,fwdmutations.vcfindels.include,SIMPLIFY=FALSE)
            
            revmutations.vcfindels.include <- lapply(includemutations.rev,function(x){x|TRUE})
            if(any(revmutations.vcfindels>0)){ #Perform filtering if any overlaps found with at least one mutation
            	invisible(apply(cbind(revmutations.granges$zmwindex[revmutations.vcfindels>0],revmutations.granges$variantindex[revmutations.vcfindels>0]),1,function(x){revmutations.vcfindels.include[[x[1]]][x[2]] <<- FALSE}))
            }
            includemutations.rev <- mapply(function(x,y){x&y},includemutations.rev,revmutations.vcfindels.include,SIMPLIFY=FALSE)
            
            sharedmutations.vcfindels.include <- lapply(includemutations.zmw,function(x){x|TRUE})
            if(any(zmwmutations.vcfindels>0)){ #Perform filtering if any overlaps found with at least one mutation
            	invisible(apply(cbind(zmwmutations.granges$zmwindex[zmwmutations.vcfindels>0],zmwmutations.granges$variantindex[zmwmutations.vcfindels>0]),1,function(x){sharedmutations.vcfindels.include[[x[1]]][x[2]] <<- FALSE}))
            }
            includemutations.zmw <- mapply(function(x,y){x&y},includemutations.zmw,sharedmutations.vcfindels.include,SIMPLIFY=FALSE)
            
            #Save regional filters for later mutation frequency calculation
            readregionalfilters.fwd <- c(readregionalfilters.fwd,vcfinsertions,vcfdeletions)
            readregionalfilters.rev <- c(readregionalfilters.rev,vcfinsertions,vcfdeletions)
            readregionalfilters.zmw <- c(readregionalfilters.zmw,vcfinsertions,vcfdeletions)
            genomeregionalfilters.fwd <- c(genomeregionalfilters.fwd,vcfinsertions,vcfdeletions)
            genomeregionalfilters.rev <- c(genomeregionalfilters.rev,vcfinsertions,vcfdeletions)
            genomeregionalfilters.zmw <- c(genomeregionalfilters.zmw,vcfinsertions,vcfdeletions)
            
            #Output filtering statistics
            sumfiltergranges <- sum(width(GenomicRanges::union(vcfinsertions,vcfdeletions)))
            
            genome.granges.filteredtemp <- GenomicRanges::setdiff(genome.granges.filtered,vcfinsertions)
            genome.granges.filteredtemp <- GenomicRanges::setdiff(genome.granges.filteredtemp,vcfdeletions)
            sumadditionalfiltergranges <- sum(width(GenomicRanges::setdiff(genome.granges.filtered,genome.granges.filteredtemp)))
            
            genome.granges.filteredtemp.ssDNA <- GenomicRanges::setdiff(genome.granges.filtered.ssDNA,vcfinsertions)
            genome.granges.filteredtemp.ssDNA <- GenomicRanges::setdiff(genome.granges.filteredtemp.ssDNA,vcfdeletions)
            sumadditionalfiltergranges.ssDNA <- sum(width(GenomicRanges::setdiff(genome.granges.filtered.ssDNA,genome.granges.filteredtemp.ssDNA)))

            cat(paste("        Filter covers:",sumfiltergranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumfiltergranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions)\n"))
            cat(paste("        Filter removes:",sumadditionalfiltergranges.ssDNA,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) beyond prior filters =",round(sumadditionalfiltergranges.ssDNA/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for ssDNA mutations\n"))           
            cat(paste("        Filter removes:",sumadditionalfiltergranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) beyond prior filters =",round(sumadditionalfiltergranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for dsDNA mutations\n"))
            
            genome.granges.filtered <- genome.granges.filteredtemp
            genome.granges.filtered.ssDNA <- genome.granges.filteredtemp.ssDNA
            
            sumremaininggranges <- sum(width(genome.granges.filtered.ssDNA))
            cat(paste("        Remaining:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions)  for ssDNA mutations\n"))
            sumremaininggranges <- sum(width(genome.granges.filtered))
            cat(paste("        Remaining:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions)  for dsDNA mutations\n"))
            
            #Delete temporary objects
            rm(vcfinsertions,vcfdeletions,includevcfindels,vcfINDELinspad,vcfINDELdelpad,padding,fwdmutations.vcfindels,revmutations.vcfindels,zmwmutations.vcfindels,fwdmutations.vcfindels.include,revmutations.vcfindels.include,sharedmutations.vcfindels.include)
            
          }else{
            cat("NOT USED\n")
          }
      }
      cat("  VCF indel filtering...DONE\n")
        
      ##Apply CCS indel filters
      if(thresholds[k,"ccsindelfilter"]){
        cat(" Applying CCS indel filters...\n")

        ##Make Granges of ccs indels, with padding specified in configuration
        #Format ccs indels
        fwdmutations.ccsindels <- stack(setNames(bamfilezmw_all[[i]][[j]]$ccsIndels.fwd[includezmws], seq_along(bamfilezmw_all[[i]][[j]]$ccsIndels.fwd[includezmws])))
        fwdmutations.ccsindels.chromannot <- stack(setNames(bamfilezmw_all[[i]][[j]]$ccsIndels.fwd[includezmws], bamfilezmw_all[[i]][[j]]$rname.fwd[includezmws]))
        mcols(fwdmutations.ccsindels)$chromannot <- mcols(fwdmutations.ccsindels.chromannot)$name
        fwdmutations.ccsindels.qname <- stack(setNames(bamfilezmw_all[[i]][[j]]$ccsIndels.fwd[includezmws], bamfilezmw_all[[i]][[j]]$qname.fwd[includezmws]))
        mcols(fwdmutations.ccsindels)$appliestoread <- mcols(fwdmutations.ccsindels.qname)$name
        
        #Delete temporary objects
        rm(fwdmutations.ccsindels.chromannot,fwdmutations.ccsindels.qname)
        
        mcols(fwdmutations.ccsindels)$insertionlengths <- stack(setNames(bamfilezmw_all[[i]][[j]]$ccsIndels.insertionlengths.fwd[includezmws], seq_along(bamfilezmw_all[[i]][[j]]$ccsIndels.insertionlengths.fwd[includezmws])))$value
        fwdccsinsertions <- fwdmutations.ccsindels[names(fwdmutations.ccsindels)=="I"]
        fwdccsdeletions <- fwdmutations.ccsindels[names(fwdmutations.ccsindels)=="D"]
        
        #Delete temporary objects
        rm(fwdmutations.ccsindels)

        revmutations.ccsindels <- stack(setNames(bamfilezmw_all[[i]][[j]]$ccsIndels.rev[includezmws], seq_along(bamfilezmw_all[[i]][[j]]$ccsIndels.rev[includezmws])))
        revmutations.ccsindels.chromannot <- stack(setNames(bamfilezmw_all[[i]][[j]]$ccsIndels.rev[includezmws], bamfilezmw_all[[i]][[j]]$rname.rev[includezmws]))
        mcols(revmutations.ccsindels)$chromannot <- mcols(revmutations.ccsindels.chromannot)$name
        revmutations.ccsindels.qname <- stack(setNames(bamfilezmw_all[[i]][[j]]$ccsIndels.rev[includezmws], bamfilezmw_all[[i]][[j]]$qname.rev[includezmws]))
        mcols(revmutations.ccsindels)$appliestoread <- mcols(revmutations.ccsindels.qname)$name
        
        #Delete temporary objects
        rm(revmutations.ccsindels.chromannot,revmutations.ccsindels.qname)
        
        mcols(revmutations.ccsindels)$insertionlengths <- stack(setNames(bamfilezmw_all[[i]][[j]]$ccsIndels.insertionlengths.rev[includezmws], seq_along(bamfilezmw_all[[i]][[j]]$ccsIndels.insertionlengths.rev[includezmws])))$value
        revccsinsertions <- revmutations.ccsindels[names(revmutations.ccsindels)=="I"]
        revccsdeletions <- revmutations.ccsindels[names(revmutations.ccsindels)=="D"]
        
        #Delete temporary objects
        rm(revmutations.ccsindels)

        #Pad ccs indel ranges with the specified padding for insertions and deletions
        ccsindelinspad <- as.numeric(unlist(str_split(thresholds[k,"ccsindelinspad"],"m|b"))[2:3])
        multpadding <- ccsindelinspad[1]
        bppadding <- ccsindelinspad[2]
        fwdpadding <- apply(cbind(round(multpadding*mcols(fwdccsinsertions)$insertionlengths),rep(bppadding,length(fwdccsinsertions))),1,max)
        revpadding <- apply(cbind(round(multpadding*mcols(revccsinsertions)$insertionlengths),rep(bppadding,length(revccsinsertions))),1,max)
        start(fwdccsinsertions) <- start(fwdccsinsertions) - fwdpadding
        start(revccsinsertions) <- start(revccsinsertions) - revpadding
        end(fwdccsinsertions) <- end(fwdccsinsertions) + fwdpadding
        end(revccsinsertions) <- end(revccsinsertions) + revpadding

        ccsindeldelpad <- as.numeric(unlist(str_split(thresholds[k,"ccsindeldelpad"],"m|b"))[2:3])
        multpadding <- ccsindeldelpad[1]
        bppadding <- ccsindeldelpad[2]
        fwdpadding <- apply(cbind(round(multpadding*width(fwdccsdeletions)),rep(bppadding,length(fwdccsdeletions))),1,max)
        revpadding <- apply(cbind(round(multpadding*width(revccsdeletions)),rep(bppadding,length(revccsdeletions))),1,max)
        start(fwdccsdeletions) <- start(fwdccsdeletions) - fwdpadding
        start(revccsdeletions) <- start(revccsdeletions) - revpadding
        end(fwdccsdeletions) <- end(fwdccsdeletions) + fwdpadding
        end(revccsdeletions) <- end(revccsdeletions) + revpadding
        
        #Delete temporary objects
        rm(fwdpadding,revpadding)

        #Make fwd strand ccs indels Granges
        fwdmutations.ccsindels <- c(fwdccsinsertions,fwdccsdeletions)
        rm(fwdccsinsertions,fwdccsdeletions)
        fwdmutations.ccsindels <- data.frame(fwdmutations.ccsindels)
        fwdmutations.ccsindels$names <- NULL
        fwdmutations.ccsindels$width <- NULL
        fwdmutations.ccsindels$insertionlengths <- NULL
        if(nrow(fwdmutations.ccsindels)==0){
          fwdmutations.ccsindels <- GRanges()
          fwdmutations.ccsindels$name <- Rle()
          fwdmutations.ccsindels$chromannot <- Rle()
          fwdmutations.ccsindels$appliestoread <- Rle()
        }else{
          fwdmutations.ccsindels <- makeGRangesFromDataFrame(fwdmutations.ccsindels,seqnames.field="name",keep.extra.columns=TRUE,ignore.strand=TRUE)
        }

        #Make rev strand ccs indels Granges
        revmutations.ccsindels <- c(revccsinsertions,revccsdeletions)
        rm(revccsinsertions,revccsdeletions)
        revmutations.ccsindels <- data.frame(revmutations.ccsindels)
        revmutations.ccsindels$names <- NULL
        revmutations.ccsindels$width <- NULL
        revmutations.ccsindels$insertionlengths <- NULL
        if(nrow(revmutations.ccsindels)==0){
          revmutations.ccsindels <- GRanges()
          revmutations.ccsindels$name <- Rle()
          revmutations.ccsindels$chromannot <- Rle()
          revmutations.ccsindels$appliestoread <- Rle()
        }else{
          revmutations.ccsindels <- makeGRangesFromDataFrame(revmutations.ccsindels,seqnames.field="name",keep.extra.columns=TRUE,ignore.strand=TRUE)
        }

        ##Make Granges objects of mutations, indexed by zmw index.
        fwdmutations.granges.zmwindex <- makeGRangesFromDataFrame(data.frame(seqnames=fwdmutations.granges$zmwindex,start=start(fwdmutations.granges),end=end(fwdmutations.granges)))
        revmutations.granges.zmwindex <- makeGRangesFromDataFrame(data.frame(seqnames=revmutations.granges$zmwindex,start=start(revmutations.granges),end=end(revmutations.granges)))
        zmwmutations.granges.zmwindex <- makeGRangesFromDataFrame(data.frame(seqnames=zmwmutations.granges$zmwindex,start=start(zmwmutations.granges),end=end(zmwmutations.granges)))
        
        ##Find overlaps of SNV mutations with the padded ccs indel ranges
          #fwd strand mutations
          fwdmutations.granges.zmwindex <- data.frame(fwdmutations.granges.zmwindex[countOverlaps(fwdmutations.granges.zmwindex,fwdmutations.ccsindels,ignore.strand=TRUE)>0])
          fwdmutations.granges.zmwindex$seqnames <- as.numeric(as.character(fwdmutations.granges.zmwindex$seqnames))

          include.fwdmutations.granges.zmwindex <- bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws]
          for(zmwindex in unique(fwdmutations.granges.zmwindex$seqnames)){
            include.fwdmutations.granges.zmwindex[[zmwindex]] <- !grepl(paste0("^",fwdmutations.granges.zmwindex[fwdmutations.granges.zmwindex$seqnames==zmwindex,"start"],"$",collapse="|"),include.fwdmutations.granges.zmwindex[[zmwindex]])
          }
          include.fwdmutations.granges.zmwindex <- mapply(function(x,y){x&y},lapply(bamfilezmw_all[[i]][[j]]$tag.fwd$rp[includezmws],function(x){x|TRUE}),include.fwdmutations.granges.zmwindex,SIMPLIFY=FALSE)

          #rev strand mutations
          revmutations.granges.zmwindex <- data.frame(revmutations.granges.zmwindex[countOverlaps(revmutations.granges.zmwindex,revmutations.ccsindels,ignore.strand=TRUE)>0])
          revmutations.granges.zmwindex$seqnames <- as.numeric(as.character(revmutations.granges.zmwindex$seqnames))

          include.revmutations.granges.zmwindex <- bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws]
          for(zmwindex in unique(revmutations.granges.zmwindex$seqnames)){
            include.revmutations.granges.zmwindex[[zmwindex]] <- !grepl(paste0("^",revmutations.granges.zmwindex[revmutations.granges.zmwindex$seqnames==zmwindex,"start"],"$",collapse="|"),include.revmutations.granges.zmwindex[[zmwindex]])
          }
          include.revmutations.granges.zmwindex <- mapply(function(x,y){x&y},lapply(bamfilezmw_all[[i]][[j]]$tag.rev$rp[includezmws],function(x){x|TRUE}),include.revmutations.granges.zmwindex,SIMPLIFY=FALSE)

          #zmw (dsDNA) mutations
          zmwmutations.granges.zmwindex <- data.frame(zmwmutations.granges.zmwindex[countOverlaps(zmwmutations.granges.zmwindex,c(fwdmutations.ccsindels,revmutations.ccsindels),ignore.strand=TRUE)>0])
          zmwmutations.granges.zmwindex$seqnames <- as.numeric(as.character(zmwmutations.granges.zmwindex$seqnames))

          include.zmwmutations.granges.zmwindex <- bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws]
          for(zmwindex in unique(zmwmutations.granges.zmwindex$seqnames)){
            include.zmwmutations.granges.zmwindex[[zmwindex]] <- !grepl(paste0("^",zmwmutations.granges.zmwindex[zmwmutations.granges.zmwindex$seqnames==zmwindex,"start"],"$",collapse="|"),include.zmwmutations.granges.zmwindex[[zmwindex]])
          }
          include.zmwmutations.granges.zmwindex <- mapply(function(x,y){x&y},lapply(bamfilezmw_all[[i]][[j]]$zmw.rp[includezmws],function(x){x|TRUE}),include.zmwmutations.granges.zmwindex,SIMPLIFY=FALSE)
          
          
        ##Make final mutation inclusion lists
        includemutations.fwd <- mapply(function(x,y){x&y},includemutations.fwd,include.fwdmutations.granges.zmwindex,SIMPLIFY=FALSE)
        includemutations.rev <- mapply(function(x,y){x&y},includemutations.rev,include.revmutations.granges.zmwindex,SIMPLIFY=FALSE)
        includemutations.zmw <- mapply(function(x,y){x&y},includemutations.zmw,include.zmwmutations.granges.zmwindex,SIMPLIFY=FALSE)

        #Save regional filters for later mutation frequency calculation
        #ccsindel filters are only stored in readregionalfilters, because it is read-specific and should not be used as a genome region filter in mutation frequency calculations.
        fwdmutations.ccsindels <- makeGRangesFromDataFrame(as.data.frame(fwdmutations.ccsindels)[,c("start","end","strand","chromannot","appliestoread")],seqnames="chromannot",keep.extra.columns=TRUE)
        revmutations.ccsindels <- makeGRangesFromDataFrame(as.data.frame(revmutations.ccsindels)[,c("start","end","strand","chromannot","appliestoread")],seqnames="chromannot",keep.extra.columns=TRUE)
        readregionalfilters.fwd <- c(readregionalfilters.fwd,fwdmutations.ccsindels)
        readregionalfilters.rev <- c(readregionalfilters.rev,revmutations.ccsindels)
        readregionalfilters.zmw <- c(readregionalfilters.zmw,fwdmutations.ccsindels,revmutations.ccsindels)

        #Output filtering statistics - fwd reads
        sumfiltergranges <- sum(width(fwdmutations.ccsindels))
        
        genome.granges.filteredtemp <- GenomicRanges::setdiff(genome.granges.filtered,fwdmutations.ccsindels)
        sumadditionalfiltergranges <- sum(width(GenomicRanges::setdiff(genome.granges.filtered,genome.granges.filteredtemp)))

        genome.granges.filteredtemp.ssDNA <- GenomicRanges::setdiff(genome.granges.filtered.ssDNA,fwdmutations.ccsindels)
        sumadditionalfiltergranges.ssDNA <- sum(width(GenomicRanges::setdiff(genome.granges.filtered.ssDNA,genome.granges.filteredtemp.ssDNA)))

        cat(paste("        Fwd strand filter covers:",sumfiltergranges,"/",sum(width(fwd.granges)),"total fwd bases sequenced =",round(sumfiltergranges/sum(width(fwd.granges))*100,1),"% of total fwd bases sequenced\n"))
        cat(paste("        Fwd strand filter removes:",sumadditionalfiltergranges.ssDNA,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) beyond prior filters =",round(sumadditionalfiltergranges.ssDNA/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for ssDNA mutations\n"))
        cat(paste("        Fwd strand filter removes:",sumadditionalfiltergranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) beyond prior filters =",round(sumadditionalfiltergranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for dsDNA mutations\n"))

        sumremaininggranges <- sum(width(genome.granges.filteredtemp.ssDNA))
        cat(paste("        Remaining after Fwd strand filter:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for ssDNA mutations\n"))
        sumremaininggranges <- sum(width(genome.granges.filteredtemp))
        cat(paste("        Remaining after Fwd strand filter:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for dsDNA mutations\n"))

        #Output filtering statistics - rev reads
        sumfiltergranges <- sum(width(revmutations.ccsindels))
        
        genome.granges.filteredtemp <- GenomicRanges::setdiff(genome.granges.filtered,revmutations.ccsindels)
        sumadditionalfiltergranges <- sum(width(GenomicRanges::setdiff(genome.granges.filtered,genome.granges.filteredtemp)))

        genome.granges.filteredtemp.ssDNA <- GenomicRanges::setdiff(genome.granges.filtered.ssDNA,revmutations.ccsindels)
        sumadditionalfiltergranges.ssDNA <- sum(width(GenomicRanges::setdiff(genome.granges.filtered.ssDNA,genome.granges.filteredtemp.ssDNA)))

        cat(paste("        Rev strand filter covers:",sumfiltergranges,"/",sum(width(rev.granges)),"total rev bases sequenced =",round(sumfiltergranges/sum(width(rev.granges))*100,1),"% of total rev bases sequenced\n"))
        cat(paste("        Rev strand filter removes:",sumadditionalfiltergranges.ssDNA,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) beyond prior filters =",round(sumadditionalfiltergranges.ssDNA/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for ssDNA mutations\n"))
        cat(paste("        Rev strand filter removes:",sumadditionalfiltergranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) beyond prior filters =",round(sumadditionalfiltergranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for dsDNA mutations\n"))

        sumremaininggranges <- sum(width(genome.granges.filteredtemp.ssDNA))
        cat(paste("        Remaining after Rev strand filter:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for ssDNA mutations\n"))
        sumremaininggranges <- sum(width(genome.granges.filteredtemp))
        cat(paste("        Remaining after Rev strand filter:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for dsDNA mutations\n"))

        #Output filtering statistics - zmw
        sumfiltergranges <- sum(width(GenomicRanges::union(fwdmutations.ccsindels,revmutations.ccsindels)))
        
        genome.granges.filteredtemp <- GenomicRanges::setdiff(genome.granges.filtered,fwdmutations.ccsindels)
        genome.granges.filteredtemp <- GenomicRanges::setdiff(genome.granges.filteredtemp,revmutations.ccsindels)
        sumadditionalfiltergranges <- sum(width(GenomicRanges::setdiff(genome.granges.filtered,genome.granges.filteredtemp)))

        genome.granges.filteredtemp.ssDNA <- GenomicRanges::setdiff(genome.granges.filtered.ssDNA,fwdmutations.ccsindels)
        genome.granges.filteredtemp.ssDNA <- GenomicRanges::setdiff(genome.granges.filteredtemp.ssDNA,revmutations.ccsindels)
        sumadditionalfiltergranges.ssDNA <- sum(width(GenomicRanges::setdiff(genome.granges.filtered.ssDNA,genome.granges.filteredtemp.ssDNA)))
        
        cat(paste("        Fwd+Rev strand filters covers:",sumfiltergranges,"/",sum(width(zmw.granges)),"total ZMW bases sequenced =",round(sumfiltergranges/sum(width(zmw.granges))*100,1),"% of total ZMW bases sequenced\n"))
        cat(paste("        Fwd+Rev strand filters removes:",sumadditionalfiltergranges.ssDNA,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) beyond prior filters =",round(sumadditionalfiltergranges.ssDNA/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for ssDNA mutations\n"))
        cat(paste("        Fwd+Rev strand filters removes:",sumadditionalfiltergranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) beyond prior filters =",round(sumadditionalfiltergranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for dsDNA mutations\n"))

        #Updating both genome.granges.filtered and genome.granges.filtered.ssDNA using Fwd+Rev strand CCS indel regions for simplicity, otherwise would need a different genome.granges.filtered for fwd, rev, and zmw
        genome.granges.filtered <- genome.granges.filteredtemp
        genome.granges.filtered.ssDNA <- genome.granges.filteredtemp.ssDNA

        sumremaininggranges <- sum(width(genome.granges.filtered.ssDNA))
        cat(paste("        Remaining after Fwd+Rev strand filters:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for ssDNA mutations\n"))
        sumremaininggranges <- sum(width(genome.granges.filtered))
        cat(paste("        Remaining after Fwd+Rev strand filters:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for dsDNA mutations\n"))

        #Delete temporary objects
        rm(fwdmutations.ccsindels,revmutations.ccsindels,fwdmutations.granges.zmwindex,revmutations.granges.zmwindex,zmwmutations.granges.zmwindex,include.fwdmutations.granges.zmwindex,include.revmutations.granges.zmwindex,include.zmwmutations.granges.zmwindex)
      }
      
      #Apply germline WGS filters
      cat(" Applying germline WGS filters...")

      #Make granges of final mutations.
      finalfwdmutations.granges <- fwdmutations.granges[unlist(includemutations.fwd)]
      finalrevmutations.granges <- revmutations.granges[unlist(includemutations.rev)]
      finalzmwmutations.granges <- zmwmutations.granges[unlist(includemutations.zmw)]
      
      #Make lists to save later in bamfilezmw_all_filtered the results of the relevant BAM/CRAM data used for filtering. These are updated with data as the filtering occurs below.
      includemutations.fwd.BAMTotalReads <- lapply(includemutations.fwd,function(x){rep(NA,length(x))})
      includemutations.fwd.BAMVariantReads <- lapply(includemutations.fwd,function(x){rep(NA,length(x))})
      includemutations.fwd.BAMVAF <- lapply(includemutations.fwd,function(x){rep(NA,length(x))})
      includemutations.rev.BAMTotalReads <- lapply(includemutations.rev,function(x){rep(NA,length(x))})
      includemutations.rev.BAMVariantReads <- lapply(includemutations.rev,function(x){rep(NA,length(x))})
      includemutations.rev.BAMVAF <- lapply(includemutations.rev,function(x){rep(NA,length(x))})
      includemutations.zmw.BAMTotalReads <- lapply(includemutations.zmw,function(x){rep(NA,length(x))})
      includemutations.zmw.BAMVariantReads <- lapply(includemutations.zmw,function(x){rep(NA,length(x))})
      includemutations.zmw.BAMVAF <- lapply(includemutations.zmw,function(x){rep(NA,length(x))})

      #Make file of final combined granges for pileup.
      finalmutations.granges <- sort(c(finalfwdmutations.granges,finalrevmutations.granges,finalzmwmutations.granges))
      
      if(length(finalmutations.granges)>0){
	      finalmutations.regionsfile <- tempfile()
	      write.table(unique(cbind(as.character(seqnames(finalmutations.granges)),start(finalmutations.granges))),finalmutations.regionsfile,quote=F,row.names=F,col.names=F,sep="\t")
	
	      #Run bcftools mpileup. For Illumina germline reference, use option -A to include orphan reads, -B to disable BAQ because this excludes some false positive alignments which we actually want in filtering, lower base quality threshold to 11 (from default of 13), and exclude PCR duplicates. For PacBio germline reference, use presets for pacbio-ccs recommended by bcftools mpileup except exclude supplementary alignments, without -D, and add -B to disable base quality (it leads to missed variants for some reason).
	      pileuptempfile <- tempfile()
	      if(bamfilters[[k]][[i]]$bamreftype=="Illumina"){
	      	system(paste("/bin/bash -c",shQuote(paste(bcftools_bin,"mpileup -I -A -B -Q 11 --ns 3328 -d 999999 -a \"INFO/AD\" -f",referenceinfo$fastaref,"-R",finalmutations.regionsfile,bamfilters[[k]][[i]]$bamfile,"2>/dev/null | grep -v ^\\#",">",pileuptempfile))))
	      }else if(bamfilters[[k]][[i]]$bamreftype=="PacBio"){
	      	system(paste("/bin/bash -c",shQuote(paste(bcftools_bin,"mpileup -I -B -Q 5 --ns 3328 -d 999999 --max-BQ 50 -F0.1 -o25 -e1 -a \"INFO/AD\" -f",referenceinfo$fastaref,"-R",finalmutations.regionsfile,bamfilters[[k]][[i]]$bamfile,"2>/dev/null | grep -v ^\\#",">",pileuptempfile))))
	      }
	      finalmutations.pileupvcf <- read.table(pileuptempfile)
	      file.remove(finalmutations.regionsfile,pileuptempfile)
	
	      #Extract pileup information
	      finalmutations.pileupvcf.qn <- str_split(finalmutations.pileupvcf[,5],",")
	      finalmutations.pileupvcf.AD <- str_split(sub(".*AD=(.*?);.*","\\1",finalmutations.pileupvcf[,8]),",")
	      finalmutations.pileupvcf.AD <- lapply(finalmutations.pileupvcf.AD,function(x){as.numeric(x)})
	      finalmutations.pileupvcf.totalreads <- unlist(lapply(finalmutations.pileupvcf.AD,function(x){sum(x)}))
	      finalmutations.pileupvcf.VAF <- lapply(finalmutations.pileupvcf.AD,function(x){x/sum(x)})
      }

      #Filter out variants that don't meet thresholds - fwd ssDNA mutations
      for(grangesindex in seq_along(finalfwdmutations.granges)){
        pileupvcf.rownum <- which(finalmutations.pileupvcf[,1]==seqnames(finalfwdmutations.granges)[grangesindex] & finalmutations.pileupvcf[,2]==start(finalfwdmutations.granges)[grangesindex])
        if(isEmpty(pileupvcf.rownum)){ #Mutation position has 0 coverage. Set as if it was detected as germline so it is not falsely considered as a somatic mutation. It is also necessary to skip the remaining code of this section if this is the case, because otherwise it leads to a script error when pileupvcf.rownum is NULL.
          includemutations.fwd[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- FALSE
          includemutations.fwd.BAMTotalReads[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- 0
          includemutations.fwd.BAMVariantReads[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- 0
          includemutations.fwd.BAMVAF[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- 0
        }else{
          pileupvcf.altnum <- grep(finalfwdmutations.granges$qn[grangesindex],finalmutations.pileupvcf.qn[[pileupvcf.rownum]])
          if(isEmpty(pileupvcf.altnum)){ #Mutation base not in pileup
            includemutations.fwd[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum] >=bamfilters[[k]][[i]]$minBAMTotalReads
            includemutations.fwd.BAMTotalReads[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum]
            includemutations.fwd.BAMVariantReads[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- 0
            includemutations.fwd.BAMVAF[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- 0
          }else{
            includemutations.fwd[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum] >=bamfilters[[k]][[i]]$minBAMTotalReads & finalmutations.pileupvcf.AD[[pileupvcf.rownum]][pileupvcf.altnum+1]<=bamfilters[[k]][[i]]$maxBAMVariantReads & finalmutations.pileupvcf.VAF[[pileupvcf.rownum]][pileupvcf.altnum+1]<=bamfilters[[k]][[i]]$maxBAMVAF
            includemutations.fwd.BAMTotalReads[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum]
            includemutations.fwd.BAMVariantReads[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.AD[[pileupvcf.rownum]][pileupvcf.altnum+1]
            includemutations.fwd.BAMVAF[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.VAF[[pileupvcf.rownum]][pileupvcf.altnum+1]
          }
        }
      }

      #Filter out variants that don't meet thresholds - rev ssDNA mutations
      for(grangesindex in seq_along(finalrevmutations.granges)){
        pileupvcf.rownum <- which(finalmutations.pileupvcf[,1]==seqnames(finalrevmutations.granges)[grangesindex] & finalmutations.pileupvcf[,2]==start(finalrevmutations.granges)[grangesindex])
        if(isEmpty(pileupvcf.rownum)){ #Mutation position has 0 coverage. Set as if it was detected as germline so it is not falsely considered as a somatic mutation. It is also necessary to skip the remaining code of this section if this is the case, because otherwise it leads to a script error when pileupvcf.rownum is NULL.
          includemutations.rev[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- FALSE
          includemutations.rev.BAMTotalReads[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- 0
          includemutations.rev.BAMVariantReads[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- 0
          includemutations.rev.BAMVAF[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- 0
        }else{
          pileupvcf.altnum <- grep(finalrevmutations.granges$qn[grangesindex],finalmutations.pileupvcf.qn[[pileupvcf.rownum]])
          if(isEmpty(pileupvcf.altnum)){ #Mutation base not in pileup
            includemutations.rev[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum] >=bamfilters[[k]][[i]]$minBAMTotalReads
            includemutations.rev.BAMTotalReads[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum]
            includemutations.rev.BAMVariantReads[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- 0
            includemutations.rev.BAMVAF[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- 0
          }else{
            includemutations.rev[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum] >=bamfilters[[k]][[i]]$minBAMTotalReads & finalmutations.pileupvcf.AD[[pileupvcf.rownum]][pileupvcf.altnum+1]<=bamfilters[[k]][[i]]$maxBAMVariantReads & finalmutations.pileupvcf.VAF[[pileupvcf.rownum]][pileupvcf.altnum+1]<=bamfilters[[k]][[i]]$maxBAMVAF
            includemutations.rev.BAMTotalReads[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum]
            includemutations.rev.BAMVariantReads[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.AD[[pileupvcf.rownum]][pileupvcf.altnum+1]
            includemutations.rev.BAMVAF[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.VAF[[pileupvcf.rownum]][pileupvcf.altnum+1]
          }
        }
      }

      #Filter out variants that don't meet thresholds - zmw dsDNA mutations
      for(grangesindex in seq_along(finalzmwmutations.granges)){
        pileupvcf.rownum <- which(finalmutations.pileupvcf[,1]==seqnames(finalzmwmutations.granges)[grangesindex] & finalmutations.pileupvcf[,2]==start(finalzmwmutations.granges)[grangesindex])
        if(isEmpty(pileupvcf.rownum)){ #Mutation position has 0 coverage. Set as if it was detected so it is not falsely considered as a somatic mutation. It is also necessary to skip the remaining code of this section if this is the case, because otherwise it leads to a script error when pileupvcf.rownum is NULL.
          includemutations.zmw[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- FALSE
          includemutations.zmw.BAMTotalReads[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- 0
          includemutations.zmw.BAMVariantReads[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- 0
          includemutations.zmw.BAMVAF[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- 0
        }else{
          pileupvcf.altnum <- grep(finalzmwmutations.granges$qn[grangesindex],finalmutations.pileupvcf.qn[[pileupvcf.rownum]])
          if(isEmpty(pileupvcf.altnum)){ #Mutation base not in pileup
            includemutations.zmw[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum] >=bamfilters[[k]][[i]]$minBAMTotalReads
            includemutations.zmw.BAMTotalReads[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum]
            includemutations.zmw.BAMVariantReads[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- 0
            includemutations.zmw.BAMVAF[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- 0
          }else{
            includemutations.zmw[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum] >=bamfilters[[k]][[i]]$minBAMTotalReads & finalmutations.pileupvcf.AD[[pileupvcf.rownum]][pileupvcf.altnum+1]<=bamfilters[[k]][[i]]$maxBAMVariantReads & finalmutations.pileupvcf.VAF[[pileupvcf.rownum]][pileupvcf.altnum+1]<=bamfilters[[k]][[i]]$maxBAMVAF
            includemutations.zmw.BAMTotalReads[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.totalreads[pileupvcf.rownum]
            includemutations.zmw.BAMVariantReads[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.AD[[pileupvcf.rownum]][pileupvcf.altnum+1]
            includemutations.zmw.BAMVAF[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalmutations.pileupvcf.VAF[[pileupvcf.rownum]][pileupvcf.altnum+1]
          }
        }
      }

      genomecovfilter.granges <- bw.granges.thresholding(bamfilters[[k]][[i]]$WGScoveragebigwig,referenceinfo$chrsizes,1,bamfilters[[k]][[i]]$minBAMTotalReads)
      genomecovfilter.granges <- genomecovfilter.granges[["lt"]]
      genomecovfilter.granges <- filter.granges.by.chroms(genomecovfilter.granges,chrs_to_analyze)

      #Save regional filters for later mutation frequency calculation
      readregionalfilters.fwd <- c(readregionalfilters.fwd,genomecovfilter.granges)
      readregionalfilters.rev <- c(readregionalfilters.rev,genomecovfilter.granges)
      readregionalfilters.zmw <- c(readregionalfilters.zmw,genomecovfilter.granges)
      genomeregionalfilters.fwd <- c(genomeregionalfilters.fwd,genomecovfilter.granges)
      genomeregionalfilters.rev <- c(genomeregionalfilters.rev,genomecovfilter.granges)
      genomeregionalfilters.zmw <- c(genomeregionalfilters.zmw,genomecovfilter.granges)

      #Output filtering statistics
      sumfiltergranges <- sum(width(genomecovfilter.granges))
      
      genome.granges.filteredtemp <- GenomicRanges::setdiff(genome.granges.filtered,genomecovfilter.granges)
      sumadditionalfiltergranges <- sum(width(GenomicRanges::setdiff(genome.granges.filtered,genome.granges.filteredtemp)))
      
      genome.granges.filteredtemp.ssDNA <- GenomicRanges::setdiff(genome.granges.filtered.ssDNA,genomecovfilter.granges)
      sumadditionalfiltergranges.ssDNA <- sum(width(GenomicRanges::setdiff(genome.granges.filtered.ssDNA,genome.granges.filteredtemp.ssDNA)))

      cat("...DONE\n")
      cat(paste("        Filter covers:",sumfiltergranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumfiltergranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions)\n"))
      cat(paste("        Filter removes:",sumadditionalfiltergranges.ssDNA,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) beyond prior filters =",round(sumadditionalfiltergranges.ssDNA/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for ssDNA mutations\n"))
      cat(paste("        Filter removes:",sumadditionalfiltergranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) beyond prior filters =",round(sumadditionalfiltergranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for dsDNA mutations\n"))

      genome.granges.filtered <- genome.granges.filteredtemp
      genome.granges.filtered.ssDNA <- genome.granges.filteredtemp.ssDNA
      sumremaininggranges <- sum(width(genome.granges.filtered.ssDNA))
      cat(paste("        Remaining:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for ssDNA mutations\n"))
      sumremaininggranges <- sum(width(genome.granges.filtered))
      cat(paste("        Remaining:",sumremaininggranges,"/",sum(width(genome.noN.granges)),"total genome bases (excluding 'N' regions) =",round(sumremaininggranges/sum(width(genome.noN.granges))*100,1),"% of total genome (excluding 'N' regions) for dsDNA mutations\n"))
      
      #Remove temporary objects.
      rm(finalmutations.granges,genomecovfilter.granges,finalfwdmutations.granges,finalrevmutations.granges,finalzmwmutations.granges)
      
      ##Apply subread pileup filters: Filters out ssDNA and ZMW (dsDNA) mutations that do not have a minimum number and fraction of subreads with the mutation (minZMWsubreadsVariantReads, minZMWsubreadsVAF, minssDNAsubreadsVariantReads, minssDNAsubreadsVAF), or that do not have a minimum overlapping coverage of subreads relative to the ZMW's ec (effective coverage) (minsubreads_cvg_fraction filter).
      cat(" Applying subread filters...\n")
      	
    	#Make lists to save later in bamfilezmw_all_filtered the results of the relevant subread data used for filtering. These are updated with data as the filtering occurs below.
      includemutations.zmw.fwdsubreadsVariantReads <- lapply(includemutations.zmw,function(x){rep(NA,length(x))})
      includemutations.zmw.revsubreadsVariantReads <- lapply(includemutations.zmw,function(x){rep(NA,length(x))})
      includemutations.zmw.fwdsubreadsVAF <- lapply(includemutations.zmw,function(x){rep(NA,length(x))})
      includemutations.zmw.revsubreadsVAF <- lapply(includemutations.zmw,function(x){rep(NA,length(x))})
      includemutations.zmw.fwdsubreadscvgfraction <- lapply(includemutations.zmw,function(x){rep(NA,length(x))})
      includemutations.zmw.revsubreadscvgfraction <- lapply(includemutations.zmw,function(x){rep(NA,length(x))})
      includemutations.fwd.ssDNAsubreadsVariantReads <- lapply(includemutations.fwd,function(x){rep(NA,length(x))})
      includemutations.fwd.ssDNAsubreadsVAF <- lapply(includemutations.fwd,function(x){rep(NA,length(x))})
      includemutations.fwd.ssDNAsubreadscvgfraction <- lapply(includemutations.fwd,function(x){rep(NA,length(x))})
      includemutations.rev.ssDNAsubreadsVariantReads <- lapply(includemutations.rev,function(x){rep(NA,length(x))})
      includemutations.rev.ssDNAsubreadsVAF <- lapply(includemutations.rev,function(x){rep(NA,length(x))})
      includemutations.rev.ssDNAsubreadscvgfraction <- lapply(includemutations.rev,function(x){rep(NA,length(x))})
	    zmwstoretrieve <- c()
	    
      #Extract and align subreads of ZMWs with mutations
      cat("  Extracting and aligning mutation subreads...")
      
      #Make granges of final Fwd and Rev ssDNA mutations.
      finalfwdmutations.granges <- fwdmutations.granges[unlist(includemutations.fwd)]
      finalrevmutations.granges <- revmutations.granges[unlist(includemutations.rev)]
      finalzmwmutations.granges <- zmwmutations.granges[unlist(includemutations.zmw)]
      
      #Make list of all ZMW holes that have a Fwd, Rev, or ZMW (dsDNA) mutation.
      zmwstoretrieve <- unique(sort(as.numeric(c(finalfwdmutations.granges$zmw,finalrevmutations.granges$zmw,finalzmwmutations.granges$zmw))))
      
      if(length(zmwstoretrieve)==0){
        bamfilezmw_filtered[[i]][[j]][[k]] <- list()
      	cat(" No mutations, skipping subread filtering!\n")
      }else{
	      zmwstoretrieve_file <- tempfile(pattern=".",tmpdir=getwd())
	      write.table(zmwstoretrieve,zmwstoretrieve_file,quote=F,row.names=F,col.names=F)
	      
	      #Extract subreads of ZMW holes and align, via SLURM job
	      tmpsubreadsbam <- paste0(tempfile(pattern=".",tmpdir=getwd()),".bam")
	      subreadsalignedbam <- paste0(tempfile(pattern=".",tmpdir=getwd()),".bam")
	      cmd <- tempfile(pattern=".",tmpdir=getwd())
	      cat(paste("source",condabase_script,"; conda activate",conda_pbbioconda_env,"; zmwfilter --include",zmwstoretrieve_file,subreadsbam,tmpsubreadsbam,"; pbmm2 align --log-level INFO -j 4",referenceinfo$genomemmiindex,tmpsubreadsbam,subreadsalignedbam,"--preset SUBREAD --sort 2>/dev/null; pbindex",subreadsalignedbam,"; rm",tmpsubreadsbam),file=cmd)
	    	Sys.chmod(cmd,mode="700")
	      subreadsjobindex <- system(paste("sbatch",yaml.config$slurm_add_options,"-t 720 --mem=64G -c 4 -o /dev/null --parsable --wrap \"$HIDEF_SINGULARITY_WRAPPER",cmd,"; rm",cmd,"\""),intern=T)
	      
	        #Wait for job to finish
	    	while(!(system(paste("sacct -o State -j",subreadsjobindex,"| tail -n 1 | awk '{print $1}'"),intern=T) %in% c("COMPLETED","FAILED"))){
	        Sys.sleep(60) #Check SLURM job status every 1 minute
	      }
	      if(system(paste("sacct -o State -j",subreadsjobindex,"| tail -n 1 | awk '{print $1}'"),intern=T)=="FAILED"){
	        cat("\n...SLURM job failed!")
	        quit()
	      }
	      
	      cat("DONE\n")
	      
	      cat("  Performing subread filtering...")
	      #Reformat aligned subreads BAM file with each ZMW hole having a different read group and sample name, so that bcftools mpileup counts reads separately for each ZMW.
	      subreadstempfile <- tempfile()
	      invisible(system(paste("/bin/bash -c",shQuote(paste(yaml.config$samtools_bin,"view -H",subreadsalignedbam,"| grep -v \"^@RG\" >",paste0(subreadstempfile,".subreads.aligned.bam.header"),"; for i in `cat",zmwstoretrieve_file,"`; do printf \"@RG\tID:$i\tSM:$i\n\"; done >>",paste0(subreadstempfile,".subreads.aligned.bam.header"),";",yaml.config$samtools_bin,"view",subreadsalignedbam,"| awk '{split($1,zmw,\"/\"); sub(\"RG:Z:[[:alnum:]]+\t\",\"RG:Z:\" zmw[2] \"\t\"); print}' >",paste0(subreadstempfile,".subreads.aligned.bam.sam"),"; cat",paste0(subreadstempfile,".subreads.aligned.bam.header"),paste0(subreadstempfile,".subreads.aligned.bam.sam"),"|",yaml.config$samtools_bin,"view -b - >",paste0(subreadstempfile,".subreads.aligned.rg.bam"),";",yaml.config$samtools_bin,"index",paste0(subreadstempfile,".subreads.aligned.rg.bam")))),intern=T))
	      
	      #Run bcftools mpileup on subreads, counting reads for each strand separately, because ssDNA subread filtering is performed using only subreads belonging to the strand on which the ssDNA mutation was called (or both strands for dsDNA mutations). Use option -A to include orphan reads, -B to disable BAQ because this excludes some false positive alignments which we actually want in filtering, and set -Q 0 in order to not apply any base quality threshold because base quality in subreads = 0. Also, exclude subreads that align as supplementary alignments (-F 2048).
	      finalmutations.granges <- sort(c(finalfwdmutations.granges,finalrevmutations.granges,finalzmwmutations.granges))
	      finalmutations.regionsfile <- tempfile()
	      write.table(unique(cbind(as.character(seqnames(finalmutations.granges)),start(finalmutations.granges))),finalmutations.regionsfile,quote=F,row.names=F,col.names=F,sep="\t")
	      
	      pileuptempfile <- tempfile()
	      invisible(system(paste("/bin/bash -c",shQuote(paste(yaml.config$bcftools_bin,"mpileup --ff 2048 -I -A -B -Q 0 -d 9999999999 -a \"INFO/AD,FORMAT/DP,FORMAT/ADF,FORMAT/ADR\" -f",referenceinfo$fastaref,"-R",finalmutations.regionsfile,paste0(subreadstempfile,".subreads.aligned.rg.bam"),"2>/dev/null | grep -v ^\\#\\#",">",pileuptempfile))),intern=T))
	      pileupvcf <- makeGRangesFromDataFrame(read.delim(pileuptempfile),seqnames.field="X.CHROM",start.field="POS",end.field="POS",keep.extra.columns=TRUE)
	      
	      #Define function to annotate granges with relevant data. Return granges unchanged if its length=0.
	      annotate.gr.withpileup <- function(gr,annotationvcf){
	      	if(length(gr)==0){
	      		return(gr)
	      	}
	        annotationvcf.qn <- strsplit(annotationvcf$ALT,",")
	      
	        FORMATsplit <- strsplit(annotationvcf$FORMAT,":")
	        ADFindex <- unlist(lapply(FORMATsplit,function(x){grep("ADF",x)}))
	        ADRindex <- unlist(lapply(FORMATsplit,function(x){grep("ADR",x)}))
	      
	        vcfhits <- findOverlaps(gr,annotationvcf,type="equal",select="first")
	        
	        #Rarely, no subreads will be aligned to site of variant called in CCS (subreads poor quality and don't align well to site of the CCS constructed by those subreads). However, bcftools mpileup will not output a record for these sites. Therefore, we identify these sites, and set their subread coverage values to 0 so they are filtered out, because their subread coverage is indeed 0.
	        vcfhits.na <- which(is.na(vcfhits))
	        vcfhits <- vcfhits[!is.na(vcfhits)]
	
	        zmwcolmatch <- match(paste0("X",gr$zmw[setdiff(seq_along(gr$zmw),vcfhits.na)]),colnames(mcols(annotationvcf)))
	        zmwcol <- as.data.frame(mcols(annotationvcf))[cbind(vcfhits,zmwcolmatch)]
	        zmwcolsplit <- strsplit(zmwcol,":")
	        
	        ADFvalues <- lapply(strsplit(mapply(function(x,y){x[y]},x=zmwcolsplit,y=ADFindex[vcfhits]),","),as.numeric)
	        ADRvalues <- lapply(strsplit(mapply(function(x,y){x[y]},x=zmwcolsplit,y=ADRindex[vcfhits]),","),as.numeric)
	        
	        altindex <- mapply(function(x,y){grep(x,y)+1},x=gr$qn[setdiff(seq_along(gr$qn),vcfhits.na)],y=annotationvcf.qn[vcfhits],USE.NAMES=F)
	        
	        #Function to add zero values at specific vector indices (postoadd), for sites without an mpileup record
	        addzerovalues <- function(x,postoadd){
	          postoadd <- postoadd-1
	          x <- rep(x, (1:2)[(1:length(x) %in% postoadd) + 1])
	          x[postoadd + 1:length(postoadd)] <- 0
	          return(x)
	        }
	        
	        gr$fwdsubreadsVariantReads <- addzerovalues(mapply(function(x,y){x[y]},x=ADFvalues,y=altindex),vcfhits.na)
	        gr$revsubreadsVariantReads <- addzerovalues(mapply(function(x,y){x[y]},x=ADRvalues,y=altindex),vcfhits.na)
	        gr$fwdsubreadsVAF <- addzerovalues(mapply(function(x,y){x[y]/sum(x)},x=ADFvalues,y=altindex),vcfhits.na)
	        gr$revsubreadsVAF <- addzerovalues(mapply(function(x,y){x[y]/sum(x)},x=ADRvalues,y=altindex),vcfhits.na)
	        
	        #A small number of variants may not be in the pileup due to differences in positions of variants in subread alignments versus their positions in ccs consensus call and ccs consensus read alignment, or they may have zero effective read coverage in the pileup (produces NA when dividing by 0 in VAF calculations). Set these to 0 variantreads and VAF (i.e. filter them out as not detected). All of the ones examined manually are artifacts mainly of ccs calling variants with minimal to no subread support.
	        for(annotatedcolumn in c("fwdsubreadsVariantReads","revsubreadsVariantReads","fwdsubreadsVAF","revsubreadsVAF")){
	          varsnotfound <- unlist(lapply(mcols(gr)[[annotatedcolumn]],length))==0
	          mcols(gr)[[annotatedcolumn]][varsnotfound] <- rep(list(0),length(which(varsnotfound==TRUE)))
	          
	          varsnotfound <- unlist(lapply(mcols(gr)[[annotatedcolumn]],is.na))
	          mcols(gr)[[annotatedcolumn]][varsnotfound] <- rep(list(0),length(which(varsnotfound==TRUE)))
	          mcols(gr)[[annotatedcolumn]] <- unlist(mcols(gr)[[annotatedcolumn]])
	        }
	        
	        return(gr)
	      }
	      
	      #Annotate final-fwd/rev/zmw-mutations.granges with pileup data
	      finalfwdmutations.granges <- annotate.gr.withpileup(finalfwdmutations.granges,pileupvcf)
	      finalrevmutations.granges <- annotate.gr.withpileup(finalrevmutations.granges,pileupvcf)
	      finalzmwmutations.granges <- annotate.gr.withpileup(finalzmwmutations.granges,pileupvcf)
	      
	      #Annotate final-fwd/rev/zmw-mutations.granges with subreads coverage fraction (specific to each ZMW), for each strand separately, for use in minsubreads_cvg_fraction filtering
	      
	        #Load subreads BAM file
	        subreadsbamfile <- scanBam(subreadsalignedbam,param=ScanBamParam(what=c("rname","pos","isize","strand"),tag=c("zm")))[[1]]
	        subreadsbamfile.gr <- makeGRangesFromDataFrame(data.frame(chrom=subreadsbamfile$rname,start=subreadsbamfile$pos,end=subreadsbamfile$pos+subreadsbamfile$isize-1,strand=subreadsbamfile$strand,zm=subreadsbamfile$tag$zm),keep.extra.columns=TRUE)
	        fwdsubreadsbamfile.gr <- subreadsbamfile.gr[strand(subreadsbamfile.gr)=="+"]
	        revsubreadsbamfile.gr <- subreadsbamfile.gr[strand(subreadsbamfile.gr)=="-"]
	      
	        #Annotate number of subreads covering each mutation. Return input granges unchanged if length=0.
	        annotatesubreadscvg <- function(gr,subreads.gr,annotationcolumn){
	        	if(length(gr)==0){
	      			return(gr)
	      		}
	          subreadscvg.temp <- join_overlap_inner(gr,subreads.gr) %>% filter(zmw==zm) %>% group_by(seqnames,start,end,zmw) %>% summarise(n=n()) %>% as_granges()
	          names(gr) <- paste(seqnames(gr),start(gr),gr$zmw,sep="_")
	          names(subreadscvg.temp) <- paste(seqnames(subreadscvg.temp),start(subreadscvg.temp),subreadscvg.temp$zmw,sep="_")
	          
	          mcols(gr)[[annotationcolumn]] <- rep(0,length(gr))
	          mcols(gr)[names(subreadscvg.temp),annotationcolumn] <- subreadscvg.temp$n
	          names(gr) <- NULL
	          mcols(gr)[[annotationcolumn]][is.na(mcols(gr)[[annotationcolumn]])] <- 0
	          return(gr)
	        }
	      
	        finalfwdmutations.granges <- annotatesubreadscvg(finalfwdmutations.granges,fwdsubreadsbamfile.gr,"fwdsubreadscvg")
	        finalrevmutations.granges <- annotatesubreadscvg(finalrevmutations.granges,revsubreadsbamfile.gr,"revsubreadscvg")
	        finalzmwmutations.granges <- annotatesubreadscvg(finalzmwmutations.granges,fwdsubreadsbamfile.gr,"fwdsubreadscvg")
	        finalzmwmutations.granges <- annotatesubreadscvg(finalzmwmutations.granges,revsubreadsbamfile.gr,"revsubreadscvg")
	      
	        #Annotate total number of subreads for each mutation's zmw
	        fwdsubreadsperzmw <- table(fwdsubreadsbamfile.gr$zm)
	        revsubreadsperzmw <- table(revsubreadsbamfile.gr$zm)
	        
	        finalfwdmutations.granges$totalfwdsubreads <- as.numeric(fwdsubreadsperzmw[match(finalfwdmutations.granges$zmw,names(fwdsubreadsperzmw))])
	        finalrevmutations.granges$totalrevsubreads <- as.numeric(revsubreadsperzmw[match(finalrevmutations.granges$zmw,names(revsubreadsperzmw))])
	        finalzmwmutations.granges$totalfwdsubreads <- as.numeric(fwdsubreadsperzmw[match(finalzmwmutations.granges$zmw,names(fwdsubreadsperzmw))])
	        finalzmwmutations.granges$totalrevsubreads <- as.numeric(revsubreadsperzmw[match(finalzmwmutations.granges$zmw,names(revsubreadsperzmw))])
	      
	        #Calculate subread coverage fraction
	        finalfwdmutations.granges$fwdsubreadscvgfraction <- finalfwdmutations.granges$fwdsubreadscvg/finalfwdmutations.granges$totalfwdsubreads
	        finalrevmutations.granges$revsubreadscvgfraction <- finalrevmutations.granges$revsubreadscvg/finalrevmutations.granges$totalrevsubreads
	        finalzmwmutations.granges$fwdsubreadscvgfraction <- finalzmwmutations.granges$fwdsubreadscvg/finalzmwmutations.granges$totalfwdsubreads
	        finalzmwmutations.granges$revsubreadscvgfraction <- finalzmwmutations.granges$revsubreadscvg/finalzmwmutations.granges$totalrevsubreads
	      
	        #Set to 0 for Na values (i.e. if total reads is 0, leads to division by zero, which leads to NA)
	        finalfwdmutations.granges$fwdsubreadscvgfraction[is.na(finalfwdmutations.granges$fwdsubreadscvgfraction)] <- 0
	        finalrevmutations.granges$revsubreadscvgfraction[is.na(finalrevmutations.granges$revsubreadscvgfraction)] <- 0
	        finalzmwmutations.granges$fwdsubreadscvgfraction[is.na(finalzmwmutations.granges$fwdsubreadscvgfraction)] <- 0
	        finalzmwmutations.granges$revsubreadscvgfraction[is.na(finalzmwmutations.granges$revsubreadscvgfraction)] <- 0
	      
	      #Remove temporary objects
	      rm(pileupvcf,finalmutations.granges,subreadsbamfile,subreadsbamfile.gr,fwdsubreadsbamfile.gr,revsubreadsbamfile.gr,fwdsubreadsperzmw,revsubreadsperzmw)
	      
	      #Perform filtering
	        #Fwd
	      for(grangesindex in seq_along(finalfwdmutations.granges)){
	        includemutations.fwd[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- finalfwdmutations.granges$fwdsubreadsVariantReads[grangesindex]>=thresholds[k,"minssDNAsubreadsVariantReads"] & finalfwdmutations.granges$fwdsubreadsVAF[grangesindex]>=thresholds[k,"minssDNAsubreadsVAF"] & finalfwdmutations.granges$fwdsubreadscvgfraction[grangesindex]>=thresholds[k,"minsubreads_cvg_fraction"]
	        
	        includemutations.fwd.ssDNAsubreadsVariantReads[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- finalfwdmutations.granges$fwdsubreadsVariantReads[grangesindex]
		  	  includemutations.fwd.ssDNAsubreadsVAF[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- finalfwdmutations.granges$fwdsubreadsVAF[grangesindex]
		  	  includemutations.fwd.ssDNAsubreadscvgfraction[[finalfwdmutations.granges$zmwindex[grangesindex]]][finalfwdmutations.granges$variantindex[grangesindex]] <- finalfwdmutations.granges$fwdsubreadscvgfraction[grangesindex]
	      }
	      
	        #Rev
	      for(grangesindex in seq_along(finalrevmutations.granges)){
	        includemutations.rev[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- finalrevmutations.granges$revsubreadsVariantReads[grangesindex]>=thresholds[k,"minssDNAsubreadsVariantReads"] & finalrevmutations.granges$revsubreadsVAF[grangesindex]>=thresholds[k,"minssDNAsubreadsVAF"] & finalrevmutations.granges$revsubreadscvgfraction[grangesindex]>=thresholds[k,"minsubreads_cvg_fraction"]
	        
	        includemutations.rev.ssDNAsubreadsVariantReads[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- finalrevmutations.granges$revsubreadsVariantReads[grangesindex]
		  	  includemutations.rev.ssDNAsubreadsVAF[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- finalrevmutations.granges$revsubreadsVAF[grangesindex]
		  	  includemutations.rev.ssDNAsubreadscvgfraction[[finalrevmutations.granges$zmwindex[grangesindex]]][finalrevmutations.granges$variantindex[grangesindex]] <- finalrevmutations.granges$revsubreadscvgfraction[grangesindex]
	      }
	
	        #ZMW
	      for(grangesindex in seq_along(finalzmwmutations.granges)){
	        includemutations.zmw[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalzmwmutations.granges$fwdsubreadsVariantReads[grangesindex]>=thresholds[k,"minZMWsubreadsVariantReads"] & finalzmwmutations.granges$fwdsubreadsVAF[grangesindex]>=thresholds[k,"minZMWsubreadsVAF"] & finalzmwmutations.granges$revsubreadsVariantReads[grangesindex]>=thresholds[k,"minZMWsubreadsVariantReads"] & finalzmwmutations.granges$revsubreadsVAF[grangesindex]>=thresholds[k,"minZMWsubreadsVAF"] & finalzmwmutations.granges$fwdsubreadscvgfraction[grangesindex]>=thresholds[k,"minsubreads_cvg_fraction"] & finalzmwmutations.granges$revsubreadscvgfraction[grangesindex]>=thresholds[k,"minsubreads_cvg_fraction"]
	        
	        includemutations.zmw.fwdsubreadsVariantReads[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalzmwmutations.granges$fwdsubreadsVariantReads[grangesindex]
	        includemutations.zmw.fwdsubreadsVAF[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalzmwmutations.granges$fwdsubreadsVAF[grangesindex]
	        includemutations.zmw.revsubreadsVariantReads[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalzmwmutations.granges$revsubreadsVariantReads[grangesindex]
	        includemutations.zmw.revsubreadsVAF[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalzmwmutations.granges$revsubreadsVAF[grangesindex]
	        includemutations.zmw.fwdsubreadscvgfraction[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalzmwmutations.granges$fwdsubreadscvgfraction[grangesindex]
	        includemutations.zmw.revsubreadscvgfraction[[finalzmwmutations.granges$zmwindex[grangesindex]]][finalzmwmutations.granges$variantindex[grangesindex]] <- finalzmwmutations.granges$revsubreadscvgfraction[grangesindex]
	      }
	      
	      #Delete temporary objects/files
	      rm(finalfwdmutations.granges,finalrevmutations.granges,finalzmwmutations.granges)
	      file.remove(zmwstoretrieve_file,paste0(subreadstempfile,".subreads.aligned.bam.header"),paste0(subreadstempfile,".subreads.aligned.bam.sam"),paste0(subreadstempfile,".subreads.aligned.rg.bam"),finalmutations.regionsfile,pileuptempfile)
	      cat("DONE\n")
	      
	      cat("  Subread filters...DONE\n")
      
	      ##Update zmwstoretrieve and extract mutation subreads BAM with final ZMW holes that have a Fwd, Rev, or ZMW (dsDNA) mutation, because some of these were filtered out by the prior subreads filters.
	      finalfwdmutations.granges <- fwdmutations.granges[unlist(includemutations.fwd)]
	      finalrevmutations.granges <- revmutations.granges[unlist(includemutations.rev)]
	      finalzmwmutations.granges <- zmwmutations.granges[unlist(includemutations.zmw)]
	      
	      zmwstoretrieve <- unique(sort(as.numeric(c(finalfwdmutations.granges$zmw,finalrevmutations.granges$zmw,finalzmwmutations.granges$zmw))))
	      zmwstoretrieve_file <- tempfile(pattern=".",tmpdir=getwd())
	      write.table(zmwstoretrieve,zmwstoretrieve_file,quote=F,row.names=F,col.names=F)
	      
	      tmpsubreadsbam <- paste0(tempfile(pattern=".",tmpdir=getwd()),".bam")
	      system(paste("/bin/bash -c",shQuote(paste("source",condabase_script,"; conda activate",conda_pbbioconda_env,"; zmwfilter --include",zmwstoretrieve_file,subreadsalignedbam,tmpsubreadsbam,"; mv",tmpsubreadsbam,paste0(outputpath,"/",analysisoutput_basename,".",i,".",j,".",k,".subreads.aligned.bam"),";",samtools_bin,"index",paste0(outputpath,"/",analysisoutput_basename,".",i,".",j,".",k,".subreads.aligned.bam"),"; pbindex",paste0(outputpath,"/",analysisoutput_basename,".",i,".",j,".",k,".subreads.aligned.bam")))),intern=T)
	      
	      rm(finalfwdmutations.granges,finalrevmutations.granges,finalzmwmutations.granges,fwdmutations.granges,revmutations.granges,zmwmutations.granges)
	      file.remove(zmwstoretrieve_file,subreadsalignedbam,paste0(subreadsalignedbam,".bai"),paste0(subreadsalignedbam,".pbi"))
      }
      
      #Make bamfilezmw_filtered that saves the yaml configuration relevant to mutation filtering, list of filtered ZMWs and mutations, region filters, BAM/CRAM germline reference data, subreads corresponding to ZMW holes with mutations, subread filtering data, and trinucleotide distributions of interrogated bases.
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["bigwig_cachedir"]] <- bigwig_cachedir
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["bamfilezmw_all_filename"]] <- bamfilezmw_all_filename
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["subreads_filename"]] <- subreads_filename
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["analysisoutput_path"]] <- analysisoutput_path
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["analysisoutput_basename"]] <- analysisoutput_basename
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["bamfilezmw_samplenames"]] <- bamfilezmw_samplenames
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["genome"]] <- genome
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["referenceinfo"]] <- referenceinfo
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["chroms_to_analyze"]] <- chroms_to_analyze
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["thresholds"]] <- yaml.config$thresholds[[k]]
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["vcffilters"]] <- vcffilters[[k]]
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["bigwigfilters"]] <- bigwigfilters[[k]]
      bamfilezmw_filtered[[i]][[j]][[k]][["yaml.config"]][["bamfilters"]] <- bamfilters[[k]]
      
      bamfilezmw_filtered[[i]][[j]][[k]][["includezmws"]] <- includezmws
      bamfilezmw_filtered[[i]][[j]][[k]][["includemutations.fwd"]] <- includemutations.fwd
      bamfilezmw_filtered[[i]][[j]][[k]][["includemutations.rev"]] <- includemutations.rev
      bamfilezmw_filtered[[i]][[j]][[k]][["includemutations.zmw"]] <- includemutations.zmw
      bamfilezmw_filtered[[i]][[j]][[k]][["readregionalfilters.fwd"]] <- readregionalfilters.fwd
      bamfilezmw_filtered[[i]][[j]][[k]][["readregionalfilters.rev"]] <- readregionalfilters.rev
      bamfilezmw_filtered[[i]][[j]][[k]][["readregionalfilters.zmw"]] <- readregionalfilters.zmw
      bamfilezmw_filtered[[i]][[j]][[k]][["genomeregionalfilters.fwd"]] <- genomeregionalfilters.fwd
      bamfilezmw_filtered[[i]][[j]][[k]][["genomeregionalfilters.rev"]] <- genomeregionalfilters.rev
      bamfilezmw_filtered[[i]][[j]][[k]][["genomeregionalfilters.zmw"]] <- genomeregionalfilters.zmw

      bamfilezmw_filtered[[i]][[j]][[k]][["BAM_CRAM_ref_data"]][["includemutations.fwd.BAMTotalReads"]] <- includemutations.fwd.BAMTotalReads
      bamfilezmw_filtered[[i]][[j]][[k]][["BAM_CRAM_ref_data"]][["includemutations.fwd.BAMVariantReads"]] <- includemutations.fwd.BAMVariantReads
      bamfilezmw_filtered[[i]][[j]][[k]][["BAM_CRAM_ref_data"]][["includemutations.fwd.BAMVAF"]] <- includemutations.fwd.BAMVAF
      bamfilezmw_filtered[[i]][[j]][[k]][["BAM_CRAM_ref_data"]][["includemutations.rev.BAMTotalReads"]] <- includemutations.rev.BAMTotalReads
      bamfilezmw_filtered[[i]][[j]][[k]][["BAM_CRAM_ref_data"]][["includemutations.rev.BAMVariantReads"]] <- includemutations.rev.BAMVariantReads
      bamfilezmw_filtered[[i]][[j]][[k]][["BAM_CRAM_ref_data"]][["includemutations.rev.BAMVAF"]] <- includemutations.rev.BAMVAF
      bamfilezmw_filtered[[i]][[j]][[k]][["BAM_CRAM_ref_data"]][["includemutations.zmw.BAMTotalReads"]] <- includemutations.zmw.BAMTotalReads
      bamfilezmw_filtered[[i]][[j]][[k]][["BAM_CRAM_ref_data"]][["includemutations.zmw.BAMVariantReads"]] <- includemutations.zmw.BAMVariantReads
      bamfilezmw_filtered[[i]][[j]][[k]][["BAM_CRAM_ref_data"]][["includemutations.zmw.BAMVAF"]] <- includemutations.zmw.BAMVAF
      
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.fwd.ssDNAsubreadsVariantReads"]] <- includemutations.fwd.ssDNAsubreadsVariantReads
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.fwd.ssDNAsubreadsVAF"]] <- includemutations.fwd.ssDNAsubreadsVAF
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.fwd.ssDNAsubreadscvgfraction"]] <- includemutations.fwd.ssDNAsubreadscvgfraction
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.rev.ssDNAsubreadsVariantReads"]] <- includemutations.rev.ssDNAsubreadsVariantReads
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.rev.ssDNAsubreadsVAF"]] <- includemutations.rev.ssDNAsubreadsVAF
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.rev.ssDNAsubreadscvgfraction"]] <- includemutations.rev.ssDNAsubreadscvgfraction
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.zmw.fwdsubreadsVariantReads"]] <- includemutations.zmw.fwdsubreadsVariantReads
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.zmw.revsubreadsVariantReads"]] <- includemutations.zmw.revsubreadsVariantReads
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.zmw.fwdsubreadsVAF"]] <- includemutations.zmw.fwdsubreadsVAF
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.zmw.revsubreadsVAF"]] <- includemutations.zmw.revsubreadsVAF
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.zmw.fwdsubreadscvgfraction"]] <- includemutations.zmw.fwdsubreadscvgfraction
      bamfilezmw_filtered[[i]][[j]][[k]][["subread_data"]][["includemutations.zmw.revsubreadscvgfraction"]] <- includemutations.zmw.revsubreadscvgfraction
      
      bamfilezmw_filtered[[i]][[j]][[k]][["zmwstoretrieve"]] <- zmwstoretrieve
      
      #Delete temporary objects
      rm(includezmws,includemutations.fwd,includemutations.rev,includemutations.zmw,readregionalfilters.fwd,readregionalfilters.rev,readregionalfilters.zmw,genomeregionalfilters.fwd,genomeregionalfilters.rev,genomeregionalfilters.zmw,includemutations.fwd.BAMTotalReads,includemutations.fwd.BAMVariantReads,includemutations.fwd.BAMVAF,includemutations.rev.BAMTotalReads,includemutations.rev.BAMVariantReads,includemutations.rev.BAMVAF,includemutations.zmw.BAMTotalReads,includemutations.zmw.BAMVariantReads,includemutations.zmw.BAMVAF,includemutations.fwd.ssDNAsubreadsVariantReads,includemutations.fwd.ssDNAsubreadsVAF,includemutations.fwd.ssDNAsubreadscvgfraction,includemutations.rev.ssDNAsubreadsVariantReads,includemutations.rev.ssDNAsubreadsVAF,includemutations.rev.ssDNAsubreadscvgfraction,includemutations.zmw.fwdsubreadsVariantReads,includemutations.zmw.revsubreadsVariantReads,includemutations.zmw.fwdsubreadsVAF,includemutations.zmw.revsubreadsVAF,includemutations.zmw.fwdsubreadscvgfraction,includemutations.zmw.revsubreadscvgfraction,zmwstoretrieve)
    }
    
  }
}
#Close output stream
sink()

return(bamfilezmw_filtered)
}

######
# Run analysis
######

# Run analysis
bamfilezmw_filtered <- analyzebamfilezmw_all(bamfilezmw_all_filename,thresholds,vcffilters,bigwigfilters,bamfilters,referenceinfo,analysisoutput_path,subreads_filename,bigwig_cachedir,chroms_to_analyze)

# Save analysis
cat("\n#### Saving analysis RDS file...")
qsave(bamfilezmw_filtered,paste0(analysisoutput_path,"/",analysisoutput_basename,".RDS"))
cat("DONE\n")
