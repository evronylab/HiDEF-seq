#!/usr/bin/env Rscript
# Usage: print_mutations.R [configuration.yaml]
#  - Uses [configuration.yaml] file per documentation
#  - Outputs to analysisoutput_path directory defined in configuration.yaml:
#     print_mutations.output.txt recording pipeline output (also saved as print_mutations.output.RDS for use by mutation_frequencies script). Non self-explanatory columns:
#				strandtype: ssDNA or dsDNA, signifies ssDNA or dsDNA mutation
#				synthesizedstrand (applies only to ssDNA mutations): the strand of the reference genome that the read aligned to. This corresponds to the strand synthesized during sequencing.
#				templatestrand (applies only to ssDNA mutations): the strand of the template being sequenced. This is the reverse of synthesizedstrand
#				ref_refplusstrand: reference base sequence relative to the reference plus strand
#				alt_refplusstrand: mutation base sequence relative to the reference plus strand
#				tnc_refplusstrand: trinucleotide context relative to the reference plus strand
#				ref32_refplusstrand (applies only to dsDNA mutations): reference base sequence relative to the reference plus strand, then transformed to strand where reference base is a pyrimidine
#				alt32_refplusstrand (applies only to dsDNA mutations): mutation base sequence on the reference plus strand, then transformed to strand where reference base is a pyrimidine
#				tnc32_refplusstrand (applies only to dsDNA mutations): trinucleotide context relative to the reference plus strand, then transformed to strand where reference base is a pyrimidine
#				ref_ssDNAsynthesizedstrand (applies only to ssDNA mutations): reference base sequence relative to the synthesized strand
#				alt_ssDNAsynthesizedstrand (applies only to ssDNA mutations): mutation base sequence relative to the synthesized strand
#				tnc_ssDNAsynthesizedstrand (applies only to ssDNA mutations): trinucleotide context relative to the synthesized strand
#				ref_ssDNAtemplatestrand (applies only to ssDNA mutations): reference base sequence relative to the template strand
#				alt_ssDNAtemplatestrand (applies only to ssDNA mutations): mutation base sequence relative to the template strand
#				tnc_ssDNAtemplatestrand (applies only to ssDNA mutations): trinucleotide context relative to the template strand
#				
# Prerequisites:
#  - R packages: configr, qs, Biostrings, BSgenome, BSgenome package for the reference genome specified in the configuration (e.g., BSgenome.Hsapiens.T2T.CHM13v1)

#Stop script for any warnings
options(warn=2)

suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(Biostrings))

cat("\n####### Print List of Mutations Script #######\n\n")

#Function to reduce 64 to 32 trinucleotide frequency with central pyrimidine. Input is vector of trinucleotide strings.
trinucleotide64to32vector <- function(x){
  suppressPackageStartupMessages(library(Biostrings))
  trinucleotides_64 <- apply(expand.grid(c("A","C","G","T"),c("A","C","G","T"),c("A","C","G","T")),1,paste,collapse="")
  trinucleotides_32_pyr <- apply(expand.grid(c("A","C","G","T"),c("C","T"),c("A","C","G","T")),1,paste,collapse="")
  trinucleotides_32_pur <- setdiff(trinucleotides_64,trinucleotides_32_pyr)
  
  x[x %in% trinucleotides_32_pur] <- as.character(reverseComplement(DNAStringSet(x[x %in% trinucleotides_32_pur])))
  return(x)
}

##Load configuration
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=1) {
  stop("Usage: Rscript print_mutations.R [configuration.yaml]\n", call.=FALSE)
}
cat("#### Loading configuration...")
yaml.config <- suppressWarnings(read.config(args[1]))

#Load BSgenome package
suppressPackageStartupMessages(library(yaml.config$BSgenomepackagename,character.only=TRUE))

cat("DONE\n")

##Load analysis files
cat("#### Loading analysis files...")
bamfilezmw_all <- qread(yaml.config$bamfilezmw_all_filename)
bamfilezmw_filtered <- qread(paste0(yaml.config$analysisoutput_path,"/",yaml.config$analysisoutput_basename,".RDS"))
cat("DONE\n")

cat("#### Printing mutations to output file...")

output.df <- c()

for(sampleid in names(bamfilezmw_filtered)){
  for(genome in names(bamfilezmw_filtered[[sampleid]])){
    for(filterid in 1:length(bamfilezmw_filtered[[sampleid]][[genome]])){
    	
    	if(length(bamfilezmw_filtered[[sampleid]][[genome]][[filterid]])==0){
    		cat(paste0("\n    Note: No ZMWs after filtering ",sampleid,".",genome,".",filterid,"\n\n"))
    		next
    	}
      
      includezmws <- bamfilezmw_filtered[[sampleid]][[genome]][[filterid]]$includezmws
      
      totalfwdrevzmwsnvs <- 0
      
      for(i in c("fwd","rev","zmw")){
        numsnvs <- unlist(lapply(bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]],sum),use.names=FALSE)

        if(i %in% c("fwd","rev")){
        	zmwstoexamine <- numsnvs > 0 & numsnvs <= yaml.config$mutratefilters$maxmutationsperssdna
        }else if(i=="zmw"){
        	zmwstoexamine <- numsnvs > 0 & numsnvs <= yaml.config$mutratefilters$maxmutationsperzmw
        }
        numsnvs <- numsnvs[zmwstoexamine]
        totalsnvs <- sum(numsnvs)
        totalfwdrevzmwsnvs <- totalfwdrevzmwsnvs + totalsnvs
        
        if(totalsnvs==0){
        	next
        } #Skip this data output code iteration if no mutations found
        
        numsnvs.output <- unlist(mapply(function(x,y){rep(x,y)},x=numsnvs,y=numsnvs),use.names=FALSE)
        
        sampleid.output <- rep(sampleid,totalsnvs)
        genome.output <- rep(genome,totalsnvs)
        filterid.output <- rep(filterid,totalsnvs)
        
        zmw.rqfwd.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$tag.fwd$rq[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
        zmw.rqrev.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$tag.rev$rq[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
        zmw.ecfwd.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$tag.fwd$ec[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
        zmw.ecrev.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$tag.rev$ec[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
        zmw.mapqfwd.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$mapq.fwd[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
        zmw.mapqrev.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$mapq.rev[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
        
        GermlineVariantReads.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["BAM_CRAM_ref_data"]][[paste0("includemutations.",i,".BAMVariantReads")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
        GermlineTotalReads.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["BAM_CRAM_ref_data"]][[paste0("includemutations.",i,".BAMTotalReads")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
        GermlineVAF.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["BAM_CRAM_ref_data"]][[paste0("includemutations.",i,".BAMVAF")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
        
        if(i=="fwd" | i=="rev"){
          strandtype.output <- rep("ssDNA",totalsnvs)
          
          zmw.name.output  <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]][[paste0("zmwname.",i)]][includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
          zmw.hole.output  <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]][[paste0("tag.",i)]]$zm[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
          zmw.isize.output <- rep(NA,totalsnvs)
          
          chrom.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]][[paste0("rname.",i)]][includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
          pos.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_all[[sampleid]][[genome]][[paste0("tag.",i)]]$rp[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          synthesizedstrand.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]][[paste0("strand.",i)]][includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
          templatestrand.output <- rep(NA,totalsnvs)
          templatestrand.output[synthesizedstrand.output=="+"] <- "-"
          templatestrand.output[synthesizedstrand.output=="-"] <- "+"
            
          ref_refplusstrand.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_all[[sampleid]][[genome]][[paste0("tag.",i)]]$rn[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          alt_refplusstrand.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_all[[sampleid]][[genome]][[paste0("tag.",i)]]$qn[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          
          tnc_refplusstrand.output <- as.character(getSeq(eval(parse(text=yaml.config$BSgenomepackagename)),resize(makeGRangesFromDataFrame(data.frame(chrom=chrom.output,start=pos.output,end=pos.output)),3,fix="center")))
          
          ref32_refplusstrand.output <- rep(NA,totalsnvs)
          alt32_refplusstrand.output <- rep(NA,totalsnvs)
          tnc32_refplusstrand.output <- rep(NA,totalsnvs)
          
          bpenddist.output <- unlist(mapply(function(x,y,a,b){pmin(x[y]-a+1,(a+b)-x[y])},x=bamfilezmw_all[[sampleid]][[genome]][[paste0("tag.",i)]]$rp[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],a=bamfilezmw_all[[sampleid]][[genome]][[paste0("pos.",i)]][includezmws][zmwstoexamine],b=bamfilezmw_all[[sampleid]][[genome]][[paste0("isize.",i)]][includezmws][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          
          if(i=="fwd"){
            isizefwd.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$isize.fwd[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
            isizerev.output <- rep(NA,totalsnvs)
            zmw.isize.output <- rep(NA,totalsnvs)
            
            ref_ssDNAsynthesizedstrand.output <- ref_refplusstrand.output
            alt_ssDNAsynthesizedstrand.output <- alt_refplusstrand.output
            ref_ssDNAtemplatestrand.output <- as.character(reverseComplement(DNAStringSet(ref_ssDNAsynthesizedstrand.output)))
            alt_ssDNAtemplatestrand.output <- as.character(reverseComplement(DNAStringSet(alt_ssDNAsynthesizedstrand.output)))
            
            tnc_ssDNAsynthesizedstrand.output <- tnc_refplusstrand.output
            tnc_ssDNAtemplatestrand.output <- as.character(reverseComplement(DNAStringSet(tnc_ssDNAsynthesizedstrand.output)))
            
            qqfwd.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_all[[sampleid]][[genome]]$tag.fwd$qq[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
            qqrev.output <- rep(NA,totalsnvs)
        
            fwdsubreadsVariantReads.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".ssDNAsubreadsVariantReads")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
            revsubreadsVariantReads.output <- rep(NA,totalsnvs)
            fwdsubreadsVAF.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".ssDNAsubreadsVAF")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
            revsubreadsVAF.output <- rep(NA,totalsnvs)
            fwdsubreadscvgfraction.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".ssDNAsubreadscvgfraction")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
            revsubreadscvgfraction.output <- rep(NA,totalsnvs)
                      
          }else if(i=="rev"){
            isizefwd.output <- rep(NA,totalsnvs)
            isizerev.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$isize.rev[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
            zmw.isize.output <- rep(NA,totalsnvs)

            ref_ssDNAsynthesizedstrand.output <- as.character(reverseComplement(DNAStringSet(ref_refplusstrand.output)))
            alt_ssDNAsynthesizedstrand.output <- as.character(reverseComplement(DNAStringSet(alt_refplusstrand.output)))
            ref_ssDNAtemplatestrand.output <- as.character(reverseComplement(DNAStringSet(ref_ssDNAsynthesizedstrand.output)))
            alt_ssDNAtemplatestrand.output <- as.character(reverseComplement(DNAStringSet(alt_ssDNAsynthesizedstrand.output)))
            
            tnc_ssDNAsynthesizedstrand.output <- as.character(reverseComplement(DNAStringSet(tnc_refplusstrand.output)))
            tnc_ssDNAtemplatestrand.output <- as.character(reverseComplement(DNAStringSet(tnc_ssDNAsynthesizedstrand.output)))
            
            qqfwd.output <- rep(NA,totalsnvs)
            qqrev.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_all[[sampleid]][[genome]]$tag.rev$qq[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
            
            revsubreadsVariantReads.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".ssDNAsubreadsVariantReads")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
            fwdsubreadsVariantReads.output <- rep(NA,totalsnvs)
            revsubreadsVAF.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".ssDNAsubreadsVAF")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
            fwdsubreadsVAF.output <- rep(NA,totalsnvs)
            revsubreadscvgfraction.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".ssDNAsubreadscvgfraction")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
            fwdsubreadscvgfraction.output <- rep(NA,totalsnvs)
          }
          
        }else if(i=="zmw"){
          strandtype.output <- rep("dsDNA",totalsnvs)
          
          zmw.name.output  <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$zmwname.fwd[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
          zmw.hole.output  <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$tag.fwd$zm[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
          
          isizefwd.output <- rep(NA,totalsnvs)
          isizerev.output <- rep(NA,totalsnvs)
          zmw.isize.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$zmw.isize[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
          
          chrom.output <- unlist(mapply(function(x,y){rep(x,y)},x=bamfilezmw_all[[sampleid]][[genome]]$rname.fwd[includezmws][zmwstoexamine],y=numsnvs),use.names=FALSE)
          pos.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_all[[sampleid]][[genome]]$zmw.rp[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]]$includemutations.zmw[zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          synthesizedstrand.output <- rep(NA,totalsnvs)
          templatestrand.output <- rep(NA,totalsnvs)
          
          ref_refplusstrand.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_all[[sampleid]][[genome]]$zmw.rn[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]]$includemutations.zmw[zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          alt_refplusstrand.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_all[[sampleid]][[genome]]$zmw.qn[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]]$includemutations.zmw[zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          ref_ssDNAsynthesizedstrand.output <- rep(NA,totalsnvs)
          alt_ssDNAsynthesizedstrand.output <- rep(NA,totalsnvs)
          ref_ssDNAtemplatestrand.output <- rep(NA,totalsnvs)
          alt_ssDNAtemplatestrand.output <- rep(NA,totalsnvs)

          tnc_refplusstrand.output <- as.character(getSeq(eval(parse(text=yaml.config$BSgenomepackagename)),resize(makeGRangesFromDataFrame(data.frame(chrom=chrom.output,start=pos.output,end=pos.output)),3,fix="center")))
          tnc_ssDNAsynthesizedstrand.output <- rep(NA,totalsnvs)
          tnc_ssDNAtemplatestrand.output <- rep(NA,totalsnvs)
          
          #ref, alt, and tnc reduced to 32 trinucleotide contexts
          ref32_refplusstrand.output <- ref_refplusstrand.output
          ref32_refplusstrand.output[ref_refplusstrand.output %in% c("A","G")] <- as.character(reverseComplement(DNAStringSet(ref32_refplusstrand.output[ref_refplusstrand.output %in% c("A","G")])))
          alt32_refplusstrand.output <- alt_refplusstrand.output
          alt32_refplusstrand.output[ref_refplusstrand.output %in% c("A","G")] <- as.character(reverseComplement(DNAStringSet(alt32_refplusstrand.output[ref_refplusstrand.output %in% c("A","G")])))
          tnc32_refplusstrand.output <- trinucleotide64to32vector(tnc_refplusstrand.output)
                    
          bpenddist.output <- unlist(mapply(function(x,y,a,b,c){pmin(x[y]-max(a,b)+1,max(a,b)+c-x[y])},x=bamfilezmw_all[[sampleid]][[genome]]$zmw.rp[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]]$includemutations.zmw[zmwstoexamine],a=bamfilezmw_all[[sampleid]][[genome]]$pos.fwd[includezmws][zmwstoexamine],b=bamfilezmw_all[[sampleid]][[genome]]$pos.rev[includezmws][zmwstoexamine],c=bamfilezmw_all[[sampleid]][[genome]]$zmw.isize[includezmws][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          
          qqfwd.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_all[[sampleid]][[genome]]$zmw.fwdqq[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          qqrev.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_all[[sampleid]][[genome]]$zmw.revqq[includezmws][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          
          fwdsubreadsVariantReads.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".fwdsubreadsVariantReads")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          revsubreadsVariantReads.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".revsubreadsVariantReads")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          fwdsubreadsVAF.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".fwdsubreadsVAF")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          revsubreadsVAF.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".revsubreadsVAF")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          fwdsubreadscvgfraction.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".fwdsubreadscvgfraction")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
          revsubreadscvgfraction.output <- unlist(mapply(function(x,y){x[y]},x=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][["subread_data"]][[paste0("includemutations.",i,".revsubreadscvgfraction")]][zmwstoexamine],y=bamfilezmw_filtered[[sampleid]][[genome]][[filterid]][[paste0("includemutations.",i)]][zmwstoexamine],SIMPLIFY=FALSE),use.names=FALSE)
        }
        
        temp.df <- data.frame(sampleid=sampleid.output,genome=genome.output,filterid=filterid.output,zmw.name=zmw.name.output,zmw.hole=zmw.hole.output,isizefwd=isizefwd.output,isizerev=isizerev.output,zmw.isize=zmw.isize.output,chrom=chrom.output,pos=pos.output,strandtype=strandtype.output,synthesizedstrand=synthesizedstrand.output,templatestrand=templatestrand.output,ref_refplusstrand=ref_refplusstrand.output,alt_refplusstrand=alt_refplusstrand.output,tnc_refplusstrand=tnc_refplusstrand.output,ref32_refplusstrand=ref32_refplusstrand.output,alt32_refplusstrand=alt32_refplusstrand.output,tnc32_refplusstrand=tnc32_refplusstrand.output,ref_ssDNAsynthesizedstrand=ref_ssDNAsynthesizedstrand.output,alt_ssDNAsynthesizedstrand=alt_ssDNAsynthesizedstrand.output,tnc_ssDNAsynthesizedstrand=tnc_ssDNAsynthesizedstrand.output,ref_ssDNAtemplatestrand=ref_ssDNAtemplatestrand.output,alt_ssDNAtemplatestrand=alt_ssDNAtemplatestrand.output,tnc_ssDNAtemplatestrand=tnc_ssDNAtemplatestrand.output,zmw.rqfwd=zmw.rqfwd.output,zmw.rqrev=zmw.rqrev.output,zmw.ecfwd=zmw.ecfwd.output,zmw.ecrev=zmw.ecrev.output,zmw.mapqfwd=zmw.mapqfwd.output,zmw.mapqrev=zmw.mapqrev.output,qqfwd=qqfwd.output,qqrev=qqrev.output,bpenddist=bpenddist.output,mutationsperstrandorzmw=numsnvs.output,GermlineVariantReads=GermlineVariantReads.output,GermlineTotalReads=GermlineTotalReads.output,GermlineVAF=GermlineVAF.output,fwdsubreadsVariantReads=fwdsubreadsVariantReads.output,fwdsubreadsVAF=fwdsubreadsVAF.output,fwdsubreadscvgfraction=fwdsubreadscvgfraction.output,revsubreadsVariantReads=revsubreadsVariantReads.output,revsubreadsVAF=revsubreadsVAF.output,revsubreadscvgfraction=revsubreadscvgfraction.output,row.names=NULL)
        output.df <- rbind(output.df,temp.df)
      }
      
      if(totalfwdrevzmwsnvs==0){
      	cat(paste0("\n    Note: No mutations after filtering ",sampleid,".",genome,".",filterid,"\n\n"))
      }
      
    }
  }
}

write.table(output.df,file=paste0(yaml.config$analysisoutput_path,"/print_mutations.output.txt"),quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
qsave(output.df,file=paste0(yaml.config$analysisoutput_path,"/print_mutations.output.RDS"))
cat("DONE\n")

if(is.null(output.df)){
	cat("  No mutations identified in any sample!\n")
}
