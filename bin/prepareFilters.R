#RUN THIS BEFORE extractVariants

#Install BSgenome reference in cache dir if needed
.libPaths(yaml.config$prepareFilters_cache_dir)
if(!yaml.config$BSgenome$BSgenome_name %in% installed.genomes()){
	if(yaml.config$BSgenome$BSgenome_name %in% available.genomes()){
		BiocManager::install(yaml.config$BSgenome$BSgenome_name,lib=yaml.config$prepareFilters_cache_dir)
	}else if(!is.null(yaml.config$BSgenome$BSgenome_file)){
		dir.create(yaml.config$prepareFilters_cache_dir, recursive = TRUE, showWarnings = FALSE)
		install.packages(yaml.config$BSgenome$BSgenome_file, repos = NULL,  lib = yaml.config$prepareFilters_cache_dir)
	}else{
		stop("ERROR: Must specify either BSgenome_name that is in available.genomes() or a BSgenome_file!", call.=FALSE)
	}
}

**Load all variants, not just analyzed chromosomes, and note that this will be used for all chunks

######################
### Load VCF variants
######################
cat("#### Loading VCF variants...")

#Get list of VCF files for the individual corresponding to the sample
individual_id <- yaml.config$samples[sample_id]$individual_id
germline_vcf_files <- yaml.config$individuals[individual_id]$germline_vcf_files

#Annotate with each VCF file (make a list of tables, each table with the same number of rows as sbses.df)
vcf_indels <- list()
vcf_snvs <- list()

for(i in germline_vcf_files){
  
  cat("     Loading VCF File:",i$germline_vcf_file,"...")
  
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
