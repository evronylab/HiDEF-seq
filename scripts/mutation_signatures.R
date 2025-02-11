#!/usr/bin/env Rscript
# Usage: mutation_signatures.R [configuration.yaml]
#  - Create [configuration.yaml] per documentation
#
# Prerequisites
#   - R packages: 
#		- mutation_frequencies script must have been previously run
#
# Outputs:
#		- Mutation signature plots
#		- Observed and corrected mutation count spectra
#		- mutation_signatures.RDS with sigfit results
#				
cat("####### Mutation Signatures Analysis #######\n\n")

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(sigfit))
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

######
#Load configuration
######
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=1) {
  stop("Usage: mutation_signatures.R [configuration.yaml]\n", call.=FALSE)
}
yaml.config <- suppressWarnings(read.config(args[1]))

#Load genomic coordinates and subtract N regions, to calculte total genome size
chrs_to_analyze <- read.table(yaml.config$chrs_to_analyze)[,1]

genome.granges <- as.data.frame(read.table(yaml.config$chrsizes))
genome.granges$V3 <- rep(1,nrow(genome.granges))
genome.granges <- makeGRangesFromDataFrame(genome.granges,seqnames.field="V1",start.field="V3",end.field="V2")
genome.granges <- genome.granges[seqnames(genome.granges) %in% chrs_to_analyze]

Nref.granges <- import(yaml.config$Nrefbed,format="BED")
Nref.granges <- Nref.granges[seqnames(Nref.granges) %in% chrs_to_analyze]
genome.noN.granges <- GenomicRanges::setdiff(genome.granges,Nref.granges)

genome.size <- sum(width(genome.noN.granges))

######
#Load mutation_frequencies RDS file created by mutation_frequencies script
######
mutation_frequencies <- qread(paste0(yaml.config$analysisoutput_path,"/mutation_frequencies.RDS"))
if(class(mutation_frequencies)!="list"){
	cat("No reads or mutations in selected chromosomes!\n")
	quit(save="no")
}

samplenames <- names(mutation_frequencies$mutations)
numsamples <- length(samplenames)

mutation_signatures <- list() #Output for results to save to mutation_signatures.RDS


######
#Sigfit: Load general data required for analysis
######
##Load COSMIC v3.2 signatures. Note: we use this human-genome-dependent version of COSMIC for signature fitting, rather than projecting these signatures to our observed trinucleotide frequencies. This is because we are inputting mutation opportunities into signature fitting that are themselves already normalized to the human genome mutation opportunities. Therefore the mutation opportunities used during signature fitting adjust only for differences between our observed trinucleotide frequencies and the trinucleotide frequencies of the human genome, such that we can use the standard human genome COSMIC signatures during fitting.
data("cosmic_signatures_v3.2")

#Order of labels for creating mutation_opportunities matrices
genome_freqs_labels <- c("ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT")

#Enable sigfit parallelization
options(mc.cores = parallel::detectCores())

######
#Sigfit: dsDNA mutational signatures analysis and plotting
######

cat("#### dsDNA mutation signature analysis...\n")

if(nrow(mutation_frequencies$SigFitcollapsed$dsDNA)==0){
	cat("       No dsDNA mutations in any samples.\n\n")
}else{
	##Make SigFit dsDNA mutation catalogue
	SigFit_dsDNA_catalogue <- build_catalogues(mutation_frequencies$SigFitcollapsed$dsDNA)
	 #Fill in samples with 0 mutation counts
	for(i in setdiff(samplenames,rownames(SigFit_dsDNA_catalogue))){
		SigFit_dsDNA_catalogue <- rbind(SigFit_dsDNA_catalogue,rep(0,ncol(SigFit_dsDNA_catalogue)))
		rownames(SigFit_dsDNA_catalogue)[nrow(SigFit_dsDNA_catalogue)] <- i
	}
	SigFit_dsDNA_catalogue <- SigFit_dsDNA_catalogue[samplenames,]
	
	##Create a matrix of mutational opportunities for each sample. This is calculated for each sample as follows:
	#  (1) Divide the trinucleotide frequencies observed in the interrogated read bases, by the sum of the trinucleotide frequencies observed in the interrogated read bases.
	#  (2) Divide the trinucleotide frequencies of the full human genome (the human genome version used by COSMIC whose trinucleotide frequencies are stored in sigfit; retrievable with human_trinuc_freqs(), by the sum of the trinucleotide frequencies of the full human genome.
	#  (3) Divide (1) by (2).
	# This produces a mutational opportunity matrix reflecting the relative over/under opportunity of our method beyond the baseline trinucleotide distribution of the human genome. For example, if a specific trinucleotide was under-observed in our reads relative to the human genome, it will have a mutation opportunity < 1, and then during signature fitting, sigfit will increase the mutation counts for that trinucleotide.
	#Note: Providing the interrogated read to genome trinucleotide frequency ratio as the mutation opportunities to sigfit produces a more accurate analysis, because signature fitting error models depend on the total number of mutations. So simply scaling the number of mutations by the interrogated read trinucleotide frequencies divided by the genome trinucleotide frequencies would change the total number of mutations and skew signature fitting results.
	mutation_opportunities <- matrix(nrow=numsamples,ncol=length(genome_freqs_labels),dimnames=list(samplenames,genome_freqs_labels))
	
	for(i in samplenames){
		mutation_opportunities[i,] <- (mutation_frequencies$trinucleotide_contexts[[i]]$dsDNA$mutationreadsfreq[genome_freqs_labels]/sum(mutation_frequencies$trinucleotide_contexts[[i]]$dsDNA$mutationreadsfreq[genome_freqs_labels]))/(human_trinuc_freqs()/sum(human_trinuc_freqs()))
	}
	
	##Plot mutation spectra with y-axis of either counts or probability, and for probability y-axis both uncorrected and corrected for mutation opportunities (3 total PDFs)
	plot_spectrum(SigFit_dsDNA_catalogue,pdf_path=paste0(yaml.config$analysisoutput_path,"/sigfit.dsDNA.counts.uncorrected.pdf"))
	
	dsDNA.probability.uncorrected <- SigFit_dsDNA_catalogue/rowSums(SigFit_dsDNA_catalogue)
		#Change samples with no mutations to all 0s instead of NA
	for(i in 1:nrow(dsDNA.probability.uncorrected)){
		if(all(is.na(dsDNA.probability.uncorrected[i,]))){
			dsDNA.probability.uncorrected[i,] <- rep(0,length(dsDNA.probability.uncorrected[i,]))
		}
	}
	plot_spectrum(dsDNA.probability.uncorrected,pdf_path=paste0(yaml.config$analysisoutput_path,"/sigfit.dsDNA.probability.uncorrected.pdf"))
	
	dsDNA.probability.corrected <- (SigFit_dsDNA_catalogue[samplenames,]/mutation_opportunities[samplenames,])/rowSums(SigFit_dsDNA_catalogue[samplenames,]/mutation_opportunities[samplenames,])
		#Change samples with no mutations to all 0s instead of NA
	for(i in 1:nrow(dsDNA.probability.corrected)){
		if(all(is.na(dsDNA.probability.corrected[i,]))){
			dsDNA.probability.corrected[i,] <- rep(0,length(dsDNA.probability.corrected[i,]))
		}
	}
	plot_spectrum(dsDNA.probability.corrected,paste0(yaml.config$analysisoutput_path,"/sigfit.dsDNA.probability.corrected.pdf"))
	
	#Output observed and corrected counts into tables
	write.table(SigFit_dsDNA_catalogue,paste0(yaml.config$analysisoutput_path,"/dsDNA.trinuc.counts.uncorrected.txt"),quote=FALSE,sep="\t",col.names=NA)
	write.table(SigFit_dsDNA_catalogue[samplenames,]/mutation_opportunities[samplenames,],paste0(yaml.config$analysisoutput_path,"/dsDNA.trinuc.counts.corrected.txt"),quote=FALSE,sep="\t",col.names=NA)
	
	##Fit signatures to known COSMIC signatures, and plot spectra, exposures, and reconstruction
	SigFit_dsDNA_fit_signatures <- list()
	for(i in samplenames){
		SigFit_dsDNA_fit_signatures[[i]] <- fit_signatures(counts=SigFit_dsDNA_catalogue[i,,drop=F],
																						signatures=cosmic_signatures_v3.2,
																						model="multinomial",
																						opportunities=mutation_opportunities[i,,drop=F],
																						iter=10000,warmup=5000,
																						chains=4,seed = 1756,control=list(adapt_delta = 0.99),
																						show_messages=FALSE, open_progress=FALSE)
	}
	mutation_signatures$SigFit_dsDNA_fit_signatures <- SigFit_dsDNA_fit_signatures
	
	SigFit_dsDNA_exposures <- list()
	for(i in samplenames){
		SigFit_dsDNA_exposures[[i]] <- retrieve_pars(SigFit_dsDNA_fit_signatures[[i]],par="exposures")
	}
	
	for(i in samplenames){
		plot_all(SigFit_dsDNA_fit_signatures[[i]],out_path=yaml.config$analysisoutput_path,prefix=paste0(i,".sigfit.dsDNA"))
	}
	
	##Redo fit_signatures only with exposures whose lower end of the Bayesian HPD interval is above a specified threshold, plus SBS1, SBS5, SBS18, and SBS40 that are commonly seen in many healthy tissues, and plot spectra, exposures, and reconstruction.
	lowerhpd_thresh <- 0.01
	healthy_signatures <- c("SBS1","SBS5","SBS18","SBS40")
	
	SigFit_dsDNA_fit_signatures_thresh <- list()
	for(i in samplenames){
		selectedsignatures <- str_sort(unique(c(rownames(cosmic_signatures_v3.2)[SigFit_dsDNA_exposures[[i]]$lower_95>lowerhpd_thresh],healthy_signatures)),numeric=TRUE)
		SigFit_dsDNA_fit_signatures_thresh[[i]] <- fit_signatures(counts=SigFit_dsDNA_catalogue[i,,drop=F],
																						signatures=cosmic_signatures_v3.2[selectedsignatures,,drop=F],
																						model="multinomial",
																						opportunities=mutation_opportunities[i,,drop=F],
																						iter=10000,chains = 4,seed = 1756,control = list(adapt_delta = 0.99),
																						show_messages=FALSE, open_progress=FALSE)
	}
	mutation_signatures$SigFit_dsDNA_fit_signatures_thresh <- SigFit_dsDNA_fit_signatures_thresh
	
	SigFit_dsDNA_exposures_thresh <- list()
	for(i in samplenames){
		SigFit_dsDNA_exposures_thresh[[i]] <- retrieve_pars(SigFit_dsDNA_fit_signatures_thresh[[i]],par="exposures")
	}
	
	for(i in samplenames){
		plot_all(SigFit_dsDNA_fit_signatures_thresh[[i]],out_path=yaml.config$analysisoutput_path,prefix=paste0(i,".sigfit_thresh.dsDNA"))
	}

}


######
#Sigfit: ssDNA UNCOLLAPSED mutational signatures plotting
######
#This uses Sigfit's capability of analyzing transcribed vs untranscribed strand mutations to separate central pyrimidine mutation context (which we arbitrarily annotate as transcribed strand) vs central purine mutation context (which we arbitrarily annotate as untranscribed strand). Important to note that transcribed and untranscribed do NOT actually correspond to transcriptional strand, but simply central pyrimidine vs central purine mutation contexts, respectively.

cat("#### ssDNA uncollapsed mutation signature analysis...\n")

if(nrow(mutation_frequencies$SigFit$ssDNA)==0){
	cat("       No ssDNA uncollapsed mutations in any samples.\n\n")
}else{
	
	##Make SigFit ssDNA mutation catalogue
	SigFit_ssDNA_catalogue <- mutation_frequencies$SigFit$ssDNA
	SigFit_ssDNA_catalogue$strand <- -1 #This assigns all mutations to the transcribed strand (= -1), but then Sigfit reverse complements and assigns to the untranscribed strand any mutations with a central purine. As a result, all central pyrimidine mutations are assigned to the transcribed strand, and all central purine mutations are assigned to the untranscribed strand, as we intend.
	
	SigFit_ssDNA_catalogue <- build_catalogues(SigFit_ssDNA_catalogue)
	 #Fill in samples with 0 mutation counts
	for(i in setdiff(samplenames,rownames(SigFit_ssDNA_catalogue))){
		SigFit_ssDNA_catalogue <- rbind(SigFit_ssDNA_catalogue,rep(0,ncol(SigFit_ssDNA_catalogue)))
		rownames(SigFit_ssDNA_catalogue)[nrow(SigFit_ssDNA_catalogue)] <- i
	}
	SigFit_ssDNA_catalogue <- SigFit_ssDNA_catalogue[samplenames,]

	##Create matrix of mutational opportunities for each sample. Since this is a ssDNA analysis separated by central pyrimidine vs central purine, this is created as follows:
	#	1) Separate trinucleotide frequencies for central pyrimidine from those for central purine, and reverse complement the labels of the central purine trinucleotide frequencies.
	#	2) Retrieve 96 central pyrimidine trinucleotide contexts in correct order from each of the 32 trinucleotide context central pyrimidine and central purine trinucleotide frequencies.
	# 3) Concatenate the 96 central pyrimidine and 96 central purine trinucleotide frequencies into a 192 trinucleotide frequencies vector, so that central pyrimidine and central purine frequencies are in the positions corresponding to transcribed and untranscribed strand mutations, respectively.
	# 4) Divide the 192 trinucleotide frequencies by their sum
	# 5) Divide the human reference genome trinucleotide frequencies (repeated twice in tandem, since each genome context is the frequency both for the central pyrimidine and central purine contexts; i.e. they are the same) by their sum (i.e. sum of the repeated twice in tandem frequencies).
	# 6) Divide (4) by (5).
	# This produces a mutational opportunity matrix reflecting the residual over/under opportunity of our method beyond the baseline trinucleotide distribution of the human genome. For example, if a specific trinucleotide was under-observed in our reads relative to the human genome, it will have a mutation opportunity < 1, and then during signature fitting, sigfit will increase the mutation counts for that trinucleotide.
	mutation_opportunities <- matrix(nrow=numsamples,ncol=length(genome_freqs_labels)*2,dimnames=list(samplenames,c(paste0("T:",genome_freqs_labels),paste0("U:",genome_freqs_labels))))
	
	for(i in samplenames){
		mutationreads_opportunities <- mutation_frequencies$trinucleotide_contexts[[i]]$fwdrev$mutationreadsfreq

		mutationreads_opportunities_pyr <- mutationreads_opportunities[substr(names(mutationreads_opportunities),2,2) %in% c("C","T")]
		mutationreads_opportunities_pur <- mutationreads_opportunities[substr(names(mutationreads_opportunities),2,2) %in% c("A","G")]
		names(mutationreads_opportunities_pur) <- reverseComplement(DNAStringSet(names(mutationreads_opportunities_pur)))
		
		mutationreads_opportunities_pyr <- mutationreads_opportunities_pyr[genome_freqs_labels]
		mutationreads_opportunities_pur <- mutationreads_opportunities_pur[genome_freqs_labels]
		
		mutation_opportunities[i,] <- c(mutationreads_opportunities_pyr,mutationreads_opportunities_pur)
		
		mutation_opportunities[i,] <- (mutation_opportunities[i,]/sum(mutation_opportunities[i,]))/(rep(human_trinuc_freqs(),2)/sum(rep(human_trinuc_freqs(),2)))

	}
	
	##Plot mutation spectra with y-axis of either counts or probability, and for probability y-axis both uncorrected and corrected for mutation opportunities (3 total PDFs)
	plot_spectrum(SigFit_ssDNA_catalogue,pdf_path=paste0(yaml.config$analysisoutput_path,"/sigfit.ssDNA.counts.uncorrected.pdf"))
	
	ssDNA.probability.uncorrected <- SigFit_ssDNA_catalogue/rowSums(SigFit_ssDNA_catalogue)
		#Change samples with no mutations to all 0s instead of NA
	for(i in 1:nrow(ssDNA.probability.uncorrected)){
		if(all(is.na(ssDNA.probability.uncorrected[i,]))){
			ssDNA.probability.uncorrected[i,] <- rep(0,length(ssDNA.probability.uncorrected[i,]))
		}
	}
	plot_spectrum(ssDNA.probability.uncorrected,pdf_path=paste0(yaml.config$analysisoutput_path,"/sigfit.ssDNA.probability.uncorrected.pdf"))
	
	ssDNA.probability.corrected <- (SigFit_ssDNA_catalogue[samplenames,]/mutation_opportunities[samplenames,])/rowSums(SigFit_ssDNA_catalogue[samplenames,]/mutation_opportunities[samplenames,])
		#Change samples with no mutations to all 0s instead of NA
	for(i in 1:nrow(ssDNA.probability.corrected)){
		if(all(is.na(ssDNA.probability.corrected[i,]))){
			ssDNA.probability.corrected[i,] <- rep(0,length(ssDNA.probability.corrected[i,]))
		}
	}
	plot_spectrum(ssDNA.probability.corrected,paste0(yaml.config$analysisoutput_path,"/sigfit.ssDNA.probability.corrected.pdf"))
	
	#Output observed and corrected counts into tables
	write.table(SigFit_ssDNA_catalogue,paste0(yaml.config$analysisoutput_path,"/ssDNA.trinuc.counts.uncorrected.txt"),quote=FALSE,sep="\t",col.names=NA)
	write.table(SigFit_ssDNA_catalogue[samplenames,]/mutation_opportunities[samplenames,],paste0(yaml.config$analysisoutput_path,"/ssDNA.trinuc.counts.corrected.txt"),quote=FALSE,sep="\t",col.names=NA)

}

######
#Sigfit: ssDNA COLLAPSED mutational signatures analysis and plotting
######

cat("#### ssDNA collapsed mutation signature analysis...\n")

if(nrow(mutation_frequencies$SigFitcollapsed$ssDNA)==0){
	cat("       No ssDNA collapsed mutations in any samples.\n\n")
}else{
	##Make SigFit ssDNA collapsed mutation catalogue
	SigFit_ssDNAcollapsed_catalogue <- build_catalogues(mutation_frequencies$SigFitcollapsed$ssDNA)
	 #Fill in samples with 0 mutation counts
	for(i in setdiff(samplenames,rownames(SigFit_ssDNAcollapsed_catalogue))){
		SigFit_ssDNAcollapsed_catalogue <- rbind(SigFit_ssDNAcollapsed_catalogue,rep(0,ncol(SigFit_ssDNAcollapsed_catalogue)))
		rownames(SigFit_ssDNAcollapsed_catalogue)[nrow(SigFit_ssDNAcollapsed_catalogue)] <- i
	}
	SigFit_ssDNAcollapsed_catalogue <- SigFit_ssDNAcollapsed_catalogue[samplenames,]
	
	##Create matrix of mutational opportunities for each sample. This is calculated for each sample as follows:
	#	1) Separate trinucleotide frequencies for central pyrimidine from those for central purine, and reverse complement the labels of the central purine trinucleotide frequencies.
	#	2) Retrieve 96 central pyrimidine trinucleotide contexts in correct order from each of the 32 trinucleotide context central pyrimidine and central purine trinucleotide frequencies.
	# 3) Sum the 96 central pyrimidine with the reverse complement 96 central purine trinucleotide frequencies.
	# 4) Divide the 96 trinucleotide frequencies by their sum
	# 5) Divide the human reference genome trinucleotide frequencies by their sum.
	# 6) Divide (4) by (5).
	# This produces a mutational opportunity matrix reflecting the residual over/under opportunity of our method beyond the baseline trinucleotide distribution of the human genome. For example, if a specific trinucleotide was under-observed in our reads relative to the human genome, it will have a mutation opportunity < 1, and then during signature fitting, sigfit will increase the mutation counts for that trinucleotide.
	mutation_opportunities <- matrix(nrow=numsamples,ncol=length(genome_freqs_labels),dimnames=list(samplenames,genome_freqs_labels))
	
	for(i in samplenames){
		mutationreads_opportunities <- mutation_frequencies$trinucleotide_contexts[[i]]$fwdrev$mutationreadsfreq

		mutationreads_opportunities_pyr <- mutationreads_opportunities[substr(names(mutationreads_opportunities),2,2) %in% c("C","T")]
		mutationreads_opportunities_pur <- mutationreads_opportunities[substr(names(mutationreads_opportunities),2,2) %in% c("A","G")]
		names(mutationreads_opportunities_pur) <- reverseComplement(DNAStringSet(names(mutationreads_opportunities_pur)))
		
		mutationreads_opportunities_pyr <- mutationreads_opportunities_pyr[genome_freqs_labels]
		mutationreads_opportunities_pur <- mutationreads_opportunities_pur[genome_freqs_labels]
		
		mutation_opportunities[i,] <- mutationreads_opportunities_pyr + mutationreads_opportunities_pur
		
		mutation_opportunities[i,] <- (mutation_opportunities[i,]/sum(mutation_opportunities[i,]))/(human_trinuc_freqs()/sum(human_trinuc_freqs()))

	}

	##Plot mutation spectra with y-axis of either counts or probability, and for probability y-axis both uncorrected and corrected for mutation opportunities (3 total PDFs)
	plot_spectrum(SigFit_ssDNAcollapsed_catalogue,pdf_path=paste0(yaml.config$analysisoutput_path,"/sigfit.ssDNAcollapsed.counts.uncorrected.pdf"))
	
	ssDNAcollapsed.probability.uncorrected <- SigFit_ssDNAcollapsed_catalogue/rowSums(SigFit_ssDNAcollapsed_catalogue)
		#Change samples with no mutations to all 0s instead of NA
	for(i in 1:nrow(ssDNAcollapsed.probability.uncorrected)){
		if(all(is.na(ssDNAcollapsed.probability.uncorrected[i,]))){
			ssDNAcollapsed.probability.uncorrected[i,] <- rep(0,length(ssDNAcollapsed.probability.uncorrected[i,]))
		}
	}
	plot_spectrum(ssDNAcollapsed.probability.uncorrected,pdf_path=paste0(yaml.config$analysisoutput_path,"/sigfit.ssDNAcollapsed.probability.uncorrected.pdf"))
	
	ssDNAcollapsed.probability.corrected <- (SigFit_ssDNAcollapsed_catalogue[samplenames,]/mutation_opportunities[samplenames,])/rowSums(SigFit_ssDNAcollapsed_catalogue[samplenames,]/mutation_opportunities[samplenames,])
		#Change samples with no mutations to all 0s instead of NA
	for(i in 1:nrow(ssDNAcollapsed.probability.corrected)){
		if(all(is.na(ssDNAcollapsed.probability.corrected[i,]))){
			ssDNAcollapsed.probability.corrected[i,] <- rep(0,length(ssDNAcollapsed.probability.corrected[i,]))
		}
	}
	plot_spectrum(ssDNAcollapsed.probability.corrected,paste0(yaml.config$analysisoutput_path,"/sigfit.ssDNAcollapsed.probability.corrected.pdf"))
	
	##Fit signatures to known COSMIC signatures, and plot spectra, exposures, and reconstruction
	SigFit_ssDNAcollapsed_fit_signatures <- list()
	for(i in samplenames){
		SigFit_ssDNAcollapsed_fit_signatures[[i]] <- fit_signatures(counts=SigFit_ssDNAcollapsed_catalogue[i,,drop=F],
																						signatures=cosmic_signatures_v3.2,
																						model="multinomial",
																						opportunities=mutation_opportunities[i,,drop=F],
																						iter=10000,warmup=5000,
																						chains=4,seed = 1756,control=list(adapt_delta = 0.99),
																						show_messages=FALSE, open_progress=FALSE)
	}
	mutation_signatures$SigFit_ssDNAcollapsed_fit_signatures <- SigFit_ssDNAcollapsed_fit_signatures
	
	SigFit_ssDNAcollapsed_exposures <- list()
	for(i in samplenames){
		SigFit_ssDNAcollapsed_exposures[[i]] <- retrieve_pars(SigFit_ssDNAcollapsed_fit_signatures[[i]],par="exposures")
	}
	
	for(i in samplenames){
		plot_all(SigFit_ssDNAcollapsed_fit_signatures[[i]],out_path=yaml.config$analysisoutput_path,prefix=paste0(i,".sigfit.ssDNAcollapsed"))
	}
	
	##Redo fit_signatures only with exposures whose lower end of the Bayesian HPD interval is above a specified threshold, plus SBS1, SBS5, SBS18, and SBS40 that are commonly seen in many healthy tissues, and plot spectra, exposures, and reconstruction.
	lowerhpd_thresh <- 0.01
	healthy_signatures <- c("SBS1","SBS5","SBS18","SBS40")
	
	SigFit_ssDNAcollapsed_fit_signatures_thresh <- list()
	for(i in samplenames){
		selectedsignatures <- str_sort(unique(c(rownames(cosmic_signatures_v3.2)[SigFit_ssDNAcollapsed_exposures[[i]]$lower_95>lowerhpd_thresh],healthy_signatures)),numeric=TRUE)
		SigFit_ssDNAcollapsed_fit_signatures_thresh[[i]] <- fit_signatures(counts=SigFit_ssDNAcollapsed_catalogue[i,,drop=F],
																						signatures=cosmic_signatures_v3.2[selectedsignatures,,drop=F],
																						model="multinomial",
																						opportunities=mutation_opportunities[i,,drop=F],
																						iter=10000,chains = 4,seed = 1756,control = list(adapt_delta = 0.99),
																						show_messages=FALSE, open_progress=FALSE)
	}
	mutation_signatures$SigFit_ssDNAcollapsed_fit_signatures_thresh <- SigFit_ssDNAcollapsed_fit_signatures_thresh
	
	SigFit_ssDNAcollapsed_exposures_thresh <- list()
	for(i in samplenames){
		SigFit_ssDNAcollapsed_exposures_thresh[[i]] <- retrieve_pars(SigFit_ssDNAcollapsed_fit_signatures_thresh[[i]],par="exposures")
	}
	
	for(i in samplenames){
		plot_all(SigFit_ssDNAcollapsed_fit_signatures_thresh[[i]],out_path=yaml.config$analysisoutput_path,prefix=paste0(i,".sigfit_thresh.ssDNAcollapsed"))
	}
	
}

qsave(mutation_signatures,paste0(yaml.config$analysisoutput_path,"/mutation_signatures.RDS"))
cat("DONE\n")
