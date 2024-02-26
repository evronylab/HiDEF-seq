########
#Configuration
########

#Load libraries
library(tidyverse)
library(sigfit)

#Load spectra, opportunities, and calls
spectra_opportunities <- readRDS("HiDEF-seq_Liu_spectra_and_opportunities.RDS")
calls <- readRDS("HiDEF-seq_Liu_calls.RDS")

#Output path
outputpath <- "~/Desktop/"

ssDNA.singlenuc.contextorder <- c("C>A","C>G","C>T","T>A","T>C","T>G","G>T","G>C","G>A","A>T","A>G","A>C")
dsDNA.singlenuc.contextorder <- c("C>A","C>G","C>T","T>A","T>C","T>G")

########
#Plots
########

#Select samples to plot with format [sample ID].[subject ID] from Suppl Table 1
samples_to_plot <- c("1.1443","2.1443","3.1443","4.1409","5.1409","6.1104") 

##
## A. dsDNA trinucleotide spectra, corrected for trinucleotide context opportunities
##
corrected_counts <- spectra_opportunities$dsDNAmutations_spectra[samples_to_plot,,drop=F]/spectra_opportunities$dsDNAmutations_opportunities[samples_to_plot,,drop=F]

plot_spectrum(corrected_counts/rowSums(corrected_counts),paste0(outputpath,"dsDNA.tri.spectra.pdf"))

##
## B. ssDNA trinucleotide spectra, corrected for trinucleotide context opportunities
##
# Note: in plots, Transcribed = central pyrimidine calls; Untranscribed = central purine calls
corrected_counts <- spectra_opportunities$ssDNAcalls_spectra[samples_to_plot,,drop=F]/spectra_opportunities$ssDNAcalls_opportunities[samples_to_plot,,drop=F]

plot_spectrum(corrected_counts/rowSums(corrected_counts),paste0(outputpath,"ssDNA.tri.spectra.pdf"))

##
## C. ssDNA and dsDNA single-nucleotide spectra (side-by-side barplots and stacked barplots) and ssDNA pyrimidine/purine fractions
##

#Format data frames
correctedcounts.ssDNA <- (spectra_opportunities$ssDNAcalls_spectra/spectra_opportunities$ssDNAcalls_opportunities) %>%
  as.data.frame %>%
  rownames_to_column("samplelabel") %>%
  pivot_longer(-samplelabel,names_to="context",values_to="counts") %>%
  mutate(
    samplelabel=fct_relevel(samplelabel,unique(samplelabel)),
    pyrpur.context = if_else(str_detect(context,"^T"),"pyrimidine","purine"),
    pyrpur.context=fct_relevel(pyrpur.context,c("pyrimidine","purine")),
    singlenuc.context = str_c(
      str_sub(context,1,2),
      str_sub(context,4,4),
      ">",
      str_sub(context,8,8)
    ),
    singlenuc.context = if_else(str_detect(singlenuc.context,"^T"),
                                str_sub(singlenuc.context,3,5),
                                str_c(
                                  str_sub(singlenuc.context,3,3) %>% DNAStringSet %>% reverseComplement %>% as.character,
                                  ">",
                                  str_sub(singlenuc.context,5,5) %>% DNAStringSet %>% reverseComplement  %>% as.character
                                )
    ),
    singlenuc.context=fct_relevel(singlenuc.context,ssDNA.singlenuc.contextorder)
  )

correctedcounts.dsDNA <- (spectra_opportunities$dsDNAmutations_spectra/spectra_opportunities$dsDNAmutations_opportunities) %>%
  as.data.frame %>%
  rownames_to_column("samplelabel") %>%
  pivot_longer(-samplelabel,names_to="context",values_to="counts") %>%
  mutate(
    samplelabel=fct_relevel(samplelabel,unique(samplelabel)),
    singlenuc.context = str_c(
      str_sub(context,2,2),
      ">",
      str_sub(context,6,6)
    ),
    singlenuc.context=fct_relevel(singlenuc.context,dsDNA.singlenuc.contextorder)
  )

#ssDNA single-nucleotide spectra (side-by-side barplots)
correctedcounts.ssDNA %>%
  filter(samplelabel %in% samples_to_plot) %>%
  group_by(samplelabel,singlenuc.context) %>%
  summarize(counts=sum(counts),.groups="drop") %>%
  group_by(samplelabel) %>%
  mutate(frac=counts/sum(counts)) %>%
  ggplot(aes(x=singlenuc.context,y=frac,fill=singlenuc.context)) +
  geom_col(position = "dodge") +
  facet_wrap(vars(samplelabel),ncol=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45,hjust=1)) +
  guides(fill=guide_legend(ncol=2)) +
  scale_fill_brewer(type="qual",palette="Paired") +
  labs(x="Context",y="Fraction of ssDNA calls",fill="Context")
ggsave(paste0(outputpath,"ssDNA.singlenuc.spectrum.pdf"))

#dsDNA single-nucleotide spectra (side-by-side barplots)
correctedcounts.dsDNA %>%
  filter(samplelabel %in% samples_to_plot) %>%
  group_by(samplelabel,singlenuc.context) %>%
  summarize(counts=sum(counts),.groups="drop") %>%
  group_by(samplelabel) %>%
  mutate(frac=counts/sum(counts)) %>%
  ggplot(aes(x=singlenuc.context,y=frac,fill=singlenuc.context)) +
  geom_col(position = "dodge") +
  facet_wrap(vars(samplelabel),ncol=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45,hjust=1)) +
  scale_fill_manual(values=c("#00bfff","#000000","#ee2c2c","#c2c2c2","#a2cd5a","#eeb4b4")) +
  labs(x="Context",y="Fraction of dsDNA mutations",fill="Context")
ggsave(paste0(outputpath,"dsDNA.singlenuc.spectrum.pdf"))

#ssDNA single-nucleotide spectra (stacked barplots)
correctedcounts.ssDNA %>%
  filter(samplelabel %in% samples_to_plot) %>%
  group_by(samplelabel,singlenuc.context) %>%
  summarize(counts=sum(counts),.groups="drop") %>%
  group_by(samplelabel) %>%
  mutate(frac=counts/sum(counts)) %>%
  ggplot(aes(x=samplelabel,y=frac,fill=singlenuc.context)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45,hjust=1)) +
  guides(fill=guide_legend(ncol=2)) +
  scale_fill_brewer(type="qual",palette="Paired") +
  labs(x="Sample",y="Fraction of ssDNA calls",fill="Context")
ggsave(paste0(outputpath,"ssDNA.singlenuc.spectrum.stacked.pdf"))

#dsDNA single-nucleotide spectra (stacked barplots)
correctedcounts.dsDNA %>%
  filter(samplelabel %in% samples_to_plot) %>%
  group_by(samplelabel,singlenuc.context) %>%
  summarize(counts=sum(counts),.groups="drop") %>%
  group_by(samplelabel) %>%
  mutate(frac=counts/sum(counts)) %>%
  ggplot(aes(x=samplelabel,y=frac,fill=singlenuc.context)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45,hjust=1)) +
  scale_fill_manual(values=c("#00bfff","#000000","#ee2c2c","#c2c2c2","#a2cd5a","#eeb4b4")) +
  labs(x="Sample",y="Fraction of dsDNA mutations",fill="Context")
ggsave(paste0(outputpath,"dsDNA.singlenuc.spectrum.stacked.pdf"))

#ssDNA pyrimidine/purine fractions
correctedcounts.ssDNA %>%
  filter(samplelabel %in% samples_to_plot) %>%
  group_by(samplelabel,pyrpur.context) %>%
  summarize(counts=sum(counts),.groups="drop") %>%
  group_by(samplelabel) %>%
  mutate(frac=counts/sum(counts)) %>%
  ggplot(aes(x=samplelabel,y=frac,fill=pyrpur.context)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45,hjust=1)) +
  labs(x="Sample",y="Fraction of ssDNA calls",fill="Context")
ggsave(paste0(outputpath,"ssDNA.pyrpur.spectrum.stacked.pdf"))