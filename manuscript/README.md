## Example analyses for the Liu, Costa, et al. manuscript
This folder contains the spectra and trinucleotide-context opportunities matrices, and also all ssDNA and dsDNA calls for all samples in the Liu, Costa, et al. manuscript.

### Files
1. HiDEF-seq_Liu_spectra_and_opportunities.RDS (R RDS format)
This list object contains 4 matrices, containing data for chromosomes 1-22,X,Y data (i.e., excludes chrM data):
   * ssDNAcalls_spectra: counts of ssDNA calls in each of the 192 trinucleotide contexts. Rows = samples (rownames: sampleID.subjectID). Columns = trinucleotide contexts (T:central pyrimidine; U:central purine).
   * dsDNAmutations_spectra: counts of dsDNA mutations in each of the 96 trinucleotide contexts. Rows = samples (rownames: sampleID.subjectID). Columns = trinucleotide contexts.
   * ssDNAcalls_opportunities: trinucleotide context opportunities for ssDNA analysis (i.e. ratio of the fraction of interrogated bases corresponding to the trinucleotide over the fraction of genome bases corresponding to the trinucleotide).
   * dsDNAmutations_opportunities: trinucleotide context opportunities for dsDNA analysis.

2. HiDEF-seq_Liu_calls.RDS (R RDS format)
This list object contains 2 data frames, containing data for chromosomes 1-22,X,Y, and M. Column definitions can be found [here](../docs/print_mutations.md).
   * ssDNAcalls: all ssDNA calls.
   * dsDNAmutations: all dsDNA mutations.
  
### Example analyses
HiDEF-seq_Liu_analysis-example.R in this folder contains example analyses.
1. First download the above RDS files.
2. Configure the script 'readRDS' commands to the path of the above downloaded RDS files.
3. Specify the samples to plot in the 'sample_to_plot' line.
4. Run the desired analyses/plots, which includes the following:
   * dsDNA trinucleotide spectra, corrected for trinucleotide context opportunities
   * ssDNA trinucleotide spectra, corrected for trinucleotide context opportunities
   * ssDNA and dsDNA single-nucleotide spectra (side-by-side barplots and stacked barplots)
   * ssDNA pyrimidine/purine fractions
