# YAML Config parameter documentation

# General tools
  1. condabase_script: $CONDA_BASE/etc/profile.d/conda.sh. This prepares the shell for activating conda environment. If using singularity, set to: /hidef/mambaforge/etc/profile.d/conda.sh
  2. conda_pbbioconda_env: path to conda environment with PacBio tools. If using singularity, set to: /hidef/bin/pbconda
  3. wigToBigWig_bin: path to wigToBigWig binary. If using singularity, set to: /hidef/bin/wigToBigWig
  4. bcftools_bin: path to bcftools binary. If using singularity, set to: bcftools
  5. samtools_bin: path to samtools >=1.12 binary. If using singularity, set to: samtools
  6. samtools_1.10_bin: path to samtools v1.10 binary for process_subreads scripts. If using singularity, set to: /hidef/bin/samtools_1.10/samtools
  7. seqkit_bin: path to seqkit binary. If using singularity, set to: /hidef/bin/seqkit
  8. kmc_bindir: path to directory containing kmc and kmc_dump binaries (https://github.com/refresh-bio/KMC). If using singularity, set to: /hidef/bin/kmc/bin
  9. gsutil_bin: Full path to gsutil binary on user's SLURM cluster account. GCloud must be installed and initialized to user's desired project.
  10. HiDEFpipeline: path to directory containing all the pipeline scripts (process_subreads.sh, postccsfilter.awk, baseposition.awk, EXTRACTREFSEQscript.awk, refseq.awk, postalign.sh, and pbmm2filter.awk). If using singularity, set to: /hidef/scripts
  11. slurm_himemjob_options: options to pass to sbatch for high memory jobs, separated by spaces as you would normally pass them to sbatch (e.g. -p himem_partition). This must also specify --mem=#G (recommended --mem=256G), where # = memory for high memory jobs (mutation_filtering, print_mutations, mutation_frequencies).
  12. slurm_add_options: additional options to pass to sbatch (can be empty)

# Genome reference information and files
  1. genome: ID of reference genome. If using singularity, set to: CHM13_v1
  2. chrsizes: tab-separated file with chromosome names (column 1) and chromosome size (column 2). If using singularity, set to: /hidef/references_CHM13v1/chrsizes.tsv
  3. chrs_to_analyze: names of chromosomes (i.e. "chr1", etc.) to analyze, one per line. If using singularity, set to: /hidef/references_CHM13v1/chr1-22X.txt
  4. Nrefbed: Compressed 3-column BED file of positions in genome with 'N' sequence. If using singularity, set to: /hidef/references_CHM13v1/chm13.draft_v1.0.N.merged.bed
  5. fastaref: reference genome fasta
  6. genomemmiindex: .mmi index of genome reference for use in pbmm2 alignment of mutation subreads
  7. BSgenomepackagename: Name of BSgenome package in R when installed, so it can be loaded with library(). If using singularity, set to: BSgenome.Hsapiens.T2T.CHM13v1
  8. gnomad_sensitivity_ref: bigwig file of gnomad variants (with bigwig score >0 for all positions) to intersect with sample's variants for germline heterozygous variant sensitivity estimate. If using singularity, set to: /hidef/references_CHM13v1/af.only.0.001.gnomad.v3.1.2.CHM13_v1.bw

# Sample configuration
 Sample names
  1. samplenames: Array of sample names. Sample names must be unique across the entire project.

 Tissues
  1. tissues: Array of tissues for each sample. Must match order of samplenames.

 Barcodes
  1. barcodes: Array of barcodes for each sample, in format bc####:XXXXXXXXXXXXXXXX, where #### = barcode ID number, and XXXXXXXXXXXXXXXX = barcode sequence. Must match order of samplenames.

# Configuration for process_subreads.sh
  1. subreads_filename: subreads file (full original SMRTcell subreads data); also requires PBI index file in same location.
  2. process_subreads_output_path: output directory for process_subreads scripts
  3. ccs_BAM_prefix: prefix for final aligned ccs BAM
  4. ccschunks: number of chunks used in pbccs
  5. minccsrq: minimum predicted quality of ccs consensus. parameter is used in pbccs.
  6. minoverlap: The minimum reciprocal overlap between the forward primary alignment and reverse primary alignment. Used in pbmm2filter.awk.
  7. remotetempdir: true or false. Whether to remove TMPoutput directory created by process_subreads.sh upon completion.

# Configuration for make_bamfilezmw_all.R
  1. bamfilezmw_all_filename: full path and file name for output RDS file created by make_bamfilezmw_all.R and the input for subsequent analysis scripts
  2. make_bamfilezmw_all_config:
      - First sub-level are samplenames (must match prior samplenames). For each, must define the below.
        - vcffile_filenames: VCF file(s) for annotation of bamfilezmw_all in make_bamfilezmw_all.R script. These are usually the same as the VCFs in vcffilters.
 
# Configuration for mutation_filtering.R
 Bigwig cache directory
  1. bigwig_cachedir: directory to store bigwig files generated during analysis, to save significant time in future analyses

 Output path and subreads bam filename
  1. analysisoutput_path: path in which to save analysis output
  2. analysisoutput_basename: base name to use for output files. This includes a log file and bamfilezmw_filtered RDS file with a list of all final mutations and other annotation generated during analysis.
  3. gcloudoutput_path: [optional], if defined, HiDEF_pipeline.sh will copy final aligned CCS and subread BAM/BAI files to this gcloud path using 'gsutil cp'.

# Mutation burden analysis filters:
For somatic mutation rate analysis, multiple distinct sets of filter configurations can be defined across the following filter configuration lists, each with a different filterset index number.
Multiple sets of filter combinations can be set up so that the analysis iterates through them and calculates a mutation rate for each. The categories of filters that are assigned to filter sets are:
  - basic filters (called 'thresholds')
  - vcffilters
  - bigwigfilters
  - bamfilters
  - mutratefilters.

Each filter set is denoted by a numeric integer index (starting at 1). Each of these categories of filters must be configured with the same number of filter set indexes; i.e. the number of filter sets in basic filters must match the number of configuration filter sets defined in vcffilters, bigwigfilters, bamfilters, and mutratefilters. However, to avoid tediously large configuration files, for vcffilters, bigwigfilters, bamfilters, mutratefilters, if ...sameforallbasicfiltersets == TRUE, then the first configuration for that type of filter is duplicated across across all filter sets for that filter with a total number of filter sets equal to the number of filter sets in basic filters. This is because basic filters are the most likely to be optimized.

When there are _fwd, _rev, _avg  options for a filter, these are filters applied to fwd, rev, and average of fwd and rev strand reads, respectively, where fwd/rev strand corresponds to the reference genome strand to which the strand sequence aligned.

<ins>Details of filters:</ins>

**thresholds:**
- First sub-level are numeric integer indexes (starting at 1) for each filter set. Under these are:
  - *General filters*
    - rqfwd, rqrev, rqavg (minimum read quality)
    - ecfwd, ecrev, ecavg (minimum effective coverage)
    - mapqfwd, mapqrev, mapqavg (minimum mapping quality)
    - numsnvsfwd, numsnvsrev, numsnvszmw (maximum number of SNVs before VCF filtering; fwd and rev include both ssDNA and dsDNA calls, and zmw are dsDNA calls; set to 999999 to disable)
    - numsnvsfwdrevdiff (maximum difference in number of SNVs between fwd and rev CCS before VCF filtering)
    - numindelsfwd, numindelsrev, numindelszmw (maximum number of indels before VCF filtering; set to 999999 to disable)
    - numsnvsfwdpostVCF, numsnvsrevpostVCF, numsnvszmwpostVCF (maximum number of SNVs right after VCF SNV filtering and basic ZMW filtering, and before other SNV filtering; set to 999999 to disable)
    - numsoftclipfwd, numsoftcliprev, numsoftclipavg (maximum number of soft clipped bases; set to 999999 to disable)
    - minqq: Minimum base quality filter. Both strands of zmw (duplex) calls must pass filter, and single-strand calls must pass filter on the strand containing the call.
    - bpends: # of basepairs to filter from 5'/3' ends of read alignments (alignment span excludes soft-clipped bases). For dsDNA calls, filter is applied to both fwd and rev strand alignments.

  - *CCS indel filters*
    - ccsindelfilter: Exclude if near ccs indel (true or false)
    - ccsindelinspad: Insertion padding: size of padding on each side of insertion: mb. 'm': multiplier relative to size of insertion (float), or 'b': basepairs. Both are calculated for each call and filtering uses the larger of the two. To exclude one of them, just set it to 0; i.e. m0 or b0.
    - ccsindeldelpad: Deletion padding: size of padding on each side of deletion. Same format as above.

  - *Subread filters:* Filters to require minimium evidence in extracted and aligned subreads. Separate filters for ssDNA and ZMW (dsDNA) calls.
    - minsubreads_cvg_fraction: Minimum fraction of subreads overlapping the calls (regardless of whether they contain the call) out of the total subreads aligned to the genome, taking into account for both the numerator and denominator terms only subreads from the same strand and ZMW in which the call was made. This filter is applied separately to each strand in which the call was made (i.e. only the call-containing strand for ssDNA calls, and to both strands for dsDNA calls so that a dsDNA call must pass this filter in both strands). This filter removes calls covered by only a fraction of a ZMW's subreads aligned to the soft-clipped region.
    - minZMWsubreadsVariantReads: Minimum number of subreads that match ZMW variant base required in both fwd and rev subreads (filtered applied to each separately, and required to pass in both).
    - minZMWsubreadsVAF: Minimum VAF among subreads that match ZMW variant base required in both fwd and rev subreads (filtered applied to each separately, and required to pass in both).
    - minssDNAsubreadsVariantReads: Minimum number of subreads that match ssDNA variant called in the consensus CCS read in order to keep the ssDNA variant. Applies only to the subreads from the same strand as the ssDNA variant.
    - minssDNAsubreadsVAF: Minimum VAF among subreads that match ssDNA variant called in the consensus CCS read in order to keep the ssDNA variant. Applies only to the subreads from the same strand as the ssDNA variant.

**vcffilters_sameforallbasicfiltersets:** true or false; whether the configuration for the first filter set is used across all basic filter filter sets
 
**vcffilters:**
- First sub-level are integer indexes (starting at 1) for each filter set.
  - Second sub-level are samplenames (must match prior samplenames)
    - Third sub-level are numeric indexes for each VCF file. Multiple VCF files can be specified.
        - vcffilename: VCF full file name path (same used during bamfilezmw_all loading)
        - SNVFILTERS: array of VCF FILTER column values to include for SNV filtering, or "." to use all
        - INDELFILTERS: array of VCF FILTER column values to include for INDEL filtering, or "." to use all
        - vcfSNVfilter: Perform SNV filtering (true or false)
        - vcfSNVDepth: Minimum total locus depth for SNVs
        - vcfSNVGQ: Minimum genotype quality for SNVs
        - vcfSNVVAF: Minimum variant allele frequency for SNVs (range 0-1)
        - vcfSNVQUAL: Minimum variant QUAL for SNVs (column 6 of VCF)
        - vcfINDELfilter: Perform INDEL filtering (true or false)
        - vcfINDELDepth: Minimum total locus depth for INDELs
        - vcfINDELGQ: Minimum genotype quality for INDELs
        - vcfINDELVAF: Minimum variant allele frequency for INDELs (range 0-1)
        - vcfINDELQUAL: Minimum variant QUAL for INDELs (column 6 of VCF)
        - vcfINDELinspad: size of padding around insertions for filtering (filter variants inside padding): mb. 'm': multiplier relative to size of insertion (float), or 'b': basepairs. Both are calculated for each call and filtering uses the larger of the two for each call. To exclude one of them, just set it to 0; i.e. m0 or b0.
        - vcfINDELdelpad: size of padding around deletions for filtering (filter variants inside padding): format per above.

**bigwigfilters_sameforallbasicfiltersets:** true or false; whether the configuration for the first filter set is used across all basic filter filter sets

**bigwigfilters:**
- First sub-level are integer indexes (starting at 1) for each filter set.
  - Second sub-level are numeric indexes for each bigwig file. Multiple bigwig files can be specified.
      - bigwigfiltertype: SNV or ZMW. This specifies if it is an SNV filter or ZMW filter.
      - bigwigfile: Bigwig full file name path
      - binsize: Size of bins across the genome to average bigwig signal (integer).
      - threshold: "[gte|lt]###".
        - If "SNV" filter: defines regions to be filtered as bins greater than or equal (gte) or less than (lt) numeric threshold ### (### = floating point between 0 to 1).
        - If "ZMW" filter, the fraction of the read span that is greater than (gte) or less than (lt) numeric threshold ### (### = floating point between 0 to 1) is calculated separately for the fwd and rev ccs reads. The average of these two fractions is then compared to zmwthreshold.
      - zmwthreshold: "[gte|lt]###". The entire ZMW is filtered out if the average fraction of the fwd and rev ccs read spans that do not pass the 'threshold' parameter above is greater than or equal (gte) or less than (lt) numeric zmwthreshold ### (### = floating point between 0 to 1), and otherwise this filter does not remove the ZMW nor any part of the ZMW. Relevant only for "ZMW" filters, blank for "SNV" filters.
      - padding:  of bp padding on both sides of each filtered region. Relevant only for "SNV" filters.
           Note: Can specify the same bigwig file twice, once as an "SNV"-filter (zmwthreshold is ignored and doesn't need to be specified), and a second time as a "ZMW" filter by setting snvorzmwfilter="ZMW".
      - applyto: "ss", "ds", or "both" -> Specifies if the filter will be applied to ssDNA or dsDNA or both types of calls. Relevant only for filters for which bigwigfiltertype = "SNV".


**bamfilters_sameforallbasicfiltersets:** true or false; whether the configuration for the first filter set is used across all basic filter filter sets
 
**bamfilters:** Configure bulk WGS BAM/CRAM file and WGS coverage files for filtering. This is used for BAM/CRAM bulk WGS pileup filtering. These are needed because GATK and DeepVariant sometimes miss low-level variants or due to wrong local haplotype assembly that calls variants in a different nearby location than in the CCS/ZMW.
- First sub-level are integer indexes (starting at 1) for each filter set.
  - Second sub-level are samplenames (must match prior samplenames)
      - bamfile: Full path of BAM/CRAM file of bulk WGS. Must be a single sample BAM/CRAM file, and index must be present.
      - bamreftype: Illumina or PacBio, depending on the type of the bulk sequencing file used for WGS pileup filtering. This changes the bcftools mpileup command used.
      - WGScoveragebigwig: Bigwig file of coverage made using samtools mpileup with the same parameters used for final SNV pileup filtering.
      - minBAMTotalReads: Minimum total read depth in bulk WGS BAM/CRAM file. If less than this, then filter the variant. This ensures that germline reads aren't called as somatic just because the bulk WGS data had low read depth. Set to 0 to disable.
      - maxBAMVariantReads: Maximum number of reads with the same variant in the bulk WGS BAM/CRAM file. If greater than this, then filter the variant. Set to 999999 to disable.
      - maxBAMVAF: Maximum VAF of the same variant in the bulk WGS BAM/CRAM file. If greater than this, then filter the variant. Set to 1 to disable.

 **mutratefilters:** These filters determine how the final mismatch and mutation burdens are calculated. These thresholds are applied identically to all samples sharing the same YAML configuration file (i.e. samples from the same analysis).
  1. maxmutationsperssdna: Maximum number of ssDNA calls in fwd or rev strand reads for that strand's ssDNA calls to be included in the final burden calculation.
  2. maxmutationsperzmw: Maximum number of dsDNA calls in a ZMW for the ZMW's dsDNA calls to be included in the final burden calculation.
