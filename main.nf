/*****************************************************************
 * Process Definitions
 *****************************************************************/

/*
  makeBarcodesFasta: Writes the barcode FASTA content to a file from a value input.
*/
process makeBarcodesFasta {
    cpus 1
    memory '2 GB'
    time '10m'
    tag "Make Barcodes Fasta"
    container "${params.hidefseq_container}"

    input:
      tuple val(run_id), val(content)
    
    output:
      tuple val(run_id), path("barcodes.fasta")
    
    script:
      """
      echo -e "${content}" > barcodes.fasta
      """
}

/*
  ccsChunk: Runs one chunk of CCS.

  #Notes:
  #1. Using ccs v8.0 (standalone version of ICS13 ccs, http://doi.org/10.5281/zenodo.10703290 ) that supports rich HiFi tags (--subread-pileup-summary-tags).
  #2. Set '--pbdc' and '--pbdc-skip-min-qv 0' to go into DeepConsensus code path to support rich HiFi tags but skip DeepConsensus polishing (by setting to skip any 100 bp window with average QV >0; i.e. every window)
  #3. Set '--binned-qvs=False' to keep full resolution quality values. Setting pptop-passes to 255 instead of 0 to match what occurs on Revio (since setting it to 0, i.e. unlimited, caused technical issues on Revio).
  #4. Setting --instrument-files-layout --min-rq -1 --movie-name <hifireads.bam> --non-hifi-prefix fail to mimic what happens on Revio, since this outputs min-rq >= 0.99 into hifireads.bam and also removes other types of failed reads (per ff tag details here: https://pacbiofileformats.readthedocs.io/en/13.0/BAM.html#use-of-read-tags-for-fail-per-read-information). We then do independent rq filtering based on pipeline configuration later in the pipeline.
  #5. Off-instrument ccs run with CPUs differs slightly from on-instrument Revio ccs run with GPUs for the parameters --max-insertion-size and window size, due to technical details, but PacBio says this should negligibly affect the output and we should not specify these.
*/
process ccsChunk {
    cpus 8
    memory '32 GB'
    time '24h'
    tag { "CCS chunk ${chunkID}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(run_id), path(bamFile), path(pbiFile), val(chunkID)

    output:
      tuple val(run_id), path("hifi_reads/${run_id}.chunk${chunkID}.hifi_reads.ccs.bam"), path("hifi_reads/${run_id}.chunk${chunkID}.hifi_reads.ccs.bam.pbi"), emit: bampbi_tuple
      path "statistics/*.ccs_report.*", emit: report
      path "statistics/*.summary.json", emit: summary

    publishDir "${params.analysis_output_dir}/logs", mode: 'copy', pattern: "statistics/*.ccs_report.*", saveAs: { filename -> new File(filename).getName() }
    publishDir "${params.analysis_output_dir}/logs", mode: 'copy', pattern: "statistics/*.summary.json", saveAs: { filename -> new File(filename).getName() }

    script:
    // Build the LD_PRELOAD command if the parameter is set.
    def ld_preload_cmd = (params.ccs_ld_preload && params.ccs_ld_preload.trim()) ? "export LD_PRELOAD=${params.ccs_ld_preload}" : ""

    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    ${ld_preload_cmd}
    ccs -j 8 --log-level INFO --by-strand --hifi-kinetics --instrument-files-layout --min-rq -1 --top-passes 255 \\
        --pbdc --pbdc-skip-min-qv 0 --subread-pileup-summary-tags --binned-qvs=False \\
        --chunk ${chunkID}/${params.ccs_chunks} \\
        --movie-name ${run_id}.chunk${chunkID} \\
        --non-hifi-prefix fail \\
        --report-file statistics/${run_id}.chunk${chunkID}.ccs_report.txt \\
        ${bamFile}
    """
}

/*
  mergeCCSchunks: Merges all CCS chunk outputs into a single BAM.
*/
process mergeCCSchunks {
    cpus 2
    memory '8 GB'
    time '6h'
    tag "Merge CCS chunks"
    container "${params.hidefseq_container}"
    
    input:
      tuple val(run_id), path(bamChunks), path(pbiChunks)

    output:
      tuple val(run_id), path("${run_id}.ccs.bam"), path("${run_id}.ccs.bam.pbi")

    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbmerge -o ${run_id}.ccs.bam ${bamChunks.join(' ')}
    """
}

/*
  filterAdapter: Filter CCS BAM file to keep only reads with ma tag == 0 (adapter detected on both ends).
*/
process filterAdapter {
    cpus 8
    memory '16 GB'
    time '10h'
    tag "Filter Bad Adapters"
    container "${params.hidefseq_container}"
    
    input:
      tuple val(run_id), path(bamFile), path(pbiFile)
    
    output:
      tuple val(run_id), path("${run_id}.ccs.filtered.bam"), path("${run_id}.ccs.filtered.bam.pbi")
    
    script:
    """
    ${params.samtools_bin} view -b -@8 -e "[ma]==0" ${bamFile} > ${run_id}.ccs.filtered.bam
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbindex ${run_id}.ccs.filtered.bam
    """
}

/*
  limaDemux: Demultiplexes the filtered BAM using lima.
*/
process limaDemux {
    cpus 8
    memory '64 GB'
    time '12h'
    tag "Lima Demultiplexing"
    container "${params.hidefseq_container}"
    
    input:
      tuple val(run_id), path(bamFile), path(pbiFile), path(barcodesFasta)

    output:
      tuple val(run_id), path("${run_id}.ccs.filtered.demux.*.bam"), emit: bam
      path "*.lima.summary", emit: lima_summary
      path "*.lima.counts", emit: lima_counts
    
    publishDir "${params.analysis_output_dir}/logs", mode: 'copy', pattern: "*.lima.summary"
    publishDir "${params.analysis_output_dir}/logs", mode: 'copy', pattern: "*.lima.counts"
    
    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    lima --ccs --split-named --min-score 80 --min-end-score 50 \\
         --min-ref-span 0.75 --same --min-scoring-regions 2 \\
         ${bamFile} \\
         ${barcodesFasta} \\
         ${run_id}.ccs.filtered.demux.bam
    """
}

/*
  pbmm2Align: Aligns a demultiplexed BAM file using pbmm2.
  and renames the output to include the sample ID.
*/
process pbmm2Align {
    cpus 8
    memory '64 GB'
    time '12h'
    tag { "pbmm2 Alignment: ${sample_id}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(run_id), val(sample_id), path(bamFile)
    
    output:
      tuple val(run_id), val(sample_id), path("${run_id}.${sample_id}.ccs.filtered.aligned.bam"), path("${run_id}.${sample_id}.ccs.filtered.aligned.bam.pbi")
    
    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbmm2 align -j 8 --preset CCS ${params.genome_mmi} ${bamFile} ${run_id}.${sample_id}.ccs.filtered.aligned.bam
    pbindex ${run_id}.${sample_id}.ccs.filtered.aligned.bam
    """
}

/*
  mergeAlignedSampleBAMs: Merges aligned BAM files from the same sample across different runs.
*/
process mergeAlignedSampleBAMs {
    cpus 4
    memory '64 GB'
    time '6h'
    tag { "Merge aligned sample BAMs: ${sample_id}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(sample_id), path(bamFiles), path(pbiFiles)

    output:
      tuple val(sample_id), path("${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.bam"), path("${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.pbi"), path("${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.bai")

    storeDir "${processReads_output_dir}"

    script:
    """
    sample_basename=${params.analysis_id}.${sample_id}.ccs.filtered.aligned

    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbmerge -o \${sample_basename}.unsorted.bam ${bamFiles.join(' ')}
    conda deactivate
    
    ${params.samtools_bin} sort -@4 -m 4G \${sample_basename}.unsorted.bam > \${sample_basename}.sorted.bam
    ${params.samtools_bin} index -@4 \${sample_basename}.sorted.bam

    conda activate ${params.conda_pbbioconda_env}
    pbindex \${sample_basename}.sorted.bam
    """
}

/*
  countZMWs: Runs zmwfilter on an input BAM file and writes the ZMW count to a file.
*/
process countZMWs {
    cpus 1
    memory '8 GB'
    time '1h'
    tag "Count ZMWs"
    container "${params.hidefseq_container}"
    
    input:
      tuple path(bamFile), path(pbiFile), val(outFileSuffix)
    
    output:
      path "*"
    
    publishDir "${params.analysis_output_dir}/logs", mode: 'copy'
    
    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    zmwfilter --show-all ${bamFile} | wc -l > \$(basename ${bamFile} .bam).${outFileSuffix}
    """
}

/*
  splitBAM: Splits BAM files into approximately equal chunks per params.analysis_chunks.
*/
process splitBAM {
    cpus 2
    memory '32 GB'
    time '4h'
    tag { "Split BAM: ${sample_id}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(sample_id), path(bamFile), path(pbiFile), path(baiFile), val(chunkID)
    
    output:
      tuple val(sample_id), path("${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam"), path("${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam.pbi"), path("${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam.bai"), val(chunkID)

    storeDir "${splitBAMs_output_dir}"

    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}

    sample_basename=\$(basename ${bamFile} .bam)
    
    zmwfilter --show-all ${bamFile} > \${sample_basename}.zmwIDs.txt
    total_zmws=\$(wc -l < \${sample_basename}.zmwIDs.txt)
    zmws_per_chunk=\$(( (total_zmws + ${params.analysis_chunks} - 1) / ${params.analysis_chunks} ))
    split -a 4 --numeric-suffixes=1 -l \$zmws_per_chunk \${sample_basename}.zmwIDs.txt \${sample_basename}.zmwIDs.chunk.

    chunk_file=\$(ls \${sample_basename}.zmwIDs.chunk.* | sort | sed -n "${chunkID}p")
    chunk_bam=\${sample_basename}.chunk${chunkID}.bam

    zmwfilter --include \$chunk_file ${bamFile} \$chunk_bam
    pbindex \$chunk_bam

    conda deactivate

    ${params.samtools_bin} index -@2 \$chunk_bam
    """
}

/*
  installBSgenome: Run installBSgenome.R
*/
process installBSgenome {
    cpus 1
    memory '16 GB'
    time '1h'
    tag { "Install BSgenome reference" }
    container "${params.hidefseq_container}"
    cache false //Always run this process because the BSgenome could have been deleted outside nextflow and because the script itself checks if the BSgenome is already installed.
      
    output:
      val(true)

    script:
    """
    installBSgenome.R -c ${params.paramsFileName}
    """
}

/*
  processGermlineVCFs: Run processGermlineVCFs.R
*/
process processGermlineVCFs {
    cpus 1
    memory '16 GB'
    time '4h'
    tag { "Process Germline VCFs" }
    container "${params.hidefseq_container}"
      
    input:
      val(individual_id)

    output:
      path "${individual_id}.germline_vcf_variants.qs2"

    storeDir "${params.cache_dir}"

    script:
    """
    processGermlineVCFs.R -c ${params.paramsFileName} -i ${individual_id}
    """
}

/*
  processGermlineBAMs: Run processGermlineBAMs.R
*/
process processGermlineBAMs {
    cpus 1
    memory '16 GB'
    time '24h'
    tag { "Process Germline BAMs" }
    container "${params.hidefseq_container}"
    
    input:
      tuple path(germline_bam_file), val(germline_bam_type)

    output:
      path("${germline_bam_file}.bw")

    storeDir "${params.cache_dir}"

    script:
    """
    #Output per-base coverage using samtools mpileup with the similar filters as used by bcftools mpileup in later call
    #filtering. This ensures this coverage data for calculating the fraction of the genome that was filtered
    #maintains correct calculation of the mutation rate.

    #Note using samtools mpileup here instead of bcftools mpileup that is later used in call filtering, because here we
    #need an output that can be re-formatted into bedgraph. Whereas in the later call filtering, we use bcftools mpileup
    #to check for low mosaicism at call positions, and the output of that is in VCF format that is easier to parse. The
    #difference between samtools mpileup and bcftools mpileup should not be significant.

    #Slightly different parameters are used for Illumina vs PacBio germline BAM to match how bcftools mpileup is run later
    #in the call filtering analysis.
    
    set -euo pipefail

    if [[ ${germline_bam_type} == Illumina ]]; then
      ${params.samtools_bin} mpileup -A -B -Q 11 -d 999999 --ff 3328 -f ${params.genome_fasta} [BAM file] 2>/dev/null | awk '{print \$1 "\t" \$2-1 "\t" \$2 "\t" \$4}' > mpileup.bg
    elif [[ ${germline_bam_type} == PacBio ]]; then
      ${params.samtools_bin} mpileup -A -B -Q 5 -d 999999 --ff 3328 -f ${params.genome_fasta} [BAM file] 2>/dev/null | awk '{print \$1 "\t" \$2-1 "\t" \$2 "\t" \$4}' > mpileup.bg
    else
      echo "ERROR: Unknown germline_bam_type: ${germline_bam_type}"
      exit 1
    fi

    sort -k1,1 -k2,2n mpileup.bg > mpileup.sorted.bg

    ${params.bedGraphToBigWig_bin} mpileup.sorted.bg <(cut -f 1,2 ${params.genome_fai}) ${germline_bam_file}.bw
    """
}

/*
  processGermlineBAMs: Run processGermlineBAMs.R
*/
process prepareRegionFilters {
    cpus 2
    memory '64 GB'
    time '24h'
    tag { "Prepare region filters: ${region_filter_file}, bin ${binsize}, threshold ${threshold}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple path(region_filter_file), val(binsize), val(threshold)

    output:
      path("${region_filter_file}.bin${binsize}.${threshold}.bw")

    storeDir "${params.cache_dir}"

    script:
    """
    #Make genome BED file to use to fill in zero values for regions not in bigwig.
    cut -f 1,2 ${params.genome_fai} | awk '{print \$1 "\t0\t" \$2}' | sort -k1,1 -k2,2n > chromsizes.bed

    if [[ ${binsize} -eq 1 ]]; then
      scale_command=""
    else
      scale_command=\$(awk "BEGIN { printf \"scale %.6f bin %d\", 1/\$binsize, \$binsize }")
    fi

    if [[ ${threshold} == gte* || ${threshold} == lt* ]]; then
      threshold_command=\$(echo ${threshold} | sed -E 's/^(gte|lt)([0-9.]+)/\\1 \\2/')
    else
      echo "ERROR: Unknown threshold type: ${threshold}"
      exit 1
    fi
    
    ${params.wiggletools_bin} \$threshold_command trim chromsizes.bed fillIn chromsizes.bed \$scale_command ${region_filter_file} \
      | ${params.wigToBigWig_bin} stdin <(cut -f 1,2 ${params.genome_fai}) ${region_filter_file}.bin${binsize}.${threshold}.bw
    """
}

/*
  extractVariantsChunk: Run extractVariants.R for an analysis chunk
*/
process extractVariantsChunk {
    cpus 2
    memory '128 GB'
    time '4h'
    tag { "Extract Variants: ${sample_id} -> chunk ${chunkID}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(sample_id), path(bamFile), path(pbiFile), path(baiFile), val(chunkID)
    
    output:
      tuple val(sample_id), path("${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.qs2"), val(chunkID)

    storeDir "${extractVariants_output_dir}"

    script:
    """
    extractVariants.R -c ${params.paramsFileName} -b ${bamFile}
    """
}


/*****************************************************************
 * Individual Workflow Definitions
 *****************************************************************/

workflow processReads {

    main:
    // create a channel of runs
    runs_ch = Channel.fromList(params.runs)

    // Create channel for the input reads file.
    reads_ch = runs_ch.map { run ->
        def reads_file = file(run.reads_file)
        def pbi_file = file("${run.reads_file}.pbi")
        return tuple(run.run_id, reads_file, pbi_file)
    }
    
    // Create barcodes FASTA for each run
    barcodes_ch = runs_ch.map { run ->
        def barcodeFastaContent = run.samples
            .collect { sample ->
              return ">${sample.barcode_id}\n${sample.barcode}"
            }
            .join("\n")
        return tuple(run.run_id, barcodeFastaContent)
    }
    
    makeBarcodesFasta( barcodes_ch )
      
    // Branch according to data type.
    if( params.reads_type == 'subreads' ) {
        // Run CCS in chunks.
        chunkIDs = Channel.of(1..params.ccs_chunks)

        ccsChunk( reads_ch | combine(chunkIDs) | map { it -> tuple(it[0], it[1], it[2], it[3]) } )
        
        // Merge all CCS chunks per run
        ccs_grouped_ch = ccsChunk.out.bampbi_tuple
            .groupTuple(by: 0) // Group by run_id
            .map { run_id, bamFiles, pbiFiles ->
              return tuple(run_id, bamFiles, pbiFiles)
            }

        mergeCCSchunks(ccs_grouped_ch)
        
        // Filter for reads with adapters on both ends.
        filterAdapter( mergeCCSchunks.out )
    }
    else if( params.reads_type == 'ccs' ) {
        // Filter for reads with adapters on both ends.
        filterAdapter( reads_ch )
    }
    else {
        error "Unsupported reads_type '${params.reads_type}'."
    }

    // Join filterAdapter and makeBarcodesFasta outputs by run_id
    limaDemux_input_ch = filterAdapter.out
        .join( makeBarcodesFasta.out, by:0 )
        .map { r_id, bamFile, pbiFile, barcodesFasta ->
          tuple(r_id, bamFile, pbiFile, barcodesFasta)
        }

    // Demultiplex with lima.
    limaDemux( limaDemux_input_ch )
    
    // Create channels for each demultiplexed sample.
    // Here, for each sample the input BAM is assumed to have the name:
    // ${run_id}.ccs.filtered.demux.${barcodeID}--${barcodeID}.bam
    demuxMap_ch = limaDemux.out.bam
    .transpose()
    .map { run_id, bamFile ->
        def m = bamFile.name =~ /${run_id}\.ccs\.filtered\.demux\.(\w+)--\1\.bam/
        if (!m) {
            error "Can't match BAM file name to run_id and barcode_id: ${bamFile.name}"
        }
        return tuple(run_id, m[0][1], bamFile)
      }

    // Create samples to align channel by matching run samples with demux BAMs
    samples_to_align_ch = Channel.fromList(params.runs)
        .flatMap { run ->
            run.samples.collect { sample ->
                return tuple(run.run_id, sample.sample_id, sample.barcode_id)
            }
        }
        .combine(demuxMap_ch, by: 0) // Combine by run_id
        .filter { run_id, sample_id, barcode_id, demux_barcode_id, demux_bam ->
            return barcode_id == demux_barcode_id
        }
        .map { run_id, sample_id, barcode_id, demux_barcode_id, demux_bam ->
            return tuple(run_id, sample_id, demux_bam)
        }

    // Run pbmm2 alignment for each sample.
    pbmm2Align( samples_to_align_ch )

    // Create channels for counting ZMWs of BAMs created during processing
    if( params.reads_type == 'subreads' ) {
      countZMWs_ch = reads_ch.map { f -> tuple(f[1], f[2], "zmwcount.txt") }
        .concat(mergeCCSchunks.out.map { f -> tuple(f[1], f[2], "zmwcount.txt") })
        .collect(flat: false)
    }
    else if( params.reads_type == 'ccs' ) {
      countZMWs_ch = reads_ch.map { f -> tuple(f[1], f[2], "zmwcount.txt") }
        .collect(flat: false)
    }

    countZMWs_ch = countZMWs_ch +
          (
          filterAdapter.out.map { f -> tuple(f[1], f[2], "zmwcount.txt") }
            .concat(
              samples_to_align_ch.map { f -> tuple(f[2], file("${f[2]}.pbi"), "zmwcount.txt") },
              pbmm2Align.out.collect(flat: false).flatMap().map { f -> tuple(f[2], f[3], "zmwcount.txt") }
            )
            .collect(flat: false)
          )

    countZMWs( countZMWs_ch.flatMap() )

    // Group pbmm2Align outputs by sample_id for merging
    pbmm2_grouped_ch = pbmm2Align.out
      .map { run_id, sample_id, bamFile, pbiFile ->
          return tuple(sample_id, bamFile, pbiFile)
      }
      .groupTuple(by: 0) // Group by sample_id

    // Merge BAM files by sample_id across runs
    mergeAlignedSampleBAMs(pbmm2_grouped_ch)

    emit:
    mergeAlignedSampleBAMs.out

}


workflow splitBAMs {

    take:
    alignedSamples_ch
    
    main:
    chunkIDs = Channel.of(1..params.analysis_chunks)

    splitBAM( alignedSamples_ch | combine(chunkIDs) | map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) } )

    emit:
    splitBAM.out

}


workflow prepareFilters {
    
    main:
    // Run installBSgenome
    installBSgenome()

    // Create channel for individual_ids for processGermlineVCFs
    individual_ids_ch = Channel
      .from(params.individuals)
      .map { individual -> individual.individual_id }

    // Run processGermlineVCFs
    processGermlineVCFs( individual_ids_ch )

    // Create channel for germline BAMs
    germlinebams_ch = Channel.fromList(params.individuals)
      .map { run ->
        def germline_bam_file = file(run.germline_bam_file)
        def germline_bam_type = file(run.germline_bam_type)
        return tuple(run.germline_bam_file, germline_bam_type)
      }

    // Run processGermlineBAMs
    processGermlineBAMs( germlinebams_ch )

    // Create channel for the region filters
    region_filters_ch = Channel.fromList(params.region_filters)
        .flatMap { region_filter ->
          def filters=[]

          region_filter.read_filters?.each { filter ->
            filters << tuple(filter.region_filter_file, filter.binsize, filter.threshold)
          }

          region_filter.bin_filters?.each { filter ->
            filters << tuple(filter.region_filter_file, filter.binsize, filter.threshold)
          }

          return filters
        }

    // Run prepareRegionFilters
    prepareRegionFilters( region_filters_ch )

    // Create a completion signal by collecting all outputs
    done_ch = Channel.mix(
            installBSgenome.out,
            processGermlineVCFs.out,
            processGermlineBAMs.out,
            prepareRegionFilters.out
          )
        .collect()
        .map { true }

    emit:
    done_ch

}


workflow extractVariants {

    take:
    splitBAMs_ch
    prepareFilters_done
    
    main:
    // Create dependency on prepareFilters_done and then remove it from the input
    splitBAMs_ch = splitBAMs_ch
        .combine(prepareFilters_done)
        .map { splitBAMs_tuple, _ -> splitBAMs_tuple }

    extractVariantsChunk( splitBAMs_ch )

    emit:
    extractVariantsChunk.out

}

/*****************************************************************
 * Main Workflow
 *****************************************************************/
// --workflow options: all, processReads, splitBAMs, prepareFilters, extractVariants

import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.DumperOptions

//Define output directories
logs_output_dir="${params.analysis_output_dir}/logs"
processReads_output_dir="${params.analysis_output_dir}/processReads"
splitBAMs_output_dir="${params.analysis_output_dir}/splitBAMs"
extractVariants_output_dir="${params.analysis_output_dir}/extractVariants"

workflow {

  // Save copy of parameters file to logs directory
  def logsDir = file("${logs_output_dir}")
  logsDir.mkdirs()

  def options = new DumperOptions()
  options.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK)
  options.setIndent(2)

  def yaml = new Yaml(options)
  def timestamp = new Date().format('yyyy_MMdd_HHmm')
  file("${logsDir}/runParams.${timestamp}.yaml").text = yaml.dump(params)

  //Get path of parameters file
  def commandLineTokens = workflow.commandLine.tokenize()
  params.paramsFileName = commandLineTokens[commandLineTokens.indexOf('-params-file') + 1]

  // Run processReads workflow
  if( params.workflow == "all" || params.workflow == "processReads" ){
    processReads()
  }

  // Run splitBAMs workflow
  if( params.workflow == "all" ){
    splitBAMs( processReads.out )
  }
  else if ( params.workflow == "splitBAMs" ){
    alignedSamples_ch = Channel.fromList(params.samples)
          .map { sample ->
              def sample_id = sample.sample_id
              def bamFile = file("${processReads_output_dir}/${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.bam")
              def pbiFile = file("${processReads_output_dir}/${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.pbi")
              def baiFile = file("${processReads_output_dir}/${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.bai")
              return tuple(sample_id, bamFile, pbiFile, baiFile)
          }
    splitBAMs( alignedSamples_ch )
  }

  // Run prepareFilters workflow
  if( params.workflow == "all" || params.workflow == "prepareFilters"){
    prepareFilters()
  }

  // Run extractVariants workflow, with dependency on prepareFilters
  if( params.workflow=="all" ){
    extractVariants( splitBAMs.out, prepareFilters.out )
  }
  else if ( params.workflow == "extractVariants" ){
    splitBAMs_ch = Channel.fromList(params.samples)
          .combine(Channel.of(1..params.analysis_chunks))
          .map { sample, chunkID ->
              def sample_id = sample.sample_id
              def bamFile = file("${splitBAMs_output_dir}/${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam")
              def pbiFile = file("${splitBAMs_output_dir}/${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam.pbi")
              def baiFile = file("${splitBAMs_output_dir}/${params.analysis_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam.bai")
              return tuple(sample_id, bamFile, pbiFile, baiFile, chunkID)
          }
    extractVariants( splitBAMs_ch, Channel.value(true) ) // 'true' parameter assumes prepareFilters was already run
  }

}