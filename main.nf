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
      val content
    
    output:
      path "barcodes.fasta"
    
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
      tuple path(bamFile), path(pbiFile), val(chunkID)

    output:
      tuple path("hifi_reads/${params.run_id}.chunk${chunkID}.hifi_reads.ccs.bam"), path("hifi_reads/${params.run_id}.chunk${chunkID}.hifi_reads.ccs.bam.pbi")
    
    publishDir "${params.processReads_output_dir}/logs", mode: 'copy', pattern: "statistics/*.ccs_report.*"
    publishDir "${params.processReads_output_dir}/logs", mode: 'copy', pattern: "statistics/*.summary.json"

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
        --movie-name ${params.run_id}.chunk${chunkID} \\
        --non-hifi-prefix fail \\
        --report-file statistics/${params.run_id}.chunk${chunkID}.ccs_report.txt \\
        ${bamFile}
    """
}

/*
  mergeCCS: Merges all CCS chunk outputs into a single BAM.
*/
process mergeCCS {
    cpus 2
    memory '8 GB'
    time '6h'
    tag "Merge CCS chunks"
    container "${params.hidefseq_container}"
    
    input:
      tuple path(bamChunks), path(pbiChunks)

    output:
      tuple path("${params.run_id}.ccs.bam"), path("${params.run_id}.ccs.bam.pbi")

    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbmerge -o ${params.run_id}.ccs.bam ${bamChunks.join(' ')}
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
      tuple path(bamFile), path(pbiFile)
    
    output:
      tuple path("${params.run_id}.ccs.filtered.bam"), path("${params.run_id}.ccs.filtered.bam.pbi")
    
    script:
    """
    ${params.samtools_bin} view -b -@8 -e "[ma]==0" ${bamFile} > ${params.run_id}.ccs.filtered.bam
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbindex ${params.run_id}.ccs.filtered.bam
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
      tuple path(bamFile), path(pbiFile)
      path barcodesFasta

    output:
      path "${params.run_id}.ccs.filtered.demux.*.bam", emit: bam
    
    publishDir "${params.processReads_output_dir}/logs", mode: 'copy', pattern: "*.lima.{summary,counts}"
    
    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    lima --ccs --split-named --min-score 80 --min-end-score 50 \\
         --min-ref-span 0.75 --same --min-scoring-regions 2 \\
         ${bamFile} \\
         ${barcodesFasta} \\
         ${params.run_id}.ccs.filtered.demux.bam
    """
}

/*
  pbmm2Align: Aligns a demultiplexed BAM file using pbmm2.
  and renames the output to include the sample name.
*/
process pbmm2Align {
    cpus 8
    memory '64 GB'
    time '12h'
    tag { "pbmm2 Alignment: ${sample_basename}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(sample_basename), val(barcodeID), path(bamFile)
    
    output:
      tuple val(sample_basename), val(barcodeID), path("${sample_basename}.aligned.sorted.bam"), path("${sample_basename}.aligned.sorted.bam.pbi"), path("${sample_basename}.aligned.sorted.bam.bai")
    
    publishDir params.processReads_output_dir, mode: 'copy', pattern: "${sample_basename}.aligned.sorted.bam*"
    
    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbmm2 align -j 8 --preset CCS ${params.genome_mmi} ${bamFile} ${sample_basename}.aligned.bam
    conda deactivate

    ${params.samtools_bin} sort -@8 -m 4G ${sample_basename}.aligned.bam > ${sample_basename}.aligned.sorted.bam
    ${params.samtools_bin} index -@8 ${sample_basename}.aligned.sorted.bam

    conda activate ${params.conda_pbbioconda_env}
    pbindex ${sample_basename}.aligned.sorted.bam
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
    
    publishDir "${params.processReads_output_dir}/logs", mode: 'copy'
    
    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    zmwfilter --show-all ${bamFile} | wc -l > \$(basename ${bamFile} .bam).${outFileSuffix}
    """
}

/*****************************************************************
 * Individual Workflow Definitions
 *****************************************************************/

workflow processReads {

    main:
    // Create channel for the input reads file.
    reads_ch = Channel
      .fromPath(params.reads_file)
      .map{ f -> tuple(f, file("${f}.pbi")) }
    
    // Make barcodes FASTA
    def barcodeFastaContent = params.samples
      .collect {
          sample ->
            def tokens = sample.barcode.tokenize(':')
            return ">${tokens[0]}\n${tokens[1]}"
      }
      .join("\n")
    
    makeBarcodesFasta( Channel.value(barcodeFastaContent) )
      
    // Branch according to data type.
    if( params.data_type == 'subreads' ) {
        // Run CCS in chunks.
        chunkIDs = Channel.of(1..params.ccs_chunks)

        ccsChunk( reads_ch | combine(chunkIDs) | map { it -> tuple(it[0], it[1], it[2]) } )
        
        // Merge all CCS chunks.
        mergeCCS( ccsChunk.out.collect(flat: false).map{ it.transpose() } )
        
        // Filter for reads with adapters on both ends.
        filterAdapter( mergeCCS.out )
    }
    else if( params.data_type == 'ccs' ) {
        // Filter for reads with adapters on both ends.
        filterAdapter( reads_ch )
    }
    else {
        error "Unsupported data_type '${params.data_type}'."
    }

    // Demultiplex with lima.
    limaDemux( filterAdapter.out, makeBarcodesFasta.out )
    
    // Create channels for each demultiplexed sample.
    // Here, for each sample the input BAM is assumed to have the name:
    // ${run_id}.ccs.filtered.demux.${barcodeID}--${barcodeID}.bam
    demuxMap_ch = limaDemux.out.bam
    .collect()
    .map {
      files ->
        def result = [:]
        files.each
        { file ->
            def m = file.name =~ /${params.run_id}\.ccs\.filtered\.demux\.(\w+)--\1\.bam/
            if (m) {
              result[m[0][1]] = file
            }
        }
        return result
      }

    samples_to_align_ch = Channel.fromList(params.samples)
    .combine(demuxMap_ch)
    .map{ sample, demuxMap ->
          def sname = sample.sample_name
          def barcodeID = sample.barcode.tokenize(':')[0]
          def sample_basename = "${params.run_id}.${sname}.ccs.filtered"
          def demux_bam = demuxMap[barcodeID]
          return tuple(sample_basename, barcodeID, demux_bam)
      }

    // Run pbmm2 alignment for each sample.
    pbmm2Align( samples_to_align_ch )

    // Create channels for counting ZMWs of BAMs created during processing
    if( params.data_type == 'subreads' ) {
      countZMWs_ch = reads_ch
        .map { f -> tuple(f[0], f[1], "subreads_zmwcount.txt") }
        .merge(mergeCCS.out | map { f -> tuple(f[0], f[1], "ccs_zmwcount.txt") })
    }
    else if( params.data_type == 'ccs' ) {
      countZMWs_ch = reads_ch | map { f -> tuple(f[0], f[1], "ccs_zmwcount.txt") }
    }

    pbmm2Align.out.collect().subscribe { println "DEBUG: pbmm2Align.out: $it" }

    countZMWs_ch = countZMWs_ch
      .merge(
        filterAdapter.out | map { f -> tuple(f[0], f[1], "filteredAdapter_zmwcount.txt") },
        limaDemux.out.bam.flatten() | map { f -> tuple(f, file("${f}.pbi"), "limaDemux_zmwcount.txt") },
        pbmm2Align.out.collect() | map { f -> tuple(f[2], f[3], "aligned_zmwcount.txt") }
      )
      .map { it.transpose() }

    countZMWs( countZMWs_ch )

    emit:
    pbmm2Align.out
}


/*****************************************************************
 * Main Workflow
 *****************************************************************/

workflow {

  // Run processReads workflow
  if( params.workflow=="all" || params.workflow == "processReads" ){
    processReads()
  }

}