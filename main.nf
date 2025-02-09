/*****************************************************************
 * Process Definitions
 *****************************************************************/

/*
  makeBarcodesFasta: Writes the barcode FASTA content to a file from a value input.
*/
process makeBarcodesFasta {
    cpus 1
    memory '2 GB'
    time '5m'
    tag "Make Barcodes Fasta"
    container '${params.hidef_container}'

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
  countZMWs: Runs zmwfilter on an input file and writes the ZMW count to a file.
*/
process countZMWs {
    cpus 1
    memory '8 GB'
    time '30m'
    tag "Count ZMWs"
    container '${params.hidef_container}'
    
    input:
      tuple path(inFile), val(outName)
    
    output:
      path outName
    
    publishDir "${params.process_reads_output_path}/logs", mode: 'copy'
    
    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    zmwfilter --show-all ${inFile} | wc -l > ${outName}
    """
}

/*
  ccsChunk: Runs one chunk of CCS.

  #Notes:
  #1. Using ccs v8.0 (standalone version of ICS13 ccs, http://doi.org/10.5281/zenodo.10703290 ) that supports rich HiFi tags (--subread-pileup-summary-tags).
  #2. Set '--pbdc' and '--pbdc-skip-min-qv 0' to go into DeepConsensus code path to support rich HiFi tags but skip DeepConsensus polishing (by setting to skip any 100 bp window with average QV >0; i.e. every window)
  #3. Set '--binned-qvs=False' to keep full resolution quality values. Setting pptop-passes to 255 instead of 0 to match what occurs on Revio (since setting it to 0, i.e. unlimited, caused technical issues on Revio).
  #4. Setting --instrument-files-layout --min-rq -1 --movie-name <hifireads.bam> --non-hifi-prefix <failreads.bam> to mimic what happens on Revio, since this outputs min-rq >= 0.99 into hifireads.bam and also removes other types of failed reads (per ff tag details here: https://pacbiofileformats.readthedocs.io/en/13.0/BAM.html#use-of-read-tags-for-fail-per-read-information). We then do independent rq filtering based on pipeline configuration later in the pipeline.
  #5. Off-instrument ccs run with CPUs differs slightly from on-instrument Revio ccs run with GPUs for the parameters --max-insertion-size and window size, due to technical details, but PacBio says this should negligibly affect the output and we should not specify these.
*/
process ccsChunk {
    cpus 8
    memory '64 GB'
    time '24h'
    tag { "CCS chunk ${chunk_id}" }
    container '${params.hidef_container}'
    
    input:
      tuple path(reads_file), val(chunk_id)

    output:
      path "${params.ccs_BAM_prefix}.ccs.chunk${chunk_id}.bam"
    
    publishDir "${params.process_reads_output_path}/logs", mode: 'copy', pattern: "*.ccsreport.*"
    publishDir "${params.process_reads_output_path}/logs", mode: 'copy', pattern: "*.ccsmetrics.*"

    script:
    // Build the LD_PRELOAD command if the parameter is set.
    def ld_preload_cmd = (params.ccs_ld_preload && params.ccs_ld_preload.trim()) ? "export LD_PRELOAD=${params.ccs_ld_preload}" : ""

    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    ${ld_preload_cmd}
    ccs -j 8 --log-level INFO --by-strand --hifi-kinetics --instrument-files-layout --min-rq -1 --top-passes 255 \\
        --pbdc --pbdc-skip-min-qv 0 --subread-pileup-summary-tags --binned-qvs=False \\
        --chunk ${chunk_id}/${params.ccschunks} \\
        --movie-name ${params.ccs_BAM_prefix}.ccs.chunk${chunk_id} \\
        --non-hifi-prefix ${params.ccs_BAM_prefix}.failreads.chunk${chunk_id} \\
        --report-file ${params.ccs_BAM_prefix}.ccsreport.chunk${chunk_id}.txt \\
        --report-json ${params.ccs_BAM_prefix}.ccsreport.chunk${chunk_id}.json \\
        --metrics-json ${params.ccs_BAM_prefix}.ccsmetrics.chunk${chunk_id}.json \\
        ${reads_file}
    """
}

/*
  mergeCCS: Merges all CCS chunk outputs into a single BAM.
*/
process mergeCCS {
    cpus 2
    memory '32 GB'
    time '4h'
    tag "Merge CCS chunks"
    container '${params.hidef_container}'
    
    input:
      path bam_chunks

    output:
      path "${params.ccs_BAM_prefix}.ccs.bam"

    publishDir "${params.process_reads_output_path}/logs", mode: 'copy'

    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbmerge -o ${params.ccs_BAM_prefix}.ccs.bam ${bam_chunks.join(' ')}
    """
}

/*
  filterAdapter: Filter CCS BAM file to keep only reads with ma tag == 0 (adapter detected on both ends).
*/
process filterAdapter {
    cpus 8
    memory '64 GB'
    time '10h'
    tag "Filter Bad Adapters"
    container '${params.hidef_container}'
    
    input:
      path inBam
    
    output:
      path "${params.ccs_BAM_prefix}.ccs.filtered.bam"
    
    script:
    """
    ${params.samtools_bin} view -b -@8 -e "[ma]==0" ${inBam} > ${params.ccs_BAM_prefix}.ccs.filtered.bam
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbindex ${params.ccs_BAM_prefix}.ccs.filtered.bam
    """
}

/*
  limaDemux: Demultiplexes the filtered BAM using lima.
  Assumes that the output for a given barcode will be named as:
      ${ccs_BAM_prefix}.ccs.demux.${barcodeId}--${barcodeId}.bam
  Resource request: 8 cores, 64GB, 12h.
*/
process limaDemux {
    cpus 8
    memory '64 GB'
    time '12h'
    tag "Lima Demultiplexing"
    container '${params.hidef_container}'
    
    input:
      path filteredBam
      path barcodesFasta

    output:
      path "${params.ccs_BAM_prefix}.ccs.filtered.demux.*.bam", emit: bam
      path "${params.ccs_BAM_prefix}.ccs.filtered.demux.*.bam.pbi", emit: bai
    
    publishDir "${params.process_reads_output_path}/logs", mode: 'copy', pattern: "*.lima.*"
    
    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    lima --ccs --same --split-named --min-score 80 --min-end-score 50 \\
         --min-ref-span 0.75 --min-scoring-regions 2 \\
         ${filteredBam} \\
         ${barcodesFasta} \\
         ${params.ccs_BAM_prefix}.ccs.filtered.demux.bam
    """
}

/*
  pbmm2Align: Aligns a demultiplexed BAM file using pbmm2.
  Resource request: 8 cores, 128GB, 30h.
  Here, for each sample the input BAM is assumed to have the name:
      ${ccs_BAM_prefix}.ccs.demux.${barcodeId}--${barcodeId}.bam
  and pbmm2Align renames the output to include the sample name.
*/
process pbmm2Align {
    cpus 8
    memory '64 GB'
    time '12h'
    tag { "pbmm2 Alignment: ${sample_basename}" }
    container '${params.hidef_container}'
    
    input:
      tuple val(sample_basename), file(demuxBam)
    
    output:
      file("${sample_basename}.aligned.bam")
    
    publishDir params.process_reads_output_path, mode: 'copy', pattern: "${sample_basename}.aligned.bam"
    
    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbmm2 align -j 8 --preset CCS \\
    ${params.genome_mmi} \\
    ${demuxBam} \\
    ${sample_basename}.aligned.bam
    """
}

/*****************************************************************
 * Individual Workflow Definitions
 *****************************************************************/

workflow processReads {

    main:

    // Create channel for the input reads file.
    reads_ch = Channel.fromPath(params.reads_filename)
    
    // Make barcodes FASTA
    def barcodeFastaContent = params.samples.collect { sample ->
        def tokens = sample.barcode.tokenize(':')
        return ">${tokens[0]}\n${tokens[1]}"
    }.join("\n")
    
    barcodesFasta = Channel.value(barcodeFastaContent) | makeBarcodesFasta
    
    // Count ZMWs on the original input.
    reads_ch.map { f -> tuple(f, "original_zmwcount.txt") } | countZMWs
    
    // Declare filteredCCS in the outer scope.
    def filteredCCS
    
    // Branch according to data type.
    if( params.data_type == 'subreads' ) {
        
        // Run CCS in chunks.
        chunk_ids = Channel.of(1..params.ccschunks)

        ccsChunks = reads_ch.combine(chunk_ids)
                        .map { r, id -> tuple(r, id) }
                        | ccsChunk
        
        // Merge all CCS chunks.
        mergedCCS = mergeCCS(ccsChunks.collect()).out
        
        // Count ZMWs after CCS merge.
        mergedCCS.map { f -> tuple(f, "ccs_zmwcount.txt") } | countZMWs
        
        // Filter for reads with adapters on both ends.
        filteredCCS = mergedCCS | filterAdapter
    }
    else if( params.data_type == 'ccs' ) {
        // Filter for reads with adapters on both ends.
        filteredCCS = reads_ch | filterAdapter
    }
    else {
        error "Unsupported data_type '${params.data_type}'."
    }

    // Count ZMWs after adapter filtering.
    filteredCCS.map { f -> tuple(f, "filtered_ccs_zmwcount.txt") } | countZMWs
        
    // Demultiplex with lima.
    demuxedCCS = limaDemux(filteredCCS, barcodesFasta)
    
    //Create channels for each demultiplexed sample.
    demuxMap = demuxedCCS.bam.collect().collect { files ->
    def result = [:]
    files.each { file ->
         def m = file.name =~ /${params.ccs_BAM_prefix}\.ccs\.filtered\.demux\.(\w+)--\1\.bam/
         if (m) {
             result[m[0][1]] = file
         }
      }
      return result
    }

    samples_to_align_ch = Channel.fromList(params.samples).map { sample ->
        def sname = sample.sample_name
        def barcodeId = sample.barcode.tokenize(':')[0]
        def sample_basename = "${params.ccs_BAM_prefix}.${sname}.ccs.filtered.demux.${barcodeId}"
        def demux_bam = demuxMap[barcodeId]
        return tuple(sample_basename, demux_bam)
    }

    // Run pbmm2 alignment for each sample.
    alignedCCS = samples_to_align_ch | pbmm2Align
    
    // Count ZMWs after pbmm2 alignment.
    alignedCCS.map { tup -> tuple(tup[1], "aligned_pbmm2_zmwcount.txt") } | countZMWs
}


/*****************************************************************
 * Main Workflow
 *****************************************************************/

workflow {
    processReads()
}