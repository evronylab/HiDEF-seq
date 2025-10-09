/*****************************************************************
 * Global variables and functions
 *****************************************************************/

//Library imports
import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.DumperOptions

//Define output directories and related helper functions
sharedLogsDir = "${params.analysis_output_dir}/${params.analysis_id}.sharedLogs"

def sampleBaseDir = { individual_id, sample_id -> "${params.analysis_output_dir}/${params.analysis_id}.${individual_id}.${sample_id}" }
def dirSampleLogs = { individual_id, sample_id -> "${sampleBaseDir(individual_id, sample_id)}/logs" }
def dirProcessReads = { individual_id, sample_id -> "${sampleBaseDir(individual_id, sample_id)}/processReads" }
def dirSplitBAMs = { individual_id, sample_id -> "${sampleBaseDir(individual_id, sample_id)}/splitBAMs" }
def dirExtractCalls = { individual_id, sample_id -> "${sampleBaseDir(individual_id, sample_id)}/extractCalls" }
def dirFilterCalls = { individual_id, sample_id -> "${sampleBaseDir(individual_id, sample_id)}/filterCalls" }
def dirCalculateBurdens = { individual_id, sample_id -> "${sampleBaseDir(individual_id, sample_id)}/calculateBurdens" }
def dirOutputResults = { individual_id, sample_id -> "${sampleBaseDir(individual_id, sample_id)}/outputResults" }

//Script to save nextflow process logs upon completion of each process
process {
    afterScript = '''
        if [[ -f ".command.log" ]]; then
            # Use environment variables or defaults
            FINAL_LOG_DIR="${NF_LOG_DIR:-''' + "${sharedLogsDir}" + '''}"
            FINAL_LOG_NAME="${NF_LOG_NAME:-${task.process}.command.log}"
            
            mkdir -p "${FINAL_LOG_DIR}"
            cp ".command.log" "${FINAL_LOG_DIR}/${FINAL_LOG_NAME}"
        fi
    '''
    
    // Set environment variables for log directories
    withName: 'processGermlineVCFs|processGermlineBAMs|prepareRegionFilters|ccsChunk' {
      env.NF_LOG_DIR = "${sharedLogsDir}"
    }
    withName: 'pbmm2Align|mergeAlignedSampleBAMs|splitBAM|extractCallsChunk|filterCallsChunk|calculateBurdensChromgroupFiltergroup|outputResultsSample|removeIntermediateFilesProcess' {
      env.NF_LOG_DIR = { dirSampleLogs(individual_id, sample_id) }
    }

    // Set environment variables for log file names
    withName: 'processGermlineVCFs' {
        env.NF_LOG_NAME = '${individual_id}.${task.process}.command.log'
    }
    
    withName: 'processGermlineBAMs' {
        env.NF_LOG_NAME = '${germline_bam_file}.${task.process}.command.log'
    }
    
    withName: 'prepareRegionFilters' {
        env.NF_LOG_NAME = '${task.process}.${region_filter_file}.bin${binsize}.${threshold}.command.log'
    }
    
    withName: 'ccsChunk' {
        env.NF_LOG_NAME = '${task.process}.chunk${chunkID}.command.log'
    }

    withName: 'filterCallsChunk' {
        env.NF_LOG_NAME = '${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.${task.process}.chunk${chunkID}.command.log'
    }
    
    withName: 'calculateBurdensChromgroupFiltergroup' {
        env.NF_LOG_NAME = '${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.${task.process}.command.log'
    }
    
    withName: 'pbmm2Align|mergeAlignedSampleBAMs|splitBAM|extractCallsChunk|outputResultsSample|removeIntermediateFilesProcess' {
        env.NF_LOG_NAME = { chunkID ? "${params.analysis_id}.${individual_id}.${sample_id}.${task.process}.chunk${chunkID}.command.log" : "${params.analysis_id}.${individual_id}.${sample_id}.${task.process}.command.log" }
    }
}

/*****************************************************************
 * Main Workflow
 *****************************************************************/
// --workflow options: all, or one or more contiguous (in the correct order) comma-separated workflows: processReads, splitBAMs, prepareFilters, extractCalls, filterCalls, calculateBurdens, outputResults, removeIntermediateFiles
// Note: the pipeline analyzes all genomic regions up to and including extractCalls. Then it restricts analysis and output files to chromgroups.

workflow {

  // Save copy of parameters file to logs directory
  def logsDir = file("${sharedLogsDir}")
  logsDir.mkdirs()

  def options = new DumperOptions()
  options.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK)
  options.setIndent(2)

  def yaml = new Yaml(options)
  def timestamp = new Date().format('yyyy_MMdd_HHmm')
  file("${logsDir}/runParams.${timestamp}.yaml").text = yaml.dump(params)

  // Save copy of run information
  file("${logsDir}/runInfo.${timestamp}.txt").text = """
  Repository: ${workflow.repository ?: 'N/A'}
  Revision/Tag: ${workflow.revision ?: 'N/A'}
  Commit ID: ${workflow.commitId ?: 'N/A'}
  Run Name: ${workflow.runName}
  Session ID: ${workflow.sessionId}
  Start Time: ${workflow.start}
  Command Line: ${workflow.commandLine}
  Project Directory: ${workflow.projectDir}
  Launch Directory: ${workflow.launchDir}
  Nextflow Version: ${workflow.nextflow.version}
  Nextflow Build: ${workflow.nextflow.build}
  """.stripIndent()

  //Get path of parameters file
  def commandLineTokens = workflow.commandLine.tokenize()
  params.paramsFileName = commandLineTokens[commandLineTokens.indexOf('-params-file') + 1]

  //Extract workflows and check that selected workflows are valid.
  def workflowsToRun = params.workflow.split(',').collect { it.trim() }
  
   //valid workflows in order that they can be run
  def validWorkflows = ["all", "processReads", "splitBAMs", "prepareFilters", "extractCalls", "filterCalls", "calculateBurdens", "outputResults", "removeIntermediateFiles"]
  
  def invalidWorkflows = workflowsToRun.findAll { !(it in validWorkflows) }
  if(invalidWorkflows.size() > 0){
    error "ERROR: Invalid workflow(s) '${invalidWorkflows.join(', ')}'. Valid options are: ${validWorkflows.join(', ')}"
  }

  if( workflowsToRun.contains("all") && workflowsToRun.size() > 1 ) {
    error "ERROR: 'all' cannot be combined with other workflows. Use either 'all' alone or specify individual workflows without 'all'."
  }

  if( !workflowsToRun.contains("all") ) {
    def workflowOrder = validWorkflows.findAll { it != "all" }
    
    // Get indices of workflows in the order specified by the user
    def specifiedIndices = workflowsToRun.collect { workflow ->
      workflowOrder.indexOf(workflow)
    }
    
    // Check if workflows are specified in correct order (indices should be ascending)
    if( specifiedIndices != specifiedIndices.sort(false) ) {
      error "ERROR: Workflows must be specified in the correct order: ${workflowOrder.join(' -> ')}. You specified: ${workflowsToRun.join(', ')}"
    }
    
    // Check if the indices are contiguous (consecutive)
    for (int i = 1; i < specifiedIndices.size(); i++) {
      if( specifiedIndices[i] != specifiedIndices[i-1] + 1 ) {
        error "ERROR: Workflows must be contiguous in this order: ${workflowOrder.join(' -> ')}. You specified: ${workflowsToRun.join(', ')}"
      }
    }
  }

  // Helper function to check if a workflow should be run
  def shouldRun = { workflowName -> 
    workflowsToRun.contains("all") || workflowsToRun.contains(workflowName)
  }

  //Initialize workflow output variables
  def processReads_out = null
  def splitBAMs_out = null
  def prepareFilters_out = null
  def extractCalls_out = null
  def filterCalls_out = null
  def calculateBurdens_out = null
  def outputResults_out = null

  // Run processReads workflow
  if (shouldRun("processReads")) {
    processReads_out = processReads()
  }

  // Run splitBAMs workflow
  if (shouldRun("splitBAMs")) {
    processReads_out = processReads_out ?:
      Channel.fromList(params.samples)
          .map { sample ->
              def individual_id = sample.individual_id
              def sample_id = sample.sample_id
              def bamFile = file("${dirProcessReads(individual_id, sample_id)}/${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam")
              def pbiFile = file("${dirProcessReads(individual_id, sample_id)}/${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.pbi")
              def baiFile = file("${dirProcessReads(individual_id, sample_id)}/${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.bai")
              return tuple(individual_id, sample_id, bamFile, pbiFile, baiFile)
          }

    splitBAMs_out = splitBAMs(processReads_out)
  }

  // Run prepareFilters workflow
  if (shouldRun("prepareFilters")) {
    prepareFilters_out = prepareFilters()
  }

  // Run extractCalls workflow
  if (shouldRun("extractCalls")) {
    splitBAMs_out = splitBAMs_out ?:
      Channel.fromList(params.samples)
          .combine(Channel.from(1..params.analysis_chunks))
          .map { sample, chunkID ->
              def individual_id = sample.individual_id
              def sample_id = sample.sample_id
              def bamFile = file("${dirSplitBAMs(individual_id, sample_id)}/${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam")
              def pbiFile = file("${dirSplitBAMs(individual_id, sample_id)}/${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam.pbi")
              def baiFile = file("${dirSplitBAMs(individual_id, sample_id)}/${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam.bai")
              return tuple(individual_id, sample_id, bamFile, pbiFile, baiFile, chunkID)
          }

    prepareFilters_out = prepareFilters_out ?: Channel.value(true)

    extractCalls_out = extractCalls(splitBAMs_out, prepareFilters_out)
  }

  //Prepare a list with all call_types.analyzein_chromgroups and call_types.SBSindel_call_types.filtergroup configured pairs.
  //Created as a list first instead of a channel to allow static calculation of its size for later downstream use
  chromgroups_filtergroups_list = params.call_types
    .collectMany { call_type ->
      def chromgroup_names
      if (call_type.analyzein_chromgroups == 'all') {
        chromgroup_names = params.chromgroups.collect { it.chromgroup }
      } else {
        chromgroup_names = call_type.analyzein_chromgroups.split(',')
      }

      call_type.SBSindel_call_types.collectMany{ SBSindel_call_type ->
        chromgroup_names.collect{ chromgroup ->
          tuple(chromgroup.trim(), SBSindel_call_type.filtergroup)
        }
      }
    }
    .unique()

  //Create a channel from the chromgroups_filtergroups_list
  chromgroups_filtergroups_ch = Channel.fromList(chromgroups_filtergroups_list)

  // Run filterCalls workflow
  if (shouldRun("filterCalls")) {
    extractCalls_out = extractCalls_out ?:
      Channel.fromList(params.samples)
          .combine(Channel.from(1..params.analysis_chunks))
          .map { sample, chunkID ->
              def individual_id = sample.individual_id
              def sample_id = sample.sample_id
              def extractCallsFile = file("${dirExtractCalls(individual_id, sample_id)}/${params.analysis_id}.${individual_id}.${sample_id}.extractCalls.chunk${chunkID}.qs2")
              return tuple(individual_id, sample_id, extractCallsFile, chunkID)
          }

    filterCalls_out = filterCalls(extractCalls_out.combine(chromgroups_filtergroups_ch))
  }

  // Run calculateBurdens workflow
  if (shouldRun("calculateBurdens")) {
    filterCalls_out = filterCalls_out ?:
      Channel.fromList(params.samples)
          .combine(Channel.from(1..params.analysis_chunks))
          .combine(chromgroups_filtergroups_ch)
          .map { sample, chunkID, chromgroup, filtergroup ->
              def individual_id = sample.individual_id
              def sample_id = sample.sample_id
              def filterCallsFile = file("${dirFilterCalls(individual_id, sample_id)}/${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.filterCalls.chunk${chunkID}.qs2")
              return tuple(individual_id, sample_id, chromgroup, filtergroup, chunkID, filterCallsFile)
          }

    def filterCalls_grouped_ch = filterCalls_out
        .groupTuple(by: [0, 1, 2, 3], size: params.analysis_chunks) // Group by individual_id, sample_id, chromgroup, filtergroup. 'size' parameter emits as soon as each group's chunks finish.
        .map { individual_id, sample_id, chromgroup, filtergroup, chunkIDs, filterCallsFiles ->
            // Sort by chunkID
            def sortedIndices = (0..<chunkIDs.size()).toList().sort { i -> chunkIDs[i] as int }
            def sortedfilterCallsFiles = sortedIndices.collect { filterCallsFiles[it] }
            return tuple(individual_id, sample_id, chromgroup, filtergroup, sortedfilterCallsFiles)
        }

    calculateBurdens_out = calculateBurdens(filterCalls_grouped_ch)
  }

  // Run outputResults workflow
  if (shouldRun("outputResults")) {
    calculateBurdens_out = calculateBurdens_out ?:
      Channel.fromList(params.samples)
          .combine(chromgroups_filtergroups_ch)
          .map { sample, chromgroup, filtergroup ->
              def individual_id = sample.individual_id
              def sample_id = sample.sample_id
              def calculateBurdensFile = file("${dirCalculateBurdens(individual_id, sample_id)}/${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.calculateBurdens.qs2")
              return tuple(individual_id, sample_id, chromgroup, filtergroup, calculateBurdensFile)
          }

    def calculateBurdens_grouped_ch = calculateBurdens_out
        .map { individual_id, sample_id, chromgroup, filtergroup, calculateBurdensFile ->
            return tuple(individual_id, sample_id, calculateBurdensFile)
        }
        .groupTuple(by: [0, 1], size: chromgroups_filtergroups_list.size()) // Group by individual_id, sample_id. Emit as soon as each sample's chromgroup/filtergroup analyses finish.

    outputResults_out = outputResults(calculateBurdens_grouped_ch).out_ch
  }

  // Run removeIntermediateFiles workflow (conditional)
  if (shouldRun("removeIntermediateFiles") && params.remove_intermediate_files) {
    outputResults_out = outputResults_out ?:
      Channel.fromList(params.samples).map { s -> tuple(s.individual_id, s.sample_id) }

    removeIntermediateFiles(outputResults_out)
  }

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
        chunkIDs = Channel.from(1..params.ccs_chunks)

        ccsChunk( reads_ch | combine(chunkIDs) | map { it -> tuple(it[0], it[1], it[2], it[3]) } )
        
        // Merge all CCS chunks per run
        ccs_grouped_ch = ccsChunk.out.bampbi_tuple
            .groupTuple(by: 0) // Group by run_id
            .map { run_id, bamFiles, pbiFiles, chunkIDs ->
              // Sort by chunkID
              def sortedIndices = (0..<chunkIDs.size()).toList().sort { i -> chunkIDs[i] as int }
              def sortedBamFiles = sortedIndices.collect { bamFiles[it] }
              def sortedPbiFiles = sortedIndices.collect { pbiFiles[it] }
              return tuple(run_id, sortedBamFiles, sortedPbiFiles)
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
    def sample_to_individual = params.samples.collectEntries { [ (it.sample_id): it.individual_id ] }

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
            return tuple(run_id, sample_to_individual[sample_id], sample_id, demux_bam)
        }

    // Run pbmm2 alignment for each sample.
    pbmm2Align( samples_to_align_ch )

    // Create channels for counting ZMWs of BAMs created during processing
    if( params.reads_type == 'subreads' ) {
      countZMWs_ch = reads_ch.map { f -> tuple(f[1], f[2], "zmwcount.txt") }
        .mix(mergeCCSchunks.out.map { f -> tuple(f[1], f[2], "zmwcount.txt") })
    }
    else if( params.reads_type == 'ccs' ) {
      countZMWs_ch = reads_ch.map { f -> tuple(f[1], f[2], "zmwcount.txt") }
    }

    countZMWs_ch = countZMWs_ch.mix(
        filterAdapter.out.map { f -> tuple(f[1], f[2], "zmwcount.txt") },
        samples_to_align_ch.map { run_id, individual_id, sample_id, demux_bam -> tuple(demux_bam, file("${demux_bam}.pbi"), "zmwcount.txt") },
        pbmm2Align.out.map { f -> tuple(f[3], f[4], "zmwcount.txt") }
      )

    // Count ZMWs
    countZMWs( countZMWs_ch )

    // Group pbmm2Align outputs by sample_id for merging
    pbmm2_grouped_ch = pbmm2Align.out
      .map { run_id, individual_id, sample_id, bamFile, pbiFile ->
          return tuple(individual_id, sample_id, bamFile, pbiFile)
      }
      .groupTuple(by: [0, 1]) // Group by individual_id, sample_id

    // Merge BAM files by sample_id across runs
    mergeAlignedSampleBAMs(pbmm2_grouped_ch)

    emit:
    mergeAlignedSampleBAMs.out

}


workflow splitBAMs {

    take:
    alignedSamples_ch
    
    main:
    chunkIDs = Channel.from(1..params.analysis_chunks)

    splitBAM( alignedSamples_ch | combine(chunkIDs) | map { it -> tuple(it[0], it[1], it[2], it[3], it[4], it[5]) } )

    emit:
    splitBAM.out

}


workflow prepareFilters {
    
    main:
    // Run installBSgenome
    installBSgenome()

    // Run extractGenomeTrinucleotides
    extractGenomeTrinucleotides()

    // Create channel for individual_ids for processGermlineVCFs
    individual_ids_ch = Channel
      .from(params.individuals)
      .map { individual -> individual.individual_id }

    // Run processGermlineVCFs
    processGermlineVCFs( individual_ids_ch )

    // Create channel for germline BAMs
    germlinebams_ch = Channel.fromList(params.individuals)
      .map { run ->
        return tuple( file(run.germline_bam_file), run.germline_bam_type )
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

          region_filter.genome_filters?.each { filter ->
            filters << tuple(filter.region_filter_file, filter.binsize, filter.threshold)
          }

          return filters
        }
        .unique() // Avoid preparing the same region filter/binsize/threshold configuration twice

    // Run prepareRegionFilters
    prepareRegionFilters( region_filters_ch )

    // Create a completion signal by collecting all outputs
    done_ch = installBSgenome.out
      .mix(
        extractGenomeTrinucleotides.out,
        processGermlineVCFs.out,
        processGermlineBAMs.out,
        prepareRegionFilters.out
      )
      .collect()
      .map { true }

    emit:
    done_ch

}


workflow extractCalls {

    take:
    splitBAMs_ch
    prepareFilters_done
    
    main:
    // Create dependency on prepareFilters_done and then remove it from the input (-2 is the second to last index)
    splitBAMs_ch = splitBAMs_ch
        .combine(prepareFilters_done)
        .map { it.dropRight(1) }

    extractCallsChunk( splitBAMs_ch )

    emit:
    extractCallsChunk.out

}


workflow filterCalls {

    take:
    extractCalls_ch

    main:

    filterCallsChunk( extractCalls_ch )

    emit:
    filterCallsChunk.out

}


workflow calculateBurdens {

    take:
    filterCalls_grouped_ch

    main:
    calculateBurdensChromgroupFiltergroup( filterCalls_grouped_ch )

    emit:
    calculateBurdensChromgroupFiltergroup.out

}

workflow outputResults {

    take:
    calculateBurdens_out

    main:
    outputResultsSample(calculateBurdens_out)

    emit:
    out_ch = outputResultsSample.out.out_ch

}

workflow removeIntermediateFiles {

    take:
    outputResults_out

    main:
    removeIntermediateFilesProcess(outputResults_out)

    emit:
    removeIntermediateFilesProcess.out

}


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
    tag "makeBarcodesFasta"
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
    tag { "ccsChunk: chunk ${chunkID}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(run_id), path(bamFile), path(pbiFile), val(chunkID)

    output:
      tuple val(run_id), path("hifi_reads/${run_id}.chunk${chunkID}.hifi_reads.ccs.bam"), path("hifi_reads/${run_id}.chunk${chunkID}.hifi_reads.ccs.bam.pbi"), val(chunkID), emit: bampbi_tuple
      path "statistics/*.ccs_report.*", emit: report
      path "statistics/*.summary.json", emit: summary

    publishDir "${sharedLogsDir}", mode: "copy", pattern: "statistics/*.ccs_report.*", saveAs: { filename -> new File(filename).getName() }
    publishDir "${sharedLogsDir}", mode: "copy", pattern: "statistics/*.summary.json", saveAs: { filename -> new File(filename).getName() }

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
    tag "mergeCCSchunks"
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
    tag "filterAdapter"
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
    tag "limaDemux"
    container "${params.hidefseq_container}"
    
    input:
      tuple val(run_id), path(bamFile), path(pbiFile), path(barcodesFasta)

    output:
      tuple val(run_id), path("${run_id}.ccs.filtered.demux.*.bam"), emit: bam
      path "*.lima.summary", emit: lima_summary
      path "*.lima.counts", emit: lima_counts
    
    publishDir "${sharedLogsDir}", mode: "copy", pattern: "*.lima.summary"
    publishDir "${sharedLogsDir}", mode: "copy", pattern: "*.lima.counts"
    
    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    lima --ccs --split-named --same --min-score ${params.lima_min_score} \\
         --min-ref-span 0.9 --min-scoring-regions 2 \\
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
    tag { "pbmm2Align: ${sample_id}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(run_id), val(individual_id), val(sample_id), path(bamFile)
    
    output:
      tuple val(run_id), val(individual_id), val(sample_id), path("${run_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.bam"), path("${run_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.bam.pbi")
    
    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbmm2 align -j 8 --preset CCS ${params.genome_mmi} ${bamFile} ${run_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.bam
    pbindex ${run_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.bam
    """
}

/*
  mergeAlignedSampleBAMs: Merges aligned BAM files from the same sample across different runs.
*/
process mergeAlignedSampleBAMs {
    cpus 4
    memory '64 GB'
    time '6h'
    tag { "mergeAlignedSampleBAMs: ${sample_id}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(individual_id), val(sample_id), path(bamFiles), path(pbiFiles)

    output:
      tuple val(individual_id), val(sample_id),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam"),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.pbi"),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.bai")

    storeDir { dirProcessReads(individual_id, sample_id) }

    script:
    """
    sample_basename=${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned

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
    tag { "countZMWs: ${bamFile}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple path(bamFile), path(pbiFile), val(outFileSuffix)
    
    output:
      path "*"
    
    publishDir "${sharedLogsDir}", mode: "copy"
    
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
    tag { "splitBAM: ${sample_id}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(individual_id), val(sample_id), path(bamFile), path(pbiFile), path(baiFile), val(chunkID)
    
    output:
      tuple val(individual_id), val(sample_id),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam"),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam.pbi"),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam.bai"),
      val(chunkID)

    storeDir { dirSplitBAMs(individual_id, sample_id) }

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
    time '8h'
    tag { "installBSgenome" }
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
  extractGenomeTrinucleotides: Extracts trinucleotides for every base in the genome
*/
process extractGenomeTrinucleotides {
    cpus 2
    memory '8 GB'
    time '6h'
    tag { "extractGenomeTrinucleotides" }
    container "${params.hidefseq_container}"

    output:
      path("${params.BSgenome.BSgenome_name}.bed.gz")
      path("${params.BSgenome.BSgenome_name}.bed.gz.tbi")

    storeDir "${params.cache_dir}"

    script:
    """
    #Extract sequences for all bases (except contig edges) from genome with seqkit,
    #convert to upper case, convert to BED format (column 2 of trinucleotides is start position of trinucleotide position),
    #and bgzip + tabix index
    ${params.seqkit_bin} sliding -S '' -s1 -W3 ${params.genome_fasta} | \
      ${params.seqkit_bin} seq -u | \
      ${params.seqkit_bin} fx2tab -Q | \
      awk -F '[:\\-\\t]' 'BEGIN {OFS="\\t"}{print \$1, \$2, \$2+1, \$4}' | \
      ${params.bgzip_bin} -c > ${params.BSgenome.BSgenome_name}.bed.gz

    ${params.tabix_bin} -@ 2 -s 1 -b 2 -e 3 ${params.BSgenome.BSgenome_name}.bed.gz
    """
}

/*
  processGermlineVCFs: Run processGermlineVCFs.R
*/
process processGermlineVCFs {
    cpus 1
    memory '16 GB'
    time '4h'
    tag { "processGermlineVCFs: ${individual_id}" }
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
    cpus 2
    memory '16 GB'
    time '24h'
    tag { "processGermlineBAMs: ${germline_bam_file}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple path(germline_bam_file), val(germline_bam_type)

    output:
      path("${germline_bam_file}.bw")
      path("${germline_bam_file}.vcf.gz*")

    storeDir "${params.cache_dir}"

    script:
    """
    #Output per-base coverage using samtools mpileup and direct BAM variant calls using bcftools mpileup
    #Use similar filters for both samtools and bcftools to ensure that samtools coverage data for calculating
    #the fraction of the genome that was filtered maintains correct calculation of the mutation rate.

    #Using samtools mpileup for genome coverage filtering, because the output can be re-formatted into bedgraph.
    #Using bcftools mpileup for call filtering, because the output is in VCF format that is easier to parse.
    #The difference between samtools mpileup and bcftools mpileup should not be significant.
    #We are not calling indels with bcftools mpileup, since that is too noisy to use for filtering.

    #Slightly different parameters are used for Illumina vs PacBio germline BAM to match how bcftools mpileup is run later
    #in the call filtering analysis.
    
    set -euo pipefail

    if [[ ${germline_bam_type} == Illumina ]]; then
      ${params.samtools_bin} mpileup -A -B -Q 11 -d 999999 --ff 3328 -f ${params.genome_fasta} ${germline_bam_file} 2>/dev/null | awk '{print \$1 "\t" \$2-1 "\t" \$2 "\t" \$4}' > mpileup.bg
      ${params.bcftools_bin} mpileup -A -B -Q 11 -d 999999 --ns 3328 -I -a "INFO/AD" -f ${params.genome_fasta} -Oz ${germline_bam_file} 2>/dev/null > ${germline_bam_file}.vcf.gz
      ${params.bcftools_bin} index -t ${germline_bam_file}.vcf.gz
    elif [[ ${germline_bam_type} == PacBio ]]; then
      ${params.samtools_bin} mpileup -A -B -Q 5 -d 999999 --ff 3328 -f ${params.genome_fasta} ${germline_bam_file} 2>/dev/null | awk '{print \$1 "\t" \$2-1 "\t" \$2 "\t" \$4}' > mpileup.bg
      ${params.bcftools_bin} mpileup -A -B -Q 5 -d 999999 --ns 3328 -I -a "INFO/AD" --max-BQ 50 -F0.1 -o25 -e1 -f ${params.genome_fasta} -Oz ${germline_bam_file} 2>/dev/null > ${germline_bam_file}.vcf.gz
      ${params.bcftools_bin} index -t ${germline_bam_file}.vcf.gz
    else
      echo "ERROR: Unknown germline_bam_type: ${germline_bam_type}"
      exit 1
    fi

    sort --parallel=2 -k1,1 -k2,2n mpileup.bg > mpileup.sorted.bg

    ${params.bedGraphToBigWig_bin} mpileup.sorted.bg <(cut -f 1,2 ${params.genome_fai}) ${germline_bam_file}.bw
    """
}

/*
  prepareRegionFilters
*/
process prepareRegionFilters {
    cpus 2
    memory '64 GB'
    time '24h'
    tag { "prepareRegionFilters: ${region_filter_file}, bin ${binsize}, threshold ${threshold}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple path(region_filter_file), val(binsize), val(threshold)

    output:
      path("${region_filter_file}.bin${binsize}.${threshold}.bw")

    storeDir "${params.cache_dir}"

    script:
    """
    #Make genome BED file to use to fill in zero values for regions not in bigwig.
    awk '{print \$1 "\t0\t" \$2}' ${params.genome_fai} | sort -k1,1 -k2,2n > chromsizes.bed

    if [[ ${binsize} -eq 1 ]]; then
      scale_command=""
    else
      scale_command=\$(awk "BEGIN { printf \\"scale %.6f bin %d\\", 1/${binsize}, ${binsize} }")
    fi

    echo "scale command: \$scale_command"

    if [[ ${threshold} == gt* || ${threshold} == lt* ]]; then
      threshold_command=\$(echo ${threshold} | sed -E 's/^(gt|gte|lt|lte)([0-9.]+)/\\1 \\2/')
    else
      echo "ERROR: Unknown threshold type: ${threshold}"
      exit 1
    fi

    echo "threshold command: \$threshold_command"
    
    ${params.wiggletools_bin} \$threshold_command trim chromsizes.bed fillIn chromsizes.bed \$scale_command ${region_filter_file} \
      | ${params.wigToBigWig_bin} stdin <(cut -f 1,2 ${params.genome_fai}) ${region_filter_file}.bin${binsize}.${threshold}.bw
    """
}

/*
  extractCallsChunk: Run extractCalls.R for an analysis chunk
*/
process extractCallsChunk {
    cpus 1
    memory { 
      def baseMemory = params.mem_extractCallsChunk as nextflow.util.MemoryUnit
      baseMemory * (1 + 0.5*(task.attempt - 1))
    }
    time { 
        def baseTime = params.time_extractCallsChunk as nextflow.util.Duration
        baseTime * (1 + (task.attempt - 1))
    }
    maxRetries params.maxRetries_extractCallsChunk
    tag { "extractCallsChunk: ${sample_id} -> chunk ${chunkID}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(individual_id), val(sample_id), path(bamFile), path(pbiFile), path(baiFile), val(chunkID)
    
    output:
      tuple val(individual_id), val(sample_id), path("${params.analysis_id}.${individual_id}.${sample_id}.extractCalls.chunk${chunkID}.qs2"), val(chunkID)

    storeDir { dirExtractCalls(individual_id, sample_id) }

    script:
    """
    extractCalls.R -c ${params.paramsFileName} -b ${bamFile} -o ${params.analysis_id}.${individual_id}.${sample_id}.extractCalls.chunk${chunkID}.qs2
    """
}

/*
  filterCallsChunk: Run filterCalls.R for each analysis chunk, chromgroup, filtergroup combination
*/
process filterCallsChunk {
    cpus 1
    memory { 
      def baseMemory = params.mem_filterCallsChunk as nextflow.util.MemoryUnit
      baseMemory * (1 + 0.5*(task.attempt - 1))
    }
    time { 
        def baseTime = params.time_filterCallsChunk as nextflow.util.Duration
        baseTime * (1 + (task.attempt - 1))
    }
    maxRetries params.maxRetries_filterCallsChunk
    tag { "filterCallsChunk: ${sample_id} -> chunk ${chunkID}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(individual_id), val(sample_id), path(extractCallsFile), val(chunkID), val(chromgroup), val(filtergroup)
    
    output:
      tuple val(individual_id), val(sample_id), val(chromgroup), val(filtergroup), val(chunkID), path("${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.filterCalls.chunk${chunkID}.qs2")

    storeDir { dirFilterCalls(individual_id, sample_id) }

    script:
    """
    filterCalls.R -c ${params.paramsFileName} -s ${sample_id} -g ${chromgroup} -v ${filtergroup} -f ${extractCallsFile} -o ${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.filterCalls.chunk${chunkID}.qs2
    """
}

/*
  calculateBurdensChromgroupFiltergroup: Run calculateBurdens.R for each sample_id x chromgroup x filtergroup combination
*/
process calculateBurdensChromgroupFiltergroup {
    cpus 2
    memory { 
      def baseMemory = params.mem_calculateBurdensChromgroupFiltergroup as nextflow.util.MemoryUnit
      baseMemory * (1 + 0.5*(task.attempt - 1))
    }
    time { 
        def baseTime = params.time_calculateBurdensChromgroupFiltergroup as nextflow.util.Duration
        baseTime * (1 + (task.attempt - 1))
    }
    maxRetries params.maxRetries_calculateBurdensChromgroupFiltergroup
    tag { "calculateBurdensChromgroupFiltergroup: ${sample_id} -> ${chromgroup} x ${filtergroup}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(individual_id), val(sample_id), val(chromgroup), val(filtergroup), path(filterCallsFiles)
    
    output:
      tuple val(individual_id), val(sample_id), val(chromgroup), val(filtergroup), path("${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.calculateBurdens.qs2")

    storeDir { dirCalculateBurdens(individual_id, sample_id) }

    script:
    """
    set -euo pipefail

    calculateBurdens.R -c ${params.paramsFileName} -s ${sample_id} -g ${chromgroup} -v ${filtergroup} -f ${filterCallsFiles.join(',')} -o ${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.calculateBurdens.qs2

    #Move bed.gz[.tbi] files safely (with checks) to output results directory.
    TARGET_DIR=${dirOutputResults(individual_id, sample_id)}/coverage_reftnc/${chromgroup}

    mkdir -p \$TARGET_DIR || {
      echo "ERROR: Failed to create output directory \$TARGET_DIR" >&2
      exit 1
    }

    for file in *.bed.gz*; do
      if [[ -f \$file ]]; then  # Check file exists
        mv \$file \$TARGET_DIR || {
          echo "ERROR: Failed to move file \$file" >&2
          exit 1
        }
      fi
    done
    """
}

/*
  outputResultsSample: Run outputResults.R for each sample_id
*/
process outputResultsSample {
    cpus 1
    memory { 
      def baseMemory = params.mem_outputResultsSample as nextflow.util.MemoryUnit
      baseMemory * (1 + 0.5*(task.attempt - 1))
    }
    time { 
        def baseTime = params.time_outputResultsSample as nextflow.util.Duration
        baseTime * (1 + (task.attempt - 1))
    }
    maxRetries params.maxRetries_outputResultsSample
    tag { "outputResultsSample: ${sample_id}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(individual_id), val(sample_id), path(calculateBurdensFiles)
    
    output:
      tuple val(individual_id), val(sample_id), emit: out_ch
      path("${params.analysis_id}.${individual_id}.${sample_id}.qs2")
      path("${params.analysis_id}.${individual_id}.${sample_id}.yaml_config.tsv")
      path("${params.analysis_id}.${individual_id}.${sample_id}.run_metadata.tsv")
      path("**/*.{tsv,vcf.bgz,vcf.bgz.tbi,pdf}")
        
    storeDir { dirOutputResults(individual_id, sample_id) }

    script:
    """
    outputResults.R -c ${params.paramsFileName} -s ${sample_id} -f ${calculateBurdensFiles.join(',')} -o ${params.analysis_id}.${individual_id}.${sample_id}
    """
}

/*
  removeIntermediateFilesProcess: Removes intermediate files to save disk space after successful completion.
*/
process removeIntermediateFilesProcess {
    cpus 1
    memory '2 GB'
    time '2h'
    tag { "removeIntermediateFilesProcess: ${sample_id}" }
    container "${params.hidefseq_container}"
    
    input:
      tuple val(individual_id), val(sample_id)
    
    output:
      val(true)
    
    script:
     """
    set -euo pipefail
    echo "Starting cleanup of intermediate files for ${params.analysis_id}.${individual_id}.${sample_id}."

    rm -rf "${dirSplitBAMs(individual_id, sample_id)}"
    rm -rf "${dirExtractCalls(individual_id, sample_id)}"
    rm -rf "${dirFilterCalls(individual_id, sample_id)}"
    rm -rf "${dirCalculateBurdens(individual_id, sample_id)}"

    echo "Cleanup completed."
    """
}