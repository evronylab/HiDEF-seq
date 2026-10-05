/*****************************************************************
 * Global variables and functions
 *****************************************************************/

//Define output directories and related helper functions
def sharedLogsDir() {
  return "${params.analysis_output_dir}/${params.analysis_id}.sharedLogs"
}

def sampleBaseDir(individual_id, sample_id) {
  return "${params.analysis_id}.${individual_id}.${sample_id}"
}

def dirSampleLogs(individual_id, sample_id) {
  return "${sampleBaseDir(individual_id, sample_id)}/logs"
}

def dirProcessReads(individual_id, sample_id) {
  return "${sampleBaseDir(individual_id, sample_id)}/processedReads"
}

def dirVerifyBAMID(individual_id, sample_id) {
  return "${sampleBaseDir(individual_id, sample_id)}/verifyBAMID"
}

def dirSplitBAMs(individual_id, sample_id) {
  return "${sampleBaseDir(individual_id, sample_id)}/splitBAMs"
}

def dirExtractCalls(individual_id, sample_id) {
  return "${sampleBaseDir(individual_id, sample_id)}/extractCalls"
}

def dirFilterCalls(individual_id, sample_id) {
  return "${sampleBaseDir(individual_id, sample_id)}/filterCalls"
}

def dirCalculateBurdens(individual_id, sample_id) {
  return "${sampleBaseDir(individual_id, sample_id)}/calculateBurdens"
}

def dirCoverage_Reftnc(individual_id, sample_id) {
  return "${sampleBaseDir(individual_id, sample_id)}/coverage_reftnc"
}

// Rebuild the original seven-field downstream tuples by exact basename, never
// glob/list order (chunk10 sorts before chunk2, and one chunk is a scalar path).
def expandAnalysisChunks(individual_id, sample_id, bamFiles, pbiFiles, baiFiles, effectiveChunks) {
  def bams = bamFiles instanceof List ? bamFiles : [bamFiles]
  def pbis = pbiFiles instanceof List ? pbiFiles : [pbiFiles]
  def bais = baiFiles instanceof List ? baiFiles : [baiFiles]
  def byName = (bams + pbis + bais).collectEntries { item -> [(item.name): item] }
  if (bams.size() != effectiveChunks || pbis.size() != effectiveChunks || bais.size() != effectiveChunks || byName.size() != 3 * effectiveChunks) {
    throw new IllegalStateException("Unexpected dispatch output set for ${individual_id}/${sample_id}")
  }
  (1..effectiveChunks).collect { chunkID ->
    def basename = "${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk${chunkID}.bam".toString()
    if (!byName.containsKey(basename) || !byName.containsKey(basename + '.pbi') || !byName.containsKey(basename + '.bai')) {
      throw new IllegalStateException("Missing dispatch chunk ${chunkID} for ${individual_id}/${sample_id}")
    }
    tuple(individual_id, sample_id, byName[basename], byName[basename + '.pbi'], byName[basename + '.bai'], chunkID, effectiveChunks)
  }
}

//Function to save nextflow process logs upon completion of each process
def generateAfterScript(logDir, logName) {
  // Strip workflow prefixes like "processReads:" from the log name
  logName = logName.replaceAll(/[\w]+:/,'')

  return """
      if [[ -f ".command.log" ]]; then
          mkdir -p "${logDir}"
          cp ".command.log" "${logDir}/${logName}"
      fi
  """
}

def signatureYaml() {
  def options = new org.yaml.snakeyaml.DumperOptions()
  options.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK)
  options.setPrettyFlow(false)
  options.setIndent(2)

  return new org.yaml.snakeyaml.Yaml(options)
}

// Function to calculate a hash for a subset of keys from a parameters file
def configHash(params_input, keys) {
  def subset = new LinkedHashMap()
  keys.each { key ->
    if (params_input.containsKey(key)) {
      subset[key] = params_input[key]
    }
  }

  def serialized = signatureYaml().dump(subset ?: [:])
  def digest = java.security.MessageDigest.getInstance('SHA-256')
  digest.update(serialized.getBytes('UTF-8'))
  digest.digest().encodeHex().toString()
}

def fileSha256(path) {
  return java.security.MessageDigest.getInstance('SHA-256')
    .digest(file(path).bytes)
    .encodeHex()
    .toString()
}

// Metadata identities avoid reading large BAMs, references or containers on the
// orchestration node. They are intentionally not described as content hashes.
def preparedInputIdentity(rawPath) {
  def source = file(rawPath).toRealPath()
  def attrs = java.nio.file.Files.readAttributes(source, java.nio.file.attribute.BasicFileAttributes)
  def changed = null
  def inode = null
  try {
    changed = java.nio.file.Files.getAttribute(source, 'unix:ctime').toString()
  } catch (UnsupportedOperationException ignored) {
    changed = null
  }
  try {
    // UnixFileKey includes st_dev, which is local to a node's mount table.
    // The same shared file can therefore have different fileKeys on Torch.
    inode = java.nio.file.Files.getAttribute(source, 'unix:ino').toString()
  } catch (UnsupportedOperationException ignored) {
    inode = null
  }
  return [kind: 'canonical-path-size-mtime-ctime-inode', path: source.toString(),
          bytes: attrs.size(), modified: attrs.lastModifiedTime().toString(),
          changed: changed, inode: inode]
}

def canonicalCacheValue(value) {
  if (value instanceof Map) {
    def result = new TreeMap()
    value.each { key, item -> result[key.toString()] = canonicalCacheValue(item) }
    return result
  }
  if (value instanceof Collection) {
    return value.collect { canonicalCacheValue(it) }
  }
  return value instanceof GString ? value.toString() : value
}

// Emit ordinary YAML scalars/collections, never tagged Nextflow/Groovy runtime
// objects such as MemoryUnit. Only original configuration keys and explicitly
// added artifact fields are passed through this conversion.
def plainConfigurationValue(value) {
  if (value == null) return null
  if (value instanceof CharSequence || value instanceof java.nio.file.Path ||
      value instanceof nextflow.util.MemoryUnit || value instanceof nextflow.util.Duration) return value.toString()
  if (value instanceof Number || value instanceof Boolean) return value
  if (value instanceof Map) {
    def result = new LinkedHashMap()
    value.each { key, item -> result[key.toString()] = plainConfigurationValue(item) }
    return result
  }
  if (value instanceof Collection) return value.collect { plainConfigurationValue(it) }
  throw new IllegalArgumentException("Unsupported configuration value type: ${value.getClass().name}")
}

def originalParamsFile(commandLine) {
  def matches = commandLine =~ /(?:^|\s)-params-file(?:\s+|=)("[^"]*"|'[^']*'|\S+)/
  if (!matches.find()) throw new IllegalArgumentException('An original -params-file YAML is required')
  def value = matches.group(1)
  if ((value.startsWith('"') && value.endsWith('"')) || (value.startsWith("'") && value.endsWith("'"))) {
    value = value.substring(1, value.length() - 1)
  }
  return file(value).toAbsolutePath()
}

def writeImmutableConfiguration(path, text) {
  def destination = file(path).toAbsolutePath()
  def expected = text.getBytes('UTF-8')
  def verifyExisting = {
    if (!java.util.Arrays.equals(java.nio.file.Files.readAllBytes(destination), expected)) {
      throw new IllegalStateException("Existing immutable configuration differs: ${destination}")
    }
  }
  if (java.nio.file.Files.exists(destination)) {
    verifyExisting.call()
    return
  }
  def temporary = java.nio.file.Files.createTempFile(destination.parent, '.effective-params-', '.tmp')
  try {
    java.nio.file.Files.write(temporary, expected)
    java.nio.channels.FileChannel.open(temporary, java.nio.file.StandardOpenOption.WRITE).withCloseable { channel ->
      channel.force(true)
    }
    try {
      // Publish one complete sibling file atomically without ever replacing an
      // existing name. Unlike ATOMIC_MOVE, createLink has no replace ambiguity.
      java.nio.file.Files.createLink(destination, temporary)
    } catch (java.nio.file.FileAlreadyExistsException ignored) {
      verifyExisting.call()
    }
  } catch (Exception failure) {
    java.nio.file.Files.deleteIfExists(temporary)
    throw failure
  }
  java.nio.file.Files.deleteIfExists(temporary)
}

def shellQuote(value) {
  return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def nearestExistingPath(path) {
  def current = file(path).toAbsolutePath()
  return java.nio.file.Files.exists(current) ? current : nearestExistingPath(current.parent)
}

def publicationMode(workDirectory, outputDirectory) {
  try {
    return java.nio.file.Files.getFileStore(nearestExistingPath(workDirectory)) ==
      java.nio.file.Files.getFileStore(nearestExistingPath(outputDirectory)) ? 'link' : 'copy'
  } catch (Exception ignored) {
    return 'copy'
  }
}

// Resolve once before process declarations: publishDir mode must be a String,
// and strict process directive scope does not resolve script helper calls.
params.publication_mode = publicationMode(workflow.workDir, params.analysis_output_dir)

def cachedBuild(entry, products, command) {
  def productOptions = products.collect { '--product ' + shellQuote(it) }.join(' ')
  return """
  set -euo pipefail
  cat > .cache.identity.json <<'HIDEF_CACHE_IDENTITY'
  ${entry.json}
  HIDEF_CACHE_IDENTITY
  cat > .cache.build.sh <<'HIDEF_CACHE_BUILD'
  set -euo pipefail
  ${command}
  HIDEF_CACHE_BUILD
  python3 ${shellQuote(params.cache_helper)} run --root ${shellQuote(params.prepared_cache_root)} --identity .cache.identity.json ${productOptions} -- bash .cache.build.sh
  """.stripIndent()
}

def canonicalBarcodePair(barcodeA, barcodeB) {
  return barcodeA == barcodeB ? barcodeA : [barcodeA, barcodeB].toSorted().join('-')
}

def parseBarcodeIds(rawBarcodeIds, fieldName, contextLabel) {
  if (!rawBarcodeIds) {
    error "${contextLabel}: ${fieldName} must be defined."
  }

  def barcodeIds = rawBarcodeIds.split('-') as List
  if (barcodeIds.size() > 2) {
    error "${contextLabel}: ${fieldName}='${rawBarcodeIds}' must contain either one barcode_id or two barcode_ids separated by '-'."
  }

  barcodeIds.each { barcodeId ->
    if (!(barcodeId ==~ /[A-Za-z0-9_]+/)) {
      error "${contextLabel}: barcode_id '${barcodeId}' contains invalid characters. Allowed pattern: [A-Za-z0-9_]+"
    }
  }

  if (barcodeIds.size() == 2 && barcodeIds[0] == barcodeIds[1]) {
    error "${contextLabel}: ${fieldName}='${rawBarcodeIds}' must contain two different barcode_ids when two are provided."
  }

  return [
    raw: rawBarcodeIds,
    ids: barcodeIds,
    mode: barcodeIds.size() == 1 ? 'same' : 'different',
    canonical: barcodeIds.size() == 1 ? canonicalBarcodePair(barcodeIds[0], barcodeIds[0]) : canonicalBarcodePair(barcodeIds[0], barcodeIds[1])
  ]
}

/*****************************************************************
 * Main Workflow
 *****************************************************************/
workflow {

  //******************
  // General configuration
  //******************

  // Save copy of parameters file to logs directory
  logsDir = file(sharedLogsDir())
  logsDir.mkdirs()

  options = new org.yaml.snakeyaml.DumperOptions()
  options.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK)
  options.setIndent(2)

  yaml = new org.yaml.snakeyaml.Yaml(options)
  originalParametersFile = originalParamsFile(workflow.commandLine)
  originalConfiguration = yaml.load(originalParametersFile.text)
  if (!(originalConfiguration instanceof Map)) error 'The original parameters YAML must contain a mapping'
  if (params.containsKey('paramsFileName')) {
    error 'paramsFileName is reserved for the generated effective YAML; use the original input configuration'
  }
  // Apply intentional command-line/config overrides to original keys, without
  // leaking unrelated injected defaults into R's scientific configuration.
  scientificConfiguration = new LinkedHashMap()
  originalConfiguration.each { name, value ->
    scientificConfiguration[name.toString()] = plainConfigurationValue(params.containsKey(name) ? params[name] : value)
  }
  timestamp = "${new Date().format('yyyy_MMdd_HHmmss_SSS')}.${workflow.sessionId}".toString()
  file("${logsDir}/runParams.${timestamp}.yaml").text = yaml.dump(scientificConfiguration)

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

  // paramsFileName is assigned once below, after resolving prepared artifacts.

  // Define parameters file components that are checked for changes to determine if a process is rerun upon resume
  signature_params = params + [sharedFunctionsHash: fileSha256("${workflow.projectDir}/bin/sharedFunctions.R")]
  config_signatures = [
    installBSgenome: configHash(signature_params + [installBSgenomeScriptHash: fileSha256("${workflow.projectDir}/bin/installBSgenome.R")], ['cache_dir', 'genome_fasta', 'genome_organism', 'circular_chromosomes', 'installBSgenomeScriptHash']),
    processGermlineVCFs: configHash(signature_params + [processGermlineVCFsScriptHash: fileSha256("${workflow.projectDir}/bin/processGermlineVCFs.R")], ['cache_dir', 'circular_chromosomes', 'individuals', 'genome_fasta', 'genome_organism', 'bcftools_bin', 'processGermlineVCFsScriptHash', 'sharedFunctionsHash']),
    extractCallsChunk: configHash(signature_params + [extractCallsScriptHash: fileSha256("${workflow.projectDir}/bin/extractCalls.R")], ['cache_dir', 'circular_chromosomes', 'genome_fasta', 'genome_organism', 'call_types', 'chromgroups', 'runs', 'barcodes', 'min_strand_overlap', 'extractCallsScriptHash', 'sharedFunctionsHash']),
    filterCallsChunkChromgroupFiltergroup: configHash(signature_params + [filterCallsScriptHash: fileSha256("${workflow.projectDir}/bin/filterCalls.R")], ['bcftools_bin', 'cache_dir', 'call_types', 'chromgroups', 'circular_chromosomes', 'filtergroups', 'genome_fai', 'genome_fasta', 'genome_organism', 'germline_vcf_types', 'individuals', 'region_filters', 'samples', 'wigToBigWig_bin', 'wiggletools_bin', 'filterCallsScriptHash', 'sharedFunctionsHash']),
    calculateBurdensChromgroupFiltergroup: configHash(signature_params + [calculateBurdensScriptHash: fileSha256("${workflow.projectDir}/bin/calculateBurdens.R"), coverageAnnotatorSourceHash: fileSha256("${workflow.projectDir}/bin/annotateCoverage.cpp")], ['analysis_id', 'bcftools_bin', 'bedtools_bin', 'bgzip_bin', 'cache_dir', 'call_types', 'chromgroups', 'circular_chromosomes', 'genome_fai', 'genome_fasta', 'genome_organism', 'individuals', 'mitochondrial_chromosome', 'samples', 'sensitivity_parameters', 'sex_chromosomes', 'tabix_bin', 'calculateBurdensScriptHash', 'coverageAnnotatorSourceHash', 'sharedFunctionsHash']),
    outputResultsSample: configHash(signature_params + [outputResultsScriptHash: fileSha256("${workflow.projectDir}/bin/outputResults.R")], ['analysis_id', 'cache_dir', 'call_types', 'chromgroups', 'circular_chromosomes', 'filtergroups', 'genome_fasta', 'genome_organism', 'region_filters', 'samples', 'outputResultsScriptHash', 'sharedFunctionsHash'])
  ]

  // Scoped prepared artifacts: changing one individual's VCF inputs does not
  // invalidate the reference, other individuals, or independent region tracks.
  params.cache_helper = "${workflow.projectDir}/bin/artifactCache.py".toString()
  params.prepared_cache_root = "${file(params.cache_dir).toAbsolutePath()}/prepared".toString()
  inputIdentities = [:]
  identifyInput = { source ->
    def name = source.toString()
    if (!inputIdentities.containsKey(name)) {
      inputIdentities[name] = preparedInputIdentity(source)
    }
    inputIdentities[name]
  }
  containerIdentity = file(params.hidefseq_container).exists() ? identifyInput.call(params.hidefseq_container) : [kind: 'container-uri', uri: params.hidefseq_container]
  workflowSource = file("${workflow.projectDir}/main.nf").text
  cacheArtifacts = [:]
  makeCacheEntry = { namespace, processName, settings, inputs, scripts, products ->
    def start = workflowSource.indexOf("\nprocess ${processName} {")
    def end = workflowSource.indexOf('\nprocess ', start + 1)
    def processSource = workflowSource.substring(start, end < 0 ? workflowSource.length() : end)
    def scriptHashes = scripts.collectEntries { script -> [(script): fileSha256("${workflow.projectDir}/bin/${script}")] }
    scriptHashes['process'] = java.security.MessageDigest.getInstance('SHA-256').digest(processSource.getBytes('UTF-8')).encodeHex().toString()
    scriptHashes['artifactCache.py'] = fileSha256(params.cache_helper)
    def toolIdentities = [container: containerIdentity]
    settings.findAll { name, value -> name in ['seqkit', 'bgzip', 'tabix', 'bcftools', 'samtools', 'bedGraphToBigWig', 'wiggletools', 'wigToBigWig'] }.each { name, executable ->
      toolIdentities[name] = file(executable).exists() ? identifyInput.call(executable) : [kind: 'container-command', command: executable]
    }
    def identity = canonicalCacheValue([schema: 1, namespace: namespace, settings: settings,
      inputs: inputs.collectEntries { name, source -> [(name): identifyInput.call(source)] },
      scripts: scriptHashes, tools: toolIdentities])
    def serialized = groovy.json.JsonOutput.toJson(identity)
    def digest = java.security.MessageDigest.getInstance('SHA-256').digest(serialized.getBytes('UTF-8')).encodeHex().toString()
    def directory = "${params.prepared_cache_root}/v1/${namespace}/${digest}".toString()
    products.each { product ->
      def destination = "${directory}/${product}".toString()
      if (cacheArtifacts.containsKey(product) && cacheArtifacts[product] != destination) {
        error "Prepared cache filename collision for '${product}'. Use distinct input basenames."
      }
      cacheArtifacts[product] = destination
    }
    def envelope = [schema: 1, namespace: namespace, serialized_identity: serialized]
    [json: groovy.json.JsonOutput.toJson(envelope), key: digest, directory: directory]
  }
  referenceEntry = makeCacheEntry.call('reference', 'installBSgenome',
    [organism: params.genome_organism, circular: params.circular_chromosomes, fastaName: file(params.genome_fasta).name],
    [fasta: params.genome_fasta], ['installBSgenome.R', 'sharedFunctions.R'], [])
  referenceSummaryEntry = makeCacheEntry.call('reference-summary', 'prepareReferenceSummary',
    [reference: referenceEntry.key], [:],
    ['prepareReferenceSummary.R', 'referenceSummaryFunctions.R', 'sharedFunctions.R'], ['referenceSummary.qs2'])
  trinucleotideEntry = makeCacheEntry.call('trinucleotides', 'extractGenomeTrinucleotides',
    [seqkit: params.seqkit_bin, bgzip: params.bgzip_bin, tabix: params.tabix_bin],
    [fasta: params.genome_fasta], [], ["${file(params.genome_fasta).name}.bed.gz".toString(), "${file(params.genome_fasta).name}.bed.gz.tbi".toString()])
  vcfEntries = [:]
  bamEntries = [:]
  params.individuals.each { individual ->
    def bamName = file(individual.germline_bam_file).name
    def vcfInputs = [fasta: params.genome_fasta]
    (individual.germline_vcf_files ?: []).eachWithIndex { vcf, i -> vcfInputs["vcf${i}"] = vcf.germline_vcf_file }
    vcfEntries[individual.individual_id] = makeCacheEntry.call('germline-vcf', 'processGermlineVCFs',
      [individual: individual.individual_id, vcfs: individual.germline_vcf_files, reference: referenceEntry.key, bcftools: params.bcftools_bin],
      vcfInputs, ['processGermlineVCFs.R', 'sharedFunctions.R'], ["${individual.individual_id}.${bamName}.germline_vcf_variants.qs2".toString()])
    bamEntries[bamName] = makeCacheEntry.call('germline-bam', 'processGermlineBAMs',
      [type: individual.germline_bam_type, samtools: params.samtools_bin, bcftools: params.bcftools_bin, bedGraphToBigWig: params.bedGraphToBigWig_bin],
      [bam: individual.germline_bam_file, fasta: params.genome_fasta, fai: params.genome_fai], [],
      ["${bamName}.bw".toString(), "${bamName}.vcf.gz".toString(), "${bamName}.vcf.gz.tbi".toString()])
  }
  regionEntries = [:]
  (params.region_filters ?: []).each { group ->
    ((group.read_filters ?: []) + (group.genome_filters ?: [])).each { region ->
      def product = "${file(region.region_filter_file).name}.bin${region.binsize}.${region.threshold}.bw".toString()
      regionEntries[product] = makeCacheEntry.call('region-filter', 'prepareRegionFilters',
        [binsize: region.binsize, threshold: region.threshold, wiggletools: params.wiggletools_bin, wigToBigWig: params.wigToBigWig_bin],
        [region: region.region_filter_file, fai: params.genome_fai], [], [product])
    }
  }
  coverageEntries = [:]
  coverageConfigurations = []
  coverageThresholds = params.filtergroups.collect { it.min_germlineBAM_TotalReads }.unique()
  params.individuals.each { individual ->
    def bamName = file(individual.germline_bam_file).name
    coverageThresholds.each { threshold ->
      def product = "${individual.individual_id}.${bamName}.minCoverage${threshold}.qs2".toString()
      def entry = makeCacheEntry.call('germline-coverage-filter', 'prepareGermlineCoverageFilters',
        [individual: individual.individual_id, threshold: threshold, rawCoverage: bamEntries[bamName].key,
         wiggletools: params.wiggletools_bin, wigToBigWig: params.wigToBigWig_bin],
        [fai: params.genome_fai], ['prepareGermlineCoverageFilters.R'], [product])
      coverageEntries[product] = entry
      coverageConfigurations << [individual_id: individual.individual_id, threshold: threshold,
        bigwig_name: "${bamName}.bw".toString(), product: product, file: "${entry.directory}/${product}".toString()]
    }
  }
  params.germline_coverage_filters = coverageConfigurations
  params.prepared_cache = [reference: referenceEntry, reference_summary: referenceSummaryEntry, trinucleotides: trinucleotideEntry, vcfs: vcfEntries, bams: bamEntries, regions: regionEntries, coverage: coverageEntries]
  params.reference_cache_dir = "${referenceEntry.directory}/library".toString()
  params.reference_summary_file = "${referenceSummaryEntry.directory}/referenceSummary.qs2".toString()
  params.cache_artifacts = cacheArtifacts
  // Keep the user's original YAML untouched; R consumers receive the resolved
  // immutable product paths through an effective configuration for this run.
  def effectiveConfiguration = new LinkedHashMap(scientificConfiguration)
  ['germline_coverage_filters', 'reference_cache_dir', 'reference_summary_file', 'cache_artifacts'].each { name ->
    effectiveConfiguration[name] = plainConfigurationValue(params[name])
  }
  def effectiveYaml = yaml.dump(effectiveConfiguration)
  def effectiveDigest = java.security.MessageDigest.getInstance('SHA-256').digest(effectiveYaml.getBytes('UTF-8')).encodeHex().toString()
  def effectiveParametersFile = "${logsDir}/effectiveParams.${effectiveDigest}.yaml".toString()
  if (file(effectiveParametersFile).toAbsolutePath() == originalParametersFile) error 'Effective YAML must not replace its input'
  writeImmutableConfiguration(effectiveParametersFile, effectiveYaml)
  params.paramsFileName = effectiveParametersFile


  // Validate barcodes section and build barcode map
  if (!params.barcodes || params.barcodes.isEmpty()) {
    error "'barcodes' section must be defined and contain at least one barcode entry."
  }

  barcode_seq_by_id = [:]
  barcode_id_by_seq = [:]
  params.barcodes.eachWithIndex { barcodeEntry, idx ->
    if (!(barcodeEntry.barcode_id ==~ /[A-Za-z0-9_]+/)) {
      error "barcodes[${idx}]: barcode_id '${barcodeEntry.barcode_id}' contains invalid characters. Allowed pattern: [A-Za-z0-9_]+"
    }
    if (barcode_seq_by_id.containsKey(barcodeEntry.barcode_id) && barcode_seq_by_id[barcodeEntry.barcode_id] != barcodeEntry.barcode) {
      error "Duplicate barcode_id '${barcodeEntry.barcode_id}' with conflicting barcode sequences in barcodes section."
    }
    if (barcode_id_by_seq.containsKey(barcodeEntry.barcode) && barcode_id_by_seq[barcodeEntry.barcode] != barcodeEntry.barcode_id) {
      error "Barcode sequence '${barcodeEntry.barcode}' is assigned to multiple barcode_id values ('${barcode_id_by_seq[barcodeEntry.barcode]}' and '${barcodeEntry.barcode_id}')."
    }
    barcode_seq_by_id[barcodeEntry.barcode_id] = barcodeEntry.barcode
    barcode_id_by_seq[barcodeEntry.barcode] = barcodeEntry.barcode_id
  }

  sample_to_individual = params.samples.collectEntries { [ (it.sample_id): it.individual_id ] }

  provisional_run_sample_configs = params.runs.collectMany { run ->
    run.samples.collect { sample ->
      def sampleContext = "run '${run.run_id}', sample '${sample.sample_id}'"
      def round1 = parseBarcodeIds(sample.barcode_ids, 'barcode_ids', sampleContext)
      round1.ids.each { barcodeId ->
        if (!barcode_seq_by_id.containsKey(barcodeId)) {
          error "${sampleContext}: barcode_id '${barcodeId}' from barcode_ids is not defined in top-level barcodes section."
        }
      }

      def round2 = null
      if (sample.barcode_ids_round2) {
        round2 = parseBarcodeIds(sample.barcode_ids_round2, 'barcode_ids_round2', sampleContext)
        round2.ids.each { barcodeId ->
          if (!barcode_seq_by_id.containsKey(barcodeId)) {
            error "${sampleContext}: barcode_id '${barcodeId}' from barcode_ids_round2 is not defined in top-level barcodes section."
          }
        }
      }

      [
        run_id: run.run_id,
        sample_id: sample.sample_id,
        individual_id: sample_to_individual[sample.sample_id],
        barcode_ids: round1.raw,
        barcode_ids_parsed: round1,
        barcode_ids_round2: round2?.raw,
        barcode_ids_round2_parsed: round2
      ]
    }
  }

  provisional_run_sample_configs
    .groupBy { sampleConfig -> sampleConfig.sample_id }
    .each { sampleId, sampleConfigs ->
      def round2Usage = sampleConfigs.collect { sampleConfig -> sampleConfig.barcode_ids_round2_parsed != null }.unique()
      if (round2Usage.size() > 1) {
        def runIds = sampleConfigs.collect { sampleConfig -> sampleConfig.run_id }.unique().sort().join(', ')
        error "sample '${sampleId}' has inconsistent barcode_ids_round2 usage across runs (${runIds}). Configure this sample with either barcode_ids only in all runs or with both barcode_ids and barcode_ids_round2 in all runs."
      }
    }

  run_sample_configs = provisional_run_sample_configs
    .groupBy { sampleConfig -> sampleConfig.run_id }
    .collectMany { runId, runConfigs ->
      def round1OnlyConfigs = runConfigs.findAll { cfg -> !cfg.barcode_ids_round2_parsed }
      def round2Configs = runConfigs.findAll { cfg -> cfg.barcode_ids_round2_parsed }

      round1OnlyConfigs.groupBy { cfg ->
        cfg.barcode_ids_parsed.canonical
      }.each { key, cfgs ->
        if (cfgs.size() > 1) {
          error "run '${runId}' has duplicate barcode_ids configuration '${key}' across round1-only samples: ${cfgs.collect { cfg -> cfg.sample_id }.join(', ')}"
        }
      }

      if (round2Configs) {
        round2Configs.groupBy { cfg ->
          "${cfg.barcode_ids_parsed.canonical}|${cfg.barcode_ids_round2_parsed.canonical}"
        }.each { key, cfgs ->
          if (cfgs.size() > 1) {
            error "run '${runId}' has duplicate barcode_ids + barcode_ids_round2 configuration '${key}' across samples: ${cfgs.collect { cfg -> cfg.sample_id }.join(', ')}. For samples using two demultiplexing rounds, each sample in a run must have a unique combination across both rounds."
          }
        }

        def overlappingRound1Keys = round1OnlyConfigs.collect { cfg -> cfg.barcode_ids_parsed.canonical }.intersect(
          round2Configs.collect { cfg -> cfg.barcode_ids_parsed.canonical }
        )
        if (overlappingRound1Keys) {
          error "run '${runId}' reuses round1 barcode_ids across round1-only and round2 samples (${overlappingRound1Keys.join(', ')}). Round1 barcode_ids must be unique for round1-only samples relative to all other samples in the same run."
        }
      }

      def round1Modes = runConfigs.collect { cfg -> cfg.barcode_ids_parsed.mode }.unique()
      if (round1Modes.size() > 1) {
        error "run '${runId}' mixes barcode_ids demultiplexing modes across samples (${round1Modes.join(', ')}). Use a single mode ('same' or 'different') per run for round1 demultiplexing."
      }

      if (round2Configs) {
        def round2Modes = round2Configs.collect { cfg -> cfg.barcode_ids_round2_parsed.mode }.unique()
        if (round2Modes.size() > 1) {
          error "run '${runId}' mixes barcode_ids_round2 demultiplexing modes across samples (${round2Modes.join(', ')}). Use a single mode ('same' or 'different') per run for round2 demultiplexing."
        }
      }

      runConfigs.collect { cfg -> cfg + [round2_enabled: cfg.barcode_ids_round2_parsed != null] }
    }

  // Create a channel of runs
  runs_ch = channel.fromList(params.runs)

  // Create channel for the input reads file.
  reads_ch = runs_ch.map { run ->
      tuple(
        run.run_id,
        file(run.reads_file),
        file("${run.reads_file}.pbi")
      )
  }

  //******************
  // makeBarcodesFasta
  //******************

  round2_run_ids = run_sample_configs
    .findAll { cfg -> cfg.round2_enabled }
    .collect { cfg -> cfg.run_id }
    .unique()

  round2_sample_keys = run_sample_configs
    .findAll { cfg -> cfg.round2_enabled }
    .collect { cfg -> "${cfg.run_id}|${cfg.individual_id}|${cfg.sample_id}|${cfg.barcode_ids}" }
    .toSet()

  // Create input channel
  makeBarcodesFasta_input_ch = runs_ch.map { run ->
    def barcodeIdsInRun = run_sample_configs
      .findAll { cfg -> cfg.run_id == run.run_id }
      .collectMany { cfg -> cfg.barcode_ids_parsed.ids }
      .unique()

    def runBarcodeFasta = barcodeIdsInRun.collect { barcodeId ->
      ">${barcodeId}\n${barcode_seq_by_id[barcodeId]}"
    }.join("\n")

    tuple(run.run_id, 'round1', runBarcodeFasta)
  }.mix(
    channel.fromList(round2_run_ids).map { run_id ->
      def barcodeIdsInRun = run_sample_configs
        .findAll { cfg -> cfg.run_id == run_id && cfg.barcode_ids_round2_parsed != null }
        .collectMany { cfg -> cfg.barcode_ids_round2_parsed.ids }
        .unique()

      def runBarcodeFasta = barcodeIdsInRun.collect { barcodeId ->
        ">${barcodeId}\n${barcode_seq_by_id[barcodeId]}"
      }.join("\n")

      tuple(run_id, 'round2', runBarcodeFasta)
    }
  )

  // Run process
  makeBarcodesFasta(makeBarcodesFasta_input_ch)


  //******************
  // ccs
  //******************

  // Run ccs if read_type is subreads, and create channels for subsequent processes

  filterAdapter_input_ch = null
  countZMWs_initial_ch = null

  if( params.reads_type == 'subreads' ) {
      // Run CCS in chunks.
      ccsChunk(
        reads_ch
          .combine(channel.from(1..params.ccs_chunks))
          .map { readTuple -> tuple(readTuple[0], readTuple[1], readTuple[2], readTuple[3]) }
      )

      mergeCCSchunks_input_ch = ccsChunk.out.bampbi_tuple
          .groupTuple(by: 0) // Group by run_id
          .map { run_id, bamFiles, pbiFiles, chunkIDs ->
            // Sort by chunkID
            def sortedIndices = (0..<chunkIDs.size()).toList().sort { i -> chunkIDs[i] as int }
            def sortedBamFiles = sortedIndices.collect { i -> bamFiles[i] }
            def sortedPbiFiles = sortedIndices.collect { i -> pbiFiles[i] }
            tuple(run_id, sortedBamFiles, sortedPbiFiles)
          }

      mergeCCSchunks(mergeCCSchunks_input_ch)

      filterAdapter_input_ch = mergeCCSchunks.out

      countZMWs_initial_ch = reads_ch.map { f -> tuple(f[1], f[2], "zmwcount.txt") }
        .mix(mergeCCSchunks.out.map { f -> tuple(f[1], f[2], "zmwcount.txt") })
  }
  else if( params.reads_type == 'ccs' ) {
      filterAdapter_input_ch = reads_ch

      countZMWs_initial_ch = reads_ch.map { f -> tuple(f[1], f[2], "zmwcount.txt") }
  }
  else {
      error "Unsupported reads_type '${params.reads_type}'."
  }

  //******************
  // filterAdapter
  //******************

  // Run process
  filterAdapter(filterAdapter_input_ch)

  //******************
  // limaDemux
  //******************

  limaDemux_round1_mode_ch = channel.fromList(run_sample_configs)
    .map { cfg -> tuple(cfg.run_id, cfg.barcode_ids_parsed.mode) }
    .unique()

  round1_barcodes_fasta_ch = makeBarcodesFasta.out
    .filter { run_id, barcodesFasta, demux_round -> demux_round == 'round1' }
    .map { run_id, barcodesFasta, demux_round -> tuple(run_id, barcodesFasta) }

  round2_barcodes_fasta_ch = makeBarcodesFasta.out
    .filter { run_id, barcodesFasta, demux_round -> demux_round == 'round2' }
    .map { run_id, barcodesFasta, demux_round -> tuple(run_id, barcodesFasta) }

  // Create input channel
  limaDemux_round1_input_ch = filterAdapter.out
      .join(round1_barcodes_fasta_ch, by: 0)
      .map { run_id, bamFile, pbiFile, barcodesFasta ->
        tuple(run_id, bamFile, pbiFile, barcodesFasta)
      }
      .join(limaDemux_round1_mode_ch, by: 0)
      .map { run_id, bamFile, pbiFile, barcodesFasta, mode ->
        tuple(run_id, null, null, null, null, null, bamFile, pbiFile, barcodesFasta, mode, params.lima_supplemental_settings)
      }

  limaDemux_round1 = limaDemuxRound1(limaDemux_round1_input_ch)

  limaDemux_round1_map_ch = limaDemux_round1.bam
    .transpose()
    .map { run_id, bamFile ->
      def m = bamFile.name =~ /.*\.demux\.([A-Za-z0-9_]+)--([A-Za-z0-9_]+)\.bam$/
      if (!m) {
          error "Can't match BAM file name to run_id and barcode_id: ${bamFile.name}"
      }
      def barcode_pair_key = canonicalBarcodePair(m[0][1], m[0][2])
      tuple(run_id, barcode_pair_key, bamFile)
    }

  mergeDemuxBams_round1_input_ch = channel.fromList(run_sample_configs)
    .map { cfg ->
      tuple(cfg.run_id, cfg.barcode_ids_parsed.canonical, cfg.individual_id, cfg.sample_id, cfg.barcode_ids)
    }
    .join(limaDemux_round1_map_ch, by: [0,1])
    .map { run_id, barcode_pair_key, individual_id, sample_id, barcode_ids, bamFile ->
      tuple(run_id, individual_id, sample_id, barcode_ids, bamFile)
    }
    .groupTuple(by: [0,1,2,3])

  mergeDemuxBams_round1 = mergeDemuxBamsRound1(mergeDemuxBams_round1_input_ch).out

  limaDemux_round2_samples_ch = channel.fromList(run_sample_configs)
    .filter { cfg -> cfg.round2_enabled }
    .map { cfg -> tuple(cfg.run_id, cfg.individual_id, cfg.sample_id, cfg.barcode_ids, cfg.barcode_ids_round2, cfg.barcode_ids_round2_parsed.mode, cfg.barcode_ids_round2_parsed.canonical) }

  limaDemux_round2_input_ch = mergeDemuxBams_round1
    .join(limaDemux_round2_samples_ch, by: [0,1,2,3])
    .combine(round2_barcodes_fasta_ch, by: 0)
    .map { run_id, individual_id, sample_id, barcode_ids, mergedBam, mergedPbi, barcode_ids_round2, mode2, barcode_pair_key_round2, barcodesFasta ->
      tuple(run_id, individual_id, sample_id, barcode_ids, barcode_ids_round2, barcode_pair_key_round2, mergedBam, mergedPbi, barcodesFasta, mode2, params.lima_round2_supplemental_settings)
    }

  limaDemux_round2 = limaDemuxRound2(limaDemux_round2_input_ch)

  limaDemux_round2_map_ch = limaDemux_round2.bam
    .transpose()
    .map { run_id, individual_id, sample_id, barcode_ids, barcode_ids_round2, barcode_pair_key_round2, bamFile ->
      def m = bamFile.name =~ /.*\.demux\.([A-Za-z0-9_]+)--([A-Za-z0-9_]+)\.bam$/
      if (!m) {
          error "Can't match BAM file name to run_id and barcode_id: ${bamFile.name}"
      }
      def barcode_pair_key = canonicalBarcodePair(m[0][1], m[0][2])
      tuple(run_id, individual_id, sample_id, barcode_ids, barcode_ids_round2, barcode_pair_key_round2, barcode_pair_key, bamFile)
    }
    .filter { run_id, individual_id, sample_id, barcode_ids, barcode_ids_round2, barcode_pair_key_round2, barcode_pair_key, bamFile ->
      barcode_pair_key == barcode_pair_key_round2
    }
    .map { run_id, individual_id, sample_id, barcode_ids, barcode_ids_round2, barcode_pair_key_round2, barcode_pair_key, bamFile ->
      tuple(run_id, individual_id, sample_id, "${barcode_ids}.${barcode_ids_round2}", bamFile)
    }
  
  mergeDemuxBams_round2_input_ch = limaDemux_round2_map_ch
    .groupTuple(by: [0,1,2,3])

  mergeDemuxBams_round2 = mergeDemuxBamsRound2(mergeDemuxBams_round2_input_ch).out

  //******************
  // pbmm2Align
  //******************

  // Create input channel
  pbmm2_round1_final_ch = mergeDemuxBams_round1
    .filter { run_id, individual_id, sample_id, barcode_ids, bamFile, pbiFile ->
      !round2_sample_keys.contains("${run_id}|${individual_id}|${sample_id}|${barcode_ids}")
    }
    .map { run_id, individual_id, sample_id, barcode_ids, bamFile, pbiFile -> tuple(run_id, individual_id, sample_id, barcode_ids, bamFile) }

  pbmm2_round2_final_ch = mergeDemuxBams_round2
    .map { run_id, individual_id, sample_id, barcode_ids, bamFile, pbiFile -> tuple(run_id, individual_id, sample_id, barcode_ids, bamFile) }

  pbmm2Align_input_ch = pbmm2_round1_final_ch
    .mix(pbmm2_round2_final_ch)

  // Run process
  pbmm2Align(pbmm2Align_input_ch)

  //******************
  // verifyBAMID
  //******************

  if (params.verifybamid_bin && params.verifybamid_resource_UD && params.verifybamid_resource_Bed && params.verifybamid_resource_Mean) {
    // Create input channel
    verifyBAMID_input_ch = pbmm2Align.out.map { run_id, individual_id, sample_id, barcode_id, bamFile, pbiFile ->
      tuple(run_id, individual_id, sample_id, barcode_id, bamFile)
    }

    // Run process
    verifyBAMID(verifyBAMID_input_ch)
  }

  //******************
  // countZMWs
  //******************

  // Create input channel
  countZMWs_input_ch = countZMWs_initial_ch.mix(
      filterAdapter.out.map { f -> tuple(f[1], f[2], "zmwcount.txt") },
      mergeDemuxBams_round1.map { run_id, individual_id, sample_id, barcode_ids, bamFile, pbiFile -> tuple(bamFile, pbiFile, "zmwcount.txt") },
      mergeDemuxBams_round2.map { run_id, individual_id, sample_id, barcode_ids, bamFile, pbiFile -> tuple(bamFile, pbiFile, "zmwcount.txt") },
      pbmm2Align.out.map { run_id, individual_id, sample_id, barcode_id, bamFile, pbiFile -> tuple(bamFile, pbiFile, "zmwcount.txt") }
    )

  // Run process
  countZMWs(countZMWs_input_ch)

  //******************
  // mergeAlignedSampleBAMs
  //******************

  // Preserve run ordering from params.runs so merge order matches the YAML configuration
  run_id_order = params.runs.withIndex().collectEntries { run, idx -> [(run.run_id): idx] }

  // Create input channel
  mergeAlignedSampleBAMs_input_ch = pbmm2Align.out
    .map { run_id, individual_id, sample_id, barcode_id, bamFile, pbiFile ->
        tuple(individual_id, sample_id, run_id_order[run_id], bamFile, pbiFile)
    }
    .groupTuple(by: [0, 1]) // Group by individual_id, sample_id
    .map { individual_id, sample_id, run_order, bamFiles, pbiFiles ->
      def ordered = [run_order, bamFiles, pbiFiles].transpose().sort { row -> row[0] }
      tuple(
        individual_id,
        sample_id,
        ordered.collect { row -> row[1] },
        ordered.collect { row -> row[2] }
      )
    }

  // Run process
  mergeAlignedSampleBAMs(mergeAlignedSampleBAMs_input_ch)

  //******************
  // splitBAM
  //******************

  countAnalysisZMWs(mergeAlignedSampleBAMs.out)
  compileBamDispatcher(channel.value(file("${projectDir}/bin/splitBamByZmw.cpp", checkIfExists: true)))

  // Keep the exact legacy enumeration once, and dispatch all chunks per sample.
  splitBAM_input_ch = mergeAlignedSampleBAMs.out
      .join(countAnalysisZMWs.out, by: [0, 1])
      .flatMap { individual_id, sample_id, bamFile, pbiFile, baiFile, zmwCountFile, zmwIdsFile ->
        def total_zmws = zmwCountFile.text.trim() as int
        if (total_zmws == 0) {
          log.warn "Skipping sample '${sample_id}' because its merged analysis BAM contains zero ZMWs."
          return []
        }
        def effective_chunks = Math.min(params.analysis_chunks as int, total_zmws)
        [tuple(individual_id, sample_id, bamFile, pbiFile, baiFile, zmwIdsFile, effective_chunks)]
      }

  splitBAM(splitBAM_input_ch, compileBamDispatcher.out)
  splitBAM_chunks_ch = splitBAM.out.flatMap { individual_id, sample_id, bamFiles, pbiFiles, baiFiles, effectiveChunks ->
    expandAnalysisChunks(individual_id, sample_id, bamFiles, pbiFiles, baiFiles, effectiveChunks)
  }

  //******************
  // installBSgenome
  //******************
  installBSgenome(channel.value(config_signatures.installBSgenome))
  BSgenome_name_ch = installBSgenome.out.map { bsgenomeNameFile -> bsgenomeNameFile.text.trim() }
  prepareReferenceSummary(BSgenome_name_ch)

  //******************
  // extractGenomeTrinucleotides
  //******************
  extractGenomeTrinucleotides()

  //******************
  // processGermlineVCFs
  //******************

  // Create input channel
  processGermlineVCFs_input_ch = channel
    .from(params.individuals)
    .map { individual -> tuple(individual.individual_id, file(individual.germline_bam_file)) }
    .combine(BSgenome_name_ch)
    .map { individual_id, germline_bam_file, BSgenome_name ->
      tuple(individual_id, germline_bam_file, config_signatures.processGermlineVCFs)
    }

  // Run process
  processGermlineVCFs(processGermlineVCFs_input_ch)

  //******************
  // processGermlineBAMs
  //******************

  // Create input channel
  processGermlineBAMs_input_ch = channel.fromList(params.individuals)
    .map { run ->
      tuple( file(run.germline_bam_file), run.germline_bam_type )
    }
    .unique() // Shared germline BAM/type inputs need only one preparation task.

  // Run process
  processGermlineBAMs(processGermlineBAMs_input_ch)

  // Reuse identical thresholds across filtergroups, samples, and chunks.
  prepareGermlineCoverageFilters_input_ch = processGermlineBAMs.out.coverage
    .flatMap { coverageFile ->
      params.germline_coverage_filters.findAll { entry -> entry.bigwig_name == coverageFile.name }
        .collect { entry -> tuple(entry.individual_id, entry.threshold, coverageFile, entry.product) }
    }
  prepareGermlineCoverageFilters(prepareGermlineCoverageFilters_input_ch)

  //******************
  // processGermlineBAMs
  //******************

  // Create input channel
  prepareRegionFilters_input_ch = channel.fromList(params.region_filters)
      .flatMap { region_filter ->
        def filters = []

        region_filter.read_filters?.each { filter ->
          filters << tuple(filter.region_filter_file, filter.binsize, filter.threshold)
        }

        region_filter.genome_filters?.each { filter ->
          filters << tuple(filter.region_filter_file, filter.binsize, filter.threshold)
        }

        filters
      }
      .unique() // Avoid preparing the same region filter/binsize/threshold configuration twice

  // Run process
  prepareRegionFilters(prepareRegionFilters_input_ch)

  //******************
  // extractCallsChunk
  //******************

  // Create a completion signal for all filter-related processes by collecting all outputs
  prepareFilters_done = BSgenome_name_ch
    .mix(
      prepareReferenceSummary.out,
      extractGenomeTrinucleotides.out,
      processGermlineVCFs.out,
      processGermlineBAMs.out,
      prepareGermlineCoverageFilters.out,
      prepareRegionFilters.out
    )
    .collect()
    .map { true }

  // Create input channel
  extractCalls_input_ch = splitBAM_chunks_ch
      .combine(prepareFilters_done)
      .map { individual_id, sample_id, bamFile, pbiFile, baiFile, chunkID, effectiveChunks, prepareFiltersReady ->
        tuple(individual_id, sample_id, bamFile, pbiFile, baiFile, chunkID, effectiveChunks, config_signatures.extractCallsChunk)
      }

  // Run process
  extractCallsChunk(extractCalls_input_ch)

  //******************
  // filterCallsChunkChromgroupFiltergroup
  //******************

  //Prepare a list with all call_types.analyzein_chromgroups and call_types.SBSindel_call_types.filtergroup configured pairs.
  //Created as a list first instead of a channel to allow static calculation of its size for later downstream use
  chromgroups_filtergroups_list = params.call_types
    .collectMany { call_type ->
      def chromgroup_names
      if (call_type.analyzein_chromgroups == 'all') {
        chromgroup_names = params.chromgroups.collect { chromgroup -> chromgroup.chromgroup }
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

  // Create input channel
  filterCallsChunkChromgroupFiltergroup_input_ch = extractCallsChunk.out
      .combine(channel.fromList(chromgroups_filtergroups_list))
      .map { individual_id, sample_id, extractCallsFile, chunkID, effectiveChunks, chromgroup, filtergroup ->
        tuple(individual_id, sample_id, extractCallsFile, chunkID, effectiveChunks, chromgroup, filtergroup, config_signatures.filterCallsChunkChromgroupFiltergroup)
      }

  // Run process
  filterCallsChunkChromgroupFiltergroup(filterCallsChunkChromgroupFiltergroup_input_ch)

  //******************
  // calculateBurdensChromgroupFiltergroup
  //******************
  // Create input channel
  calculateBurdensChromgroupFiltergroup_input_ch = filterCallsChunkChromgroupFiltergroup.out
      .map { individual_id, sample_id, chromgroup, filtergroup, chunkID, filterCallsFile, effectiveChunks ->
          tuple(groupKey([individual_id, sample_id, chromgroup, filtergroup], effectiveChunks), chunkID, filterCallsFile)
      }
      .groupTuple()
      .map { group, chunkIDs, filterCallsFiles ->
          def (individual_id, sample_id, chromgroup, filtergroup) = group.getGroupTarget()
          // Sort by chunkID
          def sortedIndices = (0..<chunkIDs.size()).toList().sort { i -> chunkIDs[i] as int }
          def sortedfilterCallsFiles = sortedIndices.collect { i -> filterCallsFiles[i] }
          return tuple(individual_id, sample_id, chromgroup, filtergroup, sortedfilterCallsFiles, config_signatures.calculateBurdensChromgroupFiltergroup)
      }

  // Run process
  compileCoverageAnnotator(channel.value(file("${projectDir}/bin/annotateCoverage.cpp", checkIfExists: true)))
  calculateBurdensChromgroupFiltergroup(calculateBurdensChromgroupFiltergroup_input_ch, compileCoverageAnnotator.out)

  //******************
  // outputResultsSample
  //******************

  // Create input channel
  outputResultsSample_input_ch = calculateBurdensChromgroupFiltergroup.out.tuple_qs2
      .map { individual_id, sample_id, chromgroup, filtergroup, calculateBurdensFile ->
          tuple(individual_id, sample_id, calculateBurdensFile)
      }
      .groupTuple(by: [0, 1], size: chromgroups_filtergroups_list.size()) // Group by individual_id, sample_id. Emit as soon as each sample's chromgroup/filtergroup analyses finish.
      .map { individual_id, sample_id, calculateBurdensFiles ->
          tuple(individual_id, sample_id, calculateBurdensFiles, config_signatures.outputResultsSample)
      }

  // Run process
  outputResultsSample(outputResultsSample_input_ch)

}


/*****************************************************************
 * Process Definitions
 *****************************************************************/

process makeBarcodesFasta {
    cpus 1
    memory '2 GB'
    time '10m'
    tag { "makeBarcodesFasta: ${run_id} ${demux_round}" }
    container "${params.hidefseq_container}"

    publishDir path: { sharedLogsDir() }, mode: "copy", pattern: "*.barcodes.fasta"

    afterScript {
      generateAfterScript(
        sharedLogsDir(),
        "${task.process}.${params.analysis_id}.${run_id}.${demux_round}.command.log"
      )
    }

    input:
      tuple val(run_id), val(demux_round), val(content)

    output:
      tuple val(run_id), path("${run_id}.${demux_round}.barcodes.fasta"), val(demux_round)

    script:
    """
    echo -e "${content}" > ${run_id}.${demux_round}.barcodes.fasta
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

    publishDir path: { sharedLogsDir() }, mode: "copy", pattern: "statistics/*.ccs_report.*", saveAs: { filename -> new File(filename).getName() }
    publishDir path: { sharedLogsDir() }, mode: "copy", pattern: "statistics/*.summary.json", saveAs: { filename -> new File(filename).getName() }

    afterScript {
      generateAfterScript(
        sharedLogsDir(),
        "${task.process}.${params.analysis_id}.${run_id}.chunk${chunkID}.command.log"
      )
    }

    input:
      tuple val(run_id), path(bamFile), path(pbiFile), val(chunkID)

    output:
      tuple val(run_id), path("hifi_reads/${run_id}.chunk${chunkID}.hifi_reads.ccs.bam"), path("hifi_reads/${run_id}.chunk${chunkID}.hifi_reads.ccs.bam.pbi"), val(chunkID), emit: bampbi_tuple
      path "statistics/*.ccs_report.*", emit: report
      path "statistics/*.summary.json", emit: summary

    script:
    // Build the LD_PRELOAD command if the parameter is set.
    def ld_preload_cmd = (params.ccs_ld_preload && params.ccs_ld_preload.trim()) ? "export LD_PRELOAD=${params.ccs_ld_preload}" : ""

    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    ${ld_preload_cmd}
    ccs -j ${task.cpus} --log-level INFO --by-strand --hifi-kinetics --instrument-files-layout --min-rq -1 --top-passes 255 \\
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

    afterScript {
      generateAfterScript(
        sharedLogsDir(),
        "${task.process}.${params.analysis_id}.${run_id}.command.log"
      )
    }

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
    memory '4 GB'
    time '4h'
    tag "filterAdapter"
    container "${params.hidefseq_container}"

    afterScript {
      generateAfterScript(
        sharedLogsDir(),
        "${task.process}.${params.analysis_id}.${run_id}.command.log"
      )
    }

    input:
      tuple val(run_id), path(bamFile), path(pbiFile)

    output:
      tuple val(run_id), path("${run_id}.ccs.filtered.bam"), path("${run_id}.ccs.filtered.bam.pbi")

    script:
    """
    ${params.samtools_bin} view -b -@ ${task.cpus} -e "[ma]==0" ${bamFile} > ${run_id}.ccs.filtered.bam
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
    memory '8 GB'
    time '4h'
    tag { "limaDemux: ${bamFile.baseName}" }
    container "${params.hidefseq_container}"

    publishDir path: { sharedLogsDir() }, mode: "copy", pattern: "*.lima.summary"
    publishDir path: { sharedLogsDir() }, mode: "copy", pattern: "*.lima.counts"

    afterScript {
      generateAfterScript(
        sharedLogsDir(),
        "${task.process}.${params.analysis_id}.${bamFile.baseName}.command.log"
      )
    }

    input:
      tuple val(run_id), val(individual_id), val(sample_id), val(barcode_ids), val(barcode_ids_round2), val(barcode_pair_key_round2), path(bamFile), path(pbiFile), path(barcodesFasta), val(mode), val(supplemental_settings)

    output:
      tuple val(run_id), val(individual_id), val(sample_id), val(barcode_ids), val(barcode_ids_round2), val(barcode_pair_key_round2), path("*.demux.*.bam"), emit: bam, optional: true
      path "*.lima.summary", emit: lima_summary
      path "*.lima.counts", emit: lima_counts

    script:
    def modeFlags = mode == 'same' ? '--same' : '--different --keep-tag-idx-order'
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    lima --ccs --split-named ${modeFlags} ${supplemental_settings} \
         ${bamFile} \
         ${barcodesFasta} \
         ${bamFile.baseName}.demux.bam
    """
}

workflow limaDemuxRound1 {
    take:
      limaDemux_input_ch

    main:
      demux = limaDemux(limaDemux_input_ch)

    emit:
      bam = demux.bam.map { run_id, individual_id, sample_id, barcode_ids, barcode_ids_round2, barcode_pair_key_round2, bamFiles ->
        tuple(run_id, bamFiles)
      }
      lima_summary = demux.lima_summary
      lima_counts = demux.lima_counts
}

workflow limaDemuxRound2 {
    take:
      limaDemux_input_ch

    main:
      demux = limaDemux(limaDemux_input_ch)

    emit:
      bam = demux.bam
      lima_summary = demux.lima_summary
      lima_counts = demux.lima_counts
}


/*
  mergeDemuxBams: Merges one or more demultiplexed BAMs for a sample into a single BAM.
*/
process mergeDemuxBams {
    cpus 2
    memory '8 GB'
    time '2h'
    tag { "mergeDemuxBams: ${run_id} ${sample_id} ${barcode_ids}" }
    container "${params.hidefseq_container}"

    afterScript {
      generateAfterScript(
        "${params.analysis_output_dir}/${dirSampleLogs(individual_id, sample_id)}",
        "${task.process}.${params.analysis_id}.${run_id}.${individual_id}.${sample_id}.${barcode_ids}.command.log"
      )
    }

    input:
      tuple val(run_id), val(individual_id), val(sample_id), val(barcode_ids), path(demuxBams)

    output:
      tuple val(run_id), val(individual_id), val(sample_id), val(barcode_ids), path("${run_id}.${individual_id}.${sample_id}.${barcode_ids}.ccs.filtered.bam"), path("${run_id}.${individual_id}.${sample_id}.${barcode_ids}.ccs.filtered.bam.pbi")

    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbmerge -o ${run_id}.${individual_id}.${sample_id}.${barcode_ids}.ccs.filtered.bam ${demuxBams.join(' ')}
    pbindex ${run_id}.${individual_id}.${sample_id}.${barcode_ids}.ccs.filtered.bam
    """
}

workflow mergeDemuxBamsRound1 {
    take:
      mergeDemuxBams_input_ch

    main:
      out = mergeDemuxBams(mergeDemuxBams_input_ch)

    emit:
      out
}

workflow mergeDemuxBamsRound2 {
    take:
      mergeDemuxBams_input_ch

    main:
      out = mergeDemuxBams(mergeDemuxBams_input_ch)

    emit:
      out
}


/*
  pbmm2Align: Aligns a demultiplexed BAM file using pbmm2.
  and renames the output to include the sample ID.
*/
process pbmm2Align {
    cpus 8
    memory '32 GB'
    time '6h'
    tag { "pbmm2Align: ${run_id} ${sample_id} ${barcode_id}" }
    container "${params.hidefseq_container}"

    afterScript {
      generateAfterScript(
        "${params.analysis_output_dir}/${dirSampleLogs(individual_id, sample_id)}",
        "${task.process}.${params.analysis_id}.${run_id}.${individual_id}.${sample_id}.${barcode_id}.command.log"
      )
    }

    input:
      tuple val(run_id), val(individual_id), val(sample_id), val(barcode_id), path(bamFile)

    output:
      tuple val(run_id), val(individual_id), val(sample_id), val(barcode_id), path("${run_id}.${individual_id}.${sample_id}.${barcode_id}.ccs.filtered.aligned.bam"), path("${run_id}.${individual_id}.${sample_id}.${barcode_id}.ccs.filtered.aligned.bam.pbi")

    script:
    """
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbmm2 align -j ${task.cpus} --preset CCS ${params.pbmm2_supplemental_settings} ${params.genome_mmi} ${bamFile} ${run_id}.${individual_id}.${sample_id}.${barcode_id}.ccs.filtered.aligned.bam
    pbindex ${run_id}.${individual_id}.${sample_id}.${barcode_id}.ccs.filtered.aligned.bam
    """
}

/*
  verifyBAMID: Runs VerifyBamID2 on aligned BAMs output by pbmm2Align.
*/
process verifyBAMID {
    cpus 2
    memory '32 GB'
    time '8h'
    tag { "verifyBAMID: ${run_id} ${sample_id} ${barcode_id}" }
    container "${params.hidefseq_container}"

    publishDir path: "${params.analysis_output_dir}",
      mode: 'copy',
      saveAs: { filename -> "${dirVerifyBAMID(individual_id, sample_id)}/${filename}" }

    afterScript {
      generateAfterScript(
        "${params.analysis_output_dir}/${dirSampleLogs(individual_id, sample_id)}",
        "${task.process}.${params.analysis_id}.${run_id}.${individual_id}.${sample_id}.${barcode_id}.command.log"
      )
    }

    input:
      tuple val(run_id), val(individual_id), val(sample_id), val(barcode_id), path(bamFile)

    output:
      tuple val(run_id), val(individual_id), val(sample_id), val(barcode_id), path("${bamFile}.verifyBAMID.selfSM"), path("${bamFile}.verifyBAMID.Ancestry")

    script:
    """
    ${params.samtools_bin} sort -@ ${task.cpus} --write-index -o ${bamFile}.sorted.bam ${bamFile}
    ${params.verifybamid_bin} --DisableSanityCheck --UDPath ${params.verifybamid_resource_UD} --BedPath ${params.verifybamid_resource_Bed} --MeanPath ${params.verifybamid_resource_Mean} --Reference ${params.genome_fasta} --BamFile ${bamFile}.sorted.bam --Output ${bamFile}.verifyBAMID
    """
}

/*
  mergeAlignedSampleBAMs: Merges aligned BAM files from the same sample across different runs.
*/
process mergeAlignedSampleBAMs {
    cpus 4
    memory '32 GB'
    time '4h'
    tag { "mergeAlignedSampleBAMs: ${sample_id}" }
    container "${params.hidefseq_container}"

    publishDir path: "${params.analysis_output_dir}",
      mode: params.publication_mode,
      saveAs: { filename -> "${dirProcessReads(individual_id, sample_id)}/${filename}" }

    afterScript {
      generateAfterScript(
        "${params.analysis_output_dir}/${dirSampleLogs(individual_id, sample_id)}",
        "${task.process}.${params.analysis_id}.${individual_id}.${sample_id}.command.log"
      )
    }

    input:
      tuple val(individual_id), val(sample_id), path(bamFiles), path(pbiFiles)

    output:
      tuple val(individual_id), val(sample_id),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam"),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.pbi"),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.bam.bai")

    script:
    """
    sample_basename=${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned

    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    pbmerge -o \${sample_basename}.unsorted.bam ${bamFiles.join(' ')}
    conda deactivate
    
    ${params.samtools_bin} sort -@ ${task.cpus} -m 4G \${sample_basename}.unsorted.bam > \${sample_basename}.sorted.bam
    ${params.samtools_bin} index -@ ${task.cpus} \${sample_basename}.sorted.bam

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
    time '10m'
    tag { "countZMWs: ${bamFile}" }
    container "${params.hidefseq_container}"

    publishDir path: { sharedLogsDir() }, mode: "copy"

    input:
      tuple path(bamFile), path(pbiFile), val(outFileSuffix)

    output:
      path "*.${outFileSuffix}"

    script:
    """
    set -euo pipefail

    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    zmwfilter --show-all ${bamFile} | wc -l > \$(basename ${bamFile} .bam).${outFileSuffix}
    """
}

/*
  countAnalysisZMWs: Counts ZMWs in merged per-sample BAM files for analysis chunk planning.
*/
process countAnalysisZMWs {
    cpus 1
    memory '8 GB'
    time '10m'
    tag { "countAnalysisZMWs: ${sample_id}" }
    container "${params.hidefseq_container}"

    input:
      tuple val(individual_id), val(sample_id), path(bamFile), path(pbiFile), path(baiFile)

    output:
      tuple val(individual_id), val(sample_id),
      path("${params.analysis_id}.${individual_id}.${sample_id}.analysis_zmwcount.txt"),
      path("${params.analysis_id}.${individual_id}.${sample_id}.analysis_zmwIDs.txt")

    script:
    """
    set -euo pipefail
    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    zmwfilter --show-all ${bamFile} > ${params.analysis_id}.${individual_id}.${sample_id}.analysis_zmwIDs.txt
    wc -l < ${params.analysis_id}.${individual_id}.${sample_id}.analysis_zmwIDs.txt > ${params.analysis_id}.${individual_id}.${sample_id}.analysis_zmwcount.txt
    """
}

/* Compile once in the pinned container; the small source is a content-hashed input. */
process compileBamDispatcher {
    cpus 1
    memory '1 GB'
    time '10m'
    container "${params.hidefseq_container}"
    cache 'deep'

    input:
      path(dispatcherSource)

    output:
      path('splitBamByZmw')

    script:
    """
    g++ -O2 -std=c++17 ${dispatcherSource} -o splitBamByZmw -lhts -Wl,-rpath,/usr/local/lib
    """
}

/*
  splitBAM: Dispatch all legacy partitions per sample, then index sequentially.
*/
process splitBAM {
    cpus 2
    memory '8 GB'
    time '4h'
    tag { "splitBAM: ${sample_id}" }
    container "${params.hidefseq_container}"

    publishDir path: "${params.analysis_output_dir}",
      mode: params.publication_mode,
      enabled: params.output_intermediate_files,
      saveAs: { filename -> "${dirSplitBAMs(individual_id, sample_id)}/${filename}" }

    afterScript {
      generateAfterScript(
        "${params.analysis_output_dir}/${dirSampleLogs(individual_id, sample_id)}",
        "${task.process}.${params.analysis_id}.${individual_id}.${sample_id}.command.log"
      )
    }

    input:
      tuple val(individual_id), val(sample_id), path(bamFile), path(pbiFile), path(baiFile), path(zmwIdsFile), val(effectiveChunks)
      path(dispatcher)

    output:
      tuple val(individual_id), val(sample_id),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk*.bam"),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk*.bam.pbi"),
      path("${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted.chunk*.bam.bai"),
      val(effectiveChunks)

    script:
    def prefix = "${params.analysis_id}.${individual_id}.${sample_id}.ccs.filtered.aligned.sorted"
    """
    set -euo pipefail
    ./${dispatcher} --input ${bamFile} --ids ${zmwIdsFile} --chunks ${effectiveChunks} \
      --output-prefix ${prefix} --threads ${task.cpus} --max-open-writers 128

    source ${params.conda_base_script}
    conda activate ${params.conda_pbbioconda_env}
    for chunk_bam in ${prefix}.chunk*.bam; do
      pbindex \$chunk_bam
    done
    conda deactivate
    for chunk_bam in ${prefix}.chunk*.bam; do
      ${params.samtools_bin} index -@ ${task.cpus} \$chunk_bam
    done
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
    cache false // Validate the complete immutable cache bundle on every launch.

    afterScript {
      generateAfterScript(
        sharedLogsDir(),
        "${task.process}.command.log"
      )
    }

    input:
      val(config_sig)

    output:
      path("BSgenome_name.txt")

    script:
    def buildCommand = """
    export HIDEF_REFERENCE_BUILD_DIR="\$PWD/library"
    installBSgenome.R -c ${params.paramsFileName}
    """
    cachedBuild(params.prepared_cache.reference, ['library', 'BSgenome_name.txt'], buildCommand)
}

/*
  prepareReferenceSummary: Reusable whole-reference N intervals and per-chromosome context counts.
*/
process prepareReferenceSummary {
    cpus 1
    memory '8 GB'
    time '2h'
    tag { "prepareReferenceSummary" }
    container "${params.hidefseq_container}"
    cache false // Validate the immutable summary bundle on each workflow launch.

    afterScript {
      generateAfterScript(sharedLogsDir(), "${task.process}.command.log")
    }

    input:
      val(BSgenome_name)

    output:
      path("referenceSummary.qs2")

    script:
    def buildCommand = """
    prepareReferenceSummary.R -c ${shellQuote(params.paramsFileName)} -o referenceSummary.qs2
    """
    cachedBuild(params.prepared_cache.reference_summary, ['referenceSummary.qs2'], buildCommand)
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
    cache false // Validate bundle integrity even when external cache files changed.


    afterScript {
      generateAfterScript(
        sharedLogsDir(),
        "${task.process}.command.log"
      )
    }

    output:
      path("${file(params.genome_fasta).name}.bed.gz")
      path("${file(params.genome_fasta).name}.bed.gz.tbi")

    script:
    def buildCommand = """
    #Convert to upper case, replace unsupported bases with N's, extract sequences for all bases
    #(except contig edges), convert to BED format (column 2 is start position of trinucleotide position),
    #and bgzip + tabix index
    ${params.seqkit_bin} seq -u ${params.genome_fasta} | \
      ${params.seqkit_bin} replace -s -p '[^ACGTN]' -r N | \
      ${params.seqkit_bin} sliding -S '' -s1 -W3 | \
      ${params.seqkit_bin} fx2tab -Q | \
      awk -F '[:\\-\\t]' 'BEGIN {OFS="\\t"}{print \$1, \$2, \$2+1, \$4}' | \
      ${params.bgzip_bin} -c > ${file(params.genome_fasta).name}.bed.gz

    ${params.tabix_bin} -@ ${task.cpus} -s 1 -b 2 -e 3 ${file(params.genome_fasta).name}.bed.gz
    """
    cachedBuild(params.prepared_cache.trinucleotides, ["${file(params.genome_fasta).name}.bed.gz", "${file(params.genome_fasta).name}.bed.gz.tbi"], buildCommand)
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
    cache false // Validate bundle integrity even when external cache files changed.


    afterScript {
      generateAfterScript(
        sharedLogsDir(),
        "${task.process}.${individual_id}.command.log"
      )
    }

    input:
      tuple val(individual_id), path(germline_bam_file), val(config_sig)

    output:
      path "${individual_id}.${germline_bam_file}.germline_vcf_variants.qs2"

    script:
    def buildCommand = """
    processGermlineVCFs.R -c ${params.paramsFileName} -i ${individual_id} -o ${individual_id}.${germline_bam_file}.germline_vcf_variants.qs2
    """
    cachedBuild(params.prepared_cache.vcfs[individual_id], ["${individual_id}.${germline_bam_file}.germline_vcf_variants.qs2"], buildCommand)
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
    cache false // Validate bundle integrity even when external cache files changed.


    afterScript {
      generateAfterScript(
        sharedLogsDir(),
        "${task.process}.${germline_bam_file}.command.log"
      )
    }

    input:
      tuple path(germline_bam_file), val(germline_bam_type)

    output:
      path("${germline_bam_file}.bw"), emit: coverage
      path("${germline_bam_file}.vcf.gz*")

    script:
    def buildCommand = """
    #Output per-base coverage using samtools mpileup and direct BAM variant calls using bcftools mpileup
    #Use similar filters for both samtools and bcftools to ensure that samtools coverage data for calculating
    #the fraction of the genome that was filtered maintains correct calculation of the mutation rate.

    #Using samtools mpileup for genome coverage filtering, because the output can be re-formatted into bedgraph.
    #Using bcftools mpileup for call filtering, because the output is in VCF format that is easier to parse.
    #The difference between samtools mpileup and bcftools mpileup should not be significant.
    #We are not calling indels with bcftools mpileup, since that is too noisy to use for filtering.
    #Treat each individual's germline BAM as one sample, even when merged read groups have different SM tags.

    #Slightly different parameters are used for Illumina vs PacBio germline BAM to match how bcftools mpileup is run later
    #in the call filtering analysis.
    
    set -euo pipefail

    if [[ ${germline_bam_type} == Illumina ]]; then
      ${params.samtools_bin} mpileup -A -B -Q 11 -d 999999 --ff 3328 -f ${params.genome_fasta} ${germline_bam_file} 2>/dev/null | awk '{print \$1 "\t" \$2-1 "\t" \$2 "\t" \$4}' > mpileup.bg
      ${params.bcftools_bin} mpileup --ignore-RG -A -B -Q 11 -d 999999 --ns 3328 -I -a "INFO/AD" -f ${params.genome_fasta} -Oz ${germline_bam_file} 2>/dev/null > ${germline_bam_file}.vcf.gz
      ${params.bcftools_bin} index -t ${germline_bam_file}.vcf.gz
    elif [[ ${germline_bam_type} == PacBio ]]; then
      ${params.samtools_bin} mpileup -A -B -Q 5 -d 999999 --ff 3328 -f ${params.genome_fasta} ${germline_bam_file} 2>/dev/null | awk '{print \$1 "\t" \$2-1 "\t" \$2 "\t" \$4}' > mpileup.bg
      ${params.bcftools_bin} mpileup --ignore-RG -A -B -Q 5 -d 999999 --ns 3328 -I -a "INFO/AD" --max-BQ 50 -F0.1 -o25 -e1 -f ${params.genome_fasta} -Oz ${germline_bam_file} 2>/dev/null > ${germline_bam_file}.vcf.gz
      ${params.bcftools_bin} index -t ${germline_bam_file}.vcf.gz
    else
      echo "ERROR: Unknown germline_bam_type: ${germline_bam_type}"
      exit 1
    fi

    sort --parallel=2 -k1,1 -k2,2n mpileup.bg > mpileup.sorted.bg

    ${params.bedGraphToBigWig_bin} mpileup.sorted.bg <(cut -f 1,2 ${params.genome_fai}) ${germline_bam_file}.bw
    """
    cachedBuild(params.prepared_cache.bams[germline_bam_file.name], ["${germline_bam_file}.bw", "${germline_bam_file}.vcf.gz", "${germline_bam_file}.vcf.gz.tbi"], buildCommand)
}

/*
  prepareGermlineCoverageFilters: Prepare each individual's whole-genome
  low-coverage intervals once per distinct threshold, independent of filtergroup.
*/
process prepareGermlineCoverageFilters {
    cpus 2
    memory '16 GB'
    time '4h'
    tag { "prepareGermlineCoverageFilters: ${individual_id} ${threshold}" }
    container "${params.hidefseq_container}"
    cache false

    afterScript {
      generateAfterScript(sharedLogsDir(), "${task.process}.${individual_id}.${threshold}.command.log")
    }

    input:
      tuple val(individual_id), val(threshold), path(coverageFile), val(product)

    output:
      path("${product}")

    script:
    def buildCommand = """
    prepareGermlineCoverageFilters.R --bigwig ${shellQuote(coverageFile)} --fai ${shellQuote(params.genome_fai)} --threshold ${shellQuote(threshold)} --wiggletools ${shellQuote(params.wiggletools_bin)} --wig_to_bigwig ${shellQuote(params.wigToBigWig_bin)} --output ${shellQuote(product)}
    """
    cachedBuild(params.prepared_cache.coverage[product], [product], buildCommand)
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
    cache false // Validate bundle integrity even when external cache files changed.


    afterScript {
      generateAfterScript(
        sharedLogsDir(),
        "${task.process}.${region_filter_file}.bin${binsize}.${threshold}.command.log"
      )
    }

    input:
      tuple path(region_filter_file), val(binsize), val(threshold)

    output:
      path("${region_filter_file}.bin${binsize}.${threshold}.bw")

    script:
    def buildCommand = """
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
    cachedBuild(params.prepared_cache.regions["${region_filter_file}.bin${binsize}.${threshold}.bw".toString()], ["${region_filter_file}.bin${binsize}.${threshold}.bw"], buildCommand)
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

    publishDir path: "${params.analysis_output_dir}",
      mode: params.publication_mode,
      enabled: params.output_intermediate_files,
      saveAs: { filename -> "${dirExtractCalls(individual_id, sample_id)}/${filename}" }

    afterScript {
      generateAfterScript(
        "${params.analysis_output_dir}/${dirSampleLogs(individual_id, sample_id)}",
        "${task.process}.${params.analysis_id}.${individual_id}.${sample_id}.chunk${chunkID}.command.log"
      )
    }

    input:
      tuple val(individual_id), val(sample_id), path(bamFile), path(pbiFile), path(baiFile), val(chunkID), val(effectiveChunks), val(config_sig)

    output:
      tuple val(individual_id), val(sample_id), path("${params.analysis_id}.${individual_id}.${sample_id}.extractCalls.chunk${chunkID}.qs2"), val(chunkID), val(effectiveChunks)

    script:
    """
    extractCalls.R -c ${params.paramsFileName} -b ${bamFile} -s ${sample_id} -o ${params.analysis_id}.${individual_id}.${sample_id}.extractCalls.chunk${chunkID}.qs2
    """
}

/*
  filterCallsChunkChromgroupFiltergroup: Run filterCalls.R for each analysis chunk, chromgroup, filtergroup combination
*/
process filterCallsChunkChromgroupFiltergroup {
    cpus 1
    memory {
      def baseMemory = params.mem_filterCallsChunkChromgroupFiltergroup as nextflow.util.MemoryUnit
      baseMemory * (1 + 0.5*(task.attempt - 1))
    }
    time {
        def baseTime = params.time_filterCallsChunkChromgroupFiltergroup as nextflow.util.Duration
        baseTime * (1 + (task.attempt - 1))
    }
    maxRetries params.maxRetries_filterCallsChunkChromgroupFiltergroup
    tag { "filterCallsChunkChromgroupFiltergroup: ${sample_id} -> chunk ${chunkID}" }
    container "${params.hidefseq_container}"

    publishDir path: "${params.analysis_output_dir}",
      mode: params.publication_mode,
      enabled: params.output_intermediate_files,
      saveAs: { filename -> "${dirFilterCalls(individual_id, sample_id)}/${filename}" }

    afterScript {
      generateAfterScript(
        "${params.analysis_output_dir}/${dirSampleLogs(individual_id, sample_id)}",
        "${task.process}.${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.chunk${chunkID}.command.log"
      )
    }

    input:
      tuple val(individual_id), val(sample_id), path(extractCallsFile), val(chunkID), val(effectiveChunks), val(chromgroup), val(filtergroup), val(config_sig)

    output:
      tuple val(individual_id), val(sample_id), val(chromgroup), val(filtergroup), val(chunkID), path("${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.filterCalls.chunk${chunkID}.qs2"), val(effectiveChunks)

    script:
    """
    filterCalls.R -c ${params.paramsFileName} -s ${sample_id} -g ${chromgroup} -v ${filtergroup} -f ${extractCallsFile} -o ${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.filterCalls.chunk${chunkID}.qs2
    """
}

/* Compile once; explicit source content identity is confined to burden tasks. */
process compileCoverageAnnotator {
    cpus 1
    memory '1 GB'
    time '10m'
    container "${params.hidefseq_container}"
    cache 'deep'

    input:
      path(annotatorSource)

    output:
      path('annotateCoverage')

    script:
    """
    g++ -O3 -std=c++17 -Wall -Wextra ${annotatorSource} -o annotateCoverage -lhts -Wl,-rpath,/usr/local/lib
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

    publishDir path: "${params.analysis_output_dir}",
      mode: params.publication_mode,
      pattern: "*.calculateBurdens.qs2",
      enabled: params.output_intermediate_files,
      saveAs: { filename -> "${dirCalculateBurdens(individual_id, sample_id)}/${filename}" }

    publishDir path: "${params.analysis_output_dir}",
      mode: params.publication_mode,
      pattern: "*.bed.gz*",
      saveAs: { filename -> "${dirCoverage_Reftnc(individual_id, sample_id)}/${chromgroup}/${filename}" }

    afterScript {
      generateAfterScript(
        "${params.analysis_output_dir}/${dirSampleLogs(individual_id, sample_id)}",
        "${task.process}.${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.command.log"
      )
    }

    input:
      tuple val(individual_id), val(sample_id), val(chromgroup), val(filtergroup), path(filterCallsFiles), val(config_sig)
      path(coverageAnnotator)

    output:
      tuple val(individual_id), val(sample_id), val(chromgroup), val(filtergroup), path("${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.calculateBurdens.qs2"), emit: tuple_qs2
      tuple val(individual_id), val(sample_id), val(chromgroup), val(filtergroup), path("*.bed.gz"), path("*.bed.gz.tbi"), emit: coverage_reftnc

    script:
    """
    calculateBurdens.R --coverage-annotator './${coverageAnnotator}' -c ${params.paramsFileName} -s ${sample_id} -g ${chromgroup} -v ${filtergroup} -f ${filterCallsFiles.join(',')} -o ${params.analysis_id}.${individual_id}.${sample_id}.${chromgroup}.${filtergroup}.calculateBurdens.qs2
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

    publishDir path: "${params.analysis_output_dir}",
      mode: params.publication_mode,
      saveAs: { filename -> "${sampleBaseDir(individual_id, sample_id)}/${filename}" }

    afterScript {
      generateAfterScript(
        "${params.analysis_output_dir}/${dirSampleLogs(individual_id, sample_id)}",
        "${task.process}.${params.analysis_id}.${individual_id}.${sample_id}.command.log"
      )
    }

    input:
      tuple val(individual_id), val(sample_id), path(calculateBurdensFiles), val(config_sig)

    output:
      tuple val(individual_id), val(sample_id), emit: out_ch
      path("${params.analysis_id}.${individual_id}.${sample_id}.outputResults.qs2")
      path("${params.analysis_id}.${individual_id}.${sample_id}.yaml_config.tsv")
      path("${params.analysis_id}.${individual_id}.${sample_id}.run_metadata.tsv")
      path("*/**/*.{tsv,vcf.bgz,vcf.bgz.tbi,pdf}")

    script:
    """
    outputResults.R -c ${params.paramsFileName} -s ${sample_id} -f ${calculateBurdensFiles.join(',')} -o ${params.analysis_id}.${individual_id}.${sample_id}
    """
}
