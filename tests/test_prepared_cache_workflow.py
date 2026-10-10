#!/usr/bin/env python3
"""Actual Nextflow prepared-cache integration on tiny tracks; compute nodes only.

Generate in the pinned container, then exercise with the host Nextflow launcher.
All cache mutation is confined to the newly created fixture directory.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import tempfile


def generate(args):
    directory = Path(args.directory).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    repo = Path(args.repo).resolve()
    inputs = directory / 'inputs'
    inputs.mkdir()
    (inputs / 'tiny.fa.fai').write_text('chr1\t200\t0\t200\t201\nchrM\t80\t0\t80\t81\n')
    (inputs / 'chrom.sizes').write_text('chr1\t200\nchrM\t80\n')
    (inputs / 'tiny.bg').write_text('chr1\t0\t20\t2\nchr1\t20\t100\t20\nchrM\t10\t40\t10\n')
    subprocess.run([args.bedgraph_to_bigwig, str(inputs / 'tiny.bg'), str(inputs / 'chrom.sizes'), str(inputs / 'coverage.bw')], check=True)
    for name in ('regionA.bw', 'regionB.bw'):
        shutil.copy2(inputs / 'coverage.bw', inputs / name)
    wrapper = directory / 'wiggletools-counted'
    wrapper.write_text('#!/usr/bin/env python3\nimport json,os,sys,time\n' +
        'log=' + repr(str(directory / 'build-events.jsonl')) + '\n' +
        'line=json.dumps(dict(cwd=os.getcwd(),args=sys.argv[1:],time=time.time(),pid=os.getpid()))+"\\n"\n' +
        'fd=os.open(log,os.O_WRONLY|os.O_APPEND|os.O_CREAT,0o600)\nos.write(fd,line.encode())\nos.close(fd)\n' +
        'if any("regionB.bw" in arg for arg in sys.argv[1:]): time.sleep(2)\n' +
        'os.execv(' + repr(args.wiggletools) + ',[' + repr(args.wiggletools) + ']+sys.argv[1:])\n')
    wrapper.chmod(0o755)
    (directory / 'bin').symlink_to(repo / 'bin', target_is_directory=True)
    source = (repo / 'main.nf').read_text()
    helpers = source[source.index('//Define output directories'):source.index('def signatureYaml()')]
    helpers += source[source.index('def fileSha256('):source.index('def canonicalBarcodePair(')]
    processes = source[source.index('process prepareGermlineCoverageFilters {'):source.index('/*\n  extractCallsChunk:')]
    identity = source[source.index('  inputIdentities = [:]'):source.index('  referenceEntry = makeCacheEntry.call(')]
    identity = identity.replace('${workflow.projectDir}', '${params.repo_dir}')
    barrier = source[source.index('  prepareFilters_done = BSgenome_name_ch'):source.index('  // Create input channel\n  extractCalls_input_ch')]
    barrier = barrier.replace('BSgenome_name_ch', "channel.value('fixture-reference-ready')")
    for name in ('prepareReferenceSummary', 'processGermlineVCFs', 'processGermlineBAMs'):
        barrier = barrier.replace(name + '.out', "channel.value('fixture-" + name + "-ready')")
    audit_process = '''process auditPreparedCache {
  cpus 1
  memory '512 MB'
  time '5m'
  container "${params.hidefseq_container}"
  cache false
  input:
    tuple val(consumer), val(ready)
    path(indexFile)
  output:
    path("consumer${consumer}.json")
  script:
    """
    python3 ${params.fixture_driver} audit --index ${indexFile} --output consumer${consumer}.json
    """
}
'''
    workflow = '''workflow {
''' + identity + '''
  regionEntries = [:]
  coverageEntries = [:]
  entries = []
  regionInputs = ['regionA.bw', 'regionB.bw'].collect { name ->
    def sourceFile = file("${params.fixture_inputs}/${name}")
    def product = "${name}.bin1.gte0.1.bw".toString()
    def entry = makeCacheEntry.call('region-filter', 'prepareRegionFilters',
      [binsize: 1, threshold: 'gte0.1', wiggletools: params.wiggletools_bin, wigToBigWig: params.wigToBigWig_bin],
      [region: sourceFile, fai: params.genome_fai], [], [product])
    regionEntries[product] = entry
    entries << [product: product, entry: entry]
    tuple(sourceFile, 1, 'gte0.1')
  }
  // Two individuals share one coverage input; duplicated filtergroup thresholds
  // still prepare only one product per individual/unique numeric threshold.
  rawIdentity = preparedInputIdentity("${params.fixture_inputs}/coverage.bw")
  thresholds = [1, 15, 0.1, 1].unique()
  coverageInputs = ['a', 'b'].collectMany { individual ->
    thresholds.collect { threshold ->
      def product = "${individual}.shared.bam.minCoverage${threshold}.qs2".toString()
      def entry = makeCacheEntry.call('germline-coverage-filter', 'prepareGermlineCoverageFilters',
        [individual: individual, threshold: threshold, rawCoverage: rawIdentity,
         wiggletools: params.wiggletools_bin, wigToBigWig: params.wigToBigWig_bin],
        [fai: params.genome_fai], ['prepareGermlineCoverageFilters.R', 'sharedFunctions.R'], [product])
      coverageEntries[product] = entry
      entries << [product: product, entry: entry]
      tuple(individual, threshold, file("${params.fixture_inputs}/coverage.bw"), product)
    }
  }
  params.prepared_cache = [regions: regionEntries, coverage: coverageEntries]
  file(params.fixture_index).text = groovy.json.JsonOutput.toJson(
    [root: params.prepared_cache_root, helper: params.cache_helper, container: params.hidefseq_container, entries: entries])
  prepareRegionFilters(channel.fromList(regionInputs))
  prepareGermlineCoverageFilters(channel.fromList(coverageInputs))
''' + barrier + '''
  consumers = channel.of(1, 2).combine(prepareFilters_done)
  auditPreparedCache(consumers, channel.value(file(params.fixture_index)))
  auditPreparedCache.out.collectFile(name: 'receipts.jsonl', storeDir: params.fixture_results, newLine: true)
}
'''
    (directory / 'fixture.nf').write_text(helpers + processes + audit_process + workflow)
    params = dict(analysis_id='cachefixture', analysis_output_dir=str(directory / 'published'),
        repo_dir=str(repo), fixture_inputs=str(inputs), genome_fai=str(inputs / 'tiny.fa.fai'),
        cache_helper=str(repo / 'bin/artifactCache.R'), prepared_cache_root=str(directory / 'cache'),
        hidefseq_container=args.container, wiggletools_bin=str(wrapper), wigToBigWig_bin=args.wig_to_bigwig,
        fixture_driver=str(Path(__file__).resolve()))
    (directory / 'base-params.json').write_text(json.dumps(params, indent=2) + '\n')
    (directory / 'fixture.config').write_text("includeConfig '" + args.cluster_config + "'\n" + '''process {
  withName: prepareRegionFilters { cpus = 1; memory = '1 GB'; time = '5m' }
  withName: prepareGermlineCoverageFilters { cpus = 1; memory = '4 GB'; time = '5m' }
}
''')
    (directory / 'compare.R').write_text('''suppressPackageStartupMessages({library(qs2); library(rtracklayer)})
a <- commandArgs(TRUE)
pairs <- read.delim(a[[1]], stringsAsFactors=FALSE)
for(i in seq_len(nrow(pairs))) {
  load <- if(grepl("\\\\.qs2$", pairs$expected[[i]])) qs2::qs_read else rtracklayer::import
  stopifnot(identical(load(pairs$expected[[i]]), load(pairs$observed[[i]])))
}
cat("Scientific objects identical:", nrow(pairs), "\\n")
''')
    print('Generated', directory)


def audit(args):
    index = json.loads(Path(args.index).read_text())
    rscript = ['Rscript', '--vanilla']
    if shutil.which('Rscript') is None:
        rscript = ['apptainer', 'exec', '--cleanenv', '-B', '/projects', index['container'], *rscript]
    paths = []
    for item in index['entries']:
        identity = json.loads(item['entry']['json'])
        with tempfile.TemporaryDirectory(dir=Path(args.index).resolve().parent) as temporary:
            identity_file = Path(temporary) / 'identity.json'
            identity_file.write_text(json.dumps(identity))
            destination = Path(subprocess.check_output([*rscript, index['helper'], 'verify',
                '--root', index['root'], '--identity', str(identity_file)], text=True).strip())
        assert destination == Path(item['entry']['directory'])
        assert (destination / item['product']).is_file()
        paths.append(str(destination / item['product']))
    assert len(paths) == 8 and len(set(paths)) == 8
    Path(args.output).write_text(json.dumps(dict(manifests_verified=len(paths), files=paths, time=time.time())) + '\n')


def events(directory):
    path = directory / 'build-events.jsonl'
    return [] if not path.exists() else [json.loads(line) for line in path.read_text().splitlines()]


def manifests(index):
    result = {}
    for item in index['entries']:
        path = Path(item['entry']['directory']) / 'manifest.complete.json'
        stat = path.stat()
        result[item['product']] = dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            inode=stat.st_ino, mtime_ns=stat.st_mtime_ns, directory=str(path.parent))
    return result


def phase_command(directory, phase, cache, work, resume=False):
    phase_dir = directory / phase
    phase_dir.mkdir()
    params = json.loads((directory / 'base-params.json').read_text())
    params.update(prepared_cache_root=str(cache), fixture_index=str(phase_dir / 'index.json'),
                  fixture_results=str(phase_dir / 'results'), analysis_output_dir=str(phase_dir / 'published'))
    (phase_dir / 'params.json').write_text(json.dumps(params, indent=2) + '\n')
    command = ['nextflow', '-config', str(directory / 'fixture.config'), 'run', str(directory / 'fixture.nf'),
        '-params-file', str(phase_dir / 'params.json'), '-work-dir', str(work),
        '-with-trace', str(phase_dir / 'trace.tsv')]
    if resume:
        command.append('-resume')
    return phase_dir, command


def check_phase(directory, phase, baseline=None):
    location = directory / phase
    rows = list(csv.DictReader((location / 'trace.tsv').open(), delimiter='\t'))
    assert len(rows) == 10 and all(row['status'] == 'COMPLETED' and row['exit'] == '0' for row in rows), rows
    index = json.loads((location / 'index.json').read_text())
    class Audit: pass
    args = Audit()
    args.index, args.output = location / 'index.json', location / 'verified.json'
    audit(args)
    receipts = [line for line in (location / 'results/receipts.jsonl').read_text().splitlines() if line.strip()]
    assert len(receipts) == 2
    assert all(json.loads(line)['manifests_verified'] == 8 for line in receipts)
    if baseline:
        by_name = {item['product']: item for item in baseline['entries']}
        pairs = location / 'scientific-pairs.tsv'
        pairs.write_text('expected\tobserved\n' + ''.join(
            str(Path(by_name[item['product']]['entry']['directory']) / item['product']) + '\t' +
            str(Path(item['entry']['directory']) / item['product']) + '\n' for item in index['entries']))
        container = json.loads((directory / 'base-params.json').read_text())['hidefseq_container']
        with (location / 'scientific.log').open('w') as output:
            subprocess.run(['apptainer', 'exec', '--cleanenv', '-B', '/projects', container,
                'Rscript', '--vanilla', str(directory / 'compare.R'), str(pairs)], check=True, stdout=output, stderr=subprocess.STDOUT)
    return index


def exercise(args):
    directory = Path(args.directory).resolve()
    launch = directory / 'launch-main'
    launch.mkdir()
    report = {}
    def sequential(phase, resume=False):
        location, command = phase_command(directory, phase, directory / 'cache', directory / 'work-main', resume)
        before = len(events(directory))
        with (location / 'nextflow.log').open('w') as output:
            subprocess.run(command, cwd=launch, stdout=output, stderr=subprocess.STDOUT, check=True)
        built = len(events(directory)) - before
        report[phase] = dict(actual_builds=built)
        return check_phase(directory, phase), built
    cold, built = sequential('cold')
    assert built == 8, built
    initial = manifests(cold)
    # Keep independent scientific snapshots, including the entry later renamed.
    baseline = json.loads(json.dumps(cold))
    snapshot = directory / 'scientific-baseline'
    snapshot.mkdir()
    for item in baseline['entries']:
        source = Path(item['entry']['directory']) / item['product']
        shutil.copy2(source, snapshot / item['product'])
        item['entry']['directory'] = str(snapshot)
    warm, built = sequential('warm-resume', True)
    assert built == 0 and manifests(warm) == initial
    check_phase(directory, 'warm-resume', baseline)
    missing = next(item for item in cold['entries'] if item['product'].startswith('regionA.'))
    Path(missing['entry']['directory']).rename(directory / 'retired-own-test-entry')
    rebuilt, built = sequential('missing-resume', True)
    assert built == 1, built
    after = manifests(rebuilt)
    assert all(after[name] == original for name, original in initial.items() if name != missing['product'])
    check_phase(directory, 'missing-resume', baseline)
    before = len(events(directory))
    concurrent = []
    for phase in ('concurrent-a', 'concurrent-b'):
        location, command = phase_command(directory, phase, directory / 'concurrent-cache', directory / ('work-' + phase))
        cwd = directory / ('launch-' + phase)
        cwd.mkdir()
        handle = (location / 'nextflow.log').open('w')
        concurrent.append((phase, subprocess.Popen(command, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT), handle))
    for phase, process, handle in concurrent:
        status = process.wait()
        handle.close()
        assert status == 0, (phase, status)
        check_phase(directory, phase, baseline)
    built = len(events(directory)) - before
    assert built == 8, built
    report['concurrent'] = dict(actual_builds=built, workflows=2, expected_without_lock=16)
    report['status'] = 'cold/warm/missing/concurrent cache and global barrier checks passed'
    (directory / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['generate', 'exercise', 'audit'])
    parser.add_argument('--directory')
    parser.add_argument('--repo', default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--container')
    parser.add_argument('--cluster-config')
    parser.add_argument('--wiggletools', default='/hidef/bin/WiggleTools/bin/wiggletools')
    parser.add_argument('--wig-to-bigwig', default='/hidef/bin/wigToBigWig')
    parser.add_argument('--bedgraph-to-bigwig', default='/hidef/bin/bedGraphToBigWig')
    parser.add_argument('--index')
    parser.add_argument('--output')
    args = parser.parse_args()
    if args.mode == 'generate' and not all((args.directory, args.container, args.cluster_config)):
        parser.error('generate requires directory, container, cluster-config')
    if args.mode == 'exercise' and not args.directory:
        parser.error('exercise requires directory')
    if args.mode == 'audit' and not all((args.index, args.output)):
        parser.error('audit requires index and output')
    {'generate': generate, 'exercise': exercise, 'audit': audit}[args.mode](args)


if __name__ == '__main__':
    main()
