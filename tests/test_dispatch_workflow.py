#!/usr/bin/env python3
"""Generate/verify a small Nextflow replay of the actual dispatch integration.

Generate inside the pinned container with PacBio tools active, run the emitted
fixture.nf with Nextflow on compute, then verify in that container. The fixture
extracts production process definitions, channel wiring and tuple expansion.
"""
import argparse
import json
from pathlib import Path
import subprocess
from test_split_bam_dispatch import fixture, record, run


def generate(args):
    directory = Path(args.directory).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    repo = Path(args.repo).resolve()
    source = (repo / 'main.nf').read_text()
    helpers = source[source.index('//Define output directories'):source.index('def signatureYaml()')]
    helpers += source[source.index('def nearestExistingPath('):source.index('def cachedBuild(')]
    processes = source[source.index('process countAnalysisZMWs {'):source.index('/*\n  installBSgenome:')]
    begin = source.index('  countAnalysisZMWs(mergeAlignedSampleBAMs.out)')
    end = source.index('  //******************\n  // installBSgenome', begin)
    wiring = source[begin:end].replace('mergeAlignedSampleBAMs.out', 'inputs').replace('${projectDir}/bin/', '${params.repo_dir}/bin/')
    inputs = directory / 'inputs'
    inputs.mkdir()
    cases = {'multi': [record(movie, hole, hole * 5) for hole in range(1, 17) for movie in ('movieA', 'movieB')],
             'single': [record('movieA', 7, 5)], 'empty': []}
    for name, records in cases.items():
        bam, _ = fixture(args.samtools, inputs, name, records)
        run([args.pbindex, bam])
        run([args.samtools, 'index', bam])
    workflow = '''workflow {
  inputs = channel.of('multi', 'single', 'empty').map { sample ->
    tuple('ind', sample, file("${params.fixture_inputs}/${sample}.bam"),
      file("${params.fixture_inputs}/${sample}.bam.pbi"), file("${params.fixture_inputs}/${sample}.bam.bai"))
  }
''' + wiring + '''
  splitBAM_chunks_ch.map { individual, sample, bam, pbi, bai, chunk, effective ->
    [individual, sample, bam, pbi, bai, chunk, effective].join('\\t')
  }.collectFile(name: 'tuples.tsv', storeDir: params.fixture_results, newLine: true)
}
'''
    (directory / 'fixture.nf').write_text(helpers + '\n' + processes + '\n' + workflow)
    configuration = dict(analysis_id='dispatchfixture', analysis_chunks=12, repo_dir=str(repo),
        fixture_inputs=str(inputs), fixture_results=str(directory / 'results'),
        analysis_output_dir=str(directory / 'published'), output_intermediate_files=True,
        hidefseq_container=args.container,
        conda_base_script=args.conda_base_script, conda_pbbioconda_env=args.conda_env,
        samtools_bin=args.samtools)
    (directory / 'params.json').write_text(json.dumps(configuration, indent=2) + '\n')
    print('Generated:', directory / 'fixture.nf')


def verify(args):
    directory = Path(args.directory).resolve()
    rows = [line.split('\t') for line in (directory / 'results/tuples.tsv').read_text().splitlines()]
    samples = {}
    for row in rows:
        if len(row) != 7 or row[0] != 'ind':
            raise AssertionError('Downstream tuple schema changed')
        samples.setdefault(row[1], []).append(row)
    if set(samples) != {'multi', 'single'}:
        raise AssertionError('Empty input not skipped or nonempty sample missing')
    report = {}
    for sample, expected_chunks in [('multi', 12), ('single', 1)]:
        selected = sorted(samples[sample], key=lambda row: int(row[5]))
        if [int(row[5]) for row in selected] != list(range(1, expected_chunks + 1)):
            raise AssertionError('Missing/duplicate chunk IDs')
        original = directory / 'inputs' / (sample + '.bam')
        ids = run([args.zmwfilter, '--show-all', original]).stdout.splitlines()
        quotient, remainder = divmod(len(ids), expected_chunks)
        offset = 0
        for row in selected:
            chunk = int(row[5])
            if int(row[6]) != expected_chunks:
                raise AssertionError('Effective chunk count changed')
            basename = f'dispatchfixture.ind.{sample}.ccs.filtered.aligned.sorted.chunk{chunk}.bam'
            if [Path(item).name for item in row[2:5]] != [basename, basename + '.pbi', basename + '.bai']:
                raise AssertionError('Scientific chunk filenames changed')
            if not all(Path(item).is_file() for item in row[2:5]):
                raise AssertionError('Output/index missing from work directory')
            published = directory / 'published' / f'dispatchfixture.ind.{sample}' / 'splitBAMs'
            for item in row[2:5]:
                work_file = Path(item)
                published_file = published / work_file.name
                if not published_file.is_file() or published_file.read_bytes() != work_file.read_bytes():
                    raise AssertionError('Published output missing/different; work output must also remain')
            size = quotient + (chunk <= remainder)
            include = directory / f'{sample}.{chunk}.include.txt'
            include.write_text('\n'.join(ids[offset:offset + size]) + '\n')
            offset += size
            legacy = directory / f'{sample}.{chunk}.legacy.bam'
            run([args.zmwfilter, '--include', include, original, legacy])
            if run([args.samtools, 'view', '--no-PG', legacy]).stdout != run([args.samtools, 'view', '--no-PG', row[2]]).stdout:
                raise AssertionError('Record/tag/order difference in ' + basename)
            for line in run([args.samtools, 'view', '--no-PG', '-H', original]).stdout.splitlines():
                if line not in run([args.samtools, 'view', '--no-PG', '-H', row[2]]).stdout.splitlines():
                    raise AssertionError('Original header line missing')
            run([args.samtools, 'quickcheck', row[2]])
            run([args.samtools, 'idxstats', row[2]])
        report[sample] = dict(chunks=expected_chunks, status='exact records, headers, indices and tuples')
    (directory / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['generate', 'verify'])
    parser.add_argument('--directory', required=True)
    parser.add_argument('--repo', default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--samtools', default='samtools')
    parser.add_argument('--pbindex', default='pbindex')
    parser.add_argument('--zmwfilter', default='zmwfilter')
    parser.add_argument('--container')
    parser.add_argument('--conda-base-script')
    parser.add_argument('--conda-env')
    args = parser.parse_args()
    if args.mode == 'generate' and not all((args.container, args.conda_base_script, args.conda_env)):
        parser.error('generate requires container, conda-base-script and conda-env')
    (generate if args.mode == 'generate' else verify)(args)


if __name__ == '__main__':
    main()
