#!/usr/bin/env python3
"""Exercise production merge-order channels with Nextflow on an allocated node.

Run: python3 tests/test_merge_input_order.py --directory NEW_FIXTURE_DIR
Requires Nextflow on PATH. No scientific workers or BAM contents are used.
The demultiplex checks exercise grouping/order only, not the upstream joins.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess


def fixture_spec(directory):
    runs = [{'run_id': name, 'samples': []} for name in ('run-z', 'run-a')]
    merge_cases, merge_expected = [], {}
    demux_cases = {'round1': [], 'round2': []}
    demux_expected = {}
    for style in ('round1', 'round2'):
        # YAML run and barcode order deliberately disagree with lexical order.
        for case, permutation in enumerate(itertools.permutations(range(4))):
            sample = f'{style}-sample-{case}'
            rows, bams, pbis = [], [], []
            for index in range(4):
                run = runs[index // 2]
                barcode = f'bc{"Z" if index % 2 == 0 else "A"}_{sample}'
                config = {'sample_id': sample, 'barcode_ids': barcode}
                if style == 'round2':
                    config['barcode_ids_round2'] = f'inner_{index}'
                    barcode += '.' + config['barcode_ids_round2']
                run['samples'].append(config)
                bam = str(directory / 'inputs' / sample / f'reverse-{3-index}.bam')
                pbi = bam + '.pbi'
                rows.append([run['run_id'], 'individual', sample, barcode, bam, pbi])
                bams.append(bam)
                pbis.append(pbi)
            merge_cases.append([rows[index] for index in permutation])
            merge_expected[sample] = ['individual', sample, bams, pbis]

        sample = f'{style}-single'
        config = {'sample_id': sample, 'barcode_ids': sample}
        barcode = sample
        if style == 'round2':
            config['barcode_ids_round2'] = 'inner_single'
            barcode += '.inner_single'
        runs[0]['samples'].append(config)
        bam = str(directory / 'inputs' / sample / 'only.bam')
        merge_cases.append([['run-z', 'individual', sample, barcode, bam, bam + '.pbi']])
        merge_expected[sample] = ['individual', sample, [bam], [bam + '.pbi']]

        for case, permutation in enumerate(itertools.permutations(range(3))):
            sample = f'{style}-demux-{case}'
            barcode = 'outer' if style == 'round1' else 'outer.inner'
            # Parent-directory order also disagrees with basename order.
            paths = [str(directory / 'inputs' / sample / f'parent-{2-index}' / name)
                     for index, name in enumerate(('a.demux.bam', 'm.demux.bam', 'z.demux.bam'))]
            rows = [['run-z', 'individual', sample, barcode, path] for path in paths]
            demux_cases[style].append([rows[index] for index in permutation])
            demux_expected[sample] = ['run-z', 'individual', sample, barcode, paths]
        sample = f'{style}-demux-single'
        barcode = 'outer' if style == 'round1' else 'outer.inner'
        path = str(directory / 'inputs' / sample / 'only.bam')
        demux_cases[style].append([['run-z', 'individual', sample, barcode, path]])
        demux_expected[sample] = ['run-z', 'individual', sample, barcode, [path]]

    def interleave(cases):
        return [rows[index] for index in range(max(map(len, cases)))
                for rows in cases if index < len(rows)]

    return {'runs': runs, 'merge_arrivals': interleave(merge_cases),
            'merge_expected': merge_expected,
            'demux_arrivals': {style: interleave(cases) for style, cases in demux_cases.items()},
            'demux_expected': demux_expected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--directory', type=Path, required=True)
    args = parser.parse_args()
    directory = args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    source_path = args.repo / 'main.nf'
    source_bytes = source_path.read_bytes()
    source = source_bytes.decode()
    start = source.index('  run_sample_order = [:]')
    end = source.index('\n  // Run process', start)
    merge_fragment = source[start:end]
    assert merge_fragment.count('pbmm2Align.out') == 1
    (directory / 'production-sample-merge.nf.txt').write_text(merge_fragment)
    fragments = [merge_fragment.replace('pbmm2Align.out', 'synthetic_alignments')]
    for style in ('round1', 'round2'):
        name = f'mergeDemuxBams_{style}_input_ch'
        start = source.index(f'  {name} = ')
        end = source.index('\n\n', start)
        block = source[start:end]
        ordering = block[block.index('    .groupTuple('):]
        (directory / f'production-demux-{style}.nf.txt').write_text(ordering)
        fragments.append(f'  {name} = synthetic_demux_{style}\n' + ordering)

    spec = fixture_spec(directory)
    (directory / 'spec.json').write_text(json.dumps(spec, indent=2) + '\n')
    prefix = r'''nextflow.enable.dsl=2
params.fixture_spec = null
params.fixture_output = null

workflow {
  def spec = new groovy.json.JsonSlurper().parseText(file(params.fixture_spec).text)
  params.runs = spec.runs
  synthetic_alignments = Channel.fromList(spec.merge_arrivals).map { row ->
    // Reproduce production's round2 interpolation; Nextflow may coerce its type.
    def barcode = row[3]
    if (row[2].startsWith('round2-')) {
      def parts = row[3].tokenize('.')
      barcode = "${parts[0]}.${parts[1]}"
    }
    tuple(row[0], row[1], row[2], barcode, file(row[4]), file(row[5]))
  }
  synthetic_demux_round1 = Channel.fromList(spec.demux_arrivals.round1).map { row ->
    tuple(row[0], row[1], row[2], row[3], file(row[4]))
  }
  synthetic_demux_round2 = Channel.fromList(spec.demux_arrivals.round2).map { row ->
    tuple(row[0], row[1], row[2], row[3], file(row[4]))
  }
'''
    suffix = r'''
  merge_checks = mergeAlignedSampleBAMs_input_ch.map { individual, sample, bams, pbis ->
    def observed = [individual, sample, bams.collect { it.toString() }, pbis.collect { it.toString() }]
    assert observed == spec.merge_expected[sample]
    [kind: 'sample', sample: sample, observed: observed]
  }
  demux_checks = mergeDemuxBams_round1_input_ch.mix(mergeDemuxBams_round2_input_ch)
    .map { run, individual, sample, barcode, bams ->
      def observed = [run, individual, sample, barcode, bams.collect { it.toString() }]
      assert observed == spec.demux_expected[sample]
      [kind: 'demux', sample: sample, observed: observed]
    }
  merge_checks.mix(demux_checks).collect().view { rows ->
    assert rows.size() == 64
    assert rows.collect { it.sample }.toSet().size() == 64
    file(params.fixture_output).text = groovy.json.JsonOutput.prettyPrint(
      groovy.json.JsonOutput.toJson(rows)) + '\n'
    "Verified 50 sample merge groups and 14 demultiplex merge groups"
  }
}
'''
    fixture = directory / 'fixture.nf'
    fixture.write_text(prefix + '\n'.join(fragments) + suffix)
    env = dict(os.environ, NXF_HOME=str(directory / 'nxf-home'),
               NXF_ASSETS=str(directory / 'assets'), NXF_ANSI_LOG='false',
               NXF_OPTS='-Xmx2g')
    subprocess.run([
        'nextflow', '-log', str(directory / 'nextflow.log'), 'run', str(fixture),
        '-work-dir', str(directory / 'work'),
        '--fixture_spec', str(directory / 'spec.json'),
        '--fixture_output', str(directory / 'observed.json'),
    ], cwd=directory, env=env, check=True)
    observations = json.loads((directory / 'observed.json').read_text())
    assert len(observations) == 64
    assert len({row['sample'] for row in observations}) == 64
    for row in observations:
        expected = spec['merge_expected' if row['kind'] == 'sample' else 'demux_expected']
        assert row['observed'] == expected[row['sample']]
    assert source_path.read_bytes() == source_bytes, 'main.nf changed during fixture'
    report = {'status': 'pass', 'sample_merge_groups': 50, 'demux_merge_groups': 14,
              'sample_arrival_permutations_per_round': 24,
              'demux_arrival_permutations_per_round': 6,
              'single_input_groups': 4, 'upstream_demux_joins_exercised': False,
              'main_nf_sha256': hashlib.sha256(source_bytes).hexdigest()}
    (directory / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
