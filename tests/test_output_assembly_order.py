#!/usr/bin/env python3
"""Check actual output channel ordering with Nextflow on an allocated node.

Run: python3 tests/test_output_assembly_order.py --directory NEW_FIXTURE_DIR
Requires Nextflow on PATH. No scientific workers or large inputs are launched.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--directory', type=Path, required=True)
    args = parser.parse_args()
    directory = args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=False)

    source_path = args.repo / 'main.nf'
    source = source_path.read_text()
    start = source.index('  outputResultsSample_input_ch = ')
    end = source.index('\n  // Run process', start)
    fragment = source[start:end]
    (directory / 'production-channel.nf.txt').write_text(fragment)
    fixture_fragment = fragment.replace(
        'calculateBurdensChromgroupFiltergroup.out.tuple_qs2', 'synthetic_inputs')
    assert fixture_fragment != fragment

    prefix = r'''nextflow.enable.dsl=2
params.chromgroups = [[chromgroup: 'nuclear-z'], [chromgroup: 'mito-a']]
params.filtergroups = [[filtergroup: 'strict-z'], [filtergroup: 'lenient-a']]
params.fixture_root = null
params.fixture_output = null

workflow {
  chromgroups_filtergroups_list = params.chromgroups.collectMany { cg ->
    params.filtergroups.collect { fg -> tuple(cg.chromgroup, fg.filtergroup) }
  }
  config_signatures = [outputResultsSample: 'fixture-signature']
  def permutations = chromgroups_filtergroups_list.permutations().toList()
  assert permutations.size() == 24
  def cases = permutations.withIndex().collect { groups, index ->
    def sample = "sample-${index}".toString()
    groups.collect { group ->
      def inputNumber = 3 - chromgroups_filtergroups_list.indexOf(group)
      tuple('individual', sample, group[0], group[1],
            file("${params.fixture_root}/${sample}/input-${inputNumber}.qs2"))
    }
  }
  // Interleave samples and deliberately reverse filenames relative to config.
  def arrivals = (0..3).collectMany { round -> cases.collect { rows -> rows[round] } }
  synthetic_inputs = Channel.fromList(arrivals)
'''
    suffix = r'''
  outputResultsSample_input_ch
    .map { individual, sample, inputs, signature ->
      assert individual == 'individual'
      assert signature == 'fixture-signature'
      def expected = (0..3).collect { i ->
        file("${params.fixture_root}/${sample}/input-${3-i}.qs2").toString()
      }
      def observed = inputs.collect { it.toString() }
      assert observed == expected
      [sample: sample, expected: expected, observed: observed]
    }
    .collect()
    .view { rows ->
      assert rows.size() == 24
      assert rows.collect { it.sample }.toSet().size() == 24
      file(params.fixture_output).text = groovy.json.JsonOutput.prettyPrint(
        groovy.json.JsonOutput.toJson([cases: rows])) + '\n'
      "Verified ${rows.size()} assembly permutations"
    }
}
'''
    fixture = directory / 'fixture.nf'
    fixture.write_text(prefix + fixture_fragment + suffix)
    env = dict(os.environ, NXF_HOME=str(directory / 'nxf-home'),
               NXF_ASSETS=str(directory / 'assets'), NXF_ANSI_LOG='false',
               NXF_OPTS='-Xmx2g')
    subprocess.run([
        'nextflow', '-log', str(directory / 'nextflow.log'), 'run', str(fixture),
        '-work-dir', str(directory / 'work'),
        '--fixture_root', str(directory / 'inputs'),
        '--fixture_output', str(directory / 'observed.json'),
    ], cwd=directory, env=env, check=True)
    observations = json.loads((directory / 'observed.json').read_text())['cases']
    assert len(observations) == 24
    assert all(row['observed'] == row['expected'] for row in observations)
    report = {'status': 'pass', 'arrival_permutations': 24, 'interleaved_samples': 24,
              'main_nf_sha256': hashlib.sha256(source.encode()).hexdigest()}
    (directory / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
