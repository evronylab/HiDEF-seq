#!/usr/bin/env python3
"""Exercise actual preparedInputIdentity() with the pinned Nextflow on compute.

Run: python3 tests/test_prepared_input_identity.py --directory NEW_FIXTURE_DIR
Optional --probe PATH records read-only metadata for cross-node comparisons.
"""
import argparse
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--probe', type=Path, action='append', default=[])
    args = parser.parse_args()
    directory = args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    source = (args.repo / 'main.nf').read_text()
    helper = source[source.index('def preparedInputIdentity('):source.index('def canonicalCacheValue(')]
    (directory / 'production-helper.nf.txt').write_text(helper)
    (directory / 'probes.json').write_text(json.dumps([str(p.resolve(strict=True)) for p in args.probe]))
    fixture = helper + r'''
workflow {
  def checks = []
  def input = file('input.txt')
  input.text = 'AAAA'
  def first = preparedInputIdentity(input)
  assert first == preparedInputIdentity(input)
  assert first.keySet() == ['kind','path','bytes','modified','changed','inode'].toSet()
  assert first.kind == 'canonical-path-size-mtime-ctime-inode'
  assert first.inode == java.nio.file.Files.getAttribute(input, 'unix:ino').toString()
  checks << 'repeat-stable-and-inode-without-node-device'
  def alias = file('alias.txt')
  java.nio.file.Files.createSymbolicLink(alias, input)
  assert preparedInputIdentity(alias) == first
  checks << 'canonical-symlink-identity'
  def modified = java.nio.file.Files.getLastModifiedTime(input)
  java.nio.file.Files.setLastModifiedTime(input,
    java.nio.file.attribute.FileTime.from(modified.toMillis()+1000, java.util.concurrent.TimeUnit.MILLISECONDS))
  assert preparedInputIdentity(input) != first
  checks << 'mtime-change-invalidates'
  def replacement = file('replacement.txt')
  replacement.text = 'BBBB'
  java.nio.file.Files.setLastModifiedTime(replacement, modified)
  java.nio.file.Files.move(replacement, input, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
  def replaced = preparedInputIdentity(input)
  assert replaced.path == first.path && replaced.bytes == first.bytes && replaced.modified == first.modified
  assert replaced.inode != first.inode && replaced != first
  checks << 'same-path-size-mtime-replacement-invalidates'
  input.text = 'longer content'
  assert preparedInputIdentity(input).bytes != replaced.bytes
  checks << 'size-change-invalidates'
  def probes = new groovy.json.JsonSlurper().parse(file('probes.json').toFile()).collect { name ->
    [identity: preparedInputIdentity(name), localFileKey: java.nio.file.Files.readAttributes(
      file(name), java.nio.file.attribute.BasicFileAttributes).fileKey().toString()]
  }
  file('results.json').text = groovy.json.JsonOutput.prettyPrint(
    groovy.json.JsonOutput.toJson([checks: checks, probes: probes]))
}
'''
    (directory / 'main.nf').write_text(fixture)
    subprocess.run(['nextflow', '-log', str(directory / 'nextflow.log'), 'run', str(directory / 'main.nf'),
                    '-work-dir', str(directory / 'work')], cwd=directory, check=True)
    report = json.loads((directory / 'results.json').read_text())
    assert len(report['checks']) == 5
    for row in report['probes']:
        identity = row['identity']
        stat = Path(identity['path']).stat()
        # Java Unix inode is a signed long; Python exposes the unsigned number.
        signed = stat.st_ino if stat.st_ino < 2**63 else stat.st_ino - 2**64
        assert identity['inode'] == str(signed)
        assert identity['bytes'] == stat.st_size
        assert 'fileKey' not in identity and 'device' not in identity
    print(json.dumps({'status': 'pass', 'checks': len(report['checks']), 'probes': len(report['probes'])}))


if __name__ == '__main__':
    main()
