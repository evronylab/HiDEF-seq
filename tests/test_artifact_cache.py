"""Prepared-artifact identity, publication, reuse and integrity contracts.

Run with Python and Rscript from the pipeline environment. Every test uses a
fresh temporary cache and leaves existing workflow caches untouched.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import random
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

HELPER = Path(__file__).resolve().parents[1] / 'bin/artifactCache.R'
SHARED = HELPER.with_name('sharedFunctions.R')
RSCRIPT = shutil.which(os.environ.get('RSCRIPT', 'Rscript'))


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False)


@unittest.skipUnless(RSCRIPT, 'Rscript is required (use the pipeline environment)')
class ArtifactCacheRTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.root = self.base / 'cache'
        self.source = self.base / 'source'
        self.source.mkdir()
        (self.source / 'result.qs2').write_bytes(b'closed product\x00')
        self.identity = {'schema': 1, 'namespace': 'unit', 'settings': {'threshold': 2}}
        self.identity_file = self.base / 'identity.json'
        self.identity_file.write_text(json.dumps(self.identity))

    def cli(self, operation, *extra, cwd=None, check=True):
        command = [RSCRIPT, str(HELPER)]
        command += [operation, '--root', str(self.root), '--identity', str(self.identity_file), *map(str, extra)]
        return subprocess.run(command, cwd=cwd, check=check, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def r(self, expression, check=True):
        return subprocess.run([RSCRIPT, '-e', 'source(' + json.dumps(str(SHARED)) + ');' + expression],
                              text=True, check=check, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def destination(self):
        return self.root / 'unit' / hashlib.sha256(canonical(self.identity).encode()).hexdigest()

    def publish(self):
        return Path(self.cli('publish', '--source', self.source).stdout.strip())

    def test_plain_and_exact_serialized_keys_match_python_bytes(self):
        identities = [self.identity]
        for number in [1.0, -0.0, 1e-7, 1e-4, 1e16, 1.2345678901234567, 9007199254740993, 10**70]:
            identities.append(dict(schema=1, namespace='unit', settings=dict(number=number, label='é/😀\x1f\x7f', empty={}, array=[])))
        for raw in ('{"namespace":"unit","schema":1,"threshold":1E-7,"label":"é"}',
                    '{"namespace":"unit","schema":1,"threshold":-0.0,"label":"\\u00e9"}'):
            identities.append(dict(schema=1, namespace='unit', serialized_identity=raw))
        payload = self.base / 'identities.json'
        payload.write_text(json.dumps(identities))
        output = self.r('x<-artifact_cache_read_json(' + json.dumps(str(payload)) + ');cat(vapply(x,artifact_cache_key,character(1)),sep="\\n")')
        expected = [hashlib.sha256(item.get('serialized_identity', canonical(item)).encode()).hexdigest() for item in identities]
        self.assertEqual(output.stdout.splitlines(), expected)

    def test_float_rendering_matches_python_for_seeded_ieee_values(self):
        rng = random.Random(9132)
        values = [0.0, -0.0, 1e-4, 1e-5, 1e15, 1e16, 1.0, 5e-324, 1.7976931348623157e308]
        while len(values) < 10009:
            value = struct.unpack('>d', rng.getrandbits(64).to_bytes(8, 'big'))[0]
            if math.isfinite(value): values.append(value)
        path = self.base / 'numbers.json'
        path.write_text(json.dumps(values))
        output = self.r('x<-artifact_cache_read_json(' + json.dumps(str(path)) + ');cat(vapply(x,artifact_cache_canonical,character(1)),sep="\\n")')
        actual = output.stdout.splitlines()
        expected = [canonical(value) for value in values]
        differences = [(i, values[i], a, b) for i, (a,b) in enumerate(zip(actual,expected)) if a != b]
        self.assertEqual(differences[:10], [])
        self.assertEqual(len(actual), len(expected))

    def test_identify_content_settings_and_complete_manifest_bytes(self):
        input_file = self.base / 'input'; input_file.write_text('reference')
        spec = dict(namespace='reference', inputs=dict(fasta=str(input_file)), settings=dict(circular=['MT']), tools={})
        path = self.base / 'spec.json'; path.write_text(json.dumps(spec))
        output = self.base / 'identified.json'
        def identify():
            return subprocess.check_output([RSCRIPT, str(HELPER), 'identify', '--spec', str(path), '--output', str(output)], text=True).strip()
        original = identify()
        expected = dict(schema=1, namespace='reference', inputs=dict(fasta=dict(bytes=9, sha256=hashlib.sha256(b'reference').hexdigest())), scripts={}, settings=spec['settings'], tools={})
        self.assertEqual(output.read_text(), canonical(expected) + '\n')
        self.assertEqual(original, hashlib.sha256(canonical(expected).encode()).hexdigest())
        input_file.write_text('Reference'); self.assertNotEqual(original, identify())
        destination = self.publish()
        manifest = json.loads((destination / 'manifest.complete.json').read_text())
        self.assertEqual((destination / 'manifest.complete.json').read_text(), canonical(manifest) + '\n')
        self.assertEqual(destination, self.destination())
        self.assertEqual(self.cli('verify').stdout.strip(), str(destination))

    def test_verification_uses_numeric_value_equality_without_rounding_large_integers(self):
        self.identity = dict(schema=1, namespace='unit', serialized_identity='{"schema":1,"namespace":"unit"}',
                             threshold=1, zero=-0.0, big=9007199254740993)
        self.identity_file.write_text(json.dumps(self.identity))
        destination = self.publish()
        manifest_path = destination/'manifest.complete.json'
        manifest = json.loads(manifest_path.read_text())
        manifest['products']['result.qs2']['bytes'] = float(manifest['products']['result.qs2']['bytes'])
        manifest_path.write_text(json.dumps(manifest))
        self.identity.update(schema=1.0, threshold=1.0, zero=0)
        self.identity_file.write_text(json.dumps(self.identity))
        self.cli('verify')
        self.identity['big'] = float(self.identity['big'])
        self.identity_file.write_text(json.dumps(self.identity))
        self.assertIn('identity mismatch', self.cli('verify', check=False).stderr)

    def test_corruption_never_overwrites(self):
        destination = self.publish(); (destination / 'result.qs2').write_bytes(b'corrupted')
        for operation, extra in [('verify', []), ('publish', ['--source', self.source])]:
            result = self.cli(operation, *extra, check=False)
            self.assertNotEqual(result.returncode, 0); self.assertIn('integrity', result.stderr)
        self.assertEqual((destination / 'result.qs2').read_bytes(), b'corrupted')

    def test_manifest_mismatch_and_incomplete_state(self):
        destination = self.publish(); path = destination / 'manifest.complete.json'
        original = json.loads(path.read_text())
        for field, value in [('state','building'),('schema',2),('key','wrong'),('identity',{})]:
            manifest = dict(original); manifest[field] = value; path.write_text(json.dumps(manifest))
            self.assertNotEqual(self.cli('verify', check=False).returncode, 0)
        path.unlink(); self.assertNotEqual(self.cli('verify', check=False).returncode, 0)

    def test_symlink_fifo_empty_and_reserved_products_rejected(self):
        (self.source / 'alias').symlink_to('result.qs2')
        self.assertIn('symlinks', self.cli('publish', '--source', self.source, check=False).stderr)
        (self.source / 'alias').unlink(); os.mkfifo(self.source / 'fifo')
        self.assertIn('regular files', self.cli('publish', '--source', self.source, check=False).stderr)
        (self.source / 'fifo').unlink(); (self.source / 'manifest.complete.json').write_text('{}')
        self.assertIn('reserved', self.cli('publish', '--source', self.source, check=False).stderr)
        (self.source / 'manifest.complete.json').unlink(); (self.source / 'result.qs2').unlink()
        self.assertIn('empty', self.cli('publish', '--source', self.source, check=False).stderr)
        self.assertFalse(self.destination().exists())

    def test_cold_creation_and_warm_reuse(self):
        counter = self.base / 'builds.txt'
        builder = ('from pathlib import Path;'
                   'Path("library/sub").mkdir(parents=True);'
                   'Path("library/sub/value").write_text("library");'
                   'Path("result.qs2").write_bytes(b"built product");'
                   'open(' + repr(str(counter)) + ',"a").write("built\\n")')
        command = ('--product', 'result.qs2', '--product', 'library', '--',
                   sys.executable, '-c', builder)
        cold = self.base / 'cold'
        warm = self.base / 'warm'
        cold.mkdir()
        warm.mkdir()
        self.cli('run', *command, cwd=cold)
        destination = self.destination()
        self.assertTrue(destination.is_dir())
        manifest = destination / 'manifest.complete.json'
        original_manifest = (manifest.read_bytes(), manifest.stat().st_ino, manifest.stat().st_mtime_ns)
        self.cli('run', *command, cwd=warm)
        self.assertEqual(counter.read_text(), 'built\n')
        self.assertEqual((manifest.read_bytes(), manifest.stat().st_ino, manifest.stat().st_mtime_ns),
                         original_manifest)
        self.assertEqual(self.cli('verify').stdout.strip(), str(destination))
        for name in ('result.qs2', 'library/sub/value'):
            self.assertEqual((warm / name).read_bytes(), (cold / name).read_bytes())
            self.assertEqual((warm / name).stat().st_ino, (destination / name).stat().st_ino)

    def test_existing_entry_is_not_overwritten(self):
        destination = self.publish()
        manifest = destination / 'manifest.complete.json'
        original = (manifest.read_bytes(), manifest.stat().st_ino, manifest.stat().st_mtime_ns)
        (self.source / 'result.qs2').write_bytes(b'another product')
        self.assertEqual(self.publish(), destination)
        self.assertEqual((destination / 'result.qs2').read_bytes(), b'closed product\x00')
        self.assertEqual((manifest.read_bytes(), manifest.stat().st_ino, manifest.stat().st_mtime_ns), original)
        self.cli('verify')

    def test_version_directory_is_not_reused_or_migrated(self):
        destination = self.publish()
        versioned = self.root / 'v1' / 'unit' / destination.name
        versioned.parent.mkdir(parents=True)
        destination.rename(versioned)
        old_manifest = (versioned / 'manifest.complete.json').read_bytes()
        work = self.base / 'task'
        work.mkdir()
        self.cli('run', '--product', 'result.qs2', '--', sys.executable, '-c',
                 'from pathlib import Path;Path("result.qs2").write_bytes(b"new path product")', cwd=work)
        self.assertEqual((destination / 'result.qs2').read_bytes(), b'new path product')
        self.assertEqual((versioned / 'result.qs2').read_bytes(), b'closed product\x00')
        self.assertEqual((versioned / 'manifest.complete.json').read_bytes(), old_manifest)
        self.assertEqual(self.cli('verify').stdout.strip(), str(destination))


if __name__ == '__main__':
    unittest.main()
