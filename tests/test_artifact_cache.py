"""R cache protocol regression tests; run inside the pinned HiDEF container.

ARTIFACT_CACHE_R_HELPER may select a frozen candidate; ARTIFACT_CACHE_PYTHON
optionally adds interoperability checks against a frozen pre-port Python helper.
"""
import concurrent.futures
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import random
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import time
import unittest

HELPER = Path(os.environ.get('ARTIFACT_CACHE_R_HELPER', Path(__file__).resolve().parents[1] / 'bin/artifactCache.R'))
SHARED = HELPER.with_name('sharedFunctions.R')
RSCRIPT = shutil.which(os.environ.get('RSCRIPT', 'Rscript'))
PYTHON_HELPER = os.environ.get('ARTIFACT_CACHE_PYTHON')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False)


@unittest.skipUnless(RSCRIPT, 'Rscript is required (use the pinned container)')
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

    def cli(self, operation, *extra, cwd=None, check=True, helper=None):
        command = ([sys.executable, helper] if helper else [RSCRIPT, str(HELPER)])
        command += [operation, '--root', str(self.root), '--identity', str(self.identity_file), *map(str, extra)]
        return subprocess.run(command, cwd=cwd, check=check, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def r(self, expression, check=True):
        return subprocess.run([RSCRIPT, '-e', 'source(' + json.dumps(str(SHARED)) + ');' + expression],
                              text=True, check=check, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def destination(self):
        return self.root / 'v1/unit' / hashlib.sha256(canonical(self.identity).encode()).hexdigest()

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
        if PYTHON_HELPER: self.cli('verify', helper=PYTHON_HELPER)
        self.identity['big'] = float(self.identity['big'])
        self.identity_file.write_text(json.dumps(self.identity))
        self.assertIn('identity mismatch', self.cli('verify', check=False).stderr)
        if PYTHON_HELPER: self.assertNotEqual(self.cli('verify', check=False, helper=PYTHON_HELPER).returncode, 0)

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

    def test_failed_copy_sync_rename_and_changed_source_never_complete(self):
        for failure in ['artifact_cache_copy<-function(...)stop("copy failed")',
                        'artifact_cache_fsync<-function(...)stop("sync failed")',
                        'file.rename<-function(...)FALSE',
                        'original<-artifact_cache_copy;artifact_cache_copy<-function(source,target){writeLines("changed",source);original(source,target)}']:
            expression = ';'.join([failure, 'id<-artifact_cache_read_json(' + json.dumps(str(self.identity_file)) + ')',
                'artifact_cache_publish(' + ','.join(map(json.dumps,[str(self.root)])) + ',id,' + json.dumps(str(self.source)) + ')'])
            self.assertNotEqual(self.r(expression, check=False).returncode, 0)
            self.assertFalse(self.destination().exists())
            self.assertFalse(list(self.root.rglob('manifest.complete.json')))
            self.assertFalse([x for x in self.root.rglob('*') if x.is_dir() and x.name.startswith('.')])

    def test_file_manifest_directory_and_parent_sync_barriers(self):
        for fail_at in (2,3,4):
            self.root=self.base/('sync-'+str(fail_at))
            expression=';'.join([
                'original<-artifact_cache_fsync;calls<-0L',
                'artifact_cache_fsync<-function(paths){calls<<-calls+1L;if(calls=='+str(fail_at)+')stop("injected sync failure");original(paths)}',
                'id<-artifact_cache_read_json('+json.dumps(str(self.identity_file))+')',
                'artifact_cache_publish('+json.dumps(str(self.root))+',id,'+json.dumps(str(self.source))+')'])
            result=self.r(expression,check=False)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('injected sync failure',result.stderr)
            # Parent-directory fsync is after atomic rename, as in the Python
            # implementation; a failure there still leaves a verifiable entry.
            self.assertEqual(self.destination().exists(), fail_at==4)
            if fail_at==4: self.cli('verify')
            else: self.assertFalse(list(self.root.rglob('manifest.complete.json')))

    def test_identify_rejects_inputs_changed_during_hashing(self):
        source=self.base/'input'; source.write_text('before')
        spec=self.base/'spec.json'; spec.write_text(json.dumps(dict(namespace='unit',inputs=dict(input=str(source)))))
        expression=';'.join([
            'original<-artifact_cache_sha256',
            'artifact_cache_sha256<-function(path){digest<-original(path);writeLines("after",path);digest}',
            'artifact_cache_identify(artifact_cache_read_json('+json.dumps(str(spec))+'))'])
        result=self.r(expression,check=False)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('changed while hashing',result.stderr)

    def test_warm_nested_restore_hardlinks_and_does_not_build(self):
        (self.source / 'library/sub').mkdir(parents=True)
        (self.source / 'library/sub/hidden').write_text('library')
        destination = self.publish(); work = self.base / 'task'; work.mkdir()
        self.cli('run','--product','result.qs2','--product','library','--',sys.executable,'-c','raise SystemExit(4)',cwd=work)
        for name in ['result.qs2','library/sub/hidden']:
            self.assertEqual((work/name).read_bytes(), (destination/name).read_bytes())
            self.assertEqual((work/name).stat().st_ino, (destination/name).stat().st_ino)
        self.assertIn('already exists', self.cli('run','--product','result.qs2','--','true',cwd=work,check=False).stderr)

    def test_cross_filesystem_copy_and_existing_target_failure(self):
        if not Path('/dev/shm').is_dir(): self.skipTest('/dev/shm unavailable')
        with tempfile.TemporaryDirectory(dir='/dev/shm') as temporary:
            target = Path(temporary)/'result.qs2'
            expression = 'artifact_cache_link_or_copy(' + json.dumps(str(self.source/'result.qs2')) + ',' + json.dumps(str(target)) + ')'
            self.r(expression)
            self.assertEqual(target.read_bytes(), (self.source/'result.qs2').read_bytes())
            self.assertNotEqual(self.r(expression, check=False).returncode, 0)

    def test_failed_missing_and_invalid_products_never_publish(self):
        work = self.base/'work'; work.mkdir()
        for product, command in [('result.qs2',[sys.executable,'-c','raise SystemExit(3)']),('absent',['true']),('../outside',['true']),('/absolute',['true']),('manifest.complete.json',['true'])]:
            result = self.cli('run','--product',product,'--',*command,cwd=work,check=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(self.destination().exists())

    def test_concurrent_publishers_and_cold_builders(self):
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            outputs = list(pool.map(lambda _: self.publish(), range(2)))
        self.assertEqual(outputs[0], outputs[1]); self.cli('verify')
        self.root = self.base/'second-cache'; counter=self.base/'count'
        works=[self.base/'work1',self.base/'work2']
        for work in works: work.mkdir()
        builder='from pathlib import Path;import time;time.sleep(0.4);Path("result.qs2").write_bytes(b"built");open('+repr(str(counter))+',"a").write("build\\n")'
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            results=list(pool.map(lambda work:self.cli('run','--product','result.qs2','--',sys.executable,'-c',builder,cwd=work),works))
        self.assertEqual(counter.read_text(),'build\n')
        self.assertEqual(results[0].stdout,results[1].stdout)

    def test_flock_interoperability_exception_and_process_death_release(self):
        lock=self.base/'lock'; marker=self.base/'acquired'
        command=[RSCRIPT,'-e','source('+json.dumps(str(SHARED))+');artifact_cache_with_lock('+json.dumps(str(lock))+',function(){file.create('+json.dumps(str(marker))+');Sys.sleep(30)})']
        with lock.open('a') as handle:
            fcntl.flock(handle,fcntl.LOCK_EX)
            process=subprocess.Popen(command,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True)
            try:
                time.sleep(.7); self.assertFalse(marker.exists()); self.assertIsNone(process.poll())
                fcntl.flock(handle,fcntl.LOCK_UN)
                deadline=time.monotonic()+10
                while not marker.exists() and time.monotonic()<deadline: time.sleep(.03)
                self.assertTrue(marker.exists())
                with self.assertRaises(BlockingIOError): fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
                process.kill(); process.wait(timeout=10)
                fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
            finally:
                if process.poll() is None: process.kill(); process.wait()
                process.stderr.close()
        result=self.r('try(artifact_cache_with_lock('+json.dumps(str(lock))+',function()stop("expected")));artifact_cache_with_lock('+json.dumps(str(lock))+',function()cat("released"))')
        self.assertIn('released',result.stdout)

    @unittest.skipUnless(PYTHON_HELPER, 'optional frozen Python baseline not specified')
    def test_python_r_publish_verify_and_mixed_concurrent_builders(self):
        self.cli('publish','--source',self.source,helper=PYTHON_HELPER); self.cli('verify')
        self.root=self.base/'r-cache'; self.publish(); self.cli('verify',helper=PYTHON_HELPER)
        self.root=self.base/'mixed-cache'; counter=self.base/'count'; works=[self.base/'one',self.base/'two']
        for work in works: work.mkdir()
        builder='from pathlib import Path;import time;time.sleep(0.4);Path("result.qs2").write_bytes(b"built");open('+repr(str(counter))+',"a").write("build\\n")'
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(self.cli,'run','--product','result.qs2','--',sys.executable,'-c',builder,cwd=work,helper=helper) for work,helper in zip(works,[None,PYTHON_HELPER])]
            for future in futures: future.result()
        self.assertEqual(counter.read_text(),'build\n')
        self.cli('verify'); self.cli('verify',helper=PYTHON_HELPER)


if __name__=='__main__': unittest.main()
