"""Tiny local tests for cache identity, completion and concurrent publication."""
import concurrent.futures
import errno
import importlib.util
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


HELPER = Path(__file__).resolve().parents[1] / "bin" / "artifactCache.py"
spec = importlib.util.spec_from_file_location("artifactCache", HELPER)
cache = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cache)


class ArtifactCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(dir=Path(__file__).parent)
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.root = self.base / "cache"
        self.source = self.base / "source"
        self.source.mkdir()
        (self.source / "result.qs2").write_bytes(b"closed product\x00")
        self.identity = {"schema": 1, "namespace": "unit", "settings": {"threshold": 2}}

    def test_content_and_settings_change_identity(self):
        source = self.base / "input"
        source.write_text("reference")
        spec = {"namespace": "reference", "inputs": {"fasta": str(source)},
                "settings": {"circular": ["MT"]}, "tools": {"container": "digest"}}
        original = cache.identify(spec)
        source.write_text("Reference")
        self.assertNotEqual(cache.key(original), cache.key(cache.identify(spec)))
        source.write_text("reference")
        spec["settings"]["circular"] = []
        self.assertNotEqual(cache.key(original), cache.key(cache.identify(spec)))

    def test_exact_serialized_identity_cross_language_contract(self):
        for raw in ('{"namespace":"unit","schema":1,"threshold":1E-7,"label":"é"}',
                    '{"namespace":"unit","schema":1,"threshold":-0.0,"label":"\\u00e9"}'):
            envelope = {"namespace": "unit", "schema": 1, "serialized_identity": raw}
            self.assertEqual(cache.key(envelope), hashlib.sha256(raw.encode("utf-8")).hexdigest())

    def test_manifest_and_corruption_detection(self):
        destination = cache.publish(self.root, self.identity, self.source)
        self.assertEqual(destination, cache.verify(self.root, self.identity))
        (destination / "result.qs2").write_bytes(b"changed product")
        with self.assertRaisesRegex(ValueError, "integrity"):
            cache.verify(self.root, self.identity)
        with self.assertRaisesRegex(ValueError, "integrity"):
            cache.publish(self.root, self.identity, self.source)

    def test_failed_copy_never_publishes_completion(self):
        with patch.object(cache.shutil, "copy2", side_effect=OSError("copy failed")):
            with self.assertRaises(OSError):
                cache.publish(self.root, self.identity, self.source)
        self.assertFalse(cache.location(self.root, self.identity).exists())
        self.assertFalse(list(self.root.rglob(cache.MANIFEST)))

    def test_concurrent_publishers_reuse_one_complete_bundle(self):
        identity_path = self.base / "identity.json"
        identity_path.write_text(json.dumps(self.identity))
        command = [sys.executable, str(HELPER), "publish", "--root", str(self.root),
                   "--identity", str(identity_path), "--source", str(self.source)]
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            outputs = list(pool.map(lambda _: subprocess.check_output(command), range(2)))
        self.assertEqual(outputs[0], outputs[1])
        cache.verify(self.root, self.identity)

    def test_run_hit_restores_files_without_running_builder(self):
        cache.publish(self.root, self.identity, self.source)
        work = self.base / "task"
        work.mkdir()
        identity_path = self.base / "identity.json"
        identity_path.write_text(json.dumps(self.identity))
        subprocess.run([sys.executable, str(HELPER), "run", "--root", str(self.root),
                        "--identity", str(identity_path), "--product", "result.qs2",
                        "--", sys.executable, "-c", "raise RuntimeError('must not run')"],
                       cwd=work, check=True, stdout=subprocess.PIPE)
        self.assertEqual((work / "result.qs2").read_bytes(), (self.source / "result.qs2").read_bytes())

    def test_concurrent_cold_runs_execute_builder_once(self):
        identity_path = self.base / "identity.json"
        identity_path.write_text(json.dumps(self.identity))
        counter = self.base / "build-count"
        work_dirs = [self.base / "task1", self.base / "task2"]
        for directory in work_dirs:
            directory.mkdir()
        builder = ("from pathlib import Path; "
                   "Path('result.qs2').write_bytes(b'complete'); "
                   "open(" + repr(str(counter)) + ", 'a').write('build\\n')")
        command = [sys.executable, str(HELPER), "run", "--root", str(self.root),
                   "--identity", str(identity_path), "--product", "result.qs2",
                   "--", sys.executable, "-c", builder]
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            results = list(pool.map(lambda directory: subprocess.run(command, cwd=directory,
                                check=True, stdout=subprocess.PIPE), work_dirs))
        self.assertEqual(counter.read_text(), "build\n")
        for directory in work_dirs:
            self.assertEqual((directory / "result.qs2").read_bytes(), b"complete")

    def test_cross_filesystem_restore_copies(self):
        target = self.base / "restored.qs2"
        with patch.object(cache.os, "link", side_effect=OSError(errno.EXDEV, "different filesystem")):
            cache.link_or_copy(self.source / "result.qs2", target)
        self.assertEqual(target.read_bytes(), (self.source / "result.qs2").read_bytes())

    def test_failed_builder_does_not_publish(self):
        identity_path = self.base / "identity.json"
        identity_path.write_text(json.dumps(self.identity))
        result = subprocess.run([sys.executable, str(HELPER), "run", "--root", str(self.root),
                                 "--identity", str(identity_path), "--product", "result.qs2",
                                 "--", sys.executable, "-c", "raise SystemExit(3)"],
                                cwd=self.source, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(cache.location(self.root, self.identity).exists())

    def test_symlink_products_are_rejected(self):
        (self.source / "alias").symlink_to("result.qs2")
        with self.assertRaisesRegex(ValueError, "symlinks"):
            cache.publish(self.root, self.identity, self.source)


if __name__ == "__main__":
    unittest.main()
