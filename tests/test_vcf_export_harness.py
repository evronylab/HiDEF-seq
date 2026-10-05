"""Protect published inputs while validating an index against its own encoding."""
import gzip
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

MODULE = Path(__file__).resolve().parents[1] / "scripts/benchmark/profile_germline_vcf.py"
spec = importlib.util.spec_from_file_location("vcf_profile", MODULE)
profile = importlib.util.module_from_spec(spec)
spec.loader.exec_module(profile)


class IndexValidationTests(unittest.TestCase):
    def check_rebuild(self, rebuilt_content):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "published.vcf.bgz"
            source.write_bytes(b"unchanged BGZF encoding")
            original_index = Path(str(source) + ".tbi")
            with gzip.open(original_index, "wb") as handle:
                handle.write(b"correct virtual offsets")
            original_bytes = original_index.read_bytes()
            def fingerprint(path):
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                stat = path.stat()
                # atime may change on read; all mutation-relevant metadata must not.
                return digest, stat.st_dev, stat.st_ino, stat.st_mode, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns
            before = {path: fingerprint(path) for path in (source, original_index)}

            def rebuild(command, **kwargs):
                staged = Path(command[-1])
                self.assertFalse(staged.is_symlink())
                self.assertNotEqual(staged.resolve(), source.resolve())
                self.assertEqual(staged.read_bytes(), source.read_bytes())
                with gzip.open(str(staged) + ".tbi", "wb") as handle:
                    handle.write(rebuilt_content)

            try:
                with patch.object(profile.subprocess, "run", side_effect=rebuild):
                    profile.validate_own_index(source)
            finally:
                self.assertEqual(original_index.read_bytes(), original_bytes)
                self.assertEqual(source.read_bytes(), b"unchanged BGZF encoding")
                self.assertEqual({path: fingerprint(path) for path in before}, before)

    def test_copy_keeps_published_data_and_index_immutable(self):
        self.check_rebuild(b"correct virtual offsets")

    def test_wrong_offsets_fail_without_touching_original(self):
        with self.assertRaisesRegex(AssertionError, "own BGZF encoding"):
            self.check_rebuild(b"wrong virtual offsets")


if __name__ == "__main__":
    unittest.main()
