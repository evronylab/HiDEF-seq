import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

MODULE = Path(__file__).resolve().parents[1] / "scripts/benchmark/compare_outputs.py"
spec = importlib.util.spec_from_file_location("compare_outputs", MODULE)
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)


class PublishedContentsTest(unittest.TestCase):
    def test_measure_fallback_waited_child_and_failure(self):
        wrapper = MODULE.parent / "measure.py"
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "metrics"
            result = subprocess.run([
                sys.executable, str(wrapper), "--out", str(prefix), "--label", "child-test",
                "--time", "/nonexistent/gnu-time", "--allocated-cpus", "2", "--",
                sys.executable, "-c",
                "import subprocess,sys; subprocess.run([sys.executable,'-c','sum(range(100000))']); sys.exit(3)",
            ])
            metrics = json.loads(Path(str(prefix) + ".json").read_text())
            self.assertEqual(result.returncode, 3)
            self.assertEqual(metrics["exit_status"], 3)
            self.assertGreater(metrics["actual_cpu_seconds"], 0)
            self.assertEqual(metrics["actual_cpu_seconds"], metrics["user_seconds"] + metrics["system_seconds"])
            self.assertEqual(metrics["allocated_cpu_hours"], 2 * metrics["wall_seconds"] / 3600)
            self.assertGreater(metrics["max_rss_kib"], 0)

    def test_vcf_encoding_metadata_and_scientific_change(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / "a.vcf.gz", Path(directory) / "b.vcf.gz"
            body = "##fileDate=20261004\n#CHROM\tPOS\nchr1\t1\n"
            with gzip.open(a, "wt") as handle:
                handle.write(body)
            with gzip.open(b, "wt") as handle:
                handle.write(body.replace("20261004", "20261005"))
            self.assertEqual(compare.compare_text(a, b)["status"], "pass")
            with gzip.open(b, "wt") as handle:
                handle.write(body.replace("chr1\t1", "chr1\t2"))
            self.assertEqual(compare.compare_text(a, b)["status"], "failure")

    def test_inventory_binary_and_missing_are_not_silently_passed(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / "a", Path(directory) / "b"
            a.mkdir(); b.mkdir()
            (a / "sample.qs2").write_bytes(b"a")
            (b / "sample.qs2").write_bytes(b"b")
            (a / "missing.tsv").write_text("x\n1\n")
            (a / "sample.run_metadata.tsv").write_text("old")
            (b / "sample.run_metadata.tsv").write_text("new")
            records = {r["path"]: r for r in compare.compare_trees(a, b, ["*.run_metadata.tsv"])}
            self.assertEqual(records["sample.qs2"]["status"], "review")
            self.assertEqual(records["missing.tsv"]["status"], "failure")
            self.assertEqual(records["sample.run_metadata.tsv"]["status"], "ignored")

    def test_bgz_vcf_inventory_and_concatenated_gzip_members(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / "a", Path(directory) / "b"
            a.mkdir(); b.mkdir()
            for suffix in (".bgz", ".bgzf"):
                name = "calls.vcf" + suffix
                (a / name).write_bytes(gzip.compress(b"##fileDate=20261004\n#CHROM\tPOS\n") +
                                      gzip.compress(b"chr1\t1\nchr2\t2\n"))
                (b / name).write_bytes(gzip.compress(b"##fileDate=20261005\n#CHROM\tPOS\nchr1\t1\nchr2\t2\n"))
                records = {r["path"]: r for r in compare.compare_trees(a, b)}
                self.assertEqual(records[name]["status"], "pass")
                (b / name).write_bytes(gzip.compress(b"##fileDate=20261005\n#CHROM\tPOS\nchr2\t2\nchr1\t1\n"))
                records = {r["path"]: r for r in compare.compare_trees(a, b)}
                self.assertEqual(records[name]["status"], "failure")

    def test_published_precision_and_row_order_are_exact(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / "a.tsv", Path(directory) / "b.tsv"
            a.write_text("value\n1.0000000000000\n2\n")
            b.write_text("value\n1.0000000000001\n2\n")
            self.assertEqual(compare.compare_text(a, b)["status"], "failure")
            b.write_text("value\n2\n1.0000000000000\n")
            self.assertEqual(compare.compare_text(a, b)["status"], "failure")

    def test_byte_fast_path_preserves_encoding_and_newline_exceptions(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / "a.bed.gz", Path(directory) / "b.bed.gz"
            body = "chr1\t0\t1\t1\t€\nchr1\t1\t2\t2\tACG\n".encode()
            a.write_bytes(gzip.compress(body))
            b.write_bytes(gzip.compress(body[:13]) + gzip.compress(body[13:]))
            self.assertTrue(compare.decompressed_bytes_equal(a, b, block_size=1))
            self.assertEqual(compare.compare_text(a, b)["reason"], "decompressed text bytes identical")
            b.write_bytes(gzip.compress(body.replace(b"\n", b"\r\n").rstrip(b"\r\n")))
            self.assertFalse(compare.decompressed_bytes_equal(a, b, block_size=7))
            self.assertEqual(compare.compare_text(a, b)["status"], "pass")

    def test_byte_fast_path_checks_utf8_crc_and_truncation(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / "a.tsv.gz", Path(directory) / "b.tsv.gz"
            for payload, error in ((gzip.compress(b"invalid\xff\n"), UnicodeDecodeError),
                                   (gzip.compress(b"valid\n")[:-4], EOFError)):
                a.write_bytes(payload); b.write_bytes(payload)
                with self.assertRaises(error):
                    compare.compare_text(a, b)
            payload = bytearray(gzip.compress(b"valid\n"))
            payload[-8] ^= 1
            a.write_bytes(payload); b.write_bytes(payload)
            with self.assertRaises(gzip.BadGzipFile):
                compare.compare_text(a, b)


if __name__ == "__main__":
    unittest.main()
