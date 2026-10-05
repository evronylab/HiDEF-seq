"""Driver tests use no R and no biological data."""
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/benchmark/benchmark_germline_formatting_fresh.py"
spec = importlib.util.spec_from_file_location("fresh_format", SCRIPT)
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)


class FreshFormattingTests(unittest.TestCase):
    def test_workers_share_packet_alternate_order_and_record_both_peaks(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for mode in ("baseline", "candidate"):
                folder = root / mode / "bin"
                folder.mkdir(parents=True)
                for name in ("outputResults.R", "sharedFunctions.R"):
                    (folder / name).write_text("# frozen fixture\n")
            (root / "input.qs2").write_bytes(b"fixture")
            commands = []

            def fake_run(command, **kwargs):
                commands.append(command)
                if "--out" in command:
                    prefix = Path(command[command.index("--out") + 1])
                    Path(str(prefix) + ".json").write_text(json.dumps(dict(actual_cpu_seconds=3, max_rss_kib=100)))
                    arm = Path(str(prefix).removesuffix(".process"))
                    Path(str(arm) + ".phase.tsv").write_text("actual_cpu_seconds\tpeak_rss_kib\n2\t80\n")
                    Path(str(arm) + ".process_peak.tsv").write_text(
                        "preformat_peak_rss_kib\tremaining_peak_rss_kib\n150\t120\n")
                return subprocess.CompletedProcess(command, 0)

            with patch.object(driver.subprocess, "run", side_effect=fake_run):
                driver.main(["--baseline-repo", str(root / "baseline"), "--candidate-repo", str(root / "candidate"),
                             "--input", str(root / "input.qs2"), "--output", str(root / "output"),
                             "--pairs", "2", "--allocated-cpus", "2"])
            workers = [c for c in commands if "format" in c]
            self.assertEqual([c[-1] for c in workers], ["baseline", "candidate", "candidate", "baseline"])
            self.assertEqual(len({c[c.index("format") + 1] for c in workers}), 1)
            self.assertIn("prepare", commands[0])
            self.assertTrue(all("compare" not in c for c in commands[:5]))
            results = json.loads((root / "output/results.json").read_text())
            self.assertEqual(len(results["comparisons"]), 2)
            self.assertTrue(all(m["process"]["max_rss_kib_across_phase_reset"] == 150 for m in results["measurements"]))


if __name__ == "__main__":
    unittest.main()
