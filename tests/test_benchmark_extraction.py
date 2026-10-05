"""Small orchestration tests; never run R or read biological data."""
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/benchmark/benchmark_extraction.py"
spec = importlib.util.spec_from_file_location("benchmark_extraction", SCRIPT)
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)


class ExtractionDriverTests(unittest.TestCase):
    def run_driver(self, comparison_status=0):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        for mode in ("baseline", "candidate"):
            folder = root / mode / "bin"
            folder.mkdir(parents=True)
            for name in ("extractCalls.R", "sharedFunctions.R"):
                (folder / name).write_text(f"# {mode} frozen source\n")
        (root / "input.bam").write_bytes(b"fixture")
        (root / "config.yaml").write_text("fixture: true\n")
        commands = []

        def fake_run(command, **kwargs):
            commands.append((command, kwargs))
            if any(str(part).endswith("/measure.py") for part in command):
                prefix = Path(command[command.index("--out") + 1])
                prefix.with_suffix(".json").write_text(json.dumps(dict(actual_cpu_seconds=1, max_rss_kib=2)))
                # Source changes mid-job must not affect the next extraction.
                (root / "candidate/bin/sharedFunctions.R").write_text("# changed live source\n")
            status = comparison_status if "compare" in command else 0
            return subprocess.CompletedProcess(command, status)

        args = ["--baseline-repo", str(root / "baseline"), "--candidate-repo", str(root / "candidate"),
                "--bam", str(root / "input.bam"), "--config", str(root / "config.yaml"),
                "--sample", "fixture", "--output", str(root / "output"), "--allocated-cpus", "2"]
        with patch.object(driver.platform, "platform", return_value="test-platform"), \
                patch.object(driver.subprocess, "run", side_effect=fake_run):
            if comparison_status:
                with self.assertRaisesRegex(RuntimeError, "Scientific comparison failed"):
                    driver.main(args)
            else:
                self.assertEqual(driver.main(args), 0)
        return root, commands

    def test_freezes_sources_alternates_pairs_and_validates_after_measurement(self):
        root, commands = self.run_driver()
        measures = [(c, k) for c, k in commands if any(str(p).endswith("/measure.py") for p in c)]
        self.assertEqual([k["cwd"].name for _, k in measures], ["baseline", "candidate", "candidate", "baseline"])
        for command, kwargs in measures:
            frozen_bin = root / "output/code" / kwargs["cwd"].name / "bin"
            self.assertEqual(kwargs["env"]["PATH"].split(":")[0], str(frozen_bin))
            self.assertIn(str(frozen_bin / "extractCalls.R"), command)
        self.assertEqual((root / "output/code/candidate/bin/sharedFunctions.R").read_text(),
                         "# candidate frozen source\n")
        self.assertTrue(all("stage" not in c and "compare" not in c for c, _ in commands[:5]))
        results = json.loads((root / "output/results.json").read_text())
        self.assertEqual(len(results["measurements"]), 4)
        self.assertEqual([x["status"] for x in results["validation"]], [0, 0])

    def test_numerical_review_is_failure_and_preserves_measurements(self):
        root, _ = self.run_driver(comparison_status=2)
        results = json.loads((root / "output/results.json").read_text())
        self.assertEqual(len(results["measurements"]), 4)
        self.assertEqual(results["validation"][0]["status"], 2)


if __name__ == "__main__":
    unittest.main()
