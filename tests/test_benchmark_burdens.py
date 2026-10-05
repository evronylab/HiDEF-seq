"""Replay orchestration test without R or biological input."""
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/benchmark/benchmark_burdens.py"
spec = importlib.util.spec_from_file_location("burdens_replay", SCRIPT)
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)


class BurdenReplayTests(unittest.TestCase):
    def test_preserves_input_order_full_script_and_all_product_comparisons(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for mode in ("baseline", "candidate"):
                folder = root / mode / "bin"
                folder.mkdir(parents=True)
                for name in ("calculateBurdens.R", "sharedFunctions.R"):
                    (folder / name).write_text("# fixture\n")
            inputs = [root / "chunk2.qs2", root / "chunk1.qs2"]
            for path in inputs:
                path.write_bytes(b"fixture")
            (root / "config.yaml").write_text("fixture: true\n")
            commands = []

            def fake_run(command, **kwargs):
                commands.append((command, kwargs))
                if "--out" in command:
                    prefix = Path(command[command.index("--out") + 1])
                    prefix.with_suffix(".json").write_text(json.dumps(dict(actual_cpu_seconds=1, max_rss_kib=2)))
                return subprocess.CompletedProcess(command, 0)

            with patch.object(driver.subprocess, "run", side_effect=fake_run):
                driver.main(["--baseline-repo", str(root / "baseline"), "--candidate-repo", str(root / "candidate"),
                             "--config", str(root / "config.yaml"), "--sample", "sample", "--chromgroup", "1-22X",
                             "--filtergroup", "lenient", "--output", str(root / "output"),
                             "--allocated-cpus", "2", "--filter-files", *map(str, inputs)])
            measured = [(c, k) for c, k in commands if "--out" in c]
            self.assertEqual(len(measured), 2)
            for command, kwargs in measured:
                self.assertEqual(command[command.index("-f") + 1], ",".join(map(str, inputs)))
                self.assertEqual(kwargs["cwd"].name, "products")
                self.assertTrue(any(p.endswith("/calculateBurdens.R") for p in command))
            file_comparisons = [c for c, _ in commands if any(p.endswith("/compare_outputs.py") for p in c)]
            self.assertEqual(len(file_comparisons), 1)
            self.assertEqual(file_comparisons[0][file_comparisons[0].index("--exclude") + 1], "calculateBurdens.qs2")
            self.assertTrue(any("compare" in c and c[-1].endswith("qs-comparison.tsv") for c, _ in commands))


if __name__ == "__main__":
    unittest.main()
