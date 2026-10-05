#!/usr/bin/env python3
"""Fresh-process paired germline formatting benchmark, run inside a compute job."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Original complete final QS2")
    parser.add_argument("--baseline-repo", required=True)
    parser.add_argument("--candidate-repo", required=True)
    parser.add_argument("--output", required=True, help="New output directory")
    parser.add_argument("--pairs", type=int, default=3)
    parser.add_argument("--allocated-cpus", type=int, required=True)
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args(argv)
    if args.pairs < 1 or args.allocated_cpus < 1:
        parser.error("pairs and allocated-cpus must be positive")
    source_input = Path(args.input).resolve(strict=True)
    repos = {mode: Path(getattr(args, mode + "_repo")).resolve(strict=True)
             for mode in ("baseline", "candidate")}
    output = Path(args.output).resolve()
    if output.exists():
        parser.error("output already exists")
    output.mkdir(parents=True)
    runner = output / "code" / "benchmark"
    runner.mkdir(parents=True)
    for name in ("benchmark_germline_formatting_fresh.py", "benchmark_germline_formatting_worker.R",
                 "germline_formatting_helpers.R", "measure.py"):
        shutil.copy2(Path(__file__).resolve().parent / name, runner / name)
    bins = {}
    for mode, repo in repos.items():
        bins[mode] = output / "code" / mode / "bin"
        bins[mode].mkdir(parents=True)
        for name in ("outputResults.R", "sharedFunctions.R"):
            shutil.copy2(repo / "bin" / name, bins[mode] / name)
    source_hashes = {str(p.relative_to(output / "code")): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in sorted((output / "code").rglob("*")) if p.is_file()}
    manifest = dict(input=str(source_input), input_bytes=source_input.stat().st_size,
                    input_mtime_ns=source_input.stat().st_mtime_ns, code_sha256=source_hashes,
                    hostname=socket.gethostname(), slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                    pairs=args.pairs, allocated_cpus=args.allocated_cpus,
                    validation="identical() of complete formatted tables, including classes/attributes/order")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if Path("/proc/cpuinfo").exists():
        shutil.copyfile("/proc/cpuinfo", output / "cpuinfo.txt")
    worker = [args.rscript, "--vanilla", str(runner / "benchmark_germline_formatting_worker.R")]
    packet = output / "raw-germline-packet.qs2"
    with (output / "prepare.log").open("w") as log:
        subprocess.run([*worker, "prepare", str(source_input), str(packet),
                        str(bins["baseline"] / "outputResults.R"), str(bins["candidate"] / "outputResults.R")],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    measurements, comparisons = [], []

    def save():
        (output / "results.json").write_text(json.dumps(dict(measurements=measurements,
                                                           comparisons=comparisons), indent=2) + "\n")

    for pair in range(1, args.pairs + 1):
        pair_dir = output / f"pair{pair:02d}"
        pair_dir.mkdir()
        for mode in (("baseline", "candidate") if pair % 2 else ("candidate", "baseline")):
            prefix = pair_dir / mode
            print(f"Pair {pair}: {mode} fresh-process formatting", flush=True)
            with prefix.with_suffix(".log").open("w") as log:
                subprocess.run([sys.executable, str(runner / "measure.py"),
                                "--label", f"format-{pair}-{mode}", "--out", str(prefix) + ".process",
                                "--allocated-cpus", str(args.allocated_cpus), "--", *worker,
                                "format", str(packet), str(bins[mode] / "outputResults.R"),
                                str(bins["candidate"] / "sharedFunctions.R"), str(prefix), mode],
                               stdout=log, stderr=subprocess.STDOUT, check=True, cwd=pair_dir)
            process = json.loads(Path(str(prefix) + ".process.json").read_text())
            with Path(str(prefix) + ".phase.tsv").open() as handle:
                phase = next(csv.DictReader(handle, delimiter="\t"))
            with Path(str(prefix) + ".process_peak.tsv").open() as handle:
                peaks = next(csv.DictReader(handle, delimiter="\t"))
            process["max_rss_kib_across_phase_reset"] = max(
                process["max_rss_kib"], *(float(value) for value in peaks.values()))
            measurements.append(dict(pair=pair, mode=mode, process=process, phase=phase))
            save()
    # No validation warm-up or retained R allocator crosses measured workers.
    for pair in range(1, args.pairs + 1):
        pair_dir = output / f"pair{pair:02d}"
        with (pair_dir / "comparison.log").open("w") as log:
            result = subprocess.run([*worker, "compare", str(pair_dir / "baseline.qs2"),
                                     str(pair_dir / "candidate.qs2"), str(pair_dir / "comparison.tsv")],
                                    stdout=log, stderr=subprocess.STDOUT)
        comparisons.append(dict(pair=pair, exit_status=result.returncode))
        save()
        if result.returncode:
            raise RuntimeError("Exact table comparison failed: " + str(pair_dir))
    print(json.dumps(dict(output=str(output), measurements=len(measurements), comparisons=len(comparisons)), indent=2))


if __name__ == "__main__":
    main()
