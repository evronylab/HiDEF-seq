#!/usr/bin/env python3
"""Freeze and compare whole extraction processes on one allocated compute node.

Run this driver inside the pipeline container, within one SLURM allocation.
All four extractions finish before validation; comparisons are outside timing.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import socket
import subprocess
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def input_identity(path):
    stat = path.stat()
    return dict(path=str(path), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                ctime_ns=stat.st_ctime_ns, device=stat.st_dev, inode=stat.st_ino)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-repo", required=True)
    parser.add_argument("--candidate-repo", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True, help="New benchmark directory")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--allocated-cpus", type=int, required=True)
    parser.add_argument("--pairs", type=int, default=2)
    args = parser.parse_args(argv)
    if args.pairs < 1 or args.allocated_cpus < 1:
        parser.error("pairs and allocated-cpus must be positive")
    repos = {mode: Path(getattr(args, mode + "_repo")).resolve(strict=True)
             for mode in ("baseline", "candidate")}
    bam = Path(args.bam).resolve(strict=True)
    config = Path(args.config).resolve(strict=True)
    output = Path(args.output).resolve()
    if output.exists():
        parser.error("benchmark output directory already exists")
    for repo in repos.values():
        for name in ("extractCalls.R", "sharedFunctions.R"):
            if not (repo / "bin" / name).is_file():
                parser.error("missing pipeline source: " + str(repo / "bin" / name))
    output.mkdir(parents=True)
    runner = output / "code" / "benchmark"
    runner.mkdir(parents=True)
    for name in ("benchmark_extraction.py", "measure.py", "compare_qs2.R"):
        shutil.copy2(Path(__file__).resolve().parent / name, runner / name)
    bins, source_hashes = {}, {}
    for mode, repo in repos.items():
        bins[mode] = output / "code" / mode / "bin"
        bins[mode].mkdir(parents=True)
        for source in sorted((repo / "bin").glob("*.R")):
            shutil.copy2(source, bins[mode] / source.name)
        source_hashes[mode] = {p.name: digest(p) for p in sorted(bins[mode].glob("*.R"))}
    frozen_config = output / "config.yaml"
    shutil.copy2(config, frozen_config)
    bam_identity = input_identity(bam)
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        shutil.copyfile(cpuinfo, output / "cpuinfo.txt")
    manifest = dict(repositories={k: str(v) for k, v in repos.items()},
                    pipeline_r_sha256=source_hashes,
                    benchmark_sha256={p.name: digest(p) for p in sorted(runner.iterdir())},
                    config_original=str(config), config_sha256=digest(frozen_config),
                    bam=bam_identity, sample=args.sample, pairs=args.pairs,
                    hostname=socket.gethostname(), platform=platform.platform(),
                    slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                    slurm_node_list=os.environ.get("SLURM_JOB_NODELIST"),
                    allocated_cpus=args.allocated_cpus,
                    validation="Exact scientific QS2 comparison; only /run_metadata ignored; numerical review fails")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    # Resolve all frozen R syntax before any expensive extraction starts.
    subprocess.run([args.rscript, "--vanilla", "-e",
                    "for (d in commandArgs(TRUE)) for (f in list.files(d, '[.]R$', full.names=TRUE)) parse(f)",
                    *map(str, bins.values())], check=True)
    measurements, validation = [], []

    def save_results():
        (output / "results.json").write_text(json.dumps(
            dict(measurements=measurements, validation=validation), indent=2) + "\n")

    for pair in range(1, args.pairs + 1):
        pair_dir = output / f"pair{pair:02d}"
        pair_dir.mkdir()
        order = ("baseline", "candidate") if pair % 2 else ("candidate", "baseline")
        for mode in order:
            if input_identity(bam) != bam_identity or digest(frozen_config) != manifest["config_sha256"]:
                raise RuntimeError("Benchmark input changed during paired run")
            arm = pair_dir / mode
            arm.mkdir()
            environment = {**os.environ, "PATH": str(bins[mode]) + os.pathsep + os.environ.get("PATH", "")}
            command = [args.rscript, "--vanilla", str(bins[mode] / "extractCalls.R"),
                       "-c", str(frozen_config), "-b", str(bam), "-s", args.sample,
                       "-o", "extractCalls.qs2"]
            print(f"pair {pair}: {mode} extraction", flush=True)
            with (arm / "extract.log").open("w") as log:
                subprocess.run([sys.executable, str(runner / "measure.py"),
                                "--label", f"extract-{pair}-{mode}", "--out", str(arm / "metrics"),
                                "--allocated-cpus", str(args.allocated_cpus), "--", *command],
                               cwd=arm, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
            measurements.append(dict(pair=pair, mode=mode,
                                     **json.loads((arm / "metrics.json").read_text())))
            save_results()
    if input_identity(bam) != bam_identity:
        raise RuntimeError("BAM changed during paired run")
    # Staging one top-level component at a time bounds comparison memory.
    comparator = [args.rscript, "--vanilla", str(runner / "compare_qs2.R")]
    for pair in range(1, args.pairs + 1):
        pair_dir = output / f"pair{pair:02d}"
        report = pair_dir / "comparison.tsv"
        scratch = pair_dir / "reference-stage"
        with (pair_dir / "comparison.log").open("w") as log:
            subprocess.run([*comparator, "stage", str(pair_dir / "baseline" / "extractCalls.qs2"),
                            str(scratch)], stdout=log, stderr=subprocess.STDOUT, check=True)
            result = subprocess.run([*comparator, "compare", str(scratch),
                                     str(pair_dir / "candidate" / "extractCalls.qs2"), str(report)],
                                    stdout=log, stderr=subprocess.STDOUT)
        validation.append(dict(pair=pair, status=result.returncode, report=str(report)))
        save_results()
        if result.returncode:
            raise RuntimeError("Scientific comparison failed or requires numerical review: " + str(report))
    print(json.dumps(dict(output=str(output), measurements=len(measurements),
                          comparisons=len(validation)), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
