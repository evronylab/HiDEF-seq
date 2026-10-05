#!/usr/bin/env python3
"""Replay complete burdens, including existing BED production, on real filtered chunks.

Use one compute allocation/container. Inputs are an explicitly ordered list of
original filterCalls QS2 files from one sample/chromgroup/filtergroup. This driver
does not synthesize chunks, batch calls, skip stages, or alter pipeline scripts.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys


def identity(path):
    stat = path.stat()
    return dict(path=str(path), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                ctime_ns=stat.st_ctime_ns, device=stat.st_dev, inode=stat.st_ino)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-repo", required=True)
    parser.add_argument("--candidate-repo", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--chromgroup", required=True)
    parser.add_argument("--filtergroup", required=True)
    parser.add_argument("--filter-files", nargs="+", required=True, help="Actual chunks, in original pipeline order")
    parser.add_argument("--output", required=True, help="New benchmark directory")
    parser.add_argument("--pairs", type=int, default=1)
    parser.add_argument("--allocated-cpus", type=int, required=True)
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args(argv)
    if args.pairs < 1 or args.allocated_cpus < 1:
        parser.error("pairs and allocated-cpus must be positive")
    repos = {mode: Path(getattr(args, mode + "_repo")).resolve(strict=True)
             for mode in ("baseline", "candidate")}
    config = Path(args.config).resolve(strict=True)
    inputs = [Path(path).resolve(strict=True) for path in args.filter_files]
    if len(set(inputs)) != len(inputs) or any("," in str(path) for path in inputs):
        parser.error("provide unique chunk files with no commas in their paths")
    output = Path(args.output).resolve()
    if output.exists():
        parser.error("output directory already exists")
    output.mkdir(parents=True)
    runner = output / "code" / "benchmark"
    runner.mkdir(parents=True)
    for name in ("benchmark_burdens.py", "measure.py", "compare_qs2.R", "compare_outputs.py"):
        shutil.copy2(Path(__file__).resolve().parent / name, runner / name)
    bins = {}
    for mode, repo in repos.items():
        bins[mode] = output / "code" / mode / "bin"
        bins[mode].mkdir(parents=True)
        for source in sorted((repo / "bin").glob("*.R")):
            shutil.copy2(source, bins[mode] / source.name)
        for name in ("calculateBurdens.R", "sharedFunctions.R"):
            if not (bins[mode] / name).is_file():
                parser.error("missing frozen source: " + str(bins[mode] / name))
    frozen_config = output / "config.yaml"
    shutil.copy2(config, frozen_config)
    input_identities = [identity(path) for path in inputs]
    manifest = dict(inputs=input_identities, sample=args.sample, chromgroup=args.chromgroup,
                    filtergroup=args.filtergroup, pairs=args.pairs, allocated_cpus=args.allocated_cpus,
                    config_original=str(config), config_sha256=hashlib.sha256(frozen_config.read_bytes()).hexdigest(),
                    code_sha256={str(p.relative_to(output / "code")): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in sorted((output / "code").rglob("*")) if p.is_file()},
                    hostname=socket.gethostname(), slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                    validation="All QS components, exact types/order (only run_metadata ignored); complete BED text and index inventory")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if Path("/proc/cpuinfo").exists():
        shutil.copyfile("/proc/cpuinfo", output / "cpuinfo.txt")
    subprocess.run([args.rscript, "--vanilla", "-e",
                    "for (d in commandArgs(TRUE)) for (f in list.files(d, '[.]R$', full.names=TRUE)) parse(f)",
                    *map(str, bins.values())], check=True)
    measurements, comparisons = [], []

    def save():
        (output / "results.json").write_text(json.dumps(dict(measurements=measurements, comparisons=comparisons), indent=2) + "\n")

    for pair in range(1, args.pairs + 1):
        for mode in (("baseline", "candidate") if pair % 2 else ("candidate", "baseline")):
            if [identity(path) for path in inputs] != input_identities:
                raise RuntimeError("Input chunks changed during replay")
            arm = output / f"pair{pair:02d}" / mode
            products = arm / "products"
            products.mkdir(parents=True)
            environment = {**os.environ, "PATH": str(bins[mode]) + os.pathsep + os.environ.get("PATH", "")}
            command = [args.rscript, "--vanilla", str(bins[mode] / "calculateBurdens.R"),
                       "-c", str(frozen_config), "-s", args.sample, "-g", args.chromgroup,
                       "-v", args.filtergroup, "-f", ",".join(map(str, inputs)), "-o", "calculateBurdens.qs2"]
            print(f"Pair {pair}: {mode}, {len(inputs)} original filtered chunks", flush=True)
            with (arm / "burdens.log").open("w") as log:
                subprocess.run([sys.executable, str(runner / "measure.py"), "--label", f"burdens-{pair}-{mode}",
                                "--out", str(arm / "metrics"), "--allocated-cpus", str(args.allocated_cpus),
                                "--", *command], cwd=products, env=environment,
                               stdout=log, stderr=subprocess.STDOUT, check=True)
            measurements.append(dict(pair=pair, mode=mode, **json.loads((arm / "metrics.json").read_text())))
            save()
    if [identity(path) for path in inputs] != input_identities:
        raise RuntimeError("Input chunks changed during replay")
    for pair in range(1, args.pairs + 1):
        pair_dir = output / f"pair{pair:02d}"
        baseline = pair_dir / "baseline" / "products"
        candidate = pair_dir / "candidate" / "products"
        scratch = pair_dir / "reference-components"
        comparator = [args.rscript, "--vanilla", str(runner / "compare_qs2.R")]
        with (pair_dir / "comparison.log").open("w") as log:
            subprocess.run([*comparator, "stage", str(baseline / "calculateBurdens.qs2"), str(scratch)],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
            qs = subprocess.run([*comparator, "compare", str(scratch), str(candidate / "calculateBurdens.qs2"),
                                 str(pair_dir / "qs-comparison.tsv")], stdout=log, stderr=subprocess.STDOUT)
            # QS is checked above; no BED/TSV/index outputs are excluded.
            files = subprocess.run([sys.executable, str(runner / "compare_outputs.py"), str(baseline),
                                    str(candidate), "--exclude", "calculateBurdens.qs2", "--report",
                                    str(pair_dir / "files-comparison.json")], stdout=log, stderr=subprocess.STDOUT)
        comparisons.append(dict(pair=pair, qs_status=qs.returncode, files_status=files.returncode))
        save()
        if qs.returncode or files.returncode:
            raise RuntimeError("Scientific comparison failed or requires review: " + str(pair_dir))
    print(json.dumps(dict(output=str(output), measurements=len(measurements), comparisons=len(comparisons)), indent=2))


if __name__ == "__main__":
    main()
