#!/usr/bin/env python3
"""Experimental paired comparison of separate filtering processes and one shared R session.

No production workflow changes. Run the entire driver inside an allocated
compute job using the pipeline container. Measurements exclude validation.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def run_separate(args):
    directory = Path(args.output)
    for index, group in enumerate(args.filtergroups, 1):
        work = directory / f"group{index:03d}"
        work.mkdir()
        command = [args.rscript, "--vanilla", str(Path(args.script_bin) / "filterCalls.R"),
                   "-c", args.config, "-f", args.extract, "-s", args.sample,
                   "-g", args.chromgroup, "-v", group, "-o", "output.qs2"]
        with (work / "filter.log").open("w") as log:
            subprocess.run(command, cwd=work, stdout=log, stderr=subprocess.STDOUT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--extract", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--chromgroup", required=True)
    parser.add_argument("--filtergroups", nargs="+", required=True)
    parser.add_argument("--output", required=True, help="New benchmark directory")
    parser.add_argument("--repo", default=str(Path(__file__).resolve().parents[2]))
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--allocated-cpus", type=int, required=True)
    parser.add_argument("--separate-worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--script-bin", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if len(args.filtergroups) < 2 or len(set(args.filtergroups)) != len(args.filtergroups):
        parser.error("provide at least two distinct filtergroups")
    if args.repeats < 1 or args.allocated_cpus < 1:
        parser.error("repeats and allocated-cpus must be positive")
    if args.separate_worker:
        run_separate(args)
        return 0
    repo = Path(args.repo).resolve(strict=True)
    output = Path(args.output).resolve()
    config = Path(args.config).resolve(strict=True)
    extract = Path(args.extract).resolve(strict=True)
    if output.exists():
        parser.error("benchmark output directory already exists")
    output.mkdir(parents=True)
    # Freeze all R pipeline source before pairing; ongoing development must not
    # make separate/shared arms run different source code. Input QS is not copied.
    script_bin = output / "code" / "bin"
    script_bin.mkdir(parents=True)
    hashes = {}
    for source in sorted((repo / "bin").glob("*.R")):
        shutil.copy2(source, script_bin / source.name)
        hashes[source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
    config_copy = output / "config.yaml"
    shutil.copy2(config, config_copy)
    runner_dir = output / "code" / "benchmark"
    runner_dir.mkdir()
    for name in ("benchmark_filter_chain.py", "benchmark_filter_chain.R", "measure.py", "compare_qs2.R"):
        shutil.copy2(repo / "scripts" / "benchmark" / name, runner_dir / name)
    environment = {**os.environ, "PATH": str(script_bin) + os.pathsep + os.environ.get("PATH", "")}
    config_hash = hashlib.sha256(config_copy.read_bytes()).hexdigest()
    metadata = dict(config_original=str(config), config_sha256=config_hash, extract=str(extract),
                    extract_bytes=extract.stat().st_size, extract_mtime_ns=extract.stat().st_mtime_ns,
                    sample=args.sample, chromgroup=args.chromgroup, filtergroups=args.filtergroups,
                    pipeline_r_sha256=hashes, repeats=args.repeats, allocated_cpus=args.allocated_cpus,
                    validation="strict QS2 comparison; no metadata or float differences automatically accepted")
    (output / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
    # Fail before expensive work if prepared reference/coverage paths are absent.
    preflight = """
        suppressPackageStartupMessages(library(configr))
        args <- commandArgs(TRUE)
        x <- suppressWarnings(read.config(args[[1]]))
        stopifnot(!is.null(x$reference_summary_file), file.exists(x$reference_summary_file),
                  length(x$germline_coverage_filters) > 0L)
        selected <- Filter(function(z) as.character(z$sample_id) == args[[2]], x$samples)
        stopifnot(length(selected) == 1L)
        groups <- Filter(function(z) as.character(z$filtergroup) %in% args[-c(1, 2)], x$filtergroups)
        stopifnot(length(groups) == length(args) - 2L)
        thresholds <- unique(vapply(groups, function(z) as.numeric(z$min_germlineBAM_TotalReads), numeric(1)))
        prepared <- Filter(function(z) as.character(z$individual_id) == as.character(selected[[1]]$individual_id) &&
                             as.numeric(z$threshold) %in% thresholds, x$germline_coverage_filters)
        stopifnot(length(prepared) == length(thresholds),
                  all(vapply(prepared, function(z) file.exists(z$file), logical(1))))
    """
    subprocess.run([args.rscript, "--vanilla", "-e", preflight, str(config_copy), args.sample, *args.filtergroups],
                   check=True, env=environment)
    measurements = []
    validation = []
    for repeat in range(1, args.repeats + 1):
        pair = output / f"pair{repeat:02d}"
        pair.mkdir()
        order = ("separate", "shared") if repeat % 2 else ("shared", "separate")
        for mode in order:
            arm = pair / mode
            arm.mkdir()
            common = [str(script_bin), str(config_copy), str(extract), args.sample,
                      args.chromgroup, str(arm), *args.filtergroups]
            if mode == "shared":
                command = [args.rscript, "--vanilla", str(runner_dir / "benchmark_filter_chain.R"), *common]
            else:
                command = [sys.executable, str(runner_dir / "benchmark_filter_chain.py"), "--separate-worker",
                           "--script-bin", str(script_bin), "--config", str(config_copy), "--extract", str(extract),
                           "--sample", args.sample, "--chromgroup", args.chromgroup, "--output", str(arm),
                           "--rscript", args.rscript, "--allocated-cpus", str(args.allocated_cpus),
                           "--filtergroups", *args.filtergroups]
            with (arm / "chain.log").open("w") as log:
                subprocess.run([sys.executable, str(runner_dir / "measure.py"), "--label", f"filter-{repeat}-{mode}",
                                "--out", str(arm / "metrics"), "--allocated-cpus", str(args.allocated_cpus),
                                "--", *command], env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
            measurements.append(dict(pair=repeat, mode=mode, **json.loads((arm / "metrics.json").read_text())))
        for index, group in enumerate(args.filtergroups, 1):
            scratch = pair / f"reference-stage-{index:03d}"
            report = pair / f"comparison-{index:03d}.tsv"
            comparator = [args.rscript, "--vanilla", str(runner_dir / "compare_qs2.R")]
            with (pair / f"comparison-{index:03d}.log").open("w") as log:
                subprocess.run([*comparator, "stage", str(pair / "separate" / f"group{index:03d}" / "output.qs2"),
                                str(scratch)], stdout=log, stderr=subprocess.STDOUT, check=True, env=environment)
                result = subprocess.run([*comparator, "compare", str(scratch),
                                         str(pair / "shared" / f"group{index:03d}" / "output.qs2"), str(report)],
                                        stdout=log, stderr=subprocess.STDOUT, env=environment)
            validation.append(dict(pair=repeat, filtergroup=group, status=result.returncode, report=str(report)))
            if result.returncode:
                (output / "results.json").write_text(json.dumps(dict(measurements=measurements, validation=validation), indent=2) + "\n")
                raise RuntimeError("Scientific output comparison failed or requires numerical review: " + str(report))
        (output / "results.json").write_text(json.dumps(dict(measurements=measurements, validation=validation), indent=2) + "\n")
    print(json.dumps(dict(output=str(output), measurements=len(measurements), comparisons=len(validation)), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
