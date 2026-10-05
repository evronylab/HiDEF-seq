#!/usr/bin/env python3
"""Paired real-data replay of complete legacy versus FASTA coverage annotation.

Run inside the pinned container in an allocated compute job. The original
chunk_runs=1e7 BED writer is AST-extracted and executed unchanged in both arms.
Common final-QS loading/row selection is prepared once and reported separately.
Each measured arm includes packet loading, writing, annotation, compression and
indexing; R also reports operation-only CPU with waited-for child processes.
"""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def process_group_rss(group):
    """Sample summed RSS, including shared pages counted in each process."""
    total = 0
    for directory in Path("/proc").iterdir():
        if not directory.name.isdigit():
            continue
        try:
            stat = (directory / "stat").read_text()
            fields = stat[stat.rfind(")") + 2:].split()
            if int(fields[2]) != group:  # Fields begin at proc stat field 3.
                continue
            pages = int((directory / "statm").read_text().split()[1])
            total += pages * os.sysconf("SC_PAGE_SIZE") // 1024
        except (FileNotFoundError, ProcessLookupError, PermissionError, IndexError, ValueError):
            pass
    return total


def measure(command, directory, runner, label, cpus, environment):
    wrapper = [sys.executable, str(runner / "measure.py"), "--label", label,
               "--out", str(directory / "metrics"), "--allocated-cpus", str(cpus), "--", *command]
    peak, samples = 0, 0
    with (directory / "operation.log").open("w") as log:
        process = subprocess.Popen(wrapper, stdout=log, stderr=subprocess.STDOUT,
                                   env=environment, start_new_session=True)
        while process.poll() is None:
            peak = max(peak, process_group_rss(process.pid))
            samples += 1
            time.sleep(0.25)
    record = json.loads((directory / "metrics.json").read_text())
    record.update(sampled_process_group_peak_rss_kib=peak, process_group_samples=samples,
                  rss_sampling_seconds=0.25,
                  rss_note="Summed process RSS sampled every 250ms; shared pages can be counted repeatedly and short peaks missed")
    (directory / "tree-metrics.json").write_text(json.dumps(record, indent=2) + "\n")
    if process.returncode:
        raise RuntimeError(f"{label} failed: see {directory / 'operation.log'}")
    return record


def read_counts(directory):
    result = {}
    for path in sorted(directory.glob("*.reftnc_plus_strand.tsv")):
        values = {}
        for line in path.read_text().splitlines():
            row, context, count = line.split("\t")
            key = (row, context)
            if key in values:
                raise ValueError("Duplicate count key: " + str(path))
            values[key] = Decimal(count)
        result[path.name] = values
    return result


def compare(left, right):
    with (left / "tool-paths.tsv").open() as handle:
        tabix = next(csv.DictReader(handle, delimiter="\t"))["tabix"]
    with (left / "output-paths.tsv").open() as handle:
        expected_rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_counts = sorted(row["annotation_row_id"] + ".reftnc_plus_strand.tsv" for row in expected_rows)
    left_counts, right_counts = read_counts(left), read_counts(right)
    if sorted(left_counts) != expected_counts or sorted(right_counts) != expected_counts:
        raise AssertionError("Context count files missing or unexpected")
    expected = sorted(path.name for path in left.glob("*.bed.gz"))
    actual = sorted(path.name for path in right.glob("*.bed.gz"))
    if not expected or expected != actual:
        raise AssertionError("Published coverage BED file sets differ or are empty")
    comparisons = []
    for name in expected:
        digest, size = hashlib.sha256(), 0
        first_line, tail = None, b""
        with gzip.open(left / name, "rb") as old, gzip.open(right / name, "rb") as new:
            while True:
                a, b = old.read(1024 * 1024), new.read(1024 * 1024)
                if a != b:
                    raise AssertionError("Decompressed BED differs: " + name)
                if not a:
                    break
                if first_line is None:
                    first_line = a.split(b"\n", 1)[0]
                tail = (tail + a)[-4096:]
                size += len(a)
                digest.update(a)
        if not (left / (name + ".tbi")).exists() or not (right / (name + ".tbi")).exists():
            raise AssertionError("Missing tabix index: " + name)
        if subprocess.check_output([tabix, "-l", str(left / name)]) != subprocess.check_output([tabix, "-l", str(right / name)]):
            raise AssertionError("Tabix contig lists differ: " + name)
        for line in (() if first_line is None else (first_line, tail.rstrip(b"\n").rsplit(b"\n", 1)[-1])):
            chromosome, start, end = line.decode().split("\t")[:3]
            region = f"{chromosome}:{int(start) + 1}-{end}"
            old_query = subprocess.check_output([tabix, str(left / name), region])
            new_query = subprocess.check_output([tabix, str(right / name), region])
            if not old_query or old_query != new_query:
                raise AssertionError("Tabix boundary query differs or is unexpectedly empty: " + name)
        comparisons.append(dict(file=name, decompressed_bytes=size, decompressed_sha256=digest.hexdigest()))
    if left_counts != right_counts:
        raise AssertionError("Observed context count maps differ")
    return dict(beds=comparisons, count_files=len(left_counts), status="exact BED bytes and numeric context counts; tabix lists and boundary queries match")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Existing final outputResults QS2")
    parser.add_argument("--baseline-script", required=True, help="Original calculateBurdens.R to freeze and AST-extract")
    parser.add_argument("--helper", required=True, help="Compiled experimental annotateCoverage")
    parser.add_argument("--chromgroup", required=True)
    parser.add_argument("--filtergroup", required=True)
    parser.add_argument("--chromosomes", default="chr22", help="Comma-separated coverage chromosomes, or all")
    parser.add_argument("--output", required=True, help="New benchmark directory")
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--allocated-cpus", type=int, required=True)
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    if args.repeats < 1 or args.allocated_cpus < 1:
        parser.error("repeats and allocated-cpus must be positive")
    output = Path(args.output).resolve()
    if output.exists():
        parser.error("output directory already exists")
    source = Path(args.input).resolve(strict=True)
    baseline = Path(args.baseline_script).resolve(strict=True)
    helper = Path(args.helper).resolve(strict=True)
    output.mkdir(parents=True)
    code = output / "code"
    code.mkdir()
    runner_directory = Path(__file__).resolve().parent
    for name in ("benchmark_coverage_annotation.R", "benchmark_coverage_annotation.py", "measure.py"):
        shutil.copy2(runner_directory / name, code / name)
    shutil.copy2(baseline, code / "calculateBurdens.original.R")
    shutil.copy2(helper, code / "annotateCoverage")
    manifest = dict(input=str(source), input_bytes=source.stat().st_size,
                    input_mtime_ns=source.stat().st_mtime_ns, baseline_script=str(baseline),
                    chromgroup=args.chromgroup, filtergroup=args.filtergroup,
                    chromosomes=args.chromosomes, repeats=args.repeats, allocated_cpus=args.allocated_cpus,
                    reference_trinucleotide_bed="Reuse existing cached BED from original yaml.config; do not rebuild",
                    original_writer_chunk_runs=10000000,
                    code_sha256={path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in code.iterdir()})
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    environment = dict(os.environ)
    runner = code / "benchmark_coverage_annotation.R"
    subprocess.run([args.rscript, "--vanilla", str(runner), "preflight", str(code / "calculateBurdens.original.R")],
                   env=environment, check=True)
    prepared = output / "prepare"
    prepared.mkdir()
    packet = output / "coverage.qs2"
    preparation = measure([args.rscript, "--vanilla", str(runner), "prepare", str(source), str(packet),
                           args.chromgroup, args.filtergroup, args.chromosomes],
                          prepared, code, "prepare-real-coverage", args.allocated_cpus, environment)
    report = dict(preparation=preparation, measurements=[], validation=[])
    for repeat in range(1, args.repeats + 1):
        pair = output / f"pair{repeat:02d}"
        pair.mkdir()
        order = ("legacy", "candidate") if repeat % 2 else ("candidate", "legacy")
        for mode in order:
            directory = pair / mode
            directory.mkdir()
            command = [args.rscript, "--vanilla", str(runner), "run", mode,
                       str(code / "calculateBurdens.original.R"), str(packet), str(code / "annotateCoverage"), str(directory)]
            metric = measure(command, directory, code, f"coverage-{repeat}-{mode}", args.allocated_cpus, environment)
            with (directory / "operation-metrics.tsv").open() as handle:
                phases = list(csv.DictReader(handle, delimiter="\t"))
            report["measurements"].append(dict(pair=repeat, mode=mode, phases=phases, **metric))
            (output / "results.json").write_text(json.dumps(report, indent=2) + "\n")
        comparison = compare(pair / "legacy", pair / "candidate")
        report["validation"].append(dict(pair=repeat, **comparison))
        (output / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(dict(output=str(output), measured_operations=2 * args.repeats,
                          comparisons=len(report["validation"])), indent=2))


if __name__ == "__main__":
    main()
