#!/usr/bin/env python3
"""Paired replay from an immutable real-coverage packet: current C++ versus R.

Execute on a compute node in the pinned HiDEF-seq container. Each new result
directory freezes all measured source. The original writer remains unchanged.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess


def verify_indexes(directory):
    """Check every saved index against an independent rebuild of its own BED."""
    with (directory / "tool-paths.tsv").open() as handle:
        tabix = next(csv.DictReader(handle, delimiter="\t"))["tabix"]
    rebuilt = directory / "rebuilt-indexes"
    rebuilt.mkdir()
    result = []
    for bed in sorted(directory.glob("*.bed.gz")):
        link = rebuilt / bed.name
        link.symlink_to(bed.resolve())
        subprocess.run([tabix, "-@2", "-s1", "-b2", "-e3", str(link)], check=True)
        with gzip.open(str(bed) + ".tbi", "rb") as old, gzip.open(str(link) + ".tbi", "rb") as new:
            payload = old.read()
            if payload != new.read():
                raise AssertionError("Tabix own-file rebuild differs: " + str(bed))
        result.append(dict(file=bed.name, status="PASS", index_payload_sha256=hashlib.sha256(payload).hexdigest()))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet", required=True)
    parser.add_argument("--baseline-script", required=True)
    parser.add_argument("--candidate-script", required=True)
    parser.add_argument("--shared-functions", required=True)
    parser.add_argument("--helper", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--chromosomes", default="all")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--allocated-cpus", type=int, required=True)
    parser.add_argument("--historical-directory", help="Optional original replay outputs for an additional exact comparison")
    parser.add_argument("--consumer-test", help="Optional test_coverage_annotation_consumer.R for ordered downstream/QS2 checks")
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    code = output / "code"
    code.mkdir()
    original = Path(__file__).resolve().parent
    for name in ("benchmark_r_coverage_annotation.py", "benchmark_r_coverage_annotation.R",
                 "benchmark_coverage_annotation.py", "measure.py"):
        shutil.copy2(original / name, code / name)
    sources = {"calculateBurdens.original.R": args.baseline_script,
               "calculateBurdens.candidate.R": args.candidate_script,
               "sharedFunctions.R": args.shared_functions, "annotateCoverage": args.helper}
    if args.consumer_test:
        sources["test_coverage_annotation_consumer.R"] = args.consumer_test
    for name, path in sources.items():
        shutil.copy2(path, code / name)
    packet = Path(args.packet).resolve(strict=True)
    stat = packet.stat()
    manifest = dict(packet=str(packet), packet_bytes=stat.st_size,
                    packet_mtime_ns=stat.st_mtime_ns, chromosomes=args.chromosomes,
                    source_paths=sources, original_writer_chunk_runs=10000000,
                    code_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in code.iterdir()},
                    repeats=args.repeats, allocated_cpus=args.allocated_cpus)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    spec = importlib.util.spec_from_file_location("paired_helpers", code / "benchmark_coverage_annotation.py")
    helpers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helpers)
    runner = code / "benchmark_r_coverage_annotation.R"
    environment = dict(os.environ)
    subprocess.run(["Rscript", "--vanilla", str(runner), "writer-contract",
                    str(code / "calculateBurdens.original.R"), str(code / "calculateBurdens.candidate.R")],
                   check=True, env=environment)
    report = dict(measurements=[], validation=[])
    if args.chromosomes != "all":
        preparation = output / "prepare"
        preparation.mkdir()
        reduced = output / "coverage.qs2"
        report["preparation"] = helpers.measure(["Rscript", "--vanilla", str(runner), "subset",
            str(packet), str(reduced), args.chromosomes], preparation, code,
            "subset-existing-real-packet", args.allocated_cpus, environment)
        packet = reduced
    for repeat in range(1, args.repeats + 1):
        pair = output / f"pair{repeat:02d}"
        pair.mkdir()
        for arm in (("candidate", "r") if repeat % 2 else ("r", "candidate")):
            directory = pair / arm
            directory.mkdir()
            implementation = code / ("sharedFunctions.R" if arm == "r" else "annotateCoverage")
            command = ["Rscript", "--vanilla", str(runner), "run", arm,
                       str(code / "calculateBurdens.original.R"), str(packet), str(implementation), str(directory)]
            metric = helpers.measure(command, directory, code, f"coverage-{repeat}-{arm}",
                                     args.allocated_cpus, environment)
            with (directory / "operation-metrics.tsv").open() as handle:
                phases = list(csv.DictReader(handle, delimiter="\t"))
            report["measurements"].append(dict(pair=repeat, mode=arm, phases=phases, **metric))
            (output / "results.json").write_text(json.dumps(report, indent=2) + "\n")
        report["validation"].append(dict(pair=repeat, **helpers.compare(pair / "candidate", pair / "r"),
            cpp_indexes=verify_indexes(pair / "candidate"), r_indexes=verify_indexes(pair / "r")))
        if args.historical_directory:
            report["historical_validation"] = helpers.compare(Path(args.historical_directory).resolve(), pair / "r")
        if args.consumer_test:
            consumer = pair / "consumer"
            subprocess.run(["Rscript", "--vanilla", str(code / "test_coverage_annotation_consumer.R"),
                str(code / "calculateBurdens.candidate.R"), str(code / "sharedFunctions.R"),
                str(pair / "candidate"), str(pair / "r"), str(consumer)], check=True, env=environment)
            report["validation"][-1]["ordered_context_consumer"] = (consumer / "PASS.txt").read_text().strip()
        (output / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    final = Path(args.packet).resolve().stat()
    if (stat.st_size, stat.st_mtime_ns) != (final.st_size, final.st_mtime_ns):
        raise AssertionError("Input coverage packet changed during replay")
    for name, expected in manifest["code_sha256"].items():
        if hashlib.sha256((code / name).read_bytes()).hexdigest() != expected:
            raise AssertionError("Frozen replay source changed: " + name)
    report["unchanged_input_and_sources"] = True
    (output / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(dict(output=str(output), status="PASS", pairs=args.repeats)))


if __name__ == "__main__":
    main()
