#!/usr/bin/env python3
"""Measure a foreground workload with GNU time; retain scheduler accounting separately."""
import argparse
import csv
import datetime
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="New output prefix (.json, .tsv, .time.txt)")
    parser.add_argument("--label", required=True)
    parser.add_argument("--time", default="/usr/bin/time", help="GNU time executable")
    parser.add_argument("--allocated-cpus", type=int, default=None)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("a command after -- is required")
    timer = shutil.which(args.time)
    prefix = Path(args.out)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = [Path(str(prefix) + suffix) for suffix in (".json", ".tsv", ".time.txt")]
    if any(path.exists() for path in outputs):
        parser.error("output prefix already exists; use a new prefix")
    allocated = args.allocated_cpus
    if allocated is None and os.environ.get("SLURM_CPUS_PER_TASK"):
        allocated = int(os.environ["SLURM_CPUS_PER_TASK"])
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    # %U + %S includes waited-for descendants. Detached children and remotely
    # submitted jobs are not included: wrap the actual task, not a submit client.
    if timer is not None:
        accounting_method = "GNU time (waited-for descendants)"
        result = subprocess.run([timer, "-o", str(outputs[2]), "-f",
                                 "HIDEF_TIME\\t%U\\t%S\\t%e\\t%M\\t%x", "--", *command])
    else:
        accounting_method = "Linux resource.getrusage(RUSAGE_CHILDREN) delta (waited-for descendants)"
        if not sys.platform.startswith("linux"):
            parser.error("GNU time fallback is supported only on Linux (RSS units differ elsewhere)")
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        wall_start = time.perf_counter()
        result = subprocess.run(command)
        elapsed = time.perf_counter() - wall_start
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        status = result.returncode if result.returncode >= 0 else 128 - result.returncode
        outputs[2].write_text("HIDEF_TIME\t{}\t{}\t{}\t{}\t{}\n".format(
            after.ru_utime - before.ru_utime, after.ru_stime - before.ru_stime,
            elapsed, after.ru_maxrss, status))
    lines = outputs[2].read_text().splitlines()
    metrics = [line.split("\t") for line in lines if line.startswith("HIDEF_TIME\t")]
    if len(metrics) != 1:
        raise RuntimeError("GNU time did not produce exactly one metric record")
    _, user, system, wall, rss, status = metrics[0]
    actual_cpu = float(user) + float(system)
    record = dict(label=args.label, command=command, accounting_method=accounting_method, started_utc=started,
                  ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  user_seconds=float(user), system_seconds=float(system),
                  actual_cpu_seconds=actual_cpu, actual_cpu_hours=actual_cpu / 3600,
                  wall_seconds=float(wall), max_rss_kib=int(rss),
                  allocated_cpus=allocated,
                  allocated_cpu_hours=None if allocated is None else allocated * float(wall) / 3600,
                  exit_status=int(status), wrapper_returncode=result.returncode,
                  slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                  slurm_step_id=os.environ.get("SLURM_STEP_ID"))
    outputs[0].write_text(json.dumps(record, indent=2) + "\n")
    with outputs[1].open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=record, delimiter="\t")
        writer.writeheader()
        writer.writerow({**record, "command": json.dumps(command)})
    return result.returncode if result.returncode >= 0 else 128 - result.returncode


if __name__ == "__main__":
    sys.exit(main())
