#!/usr/bin/env python3
"""Collect Slurm allocation and step accounting without double-counting them."""
import argparse
import csv
import datetime
import json
from pathlib import Path
import subprocess

FIELDS = ["JobIDRaw", "JobName", "State", "ExitCode", "ElapsedRaw", "AllocCPUS",
          "CPUTimeRAW", "TotalCPU", "UserCPU", "SystemCPU", "MaxRSS", "AveRSS",
          "ReqMem", "AllocTRES", "NodeList"]


def cpu_seconds(value):
    if not value or value in ("Unknown", "INVALID", "N/A"):
        return None
    days, clock = value.split("-", 1) if "-" in value else ("0", value)
    parts = [float(part) for part in clock.split(":")]
    seconds = 0.0
    for part in parts:
        seconds = seconds * 60 + part
    return int(days) * 86400 + seconds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", required=True, help="Comma-separated job IDs")
    parser.add_argument("--out", required=True, help="New output prefix")
    args = parser.parse_args()
    prefix = Path(args.out)
    outputs = [Path(str(prefix) + suffix) for suffix in (".json", ".tsv", ".sacct.txt")]
    if any(path.exists() for path in outputs):
        parser.error("output prefix already exists")
    command = ["sacct", "-j", args.jobs, "--parsable2", "--noheader", "--units=K",
               "--format=" + ",".join(FIELDS)]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    records = []
    for line in result.stdout.splitlines():
        values = line.split("|")
        if len(values) != len(FIELDS):
            raise ValueError("Unexpected sacct field count: " + repr(line))
        record = dict(zip(FIELDS, values))
        record["actual_cpu_seconds"] = cpu_seconds(record["TotalCPU"])
        record["allocated_cpu_hours"] = float(record["CPUTimeRAW"]) / 3600 if record["CPUTimeRAW"] else None
        record["record_kind"] = "step" if "." in record["JobIDRaw"] else "allocation"
        records.append(record)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs[2].write_text(result.stdout)
    outputs[0].write_text(json.dumps(dict(collected_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                         command=command, records=records), indent=2) + "\n")
    with outputs[1].open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS + ["actual_cpu_seconds", "allocated_cpu_hours", "record_kind"], delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    main()
