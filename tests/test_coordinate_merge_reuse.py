#!/usr/bin/env python3
"""Probe shared per-run coordinate-sort reuse with the pinned pbmerge.

No workflow changes. Compare current pbmerge(raw inputs)+samtools sort with
samtools sort(each raw input)+pbmerge, including exact SAM record order. A
coordinate-indexable BAM alone is insufficient: any ordered-record mismatch
rejects removal of the final sort. Results JSON is a probe report, not approval.
"""
import argparse
import json
from pathlib import Path
import subprocess

from test_split_bam_dispatch import record


def run(command):
    result = subprocess.run([str(value) for value in command], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    if result.returncode:
        raise RuntimeError(" ".join(str(value) for value in command) + "\n" + result.stderr)
    return result.stdout


def header_without_programs(samtools, bam):
    return [line for line in run([samtools, "view", "--no-PG", "-H", bam]).splitlines()
            if not line.startswith("@PG\t")]


def first_mismatch(expected, actual):
    for index, (left, right) in enumerate(zip(expected, actual)):
        if left != right:
            return {"record": index + 1, "legacy": left, "candidate": right}
    if len(expected) != len(actual):
        return {"record_count": [len(expected), len(actual)]}
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--pbmerge", default="pbmerge")
    parser.add_argument("--pbindex", default="pbindex")
    args = parser.parse_args()
    directory = Path(args.output_dir).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    header = ["@HD\tVN:1.6\tSO:unknown\tpb:3.0.1", "@SQ\tSN:chr1\tLN:1000",
              "@SQ\tSN:chr2\tLN:1000",
              "@RG\tID:00000001\tPL:PACBIO\tPU:movieA\tDS:READTYPE=CCS",
              "@RG\tID:00000002\tPL:PACBIO\tPU:movieB\tDS:READTYPE=CCS",
              "@CO\tCoordinate ties and unmapped tails must retain legacy record order"]
    rows = {
        "A": [record("movieA", 10, 10, 16), record("movieA", 20, 5),
              record("movieA", 30, 10), record("movieA", 40, 10, 256),
              record("movieA", 50, 10, 2064), record("movieA", 60, 20).replace("\tchr1\t", "\tchr2\t"),
              record("movieA", 70, 0, 4), record("movieA", 80, 0, 4)],
        "B": [record("movieB", 10, 10), record("movieB", 20, 10, 16),
              record("movieB", 30, 5), record("movieB", 40, 10, 2048),
              record("movieB", 50, 10, 272), record("movieB", 60, 20).replace("\tchr1\t", "\tchr2\t"),
              record("movieB", 70, 0, 4), record("movieB", 80, 0, 4)],
        "C": [record("movieA", 90, 10, 16), record("movieA", 100, 5),
              record("movieA", 110, 10), record("movieA", 120, 0, 4)],
    }
    raw, sorted_inputs = {}, {}
    for name, records in rows.items():
        sam = directory / f"{name}.raw.sam"
        sam.write_text("\n".join(header + records) + "\n")
        raw[name] = directory / f"{name}.raw.bam"
        sorted_inputs[name] = directory / f"{name}.verify.sorted.bam"
        run([args.samtools, "view", "--no-PG", "-b", "-o", raw[name], sam])
        run([args.pbindex, raw[name]])
        # Exactly the current VerifyBAMID sorting options.
        run([args.samtools, "sort", "-@", "2", "--write-index", "-o", sorted_inputs[name], raw[name]])
        run([args.pbindex, sorted_inputs[name]])
    results = {}
    for label, ordering in [("single", ["A"]), ("AB", ["A", "B"]), ("BA", ["B", "A"]),
                            ("ABC", ["A", "B", "C"]), ("CBA", ["C", "B", "A"])]:
        merged = directory / f"{label}.legacy.unsorted.bam"
        legacy = directory / f"{label}.legacy.sorted.bam"
        candidate = directory / f"{label}.candidate.sorted.bam"
        run([args.pbmerge, "-o", merged] + [raw[name] for name in ordering])
        # Exactly the current postmerge sort options, with explicit output name.
        run([args.samtools, "sort", "-@", "4", "-m", "4G", "-o", legacy, merged])
        run([args.pbmerge, "-o", candidate] + [sorted_inputs[name] for name in ordering])
        for bam in (legacy, candidate):
            run([args.samtools, "quickcheck", bam])
            run([args.samtools, "index", "-@", "4", bam])
            run([args.pbindex, bam])
        expected = run([args.samtools, "view", "--no-PG", legacy]).splitlines()
        actual = run([args.samtools, "view", "--no-PG", candidate]).splitlines()
        expected_header = header_without_programs(args.samtools, legacy)
        actual_header = header_without_programs(args.samtools, candidate)
        results[label] = {"input_order": ordering, "records": len(expected),
                          "exact_ordered_records": expected == actual,
                          "record_multiset_equal": sorted(expected) == sorted(actual),
                          "header_equal_except_PG": expected_header == actual_header,
                          "first_mismatch": first_mismatch(expected, actual),
                          "legacy_header": expected_header, "candidate_header": actual_header,
                          "BAI_PBI_quickcheck": "passed"}
    equivalent = all(case["exact_ordered_records"] and case["header_equal_except_PG"]
                     for case in results.values())
    report = {"cases": results, "all_fixtures_equivalent": equivalent,
              "recommendation": ("Real unsorted data and VerifyBAMID replay still required"
                                 if equivalent else "Reject final-sort removal; retain current ordering")}
    (directory / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
