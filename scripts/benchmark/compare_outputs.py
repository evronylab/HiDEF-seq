#!/usr/bin/env python3
"""Stream-compare published text; inventory every file and report unvalidated binaries."""
import argparse
import fnmatch
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sys

TEXT_SUFFIXES = {".tsv", ".csv", ".txt", ".bed", ".vcf", ".yaml", ".yml"}
GZIP_SUFFIXES = {".gz", ".bgz", ".bgzf"}


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def scientific_lines(path, vcf_ignore_keys=("fileDate",)):
    compressed = path.suffix in GZIP_SUFFIXES
    opener = gzip.open if compressed else open
    is_vcf = (path.with_suffix("") if compressed else path).suffix == ".vcf"
    with opener(path, "rt", encoding="utf-8", newline=None) as handle:
        for line_number, line in enumerate(handle, 1):
            if is_vcf and any(line.startswith("##" + key + "=") for key in vcf_ignore_keys):
                continue
            yield line_number, line.rstrip("\r\n")


def compare_text(reference, candidate, vcf_ignore_keys=("fileDate",)):
    for left, right in itertools.zip_longest(scientific_lines(reference, vcf_ignore_keys),
                                            scientific_lines(candidate, vcf_ignore_keys)):
        if left is None or right is None:
            return {"status": "failure", "reason": "scientific line count differs"}
        if left[1] != right[1]:
            return {"status": "failure", "reason": "scientific text or order differs",
                    "reference_line": left[0], "candidate_line": right[0]}
    return {"status": "pass", "reason": "scientific text identical"}


def compare_trees(reference, candidate, excludes=(), vcf_ignore_keys=("fileDate",)):
    files = lambda root: {str(p.relative_to(root)): p for p in root.rglob("*") if p.is_file()}
    left, right = files(reference), files(candidate)
    records = []
    for name in sorted(left.keys() | right.keys()):
        record = {"path": name}
        if any(fnmatch.fnmatchcase(name, pattern) for pattern in excludes):
            record.update(status="ignored", reason="explicit metadata whitelist")
        elif name not in left or name not in right:
            record.update(status="failure", reason="missing in " + ("reference" if name not in left else "candidate"))
        elif digest(left[name]) == digest(right[name]):
            record.update(status="pass", reason="byte identical")
        else:
            filename = Path(name)
            suffix = (filename.with_suffix("") if filename.suffix in GZIP_SUFFIXES else filename).suffix
            if suffix in TEXT_SUFFIXES:
                record.update(compare_text(left[name], right[name], vcf_ignore_keys))
            else:
                record.update(status="review", reason="binary differs; requires format-specific validation")
        records.append(record)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--exclude", action="append", default=["*.run_metadata.tsv"],
                        help="Explicit run-metadata glob; repeatable, recorded in report")
    parser.add_argument("--vcf-ignore-key", action="append", default=["fileDate"],
                        help="Exact VCF metadata key; repeatable, recorded in report")
    args = parser.parse_args()
    if not args.reference.is_dir() or not args.candidate.is_dir():
        parser.error("both input directories must exist")
    if args.report.exists():
        parser.error("report already exists")
    records = compare_trees(args.reference, args.candidate, args.exclude, args.vcf_ignore_key)
    counts = {status: sum(row["status"] == status for row in records)
              for status in ("pass", "failure", "review", "ignored")}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(dict(reference=str(args.reference.resolve()),
                                          candidate=str(args.candidate.resolve()),
                                          excludes=args.exclude, vcf_ignore_keys=args.vcf_ignore_key,
                                          counts=counts, files=records), indent=2) + "\n")
    print(json.dumps(counts))
    return 1 if counts["failure"] else 2 if counts["review"] else 0


if __name__ == "__main__":
    sys.exit(main())
