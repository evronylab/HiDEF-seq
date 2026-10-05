#!/usr/bin/env python3
"""Compare experimental FASTA annotation with the complete legacy BED chain.

Run inside the pinned container after compiling bin/annotateCoverage.cpp.
This is a semantic fixture, not a production CPU/memory benchmark. The legacy
reference BED, merged union, one-base expansion, intersections, awk counting,
bgzip and tabix steps are all executed. The R BED writer is not modified.
"""
import argparse
import decimal
import gzip
import itertools
import json
from pathlib import Path
import shlex
import subprocess
import sys


def quote(value):
    return shlex.quote(str(value))


def shell(script, directory):
    result = subprocess.run(["/bin/bash", "-c", "set -euo pipefail\n" + script], cwd=directory,
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        sys.stderr.write(result.stderr)
        result.check_returncode()
    return result


def counts(path):
    result = {}
    for line in path.read_text().splitlines():
        row, context, value = line.split("\t")
        result[(row, context)] = decimal.Decimal(value)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--helper", required=True)
    parser.add_argument("--output-dir", required=True)
    for tool in ("samtools", "seqkit", "bedtools", "bgzip", "tabix"):
        parser.add_argument("--" + tool, default=tool)
    args = parser.parse_args()
    args.helper = str(Path(args.helper).resolve())
    directory = Path(args.output_dir).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    # Every ACGTN triplet occurs in chrAll. Other contigs probe normalization,
    # FASTA order rather than lexical order, lengths 1/2, and linear chrM edges.
    sequences = {
        "chr2": "acgtnRYSWKMBDHVacgtNNACgt",
        "chr1": "ACGTACGTNNACGT",
        "chrOne": "n",
        "chrTwo": "aR",
        "chrM": "tACgN",
        "chrAll": "".join("".join(triplet) for triplet in itertools.product("ACGTN", repeat=3)),
    }
    reference = directory / "reference.fa"
    reference.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in sequences.items()))
    shell(f"{quote(args.samtools)} faidx reference.fa", directory)
    rows = {
        1: ["chr2\t0\t1\t02", "chr2\t1\t6\t2", "chr2\t6\t12\t3",
            f"chr2\t12\t{len(sequences['chr2'])}\t4", "chr1\t0\t4\t2e+00",
            "chr1\t4\t7\t5", "chr1\t10\t14\t1", "chrOne\t0\t1\t9",
            "chrTwo\t0\t2\t7", "chrM\t0\t5\t1", f"chrAll\t0\t{len(sequences['chrAll'])}\t1"],
        2: ["chr2\t2\t5\t1.5", "chr2\t5\t10\t0.5", "chr1\t1\t3\t1000001",
            "chrM\t1\t4\t0.1"],
        3: [],
        4: ["chr1\t1\t3\t2147483648", "chrM\t1\t2\t1000000000000001.5"],
    }
    for row, lines in rows.items():
        (directory / f"{row}.bed").write_text("".join(line + "\n" for line in lines))
    # Exact reference preprocessing from extractGenomeTrinucleotides.
    reference_command = (
        f"{quote(args.seqkit)} seq -u reference.fa | "
        f"{quote(args.seqkit)} replace -s -p '[^ACGTN]' -r N | "
        f"{quote(args.seqkit)} sliding -S '' -s1 -W3 | "
        f"{quote(args.seqkit)} fx2tab -Q | "
        "awk -F '[:\\-\\t]' 'BEGIN {OFS=\"\\t\"}{print $1, $2, $2+1, $4}' | "
        f"{quote(args.bgzip)} -c > reference.trinuc.bed.gz\n"
        f"{quote(args.tabix)} -s1 -b2 -e3 reference.trinuc.bed.gz")
    shell(reference_command, directory)
    union_inputs = " ".join(
        "<(awk -v OFS='\\t' 'FNR==NR{ord[$1]=NR-1;next}{print ord[$1],$1,$2,$3}' "
        f"reference.fa.fai {row}.bed)" for row in rows)
    union_command = (
        f"sort -m -k1,1n -k3,3n -k4,4n {union_inputs} | cut -f2- | "
        f"{quote(args.bedtools)} merge -i stdin > all.bed\n"
        f"{quote(args.bedtools)} makewindows -w 1 -b all.bed | "
        f"{quote(args.bedtools)} intersect -sorted -loj -wa -wb -a stdin -b reference.trinuc.bed.gz -g reference.fa.fai | "
        "cut -f 1-3,7 > all.trinuc.bed")
    shell(union_command, directory)
    results = {}
    observed_contexts = set()
    for row in rows:
        legacy_counts = directory / f"{row}.legacy.counts.tsv"
        candidate_counts = directory / f"{row}.candidate.counts.tsv"
        legacy_bed = directory / f"{row}.legacy.bed.gz"
        candidate_bed = directory / f"{row}.candidate.bed.gz"
        legacy_command = (
            f"{quote(args.bedtools)} intersect -sorted -wa -wb -a {row}.bed -b all.trinuc.bed -g reference.fa.fai | "
            "awk -v OFS='\\t' " + f"-v row_id={row} -v counts={quote(legacy_counts)} " +
            "'{print $5,$6,$7,$4,$8; sum[$8]+=$4} END{if(length(sum)==0){print row_id,\"NA\",0 > counts}"
            "else{for(k in sum){print row_id,k,sum[k] > counts}}}' | " +
            f"{quote(args.bgzip)} -c > {quote(legacy_bed)}\n"
            f"{quote(args.tabix)} -@2 -s1 -b2 -e3 {quote(legacy_bed)}")
        shell(legacy_command, directory)
        helper_command = (
            f"{quote(args.helper)} --bed {row}.bed --fasta reference.fa --fai reference.fa.fai "
            f"--row-id {row} --counts {quote(candidate_counts)} --bed-output - | "
            f"{quote(args.bgzip)} -c > {quote(candidate_bed)}\n"
            f"{quote(args.tabix)} -@2 -s1 -b2 -e3 {quote(candidate_bed)}")
        shell(helper_command, directory)
        with gzip.open(legacy_bed, "rb") as expected, gzip.open(candidate_bed, "rb") as actual:
            expected_bytes, actual_bytes = expected.read(), actual.read()
        if expected_bytes != actual_bytes:
            raise AssertionError(f"row {row}: annotated BED bytes differ (coordinates, order, depth tokens or context)")
        if counts(legacy_counts) != counts(candidate_counts):
            raise AssertionError(f"row {row}: context counts differ")
        counts_only = directory / f"{row}.counts-only.tsv"
        shell(f"{quote(args.helper)} --bed {row}.bed --fasta reference.fa --fai reference.fa.fai "
              f"--row-id {row} --counts {quote(counts_only)}", directory)
        if counts(counts_only) != counts(legacy_counts):
            raise AssertionError(f"row {row}: counts-only mode differs")
        observed_contexts.update(context for _, context in counts(candidate_counts))
        results[str(row)] = {"annotated_bytes": len(actual_bytes), "counts": len(counts(candidate_counts))}
    expected_contexts = {"".join(t) for t in itertools.product("ACGTN", repeat=3)} | {".", "NA"}
    if observed_contexts != expected_contexts:
        raise AssertionError("Not all 125 normalized contexts, dot and empty NA were tested")
    # The legacy reference BED parser mishandles these contig names. Fail safely
    # so a caller can retain its legacy path instead of silently correcting it.
    for name in ("chr-1", "chr:1"):
        unsafe = directory / ("unsafe-" + name.replace(":", "colon") + ".fa")
        unsafe.write_text(f">{name}\nACGT\n")
        shell(f"{quote(args.samtools)} faidx {quote(unsafe)}", directory)
        result = subprocess.run([args.helper, "--bed", str(directory / "3.bed"), "--fasta", str(unsafe),
                                 "--fai", str(unsafe) + ".fai", "--row-id", "3", "--counts", str(directory / "unsafe.tsv")],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0 or "legacy-fallback-required" not in result.stderr:
            raise AssertionError("Legacy-ambiguous reference contig was not rejected")
    results["legacy_ambiguous_contigs"] = "rejected for fallback"
    # A valid first row followed by a bad row must fail after output starts.
    # The caller must detect the failure despite bgzip successfully closing a
    # partial stream; it must never treat absent fresh counts as a valid row.
    late_bad = directory / "late-invalid.bed"
    late_bad.write_text("chr1\t1\t3\t2\nchr1\t2\t4\t1\n")
    late_counts = directory / "late-invalid.counts.tsv"
    helper_args = [args.helper, "--bed", str(late_bad), "--fasta", str(reference),
                   "--fai", str(reference) + ".fai", "--row-id", "late", "--counts", str(late_counts),
                   "--bed-output", "-"]
    rejected = subprocess.run(helper_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    expected_partial = "chr1\t1\t2\t2\tACG\nchr1\t2\t3\t2\tCGT\n"
    if rejected.returncode == 0 or rejected.stdout != expected_partial or late_counts.exists():
        raise AssertionError("Late-invalid row did not fail with partial BED and absent fresh counts")
    if "non-overlapping" not in rejected.stderr:
        raise AssertionError("Late-invalid fixture failed for an unrelated reason")
    partial_gzip = directory / "late-invalid.partial.bed.gz"
    success_marker = directory / "late-invalid.unexpected-success"
    pipeline = " ".join(quote(item) for item in helper_args) + " | " + quote(args.bgzip) + " -c > " + quote(partial_gzip)
    pipeline += "\nprintf 'unexpected success' > " + quote(success_marker)
    result = subprocess.run(["/bin/bash", "-c", "set -euo pipefail\n" + pipeline],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode == 0 or late_counts.exists() or success_marker.exists():
        raise AssertionError("pipefail did not stop the caller after annotation failed")
    with gzip.open(partial_gzip, "rt") as partial:
        if partial.read() != expected_partial:
            raise AssertionError("Late failure fixture did not exercise a valid partial bgzip output")
    results["late_invalid_row"] = "helper and pipefail caller failed; partial BED present, fresh counts and success marker absent"
    (directory / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
