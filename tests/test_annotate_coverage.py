#!/usr/bin/env python3
"""Build tiny independent coverage annotation fixtures using interval tools.

Reference contexts, interval intersections, depth counting and indexed BED
outputs provide expected values for test_coverage_annotation.R.
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
    parser.add_argument("--output-dir", required=True)
    for tool in ("samtools", "seqkit", "bedtools", "bgzip", "tabix"):
        parser.add_argument("--" + tool, default=tool)
    args = parser.parse_args()
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
    # Normalize reference bases and enumerate centered trinucleotide contexts.
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
        expected_counts = directory / f"{row}.expected.counts.tsv"
        expected_bed = directory / f"{row}.expected.bed.gz"
        expected_command = (
            f"{quote(args.bedtools)} intersect -sorted -wa -wb -a {row}.bed -b all.trinuc.bed -g reference.fa.fai | "
            "awk -v OFS='\\t' " + f"-v row_id={row} -v counts={quote(expected_counts)} " +
            "'{print $5,$6,$7,$4,$8; sum[$8]+=$4} END{if(length(sum)==0){print row_id,\"NA\",0 > counts}"
            "else{for(k in sum){print row_id,k,sum[k] > counts}}}' | " +
            f"{quote(args.bgzip)} -c > {quote(expected_bed)}\n"
            f"{quote(args.tabix)} -@2 -s1 -b2 -e3 {quote(expected_bed)}")
        shell(expected_command, directory)
        with gzip.open(expected_bed, "rb") as stream:
            annotated_bytes = stream.read()
        observed_contexts.update(context for _, context in counts(expected_counts))
        results[str(row)] = {"annotated_bytes": len(annotated_bytes), "counts": len(counts(expected_counts))}
    expected_contexts = {"".join(t) for t in itertools.product("ACGTN", repeat=3)} | {".", "NA"}
    if observed_contexts != expected_contexts:
        raise AssertionError("Not all 125 normalized contexts, dot and empty NA were tested")
    (directory / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
