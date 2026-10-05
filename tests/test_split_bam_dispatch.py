#!/usr/bin/env python3
"""Run inside the pinned container; compare the dispatcher with zmwfilter.

The caller compiles splitBamByZmw.cpp, activates the PacBio environment, and
provides a new output directory. No Python packages beyond stdlib are needed.
Missing-tag/QNAME disagreement probes verify the observed legacy behavior:
missing zm is indexed as zero and an existing zm takes precedence over QNAME.
"""
import argparse
import gzip
import json
from pathlib import Path
import struct
import subprocess


def run(command, **kwargs):
    return subprocess.run([str(item) for item in command], check=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kwargs)


def fixture(samtools, directory, name, records, reference_length=1000):
    sam = directory / (name + ".sam")
    header = ["@HD\tVN:1.6\tSO:coordinate", f"@SQ\tSN:chr1\tLN:{reference_length}",
              "@RG\tID:00000001\tPL:PACBIO\tPU:movieA\tDS:READTYPE=CCS",
              "@RG\tID:00000002\tPL:PACBIO\tPU:movieB\tDS:READTYPE=CCS",
              "@CO\tHeader comment survives dispatch"]
    sam.write_text("\n".join(header + records) + "\n")
    bam = directory / (name + ".bam")
    run([samtools, "view", "--no-PG", "-b", "-o", bam, sam])
    return bam, header


def record(movie, hole, position, flag=0, zm=True, qname_hole=None):
    qname = f"{movie}/{hole if qname_hole is None else qname_hole}/ccs"
    rg = "00000001" if movie == "movieA" else "00000002"
    fields = [qname, str(flag), "chr1", str(position), "60", "2M1I2M", "=", "30", "12", "ACGTA", "IIIII"]
    if flag & 4:
        fields[2:9] = ["*", "0", "0", "*", "*", "0", "0"]
    tags = ([f"zm:i:{hole}"] if zm else []) + [f"RG:Z:{rg}", "rq:f:0.998", "np:i:7",
             "bc:B:s,1,2", "ZZ:Z:retain_all_tags", "HH:H:0AFFFF", "aa:A:T", "xx:B:f,0.5,1.25"]
    return "\t".join(fields + tags)


def compare(args, directory, bam, header, ids, label, chunks, max_writers=16):
    ids_file = directory / (label + ".ids.txt")
    ids_file.write_text("".join(str(item) + "\n" for item in ids))
    prefix = directory / label
    result = run([args.dispatcher, "--input", bam, "--ids", ids_file,
                  "--chunks", chunks, "--output-prefix", prefix, "--threads", 2,
                  "--max-open-writers", max_writers])
    quotient, remainder = divmod(len(ids), chunks)
    offset, counts = 0, []
    for chunk in range(1, chunks + 1):
        size = quotient + (chunk <= remainder)
        selected = ids[offset:offset + size]
        offset += size
        include = directory / f"{label}.chunk{chunk}.include.txt"
        include.write_text("".join(str(item) + "\n" for item in selected))
        legacy = directory / f"{label}.chunk{chunk}.legacy.bam"
        candidate = directory / f"{label}.chunk{chunk}.bam"
        run([args.zmwfilter, "--include", include, bam, legacy])
        expected = run([args.samtools, "view", "--no-PG", legacy]).stdout
        actual = run([args.samtools, "view", "--no-PG", candidate]).stdout
        if actual != expected:
            raise AssertionError(f"{label} chunk {chunk}: SAM records/tags/order differ")
        actual_header = run([args.samtools, "view", "--no-PG", "-H", candidate]).stdout.splitlines()
        if any(line not in actual_header for line in header):
            raise AssertionError(f"{label} chunk {chunk}: original header line missing")
        run([args.samtools, "quickcheck", candidate])
        run([args.samtools, "index", candidate])
        run([args.pbindex, candidate])
        counts.append(len(actual.splitlines()))
    return {"chunks": chunks, "records_per_chunk": counts, "dispatcher_stdout": result.stdout}


def has_long_cigar_transport(bam, operation_count):
    """Confirm this fixture exercises BAM's CG:B:I transport representation."""
    data = gzip.decompress(bam.read_bytes())
    if data[:4] != b"BAM\x01":
        raise AssertionError("Fixture is not BAM")
    offset = 8 + struct.unpack_from("<i", data, 4)[0]
    references = struct.unpack_from("<i", data, offset)[0]
    offset += 4
    for _ in range(references):
        name_length = struct.unpack_from("<i", data, offset)[0]
        offset += 4 + name_length + 4
    while offset < len(data):
        block_length = struct.unpack_from("<i", data, offset)[0]
        record_bytes = data[offset + 4:offset + 4 + block_length]
        name_length = struct.unpack_from("<I", record_bytes, 8)[0] & 255
        cigar_count = struct.unpack_from("<I", record_bytes, 12)[0] & 65535
        sequence_length = struct.unpack_from("<i", record_bytes, 16)[0]
        aux_offset = 32 + name_length + cigar_count * 4 + (sequence_length + 1) // 2 + sequence_length
        aux = record_bytes[aux_offset:]
        tag = aux.find(b"CGBI")
        if cigar_count == 2 and tag >= 0 and struct.unpack_from("<I", aux, tag + 4)[0] == operation_count:
            return True
        offset += 4 + block_length
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dispatcher", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--zmwfilter", default="zmwfilter")
    parser.add_argument("--pbindex", default="pbindex")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    args.dispatcher = str(Path(args.dispatcher).resolve())
    directory = Path(args.output_dir).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    bam, header = fixture(args.samtools, directory, "adversarial", [
        record("movieA", 11, 1), record("movieB", 11, 1),
        record("movieA", 7, 5), record("movieA", 11, 10, 256),
        record("movieA", 13, 15, 2048), record("movieB", 7, 15, 16),
        record("movieA", 13, 20), record("movieB", 22, 0, 4)])
    run([args.pbindex, bam])
    observed = run([args.zmwfilter, "--show-all", bam]).stdout
    (directory / "legacy-show-all.txt").write_text(observed)
    ids = [int(item) for item in observed.splitlines()]
    results = {"legacy_enumeration": ids}
    results["enumerated"] = compare(args, directory, bam, header, ids, "enumerated", min(3, len(ids)))
    results["duplicate_ids_cross_chunks"] = compare(args, directory, bam, header,
                                                     [11, 7, 11, 13, 22], "duplicate_ids", 3)
    results["one_chunk"] = compare(args, directory, bam, header, ids, "one_chunk", 1)
    results["one_id_per_chunk"] = compare(args, directory, bam, header, ids, "one_id_per_chunk", len(ids))
    results["writer_groups_one"] = compare(args, directory, bam, header,
                                            [11, 7, 11, 13, 22], "writer_groups_one", 3, max_writers=1)
    results["writer_groups_two"] = compare(args, directory, bam, header,
                                            [11, 7, 11, 13, 22], "writer_groups_two", 3, max_writers=2)
    # More than 65,535 CIGAR operations are represented using CG:B:I in BAM.
    # Compare interpreted records/tags and order, never compressed BAM bytes.
    long_fields = record("movieA", 77, 10).split("\t")
    long_fields[5] = "1M1I" * 32768
    long_fields[9] = "AC" * 32768
    long_fields[10] = "I" * 65536
    long_bam, long_header = fixture(args.samtools, directory, "long_cigar", [
        record("movieA", 1, 1), "\t".join(long_fields), record("movieB", 2, 50000)], reference_length=100000)
    if not has_long_cigar_transport(long_bam, 65536):
        raise AssertionError("Long-CIGAR fixture did not produce CG:B:I transport encoding")
    run([args.pbindex, long_bam])
    long_ids = [int(item) for item in run([args.zmwfilter, "--show-all", long_bam]).stdout.splitlines()]
    results["long_cigar"] = compare(args, directory, long_bam, long_header, long_ids, "long_cigar_compare", 1)
    original_records = run([args.samtools, "view", "--no-PG", long_bam]).stdout
    dispatched_records = run([args.samtools, "view", "--no-PG", directory / "long_cigar_compare.chunk1.bam"]).stdout
    if original_records != dispatched_records or long_fields[5] not in dispatched_records:
        raise AssertionError("Long CIGAR, source tags, or ordered records changed")
    results["long_cigar"]["transport_operations"] = 65536
    compatibility_failures = []
    # Deliberately probe records that distinguish tag and QNAME semantics.
    for label, item in [("missing_zm", record("movieA", 11, 1, zm=False)),
                        ("tag_qname_disagree", record("movieA", 11, 1, qname_hole=44))]:
        probe_bam, probe_header = fixture(args.samtools, directory, label, [item])
        run([args.pbindex, probe_bam])
        legacy = subprocess.run([args.zmwfilter, "--show-all", str(probe_bam)],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        probe = {"show_all_status": legacy.returncode, "show_all_stdout": legacy.stdout,
                 "show_all_stderr": legacy.stderr}
        if legacy.returncode == 0:
            try:
                probe["comparison"] = compare(args, directory, probe_bam, probe_header,
                                               [int(line) for line in legacy.stdout.splitlines()], label + "_compare", 1)
            except subprocess.CalledProcessError as error:
                probe["candidate_failure"] = {"status": error.returncode, "stderr": error.stderr}
                compatibility_failures.append(label)
            except AssertionError as error:
                probe["candidate_mismatch"] = str(error)
                compatibility_failures.append(label)
        results[label] = probe
    # An invalid writer limit must fail before creating any output stream.
    bounded = subprocess.run([args.dispatcher, "--input", str(bam), "--ids", str(directory / "enumerated.ids.txt"),
                              "--chunks", "3", "--max-open-writers", "0", "--output-prefix", str(directory / "bounded")],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if bounded.returncode == 0 or list(directory.glob("bounded.chunk*.bam")):
        raise AssertionError("invalid writer bound was not rejected before output creation")
    results["invalid_writer_bound"] = "rejected before output creation"
    truncated = directory / "missing_eof.bam"
    truncated.write_bytes(bam.read_bytes()[:-28])
    rejected = subprocess.run([args.dispatcher, "--input", str(truncated), "--ids", str(directory / "enumerated.ids.txt"),
                               "--chunks", "1", "--output-prefix", str(directory / "truncated")],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if rejected.returncode == 0 or list(directory.glob("truncated.chunk*.bam")):
        raise AssertionError("truncated input was not rejected before output creation")
    results["missing_eof"] = "rejected before output creation"
    (directory / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))
    if compatibility_failures:
        raise AssertionError("Legacy-compatible probe failed: " + ", ".join(compatibility_failures))


if __name__ == "__main__":
    main()
