#!/usr/bin/env python3
"""Freeze and profile the unchanged germline VCF export with final QS resident."""
import argparse
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    stat = path.stat()
    # Reads may update atime; retain metadata that detects content/inode changes.
    return dict(sha256=digest.hexdigest(), device=stat.st_dev, inode=stat.st_ino,
                mode=stat.st_mode, size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                ctime_ns=stat.st_ctime_ns)


def validate_own_index(path):
    # Rebuild against the same BGZF encoding, using the pinned index writer.
    # Comparing indices from different BGZF encodings would be invalid.
    with tempfile.TemporaryDirectory(prefix="vcf-index-") as folder:
        staged = Path(folder) / "input.vcf.bgz"
        # Rsamtools normalizes paths; symlinks would resolve back to inputs.
        shutil.copyfile(path, staged)
        subprocess.run(["Rscript", "--vanilla", "-e",
                        'Rsamtools::indexTabix(commandArgs(TRUE)[[1]], format="vcf")',
                        str(staged)], check=True, stdout=subprocess.DEVNULL)
        with gzip.open(str(staged) + ".tbi", "rb") as rebuilt, gzip.open(str(path) + ".tbi", "rb") as original:
            if rebuilt.read() != original.read():
                raise AssertionError(f"Index differs from rebuild for its own BGZF encoding: {path}")


def compare_vcf(reference, candidate, own_indices=False):
    def lines(path):
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if not line.startswith("##fileDate="):
                    yield line
    records = 0
    boundaries = {}
    for number, (left, right) in enumerate(itertools.zip_longest(lines(reference), lines(candidate)), 1):
        if left != right:
            raise AssertionError(f"VCF mismatch at scientific line {number}")
        if not left.startswith("#"):
            chrom, pos = left.split("\t", 2)[:2]
            boundaries.setdefault(chrom, [pos, pos])[1] = pos
            records += 1
    tabix = shutil.which("tabix")
    if not tabix:
        raise RuntimeError("tabix required for index validation")
    def query(path, *args):
        return subprocess.check_output([tabix, *args, str(path)] if args == ("-l",)
                                       else [tabix, str(path), *args])
    if query(reference, "-l") != query(candidate, "-l"):
        raise AssertionError("Index sequence names differ")
    for chrom, positions in boundaries.items():
        for pos in set(positions):
            interval = f"{chrom}:{pos}-{pos}"
            if query(reference, interval) != query(candidate, interval):
                raise AssertionError(f"Index query differs: {interval}")
    if own_indices:
        validate_own_index(reference)
        validate_own_index(candidate)
    return dict(scientific_text_identical=True, ignored_header_keys=["fileDate"],
                own_encoding_index_rebuild_identical=own_indices,
                records=records, indexed_chromosomes=len(boundaries),
                index_check="sequence names and first/last record positions on each chromosome")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--reference-vcf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chromgroup", default="1-22X")
    parser.add_argument("--filtergroup", default="lenient")
    parser.add_argument("--pairs", type=int, default=0,
                        help="0: unchanged diagnostic profile; positive: alternating fresh default/native buffer pairs")
    args = parser.parse_args()
    if args.pairs < 0:
        parser.error("pairs must be nonnegative")
    for path in (args.input, args.reference_vcf, Path(str(args.reference_vcf) + ".tbi")):
        if not path.is_file():
            parser.error(f"Missing input: {path}")
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output.resolve()
    code = output / "code"
    code.mkdir()
    sources = [args.repo / "bin" / name for name in ("outputResults.R", "sharedFunctions.R")]
    sources += [Path(__file__).parent / name for name in (
        "profile_germline_vcf.R", "profile_germline_vcf.py", "germline_formatting_helpers.R", "measure.py")]
    if args.pairs:
        sources.append(args.repo / "tests" / "test_vcf_native_buffer.R")
    hashes = {}
    for source in sources:
        target = code / source.name
        shutil.copyfile(source, target)
        hashes[source.name] = hashlib.sha256(target.read_bytes()).hexdigest()
    manifest = dict(input=str(args.input.resolve()), reference_vcf=str(args.reference_vcf.resolve()),
                    input_size=args.input.stat().st_size, source_sha256=hashes,
                    hostname=os.uname().nodename, slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                    chromgroup=args.chromgroup, filtergroup=args.filtergroup, pairs=args.pairs,
                    native_nchunk=100000 if args.pairs else None)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    worker = str(code / "profile_germline_vcf.R")
    preflight = ["Rscript", "--vanilla", worker, str(code), "--preflight",
                 args.chromgroup, args.filtergroup, str(output / "preflight")]
    if args.pairs:
        preflight.append("default")
    subprocess.run(preflight, check=True)
    if args.pairs:
        # Actual pinned BGZF/index fixtures establish scratch-copy immutability
        # before any large QS is loaded or any published file is compared.
        fixture_paths = sorted(output.glob("preflight.fixture.FALSE.*.vcf.bgz"))
        if len(fixture_paths) != 3:
            raise AssertionError("Expected three native-buffer BGZF fixtures")
        fixture_inputs = fixture_paths + [Path(str(path) + ".tbi") for path in fixture_paths]
        before = {str(path): fingerprint(path) for path in fixture_inputs}
        for path in fixture_paths:
            validate_own_index(path)
        after = {str(path): fingerprint(path) for path in fixture_inputs}
        if before != after:
            raise AssertionError("Index verification mutated an original BGZF/index fixture")
        (output / "preflight-input-immutability.json").write_text(json.dumps(
            dict(unchanged=True, before=before, after=after), indent=2) + "\n")
        print("Pinned index rebuild and original BGZF/index hash/stat immutability passed", flush=True)

    def run(prefix, mode):
        subprocess.run([sys.executable, str(code / "measure.py"), "--out", str(prefix) + ".metrics",
                        "--label", "germline-vcf-" + mode, "--allocated-cpus", "2", "--",
                        "Rscript", "--vanilla", worker, str(code), str(args.input.resolve()),
                        args.chromgroup, args.filtergroup, str(prefix), mode], check=True)

    if not args.pairs:
        prefix = output / "germline"
        run(prefix, "profile")
        result = compare_vcf(args.reference_vcf.resolve(), Path(str(prefix) + ".vcf.bgz"))
        (output / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
        return
    for pair in range(1, args.pairs + 1):
        folder = output / f"pair{pair:02d}"
        folder.mkdir()
        for mode in (("default", "100000") if pair % 2 else ("100000", "default")):
            print(f"Pair {pair}: fresh worker {mode}", flush=True)
            run(folder / mode, mode)
        default = folder / "default.vcf.bgz"
        candidate = folder / "100000.vcf.bgz"
        published = compare_vcf(args.reference_vcf.resolve(), default, own_indices=True)
        result = compare_vcf(default, candidate, own_indices=True)
        (folder / "comparison.json").write_text(json.dumps(
            dict(published_vs_default=published, default_vs_native=result), indent=2) + "\n")
        print(f"Pair {pair}: all scientific text and own-encoding index comparisons passed", flush=True)


if __name__ == "__main__":
    main()
