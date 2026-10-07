#!/usr/bin/env python3
"""Dispatch binary BAM records to the legacy ordered ZMW ID partitions.

Requires pysam from the pinned PacBio environment. Input decompression uses
--threads; each output compresses synchronously so writer count does not multiply
the thread budget. More chunks than --max-open-writers require additional reads.
"""
import argparse
from contextlib import ExitStack
import os
import re
import sys

import pysam


def positive_integer(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def assignments(ids_path, chunks):
    """Partition lines before deduplicating: duplicate IDs can span chunks."""
    ids = []
    with open(ids_path) as source:
        for line in source:
            text = line.strip(" \t\r\n")
            if not re.fullmatch(r"[+-]?[0-9]+", text):
                raise ValueError("invalid ZMW ID: " + repr(text))
            value = int(text)
            if not 0 <= value <= 2147483647:
                raise ValueError("ZMW ID outside nonnegative int32 range")
            ids.append(value)
    if len(ids) < chunks:
        raise ValueError("ID count is smaller than chunk count; legacy chunks would be empty")
    quotient, remainder = divmod(len(ids), chunks)
    destinations = {}
    offset = 0
    for chunk in range(chunks):
        size = quotient + (chunk < remainder)
        for zmw in ids[offset:offset + size]:
            matches = destinations.setdefault(zmw, [])
            if not matches or matches[-1] != chunk:
                matches.append(chunk)
        offset += size
    return destinations


def hole_number(record):
    # pbindex/zmwfilter index an absent zm as zero, even when QNAME disagrees.
    try:
        value, kind = record.get_tag("zm", with_value_type=True)
    except KeyError:
        return 0
    if kind not in "cCsSiI":
        raise ValueError("zm tag must have an integer BAM type")
    if not 0 <= value <= 2147483647:
        raise ValueError("invalid numeric zm tag")
    return value


def split_group(args, destinations, paths, first, last):
    with ExitStack() as streams:
        reader = streams.enter_context(pysam.AlignmentFile(
            args.input, "rb", threads=args.threads, check_sq=False))
        if not reader.is_bam:
            raise ValueError("input must be BAM")
        reader.check_truncation()
        writers = [streams.enter_context(pysam.AlignmentFile(
            path, "wb", template=reader)) for path in paths[first:last]]
        counts = [0] * len(writers)
        for record in reader:
            matches = destinations.get(hole_number(record))
            if matches is None:
                raise ValueError("BAM ZMW is missing from ID enumeration: " + record.query_name)
            for chunk in matches:
                if first <= chunk < last:
                    writers[chunk - first].write(record)
                    counts[chunk - first] += 1
    # Closing through ExitStack checks final writes before reporting success.
    if not all(counts):
        raise ValueError("output chunk contains no BAM records")
    for chunk, count in enumerate(counts, first):
        print("{}\t{}\t{}".format(chunk + 1, count, paths[chunk]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--ids", required=True)
    parser.add_argument("--chunks", required=True, type=positive_integer)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--threads", type=positive_integer, default=1,
                        help="input decompression threads (default: 1)")
    parser.add_argument("--max-open-writers", type=positive_integer, default=128)
    args = parser.parse_args()
    if args.threads > 4096 or args.max_open_writers > 4096:
        parser.error("--threads and --max-open-writers must not exceed 4096")
    try:
        destinations = assignments(args.ids, args.chunks)
        paths = ["{}.chunk{}.bam".format(args.output_prefix, chunk + 1)
                 for chunk in range(args.chunks)]
        for path in paths:
            if os.path.lexists(path):
                raise ValueError("refusing to overwrite existing output " + path)
        for first in range(0, args.chunks, args.max_open_writers):
            split_group(args, destinations, paths, first,
                        min(args.chunks, first + args.max_open_writers))
    except (OSError, ValueError) as error:
        print("splitBamByZmw: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
