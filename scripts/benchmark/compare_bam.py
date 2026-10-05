#!/usr/bin/env python3
"""Validate ordered BAM contents and optionally rebuild indices in scratch only."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import struct
import zlib

MAX_FRAME_BYTES=64*1024*1024


class UnsupportedFrame(Exception):
    pass


def payload_digest(path, compressed=False):
    digest = hashlib.sha256()
    with (gzip.open(path,"rb") if compressed else path.open("rb")) as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_pg_command(line):
    if not line.startswith(b"@PG\t"):
        return line
    fields = line.rstrip(b"\r\n").split(b"\t")
    commands = [field for field in fields[1:] if field.startswith(b"CL:")]
    if len(commands)>1:
        raise ValueError("Duplicate CL fields in @PG header require review")
    return b"\t".join(field for field in fields if not field.startswith(b"CL:")) + b"\n"


def header_and_records(stream):
    header=[]
    for line in stream:
        if line.startswith((b"@HD\t",b"@SQ\t",b"@RG\t",b"@PG\t",b"@CO\t")):
            header.append(line)
        else:
            return header,itertools.chain((line,),stream)
    return header,iter(())


def compare_headers(left_header,right_header):
    try:
        left_normal=[normalize_pg_command(line) for line in left_header]
        right_normal=[normalize_pg_command(line) for line in right_header]
        if left_normal==right_normal:
            header=dict(status="pass",reason="all header fields/order exact except explicit @PG CL metadata")
        elif ([line for line in left_header if not line.startswith(b"@PG\t")]==
              [line for line in right_header if not line.startswith(b"@PG\t")]):
            header=dict(status="review",reason="other @PG provenance/count/order differs; no automatic whitelist")
        else:
            header=dict(status="failure",reason="non-PG header fields or order differ")
    except ValueError as error:
        header=dict(status="review",reason=str(error))
    header.update(reference_lines=len(left_header),candidate_lines=len(right_header),
                  reference_sha256=hashlib.sha256(b"".join(left_header)).hexdigest(),
                  candidate_sha256=hashlib.sha256(b"".join(right_header)).hexdigest())
    return header


def compare_records(left_records,right_records,encoding):
    counts=[0,0]
    first=None
    for index,(left,right) in enumerate(itertools.zip_longest(left_records,right_records),1):
        counts[0]+=left is not None
        counts[1]+=right is not None
        if left!=right and first is None:
            first=dict(record=index,
                       reference_sha256=hashlib.sha256(left).hexdigest() if left is not None else None,
                       candidate_sha256=hashlib.sha256(right).hexdigest() if right is not None else None)
    records=dict(status="pass" if first is None else "failure",reference_records=counts[0],
                 candidate_records=counts[1],first_difference=first,
                 encoding=encoding,reason="ordered record bytes exact" if first is None else "ordered record bytes differ")
    return records


def compare_sam_streams(reference,candidate):
    left_header,left_records=header_and_records(reference)
    right_header,right_records=header_and_records(candidate)
    return dict(headers=compare_headers(left_header,right_header),
                records=compare_records(left_records,right_records,"SAM text"))


def read_exact(stream,size):
    if size>MAX_FRAME_BYTES:
        raise UnsupportedFrame(f"BAM frame exceeds bounded reader limit {MAX_FRAME_BYTES}: {size}")
    data=stream.read(size)
    if len(data)!=size: raise EOFError("Truncated decompressed BAM frame")
    return data


def bam_header(stream):
    if read_exact(stream,4)!=b"BAM\x01": raise ValueError("Not BAM version 1")
    text=read_exact(stream,struct.unpack("<I",read_exact(stream,4))[0])
    count=struct.unpack("<I",read_exact(stream,4))[0]
    if count>1000000: raise UnsupportedFrame("BAM reference dictionary exceeds bounded reader limit")
    dictionary=[]
    dictionary_bytes=0
    for _ in range(count):
        size=struct.unpack("<I",read_exact(stream,4))[0]
        dictionary_bytes+=size+8
        if dictionary_bytes>MAX_FRAME_BYTES: raise UnsupportedFrame("BAM reference dictionary is too large")
        name=read_exact(stream,size)
        if not name or name[-1:]!=b"\0": raise ValueError("Invalid BAM reference name framing")
        length=read_exact(stream,4)
        dictionary.append((name,length))
    return text.splitlines(keepends=True),dictionary


def bam_records(stream):
    while True:
        prefix=stream.read(4)
        if not prefix: return
        if len(prefix)!=4: raise EOFError("Truncated BAM record length")
        size=struct.unpack("<I",prefix)[0]
        if size<32: raise ValueError("BAM record shorter than fixed fields")
        yield read_exact(stream,size)


def compare_bam_streams(reference,candidate):
    left_header,left_dictionary=bam_header(reference)
    right_header,right_dictionary=bam_header(candidate)
    return dict(headers=compare_headers(left_header,right_header),
                reference_dictionary=dict(status="pass" if left_dictionary==right_dictionary else "failure",
                                          reason="binary reference names/lengths/order compared exactly"),
                records=compare_records(bam_records(reference),bam_records(candidate),"raw decompressed BAM blocks"))


def command_result(command):
    process=subprocess.run([str(x) for x in command],capture_output=True,text=True)
    return dict(command=[str(x) for x in command],returncode=process.returncode,
                stdout=process.stdout,stderr=process.stderr)


def rebuild_indices(bam,scratch,samtools,pbindex):
    """Regenerate indices for this same BAM; never compare offsets across BAMs."""
    scratch.mkdir()
    link=scratch/"input.bam"
    link.symlink_to(bam.resolve())
    results=[]
    for suffix,command,compressed in (
        (".bai",[samtools,"index",link,Path(str(link)+".bai")],False),
        (".pbi",[pbindex,link],True)):
        original=Path(str(bam)+suffix)
        generated=Path(str(link)+suffix)
        if not original.is_file():
            results.append(dict(index=suffix,status="failure",reason="published index missing",path=str(original)))
            continue
        before=original.stat()
        original_digest=payload_digest(original,compressed)
        execution=command_result(command)
        after=original.stat()
        if (before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=(after.st_size,after.st_mtime_ns,after.st_ctime_ns):
            raise RuntimeError("Input index changed during validation: "+str(original))
        if execution["returncode"] or not generated.is_file():
            results.append(dict(index=suffix,status="failure",reason="scratch index rebuild failed",execution=execution))
        elif execution["stderr"]:
            results.append(dict(index=suffix,status="review",reason="index builder diagnostics require review",execution=execution))
        else:
            digest=payload_digest(generated,compressed)
            same=original_digest==digest
            results.append(dict(index=suffix,status="pass" if same else "review",
                                reason="index payload matches rebuild against its own BAM" if same else
                                       "index differs from rebuild; may be invalid or use a different valid encoding",
                                original_payload_sha256=original_digest,rebuilt_payload_sha256=digest,
                                rebuilt=str(generated),execution=execution))
    return results


def compare_bams(reference,candidate,output,samtools="samtools",pbindex="pbindex",rebuild=False):
    output.mkdir(parents=True)
    reference,candidate=reference.resolve(),candidate.resolve()
    before=[path.stat() for path in (reference,candidate)]
    report=dict(reference=str(reference),candidate=str(candidate),ignore_metadata=["@PG:CL"],
                index_rebuild_requested=rebuild,tools=dict(samtools=samtools,pbindex=pbindex))
    quick=[command_result([samtools,"quickcheck","-v",path]) for path in (reference,candidate)]
    report["quickcheck"]=quick
    if any(result["returncode"] for result in quick):
        return dict(report,status="failure",reason="BAM quickcheck failed")
    try:
        # BGZF is concatenated gzip. CRC/length errors are checked while reading;
        # record bytes are never reserialized through SAM or another BAM writer.
        with gzip.open(reference,"rb") as left,gzip.open(candidate,"rb") as right:
            report.update(compare_bam_streams(left,right))
    except UnsupportedFrame as error:
        return dict(report,status="review",reason=str(error))
    except (EOFError,ValueError,OSError,zlib.error) as error:
        return dict(report,status="failure",reason=f"BAM framing/decompression failed: {error}")
    statuses=[report["headers"]["status"],report["records"]["status"],report["reference_dictionary"]["status"]]
    if any(result["stderr"] for result in quick): statuses.append("review")
    if rebuild and "failure" not in statuses:
        report["indices"]={name:rebuild_indices(path,output/name,samtools,pbindex)
                           for name,path in (("reference",reference),("candidate",candidate))}
        statuses.extend(item["status"] for checks in report["indices"].values() for item in checks)
    else:
        report["indices"]={"status":"review","reason":"index validation not performed"}
        statuses.append("review")
    after=[path.stat() for path in (reference,candidate)]
    if any((a.st_size,a.st_mtime_ns,a.st_ctime_ns)!=(b.st_size,b.st_mtime_ns,b.st_ctime_ns)
           for a,b in zip(before,after)):
        statuses.append("failure")
        report["input_changed_during_validation"]=True
    return dict(report,status="failure" if "failure" in statuses else "review" if "review" in statuses else "pass")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference",type=Path)
    parser.add_argument("candidate",type=Path)
    parser.add_argument("--output",type=Path,required=True,help="New scratch/report directory")
    parser.add_argument("--samtools",default="samtools")
    parser.add_argument("--pbindex",default="pbindex")
    parser.add_argument("--rebuild-indexes",action="store_true",help="Extra BAM reads, validation cost only")
    args=parser.parse_args()
    if args.output.exists(): parser.error("output directory already exists")
    try:
        report=compare_bams(args.reference,args.candidate,args.output,args.samtools,args.pbindex,args.rebuild_indexes)
    except Exception as error:
        report=dict(status="review",reason=f"{type(error).__name__}: {error}")
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"report.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report))
    return {"pass":0,"failure":1,"review":2}[report["status"]]


if __name__=="__main__":
    raise SystemExit(main())
