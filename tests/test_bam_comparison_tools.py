#!/usr/bin/env python3
"""Tiny pinned-tool integration fixtures; run on compute, never production BAMs."""
import argparse
import hashlib
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import struct


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--samtools",default="samtools")
    parser.add_argument("--pbindex",default="pbindex")
    parser.add_argument("--bgzip",default="bgzip")
    args=parser.parse_args()
    args.output.mkdir(parents=True)
    spec=importlib.util.spec_from_file_location("compare_bam",Path(__file__).resolve().parents[1]/"scripts/benchmark/compare_bam.py")
    compare=importlib.util.module_from_spec(spec);spec.loader.exec_module(compare)
    def run(command): subprocess.run([str(x) for x in command],check=True,capture_output=True)
    def fixture(name,command="old /path",version="1",zm=2,reverse=False):
        header=["@HD\tVN:1.6\tSO:coordinate","@SQ\tSN:chr1\tLN:1000",
                "@RG\tID:00000001\tPL:PACBIO\tDS:READTYPE=CCS\tPU:movie\tPM:SEQUEL",
                f"@PG\tID:align\tPN:align\tVN:{version}\tCL:{command}"]
        records=["movie/1/ccs\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:00000001\tzm:i:1\trq:f:0.99\tnp:i:5",
                 f"movie/{zm}/ccs\t0\tchr1\t10\t60\t4M\t*\t0\t0\tTGCA\tIIII\tRG:Z:00000001\tzm:i:{zm}\trq:f:0.98\tnp:i:6"]
        if reverse: records.reverse()
        sam=args.output/(name+".sam");bam=args.output/(name+".bam")
        sam.write_text("\n".join(header+records)+"\n")
        run([args.samtools,"view","--no-PG","-b","-o",bam,sam])
        run([args.samtools,"index",bam]);run([args.pbindex,bam])
        return bam
    reference=fixture("reference")
    candidate=fixture("candidate",command="new /different/path")
    paths=[Path(str(path)+suffix) for path in (reference,candidate) for suffix in ("", ".bai", ".pbi")]
    before={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    results={}
    results["metadata_only"]=compare.compare_bams(reference,candidate,args.output/"metadata-only",args.samtools,args.pbindex,True)
    assert results["metadata_only"]["status"]=="pass",results["metadata_only"]
    assert before=={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},"Scratch rebuild changed source files"
    results["unchecked_indices"]=compare.compare_bams(reference,candidate,args.output/"unchecked",args.samtools,args.pbindex,False)
    assert results["unchecked_indices"]["status"]=="review"
    provenance=fixture("provenance",version="2")
    results["other_provenance"]=compare.compare_bams(reference,provenance,args.output/"provenance-check",args.samtools,args.pbindex,True)
    assert results["other_provenance"]["headers"]["status"]=="review"
    reordered=fixture("reordered",reverse=True)
    results["order"]=compare.compare_bams(reference,reordered,args.output/"order-check",args.samtools,args.pbindex,True)
    assert results["order"]["records"]["status"]=="failure"
    different=fixture("different",zm=3)
    results["tags"]=compare.compare_bams(reference,different,args.output/"tags-check",args.samtools,args.pbindex,True)
    assert results["tags"]["records"]["status"]=="failure"
    raw=gzip.decompress(candidate.read_bytes())
    value=struct.pack("<f",0.99)
    marker=b"rqf"+value
    assert raw.count(marker)==1
    next_value=struct.pack("<I",struct.unpack("<I",value)[0]+1)
    mutated=args.output/"float-ulp.bam"
    with mutated.open("wb") as handle:
        subprocess.run([args.bgzip,"-c"],input=raw.replace(marker,b"rqf"+next_value,1),stdout=handle,check=True)
    results["float_ulp"]=compare.compare_bams(reference,mutated,args.output/"float-check",args.samtools,args.pbindex,False)
    assert results["float_ulp"]["records"]["status"]=="failure"
    sam_left=subprocess.check_output([args.samtools,"view","--no-PG",str(reference)])
    sam_right=subprocess.check_output([args.samtools,"view","--no-PG",str(mutated)])
    results["float_ulp"]["sam_text_identical_despite_ulp"]=sam_left==sam_right
    original=Path(str(candidate)+".pbi").read_bytes()
    Path(str(candidate)+".pbi").write_bytes(Path(str(different)+".pbi").read_bytes())
    results["stale_pbi"]=compare.compare_bams(reference,candidate,args.output/"stale-check",args.samtools,args.pbindex,True)
    assert results["stale_pbi"]["status"]=="review"
    assert results["stale_pbi"]["indices"]["candidate"][1]["status"]=="review"
    Path(str(candidate)+".pbi").write_bytes(original)
    truncated=args.output/"truncated.bam"
    truncated.write_bytes(candidate.read_bytes()[:-40])
    results["truncated"]=compare.compare_bams(reference,truncated,args.output/"truncated-check",args.samtools,args.pbindex,False)
    assert results["truncated"]["status"]=="failure"
    (args.output/"results.json").write_text(json.dumps(results,indent=2)+"\n")
    print("PASS: exact raw records/header/dictionary, one-ULP floats, scratch-only BAI/PBI rebuilds, source immutability, ordering, tags, stale PBI and truncation")


if __name__=="__main__": main()
