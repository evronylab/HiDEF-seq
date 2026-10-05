import importlib.util
import io
import struct
from pathlib import Path
import unittest

SPEC=importlib.util.spec_from_file_location("compare_bam",Path(__file__).resolve().parents[1]/"scripts/benchmark/compare_bam.py")
MODULE=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class BamComparisonTests(unittest.TestCase):
    header=b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n@PG\tID:align\tPN:align\tVN:1\tCL:old /path\n"
    records=b"read1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\tzm:i:1\nread2\t0\tchr1\t1\t60\t1M\t*\t0\t0\tC\tI\tzm:i:2\n"

    def compare(self,left,right):
        return MODULE.compare_sam_streams(io.BytesIO(left),io.BytesIO(right))

    def test_only_pg_command_is_normalized(self):
        result=self.compare(self.header+self.records,self.header.replace(b"old /path",b"new /other")+self.records)
        self.assertEqual(result["headers"]["status"],"pass")
        self.assertEqual(result["records"]["status"],"pass")
        result=self.compare(self.header+self.records,self.header.replace(b"VN:1\tCL",b"VN:2\tCL")+self.records)
        self.assertEqual(result["headers"]["status"],"review")
        result=self.compare(self.header+self.records,self.header.replace(b"LN:100",b"LN:101")+self.records)
        self.assertEqual(result["headers"]["status"],"failure")

    def test_pg_count_order_and_duplicate_fields_require_review(self):
        for header in (self.header+b"@PG\tID:extra\tPN:extra\n",self.header.replace(b"CL:old /path",b"CL:one\tCL:two")):
            result=self.compare(self.header+self.records,header+self.records)
            self.assertEqual(result["headers"]["status"],"review")
            self.assertEqual(result["records"]["status"],"pass")

    def test_record_order_tags_and_truncation_remain_exact(self):
        for records in (b"\n".join(reversed(self.records.rstrip(b"\n").split(b"\n")))+b"\n",
                        self.records.replace(b"zm:i:2",b"zm:i:3"),self.records.splitlines(keepends=True)[0]):
            result=self.compare(self.header+self.records,self.header+records)
            self.assertEqual(result["records"]["status"],"failure")

    def test_empty_and_non_pg_header_order(self):
        self.assertEqual(self.compare(self.header,self.header)["records"]["reference_records"],0)
        lines=self.header.splitlines(keepends=True)
        result=self.compare(self.header,b"".join([lines[1],lines[0],lines[2]]))
        self.assertEqual(result["headers"]["status"],"failure")

    def raw_bam(self,records,dictionary=((b"chr1\0",100),)):
        return (b"BAM\1"+struct.pack("<I",len(self.header))+self.header+
                struct.pack("<I",len(dictionary))+
                b"".join(struct.pack("<I",len(name))+name+struct.pack("<I",length)
                         for name,length in dictionary)+
                b"".join(struct.pack("<I",len(record))+record for record in records))

    def test_raw_records_keep_one_ulp_float_difference(self):
        value=struct.pack("<f",0.99)
        next_value=struct.pack("<I",struct.unpack("<I",value)[0]+1)
        left=self.raw_bam([b"\0"*32+b"rqf"+value])
        right=self.raw_bam([b"\0"*32+b"rqf"+next_value])
        result=MODULE.compare_bam_streams(io.BytesIO(left),io.BytesIO(right))
        self.assertEqual(result["records"]["status"],"failure")
        self.assertEqual(result["records"]["encoding"],"raw decompressed BAM blocks")

    def test_reference_dictionary_and_bounded_record_framing(self):
        left=self.raw_bam([])
        right=self.raw_bam([],dictionary=((b"chr1\0",101),))
        result=MODULE.compare_bam_streams(io.BytesIO(left),io.BytesIO(right))
        self.assertEqual(result["reference_dictionary"]["status"],"failure")
        with self.assertRaises(EOFError): list(MODULE.bam_records(io.BytesIO(struct.pack("<I",32)+b"short")))
        with self.assertRaises(MODULE.UnsupportedFrame):
            list(MODULE.bam_records(io.BytesIO(struct.pack("<I",MODULE.MAX_FRAME_BYTES+1))))


if __name__=="__main__": unittest.main()
