import importlib.util
from pathlib import Path
import tempfile
import unittest
import zlib

MODULE = Path(__file__).resolve().parents[1] / "scripts/benchmark/compare_pdf.py"
spec = importlib.util.spec_from_file_location("compare_pdf", MODULE)
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)
AVAILABLE = importlib.util.find_spec("pypdf") is not None


def fixture(date="20261004000000", drawing=b"BT /F1 12 Tf 10 10 Td (Science) Tj ET",
            font=b"fixture-font-program", width="100.000000000000000001", level=6,
            producer="fixture", alias=False, unsupported=False):
    def stream(data):
        encoded = zlib.compress(data, level)
        filter_name = b"/Unknown" if unsupported else b"/FlateDecode"
        return b"<< /Length " + str(len(encoded)).encode() + b" /Filter " + filter_name + b" >>\nstream\n" + encoded + b"\nendstream"
    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width} 100] /Resources << /Font << /F1 6 0 R /F2 {'9' if alias else '6'} 0 R >> >> /Contents 4 0 R >>".encode(),
        stream(drawing),
        f"<< /CreationDate (D:{date}) /ModDate (D:{date}) /Producer ({producer}) >>".encode(),
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /FontDescriptor 7 0 R >>",
        b"<< /Type /FontDescriptor /FontName /Helvetica /FontFile 8 0 R >>",
        stream(font),
    ]
    if alias:
        objects.append(objects[5])
    data = b"%PDF-1.7\n"
    offsets = [0]
    for index, obj in enumerate(objects, 1):
        offsets.append(len(data))
        data += f"{index} 0 obj\n".encode() + obj + b"\nendobj\n"
    startxref = len(data)
    data += f"xref\n0 {len(offsets)}\n0000000000 65535 f \n".encode()
    for offset in offsets[1:]:
        data += f"{offset:010} 00000 n \n".encode()
    data += f"trailer\n<< /Size {len(offsets)} /Root 1 0 R /Info 5 0 R >>\nstartxref\n{startxref}\n%%EOF\n".encode()
    return data


@unittest.skipUnless(AVAILABLE, "Validation-only pinned pypdf is not on PYTHONPATH")
class PdfComparisonTest(unittest.TestCase):
    def check(self, left, right, expected, **options):
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / "a.pdf", Path(directory) / "b.pdf"
            a.write_bytes(left); b.write_bytes(right)
            result = compare.compare_pdf(a, b, **options)
            self.assertEqual(result["status"], expected, result)
            return result

    def test_dates_and_compression_are_only_default_exceptions(self):
        result = self.check(fixture(), fixture(date="20261005000000", level=1), "pass")
        self.assertEqual(result["ignored_candidate"], ["/Info/CreationDate", "/Info/ModDate"])
        self.check(fixture(), fixture(producer="different"), "failure")
        self.check(fixture(), fixture(producer="different"), "pass",
                   ignore_info=("CreationDate", "ModDate", "Producer"))

    def test_drawing_text_and_font_bytes_remain_exact(self):
        self.check(fixture(), fixture(drawing=b"BT (Different science) Tj ET"), "failure")
        self.check(fixture(), fixture(font=b"changed-font-program"), "failure")
        self.check(fixture(drawing=b"(/CreationDate old) Tj"),
                   fixture(drawing=b"(/CreationDate new) Tj"), "failure")

    def test_decimal_differences_below_float_precision_are_not_lost(self):
        self.check(fixture(), fixture(width="100.000000000000000002"), "failure")

    def test_indirect_alias_structure_is_preserved(self):
        self.check(fixture(), fixture(alias=True), "failure")
        result = self.check(fixture(), fixture().replace(b"%PDF-1.7",b"%PDF-1.3"), "failure")
        self.assertEqual(result["first_difference"]["canonical_path"], "/header")

    def test_unsupported_and_malformed_documents_require_review(self):
        self.check(fixture(), fixture(unsupported=True), "review")
        self.check(fixture(), b"not a PDF", "review")


if __name__ == "__main__":
    unittest.main()
