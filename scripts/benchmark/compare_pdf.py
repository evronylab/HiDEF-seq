#!/usr/bin/env python3
"""Conservative exact comparison of rooted PDF object graphs and decoded streams.

Requires validation-only pypdf==6.19.0. Unrecognized features remain review-required.
This is not a visual similarity test or a general PDF conformance validator.
"""
import argparse
from contextlib import contextmanager
from decimal import Decimal
import hashlib
import json
import logging
from pathlib import Path

PINNED_PYPDF = "6.19.0"
DEFAULT_INFO_IGNORES = ("CreationDate", "ModDate")


class NeedsReview(Exception):
    pass


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def first_difference(left, right, path=""):
    if type(left) is type(right) and isinstance(left, dict) and left.keys() == right.keys():
        for key in left:
            if left[key] != right[key]:
                return first_difference(left[key], right[key], path + "/" + key)
    elif type(left) is type(right) and isinstance(left, list) and len(left) == len(right):
        for index, (a, b) in enumerate(zip(left, right)):
            if a != b:
                return first_difference(a, b, path + "/" + str(index))
    return dict(canonical_path=path, reference=repr(left)[:500], candidate=repr(right)[:500])


@contextmanager
def exact_pdf_decimals():
    # The pinned reader constructs FloatObject from the original numeric token.
    # Retain that decimal token's value before binary-float rounding, so distinct
    # page/font numeric values cannot accidentally compare equal after parsing.
    import pypdf.generic._base as base
    original = base.FloatObject

    class ExactFloat(original):
        def __new__(cls, value="0.0", context=None):
            result = super().__new__(cls, value, context)
            result.exact_decimal = Decimal(value.decode("ascii") if isinstance(value, bytes) else str(value))
            if not result.exact_decimal.is_finite():
                raise NeedsReview("Non-finite PDF number")
            return result

    base.FloatObject = ExactFloat
    try:
        yield
    finally:
        base.FloatObject = original


def canonical_pdf(path, ignore_info=DEFAULT_INFO_IGNORES, ignore_document_id=False):
    import pypdf
    from pypdf import PdfReader
    from pypdf.generic import (ArrayObject, BooleanObject, ByteStringObject, DictionaryObject,
                               FloatObject, IndirectObject, NameObject, NullObject,
                               NumberObject, StreamObject, TextStringObject)
    if pypdf.__version__ != PINNED_PYPDF:
        raise NeedsReview(f"Expected pypdf=={PINNED_PYPDF}; found {pypdf.__version__}")
    if any(not name or "/" in name for name in ignore_info):
        raise NeedsReview("Info whitelist entries must be exact top-level metadata key names")
    nodes, seen, ignored = [], {}, []
    warnings_seen = []

    class Capture(logging.Handler):
        def emit(self, record):
            warnings_seen.append(record.getMessage())

    handler = Capture(level=logging.WARNING)
    logger = logging.getLogger("pypdf")
    logger.addHandler(handler)
    try:
        with exact_pdf_decimals(), open(path, "rb") as source:
            reader = PdfReader(source, strict=True)
            if reader.is_encrypted or "/Encrypt" in reader.trailer:
                raise NeedsReview("Encrypted PDF is unsupported")
            if "/Prev" in reader.trailer:
                raise NeedsReview("Incremental PDF revisions require separate review")

            def walk(obj, location):
                if isinstance(obj, IndirectObject):
                    key = (obj.idnum, obj.generation)
                    if key not in seen:
                        if len(nodes) >= 100000:
                            raise NeedsReview("Object graph limit exceeded")
                        index = len(nodes)
                        seen[key] = index
                        nodes.append(None)
                        nodes[index] = walk(obj.get_object(), location)
                    return ["ref", seen[key]]
                if isinstance(obj, StreamObject):
                    if any(key in obj for key in ("/F", "/FFilter", "/FDecodeParms")):
                        raise NeedsReview("External stream data is unsupported")
                    filters = obj.get("/Filter", [])
                    if isinstance(filters, IndirectObject):
                        filters = filters.get_object()
                    if isinstance(filters, NameObject):
                        filters = [filters]
                    if any(str(item) not in ("/FlateDecode", "/Fl") for item in filters):
                        raise NeedsReview(f"Unsupported stream filter at {location}: {filters}")
                    decoded = obj.get_data()
                    attrs = {key: value for key, value in obj.items()
                             if str(key) not in ("/Length", "/Filter", "/DecodeParms")}
                    return ["stream", mapping(attrs, location), len(decoded), sha256(decoded)]
                if isinstance(obj, DictionaryObject):
                    return ["dict", mapping(obj, location)]
                if isinstance(obj, ArrayObject):
                    return ["array", [walk(item, location + f"/[{index}]") for index, item in enumerate(obj)]]
                if isinstance(obj, NullObject):
                    return ["null"]
                if isinstance(obj, BooleanObject):
                    return ["boolean", obj.value]
                if isinstance(obj, FloatObject):
                    if not hasattr(obj, "exact_decimal"):
                        raise NeedsReview("PDF float lacks its exact original decimal value")
                    # Decimal equality is exact; as_tuple avoids context-dependent normalize rounding.
                    value = obj.exact_decimal
                    parts = value.as_tuple()
                    digits, exponent = list(parts.digits), parts.exponent
                    while len(digits) > 1 and digits[-1] == 0:
                        digits.pop(); exponent += 1
                    return ["decimal", parts.sign, digits, exponent]
                if isinstance(obj, NumberObject):
                    return ["integer", int(obj)]
                if isinstance(obj, NameObject):
                    return ["name", str(obj)]
                if isinstance(obj, TextStringObject):
                    return ["text", str(obj)]
                if isinstance(obj, ByteStringObject):
                    return ["bytes", bytes(obj).hex()]
                raise NeedsReview(f"Unsupported PDF object type at {location}: {type(obj).__name__}")

            def mapping(obj, location):
                result = []
                for key in sorted(obj, key=str):
                    name = str(key)
                    child = location + name
                    if name in ("/OpenAction", "/AA", "/JavaScript", "/JS", "/AcroForm", "/XFA"):
                        raise NeedsReview(f"Interactive/active PDF feature requires review: {child}")
                    if location == "/Info" and name[1:] in ignore_info:
                        value = obj[key]
                        if isinstance(value, IndirectObject):
                            value = value.get_object()
                        if not isinstance(value, (TextStringObject, ByteStringObject, NullObject)):
                            raise NeedsReview(f"Ignored Info metadata is not a scalar string: {child}")
                        ignored.append(child)
                        continue
                    result.append([name, walk(obj.raw_get(key) if hasattr(obj, "raw_get") else obj[key], child)])
                return result

            # These entries describe serialization/xref layout, not document content.
            trailer = {key: value for key, value in reader.trailer.items()
                       if str(key) not in ("/Size", "/XRefStm")}
            if ignore_document_id and "/ID" in trailer:
                ignored.append("/ID")
                del trailer["/ID"]
            root = mapping(trailer, "")
            result = {"header": reader.pdf_header, "trailer": root, "objects": nodes}
        if warnings_seen:
            raise NeedsReview("Parser warnings: " + "; ".join(warnings_seen))
        return result, ignored
    finally:
        logger.removeHandler(handler)


def compare_pdf(reference, candidate, ignore_info=DEFAULT_INFO_IGNORES, ignore_document_id=False):
    left_bytes, right_bytes = Path(reference).read_bytes(), Path(candidate).read_bytes()
    result = dict(reference=str(reference), candidate=str(candidate),
                  reference_sha256=sha256(left_bytes), candidate_sha256=sha256(right_bytes),
                  pypdf_version=PINNED_PYPDF, ignored_info_keys=list(ignore_info),
                  ignore_document_id=ignore_document_id)
    if left_bytes == right_bytes:
        return dict(result, status="pass", reason="byte identical")
    try:
        left, ignored_left = canonical_pdf(reference, ignore_info, ignore_document_id)
        right, ignored_right = canonical_pdf(candidate, ignore_info, ignore_document_id)
    except Exception as error:
        return dict(result, status="review", reason=f"{type(error).__name__}: {error}")
    if left == right:
        return dict(result, status="pass", reason="rooted PDF graph and all decoded streams identical",
                    ignored_reference=ignored_left, ignored_candidate=ignored_right,
                    indirect_objects=len(left["objects"]))
    return dict(result, status="failure", reason="parsed document graph or decoded content differs",
                ignored_reference=ignored_left, ignored_candidate=ignored_right,
                first_difference=first_difference(left,right))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--ignore-info-key", action="append", default=list(DEFAULT_INFO_IGNORES),
                        help="Explicit additional Info metadata key only; never replaces content text")
    parser.add_argument("--ignore-document-id", action="store_true")
    args = parser.parse_args()
    if args.report.exists():
        parser.error("report already exists")
    result = compare_pdf(args.reference, args.candidate, args.ignore_info_key, args.ignore_document_id)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
    return {"pass": 0, "failure": 1, "review": 2}[result["status"]]


if __name__ == "__main__":
    raise SystemExit(main())
