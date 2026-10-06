"""Small evidence fixtures; these tests never parse BAMs or launch index tools."""
from collections import Counter
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "accept_bam_coordinate_diagnoses",
    Path(__file__).resolve().parents[1] / "scripts/benchmark/accept_bam_coordinate_diagnoses.py",
)
GATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GATE)


class EvidenceFixture:
    def __init__(self, root):
        self.root = Path(root)
        self.plan_path = self.root / "plan.json"
        self.plan = {
            "approval": {"rule": GATE.RULE, "source": "explicit user message", "date": "2026-10-06"},
            "evidence_root": str(self.root), "evidence_sha256": {},
            "original_results": "original/results.json",
            "original_summary": "original/summary.json",
            "original_manifest": "original/manifest.json",
            "publication_source_audit": "original/publication-source-audit.json",
            "diagnoses": [], "input_identities": {},
        }
        self.results = [{"path": f"science/{i:03d}.tsv", "kind": "text", "status": "pass"}
                        for i in range(705)]
        self.manifest = {"reference": str(self.root / "baseline"),
                         "candidate": str(self.root / "candidate")}
        self.documents = {}
        for lib in (1, 2):
            path = f"sample-LIB{lib}/processedReads/sample-LIB{lib}.bam"
            directory = f"diagnosis-lib{lib}"
            self.plan["diagnoses"].append({"path": path, "directory": directory, "job_id": str(100 + lib)})
            original = {
                "path": path, "kind": "bam", "status": "failure",
                "headers": {"status": "pass", "reference_lines": 2, "candidate_lines": 2},
                "reference_dictionary": {"status": "pass"},
                "quickcheck": [{"returncode": 0, "stderr": ""}] * 2,
                "records": {"status": "failure", "reason": "ordered record bytes differ",
                            "reference_records": 100 + lib, "candidate_records": 100 + lib},
                "tools": {"samtools": "samtools", "pbindex": "/pinned/pbindex"},
            }
            indices = {}
            for side in ("reference", "candidate"):
                bam = Path(self.manifest[side]) / path
                original[side] = str(bam)
                for suffix in ("", ".bai", ".pbi"):
                    source = Path(str(bam) + suffix)
                    source.parent.mkdir(parents=True, exist_ok=True)
                    source.write_bytes(f"synthetic {lib} {side} {suffix}".encode())
                    self.plan["input_identities"][str(source)] = GATE.identity(source)
                link = self.root / directory / side / "input.bam"
                link.parent.mkdir(parents=True, exist_ok=True)
                link.symlink_to(bam)
                indices[side] = []
                for suffix in (".bai", ".pbi"):
                    rebuilt = Path(str(link) + suffix)
                    rebuilt.write_bytes(b"synthetic index evidence")
                    command = (["samtools", "index", str(link), str(rebuilt)] if suffix == ".bai"
                               else ["/pinned/pbindex", str(link)])
                    indices[side].append({
                        "index": suffix, "status": "pass", "rebuilt": str(rebuilt),
                        "original_payload_sha256": "a" * 64, "rebuilt_payload_sha256": "a" * 64,
                        "execution": {"command": command, "returncode": 0, "stderr": ""},
                    })
            source_hashes = {}
            for label, name in (("helper", "compare_bam.py"), ("diagnostic", "diagnose.py")):
                relative = directory + "/" + name
                source = self.root / relative
                source.write_text("# frozen synthetic " + label + "\n")
                self.plan["evidence_sha256"][relative] = GATE.sha(source)
                source_hashes[label] = GATE.sha(source)
            self.results.append(original)
            self.results.extend({"path": path + suffix, "kind": "bam_index", "status": "covered", "owner": path}
                                for suffix in (".bai", ".pbi"))
            self.documents[directory + "/failed-comparison.json"] = copy.deepcopy(original)
            self.documents[directory + "/record-diagnosis.json"] = {
                "job_id": str(100 + lib), "headers": copy.deepcopy(original["headers"]),
                "reference_dictionary_exact": True, "source_sha256": source_hashes,
                "totals": {"reference_records": 100 + lib, "candidate_records": 100 + lib,
                           "groups": 80, "groups_with_order_only_difference": 2,
                           "different_group_coordinates": 0, "different_record_multisets": 0},
            }
            self.documents[directory + "/index-diagnosis.json"] = indices
            self.documents[directory + "/completion.json"] = {
                "diagnostic_complete": True, "inputs_unchanged": True,
                "index_statuses": {"reference": ["pass", "pass"], "candidate": ["pass", "pass"]},
            }
        required = [{"path": x["path"], "kind": x["kind"]} for x in self.results]
        self.manifest["required_publications"] = {"required": required}
        self.manifest["inventory"] = copy.deepcopy(required)
        self.results.extend({"path": p, "status": "pass"}
                            for p in ("@source-config", "@publication-source-audit"))
        self.documents[self.plan["original_manifest"]] = self.manifest
        self.documents[self.plan["original_results"]] = self.results
        self.documents[self.plan["original_summary"]] = {
            "scope": "full-publication", "status": "failure", "required_publication_count": 711,
            "counts": dict(Counter(x["status"] for x in self.results)),
        }
        checks = []
        for row in self.results:
            if row.get("kind") in ("bam", "bam_index"):
                stat = (Path(self.manifest["candidate"]) / row["path"]).stat()
                checks.append({"path": row["path"], "status": "pass", "publication_stat": [
                    stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]})
        self.documents[self.plan["publication_source_audit"]] = {"status": "pass", "checks": checks}
        self.save()

    def save(self):
        for relative, document in self.documents.items():
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(document))
            self.plan["evidence_sha256"][relative] = GATE.sha(path)
        for item in self.plan["diagnoses"]:
            item["diagnostic_start_ns"] = min(
                (self.root / item["directory"] / name).stat().st_mtime_ns
                for name in ("compare_bam.py", "diagnose.py", "failed-comparison.json"))
        self.plan_path.write_text(json.dumps(self.plan))

    def doc(self, name, lib=1):
        return self.documents[f"diagnosis-lib{lib}/{name}.json"]

    def row(self, suffix="", lib=1):
        path = self.plan["diagnoses"][lib - 1]["path"] + suffix
        return next(x for x in self.results if x["path"] == path)


class BamCoordinateAcceptanceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.fixture = EvidenceFixture(self.temporary.name)

    def reject(self, message):
        self.fixture.save()
        with self.assertRaisesRegex(ValueError, message):
            GATE.accept(self.fixture.plan_path)

    def test_accepts_only_six_paths_and_preserves_original_failure_bytes(self):
        f = self.fixture
        before = {p: (f.root / p).read_bytes() for p in f.plan["evidence_sha256"]}
        result = GATE.accept(f.plan_path)
        self.assertEqual(result["status"], "pass")
        self.assertEqual(result["required_publications_accepted"], 711)
        self.assertEqual(result["non_bam_original_passes"], 705)
        self.assertEqual(Counter(x["kind"] for x in result["accepted"]), {"bam": 2, "bam_index": 4})
        self.assertEqual(result["original_report_status"], "failure")
        self.assertFalse(result["original_report_modified"])
        self.assertEqual(before, {p: (f.root / p).read_bytes() for p in before})

    def test_wrong_sample_diagnosis_cannot_replace_required_bam(self):
        self.fixture.plan["diagnoses"][0]["path"] = "wrong-sample.bam"
        self.reject("exactly the two required BAMs")

    def test_swapped_failed_sample_report_is_rejected_even_when_repinned(self):
        self.fixture.documents["diagnosis-lib1/failed-comparison.json"] = copy.deepcopy(
            self.fixture.doc("failed-comparison", 2))
        self.reject("different original BAM result")

    def test_index_input_symlink_cannot_point_at_other_sample(self):
        f = self.fixture
        link = f.root / "diagnosis-lib1/candidate/input.bam"
        link.unlink()
        link.symlink_to(f.row(lib=2)["candidate"])
        self.reject("scratch symlink targets a different BAM")

    def test_altered_pinned_source_and_report_are_rejected(self):
        f = self.fixture
        for relative in ("diagnosis-lib1/compare_bam.py", "diagnosis-lib1/record-diagnosis.json"):
            with self.subTest(relative=relative):
                path = f.root / relative
                original = path.read_bytes()
                path.write_bytes(original + b" ")
                with self.assertRaisesRegex(ValueError, "Evidence checksum mismatch"):
                    GATE.accept(f.plan_path)
                path.write_bytes(original)

    def test_repinned_source_must_match_recorded_source_hash(self):
        self.fixture.doc("record-diagnosis")["source_sha256"]["diagnostic"] = "b" * 64
        self.reject("source pin mismatch")

    def test_coordinate_and_record_multiset_differences_remain_blocking(self):
        totals = self.fixture.doc("record-diagnosis")["totals"]
        for key, message in (("different_group_coordinates", "Coordinates"),
                             ("different_record_multisets", "Record contents")):
            with self.subTest(key=key):
                totals[key] = 1
                self.reject(message)
                totals[key] = 0

    def test_count_mismatch_is_rejected(self):
        self.fixture.doc("record-diagnosis")["totals"]["candidate_records"] -= 1
        self.reject("record count mismatch")

    def test_missing_index_is_rejected(self):
        self.fixture.doc("index-diagnosis")["candidate"].pop()
        self.reject("Missing or duplicate BAI/PBI")

    def test_failed_index_is_rejected_even_with_matching_completion(self):
        self.fixture.doc("index-diagnosis")["candidate"][0]["status"] = "failure"
        self.fixture.doc("completion")["index_statuses"]["candidate"][0] = "failure"
        self.reject("Index rebuild did not pass")

    def test_index_stderr_is_not_silently_accepted(self):
        self.fixture.doc("index-diagnosis")["reference"][1]["execution"]["stderr"] = "warning"
        self.reject("Index rebuild did not pass")

    def test_wrong_index_command_input_is_rejected(self):
        self.fixture.doc("index-diagnosis")["candidate"][1]["execution"]["command"][-1] = "other.bam"
        self.reject("Index command used a different BAM")

    def test_index_payload_mismatch_is_rejected(self):
        self.fixture.doc("index-diagnosis")["reference"][0]["rebuilt_payload_sha256"] = "b" * 64
        self.reject("Index payload differs")

    def test_companion_owner_mismatch_is_rejected(self):
        self.fixture.row(".bai")["owner"] = self.fixture.row(lib=2)["path"]
        self.reject("companion ownership differs")

    def test_incomplete_or_changed_input_diagnostic_is_rejected(self):
        completion = self.fixture.doc("completion")
        for key in ("diagnostic_complete", "inputs_unchanged"):
            with self.subTest(key=key):
                completion[key] = False
                self.reject("Incomplete or changed-input")
                completion[key] = True

    def test_extra_scientific_failure_is_not_waived(self):
        f = self.fixture
        f.results[0]["status"] = "failure"
        f.documents[f.plan["original_summary"]]["counts"] = dict(Counter(x["status"] for x in f.results))
        self.reject("Unresolved failure/review outside")

    def test_current_bam_or_index_identity_change_is_rejected(self):
        for suffix in ("", ".bai", ".pbi"):
            with self.subTest(suffix=suffix):
                f = EvidenceFixture(self.fixture.root / ("changed" + (suffix or ".bam")))
                path = Path(f.row()["candidate"] + suffix)
                path.write_bytes(path.read_bytes() + b"changed")
                with self.assertRaisesRegex(ValueError, "Input identity changed"):
                    GATE.accept(f.plan_path)

    def test_missing_user_approval_is_rejected(self):
        self.fixture.plan["approval"]["rule"] = "Any record ordering is permitted"
        self.reject("Missing approved rule")

    def test_same_count_with_replaced_result_path_is_rejected(self):
        self.fixture.results[0]["path"] = "unreviewed/replacement.tsv"
        self.reject("[Ii]nventory|[Pp]ath")

    def test_candidate_provenance_stat_mismatch_is_rejected(self):
        audit = self.fixture.documents[self.fixture.plan["publication_source_audit"]]
        audit["checks"][0]["publication_stat"][1] += 1
        self.reject("[Pp]rovenance|[Ss]tat|[Ii]dentity")

    def test_repinned_input_change_after_startup_is_rejected(self):
        f = self.fixture
        path = Path(f.row()["reference"])
        path.write_bytes(path.read_bytes() + b"changed")
        f.plan["input_identities"][str(path)] = GATE.identity(path)
        self.reject("modified after diagnostic startup")


if __name__ == "__main__":
    unittest.main()
