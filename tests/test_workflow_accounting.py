import importlib.util
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest

SPEC = importlib.util.spec_from_file_location("workflow_accounting", Path(__file__).resolve().parents[1] /
                                            "scripts/benchmark/summarize_workflow.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class WorkflowAccountingTests(unittest.TestCase):
    def allocation(self, job="1", **updates):
        return dict(dict(JobIDRaw=job,State="COMPLETED",ExitCode="0:0",actual_cpu_seconds=100,
                         allocated_cpu_hours=1,MaxRSS="",NodeList="node"),**updates)

    def task(self, **updates):
        return dict(dict(name="extractCallsChunk (sample)",process="extractCallsChunk",stage="extraction",
                         status="COMPLETED",exit="0"),**updates)

    def test_annotation_compiler_cost_is_retained_in_burdens(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as temporary:
            trace = Path(temporary) / "compiler.tsv"
            trace.write_text("native_id\tname\tstatus\texit\tpeak_rss\n"
                             "10\tcompileCoverageAnnotator\tCOMPLETED\t0\t1 GB\n"
                             "11\tcompileCoverageAnnotator\tFAILED\t1\t1 GB\n")
            jobs = MODULE.trace_jobs([trace])
            report = MODULE.summarize(jobs, [dict(records=[
                self.allocation("10", actual_cpu_seconds=4),
                self.allocation("11", State="FAILED", ExitCode="1:0", actual_cpu_seconds=1)])])
            self.assertEqual(report["stages"]["burdens"]["observed_actual_cpu_seconds"], 5)
            self.assertEqual(report["stages"]["burdens"]["observed_unsuccessful_cpu_seconds"], 1)
            self.assertEqual(report["scopes"]["common_downstream"]["observed_actual_cpu_seconds"], 5)

    def test_resumed_job_is_counted_once_and_steps_are_not_added_twice(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as temporary:
            one, two = Path(temporary) / "one.tsv", Path(temporary) / "two.tsv"
            header = "native_id\tname\tstatus\texit\tpeak_rss\n"
            one.write_text(header + "10\textractCallsChunk (sample)\tCOMPLETED\t0\t1 GB\n")
            two.write_text(header + "10\textractCallsChunk (sample)\tCACHED\t0\t1 GB\n")
            jobs = MODULE.trace_jobs([one, two])
            allocation = dict(JobIDRaw="10", State="COMPLETED", ExitCode="0:0", actual_cpu_seconds=100,
                              allocated_cpu_hours=1, MaxRSS="", NodeList="node")
            step = {**allocation, "JobIDRaw": "10.batch", "MaxRSS": "2048K"}
            report = MODULE.summarize(jobs, [dict(records=[allocation, step])])
            self.assertEqual(len(report["tasks"]), 1)
            self.assertEqual(report["stages"]["extraction"]["observed_actual_cpu_seconds"], 100)
            self.assertEqual(report["stages"]["extraction"]["slurm_peak_rss_kib"], 2048)
            self.assertEqual(report["unresolved_jobs"], [])

    def test_failed_work_counts_as_overhead_and_missing_cpu_stays_unresolved(self):
        base = dict(name="splitBAM (sample)", process="splitBAM", stage="bam_dispatch", status="FAILED", exit="1")
        jobs = {"1": base, "2": {**base, "status": "COMPLETED", "exit": "0"}, "3": base}
        failed = dict(JobIDRaw="1", State="FAILED", ExitCode="1:0", actual_cpu_seconds=20,
                      allocated_cpu_hours=.5, MaxRSS="100K", NodeList="node")
        incomplete = {**failed, "JobIDRaw": "2", "State": "COMPLETED", "ExitCode": "0:0",
                      "actual_cpu_seconds": None}
        result = MODULE.summarize(jobs, [dict(records=[failed, incomplete])])
        self.assertEqual(result["unresolved_jobs"], ["2", "3"])
        self.assertEqual(result["observed_actual_cpu_hours"], 20 / 3600)
        self.assertEqual(result["observed_successful_cpu_hours"], 0)
        self.assertIsNone(result["tasks"][2]["actual_cpu_seconds"])

    def test_signal_in_step_marks_allocation_cpu_as_lower_bound(self):
        allocation = self.allocation()
        step = self.allocation("1.batch",State="CANCELLED by 123",ExitCode="0:15",actual_cpu_seconds=80)
        result = MODULE.summarize({"1":self.task()},[dict(records=[allocation,step])])
        self.assertEqual(result["observed_actual_cpu_hours"],100/3600)
        self.assertEqual(result["unresolved_jobs"],["1"])
        self.assertEqual(result["cpu_lower_bound_jobs"],["1"])
        self.assertIn("signal_interruption_may_omit_child_cpu",result["tasks"][0]["cpu_incomplete_reasons"])
        self.assertEqual(result["tasks"][0]["interrupted_accounting_records"],["1.batch"])

    def test_signal_only_exit_code_and_truncated_oom_state(self):
        for updates in (dict(State="FAILED",ExitCode="0:9"),dict(State="OUT_OF_ME+",ExitCode="1:0")):
            result = MODULE.summarize({"1":self.task(status="FAILED",exit="1")},
                                     [dict(records=[self.allocation(**updates)])])
            self.assertTrue(result["tasks"][0]["terminal_accounting"])
            self.assertFalse(result["tasks"][0]["cpu_measurement_complete"])
            self.assertEqual(result["cpu_lower_bound_jobs"],["1"])

    def test_lagging_allocation_and_live_steps_remain_unresolved(self):
        result = MODULE.summarize({"1":self.task()},[dict(records=[self.allocation(actual_cpu_seconds=0),
                                self.allocation("1.batch",actual_cpu_seconds=50)])])
        self.assertEqual(result["observed_actual_cpu_hours"],0)
        self.assertEqual(result["unresolved_jobs"],["1"])
        self.assertIn("allocation_cpu_below_reported_step_cpu",result["tasks"][0]["cpu_incomplete_reasons"])
        result = MODULE.summarize({"1":self.task()},[dict(records=[self.allocation(),
                                self.allocation("1.batch",State="RUNNING")])])
        self.assertIn("nonterminal_accounting_steps",result["tasks"][0]["cpu_incomplete_reasons"])

    def test_later_snapshot_resolves_missing_cpu_but_keeps_interrupt_history(self):
        result = MODULE.summarize({"1":self.task()},[
            dict(records=[self.allocation(actual_cpu_seconds=None)]),dict(records=[self.allocation()])])
        self.assertEqual(result["unresolved_jobs"],[])
        result = MODULE.summarize({"1":self.task()},[
            dict(records=[self.allocation(State="PREEMPTED")]),dict(records=[self.allocation()])])
        self.assertEqual(result["cpu_lower_bound_jobs"],["1"])

    def test_allocation_below_combined_steps_stays_unresolved_without_imputation(self):
        records = [self.allocation(), self.allocation("1.batch",actual_cpu_seconds=60),
                   self.allocation("1.0",actual_cpu_seconds=60),
                   self.allocation("1.extern",actual_cpu_seconds=0)]
        # Repeated snapshots must not duplicate any step in the consistency sum.
        result = MODULE.summarize({"1":self.task()},[dict(records=records),dict(records=records)])
        self.assertEqual(result["unresolved_jobs"],["1"])
        self.assertEqual(result["cpu_lower_bound_jobs"],["1"])
        self.assertEqual(result["tasks"][0]["reported_step_cpu_seconds"],120)
        self.assertEqual(result["observed_actual_cpu_hours"],100/3600)
        result = MODULE.summarize({"1":self.task()},[dict(records=records),
            dict(records=[self.allocation(actual_cpu_seconds=120)])])
        self.assertEqual(result["unresolved_jobs"],[])
        self.assertEqual(result["observed_actual_cpu_hours"],120/3600)

    def test_long_cpu_display_truncation_is_allowed_but_real_gap_is_not(self):
        allocation = self.allocation(actual_cpu_seconds=3600,TotalCPU="01:00:00")
        for last_step, unresolved in ((.8,[]),(2.5,["1"])):
            result = MODULE.summarize({"1":self.task()},[dict(records=[allocation,
                self.allocation("1.batch",actual_cpu_seconds=3599.8),
                self.allocation("1.0",actual_cpu_seconds=last_step)])])
            self.assertEqual(result["unresolved_jobs"],unresolved)
            self.assertEqual(result["observed_actual_cpu_hours"],1)

    def test_preparation_downstream_and_unsuccessful_attempts_separate(self):
        jobs = {"1":self.task(process="installBSgenome",stage="preparation"),
                "2":self.task(status="FAILED",exit="1"),"3":self.task()}
        result = MODULE.summarize(jobs,[dict(records=[self.allocation(actual_cpu_seconds=300),
            self.allocation("2",State="FAILED",ExitCode="1:0",actual_cpu_seconds=20),
            self.allocation("3",actual_cpu_seconds=100)])])
        self.assertEqual(result["scopes"]["preparation"]["observed_actual_cpu_seconds"],300)
        downstream = result["scopes"]["common_downstream"]
        self.assertEqual(downstream["observed_actual_cpu_seconds"],120)
        self.assertEqual(downstream["observed_successful_cpu_seconds"],100)
        self.assertEqual(downstream["observed_unsuccessful_cpu_seconds"],20)

    def test_exclude_cached_requires_one_invocation(self):
        result = subprocess.run([sys.executable,str(SPEC.origin),"--trace","one.tsv","--trace","two.tsv",
                                 "--exclude-cached","--job-ids-only"],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn("requires exactly one --trace",result.stderr)
        result = MODULE.summarize({"1":self.task(status="CACHED")},
                                 [dict(records=[self.allocation()])],include_cached=False)
        self.assertEqual(result["tasks"],[])


if __name__ == "__main__":
    unittest.main()
