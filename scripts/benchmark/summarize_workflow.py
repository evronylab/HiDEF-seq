#!/usr/bin/env python3
"""Join Nextflow traces to saved Slurm accounting without counting steps twice.

Cached tasks retain their original job IDs: include their recorded workload cost
by default, or exclude it explicitly for one supplied resumed invocation.
Preparation reused through storeDir has no trace row and no measured cost here.
"""
import argparse
import csv
import datetime
import json
from pathlib import Path
import re


STAGES = {
    **dict.fromkeys(("installBSgenome", "prepareReferenceSummary", "extractGenomeTrinucleotides",
                    "processGermlineVCFs", "processGermlineBAMs", "prepareGermlineCoverageFilters",
                    "prepareRegionFilters"), "preparation"),
    **dict.fromkeys(("ccsChunk", "mergeCCSchunks", "filterAdapter", "makeBarcodesFasta", "countZMWs",
                    "limaDemux", "mergeDemuxBams", "pbmm2Align", "verifyBAMID",
                    "mergeAlignedSampleBAMs"), "read_processing"),
    **dict.fromkeys(("countAnalysisZMWs", "compileBamDispatcher", "splitBAM"), "bam_dispatch"),
    "extractCallsChunk": "extraction",
    "filterCallsChunkChromgroupFiltergroup": "filtering",
    "calculateBurdensChromgroupFiltergroup": "burdens",
    "outputResultsSample": "output",
}
TERMINAL_STATES = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
                   "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED"}
INTERRUPTED_STATES = {"CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED",
                      "BOOT_FAIL", "DEADLINE", "REVOKED"}


def state_name(row):
    state = row.get("State", "").split()[0] if row.get("State") else ""
    if state.endswith("+"):
        # sacct can truncate long state names, e.g. OUT_OF_ME+.
        matches = [item for item in TERMINAL_STATES if item.startswith(state[:-1])]
        return matches[0] if len(matches) == 1 else state
    return state


def signal_interrupted(row):
    signal = row.get("ExitCode", "").partition(":")[2]
    return state_name(row) in INTERRUPTED_STATES or (signal.isdigit() and int(signal) != 0)


def trace_jobs(paths):
    jobs = {}
    for path in paths:
        with Path(path).open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                job = row["native_id"]
                if not re.fullmatch(r"\d+(?:_\d+)?", job):
                    raise ValueError("Expected a Slurm job ID: " + repr(job))
                process = row["name"].split(" (", 1)[0].rsplit(":", 1)[-1]
                if process not in STAGES:
                    raise ValueError("Add an explicit stage mapping for process " + process)
                if job in jobs and jobs[job]["name"] != row["name"]:
                    raise ValueError("Conflicting trace task names for job " + job)
                # An executed record wins over its later cached appearance.
                if job not in jobs or row["status"] != "CACHED":
                    jobs[job] = {**row, "process": process, "stage": STAGES[process]}
    return jobs


def rss_kib(value):
    if value in (None, "", "Unknown", "N/A"):
        return None
    match = re.fullmatch(r"([0-9.]+)([KMGT]?)", value)
    if not match:
        raise ValueError("Unexpected Slurm RSS (collect with --units=K): " + value)
    return float(match[1]) * {"": 1, "K": 1, "M": 1024, "G": 1024**2, "T": 1024**3}[match[2]]


def summarize(jobs, accounting, include_cached=True):
    records = {}
    interrupted_records = set()
    for snapshot in accounting:
        # Later supplied accounting snapshots replace earlier incomplete ones.
        for row in snapshot["records"]:
            records[row["JobIDRaw"]] = row
            if signal_interrupted(row):
                interrupted_records.add(row["JobIDRaw"])
    tasks, stages, unresolved, lower_bound = [], {}, [], []
    for job, trace in jobs.items():
        if not include_cached and trace["status"] == "CACHED":
            continue
        allocation = records.get(job)
        steps = [r for key, r in records.items() if key.startswith(job + ".")]
        rss = [rss_kib(r.get("MaxRSS")) for r in ([allocation] if allocation else []) + steps]
        rss = [value for value in rss if value is not None]
        complete = allocation is not None and state_name(allocation) in TERMINAL_STATES
        cpu = allocation.get("actual_cpu_seconds") if allocation else None
        allocated = allocation.get("allocated_cpu_hours") if allocation else None
        interrupted = sorted(key for key in interrupted_records if key == job or key.startswith(job + "."))
        step_cpus = [r["actual_cpu_seconds"] for r in steps if r.get("actual_cpu_seconds") is not None]
        # Small allowance for separately rounded sacct fields; no CPU is imputed.
        lagging = cpu is not None and any(value > cpu + 0.01 for value in step_cpus)
        active_steps = [r["JobIDRaw"] for r in steps if state_name(r) not in TERMINAL_STATES]
        reasons = []
        if not complete:
            reasons.append("missing_or_nonterminal_allocation")
        if cpu is None:
            reasons.append("missing_actual_cpu")
        if interrupted:
            reasons.append("signal_interruption_may_omit_child_cpu")
        if lagging:
            reasons.append("allocation_cpu_below_reported_step_cpu")
        if active_steps:
            reasons.append("nonterminal_accounting_steps")
        cpu_complete = not reasons
        if not cpu_complete or allocated is None:
            unresolved.append(job)
        if interrupted or lagging:
            lower_bound.append(job)
        success = bool(allocation and state_name(allocation) == "COMPLETED" and
                       allocation["ExitCode"] == "0:0" and trace["exit"] == "0" and
                       trace["status"] in {"COMPLETED", "CACHED"})
        task = dict(job_id=job, name=trace["name"], process=trace["process"], stage=trace["stage"],
                    trace_status=trace["status"], accounting_state=allocation.get("State") if allocation else None,
                    successful=success, terminal_accounting=complete, actual_cpu_seconds=cpu,
                    cpu_measurement_complete=cpu_complete, cpu_incomplete_reasons=reasons,
                    interrupted_accounting_records=interrupted, nonterminal_steps=active_steps,
                    allocated_cpu_hours=allocated, slurm_peak_rss_kib=max(rss) if rss else None,
                    trace_peak_rss=trace.get("peak_rss"), node=allocation.get("NodeList") if allocation else None)
        tasks.append(task)
        stage = stages.setdefault(trace["stage"], dict(tasks=0, successful_tasks=0, unresolved_jobs=[],
                                  observed_actual_cpu_seconds=0.0, observed_successful_cpu_seconds=0.0,
                                  observed_unsuccessful_cpu_seconds=0.0,
                                  observed_allocated_cpu_hours=0.0, slurm_peak_rss_kib=None))
        stage["tasks"] += 1
        stage["successful_tasks"] += success
        if job in unresolved:
            stage["unresolved_jobs"].append(job)
        if cpu is not None:
            stage["observed_actual_cpu_seconds"] += cpu
            if success:
                stage["observed_successful_cpu_seconds"] += cpu
            else:
                stage["observed_unsuccessful_cpu_seconds"] += cpu
        if allocated is not None:
            stage["observed_allocated_cpu_hours"] += allocated
        if rss:
            stage["slurm_peak_rss_kib"] = max([stage["slurm_peak_rss_kib"] or 0, *rss])
    scopes = {}
    for name, selected in (("preparation", [s for key, s in stages.items() if key == "preparation"]),
                           ("common_downstream", [s for key, s in stages.items() if key != "preparation"])):
        scopes[name] = dict(tasks=sum(s["tasks"] for s in selected),
                           observed_actual_cpu_seconds=sum(s["observed_actual_cpu_seconds"] for s in selected),
                           observed_successful_cpu_seconds=sum(s["observed_successful_cpu_seconds"] for s in selected),
                           observed_unsuccessful_cpu_seconds=sum(s["observed_unsuccessful_cpu_seconds"] for s in selected),
                           observed_allocated_cpu_hours=sum(s["observed_allocated_cpu_hours"] for s in selected),
                           unresolved_jobs=[job for s in selected for job in s["unresolved_jobs"]])
    return dict(tasks=tasks, stages=stages, scopes=scopes, unresolved_jobs=unresolved,
                cpu_lower_bound_jobs=lower_bound,
                observed_actual_cpu_hours=sum(s["observed_actual_cpu_seconds"] for s in stages.values()) / 3600,
                observed_successful_cpu_hours=sum(s["observed_successful_cpu_seconds"] for s in stages.values()) / 3600,
                observed_unsuccessful_cpu_hours=sum(s["observed_unsuccessful_cpu_seconds"] for s in stages.values()) / 3600,
                observed_allocated_cpu_hours=sum(s["observed_allocated_cpu_hours"] for s in stages.values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace", action="append", required=True, type=Path)
    parser.add_argument("--accounting", action="append", default=[], type=Path,
                        help="collect_sacct.py JSON; later snapshots supersede earlier ones")
    parser.add_argument("--job-ids-only", action="store_true")
    parser.add_argument("--exclude-cached", action="store_true")
    parser.add_argument("--preparation-scope", default="unknown",
                        choices=("cold", "verified-cache", "historical-unverified", "mixed", "unknown"),
                        help="Required provenance classification for interpreting preparation; default unknown")
    parser.add_argument("--workflow-complete", action="store_true",
                        help="Assert independently verified workflow termination; does not establish science/cache equivalence")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.exclude_cached and len(args.trace) != 1:
        parser.error("--exclude-cached requires exactly one --trace to define the invocation")
    jobs = trace_jobs(args.trace)
    if args.job_ids_only:
        print(",".join(jobs))
        return
    if not args.accounting or args.report is None:
        parser.error("--accounting and --report are required unless --job-ids-only")
    if args.report.exists():
        parser.error("report already exists")
    report = summarize(jobs, [json.loads(path.read_text()) for path in args.accounting], not args.exclude_cached)
    report.update(collected_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  traces=[str(p.resolve()) for p in args.trace], accounting=[str(p.resolve()) for p in args.accounting],
                  include_cached_workload=not args.exclude_cached,
                  preparation_scope=args.preparation_scope,
                  comparability_established=False,
                  workflow_terminal=args.workflow_complete,
                  partial=not args.workflow_complete or bool(report["unresolved_jobs"]),
                  notes=["Actual CPU uses allocation records only; never add their step records again.",
                         "Observed totals with missing/live accounting are partial, not complete workload costs.",
                         "storeDir reuse without a trace row and prior cache construction are unmeasured.",
                         "Preparation and common downstream costs are separate; scope labels alone do not establish comparability.",
                         "Actual totals include all supplied attempts; unsuccessful CPU includes failed retries and live/incomplete work.",
                         "Signal-interrupted steps may omit child CPU; reported totals remain observed lower bounds.",
                         "workflow-complete asserts termination only, not scientific or cache equivalence.",
                         "RSS is the maximum recorded task/step RSS, not simultaneous pipeline memory.",
                         "Allocated CPU-hours and actual CPU-hours are distinct; wall time is not summed here."])
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("partial", "unresolved_jobs", "observed_actual_cpu_hours",
                                            "observed_successful_cpu_hours", "observed_allocated_cpu_hours")}))


if __name__ == "__main__":
    main()
