#!/usr/bin/env python3
"""Sample Linux host load during one benchmark, without adjusting measured times."""

import argparse
import json
import math
import os
import signal
import statistics
import time
from pathlib import Path


def process_ticks(pid):
    """Count live descendant CPU ticks once per process, not once per thread."""
    ticks, descendants = 0, []
    pending = [pid]
    while pending:
        current = pending.pop()
        try:
            stat = Path(f"/proc/{current}/stat").read_text().rsplit(")", 1)[1].split()
            ticks += int(stat[11]) + int(stat[12])
            children = Path(f"/proc/{current}/task/{current}/children").read_text().split()
            pending.extend(int(child) for child in children)
            if current != pid:
                descendants.append(current)
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return ticks, descendants


def pressure(kind):
    try:
        first = Path(f"/proc/pressure/{kind}").read_text().splitlines()[0]
        return float(dict(item.split("=", 1) for item in first.split()[1:])["avg10"])
    except (OSError, KeyError, IndexError):
        return None


def host_sample(pid, output):
    cpu = [int(value) for value in Path("/proc/stat").read_text().splitlines()[0].split()[1:9]]
    memory = {}
    for line in Path("/proc/meminfo").read_text().splitlines():
        key, value = line.split(":", 1)
        memory[key] = int(value.split()[0])
    ticks, descendants = process_ticks(pid)
    sample = {
        "timestamp": time.time(),
        "load1": os.getloadavg()[0],
        "mem_available_mib": memory["MemAvailable"] / 1024,
        "swap_used_mib": (memory["SwapTotal"] - memory["SwapFree"]) / 1024,
        "cpu_pressure_some_percent": pressure("cpu"),
        "memory_pressure_some_percent": pressure("memory"),
        "io_pressure_some_percent": pressure("io"),
        "disk_free_gib": os.statvfs(output).f_bavail * os.statvfs(output).f_frsize / 1024**3,
    }
    return sample, cpu, ticks, descendants


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pid", type=int)
    parser.add_argument("output", type=Path)
    parser.add_argument("--interval", type=float, default=1)
    parser.add_argument("--min-free-disk-gib", type=float, default=3)
    args = parser.parse_args()
    if args.pid <= 0 or not math.isfinite(args.interval) or not math.isfinite(args.min_free_disk_gib) or args.interval <= 0 or args.min_free_disk_gib < 0:
        parser.error("interval must be positive and disk reserve must not be negative")
    samples = []
    previous_cpu = previous_ticks = previous_time = None
    cores = os.cpu_count() or 1
    abort_reason = None
    with (args.output / "host_workload.jsonl").open("w", encoding="utf-8") as handle:
        while Path(f"/proc/{args.pid}").exists():
            sample, cpu, ticks, descendants = host_sample(args.pid, args.output)
            now = time.monotonic()
            if previous_cpu is not None:
                elapsed = now - previous_time
                total = sum(cpu) - sum(previous_cpu)
                idle = cpu[3] + cpu[4] - previous_cpu[3] - previous_cpu[4]
                sample["host_cpu_busy_percent"] = 100 * (1 - idle / total) if total > 0 else None
                # This is a workload estimate, not a correction to elapsed time.
                job_cpu = 100 * (ticks - previous_ticks) / os.sysconf("SC_CLK_TCK") / elapsed / cores
                sample["job_cpu_host_percent"] = job_cpu if job_cpu >= 0 else None
                sample["other_cpu_busy_estimate_percent"] = (
                    max(0, sample["host_cpu_busy_percent"] - job_cpu)
                    if job_cpu >= 0 and total > 0 else None
                )
            samples.append(sample)
            if len(samples) % 60 == 0:
                print("JCODE_PROGRESS " + json.dumps({
                    "message": f"{args.output.parent.name}: {len(samples)} workload samples, {sample['disk_free_gib']:.1f} GiB disk free",
                }), flush=True)
            handle.write(json.dumps(sample, allow_nan=False) + "\n")
            handle.flush()
            if sample["disk_free_gib"] < args.min_free_disk_gib:
                abort_reason = f"Free disk space fell below {args.min_free_disk_gib:g} GiB"
                if not descendants:
                    time.sleep(0.1)
                    continue
                # Stop only this measured job's descendants. Keep GNU time alive to record failure.
                for child in reversed(descendants):
                    try:
                        os.kill(child, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                break
            previous_cpu, previous_ticks, previous_time = cpu, ticks, now
            time.sleep(args.interval)
    summary = {
        "sample_count": len(samples),
        "logical_cpus": cores,
        "interval_seconds": args.interval,
        "abort_reason": abort_reason,
        "note": "Host CPU includes the benchmark. Other CPU is an estimate. Elapsed times are not corrected for workload.",
    }
    for key in samples[0] if samples else ():
        if key == "timestamp":
            continue
        values = [sample[key] for sample in samples if sample.get(key) is not None]
        summary[key] = {"mean": statistics.fmean(values), "min": min(values), "max": max(values)} if values else None
    # CPU fields are absent in the first sample, so include them explicitly.
    for key in ("host_cpu_busy_percent", "job_cpu_host_percent", "other_cpu_busy_estimate_percent"):
        values = [sample[key] for sample in samples if sample.get(key) is not None]
        summary[key] = {"mean": statistics.fmean(values), "min": min(values), "max": max(values)} if values else None
    (args.output / "host_workload.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
