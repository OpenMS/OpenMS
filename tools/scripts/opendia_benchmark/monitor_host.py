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


def process_tree(pid):
    """Snapshot processes, including children created by non-main threads."""
    processes = {}
    pending = [pid]
    while pending:
        current = pending.pop()
        if current in processes:
            continue
        try:
            stat = Path(f"/proc/{current}/stat").read_text().rsplit(")", 1)[1].split()
            processes[current] = stat
            for task in Path(f"/proc/{current}/task").iterdir():
                try:
                    pending.extend(int(child) for child in (task / "children").read_text().split())
                except (FileNotFoundError, ProcessLookupError, PermissionError):
                    continue
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return processes


def process_ticks(pid):
    """Estimate cumulative tree CPU, retaining CPU of waited children."""
    processes = process_tree(pid)
    # cutime/cstime retain reaped child work. Snapshot races still make this an estimate.
    ticks = sum(sum(int(stat[index]) for index in (11, 12, 13, 14)) for stat in processes.values())
    descendants = [child for child in processes if child != pid]
    return ticks, descendants


def stop_descendants(pid, grace_seconds=2):
    """Stop only the measured tree, retaining pidfds if children become orphans."""
    tracked = {}
    deadline = time.monotonic() + grace_seconds
    initial = process_tree(pid)
    root_started = initial[pid][19] if pid in initial else None
    try:
        while True:
            # Pin all identities before signalling a parent that could orphan its children.
            processes = process_tree(pid)
            if pid not in processes or processes[pid][19] != root_started:
                processes = {}
            # Orphans may themselves create children during their TERM cleanup.
            for child, started in list(tracked):
                subtree = process_tree(child)
                if child in subtree and subtree[child][19] == started:
                    processes.update(subtree)
            new_descriptors = []
            for child, stat in processes.items():
                identity = (child, stat[19])  # PID plus kernel start time
                if child == pid or stat[0] == "Z" or identity in tracked:
                    continue
                try:
                    descriptor = os.pidfd_open(child)
                except ProcessLookupError:
                    continue
                try:
                    current = Path(f"/proc/{child}/stat").read_text().rsplit(")", 1)[1].split()
                    if current[19] != stat[19]:
                        os.close(descriptor)
                        continue
                except (FileNotFoundError, ProcessLookupError):
                    os.close(descriptor)
                    continue
                tracked[identity] = descriptor
                new_descriptors.append(descriptor)
            for descriptor in reversed(new_descriptors):
                try:
                    signal.pidfd_send_signal(descriptor, signal.SIGTERM)
                except ProcessLookupError:
                    pass

            alive = []
            for (child, started), descriptor in tracked.items():
                try:
                    stat = Path(f"/proc/{child}/stat").read_text().rsplit(")", 1)[1].split()
                    if stat[19] == started and stat[0] != "Z":
                        alive.append(descriptor)
                except (FileNotFoundError, ProcessLookupError):
                    continue
            if not alive and tracked:
                return True
            if not alive and not Path(f"/proc/{pid}").exists():
                return True
            now = time.monotonic()
            if now >= deadline:
                for descriptor in alive:
                    try:
                        signal.pidfd_send_signal(descriptor, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                if now >= deadline + 2:
                    return not alive
            time.sleep(0.05)
    finally:
        for descriptor in tracked.values():
            os.close(descriptor)


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
                # Leave GNU time alive to record failure, but do not abandon resistant children.
                if not stop_descendants(args.pid):
                    abort_reason += "; some descendants did not stop after SIGKILL"
                break
            previous_cpu, previous_ticks, previous_time = cpu, ticks, now
            time.sleep(args.interval)
    summary = {
        "sample_count": len(samples),
        "logical_cpus": cores,
        "interval_seconds": args.interval,
        "abort_reason": abort_reason,
        "note": "Host CPU includes the benchmark. Job CPU and other CPU are estimates from process snapshots. Elapsed times are not corrected for workload.",
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
    return 2 if abort_reason else 0


if __name__ == "__main__":
    raise SystemExit(main())
