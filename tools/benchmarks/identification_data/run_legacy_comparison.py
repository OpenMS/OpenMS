"""Measure identification I/O in fresh Linux processes and verify roundtrips."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import time


def dimensions(value):
    try:
        rows, runs = (int(part) for part in value.split(":"))
        if rows < runs or runs < 1:
            raise ValueError
        return rows, runs
    except ValueError as error:
        raise argparse.ArgumentTypeError("Expected ROWS:RUNS with ROWS >= RUNS > 0") from error


def positive(value):
    result = int(value)
    if result < 1 or result > 2147483647:
        raise argparse.ArgumentTypeError("Expected a positive int32 value")
    return result


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("binary", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--formats", nargs="+", choices=["idxml", "parquet", "oms", "native"],
                    default=["idxml", "parquet", "oms", "native"])
parser.add_argument("--cases", nargs="+", type=dimensions,
                    default=[(1000, 2), (100000, 1), (1000000, 1), (1000000, 1000)])
parser.add_argument("--native-threads", nargs="+", type=positive, default=[1])
parser.add_argument("--read-repetitions", type=positive, default=3)
args = parser.parse_args()
if len(set(args.formats)) != len(args.formats) or len(set(args.cases)) != len(args.cases):
    parser.error("Formats and cases must be unique")
if len(set(args.native_threads)) != len(args.native_threads):
    parser.error("Thread counts must be unique")
binary = args.binary.resolve()
out = args.output.resolve()
out.mkdir(parents=True, exist_ok=False)
env = os.environ.copy()
env["OMP_NUM_THREADS"] = "1"
records = []
(out / "method.json").write_text(json.dumps({
    "formats": args.formats, "cases": args.cases, "native_threads": args.native_threads,
    "read_warmups": 1, "read_repetitions": args.read_repetitions, "write_repetitions": 1,
    "OMP_NUM_THREADS": "1", "cache": "warm; caches not dropped",
}, indent=2) + "\n")


def measure(label, command, threads):
    print(label, flush=True)
    start = time.monotonic()
    with (out / (label + ".log")).open("w") as log:
        child = subprocess.Popen([str(binary)] + [str(x) for x in command], env=env,
                                 stdout=log, stderr=subprocess.STDOUT, text=True)
        _, status, usage = os.wait4(child.pid, 0)
        child.returncode = os.waitstatus_to_exitcode(status)
    output = (out / (label + ".log")).read_text()
    phases, digests = {}, []
    for line in output.splitlines():
        if line.startswith("TIME\t"):
            _, phase, value = line.split("\t")
            phases[phase] = float(value)
        if line.startswith("DIGEST\t"):
            digests.append(line.split("\t")[1:])
    record = {"label": label, "args": [str(x) for x in command], "threads": threads,
              "code": child.returncode, "elapsed": time.monotonic() - start,
              "rss_kib": usage.ru_maxrss, "phases": phases, "digests": digests}
    record["user_cpu_s"] = usage.ru_utime
    record["system_cpu_s"] = usage.ru_stime
    records.append(record)
    (out / "measurements.json").write_text(json.dumps(records, indent=2) + "\n")
    print(json.dumps(record), flush=True)
    if child.returncode:
        raise RuntimeError(label + " failed:\n" + output[-3000:])
    if not digests or any(d != digests[0] for d in digests):
        raise RuntimeError("Conversion digest mismatch " + label)
    return digests[0]


variants = []
for fmt in args.formats:
    for threads in args.native_threads if fmt == "native" else [1]:
        label = fmt if fmt != "native" or args.native_threads == [1] else f"native-t{threads}"
        variants.append((label, fmt, threads))
for rows, runs in args.cases:
    case = f"{rows}-{runs}"
    expected = None
    paths = {label: out / (case + "-" + label) for label, _, _ in variants}
    for label, fmt, threads in variants:
        command = ["write", fmt, paths[label], rows, runs]
        if fmt == "native" and threads != 1:
            command.append(threads)
        got = measure(case + "-" + label + "-write", command, threads)
        if expected is None:
            expected = got
        if got != expected:
            raise RuntimeError("Input mismatch " + label)
    for repetition in range(args.read_repetitions + 1):
        shift = repetition % len(variants)
        for label, fmt, threads in variants[shift:] + variants[:shift]:
            command = ["read", fmt, paths[label]]
            if fmt == "native" and threads != 1:
                command.append(threads)
            got = measure(case + "-" + label + "-read-" + str(repetition), command, threads)
            if got != expected:
                raise RuntimeError(f"Roundtrip mismatch {label}: {got} versus {expected}")
    sizes = {label: sum(p.stat().st_size for p in path.rglob("*") if p.is_file())
             if path.is_dir() else path.stat().st_size for label, path in paths.items()}
    (out / (case + "-sizes.json")).write_text(json.dumps(sizes, indent=2) + "\n")
print("Comparison completed", flush=True)
