#!/usr/bin/env python3
"""Exercise real Linux workload sampling and benchmark-only disk termination."""

import json
import os
import resource
import signal
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

from monitor_host import process_ticks, process_tree


def running(pid):
    tree = process_tree(pid)
    return pid in tree and tree[pid][0] != "Z"


def test_cpu_accounting():
    before = process_ticks(os.getpid())[0]
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    subprocess.run([
        sys.executable, "-c",
        "import time; end = time.process_time() + 0.25\nwhile time.process_time() < end: pass",
    ], check=True)
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    actual_cpu = usage_after.ru_utime + usage_after.ru_stime - usage_before.ru_utime - usage_before.ru_stime
    observed_cpu = (process_ticks(os.getpid())[0] - before) / os.sysconf("SC_CLK_TCK")
    assert observed_cpu >= actual_cpu - 0.04, (observed_cpu, actual_cpu)

    # A subprocess belongs to the creating thread's children list, not the main thread's.
    ready, done = threading.Event(), threading.Event()
    children = []

    def worker():
        child = subprocess.Popen(["sleep", "30"])
        children.append(child)
        ready.set()
        try:
            done.wait(timeout=5)
        finally:
            child.terminate()
            child.wait(timeout=5)

    thread = threading.Thread(target=worker)
    thread.start()
    try:
        assert ready.wait(timeout=5)
        _, descendants = process_ticks(os.getpid())
        assert descendants.count(children[0].pid) == 1
    finally:
        done.set()
        thread.join(timeout=5)
        assert not thread.is_alive()


def test_sampling(root, mode):
    monitor = Path(__file__).with_name("monitor_host.py")
    output = root / mode
    output.mkdir()
    guard = mode != "normal"
    child = None
    if mode == "orphan_thread":
        # TERM kills this parent but not its worker-thread child. The guard must retain it.
        program = """import subprocess, sys, threading, time

def worker():
    child = subprocess.Popen([sys.executable, '-c',
        'import os,signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); print(os.getpid(),flush=True); time.sleep(30)'])
    child.wait()
threading.Thread(target=worker).start()
time.sleep(30)
"""
    elif mode == "resistant":
        program = "import os,signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); print(os.getpid(),flush=True); time.sleep(30)"
    else:
        program = f"import os,time; print(os.getpid(),flush=True); time.sleep({30 if guard else 0.5})"
    job = subprocess.Popen([
        "/usr/bin/time", "-v", "-o", str(output / "resources.txt"), sys.executable, "-c", program,
    ], stdout=subprocess.PIPE, text=True)
    sampler = None
    try:
        child = int(job.stdout.readline())
        sampler = subprocess.Popen([
            sys.executable, monitor, str(job.pid), output, "--interval", "0.1",
            "--min-free-disk-gib", "1000000000" if guard else "0",
        ])
        status = job.wait(timeout=8)
        assert sampler.wait(timeout=8) == (2 if guard else 0)
        assert not running(child), (mode, child)
        data = json.loads((output / "host_workload.json").read_text())
        assert data["sample_count"] > 0
        assert data["logical_cpus"] >= 1
        assert data["mem_available_mib"]["min"] > 0
        assert bool(data["abort_reason"]) == guard
        assert (status != 0) == guard
        raw = [json.loads(line) for line in (output / "host_workload.jsonl").read_text().splitlines()]
        assert len(raw) == data["sample_count"]
        # GNU time itself must survive to write its resource/failure record.
        assert "Maximum resident set size" in (output / "resources.txt").read_text()
        if not guard:
            assert 0 <= data["host_cpu_busy_percent"]["mean"] <= 100
    finally:
        if sampler is not None and sampler.poll() is None:
            sampler.kill()
            sampler.wait(timeout=5)
        # Only these synthetic processes are candidates for test cleanup.
        descendants = process_tree(job.pid)
        if child is not None and running(child):
            descendants.update(process_tree(child))
        for pid in reversed(descendants):
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        job.wait(timeout=5)
        job.stdout.close()


def main():
    test_cpu_accounting()
    with tempfile.TemporaryDirectory(dir=os.environ.get("JCODE_SCRATCH_DIR")) as temporary:
        # An unrelated process must remain alive through every guard test.
        unrelated = subprocess.Popen(["sleep", "30"])
        try:
            for mode in ("normal", "cooperative", "resistant", "orphan_thread"):
                test_sampling(Path(temporary), mode)
                assert unrelated.poll() is None
        finally:
            unrelated.terminate()
            unrelated.wait(timeout=5)
    print("PASS: host sampling, waited-child CPU, threaded descendants, TERM/KILL and orphan guards, GNU time records, unrelated process safety")


if __name__ == "__main__":
    main()
