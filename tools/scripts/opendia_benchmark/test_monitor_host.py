#!/usr/bin/env python3
"""Exercise real Linux workload sampling and the benchmark-only disk guard."""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    monitor = Path(__file__).with_name("monitor_host.py")
    with tempfile.TemporaryDirectory(dir=os.environ.get("JCODE_SCRATCH_DIR")) as temporary:
        root = Path(temporary)
        for guard in (False, True):
            output = root / str(guard)
            output.mkdir()
            job = subprocess.Popen(["bash", "-c", "sleep 30 & wait $!" if guard else "sleep 1 & wait $!"])
            try:
                sampler = subprocess.Popen([
                    sys.executable, monitor, str(job.pid), output, "--interval", "0.1",
                    "--min-free-disk-gib", "1000000000" if guard else "0",
                ])
                status = job.wait(timeout=5)
                assert sampler.wait(timeout=5) == 0
                data = json.loads((output / "host_workload.json").read_text())
                assert data["sample_count"] > 0
                assert data["logical_cpus"] >= 1
                assert data["mem_available_mib"]["min"] > 0
                assert bool(data["abort_reason"]) == guard
                assert (status != 0) == guard
                raw = [json.loads(line) for line in (output / "host_workload.jsonl").read_text().splitlines()]
                assert len(raw) == data["sample_count"]
                if not guard:
                    assert 0 <= data["host_cpu_busy_percent"]["mean"] <= 100
            finally:
                if job.poll() is None:
                    job.terminate()
                    job.wait(timeout=5)
    print("PASS: real host sampling, summary/raw agreement, clean shutdown, and low-disk guard")


if __name__ == "__main__":
    main()
