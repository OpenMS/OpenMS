#!/usr/bin/env python3
"""
Top-level runner for OpenMS Issue #10425:
Benchmark ModifiedSincSmoother against SavitzkyGolayFilter on Chromatograms.

Usage:
    python src/pyOpenMS/tests/benchmark_chromatogram_smoothing.py --quick
    python src/pyOpenMS/tests/benchmark_chromatogram_smoothing.py --help
"""

import sys
from pathlib import Path

# Add benchmark_smoothing to path
pkg_dir = Path(__file__).resolve().parent / "benchmark_smoothing"
if str(pkg_dir) not in sys.path:
    sys.path.insert(0, str(pkg_dir))

from run_benchmark import main

if __name__ == "__main__":
    main()
