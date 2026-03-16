#!/usr/bin/env python3
"""Compare two benchmark JSON result files and report regressions.

Usage:
    python compare_benchmarks.py baseline.json stable_abi.json
    python compare_benchmarks.py baseline.json stable_abi.json --threshold 5.0
"""
import argparse
import json
import sys
from collections import defaultdict


def load_results(path):
    """Load benchmark JSON and return {category: {name: median_s}}."""
    with open(path) as f:
        data = json.load(f)
    by_category = defaultdict(dict)
    for result in data.get("results", []):
        by_category[result["category"]][result["name"]] = result["median_s"]
    return dict(by_category)


def compare(baseline, current, warn_threshold=5.0, fail_threshold=15.0):
    """Compare two result dicts. Returns (rows, overall_delta, passed)."""
    rows = []
    all_deltas = []

    all_categories = sorted(set(baseline) | set(current))
    for cat in all_categories:
        base_benchmarks = baseline.get(cat, {})
        curr_benchmarks = current.get(cat, {})
        all_names = sorted(set(base_benchmarks) | set(curr_benchmarks))

        for name in all_names:
            base_val = base_benchmarks.get(name)
            curr_val = curr_benchmarks.get(name)

            if base_val is None:
                rows.append((cat, name, None, curr_val, None, "NEW"))
                continue
            if curr_val is None:
                rows.append((cat, name, base_val, None, None, "REMOVED"))
                continue

            if base_val > 0:
                delta_pct = ((curr_val - base_val) / base_val) * 100
            else:
                delta_pct = 0.0

            all_deltas.append(delta_pct)

            if delta_pct > fail_threshold:
                status = "FAIL"
            elif delta_pct > warn_threshold:
                status = "WARN"
            else:
                status = "PASS"

            rows.append((cat, name, base_val, curr_val, delta_pct, status))

    overall_delta = sum(all_deltas) / len(all_deltas) if all_deltas else 0.0
    passed = all(r[5] in ("PASS", "NEW") for r in rows)
    return rows, overall_delta, passed


def print_report(rows, overall_delta, passed):
    """Print color-coded comparison table."""
    colors = {
        "PASS": "\033[32m",  # green
        "WARN": "\033[33m",  # yellow
        "FAIL": "\033[31m",  # red
        "NEW": "\033[36m",   # cyan
        "REMOVED": "\033[90m",  # gray
    }
    reset = "\033[0m"

    print(f"\n{'Category':<25} {'Benchmark':<50} {'Base (ms)':>10} {'Curr (ms)':>10} {'Delta':>8} {'Status':>8}")
    print("-" * 115)

    for cat, name, base_val, curr_val, delta_pct, status in rows:
        color = colors.get(status, "")
        base_str = f"{base_val * 1000:.3f}" if base_val is not None else "N/A"
        curr_str = f"{curr_val * 1000:.3f}" if curr_val is not None else "N/A"
        delta_str = f"{delta_pct:+.1f}%" if delta_pct is not None else "---"
        print(f"{cat:<25} {name:<50} {base_str:>10} {curr_str:>10} {color}{delta_str:>8} {status:>8}{reset}")

    print("-" * 115)
    overall_color = colors["PASS"] if passed else colors["FAIL"]
    print(f"Overall mean delta: {overall_color}{overall_delta:+.1f}%{reset}")
    print(f"Result: {overall_color}{'PASSED' if passed else 'FAILED'}{reset}\n")


def main():
    parser = argparse.ArgumentParser(description="Compare benchmark results")
    parser.add_argument("baseline", help="Path to baseline JSON results")
    parser.add_argument("current", help="Path to current JSON results")
    parser.add_argument("--warn-threshold", type=float, default=5.0,
                        help="Percent regression to trigger warning (default: 5.0)")
    parser.add_argument("--fail-threshold", type=float, default=15.0,
                        help="Percent regression to trigger failure (default: 15.0)")
    args = parser.parse_args()

    baseline = load_results(args.baseline)
    current = load_results(args.current)
    rows, overall_delta, passed = compare(baseline, current, args.warn_threshold, args.fail_threshold)
    print_report(rows, overall_delta, passed)

    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
