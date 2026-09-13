#!/usr/bin/env python3
"""Compare two benchmark JSON result files and report regressions.

Both files are produced by ``benchmark_pyopenms.py --json <file>``. The
comparison is per benchmark name within a category, so adding or removing
benchmarks between the two runs is tolerated (reported as NEW/REMOVED).

Usage:
    python compare_benchmarks.py baseline.json candidate.json
    python compare_benchmarks.py baseline.json candidate.json --fail-threshold 10
    python compare_benchmarks.py baseline.json candidate.json --markdown
"""
import argparse
import json
import statistics
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
    """Compare two result dicts. Returns (rows, summary)."""
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

            delta_pct = ((curr_val - base_val) / base_val) * 100 if base_val > 0 else 0.0
            all_deltas.append(delta_pct)

            if delta_pct > fail_threshold:
                status = "FAIL"
            elif delta_pct > warn_threshold:
                status = "WARN"
            else:
                status = "PASS"

            rows.append((cat, name, base_val, curr_val, delta_pct, status))

    summary = {
        "mean_delta": statistics.mean(all_deltas) if all_deltas else 0.0,
        # The median is the robust headline number: a handful of microbenchmarks
        # with huge relative deltas would otherwise dominate the mean.
        "median_delta": statistics.median(all_deltas) if all_deltas else 0.0,
        "n": len(all_deltas),
        "passed": all(r[5] in ("PASS", "NEW") for r in rows),
    }
    return rows, summary


def category_summary(rows):
    """Aggregate per-category median deltas, in the row order of first appearance."""
    per_cat = defaultdict(list)
    for cat, _name, _b, _c, delta, _status in rows:
        if delta is not None:
            per_cat[cat].append(delta)
    return [(cat, statistics.median(d), len(d)) for cat, d in per_cat.items()]


def print_report(rows, summary):
    """Print color-coded comparison table."""
    colors = {
        "PASS": "\033[32m",     # green
        "WARN": "\033[33m",     # yellow
        "FAIL": "\033[31m",     # red
        "NEW": "\033[36m",      # cyan
        "REMOVED": "\033[90m",  # gray
    }
    reset = "\033[0m"

    print(f"\n{'Category':<22} {'Benchmark':<52} {'Base (ms)':>11} {'Curr (ms)':>11} {'Delta':>8} {'Status':>8}")
    print("-" * 116)

    for cat, name, base_val, curr_val, delta_pct, status in rows:
        color = colors.get(status, "")
        base_str = f"{base_val * 1000:.3f}" if base_val is not None else "N/A"
        curr_str = f"{curr_val * 1000:.3f}" if curr_val is not None else "N/A"
        delta_str = f"{delta_pct:+.1f}%" if delta_pct is not None else "---"
        print(f"{cat:<22} {name:<52} {base_str:>11} {curr_str:>11} {color}{delta_str:>8} {status:>8}{reset}")

    print("-" * 116)
    print(f"\n{'Category':<22} {'Median delta':>13} {'N':>5}")
    print("-" * 42)
    for cat, med, n in category_summary(rows):
        color = colors["FAIL"] if med > 15 else colors["WARN"] if med > 5 else colors["PASS"]
        print(f"{cat:<22} {color}{med:>+12.1f}%{reset} {n:>5}")
    print("-" * 42)

    overall_color = colors["PASS"] if summary["passed"] else colors["FAIL"]
    print(f"\nOverall median delta: {overall_color}{summary['median_delta']:+.1f}%{reset} "
          f"(mean {summary['mean_delta']:+.1f}%, n={summary['n']})")
    print(f"Result: {overall_color}{'PASSED' if summary['passed'] else 'FAILED'}{reset}\n")


def print_markdown(rows, summary):
    """Print a GitHub-flavoured markdown report (per-category medians + outliers)."""
    print("| Category | Median delta | Benchmarks |")
    print("|---|---:|---:|")
    for cat, med, n in category_summary(rows):
        print(f"| {cat} | {med:+.1f}% | {n} |")
    print()
    print(f"Overall median delta: **{summary['median_delta']:+.1f}%** "
          f"(mean {summary['mean_delta']:+.1f}%, n={summary['n']})")
    print()
    print("| Category | Benchmark | Base (ms) | Curr (ms) | Delta |")
    print("|---|---|---:|---:|---:|")
    for cat, name, base_val, curr_val, delta_pct, status in rows:
        base_str = f"{base_val * 1000:.3f}" if base_val is not None else "N/A"
        curr_str = f"{curr_val * 1000:.3f}" if curr_val is not None else "N/A"
        delta_str = f"{delta_pct:+.1f}%" if delta_pct is not None else status
        print(f"| {cat} | {name} | {base_str} | {curr_str} | {delta_str} |")


def main():
    parser = argparse.ArgumentParser(description="Compare benchmark results")
    parser.add_argument("baseline", help="Path to baseline JSON results")
    parser.add_argument("current", help="Path to current JSON results")
    parser.add_argument("--warn-threshold", type=float, default=5.0,
                        help="Percent regression to trigger warning (default: 5.0)")
    parser.add_argument("--fail-threshold", type=float, default=15.0,
                        help="Percent regression to trigger failure (default: 15.0)")
    parser.add_argument("--markdown", action="store_true",
                        help="Emit a markdown report instead of the colored table")
    args = parser.parse_args()

    baseline = load_results(args.baseline)
    current = load_results(args.current)
    rows, summary = compare(baseline, current, args.warn_threshold, args.fail_threshold)
    if args.markdown:
        print_markdown(rows, summary)
    else:
        print_report(rows, summary)

    sys.exit(0 if summary["passed"] else 1)


if __name__ == "__main__":
    main()
