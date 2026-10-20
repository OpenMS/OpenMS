#!/usr/bin/env python3
"""Make a short OpenDIA benchmark report and two SVG plots."""

import argparse
import html
import json
import re
import subprocess
from pathlib import Path


def load_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def parse_seconds(value):
    parts = value.split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    if len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    return float(parts[0])


def read_resources(path):
    text = path.read_text(encoding="utf-8")
    elapsed = re.search(r"Elapsed \(wall clock\).*?: ([0-9:.]+)", text)
    memory = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    status = re.search(r"Exit status: (\d+)", text)
    if not elapsed or not memory or not status:
        raise SystemExit(f"Cannot read resource data from {path}")
    resource = {
        "seconds": parse_seconds(elapsed.group(1)),
        "memory_mib": int(memory.group(1)) / 1024,
        "status": int(status.group(1)),
    }
    if resource["status"] != 0:
        raise SystemExit(f"Benchmark failed in {path}")
    return resource


def read_timsim(root):
    return {
        "metrics": load_json(root / "validation.json"),
        "target": read_resources(root / "target_only/logs/resources.txt"),
        "entrapment": read_resources(root / "entrapment/logs/resources.txt"),
        "revision": (root / "target_only/meta/source_revision.txt").read_text().strip(),
        "binary": (root / "target_only/meta/OpenDIA.sha256").read_text().split()[0],
    }


def percent_drop(old, new):
    return 100 * (old - new) / old


def change_text(name, old, new, measure):
    change = abs(percent_drop(old, new))
    if new <= old:
        return f"The {name} used {change:.0f}% less {measure}."
    return f"The {name} used {change:.0f}% more {measure}."


def format_duration(seconds):
    hours, remainder = divmod(round(seconds), 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours} h {minutes} min {seconds} s"
    if minutes:
        return f"{minutes} min {seconds} s"
    return f"{seconds} s"


def write_chart(path, title, unit, labels, old_values, new_values):
    width, height = 680, 320
    left, top, bottom = 80, 55, 55
    chart_height = height - top - bottom
    maximum = max(old_values + new_values) * 1.15
    group_width = (width - left - 30) / len(labels)
    bar_width = 58
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="18">{html.escape(title)}</text>',
        f'<line x1="{left}" y1="{top + chart_height}" x2="{width - 20}" y2="{top + chart_height}" stroke="#333"/>',
    ]
    for index, label in enumerate(labels):
        center = left + group_width * (index + 0.5)
        for offset, value, color in (
            (-bar_width / 2, old_values[index], "#999999"),
            (bar_width / 2, new_values[index], "#2b6cb0"),
        ):
            bar_height = chart_height * value / maximum
            x = center + offset - bar_width / 2
            y = top + chart_height - bar_height
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width}" height="{bar_height:.1f}" fill="{color}"/>')
            parts.append(f'<text x="{x + bar_width / 2:.1f}" y="{y - 7:.1f}" text-anchor="middle" font-family="sans-serif" font-size="12">{value:.1f}</text>')
        parts.append(f'<text x="{center:.1f}" y="{height - 28}" text-anchor="middle" font-family="sans-serif" font-size="13">{html.escape(label)}</text>')
    parts.extend(
        [
            f'<text x="18" y="{top + chart_height / 2}" transform="rotate(-90 18 {top + chart_height / 2})" text-anchor="middle" font-family="sans-serif" font-size="13">{html.escape(unit)}</text>',
            '<rect x="420" y="42" width="14" height="14" fill="#999999"/>',
            '<text x="440" y="54" font-family="sans-serif" font-size="12">Before update</text>',
            '<rect x="540" y="42" width="14" height="14" fill="#2b6cb0"/>',
            '<text x="560" y="54" font-family="sans-serif" font-size="12">Current</text>',
            "</svg>",
        ]
    )
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("current_timsim", type=Path)
    parser.add_argument("--previous-timsim", type=Path, required=True)
    parser.add_argument("--develop-revision")
    parser.add_argument("--astral-metrics", type=Path)
    parser.add_argument("--astral-time", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    if bool(args.astral_metrics) != bool(args.astral_time):
        parser.error("Use --astral-metrics and --astral-time together")

    current = read_timsim(args.current_timsim)
    previous = read_timsim(args.previous_timsim)
    if args.develop_revision:
        source_dir = Path(__file__).resolve().parents[3]
        result = subprocess.run(
            [
                "git",
                "-C",
                source_dir,
                "merge-base",
                "--is-ancestor",
                args.develop_revision,
                current["revision"],
            ],
            check=False,
        )
        if result.returncode != 0:
            raise SystemExit("The develop revision is not in the benchmark source")
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    labels = ["Target only", "Target and entrapment"]
    old_time = [previous["target"]["seconds"], previous["entrapment"]["seconds"]]
    new_time = [current["target"]["seconds"], current["entrapment"]["seconds"]]
    old_memory = [previous["target"]["memory_mib"], previous["entrapment"]["memory_mib"]]
    new_memory = [current["target"]["memory_mib"], current["entrapment"]["memory_mib"]]
    write_chart(output / "runtime.svg", "OpenDIA run time", "Seconds", labels, old_time, new_time)
    write_chart(output / "peak_memory.svg", "OpenDIA peak memory", "MiB", labels, old_memory, new_memory)

    metrics = current["metrics"]
    previous_metrics = previous["metrics"]
    combined_time = change_text("two searches together", sum(old_time), sum(new_time), "time")
    target_time = change_text("target search", old_time[0], new_time[0], "time")
    entrapment_time = change_text("mixed search", old_time[1], new_time[1], "time")
    target_memory = change_text("target search", old_memory[0], new_memory[0], "peak memory")
    entrapment_memory = change_text("mixed search", old_memory[1], new_memory[1], "peak memory")

    lines = [
        "# OpenDIA benchmark report",
        "",
        f"The source before the update is `{previous['revision'][:12]}`.",
        f"The current source is `{current['revision'][:12]}`. The OpenDIA binary is `{current['binary'][:12]}`.",
    ]
    if args.develop_revision:
        lines.append(f"This branch includes develop `{args.develop_revision[:12]}`.")
    lines.extend(
        [
            "Each TimSim value is from one run with four threads.",
            "",
            "## Main result",
            "",
            f"{combined_time}",
            f"{target_time} {entrapment_time}",
            f"{target_memory} {entrapment_memory}",
            "",
            "| Test | Before update | Current branch |",
            "|---|---:|---:|",
            f"| Target run time | {old_time[0]:.2f} s | {new_time[0]:.2f} s |",
            f"| Target peak memory | {old_memory[0]:.0f} MiB | {new_memory[0]:.0f} MiB |",
            f"| Target and entrapment run time | {old_time[1]:.2f} s | {new_time[1]:.2f} s |",
            f"| Target and entrapment peak memory | {old_memory[1]:.0f} MiB | {new_memory[1]:.0f} MiB |",
            "",
            "![Run time](runtime.svg)",
            "",
            "![Peak memory](peak_memory.svg)",
            "",
            "## Result quality",
            "",
            "| Measure | Before update | Current branch |",
            "|---|---:|---:|",
            f"| Target recovery | {previous_metrics['target_identification']['overall_recovery']:.1%} | {metrics['target_identification']['overall_recovery']:.1%} |",
            f"| Precursor match | {previous_metrics['quantification']['precursor_realized_vs_observed_pearson']:.4f} | {metrics['quantification']['precursor_realized_vs_observed_pearson']:.4f} |",
            f"| Protein match | {previous_metrics['quantification']['protein_realized_vs_observed_pearson']:.4f} | {metrics['quantification']['protein_realized_vs_observed_pearson']:.4f} |",
            f"| Entrapment rate | {previous_metrics['entrapment_export']['pooled_empirical_fdp']:.3%} | {metrics['entrapment_export']['pooled_empirical_fdp']:.3%} |",
            f"| Target rows in mixed search | {previous_metrics['entrapment_export']['target_rows']:,} | {metrics['entrapment_export']['target_rows']:,} |",
        ]
    )

    if args.astral_metrics and args.astral_time:
        astral = load_json(args.astral_metrics)
        astral_resources = read_resources(args.astral_time)
        identification = astral["identification"]
        lines.extend(
            [
                "",
                "## Astral reference",
                "",
                f"The reference run found {identification['target_rows']:,} target rows and {identification['entrapment_rows']:,} entrapment rows.",
                f"The run time was {format_duration(astral_resources['seconds'])}. Peak memory was {astral_resources['memory_mib'] / 1024:.1f} GiB.",
                "This result is from older code. This test did not run the full Astral data again.",
            ]
        )

    (output / "FINAL_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
