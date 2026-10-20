#!/usr/bin/env python3
"""Make a short OpenDIA benchmark report and SVG plots."""

import argparse
import html
import json
import os
import re
import shlex
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
    threads = set()
    revisions = set()
    binaries = set()
    for name in ("target_only", "entrapment"):
        meta = root / name / "meta"
        command = shlex.split((meta / "command.sh").read_text())
        threads.add(command[command.index("-threads") + 1])
        revisions.add((meta / "source_revision.txt").read_text().strip())
        binaries.add((meta / "OpenDIA.sha256").read_text().split()[0])
    if len(threads) != 1 or len(revisions) != 1 or len(binaries) != 1:
        raise SystemExit(f"The two TimSim searches have different provenance: {root}")
    return {
        "metrics": load_json(root / "validation.json"),
        "target": read_resources(root / "target_only/logs/resources.txt"),
        "entrapment": read_resources(root / "entrapment/logs/resources.txt"),
        "revision": revisions.pop(),
        "binary": binaries.pop(),
        "threads": threads.pop(),
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


def write_chart(path, title, unit, labels, series, decimals=1):
    width, height = 680, 320
    left, top, bottom = 80, 55, 55
    chart_height = height - top - bottom
    finite_values = [value for _, _, values in series for value in values if value is not None]
    maximum = max([0, *finite_values]) * 1.15 or 1
    minimum = min([0, *finite_values]) * 1.15
    scale = chart_height / (maximum - minimum)
    baseline = top + maximum * scale
    group_width = (width - left - 30) / len(labels)
    bar_width = min(58, group_width / (len(series) + 1))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="18">{html.escape(title)}</text>',
        f'<line x1="{left}" y1="{baseline}" x2="{width - 20}" y2="{baseline}" stroke="#333"/>',
    ]
    for index, label in enumerate(labels):
        center = left + group_width * (index + 0.5)
        for series_index, (_, color, values) in enumerate(series):
            offset = (series_index - (len(series) - 1) / 2) * bar_width * 1.15
            value = values[index]
            x = center + offset - bar_width / 2
            if value is None:
                parts.append(f'<text x="{x + bar_width / 2:.1f}" y="{baseline - 7:.1f}" text-anchor="middle" font-family="sans-serif" font-size="12">n/a</text>')
                continue
            bar_height = abs(value) * scale
            y = baseline - max(value, 0) * scale
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width}" height="{bar_height:.1f}" fill="{color}"/>')
            parts.append(f'<text x="{x + bar_width / 2:.1f}" y="{y - 7:.1f}" text-anchor="middle" font-family="sans-serif" font-size="12">{value:.{decimals}f}</text>')
        parts.append(f'<text x="{center:.1f}" y="{height - 28}" text-anchor="middle" font-family="sans-serif" font-size="13">{html.escape(label)}</text>')
    parts.append(f'<text x="18" y="{top + chart_height / 2}" transform="rotate(-90 18 {top + chart_height / 2})" text-anchor="middle" font-family="sans-serif" font-size="13">{html.escape(unit)}</text>')
    if len(series) > 1:
        legend_x = width - 130 * len(series)
        for index, (name, color, _) in enumerate(series):
            x = legend_x + index * 130
            parts.append(f'<rect x="{x}" y="42" width="14" height="14" fill="{color}"/>')
            parts.append(f'<text x="{x + 20}" y="54" font-family="sans-serif" font-size="12">{html.escape(name)}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def comparison(old, new):
    return [("Before", "#999999", old), ("Current", "#2b6cb0", new)]


def format_match(value):
    return f"{value:.3f}" if value is not None else "n/a"


def quality_report(output, metrics, astral):
    identification = metrics["target_identification"]
    mixed = metrics["entrapment_export"]
    quantification = metrics["quantification"]
    id_labels = ["TimSim target", "TimSim mixed"]
    id_values = [identification["identified_target_run_pairs"], mixed["target_rows"]]
    quant_labels = ["TimSim precursor", "TimSim protein"]
    quant_values = [quantification[f"{level}_realized_vs_observed_pearson"] for level in ("precursor", "protein")]
    lines = [
        "## Identification",
        "",
        "A row is one precursor in one run. TimSim has six runs. Astral has three runs. Do not compare their total counts as equal tests.",
        f"TimSim found {identification['identified_target_run_pairs']:,} of {identification['expected_target_run_pairs']:,} expected results ({identification['overall_recovery']:.1%}).",
        f"The mixed search found {mixed['target_rows']:,} target rows and {mixed['entrapment_rows']:,} entrapment rows.",
    ]
    if astral:
        ids = astral["identification"]
        id_labels.append("Astral reference")
        id_values.append(ids["target_rows"])
        lines.append(f"The older Astral run found {ids['target_rows']:,} target rows ({ids['unique_target_precursors']:,} distinct precursors) and {ids['entrapment_rows']:,} entrapment rows.")
        for level in ("precursor", "protein"):
            values = [item["pearson_log2"] for item in astral["quantification_repeatability"][level]["pairwise_log2_correlations"].values()]
            values = [value for value in values if value is not None]
            quant_labels.append(f"Astral {level}")
            quant_values.append(sum(values) / len(values) if values else None)
    write_chart(output / "identification.svg", "Target results", "Rows", id_labels, [("Results", "#2b6cb0", id_values)], decimals=0)
    write_chart(output / "quantification.svg", "Quantification match", "Correlation", quant_labels, [("Match", "#2b6cb0", quant_values)], decimals=3)
    lines.extend([
        "",
        "![Identification](identification.svg)",
        "",
        "## Quantification",
        "",
        "| Data | Precursor match | Protein match |",
        "|---|---:|---:|",
        f"| TimSim truth | {format_match(quant_values[0])} | {format_match(quant_values[1])} |",
    ])
    if astral:
        lines.append(f"| Astral reference | {format_match(quant_values[2])} | {format_match(quant_values[3])} |")
    lines.extend([
        "",
        "TimSim compares measured fold changes with known fold changes. A value of 1 means the same linear pattern, not exact amounts.",
        "",
        "![Quantification](quantification.svg)",
    ])
    if astral:
        lines.extend(["", "Astral uses the mean correlation between repeat runs. This measures repeatability, not accuracy.", "", "| Astral reference | Present in all 3 runs | Median variation |", "|---|---:|---:|"])
        for level in ("precursor", "protein"):
            group = astral["quantification_repeatability"][level]
            cv = f"{group['median_cv_complete']:.2%}" if group['median_cv_complete'] is not None else "n/a"
            lines.append(f"| {level.capitalize()} | {group['complete_entities']:,} / {group['entities']:,} | {cv} |")
    lines.extend(["", "## FDR estimate and calibration check", "", "| Data | q-value limit | Entrapment share |", "|---|---:|---:|", f"| TimSim mixed | {metrics['q_threshold']:.2%} | {mixed['pooled_empirical_fdp']:.3%} |"])
    datasets = [("TimSim", "#2b6cb0", metrics)]
    if astral:
        lines.append(f"| Astral reference | {astral['q_threshold']:.2%} | {astral['identification']['empirical_fdp']:.3%} |")
        datasets.append(("Astral reference", "#dd6b20", astral))
    # Final exports support only a diagnostic check, not full FDR calibration.
    curve_labels = [f"{point['q_threshold']:.2%}" for point in metrics["fdr_curve"]]
    series = [("q limit", "#999999", [100 * point["q_threshold"] for point in metrics["fdr_curve"]])]
    for name, color, data in datasets:
        if [point["q_threshold"] for point in data["fdr_curve"]] != [point["q_threshold"] for point in metrics["fdr_curve"]]:
            raise SystemExit("FDR curves must use the same q-value limits")
        series.append((name, color, [100 * point["entrapment_share"] if point["entrapment_share"] is not None else None for point in data["fdr_curve"]]))
    write_chart(output / "fdr.svg", "Entrapment share at each q-value limit", "Percent", curve_labels, series, decimals=3)
    lines.extend([
        "",
        "![FDR check](fdr.svg)",
        "",
        "The q-value is the reported error estimate. The bars show the fraction of accepted rows from absent sequences.",
        "FDR means false discovery rate: the fraction of accepted results that are wrong.",
        "An n/a label means that no rows passed the limit.",
        f"TimSim uses `{metrics['entrapment_q_column']}`. " + (f"Astral uses `{astral['q_column']}`." if astral else ""),
        "These shares are not adjusted for search-space size. They do not prove that FDR is calibrated. Zero observed errors does not prove zero error risk.",
        "The input files contain only exported results. Rejected results are not available for a full calibration test.",
    ])
    return lines


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
    write_chart(output / "runtime.svg", "OpenDIA run time", "Seconds", labels, comparison(old_time, new_time))
    write_chart(output / "peak_memory.svg", "OpenDIA peak memory", "MiB", labels, comparison(old_memory, new_memory))

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
    lines.extend(["Both runners use `batchSize 0` and `cacheWorkingInMemory`. The scheduler selects its batch sizes and concurrency automatically.", "The shared settings are not proof of the fastest possible run."])
    lines.extend(
        [
            f"Each TimSim value is from one run. The previous runs used {previous['threads']} threads. The current runs used {current['threads']} threads.",
            "Run times can change with host load. Repeat these tests on an idle host before drawing a speed conclusion.",
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

    astral = None
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
                "The new shared settings have not been tested on the full Astral data.",
            ]
        )
        host_memory = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 1024**3
        if astral_resources["memory_mib"] / 1024 > host_memory:
            lines.append(f"The reference run used more memory than this host has ({host_memory:.1f} GiB).")

    lines.extend(["", *quality_report(output, metrics, astral)])
    (output / "FINAL_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
