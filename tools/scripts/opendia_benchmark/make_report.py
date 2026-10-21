#!/usr/bin/env python3
"""Make a research report from completed OpenDIA benchmark measurements.

This command writes Markdown and native SVG only. It does not run a benchmark.
"""

import argparse
import html
import json
import math
import re
import shlex
import statistics
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
    if "Command terminated by signal" in text:
        raise SystemExit(f"Benchmark was terminated by a signal in {path}")
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


def format_duration(seconds):
    hours, remainder = divmod(round(seconds), 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours} h {minutes} min {seconds} s"
    if minutes:
        return f"{minutes} min {seconds} s"
    return f"{seconds} s"


def write_chart(path, title, unit, labels, series, decimals=1):
    """Write grouped bars. None denotes unavailable, never zero."""
    width, height = max(900, 160 * len(labels)), 440
    left, top, bottom = 90, 85, 90
    chart_height = height - top - bottom
    finite_values = [v for _, _, values in series for v in values if v is not None and math.isfinite(v)]
    maximum = max([0, *finite_values]) * 1.18 or 1
    minimum = min([0, *finite_values]) * 1.18
    scale = chart_height / (maximum - minimum)
    baseline = top + maximum * scale
    group_width = (width - left - 30) / max(len(labels), 1)
    bar_width = min(58, group_width / (len(series) + 1))
    parts = svg_start(width, height, title)
    if not labels:
        parts.append(svg_text(width / 2, height / 2, "n/a: no available measurements"))
    parts.append(f'<line x1="{left}" y1="{baseline}" x2="{width - 20}" y2="{baseline}" stroke="#333"/>')
    for index, label in enumerate(labels):
        center = left + group_width * (index + 0.5)
        for series_index, (_, color, values) in enumerate(series):
            offset = (series_index - (len(series) - 1) / 2) * bar_width * 1.15
            value = values[index]
            x = center + offset - bar_width / 2
            if value is None or not math.isfinite(value):
                parts.append(svg_text(x + bar_width / 2, baseline - 7, "n/a", 12))
                continue
            bar_height = abs(value) * scale
            y = baseline - max(value, 0) * scale
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width}" height="{bar_height:.1f}" fill="{color}"/>')
            parts.append(svg_text(x + bar_width / 2, y - 7, f"{value:.{decimals}f}", 12))
        parts.append(svg_label(center, height - 65, label))
    parts.append(f'<text x="20" y="{top + chart_height / 2}" transform="rotate(-90 20 {top + chart_height / 2})" text-anchor="middle" font-family="sans-serif" font-size="13">{html.escape(unit)}</text>')
    for index, (name, color, _) in enumerate(series):
        x = 100 + index * 250
        parts.append(f'<rect x="{x}" y="48" width="14" height="14" fill="{color}"/>')
        parts.append(f'<text x="{x + 20}" y="60" font-family="sans-serif" font-size="12">{html.escape(name)}</text>')
    save_svg(path, parts)


def svg_text(x, y, text, size=13):
    return f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="middle" font-family="sans-serif" font-size="{size}">{html.escape(str(text))}</text>'


def svg_label(x, y, text):
    return "\n".join(svg_text(x, y + index * 15, part.strip(), 11)
                     for index, part in enumerate(text.split(" / ")))


def svg_start(width, height, title):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(title, quote=True)}">', '<rect width="100%" height="100%" fill="white"/>', svg_text(width / 2, 28, title, 18)]


def save_svg(path, parts):
    path.write_text("\n".join([*parts, "</svg>"]) + "\n", encoding="utf-8")


def write_distribution(path, title, unit, groups):
    """Plot individual observations, with mean and sample SD for n >= 2."""
    width, height = max(900, 160 * len(groups)), 410
    left, top, bottom = 90, 75, 105
    summaries = []
    for label, values in groups:
        values = [v for v in values if v is not None and math.isfinite(v)]
        mean = statistics.mean(values) if values else None
        sd = statistics.stdev(values) if len(values) > 1 else None
        summaries.append((label, values, mean, sd))
    upper = [max(values + [mean + (sd or 0)]) for _, values, mean, sd in summaries if values]
    lower = [mean - (sd or 0) for _, values, mean, sd in summaries if values]
    maximum = max([0, *upper]) * 1.15 or 1
    minimum = min([0, *lower]) * 1.15
    scale = (height - top - bottom) / (maximum - minimum)
    base = top + maximum * scale
    parts = svg_start(width, height, title)
    parts.append(svg_text(width / 2, 51, "Dots: individual observations. Crossbar: mean. Whiskers: sample SD.", 12))
    if not groups:
        parts.append(svg_text(width / 2, height / 2, "n/a: no available observations"))
    for tick in range(5):
        value = minimum + (maximum - minimum) * tick / 4
        y = base - value * scale
        parts.append(f'<line x1="{left}" y1="{y}" x2="{width - 25}" y2="{y}" stroke="#e5e7eb"/>')
        parts.append(svg_text(left - 35, y + 4, f"{value:.3g}", 11))
    for index, (label, values, mean, sd) in enumerate(summaries):
        x = left + (width - left - 30) * (index + 0.5) / max(len(groups), 1)
        if mean is None:
            parts.append(svg_text(x, base - 10, "n/a"))
        else:
            y = base - mean * scale
            parts.append(f'<line x1="{x - 20}" y1="{y}" x2="{x + 20}" y2="{y}" stroke="#111" stroke-width="3"/>')
            if sd is not None:
                low, high = base - (mean - sd) * scale, base - (mean + sd) * scale
                parts.append(f'<path d="M {x} {low} V {high} M {x - 8} {low} H {x + 8} M {x - 8} {high} H {x + 8}" fill="none" stroke="#111"/>')
            for dot, value in enumerate(values):
                jitter = (dot - (len(values) - 1) / 2) * min(9, 60 / max(len(values), 1))
                parts.append(f'<circle cx="{x + jitter}" cy="{base - value * scale}" r="4" fill="#2b6cb0" fill-opacity="0.65"><title>{html.escape(label)}: {value:g}</title></circle>')
        parts.append(svg_label(x, height - bottom + 23, label))
        parts.append(svg_text(x, height - bottom + 83, f"n={len(values)}", 12))
    parts.append(f'<text x="18" y="190" transform="rotate(-90 18 190)" text-anchor="middle" font-family="sans-serif" font-size="13">{html.escape(unit)}</text>')
    save_svg(path, parts)


def quantile(values, fraction):
    values = sorted(values)
    index = (len(values) - 1) * fraction
    lower = int(index)
    upper = min(lower + 1, len(values) - 1)
    return values[lower] + (values[upper] - values[lower]) * (index - lower)


def write_cv_chart(path, groups, title="Within-condition entity CV distribution"):
    """Entity CVs are a distribution, not independent biological replicates."""
    width, height = max(900, 150 * len(groups)), 430
    summaries = []
    for label, values in groups:
        values = [v for v in values if v is not None and math.isfinite(v)]
        summaries.append((label, values, [quantile(values, q) for q in (.1, .25, .5, .75, .9)] if values else []))
    maximum = max([1, *[points[-1] for _, _, points in summaries if points]]) * 1.15
    parts = svg_start(width, height, title)
    parts.append(svg_text(width / 2, 50, "Line: median. Box: IQR. Whiskers: 10th to 90th percentiles. n: entities, not replicate runs.", 12))
    if not groups:
        parts.append(svg_text(width / 2, height / 2, "n/a: no available CV distribution"))
    base, top = 310, 80
    def y(value):
        return base - value / maximum * (base - top)
    for tick in range(5):
        value = maximum * tick / 4
        parts.append(f'<line x1="80" y1="{y(value)}" x2="{width - 20}" y2="{y(value)}" stroke="#e5e7eb"/>')
        parts.append(svg_text(45, y(value) + 4, f"{value:.2g}%", 11))
    for index, (label, values, points) in enumerate(summaries):
        x = 80 + (width - 100) * (index + .5) / max(len(groups), 1)
        if points:
            p10, p25, median, p75, p90 = points
            parts.append(f'<path d="M {x} {y(p10)} V {y(p90)} M {x - 10} {y(p10)} H {x + 10} M {x - 10} {y(p90)} H {x + 10}" fill="none" stroke="#333"/>')
            parts.append(f'<rect x="{x - 23}" y="{y(p75)}" width="46" height="{y(p25) - y(p75)}" fill="#c3d9ef" stroke="#2b6cb0"/>')
            parts.append(f'<line x1="{x - 23}" y1="{y(median)}" x2="{x + 23}" y2="{y(median)}" stroke="#111" stroke-width="2"/>')
        else:
            parts.append(svg_text(x, base - 10, "n/a"))
        parts.append(svg_label(x, base + 23, label))
        parts.append(svg_text(x, base + 85, f"n={len(values)}", 11))
    save_svg(path, parts)


def write_truth_chart(path, title, effects):
    pairs = [(p.get("realized_log2_fc"), p.get("observed_log2_fc")) for p in effects]
    pairs = [(x, y) for x, y in pairs if x is not None and y is not None and math.isfinite(x) and math.isfinite(y)]
    parts = svg_start(650, 470, title)
    if not pairs:
        parts.append(svg_text(325, 230, "n/a: no matched truth effects"))
    else:
        low = min(0, *[min(pair) for pair in pairs]) - .2
        high = max(0, *[max(pair) for pair in pairs]) + .2
        def x(value):
            return 85 + (value - low) / (high - low) * 500
        def y(value):
            return 390 - (value - low) / (high - low) * 320
        parts.append(f'<path d="M 85 70 V 390 H 585 M {x(low)} {y(low)} L {x(high)} {y(high)}" fill="none" stroke="#999"/>')
        for truth, observed in pairs:
            parts.append(f'<circle cx="{x(truth):.2f}" cy="{y(observed):.2f}" r="2" fill="#2b6cb0" fill-opacity="0.25"/>')
        for tick in range(5):
            value = low + (high - low) * tick / 4
            parts.append(svg_text(x(value), 411, f"{value:.2g}", 11))
            parts.append(svg_text(55, y(value) + 4, f"{value:.2g}", 11))
    parts.extend([svg_text(325, 450, "Realized log2 fold change"), '<text x="20" y="235" transform="rotate(-90 20 235)" text-anchor="middle" font-family="sans-serif" font-size="13">Observed log2 fold change</text>', svg_text(325, 51, f"Matched entities: {len(pairs)}. Grey diagonal: exact agreement.", 12)])
    save_svg(path, parts)


def optional_text(path):
    return path.read_text(encoding="utf-8").strip() if path.exists() else "not recorded"


def read_job(root):
    meta = root / "meta"
    command_text = optional_text(meta / "command.sh")
    command = shlex.split(command_text)
    exit_code = meta / "exit_code.txt"
    if exit_code.exists() and exit_code.read_text().strip() != "0":
        raise SystemExit(f"Benchmark did not finish successfully in {root}: {exit_code.read_text().strip()}")
    workload_path = root / "logs/host_workload.json"
    workload = (load_json(workload_path) or {}) if workload_path.exists() else {}
    if workload.get("abort_reason"):
        raise SystemExit(f"Benchmark aborted in {root}: {workload['abort_reason']}")
    return {
        "root": root,
        "resource": read_resources(root / "logs/resources.txt"),
        "workload": workload,
        "revision": optional_text(meta / "source_revision.txt"),
        "binary": optional_text(meta / "OpenDIA.sha256"),
        "source_status": optional_text(meta / "source_status.txt"),
        "limits": optional_text(meta / "resource_limits.txt"),
        "command": command_text,
        "threads": command[command.index("-threads") + 1] if "-threads" in command else "not recorded",
    }


def read_timsim(root):
    jobs = {name: read_job(root / name) for name in ("target_only", "entrapment")}
    return {"metrics": load_json(root / "validation.json"), "jobs": jobs}


def number_text(value, digits=3):
    return f"{value:.{digits}f}" if value is not None and math.isfinite(value) else "n/a"


def percent_text(value):
    return f"{value:.2%}" if value is not None and math.isfinite(value) else "n/a"


def mean_sd_text(values):
    if not values:
        return "n/a"
    mean = statistics.mean(values)
    sd = statistics.stdev(values) if len(values) > 1 else None
    return f"{number_text(mean, 1)} ± {number_text(sd, 1)}"


def identification_report(output, datasets):
    lines = ["## Identification", "", "Counts are distinct target precursor IDs, not rows. TimSim has three control and three treatment runs. Astral has three repeat runs.", "Precursor-only uses m_score <= 1%. Standard also requires global peptide and protein q <= 1%.", "Dots show individual runs. Error bars show sample SD, the run-to-run spread.", ""]
    for view_name, title in (("precursor_only", "Precursor-only 1% FDR"), ("standard", "Standard all-three-filter 1% view")):
        header = "| Search / condition | n | Precursor mean ± SD | Peptide mean ± SD | Protein / FASTA groups mean ± SD |" if view_name == "standard" else "| Search / condition | Runs | Mean precursor IDs | Sample SD | Individual run counts |"
        lines.extend([f"### {title}", "", header, "|---|---:|---:|---:|---|"])
        groups = []
        for name, metrics in datasets:
            view = metrics.get("views", {}).get(view_name, {})
            if view.get("status") != "ok":
                lines.append(f"| {name} | n/a | n/a | n/a | {view.get('reason', 'Metrics unavailable')} |")
                continue
            for condition, levels in view.get("count_summary", {}).items():
                values = levels.get("precursor", {}).get("values", [])
                label = f"{name} / {condition}"
                groups.append((label, values))
                mean = statistics.mean(values) if values else None
                sd = statistics.stdev(values) if len(values) > 1 else None
                if view_name == "standard":
                    summaries = [mean_sd_text(levels.get(level, {}).get("values", [])) for level in ("precursor", "peptide", "protein_label")]
                    lines.append(f"| {label} | {len(values)} | " + " | ".join(summaries) + " |")
                else:
                    lines.append(f"| {label} | {len(values)} | {number_text(mean, 1)} | {number_text(sd, 1)} | {', '.join(str(v) for v in values)} |")
        if view_name == "standard":
            filename = "identification_standard.svg"
            write_distribution(output / filename, title, "Precursor IDs (1% FDR)", groups)
            lines.extend(["", f"![{title}]({filename})", ""])
        else:
            lines.append("")
    lines.extend(["Astral FASTA entries are peptide pseudo-proteins, not biological protein accessions. Its group counts refer to these peptide entries.", "Historical Astral with global inference disabled has no standard view, not zero IDs.", ""])
    return lines


def quantification_report(output, datasets):
    lines = ["## Quantitative variation and completeness", "", "CV is sample SD / mean abundance for entities present in all three repeat runs, using standard-filter native matrices.", "Completeness is the fraction present in all three runs, against the target union across all conditions (including condition-missing entities).", "Mixed-search detail remains in metrics JSON.", "", "| Search / condition / level | Complete / union entities | Completeness | Median CV | CV entities |", "|---|---:|---:|---:|---:|"]
    quant_plots = []
    for dataset_index, (name, metrics) in enumerate(datasets):
        if name == "TimSim mixed":
            continue
        cv_groups = []
        standard = metrics.get("views", {}).get("standard", {})
        quantification = standard.get("quantification", {})
        if standard.get("status") != "ok" or not quantification:
            lines.append(f"| {name} | n/a | n/a | n/a | n/a |")
            continue
        for level, data in quantification.items():
            if not isinstance(data, dict):
                continue
            if data.get("status") == "NA" or not data.get("conditions"):
                lines.append(f"| {name} / {level}: {data.get('reason', 'Native matrix unavailable')} | n/a | n/a | n/a | n/a |")
                continue
            for condition, group in data.get("conditions", {}).items():
                level_label = "FASTA group (peptide entries)" if level == "protein_label" and name.startswith("Astral") else "protein-inference group" if level == "protein_label" else level
                label = f"{name} / {condition} / {level_label}"
                entities = group.get('entities', 0)
                complete = group.get('complete_entities', 0)
                fraction = complete / entities if entities else None
                lines.append(f"| {label} | {complete} / {entities} | {percent_text(fraction)} | {percent_text(group.get('median_cv_complete'))} | {group.get('cv_n', 0)} |")
                if level == "peptide":
                    continue
                plot_label = f"{condition} / {level_label}"
                cv_groups.append((plot_label, [100 * v for v in group.get("cv_values", [])]))
        cv_filename = f"quantification_cv_{dataset_index}.svg"
        write_cv_chart(output / cv_filename, cv_groups, f"{name}: within-condition entity CV")
        quant_plots.extend(["", f"![{name} CV]({cv_filename})"])
    lines.extend(["", "Boxes show median and IQR (the middle half of entities). Whiskers span the 10th to 90th percentiles. Entity counts are not replicate counts.", *quant_plots, "", "## TimSim truth accuracy", "", "Observed treatment/control log2 fold changes are compared with known changes. Correlation alone does not prove exact agreement.", "Native standard matrices are preferred. Export-based protein sums are descriptive only.", "", "| TimSim search | Level / source | Matched effects | Truth correlation | Median absolute log2 error |", "|---|---|---:|---:|---:|"])
    truth_plots = []
    for dataset_index, (name, metrics) in enumerate(datasets):
        if "target_identification" not in metrics:
            continue
        truth = metrics.get("quantification", {})
        native = truth.get("native_standard_effects", {})
        if native:
            effect_groups = [(level, "native standard matrix", data.get("effects", []), data.get("pearson")) for level, data in native.items()]
        else:
            effect_groups = [("precursor", "export", truth.get("precursor_effect_values", []), truth.get("precursor_realized_vs_observed_pearson")), ("protein descriptive sum", "export", truth.get("protein_descriptive_sum_effect_values", []), truth.get("protein_realized_vs_observed_pearson"))]
        for level_index, (level, source, effects, correlation) in enumerate(effect_groups):
            errors = [abs(p["observed_log2_fc"] - p["realized_log2_fc"]) for p in effects if p.get("observed_log2_fc") is not None and p.get("realized_log2_fc") is not None and math.isfinite(p["observed_log2_fc"]) and math.isfinite(p["realized_log2_fc"])]
            lines.append(f"| {name} | {level} / {source} | {len(errors)} | {number_text(correlation)} | {number_text(statistics.median(errors) if errors else None)} |")
            if level == "peptide":
                continue
            filename = f"truth_{dataset_index}_{level_index}.svg"
            write_truth_chart(output / filename, f"{name}: {level} truth accuracy", effects)
            truth_plots.extend(["", f"![{name} {level} truth accuracy]({filename})"])
    lines.extend(truth_plots)
    lines.extend(["", "Astral has no quantitative ground truth here. Its three-run CV and completeness measure repeatability only.", ""])
    return lines


def short_run_name(run):
    stem = Path(run).stem
    timsim = re.fullmatch(r"OpenSwathTimSim_(control|treatment)_(\d+)", stem)
    if timsim:
        return f"{timsim[1]}{int(timsim[2]):02d}"
    astral = re.fullmatch(r"LFQ_Astral_DIA_15min_50ng_Human_(\d+)", stem)
    return f"Astral{int(astral[1]):02d}" if astral else run


def fdr_report(output, datasets):
    lines = ["## Per-run entrapment diagnostic", "", "Tables use m_score <= 1%. Plots and JSON retain all cutoffs. Full run names remain in JSON.", "Raw share counts absent entrapment sequences among accepted non-decoy precursor IDs. It is not search-space adjusted or calibrated FDR.", "Dots are runs. Mean and sample SD show descriptive spread across the study, not replicate uncertainty. TimSim combines both conditions here.", "No accepted IDs means n/a, not zero share. Zero entrapment does not prove zero error risk. Export limits bound the curves."]
    for dataset_index, (name, metrics) in enumerate(datasets):
        lines.extend(["", f"### {name}", "", "| Search | Run | m_score cutoff | Accepted | Entrapment | Raw share |", "|---|---|---:|---:|---:|---:|"])
        curves = metrics.get("fdr_curve_per_run", {})
        if not curves:
            lines.append(f"| {name} | n/a | n/a | n/a | n/a | n/a |")
            continue
        labels = sorted({point["q_threshold"] for points in curves.values() for point in points})
        groups = []
        for cutoff in labels:
            values = [100 * p["entrapment_share"] for points in curves.values() for p in points if p["q_threshold"] == cutoff and p["entrapment_share"] is not None]
            groups.append((f"m_score <= {cutoff:.2%}", values))
        filename = f"entrapment_{dataset_index}.svg"
        write_distribution(output / filename, f"{name}: raw entrapment share (descriptive across {len(curves)} runs)", "Percent", groups)
        for run, points in curves.items():
            for point in points:
                if math.isclose(point["q_threshold"], .01):
                    lines.append(f"| {name} | {short_run_name(run)} | {point['q_threshold']:.2%} | {point['rows']} | {point['entrapment_rows']} | {percent_text(point['entrapment_share'])} |")
        lines.extend(["", f"![{name} entrapment]({filename})", ""])
    return lines


def workload_cell(workload, field, keys, divisor=1):
    data = workload.get(field) or {}
    values = [data.get(key) for key in keys]
    if not all(isinstance(value, (int, float)) and math.isfinite(value) for value in values):
        return "n/a"
    values = [f"{value / divisor:.2f}" for value in values]
    return f"{values[0]} [{values[1]}, {values[2]}]" if len(values) == 3 else " / ".join(values)


def resource_report(output, timsim, previous, astral_jobs):
    lines = ["## Runtime and peak memory", "", "Each search uses all input files. We measured each search once, so time and memory have no error bars. The two TimSim searches are not repeats.", "Peak RSS is peak process RAM from GNU time. It does not add RAM from simultaneous child processes or file caches.", "Old filters, memory limits and host load can differ. These times do not measure an isolated code speed gain.", ""]
    timsim_jobs = []
    for search, label in (("target_only", "target"), ("entrapment", "mixed")):
        if previous:
            timsim_jobs.append((f"Previous {label}", previous["jobs"][search]))
        timsim_jobs.append((f"Current {label}", timsim["jobs"][search]))
    datasets = [("timsim", "TimSim diaPASEF", timsim_jobs)]
    if astral_jobs:
        datasets.append(("astral", "Astral", astral_jobs))
    all_jobs = []
    for stem, title, jobs in datasets:
        all_jobs.extend((f"{title} / {label}", job) for label, job in jobs)
        lines.extend([f"### {title}", "", "| Search | Elapsed time | Peak RSS |", "|---|---:|---:|"])
        for label, job in jobs:
            resource = job["resource"]
            lines.append(f"| {label} | {format_duration(resource['seconds'])} | {resource['memory_mib'] / 1024:.2f} GiB |")
        for key, unit, divisor, suffix in (("seconds", "Seconds", 1, "runtime"), ("memory_mib", "GiB", 1024, "peak_memory")):
            filename = f"{stem}_{suffix}.svg"
            write_chart(output / filename, f"{title}: {suffix.replace('_', ' ')} (n=1)", unit, [label for label, _ in jobs], [("Measured search", "#2b6cb0", [job["resource"][key] / divisor for _, job in jobs])])
            lines.extend(["", f"![{title} {suffix}]({filename})", ""])
    lines.extend(["### Host workload", "", "Host measurements include the benchmark. We did not adjust times for load. We do not use the old estimate of other-process CPU use. Missing data is n/a.", "", "| Search | CPU % mean [range] | Load mean/max | I/O pressure % mean/max | Min available GiB | Max swap GiB |", "|---|---:|---:|---:|---:|---:|"])
    for label, job in all_jobs:
        cells = [workload_cell(job["workload"], field, keys, divisor) for field, keys, divisor in (
            ("host_cpu_busy_percent", ("mean", "min", "max"), 1),
            ("load1", ("mean", "max"), 1),
            ("io_pressure_some_percent", ("mean", "max"), 1),
            ("mem_available_mib", ("min",), 1024),
            ("swap_used_mib", ("max",), 1024),
        )]
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    lines.extend(["", "### Provenance and limits", "", "Binary SHA-256 hashes are shortened. Dirty means uncommitted changes. Full provenance and commands remain in job metadata.", "", "| Search | Source revision | Binary hash | Threads | Memory / swap limits |", "|---|---|---|---:|---|"])
    for label, job in all_jobs:
        status = job["source_status"]
        revision = job["revision"][:12] if job["revision"] != "not recorded" else "not recorded"
        if status and status != "not recorded":
            revision += " (dirty)"
        binary_hash = job["binary"].split()[0][:12] if job["binary"] != "not recorded" else "not recorded"
        limits = dict(line.split("=", 1) for line in job["limits"].splitlines() if "=" in line)
        limits_text = f"{limits.get('memory_max', 'n/a')} / {limits.get('memory_swap_max', 'n/a')}"
        lines.append(f"| {label} | {revision} | {binary_hash} | {job['threads']} | {limits_text} |")
    lines.append("")
    for stem, title, jobs in datasets:
        paths = dict.fromkeys(str(job["root"].parent if stem == "timsim" else job["root"]) for _, job in jobs)
        lines.append(f"{title} data: " + ", ".join(f"`{path}`" for path in paths) + ".")
    lines.append("")
    return lines


def find_metrics(root):
    for name in ("metrics.json", "validation.json"):
        if (root / name).exists():
            return root / name
    raise SystemExit(f"No metrics.json or validation.json in {root}. Generate metrics after the full run completes.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("current_timsim", type=Path)
    parser.add_argument("--previous-timsim", type=Path, help="Optional preserved TimSim baseline")
    parser.add_argument("--develop-revision")
    parser.add_argument("--astral", type=Path, help="Completed current Astral job directory")
    parser.add_argument("--previous-astral", type=Path, help="Optional completed historical Astral job directory")
    parser.add_argument("--astral-metrics", type=Path, help="Explicit Astral metrics JSON (also supports historical jobs)")
    parser.add_argument("--astral-time", type=Path, help="Historical Astral GNU time resource log")
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    if args.astral and args.astral_time:
        parser.error("Use --astral for the current job or --astral-time for a historical job, not both")
    if args.astral_time and not args.astral_metrics:
        parser.error("--astral-time requires --astral-metrics")
    if args.astral_metrics and not (args.astral or args.astral_time):
        parser.error("--astral-metrics requires --astral or --astral-time")
    current = read_timsim(args.current_timsim)
    previous = read_timsim(args.previous_timsim) if args.previous_timsim else None
    if args.develop_revision:
        result = subprocess.run(["git", "-C", Path(__file__).resolve().parents[3], "merge-base", "--is-ancestor", args.develop_revision, current["jobs"]["target_only"]["revision"]], check=False)
        if result.returncode:
            raise SystemExit("The develop revision is not in the benchmark source")
    datasets = [("TimSim target", current["metrics"])]
    mixed = {"views": current["metrics"].get("entrapment_views", {})}
    datasets.append(("TimSim mixed", mixed))
    # The metrics file describes target-only IDs and a separate mixed-search diagnostic.
    astral_jobs = []
    if args.previous_astral:
        astral_jobs.append(("Historical", read_job(args.previous_astral)))
    if args.astral:
        astral_jobs.append(("Current", read_job(args.astral)))
        datasets.append(("Astral", load_json(args.astral_metrics or find_metrics(args.astral))))
    elif args.astral_metrics:
        resource = read_resources(args.astral_time)
        job = {"root": args.astral_time.parent, "resource": resource, "workload": {}, "revision": "not recorded", "binary": "not recorded", "threads": "not recorded", "limits": "not recorded", "source_status": "not recorded", "command": "not recorded"}
        astral_jobs.append(("Historical", job))
        datasets.append(("Astral historical", load_json(args.astral_metrics)))
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    lines = ["# OpenDIA research benchmark report", "", "TimSim has six simulated diaPASEF runs. Astral has three experimental repeat runs.", "The report separates identification, abundance accuracy, repeat variation, time and memory.", ""]
    lines.extend(resource_report(output, current, previous, astral_jobs))
    lines.extend(identification_report(output, datasets))
    lines.extend(quantification_report(output, datasets))
    fdr_datasets = [("TimSim mixed", current["metrics"]), *[(name, metrics) for name, metrics in datasets if name.startswith("Astral")]]
    lines.extend(fdr_report(output, fdr_datasets))
    (output / "FINAL_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
