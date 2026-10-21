#!/usr/bin/env python3
"""Synthetic report regression tests. Never regenerate benchmark artifacts."""

import copy
import json
import math
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from make_report import (
    identification_report,
    quantification_report,
    quantile,
    read_job,
    read_resources,
    short_run_name,
    workload_cell,
    write_distribution,
)


SVG = "{http://www.w3.org/2000/svg}"
SCRIPT = Path(__file__).with_name("make_report.py")


def metrics(astral=False):
    conditions = {f"run{i}": "repeat" if astral else "control" if i < 3 else "treatment"
                  for i in range(3 if astral else 6)}
    groups = sorted(set(conditions.values()))
    summary = {group: {"precursor": {"values": [8, 10, 12], "n": 3, "mean": 10, "sample_sd": 2}}
               for group in groups}
    for levels in summary.values():
        levels["peptide"] = {"values": [4, 6, 8], "n": 3, "mean": 6, "sample_sd": 2}
        levels["protein_label"] = {"values": [1, 2, 3], "n": 3, "mean": 2, "sample_sd": 1}
    quant = {"precursor": {"source": "native_matrix", "status": "ok", "conditions": {
        group: {"entities": 4, "complete_entities": 3, "completeness": {"0": 1, "3": 3},
                "median_cv_complete": .1, "cv_n": 3, "cv_values": [.05, .1, .2]}
        for group in groups}}}
    quant["protein_label"] = copy.deepcopy(quant["precursor"])
    quant["peptide"] = copy.deepcopy(quant["precursor"])
    views = {"precursor_only": {"status": "ok", "count_summary": summary},
             "standard": {"status": "ok", "count_summary": summary, "quantification": quant}}
    result = {"schema_version": 2, "runs": list(conditions), "conditions": conditions,
              "views": views, "fdr_curve_per_run": {
                  run: [{"q_threshold": .001, "rows": 0, "entrapment_rows": 0, "entrapment_share": None},
                        {"q_threshold": .01, "rows": 100, "entrapment_rows": i + 1, "entrapment_share": (i + 1) / 100}]
                  for i, run in enumerate(conditions)}}
    if not astral:
        result["entrapment_views"] = copy.deepcopy(views)
        result["target_identification"] = {"overall_recovery": 1}
        result["quantification"] = {"native_standard_effects": {"precursor": {
            "n": 2, "pearson": .9, "effects": [
                {"entity": "a", "realized_log2_fc": -1, "observed_log2_fc": -.75},
                {"entity": "b", "realized_log2_fc": 1, "observed_log2_fc": 1.5}]}}}
        effects = result["quantification"]["native_standard_effects"]
        effects["peptide"] = copy.deepcopy(effects["precursor"])
        effects["protein_label"] = copy.deepcopy(effects["precursor"])
    return result


def job(root, complete=True):
    (root / "meta").mkdir(parents=True)
    (root / "logs").mkdir()
    (root / "meta/command.sh").write_text("OpenDIA -threads 4\n")
    (root / "meta/source_revision.txt").write_text("abc123\n")
    (root / "meta/exit_code.txt").write_text("0\n")
    (root / "meta/source_status.txt").write_text(" M source.cpp\n")
    (root / "meta/OpenDIA.sha256").write_text("f" * 64 + "  /build/OpenDIA\n")
    (root / "meta/resource_limits.txt").write_text("memory_max=20G\nmemory_swap_max=6G\n")
    if complete:
        (root / "logs/resources.txt").write_text(
            "Elapsed (wall clock) time (h:mm:ss or m:ss): 1:14.09\n"
            "Maximum resident set size (kbytes): 2048\nExit status: 0\n")
    workload = {key: {"mean": 4, "min": 2, "max": 8} for key in (
        "host_cpu_busy_percent", "load1", "io_pressure_some_percent", "mem_available_mib", "swap_used_mib")}
    workload["other_cpu_busy_estimate_percent"] = {"mean": 999, "min": 999, "max": 999}
    (root / "logs/host_workload.json").write_text(json.dumps(workload))


class ReportTests(unittest.TestCase):
    def test_sample_sd_and_individual_dots(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            lines = identification_report(output, [("TimSim target", metrics())])
            self.assertIn("| TimSim target / control | 3 | 10.0 | 2.0 | 8, 10, 12 |", lines)
            self.assertIn("| TimSim target / control | 3 | 10.0 ± 2.0 | 6.0 ± 2.0 | 2.0 ± 1.0 |", lines)
            tree = ET.parse(output / "identification_standard.svg")
            self.assertEqual(len(tree.findall(f"{SVG}circle")), 6)
            self.assertEqual(len(tree.findall(f"{SVG}path")), 2)
            self.assertFalse((output / "identification_precursor_only.svg").exists())
            self.assertIn("Precursor IDs (1% FDR)", (output / "identification_standard.svg").read_text())

    def test_zero_ids_and_na_standard_are_distinct(self):
        data = metrics(True)
        data["views"]["precursor_only"]["count_summary"]["repeat"]["precursor"]["values"] = [0, 0, 0]
        data["views"]["standard"] = {"status": "NA", "reason": "Global inference disabled", "count_summary": None}
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            text = "\n".join(identification_report(output, [("Astral historical", data)]))
            self.assertIn("| 3 | 0.0 | 0.0 | 0, 0, 0 |", text)
            self.assertIn("Global inference disabled", text)
            self.assertIn("not biological protein accessions", text)
            self.assertEqual(len(ET.parse(output / "identification_standard.svg").findall(f"{SVG}circle")), 0)

    def test_cv_distribution_and_native_truth(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            text = "\n".join(quantification_report(output, [("TimSim target", metrics()), ("Astral", metrics(True))]))
            self.assertIn("| 3 / 4 | 75.00% | 10.00% | 3 |", text)
            self.assertIn("| precursor / native standard matrix | 2 | 0.900 | 0.375 |", text)
            self.assertIn("FASTA group (peptide entries)", text)
            cv = (output / "quantification_cv_0.svg").read_text()
            self.assertIn("Box: IQR", cv)
            self.assertIn("not replicate runs", cv)
            self.assertNotIn("sample SD", cv)
            self.assertIn(" / peptide |", text)
            self.assertNotIn(">peptide</text>", cv)
            self.assertEqual(len(list(output.glob("truth_*.svg"))), 2)
            self.assertFalse((output / "truth_0_1.svg").exists())
            self.assertEqual(quantile([1, 2, 3], .25), 1.5)

    def test_reject_interrupted_or_aborted_jobs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "job"
            job(root)
            resources = root / "logs/resources.txt"
            original = resources.read_text()
            resources.write_text("Command terminated by signal 9\n" + original)
            with self.assertRaisesRegex(SystemExit, "signal"):
                read_resources(resources)
            resources.write_text(original)
            (root / "meta/exit_code.txt").write_text("137\n")
            with self.assertRaisesRegex(SystemExit, "successfully"):
                read_job(root)
            (root / "meta/exit_code.txt").write_text("0\n")
            workload_path = root / "logs/host_workload.json"
            workload = json.loads(workload_path.read_text())
            workload["abort_reason"] = "memory cap"
            workload_path.write_text(json.dumps(workload))
            with self.assertRaisesRegex(SystemExit, "memory cap"):
                read_job(root)
            workload_path.write_text("null")
            self.assertEqual(read_job(root)["workload"], {})
            # An unrelated ancestor manifest must not be assigned to this job.
            (root.parent / "build-manifest.txt").write_text("wrong historical provenance")
            self.assertNotIn("build_manifest", read_job(root))

    def test_short_run_names_and_workload_cells(self):
        self.assertEqual(short_run_name("OpenSwathTimSim_control_01"), "control01")
        self.assertEqual(short_run_name("OpenSwathTimSim_treatment_03"), "treatment03")
        self.assertEqual(short_run_name("LFQ_Astral_DIA_15min_50ng_Human_02.raw"), "Astral02")
        self.assertEqual(short_run_name("unrecognized.run"), "unrecognized.run")
        workload = {"mem_available_mib": {"min": 2048}, "swap_used_mib": {"max": None}}
        self.assertEqual(workload_cell(workload, "mem_available_mib", ("min",), 1024), "2.00")
        self.assertEqual(workload_cell(workload, "swap_used_mib", ("max",), 1024), "n/a")

    def test_plot_edge_cases(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "plot.svg"
            write_distribution(path, "A & B", "IDs", [("empty", []), ("one", [0]), ("vary", [0, 0, 10])])
            tree = ET.parse(path)
            self.assertEqual(len(tree.findall(f"{SVG}circle")), 4)
            self.assertEqual(len(tree.findall(f"{SVG}path")), 1)
            for circle in tree.findall(f"{SVG}circle"):
                self.assertTrue(math.isfinite(float(circle.attrib["cy"])))
                self.assertGreaterEqual(float(circle.attrib["cy"]), 75)
                self.assertLessEqual(float(circle.attrib["cy"]), 305)

    def test_cli_completed_and_unfinished_jobs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ("current", "previous"):
                for search in ("target_only", "entrapment"):
                    job(root / name / search)
                (root / name / "validation.json").write_text(json.dumps(metrics()))
            job(root / "astral")
            (root / "astral/metrics.json").write_text(json.dumps(metrics(True)))
            (root / "build-manifest.txt").write_text("source_dirty=true\nCMAKE_BUILD_TYPE:STRING=Release\n")
            output = root / "output"
            command = [sys.executable, str(SCRIPT), str(root / "current"), "--previous-timsim", str(root / "previous"), "--astral", str(root / "astral"), "--output-dir", str(output)]
            completed = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            text = (output / "FINAL_REPORT.md").read_text()
            for filename in ("timsim_runtime.svg", "timsim_peak_memory.svg", "astral_runtime.svg", "astral_peak_memory.svg"):
                tree = ET.parse(output / filename)
                self.assertFalse(tree.findall(f"{SVG}circle"))
                self.assertFalse(tree.findall(f"{SVG}path"))
            self.assertIn("The two TimSim searches are not repeats", text)
            self.assertIn("These times do not measure an isolated code speed gain", text)
            self.assertIn("not search-space adjusted or calibrated FDR", text)
            self.assertIn("20G / 6G", text)
            self.assertNotIn("CMAKE_BUILD_TYPE:STRING=Release", text)
            self.assertIn("f" * 12, text)
            self.assertNotIn("f" * 64, text)
            self.assertNotIn("999.00", text)
            self.assertIn("| TimSim mixed / control |", text)
            self.assertIn("| TimSim mixed | run0 | 1.00% | 100 | 1 | 1.00% |", text)
            self.assertNotIn("| TimSim mixed | run0 | 0.10% |", text)
            entrapment = ET.parse(output / "entrapment_0.svg")
            self.assertEqual(len(entrapment.findall(f"{SVG}circle")), 6)
            self.assertEqual(len(entrapment.findall(f"{SVG}path")), 1)
            self.assertTrue((output / "entrapment_1.svg").exists())
            self.assertEqual(text.count("!["), 11)
            self.assertFalse(list(output.glob("quantification_completeness*.svg")))
            self.assertFalse((output / "quantification_cv_1.svg").exists())
            self.assertIn("not replicate uncertainty", text)
            for path in output.glob("*.svg"):
                ET.parse(path)
            historical = metrics(True)
            historical["views"]["standard"] = {"status": "NA", "reason": "Global inference disabled"}
            historical_metrics = root / "historical.json"
            historical_metrics.write_text(json.dumps(historical))
            legacy_output = root / "legacy"
            legacy = subprocess.run([
                sys.executable, str(SCRIPT), str(root / "current"),
                "--previous-timsim", str(root / "previous"),
                "--astral-metrics", str(historical_metrics),
                "--astral-time", str(root / "astral/logs/resources.txt"),
                "--output-dir", str(legacy_output),
            ], capture_output=True, text=True)
            self.assertEqual(legacy.returncode, 0, legacy.stderr)
            legacy_text = (legacy_output / "FINAL_REPORT.md").read_text()
            self.assertIn("Global inference disabled", legacy_text)
            self.assertIn("| not recorded | not recorded | not recorded |", legacy_text)
            # A first report needs no previous baseline. Null monitoring is n/a.
            (root / "current/target_only/logs/host_workload.json").write_text(json.dumps({"load1": None}))
            first_output = root / "first"
            first = subprocess.run([
                sys.executable, str(SCRIPT), str(root / "current"),
                "--output-dir", str(first_output),
            ], capture_output=True, text=True)
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertNotIn("| Previous target |", (first_output / "FINAL_REPORT.md").read_text())
            # Incomplete full Astral must fail before creating any report output.
            (root / "astral/logs/resources.txt").unlink()
            blocked = root / "blocked"
            command[-1] = str(blocked)
            failed = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertFalse(blocked.exists())


if __name__ == "__main__":
    unittest.main()
