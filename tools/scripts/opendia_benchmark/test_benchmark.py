#!/usr/bin/env python3
"""Small regression check. Optional arguments are captured meta/command.sh files."""

import csv
import json
import os
import shlex
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from make_report import parse_seconds, read_resources, write_chart
from validate_timsim import best_results, fdr_curve


def main():
    assert parse_seconds("1:27:54") == 5274
    assert parse_seconds("1:14.09") == 74.09
    points = fdr_curve([(0.0005, False), (0.001, True), (0.02, True)], 0.01)
    assert points[0]["entrapment_share"] is None
    assert points[1]["entrapment_share"] == 0
    assert points[2]["entrapment_share"] == 0.5
    assert points[-1]["rows"] == 2
    assert fdr_curve([], 0.01)[-1]["entrapment_share"] is None

    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        resources = root / "resources.txt"
        resources.write_text("Elapsed (wall clock) time (h:mm:ss or m:ss): 1:14.09\nMaximum resident set size (kbytes): 2048\nExit status: 0\n")
        assert read_resources(resources) == {"seconds": 74.09, "memory_mib": 2, "status": 0}
        resources.write_text(resources.read_text().replace("status: 0", "status: 1"))
        try:
            read_resources(resources)
        except SystemExit:
            pass
        else:
            raise AssertionError("A failed run was accepted")

        result = root / "results.tsv"
        columns = ["run_name", "Sequence", "Charge", "Intensity", "RT", "EXP_IM", "m_score", "FullPeptideName", "ProteinName"]
        with result.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(columns)
            for index, q in enumerate(("nan", "", "inf", "-0.1", "0.02", "0.001")):
                writer.writerow(["run", f"PEPTIDE{index}", 2, 10, 1, 1, q, f"PEPTIDE{index}", "TARGET"])
        rows, q_column = best_results(result, 0.01, ("m_score",))
        assert len(rows) == 1 and q_column == "m_score"
        metrics_script = Path(__file__).with_name("metrics_astral.py")
        completed = subprocess.run([sys.executable, metrics_script, result], capture_output=True, text=True, check=True)
        assert json.loads(completed.stdout)["identification"]["target_rows"] == 1
        rejected = subprocess.run([sys.executable, metrics_script, result, "--q-threshold", "nan"], capture_output=True, text=True)
        assert rejected.returncode == 2 and "between 0 and 1" in rejected.stderr

        # Chrome can exit zero without writing a PDF. Keep the old report intact.
        fake_browser = root / "google-chrome"
        fake_browser.write_text("#!/bin/sh\nexit 0\n")
        fake_browser.chmod(0o755)
        report = root / "report.md"
        report.write_text("# Report\n")
        pdf = report.with_suffix(".pdf")
        pdf.write_bytes(b"%PDF-old report")
        environment = dict(os.environ, PATH=f"{root}:{os.environ.get('PATH', '')}")
        rendered = subprocess.run([
            sys.executable, Path(__file__).with_name("render_report.py"), report,
        ], env=environment, capture_output=True, text=True)
        assert rendered.returncode != 0
        assert "did not create a PDF" in rendered.stderr
        assert pdf.read_bytes() == b"%PDF-old report"

        for values in ([0, 1, None], [-1, 0, 1], [None, None, None]):
            chart = root / "chart.svg"
            write_chart(chart, "A & B", "Value", ["a", "b", "c"], [("Data", "blue", values)])
            tree = ET.parse(chart)
            for rectangle in tree.findall("{http://www.w3.org/2000/svg}rect"):
                if not rectangle.get("height", "").endswith("%"):
                    assert float(rectangle.get("height", "0")) >= 0
            if None in values:
                assert "n/a" in chart.read_text()

    # Verify actual captured commands, not just script text.
    expected = {
        "readOptions": "cacheWorkingInMemory",
        "batchSize": "0",
        "innerBatchSize": "-1",
        "maxConcurrentSwaths": "-1",
        "outer_loop_threads": "-1",
    }
    analysis_expected = {
        "-Inference:peptide:global": "true",
        "-Inference:protein:global": "true",
        "-Export:results:max_rs_peakgroup_qvalue": "0.05",
        "-Export:results:use_peptide_scores": "false",
        "-Export:results:use_protein_scores": "false",
        "-Export:results:exclude_decoys": "false",
        "-Export:matrix:max_rs_peakgroup_qvalue": "0.01",
        "-Export:matrix:use_peptide_scores": "true",
        "-Export:matrix:use_protein_scores": "true",
        "-Export:matrix:max_global_peptide_qvalue": "0.01",
        "-Export:matrix:max_global_protein_qvalue": "0.01",
    }
    for level in ("peptide", "protein", "gene"):
        for context in ("run_specific", "experiment_wide"):
            analysis_expected[f"-Inference:{level}:{context}"] = "false"
    analysis_expected["-Inference:gene:global"] = "false"
    for argument in sys.argv[1:]:
        command = shlex.split(Path(argument).read_text())
        for option, value in analysis_expected.items():
            assert command.count(option) == 1, (argument, option)
            assert command[command.index(option) + 1] == value, (argument, option)
        for name, value in expected.items():
            option = f"-TargetedDataExtraction:{name}"
            assert command.count(option) == 1, (argument, option)
            assert command[command.index(option) + 1] == value, (argument, option)
        if os.environ.get("OPEN_DIA"):
            # Run the real parser, then stop at a missing library before costly I/O.
            # -write_ini is not a dry run: it writes defaults, not resolved values.
            with tempfile.TemporaryDirectory() as temporary:
                command[2] = str(Path(os.environ["OPEN_DIA"]).resolve())
                library = Path(command[command.index("-tr") + 1])
                command[command.index("-tr") + 1] = str(Path(temporary) / ("missing" + library.suffix))
                command.extend(["-debug", "2"])
                accepted = subprocess.run(command, capture_output=True, text=True, timeout=30,
                                          env=dict(os.environ, OPENMS_DISABLE_UPDATE_CHECK="ON"))
                log = accepted.stdout + accepted.stderr
                assert accepted.returncode != 0 and "Error: File not found" in log
                assert "Parameters passed to" not in log and "Unknown option(s)" not in log
                for name, value in expected.items():
                    assert f'"TargetedDataExtraction|{name}" -> "{value}"' in log, (argument, name)
                for option, value in analysis_expected.items():
                    prefix, name = option[1:].rsplit(":", 1)
                    key = f"{prefix}|{name}"
                    assert f'"{key}" -> "{value}"' in log, (argument, option)
    print(f"PASS: metrics, invalid q-values, resource failures, SVG edge cases, PDF failures, {len(sys.argv) - 1} captured commands")
    if os.environ.get("OPEN_DIA") and len(sys.argv) > 1:
        print("PASS: real OpenDIA resolved shared speed, inference and export options; this parser check does not process the data")


if __name__ == "__main__":
    main()
