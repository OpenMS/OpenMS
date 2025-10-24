#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.9"
# dependencies = [
#   "pyprophet>=2",
#   "pyopenms>=3",
#   "CTDopts @ git+https://github.com/WorkflowConversion/CTDopts.git",
# ]
# ///


from __future__ import annotations

import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable, List

import pyopenms as pms

# CTDopts / pyTOPP helpers
import CTDopts  # noqa: F401
from CTDopts.CTDopts import CTDModel
from CTDsupport import HelpConfig, print_model_help, addParamToCTDopts, parseCTDCommandLine, parseCTDCommandLinePure, parse_args_ctdopts
from CTDopts.CTDopts import flatten_dict


# -------------------------- small helpers ---------------------------------
def _as_bool(v) -> bool:
    if isinstance(v, bool):
        return v
    if v is None:
        return False
    s = str(v).strip().lower()
    return s in ("true", "1", "yes", "on")


def _tok_extra(s: str) -> List[str]:
    """Split shell-like string into tokens; undo legacy '+--' or '\\-' prefixes."""
    if not s:
        return []
    toks = shlex.split(s)
    out: List[str] = []
    for t in toks:
        if t.startswith("+--") or t.startswith("+-") or t.startswith(r"\-"):
            out.append(t[1:])
        else:
            out.append(t)
    return out


def _canonical_levels(csv: str) -> List[str]:
    allowed = {"ms1", "ms2", "ms1ms2", "transition", "alignment"}
    raw = [t.strip().lower() for t in (csv or "").split(",") if t.strip()]
    seen = set(raw)
    if not seen:
        seen = {"ms2"}
    bad = [t for t in raw if t not in allowed]
    if bad:
        raise ValueError(
            f"Unknown level(s): {','.join(bad)}. Allowed: {','.join(sorted(allowed))}"
        )
    if "ms1ms2" in seen and "ms2" in seen:
        raise ValueError("Invalid combination: 'ms1ms2' cannot be combined with 'ms2'.")
    bases = {"ms1", "ms2", "ms1ms2"}
    if "transition" in seen and not (seen & bases):
        raise ValueError("Invalid combination: 'transition' requires ms1|ms2|ms1ms2.")
    order = ["ms1", "ms2", "ms1ms2", "transition", "alignment"]
    return [l for l in order if l in seen]


def _parse_contexts(csv: str) -> List[str]:
    allowed = ["run-specific", "experiment-wide", "global"]
    if not csv:
        return ["run-specific"]
    out = []
    for t in [x.strip().lower() for x in csv.split(",") if x.strip()]:
        if t == "runspecific":
            t = "run-specific"
        if t == "experimentwide":
            t = "experiment-wide"
        if t not in allowed:
            raise ValueError(f"Unknown context '{t}'. Allowed: {', '.join(allowed)}")
        if t not in out:
            out.append(t)
    return out or ["run-specific"]


def _call_pyprophet(exe: str, args: Iterable[str], env=None) -> None:
    # Prefer module mode so we use the current (uv) environment
    cmd = [sys.executable, "-m", "pyprophet", *args] if exe == "pyprophet" else [exe, *args]
    p = subprocess.run(cmd, env=env)
    if p.returncode != 0:
        raise RuntimeError(f"PyProphet failed with exit code {p.returncode}")


def _derive_tsv(out_osw: Path) -> Path:
    return out_osw.with_suffix(".tsv")


def _derive_matrix(out_osw: Path) -> Path:
    return out_osw.with_name(out_osw.stem + ".matrix.tsv")


def _derive_compound(out_osw: Path) -> Path:
    return out_osw.with_name(out_osw.stem + ".compound.tsv")


# -------------------------- CTD model -------------------------------------
def build_model() -> CTDModel:
    model = CTDModel(
        name="PyProphet",
        version="1.0",
        description=(
            "pyTOPP tool to run PyProphet scoring, optional peptide/protein/gene inference, "
            "peptidoform inference (IPF), and exports on OSW files."
        ),
        manual="RTF",
        docurl="https://openswath.org/",
        category="Quantitation",
        executableName="pyprophet",
        executablePath=""
    )

    # I/O
    model.add(
        "in",
        required=True,
        type="input-file",
        is_list=True,
        file_formats=["OSW"],
        description="Input OSW file(s). If >1, a merge step is executed first."
    )
    model.add(
        "out",
        required=False,
        type="output-file",
        is_list=False,
        file_formats=["OSW"],
        description="Output OSW file. If omitted and a single input is given, scoring is in-place."
    )

    # logging / threading / exe
    model.add(
        "threads",
        required=False,
        type="int",
        description="Worker threads used by PyProphet scoring.",
        default_value="1"
    )
    model.add(
        "log_level",
        required=False,
        type="string",
        description="PyProphet log level.",
        default_value="INFO",
        restrictions=["TRACE","DEBUG","INFO","SUCCESS","WARNING","ERROR","CRITICAL"]
    )
    model.add(
        "pyprophet_executable",
        required=False,
        type="string",
        description="PyProphet executable path or 'pyprophet' to run as a Python module.",
        default_value="pyprophet"
    )
    model.add(
        "dry_run",
        required=False,
        type="bool",
        description="Print commands without executing.",
        default_value="false",
        restrictions=["true","false"]
    )

    # scoring
    model.add(
        "scoring:run_score",
        required=False,
        type="bool",
        description="Run semi-supervised scoring.",
        default_value="true",
        restrictions=["true","false"]
    )
    model.add(
        "scoring:level",
        required=False,
        type="string",
        description="OSW data level(s) for scoring (CSV).",
        default_value="ms1ms2",
        restrictions=["ms1","ms2","ms1ms2","transition","alignment"]
    )
    model.add(
        "scoring:classifier",
        required=False,
        type="string",
        description="Classifier for semi-supervised scoring.",
        default_value="SVM",
        restrictions=["LDA","SVM","XGBoost","HistGradientBoosting"]
    )
    model.add(
        "scoring:subsample_ratio",
        required=False,
        type="float",
        description="Subsample ratio for training; weights applied to full set.",
        default_value="1.0"
    )
    model.add(
        "scoring:apply_weights",
        required=False,
        type="input-file",
        is_list=False,
        description="Pre-trained weights file to apply (skip training)."
        # file_formats can be omitted; pyprophet accepts various formats
    )
    model.add(
        "scoring:autotune",
        required=False,
        type="bool",
        description="Autotune classifier hyperparameters.",
        default_value="false",
        restrictions=["true","false"]
    )
    model.add(
        "scoring:xeval_num_iter",
        required=False,
        type="int",
        description="Cross-validation iterations for semi-supervised learning.",
        default_value="10"
    )
    model.add(
        "scoring:ss_num_iter",
        required=False,
        type="int",
        description="Semi-supervised iterations.",
        default_value="10"
    )
    model.add(
        "scoring:ss_initial_fdr",
        required=False,
        type="float",
        description="Initial FDR cutoff for best scoring targets.",
        default_value="0.15"
    )
    model.add(
        "scoring:ss_iteration_fdr",
        required=False,
        type="float",
        description="Iteration FDR cutoff for best scoring targets.",
        default_value="0.05"
    )
    model.add(
        "scoring:ss_main_score",
        required=False,
        type="string",
        description='Main score to start SSL (e.g., "var_xcorr_shape" or "auto").',
        default_value="auto"
    )
    model.add(
        "scoring:ss_scale_features",
        required=False,
        type="bool",
        description="Scale / standardize features to unit variance.",
        default_value="false",
        restrictions=["true","false"]
    )
    model.add(
        "scoring:extra",
        required=False,
        type="string",
        description='Extra args passed verbatim to "pyprophet score" '
                    '(e.g. "--parametric --pi0_lambda 0.4 0 0").',
        default_value=""
    )

    # infer
    model.add(
        "infer:context",
        required=False,
        type="string",
        description="Context(s) for peptide/protein/gene inference (CSV).",
        default_value="global",
        restrictions=["run-specific","experiment-wide","global"]
    )
    model.add(
        "infer:run_peptide",
        required=False,
        type="bool",
        description="Run peptide inference.",
        default_value="true",
        restrictions=["true","false"]
    )
    model.add(
        "infer:run_protein",
        required=False,
        type="bool",
        description="Run protein inference.",
        default_value="true",
        restrictions=["true","false"]
    )
    model.add(
        "infer:run_gene",
        required=False,
        type="bool",
        description="Run gene inference.",
        default_value="false",
        restrictions=["true","false"]
    )
    model.add(
        "infer:peptide:extra",
        required=False,
        type="string",
        description='Extra args for "pyprophet infer peptide".',
        default_value=""
    )
    model.add(
        "infer:protein:extra",
        required=False,
        type="string",
        description='Extra args for "pyprophet infer protein".',
        default_value=""
    )
    model.add(
        "infer:gene:extra",
        required=False,
        type="string",
        description='Extra args for "pyprophet infer gene".',
        default_value=""
    )

    # IPF
    model.add(
        "infer:run_ipf",
        required=False,
        type="bool",
        description="Run peptidoform inference (IPF).",
        default_value="false",
        restrictions=["true","false"]
    )
    model.add(
        "infer:ipf_ms1_scoring",
        required=False,
        type="bool",
        description="Use MS1 precursor information in IPF scoring.",
        default_value="false",
        restrictions=["true","false"]
    )
    model.add(
        "infer:ipf_ms2_scoring",
        required=False,
        type="bool",
        description="Use MS2 precursor information in IPF scoring.",
        default_value="false",
        restrictions=["true","false"]
    )
    model.add(
        "infer:ipf_h0",
        required=False,
        type="bool",
        description="Include possibility that peak groups are not covered by peptidoform space.",
        default_value="false",
        restrictions=["true","false"]
    )
    model.add(
        "infer:ipf_grouped_fdr",
        required=False,
        type="bool",
        description="Compute grouped FDR instead of pooled FDR.",
        default_value="false",
        restrictions=["true","false"]
    )
    model.add(
        "infer:ipf_max_precursor_pep",
        required=False,
        type="float",
        description="Max PEP for precursors considered in IPF.",
        default_value="0.7"
    )
    model.add(
        "infer:ipf_max_peakgroup_pep",
        required=False,
        type="float",
        description="Max PEP for peak-groups considered in IPF.",
        default_value="0.7"
    )
    model.add(
        "infer:ipf_max_precursor_peakgroup_pep",
        required=False,
        type="float",
        description="Max integrated precursor-peakgroup PEP.",
        default_value="0.4"
    )
    model.add(
        "infer:ipf_max_transition_pep",
        required=False,
        type="float",
        description="Max PEP for transitions considered in IPF.",
        default_value="0.6"
    )

    # export
    model.add(
        "export:score_report",
        required=False,
        type="bool",
        description="Generate a summary score report (PDF).",
        default_value="true",
        restrictions=["true","false"]
    )
    model.add(
        "export:score_plots",
        required=False,
        type="bool",
        description="Export score plots (PDF).",
        default_value="false",
        restrictions=["true","false"]
    )
    model.add(
        "export:run_tsv",
        required=False,
        type="bool",
        description="Export Proteomics TSV (long format).",
        default_value="true",
        restrictions=["true","false"]
    )
    model.add(
        "export:tsv:extra",
        required=False,
        type="string",
        description='Extra args for "pyprophet export tsv".',
        default_value=""
    )
    model.add(
        "export:run_matrix",
        required=False,
        type="bool",
        description="Export Proteomics matrix TSV.",
        default_value="true",
        restrictions=["true","false"]
    )
    model.add(
        "export:matrix:extra",
        required=False,
        type="string",
        description='Extra args for "pyprophet export matrix".',
        default_value=""
    )
    model.add(
        "export:run_compound",
        required=False,
        type="bool",
        description="Export Small Molecules TSV.",
        default_value="false",
        restrictions=["true","false"]
    )
    model.add(
        "export:compound:format",
        required=False,
        type="string",
        description="Small-molecule export format.",
        default_value="legacy_merged",
        restrictions=["matrix","legacy_merged"]
    )
    model.add(
        "export:compound:max_rs_peakgroup_qvalue",
        required=False,
        type="float",
        description="Filter by max run-specific peakgroup q-value.",
        default_value="0.05"
    )

    return model

# -------------------------- main logic ------------------------------------
def main() -> int:
    model = build_model()
    
    cfg = HelpConfig(
    tool_name=f"{model.name} (pyTOPP)",
    binary_name="PyProphet.py",
    show_description=True,
    advanced_by_name={
            "scoring:extra", "scoring:subsample_ratio", "scoring:apply_weights",
            "scoring:xeval_num_iter", "scoring:ss_num_iter",
            "scoring:ss_initial_fdr", "scoring:ss_iteration_fdr",
            "scoring:ss_main_score", "scoring:ss_scale_features", "scoring:autotune",
            "infer:peptide:extra", "infer:protein:extra", "infer:gene:extra",
            "infer:ipf_ms1_scoring", "infer:ipf_ms2_scoring", "infer:ipf_h0",
            "infer:ipf_grouped_fdr", "infer:ipf_max_precursor_pep", "infer:ipf_max_peakgroup_pep",
            "infer:ipf_max_precursor_peakgroup_pep", "infer:ipf_max_transition_pep",
            "export:tsv:extra", "export:matrix:extra",
        },
    )

    argv = sys.argv[1:]
    if "-h" in argv or "--help" in argv:
        print_model_help(model, advanced=False, cfg=cfg)
        sys.exit(0)
    if "--help-advanced" in argv or "--helphelp" in argv:
        print_model_help(model, advanced=True, cfg=cfg)
        sys.exit(0)

    # No algorithm defaults (we're just wrapping PyProphet); still provide an empty Param
    defaults = pms.Param()
    addParamToCTDopts(defaults, model)

    # parse command line / handle -write_ini/-ini
    arg_dict, _ = parseCTDCommandLinePure(sys.argv, model)

    in_list = arg_dict.get("in", [])
    if isinstance(in_list, str):
        in_list = [in_list]
    in_list = [Path(x) for x in in_list]
    if not in_list:
        print("No input (-in) given.", file=sys.stderr)
        return 2

    out_osw = arg_dict.get("out", "")
    out_osw = Path(out_osw) if out_osw else (in_list[0] if len(in_list) == 1 else None)
    if out_osw is None:
        print("Multiple inputs provided but no -out given.", file=sys.stderr)
        return 2

    pyexe = arg_dict.get("pyprophet_executable", "pyprophet")
    loglevel = arg_dict.get("log_level", "INFO")
    dry = _as_bool(arg_dict.get("dry_run", False))
    threads = int(arg_dict.get("threads", 1))

    # Merge (if needed)
    current_in = in_list[0]
    if len(in_list) > 1:
        merged = Path(tempfile.gettempdir()) / f"{next(tempfile._get_candidate_names())}_merged.osw"
        merge_cmd = [
            "--log-level", loglevel, "--no-log-colorize",
            "merge", "osw",
            "--out", str(merged),
            "--template", str(in_list[0]),
            *[str(p) for p in in_list],
        ]
        print("====== Merging Input OSWs ======")
        print("CMD:", sys.executable, "-m", "pyprophet", *merge_cmd)
        if not dry:
            _call_pyprophet(pyexe, merge_cmd)
        current_in = merged

    # Score
    if _as_bool(arg_dict.get("scoring:run_score", True)):
        try:
            levels = _canonical_levels(arg_dict.get("scoring:level", "ms1ms2"))
        except ValueError as e:
            print(str(e), file=sys.stderr)
            return 2

        for lvl in levels:
            score_cmd = [
                "--log-level", loglevel, "--no-log-colorize",
                "score",
                "--in", str(current_in),
                "--out", str(out_osw),
                "--level", lvl,
                "--classifier", arg_dict.get("scoring:classifier", "SVM"),
                "--subsample_ratio", str(arg_dict.get("scoring:subsample_ratio", 1.0)),
                "--threads", str(threads),
                "--xeval_num_iter", str(arg_dict.get("scoring:xeval_num_iter", 10)),
                "--ss_num_iter", str(arg_dict.get("scoring:ss_num_iter", 10)),
                "--ss_initial_fdr", str(arg_dict.get("scoring:ss_initial_fdr", 0.15)),
                "--ss_iteration_fdr", str(arg_dict.get("scoring:ss_iteration_fdr", 0.05)),
            ]
            if arg_dict.get("scoring:apply_weights"):
                score_cmd += ["--apply_weights", arg_dict["scoring:apply_weights"]]
            if _as_bool(arg_dict.get("scoring:autotune", False)):
                score_cmd += ["--autotune"]
            if arg_dict.get("scoring:ss_main_score", "auto"):
                score_cmd += ["--ss_main_score", arg_dict.get("scoring:ss_main_score", "auto")]
            if _as_bool(arg_dict.get("scoring:ss_scale_features", False)):
                score_cmd += ["--ss_scale_features"]
            score_cmd += _tok_extra(arg_dict.get("scoring:extra", ""))

            env = os.environ.copy()
            if arg_dict.get("scoring:classifier", "SVM") == "HistGradientBoosting":
                env.setdefault("OMP_NUM_THREADS", str(max(1, (os.cpu_count() or 1) // max(1, threads))))

            print(f"====== Scoring level: {lvl} ======")
            print("CMD:", sys.executable, "-m", "pyprophet", *score_cmd)
            if not dry:
                _call_pyprophet(pyexe, score_cmd, env=env)
            current_in = out_osw

    # Inference (peptide/protein/gene)
    run_pep = _as_bool(arg_dict.get("infer:run_peptide", True))
    run_pro = _as_bool(arg_dict.get("infer:run_protein", True))
    run_gene = _as_bool(arg_dict.get("infer:run_gene", False))
    if run_pep or run_pro or run_gene:
        try:
            contexts = _parse_contexts(arg_dict.get("infer:context", "global"))
        except ValueError as e:
            print(str(e), file=sys.stderr)
            return 2

        def infer_one(what: str, enabled: bool, ctx: str, extra_key: str):
            if not enabled:
                return
            cmd = [
                "--log-level", loglevel, "--no-log-colorize",
                "infer", what,
                "--in", str(out_osw),
                "--out", str(out_osw),
                "--context", ctx,
            ]
            cmd += _tok_extra(arg_dict.get(extra_key, ""))
            print(f"====== {what} inference ({ctx}) ======")
            print("CMD:", sys.executable, "-m", "pyprophet", *cmd)
            if not dry:
                _call_pyprophet(pyexe, cmd)

        for ctx in contexts:
            infer_one("peptide", run_pep, ctx, "infer:peptide:extra")
            infer_one("protein", run_pro, ctx, "infer:protein:extra")
            infer_one("gene", run_gene, ctx, "infer:gene:extra")

    # IPF
    if _as_bool(arg_dict.get("infer:run_ipf", False)):
        def _ff(k: str) -> List[str]:
            return [f"--{k}"] if _as_bool(arg_dict.get(k, False)) else [f"--no-{k}"]

        cmd = [
            "--log-level", loglevel, "--no-log-colorize",
            "infer", "peptidoform",
            "--in", str(out_osw),
            "--out", str(out_osw),
            *_ff("infer:ipf_ms1_scoring"),
            *_ff("infer:ipf_ms2_scoring"),
            *_ff("infer:ipf_h0"),
            *_ff("infer:ipf_grouped_fdr"),
            "--ipf_max_precursor_pep", str(arg_dict.get("infer:ipf_max_precursor_pep", 0.7)),
            "--ipf_max_peakgroup_pep", str(arg_dict.get("infer:ipf_max_peakgroup_pep", 0.7)),
            "--ipf_max_precursor_peakgroup_pep", str(arg_dict.get("infer:ipf_max_precursor_peakgroup_pep", 0.4)),
            "--ipf_max_transition_pep", str(arg_dict.get("infer:ipf_max_transition_pep", 0.6)),
        ]
        print("====== Peptidoform inference ======")
        print("CMD:", sys.executable, "-m", "pyprophet", *cmd)
        if not dry:
            _call_pyprophet(pyexe, cmd)

    # Exports
    do_report = _as_bool(arg_dict.get("export:score_report", True))
    do_plots = _as_bool(arg_dict.get("export:score_plots", False))
    do_tsv = _as_bool(arg_dict.get("export:run_tsv", True))
    do_matrix = _as_bool(arg_dict.get("export:run_matrix", True))
    do_comp = _as_bool(arg_dict.get("export:run_compound", False))

    if do_comp and (do_tsv or do_matrix):
        print("Options export:run_compound and (export:run_tsv/export:run_matrix) are mutually exclusive.",
              file=sys.stderr)
        return 2

    def run_export(sub: str, out_path: Path | None = None, extra_key: str | None = None):
        cmd = ["--log-level", loglevel, "--no-log-colorize", "export", sub, "--in", str(out_osw)]
        if out_path:
            cmd += ["--out", str(out_path)]
        if extra_key:
            cmd += _tok_extra(arg_dict.get(extra_key, ""))
        print(f"====== Export {sub} ======")
        print("CMD:", sys.executable, "-m", "pyprophet", *cmd)
        if not dry:
            _call_pyprophet(pyexe, cmd)

    if do_report:
        run_export("score-report")
    if do_plots:
        run_export("score-plots")
    if do_comp:
        cmd = [
            "--log-level", loglevel, "--no-log-colorize",
            "export", "compound",
            "--format", arg_dict.get("export:compound:format", "legacy_merged"),
            "--in", str(out_osw),
            "--out", str(_derive_compound(out_osw)),
            "--max_rs_peakgroup_qvalue", str(arg_dict.get("export:compound:max_rs_peakgroup_qvalue", 0.05)),
        ]
        print(f"====== Export compound ({arg_dict.get('export:compound:format', 'legacy_merged')}) ======")
        print("CMD:", sys.executable, "-m", "pyprophet", *cmd)
        if not dry:
            _call_pyprophet(pyexe, cmd)
    if do_tsv:
        run_export("tsv", _derive_tsv(out_osw), "export:tsv:extra")
    if do_matrix:
        run_export("matrix", _derive_matrix(out_osw), "export:matrix:extra")

    print("PyProphetAdapter (pyTOPP) finished successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
