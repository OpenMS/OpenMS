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
_TRUE_TOKENS  = {"true", "1", "yes", "on"}
_FALSE_TOKENS = {"false", "0", "no", "off"}

def _as_bool(v) -> bool:
    if isinstance(v, bool):
        return v
    if v is None:
        return False
    s = str(v).strip().lower()
    return s in ("true", "1", "yes", "on")

def _is_bool_param(p) -> bool:
    t = getattr(p, "type", None)
    return (t is bool) or (isinstance(t, str) and t.lower() == "bool")


def _coerce_bool_from_argv(argv: list[str], name: str, current: bool) -> bool:
    """
    If the user explicitly set -<name> true/false (or =true/=false), or used
    -no-<name>, override the current boolean. Otherwise keep current.
    """
    flag    = f"-{name}"
    noflag  = f"-no-{name}"

    it = enumerate(argv)
    for i, tok in it:
        if tok == noflag:
            return False

        if tok == flag:
            # Value in next token?
            if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                v = argv[i + 1].strip().lower()
                if v in _TRUE_TOKENS:  return True
                if v in _FALSE_TOKENS: return False
            # No explicit value -> treat as “set True” and continue scanning
            continue

        if tok.startswith(flag + "="):
            v = tok.split("=", 1)[1].strip().lower()
            if v in _TRUE_TOKENS:  return True
            if v in _FALSE_TOKENS: return False

    return current


def _is_nullish(v) -> bool:
    # CTDopts uses a private sentinel class named "_Null" for “no value”
    return v is None or type(v).__name__ in {"_Null", "Null"} or v == ""

def _first_value_for_flag(argv: list[str], name: str) -> str | None:
    """
    Return the raw token following -<name>, or a value from -<name>=VALUE if used.
    Does NOT interpret the value; it’s exactly what the shell passed.
    """
    flag = f"-{name}"
    for i, t in enumerate(argv):
        if t == flag:
            if i + 1 < len(argv):
                return argv[i + 1]
            return ""
        if t.startswith(flag + "="):
            return t.split("=", 1)[1]
    return None

def _tok_extra(s: str) -> list[str]:
    """Split shell-like string into tokens; undo '+--', '+-' and '\\-' escapes."""
    if _is_nullish(s) or not isinstance(s, str) or not s.strip():
        return []
    # allow an explicit '-- ' sentinel inside the value (optional nicety)
    s = s.strip()
    if s.startswith("-- "):
        s = s[3:]
    toks = shlex.split(s)
    out: list[str] = []
    for t in toks:
        if t.startswith("+--") or t.startswith("+-") or t.startswith("\\-"):
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
    cmd = ["pyprophet", *args] if exe == "pyprophet" else [exe, *args]
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
        type=int,
        description="Worker threads used by PyProphet scoring.",
        default=1,
    )
    model.add(
        "log_level",
        required=False,
        type=str,
        description="PyProphet log level.",
        default="INFO",
        choices=["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"],
    )
    model.add(
        "pyprophet_executable",
        required=False,
        type=str,
        description="PyProphet executable path or 'pyprophet' to run as a Python module.",
        default="pyprophet",
    )
    model.add(
        "dry_run",
        required=False,
        type=bool,
        description="Print commands without executing.",
        default=False,
    )

    # scoring
    model.add(
        "scoring:run_score",
        required=False,
        type=bool,
        description="Run semi-supervised scoring.",
        default=True,
    )
    model.add(
        "scoring:level",
        required=False,
        type=str,
        description="OSW data level(s) for scoring (CSV).",
        default="ms1ms2",
        choices=["ms1", "ms2", "ms1ms2", "transition", "alignment"],
    )
    model.add(
        "scoring:classifier",
        required=False,
        type=str,
        description="Classifier for semi-supervised scoring.",
        default="SVM",
        choices=["LDA", "SVM", "XGBoost", "HistGradientBoosting"],
    )
    model.add(
        "scoring:subsample_ratio",
        required=False,
        type=float,
        description="Subsample ratio for training; weights applied to full set.",
        default=1.0,
    )
    model.add(
        "scoring:apply_weights",
        required=False,
        type="input-file",
        is_list=False,
        description="Pre-trained weights file to apply (skip training)."
        # file_formats intentionally omitted; pyprophet accepts several formats
    )
    model.add(
        "scoring:autotune",
        required=False,
        type=bool,
        description="Autotune classifier hyperparameters.",
        default=False,
    )
    model.add(
        "scoring:xeval_num_iter",
        required=False,
        type=int,
        description="Cross-validation iterations for semi-supervised learning.",
        default=10,
    )
    model.add(
        "scoring:ss_num_iter",
        required=False,
        type=int,
        description="Semi-supervised iterations.",
        default=10,
    )
    model.add(
        "scoring:ss_initial_fdr",
        required=False,
        type=float,
        description="Initial FDR cutoff for best scoring targets.",
        default=0.15,
    )
    model.add(
        "scoring:ss_iteration_fdr",
        required=False,
        type=float,
        description="Iteration FDR cutoff for best scoring targets.",
        default=0.05,
    )
    model.add(
        "scoring:ss_main_score",
        required=False,
        type=str,
        description='Main score to start SSL (e.g., "var_xcorr_shape" or "auto").',
        default="auto",
    )
    model.add(
        "scoring:ss_scale_features",
        required=False,
        type=bool,
        description="Scale / standardize features to unit variance.",
        default=False,
    )
    model.add(
        "scoring:extra",
        required=False,
        type=str,
        description='Extra args passed verbatim to "pyprophet score" '
                    '(e.g. "--parametric --pi0_lambda 0.4 0 0").',
        default="",
    )

    # infer
    model.add(
        "infer:context",
        required=False,
        type=str,
        description="Context(s) for peptide/protein/gene inference (CSV).",
        default="global",
        choices=["run-specific", "experiment-wide", "global"],
    )
    model.add(
        "infer:run_peptide",
        required=False,
        type=bool,
        description="Run peptide inference.",
        default=True,
    )
    model.add(
        "infer:run_protein",
        required=False,
        type=bool,
        description="Run protein inference.",
        default=True,
    )
    model.add(
        "infer:run_gene",
        required=False,
        type=bool,
        description="Run gene inference.",
        default=False,
    )
    model.add(
        "infer:peptide:extra",
        required=False,
        type=str,
        description='Extra args for "pyprophet infer peptide".',
        default="",
    )
    model.add(
        "infer:protein:extra",
        required=False,
        type=str,
        description='Extra args for "pyprophet infer protein".',
        default="",
    )
    model.add(
        "infer:gene:extra",
        required=False,
        type=str,
        description='Extra args for "pyprophet infer gene".',
        default="",
    )

    # IPF
    model.add(
        "infer:run_ipf",
        required=False,
        type=bool,
        description="Run peptidoform inference (IPF).",
        default=False,
    )
    model.add(
        "infer:ipf_ms1_scoring",
        required=False,
        type=bool,
        description="Use MS1 precursor information in IPF scoring.",
        default=False,
    )
    model.add(
        "infer:ipf_ms2_scoring",
        required=False,
        type=bool,
        description="Use MS2 precursor information in IPF scoring.",
        default=False,
    )
    model.add(
        "infer:ipf_h0",
        required=False,
        type=bool,
        description="Include possibility that peak groups are not covered by peptidoform space.",
        default=False,
    )
    model.add(
        "infer:ipf_grouped_fdr",
        required=False,
        type=bool,
        description="Compute grouped FDR instead of pooled FDR.",
        default=False,
    )
    model.add(
        "infer:ipf_max_precursor_pep",
        required=False,
        type=float,
        description="Max PEP for precursors considered in IPF.",
        default=0.7,
    )
    model.add(
        "infer:ipf_max_peakgroup_pep",
        required=False,
        type=float,
        description="Max PEP for peak-groups considered in IPF.",
        default=0.7,
    )
    model.add(
        "infer:ipf_max_precursor_peakgroup_pep",
        required=False,
        type=float,
        description="Max integrated precursor-peakgroup PEP.",
        default=0.4,
    )
    model.add(
        "infer:ipf_max_transition_pep",
        required=False,
        type=float,
        description="Max PEP for transitions considered in IPF.",
        default=0.6,
    )

    # export
    model.add(
        "export:score_report",
        required=False,
        type=bool,
        description="Generate a summary score report (PDF).",
        default=True,
    )
    model.add(
        "export:score_plots",
        required=False,
        type=bool,
        description="Export score plots (PDF).",
        default=False,
    )
    model.add(
        "export:run_tsv",
        required=False,
        type=bool,
        description="Export Proteomics TSV (long format).",
        default=True,
    )
    model.add(
        "export:tsv:extra",
        required=False,
        type=str,
        description='Extra args for "pyprophet export tsv".',
        default="",
    )
    model.add(
        "export:run_matrix",
        required=False,
        type=bool,
        description="Export Proteomics matrix TSV.",
        default=True,
    )
    model.add(
        "export:matrix:extra",
        required=False,
        type=str,
        description='Extra args for "pyprophet export matrix".',
        default="",
    )
    model.add(
        "export:run_compound",
        required=False,
        type=bool,
        description="Export Small Molecules TSV.",
        default=False,
    )
    model.add(
        "export:compound:format",
        required=False,
        type=str,
        description="Small-molecule export format.",
        default="legacy_merged",
        choices=["matrix", "legacy_merged"],
    )
    model.add(
        "export:compound:max_rs_peakgroup_qvalue",
        required=False,
        type=float,
        description="Filter by max run-specific peakgroup q-value.",
        default=0.05,
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
    arg_dict = parseCTDCommandLinePure(sys.argv[1:], model)
    
    # Recover any :extra raw strings 
    extra_keys = ["scoring:extra", "infer:peptide:extra", "infer:protein:extra",
                "infer:gene:extra", "export:tsv:extra", "export:matrix:extra"]
    for k in extra_keys:
        if _is_nullish(arg_dict.get(k, None)):
            raw = _first_value_for_flag(sys.argv[1:], k)
            if raw is not None:
                arg_dict[k] = raw

    # **Boolean overrides from argv** (`-infer:run_protein false`, etc.)
    for p in model.list_parameters():
        if _is_bool_param(p):
            cur = _as_bool(arg_dict.get(p.name, getattr(p, "default", False)))
            arg_dict[p.name] = _coerce_bool_from_argv(sys.argv[1:], p.name, cur)

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
            w = arg_dict.get("scoring:apply_weights", None)
            if not _is_nullish(w):
                score_cmd += ["--apply_weights", str(w)]
            if _as_bool(arg_dict.get("scoring:autotune", False)):
                score_cmd += ["--autotune"]
            ms = arg_dict.get("scoring:ss_main_score", "auto")
            if not _is_nullish(ms):
                score_cmd += ["--ss_main_score", str(ms)]
            if _as_bool(arg_dict.get("scoring:ss_scale_features", False)):
                score_cmd += ["--ss_scale_features"]
            score_cmd += _tok_extra(arg_dict.get("scoring:extra", ""))

            env = os.environ.copy()
            if arg_dict.get("scoring:classifier", "SVM") == "HistGradientBoosting":
                env.setdefault("OMP_NUM_THREADS", str(max(1, (os.cpu_count() or 1) // max(1, threads))))

            print(f"====== Scoring level: {lvl} ======")
            print("CMD:", "pyprophet", *score_cmd)
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
