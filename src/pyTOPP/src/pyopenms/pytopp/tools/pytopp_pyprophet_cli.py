#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.9"
# dependencies = [
#   "pyprophet>=3",
#   "pyopenms>=3",
#   "CTDopts @ git+https://github.com/WorkflowConversion/CTDopts.git",
# ]
# ///
from __future__ import print_function


# --- make sure our vendored pyopenms/pytopp is on the pyopenms namespace -----
def _ensure_pyopenms_namespace():
    import importlib
    import importlib.machinery
    import os
    import sys
    import types

    # find the "site" root that contains pyopenms/
    here = os.path.abspath(__file__)
    d = os.path.dirname(here)
    site = None
    for _ in range(10):
        if os.path.isdir(os.path.join(d, "pyopenms")):
            site = d
            break
        nd = os.path.dirname(d)
        if nd == d:
            break
        d = nd
    if not site:
        return

    if site not in sys.path:
        sys.path.insert(0, site)

    pyo_dir = os.path.join(site, "pyopenms")

    # Prefer a real installed pyopenms if available
    try:
        m = importlib.import_module("pyopenms")
    except Exception:
        m = None

    if m is not None:
        # extend its namespace to also search our vendored path
        try:
            path_list = list(m.__path__)
        except Exception:
            path_list = []
        if pyo_dir not in path_list:
            try:
                m.__path__.append(pyo_dir)
            except Exception:
                m.__path__ = path_list + [pyo_dir]
        return

    # otherwise create a namespace stub pointing at our vendored path
    m = types.ModuleType("pyopenms")
    m.__package__ = "pyopenms"
    m.__path__ = [pyo_dir]
    m.__spec__ = importlib.machinery.ModuleSpec("pyopenms", loader=None, is_package=True)
    sys.modules["pyopenms"] = m

_ensure_pyopenms_namespace()
# -----------------------------------------------------------------------------

import os
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path
from typing import Iterable, List, Optional

import CTDopts  # noqa: F401
import pyopenms as poms
from CTDopts.CTDopts import CTDModel
from pyopenms.pytopp.util import (
    as_bool,
    coerce_all_bools_from_argv,
    is_not_nullish,
    is_nullish,
    recover_missing_extras,
    run_tool,
    tok_extra,
)
from pyopenms.pytopp.CTDsupport import (
    CTDHelpConfig,
    CTDHelpPrinter,
    PyTOPPArgParser,
)

PYPROPHET_VERSION = version("pyprophet")

# -------------------------- PyProphet-specific helpers ---------------------
def _canonical_levels(csv: str) -> List[str]:
    allowed = {"ms1", "ms2", "ms1ms2", "transition", "alignment"}
    raw = [t.strip().lower() for t in (csv or "").split(",") if t.strip()]
    seen = set(raw) or {"ms2"}
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
        t = {"runspecific": "run-specific", "experimentwide": "experiment-wide"}.get(t, t)
        if t not in allowed:
            raise ValueError(f"Unknown context '{t}'. Allowed: {', '.join(allowed)}")
        if t not in out:
            out.append(t)
    return out or ["run-specific"]


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
        version=PYPROPHET_VERSION,
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
    model.add("in",  required=True,  type="input-file", is_list=True,  file_formats=["osw", "oswpq", "oswpqd", "parquet"],
              description="Input OSW file(s). If input is sqlite-based OSW and more than 1 file is passed, a merge step is executed first.")
    model.add("out", required=False, type="output-file", is_list=False, file_formats=["osw", "oswpq", "oswpqd", "parquet"],
              description="Output OSW file. If omitted and a single input is given, scoring is done in-place.")

    # logging / threading / exe
    model.add("threads",             required=False, type=int,  default=1,
              description="Worker threads used by PyProphet scoring.")
    model.add("log_level",           required=False, type=str,  default="INFO",
              choices=["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"],
              description="PyProphet log level.")
    model.add("pyprophet_executable", required=False, type=str, default="pyprophet",
              description="PyProphet executable path or 'pyprophet' to run as a Python module.")
    model.add("dry_run",             required=False, type=bool, default=False,
              description="Print commands without executing.")

    # scoring
    model.add("scoring:run_score",       required=False, type=bool,  default=True,  description="Run semi-supervised scoring.")
    model.add("scoring:level",           required=False, type=str,   default="ms1ms2",
              description="Either 'ms1', 'ms2', 'ms1ms2', 'transition', or 'alignment'; the data level selected for scoring. 'ms1ms2 integrates both MS1- and MS2-level scores and can be used instead of 'ms2'-level results. Can be comma-separated to perform scoring per level, i.e., ms1ms2,transition") # We do not use choices here to allow for comma-separated multi-values. We use a custom parser and validator instead.
    model.add("scoring:classifier",      required=False, type=str,   default="SVM",
              choices=["LDA", "SVM", "XGBoost", "HistGradientBoosting"], description="Classifier for semi-supervised scoring.")
    model.add("scoring:subsample_ratio", required=False, type=float, default=1.0,
              description="Subsample ratio for training; weights applied to full set.")
    model.add("scoring:apply_weights",   required=False, type="input-file", is_list=False,
              description="Pre-trained weights file to apply (skip training).")
    model.add("scoring:autotune",        required=False, type=bool,  default=False,
              description="Autotune classifier hyperparameters.")
    model.add("scoring:xeval_num_iter",  required=False, type=int,   default=10,
              description="Cross-validation iterations for semi-supervised learning.")
    model.add("scoring:ss_num_iter",     required=False, type=int,   default=10,
              description="Semi-supervised iterations.")
    model.add("scoring:ss_initial_fdr",  required=False, type=float, default=0.15,
              description="Initial FDR cutoff for best scoring targets.")
    model.add("scoring:ss_iteration_fdr",required=False, type=float, default=0.05,
              description="Iteration FDR cutoff for best scoring targets.")
    model.add("scoring:ss_main_score",   required=False, type=str,   default="auto",
              description='Main score to start SSL (e.g., "var_xcorr_shape" or "auto").')
    model.add("scoring:ss_scale_features",required=False, type=bool, default=False,
              description="Scale / standardize features to unit variance.")
    model.add("scoring:extra",           required=False, type=str,   default="",
              description='Extra args passed verbatim to "pyprophet score" (e.g. "+--parametric +--pi0_lambda 0.4 0 0").')

    # infer
    model.add("infer:context",        required=False, type=str, default="global",
              description="Context(s) for peptide/protein/gene inference ('run-specific', 'experiment-wide', 'global'). Can be comma-separated to perform scoring per level, i.e., global,run-specific") # We do not use choices here to allow for comma-separated multi-values. We use a custom parser and validator instead.
    model.add("infer:run_peptide",    required=False, type=bool, default=True,  description="Run peptide inference.")
    model.add("infer:run_protein",    required=False, type=bool, default=True,  description="Run protein inference.")
    model.add("infer:run_gene",       required=False, type=bool, default=False, description="Run gene inference.")
    model.add("infer:peptide:extra",  required=False, type=str,  default="",    description='Extra args for "pyprophet infer peptide".')
    model.add("infer:protein:extra",  required=False, type=str,  default="",    description='Extra args for "pyprophet infer protein".')
    model.add("infer:gene:extra",     required=False, type=str,  default="",    description='Extra args for "pyprophet infer gene".')

    # IPF
    model.add("infer:run_ipf",                    required=False, type=bool,  default=False, description="Run peptidoform inference (IPF).")
    model.add("infer:ipf_ms1_scoring",            required=False, type=bool,  default=False, description="Use MS1 precursor information in IPF scoring.")
    model.add("infer:ipf_ms2_scoring",            required=False, type=bool,  default=False, description="Use MS2 precursor information in IPF scoring.")
    model.add("infer:ipf_h0",                     required=False, type=bool,  default=False, description="Include hypothesis that peak groups not covered by peptidoform space.")
    model.add("infer:ipf_grouped_fdr",            required=False, type=bool,  default=False, description="Compute grouped FDR instead of pooled FDR.")
    model.add("infer:ipf_max_precursor_pep",      required=False, type=float, default=0.7,   description="Max PEP for precursors considered in IPF.")
    model.add("infer:ipf_max_peakgroup_pep",      required=False, type=float, default=0.7,   description="Max PEP for peak-groups considered in IPF.")
    model.add("infer:ipf_max_precursor_peakgroup_pep", required=False, type=float, default=0.4, description="Max integrated precursor-peakgroup PEP.")
    model.add("infer:ipf_max_transition_pep",     required=False, type=float, default=0.6,   description="Max PEP for transitions considered in IPF.")

    # export
    model.add("export:score_report",   required=False, type=bool, default=True,  description="Generate a summary score report (PDF).")
    model.add("export:score_plots",    required=False, type=bool, default=False, description="Export score plots (PDF).")
    model.add("export:run_tsv",        required=False, type=bool, default=True,  description="Export Proteomics TSV (long format).")
    model.add("export:tsv:extra",      required=False, type=str,  default="",    description='Extra args for "pyprophet export tsv".')
    model.add("export:run_matrix",     required=False, type=bool, default=True,  description="Export Proteomics matrix TSV.")
    model.add("export:matrix:extra",   required=False, type=str,  default="",    description='Extra args for "pyprophet export matrix".')
    model.add("export:run_compound",   required=False, type=bool, default=False, description="Export Small Molecules TSV.")
    model.add("export:compound:format",required=False, type=str,  default="legacy_merged",
              choices=["matrix", "legacy_merged"], description="Small-molecule export format.")
    model.add("export:compound:max_rs_peakgroup_qvalue", required=False, type=float, default=0.05,
              description="Filter by max run-specific peakgroup q-value.")

    return model


# -------------------------- main logic ------------------------------------
def main() -> int:
    model = build_model()

    cfg = CTDHelpConfig(
        tool_name=f"{model.name} v{model.version} (pyTOPP)",
        binary_name="pyprophet",
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
            "export:compound:format", "export:compound:max_rs_peakgroup_qvalue",
        },
    )
    cli_opt_printer = CTDHelpPrinter(model, cfg=cfg)

    argv = sys.argv[1:]
    if "-h" in argv or "--help" in argv:
        cli_opt_printer.print(advanced=False); return 0
    if "--help-advanced" in argv or "--helphelp" in argv:
        cli_opt_printer.print(advanced=True); return 0
    
    # Parse CTD/CLI; this returns validated & type-cast values
    parser = PyTOPPArgParser(model)
    arg_dict = parser.parse(sys.argv[1:], mode="pure")
    
    # Repair dropped :extra strings & coerce booleans based on raw argv
    recover_missing_extras(arg_dict, argv, [
        "scoring:extra", "infer:peptide:extra", "infer:protein:extra",
        "infer:gene:extra", "export:tsv:extra", "export:matrix:extra",
    ])
    coerce_all_bools_from_argv(model, arg_dict, argv)

    # Inputs / outputs
    in_list = arg_dict.get("in", [])
    if isinstance(in_list, str):
        in_list = [in_list]
    in_list = [Path(x) for x in in_list]
    if not in_list:
        print("No input (-in) given.", file=sys.stderr); return 2

    out_osw = arg_dict.get("out", "")
    out_osw = Path(out_osw) if is_not_nullish(out_osw) else (in_list[0] if len(in_list) == 1 else None)
    if out_osw is None:
        print("Multiple inputs provided but no -out given.", file=sys.stderr); return 2

    pyexe    = arg_dict.get("pyprophet_executable", "pyprophet")
    loglevel = arg_dict.get("log_level", "INFO")
    dry      = as_bool(arg_dict.get("dry_run", False))
    threads  = int(arg_dict.get("threads", 1))

    # Merge if needed
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
        run_tool(pyexe, merge_cmd, dry=dry)
        current_in = merged

    # Score
    if as_bool(arg_dict.get("scoring:run_score", True)):
        try:
            levels = _canonical_levels(arg_dict.get("scoring:level", "ms1ms2"))
        except ValueError as e:
            print(str(e), file=sys.stderr); return 2

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
            if not is_nullish(w):
                score_cmd += ["--apply_weights", str(w)]
            if as_bool(arg_dict.get("scoring:autotune", False)):
                score_cmd += ["--autotune"]
            ms = arg_dict.get("scoring:ss_main_score", "auto")
            if not is_nullish(ms):
                score_cmd += ["--ss_main_score", str(ms)]
            if as_bool(arg_dict.get("scoring:ss_scale_features", False)):
                score_cmd += ["--ss_scale_features"]
            score_cmd += tok_extra(arg_dict.get("scoring:extra", ""))

            env = os.environ.copy()
            if arg_dict.get("scoring:classifier", "SVM") == "HistGradientBoosting":
                env.setdefault("OMP_NUM_THREADS", str(max(1, (os.cpu_count() or 1) // max(1, threads))))

            print(f"====== Scoring level: {lvl} ======")
            run_tool(pyexe, score_cmd, env=env, dry=dry)
            current_in = out_osw

    # Inference (peptide/protein/gene)
    run_pep  = as_bool(arg_dict.get("infer:run_peptide", True))
    run_pro  = as_bool(arg_dict.get("infer:run_protein", True))
    run_gene = as_bool(arg_dict.get("infer:run_gene", False))
    if run_pep or run_pro or run_gene:
        try:
            contexts = _parse_contexts(arg_dict.get("infer:context", "global"))
        except ValueError as e:
            print(str(e), file=sys.stderr); return 2

        def infer_one(what: str, enabled: bool, ctx: str, extra_key: str):
            if not enabled:
                return
            cmd = [
                "--log-level", loglevel, "--no-log-colorize",
                "infer", what,
                "--in", str(out_osw),
                "--out", str(out_osw),
                "--context", ctx,
            ] + tok_extra(arg_dict.get(extra_key, ""))
            print(f"====== {what} inference ({ctx}) ======")
            run_tool(pyexe, cmd, dry=dry)

        for ctx in contexts:
            infer_one("peptide", run_pep,  ctx, "infer:peptide:extra")
            infer_one("protein", run_pro,  ctx, "infer:protein:extra")
            infer_one("gene",    run_gene, ctx, "infer:gene:extra")

    # IPF
    if as_bool(arg_dict.get("infer:run_ipf", False)):
        def _ff(k: str) -> List[str]:
            return [f"--{k}"] if as_bool(arg_dict.get(k, False)) else [f"--no-{k}"]

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
        run_tool(pyexe, cmd, dry=dry)

    # Exports
    do_report = as_bool(arg_dict.get("export:score_report", True))
    do_plots  = as_bool(arg_dict.get("export:score_plots", False))
    do_tsv    = as_bool(arg_dict.get("export:run_tsv", True))
    do_matrix = as_bool(arg_dict.get("export:run_matrix", True))
    do_comp   = as_bool(arg_dict.get("export:run_compound", False))

    if do_comp and (do_tsv or do_matrix):
        print("Options export:run_compound and (export:run_tsv/export:run_matrix) are mutually exclusive.", file=sys.stderr)
        return 2

    def run_export(sub: str, out_path: Optional[Path] = None, extra_key: Optional[str] = None):
        cmd = ["--log-level", loglevel, "--no-log-colorize", "export", sub, "--in", str(out_osw)]
        if out_path:
            cmd += ["--out", str(out_path)]
        if extra_key:
            cmd += tok_extra(arg_dict.get(extra_key, ""))
        print(f"====== Export {sub} ======")
        run_tool(pyexe, cmd, dry=dry)

    if do_report:  run_export("score-report")
    if do_plots:   run_export("score-plots")
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
        run_tool(pyexe, cmd, dry=dry)
    if do_tsv:     run_export("tsv",    _derive_tsv(out_osw),    "export:tsv:extra")
    if do_matrix:  run_export("matrix", _derive_matrix(out_osw), "export:matrix:extra")

    print("PyProphet (pyTOPP) finished successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
