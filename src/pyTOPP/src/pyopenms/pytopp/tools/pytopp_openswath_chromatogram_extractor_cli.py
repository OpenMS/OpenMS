#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.9"
# dependencies = [
#   "pyopenms>=3",
#   "CTDopts @ git+https://github.com/WorkflowConversion/CTDopts.git",
# ]
# ///
from __future__ import annotations

import sys
from pathlib import Path
from typing import List

import pyopenms as poms
from CTDopts.CTDopts import CTDModel


from pyopenms.pytopp.ctdsupport import (
    CTDHelpConfig,
    CTDHelpPrinter,
    PyTOPPArgParser,
)



def build_model() -> CTDModel:
    """
    CTD model for OpenSwathChromatogramExtractor (pyTOPP wrapper).

    Options mirror the legacy script's arguments, with types/defaults exposed via CTDopts.
    """
    model = CTDModel(
        name="OpenSwathChromatogramExtractor",
        version=getattr(poms, "__version__", "3"),
        description="Extract chromatograms (XICs) from MS data given a TraML transition list.",
        manual="RTF",
        docurl="https://pyopenms.readthedocs.io/en/latest/apidocs/_autosummary/pyopenms/pyopenms.ChromatogramExtractor.html",
        category="Quantitation",
        executableName="OpenSwathChromatogramExtractor",
        executablePath="",
    )

    # I/O
    model.add(
        "in",
        required=True,
        type="input-file",
        is_list=True,
        file_formats=["mzML", "mzXML", "mzData"],
        description="Input spectrum file(s). Multiple files are concatenated into one output.",
    )
    model.add(
        "tr",
        required=True,
        type="input-file",
        is_list=False,
        file_formats=["TraML"],
        description="TraML file with transitions.",
    )
    model.add(
        "out",
        required=True,
        type="output-file",
        is_list=False,
        file_formats=["mzML", "mzXML", "mzData"],
        description="Output chromatogram mzML.",
    )

    # Extraction options
    model.add(
        "extraction_window",
        required=False,
        type=float,
        default=0.05,
        description="m/z extraction window (Th unless --ppm is set).",
    )
    model.add(
        "ppm",
        required=False,
        type=bool,
        default=False,
        description="Interpret extraction_window in ppm instead of Th.",
    )
    model.add(
        "rt_extraction_window",
        required=False,
        type=float,
        default=-1.0,
        description="RT extraction window (seconds). Use -1 to disable.",
    )
    model.add(
        "extraction_function",
        required=False,
        type=str,
        default="tophat",
        choices=["tophat", "bartlett"],
        description="Extraction function.",
    )

    # SWATH-specific helpers
    model.add(
        "is_swath",
        required=False,
        type=bool,
        default=False,
        description="Treat input as SWATH map and preselect transitions per window.",
    )
    model.add(
        "min_upper_edge_dist",
        required=False,
        type=float,
        default=0.0,
        description="Minimal distance to upper SWATH edge (Th) to consider a precursor.",
    )

    return model


def _load_targeted(tr_path: Path) -> poms.TargetedExperiment:
    targeted = poms.TargetedExperiment()
    poms.TraMLFile().load(str(tr_path), targeted)
    return targeted


def _load_ms(in_path: Path) -> poms.MSExperiment:
    exp = poms.MSExperiment()
    poms.FileHandler().loadExperiment(str(in_path), exp)
    return exp


def _append_chroms(dst: poms.MSExperiment, src: List[poms.OSChromatogram]) -> None:
    for chrom in src:
        # Convert OSChromatogram to MSChromatogram
        # TODO: Maybe OSChromatogram should have a conversion method, if one doesn't already exist?
        int = chrom.getIntensityArray()
        rt = chrom.getTimeArray()
        ms_chrom = poms.MSChromatogram()
        ms_chrom.set_peaks([rt, int])
        dst.addChromatogram(ms_chrom)


def _tag_processing_smoothing(exp: poms.MSExperiment) -> None:
    dp = poms.DataProcessing()
    pa = poms.DataProcessing().ProcessingAction().SMOOTHING
    dp.setProcessingActions(set([pa]))
    chroms = exp.getChromatograms()
    for c in chroms:
        lst = c.getDataProcessing()
        lst.append(dp)
        c.setDataProcessing(lst)
    exp.setChromatograms(chroms)


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    model = build_model()
    cfg = CTDHelpConfig(
        tool_name=f"{model.name} (pyTOPP)",
        binary_name="OpenSwathChromatogramExtractor",
        show_description=True,
    )
    printer = CTDHelpPrinter(model, cfg=cfg)

    # Help
    if "-h" in argv or "--help" in argv:
        printer.print(advanced=False)
        return 0
    if "--help-advanced" in argv or "--helphelp" in argv:
        printer.print(advanced=True)
        return 0

    # Parse via CTDopts (pure mode; no Param synchronization needed here)
    parser = PyTOPPArgParser(model)
    try:
        args = parser.parse(argv, mode="pure")
    except SystemExit as e:
        return int(getattr(e, "code", 2))

    # Resolve I/O
    in_list = args.get("in", [])
    if isinstance(in_list, str):
        in_list = [in_list]
    in_paths = [Path(p) for p in in_list]
    tr_path = Path(args.get("tr", ""))
    out_path = Path(args.get("out", ""))

    if not in_paths:
        print("No input (-in) provided.", file=sys.stderr)
        return 2
    if not tr_path:
        print("No transition file (-tr) provided.", file=sys.stderr)
        return 2
    if not out_path:
        print("No output (-out) provided.", file=sys.stderr)
        return 2

    # Options
    mz_win = float(args.get("extraction_window", 0.05))
    ppm = bool(args.get("ppm", False))
    rt_win = float(args.get("rt_extraction_window", -1.0))
    func = str(args.get("extraction_function", "tophat"))
    is_swath = bool(args.get("is_swath", False))
    min_upper_edge = float(args.get("min_upper_edge_dist", 0.0))

    # Load transitions
    targeted = _load_targeted(tr_path)

    # Aggregate chromatograms from all inputs
    output = poms.MSExperiment()
    trafo = poms.TransformationDescription()

    for ip in in_paths:
        exp = _load_ms(ip)

        # If requested, restrict transitions to the SWATH window
        transition_exp_used = poms.TargetedExperiment()
        proceed = True
        if is_swath:
            helper = poms.OpenSwathHelper()
            proceed = helper.checkSwathMapAndSelectTransitions(
                exp, targeted, transition_exp_used, min_upper_edge
            )
        else:
            transition_exp_used = targeted

        if not proceed:
            continue
        
        # Prepare extraction coordinates
        output_chromatograms = []
        extraction_coordinates = []
        poms.ChromatogramExtractor.prepare_coordinates(
            output_chromatograms,
            extraction_coordinates,
            transition_exp_used,
            rt_win,
            False,  # ms1
            0,      # ms1_isotopes
        )
        
        # Create input map
        input_map = poms.SpectrumAccessOpenMS(exp)

        extractor = poms.ChromatogramExtractor()
        extractor.extractChromatograms(
            input_map,
            output_chromatograms,
            extraction_coordinates,
            mz_win,
            ppm,
            -1, # Assumes no IonMobility dimension
            func,
        )
        _append_chroms(output, output_chromatograms)

    # annotate & store
    _tag_processing_smoothing(output)
    poms.FileHandler().storeExperiment(str(out_path), output)

    print("OpenSwathChromatogramExtractor (pyTOPP) finished successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
