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
import logging
from pathlib import Path

import pyopenms as poms
from CTDopts.CTDopts import CTDModel

from pyopenms.pytopp.ctdsupport import (
    CTDHelpConfig,
    CTDHelpPrinter,
    PyTOPPArgParser,
    addParamToCTDopts,
)


from pyopenms.pytopp.common import addDataProcessing as _addDP



# -------------------------- CTD model -------------------------------------
def build_model() -> CTDModel:
    """
    Create a CTD model that contains:
      - I/O flags: -in/-out
      - All PeakPickerHiRes algorithm parameters (converted from OpenMS Param)
    """
    model = CTDModel(
        name="PeakPickerHiRes",
        version=getattr(poms, "__version__", "3"),
        description=(
            "High-resolution peak picking on profile MS data. "
            "Wraps OpenMS PeakPickerHiRes."
        ),
        manual="RTF",
        docurl="https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/",
        category="Signal Processing",
        executableName="PeakPickerHiRes",
        executablePath="",
    )

    # I/O
    model.add(
        "in",
        required=True,
        type="input-file",
        is_list=False,
        file_formats=["mzML", "mzXML", "mzData"],
        description="Input profile spectrum file.",
    )
    model.add(
        "out",
        required=True,
        type="output-file",
        is_list=False,
        file_formats=["mzML", "mzXML", "mzData"],
        description="Output picked (centroided) spectrum file.",
    )

    # Convert algorithm defaults from OpenMS Param to CTD model options
    defaults_param = poms.PeakPickerHiRes().getDefaults()
    addParamToCTDopts(defaults_param, model)

    return model


# -------------------------- main logic ------------------------------------
def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    model = build_model()
    cfg = CTDHelpConfig(
        tool_name=f"{model.name} (pyTOPP)",
        binary_name="PeakPickerHiRes",
        show_description=True,
        # You can move verbose/rarely-used params to advanced by name if desired:
        # advanced_by_name={"signal_to_noise", "ms_levels", ...},
    )
    printer = CTDHelpPrinter(model, cfg=cfg)

    # Help first (matches your other tools)
    if "-h" in argv or "--help" in argv:
        printer.print(advanced=False)
        return 0
    if "--help-advanced" in argv or "--helphelp" in argv:
        printer.print(advanced=True)
        return 0

    # Parse CLI/CTD in "pyopenms" mode, so we get an OpenMS Param in sync
    om_param = poms.Param()  # will be filled from parse() via CTD roundtrip
    parser = PyTOPPArgParser(model)
    try:
        validated, om_param = parser.parse(argv, mode="pyopenms", openms_param=om_param)
    except SystemExit as e:
        # parser prints its own message (unknown flags, etc.)
        return int(getattr(e, "code", 2))

    # Resolve I/O
    in_path = Path(validated.get("in", ""))
    out_path = Path(validated.get("out", ""))

    if not in_path:
        print("Missing -in.", file=sys.stderr)
        return 2
    if not out_path:
        print("Missing -out.", file=sys.stderr)
        return 2

    # Load input
    exp_in = poms.MSExperiment()
    mzml = poms.MzMLFile()
    mzml.setLogType(poms.LogType.CMD)
    mzml.load(str(in_path), exp_in)

    # Quick sanity: warn if data already looks centroided
    if exp_in.size() > 0:
        est = poms.PeakTypeEstimator()
        try:
            s0 = exp_in[0]
            if est.estimateType(s0) == poms.SpectrumSettings.SpectrumType.PEAKS:
                logging.warning("Input appears already centroided (PEAKS).")
        except Exception:
            pass

    # Require sorted spectra (OpenMS algorithm expects this)
    if any(not s.isSorted() for s in exp_in):
        print("Not all spectra are sorted by m/z.", file=sys.stderr)
        return 2

    # Configure and run picker
    picker = poms.PeakPickerHiRes()
    if om_param is not None and isinstance(om_param, poms.Param):
        picker.setParameters(om_param)

    exp_out = poms.MSExperiment()
    picker.pickExperiment(exp_in, exp_out, True)

    # Annotate processing (if helper available)
    if _addDP is not None:
        try:
            exp_out = _addDP(
                exp_out,
                picker.getParameters() if hasattr(picker, "getParameters") else om_param,
                poms.DataProcessing.ProcessingAction.PEAK_PICKING,
            )
        except Exception:
            # If helper shape changed, fall back silently
            pass

    # Store output (FileHandler handles multiple formats by extension)
    fh = poms.FileHandler()
    fh.storeExperiment(str(out_path), exp_out)

    print("PeakPickerHiRes (pyTOPP) finished successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
