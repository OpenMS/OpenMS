"""Optional confirmation with ProSE itself on a subset of the spectra.

Two runs, both on the top-TIC MS2 subset:

1. closed search with the proposed settings and the calibration pass on;
   ProSE's -summary_out YAML reports per file the precursor and fragment
   error (median, MAD, recommended ppm) and the PSM/peptide/protein counts;
2. optional open search (precursor window of +-500 Da turns ProSE's open
   mode on) with -out_mod_analysis_dir; the delta-mass tables show the mass
   shifts that the closed search is missing (a Cys alkylation that is not
   set, labels, frequent variable modifications).

This needs a ProSE binary (OpenMS build or installer); the Python
mini-search in hypotheses.py covers the same questions when there is none.
"""

from __future__ import annotations

import csv
import glob
import os
import shutil
import subprocess
import tempfile

from .evidence import Evidence
from .spectra import write_subset_mzml


def find_prose(path=None):
    return path or shutil.which("ProSE")


def _run(cmd, log):
    log("  $ " + " ".join(cmd[:6]) + (" ..." if len(cmd) > 6 else ""))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:] or r.stdout[-2000:])
    return r


def confirm_with_prose(prose, mzml, fasta, cli_args, n_spectra=5000, open_search=False, log=print):
    import yaml
    ev = []
    with tempfile.TemporaryDirectory(prefix="prose_auto_") as tmp:
        sub = os.path.join(tmp, "subset.mzML")
        n = write_subset_mzml(mzml, sub, n_spectra)
        summ = os.path.join(tmp, "summary.yaml")
        args = [a for a in cli_args]
        # make sure the calibration pass runs and reports the errors
        if "-Search:calibration:enabled" in args:
            i = args.index("-Search:calibration:enabled")
            args[i + 1] = "true"
        else:
            args += ["-Search:calibration:enabled", "true"]
        _run([prose, "-in", sub, "-database", fasta, "-out", os.path.join(tmp, "out.idXML"),
              "-summary_out", summ, "-no_summary"] + args, log)
        with open(summ) as fh:
            s = yaml.safe_load(fh)
        pf = (s.get("per_file") or [{}])[0]
        stats = {k: pf.get(k) for k in ("ms2_spectra", "target_psms", "unique_peptides", "unique_proteins",
                                         "achieved_psm_fdr")}
        frac = (pf.get("target_psms") or 0) / max(1, n)
        ev.append(Evidence("organism", os.path.basename(fasta), min(0.9, 0.3 + frac), "search", "ProSE subset search",
                           f"ProSE identified {pf.get('target_psms')} of {n} subset spectra "
                           f"({frac:.0%}) against {os.path.basename(fasta)}", stats))
        for kind in ("precursor", "fragment"):
            e = pf.get(f"{kind}_error") or {}
            if e.get("valid"):
                ev.append(Evidence(f"{kind}_offset", {"value": e["median_ppm"], "unit": "ppm"}, 0.9, "search",
                                   "ProSE calibration", f"ProSE {kind} error median {e['median_ppm']} ppm, "
                                   f"MAD {e['mad_ppm']} ppm"))
                ev.append(Evidence(f"{kind}_tolerance", {"value": e["recommended_ppm"], "unit": "ppm"}, 0.9,
                                   "search", "ProSE calibration",
                                   f"ProSE recommends {e['recommended_ppm']} ppm from its calibration PSMs"))
        if open_search:
            mods = os.path.join(tmp, "mods")
            os.makedirs(mods)
            oargs = [a for a in args]
            for key, val in (("-Search:precursor:mass_tolerance_lower", "500"),
                             ("-Search:precursor:mass_tolerance_upper", "500"),
                             ("-Search:precursor:mass_tolerance_unit", "Da"),
                             ("-Search:calibration:enabled", "false")):
                if key in oargs:
                    oargs[oargs.index(key) + 1] = val
                else:
                    oargs += [key, val]
            _run([prose, "-in", sub, "-database", fasta, "-out", os.path.join(tmp, "open.idXML"),
                  "-out_mod_analysis_dir", mods, "-no_summary"] + oargs, log)
            shifts = []
            for f in glob.glob(os.path.join(mods, "*")):
                if "aggregate" in os.path.basename(f).lower() or len(glob.glob(os.path.join(mods, "*"))) == 1:
                    with open(f) as fh:
                        rows = list(csv.DictReader(fh, delimiter="\t"))
                    shifts.extend(rows[:15])
            if shifts:
                ev.append(Evidence("variable_mods", [r for r in shifts], 0.5, "search", "ProSE open search",
                                   "most frequent mass shifts of the open search (see the table)"))
    return ev
