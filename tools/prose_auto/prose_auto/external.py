"""Adapters for external tools: FASTag and DIAspeXtract.

FASTag (github.com/okohlbacher/FASTag) reads sequence tags from MS2 spectra.
prose-auto uses two of its outputs:

* `-species`: a genus call from tag 7-mers against its bundled index of ~50
  reference proteomes (model organisms, contaminants, pathogens, microbiome).
  Database free, so it also works when the organism is not in the panel.
  Genus level only: close relatives (mammals) share most 7-mers, and reagent
  genera (Bos, Sus, Oryctolagus, ...) show up in almost every run.
* `-recon_fasta` + `-delta_out`: tags placed on a proteome with NO fixed
  modification; the histogram of flank mass gaps shows the Cys alkylation
  (+57.021 carbamidomethyl, +45.988 methylthio, ...) and labels as peaks.

DIAspeXtract (github.com/okohlbacher/diaspextract) turns diaPASEF (Bruker .d,
or DIA mzML) into pseudo-MS/MS spectra with one precursor each, so the DDA
detectors and the mini-search can run on DIA data too.
"""

from __future__ import annotations

import csv
import os
import subprocess
import tempfile

from .evidence import Evidence

REAGENT_GENERA = {"Bos", "Sus", "Oryctolagus", "Gallus", "Capra", "Ovis", "Equus", "Lysobacter",
                  "Achromobacter", "Staphylococcus"}

# Cys alkylation and isobaric labels as mass gaps when the search sets no fixed modification
DELTA_CANDIDATES = {
    "Carbamidomethyl (C)": 57.021464,
    "Propionamide (C)": 71.037114,
    "Nethylmaleimide (C)": 125.047679,
    "Methylthio (C)": 45.987721,
    "Carboxymethyl (C)": 58.005479,
}


def _run(cmd, log, timeout=6 * 3600):
    log("  $ " + " ".join(str(c) for c in cmd))
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError((r.stderr or r.stdout)[-1500:])
    return r


def _tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _frag_args(frag_tol, frag_unit):
    return ["-fragment_tolerance", frag_tol, "-fragment_tolerance_unit", frag_unit]


def fastag_species(fastag, mzml, frag_tol=20.0, frag_unit="ppm", fixed=("Carbamidomethyl (C)",),
                   subsample=0.0, threads=0, log=print):
    """Genus call from FASTag -species."""
    with tempfile.TemporaryDirectory(prefix="fastag_") as tmp:
        out = os.path.join(tmp, "tags.tsv")
        sp = os.path.join(tmp, "species.tsv")
        cmd = [fastag, "-in", mzml, "-out", out, "-species", "-tag_length", 7, "-species_out", sp,
               "-threads", threads] + _frag_args(frag_tol, frag_unit)
        cmd += ["-fixed_modifications"] + list(fixed)
        if subsample:
            cmd += ["-subsample_fraction", subsample]
        _run(cmd, log)
        rows = _tsv(sp) if os.path.exists(sp) else []
    if not rows:
        return [Evidence("organism", "unknown", 0.1, "data", "FASTag species", "FASTag reported no taxon")]
    for r in rows:
        r["adjusted"] = float(r.get("adjusted") or r.get("observed") or 0)
        r["observed"] = float(r.get("observed") or 0)
    rows.sort(key=lambda r: -r["adjusted"])
    table = [{k: r.get(k) for k in ("name", "observed", "adjusted", "expected", "qvalue")} for r in rows[:15]]
    # the top non-reagent genus is the call; a reagent genus on top is reported as such
    top = rows[0]
    real = [r for r in rows if r["name"] not in REAGENT_GENERA]
    call = real[0] if real and (top["name"] in REAGENT_GENERA and real[0]["adjusted"] >= 0.5 * top["adjusted"]) else top
    nxt = next((r for r in rows if r is not call), None)
    margin = (call["adjusted"] - (nxt["adjusted"] if nxt else 0)) / max(call["adjusted"], 1)
    n = call["observed"]
    conf = 0.25 + min(0.35, 0.12 * (n ** 0.5) / 3) + min(0.3, margin)
    conf = min(0.85, conf if n >= 20 else min(conf, 0.35))
    detail = (f"FASTag genus call {call['name']} from {int(n)} tag-carrying spectra"
              + (f"; next {nxt['name']} ({nxt['adjusted']:.0f} vs {call['adjusted']:.0f})" if nxt else "")
              + ("; reagent genus on top skipped" if call is not top else ""))
    return [Evidence("organism", call["name"], conf, "data", "FASTag species", detail, {"taxa": table})]


def fastag_alkylation(fastag, mzml, fasta, frag_tol=20.0, frag_unit="ppm", threads=0, log=print):
    """Cys alkylation from the FASTag mass-gap histogram without fixed mods."""
    with tempfile.TemporaryDirectory(prefix="fastag_") as tmp:
        out = os.path.join(tmp, "tags.tsv")
        delta = os.path.join(tmp, "delta.tsv")
        cmd = [fastag, "-in", mzml, "-out", out, "-tag_length", 4, "-extension", 6, "-threads", threads,
               "-recon_fasta", fasta, "-recon_out", os.path.join(tmp, "recon.tsv"), "-delta_out", delta,
               "-fixed_modifications"] + _frag_args(frag_tol, frag_unit)
        _run(cmd, log)
        rows = _tsv(delta) if os.path.exists(delta) else []
    hist = []
    for r in rows:
        try:
            hist.append((float(r.get("delta")), int(float(r.get("spectra") or 0)), r.get("top_interps", "")))
        except (TypeError, ValueError):
            continue
    total = sum(c for _, c, _ in hist)
    tol = 0.02 if frag_unit == "ppm" else 0.5
    counts = {name: sum(c for d, c, _ in hist if abs(d - m) <= tol) for name, m in DELTA_CANDIDATES.items()}
    zero = sum(c for d, c, _ in hist if abs(d) <= tol)
    stats = {"counts": counts, "zero_bin": zero, "spectra": total,
             "top_bins": sorted(hist, key=lambda x: -x[1])[:15]}
    best = max(counts, key=counts.get)
    second = sorted(counts.values())[-2] if len(counts) > 1 else 0
    if counts[best] < 10:
        return [Evidence("cys_alkylation", "unknown", 0.15, "data", "FASTag mass gaps",
                         f"no Cys alkylation peak in the FASTag gap histogram ({total} spectra placed)", stats)]
    ratio = (counts[best] + 1) / (second + 1)
    conf = min(0.85, 0.4 + 0.15 * ratio ** 0.5)
    return [Evidence("cys_alkylation", best, conf, "data", "FASTag mass gaps",
                     f"{counts[best]} placed spectra with a {DELTA_CANDIDATES[best]:+.3f} Da gap "
                     f"(next candidate {second}); searched without fixed modifications", stats)]


# DIAspeXtract on a constant-mobility copy: mobility tracing, weighting, charge veto and twin folding off
NO_MOBILITY_OPTIONS = ["-trace:detector", "openms", "-gate:delta_im", "1.0", "-assembly:im_weight_sigma", "0",
                       "-charge:im_charge_veto", "false", "-assembly:twin_im_tolerance", "-1",
                       "-trace:im_slice_merge", "-1"]


def with_constant_mobility(src, dst, k0=1.0):
    """Copy of a DIA mzML without ion mobility in which every peak carries 1/K0 = k0, so DIAspeXtract accepts it."""
    import numpy as np
    import pyopenms as oms
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(src), exp)
    spectra = exp.getSpectra()  # iterating an MSExperiment yields copies
    for s in spectra:
        fda = oms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.set_data(np.full(s.size(), k0, dtype=np.float32))
        s.setFloatDataArrays([f for f in s.getFloatDataArrays() if f.getName() != "Ion Mobility"] + [fda])
        s.setDriftTimeUnit(oms.DriftTimeUnit.VSSC)
    exp.setSpectra(spectra)
    exp.updateRanges()
    oms.MzMLFile().store(str(dst), exp)
    return dst


def diaspextract(exe, inp, out_mzml, threads=0, extra=(), log=print):
    """diaPASEF .d (or DIA mzML) -> pseudo-MS/MS mzML."""
    cmd = [exe, "-in", inp, "-out", out_mzml] + (["-threads", threads] if threads else []) + list(extra)
    _run(cmd, log, timeout=12 * 3600)
    return out_mzml
