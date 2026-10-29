"""Instrument, analyzer, fragmentation and acquisition mode from file metadata.

Sources: the mzML instrument configuration (model CV term, analyzers), the
precursor activation CV terms, Thermo filter strings ('FTMS + c NSI Full ms2
582.32@hcd27.00 [100-1205]'), isolation windows and ion mobility arrays.
"""

from __future__ import annotations

import re
from collections import Counter

import numpy as np

from .evidence import Evidence

# Instrument model -> (instrument family, PeptDeep instrument class,
#                      default precursor ppm, default fragment tol, unit)
# The defaults are only priors: tolerance.py measures the data.
_INSTRUMENTS = [
    (r"astral", ("Orbitrap Astral", "ThermoTOF", 5.0, 10.0, "ppm")),
    (r"exploris", ("Orbitrap Exploris", "QE", 10.0, 20.0, "ppm")),
    (r"q ?exactive|\bqe\b", ("Q Exactive", "QE", 10.0, 20.0, "ppm")),
    (r"lumos|eclipse|fusion|ascend|tribrid", ("Orbitrap Tribrid", "Lumos", 10.0, 20.0, "ppm")),
    (r"orbitrap (xl|velos|elite|classic)|ltq orbitrap", ("LTQ Orbitrap", "QE", 10.0, 0.5, "Da")),
    (r"ltq|lcq|velos pro|linear ion trap", ("Ion trap", "QE", 0.5, 0.5, "Da")),
    (r"timstof|tims tof", ("timsTOF", "timsTOF", 15.0, 20.0, "ppm")),
    (r"zenotof|triple ?tof|x500|tripletof", ("SCIEX TOF", "SciexTOF", 20.0, 30.0, "ppm")),
    (r"synapt|xevo|vion|select series", ("Waters Q-TOF", "SciexTOF", 20.0, 30.0, "ppm")),
    (r"impact|maxis|compact", ("Bruker Q-TOF", "timsTOF", 20.0, 30.0, "ppm")),
]

HIGH_RES = {"ORBITRAP", "FOURIERTRANSFORM", "CYCLOTRON", "TOF", "ASTRAL"}


def classify_instrument(name: str):
    n = (name or "").lower()
    for pat, info in _INSTRUMENTS:
        if re.search(pat, n):
            return info
    return None


def _analyzer_from_filter(fs: str):
    fs = fs.strip()
    for prefix, name in (("FTMS", "ORBITRAP"), ("ITMS", "IT"), ("ASTMS", "ASTRAL"),
                         ("TOFMS", "TOF")):
        if fs.startswith(prefix):
            return name
    return None


def detect_instrument(run) -> list:
    ev = []
    name = run.instrument_name
    info = classify_instrument(name)
    if name:
        ev.append(Evidence("instrument", name, 0.95, "metadata", "mzML instrument CV",
                           f"Instrument model in the mzML header: {name}"))
    if info:
        ev.append(Evidence("instrument_class", info[1], 0.9, "metadata", "model lookup",
                           f"{name} belongs to the {info[0]} family"))
    # Per-level analyzers: filter strings first, file-level analyzers second.
    for level, scans, key in ((1, run.ms1, "ms1_analyzer"), (2, run.ms2, "ms2_analyzer")):
        c = Counter(a for a in (_analyzer_from_filter(s.filter_string) for s in scans) if a)
        if c:
            a, n = c.most_common(1)[0]
            ev.append(Evidence(key, a, 0.95, "metadata", "filter string",
                               f"{n} of {len(scans)} MS{level} filter strings start with "
                               f"{'FTMS' if a == 'ORBITRAP' else a}", {"counts": dict(c)}))
        elif run.analyzers:
            # A single analyzer in the header covers both levels; with several
            # (e.g. Orbitrap + ion trap) the level is ambiguous.
            # Quadrupoles in the header are mass filters, not the analyzer
            # that records the spectrum.
            uniq = sorted(set(run.analyzers) - {"QUADRUPOLE"}) or sorted(set(run.analyzers))
            conf = 0.8 if len(uniq) == 1 else 0.4
            hr = [a for a in uniq if a in HIGH_RES]
            a = uniq[0] if len(uniq) == 1 else (hr[0] if hr and level == 1 else
                                                 next((a for a in uniq if a not in HIGH_RES), uniq[0]))
            ev.append(Evidence(key, a, conf, "metadata", "mzML analyzers",
                               f"Analyzers in the mzML header: {', '.join(uniq)}"))
    # Fragmentation
    acts = Counter("+".join(s.activation) for s in run.ms2 if s.activation)
    if acts:
        a, n = acts.most_common(1)[0]
        frac = n / max(1, sum(acts.values()))
        ev.append(Evidence("fragmentation", a, 0.6 + 0.35 * frac, "metadata", "activation CV",
                           f"{n} of {sum(acts.values())} MS2 spectra use {a}",
                           {"counts": dict(acts)}))
    energies = np.array([s.energy for s in run.ms2 if s.energy > 0])
    if energies.size:
        med = float(np.median(energies))
        stepped = len(set(np.round(energies, 1))) > 1
        ev.append(Evidence("collision_energy", round(med, 1), 0.8, "metadata", "activation energy",
                           f"median collision energy {med:.1f}" + (" (several values: stepped or per charge)" if stepped else ""),
                           {"distinct": sorted(set(np.round(energies, 1).tolist()))[:20]}))
    im = sum(1 for s in run.ms2 + run.ms1 if s.ion_mobility)
    total = len(run.ms2) + len(run.ms1)
    ev.append(Evidence("ion_mobility", im > 0.5 * max(1, total), 0.9, "metadata", "IM arrays",
                       f"{im} of {total} spectra carry ion mobility information"))
    return ev


def detect_acquisition(run) -> list:
    """DDA vs DIA vs PRM from the isolation windows.

    DDA: narrow windows (<= ~3 Th) centred on many different precursors, which
    carry a charge. DIA: wide windows (typically 2-50 Th), a small set of
    window centres that repeats every cycle, no charge. PRM: narrow windows,
    a small repeating set of targets."""
    ms2 = run.ms2
    if not ms2:
        return [Evidence("acquisition", "MS1 only", 0.9, "metadata", "isolation windows",
                         "no MS2 spectra in the file")]
    widths = np.array([s.iso_width for s in ms2])
    have_w = widths > 0
    centers = np.round(np.array([s.iso_center for s in ms2]), 1)
    n_unique = len(set(centers.tolist()))
    uniq_ratio = n_unique / len(ms2)
    charged = float(np.mean([s.charge > 0 for s in ms2]))
    med_w = float(np.median(widths[have_w])) if have_w.any() else 0.0
    stats = {"median_isolation_width": round(med_w, 3), "unique_window_centres": n_unique,
             "ms2": len(ms2), "fraction_with_charge": round(charged, 3)}
    # Repetition: in DIA every centre is hit about once per cycle.
    repeats = Counter(centers.tolist())
    rep_med = float(np.median(list(repeats.values())))
    stats["median_repeats_per_centre"] = rep_med
    ev = []
    if not have_w.any():
        mode, conf, why = ("DDA" if uniq_ratio > 0.5 else "DIA"), 0.5, "no isolation windows recorded"
    elif med_w >= 4.0 and uniq_ratio < 0.2:
        mode, conf, why = "DIA", 0.97, f"wide windows (median {med_w:.1f} Th) on {n_unique} repeating centres"
    elif med_w >= 4.0 and uniq_ratio > 0.5:
        # older ion traps (LTQ, LTQ FT) isolate DDA precursors with 3-4 Th windows
        mode, conf, why = "DDA", 0.7, \
            f"wide windows (median {med_w:.1f} Th) on {n_unique} mostly distinct precursors"
    elif med_w >= 4.0:
        mode, conf, why = "DIA", 0.75, f"wide windows (median {med_w:.1f} Th), but many distinct centres"
    elif uniq_ratio < 0.05 and rep_med >= 5:
        mode, conf, why = "PRM", 0.8, f"narrow windows on {n_unique} repeating targets"
    else:
        mode, conf, why = "DDA", 0.95 if uniq_ratio > 0.3 else 0.7, \
            f"narrow windows (median {med_w:.2f} Th) on {n_unique} distinct precursors"
    if mode == "DIA" and any(s.ion_mobility for s in ms2):
        mode = "DIA (diaPASEF)"
    ev.append(Evidence("acquisition", mode, conf, "metadata", "isolation windows", why, stats))
    if mode.startswith("DIA") and have_w.any():
        lo = np.array([s.iso_center - s.iso_lower for s in ms2])
        hi = np.array([s.iso_center + s.iso_upper for s in ms2])
        win = sorted(set(zip(np.round(lo, 1).tolist(), np.round(hi, 1).tolist())))
        overlap = any(win[i + 1][0] < win[i][1] - 0.5 for i in range(len(win) - 1))
        ev.append(Evidence("dia_windows", {"n_windows": len(win), "min_mz": win[0][0],
                                           "max_mz": win[-1][1], "median_width": med_w,
                                           "overlapping": overlap,
                                           "windows": win[:200]},
                           0.95, "metadata", "isolation windows",
                           f"{len(win)} windows from {win[0][0]} to {win[-1][1]} Th"
                           + (" (overlapping/staggered)" if overlap else "")))
    # Charge range from DDA precursors. DIA windows carry no real charge (converters write a placeholder 2+, which
    # made every Orbitrap DIA run 2..2 and halved the PSMs on its pseudo-spectra); cli re-measures the charges on
    # DIAspeXtract pseudo-spectra with charge_range()
    if mode.startswith("DIA"):
        ev.append(Evidence("charge_range", [2, 4], 0.5, "default", "DIA default",
                           "DIA spectra carry no charge; 2..4 is the usual range"))
    else:
        ev += charge_range(ms2)
    return ev


def charge_range(ms2) -> list:
    """precursor charge range covering 99% of the assigned charges (DDA or pseudo-MS/MS spectra)"""
    charges = np.array([s.charge for s in ms2 if s.charge > 0])
    if charges.size < 20:
        return []
    lo, hi = int(np.percentile(charges, 0.5)), int(np.percentile(charges, 99.5))
    return [Evidence("charge_range", [max(1, lo), max(lo, hi)], 0.85, "data",
                     "precursor charges", f"99% of assigned precursor charges lie in {lo}..{hi}",
                     {"histogram": {int(k): int(v) for k, v in Counter(charges.tolist()).items()}})]
