"""Mass accuracy from the data alone, before any database search.

Two kinds of signal, both database free:

* Precision (spread) from repeated measurements of the same ion
  (Param-Medic, May et al. J Proteome Res 2017): pairs of MS2 spectra of the
  same precursor (same charge, precursor m/z within 50 ppm, >= 20 of the top
  40 fragments shared) give precursor and fragment m/z differences; a
  Gaussian + uniform mixture fitted by EM gives sigma_pair, and the
  per-measurement sigma is sigma_pair / sqrt(2). For DIA the same window in
  neighbouring cycles plays the role of the repeated precursor, and
  neighbouring MS1 scans give the MS1 precision in every mode.

* Accuracy (offset) from ions of known mass that are present in most runs:
  the polysiloxane lock-mass ions in MS1 (445.12003, 371.10123, 519.13882)
  and, in HCD MS2 spectra, immonium ions (L/I 86.09643, F 120.08078,
  Y 136.07569, P 70.06512, H 110.07127) and tryptic y1 ions (K 147.11280,
  R 175.11895). Param-Medic cannot see the offset; these probes can, and
  a ProSE calibration pre-search (prose_runner.py) confirms it.

The recommended tolerance is |offset| + 4 sigma (the 99.99% range of a
normal error), rounded up and kept above an instrument-dependent floor.
"""

from __future__ import annotations

import math

import numpy as np

from .evidence import Evidence

LOCK_MASSES = (445.120025, 371.101233, 519.138815)
IMMONIUM = {"L/I": 86.096426, "F": 120.080776, "Y": 136.075690, "P": 70.065126,
            "H": 110.071273, "W": 159.091675}
Y1 = {"y1-K": 147.112804, "y1-R": 175.118952}


def fit_gauss_uniform(x: np.ndarray, half_width: float, iters: int = 300):
    """EM for a Gaussian + uniform(-half_width, half_width) mixture.
    Returns (mu, sigma, gaussian_fraction) or None."""
    x = np.asarray(x, dtype=float)
    x = x[np.abs(x) <= half_width]
    if x.size < 30:
        return None
    mu = float(np.median(x))
    sigma = float(1.4826 * np.median(np.abs(x - mu))) or half_width / 10
    pi = 0.5
    u = 1.0 / (2 * half_width)
    for _ in range(iters):
        g = pi * np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))
        r = g / (g + (1 - pi) * u)
        rs = r.sum()
        if rs < 5:
            return None
        mu_n = float((r * x).sum() / rs)
        sigma_n = float(math.sqrt((r * (x - mu_n) ** 2).sum() / rs)) or 1e-9
        pi_n = float(rs / x.size)
        if abs(mu_n - mu) < 1e-9 and abs(sigma_n - sigma) < 1e-9:
            break
        mu, sigma, pi = mu_n, sigma_n, pi_n
    return mu, sigma, pi


def _top_peaks(scan, n=40):
    if scan.mz.size == 0:
        return scan.mz
    idx = np.argsort(scan.intensity)[::-1][:n]
    return np.sort(scan.mz[idx])


def _match(a, b, tol, ppm):
    """Greedy nearest matching of two sorted m/z arrays; returns diffs of
    matched pairs (ppm or Th)."""
    if a.size == 0 or b.size == 0:
        return np.empty(0)
    j = np.searchsorted(b, a)
    j = np.clip(j, 1, b.size - 1) if b.size > 1 else np.zeros_like(j)
    cand = np.stack([b[np.maximum(j - 1, 0)], b[j]], axis=1)
    near = cand[np.arange(a.size), np.argmin(np.abs(cand - a[:, None]), axis=1)]
    d = near - a
    if ppm:
        d = d / a * 1e6
    return d[np.abs(d) <= tol]


def _pairs_dda(ms2, prec_ppm=50.0, rt_window=600.0, max_pairs=5000):
    """Candidate pairs of MS2 spectra of the same precursor ion."""
    s = sorted((x for x in ms2 if x.precursor_mz > 0), key=lambda x: x.precursor_mz)
    mzs = np.array([x.precursor_mz for x in s])
    pairs = []
    for i in range(len(s)):
        hi = mzs[i] * (1 + prec_ppm * 1e-6)
        j = i + 1
        while j < len(s) and mzs[j] <= hi:
            a, b = s[i], s[j]
            if (a.charge == b.charge or 0 in (a.charge, b.charge)) and abs(a.rt - b.rt) <= rt_window:
                pairs.append((a, b))
            j += 1
        if len(pairs) >= max_pairs * 4:
            break
    return pairs


def _pairs_dia(ms2):
    """The same isolation window in consecutive cycles."""
    by_win = {}
    for x in ms2:
        by_win.setdefault(round(x.iso_center, 1), []).append(x)
    pairs = []
    for scans in by_win.values():
        scans.sort(key=lambda x: x.rt)
        pairs.extend(zip(scans[:-1], scans[1:]))
    return pairs


def _recommend(offset, sigma, floor, step):
    tol = abs(offset) + 4.0 * sigma
    tol = max(floor, math.ceil(tol / step) * step)
    return round(tol, 4)


def detect_tolerances(run, acquisition="DDA", ms2_high_res=True, ms1_high_res=True) -> list:
    ev = []
    dia = acquisition.startswith("DIA")
    frag_unit = "ppm" if ms2_high_res else "Da"
    frag_tol_match = 30.0 if ms2_high_res else 0.6     # pairing tolerance
    stats = {}

    # --- fragment precision (and precursor precision for DDA) from repeated ions
    pairs = _pairs_dia(run.ms2) if dia else _pairs_dda(run.ms2)
    prec_d, frag_d, n_good = [], [], 0
    for a, b in pairs:
        pa, pb = _top_peaks(a), _top_peaks(b)
        d = _match(pa, pb, frag_tol_match, ms2_high_res)
        if d.size >= 20:
            n_good += 1
            frag_d.append(d)
            if not dia:
                prec_d.append((b.precursor_mz - a.precursor_mz) / a.precursor_mz * 1e6)
            if n_good >= 5000:
                break
    stats["pairs_candidate"] = len(pairs)
    stats["pairs_same_ion"] = n_good

    if frag_d:
        fd = np.concatenate(frag_d)
        fit = fit_gauss_uniform(fd, frag_tol_match)
        if fit:
            sigma = fit[1] / math.sqrt(2)
            stats["fragment_sigma"] = round(sigma, 4)
            rec = _recommend(0.0, sigma, 5.0 if ms2_high_res else 0.2, 1.0 if ms2_high_res else 0.05)
            conf = min(0.9, 0.4 + 0.1 * math.log10(max(10, n_good)))
            ev.append(Evidence("fragment_tolerance", {"value": rec, "unit": frag_unit}, conf, "data",
                               "repeated-ion pairs",
                               f"fragment precision sigma {sigma:.3g} {frag_unit} from {n_good} "
                               f"{'window repeats' if dia else 'same-precursor spectrum pairs'}; "
                               f"tolerance 4 sigma = {rec} {frag_unit} (offset not included)",
                               {"sigma": sigma, "pairs": n_good, "gauss_fraction": fit[2]}))
    if prec_d:
        fit = fit_gauss_uniform(np.array(prec_d), 50.0)
        if fit:
            sigma = fit[1] / math.sqrt(2)
            rec = _recommend(0.0, sigma, 3.0 if ms1_high_res else 0.3, 1.0)
            conf = min(0.85, 0.35 + 0.1 * math.log10(max(10, len(prec_d))))
            ev.append(Evidence("precursor_tolerance", {"value": rec, "unit": "ppm"}, conf, "data",
                               "repeated-ion pairs",
                               f"precursor precision sigma {sigma:.2f} ppm from {len(prec_d)} pairs; "
                               f"tolerance 4 sigma = {rec} ppm (offset not included)",
                               {"sigma": sigma, "pairs": len(prec_d), "gauss_fraction": fit[2]}))

    # --- MS1 precision from neighbouring survey scans (every acquisition mode)
    if ms1_high_res and len(run.ms1) >= 10:
        ms1 = sorted(run.ms1, key=lambda x: x.rt)
        d = [_match(_top_peaks(a, 100), _top_peaks(b, 100), 20.0, True)
             for a, b in zip(ms1[:-1], ms1[1:])]
        d = np.concatenate(d) if d else np.empty(0)
        fit = fit_gauss_uniform(d, 20.0)
        if fit:
            sigma = fit[1] / math.sqrt(2)
            stats["ms1_sigma_ppm"] = round(sigma, 3)
            ev.append(Evidence("precursor_tolerance",
                               {"value": _recommend(0.0, sigma, 3.0, 1.0), "unit": "ppm"},
                               0.55 if dia else 0.45, "data", "MS1 scan-to-scan",
                               f"MS1 peak m/z precision sigma {sigma:.2f} ppm between neighbouring "
                               f"survey scans ({d.size} matched peaks)",
                               {"sigma": sigma, "matched_peaks": int(d.size)}))

    # --- offsets from ions of known mass
    if ms1_high_res and run.ms1:
        errs = []
        for s in run.ms1:
            for lm in LOCK_MASSES:
                if s.mz.size and s.scan_low <= lm <= s.scan_high:
                    i = np.searchsorted(s.mz, lm)
                    for k in (i - 1, i):
                        if 0 <= k < s.mz.size:
                            e = (s.mz[k] - lm) / lm * 1e6
                            if abs(e) <= 15:
                                errs.append(e)
        frac = len(errs) / max(1, len(run.ms1))
        if len(errs) >= 30 and frac >= 0.2:
            fit = fit_gauss_uniform(np.array(errs), 15.0)
            if fit and fit[2] > 0.3:
                ev.append(Evidence("precursor_offset", {"value": round(fit[0], 2), "unit": "ppm"},
                                   min(0.8, 0.4 + frac / 2), "data", "lock-mass probe",
                                   f"polysiloxane ions found in {frac:.0%} of MS1 scans, "
                                   f"median error {fit[0]:+.2f} ppm, sigma {fit[1]:.2f} ppm",
                                   {"n": len(errs), "mu": fit[0], "sigma": fit[1]}))
    if ms2_high_res and run.ms2:
        errs = []
        probes = list(IMMONIUM.values()) + list(Y1.values())
        for s in run.ms2:
            if s.mz.size == 0 or s.scan_low > 140:
                continue
            base = s.intensity.max()
            for p in probes:
                i = np.searchsorted(s.mz, p)
                for k in (i - 1, i):
                    if 0 <= k < s.mz.size and s.intensity[k] >= 0.01 * base:
                        e = (s.mz[k] - p) / p * 1e6
                        if abs(e) <= 30:
                            errs.append(e)
        if len(errs) >= 50:
            fit = fit_gauss_uniform(np.array(errs), 30.0)
            if fit and fit[2] > 0.3:
                rec = _recommend(fit[0], fit[1], 5.0, 1.0)
                ev.append(Evidence("fragment_offset", {"value": round(fit[0], 2), "unit": "ppm"},
                                   0.7, "data", "immonium/y1 probe",
                                   f"{len(errs)} immonium and y1 peaks: offset {fit[0]:+.2f} ppm, "
                                   f"sigma {fit[1]:.2f} ppm", {"n": len(errs), "mu": fit[0], "sigma": fit[1]}))
                ev.append(Evidence("fragment_tolerance", {"value": rec, "unit": "ppm"}, 0.75, "data",
                                   "immonium/y1 probe",
                                   f"absolute fragment error from ions of known mass: "
                                   f"|offset| + 4 sigma = {rec} ppm", {"n": len(errs)}))
    return ev
