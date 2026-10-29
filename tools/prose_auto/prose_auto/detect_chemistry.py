"""Labels, enrichment and alkylation from diagnostic ions, before any search.

Each test counts the fraction of MSn spectra with a peak at a diagnostic m/z
and compares it with the same count at control m/z values that no peptide ion
prefers. Enrichment = f(diagnostic) / f(control); a label or enrichment is
called when the diagnostic ions are frequent and strongly enriched.

* Isobaric labels: TMT (0/2/6/10/11), TMTpro (16/18), iTRAQ 4/8 reporter ions.
  TMT 10/11 vs 6 needs the 6 mDa N/C isotopologue split (high-res MS2 only).
* SILAC: MS1 peak pairs spaced by Lys8 (8.0142 Da), Arg10 (10.0083 Da),
  Arg6 (6.0201 Da), Lys4 (4.0251 Da, same as dimethyl d4) at charge 2 and 3,
  against the same spacing shifted by whole 13C steps, on many distinct
  precursors, with the heavy peak separate from the light isotope envelope.
* Phospho enrichment: neutral loss of H3PO4 (97.9769 Da) from the precursor,
  phosphotyrosine immonium 216.0420.
* Glyco: HexNAc oxonium ions 204.0867, 138.0545, 366.1395.
* Acetyl-lysine: 126.0913 / 143.1179 immonium ions.
* Cys alkylation: the cysteine immonium ion shifted by the alkylating group
  (CAM 133.0430, propionamide 147.0587, NEM 201.0692, MMTS 122.0093,
  carboxymethyl 134.0270, free Cys 76.0216). This is weak evidence (Cys
  immonium ions are small); search.py's mini-search competition is the
  stronger test.
"""

from __future__ import annotations

import numpy as np

from .evidence import Evidence

H3PO4 = 97.976896

REPORTERS = {
    "TMT126": 126.127726, "TMT127N": 127.124761, "TMT127C": 127.131081,
    "TMT128N": 128.128116, "TMT128C": 128.134436, "TMT129N": 129.131471,
    "TMT129C": 129.137790, "TMT130N": 130.134825, "TMT130C": 130.141145,
    "TMT131N": 131.138180, "TMT131C": 131.144499, "TMT132N": 132.141535,
    "TMT132C": 132.147855, "TMT133N": 133.144890, "TMT133C": 133.151210,
    "TMT134N": 134.148245, "TMT134C": 134.154565, "TMT135N": 135.151600,
    "iTRAQ113": 113.107873, "iTRAQ114": 114.111228, "iTRAQ115": 115.108263,
    "iTRAQ116": 116.111618, "iTRAQ117": 117.114973, "iTRAQ118": 118.112008,
    "iTRAQ119": 119.115363, "iTRAQ121": 121.122072,
}

DIAG = {
    "pY_immonium": 216.042040,
    "HexNAc_204": 204.086649, "HexNAc_138": 138.054500, "HexNAcHex_366": 366.139472,
    "AcK_126": 126.091340, "AcK_143": 143.117890,
}

CYS_IMMONIUM = {
    "Carbamidomethyl (C)": 133.043010,
    "Propionamide (C)": 147.058660,
    "Nethylmaleimide (C)": 201.069240,
    "Methylthio (C)": 122.009280,
    "Carboxymethyl (C)": 134.027030,
    "none (free Cys)": 76.021560,
}

SILAC_DELTAS = {"Lys8": 8.014199, "Arg10": 10.008269, "Arg6": 6.020129, "Lys4/dimethyl-d4": 4.025107}
CONTROL_SHIFTS = (0.27, 0.41, -0.33, 0.58, -0.46)


def _has_peak(scan, mz, tol, ppm, min_rel=0.0):
    if scan.mz.size == 0 or not (scan.scan_low - 0.5 <= mz <= scan.scan_high + 0.5):
        return None  # outside the scan range: no information
    t = mz * tol * 1e-6 if ppm else tol
    i = np.searchsorted(scan.mz, mz - t)
    j = np.searchsorted(scan.mz, mz + t)
    if j <= i:
        return False
    if min_rel > 0:
        return scan.intensity[i:j].max() >= min_rel * scan.intensity.max()
    return True


def _frac(scans, mz, tol, ppm, min_rel=0.0):
    hits = [h for h in (_has_peak(s, mz, tol, ppm, min_rel) for s in scans) if h is not None]
    if not hits:
        return None, 0
    return float(np.mean(hits)), len(hits)


def _enrichment(scans, mz, tol, ppm, min_rel=0.0):
    f, n = _frac(scans, mz, tol, ppm, min_rel)
    if f is None:
        return None
    ctrl = [c for c in (_frac(scans, mz + d, tol, ppm, min_rel)[0] for d in CONTROL_SHIFTS) if c is not None]
    bg = float(np.mean(ctrl)) if ctrl else 0.0
    return {"fraction": round(f, 4), "background": round(bg, 4),
            "enrichment": round(f / max(bg, 0.002), 2), "n": n}


def detect_labels(run, ms2_high_res=True) -> list:
    # Reporter ions are in MS3 for SPS-MS3 methods, in MS2 otherwise.
    scans = run.ms3 if len(run.ms3) > 0.2 * max(1, len(run.ms2)) else run.ms2
    level = "MS3" if scans is run.ms3 else "MS2"
    hr = ms2_high_res if level == "MS2" else True
    tol, ppm = (15.0, True) if hr else (0.3, False)
    res = {k: _enrichment(scans, v, tol, ppm) for k, v in REPORTERS.items()}
    covered = {k: r for k, r in res.items() if r is not None}
    if not scans:
        return detect_silac(run)
    if not covered:
        low = float(np.median([s.scan_low for s in scans]))
        return [Evidence("label", "unknown", 0.2, "data", "reporter ions",
                         f"reporter region (113-135 Th) is outside the {level} scan range "
                         f"(scans start at {low:.0f} Th); isobaric labels cannot be "
                         "seen in these spectra")] + detect_silac(run)

    def present(keys):
        rs = [covered.get(k) for k in keys]
        rs = [r for r in rs if r]
        if len(rs) < len(keys):
            return False, rs
        return all(r["fraction"] >= 0.25 and r["enrichment"] >= 4 for r in rs), rs

    tmt6 = ["TMT126", "TMT127N", "TMT128C", "TMT129N", "TMT130C", "TMT131N"]
    tmt10_extra = ["TMT127C", "TMT128N", "TMT129C", "TMT130N"]
    tmtpro_extra = ["TMT132N", "TMT133C", "TMT134N"]
    itraq4 = ["iTRAQ114", "iTRAQ115", "iTRAQ116", "iTRAQ117"]
    itraq8_extra = ["iTRAQ113", "iTRAQ118", "iTRAQ119", "iTRAQ121"]
    stats = {"level": level, "reporters": covered}

    call, conf, why = "label free", 0.0, ""
    ok_tmt, rs = present(tmt6)
    ok_126, _ = present(["TMT126"])
    if ok_tmt or ok_126:
        n_pro = sum(1 for k in tmtpro_extra if covered.get(k) and covered[k]["fraction"] >= 0.2
                    and covered[k]["enrichment"] >= 4)
        if n_pro >= 2:
            call = "TMTpro (16/18plex)"
        elif ok_tmt:
            # 10/11plex needs both members of the N/C pairs; at low res they merge
            n10 = sum(1 for k in tmt10_extra if covered.get(k) and covered[k]["enrichment"] >= 4
                      and covered[k]["fraction"] >= 0.2)
            call = "TMT10/11plex" if (hr and n10 >= 3) else ("TMT6plex" if hr else "TMT6/10/11plex")
        else:
            # 127C and 129C are also the 13C isotopes of 126 and 128C, so 126 + 127C alone is not TMT2plex
            # when channels only TMT6 and up carry (128C, 130C, 131N) are there (PXD000001: 6plex, 127N/129N
            # unused, called TMT2plex and searched with the wrong reagent mass, 0 PSMs)
            n6 = sum(present([k])[0] for k in ("TMT128C", "TMT130C", "TMT131N"))
            ok_2, _ = present(["TMT126", "TMT127C"])
            call = (("TMT6plex" if hr else "TMT6/10/11plex") if n6 >= 2 else
                    "TMT2plex" if ok_2 else "TMT0 / TMT (126 only)")
        top = covered["TMT126"]
        conf = min(0.97, 0.6 + 0.4 * top["fraction"])
        why = f"TMT 126 reporter in {top['fraction']:.0%} of {level} spectra ({top['enrichment']}x background)"
    else:
        ok4, rs4 = present(itraq4)
        if ok4:
            ok8 = sum(1 for k in itraq8_extra if covered.get(k) and covered[k]["enrichment"] >= 4) >= 3
            call = "iTRAQ8plex" if ok8 else "iTRAQ4plex"
            f = float(np.mean([r["fraction"] for r in rs4]))
            conf = min(0.95, 0.6 + 0.4 * f)
            why = f"iTRAQ reporters in {f:.0%} of {level} spectra"
    ev = []
    if call != "label free":
        ev.append(Evidence("label", call, conf, "data", "reporter ions", why, stats))
    else:
        mx = max(covered.values(), key=lambda r: r["fraction"])
        ev.append(Evidence("label", "no isobaric label", 0.8, "data", "reporter ions",
                           f"no reporter ion set is enriched (most frequent reporter in "
                           f"{mx['fraction']:.0%} of spectra)", stats))
    ev.extend(detect_silac(run))
    return ev


C13 = 1.0033548


def _silac_controls(d):
    """Control spacings for a SILAC delta: the same delta shifted by whole 13C spacings, so the control has the
    mass defect of the delta. Peptide masses cluster near integer spacings, so fractional controls (7.3, 9.27 Da)
    make every near-integer delta look enriched; on label-free runs Arg6 came out 5-22x 'enriched' that way.
    Shifts that land on another SILAC delta are left out."""
    out = []
    for k in (-3, -2, -1, 1, 2, 3):
        c = d + k * C13
        if c > 1.5 and all(abs(c - s) > 0.03 for s in SILAC_DELTAS.values()):
            out.append(c)
    return out


def _pairs(mz, it, dm, z):
    """Indices of light peaks with a heavy partner dm Th above: comparable intensity (0.15-7x; Arg6 spaces like the
    light peptide's M+6 isotope, which is far weaker), a light M+1 isotope (a peptide, not noise), and a gap below
    the heavy peak (the heavy is not just the next peak of the light isotope envelope)."""
    def find(target):
        j = np.clip(np.searchsorted(mz, target), 1, mz.size - 1)
        j = np.where(np.abs(mz[j - 1] - target) < np.abs(mz[j] - target), j - 1, j)
        return j, np.abs(mz[j] - target) <= target * 5e-6
    j, ok = find(mz + dm)
    ratio = it[j] / it
    k1, ok1 = find(mz + C13 / z)
    g, okg = find(mz + dm - C13 / z)
    gap = ~okg | (it[g] < 0.3 * it[j]) | (np.abs(dm - C13 / z) < 1e-3)
    return np.nonzero(ok & (ratio >= 0.15) & (ratio <= 7) & ok1 & gap)[0]


def detect_silac(run) -> list:
    if len(run.ms1) < 20:
        return []
    counts = {k: 0 for k in SILAC_DELTAS}
    ctrl = {k: [0] * len(_silac_controls(d)) for k, d in SILAC_DELTAS.items()}
    light = {k: set() for k in SILAC_DELTAS}   # distinct light precursors (m/z, z) carrying a pair
    for s in run.ms1:
        if s.mz.size < 10:
            continue
        idx = np.argsort(s.intensity)[::-1][:150]
        order = np.argsort(s.mz[idx])
        mz, it = s.mz[idx][order], s.intensity[idx][order]
        for z in (2, 3):
            for k, d in SILAC_DELTAS.items():
                hit = _pairs(mz, it, d / z, z)
                counts[k] += hit.size
                light[k].update((round(float(m), 2), z) for m in mz[hit])
                for i, c in enumerate(_silac_controls(d)):
                    ctrl[k][i] += _pairs(mz, it, c / z, z).size
    bg = {k: max(1.0, float(np.mean(v))) for k, v in ctrl.items()}
    enr = {k: round(counts[k] / bg[k], 2) for k in counts}
    distinct = {k: len(v) for k, v in light.items()}
    best = max(enr, key=enr.get)
    stats = {"pair_counts": counts, "control_mean": bg, "enrichment": enr, "distinct_light": distinct}
    # a real SILAC sample: the delta beats its own isotope-shifted neighbours, on many distinct peptides
    if enr[best] >= 3 and counts[best] >= 50 and distinct[best] >= 30:
        name = {"Lys8": "SILAC (Lys8/Arg10)", "Arg10": "SILAC (Lys8/Arg10)",
                "Arg6": "SILAC (Arg6)", "Lys4/dimethyl-d4": "SILAC Lys4 or dimethyl (d4)"}[best]
        return [Evidence("label", name, min(0.9, 0.5 + enr[best] / 40), "data", "MS1 heavy/light pairs",
                         f"MS1 peak pairs {SILAC_DELTAS[best]:.4f} Da apart on {distinct[best]} precursors are "
                         f"{enr[best]}x more frequent than the same spacing shifted by whole 13C steps", stats)]
    return [Evidence("label", "no SILAC", 0.6, "data", "MS1 heavy/light pairs",
                     f"no SILAC spacing enriched (best {best}: {enr[best]}x)", stats)]


def detect_enrichment(run, ms2_high_res=True, dia=False) -> list:
    """dia: wide isolation windows co-fragment many precursors, so oxonium ions are everywhere (HexNAc 204 in
    34-43% of the spectra of two phospho DIA runs); no glyco call is made from DIA spectra."""
    tol, ppm = (15.0, True) if ms2_high_res else (0.4, False)
    ev = []
    # Phospho neutral loss from the precursor: needs the charge (DDA)
    nl, ctrl, n = 0, [0, 0], 0
    for s in run.ms2:
        if s.charge <= 0 or s.precursor_mz <= 0 or s.mz.size == 0:
            continue
        n += 1
        if _has_peak(s, s.precursor_mz - H3PO4 / s.charge, tol, ppm, 0.2):
            nl += 1
        for i, c in enumerate((91.5, 104.5)):
            if _has_peak(s, s.precursor_mz - c / s.charge, tol, ppm, 0.2):
                ctrl[i] += 1
    diag = {k: _enrichment(run.ms2, v, tol, ppm, 0.0) for k, v in DIAG.items()}
    stats = {"neutral_loss_fraction": round(nl / max(1, n), 4),
             "neutral_loss_control": round(np.mean(ctrl) / max(1, n), 4), "diagnostic": diag}
    calls = []
    if n >= 100:
        f, b = nl / n, np.mean(ctrl) / n
        # HCD loses H3PO4 less often than CID (5% of the spectra of an enriched Orbitrap Elite run), but the
        # control losses stay near zero on high-res MS2; unenriched runs reach the same fraction only with an
        # equally high control (ion-trap MS2, 0.4 Da windows)
        if f >= 0.03 and f >= 8 * max(b, 0.003):
            calls.append(("Phospho (STY)", min(0.9, 0.6 + 2 * f), f"{f:.1%} of MS2 spectra show a strong "
                          f"-98 Da neutral loss from the precursor (control {b:.2%})"))
    py = diag.get("pY_immonium")
    if py and py["fraction"] >= 0.1 and py["enrichment"] >= 5:
        calls.append(("Phospho (Y)", 0.6, f"pY immonium ion in {py['fraction']:.0%} of spectra"))
    g = diag.get("HexNAc_204") if not dia else None
    if g and g["fraction"] >= 0.15 and g["enrichment"] >= 5:
        calls.append(("Glyco (HexNAc oxonium)", min(0.9, 0.5 + g["fraction"]),
                      f"HexNAc oxonium 204.087 in {g['fraction']:.0%} of spectra ({g['enrichment']}x)"))
    a = diag.get("AcK_126") if ms2_high_res else None
    if a and a["fraction"] >= 0.15 and a["enrichment"] >= 5:
        calls.append(("Acetyl (K)", 0.6, f"acetyl-lysine immonium ion in {a['fraction']:.0%} of spectra"))
    if calls:
        for name, conf, why in calls:
            ev.append(Evidence("enrichment", name, conf, "data", "diagnostic ions", why, stats))
    else:
        # with >= 100 charged DDA spectra the neutral-loss test could have seen an enrichment; DIA spectra
        # (no precursor charge) leave the question to the records and the text
        ev.append(Evidence("enrichment", "none", 0.6 if n >= 100 else 0.3, "data", "diagnostic ions",
                           "no phospho, glyco or acetyl signature above background"
                           + ("" if n >= 100 else " (DIA or too few charged spectra: weak)"), stats))
    return ev


def detect_cys_alkylation(run, ms2_high_res=True) -> list:
    if not ms2_high_res:
        return [Evidence("cys_alkylation", "unknown", 0.1, "data", "Cys immonium",
                         "low-resolution MS2: Cys immonium ions cannot be resolved")]
    res = {k: _enrichment(run.ms2, v, 10.0, True, 0.0) for k, v in CYS_IMMONIUM.items()}
    res = {k: r for k, r in res.items() if r}
    if not res:
        return [Evidence("cys_alkylation", "unknown", 0.1, "data", "Cys immonium",
                         "Cys immonium region outside the MS2 scan range")]
    best = max(res, key=lambda k: res[k]["fraction"] * min(res[k]["enrichment"], 20))
    r = res[best]
    if r["fraction"] >= 0.03 and r["enrichment"] >= 3:
        return [Evidence("cys_alkylation", best, min(0.6, 0.3 + r["fraction"] * 3), "data",
                         "Cys immonium", f"{best} Cys immonium ion in {r['fraction']:.1%} of MS2 "
                         f"spectra ({r['enrichment']}x background); weak evidence, confirm by mini-search",
                         {"immonium": res})]
    return [Evidence("cys_alkylation", "unknown", 0.15, "data", "Cys immonium",
                     "no Cys immonium ion above background", {"immonium": res})]
