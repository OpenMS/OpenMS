"""Hypothesis tests by mini-search (tagsearch.score_search).

Each test searches the same subset of spectra under competing settings and
compares the number of PSMs at 1% FDR. The setting that explains clearly more
spectra wins. These are the database-backed counterparts of the
diagnostic-ion detectors, and they work for things no diagnostic ion shows
(organism, protease, Cys alkylation on low-resolution data).
"""

from __future__ import annotations

import math
import time

import numpy as np

from . import tagsearch as ts
from .evidence import Evidence


def _psms(scans, peaks, index, prec_ppm, frag_tol, frag_ppm, mods, fdr=0.01):
    ps = ts.score_search(scans, index, prec_ppm, frag_tol, frag_ppm, mods, peaks=peaks)
    ok, thr = ts.fdr_filter(ps, fdr)
    return ok, len(ps)


def _margin_conf(best, second):
    """Confidence from the PSM count margin of the winner over the runner-up."""
    if best <= 0:
        return 0.0
    ratio = (best + 1) / (second + 1)
    return float(min(0.95, 0.4 + 0.25 * math.log2(ratio))) if ratio > 1 else 0.3


def fragment_errors(scans, peaks, psms, mods, frag_tol, frag_ppm):
    errs = []
    for k, _sc, _d, seq, _e in psms:
        th = ts._fragments(seq, mods, 1)
        pk = peaks[k]
        j = np.clip(np.searchsorted(pk, th), 1, pk.size - 1)
        cand = np.stack([pk[j - 1], pk[j]], 1)
        near = cand[np.arange(th.size), np.argmin(np.abs(cand - th[:, None]), 1)]
        d = (near - th) / th * 1e6 if frag_ppm else near - th
        errs.extend(d[np.abs(d) <= frag_tol].tolist())
    return np.asarray(errs)


class HypothesisTester:
    def __init__(self, run, frag_tol, frag_ppm, prec_ppm=20.0, n_spectra=2000,
                 mods: ts.Mods | None = None, log=print):
        self.run = run
        self.scans = ts.select_scans(run, n_spectra)
        self.peaks = [ts.prepare_peaks(s) for s in self.scans]
        self.frag_tol, self.frag_ppm, self.prec_ppm = frag_tol, frag_ppm, prec_ppm
        self.mods = mods or ts.Mods()
        self.log = log
        self.best_db = None       # (name, path)
        self.best_psms = []

    # -- organism -------------------------------------------------------------
    def organism(self, panel: dict, missed=1) -> list:
        res = {}
        for name, path in panel.items():
            t = time.time()
            idx = ts.build_index(name, ts.read_fasta(path), self.mods, missed=missed)
            ok, n = _psms(self.scans, self.peaks, idx, self.prec_ppm, self.frag_tol,
                          self.frag_ppm, self.mods)
            res[name] = {"psms": len(ok), "scans": {p[0] for p in ok}, "peptides": len(idx.seqs),
                         "seconds": round(time.time() - t, 1), "_ok": ok, "path": path}
            self.log(f"  organism {name}: {len(ok)} PSMs at 1% FDR ({len(idx.seqs)} peptides, "
                     f"{res[name]['seconds']} s)")
        if not res:
            return []
        ranked = sorted(res, key=lambda k: -res[k]["psms"])
        best = ranked[0]
        stats = {k: {kk: vv for kk, vv in v.items() if kk not in ("scans", "_ok", "path")}
                 for k, v in res.items()}
        # spectra a database explains that the winner does not (shared spike-ins,
        # contaminants and conserved proteins make the raw counts overlap)
        others = set().union(*(res[o]["scans"] for o in ranked[1:])) if len(ranked) > 1 else set()
        for k in ranked:
            if k == best:
                stats[k]["unique_vs_best"] = len(res[best]["scans"] - others)
            else:
                stats[k]["unique_vs_best"] = len(res[k]["scans"] - res[best]["scans"])
        n = len(self.scans)
        top = res[best]["psms"]
        second = res[ranked[1]]["psms"] if len(ranked) > 1 else 0
        ev = []
        if top < max(10, 0.01 * n):
            ev.append(Evidence("organism", "not in panel", 0.5, "search", "panel mini-search",
                               f"no panel proteome explains more than {top} of {n} spectra at 1% FDR; "
                               "the organism is probably not in the panel (or the settings are far off)",
                               stats))
            return ev
        conf = _margin_conf(top, second)
        detail = (f"{best} explains {top} of {n} spectra at 1% FDR "
                  f"({top / n:.0%}); next: {ranked[1]} with {second}" if len(ranked) > 1 else
                  f"{best} explains {top} of {n} spectra at 1% FDR")
        ev.append(Evidence("organism", best, conf, "search", "panel mini-search", detail, stats))
        # a second organism with many spectra of its own: a mixed sample
        for k in ranked[1:]:
            u = stats[k]["unique_vs_best"]
            if u >= max(10, 0.2 * top):
                ev.append(Evidence("organism", f"{best} + {k}", 0.5, "search", "panel mini-search",
                                   f"{k} explains {u} spectra that {best} does not: mixed sample, "
                                   "spike-ins or xenograft?", stats))
        self.best_db = (best, res[best]["path"])
        self.best_psms = res[best]["_ok"]
        return ev

    def _db(self):
        return self.best_db

    # -- accuracy -------------------------------------------------------------
    def accuracy(self) -> list:
        ev = []
        ok = self.best_psms
        if len(ok) < 20:
            return ev
        errs = np.array([p[4] for p in ok])
        med = float(np.median(errs))
        sd = float(1.4826 * np.median(np.abs(errs - med)))
        rec = max(3.0, math.ceil(abs(med) + 4 * sd))
        conf = min(0.95, 0.5 + 0.1 * math.log10(len(ok)) * 2)
        ev.append(Evidence("precursor_offset", {"value": round(med, 2), "unit": "ppm"}, conf, "search",
                           "mini-search PSMs", f"median precursor error of {len(ok)} PSMs: {med:+.2f} ppm"))
        ev.append(Evidence("precursor_tolerance", {"value": rec, "unit": "ppm"}, conf, "search",
                           "mini-search PSMs",
                           f"{len(ok)} PSMs at 1% FDR: offset {med:+.2f} ppm, robust SD {sd:.2f} ppm; "
                           f"|offset| + 4 SD = {rec} ppm", {"median": med, "sd": sd, "n": len(ok)}))
        fe = fragment_errors(self.scans, self.peaks, ok, self.mods, self.frag_tol, self.frag_ppm)
        if fe.size >= 100:
            fmed = float(np.median(fe))
            fsd = float(1.4826 * np.median(np.abs(fe - fmed)))
            unit = "ppm" if self.frag_ppm else "Da"
            step = 1.0 if self.frag_ppm else 0.05
            frec = abs(fmed) + 4 * fsd
            frec = max(5.0 if self.frag_ppm else 0.2, math.ceil(frec / step) * step)
            ev.append(Evidence("fragment_offset", {"value": round(fmed, 4), "unit": unit}, conf, "search",
                               "mini-search PSMs", f"median fragment error {fmed:+.3g} {unit}"))
            ev.append(Evidence("fragment_tolerance", {"value": round(frec, 3), "unit": unit}, conf * 0.9,
                               "search", "mini-search PSMs",
                               f"{fe.size} matched fragments: offset {fmed:+.3g}, SD {fsd:.3g} {unit}; "
                               f"|offset| + 4 SD = {frec:.3g} {unit} (bounded by the {self.frag_tol} "
                               f"{unit} pre-search tolerance)", {"median": fmed, "sd": fsd, "n": int(fe.size)}))
        return ev

    # -- competing settings on the best database -------------------------------
    def _compete(self, parameter, variants, missed=1, enzyme="Trypsin", method=""):
        if not self.best_db:
            return []
        name, path = self.best_db
        prots = list(ts.read_fasta(path))
        counts = {}
        for label, mods, enz in variants:
            idx = ts.build_index(name, prots, mods, enzyme=enz or enzyme, missed=missed)
            ok, _ = _psms(self.scans, self.peaks, idx, self.prec_ppm, self.frag_tol, self.frag_ppm, mods)
            counts[label] = {"psms": len(ok), "with_C": sum(1 for p in ok if "C" in p[3])}
            self.log(f"  {parameter} {label}: {counts[label]}")
        return counts

    def alkylation(self) -> list:
        variants = []
        for name, dm in ts.CYS_MODS.items():
            m = ts.Mods(cys=dm, k=self.mods.k, nterm=self.mods.nterm)
            variants.append((name, m, None))
        counts = self._compete("cys_alkylation", variants)
        if not counts:
            return []
        ranked = sorted(counts, key=lambda k: -counts[k]["with_C"])
        best, second = ranked[0], ranked[1]
        b, s = counts[best]["with_C"], counts[second]["with_C"]
        if b < 5:
            return [Evidence("cys_alkylation", "unknown", 0.2, "search", "alkylation competition",
                             f"too few Cys peptides identified ({b}) to decide", {"counts": counts})]
        return [Evidence("cys_alkylation", best, _margin_conf(b, s), "search", "alkylation competition",
                         f"{best}: {b} Cys-containing PSMs at 1% FDR; next {second}: {s}",
                         {"counts": counts})]

    def enzyme(self) -> list:
        variants = [(e, self.mods, e) for e in ("Trypsin", "Lys-C", "Arg-C", "Chymotrypsin", "Glu-C", "Asp-N")]
        counts = self._compete("enzyme", variants)
        if not counts:
            return []
        ranked = sorted(counts, key=lambda k: -counts[k]["psms"])
        b, s = counts[ranked[0]]["psms"], counts[ranked[1]]["psms"]
        val = ranked[0]
        # Trypsin and Lys-C/Arg-C share many peptides; call the broader rule
        # unless the narrower one is clearly better
        return [Evidence("enzyme", val, _margin_conf(b, s), "search", "protease competition",
                         f"{val}: {b} PSMs at 1% FDR; next {ranked[1]}: {s}",
                         {"counts": {k: v["psms"] for k, v in counts.items()}})]
