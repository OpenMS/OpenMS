"""Iterative tuning: start from prose-auto's settings and improve them with ProSE.

Three tiers, each cheaper than the next is expensive:

1. prose-auto proposes settings (detectors, mini-search, texts, consistency rules).
2. ProSE searches a top-TIC subset of the run with those settings. Its
   -summary_out YAML gives the target PSMs at the FDR threshold and the
   precursor and fragment error of the calibration PSMs; an optional open
   search (-out_mod_analysis_dir) gives the most frequent unexplained mass
   shifts. These diagnostics turn into a shortlist of candidate changes
   (widen or narrow a tolerance, add a variable modification seen as a mass
   shift, Trypsin/P, more missed cleavages, the runner-up organism, ...),
   each with a reason. Changes that would break a hard consistency rule
   (schema.RULES) are not proposed.
3. A chooser orders the shortlist: Jev (one `choice` question over the
   shortlist, with the diagnostics and the results of earlier steps as
   state), or without Jev the rule order. The first candidates are searched;
   a change is kept when it gains at least `min_gain` PSMs (or entrapment
   PSMs at 1% FDP when a scorer is given), and the loop repeats from the new
   best setting until nothing gains or the budget is used.

Every step is logged (situation, candidates, choice, result) to a JSONL file,
so the chooser can later be calibrated against what actually helped.

The scorer is pluggable: by default the PSM count from ProSE's summary; with
--entrap, the benchmark's `entrap.py eval` (PSMs at 1% entrapment FDP), run
against a database that contains the entrapment twins.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import subprocess
import sys
import tempfile
import time

from . import schema as S
from .prose_params import to_cli
from .spectra import write_subset_mzml

COMMON = {"Search:peptide:missed_cleavages": 2, "Search:modifications:variable_max_per_peptide": 2,
          "Search:rescore": "true"}

# mass shifts from the open search -> the ProSE variable modification they suggest
SHIFT_MODS = [
    (79.966, ["Phospho (S)", "Phospho (T)", "Phospho (Y)"]), (0.984, ["Deamidated (N)", "Deamidated (Q)"]),
    (42.011, ["Acetyl (Protein N-term)"]), (-17.027, ["Gln->pyro-Glu (N-term Q)"]),
    (43.006, ["Carbamyl (N-term)"]), (14.016, ["Methyl (K)"]), (114.043, ["GG (K)"]),
    (15.995, ["Oxidation (M)"]), (57.021, ["Carbamidomethyl (C)"]), (-57.021, []),
    (8.014, ["Label:13C(6)15N(2) (K)"]), (10.008, ["Label:13C(6)15N(4) (R)"]), (6.020, ["Label:13C(6) (R)"]),
]


def _get(P, k, d=None):
    return P.get("Search:" + k, d)


class ProSE:
    """Runs ProSE on one (subset) mzML and returns score + diagnostics."""

    def __init__(self, exe, database, threads=8, extra=(), entrap=None, log=print):
        self.exe, self.db, self.threads, self.extra, self.entrap, self.log = exe, database, threads, list(extra), entrap, log
        self.n = 0

    def __call__(self, mzml, P, open_search=False):
        import yaml
        self.n += 1
        with tempfile.TemporaryDirectory(prefix="prose_tune_") as tmp:
            out, summ = os.path.join(tmp, "out.idXML"), os.path.join(tmp, "summary.yaml")
            params = dict(P)
            if open_search:
                params.update({"Search:precursor:mass_tolerance_lower": 500.0,
                               "Search:precursor:mass_tolerance_upper": 500.0,
                               "Search:precursor:mass_tolerance_unit": "Da",
                               "Search:calibration:enabled": "false"})
            cmd = [self.exe, "-in", mzml, "-database", self.db, "-out_idxml", out, "-summary_out", summ,
                   "-threads", str(self.threads), "-no_summary"] + to_cli(params) + self.extra
            mods = None
            if open_search:
                mods = os.path.join(tmp, "mods")
                os.makedirs(mods)
                cmd += ["-out_mod_analysis_dir", mods]
            t = time.time()
            r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError((r.stderr or r.stdout)[-1500:])
            with open(summ) as fh:
                pf = ((yaml.safe_load(fh) or {}).get("per_file") or [{}])[0]
            # the summary's target_psms is unfiltered unless -Search:FDR:PSM is set, and filtering would cut the
            # q-value curve the entrapment estimate needs: count rank-1 targets at q <= 0.01 in the idXML instead
            res = {"psms": _psms_at(out, 0.01), "peptides": pf.get("unique_peptides"),
                   "ms2": pf.get("ms2_spectra"), "precursor_error": pf.get("precursor_error") or {},
                   "fragment_error": pf.get("fragment_error") or {}, "seconds": round(time.time() - t, 1)}
            if self.entrap and not open_search:
                ej = os.path.join(tmp, "entrap.json")
                subprocess.run([sys.executable, self.entrap, "eval", ej, self.db, out], check=True,
                               capture_output=True, text=True)
                with open(ej) as fh:
                    res["entrap_psms"] = json.load(fh).get("psm", {}).get("N_T_at_equal_fdp_1pct")
            res["score"] = res.get("entrap_psms") if res.get("entrap_psms") is not None else res["psms"]
            if mods:
                res["shifts"] = _read_shifts(mods)
            return res


def _psms_at(idxml, q_max):
    """rank-1 target PSMs with q-value <= q_max (ProSE's main score is the q-value)"""
    import re
    n = 0
    txt = open(idxml).read()
    for m in re.finditer(r"<PeptideIdentification\b(.*?)</PeptideIdentification>", txt, re.S):
        h = re.search(r'<PeptideHit score="([^"]+)"(.*?)</PeptideHit>', m.group(1), re.S)
        if h and 'value="decoy"' not in h.group(2) and float(h.group(1)) <= q_max:
            n += 1
    return n


def _read_shifts(d, top=15):
    """(delta mass, target PSMs) from ProSE's <run>_ModificationAnalysis_DeltaMassStats.tsv: '#' comment lines,
    then columns DeltaMass, Count, UniquePeptides, ..."""
    import csv
    import glob
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "*DeltaMassStats.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                try:
                    m, n = float(r["DeltaMass"]), float(r["Count"])
                except (KeyError, TypeError, ValueError):
                    continue
                rows.append((m, n))
    rows.sort(key=lambda x: -x[1])
    return rows[:top]


# ---------------------------------------------------------------- candidates

def candidates(P, res, history, facts, alt_dbs=()):
    """Shortlist of changes from the diagnostics. Each: (key, description, new params)."""
    out = []

    def add(key, desc, **changes):
        Q = copy.deepcopy(P)
        for k, v in changes.items():
            Q["Search:" + k.replace("__", ":")] = v
        if Q != P and key not in {h["action"] for h in history}:
            out.append((key, desc, Q))

    pe, fe = res.get("precursor_error") or {}, res.get("fragment_error") or {}
    pu, pt = _get(P, "precursor:mass_tolerance_unit", "ppm"), float(_get(P, "precursor:mass_tolerance_upper", 10))
    fu, ft = _get(P, "fragment:mass_tolerance_unit", "ppm"), float(_get(P, "fragment:mass_tolerance", 20))
    if pu == "ppm":
        rec = pe.get("recommended_ppm") if pe.get("valid") else None
        if rec and rec > 0.8 * pt:
            add("widen_precursor", f"calibration PSMs need {rec} ppm, close to the {pt} ppm window",
                precursor__mass_tolerance_lower=round(max(rec * 1.5, pt * 1.5), 1),
                precursor__mass_tolerance_upper=round(max(rec * 1.5, pt * 1.5), 1))
        if rec and rec < 0.4 * pt:
            add("narrow_precursor", f"calibration PSMs fit within {rec} ppm of a {pt} ppm window",
                precursor__mass_tolerance_lower=round(max(rec * 1.5, 3), 1),
                precursor__mass_tolerance_upper=round(max(rec * 1.5, 3), 1))
        if not rec:
            add("widen_precursor", "no calibration result: try a wider precursor window",
                precursor__mass_tolerance_lower=pt * 2, precursor__mass_tolerance_upper=pt * 2)
    if fu == "ppm":
        rec = fe.get("recommended_ppm") if fe.get("valid") else None
        if rec and rec > 0.8 * ft:
            add("widen_fragment", f"fragment errors need {rec} ppm, close to the {ft} ppm window",
                fragment__mass_tolerance=round(max(rec * 1.5, ft * 1.5), 1))
        if rec and rec < 0.4 * ft:
            add("narrow_fragment", f"fragment errors fit within {rec} ppm of a {ft} ppm window",
                fragment__mass_tolerance=round(max(rec * 1.5, 5), 1))
        if not rec and facts.get("ms2_analyzer") != "IT":
            add("widen_fragment", "no fragment calibration: try a wider fragment window", fragment__mass_tolerance=ft * 2)
    if _get(P, "calibration:enabled") != "true":
        add("calibration_on", "let ProSE calibrate the masses first", calibration__enabled="true")
    var = list(_get(P, "modifications:variable", []) or [])
    fixed = list(_get(P, "modifications:fixed", []) or [])
    for m, n in res.get("shifts") or []:
        for mass, mods in SHIFT_MODS:
            if abs(m - mass) < 0.01 and mods and not set(mods) <= set(var + fixed):
                add("add_" + mods[0].split(" ")[0].lower(), f"open search: {int(n)} PSMs carry a {m:+.3f} Da shift",
                    modifications__variable=var + [x for x in mods if x not in var])
    if len(var) > 3:
        for m in var:
            if m != "Oxidation (M)":
                add("drop_" + m.split(" ")[0].lower(), f"{len(var)} variable modifications inflate the search space",
                    modifications__variable=[x for x in var if x != m])
    enz = _get(P, "enzyme")
    if enz == "Trypsin":
        add("trypsin_p", "allow cleavage before proline", enzyme="Trypsin/P")
    if (_get(P, "peptide:missed_cleavages", 2) or 2) < 3:
        add("missed_cleavages_3", "allow one more missed cleavage", peptide__missed_cleavages=3)
    # hard rules: never propose a combination that cannot be right
    keep = []
    for key, desc, Q in out:
        f = dict(facts, fragment_tolerance={"value": _get(Q, "fragment:mass_tolerance"),
                                            "unit": _get(Q, "fragment:mass_tolerance_unit", "ppm")},
                 enzyme=_get(Q, "enzyme"))
        if not S.check(f, kinds=("hard",)):
            keep.append((key, desc, Q))
    return keep


def situation(P, res, history):
    lines = [f"current: {res['score']} PSMs ({res.get('psms')} at the FDR threshold, "
             f"{res.get('ms2')} spectra searched)",
             f"precursor window {_get(P, 'precursor:mass_tolerance_upper')} {_get(P, 'precursor:mass_tolerance_unit')}, "
             f"error: {res.get('precursor_error')}",
             f"fragment window {_get(P, 'fragment:mass_tolerance')} {_get(P, 'fragment:mass_tolerance_unit')}, "
             f"error: {res.get('fragment_error')}",
             f"enzyme {_get(P, 'enzyme')}, fixed {_get(P, 'modifications:fixed')}, "
             f"variable {_get(P, 'modifications:variable')}"]
    if res.get("shifts"):
        lines.append("open-search mass shifts (Da, PSMs): " + ", ".join(f"{m:+.3f}:{int(n)}" for m, n in res["shifts"][:8]))
    for h in history[-8:]:
        lines.append(f"tried {h['action']}: {h['gain']:+.1%} ({'kept' if h['kept'] else 'dropped'})")
    return "\n".join(lines)


def order(cands, P, res, history, jev=None, log=print):
    if not jev or len(cands) < 2:
        return cands
    from .jev import choose
    actions = {k: d for k, d, _ in cands}
    actions["stop"] = "no change is likely to help; keep the current settings"
    try:
        pick, probs = choose(jev, situation(P, res, history), actions,
                             "You tune a proteomics database search. Pick the change most likely to increase "
                             "the identified PSMs at 1% FDR without inflating false identifications.")
    except Exception as e:
        log(f"  Jev chooser failed ({e}); using the rule order")
        return cands
    log(f"  Jev picks {pick} (p={probs.get(pick, 0):.2f})")
    rank = {k: -probs.get(k, 0) for k in actions}
    if pick == "stop" and probs.get("stop", 0) > 0.8:
        return []
    return sorted(cands, key=lambda c: rank.get(c[0], 0))


def tune(mzml, P0, prose: ProSE, facts=None, jev=None, n_spectra=5000, rounds=4, per_round=3,
         min_gain=0.01, open_search=True, log=print, trace=None):
    facts = facts or {}
    with tempfile.TemporaryDirectory(prefix="prose_tune_sub_") as tmp:
        sub = os.path.join(tmp, "subset.mzML")
        n = write_subset_mzml(mzml, sub, n_spectra)
        log(f"  subset: {n} MS2 spectra")
        P = {**{"Search:" + k if not k.startswith("Search:") else k: v for k, v in COMMON.items()}, **P0}
        best = prose(sub, P)
        if open_search:
            try:
                best["shifts"] = prose(sub, P, open_search=True).get("shifts")
            except Exception as e:
                log(f"  open search failed: {e}")
        start = best["score"]
        log(f"  tier 2 start: {start} PSMs")
        history = []
        for rnd in range(rounds):
            cands = order(candidates(P, best, history, facts), P, best, history, jev, log)
            if not cands:
                break
            improved = False
            for key, desc, Q in cands[:per_round]:
                try:
                    r = prose(sub, Q)
                except Exception as e:
                    log(f"  {key}: failed ({e})")
                    history.append({"round": rnd, "action": key, "reason": desc, "gain": 0.0, "kept": False})
                    continue
                gain = (r["score"] - best["score"]) / max(1, best["score"])
                kept = gain >= min_gain
                step = {"round": rnd, "action": key, "reason": desc, "score": r["score"], "gain": gain, "kept": kept,
                        "situation": situation(P, best, history)}
                history.append(step)
                if trace:
                    with open(trace, "a") as fh:
                        fh.write(json.dumps(step, default=str) + "\n")
                log(f"  {key}: {r['score']} PSMs ({gain:+.1%}){' kept' if kept else ''}")
                if kept:
                    r.setdefault("shifts", best.get("shifts"))
                    P, best, improved = Q, r, True
                    break
            if not improved:
                break
        return {"start_score": start, "final_score": best["score"], "gain": (best["score"] - start) / max(1, start),
                "params": P, "history": history, "searches": prose.n, "subset_spectra": n}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="prose-auto-tune", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("report", help="<run>.auto.json from prose-auto analyze")
    ap.add_argument("--mzml", help="spectra (default: the report's input)")
    ap.add_argument("--prose", required=True, help="ProSE executable")
    ap.add_argument("--database", required=True, help="FASTA (with entrapment twins when --entrap is used)")
    ap.add_argument("--entrap", help="entrap.py of the benchmark: score by PSMs at 1%% entrapment FDP")
    ap.add_argument("--jev", action="store_true", help="let Jev order the candidate changes (JEV_API_KEY)")
    ap.add_argument("--spectra", type=int, default=5000)
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--per-round", type=int, default=3)
    ap.add_argument("--min-gain", type=float, default=0.01)
    ap.add_argument("--no-open-search", action="store_true")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--extra", default="", help="extra ProSE arguments, e.g. '-Search:report:top_hits 10'")
    ap.add_argument("-o", "--out", default=None, help="output JSON (default <report>.tune.json)")
    a = ap.parse_args(argv)
    with open(a.report) as fh:
        rep = json.load(fh)
    from .evidence import Evidence  # noqa: F401  (decisions are plain dicts here)
    facts = {k: d.get("value") for k, d in rep["decisions"].items()}
    jev = None
    if a.jev:
        from .jev import Jev
        jev = Jev()
    prose = ProSE(a.prose, a.database, a.threads, a.extra.split(), a.entrap, log=print)
    out = a.out or a.report.replace(".auto.json", "") + ".tune.json"
    res = tune(a.mzml or rep["input"]["given"], rep["prose"]["params"], prose, facts, jev, a.spectra, a.rounds,
               a.per_round, a.min_gain, not a.no_open_search, trace=out.replace(".json", ".trace.jsonl"))
    if jev:
        res["jev"] = {"calls": jev.calls, "cost_usd": round(jev.cost_usd, 5)}
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    print(f"{res['start_score']} -> {res['final_score']} PSMs ({res['gain']:+.1%}) in {res['searches']} searches; {out}")


if __name__ == "__main__":
    main()
