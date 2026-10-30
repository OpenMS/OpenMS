"""Turn all evidence into one decision per parameter.

Policy
------
* Parameters the instrument and the spectra determine (tolerances, offsets,
  acquisition, analyzers, fragmentation, labels): measurements win
  (search > data > metadata), then repository, user, text, llm, default.
* Sample-preparation parameters (organism, protease, Cys alkylation,
  enrichment): a clear mini-search result wins; otherwise the user's own
  description, then the repository record, then methods text, then the LLM,
  then weak diagnostic-ion evidence, then the default.
* Independent sources that agree raise the confidence
  (1 - prod(1 - c_i)); a disagreeing source with confidence >= 0.4 is a
  conflict, and a conflict between a measurement and what the user or the
  repository says puts the parameter on 'review'.
* Tolerances combine: |offset| + 4 sigma from a mini-search when there is
  one; otherwise the data precision plus a measured offset; otherwise the
  instrument default, flagged for review and for ProSE's calibration pass.
"""

from __future__ import annotations

import math

from .detect_instrument import classify_instrument
from .evidence import Decision, Evidence

MEASURED = {"precursor_tolerance", "fragment_tolerance", "precursor_offset", "fragment_offset",
            "acquisition", "dia_windows", "ms1_analyzer", "ms2_analyzer", "fragmentation",
            "collision_energy", "ion_mobility", "charge_range", "instrument", "instrument_class", "label"}
ORDER_MEASURED = ["search", "data", "metadata", "repository", "user", "text", "llm", "default"]
ORDER_SAMPLE = ["user", "repository", "text", "llm", "search", "data", "metadata", "default"]

NEGATIVE = {"no isobaric label", "no SILAC", "none", "unknown", "label free"}


def _same(a, b):
    if isinstance(a, dict) and isinstance(b, dict) and "value" in a and "value" in b:
        if a.get("unit", "").lower() != b.get("unit", "").lower():
            return False
        x, y = abs(float(a["value"])), abs(float(b["value"]))
        return max(x, y) <= 2.0 * max(min(x, y), 1e-9) or abs(x - y) <= 2.0
    na, nb = str(a).lower(), str(b).lower()
    if na == nb:
        return True
    # loose matches: 'TMT10/11plex' vs 'TMT', 'Orbitrap Fusion' vs 'Orbitrap Fusion ETD',
    # 'Escherichia coli K-12' vs 'Escherichia coli', 'DIA (diaPASEF)' vs 'DIA'
    return na.split(" (")[0] in nb or nb.split(" (")[0] in na or \
        (na[:3] == nb[:3] and na[:3] in ("tmt", "itr", "sil", "dia", "dda"))


def _organism_key(v):
    return " ".join(str(v).replace("_", " ").split()[:2]).lower()


def _stated_untested(d, evs):
    """A search can only pick among the databases it searched. When the user or a repository states an organism of
    another genus that no search tested (outside the panel, proteome fetch failed or DIA without a search), the
    nearest panel organism winning says nothing against it: the stated organism wins, flagged for review."""
    if d.source != "search":
        return d
    g = lambda v: str(v).replace("_", " ").split(" ")[0].lower()
    searched = {g(k) for e in evs if e.source == "search" for k in (e.stats or {})}
    stated = [e for e in evs if e.source in ("user", "repository") and e.confidence >= 0.5
              and g(e.value) != g(d.value) and g(e.value) not in searched]
    if not stated:
        return d
    s = max(stated, key=lambda e: e.confidence)
    note = f"the mini-search prefers {d.value}, but {s.source} states {s.value}, which no searched database covers"
    return Decision("organism", s.value, round(s.confidence * 0.9, 3), s.source, evs, d.conflicts + [note], "review")


def _combine(winner: Evidence, evs: list):
    agree = [e for e in evs if e is not winner and _same(e.value, winner.value) and e.source != winner.source]
    p = 1 - winner.confidence
    for e in agree:
        p *= 1 - 0.5 * e.confidence   # agreement counts half: sources are not fully independent
    return min(0.99, 1 - p), agree


def _pick(evs, order, strong_first=None):
    if strong_first:
        strong = [e for e in evs if e.source == strong_first and e.confidence >= 0.7
                  and str(e.value) not in NEGATIVE]
        if strong:
            return max(strong, key=lambda e: e.confidence)
    for src in order:
        c = [e for e in evs if e.source == src]
        if c:
            return max(c, key=lambda e: e.confidence)
    return None


def _decide(param, evs, order=None):
    measured = order is ORDER_MEASURED or param in MEASURED
    order = order or (ORDER_MEASURED if measured else ORDER_SAMPLE)
    w = _pick(evs, order, strong_first=None if measured else "search")
    if w is None:
        return None
    conf, agree = _combine(w, evs)
    conflicts = []
    for e in evs:
        if e is w or e.confidence < 0.4 or _same(e.value, w.value):
            continue
        if param == "organism" and _organism_key(e.value) == _organism_key(w.value):
            continue
        conflicts.append(f"{e.source} ({e.method}) says {e.value}: {e.detail}")
        conf *= 0.85
    status = "auto" if conf >= 0.7 and not conflicts else ("default" if w.source == "default" else "review")
    if conflicts and any(e.source in ("user", "repository") for e in evs if not _same(e.value, w.value)) \
            and w.source in ("data", "search", "metadata"):
        status = "review"
    return Decision(param, w.value, round(conf, 3), w.source, evs, conflicts, status)


def _label(evs):
    pos = [e for e in evs if str(e.value) not in NEGATIVE]
    measured_neg = [e for e in evs if str(e.value) in ("no isobaric label",) and e.source == "data"]
    silac_neg = [e for e in evs if str(e.value) == "no SILAC" and e.source == "data"]
    data_pos = [e for e in pos if e.source == "data" and e.confidence >= 0.6]
    if data_pos:
        d = _decide("label", [e for e in evs if e in data_pos or e.source != "data"])
        d.value = max(data_pos, key=lambda e: e.confidence).value
        d.evidence = evs
        return d
    if measured_neg and silac_neg and not [e for e in pos if e.source in ("user", "repository") and "silac" in str(e.value).lower()]:
        lf = Evidence("label", "label free", min(measured_neg[0].confidence, 0.9), "data", "reporter ions + MS1 pairs",
                      "no reporter ions and no heavy/light MS1 pairs")
        d = _decide("label", [lf] + [e for e in evs if e.source != "data"])
        d.evidence = evs + [lf]
        return d
    return _decide("label", evs)


def _tolerance(kind, evs, offsets, instrument, ms_analyzer):
    """kind: precursor | fragment"""
    param = f"{kind}_tolerance"
    search = [e for e in evs if e.source == "search"]
    if search:
        d = _decide(param, evs)
        return d
    data = [e for e in evs if e.source == "data"]
    info = classify_instrument(instrument or "")
    if kind == "precursor":
        dflt = {"value": info[2] if info else 10.0, "unit": "ppm"}
    else:
        lowres = ms_analyzer in ("IT", "LIT", "PAULIONTRAP")
        dflt = {"value": 0.5, "unit": "Da"} if lowres else {"value": info[3] if info and info[4] == "ppm" else 20.0, "unit": "ppm"}
    default = Evidence(param, dflt, 0.4, "default", "instrument class",
                       f"default for {info[0] if info else 'an unknown instrument'}"
                       + (" with ion-trap MS2" if kind == "fragment" and dflt["unit"] == "Da" else ""))
    if not data:
        d = _decide(param, evs + [default])
        d.evidence = evs + [default]
        if d.source in ("default", "llm"):
            d.status = "review"
            d.conflicts.append("not measured: enable ProSE calibration (Search:calibration:enabled) "
                               "or run the mini-search with a FASTA")
        return d
    best = max(data, key=lambda e: e.confidence)
    v = dict(best.value)
    off = [e for e in offsets if e.source in ("data", "search")]
    if off:
        o = max(off, key=lambda e: e.confidence)
        if o.value.get("unit") == v.get("unit"):
            v["value"] = round(v["value"] + abs(o.value["value"]), 3)
        note = f"precision from {best.method} plus the {o.method} offset {o.value['value']:+} {o.value['unit']}"
        conf = min(best.confidence, o.confidence) + 0.1
    else:
        # precision only: never go below the instrument default, the offset is unknown
        if v.get("unit") == dflt["unit"]:
            v["value"] = max(v["value"], dflt["value"])
        note = f"precision from {best.method}; offset unknown, kept at least at the instrument default"
        conf = best.confidence * 0.8
    comb = Evidence(param, v, conf, "data", "combined", note)
    # the combined value is the decision; the raw precision estimates are its inputs
    conflicts = [f"{e.source} ({e.method}) says {e.value}: {e.detail}" for e in evs
                 if e.source in ("repository", "user", "text") and e.confidence >= 0.4 and not _same(e.value, v)]
    status = "auto" if off and conf >= 0.7 and not conflicts else "review"
    if not off:
        conflicts.append("systematic offset not measured: ProSE's calibration pass is recommended")
    return Decision(param, v, round(min(conf, 0.99), 3), "data", evs + [default, comb], conflicts, status)


def reconcile(evidence: list) -> dict:
    by = {}
    for e in evidence:
        by.setdefault(e.parameter, []).append(e)
    out = {}
    get = lambda p: out.get(p).value if out.get(p) else None
    for p in ("instrument", "instrument_class", "ms1_analyzer", "ms2_analyzer", "fragmentation",
              "collision_energy", "ion_mobility", "acquisition", "dia_windows", "charge_range",
              "precursor_offset", "fragment_offset", "enzyme", "cys_alkylation", "enrichment",
              "variable_mods"):
        if p in by:
            # enrichment: when the spectra could test it (>= 100 charged DDA spectra, confidence >= 0.6) the
            # measurement wins; project texts mention phospho for atlases whose other runs were enriched
            tested = p == "enrichment" and any(e.source == "data" and e.confidence >= 0.6 for e in by[p])
            d = _decide(p, by[p], ORDER_MEASURED if tested else None)
            if d:
                out[p] = d
    if "label" in by:
        d = _label(by["label"])
        if d:
            out["label"] = d
    # organism: the search result names a panel file; match it to the stated species
    if "organism" in by:
        # "not in panel" only says the panel search could not decide; it must not outrank a real organism call
        # (FASTag species, records, text): on DIA runs and TMT runs with reporters in MS3 the mini-search finds too
        # few PSMs while FASTag names the species
        real = [e for e in by["organism"] if e.value not in ("not in panel", "unknown") and e.confidence >= 0.5]
        if real:
            by["organism"] = [e for e in by["organism"] if e.value != "not in panel"]
        d = _decide("organism", by["organism"])
        if d:
            out["organism"] = _stated_untested(d, by["organism"])
    instrument = get("instrument") or ""
    for kind in ("precursor", "fragment"):
        param = f"{kind}_tolerance"
        d = _tolerance(kind, by.get(param, []), by.get(f"{kind}_offset", []), instrument,
                       get("ms2_analyzer") if kind == "fragment" else get("ms1_analyzer"))
        if d:
            out[param] = d
    # defaults for what nothing determined
    for p, v, why in (("enzyme", "Trypsin", "most common protease"),
                      ("cys_alkylation", "Carbamidomethyl (C)", "most common alkylation"),
                      ("label", "label free", "no label evidence"),
                      ("enrichment", "none", "no enrichment evidence")):
        if p not in out or (out[p].value == "unknown" and p in ("cys_alkylation", "enzyme")):
            e = Evidence(p, v, 0.3, "default", "default", why)
            prev = out[p].evidence if p in out else []
            out[p] = Decision(p, v, 0.3, "default", prev + [e], [], "default")
    if "organism" not in out:
        out["organism"] = Decision("organism", "unknown", 0.0, "default", [], ["give a FASTA, a "
                                   "description, a PXD accession or an organism panel"], "unknown")
    return out
