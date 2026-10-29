"""Benchmark: how well do the study texts alone predict the search settings?

For PRIDE datasets that have a curated SDRF (bigbio/sdrf-annotated-datasets),
the texts a reader would see before downloading anything (PRIDE title,
project description, sample and data processing protocols, optionally the
paper's methods from Europe PMC) are read by

* jev     one multi-head Jev request (schema.HEADS + flags, joint decoding)
* llm     one multi-head LLM request (llm.extract_study), optional
* regex   the keyword extractor (text_extract.py)
* pride   PRIDE's own structured fields (instrument, organism, quantification),
          as the "just ask the repository" baseline

and scored against the SDRF, which none of the readers sees. Truth is
dataset level: a prediction counts as correct when it matches any value the
SDRF lists for that column (multi-instrument or multi-species studies list
several). Heads the SDRF does not fill are not scored.

Outputs: <out>/rows.tsv (dataset x head x predictor), <out>/summary.md
(accuracy, coverage, accuracy when answered, Jev calibration by probability),
<out>/raw/<PXD>.json (texts used, answers).

    python -m prose_auto.bench_text --sdrf-repo sdrf-annotated-datasets -n 200 -o bench_text/ [--llm codex] [--paper]
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import random
import re
import time
from collections import defaultdict

from . import repository as repo
from . import schema as S
from . import text_extract as tx

HEADS = ["organism", "instrument", "instrument_family", "label", "label_class", "enzyme", "cys_alkylation",
         "acquisition", "fragmentation", "precursor_tolerance", "fragment_tolerance"]
MOD_FLAGS = {"Oxidation (M)": "oxidation", "Acetyl (Protein N-term)": "acetyl", "Deamidated (N)": "deamidated",
             "Gln->pyro-Glu (N-term Q)": "pyro", "Phospho (S)": "phospho", "Carbamyl (N-term)": "carbamyl",
             "Methyl (K)": "methyl"}

# ---------------------------------------------------------------- normalisation

FAMILIES = [("astral", "Orbitrap Astral"), ("exploris", "Orbitrap Exploris"), ("ascend", "Orbitrap Tribrid"),
            ("eclipse", "Orbitrap Tribrid"), ("lumos", "Orbitrap Tribrid"), ("fusion", "Orbitrap Tribrid"),
            ("exactive", "Q Exactive"), ("elite", "LTQ Orbitrap"), ("velos", "LTQ Orbitrap"),
            ("ltq orbitrap", "LTQ Orbitrap"), ("orbitrap", "LTQ Orbitrap"), ("timstof", "timsTOF"),
            ("tripletof", "SCIEX TOF"), ("zenotof", "SCIEX TOF"), ("6600", "SCIEX TOF"), ("5600", "SCIEX TOF"),
            ("synapt", "Waters"), ("xevo", "Waters"), ("ltq", "ion trap"), ("velos pro", "ion trap"),
            ("maxis", "Bruker Q-TOF"), ("impact", "Bruker Q-TOF")]


def _nt(cell):
    kv = repo._kv(cell or "")
    return (kv.get("NT") or cell or "").strip()


def _clean(s):
    s = str(s or "").lower()
    s = re.sub(r"thermo( fisher)?( scientific)?|bruker( daltonics)?|sciex|ab sciex|waters|instrument model", " ", s)
    return re.sub(r"[^a-z0-9+]+", " ", s).strip()


def canon_instrument(s):
    c = _clean(s)
    if not c:
        return None
    key = lambda x: _clean(x).replace("ltq ", "").replace(" ", "")
    for name in sorted(tx.INSTRUMENTS, key=len, reverse=True):
        if key(name) and key(name) in key(c):
            return name
    return c


def instrument_family(s):
    c = _clean(s)
    if not c:
        return None
    for key, fam in FAMILIES:
        if key in c:
            return fam
    return c


def label_from_sdrf(values):
    """comment[label] values of a dataset -> (exact prose-auto label, class)."""
    vals = [_nt(v).lower() for v in values if v]
    if not vals:
        return None, None
    if all("label free" in v or v in ("lfq", "label-free") for v in vals):
        return "label free", "label free"
    tmt = {v for v in vals if "tmt" in v}
    itraq = {v for v in vals if "itraq" in v}
    if tmt:
        n = len(tmt)
        if any("pro" in v or re.search(r"13[2-5]", v) for v in tmt) or n > 11:
            return "TMTpro (16/18plex)", "TMT"
        return ("TMT2plex" if n <= 2 else "TMT6plex" if n <= 6 else "TMT10/11plex"), "TMT"
    if itraq:
        return ("iTRAQ4plex" if len(itraq) <= 4 else "iTRAQ8plex"), "iTRAQ"
    if any("silac" in v or "heavy" in v or "medium" in v for v in vals):
        return None, "SILAC"
    if any("dimethyl" in v for v in vals):
        return None, "dimethyl"
    return None, vals[0]


def label_class(v):
    v = str(v or "")
    return ("TMT" if v.startswith("TMT") else "iTRAQ" if v.startswith("iTRAQ") else
            "SILAC/dimethyl" if v.startswith("SILAC") and "dimethyl" in v else   # Lys4 and dimethyl-d4: same mass
            "SILAC" if v.startswith("SILAC") else "dimethyl" if "dimethyl" in v else
            "label free" if v == "label free" else v or None)


def enzyme_norm(v):
    c = _clean(_nt(v))
    if not c:
        return None
    if "chymotrypsin" in c:         # contains "trypsin": test it first
        return "Chymotrypsin"
    if "trypsin" in c:
        return "Trypsin"            # Trypsin/P and Lys-C + trypsin count as trypsin
    if "glutamyl" in c:              # glutamyl endopeptidase = Glu-C
        return "Glu-C"
    for k in ("lys c", "lysc"):
        if k in c:
            return "Lys-C"
    return {"chymotrypsin": "Chymotrypsin", "glu c": "Glu-C", "asp n": "Asp-N", "arg c": "Arg-C",
            "no cleavage": "unspecific cleavage", "unspecific": "unspecific cleavage",
            "non specific": "unspecific cleavage"}.get(c, c)


def frag_norm(v):
    c = _clean(_nt(v))
    if not c:
        return None
    if "beam" in c or "hcd" in c or "higher energy" in c:
        return "HCD"
    if "ethcd" in c:
        return "EThcD"
    if "electron transfer" in c or c == "etd":
        return "ETD"
    if "collision" in c or c == "cid":
        return "CID"
    return c


def acq_norm(v):
    c = _clean(_nt(v))
    if not c:
        return None
    if "independent" in c or "dia" in c.split() or "swath" in c or c.startswith("diapasef") or "dia pasef" in c:
        return "DIA"
    if "dependent" in c or "dda" in c.split() or "pasef" in c:
        return "DDA"
    if "parallel" in c or "prm" in c.split() or "srm" in c or "selected reaction" in c:
        return "PRM"
    return c


def tol_bucket(t, kind):
    if not isinstance(t, dict):
        return None
    ppm = t["value"] * (1.0 if t["unit"] == "ppm" else 1000.0)
    if kind == "precursor":
        return "<=5 ppm" if ppm <= 5.5 else "10 ppm" if ppm <= 15 else "20 ppm" if ppm <= 30 else \
            "50 ppm" if ppm <= 100 else ">=0.5 Da"
    if t["unit"] != "ppm":  # Da: the option cut-offs are in Da (0.03 Da is "0.05 Da", not "20 ppm")
        da = float(t["value"])
        return "<=10 ppm" if da <= 0.0105 else "20 ppm" if da <= 0.025 else "0.05 Da" if da <= 0.1 else \
            "0.5 Da" if da <= 0.65 else ">=0.8 Da"
    return "<=10 ppm" if ppm <= 10.5 else "20 ppm" if ppm <= 30 else "0.05 Da" if ppm <= 100 else \
        "0.5 Da" if ppm <= 650 else ">=0.8 Da"


def org_norm(v):
    c = str(_nt(v) or "").strip().lower()
    c = re.sub(r"\s*\(.*\)$", "", c)
    return " ".join(c.split()[:2]) if c and c not in ("not available", "not applicable") else None


# ---------------------------------------------------------------- truth

def sdrf_truth(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if len(rows) < 2:
        return {}
    hdr = [h.strip().lower() for h in rows[0]]

    def col(*names):
        vals = []
        for i, h in enumerate(hdr):
            if any(h == n or h.startswith(n) for n in names):
                vals += [r[i] for r in rows[1:] if i < len(r) and r[i].strip()
                         and _nt(r[i]).lower() not in ("not available", "not applicable", "na", "n/a")]
        return vals

    T = {}
    orgs = {org_norm(v) for v in col("characteristics[organism]")} - {None}
    if orgs:
        T["organism"] = orgs
    ins = col("comment[instrument]")
    if ins:
        T["instrument"] = {canon_instrument(_nt(v)) for v in ins} - {None}
        T["instrument_family"] = {instrument_family(_nt(v)) for v in ins} - {None}
    lab, cls = label_from_sdrf(col("comment[label]"))
    if lab:
        T["label"] = {lab}
    if cls:
        T["label_class"] = {cls}
    enz = {enzyme_norm(v) for v in col("comment[cleavage agent details]")} - {None}
    if enz:
        T["enzyme"] = enz
    fixed_c, var = set(), set()
    for v in col("comment[modification parameters]"):
        kv = repo._kv(v)
        name, mt, ta = kv.get("NT", ""), kv.get("MT", "").lower(), kv.get("TA", "")
        if mt == "fixed" and "C" in ta:
            fixed_c.add(name)
        elif mt == "variable":
            var.add(name.lower())
    alk_map = {"carbamidomethyl": "Carbamidomethyl (C)", "methylthio": "Methylthio (C)",
               "propionamide": "Propionamide (C)", "nethylmaleimide": "Nethylmaleimide (C)",
               "carboxymethyl": "Carboxymethyl (C)"}
    alk = {alk_map.get(re.sub(r"[^a-z]", "", n.lower()), n) for n in fixed_c}
    if alk:
        T["cys_alkylation"] = alk
    if var:
        T["variable_mods"] = var
    acq = {acq_norm(v) for v in col("comment[proteomics data acquisition method]")} - {None}
    if acq:
        T["acquisition"] = acq
    fr = {frag_norm(v) for v in col("comment[dissociation method]")} - {None}
    if fr:
        T["fragmentation"] = fr
    for kind in ("precursor", "fragment"):
        b = {tol_bucket(repo._tol(v), kind) for v in col(f"comment[{kind} mass tolerance]")} - {None}
        if b:
            T[f"{kind}_tolerance"] = b
    return T


# ---------------------------------------------------------------- predictions

def _norm_pred(head, v):
    if v in (None, "", S.NOT_STATED, "unknown"):
        return None
    if head == "organism":
        return "other" if v == "other species" else "mixed species" if v == "mixed species" else org_norm(v)
    if head == "instrument":
        return "other" if v == "other instrument" else canon_instrument(v)
    if head == "instrument_family":
        return instrument_family(v)
    if head == "enzyme":
        return enzyme_norm(v)
    if head == "acquisition":
        return acq_norm({"DDA SPS-MS3": "DDA", "ddaPASEF": "DDA", "DIA (diaPASEF)": "DIA"}.get(v, v))
    if head == "fragmentation":
        return frag_norm(v)
    if head == "label_class":
        return label_class(v)
    if head.endswith("_tolerance"):
        return v if isinstance(v, str) else tol_bucket(v, head.split("_")[0])
    return v


def predictions_from_heads(heads: dict):
    P = {}
    for h in HEADS:
        src = {"instrument_family": "instrument", "label_class": "label"}.get(h, h)
        P[h] = _norm_pred(h, heads.get(src))
    return P


def regex_predict(text):
    ev = tx.extract(text, source="text")
    heads = {}
    for e in ev:
        if e.parameter in ("organism", "instrument", "label", "enzyme", "cys_alkylation", "acquisition",
                           "fragmentation") and e.parameter not in heads:
            v = e.value
            if e.parameter == "label" and v == "SILAC":
                v = "SILAC (Lys8/Arg10)"
            heads[e.parameter] = v
        if e.parameter.endswith("_tolerance") and e.parameter not in heads:
            heads[e.parameter] = tol_bucket(e.value, e.parameter.split("_")[0])
    return predictions_from_heads(heads)


def correct(head, pred, truth):
    if pred is None or truth is None:
        return None
    if pred == "other":
        known = {org_norm(o) for o in list(tx.ORGANISMS) + S.PRIDE_SPECIES} if head == "organism" else {canon_instrument(i) for i in tx.INSTRUMENTS}
        return int(not (truth & known))
    if head == "organism" and pred == "mixed species":    # several species in the (dataset-level) truth
        return int(len(truth) >= 2)
    if head == "label_class" and pred == "SILAC/dimethyl":
        return int(bool(truth & {"SILAC", "dimethyl"}))
    if head == "instrument" and pred not in truth:
        # an SDRF that names no model ("ltq orbitrap", "timstof series", "synapt ms") is right for any model
        # of that family
        models = {canon_instrument(i) for i in tx.INSTRUMENTS}
        generic = [x for x in truth if x not in models]
        return int(any(instrument_family(x) == instrument_family(pred) for x in generic))
    return int(pred in truth)


# ---------------------------------------------------------------- run

def texts_for(acc, paper=False, log=print):
    _, t = repo.pride_evidence(acc, log=lambda *a: None)
    if not t:
        return None, {}
    texts = {"PRIDE title": t.get("title", ""), "PRIDE description": t.get("description", ""),
             "PRIDE sample protocol": t.get("sample_protocol", ""),
             "PRIDE data protocol": t.get("data_protocol", "")}
    if paper:
        for ref in t.get("references") or []:
            key = (ref.get("pubmedId") or ref.get("doi")) if isinstance(ref, dict) else None
            if key:
                m = repo.methods_from_europepmc(key, log=lambda *a: None)
                if m:
                    texts["paper methods"] = m[:40000]
                    break
    import json as _j
    raw = repo.fetch(repo.PRIDE.format(acc=acc), log=lambda *a: None)
    pride = _j.loads(raw) if raw else {}
    return {k: v for k, v in texts.items() if v}, pride


def pride_predict(pride):
    heads = {}
    if pride.get("organisms"):
        heads["organism"] = pride["organisms"][0].get("name")
    if pride.get("instruments"):
        heads["instrument"] = pride["instruments"][0].get("name")
    q = " ".join(x.get("name", "") for x in pride.get("quantificationMethods", [])).lower()
    if q:
        heads["label"] = ("label free" if "label free" in q or "label-free" in q else
                          "TMT6plex" if "tmt" in q else "iTRAQ4plex" if "itraq" in q else
                          "SILAC (Lys8/Arg10)" if "silac" in q else None)
    P = predictions_from_heads(heads)
    for h in ("label",):
        if P.get(h) and P[h] != "label free":
            P[h] = None   # PRIDE gives the method, not the plex
    return P


def select(sdrf_repo, n, seed, only=None):
    files = sorted(glob.glob(os.path.join(sdrf_repo, "datasets", "PXD*", "*.sdrf.tsv")))
    by_acc = {}
    for f in files:
        acc = os.path.basename(os.path.dirname(f))
        if only is None or acc in only:
            by_acc.setdefault(acc, f)
    accs = sorted(by_acc)
    random.Random(seed).shuffle(accs)
    return [(a, by_acc[a]) for a in accs[:n]]


def run(a):
    from . import jev as J
    from . import llm as L
    os.makedirs(os.path.join(a.out, "raw"), exist_ok=True)
    jev = J.Jev() if not a.no_jev else None
    backend = L.make_backend(a.llm) if a.llm else None
    rows, done = [], 0
    accs = getattr(a, "accs", None)
    only = {l.split()[0] for l in open(accs) if l.startswith("PXD")} if accs else None
    for acc, sdrf in select(a.sdrf_repo, a.n * 3, a.seed, only):
        if done >= a.n:
            break
        T = sdrf_truth(sdrf)
        if not T:
            continue
        texts, pride = texts_for(acc, a.paper)
        if not texts or sum(len(v) for v in texts.values()) < 200:
            continue
        rec = {"acc": acc, "texts": {k: len(v) for k, v in texts.items()}, "truth": {k: sorted(map(str, v)) for k, v in T.items()}}
        preds = {"regex": regex_predict("\n\n".join(texts.values())), "pride": pride_predict(pride)}
        probs = {}
        if jev:
            t = time.time()
            ev, rep = J.read_study(jev, texts)
            if rep.get("error"):
                rec["jev_error"] = rep["error"]
            else:
                preds["jev"] = predictions_from_heads(rep["heads"])
                probs = {h: rep["probabilities"].get({"instrument_family": "instrument", "label_class": "label"}.get(h, h), {})
                         .get(rep["heads"].get({"instrument_family": "instrument", "label_class": "label"}.get(h, h)), None)
                         for h in HEADS}
                rec["jev"] = {"heads": rep["heads"], "flags": rep["flags"], "changed": rep["changed_by_rules"],
                              "text_rules": rep["text_rules"], "seconds": round(time.time() - t, 2)}
                for mod, key in MOD_FLAGS.items():
                    if "variable_mods" in T:
                        truth = int(any(key in m for m in T["variable_mods"]))
                        rows.append({"acc": acc, "head": "var:" + key, "predictor": "jev",
                                     "pred": round(rep["flags"].get(mod, 0), 3), "truth": truth, "correct": "",
                                     "p": rep["flags"].get(mod, 0)})
        if backend:
            ev, rej, conf = L.extract_study(backend, texts)
            heads = {}
            for e in ev:
                heads.setdefault({"variable_mods": "_vm"}.get(e.parameter, e.parameter), e.value)
            preds["llm"] = predictions_from_heads(heads)
            rec["llm"] = {"rejected": len(rej), "conflicts": conf}
        for pred_name, P in preds.items():
            for h in HEADS:
                truth = T.get(h)
                if truth is None:
                    continue
                c = correct(h, P.get(h), truth)
                rows.append({"acc": acc, "head": h, "predictor": pred_name, "pred": P.get(h) or "",
                             "truth": "|".join(sorted(map(str, truth))), "correct": "" if c is None else c,
                             "p": probs.get(h) if pred_name == "jev" else ""})
        rec["pred"] = preds
        with open(os.path.join(a.out, "raw", f"{acc}.json"), "w", encoding="utf-8") as fh:
            json.dump(rec, fh, indent=1, default=str)
        done += 1
        print(f"{done}/{a.n} {acc}", flush=True)
    with open(os.path.join(a.out, "rows.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["acc", "head", "predictor", "pred", "truth", "correct", "p"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    summary = summarize(rows, done, jev)
    with open(os.path.join(a.out, "summary.md"), "w", encoding="utf-8") as fh:
        fh.write(summary)
    print(summary)


def summarize(rows, n, jev=None):
    acc = defaultdict(lambda: [0, 0, 0])      # scored, answered, correct
    for r in rows:
        if r["head"].startswith("var:"):
            continue
        k = (r["head"], r["predictor"])
        acc[k][0] += 1
        if r["correct"] != "":
            acc[k][1] += 1
            acc[k][2] += int(r["correct"])
    preds = sorted({r["predictor"] for r in rows})
    L_ = [f"# Settings from study texts vs curated SDRF ({n} PRIDE datasets)", "",
          "Readers see the PRIDE title, description and protocols (and the paper's methods with --paper); "
          "none sees the SDRF. Cells: correct / scored (correct / answered).", "",
          "| Head | " + " | ".join(preds) + " |", "|---|" + "---|" * len(preds)]
    for h in HEADS:
        cells = []
        for p in preds:
            s, an, c = acc[(h, p)]
            cells.append(f"{c}/{s} ({c}/{an})" if s else "–")
        L_.append(f"| {h} | " + " | ".join(cells) + " |")
    # Jev calibration: accuracy by probability bin
    bins = defaultdict(lambda: [0, 0])
    for r in rows:
        if r["predictor"] == "jev" and r["correct"] != "" and r["p"] not in ("", None) and not r["head"].startswith("var:"):
            b = min(4, int(float(r["p"]) * 5))
            bins[b][0] += 1
            bins[b][1] += int(r["correct"])
    if bins:
        L_ += ["", "Jev calibration (answered heads): probability bin, accuracy", ""]
        for b in sorted(bins):
            t, c = bins[b]
            L_.append(f"- p {b / 5:.1f}–{(b + 1) / 5:.1f}: {c}/{t} ({c / t:.0%})")
    # variable-mod flags
    vm = defaultdict(lambda: [0, 0, 0, 0])   # tp fp fn tn at 0.5
    for r in rows:
        if r["head"].startswith("var:"):
            y, p = int(r["truth"]), float(r["pred"]) >= 0.5
            vm[r["head"]][(0 if p else 2) + (0 if y else 1)] += 1
    if vm:
        L_ += ["", "Jev variable-modification flags vs SDRF variable mods (p ≥ 0.5): TP / FP / FN / TN", ""]
        for k, (tp, fp, fn, tn) in sorted(vm.items()):
            L_.append(f"- {k[4:]}: {tp} / {fp} / {fn} / {tn}")
    if jev:
        L_ += ["", f"Jev: {jev.calls} requests, about ${jev.cost_usd:.3f}."]
    return "\n".join(L_) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(prog="prose-auto-bench-text", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sdrf-repo", required=True, help="local clone of bigbio/sdrf-annotated-datasets")
    ap.add_argument("-n", type=int, default=100, help="datasets to score")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--paper", action="store_true", help="add the paper's methods from Europe PMC")
    ap.add_argument("--llm", help="also run the generative LLM reader (e.g. codex, claude)")
    ap.add_argument("--no-jev", action="store_true")
    ap.add_argument("--accs", help="file with the accessions to score (one per line); default: a random draw")
    ap.add_argument("-o", "--out", default="bench_text")
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()
