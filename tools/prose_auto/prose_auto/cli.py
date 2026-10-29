"""prose-auto: propose ProSE search settings from the data itself.

    prose-auto panel proteomes/        (download the model-organism panel once)
    prose-auto analyze run.mzML [--panel proteomes/ | --fasta db.fasta]
                                [--description "HeLa, TMTpro, IAA" | --description @notes.txt]
                                [--methods methods.txt] [--accession PXD012345]
                                [--llm none|claude|codex|cli:<cmd>|openai:<model>] [--jev]
                                [--prose /path/to/ProSE] [--open-search] [-o outdir]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

from . import detect_chemistry as chem
from . import detect_instrument as inst
from . import detect_tolerance as tol
from . import external as ext
from . import jev as J
from . import llm as L
from . import prose_params as pp
from . import repository as repo
from . import tagsearch as ts
from . import text_extract as tx
from .hypotheses import HypothesisTester
from .evidence import Evidence
from .reconcile import reconcile
from .spectra import load_mzml

LOW_RES = {"IT", "LIT", "PAULIONTRAP", "RADIALEJECTIONLINEARIONTRAP", "AXIALEJECTIONLINEARIONTRAP"}


def _arg_text(v):
    if not v:
        return ""
    if v.startswith("@") and os.path.exists(v[1:]):
        return open(v[1:], encoding="utf-8", errors="replace").read()
    if os.path.exists(v):
        return open(v, encoding="utf-8", errors="replace").read()
    return v


def _first(evs, param):
    c = [e for e in evs if e.parameter == param]
    return max(c, key=lambda e: e.confidence) if c else None


def analyze(a, log=print):
    t0 = time.time()
    timings = {}
    ev = []
    os.makedirs(a.out, exist_ok=True)
    base = os.path.join(a.out, os.path.splitext(os.path.basename(a.mzml.rstrip("/")))[0])
    spectra_file = a.mzml          # the MS2 spectra the detectors and searches see
    if a.mzml.rstrip("/").lower().endswith(".d"):
        # Bruker diaPASEF: pseudo-MS/MS spectra first (pyopenms cannot read .d)
        if not a.diaspextract:
            raise SystemExit("a Bruker .d input needs --diaspextract")
        log(f"[0/7] DIAspeXtract: {a.mzml} -> pseudo-MS/MS spectra")
        try:
            spectra_file = ext.diaspextract(a.diaspextract, a.mzml, base + ".pseudo.mzML", a.threads, log=log)
        except RuntimeError as e:
            raise SystemExit(f"DIAspeXtract cannot read {a.mzml} (only diaPASEF .d is supported; "
                             f"DDA-PASEF needs a conversion to mzML first): {str(e).strip()[-300:]}")
        ev.append(Evidence("acquisition", "DIA (diaPASEF)", 0.9, "metadata", "Bruker .d input",
                           "a Bruker .d run converted to pseudo-spectra by DIAspeXtract"))
        ev.append(Evidence("ion_mobility", True, 0.95, "metadata", "Bruker .d input", "timsTOF data"))
    log(f"[1/7] reading {spectra_file}")
    run = load_mzml(spectra_file, max_ms2=a.max_ms2)
    timings["read"] = round(time.time() - t0, 1)

    log("[2/7] instrument, fragmentation, acquisition mode")
    ev += inst.detect_instrument(run)
    if spectra_file == a.mzml:
        ev += inst.detect_acquisition(run)
    acq = _first(ev, "acquisition").value
    im = _first(ev, "ion_mobility")
    if str(acq).startswith("DIA") and a.diaspextract and spectra_file == a.mzml and im and im.value is True:
        # DIAspeXtract needs ion mobility (diaPASEF); Orbitrap DIA stays on the DIA path
        log("  DIA mzML with ion mobility: DIAspeXtract pseudo-MS/MS spectra for the DDA-style stages")
        try:
            spectra_file = ext.diaspextract(a.diaspextract, a.mzml, base + ".pseudo.mzML", a.threads, log=log)
            run = load_mzml(spectra_file, max_ms2=a.max_ms2)
        except RuntimeError as e:
            log(f"  DIAspeXtract failed, continuing without pseudo-spectra: {str(e).strip()[-300:]}")
    elif str(acq).startswith("DIA") and a.diaspextract and spectra_file == a.mzml:
        log("  DIA mzML without ion mobility: DIAspeXtract skipped (it needs diaPASEF)")
    pseudo = spectra_file != a.mzml
    search_acq = "DDA" if pseudo else acq      # pseudo-spectra have one precursor each
    ms1a = (_first(ev, "ms1_analyzer").value if _first(ev, "ms1_analyzer") else "ORBITRAP")
    ms2a = (_first(ev, "ms2_analyzer").value if _first(ev, "ms2_analyzer") else "ORBITRAP")
    hr1, hr2 = ms1a not in LOW_RES, ms2a not in LOW_RES

    log("[3/7] mass accuracy from repeated ions and known ions")
    t = time.time()
    ev += tol.detect_tolerances(run, search_acq, ms2_high_res=hr2, ms1_high_res=hr1)
    timings["tolerance"] = round(time.time() - t, 1)

    log("[4/7] labels, enrichment, Cys alkylation (diagnostic ions)")
    t = time.time()
    ev += chem.detect_labels(run, hr2) + chem.detect_enrichment(run, hr2) + chem.detect_cys_alkylation(run, hr2)
    timings["chemistry"] = round(time.time() - t, 1)

    # mini-search hypotheses (DDA only: DIA spectra have no single precursor)
    panel = ts.load_panel(a.panel) if a.panel else {}
    if a.fasta:
        panel = {**panel, os.path.splitext(os.path.basename(a.fasta))[0]: a.fasta}
    h = None
    if panel and not str(search_acq).startswith("DIA") and not a.no_search:
        log(f"[5/7] mini-search on {min(a.spectra, len(run.ms2))} spectra against {len(panel)} database(s)")
        t = time.time()
        lab = _first([e for e in ev if e.parameter == "label" and e.source == "data"
                      and e.value not in ("no isobaric label", "no SILAC", "unknown")], "label")
        mods = ts.Mods()
        if lab and str(lab.value).startswith("TMTpro"):
            mods.k = mods.nterm = ts.LABEL_MASS["TMTpro"]
        elif lab and str(lab.value).startswith("TMT"):
            mods.k = mods.nterm = ts.LABEL_MASS["TMT"]
        elif lab and str(lab.value).startswith("iTRAQ8"):
            mods.k = mods.nterm = ts.LABEL_MASS["iTRAQ8plex"]
        elif lab and str(lab.value).startswith("iTRAQ"):
            mods.k = mods.nterm = ts.LABEL_MASS["iTRAQ4plex"]
        frag_tol = 20.0 if hr2 else 0.4
        h = HypothesisTester(run, frag_tol, hr2, prec_ppm=20.0 if hr1 else 500.0,
                             n_spectra=a.spectra, mods=mods, log=log)
        ev += h.organism(panel)
        if h.best_db:
            ev += h.accuracy() + h.alkylation()
            if not a.skip_enzyme:
                ev += h.enzyme()
        timings["mini_search"] = round(time.time() - t, 1)
    else:
        log("[5/7] mini-search skipped" + (" (DIA without pseudo-spectra)" if str(search_acq).startswith("DIA") else
                                           " (no --panel/--fasta)" if not panel else ""))
    if a.fastag:
        t = time.time()
        ftol, funit = (20.0, "ppm") if hr2 else (0.3, "Da")
        log("  FASTag: genus call from sequence tags")
        try:
            ev += ext.fastag_species(a.fastag, spectra_file, ftol, funit, threads=a.threads,
                                     subsample=a.fastag_subsample, log=log)
            db = (h.best_db[1] if h is not None and h.best_db else None) or a.fasta
            if db:
                log(f"  FASTag: mass-gap histogram against {os.path.basename(db)} without fixed modifications")
                ev += ext.fastag_alkylation(a.fastag, spectra_file, db, ftol, funit, threads=a.threads, log=log)
        except Exception as e:
            log(f"  FASTag failed: {e}")
        timings["fastag"] = round(time.time() - t, 1)

    log("[6/7] repository records and text")
    texts = {}
    if a.accession:
        e, texts = repo.repository_evidence(a.accession, log=log, fetch_paper=not a.no_paper)
        ev += e
    desc = _arg_text(a.description)
    methods = _arg_text(a.methods)
    if desc:
        ev += tx.extract(desc, source="user", base_conf=0.75)
    if methods:
        ev += tx.extract(methods, source="text", base_conf=0.55)
    for k in ("sample_protocol", "data_protocol", "description", "methods"):
        if texts.get(k):
            ev += tx.extract(texts[k], source="text", base_conf=0.5)
    study = {"user's description": desc, "methods section": methods,
             "PRIDE title": texts.get("title", ""), "PRIDE description": texts.get("description", ""),
             "PRIDE sample protocol": texts.get("sample_protocol", ""),
             "PRIDE data protocol": texts.get("data_protocol", ""),
             "paper methods": texts.get("methods", ""), "SDRF": texts.get("sdrf", "")[:8000]}
    study = {k: v for k, v in study.items() if v}
    backend = L.make_backend(a.llm)
    llm_rejected, llm_text_conflicts = [], []
    if backend and study:
        log(f"  LLM: {backend.name} (one multi-head request over {len(study)} texts)")
        e, rej, llm_text_conflicts = L.extract_study(backend, study, user=bool(desc))
        ev += e
        llm_rejected += rej
        if a.llm_priors:
            ins = _first(ev, "instrument")
            frg = _first(ev, "fragmentation")
            ev += L.literature_priors(backend, ins.value if ins else "", acq, frg.value if frg else "", ms2a)

    jev = J.Jev() if a.jev else None
    jev_report, jev_issues = {}, []
    if jev and study:
        log(f"  Jev: {jev.name} (one multi-head request over {len(study)} texts)")
        e, jev_report = J.read_study(jev, study, user=bool(desc))
        ev += e

    log("[7/7] reconciling")
    decisions = reconcile(ev)
    prose = pp.to_prose(decisions)
    if a.prose:
        from .prose_runner import confirm_with_prose, find_prose
        exe = find_prose(a.prose)
        db = a.fasta or (panel.get(decisions["organism"].value) if panel else None)
        if exe and db:
            log(f"  confirming with ProSE ({exe})")
            try:
                ev += confirm_with_prose(exe, a.mzml, db, pp.to_cli(prose["params"]),
                                         open_search=a.open_search, log=log)
                decisions = reconcile(ev)
                prose = pp.to_prose(decisions)
            except Exception as e:  # keep the data-only result
                log(f"  ProSE confirmation failed: {e}")
    review = L.review(backend, decisions) if backend and a.llm_review else []
    if jev and study:
        jev_issues = J.cross_check(jev, study, decisions)
    consistency = J.consistency(decisions)
    timings["total"] = round(time.time() - t0, 1)
    return {
        "tool": "prose-auto 0.1 (prototype)",
        "input": {**run.summary(), "given": a.mzml, "pseudo_spectra": spectra_file if pseudo else None},
        "decisions": {k: d.to_dict() for k, d in decisions.items()},
        "prose": prose,
        "prose_cli": pp.to_cli(prose["params"]),
        "sdrf": pp.to_sdrf_technical(decisions),
        "llm": {"backend": backend.name if backend else None, "rejected": llm_rejected, "review": review},
        "jev": {"model": jev.name if jev else None, "cross_check": jev_issues, "study": jev_report,
                "calls": jev.calls if jev else 0, "cost_usd": round(jev.cost_usd, 5) if jev else 0},
        "consistency": consistency,
        "text_conflicts": llm_text_conflicts + (jev_report.get("text_rules", []) if jev_report else []),
        "timings_s": timings,
    }


STATUS = {"auto": "auto", "review": "REVIEW", "default": "default", "unknown": "unknown"}


def _fmt(v):
    if isinstance(v, dict) and "value" in v:
        return f"{v['value']} {v.get('unit', '')}".strip()
    if isinstance(v, dict):
        return ", ".join(f"{k}={v[k]}" for k in list(v)[:4] if k != "windows")
    if isinstance(v, list):
        return ", ".join(str(x) for x in v[:6])
    return str(v)


def markdown(rep):
    i = rep["input"]
    L_ = [f"# prose-auto report: {os.path.basename(i['path'])}", "",
          f"{i['instrument'] or 'unknown instrument'}, {i['ms1']} MS1 / {i['ms2']} MS2 / {i['ms3']} MS3 spectra. "
          f"Run time {rep['timings_s'].get('total')} s.", "",
          "| Parameter | Value | Confidence | Status | From |", "|---|---|---|---|---|"]
    order = ["acquisition", "instrument", "ms1_analyzer", "ms2_analyzer", "fragmentation", "precursor_tolerance",
             "precursor_offset", "fragment_tolerance", "fragment_offset", "charge_range", "label", "cys_alkylation",
             "enrichment", "enzyme", "organism", "variable_mods", "dia_windows", "collision_energy",
             "instrument_class", "ion_mobility"]
    ds = rep["decisions"]
    for k in order + [k for k in ds if k not in order]:
        if k in ds:
            d = ds[k]
            L_.append(f"| {d['label']} | {_fmt(d['value'])} | {d['confidence']:.2f} | {STATUS[d['status']]} | {d['source']} |")
    L_ += ["", "## Conflicts and notes", ""]
    any_c = False
    for k, d in ds.items():
        for c in d["conflicts"]:
            L_.append(f"- **{d['label']}**: {c}")
            any_c = True
    for n in rep["prose"]["notes"]:
        L_.append(f"- {n}")
        any_c = True
    for r in rep.get("consistency", []):
        L_.append(f"- Consistency ({r.get('severity')}): **{r.get('parameter')}** {r.get('issue')}")
        any_c = True
    for r in rep.get("text_conflicts", []):
        L_.append(f"- Texts disagree with each other ({r.get('kind', '')}): {r.get('message') or r.get('issue')}")
        any_c = True
    for r in rep.get("jev", {}).get("cross_check", []):
        L_.append(f"- Jev cross-check ({r.get('severity')}): **{r.get('parameter')}** {r.get('issue')}")
        any_c = True
    for r in rep["llm"]["review"]:
        L_.append(f"- LLM review ({r.get('severity')}): **{r.get('parameter')}** {r.get('issue')}")
        any_c = True
    if not any_c:
        L_.append("- none")
    L_ += ["", "## ProSE parameters", "", "```", "ProSE -in <mzML> -database <fasta> -out <idXML> \\",
           "  " + " ".join(f'"{x}"' if " " in x or not x else x for x in rep["prose_cli"]), "```", "",
           "## Evidence", ""]
    for k in order + [k for k in ds if k not in order]:
        if k not in ds:
            continue
        d = ds[k]
        L_.append(f"**{d['label']}** = {_fmt(d['value'])}")
        for e in d["evidence"]:
            q = f' — "{e["quote"]}"' if e.get("quote") else ""
            L_.append(f"- [{e['source']}/{e['method']}] {_fmt(e['value'])} ({e['confidence']:.2f}): {e['detail']}{q}")
        L_.append("")
    if rep["llm"]["rejected"]:
        L_ += ["## LLM statements rejected by the quote check", ""]
        for r in rep["llm"]["rejected"][:20]:
            L_.append(f"- {r.get('source', '')}: {r.get('parameter')} = {r.get('value')} ({r.get('why', r.get('error'))})")
    return "\n".join(L_) + "\n"


MODEL_PROTEOMES = {
    "Homo_sapiens": "UP000005640", "Mus_musculus": "UP000000589", "Rattus_norvegicus": "UP000002494",
    "Bos_taurus": "UP000009136", "Sus_scrofa": "UP000008227", "Danio_rerio": "UP000000437",
    "Drosophila_melanogaster": "UP000000803", "Caenorhabditis_elegans": "UP000001940",
    "Saccharomyces_cerevisiae": "UP000002311", "Arabidopsis_thaliana": "UP000006548",
    "Escherichia_coli": "UP000000625", "Bacillus_subtilis": "UP000001570",
}


def fetch_panel(out, reviewed=True, log=print):
    """Download the model-organism panel (UniProt reference proteomes)."""
    import urllib.request
    os.makedirs(out, exist_ok=True)
    for name, up in MODEL_PROTEOMES.items():
        q = f"(proteome:{up})" + ("+AND+(reviewed:true)" if reviewed else "")
        url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query={q}"
        dst = os.path.join(out, f"{name}.fasta")
        if os.path.getsize(dst) if os.path.exists(dst) else 0:
            log(f"  {name}: present, skipped")
            continue
        log(f"  {name} ({up}) -> {dst}")
        # UniProt's stream endpoint returns sporadic HTTP 500s: retry, and never leave a partial file
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=600) as r:
                    data = r.read()
                break
            except OSError as e:  # URLError and HTTPError are OSErrors
                log(f"    attempt {attempt + 1} failed: {e}")
                time.sleep(10 * 2 ** attempt)
        else:
            log(f"  {name}: giving up; the panel lacks this organism")
            continue
        with open(dst + ".part", "wb") as fh:
            fh.write(data)
        os.replace(dst + ".part", dst)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="prose-auto", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    an = sub.add_parser("analyze", help="propose search settings for one mzML file")
    an.add_argument("mzml")
    an.add_argument("--panel", help="directory of model-organism FASTA files (one per organism)")
    an.add_argument("--fasta", help="the user's database; searched in addition to the panel")
    an.add_argument("--description", help="free-text description of the experiment, or @file")
    an.add_argument("--methods", help="methods section (file or text)")
    an.add_argument("--accession", help="PXD/MSV accession: PRIDE, MassIVE and curated SDRF records")
    an.add_argument("--no-paper", action="store_true", help="do not fetch the paper from Europe PMC")
    an.add_argument("--llm", default="none", help="none | claude | codex | cli:<cmd> | openai:<model>")
    an.add_argument("--jev", action="store_true",
                    help="classify the texts with Jev and cross-check the settings (JEV_API_KEY)")
    an.add_argument("--llm-priors", action="store_true", help="ask the LLM for literature tolerance priors")
    an.add_argument("--llm-review", action="store_true", help="let the LLM review the final settings")
    an.add_argument("--prose", nargs="?", const="", help="confirm with a ProSE subset search (path or PATH)")
    an.add_argument("--open-search", action="store_true", help="with --prose: also run an open search")
    an.add_argument("--fastag", help="FASTag executable: genus call and Cys alkylation from sequence tags")
    an.add_argument("--fastag-subsample", type=float, default=0.0, help="fraction of spectra FASTag tags (0 = all)")
    an.add_argument("--diaspextract", help="DIAspeXtract executable: pseudo-MS/MS spectra for DIA / .d input")
    an.add_argument("--threads", type=int, default=0, help="threads for external tools (0 = their default)")
    an.add_argument("--spectra", type=int, default=2000, help="MS2 spectra for the mini-search")
    an.add_argument("--max-ms2", type=int, default=20000, help="MS2 spectra read for the detectors")
    an.add_argument("--no-search", action="store_true", help="skip the mini-search")
    an.add_argument("--skip-enzyme", action="store_true", help="skip the protease competition")
    an.add_argument("-o", "--out", default="prose_auto_out")
    an.add_argument("--json", action="store_true", help="print the report JSON to stdout (for the GUI)")
    pn = sub.add_parser("panel", help="download the model-organism FASTA panel from UniProt")
    pn.add_argument("out")
    pn.add_argument("--all-entries", action="store_true", help="TrEMBL too (default: Swiss-Prot only)")
    a = ap.parse_args(argv)
    if a.cmd == "panel":
        fetch_panel(a.out, reviewed=not a.all_entries)
    if a.cmd == "analyze":
        log = (lambda *x: print(*x, file=sys.stderr)) if a.json else print
        rep = analyze(a, log)
        os.makedirs(a.out, exist_ok=True)
        base = os.path.join(a.out, os.path.splitext(os.path.basename(a.mzml))[0])
        with open(base + ".auto.json", "w") as fh:
            json.dump(rep, fh, indent=1, default=str)
        with open(base + ".auto.md", "w") as fh:
            fh.write(markdown(rep))
        with open(base + ".prose.ini", "w") as fh:
            fh.write(pp.to_ini_snippet(rep["prose"]["params"]))
        with open(base + ".sdrf_technical.tsv", "w") as fh:
            fh.write("\t".join(rep["sdrf"]) + "\n" + "\t".join(rep["sdrf"].values()) + "\n")
        if a.json:
            print(json.dumps(rep, default=str))
        else:
            print(markdown(rep).split("## Evidence")[0])
            print(f"written: {base}.auto.json, .auto.md, .prose.ini, .sdrf_technical.tsv")


if __name__ == "__main__":
    main()
