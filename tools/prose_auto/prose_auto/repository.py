"""Evidence from public repositories: PRIDE, MassIVE (PROXI), curated SDRF.

* PRIDE Archive API v3 (project record): organisms, instruments, identified
  PTMs, quantification methods, the free-text sample and data processing
  protocols, and the publication (PubMed id / DOI).
* MassIVE / ProteomeXchange via PROXI: species, instruments, modifications.
* bigbio/sdrf-annotated-datasets on GitHub: the curated SDRF of the dataset
  (13k accessions in 2026), which carries the search-relevant comment
  columns: label, instrument, modification parameters, cleavage agent,
  dissociation method, precursor/fragment mass tolerance, acquisition method.
* Europe PMC full text of the publication, for its methods section.

Repository metadata is deposited by people and is sometimes wrong (the
bigbio curators keep a list of PRIDE instrument fields that contradict the
raw files). It is evidence to check against the data, not ground truth.
"""

from __future__ import annotations

import csv
import io
import json
import re
import urllib.request

from .evidence import Evidence

PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{acc}"
PROXI = "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/{acc}"
SDRF = "https://raw.githubusercontent.com/bigbio/sdrf-annotated-datasets/master/datasets/{acc}/{acc}.sdrf.tsv"
EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&resultType=lite"
EPMC_FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"


def _get(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "prose-auto/0.1",
                                               "Accept": "application/json, text/plain, */*"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def fetch(url, log=print):
    try:
        return _get(url)
    except Exception as e:  # network policies, 404s
        log(f"  could not fetch {url}: {e}")
        return None


# --- SDRF -------------------------------------------------------------------

def _kv(cell):
    """'NT=Oxidation;MT=Variable;TA=M;AC=Unimod:35' -> dict"""
    out = {}
    for part in cell.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip().upper()] = v.strip()
    return out


def _tol(cell):
    m = re.match(r"\s*([\d.]+)\s*(ppm|da|mmu|th)?", cell or "", re.I)
    if not m:
        return None
    v, u = float(m.group(1)), (m.group(2) or "ppm").lower()
    if u == "mmu":
        return {"value": v / 1000.0, "unit": "Da"}
    return {"value": v, "unit": "ppm" if u == "ppm" else "Da"}


def sdrf_evidence(text, source_name="SDRF") -> list:
    rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
    if len(rows) < 2:
        return []
    head = [h.strip().lower() for h in rows[0]]
    body = rows[1:]

    def col(name):
        idx = [i for i, h in enumerate(head) if h == name]
        vals = []
        for r in body:
            for i in idx:
                if i < len(r) and r[i].strip() and r[i].strip().lower() not in ("not available", "not applicable", "na"):
                    vals.append(r[i].strip())
        return vals

    ev = []
    conf = 0.7
    src = f"{source_name} (curated SDRF)"

    def most(vals):
        return max(set(vals), key=vals.count) if vals else None

    org = most(col("characteristics[organism]"))
    if org:
        ev.append(Evidence("organism", org, conf, "repository", "sdrf", f"{src}: organism {org}"))
    labels = {_kv(l).get("NT", l) for l in col("comment[label]")}
    if labels:
        lab = "label free" if any("label free" in l.lower() for l in labels) else \
            ("TMT" if any(l.upper().startswith("TMT") for l in labels) else
             ("iTRAQ" if any("itraq" in l.lower() for l in labels) else
              ("SILAC" if any("silac" in l.lower() for l in labels) else ", ".join(sorted(labels)))))
        if lab.startswith("TMT"):
            n = len(labels)
            lab = "TMTpro (16/18plex)" if n > 11 else ("TMT10/11plex" if n > 6 else f"TMT{n}plex")
        ev.append(Evidence("label", lab, conf, "repository", "sdrf", f"{src}: labels {sorted(labels)[:6]}"))
    ins = most(col("comment[instrument]"))
    if ins:
        ev.append(Evidence("instrument", _kv(ins).get("NT", ins), conf, "repository", "sdrf", f"{src}: {ins}"))
    acq = most(col("comment[proteomics data acquisition method]"))
    if acq:
        nt = _kv(acq).get("NT", acq).lower()
        mode = "DIA" if "independent" in nt or "dia" in nt else ("DDA" if "dependent" in nt or "dda" in nt else nt)
        ev.append(Evidence("acquisition", mode, conf, "repository", "sdrf", f"{src}: {acq}"))
    dis = most(col("comment[dissociation method]"))
    if dis:
        ev.append(Evidence("fragmentation", _kv(dis).get("NT", dis).upper(), conf, "repository", "sdrf", f"{src}: {dis}"))
    enz = most(col("comment[cleavage agent details]"))
    if enz:
        ev.append(Evidence("enzyme", _kv(enz).get("NT", enz), conf, "repository", "sdrf", f"{src}: {enz}"))
    for name, param in (("comment[precursor mass tolerance]", "precursor_tolerance"),
                        ("comment[fragment mass tolerance]", "fragment_tolerance")):
        t = _tol(most(col(name)) or "")
        if t:
            ev.append(Evidence(param, t, 0.5, "repository", "sdrf",
                               f"{src}: {most(col(name))} (what the submitters searched with)"))
    mods = sorted(set(col("comment[modification parameters]")))
    fixed, var = [], []
    for m in mods:
        kv = _kv(m)
        name, ta, mt = kv.get("NT", ""), kv.get("TA", ""), kv.get("MT", "").lower()
        pp = kv.get("PP", "")
        site = ta if ta else ("N-term" if "n-term" in pp.lower() else pp)
        label = f"{name} ({site})" if site else name
        (fixed if mt == "fixed" else var).append(label)
    cys = [f for f in fixed if f.endswith("(C)")]
    if cys:
        ev.append(Evidence("cys_alkylation", cys[0], conf, "repository", "sdrf", f"{src}: fixed {cys[0]}"))
    if var:
        ev.append(Evidence("variable_mods", var, 0.5, "repository", "sdrf", f"{src}: variable {var}"))
    enr = [v for v in var if v.lower().startswith(("phospho", "glygly", "gg", "acetyl (k)"))]
    if enr:
        ev.append(Evidence("enrichment", enr[0], 0.4, "repository", "sdrf",
                           f"{src}: {enr[0]} among the variable modifications"))
    return ev


# --- PRIDE / MassIVE -----------------------------------------------------------

def pride_evidence(acc, log=print):
    txt = fetch(PRIDE.format(acc=acc), log)
    if not txt:
        return [], {}
    try:
        p = json.loads(txt)
    except json.JSONDecodeError:
        return [], {}
    ev = []
    src = f"PRIDE {acc}"
    for o in p.get("organisms", [])[:3]:
        ev.append(Evidence("organism", o.get("name"), 0.6, "repository", "pride", f"{src}: organism {o.get('name')}"))
    for i in p.get("instruments", [])[:2]:
        ev.append(Evidence("instrument", i.get("name"), 0.6, "repository", "pride", f"{src}: instrument {i.get('name')}"))
    for q in p.get("quantificationMethods", [])[:2]:
        name = q.get("name", "")
        ev.append(Evidence("label", name, 0.4, "repository", "pride", f"{src}: quantification {name}"))
    ptms = [m.get("name") for m in p.get("identifiedPTMStrings", [])]
    if ptms:
        ev.append(Evidence("variable_mods", ptms, 0.3, "repository", "pride", f"{src}: identified PTMs {ptms}"))
    texts = {
        "title": p.get("title", ""),
        "description": p.get("projectDescription", ""),
        "sample_protocol": p.get("sampleProcessingProtocol", ""),
        "data_protocol": p.get("dataProcessingProtocol", ""),
        "references": p.get("references", []),
    }
    return ev, texts


def proxi_evidence(acc, log=print):
    txt = fetch(PROXI.format(acc=acc), log)
    if not txt:
        return [], {}
    try:
        p = json.loads(txt)
    except json.JSONDecodeError:
        return [], {}
    ev = []
    src = f"PROXI {acc}"
    for sp in p.get("species", []):
        # [{"terms": [{"name": "taxonomy: scientific name", "value": "Homo sapiens"}, ...]}]
        # or, from some PROXI servers, a bare list of terms per species
        terms = sp if isinstance(sp, list) else (sp.get("terms") or [sp])
        for term in terms:
            if not isinstance(term, dict):
                continue
            if "scientific name" in str(term.get("name", "")).lower() and term.get("value"):
                ev.append(Evidence("organism", term["value"], 0.6, "repository", "proxi", f"{src}: species"))
                break
    for ins in p.get("instruments", []):
        if not isinstance(ins, dict):
            continue
        ev.append(Evidence("instrument", ins.get("name"), 0.6, "repository", "proxi", f"{src}: instrument"))
    mods = [m.get("name") for m in p.get("modifications", []) if isinstance(m, dict) and m.get("name")]
    if mods:
        ev.append(Evidence("variable_mods", mods, 0.3, "repository", "proxi", f"{src}: modifications {mods}"))
    return ev, {"title": p.get("title", ""), "description": p.get("summary", ""),
                "references": p.get("publications", [])}


def methods_from_europepmc(pmid_or_doi, log=print):
    """The methods section of an open-access paper (Europe PMC JATS XML)."""
    import urllib.parse
    import xml.etree.ElementTree as ET
    q = urllib.parse.quote(f'EXT_ID:{pmid_or_doi}' if str(pmid_or_doi).isdigit() else f'DOI:"{pmid_or_doi}"')
    txt = fetch(EPMC_SEARCH.format(q=q), log)
    if not txt:
        return ""
    hits = json.loads(txt).get("resultList", {}).get("result", [])
    pmcid = next((h.get("pmcid") for h in hits if h.get("pmcid")), None)
    if not pmcid:
        return ""
    xml = fetch(EPMC_FULLTEXT.format(pmcid=pmcid), log)
    if not xml:
        return ""
    root = ET.fromstring(xml)
    out = []
    for sec in root.iter("sec"):
        title = "".join(sec.findtext("title") or "").lower()
        if sec.get("sec-type", "").startswith("method") or re.search(r"method|experimental|procedure", title):
            out.append(" ".join("".join(sec.itertext()).split()))
    return "\n\n".join(out)


def repository_evidence(acc, log=print, fetch_paper=True):
    """All repository evidence for a PXD/MSV accession, plus the texts
    (protocols, methods) for text_extract / llm."""
    ev, texts = [], {}
    if acc.upper().startswith("PXD"):
        e, t = pride_evidence(acc, log)
        ev += e
        texts.update(t)
    if acc.upper().startswith(("MSV", "PXD")) and not ev:
        e, t = proxi_evidence(acc, log)
        ev += e
        texts.update(t)
    sdrf = fetch(SDRF.format(acc=acc.upper()), log)
    if sdrf:
        ev += sdrf_evidence(sdrf, f"bigbio {acc}")
        texts["sdrf"] = sdrf
    if fetch_paper:
        for ref in texts.get("references", []) or []:
            key = ref.get("pubmedId") or ref.get("doi") if isinstance(ref, dict) else None
            if key:
                m = methods_from_europepmc(key, log)
                if m:
                    texts["methods"] = m
                    break
    return ev, texts
