"""Rule-based extraction of search settings from free text.

Works without any LLM on a user's description ("HeLa digest, TMTpro 16plex,
IAA alkylation, Exploris 480, DIA"), a methods section, or the PRIDE
sample/data processing protocols. Every hit keeps the verbatim span it came
from, so the GUI can show it and the LLM layer can be checked against it.
"""

from __future__ import annotations

import re

from .evidence import Evidence

ORGANISMS = {
    "Homo sapiens": r"\bhomo sapiens\b|\bhuman\b|\bhela\b|\bhek ?293|\bjurkat\b|\ba549\b|\bmcf-?7\b|\bk562\b|\bu2os\b|\bpatients?\b|\bhuman (plasma|serum)\b",
    "Mus musculus": r"\bmus musculus\b|\bmouse\b|\bmice\b|\bmurine\b",
    "Rattus norvegicus": r"\brattus\b|\brats?\b",
    "Saccharomyces cerevisiae": r"\bsaccharomyces\b|\byeast\b|\bs\. ?cerevisiae\b",
    "Escherichia coli": r"\bescherichia\b|\be\. ?coli\b",
    "Arabidopsis thaliana": r"\barabidopsis\b",
    "Drosophila melanogaster": r"\bdrosophila\b",
    "Caenorhabditis elegans": r"\bc\. ?elegans\b|\bcaenorhabditis\b",
    "Danio rerio": r"\bzebrafish\b|\bdanio\b",
    "Bos taurus": r"\bbos taurus\b|\bbovine\b|\bbsa\b",
    "Sus scrofa": r"\bporcine\b|\bpigs?\b|\bsus scrofa\b",
    "Plasmodium falciparum": r"\bplasmodium\b",
    "Mycobacterium tuberculosis": r"\bmycobacterium tuberculosis\b|\bm\. ?tuberculosis\b",
}

ENZYMES = {
    "Trypsin/P": r"\btrypsin/p\b",
    "Lys-C": r"\blys-?c\b(?!.*trypsin)",
    "Trypsin": r"\btrypsin\b|\btryptic\b",
    "Chymotrypsin": r"\bchymotrypsin\b",
    "Glu-C": r"\bglu-?c\b",
    "Asp-N": r"\basp-?n\b",
    "Arg-C": r"\barg-?c\b",
    "unspecific cleavage": r"\bno enzyme\b|\bnon-?specific\b|\bunspecific\b|\bhla\b|\bimmunopeptid",
}

ALKYLATION = {
    "Carbamidomethyl (C)": r"\biodoacetamide\b|\bchloroacetamide\b|\biaa\b|\bcaa\b|\bcarbamidomethyl",
    "Propionamide (C)": r"\bacrylamide\b|\bpropionamide\b",
    "Nethylmaleimide (C)": r"\bn-?ethylmaleimide\b|\bnem\b",
    "Methylthio (C)": r"\bmmts\b|\bmethyl methanethiosulfonate\b|\bmethylthio\b",
    "Carboxymethyl (C)": r"\biodoacetic acid\b|\bcarboxymethyl",
    "none (free Cys)": r"\bno alkylation\b|\bnot alkylated\b|\bwithout alkylation\b|\bnon-?alkylated\b",
}

LABELS = {
    "TMTpro (16/18plex)": r"\btmt ?pro\b|\btmt ?1[68] ?-?plex\b",
    "TMT10/11plex": r"\btmt ?1[01] ?-?plex\b",
    "TMT6plex": r"\btmt ?6 ?-?plex\b|\btmt sixplex\b",
    "TMT2plex": r"\btmt ?2 ?-?plex\b|\btmt duplex\b",
    "iTRAQ8plex": r"\bitraq ?8 ?-?plex\b|\b8-?plex itraq\b",
    "iTRAQ4plex": r"\bitraq ?4 ?-?plex\b|\b4-?plex itraq\b|\bitraq\b",
    "SILAC": r"\bsilac\b|\blys ?8\b|\barg ?10\b|\bheavy (lysine|arginine)\b",
    "dimethyl": r"\bdimethyl label",
    "label free": r"\blabel[- ]free\b|\blfq\b",
    "TMT": r"\btmt\b|\btandem mass tag",
}

ENRICHMENT = {
    "Phospho (STY)": r"\bphospho|\btio2\b|\bimac\b|\bfe-?nta\b|\bti4\+",
    "GlyGly (K)": r"\bk-?[eε]?-?gg\b|\bdi-?gly|\bubiquitin(ome|ation)\b",
    "Acetyl (K)": r"\bacetyl-?lysine\b|\bacetylome\b",
    "Glyco": r"\bglyco(peptide|proteom)",
}

ACQUISITION = {
    "DIA (diaPASEF)": r"\bdia-?pasef\b",
    "DIA": r"\bdata[- ]independent\b|\bdia\b|\bswath\b",
    "PRM": r"\bprm\b|\bparallel reaction monitoring\b",
    "DDA": r"\bdata[- ]dependent\b|\bdda\b|\btop ?\d+\b|\bdda-?pasef\b",
}

FRAGMENTATION = {
    "EThcD": r"\bethcd\b", "ETD": r"\betd\b|\belectron transfer", "HCD": r"\bhcd\b|\bhigher-energy collision",
    "CID": r"\bcid\b|\bcollision-induced\b", "UVPD": r"\buvpd\b",
}

INSTRUMENTS = [
    "Orbitrap Astral", "Orbitrap Exploris 480", "Orbitrap Exploris 240", "Orbitrap Exploris 120",
    "Orbitrap Eclipse", "Orbitrap Fusion Lumos", "Orbitrap Fusion", "Orbitrap Ascend",
    "Q Exactive HF-X", "Q Exactive HF", "Q Exactive Plus", "Q Exactive", "LTQ Orbitrap Velos",
    "LTQ Orbitrap Elite", "LTQ Orbitrap XL", "timsTOF Pro 2", "timsTOF Pro", "timsTOF HT", "timsTOF SCP",
    "timsTOF Ultra", "ZenoTOF 7600", "TripleTOF 6600", "TripleTOF 5600", "Synapt G2", "Xevo G2",
]

_TOL = re.compile(r"(precursor|ms1|parent|peptide|fragment|ms/?ms|ms2|product)[^.;]{0,60}?"
                  r"(?:tolerance|accuracy|mass deviation|window)[^.;]{0,30}?(\d+(?:\.\d+)?)\s*(ppm|da|mmu|th)",
                  re.I)


def _span(text, m, pad=40):
    a, b = max(0, m.start() - pad), min(len(text), m.end() + pad)
    return text[a:b].replace("\n", " ").strip()


def extract(text: str, source="text", base_conf=0.55) -> list:
    """Evidence from free text. source is 'user' for the user's own
    description (more trusted) or 'text' for methods/protocols."""
    if not text:
        return []
    ev = []
    low = text.lower()

    def first(table, parameter, conf, single=True):
        hits = []
        for value, pat in table.items():
            m = re.search(pat, low)
            if m:
                hits.append((m.start(), value, m))
        hits.sort()
        for _, value, m in hits[:1 if single else len(hits)]:
            ev.append(Evidence(parameter, value, conf, source, "keywords",
                               f"text mentions '{text[m.start():m.end()]}'", quote=_span(text, m)))
        return hits

    first(ORGANISMS, "organism", base_conf, single=False)
    first(ENZYMES, "enzyme", base_conf)
    first(ALKYLATION, "cys_alkylation", base_conf)
    lab = first(LABELS, "label", base_conf)
    first(ENRICHMENT, "enrichment", base_conf * 0.9, single=False)
    first(ACQUISITION, "acquisition", base_conf)
    first(FRAGMENTATION, "fragmentation", base_conf * 0.9)
    for ins in INSTRUMENTS:
        m = re.search(re.escape(ins.lower()).replace(r"\ ", r"[\s-]?"), low)
        if m:
            ev.append(Evidence("instrument", ins, base_conf, source, "keywords",
                               f"text names the instrument {ins}", quote=_span(text, m)))
            break
    for m in _TOL.finditer(text):
        kind, val, unit = m.group(1).lower(), float(m.group(2)), m.group(3).lower()
        if unit == "mmu":
            val, unit = val / 1000, "da"
        param = "precursor_tolerance" if kind in ("precursor", "ms1", "parent", "peptide") else "fragment_tolerance"
        ev.append(Evidence(param, {"value": val, "unit": "ppm" if unit == "ppm" else "Da"},
                           base_conf * 0.8, source, "tolerance pattern",
                           f"text gives a {kind} tolerance of {val} {unit}", quote=_span(text, m, 10)))
    return ev
