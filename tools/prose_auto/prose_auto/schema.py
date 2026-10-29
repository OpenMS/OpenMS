"""The ProSE parameter heads read from text, and the rules that tie them together.

One schema serves every text reader: Jev gets each head as a typed question
(one request answers all heads), the generative LLM gets the same heads and
values in its prompt, and the reconciled decisions are checked against the
same rules. A head's options carry the prose-auto value they map to, so text
answers become ordinary Evidence.

Heads are mutually exclusive inside (one choice per head). Rules capture what
is exclusive or implausible *across* heads, e.g. TMT with ion-trap CID MS2
(reporter ions fall below the low-mass cut-off), diaPASEF on an Orbitrap, an
immunopeptidome digested with trypsin. A rule is

    hard   the combination cannot be right; one of the values is wrong
    soft   possible but unusual; worth a look before a long search
    imply  one value requires another (filled in, or flagged when it contradicts)

`decode()` picks the most probable combination of head answers that breaks no
hard rule, so a near-tie in one head is settled by the others.
"""

from __future__ import annotations

import itertools
import math
import re
from dataclasses import dataclass, field

from . import text_extract as tx

NOT_STATED = "not_stated"


@dataclass
class Head:
    name: str
    question: str
    options: dict                      # key -> (value, description)
    parameter: str | None = None       # Evidence parameter it feeds; None = context only
    kind: str = "choice"               # choice | flag (yes/no)
    hint: str = ""                     # extra instruction (normalisation, priorities)
    cap: float = 0.6                   # max evidence confidence from text


def _opts(pairs):
    return {re.sub(r"[^a-z0-9]+", "_", v.lower()).strip("_")[:60]: (v, d) for v, d in pairs}


ISOBARIC = {"TMT6plex", "TMT10/11plex", "TMTpro (16/18plex)", "TMT2plex", "iTRAQ4plex", "iTRAQ8plex"}
SILAC = {"SILAC (Lys8/Arg10)", "SILAC (Arg6)", "SILAC Lys4 or dimethyl (d4)"}
METABOLIC = {"15N metabolic labelling", "heavy water (2H2O) labelling"}
TIMS = ("timsTOF",)
NO_ETD = ("Q Exactive", "Exploris", "Astral", "TripleTOF", "ZenoTOF", "timsTOF", "LTQ Orbitrap XL", "Synapt", "Xevo")

def _pride_species():
    """Species list covering 95% of PRIDE raw files (data/organisms_pride.tsv)."""
    import os
    f = os.path.join(os.path.dirname(__file__), "data", "organisms_pride.tsv")
    try:
        rows = [l.split("\t") for l in open(f) if l.strip() and not l.startswith(("#", "rank"))]
        return [r[1] for r in rows]
    except OSError:
        return []


PRIDE_SPECIES = _pride_species()

# instrument options: the models of all bigbio SDRFs (>= 5 datasets, spellings merged) with the aliases methods
# texts use; tuned on the text-only benchmark (dev split, rounds 1-4)
INSTRUMENT_OPTIONS = [
    ("Q Exactive", "Q Exactive: QE, Q-Exactive, QExactive (not Plus/HF/HF-X)"),
    ("Orbitrap Fusion Lumos", "Orbitrap Fusion Lumos: Lumos, Fusion Lumos Tribrid"),
    ("Q Exactive HF", "Q Exactive HF: QE HF, Q-Exactive HF, QE-HF (not HF-X)"),
    ("Orbitrap Exploris 480", "Orbitrap Exploris 480: Exploris 480"),
    ("Orbitrap Fusion", "Orbitrap Fusion: Fusion Tribrid (not Lumos)"),
    ("LTQ Orbitrap Velos", "LTQ Orbitrap Velos: Orbitrap Velos, LTQ-Orbitrap Velos (not Velos Pro)"),
    ("Q Exactive Plus", "Q Exactive Plus: QE Plus, Q-Exactive Plus, QE+"),
    ("Q Exactive HF-X", "Q Exactive HF-X: QE HF-X, Q-Exactive HFX, QE-HFX"),
    ("LTQ Orbitrap", "LTQ Orbitrap: LTQ-Orbitrap, LTQ Orbitrap Classic or Discovery"),
    ("Orbitrap Eclipse", "Orbitrap Eclipse: Eclipse Tribrid"),
    ("LTQ Orbitrap Elite", "LTQ Orbitrap Elite: Orbitrap Elite"),
    ("Orbitrap Elite", "Orbitrap Elite: Orbitrap Elite"),
    ("timsTOF Pro", "timsTOF Pro: timsTOF Pro (not Pro 2)"),
    ("Orbitrap Astral", "Orbitrap Astral: Astral (not Astral Zoom)"),
    ("TripleTOF 5600", "TripleTOF 5600: TripleTOF 5600, 5600+, AB SCIEX 5600"),
    ("LTQ Orbitrap XL", "LTQ Orbitrap XL: Orbitrap XL, LTQ-Orbitrap XL"),
    ("Orbitrap Exploris 240", "Orbitrap Exploris 240: Exploris 240"),
    ("timsTOF Pro 2", "timsTOF Pro 2: timsTOF Pro 2"),
    ("TripleTOF 6600", "TripleTOF 6600: TripleTOF 6600, 6600+"),
    ("Orbitrap Velos Pro", "Orbitrap Velos Pro: LTQ Orbitrap Velos Pro, Velos Pro"),
    ("LTQ", "LTQ: LTQ linear ion trap, LTQ XL, Velos Pro ion trap alone"),
    ("timsTOF HT", "timsTOF HT: timsTOF HT"),
    ("LTQ FT", "LTQ FT: LTQ-FT, LTQ FT Ultra, FT-ICR"),
    ("Orbitrap Ascend", "Orbitrap Ascend: Ascend Tribrid"),
    ("timsTOF SCP", "timsTOF SCP: timsTOF SCP"),
    ("Orbitrap Fusion ETD", "Orbitrap Fusion ETD: Orbitrap Fusion ETD"),
    ("maXis", "maXis: Bruker maXis, maXis II, impact II Q-TOF"),
    ("timsTOF fleX", "timsTOF fleX: timsTOF fleX (MALDI-2)"),
    ("ZenoTOF 7600", "ZenoTOF 7600: ZenoTOF 7600"),
    ("Orbitrap Astral Zoom", "Orbitrap Astral Zoom: Astral Zoom"),
    ("TripleTOF 5600+", "TripleTOF 5600+: TripleTOF 5600+"),
    ("autoflex", "autoflex: autoflex"),
    ("qstar", "qstar: qstar"),
    ("impact II", "impact II: impact II"),
    ("LTQ Velos", "LTQ Velos: LTQ Velos"),
    ("Orbitrap", "Orbitrap: Orbitrap"),
    ("timsTOF", "timsTOF: timsTOF"),
    ("6520A Quadrupole Time-of-Flight LC/MS", "6520A Quadrupole Time-of-Flight LC/MS: 6520A Quadrupole Time-of-Flight LC/MS"),
    ("ultraflex", "ultraflex: ultraflex"),
    ("timsTOF Ultra 2", "timsTOF Ultra 2: timsTOF Ultra 2"),
    ("timsTOF Ultra", "timsTOF Ultra: timsTOF Ultra (and Ultra 2)"),
    ("SYNAPT G2-Si", "SYNAPT G2-Si: SYNAPT G2-Si"),
    ("Synapt G2 HDMS", "Synapt G2 HDMS: Synapt G2 HDMS"),
    ("amazon speed etd", "amazon speed etd: amazon speed etd"),
    ("LTQ FT Ultra", "LTQ FT Ultra: LTQ FT Ultra"),
    ("LTQ Orbitrap XL ETD", "LTQ Orbitrap XL ETD: LTQ Orbitrap XL ETD"),
    ("Synapt G2 MS", "Synapt G2 MS: Synapt G2 MS"),
    ("LTQ Orbitrap Discovery", "LTQ Orbitrap Discovery: LTQ Orbitrap Discovery"),
    ("Orbitrap Exploris 120", "Orbitrap Exploris 120: Exploris 120"),
    ("Synapt G2", "Synapt G2: Synapt G2, G2-S, G2-Si (Waters)"),
    ("Xevo G2", "Xevo G2: Xevo G2-XS QTof"),
    ("other instrument", "a mass spectrometer not in this list"),
]


HEADS = [
    Head("organism", "Which species do the analysed proteins come from (the database to search)?",
         _opts([(o, o) for o in dict.fromkeys(list(tx.ORGANISMS) + PRIDE_SPECIES)] + [
             ("other species", "a species not in this list"),
             ("mixed species", "proteins of several species are in the SAME sample and searched together on purpose: a human/yeast/E. coli benchmark mix, xenograft, host-pathogen infection, microbiome. NOT for a study whose different samples come from different species, NOT for contaminants, BSA, trypsin or spike-in standard peptides")]),
         "organism", hint="Cell lines imply their species (HeLa, HEK293 = human). Ignore reagent proteins "
                          "(BSA, trypsin) and contaminant databases. " + "If one species dominates the samples, answer that species."),
    Head("sample_type", "What kind of sample was analysed?",
         _opts([("cell line", "cultured cells"), ("tissue", "tissue or organ"),
                ("body fluid", "plasma, serum, urine, CSF, saliva"),
                ("microbial culture", "a single bacterial, yeast or other microbial culture"),
                ("microbiome", "a microbial community: gut, soil, water"),
                ("immunopeptidome", "HLA/MHC-bound peptides, no protease digestion"),
                ("purified protein", "a purified protein, complex or standard mix")])),
    Head("enzyme", "Which protease produced the peptides?",
         _opts([("Trypsin", "trypsin, including Lys-C + trypsin double digests"),
                ("Trypsin/P", "trypsin explicitly allowing cleavage before proline"),
                ("Lys-C", "Lys-C alone, without trypsin"), ("Arg-C", "Arg-C"),
                ("Chymotrypsin", "chymotrypsin"), ("Glu-C", "Glu-C (V8)"), ("Asp-N", "Asp-N"),
                ("unspecific cleavage", "no protease or unspecific search: immunopeptides, endogenous peptides")]),
         "enzyme", hint="Lys-C followed by trypsin, or a Lys-C/trypsin mix, is Trypsin."),
    Head("missed_cleavages", "How many missed cleavages did the authors allow?",
         _opts([("0", "none"), ("1", "one"), ("2", "two"), ("3+", "three or more")])),
    Head("cys_alkylation", "How were cysteines alkylated (reagent or the resulting modification)?",
         _opts([("Carbamidomethyl (C)", "iodoacetamide (IAA) or chloroacetamide (CAA), carbamidomethylation"),
                ("Propionamide (C)", "acrylamide, propionamide"),
                ("Nethylmaleimide (C)", "N-ethylmaleimide (NEM)"),
                ("Methylthio (C)", "MMTS, methyl methanethiosulfonate"),
                ("Carboxymethyl (C)", "iodoacetic acid, carboxymethylation"),
                ("none (free Cys)", "explicitly not reduced or not alkylated")]),
         "cys_alkylation", hint="Reduction alone (DTT, TCEP) says nothing about alkylation."),
    Head("label", "Which quantification labelling was used?",
         _opts([("label free", "no quantification label: label-free quantification, spectral counting, or no quantification. Proximity labelling (biotin, APEX, BioID), crosslinking and isotope-labelled standard peptides are label free"),
                ("TMT2plex", "TMT duplex"), ("TMT6plex", "TMT sixplex"),
                ("TMT10/11plex", "TMT 10plex or 11plex"), ("TMTpro (16/18plex)", "TMTpro 16plex or 18plex"),
                ("iTRAQ4plex", "iTRAQ 4plex"), ("iTRAQ8plex", "iTRAQ 8plex"),
                ("SILAC (Lys8/Arg10)", "SILAC with heavy lysine +8 and arginine +10"),
                ("SILAC (Arg6)", "SILAC with arginine +6 (and/or lysine +6/+4, medium channels)"),
                ("SILAC Lys4 or dimethyl (d4)", "dimethyl labelling or lysine +4"),
                ("15N metabolic labelling", "15N-labelled organism or medium"),
                ("heavy water (2H2O) labelling", "D2O / heavy water turnover labelling")]),
         "label", hint="Label-free spike-ins of heavy standard peptides are still label free. " + "'TMT6plex' is also the modification name used for TMT 10/11plex: without a plex count or channel names, TMT10/11plex and TMT6plex are equally likely; prefer the plex the texts count."),
    Head("enrichment", "Which PTM enrichment was performed before MS?",
         _opts([("none", "whole proteome, no PTM enrichment"),
                ("Phospho (STY)", "phosphopeptide enrichment: TiO2, IMAC, Fe-NTA, Ti4+"),
                ("Phospho (Y)", "phosphotyrosine immunoprecipitation"),
                ("GlyGly (K)", "K-GG / di-glycine remnant (ubiquitin) enrichment"),
                ("Acetyl (K)", "acetyl-lysine immunoprecipitation"),
                ("Glyco", "glycopeptide enrichment: HILIC, lectins, hydrazide"),
                ("crosslinking", "chemical cross-linking (DSSO, BS3, DSS)")]),
         "enrichment", hint="A modification that was only searched for is not an enrichment."),
    Head("acquisition", "Which acquisition mode was used?",
         _opts([("DDA", "data-dependent acquisition (top-N, top speed)"),
                ("DDA SPS-MS3", "DDA with synchronous precursor selection MS3 for TMT reporters"),
                ("ddaPASEF", "DDA on a timsTOF (PASEF)"),
                ("DIA", "data-independent acquisition, SWATH"),
                ("DIA (diaPASEF)", "DIA on a timsTOF"),
                ("PRM", "parallel or selected reaction monitoring of targets")]),
         "acquisition"),
    Head("fragmentation", "Which MS/MS fragmentation was used?",
         _opts([("HCD", "HCD, beam-type CID"), ("CID", "ion-trap CID"), ("ETD", "ETD"),
                ("EThcD", "EThcD"), ("UVPD", "UVPD")]), "fragmentation"),
    Head("ms2_analyzer", "Which analyzer recorded the MS/MS spectra?",
         _opts([("ORBITRAP", "Orbitrap"), ("IT", "linear or 3D ion trap"), ("TOF", "time of flight")]),
         "ms2_analyzer"),
    Head("instrument", "Which mass spectrometer was used?", _opts(INSTRUMENT_OPTIONS), "instrument",
         hint="Name a model only if a model name appears in the texts; a vendor or 'Orbitrap mass spectrometer' alone is not_stated. 'CID in the LTQ' means an LTQ Orbitrap, never a Q Exactive."),
    Head("precursor_tolerance", "Which precursor mass tolerance did the authors search with?",
         _opts([("<=5 ppm", "5 ppm or less"), ("10 ppm", "6 to 15 ppm"), ("20 ppm", "16 to 30 ppm"),
                ("50 ppm", "31 to 100 ppm"), (">=0.5 Da", "0.5 Da or wider (low resolution or open search)")]),
         "precursor_tolerance", cap=0.4, hint="MaxQuant reports a first-search and a main-search tolerance; answer the main search."),
    Head("fragment_tolerance", "Which fragment mass tolerance did the authors search with?",
         _opts([(k, v) for k, v in {"<=10 ppm": "0.01 Da or less / up to 10 ppm", "20 ppm": "0.011-0.025 Da / 11-30 ppm", "0.05 Da": "0.03-0.1 Da / 31-100 ppm", "0.5 Da": "0.2-0.6 Da", ">=0.8 Da": "0.7 Da or wider"}.items()]),
         "fragment_tolerance", cap=0.4, hint="Convert a Da value with the Da limits of the options; do not compare Da numbers with ppm."),
]

# yes/no heads: variable modifications the authors searched, and search types ProSE does not cover
VARIABLE_MOD_FLAGS = {
    "Oxidation (M)": "methionine oxidation",
    "Acetyl (Protein N-term)": "protein N-terminal acetylation",
    "Deamidated (N)": "deamidation of asparagine/glutamine",
    "Gln->pyro-Glu (N-term Q)": "pyro-glutamate from N-terminal glutamine",
    "Carbamyl (N-term)": "carbamylation (urea artefact)",
    "Phospho (S)": "phosphorylation of serine/threonine/tyrosine",
    "Methyl (K)": "lysine/arginine methylation",
}
SCOPE_FLAGS = {
    "open_search": "an open or mass-tolerant search for unknown modifications",
    "de_novo": "de novo sequencing",
    "top_down": "intact proteins (top-down), not peptides",
    "crosslink_search": "a cross-link search engine (XlinkX, pLink, MeroX, xiSEARCH)",
    "glycan_search": "an intact glycopeptide search (Byonic, pGlyco, MSFragger-Glyco)",
    "sample_multiplexed_species": "samples combining several species (spike-in, xenograft, two-proteome benchmark)",
}

HEAD_BY_NAME = {h.name: h for h in HEADS}

# representative values for tolerance buckets, used for evidence and rule checks
TOL_VALUE = {
    "<=5 ppm": {"value": 5.0, "unit": "ppm"}, "10 ppm": {"value": 10.0, "unit": "ppm"},
    "20 ppm": {"value": 20.0, "unit": "ppm"}, "50 ppm": {"value": 50.0, "unit": "ppm"},
    ">=0.5 Da": {"value": 0.5, "unit": "Da"},
    "<=10 ppm": {"value": 10.0, "unit": "ppm"}, "0.05 Da": {"value": 0.05, "unit": "Da"},
    "0.5 Da": {"value": 0.5, "unit": "Da"}, ">=0.8 Da": {"value": 0.8, "unit": "Da"},
}


# ---------------------------------------------------------------- rules

@dataclass
class Rule:
    id: str
    kind: str              # hard | soft | imply
    heads: tuple           # heads it reads (for decoding and reporting)
    test: object           # f(facts) -> bool: True = the rule fires
    message: str
    fix: dict = field(default_factory=dict)   # imply: {head: value} to fill in


def _ppm(t):
    """Tolerance dict -> ppm at m/z 1000 (rough, for rule checks only)."""
    if not isinstance(t, dict) or "value" not in t:
        return None
    return float(t["value"]) * (1.0 if t.get("unit", "ppm") == "ppm" else 1000.0)


def _has(v, *needles):
    return v is not None and any(n.lower() in str(v).lower() for n in needles)


RULES = [
    Rule("isobaric_it_cid", "hard", ("label", "ms2_analyzer", "fragmentation", "acquisition"),
         lambda f: f.get("label") in ISOBARIC and f.get("ms2_analyzer") == "IT" and f.get("fragmentation") == "CID"
         and f.get("acquisition") != "DDA SPS-MS3",
         "isobaric reporter ions (m/z 114-135) fall below the low-mass cut-off of ion-trap CID; "
         "TMT/iTRAQ needs HCD or SPS-MS3"),
    Rule("isobaric_dia", "soft", ("label", "acquisition"),
         lambda f: f.get("label") in ISOBARIC and _has(f.get("acquisition"), "DIA"),
         "isobaric labels with DIA are unusual: co-isolated precursors mix the reporter ions"),
    Rule("sps_needs_isobaric", "soft", ("label", "acquisition"),
         lambda f: f.get("acquisition") == "DDA SPS-MS3" and f.get("label") not in ISOBARIC | {None},
         "SPS-MS3 is used for TMT/iTRAQ reporters; the label head says otherwise"),
    Rule("pasef_needs_tims", "hard", ("acquisition", "instrument"),
         lambda f: _has(f.get("acquisition"), "PASEF") and f.get("instrument") is not None
         and not _has(f.get("instrument"), *TIMS) and f.get("instrument") != "other instrument",
         "PASEF acquisition exists only on timsTOF instruments"),
    Rule("tims_is_tof", "hard", ("instrument", "ms2_analyzer"),
         lambda f: _has(f.get("instrument"), "timsTOF", "TripleTOF", "ZenoTOF", "Synapt", "Xevo")
         and f.get("ms2_analyzer") in ("ORBITRAP", "IT"),
         "this instrument records MS/MS in a TOF, not an Orbitrap or ion trap"),
    Rule("orbitrap_not_tof", "hard", ("instrument", "ms2_analyzer"),
         lambda f: _has(f.get("instrument"), "Orbitrap", "Exactive", "Exploris", "Astral") and f.get("ms2_analyzer") == "TOF",
         "this instrument has no TOF analyzer (the Astral analyzer counts as Orbitrap-class here)"),
    Rule("etd_instrument", "hard", ("fragmentation", "instrument"),
         lambda f: _has(f.get("fragmentation"), "ETD") and _has(f.get("instrument"), *NO_ETD),
         "this instrument has no ETD source"),
    Rule("immuno_unspecific", "hard", ("sample_type", "enzyme"),
         lambda f: f.get("sample_type") == "immunopeptidome" and f.get("enzyme") not in (None, "unspecific cleavage"),
         "immunopeptides are not digested: search them unspecific"),
    Rule("unspecific_needs_reason", "soft", ("sample_type", "enzyme"),
         lambda f: f.get("enzyme") == "unspecific cleavage" and f.get("sample_type") not in (None, "immunopeptidome", "body fluid"),
         "an unspecific search on a digested sample multiplies the search space"),
    Rule("diglycine_needs_trypsin", "soft", ("enrichment", "enzyme"),
         lambda f: f.get("enrichment") == "GlyGly (K)" and f.get("enzyme") not in (None, "Trypsin", "Trypsin/P"),
         "the K-GG remnant is a tryptic product; another protease leaves a different remnant"),
    Rule("fragment_tol_it", "hard", ("ms2_analyzer", "fragment_tolerance"),
         lambda f: f.get("ms2_analyzer") == "IT" and (_ppm(f.get("fragment_tolerance")) or 1e9) < 100,
         "ion-trap MS/MS cannot support a ppm-level fragment tolerance (use about 0.5 Da)"),
    Rule("fragment_tol_wide", "soft", ("ms2_analyzer", "fragment_tolerance"),
         lambda f: f.get("ms2_analyzer") in ("ORBITRAP", "TOF") and (_ppm(f.get("fragment_tolerance")) or 0) >= 300,
         "a fragment tolerance of 0.3 Da or more wastes high-resolution MS/MS (costs identifications)"),
    Rule("precursor_tol_wide", "soft", ("instrument", "precursor_tolerance"),
         lambda f: _has(f.get("instrument"), "Orbitrap", "Exactive", "Exploris", "Astral", "timsTOF")
         and (_ppm(f.get("precursor_tolerance")) or 0) >= 100 and not f.get("open_search"),
         "a precursor window of 0.1 Da or more on a high-resolution instrument only makes sense for an open search"),
    Rule("metabolic_not_supported", "soft", ("label",),
         lambda f: f.get("label") in METABOLIC,
         "15N and heavy-water labelling shift every peptide by a variable amount; ProSE cannot search them as fixed or variable mods"),
    Rule("mixed_species", "soft", ("organism", "sample_type"),
         lambda f: f.get("organism") == "mixed species" or f.get("sample_type") == "microbiome",
         "several species: search a combined database (or a metaproteomics workflow), not one model proteome"),
    Rule("unsupported_scope", "soft", ("top_down", "crosslink_search", "glycan_search"),
         lambda f: any(f.get(k) for k in ("top_down", "crosslink_search", "glycan_search")),
         "the study used a search type ProSE does not cover (top-down, cross-link or intact glycopeptide)"),
    Rule("phospho_needs_var", "imply", ("enrichment",),
         lambda f: _has(f.get("enrichment"), "Phospho"), "phospho enrichment: search Phospho (STY) as variable",
         {"variable_mods": ["Phospho (S)", "Phospho (T)", "Phospho (Y)"]}),
    Rule("pasef_implies_tof", "imply", ("acquisition",),
         lambda f: _has(f.get("acquisition"), "PASEF"), "PASEF implies TOF MS/MS and ion mobility",
         {"ms2_analyzer": "TOF", "ion_mobility": True}),
]


def check(facts: dict, kinds=("hard", "soft")) -> list:
    """Rules that fire on a set of facts (head -> value; None or absent = unknown)."""
    out = []
    for r in RULES:
        if r.kind in kinds:
            try:
                if r.test(facts):
                    out.append({"rule": r.id, "kind": r.kind, "heads": list(r.heads), "message": r.message,
                                "values": {h: facts.get(h) for h in r.heads if facts.get(h) is not None}})
            except Exception:
                continue
    return out


def facts_from_decisions(decisions: dict) -> dict:
    """Map reconciled decisions onto head names, for checking the final settings."""
    f = {}
    for k in ("organism", "enzyme", "cys_alkylation", "label", "enrichment", "acquisition", "fragmentation",
              "ms2_analyzer", "instrument", "precursor_tolerance", "fragment_tolerance"):
        d = decisions.get(k)
        if d is not None and d.value not in (None, "", "unknown"):
            f[k] = d.value
    a = str(f.get("ms2_analyzer", "")).upper()
    if a:
        f["ms2_analyzer"] = "IT" if ("TRAP" in a and "ORBI" not in a) or a in ("IT", "LIT") else \
            "TOF" if "TOF" in a else "ORBITRAP" if "ORBI" in a or "ASTRAL" in a else a
    return f


# ---------------------------------------------------------------- decoding

def decode(answers: dict, top=3, min_p=0.05):
    """Most probable joint assignment of the choice heads that breaks no hard rule.

    answers: head -> {"choice": value, "probabilities": {value: p}} with values
    already mapped to prose-auto values (NOT_STATED allowed). Returns
    (assignment, changed) where changed lists heads whose value differs from
    the per-head argmax because a hard rule forced it."""
    heads = [h for h in answers if h in HEAD_BY_NAME]
    cands = {}
    for h in heads:
        probs = answers[h].get("probabilities") or {answers[h].get("choice"): 1.0}
        ranked = sorted(probs.items(), key=lambda kv: -kv[1])
        cands[h] = [(v, p) for v, p in ranked[:top] if p >= min_p] or [(NOT_STATED, 1.0)]
    # only heads that take part in a hard rule are searched jointly; the rest keep their argmax
    linked = {h for r in RULES if r.kind == "hard" for h in r.heads}
    joint = [h for h in heads if h in linked and len(cands[h]) > 1]
    fixed = {h: cands[h][0][0] for h in heads if h not in joint}
    best, best_score = None, -math.inf
    for combo in itertools.product(*[cands[h] for h in joint]) if joint else [()]:
        a = dict(fixed)
        score = 0.0
        for h, (v, p) in zip(joint, combo):
            a[h] = v
            score += math.log(p + 1e-6)
        facts = {k: (None if v == NOT_STATED else v) for k, v in a.items()}
        score -= 20.0 * len(check(facts, kinds=("hard",)))
        if score > best_score:
            best, best_score = a, score
    changed = [h for h in joint if best[h] != cands[h][0][0]]
    return best, changed


# ---------------------------------------------------------------- text for prompts

def study_state(texts: dict, limit=60_000) -> dict:
    """One structured state from all texts of a study, labelled by source, so a
    single request reads them together (Jev 'state', or the LLM prompt)."""
    order = ["user's description", "SDRF", "PRIDE title", "PRIDE description", "PRIDE sample protocol",
             "PRIDE data protocol", "methods section", "paper methods"]
    state, used = {}, 0
    for k in order + [k for k in texts if k not in order]:
        t = (texts.get(k) or "").strip()
        if t and k not in state:
            t = t[: max(0, limit - used)]
            state[k] = t
            used += len(t)
    return state


def prompt_schema() -> str:
    """The heads, values and rules as text for a generative LLM prompt."""
    lines = []
    for h in HEADS:
        vals = "; ".join(f'"{v}"' for v, _ in h.options.values())
        lines.append(f"- {h.name}: {h.question} One of: {vals}." + (f" {h.hint}" if h.hint else ""))
    lines.append("- variable_mods (list): any of " + "; ".join(f'"{m}"' for m in VARIABLE_MOD_FLAGS) + ".")
    lines.append("- scope flags (true/false): " + "; ".join(f"{k} ({d})" for k, d in SCOPE_FLAGS.items()) + ".")
    rules = [f"- [{r.kind}] {r.message}" for r in RULES if r.kind != "imply"]
    return "HEADS:\n" + "\n".join(lines) + "\n\nCONSISTENCY RULES:\n" + "\n".join(rules)
