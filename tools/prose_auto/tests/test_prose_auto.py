"""Tests on OpenMS test data (run from tools/prose_auto: python -m pytest tests)."""

import json
import os

import pytest

from prose_auto import detect_chemistry as chem
from prose_auto import detect_instrument as inst
from prose_auto import detect_tolerance as tol
from prose_auto import llm as L
from prose_auto import prose_params as pp
from prose_auto import tagsearch as ts
from prose_auto import text_extract as tx
from prose_auto.evidence import Evidence
from prose_auto.hypotheses import HypothesisTester
from prose_auto.reconcile import reconcile
from prose_auto.repository import sdrf_evidence
from prose_auto.spectra import load_mzml

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
TOPP = os.path.join(ROOT, "src", "tests", "topp")
BSA = os.path.join(ROOT, "share", "OpenMS", "examples", "BSA")


def first(evs, p):
    c = [e for e in evs if e.parameter == p]
    return max(c, key=lambda e: e.confidence) if c else None


@pytest.fixture(scope="module")
def bsa():
    return load_mzml(os.path.join(BSA, "BSA1.mzML"))


def test_instrument_and_dda(bsa):
    ev = inst.detect_instrument(bsa) + inst.detect_acquisition(bsa)
    assert first(ev, "instrument").value == "LTQ Orbitrap XL"
    assert first(ev, "ms1_analyzer").value == "ORBITRAP"
    assert first(ev, "ms2_analyzer").value == "IT"
    assert first(ev, "fragmentation").value == "CID"
    assert first(ev, "acquisition").value == "DDA"


def test_dia_windows():
    run = load_mzml(os.path.join(TOPP, "DIAuditor_1_input.mzML"))
    ev = inst.detect_acquisition(run)
    assert first(ev, "acquisition").value == "DIA"
    w = first(ev, "dia_windows").value
    assert w["n_windows"] >= 4 and w["min_mz"] == 400.0


def test_itraq_reporters():
    run = load_mzml(os.path.join(TOPP, "IsobaricAnalyzer_input_1.mzML"))
    lab = [e for e in chem.detect_labels(run, True) if e.method == "reporter ions"][0]
    assert lab.value == "iTRAQ4plex"


def test_low_res_tolerances(bsa):
    ev = tol.detect_tolerances(bsa, "DDA", ms2_high_res=False, ms1_high_res=True)
    f = first(ev, "fragment_tolerance")
    assert f.value["unit"] == "Da" and 0.2 <= f.value["value"] <= 0.8
    p = first(ev, "precursor_tolerance")
    assert p.value["unit"] == "ppm" and p.value["value"] <= 10


def test_gauss_uniform_fit():
    import numpy as np
    rng = np.random.default_rng(1)
    x = np.concatenate([rng.normal(2.0, 1.5, 2000), rng.uniform(-30, 30, 1000)])
    mu, sigma, pi = tol.fit_gauss_uniform(x, 30)
    assert abs(mu - 2.0) < 0.2 and abs(sigma - 1.5) < 0.2 and 0.6 < pi < 0.75


def _write(path, recs):
    with open(path, "w") as fh:
        for h, s in recs:
            fh.write(f">{h}\n{s}\n")


def test_mini_search_alkylation(bsa, tmp_path):
    # BSA plus an unrelated proteome-like decoy organism (shuffled BSA-free mix)
    mix = os.path.join(ROOT, "share", "OpenMS", "examples", "TOPPAS", "data", "BSA_Identification",
                       "18Protein_SoCe_Tr_detergents_trace.fasta")
    recs = [(a, s) for a, s in ts.read_fasta(mix) if "SORC5" not in a]
    ecoli = os.path.join(ROOT, "share", "OpenMS", "examples", "TOPPAS", "data", "Identification",
                         "target_decoy_Ecoli_K12_TaxID_83333.proteomes.fasta")
    _write(tmp_path / "mix.fasta", recs)
    h = HypothesisTester(bsa, 0.4, False, prec_ppm=20, n_spectra=1200)
    ev = h.organism({"mix": str(tmp_path / "mix.fasta"), "ecoli": ecoli})
    assert first(ev, "organism").value == "mix"
    alk = h.alkylation()[0]
    assert alk.value == "Carbamidomethyl (C)" and alk.confidence > 0.8
    acc = h.accuracy()
    assert abs(first(acc, "precursor_offset").value["value"]) < 3


def test_text_extract():
    ev = tx.extract("HeLa cells, alkylated with iodoacetamide, TMT10plex, Q Exactive HF, DIA", "user")
    vals = {e.parameter: e.value for e in ev}
    assert vals["organism"] == "Homo sapiens"
    assert vals["cys_alkylation"] == "Carbamidomethyl (C)"
    assert vals["label"] == "TMT10/11plex"
    assert vals["instrument"] == "Q Exactive HF"
    assert vals["acquisition"] == "DIA"


class FakeLLM(L.Backend):
    name = "fake"

    def __init__(self, answer):
        self.answer = answer

    def complete(self, prompt):
        return self.answer


def test_llm_quote_guard():
    text = "Proteins were alkylated with iodoacetamide and digested with trypsin."
    ans = json.dumps([
        {"parameter": "cys_alkylation", "value": "carbamidomethyl", "quote": "alkylated with iodoacetamide", "confidence": 0.9},
        {"parameter": "label", "value": "TMT", "quote": "labelled with TMT", "confidence": 0.9},   # not in text
        {"parameter": "colour", "value": "blue", "quote": "trypsin", "confidence": 0.9},
    ])
    ev, rej = L.extract_settings(FakeLLM("Here you go:\n```json\n" + ans + "\n```"), text, "methods")
    assert [(e.parameter, e.value) for e in ev] == [("cys_alkylation", "Carbamidomethyl (C)")]
    assert len(rej) == 2


def test_sdrf_parse():
    sdrf = ("source name\tcharacteristics[organism]\tcomment[label]\tcomment[instrument]\t"
            "comment[modification parameters]\tcomment[cleavage agent details]\tcomment[precursor mass tolerance]\n"
            "s1\tHomo sapiens\tTMT126\tNT=Q Exactive HF;AC=MS:1002523\tNT=Carbamidomethyl;TA=C;MT=fixed\t"
            "NT=Trypsin;AC=MS:1001251\t10 ppm\n")
    vals = {e.parameter: e.value for e in sdrf_evidence(sdrf)}
    assert vals["organism"] == "Homo sapiens" and vals["cys_alkylation"] == "Carbamidomethyl (C)"
    assert vals["precursor_tolerance"] == {"value": 10.0, "unit": "ppm"}


def test_reconcile_conflict_and_params():
    ev = [
        Evidence("label", "TMT10/11plex", 0.95, "data", "reporter ions"),
        Evidence("label", "label free", 0.75, "user", "keywords"),
        Evidence("precursor_tolerance", {"value": 4, "unit": "ppm"}, 0.9, "search", "mini-search PSMs"),
        Evidence("fragment_tolerance", {"value": 0.0, "unit": "ppm"}, 0.0, "data", "x"),
        Evidence("cys_alkylation", "Carbamidomethyl (C)", 0.9, "search", "alkylation competition"),
        Evidence("cys_alkylation", "none (free Cys)", 0.75, "user", "keywords"),
    ]
    d = reconcile(ev)
    assert d["label"].value == "TMT10/11plex" and d["label"].status == "review" and d["label"].conflicts
    assert d["cys_alkylation"].value == "Carbamidomethyl (C)" and d["cys_alkylation"].status == "review"
    P = pp.to_prose(d)["params"]
    assert "TMT6plex (K)" in P["Search:modifications:fixed"]
    assert P["Search:precursor:mass_tolerance_lower"] == 4
    assert "<ITEMLIST" in pp.to_ini_snippet(P)
