"""Jev (TypeSafe AI) as a typed classifier for the text evidence.

Jev does not generate text: it picks one label out of a closed set (up to 255)
and returns a probability for every label, so its answers cannot leave the
prose-auto vocabularies. That makes it a good fit for two jobs:

* classify: read a description, methods section or PRIDE protocol and pick
  the stated organism, protease, Cys alkylation, label, enrichment,
  acquisition, fragmentation and instrument, each with a "not stated" option.
  All questions go out in one /v1/decide request.
* cross-check: for every final decision, ask whether the text supports it,
  contradicts it or does not mention it. A contradiction is reported like an
  LLM review issue.

Numbers (tolerances) are not classification questions; they stay with the
regex extractor and the generative LLM. Jev evidence has source 'llm' and
confidence <= 0.7, like the LLM layer.

API: POST {JEV_URL or https://jevtypesafeai.com/api/v1/decide}, header
'Authorization: Bearer $JEV_API_KEY', body {model, state, questions};
questions are {"type": "choice", "instructions": ..., "criteria": {label: description}}
and answers {"choice", "confidence", "probabilities"}.
"""

from __future__ import annotations

import json
import os
import re
import urllib.request

from . import text_extract as tx
from .evidence import Evidence, PARAMETERS

URL = "https://jevtypesafeai.com/api/v1/decide"
NOT_STATED = "not_stated"


def _slug(s):
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")[:60]


def _vocab(table, describe=None):
    return {_slug(v): (v, describe(v) if describe else v) for v in table}


QUESTIONS = {
    "organism": ("Which organism were the analysed proteins from (the species to search against)?",
                 _vocab(tx.ORGANISMS)),
    "enzyme": ("Which protease was used for digestion?", _vocab(tx.ENZYMES)),
    "cys_alkylation": ("Which cysteine alkylation (reagent or resulting modification) was used?",
                       _vocab(tx.ALKYLATION, lambda v: {
                           "Carbamidomethyl (C)": "iodoacetamide or chloroacetamide, carbamidomethyl",
                           "Propionamide (C)": "acrylamide, propionamide",
                           "Nethylmaleimide (C)": "N-ethylmaleimide (NEM)",
                           "Methylthio (C)": "MMTS, methyl methanethiosulfonate",
                           "Carboxymethyl (C)": "iodoacetic acid, carboxymethyl",
                           "none (free Cys)": "cysteines were explicitly not alkylated"}.get(v, v))),
    "label": ("Which quantification labelling was used?", _vocab(tx.LABELS)),
    "enrichment": ("Which PTM enrichment was performed, if any?",
                   {**_vocab(tx.ENRICHMENT), "none": ("none", "no PTM enrichment, whole proteome")}),
    "acquisition": ("Which MS acquisition mode was used?", _vocab(tx.ACQUISITION)),
    "fragmentation": ("Which fragmentation method was used for MS/MS?", _vocab(tx.FRAGMENTATION)),
    "instrument": ("Which mass spectrometer model was used?", _vocab(tx.INSTRUMENTS)),
}


class Jev:
    def __init__(self, key=None, url=None, model=None, transport=None, timeout=60):
        self.key = key or os.environ.get("JEV_API_KEY", "")
        self.url = url or os.environ.get("JEV_URL", URL)
        self.model = model or os.environ.get("JEV_MODEL", "jev-latest")
        self.timeout = timeout
        self.transport = transport or self._http   # tests inject a fake
        self.name = f"jev:{self.model}"

    def _http(self, body):
        headers = {"Content-Type": "application/json", "Authorization": f"Bearer {self.key}"}
        if len(json.dumps(body["state"])) > 400_000:
            headers["TypeSafe-Danger-Zone-Long-Context"] = "true"
        req = urllib.request.Request(self.url, data=json.dumps(body).encode(), headers=headers)
        with urllib.request.urlopen(req, timeout=self.timeout) as r:
            return json.loads(r.read())

    def decide(self, state, questions):
        return self.transport({"model": self.model, "state": state, "questions": questions}).get("answers", {})


def _choice(instructions, vocab):
    criteria = {k: d for k, (_, d) in vocab.items()}
    criteria[NOT_STATED] = "the text does not say"
    return {"type": "choice", "instructions": instructions + " Answer not_stated unless the text says so.",
            "criteria": criteria}


def classify(jev: Jev, text: str, source_label: str, user=False, min_p=0.5):
    """Evidence for the categorical parameters. Returns (evidence, raw answers)."""
    if not jev or not text.strip():
        return [], {}
    qs = {p: _choice(instr, vocab) for p, (instr, vocab) in QUESTIONS.items()}
    try:
        ans = jev.decide({"source": source_label, "text": text[:200_000]}, qs)
    except Exception as e:
        return [], {"error": str(e)}
    ev, cap = [], 0.7 if user else 0.6
    for p, a in ans.items():
        key = a.get("choice")
        if p not in QUESTIONS or key in (None, NOT_STATED):
            continue
        prob = float((a.get("probabilities") or {}).get(key, a.get("confidence", 0)))
        if prob < min_p or key not in QUESTIONS[p][1]:
            continue
        value = QUESTIONS[p][1][key][0]
        if p == "enrichment" and value == "none":
            continue
        ev.append(Evidence(p, value, min(cap, 0.7 * prob), "llm", jev.name,
                           f"Jev classified the {source_label}: {value} (p={prob:.2f})"))
    return ev, ans


def cross_check(jev: Jev, text: str, decisions: dict, source_label="description"):
    """Ask, for every decided setting, whether the text supports or contradicts it.
    Returns review issues in the LLM review format."""
    if not jev or not text.strip():
        return []
    qs = {}
    for p, d in decisions.items():
        if d.value in (None, "", "unknown") or p in ("charge_range", "dia_windows", "collision_energy"):
            continue
        qs[p] = {"type": "choice",
                 "instructions": f"The search will use {PARAMETERS[p]} = {d.value}. Does the text agree?",
                 "criteria": {"supports": "the text states this or something equivalent",
                              "contradicts": "the text states a different value",
                              "not_mentioned": "the text says nothing about it"}}
    try:
        ans = jev.decide({"source": source_label, "text": text[:200_000]}, qs)
    except Exception as e:
        return [{"parameter": "-", "issue": f"Jev cross-check failed: {e}", "severity": "low"}]
    out = []
    for p, a in ans.items():
        prob = float((a.get("probabilities") or {}).get("contradicts", 0))
        if a.get("choice") == "contradicts" and prob >= 0.5:
            out.append({"parameter": p, "severity": "high" if prob >= 0.8 else "medium",
                        "issue": f"Jev: the {source_label} contradicts {decisions[p].value} (p={prob:.2f})"})
    return out
