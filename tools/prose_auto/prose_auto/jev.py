"""Jev (TypeSafe AI) as a typed, multi-head classifier for the text evidence.

Jev does not generate text: it picks one label out of a closed set (up to 255)
and returns a probability for every label, so its answers cannot leave the
prose-auto vocabularies. prose-auto uses it in three ways:

* read_study: all texts of a study (the user's description, SDRF, PRIDE
  title/description/protocols, methods section) go into one structured
  `state`, and one /v1/decide request asks every head of schema.HEADS as a
  `choice` question (each with a "not stated" option), plus a `noul` (yes/no)
  question per variable modification and per out-of-scope search type. The
  head answers are decoded jointly (schema.decode): the most probable
  combination that breaks no hard consistency rule wins, so TMT + ion-trap
  CID or diaPASEF + Orbitrap cannot come out of the text together. Soft rules
  and decoding changes are reported as text conflicts.
* cross_check: for every final decision, ask whether the texts support it,
  contradict it or do not mention it.
* choose: pick the next step of the tuning loop (tune.py) from a closed list,
  given a description of the last search's diagnostics.

Numbers stay with the regex extractor and the generative LLM; Jev answers the
tolerance heads only as coarse buckets, with low confidence.

API: two routes with the same body and answers. A TypeSafe console key goes to
the official endpoint POST https://api.typesafe.ai/v1/systemone; a jv_live_
key from jevtypesafeai.com goes to its hosted proxy POST
https://jevtypesafeai.com/api/v1/decide. Both take 'Authorization: Bearer <key>'
(JEV_API_KEY or TYPESAFE_API_KEY; JEV_URL overrides the endpoint) and a body
{model, state, questions}; choice questions are {"type": "choice", "instructions": ...,
"criteria": {label: description}} and answer {"choice", "confidence",
"probabilities"}; noul questions answer {"noul": p(yes)}.
"""

from __future__ import annotations

import json
import os
import urllib.request

from . import schema as S
from .evidence import Evidence, PARAMETERS

URL_OFFICIAL = "https://api.typesafe.ai/v1/systemone"
URL_HOSTED = "https://jevtypesafeai.com/api/v1/decide"
USD_PER_INPUT_TOKEN = 0.42e-6      # list price; output tokens are free
NOT_STATED = S.NOT_STATED

READ_INSTRUCTIONS = ("You read the description of a proteomics experiment. The state holds the texts by "
                     "source; the user's own description and the SDRF are the most specific. ")


class Jev:
    def __init__(self, key=None, url=None, model=None, transport=None, timeout=60):
        self.key = (key or os.environ.get("JEV_API_KEY") or os.environ.get("TYPESAFE_API_KEY") or "").strip()
        self.url = url or os.environ.get("JEV_URL") or (URL_HOSTED if self.key.startswith("jv_") else URL_OFFICIAL)
        self.model = model or os.environ.get("JEV_MODEL", "jev-latest")
        self.timeout = timeout
        self.transport = transport or self._http   # tests inject a fake
        self.name = f"jev:{self.model}"
        self.cost_usd = 0.0
        self.calls = 0

    def _http(self, body):
        headers = {"Content-Type": "application/json", "Authorization": f"Bearer {self.key}"}
        if len(json.dumps(body["state"])) > 400_000:
            headers["TypeSafe-Danger-Zone-Long-Context"] = "true"
        req = urllib.request.Request(self.url, data=json.dumps(body).encode(), headers=headers)
        with urllib.request.urlopen(req, timeout=self.timeout) as r:
            return json.loads(r.read())

    def decide(self, state, questions):
        out = self.transport({"model": self.model, "state": state, "questions": questions})
        self.calls += 1
        u = out.get("usage") or {}
        self.cost_usd += float(u.get("cost_usd") or (u.get("input_tokens") or 0) * USD_PER_INPUT_TOKEN)
        return out.get("answers", {})


def _choice(instructions, criteria):
    return {"type": "choice", "instructions": instructions, "criteria": criteria}


def _head_question(h: S.Head):
    criteria = {k: d for k, (_, d) in h.options.items()}
    criteria[NOT_STATED] = "the texts do not say"
    return _choice(READ_INSTRUCTIONS + h.question + (" " + h.hint if h.hint else "")
                   + " Answer not_stated unless the texts say so or clearly imply it.", criteria)


def study_questions():
    qs = {h.name: _head_question(h) for h in S.HEADS}
    for mod, desc in S.VARIABLE_MOD_FLAGS.items():
        qs["var:" + mod] = {"type": "noul", "instructions": READ_INSTRUCTIONS +
                            f"Did the authors search {desc} as a variable modification?"}
    for flag, desc in S.SCOPE_FLAGS.items():
        qs["scope:" + flag] = {"type": "noul", "instructions": READ_INSTRUCTIONS + f"Did the study use {desc}?"}
    return qs


def _map_answer(h: S.Head, a: dict):
    """Jev answer with option keys -> {choice, probabilities} with prose-auto values."""
    probs = {}
    for k, p in (a.get("probabilities") or {}).items():
        v = NOT_STATED if k == NOT_STATED else h.options.get(k, (None,))[0]
        if v is not None:
            probs[v] = probs.get(v, 0.0) + float(p)
    key = a.get("choice")
    choice = NOT_STATED if key == NOT_STATED else h.options.get(key, (None,))[0]
    if choice is not None and choice not in probs:
        probs[choice] = float(a.get("confidence", 1.0))
    return {"choice": choice, "probabilities": probs}


def read_study(jev: Jev, texts: dict, user=False, min_p=0.5):
    """All texts of a study, all heads, one request.

    Returns (evidence, report) where report holds the decoded heads, the
    per-head probabilities, flags, rules that fire on the text alone, and the
    heads the joint decoding changed."""
    state = S.study_state(texts)
    if not jev or not state:
        return [], {}
    try:
        ans = jev.decide(state, study_questions())
    except Exception as e:
        return [], {"error": str(e)}
    mapped = {h.name: _map_answer(h, ans[h.name]) for h in S.HEADS if h.name in ans}
    decoded, changed = S.decode(mapped)
    flags = {k.split(":", 1)[1]: float(a.get("noul", 0)) for k, a in ans.items() if k.startswith(("var:", "scope:"))}
    facts = {h: (None if v == NOT_STATED else v) for h, v in decoded.items()}
    facts.update({k: p >= 0.5 for k, p in flags.items() if k in S.SCOPE_FLAGS})
    ev, cap_user = [], 0.7 if "user's description" in state else 0.6
    for h in S.HEADS:
        v = decoded.get(h.name)
        if not h.parameter or v in (None, NOT_STATED) or h.parameter not in PARAMETERS:
            continue
        p = mapped[h.name]["probabilities"].get(v, 0.0)
        if p < min_p and h.name not in changed:
            continue
        if h.name == "enrichment" and v == "none":
            continue
        if h.name in ("precursor_tolerance", "fragment_tolerance"):
            v = S.TOL_VALUE.get(v)
            if not v:
                continue
        if h.name == "acquisition":
            v = {"DDA SPS-MS3": "DDA", "ddaPASEF": "DDA"}.get(v, v)
        if h.name == "organism" and v in ("other species", "mixed species"):
            continue
        conf = min(h.cap, cap_user, 0.7 * p) * (0.8 if h.name in changed else 1.0)
        ev.append(Evidence(h.parameter, v, conf, "llm", jev.name,
                           f"Jev read the study texts: {h.name} = {decoded[h.name]} (p={p:.2f})"
                           + ("; chosen by the consistency rules over the top answer" if h.name in changed else "")))
    mods = [m for m in S.VARIABLE_MOD_FLAGS if flags.get(m, 0) >= 0.6]
    if mods:
        ev.append(Evidence("variable_mods", mods, 0.5, "llm", jev.name,
                           "Jev: the texts say these were searched as variable modifications"))
    for r in S.RULES:
        if r.kind == "imply" and r.test(facts):
            for k, val in r.fix.items():
                if k in PARAMETERS and k != "ion_mobility":
                    ev.append(Evidence(k, val, 0.45, "llm", jev.name, f"implied by the text: {r.message}"))
    report = {"heads": decoded, "probabilities": {k: v["probabilities"] for k, v in mapped.items()},
              "changed_by_rules": changed, "flags": flags, "text_rules": S.check(facts),
              "sources": list(state)}
    return ev, report


def classify(jev: Jev, text: str, source_label: str, user=False, min_p=0.5):
    """One text on its own (kept for single-source use)."""
    return read_study(jev, {source_label: text}, user=user, min_p=min_p)


def cross_check(jev: Jev, texts, decisions: dict, source_label="description"):
    """Ask, for every decided setting, whether the texts support or contradict it.
    texts is a str or a {source: text} dict. Returns review issues."""
    state = S.study_state(texts if isinstance(texts, dict) else {source_label: texts})
    if not jev or not state:
        return []
    qs = {}
    for p, d in decisions.items():
        if d.value in (None, "", "unknown") or p in ("charge_range", "dia_windows", "collision_energy"):
            continue
        qs[p] = _choice(READ_INSTRUCTIONS + f"The search will use {PARAMETERS[p]} = {d.value}. Do the texts agree?",
                        {"supports": "the texts state this or something equivalent",
                         "contradicts": "the texts state a different value",
                         "not_mentioned": "the texts say nothing about it"})
    try:
        ans = jev.decide(state, qs)
    except Exception as e:
        return [{"parameter": "-", "issue": f"Jev cross-check failed: {e}", "severity": "low"}]
    out = []
    for p, a in ans.items():
        prob = float((a.get("probabilities") or {}).get("contradicts", 0))
        if a.get("choice") == "contradicts" and prob >= 0.5:
            out.append({"parameter": p, "severity": "high" if prob >= 0.8 else "medium",
                        "issue": f"Jev: the texts contradict {decisions[p].value} (p={prob:.2f})"})
    return out


def consistency(decisions: dict):
    """Consistency rules on the final settings (no Jev call): review issues."""
    return [{"parameter": ",".join(r["heads"]), "severity": "high" if r["kind"] == "hard" else "low",
             "issue": f"{r['message']} ({', '.join(f'{k}={v}' for k, v in r['values'].items())})"}
            for r in S.check(S.facts_from_decisions(decisions))]


def choose(jev: Jev, situation: str, actions: dict, instructions: str):
    """Pick one action from a closed set. Returns (action, probabilities)."""
    a = jev.decide({"situation": situation}, {"next": _choice(instructions, actions)}).get("next", {})
    return a.get("choice"), a.get("probabilities") or {}
