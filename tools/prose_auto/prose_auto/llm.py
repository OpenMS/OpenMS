"""LLM layer: read descriptions and methods, look up literature priors, and
review the final settings. Typed classification of the same texts with Jev
lives in jev.py.

The LLM never measures anything. Its output is evidence of source 'llm' with
modest confidence, and every extracted setting must quote the text it came
from: a quote that is not in the input is dropped (hallucination guard).
Values are normalised through the same vocabularies as text_extract.py.

Backends
--------
cli:<command>      run a CLI that reads the prompt on stdin and answers on
                   stdout, e.g. 'cli:codex exec -', 'cli:claude -p',
                   'cli:llm -m gpt-4o', 'cli:ollama run llama3.1'
openai:<model>     an OpenAI-compatible /chat/completions endpoint
                   (OPENAI_BASE_URL, OPENAI_API_KEY)
none               no LLM (the rule-based extractor still runs)
"""

from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import urllib.request

from . import schema as S
from . import text_extract as tx
from .evidence import Evidence, PARAMETERS

PRESETS = {
    "codex": "cli:codex exec --skip-git-repo-check -",  # outside a git repo codex refuses without it
    "claude": "cli:claude -p",
    "llm": "cli:llm",
}


class Backend:
    name = "none"

    def complete(self, prompt: str) -> str:
        raise NotImplementedError


class CommandBackend(Backend):
    def __init__(self, command: str, timeout=300):
        self.cmd = shlex.split(command)
        self.name = f"cli:{command}"
        self.timeout = timeout

    def complete(self, prompt):
        r = subprocess.run(self.cmd, input=prompt, capture_output=True, text=True, timeout=self.timeout)
        if r.returncode != 0:
            raise RuntimeError(f"{self.cmd[0]} failed: {r.stderr.strip()[:300]}")
        return r.stdout


class OpenAIBackend(Backend):
    def __init__(self, model):
        self.model = model
        self.name = f"openai:{model}"
        self.url = os.environ.get("OPENAI_BASE_URL", "https://api.openai.com/v1").rstrip("/") + "/chat/completions"
        self.key = os.environ.get("OPENAI_API_KEY", "")

    def complete(self, prompt):
        body = json.dumps({"model": self.model, "temperature": 0,
                           "messages": [{"role": "user", "content": prompt}]}).encode()
        req = urllib.request.Request(self.url, data=body, headers={
            "Content-Type": "application/json", "Authorization": f"Bearer {self.key}"})
        with urllib.request.urlopen(req, timeout=300) as r:
            return json.loads(r.read())["choices"][0]["message"]["content"]


def make_backend(spec: str | None) -> Backend | None:
    if not spec or spec == "none":
        return None
    spec = PRESETS.get(spec, spec)
    if spec.startswith("cli:"):
        return CommandBackend(spec[4:])
    if spec.startswith("openai:"):
        return OpenAIBackend(spec[7:])
    raise ValueError(f"unknown LLM backend {spec!r}")


def _json_block(s: str):
    m = re.search(r"```(?:json)?\s*(.*?)```", s, re.S)
    if m:
        s = m.group(1)
    start = min([i for i in (s.find("["), s.find("{")) if i >= 0], default=-1)
    if start < 0:
        raise ValueError("no JSON in LLM answer")
    return json.loads(s[start:s.rfind("]" if s[start] == "[" else "}") + 1])


_VOCAB = {
    "label": tx.LABELS, "cys_alkylation": tx.ALKYLATION, "enzyme": tx.ENZYMES,
    "enrichment": tx.ENRICHMENT, "acquisition": tx.ACQUISITION,
    "fragmentation": tx.FRAGMENTATION, "organism": tx.ORGANISMS,
}

EXTRACT_PROMPT = """You extract mass-spectrometry search settings from text for a proteomics search engine.
Return ONLY a JSON array. Each element: {{"parameter": one of {params}, "value": string or
{{"value": number, "unit": "ppm"|"Da"}} for tolerances, "quote": the exact words from the text that
support it (copied verbatim, at most 30 words), "confidence": 0..1}}.
Use these value conventions: label like "TMT10/11plex", "TMTpro (16/18plex)", "iTRAQ4plex", "SILAC",
"label free"; cys_alkylation as a Unimod name with site, e.g. "Carbamidomethyl (C)", or
"none (free Cys)"; enzyme like "Trypsin", "Lys-C"; organism as the Latin species name; acquisition
"DDA", "DIA", "PRM"; enrichment like "Phospho (STY)". Report only what the text states or clearly
implies (e.g. iodoacetamide implies Carbamidomethyl (C); HeLa implies Homo sapiens). Do not guess.

TEXT ({source}):
\"\"\"{text}\"\"\"
"""


def _norm(s):
    return re.sub(r"\s+", " ", s.lower()).strip()


def extract_settings(backend: Backend, text: str, source_label: str, user=False) -> tuple:
    """LLM extraction with quote verification. Returns (evidence, rejected)."""
    if not backend or not text.strip():
        return [], []
    prompt = EXTRACT_PROMPT.format(params=list(PARAMETERS), source=source_label, text=text[:20000])
    try:
        items = _json_block(backend.complete(prompt))
    except Exception as e:
        return [], [{"error": str(e)}]
    ev, rejected = [], []
    ntext = _norm(text)
    for it in items if isinstance(items, list) else []:
        p, v, q = it.get("parameter"), it.get("value"), str(it.get("quote", ""))
        if p not in PARAMETERS or v in (None, ""):
            rejected.append({**it, "why": "unknown parameter"})
            continue
        if not q or _norm(q) not in ntext:
            rejected.append({**it, "why": "quote not found in the text"})
            continue
        conf = float(it.get("confidence", 0.5))
        if p in _VOCAB and isinstance(v, str):
            hit = next((name for name, pat in _VOCAB[p].items() if re.search(pat, v.lower())), None)
            if hit:
                v = hit
            else:
                conf *= 0.6   # outside the vocabulary: keep, but trust less
        ev.append(Evidence(p, v, min(conf, 0.75 if user else 0.6), "llm", backend.name,
                           f"LLM read the {source_label}", quote=q))
    return ev, rejected


STUDY_PROMPT = """You set up a proteomics database search (search engine ProSE) from what the authors and
the user wrote about a study. The texts are given below by source; the user's own description and the
SDRF are the most specific, the paper's methods the most detailed.

Answer every head of the schema that the texts state or clearly imply (HeLa implies Homo sapiens,
iodoacetamide implies Carbamidomethyl (C), "TMT 16plex" implies TMTpro (16/18plex)). Leave a head out
when the texts do not say. Use exactly the listed values, except for tolerances, which are numbers.
Every answer needs a verbatim quote (at most 30 words) copied from the texts and the source it is from.
Then check the answers against the consistency rules and against each other: list every place where
the texts contradict each other or a rule (e.g. the abstract says TMT but the protocol says label free;
a TMT study measured with ion-trap CID), quoting both sides.

{schema}

Answer ONLY with JSON:
{{"settings": [{{"head": name, "value": value or {{"value": number, "unit": "ppm"|"Da"}},
                 "quote": "...", "source": "...", "confidence": 0..1}}],
  "variable_mods": [{{"value": modification, "quote": "..."}}],
  "scope": {{flag: true}},
  "contradictions": [{{"heads": [names], "issue": one sentence, "quotes": ["...", "..."]}}]}}

TEXTS:
{texts}
"""


def extract_study(backend: Backend, texts: dict, user=False) -> tuple:
    """One multi-head request over all texts of a study.

    Returns (evidence, rejected, conflicts): conflicts are the contradictions
    the LLM found between texts, plus the consistency rules that fire on its
    answers. Quotes are verified against the texts as in extract_settings."""
    state = S.study_state(texts)
    if not backend or not state:
        return [], [], []
    body = "\n\n".join(f"[{k}]\n{v}" for k, v in state.items())
    try:
        out = _json_block(backend.complete(STUDY_PROMPT.format(schema=S.prompt_schema(), texts=body)))
    except Exception as e:
        return [], [{"error": str(e)}], []
    if not isinstance(out, dict):
        return [], [{"error": "answer is not a JSON object"}], []
    ntext = _norm(body)
    ok = lambda q: bool(q) and _norm(q) in ntext
    ev, rejected, facts = [], [], {}
    cap = 0.75 if user else 0.6
    for it in out.get("settings") or []:
        h = S.HEAD_BY_NAME.get(it.get("head"))
        v, q = it.get("value"), str(it.get("quote", ""))
        if h is None or v in (None, ""):
            rejected.append({**it, "why": "unknown head"})
            continue
        if not ok(q):
            rejected.append({**it, "why": "quote not found in the texts"})
            continue
        if h.name.endswith("_tolerance"):
            if not (isinstance(v, dict) and "value" in v):
                v = S.TOL_VALUE.get(str(v))
            if not v:
                rejected.append({**it, "why": "tolerance without a number"})
                continue
            v = {"value": float(v["value"]), "unit": "Da" if str(v.get("unit", "ppm")).lower() == "da" else "ppm"}
        else:
            allowed = {str(val).lower(): val for val, _ in h.options.values()}
            if str(v).lower() not in allowed:
                rejected.append({**it, "why": f"value not in the {h.name} vocabulary"})
                continue
            v = allowed[str(v).lower()]
        facts[h.name] = v
        if not h.parameter or v in ("other species", "mixed species", "other instrument") or \
                (h.name == "enrichment" and v == "none"):
            continue
        if h.name == "acquisition":
            v = {"DDA SPS-MS3": "DDA", "ddaPASEF": "DDA"}.get(v, v)
        conf = min(float(it.get("confidence", 0.5)), cap, h.cap if h.name.endswith("_tolerance") else cap)
        ev.append(Evidence(h.parameter, v, conf, "llm", backend.name,
                           f"LLM read the {it.get('source') or 'texts'}", quote=q))
    mods = [m.get("value") for m in out.get("variable_mods") or [] if isinstance(m, dict) and ok(str(m.get("quote", "")))]
    mods = [m for m in mods if m in S.VARIABLE_MOD_FLAGS]
    if mods:
        ev.append(Evidence("variable_mods", mods, 0.5, "llm", backend.name, "LLM: searched as variable modifications"))
    facts.update({k: bool(v) for k, v in (out.get("scope") or {}).items() if k in S.SCOPE_FLAGS})
    conflicts = [{"kind": "text", "heads": c.get("heads", []), "message": c.get("issue", "")}
                 for c in out.get("contradictions") or [] if isinstance(c, dict)
                 and any(ok(str(q)) for q in c.get("quotes") or [])]
    conflicts += S.check(facts)
    return ev, rejected, conflicts


PRIOR_PROMPT = """You are a proteomics methods expert. For a {acq} run on a {instrument} with {frag}
fragmentation{extra}, what precursor and fragment mass tolerances do published methods sections and
tool defaults typically use for a database search? Answer ONLY with JSON:
{{"precursor": {{"value": number, "unit": "ppm"|"Da"}}, "fragment": {{"value": number, "unit": "ppm"|"Da"}},
"basis": "one sentence on where these numbers come from (tools, papers)"}}"""


def literature_priors(backend: Backend, instrument, acquisition, fragmentation, ms2_analyzer=""):
    if not backend:
        return []
    extra = f" and a {ms2_analyzer} MS2 analyzer" if ms2_analyzer else ""
    try:
        a = _json_block(backend.complete(PRIOR_PROMPT.format(
            acq=acquisition or "DDA", instrument=instrument or "unknown instrument",
            frag=fragmentation or "HCD", extra=extra)))
    except Exception:
        return []
    ev = []
    for key, param in (("precursor", "precursor_tolerance"), ("fragment", "fragment_tolerance")):
        t = a.get(key)
        if isinstance(t, dict) and "value" in t:
            ev.append(Evidence(param, {"value": float(t["value"]), "unit": t.get("unit", "ppm")}, 0.3,
                               "llm", backend.name, f"literature prior: {a.get('basis', '')}"))
    return ev


REVIEW_PROMPT = """You review the automatically chosen settings of a proteomics database search.
Below are the chosen settings and the evidence behind each (source: data, search, metadata,
repository, user, text, llm). Point out contradictions or settings that look implausible together
(for example TMT on an instrument whose MS2 cannot see reporter ions, a 3 ppm precursor tolerance on
an ion trap, an organism that contradicts the user's description). Do not restate settings that are
fine. Answer ONLY with a JSON array of {{"parameter": name, "issue": one sentence,
"severity": "high"|"medium"|"low"}}; an empty array if nothing is wrong.

SETTINGS:
{settings}
"""


def review(backend: Backend, decisions: dict):
    if not backend:
        return []
    compact = {k: {"value": d.value, "confidence": round(d.confidence, 2), "source": d.source,
                   "evidence": [f"{e.source}/{e.method}: {e.value} ({e.detail[:120]})" for e in d.evidence][:6]}
               for k, d in decisions.items()}
    try:
        out = _json_block(backend.complete(REVIEW_PROMPT.format(settings=json.dumps(compact, indent=1, default=str))))
        return [o for o in out if isinstance(o, dict) and o.get("parameter")]
    except Exception as e:
        return [{"parameter": "-", "issue": f"LLM review failed: {e}", "severity": "low"}]
