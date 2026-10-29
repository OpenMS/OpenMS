# prose-auto: an Auto Mode prototype for ProSE

`prose-auto` proposes ProSE search settings from an mzML file. It uses the
spectra and the file metadata, and optionally a FASTA or a model-organism
panel, a free-text description, a methods section, and the PRIDE/MassIVE/SDRF
records of a dataset. Every proposed setting comes with a confidence, its
source and the evidence behind it, and a setting is flagged for review
when the sources disagree.

This is a temporary exploration tool. It shows which detectors work and what
the GUI needs. It is not meant to be the final implementation: the detectors
that work belong in ProSE (C++) or in a prose-gui sidecar, as described in the
Auto Mode report.

## Install and run

```bash
pip install pyopenms numpy pyyaml          # or: pip install -e tools/prose_auto
cd tools/prose_auto
python -m prose_auto.cli panel proteomes/  # once: UniProt Swiss-Prot of 12 model organisms
python -m prose_auto.cli analyze run.mzML --panel proteomes/ \
    --description "HeLa, TMTpro 16plex, IAA, Exploris 480" \
    --accession PXD012345 --llm claude -o out/
```

The outputs, written per input to `out/`:

| File | Content |
|---|---|
| `<run>.auto.json` | everything: decisions, evidence, ProSE parameters, SDRF columns, LLM review (the GUI contract) |
| `<run>.auto.md` | readable report: settings table, conflicts, ProSE command line, evidence |
| `<run>.prose.ini` | ProSE INI fragment (`ProSE:1:Search:...`) |
| `<run>.sdrf_technical.tsv` | SDRF technical columns (instrument, label, tolerances, modifications) |

`--json` prints the report JSON to stdout and the log to stderr, so prose-gui
can run it as a sidecar.

## External tools (optional)

| Flag | Tool | What it adds |
|---|---|---|
| `--fastag <exe>` | [FASTag](https://github.com/okohlbacher/FASTag) | a database-free genus call (`-species`, reagent genera such as Bos/Sus are skipped when a real genus is close) and, with a FASTA, the Cys alkylation from the tag mass-gap histogram searched without fixed modifications (`-recon_fasta`, `-delta_out`) |
| `--fastag-subsample <f>` | FASTag | analyse only a fraction of the MS2 spectra |
| `--diaspextract <exe>` | [DIAspeXtract](https://github.com/okohlbacher/diaspextract) | converts a diaPASEF `.d` folder (required for `.d` input) or a DIA mzML with ion mobility into pseudo-MS/MS spectra (Orbitrap DIA without ion mobility skips it and stays on the DIA path; DDA-PASEF `.d` is not supported and stops with an error), so the mini-search and FASTag also run on DIA data |
| `--threads <n>` | both | threads for the external tools |

FASTag needs high-resolution MS2 to be useful: on the ion-trap BSA example it
places only a handful of tags, and its evidence stays below the mini-search.

## Pipeline

| Step | Module | What it decides | Basis |
|---|---|---|---|
| 1 metadata | `detect_instrument.py` | instrument model and class, MS1/MS2 analyzer (filter string `FTMS`/`ITMS`, mzML CV), activation, NCE, ion mobility | mzML CV terms, Thermo filter strings |
| 2 acquisition | `detect_instrument.py` | DDA / DIA / diaPASEF / PRM, DIA window scheme, charge range | isolation-window width, number of distinct window centres, repeats per cycle, charge annotation |
| 3 mass precision | `detect_tolerance.py` | fragment and precursor precision without a database | Param-Medic pairs (May et al. 2017), extended to DIA window repeats and MS1 scan-to-scan; Gaussian + uniform EM |
| 3 mass offset | `detect_tolerance.py` | systematic offset without a database | lock-mass polysiloxane ions (MS1), immonium and y1 ions (HCD MS2) |
| 4 chemistry | `detect_chemistry.py` | TMT/TMTpro/iTRAQ plex, SILAC, phospho/glyco/acetyl enrichment, Cys alkylation (weak) | reporter and diagnostic ions against control m/z values (Param-Medic 2019 idea); MS1 heavy/light pairs with an intensity check |
| 5 mini-search | `tagsearch.py`, `hypotheses.py` | organism, Cys alkylation, protease, calibrated offsets and tolerances | binomial b/y scoring with target/reversed-decoy FDR on the top-TIC subset; competing hypotheses compared by PSMs at 1% FDR |
| 6 records and text | `repository.py`, `text_extract.py` | what the submitters and the user say | PRIDE API v3, PROXI, bigbio/sdrf-annotated-datasets, Europe PMC methods; keyword rules with verbatim quotes |
| 6 LLM | `llm.py` | the same from free text, literature priors, a final review | any CLI or OpenAI-compatible endpoint; quotes must occur in the text |
| 7 reconcile | `reconcile.py` | one decision per parameter, confidence, conflicts, status | measured beats stated for instrument-determined settings; a clear mini-search result beats text for sample prep |
| (opt) ProSE | `prose_runner.py` | confirmation with ProSE's own calibration pass and an open search | `-summary_out` YAML (`precursor_error`, `fragment_error`), `-out_mod_analysis_dir` |

## LLM backends

`--llm` takes `none` (the default; the rule-based extractor still runs),
`claude` (`claude -p`), `codex` (`codex exec -`), `cli:<any command reading
stdin>`, or `openai:<model>` (`OPENAI_BASE_URL`, `OPENAI_API_KEY`). It
answers three prompts: extract settings from text, typical tolerances for an
instrument/acquisition (`--llm-priors`), and a review of the final settings
(`--llm-review`).

## Jev (typed classifier)

`--jev` sends the same texts to [Jev](https://jevtypesafeai.com/jev/classifier)
(`POST /api/v1/decide`, key in `JEV_API_KEY`, optional `JEV_URL` and
`JEV_MODEL`). Jev picks one label from a closed set and returns a probability
for each, so its answers stay inside the prose-auto vocabularies and need no
quote guard. `jev.py` asks, in one request per text, for the organism,
protease, Cys alkylation, label, enrichment, acquisition, fragmentation and
instrument, each with a `not_stated` option; answers with p ≥ 0.5 become
`llm` evidence with confidence 0.7·p (≤ 0.6 for methods text, ≤ 0.7 for the
user's description). After reconciling, it asks per setting whether each text
supports, contradicts or does not mention it, and lists contradictions under
"Conflicts and notes". Tolerances are numbers, not labels, so they stay with
the regex extractor and `--llm`. Cost is about $0.001 per request.

The LLM never overrides a measurement. Its findings are evidence of source
`llm` with confidence ≤ 0.75. Any extracted setting whose quote does not
occur verbatim in the input is rejected and listed in the report.

## Tests

```bash
cd tools/prose_auto && python -m pytest tests -q
```

The tests use OpenMS test data: BSA1 (LTQ Orbitrap XL, CID/ion trap), the
DIAuditor and OpenSWATH DIA inputs, and the IsobaricAnalyzer iTRAQ input. A
fake LLM covers the quote guard.

## Known limits (prototype)

* The mini-search is pure Python. On the BSA example (1,120 MS2 spectra, four
  databases with 6.6k to 870k peptides) the whole analysis takes about 20 s;
  a full human proteome panel needs ProSE's fragment index instead.
* The mini-search has no variable modifications and covers DDA only (DIA
  spectra have no single precursor) unless `--diaspextract` converts the DIA
  run to pseudo-spectra first. Without it, DIA tolerances come from window
  repeats and MS1 precision only.
* The diagnostic-ion Cys test is weak by design; the alkylation competition
  in the mini-search is the real test.
* Real-data validation so far covers only the small OpenMS test files: this
  container cannot reach PRIDE, UniProt or MassIVE (see the report).
