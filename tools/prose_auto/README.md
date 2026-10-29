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
stdin>`, `openai:<model>` (`OPENAI_BASE_URL`, `OPENAI_API_KEY`), or `djev`.

`djev` is a placeholder for a dedicated literature service: it reads
`DJEV_URL`/`DJEV_TOKEN` and, until its API is wired in, passes the same
prompts to `PROSE_AUTO_DJEV_FALLBACK` (default `claude -p`). The three
prompts (extract settings from text, typical tolerances for an
instrument/acquisition, review the final settings) are the questions such a
service would answer, so only `DjevBackend.complete()` needs to change.

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
  spectra have no single precursor). DIA tolerances come from window repeats
  and MS1 precision only.
* The diagnostic-ion Cys test is weak by design; the alkylation competition
  in the mini-search is the real test.
* Real-data validation so far covers only the small OpenMS test files: this
  container cannot reach PRIDE, UniProt or MassIVE (see the report).
