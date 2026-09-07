Implemented and measured in a sandbox: MS1 envelope precursor mono-correction, opt-in as `--precursor-mono {off,auto}`, on branch `claude/issue-64-precursor-mono` (commit applies on `main` `1254cb7`). Patch file, per-scan dumps, eval outputs and scripts: https://github.com/OpenMS/OpenMS/tree/claude/github-issue-64-measure-38qg2r/andes-issue-64 (transfer location; the branch here carries the same commit).

### What it does

- For each MS2, the preceding MS1 is read (same `Ms1Link` streaming as `--chimeric`), and the observed envelope at the reported charge is fitted against an **exact multi-element isotope model of a glycopeptide** (peptide averagine blended 50/50 by mass with C79H129N5O57, an average complex N-glycan) under the hypotheses "recorded precursor = M+k", k = 0..6. The fit is a cosine over M−1..M+n with a **zero-weight pre-monoisotope slot**, so the recorded hypothesis loses on a mis-picked scan (it leaves a real peak one isotope below its supposed mono); a KL on the observed distribution cannot see that.
- A shift is applied when the best k > 0 has fit ≥ 0.90, beats the recorded hypothesis by ≥ 0.15, and its mono is ≥ 3× the MS1 median. The applied shift is **k − 1** (back-off) and the search keeps the default `0..2` window, so the sweep's +1 takes the last step: an overshoot by one can never lose the true mass, an undershoot costs nothing. Nothing is widened.
- Glyco PIN gains `MonoShift`/`MonoFit`/`MonoFitGain`/`MonoSNR` only when something was fitted; `--precursor-mono-dump` writes every hypothesis' fit per MS2. `off` (default), MGF, and mzML without MS1 are byte-identical (goldens unchanged; `main`'s binary vs this one flag-off: identical on a subset run).
- Also `--glyco-index-sequon-only` (hidden): sequon-bearing peptides only in the candidate index, 3.4 GB instead of ~27 GB. Needed because the 16 GB sandbox killed the full index twice. The glyco scorer never reads a non-sequon candidate; measured on a 1,500-protein subset: same rows, peptides and labels, 16 of 7,113 rows differ in `RawScore`/`CandidateRankEntropy` only.

### Step 1 (isolation target vs trailer): settled

Probe over `thermorawfilereader` 0.7.0 on `MouseLiver-Z-T-1.raw`: the reader's precursor m/z equals **both** the isolation window target and the trailer `Monoisotopic M/Z` on all 45,905 MS2. Nothing to gain from the trailer; the +4 scans are firmware picks.

### Offline validation (dump, 23 s for the file)

True offset per pGlyco2 reference scan = the k making recorded − (peptide + Cam-C + Ox-M + glycan) an isotope multiple. Best-fitting hypothesis vs true offset (3,824 resolvable scans):

| true k | n | picks 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| −1 | 30 | 25 | 3 | 2 | 0 | 0 | 0 | 0 |
| 0 | 3,488 | 3,459 | 26 | 1 | 0 | 1 | 0 | 1 |
| 1 | 158 | 9 | 137 | 11 | 0 | 1 | 0 | 0 |
| 2 | 25 | 0 | 1 | 22 | 1 | 1 | 0 | 0 |
| 3 | 10 | 2 | 0 | 0 | 7 | 1 | 0 | 0 |
| 4 | 88 | 1 | 0 | 0 | 1 | 86 | 0 | 0 |
| 5 | 17 | 2 | 0 | 0 | 0 | 1 | 14 | 0 |
| 6 | 8 | 0 | 1 | 0 | 0 | 0 | 1 | 6 |

Scans reachable by the `0..2` sweep: 3,671 → **3,783** (+112) with **1** wrong shift (a true −1 scan). Threshold sweeps move this by ±4; back-off 0 gives +106 with 8 wrong shifts, hence back-off 1. 447 of 45,905 MS2 are shifted in the run (1:106 2:37 3:202 4:74 5:28).

### Full A/B on T-1 (one session, one binary, one database, Percolator 3.7.1 `--seed 42 -Y`)

4-thread 16 GB VM, native `.raw`, `mouse_entrap.fasta` rebuilt with the repo script, gated NeuGc default (852 compositions), both arms `--glyco-index-sequon-only`. Arm B reproduces the README quick tier within seed noise (42,108 rows both; 7,109 vs 7,122 PSMs; 3,361 vs 3,362 pGlyco2 confirmed; 2,669 = 2,669 MSFragger; 96.3% = 96.3% peptidoform; 37 vs 38 entrapment hits).

| | B: off | E: auto | acceptance |
|---|---:|---:|---|
| search wall | 9,561 s | 9,730 s (+1.8%) | — |
| glycoPSMs @1% | 7,109 | **7,225** | — |
| entrapment FDP (1:1) | 1.10% (37; CI 0.78–1.52) | **0.91%** (31; CI 0.62–1.29) | inside B's CI ✔ |
| pGlyco2 confirmed / 3,877 | 3,361 | **3,444** | ≥ 3,386 ✔ |
| gained / lost vs B | — | **+89 / −6** (all 6 lost are FDR-rejected, right answer emitted) | < 64 flips ✔ |
| MSFragger confirmed / 3,040 | 2,669 | **2,746** | — |
| peptidoform agreement pGlyco2 / MSFragger | 96.3% / 95.7% | **96.8% / 95.9%** | ≥ 96.3% ✔ |
| the 84 targets | 0 confirmed | **69 confirmed**, 5 FDR-rejected, 6 wrong target, 4 decoy | — |
| the 62 targets at +4 | 0 | **55 confirmed** + 4 FDR-rejected (59 emitted right) | ≥ 58 confirmed: **short by 3** |
| accepted by effective offset (`MonoShift` + `isotope_error`) | +0 5,332 · +1 1,347 · +2 430 | +0 5,325 · +1 1,337 · +2 402 · **+3 16 · +4 80 · +5 50 · +6 15** | no pile-up at the largest offset ✔ |
| HexNAc3 share | 5.4% at +0 | 5.5% at +0; **1.2% at +4**, 0% at +5/+6 (arm D: 79%) | matches offset-0 ✔ |
| reference coverage of +4 / +5 / +6 | — | 71 of 80 / 15 of 50 / 4 of 15 in pGlyco2; **0 entrapment hits** in +3..+6 | — |

The +4 tier is the real firmware-failure population (89% reference coverage against 42% for the run, no entrapment, 1.2% HexNAc3), not the offset-0 population shifted; the largest shift carries 15 PSMs against 613/735 in arms C/D.

**Where it falls short.** The strict "≥ 58 of 62 confirmed at +4" lands at 55: four more have the right peptidoform emitted but q just above 1%, three lose the collapse. Two levers without widening the window, not yet run (each is another ~2.7 h arm on that host): `--precursor-mono-backoff 0` (8 wrong shifts offline instead of 1), or a `0..1` window once the correction is trusted, which would also remove the (k, X) ≡ (k−1, X') degeneracy inside the sweep.

Verification: `cargo test --workspace`, `cargo clippy -D warnings` (1.87), `cargo fmt --check`, and the `thermo_raw` test against the real file, all green. Docs: `DOCS.md` flag tables and a new section in `docs/benchmarks/README.md` with the full tables and provenance.

---
_Generated by [Claude Code](https://claude.ai/code)_
