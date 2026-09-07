# bigbio/andes#64 — precursor mono-correction from the MS1 envelope: implementation and measurement

Branch `claude/issue-64-precursor-mono`, commit `bc50a35` on top of `main` `1254cb7`
(patch: `0001-glyco-correct-precursors-to-the-monoisotopic-peak-fr.patch`, 19 files,
+1,424/−27). Not pushed: this session has no write access to bigbio/andes (the git proxy
refuses to inject a credential; `add_repo` was denied). Apply with `git am`.

## What was built

* `--precursor-mono {off,auto}` (glyco; mzML or Thermo `.raw`). For each MS2 the preceding
  MS1 is read, and the observed envelope at the reported charge is fitted against an exact
  multi-element isotope model of a glycopeptide (peptide averagine blended 50/50 by mass
  with an average complex N-glycan, C79H129N5O57) under the hypotheses "recorded = M+k",
  k = 0..6. The fit is a cosine over M−1..M+n with a zero-weight pre-monoisotope slot, so a
  hypothesis that leaves a real peak one isotope below its monoisotope (what the recorded
  hypothesis does on a mis-picked scan) loses. A shift is applied when the best k > 0 has
  fit ≥ 0.90, beats the recorded hypothesis by ≥ 0.15, and its monoisotope is ≥ 3× the
  MS1 median intensity; the applied shift is **k − 1** (back-off), and the search keeps
  its default `0..2` window so the sweep's +1 takes the last step (an overshoot by one can
  never lose the true mass; an undershoot costs nothing). Nothing is widened.
* Glyco PIN columns `MonoShift`, `MonoFit`, `MonoFitGain`, `MonoSNR`, present only when the
  run fitted at least one MS2; `--precursor-mono-dump` writes every hypothesis' fit per MS2.
  Hidden tuning flags for every threshold.
* `--glyco-index-sequon-only` (hidden): enumerate only N-X-S/T-bearing peptides into the
  candidate index. Needed here because the full mouse-entrapment index (~27 GB) cannot fit
  the 16 GB sandbox (killed twice at the cgroup limit). The glyco scorer never reads a
  non-sequon candidate; measured on a 1,500-protein subset the PIN has the same rows,
  peptides and labels, with 16 of 7,113 rows differing in `RawScore`/`CandidateRankEntropy`.
* Byte-identity: `main`'s binary and this binary with the flag off produce identical
  `.glyco.pin` on the subset run; both goldens pass unchanged; `--precursor-mono auto` on
  MGF or on an mzML without MS1 warns and is byte-identical to off (new tests).
* Verification: `cargo test --workspace` (all 71 test binaries), `cargo clippy -D warnings`
  (pinned 1.87), `cargo fmt --check`, the `thermo_raw` integration test against the real
  file, all green.

## Step 1 of the issue (isolation target vs trailer)

Probe over `thermorawfilereader` 0.7.0 on `MouseLiver-Z-T-1.raw`: the reader's precursor
m/z equals BOTH the isolation window target and the trailer `Monoisotopic M/Z` on all
45,905 MS2 (`thermo_trailer_probe.tsv`). The trailer holds nothing the reader is not
already using; the +4 scans are firmware picks.

## Offline validation of the corrector (`precursor_mono_dump.tsv`, 23 s for the file)

True offset per pGlyco2 reference scan = the integer k making
recorded − (peptide + Cam-C + Ox-M + glycan) an isotope multiple (53 of 3,877 unresolved).

| true offset | n | fit picks 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| −1 | 30 | 25 | 3 | 2 | 0 | 0 | 0 | 0 |
| 0 | 3,488 | 3,459 | 26 | 1 | 0 | 1 | 0 | 1 |
| 1 | 158 | 9 | 137 | 11 | 0 | 1 | 0 | 0 |
| 2 | 25 | 0 | 1 | 22 | 1 | 1 | 0 | 0 |
| 3 | 10 | 2 | 0 | 0 | 7 | 1 | 0 | 0 |
| 4 | 88 | 1 | 0 | 0 | 1 | 86 | 0 | 0 |
| 5 | 17 | 2 | 0 | 0 | 0 | 1 | 14 | 0 |
| 6 | 8 | 0 | 1 | 0 | 0 | 0 | 1 | 6 |

Scans reachable by the `0..2` sweep: 3,671 without correction → **3,783 with** (+112),
**1 wrong shift** (a true −1 scan). Threshold sweep moved this by ±4; back-off 0 gives
+106 with 8 wrong shifts (why back-off 1 is the default). 447 of 45,905 MS2 are shifted
in the whole run (1:106, 2:37, 3:202, 4:74, 5:28 after back-off).

## Full A/B (one session, one binary, one database, Percolator 3.7.1 `--seed 42 -Y`)

4-thread 16 GB VM, native `.raw` (sha256 `2f0142b7…`), `mouse_entrap.fasta` rebuilt from
UniProt with the repo script (34,554 sequences), gated NeuGc default (852 compositions),
both arms `--glyco-index-sequon-only`. Arm B reproduces the README quick tier within
Percolator seed noise (42,108 rows both; 7,109 vs 7,122 PSMs; 3,361 vs 3,362 pGlyco2
confirmed; 2,669 = 2,669 MSFragger; 96.3% = 96.3% peptidoform; 37 vs 38 entrapment hits).

| | B: off | E: auto | acceptance (#64) |
|---|---:|---:|---|
| search wall | 9,561 s | 9,730 s (+1.8%) | — |
| glycoPSMs @1% | 7,109 | **7,225** | — |
| entrapment FDP (1:1) | 1.10% (37; CI 0.78–1.52) | **0.91%** (31; CI 0.62–1.29) | inside B's CI ✔ |
| pGlyco2 confirmed / 3,877 | 3,361 | **3,444** | ≥ 3,386 ✔ |
| gained / lost vs B | — | **+89 / −6** (all 6 FDR-rejected, right answer emitted) | < 64 flips ✔ |
| MSFragger confirmed / 3,040 | 2,669 | **2,746** | — |
| peptidoform agreement pGlyco2 / MSFragger | 96.3% / 95.7% | **96.8% / 95.9%** | ≥ 96.3% ✔ |
| the 84 targets | 0 confirmed | **69 confirmed**, 5 FDR-rejected, 6 wrong target, 4 decoy | — |
| the 62 targets at +4 | 0 | **55 confirmed** + 4 FDR-rejected | ≥ 58 confirmed: **short by 3** |
| accepted by effective offset (shift + isotope_error) | +0 5,332 · +1 1,347 · +2 430 | +0 5,325 · +1 1,337 · +2 402 · +3 16 · +4 80 · +5 50 · +6 15 | no pile-up at the largest offset ✔ |
| HexNAc3 share | 5.4% at +0 | 5.5% at +0; **1.2% at +4**, 0% at +5/+6 (arm D: 79%) | matches offset-0 ✔ |
| reference coverage of +4 / +5 / +6 | — | 71 of 80 / 15 of 50 / 4 of 15 in pGlyco2; **0 entrapment** | — |
| goldens with flag off | unchanged | unchanged | ✔ |

The +4 tier is the real firmware-failure population (89% reference coverage against 42%
for the run, no entrapment hit, 1.2% HexNAc3), not the offset-0 population shifted. There
is no pile-up: the largest shift carries 15 PSMs against 613 and 735 in arms C and D.

Short of the strict "≥ 58 of 62 confirmed at +4": 55 confirmed; four more have the right
peptidoform emitted but q slightly above 1% (`testset_outcomes.tsv`); three lose the
collapse. Arms C/D reached 58 by widening the window at the cost of the degenerate tiers.
Two levers are available without widening: `--precursor-mono-backoff 0` (recovers the
+1-isotope step but measured 8 wrong shifts offline instead of 1), or a `0..1` window once
the correction is trusted, which would remove the (k, X) ≡ (k−1, X') degeneracy inside
the sweep. Neither was run: each is another 2.7 h arm on this host.

## Artifacts (this directory)

`armB_eval.txt`, `armE_eval.txt` (eval_yield / eval_entrap / score_vs_truth ×2 /
agreement ×2), `*_provenance.txt`, `offset_analysis.txt`, `testset_outcomes.tsv` (the 84
scans with both arms' buckets and the fit values), `offline_validation.txt`,
`precursor_mono_dump.tsv` (default-run dump), `armE_precursor_mono_dump.tsv`,
`thermo_trailer_probe.tsv`, `index_equivalence.txt`, and the scripts used
(`run_arm.sh`, `eval_arm.sh`, `offset_analysis.py`, `score_dump.py`).
