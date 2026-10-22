# #10693, point 1: why PIP-ECHO produces so few candidate transfers

Data and pipeline as in `../pr10637` (ProteoBench HF-X, #10585 recipe, plain HYE database, `-proteinFDR 0.01`,
PIP-ECHO `fdr 0.05`, shape features off). The arm is the `off_q05` arm of `../pr10637`, rerun from the same `-feat_dir`
checkpoints with a diagnostic build. Its MSstats output is byte-identical to that arm's.

## Method

- `diag_link.patch` (local investigation only, not for a PR): with `PIPECHO_DIAG=<file>`,
  `Impl::link_donors_and_acceptors()` writes one line per donor feature and acceptor run: whether it found a
  candidate, the number of unidentified features of that charge within the m/z tolerance anywhere in the run, how many
  of them lie inside the local RT window, and the RT distance of the nearest one. Applied to the PR head `ce6eb11b`
  without reconfiguring, so the build stamp matches the checkpoints. (The patch writes m/z with 6 significant digits;
  the scripts below use the theoretical m/z of the sequence instead.)
- `candidates.py`: for every reported precursor (sequence + charge) and every run without an identified feature of
  it, the status of that pair:
  - `accepted`: the transfer was accepted;
  - `candidate`: there was a candidate, but it was rejected by the transfer FDR;
  - `outside_rt`: unidentified features at that m/z exist, but none inside the RT window;
  - `no_feature`: there is no unidentified feature at that m/z.
- `candidates_detail.py`: how far outside the window the nearest feature lies, whether an identified feature of another
  precursor occupies the spot, and the species and condition of the missing run.
- `xic_check.py`: MS1 signal in the alignment reference run (A_02, where the predicted RT is the run's own RT). It
  extracts M0 and M+1 (±10 ppm) inside the local window, smooths both over 9 scans, and calls the signal present if
  ≥ 10 scans within ±15 scans of the M0 apex carry both isotopes and the traces correlate (Pearson ≥ 0.6). Run for
  1,500 pairs per status; the null is the same test with the precursor mass shifted by +5 or −11 Da.
- `run_seeds.sh` (detection of A_01 alone with `-debug 667`), `seed_check.py`, `seed_skip.py`: check whether a
  FeatureFinderMultiplex seed existed at the missing spot and whether FFID skipped it (seed within ±30 s of a PSM at
  its m/z or one of 5 lower isotopes, `addSeeds_`) or dropped it later.
- `upper_bound.py`: Monte Carlo estimate of how many precursors would reach ≥ 3 runs if every pair with real signal got
  a candidate.

## Results

Reported precursors × runs without an identification (`results/cand_status.tsv`), by the number of runs the precursor
is identified in:

| identified in | pairs | accepted | candidate, rejected | outside RT window | no feature |
|---|---|---|---|---|---|
| 1 run | 50,010 | 11.2 % | 1.5 % | 50.6 % | 36.8 % |
| 2 runs | 38,296 | 15.5 % | 1.2 % | 48.1 % | 35.2 % |
| all | 141,706 | 17.3 % | 1.2 % | 47.6 % | 33.9 % |

1. **Scoring is not the limit.** 93 % of the candidates are accepted (24,460 of 26,186). Of the 19,576 precursors
   identified in 1–2 runs, 4,705 reach ≥ 3 runs; accepting every candidate would make it 5,177.
2. **The RT window is not the limit.** For `outside_rt`, the nearest feature lies more than 10 half-windows away in
   82 % of the pairs, and within 2 half-windows in only 2 %. These are unrelated features at the same m/z. The median
   half-window is 22.8 s.
3. **The acceptor feature is missing.** In 82 % of the pairs there is no unidentified feature at the spot. Only 5–8 %
   of the spots are taken by an identified feature of another precursor. Missing candidates follow abundance: for
   E. coli, 92 % of the pairs have no candidate in condition A (4× lower) against 79 % in B; for yeast, 87 % in B
   (2× lower) against 79 % in A.
4. **The signal is there for about half of them.** In A_02, the M0/M+1 test finds signal in:
   - 77 % of the accepted pairs (sensitivity);
   - 43 % of `no_feature` and 49 % of `outside_rt` pairs;
   - 8–9 % and 14 % at shifted mass (null).

   Corrected for sensitivity and null, about 51 % and 56 % of the pairs without a candidate have real co-eluting
   signal. Its median apex intensity is 1.5–2.5e6, as for accepted transfers and far above the seed threshold of 1e4.
5. **Where the feature is lost (A_01, 100,326 seeds; 29,418 unidentified FFID features in the run):** of the pairs
   without a candidate,
   - 72.8 % have no seed;
   - 20.3 % have a seed that FFID skipped next to a PSM. The blocking PSM matches through an isotope offset of 1–5 in
     47 % of the cases; it is a decoy in 16 % and above q 0.01 in 39 %; it is the same precursor in only 3 %;
   - 6.8 % have a seed that FFID dropped later.
6. **Upper bound.** If every pair with real signal got a candidate (signal share 0.4–0.6, acceptance as today), about
   15.8k–18.4k instead of 4.7k of the 19,576 precursors would reach ≥ 3 runs, i.e. +11k to +14k precursors. This is
   rough: it rests on one run's signal test, and signal at the right m/z and RT is not proof of the right peptide.

## Reproduce

```bash
# diagnostic arm (from ../pr10637 checkpoints and IDs)
git apply diag_link.patch   # on ce6eb11b, build ProteomicsLFQ without reconfiguring, install to /opt/openms-diag
./run_diag.sh               # writes link.tsv and the arm's outputs
python3 candidates.py link.tsv hye_off_q05/out.consensusXML.gz cand
python3 candidates_detail.py link.tsv cand_pairs.tsv.gz hye_off_q05/out.consensusXML.gz
python3 xic_check.py LFQ_..._A_Sample_Alpha_02.mzML cand_pairs.tsv.gz 1 1500 1 [5.0|-11.0]
./run_seeds.sh              # A_01 alone, -debug 667: debug_seeds_fraction_1_1.featureXML, debug_fraction_1_1.featureXML
python3 seed_check.py link.tsv cand_pairs.tsv.gz 0 debug_seeds_fraction_1_1.featureXML debug_fraction_1_1.featureXML
python3 seed_skip.py        # paths inside; reuses seed_check.py
python3 upper_bound.py cand_pairs.tsv.gz 0.5 0.93
```
Python: `pandas`, `pyopenms` 3.6.0 (the proteobench venv of `../pr10637` plus `pip install pyopenms`).
