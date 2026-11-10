# Round 7 result: refined peak picking (P1) and refined quantification

`PLAN7.md` was written before the run (sha256 in `PLAN7.md.sha256`, time in `PLAN7.written_at`). Code: `dfbe38b470`. Evaluator: `evaluate7.py`, output in `EVALUATION7.txt`.

## Equivalence

- HeLa default output: IDENTICAL (fresh directory `hela-default-cmp7`).
- Internal check on both raws: `crossfit.tsv` columns 1–6 and 9, and `features.tsv` columns 1–9 and 15–16, are IDENTICAL to base (`CMP7.txt`).

## ID rule: FAIL on both raws

| | N at FDP_c 0.5 / 1 / 1.5 / 2 % | Rule 1 (N@1 %) | q ≤ 0.01: N / FDP_c | q ≤ 0.05: N / FDP_c | Time, RSS |
|---|---|---|---|---|---|
| raw 01 base | 899 / 1029 / 1270 / 1321 | | 1068 / 0.013 | 1638 / 0.065 | 12:58, 24.9 GB |
| raw 01 P1 | 916 / 1024 / 1252 / 1305 | 1024 < 1050, fail | 1106 / 0.013 | 1642 / 0.065 | 18.2 min, 23.7 GB |
| raw 02 base | 1230 / 1298 / 1331 / 1372 | | 1236 / 0.007 | 1788 / 0.074 | 14:17, 24.7 GB |
| raw 02 P1 | 1243 / 1323 / 1346 / 1415 | 1323 < 1324, fail | 1234 / 0.005 | 1820 / 0.075 | 18.8 min, 23.6 GB |

- Rules 2, 3 and 4 pass on both raws. Rule 1 fails on both raws. On raw 02 it fails by 1 ID.
- Paired bootstrap, P1 minus base, 90 % interval:
  - raw 01: −24..+36 at 1 %, −39..+31 at 2 %.
  - raw 02: −1..+38 at 1 %, −26..+45 at 2 %.
- There is no clear loss and no clear gain. The difference is in the noise of this metric.
- Mechanism (reported only): the entrapment-over-decoy AUC of the boundary features is 0.508–0.510 with P1, against 0.509–0.512 for base. The stage-2 score stays at 0.500–0.501.
- Wall time increased by 4–5 min (18:12 and 18:47 against 12:58 and 14:17). Both P1 runs started after the machine suspend (about 08:09–09:02), so the suspend does not explain this. During the raw 02 run, small Python scripts also ran on the machine. The 35 min limit holds. The extra time of the re-pick is not measured separately.

## Quantification rule: FAIL

Population: base 1 % non-entrapment precursors in both raws, at the base winning windows (953 precursors).

| | Value |
|---|---|
| Robust SD of log2(01/02), naive | 0.1398 |
| Robust SD of log2(01/02), refined | 0.1475 (ratio 1.055) |
| Bootstrap 90 % of refined − naive | −0.0026..+0.0167 |
| Coverage of refined | 922 / 953 (0.967) |

- Criterion 1 needed a ratio ≤ 0.95. The ratio is 1.055. FAIL.
- Criterion 2 needed the interval below 0. It is not. FAIL.
- Criterion 3 (coverage ≥ 80 %) passes.

Reported only:
- Refined boundaries differ from naive in 26 % (raw 01) and 32 % (raw 02) of windows. The median width does not change (15.0 s and 14.1 s).
- 69 % (raw 01) and 57 % (raw 02) of windows have 5 or 6 fragments with r ≥ 0.9.
- Spearman with DIA-NN `Precursor.Quantity` on raw 01 (907 shared): 0.892 for naive and 0.892 for refined.
- The DIA-NN 01/02 robust SD is 0.335 over 72k precursors. Our 0.14 is over a much smaller and brighter set, so the two numbers are not comparable.

## Reading

- On the precursors that we identify at 1 %, most fragments already agree with the consensus shape. Removing the few fragments that do not agree does not make the replicate ratio more stable. It makes it slightly less stable, possibly because fewer fragments give a noisier sum (not tested).
- Refined picking changes the boundary in one quarter to one third of windows, but it does not change the IDs by more than the noise.
- Both options stay opt-in. The defaults do not change.
- With rounds 6 and 7, transition refinement does not help detection or quantification on these runs. The work moves to the single voxel scoring model (`../../voxel-scoring-model/`).
