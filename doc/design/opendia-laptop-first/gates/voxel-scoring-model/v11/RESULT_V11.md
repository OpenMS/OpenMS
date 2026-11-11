# Voxel v11 result: native-sampling fragment evidence on unique fragments only (D11FUWT). FAIL

**Plan.** `PLAN_V11.md`, written 2026-11-11T21:29:48Z. `chain_v11.sh` checked the sha (OK) before training. There were no deviations.

## Result (`EVALUATION_V11.txt`, `NULL_V11.txt`, `VS_D11FWT_V11.txt`)

| Raw | Arm | N@FDP_c 0.5 / 1 / 1.5 / 2 % | q ≤ 0.01 N, FDP_c | AUC (entrapment > decoy) |
|---|---|---|---|---|
| 01 | D11FWT | 971 / 1496 / 1616 / 1716 | 1614, 0.014 | 0.4907 |
| 01 | D11FUWT | 944 / 1133 / 1483 / **1552** | 1725, **0.030** | 0.4965 |
| 02 | D11FWT | 1127 / 1426 / 1625 / 1736 | 1536, 0.013 | 0.4958 |
| 02 | D11FUWT | 1242 / 1307 / 1460 / **1518** | 1651, **0.032** | 0.5009 |

- **Rule 1.3 fails on both raws:** FDP_c is 0.030 / 0.032. This is worse than V10 (0.023 / 0.018).
- **Rule 2 fails on both raws:** N@2 % is 1552 / 1518 against 1716 / 1736.
- **Rule 3 passes:** null-winner AUC is 0.4950 / 0.4950, against 0.4933 / 0.4936 for D11FWT (+0.0017 / +0.0014).

**Funnel (raw 01):** D went from 866 to 1207 and E from 1413 to 1100.

**Paired bootstrap against D11FWT, @2 %:** −340..+17 on raw 01 and −332..+10 on raw 02.

## What the mask did (`NATU_SEPARATION_raw01.txt`, `DIAG_V11_ENTRAP.txt`; descriptive, labels only to count)

- **Fraction masked:** 40.6 % / 40.5 % of the 12-fragment slots are shared under the rule. The median number of unique fragments is 7 for both targets and decoys. The mean is lower for targets (6.78) than for decoys (7.28). This is consistent with target fragments colliding with other library target fragments more often than shuffled decoy fragments do.
- **Pre-plan look (raw 01):**
  - The key check was not met. Entrapment against decoys, score-matched, was 0.538 on `apex_logp` against 0.548 for E.
  - Masking also removed evidence from true human peaks: all human against decoys fell from 0.587 to 0.548.
- **Accepted-set change at q ≤ 0.01, D11FUWT against D11FWT:**

| | Gained human / entrapment / decoys | Lost human / entrapment / decoys |
|---|---|---|
| Raw 01 | 284 / 22 / 14 | 188 / 7 / 13 |
| Raw 02 | 253 / 24 / 12 | 154 / 8 / 11 |

- **Gained precursors** have many unique fragments (median 8.5–9.5) and strong evidence. Entrapment `apex_logp` on unique fragments is 5.1–5.4, close to human at 5.7–6.1.
- **Lost precursors** have few unique fragments (median 2.5–5).

**Reading (softened).**
- Masking shared fragments did not remove the entrapment advantage in native evidence. Entrapment winners still match their unique fragments in real spectra.
- So the V10 excess is not mainly explained by library-visible fragment sharing. It may come from:
  - co-isolated signal that is not in the library;
  - near-isobaric fragments beyond 2 × tolerance;
  - entrapment peptides that share sequence features with real human peptides.
- The mask also penalises targets in crowded windows, because they lose fragments. The network then sorts by "many unique fragments", and that trait is no more specific to true peptides than to entrapment.

## Outcome

- **FAIL.** The base stays D11FWT.
- As agreed before the run, this ends the native-sampling fragment-evidence line for stage D (V9, V10, V11). The next step goes back to the drawing board, with the user.
