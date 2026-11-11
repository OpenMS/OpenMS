# Voxel v10 result: native-sampling fragment evidence, empirical per-m/z null (D11FEWT). FAIL

**Plan.** `PLAN_V10.md`, written 2026-11-11T20:51Z. The sha was checked OK by `chain_v10.sh` before training.

## Result (`EVALUATION_V10.txt`, `NULL_V10.txt`, `VS_D11FWT_V10.txt`)

| Raw | Arm | N@FDP_c 0.5 / 1 / 1.5 / 2 % | q ≤ 0.01 N, FDP_c | AUC (entrapment > decoy) |
|---|---|---|---|---|
| 01 | D11FWT | 971 / 1496 / 1616 / 1716 | 1614, 0.014 | 0.4907 |
| 01 | D11FEWT | 1267 / 1347 / 1550 / **1676** | 1730, **0.023** | 0.4964 |
| 02 | D11FWT | 1127 / 1426 / 1625 / 1736 | 1536, 0.013 | 0.4958 |
| 02 | D11FEWT | 1377 / 1448 / 1554 / **1695** | 1678, **0.018** | 0.5012 |

- **Rule 1.3 fails on both raws:** FDP_c is 0.023 / 0.018.
- **Rule 2 fails on both raws:** N@2 % is 1676 / 1695 against 1716 / 1736.
- **Rule 3 passes:** null-winner AUC 0.4937 / 0.4940.

**Funnel:** D 866 → 1019 (raw 01) and 929 → 943 (raw 02); E 1413 → 1293 and 1362 → 1387.

**Paired bootstrap against D11FWT, @2 %:** −199..+135 on raw 01 and −217..+168 on raw 02.

## What the new features did (`DIAG_V10_ENTRAP.txt`; descriptive, labels only to count)

At q ≤ 0.01, D11FEWT changes the accepted set by a few hundred targets per raw:

| | Gained human / entrapment / decoys | Lost human / entrapment / decoys |
|---|---|---|
| Raw 01 | 262 / 16 / 12 | 155 / 7 / 11 |
| Raw 02 | 258 / 11 / 13 | 121 / 6 / 12 |

- **Gained targets** have strong native evidence: median `apex_logp` 6.4–6.7 for human and 5.8–6.0 for entrapment.
- **Accepted decoys** sit at a similar level: 5.5–5.9.
- **Lost human targets** have weak native evidence: 3.7–3.9.

**Reading (softened).**
- The empirical-null evidence is high for entrapment targets as well as for human ones. Entrapment peptides come from another proteome and are present in the library at equal size. A few can match real co-eluting signal in a 2 Th window: shared fragments with a co-isolated human peptide.
- Decoys rarely do that, because they are shuffled sequences. So the feature carries a target-against-decoy asymmetry, not only true-against-false information. It also moves the score toward strong native evidence and away from the weak true peaks this plan aimed to recover (D grew on raw 01).
- This matches the pre-plan risk: entrapment against decoys was 0.55 on `apex_logp`.

## Outcome

- **Fail.** D11FWT stays the base.
- **Native-sampling fragment evidence did not recover stage-D misses under either chance model.**
  - The uniform null was mis-calibrated (V9).
  - The empirical null is calibrated on decoys (P(L ≥ 2) = 0.033). But it rewards entrapment about as much as weak human peaks, so the FDP_c at q ≤ 0.01 rose.
- **Next step: the user's decision.** No further native-evidence variant is run without a new idea that addresses the entrapment asymmetry. One example would be scoring only fragments not shared with any co-isolated library precursor; that is untested.
