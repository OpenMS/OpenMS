# Voxel v7 result: cross-run evidence (D11FXWT): FAIL

Plan: `PLAN_V7.md` (sha256 `cdcaac8cbb…`, written 2026-11-11T14:32:01Z, sha check OK). Deviation in `PLAN_V7.notes`: a restart at 4 torch threads because the shared CPU was oversubscribed. Base: D11FWT. Raws 01 and 02. Entrapment labels are used for evaluation only.

## Numbers (`EVALUATION_V7.txt`, `NULL_V7.txt`, `VS_D11FWT_V7.txt`)

| Raw | Arm | N@0.5 / 1 / 1.5 / 2 % | q ≤ 0.01 N (FDP_c) | Entrap AUC | Null-winner AUC |
|---|---|---|---|---|---|
| 01 | BASE | 900 / 1024 / 1267 / 1314 | 1081 (0.013) | 0.531 | |
| 01 | D11FWT | 971 / 1496 / 1616 / 1716 | 1614 (0.014) | 0.491 | 0.4933 |
| 01 | D11FXWT | 1291 / 1489 / 1617 / **1782** | 1739 (**0.017**) | 0.498 | 0.4939 |
| 02 | BASE | 1215 / 1288 / 1333 / 1376 | 1241 (0.008) | 0.535 | |
| 02 | D11FWT | 1127 / 1426 / 1625 / 1736 | 1536 (0.013) | 0.496 | 0.4936 |
| 02 | D11FXWT | 936 / 1507 / 1659 / **1746** | 1721 (**0.019**) | 0.495 | 0.4943 |

## Verdict: FAIL

**Rule 1 fails.** PLAN_V2 rule 3 requires FDP_c at q ≤ 0.01 ≤ 0.015. D11FXWT gives 0.017 (raw 01) and 0.019 (raw 02).

The other rules pass on both raws:
- Bootstrap against BASE: N@2 % interval above 0; N@1 % interval not entirely below 0.
- Entrapment AUC ≤ BASE.
- N@2 % above D11FWT: 1782 > 1716 and 1746 > 1736.
- Null-winner AUC ≤ D11FWT + 0.005.

The paired bootstrap against D11FWT includes 0 at 1 % and 2 % on both raws.

## Readings (descriptive)

- **Stage D is unchanged.** Stage D (right window, below the 1 % cut) has 870 / 844 against 866 / 929 for D11FWT. Stage C rose to 124 / 140 from 113 / 109. The diagnostic AUC 0.74 of the other run's *score* did not carry over to raw other-run features in this model.
- **Entrapment beat the decoys on cross-run consistency.** FDP_c at q ≤ 0.01 rose on both raws, so entrapment targets gained more from cross-run consistency than decoys did. One possible cause, not tested: an entrapment precursor that matches a real interfering signal in run a matches the same signal in run b (same sample, same RT), while a decoy matches it less often. If so, cross-run agreement is not symmetric between decoys and entrapment, and the decoy-based q-value underestimates the error.
- **Base unchanged.** D11FWT stays the base.
