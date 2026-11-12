# V18 result: charge-sibling co-elution feature (D11FSWT): FAIL

12 November 2026. Plan `PLAN_V18.md`, frozen and hashed at 2026-11-12T10:03:59Z before any V18 feature or score was computed. All hashes re-checked OK at run start. D11FWT stays the base.

## Result

| | raw01 D11FWT | raw01 D11FSWT | raw02 D11FWT | raw02 D11FSWT |
|---|---:|---:|---:|---:|
| N@FDP_c 0.5 / 1 / 1.5 / 2 % | 971 / 1496 / 1616 / 1716 | 1120 / 1258 / 1484 / 1657 | 1127 / 1426 / 1625 / 1736 | 1137 / 1417 / 1554 / 1585 |
| q ≤ 0.01 targets (human / entrapment) | 1614 (1603 / 11) | 1151 (1148 / 3) | 1536 (1526 / 10) | 1477 (1469 / 8) |
| q ≤ 0.01 FDP_c | 0.014 | 0.005 | 0.013 | 0.011 |
| Paired decoys above q ≤ 0.01 | 15 | 10 | 14 | 13 |
| q ≤ 0.01 human peptides | 1503 | 947 | 1431 | 1267 |
| Entrapment-over-decoy winner AUC | 0.4907 | 0.5036 | 0.4958 | 0.4999 |
| Null-window winner AUC | 0.4933 | 0.4950 | 0.4936 | 0.4937 |

Pass rule:
1. PLAN_V2 rules 1-4 against BASE: pass on both raws (N@2 % interval +151..+621 / +102..+408, N@1 % −44..+615 / −115..+335, FDP_c 0.005 / 0.011, AUC 0.5036 / 0.4999 ≤ 0.531 / 0.535).
2. **N@2 % above D11FWT: FAIL on both raws** (1657 vs 1716, 1585 vs 1736; intervals −267..+191, −236..+75).
3. Null-window AUC ≤ D11FWT + 0.005: pass (0.4950 / 0.4937).

At q ≤ 0.01 V18 gains 210 / 169 precursors (entrapment 2 / 5) but loses 673 / 228 (entrapment 10 / 7). Losses include 506 / 119 previously accepted E and 77 / 79 old D. Old D gained: 65 / 82. Net human peptides: −556 / −164. Fewer entrapments, but at a much larger cost in true IDs.

Integrity: the in-runner D11F control reproduced `scores_D11F` and `lin_D11F` exactly (max abs error 0) on both raws before the arm was trained. Sibling-feature tests passed (4). Runtime 2:38 / 2:24 per raw.

## Why it failed (post-failure diagnostic, `DIAG_V18.txt`)

The model learns the feature from its training positives, which are confident targets. Among winners with a sibling, a sibling window within 3 s exists for 60 % / 58 % of accepted E, but only 27-28 % of decoys and 29 % of entrapment and all human winners (look tables). So training puts a large weight on having a nearby sibling window. Among weak candidates, though, that indicator barely separates anything: the all-human vs decoy AUC was only 0.510 / 0.506. Decoys and entrapments with a nearby sibling window are common (about 11 % of decoy winners), and the model promotes them:

- Fraction above the decoy p99, before → after: **decoys with a near sibling 0.92 % → 2.66 %** (raw01), 1.11 % → 2.68 % (raw02). **Entrapment with a near sibling 1.12 % → 3.00 %**, 1.31 % → 3.37 %. **Decoys without a sibling 1.07 % → 0.78 %**.
- The fold-wise decoy tail therefore fills with sibling-supported decoys. Everything without strong sibling support (about 60 % of IDs have no sibling at all) is pushed down relative to that tail. Accepted E with no sibling stayed near the top, but enough fell below the q ≤ 0.01 cut to lose 506 on raw01.
- D with a near sibling did improve (fraction ≥ p99.9: 0.387 → 0.562 raw01, 0.505 → 0.562 raw02). D without a near sibling lost ground.

So charge-sibling support is strong evidence for confident targets, but for weak candidates it is mostly chance co-elution. The look's band counts (78 human / 0 entrapment / 1 decoy) were for a strong sibling (z > 2) on top of an already-near-cut own score. A trained model does not limit itself to that band. It rewards the much more common weak nearby-sibling cases, where decoys and entrapments are as frequent as true targets. A post-hoc rule restricted to that band would be a tuned threshold chosen from the look, which the plan excluded.

## Decision

Retain D11FWT. Stop the charge-sibling feature line. Do not retune (threshold, clip, decoy-only exclusion), because any variant that treats decoys differently from targets breaks target-decoy competition.

What this adds to V16: the remaining Stage D signal we have found is either within the precursor's own window (rewards entrapment), from other runs (V7, rewards entrapment), or from related precursors (V18, rewards decoys and entrapment). The deep Stage D group still has no measured evidence that separates it from false targets.
