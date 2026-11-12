# V19 result: native-sampling peak geometry replaces search-run summaries (D11FGWT): FAIL

12 November 2026. `PLAN_V19.md` was frozen and hashed at 2026-11-12T11:34:02Z, after the raw01 feature look and before raw02 features or any V19 score. All hashes were OK at run start. The in-runner D11F control reproduced exactly on both raws (max abs error 0). D11FWT stays the base.

## Measurement finding (this holds regardless of the arm)

Three D11F inputs come from the search's streaming signal run, `OpenDIASearch::updateCandidate_`. A run only closes after a 30 s gap with no matched signal (`peak_gap_seconds = 30`):

| Input | What it actually measures | Raw01 median (accepted / DIA-NN-only / decoy) |
|---|---|---|
| `log_peak_width` | duration of the whole signal run | 46 / 54 / 66 s (DIA-NN FWHM 3.6 / 3.2 s) |
| `log_matched_scans` | scans in the run | 14 / 16 / 19 |
| `transition_support` | max fragments matched in **any** scan of the run | 6 / 4 / 2 |

Spearman of our width against DIA-NN FWHM is 0.02. Native FWHM (`build_peak.py`) gets 0.56. Width and scan count are run-length descriptors. Decoys have longer runs than real peptides, because random matches keep a run open (AUC for width, DN vs decoy: 0.41).

## Result

| | raw01 D11FWT | raw01 D11FGWT | raw02 D11FWT | raw02 D11FGWT |
|---|---:|---:|---:|---:|
| N@FDP_c 0.5 / 1 / 1.5 / 2 % | 971 / 1496 / 1616 / 1716 | 1328 / 1445 / 1458 / 1570 | 1127 / 1426 / 1625 / 1736 | 1265 / 1456 / 1469 / 1580 |
| q ≤ 0.01 targets (human / entrapment) | 1614 (1603 / 11) | 1456 (1447 / 9) | 1536 (1526 / 10) | 1370 (1364 / 6) |
| q ≤ 0.01 FDP_c | 0.014 | 0.012 | 0.013 | 0.009 |
| Paired decoys above q ≤ 0.01 | 15 | 13 | 14 | 12 |
| Null-window winner AUC | 0.4933 | 0.4945 | 0.4936 | 0.4939 |

Pass rule:
1. PLAN_V2 rules 1-4: pass on both raws.
2. **N@2 % above D11FWT: FAIL on both raws.** 1570 vs 1716 and 1580 vs 1736. The intervals −249..−55 and −235..−7 are entirely below zero.
3. Null AUC: pass.

At q ≤ 0.01 the arm gains 50 / 51 precursors (32 / 38 old D) and loses 208 / 217 (109 / 118 accepted E). Entrapment counts fall slightly (9 / 6 against 11 / 10), but only together with the lost IDs.

**The curve shape changed.** N@0.5 % rises from 971 to 1328 on raw01 and from 1127 to 1265 on raw02, so the native features make the top of the list cleaner. N@1 % is about flat (1445 vs 1496, 1456 vs 1426). The loss is concentrated at 1.5-2 %. The old run-length inputs carried ranking information for the middle of the list that the native replacements do not provide.

## Reading

- Replacing the three inputs removed the run-length descriptors, and the model lost ranking power in the middle of the list. Those descriptors were therefore informative in combination with the other inputs, even though they do not measure the peak. Per window, real peptides still have shorter runs than decoys (DN 54 s, decoys 66 s). Why the combination helps was not measured here.
- Native support (fragments matched in at least 2 of 3 apex scans) is a stricter count than the run maximum. It separates confident targets as well (AUC 0.98 vs 0.99), but weak DN targets less well (0.84 vs 0.91). Weak targets lose support they had under the looser count.
- So the measurement problem is real, but swapping one set of inputs for the other loses information. Keeping both sets would add 3 inputs, which this plan excluded on purpose (a "replace, not add" test).

## Decision

Retain D11FWT. Do not retune this arm (offsets, ±2-scan search, half-maximum rule). The measurement finding stands as an input for the C++ scorer design: OpenDIA's `peak_width_seconds` and `matched_scans` are not peak descriptors and should not be named as such.
