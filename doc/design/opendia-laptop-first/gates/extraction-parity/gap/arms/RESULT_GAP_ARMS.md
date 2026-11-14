# RESULT_GAP: three label-free evidence blocks for the right-peak-but-rejected IDs

Plan: `PLAN_GAP.md`.
- Frozen 2026-11-14T11:33:34Z, with the hashes in `GAP_FROZEN.sha256` checked OK at chain start.
- Amendment 1 is at 2026-11-14T12:32Z (`GAP_AMEND1.sha256`).
- Codon jobs: 28191789 (features, training, EVALUATION) and 28246467 (amendment 1, remaining evaluation steps).
- Outputs are in `run/score/`.

**Control.** The D11FD retrained in this run reproduces arm E's D11FDWT exactly: max abs diff 0.0 on both raws (`CONTROL.txt`). So every difference below comes from the added blocks.

**Amendment 1** (see plan): `compare_gap.py` lost its arm list because it overwrote a variable when it exec'd evaluate_v2's header. Only the evaluation steps were rerun. EVALUATION.txt had been read before the fix. The rules were not changed.

## Verdict: D11FDI and D11FDAGI pass; D11FDA and D11FDG fail (rule 2)

| Rule | A | G | I | AGI |
|---|---|---|---|---|
| 1 FDP_c at q ≤ 0.01 ≤ 0.015 (raw01 / raw02) | 0.0070 / 0.0105 | 0.0062 / 0.0041 | 0.0070 / 0.0083 | 0.0090 / 0.0046 |
| 1 AUC(entrapment > decoy winner) ≤ control + 0.005 (0.5035 / 0.5035) | 0.4994 / 0.4989 | 0.4997 / 0.5003 | 0.4994 / 0.4996 | 0.5007 / 0.5003 |
| 1 null-winner AUC ≤ control + 0.005 (0.5014 / 0.5034) | 0.4962 / 0.4976 | 0.4964 / 0.4976 | 0.4967 / 0.4976 | 0.4967 / 0.4975 |
| 2 ΔN@FDP_c 2 % 90 % interval > 0 (raw01 / raw02) | **−70..+111 / −85..+133 FAIL** | **−107..+56 / −110..+80 FAIL** | +4..+183 / +150..+397 | +59..+272 / +102..+373 |
| 2 ΔN@FDP_c 1 % point > 0, interval not all < 0 | +117 (−54..+373) / +13 (−162..+166) | +33 / +72 | +185 (+30..+519) / +293 (−101..+490) | +152 (−321..+459) / +354 (+160..+561) |
| 3 pooled rho at q ≤ 0.10 ≤ 1.3 | 1.01 | 1.01 | 0.95 | 0.97 |
| 3 pooled rho 90 % interval at q ≤ 0.05 includes 1 | [0.68, 1.05] | [0.65, 1.02] | [0.67, 1.03] | [0.82, 1.21] |
| **Result** | FAIL | FAIL | **PASS** | **PASS** |

**Choice (plan rule).**
- AGI has the highest summed N@FDP_c 1 %: 1846 + 1982 = 3828, against 3800 for I.
- The margin of 28 is far below AGI's own ΔN@1 % interval half-width (390 / 200). So **D11FDI**, the simplest passing block, is named as well.
- In practice the two are equivalent. I is the cheaper of the two, with 4 columns.

**Expectation (written in the plan).** I expected +50 to +150 per block, with I ≥ A ≥ G. I came out larger than that (+185 / +293). A against G is mixed: A was larger on raw01, G on raw02.

## Numbers

| Raw | Arm | N@FDP_c 0.5 / 1 / 2 % | q ≤ 0.01 N (FDP_c) |
|---|---|---|---|
| 01 | control D11FDWT | 1251 / 1694 / 1908 | 1667 (0.0096) |
| 01 | D11FDA | 1441 / 1811 / 1914 | 1712 (0.0070) |
| 01 | D11FDG | 1265 / 1727 / 1866 | 1608 (0.0062) |
| 01 | **D11FDI** | 1505 / **1879** / 1984 | 1707 (0.0070) |
| 01 | D11FDAGI | 1289 / 1846 / 2079 | 1781 (0.0090) |
| 02 | control D11FDWT | 1367 / 1628 / 1903 | 1457 (0.0069) |
| 02 | D11FDA | 1367 / 1641 / 1875 | 1527 (0.0105) |
| 02 | D11FDG | 1474 / 1700 / 1837 | 1466 (0.0041) |
| 02 | **D11FDI** | 1303 / **1921** / 2127 | 1690 (0.0083) |
| 02 | D11FDAGI | 1813 / 1982 / 2072 | 1754 (0.0046) |

**Reference.** DIA-NN 2.2 with its own decoys gives N@FDP_c 1 % of 2067 / 2301 (`diann-decoys-20261114/RESULT.md`). For D11FDI the gap to DIA-NN goes from 373 / 673 to 188 / 380.

## Funnel recovery (recorded, `RECOVERY.txt`; DIA-NN human q ≤ 0.01 IDs whose OpenDIA winner is DIA-NN's peak)

| Raw | Arm | recovered of the 626 / 644 control-rejected | lost of the control-accepted | accepted that are not DIA-NN IDs |
|---|---|---|---|---|
| 01 | A | 113 | 33 | 149 (control 114) |
| 01 | G | 49 | 30 | 127 |
| 01 | I | **149** | 17 | 164 |
| 01 | AGI | 148 | 31 | 146 |
| 02 | A | 48 | 41 | 100 (control 95) |
| 02 | G | 73 | 16 | 108 |
| 02 | I | **234** | 21 | 171 |
| 02 | AGI | 269 | 16 | 191 |

The diagnostic ceiling in `RESULT_GAP.md` estimated about 100–140 of the 626 recoverable at today's decoy rate from the old feature set. **I recovers 149 / 234**, which exceeds that ceiling. The explaining-away block therefore carries evidence that the 20 scorer inputs did not contain. It also adds non-DIA-NN accepts (+50 / +76). Those are within the FDP_c budget, since entrapment FDP_c stays at or below 1 %.

## Reading

1. **Shared-signal explaining away (step 3) is what recovers IDs.**
   - Removing what co-eluting claimants are expected to put on a shared fragment makes weak true peptides score above decoys.
   - It does this without raising FDP_c, entrapment AUC or the null-window AUC, and pooled rho stays at 0.83–0.95.
   - This is the first evidence for the `DESIGN_INTERFERENCE_GRAPH.md` direction. It comes from a one-pass approximation, so the full joint fit may do better. That is not shown.
2. **Mass consistency and dominance (step 1).**
   - These help on raw01 (+117 at 1 %) but not on raw02 (+13).
   - The raw01 gain is the raw that was looked at before the plan, so it should not be over-read.
3. **Anchor-corrected fragment intensities (step 2) carry little new evidence.**
   - The anchor model removes only 5–6 % of the intensity error.
   - G lowers FDP_c (0.006 / 0.004) but does not raise N@2 %.
4. **Combined.** AGI is not better than I alone at 1 %, and it is better at 0.5 % on raw02 only. The blocks overlap.

## Limits
- Each arm was run once.
- Raw01 was seen descriptively before the plan, so raw02 is the cleaner replication. I passes on raw02 with the larger gain (+293).
- The bootstrap 90 % intervals of ΔN@1 % for I include 0 on raw02. The pass rule only required that the interval not lie entirely below 0, and the 2 % interval is clearly above 0 on both raws.

## Next (proposed, not run)
1. **Port D11FDI's block to C++**, in OpenDIA next to the picker:
   - the claimant graph over the isolation window plus RT overlap plus fragment tolerance;
   - each claimant's abundance from its unshared fragments;
   - residual evidence.
   Wall time in Python was 2:02 / 11 GB on raw01, so it will be cheap in C++.
2. **Full joint fit** (the iterative noisy-OR / EM of the design doc) as a new arm against D11FDI.

## Post-hoc validation (2026-11-14, after the ruling; descriptive, does not change the verdict)

**1. Chain wiring.**
- All 28 frozen hashes and the 3 amendment hashes were checked OK.
- make_w and calib_tail ran for all 10 arm × raw combinations, and calib_tail reports "window choice unchanged: True" for each.
- The control reproduces arm E's D11FDWT exactly (max abs diff 0.0).

**2. Synthetic test of features_I** (`scripts/test_features_I.py`):
- A planted claimant on one shared fragment gives the analytic explained share (0.5952, matched to 1e-4).
- The residual library cosine rises (+0.063).
- A precursor in another isolation window, and a pair without RT overlap, get 0 links.
- **Bug found:** the script crashed when no link exists at all, which happens with 0 links in a run. It was fixed by guarding the empty concatenation and the de-duplication mask.
  - The real runs never hit that path.
  - The fixed script reproduces the codon raw01 featI.npy byte for byte (sha256 382896e7…).

**3. Null-window tell check** (`NULL_FEAT_raw01.txt`, label-free; windows with |library RT error| > 3σ, where no peptide is expected):
- AUC(target > decoy) per featA/G/I column lies between 0.497 and 0.506 over windows and over null winners.
- The largest is i_explained at 0.5058. So targets are slightly more often explained by claimants in null windows, but the effect is small.
- No column carries a library or sequence tell of practical size.

**4. Seed sensitivity** (`SEED_COMPARE_<1..4>.txt`, `seed_check.py`):
- D11FD and D11FDI were retrained with every training seed shifted by 1–4, with the same W/T and evaluation.
- ΔN@FDP_c 1 % / 2 % of D11FDI over its own-seed control:

| Raw | seed 0 (frozen run) | 1 | 2 | 3 | 4 | range |
|---|---|---|---|---|---|---|
| 01 @1 % | +185 | +126 | +124 | +136 | +118 | +118..+185 |
| 01 @2 % | +76 | +39 | +81 | +135 | +56 | +39..+135 |
| 02 @1 % | +293 | +342 | +273 | +199 | +351 | +199..+351 |
| 02 @2 % | +224 | +256 | +271 | +189 | +291 | +189..+291 |

- Across seeds, FDP_c at q ≤ 0.01 for D11FDI is 0.0060–0.0071 (raw01) and 0.0012–0.0083 (raw02). The control's range is 0.0073–0.0096 and 0.0042–0.0082, so D11FDI never exceeds the rule-1c limit.
- **The gain is positive at every seed on both raws, and on raw02 it is larger than the control's own seed spread** (control N@1 % 1609–1665).
- The frozen run sits at the top of the raw01 range. On raw01 the typical gain is about +125 at 1 % rather than +185, and +40 to +135 at 2 %. Several raw01 seeds have 2 % intervals that touch 0 (−17..+150, −19..+149), so the raw01 2 % pass in the frozen run is not robust to the seed.
