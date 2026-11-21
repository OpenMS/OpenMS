# RESULT_JOINT: charge-state sibling evidence is real new information; interference marks the missed peaks (21 Nov 2026, diagnostic)

The gate is in `PRECHECK_JOINT.md` (`8fbde5598b`). The script is `precheck_joint.py` (sha256 `be1934a7…`), run on the PRECHECK_CEIL0 feature tables at seed 0. The full output is in `EVAL_JOINT.txt`.

## Sibling (charge-state) features

Windows with a co-eluting sibling: D r0 target 42.9 % vs decoy 42.4 %, H r0 15.0 % vs 13.9 %. The features are symmetric between targets and decoys, as designed.

**HeLa IDs at own q ≤ 0.01:**

| Arm | D r0 | H r0 | D r7 | H r7 |
|---|---:|---:|---:|---:|
| A tool | 19,097 | 34,702 | 15,859 | 30,888 |
| B label-free GBM | 20,261 | 37,178 | 15,961 | 31,935 |
| **Bs = B + siblings** | **22,761 (+12.3 % vs B)** | **41,941 (+12.8 %)** | 15,571 (−2.4 %) | **37,152 (+16.3 %)** |
| C DIA-NN ceiling | 20,644 | 37,444 | 15,583 | 32,462 |
| **Cs = C + siblings** | **25,077 (+21.5 % vs C)** | **41,100 (+9.8 %)** | 16,365 (+5.0 %) | 33,964 (+4.6 %) |
| DIA-NN | 27,901 | 44,866 | 25,881 | 41,231 |

**At r7, N at FDP_c 1 % (lower bound):**

| Arm | D | H |
|---|---|---|
| A | 15,465 (1.12 %) | 30,543 (1.20 %) |
| B | 15,387 (1.21 %) | 31,250 (1.26 %) |
| Bs | 14,780 (1.27 %), **−3.9 %** vs B | **36,588 (1.18 %), +17.1 %** vs B |
| C | 15,330 (1.01 %) | 31,666 (1.24 %) |
| Cs | 16,031 (0.98 %) | 33,362 (1.16 %) |

## Interference profile (tool score, r0)

| Class | D: `i_explained` median / share ≥ 0.5 | H: median / share ≥ 0.5 |
|---|---|---|
| A_low (right window best, not accepted) | 0.454 / **45.2 %** | 0.468 / **47.7 %** |
| Accepted HeLa targets | 0.298 / 28.0 % | 0.262 / 25.6 % |
| Decoys within ±2 of the 1 % cut | 0.454 / 46.5 % | 0.402 / 43.2 % |

## Decision rule

| Rule | Result | Verdict |
|---|---|---|
| 1. Bs/B ≥ +3 % at r0 on D and H; ≥ +2 % N@FDP_c at r7 on D and H; lower bound ≤ +0.2 pp | r0 +12.3 / +12.8 %, r7 H +17.1 % (lower −0.08 pp), **but D r7 −3.9 % (lower +0.06 pp)** | **Fails on D r7** |
| 2. A_low `i_explained` ≥ 0.5 share ≥ 25 % and ≥ 1.5× accepted, on D or H at r0 | D 45.2 % (1.62×), H 47.7 % (1.86×) | **Passes on both**: joint deconvolution is a lever |
| 3. Cs/C < +3 % on both means the sibling information is already implicit | Cs/C +21.5 / +9.8 % (r0), +5.0 / +4.6 % (r7) | **The sibling information is new** |

## Reading

- **Sibling evidence is the largest new signal found so far.**
  - With DIA-NN labels, adding four sibling features lifts the ceiling by +21.5 % (D) and +9.8 % (H) at r0. Every scorer change tested before moved it by less than 1 %.
  - It also works label-free: +12-13 % at r0, and +17 % at true FDP 1 % on H r7 with *better* calibration. That would take H from 74 % to about 87 % of DIA-NN.
  - Why it works: a true peptide elutes at the same time in every charge state it is seen in. A decoy or noise peak does not.
- **D r7 is the exception.** The label-free learner already fails there (B ≈ A, while Bt gained +4.3 %). The ceiling still gains with siblings on D r7 (Cs: N@FDP_c 1 % +4.6 % vs C, lower bound 0.98 %), so the information is present. The plain self-trained GBM does not use it on the crowded, large-library run.
  - A plausible cause: at r7 most of D's windows are entrapment or decoy windows with co-eluting siblings by chance (33 % of windows). Self-training with all decoys as negatives then fits that chance structure.
  - Whether the hard-decoy variant (Bt + siblings) fixes this is the follow-up question.
- **Deconvolution is also indicated.** About 46 % of A_low peaks are at least half explained by co-eluting claimants, 1.6-1.9× the rate of accepted targets. So the missed peaks are mostly interfered peaks.
  - Caveat: decoys near the cut are equally interfered (43-47 %). Deconvolution helps only if it uncovers target-specific residual signal, not just more interference.
  - The current resolver is one pass with a median scale, not a joint fit. A joint NNLS fit of the co-eluting library patterns is the next test.
- **Caveat:** this is one seed on one run per type, and at r0 there is no entrapment check. The r7 H result has one.

## Follow-up fixed before running (PRECHECK_JOINT2)

Seed robustness and the D r7 failure:

- Add arms Bt (hard decoys) and **Bts** (hard decoys plus siblings) to the same script, `precheck_joint2.py`. The method is otherwise unchanged.
- Run on:
  - the r ≈ 1.56 feature tables of RESULT_BT156 (`bt156-20261121/{D,H}_s{0,1}`, r = 1.1225);
  - D r7 and H r7 seed 0.
- **Rule:** the better of Bs and Bts goes to a C++ opt-in (sibling features in the scorer input) if:
  - at r ≈ 1.56 it gains ≥ +3 % N@FDP_c 1 % over A on all 4 runs, with the mean lower bound rising by ≤ +0.2 pp per run;
  - and at r7 it gains ≥ +2 % over A on both D and H.

  Otherwise it is recorded, with the per-run breakdown.
