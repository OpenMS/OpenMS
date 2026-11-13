# q-value calibration diagnosis (step 1 after PLAN_XIC_IDS)

13 November 2026. Evaluation only: it reads entrapment labels and the ProteoBench DIA-NN reference. No model, score,
population or rule was changed. Scripts and outputs are in `qcal-20261113/`: `diagnose.py` → `DIAG_xic.txt`,
`DIAG_old.txt`; `on_real.py` → `ON_REAL.txt`; `parent.py` → `PARENT.txt`; `parent2.py` → `PARENT2.txt`;
`entrap_origin.py`.

## Question
At q ≤ 0.01 the xic D11FWT arm has entrapment FDP_c 0.016 (raw01) and 0.017 (raw02). Rule 1c allows at most 0.015.
Why do the decoys undercount false targets?

The measure used throughout is the calibration ratio rho = E (1 + 1/r) / D:
- E is the number of entrapment target winners above a score cut.
- D is the number of paired decoy winners above the same cut, using the same TDC as `evaluate_v2.py`.
- rho = 1 means the decoys count false targets correctly. rho > 1 means the q-values are optimistic.
- The interval is a 90 % binomial interval on the E-versus-D split.

## Findings

**1. The xic extraction did not cause this.** rho is above 1 at every q in both rebuilds, and also for DIA-NN given
our library:

| | q ≤ 0.01 | q ≤ 0.05 | q ≤ 0.10 | q ≤ 0.20 |
|---|---|---|---|---|
| xic raw01 | 1.76 (E 14, D 16) | 1.86 | 1.47 | 1.43 |
| xic raw02 | 1.88 (E 14, D 15) | 1.65 | 1.62 | 1.51 |
| old base raw01 | 1.48 (E 11, D 15) | 1.68 | 1.60 | 1.45 |
| old base raw02 | 1.44 (E 10, D 14) | 1.50 | 1.31 | 1.36 |
| DIA-NN 2.2, our library, raw01 (D = q·N) | 1.88 (FDP_c 0.0188) | 1.20 | 1.28 | |
| DIA-NN 2.3.1, ProteoBench 8.65 M library | 0.73 (FDP_c 0.0073) | 0.66 | | |

- The old base passed rule 1c (0.014 / 0.013) with 11 and 10 entrapments. The xic arm failed with 14 and 14. At
  these counts that difference is noise: the intervals overlap widely, and at q 0.05 to 0.20 the two rebuilds have
  the same rho.
- DIA-NN on the same library has FDP_c 0.0188 at its q ≤ 0.01. That fails rule 1c too.

**2. Candidate windows are not unequal between classes.** Human targets, entrapment targets and both kinds of decoy
each have 4.95 to 4.97 windows per precursor, with about 98 % having 5. The xic population has no run-edge windows.
No stratum points to an extraction artefact: peak rank, transition support, RT error, charge, m/z, RT tertile and
fold all show rho of about 1.1 to 2.

**3. Part of the excess comes from how the entrapments were built.**
- The entrapment peptides are FDRBench `_p_target` shuffles of human peptides. They keep the amino-acid composition,
  and therefore the precursor mass.
- The 100K-sampled library contains only about 4.6 % of their parents. The ProteoBench reference identified many of
  the missing parents in these raws.
- In the population, 3.7 % of entrapment precursors have a reference-identified peptide that is a different sequence
  with the same composition and charge. For the other classes the share is 0.7–0.8 %.
- Among accepted entrapment winners, many sit on such a parent eluting within 10 s of the apex: 5 of 14 (raw01) and
  4 of 14 (raw02) at q ≤ 0.01, and 28 of 194 and 27 of 211 at q ≤ 0.10.
- These parents share on average 2.9 of the 6 library fragments (b/y ions, 10 ppm), and 59 % share 3 or more.
- A decoy cannot copy this. Its precursor m/z is shifted (median |Δ| 28 Da), so its signal never comes from a real
  co-eluting peptide of the same mass.
- These are real false discoveries, caused by a real peptide that is missing from the library. Decoys cannot model
  them.
- DIA-NN on our library shows the same pattern: 9 of its 22 accepted entrapments at q ≤ 0.01 have a co-eluting
  parent.
- With the full predicted library, the parents are scored themselves and take their signal, which fits rho 0.7 there.

Decomposition, counting only parents that the reference identified, so the parent count is a lower bound:

| q ≤ 0.01 | FDP_c | FDP_c without parent-explained entrapments |
|---|---|---|
| xic raw01 | 0.016 | 0.010 |
| xic raw02 | 0.017 | 0.012 |
| DIA-NN, our library | 0.019 | 0.011 |

**4. Some excess is left over.** Winners without a known parent still have rho of 1.4 to 1.6 at q 0.10 to 0.20. DIA-NN
on our library is at 1.2 to 1.3 there. Part of this is probably the same mechanism with parents the reference did not
identify. The rest is unexplained.

## What this means
- The rule 1c failure is a property of the benchmark library (shuffled entrapments whose parents are missing)
  combined with any decoy scheme that shifts mass. It is not a bug in the new extraction.
- Under this rule the reference tool also fails on the same library.
- No label-free q-value correction can see this from decoys alone. The fix has to be in scoring: evidence that the
  signal belongs to a different, co-eluting species.

## Options (for a new frozen plan; nothing has been run)
- **A. Keep rule 1c and add label-free interference evidence to scoring.** For example:
  - the fraction of the apex-scan intensity near the precursor that the library fragments do not explain;
  - a fragment-set agreement test between shared and unshared fragments, since a shuffle partner explains only some
    of the 6;
  - the MS1 isotope pattern against the library composition.

  Judge on the unchanged rules 1–4 on both raws. This merges with step 3 (scoring).
- **B. Restate rule 1c against the reference tool on the same library:** pass if FDP_c at q ≤ 0.01 is no worse than
  DIA-NN's 0.0188 on our library. This changes a frozen rule and needs your decision.
- **C. Report both quantities** in every evaluation from now on: FDP_c and the parent-explained share. This is a
  diagnostic only and must not be used for any decision.
