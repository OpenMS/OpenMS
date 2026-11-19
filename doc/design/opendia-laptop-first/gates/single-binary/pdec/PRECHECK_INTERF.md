# PRECHECK_INTERF: OpenDIA's interference features already flag the shared-signal false targets (diagnostic)

Written 19 Nov 07:20. Runs: `entlib-odpfeat-20261119`, D and H seed 0, CUR, `-decoys predict -features_out`, r ≈ 7. The
predicted library is the same as PLAN_ODP (sha256 76f2a6c2…). Labels are used for evaluation only, and no verdict
changes. Script `interfdiag.py`, output `INTERFDIAG.txt`.

**Classes** (accepted precursor winners at q ≤ 0.01, D and H pooled):
- **ES:** entrapments sharing ≥ 3 fragments with a co-eluting identified HeLa target (n 45).
- **EO:** other entrapments (278).
- **DE:** decoys of entrapment pairs (235).
- **HT:** HeLa targets (sample of 2,000).

| feature | median ES / EO / DE / HT | AUC(ES > DE) | AUC(EO > DE) | AUC(HT > DE) |
|---|---|---|---|---|
| `i_explained` (share of area explained by claimants) | **1.000** / 0.831 / 0.857 / 0.764 | **0.744** | 0.483 | 0.416 |
| `i_libcos_resid` (library cosine of the residual) | **0.000** / 0.510 / 0.496 / 0.617 | **0.292** | 0.510 | 0.601 |
| `i_dlibcos` (residual minus observed cosine) | **−0.927** / −0.436 / −0.466 / −0.365 | **0.264** | 0.497 | 0.562 |
| `i_nclaim` (log 1 + claimants) | **5.32** / 4.09 / 4.11 / 3.99 | **0.747** | 0.504 | 0.469 |
| final score | 33.0 / 26.7 / 23.8 / 43.4 | 0.707 | 0.590 | 0.881 |

**Share of each class that is fully explained (`i_explained` ≥ 0.999, residual cosine 0):**

| class | share |
|---|---|
| ES | **53 %** |
| EO | 26 % |
| DE | 23 % |
| HT | 19 % |

## Reading
1. **The information is there.** The interference graph already sees the shared-signal false targets.
   - Half of them are entirely explained by co-eluting claimants, with a zero residual library cosine. Decoys and true
     HeLa targets reach that state only 19-26 % of the time.
   - On ES against DE, all four interference features separate in the "more interfered" direction (AUC 0.26-0.29 or
     0.74-0.75).
2. **The scorer does not use it, or uses it the wrong way.**
   - ES still scores *higher* than decoys (AUC 0.707) and is accepted.
   - Among all accepted winners, 19.9 % (6,011) are fully explained. Most of them are real HeLa targets (5,843). So
     "fully explained" alone cannot be a veto: true targets that share a fragment set with a co-eluting HeLa precursor
     also look like this.
   - The features are weak overall (AUC HT > DE 0.42-0.60). A linear scorer trained on target-vs-decoy gives them little
     weight, because for most windows they do not separate targets from decoys.
3. **What would help, label-free:** use the explained share **conditional on evidence strength**.
   - For a true target, the claimant explains only part of its signal and its own unshared fragments carry the rest.
   - For a shared-signal false target, the shared fragments are all it has.
   - A feature such as *unexplained supporting fragments* (the count of fragments with residual area above noise after
     claimants) or *score after removing claimed area* separates the two. A linear scorer can use either directly.
4. **Size of the effect.** ES is 45 of 323 accepted entrapments (14 %). Removing their excess alone would close roughly
   a third of the remaining gap. That is an estimate from these counts, not a measured result.

## Proposed next step (design only, not started)
- **What:** a residual-evidence feature in the interference resolver. Per window, compute the number of fragments whose
  observed central area exceeds the claimed expectation by more than k × noise. Optionally add a residual-only
  coelution score.
- **Properties:** run-derived, identical for targets and decoys.
- **How to test:** a frozen plan on the r ≈ 7 library with the PLAN_ODP criteria. D is design, H is held out.
