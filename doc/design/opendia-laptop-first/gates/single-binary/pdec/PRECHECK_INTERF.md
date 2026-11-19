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

## Addendum (07:40): what a label-free shared-signal rule could gain
Three follow-up checks on the same two runs. Labels are used for evaluation only.

1. **Siblings are not the explaining claimant** (`SIBDIAG.txt`).
   - Of the 5,843 fully explained HeLa targets, only 9 % have a co-eluting sibling, meaning the same sequence at another
     charge or modification state. 5,317 have none.
   - So the "fully explained" HeLa targets are mostly explained by *other* co-eluting HeLa peptides. That is a
     limitation of `i_explained`, which overattributes in dense regions. A sibling exemption would not fix it.
2. **Stronger-claimant fragment sharing** (`STRONGCLAIM.txt`). The rule: ≥ k of a winner's fragments match a co-eluting,
   higher-scoring accepted winner.
   - At k ≥ 3, the share of each class that matches is: entrapments 9.6 / 11.4 %, decoys of entrapment pairs
     7.9 / 5.8 %, HeLa targets 1.5 / 1.1 % (D / H).
   - That is an enrichment, but with an overlap.
3. **Fair simulation** (`RULESIM.txt`). The rule is applied to targets and decoys alike among winners with q ≤ 0.05, and
   q is recomputed by plain TDC. The baseline with plain TDC reproduces the tool exactly: N 8,483 against 8,484 on D, and
   21,386 on H.

| run | baseline lower / combined | rule k = 3: N change, lower / combined | rule k = 4: N change, lower / combined |
|---|---|---|---|
| D s0 | 0.0111 / 0.0127 | −0.26 %, **0.0106 / 0.0121** | −0.94 %, 0.0106 / 0.0121 |
| H s0 | 0.0107 / 0.0122 | −0.37 %, **0.0102 / 0.0117** | −0.77 %, 0.0102 / 0.0117 |

**Reading.**
- A label-free, run-derived rule (drop a winner whose fragments are claimed by a stronger co-eluting winner) moves the
  true FDP at q ≤ 0.01 by about −0.05 percentage points at a cost of 0.3-0.9 % of IDs.
- That is real but small: roughly a quarter to a third of the remaining excess, in line with the 14 % estimate above. It
  would bring H to about 1.02-1.17 %. D would still exceed the 1.0 / 1.2 % limits (1.06 / 1.21 %).
- **This is one seed per run, and the effect (about 20 entrapments) is within seed noise.** A frozen plan would need
  5 seeds. Even then it is unlikely to meet the PLAN_ODP limits on its own.

**Recommendation.**
- **Not yet:** a frozen plan for this rule alone. The expected gain does not close the gap.
- **Better next step:** accept that predicted decoys bring q to about 1.1× at r ≈ 7. For realistic libraries, the
  r = 1.56 results were calibrated at 0.53-1.02 %.
- **Or treat this as a residual to model rather than remove.** For example, add a small, label-free q correction from
  the decoy-vs-target fragment-sharing rate. That is a separate design question for the user.
