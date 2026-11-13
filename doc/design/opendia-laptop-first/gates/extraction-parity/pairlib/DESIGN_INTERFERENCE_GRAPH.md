# Design draft: label-free interference resolution as a hierarchical Bayesian graph (EPIFANY-style)

Status: draft for discussion. Nothing here is frozen or run.

## Problem
Two library precursors can explain the same physical signal:
- an entrapment and its co-eluting parent (same composition, same m/z, isolated in the same DIA window), and
- more generally any precursors in one isolation window whose fragments fall on the same m/z at the same RT.

Today every precursor is scored alone, so both can win on one peak. DIA-NN removes this with a greedy step that drops
the weaker precursor when two explain the same signal ("Removing interfering precursors"). We want the same effect
from a probabilistic model that shares evidence instead of a hard greedy cut, with no labels and no fixed constants.

## Graph (three levels, as in EPIFANY's protein / peptide / PSM graph)
1. **Precursor presence** P_p ∈ {0, 1}, one per library precursor. Prior π, estimated from the data: the target-decoy
   mixture of the current score, as the EPIFANY prior is set from the search.
2. **Peak choice** C_p ∈ {none, w1..wK}, the candidate peak groups of p (the up to 5 xic peaks). P_p = 0 forces none.
   The evidence of a peak group is the scorer's calibrated likelihood ratio (the D11FDWT logit → LR via the decoy
   null), so the graph uses what the scorer already knows.
3. **Fragment signals** S_j: a physical XIC peak in one isolation window, i.e. one fragment m/z bin (fitted ppm
   tolerance) × one apex scan range. Peak group (p, w) connects to S_j when one of p's fragments matches S_j at w's
   apex.

The likelihood of a shared signal is the explaining-away part.
- An intensity-aware noisy-OR: S_j's observed area must be explained by the sum of library-scaled abundances of the
  active parents plus background.
- If p and q both claim S_j but S_j's intensity fits one abundance, the pair is penalised. The same holds for
  fragment patterns where only p's other fragments are present.
- This is where the fragment-deficit signal (featD) enters properly. It is no longer a per-precursor feature but a
  joint fit.

## Inference
- Connected components: precursors in the same isolation window whose candidate peaks overlap in RT and share at
  least one fragment signal. Most components are small, and singletons pass through unchanged.
- Loopy belief propagation with the Evergreen library already in OpenMS (`MessagePasserFactory`,
  `BayesianProteinInferenceAlgorithm`). Exact enumeration is used for small components.
- Model parameters (noise floor, α/β of the noisy-OR, intensity tolerance) are estimated per run by EM over all
  components. No fixed constants.

## Output and FDR
- Posterior P(C_p = w) per peak group. A precursor's winner and its score come from the posterior instead of the
  isolated score.
- Decoys go through the same graph. Decoys are mass-shifted, so they rarely share signals, and their posteriors stay
  close to their isolated scores.
- Risk: this asymmetry is exactly what made decoys miss the entrapment null. Shared-signal penalties then mainly hit
  targets, which is the intended effect. FDR is still estimated by TDC on the posterior scores, and the entrapment
  rule checks it.

## Measurable first steps
1. After PLAN_PAIRLIB, on its search output (no labels):
   - the size distribution of components;
   - the share of target winners at q ≤ 0.05 in non-singleton components.
   Then, with labels and descriptively, the share of entrapment winners whose component contains their parent.
2. A Python prototype on those components, with a frozen plan and the same rules, as a new arm.
3. Only if it passes: a C++ port into OpenDIA, next to the picker, using Evergreen.

## Step 1 result on the old 100K library (xicids raw01, D11FDWT; `interf_pre/COMPONENTS_xicids_raw01.txt`)
Script: `scripts/interference_components.py`. Edges join candidate windows of different precursors in the same
isolation window whose peak ranges overlap, weighted by shared library fragments within the fitted 7.1 ppm tolerance.
Run time is 2:49 and peak memory 2.2 GB.

- **Graph size.** 2.75 M windows give 4.3 M edges.
  - 95 % of the edges share a single fragment. Mostly common low-mass ions, and only weakly informative.
  - At ≥ 2 shared fragments, 11 % of windows sit in non-singleton components. Components stay small (p99 2
    precursors, max 110), so exact or loopy inference is cheap.
- **Winners that share ≥ 2 fragments with another precursor's overlapping window** (labels used descriptively):

  | class | q ≤ 0.01 | q ≤ 0.05 |
  |---|---|---|
  | human | 12 % (of 1984) | 13 % (of 2359) |
  | entrapment | 39 % (of 23) | 32 % (of 109) |
  | decoy | 58 % (of 19) | 40 % (of 122) |

- **What this means for the design:**
  1. Shared signal is strong label-free evidence of a false winner, for both entrapments and decoys. Decoys use it
     even more, so penalising it should not by itself skew TDC calibration against targets.
  2. In under 5 % of cases is the claimant a stronger winner of its own. A DIA-NN-style greedy rule (drop the weaker
     of two precursors that explain the same signal) would therefore remove almost none of these false winners: they
     beat the claimants.
     - The model must explain away by **fit**: does the shared fragment's intensity, together with the window's other
       fragments, fit p's library pattern better than the claimant's? It cannot use score order alone.
     - This is the intensity-aware noisy-OR above. featD captures the single-precursor half of it.
  3. On this library, 5 of 23 entrapment winners at q ≤ 0.01 have a same-composition claimant, which is their
     parent. The rest share signal with unrelated precursors.
     - The PAIRLIB rerun will show how much this changes when every parent is present.

