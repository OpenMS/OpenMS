# PRECHECK_SIBCPP: C++ charge-state sibling rescoring with the network scorer (21 Nov 2026)

**Why.** RESULT_JOINT2 showed that hard decoys plus siblings (Bts) give +20-26 % at FDP_c 1 % with a GBM. Siblings alone with all decoys (Bs) gave H +18 %, but D only −0.3 / +3.5 %. The C++ scorer cannot use hard decoys (RESULT_BTCPP), so this precheck tests siblings alone in C++. It also confirms that the two-stage scheme is valid.

**Fold check.** Folds are assigned by the unmodified-sequence family (`fnv1a(family) % 3`). On the r ≈ 1.56 library, 0 of 245,711 multi-precursor peptides span two folds. Siblings always share a fold, so the first-stage scores used for the sibling features are out-of-fold for the precursor. That settles the cross-fitting caveat in RESULT_JOINT2.

**Code:** an opt-in `-sibling_features`. The workflow scores once, adds `s_n`, `s_co`, `s_best` and `s_best_d` from the calibrated first scores, and scores again with the same options. Commit and container are recorded in the result.

## Setup (fixed before any result)

- **Library:** r ≈ 1.56. Runs D and H, seeds 0 and 1.
- **Arms:**
  - **cur:** the new container with no flag. It must be byte-identical to `bt156-20261121/<R>_s<seed>` on all 4 runs, or the precheck stops.
  - **sib:** `-sibling_features`.
- **Evaluator:** frozen `eval_pdec.py`.

## Decision rule (written before any result)

- **Pass** (frozen r7 plan PLAN_SIB, D seeds 0-4, +3 % bar, H after the D verdict) if all hold:
  - mean N@FDP_c 1 % gain ≥ +5 % on H;
  - mean gain ≥ 0 on D;
  - every H seed gains;
  - mean lower bound up by ≤ +0.2 pp on each run.
- **Partial:** if H passes but D loses, siblings stay opt-in and wait for the hard-negative scorer redesign, which D needs (RESULT_JOINT2: the D gain comes only from Bts).
- **Fail** otherwise.
