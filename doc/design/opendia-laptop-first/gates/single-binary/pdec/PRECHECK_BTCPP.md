# PRECHECK_BTCPP: C++ hard-decoy self-training, default byte check and two-seed r ≈ 1.56 check (21 Nov 2026)

**Code:** OpenMS `f3005a4b14` adds the opt-in `-hard_decoys <fraction>` and `-scorer_self_q <list>`. It passes 26/26 ctests, including a new unit test showing that fraction 0 equals the option unset. The container is `opendia-pd10-btopt.sif` (sha256 `a3a75146…`).

**What changes.** In each `CrossFitScorer` self-training round, the negatives become all training-fold windows of the hardest `fraction` of decoy precursor winners (at least 10). By default they are the winners of all decoy precursors. The learner is still the network with the decoy-whitened readout, not the GBM used in RESULT_BT156, so the gain may not carry over.

## Setup (fixed before any result)

- **Library:** r ≈ 1.56 (`/homes/sing/hela_entrap_codon.predicted_decoys.pqp`). Runs D and H, seeds 0 and 1. Run dirs are `$C/btcpp-20261121/runs/<R>_<arm>_s<seed>`.
- **Arms:**
  - **cur:** no new flags.
  - **bt1:** `-hard_decoys 0.02`, keeping the default self_q 0.15 and 0.05.
  - **bt2:** `-hard_decoys 0.02 -scorer_self_q "0.01 0.01 0.01"`. This is the design of `precheck_ceil0.py` arm Bt.
- **Byte check:** cur `precursors.tsv` must be byte-identical to `bt156-20261121/<R>_s<seed>` (container rtrdef, same library and seed) for all four runs. If any differs, stop.
- **Evaluator:** frozen `eval_pdec.py` (sha256 `bb03387d…`). The measure is N at combined FDP 1 % with the lower and combined bounds.

## Decision rule (written before any result)

Take the better of bt1 and bt2 by mean D gain.

- **Pass** if all of these hold:
  - the mean gain on D at N@FDP_c 1 % is ≥ +2 %, and both D seeds gain;
  - the mean gain on H is ≥ 0;
  - the mean lower bound rises by ≤ +0.2 pp on D and on H.

  The arm then goes to a frozen r7 plan (PLAN_BT, D seeds 0-4, a +3 % bar, H after the D verdict).
- **Fail** otherwise. The flag stays opt-in. Possible next steps would be a GBM-like learner, or using hard negatives only in the last round. Neither is part of this precheck.
