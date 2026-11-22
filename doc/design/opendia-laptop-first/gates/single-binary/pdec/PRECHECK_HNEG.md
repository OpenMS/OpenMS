# PRECHECK_HNEG: a hard-negative GBM final stage in C++ (22 Nov 2026)

**Why.** RESULT_GAP2 branch 1 applies. On the sibling default, the label-free hard-decoy GBM (Bt) gains D +12.6 % and H +7.4 % at r7 N @ FDP_c 1 %. Plain label-free GBM (B) gains +3.9 % and +2.8 %. The network scorer cannot use hard negatives (RESULT_BTCPP), so this precheck ports the Python Bt procedure as a C++ final stage after the network and sibling rescoring.

**Code.** An opt-in `-gbm_rescore <frac>`, from a new commit with the commit and container recorded in the result. It runs 3 relabel rounds. In each round:
- **Positives:** the target winners at pooled q ≤ 0.01.
- **Negatives:** all windows of the top `<frac>` decoy precursors by best-window score. With `<frac>` = 0, all decoy winners are used instead.
- **Model:** a histogram GBM is cross-fitted by the library fold, using the same 49 features: 300 rounds, 31 leaves, learning rate 0.05, balanced weights, and early stopping on a 10 % validation share.
- **Score:** each window is scored by its held-out fold's model. The best window per precursor is calibrated per fold, as now, and q comes from the same paired TDC.

## Setup (fixed before any result)

- **Library:** r ≈ 1.56. Runs D and H, seeds 0 and 1.
- **Arms:**
  - **cur:** the new container with no flag. It must be byte-identical to `sibcpp-20261121/runs/<R>_sib_s<seed>` (the sibling default) on all 4 runs, or the precheck stops.
  - **gbt:** `-gbm_rescore 0.02`, the primary arm.
  - **gbm:** `-gbm_rescore 0`, a control for how much comes from the hard negatives.
- **Evaluator:** frozen `eval_pdec.py`.

## Decision rule (written before any result, for gbt)

- **Pass** if all hold. A pass goes to a frozen r7 plan PLAN_GBT (D seeds 0-4, +3 % bar, H after the D verdict).
  - mean N @ FDP_c 1 % gain ≥ +3 % on D and on H;
  - every seed gains;
  - mean lower bound up by ≤ +0.2 pp on each run;
  - median wall ≤ 1.5× cur.
- **Fail** otherwise. The gbm arm is reported but decides nothing.
