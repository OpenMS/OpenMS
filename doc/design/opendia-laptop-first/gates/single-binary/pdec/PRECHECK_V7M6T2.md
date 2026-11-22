# PRECHECK_V7M6T2: V7M6 re-run with the faster interference resolver (22 Nov 2026)

**Why.** RESULT_V7M6T failed rule B only on D r7 wall time: 1.504× against 1.5×. The interference step was 62 % of V7M6's extra time and grew 2.44× for 1.39× windows. With the user's go-ahead, the resolver was made faster without changing its output: OpenMS `955516d84c`, container `opendia-pd10-fastint.sif`, sha256 `8de3a968`.
- **Resolver change:** window tests run on the index entries, the fragment test runs once per claimant precursor (not once per window), and claimant data is packed in RT order.
- **Bench:** on the dumped D r7 V7M6 first-pass inputs (14.1 M windows, 8 threads), the resolver went from 1,562 s to 878 s, with byte-identical output.
- **Tests:** 25 of 25 unit tests pass, including the all-pairs brute force and the thread-independence test.

## Setup (unchanged from PRECHECK_V7M6T, fixed before any result)

Everything in PRECHECK_V7M6T stays: cells, arms, whole-node Ice Lake jobs in ABBA order, the wall-ratio definition, the IDs from rep 1, the rep-2 determinism check, the evaluator `precheck_v7m6t.py`, and all four rule B bars.

The only changes are the container (`opendia-pd10-fastint.sif`) and the directory (`v7m6t2-20261122`, submit script `v7m6t2_submit.sh`).

**Byte-identity precondition.** Before the cells run, `fastint-20261122` must show the new container byte-identical on three runs. If any differs, this precheck stops.
- D r≈1.56 default vs `gbmdef-20261122/D_def_s0`;
- H r7 default vs `entlib-gbt-20261122/runs/H_gbt_s0`;
- D r7 V7M6 vs `v7m6t-20261122/D_r7_V7M6_1`.

Since the IDs are then identical to RESULT_V7M6T by construction, this re-run decides on wall time. The ID bars are still evaluated.

## Decision rule (rule B, bars unchanged)

V7M6 **passes** if all hold: r0 HeLa ≥ +3 % on D and H, r7 N @ FDP_c 1 % ≥ +2 % on D and H, r7 lower bound up by ≤ +0.2 pp, and matched-node wall ratio ≤ 1.5× on every cell. A pass goes to a frozen 5-seed r7 plan (PLAN_V7M6: D seeds 0-4 first, H after the D verdict). Anything else is a fail.
