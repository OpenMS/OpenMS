# PRECHECK_V7M6T: V7M6 wall time re-measured on matched nodes (22 Nov 2026)

**Why.** In RESULT_GAP2, V7M6 (`-valley_fraction 0.7 -max_peaks 6`) passed every rule B bar except wall time on one cell: D r7 at 1.52× against 1.5×. The two arms of that cell ran on different nodes (`hl-codon-cn-009` and `-010`), and the same arm was 0.79× on H r7, so node noise is the likely cause. The user asked to re-time it on matched nodes. The bars are not changed.

**What changed since GAP2.** The default now has the hard-negative GBM stage (PLAN_GBT) and the cleanup (OpenMS `9ba060f978`, container `opendia-pd10-best.sif`, sha256 `1bfd94e4`, byte-identical to the GBM default). GAP2 measured V7M6's ID gains on the sibling default before the GBM. So this precheck applies **all of rule B again, on the current default**, with the same bars. That is stricter than re-timing alone, not looser.

## Setup (fixed before any result)

- **Cells:** D and H × r0 and r7, the GAP2 libraries (r0 `libsize-20261121/r0/lib.pqp`, r7 `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`). Seed 0, 8 threads, no `-features_out`.
- **Arms:** def (no flag) and V7M6 (`-valley_fraction 0.7 -max_peaks 6`), both from `opendia-pd10-best.sif`.
- **Matched nodes:** one job per cell on an Intel Ice Lake node (`-C icelake`), holding all 48 of its 48 cores so no other job shares it (codon refuses `--exclusive`; amended before any job ran). The tool still runs with 8 threads. The job first reads the raw data and library once to warm the cache, then runs def, V7M6, V7M6, def in that order (ABBA), each under `/usr/bin/time -v`.
- **Wall ratio per cell:** (V7M6 rep 1 + V7M6 rep 2) / (def rep 1 + def rep 2).
- **IDs:** from rep 1 of each arm, with the GAP2 part B definitions: HeLa IDs at own q ≤ 0.01, and at r7 N @ FDP_c 1 % and the lower bound. Rep 2 must be byte-identical to rep 1 (a determinism check, reported).
- **Sanity check:** H r7 def rep 1 must be byte-identical to `entlib-gbt-20261122/runs/H_gbt_s0` (the current default on that cell). If it is not, the precheck stops.
- **Evaluator:** `precheck_v7m6t.py`, committed with this file.

## Decision rule (rule B, bars unchanged from PRECHECK_CANDREC / GAP2)

V7M6 **passes** if all hold on the current default:
- r0 HeLa ≥ +3 % on D and H;
- r7 N @ FDP_c 1 % ≥ +2 % on D and H;
- r7 lower bound rises by ≤ +0.2 pp on D and H;
- matched-node wall ratio ≤ 1.5× on every cell.

A pass goes to a frozen 5-seed r7 plan (PLAN_V7M6: D seeds 0-4 first, H after the D verdict). A fail on any bar is a fail and is reported as such, with each bar's value.
