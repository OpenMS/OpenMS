# PRECHECK_BLEND: does the network/GBM score blend keep IDs and restore own-q calibration? (frozen 25 Nov 2026)

**Why.** DIAG_FDR (local; Astral raw 01/02 × 5 seeds, the dev set) found the following:
- The GBM rescore's own-q excess on Astral is confined to q ≤ 0.01.
- It sits mostly in entrapments that share library fragments with an isobaric human peptide, a class the pseudo-reverse decoys cannot represent.
- An offline per-precursor blend of the network and GBM calibrated scores (w = 0.5, recalibrated, recompeted) brought the pooled bound from 1.17 to 0.75 % on the dev set, and median N @ FDP_c 1 % from 2,136 to 2,165.
- The opt-in `-gbm_blend <w>` (OpenMS `b3391c92e8`) reproduces the simulation exactly on raw01_s0.

The dev set was used to find the mechanism and to compare weights, so it cannot confirm the effect. This precheck is the test on data not used to choose it.

**Written** before any remote run with the blend, before any held-out result, and before the container exists.
- **Seen:** the dev-set offline numbers (DIAG_FDR), and the outcomes of the current default on all sets (RESULT_ASTRAL_LAPTOP, RESULT_GBT_D/H, RESULT_GBMNULL).
- **Not seen:** any blend output on raw 03, HeLa D or HeLa H.

## Arm (one, fixed now)

| Arm | Flags | Reference (existing outputs, read only) |
|---|---|---|
| **B5** | `-gbm_blend 0.5` | G0 = the current default |

- **Weight:** w = 0.5 is the equal weight. No other weight is run or reported in this precheck.
- **G0 references** (all produced by binaries whose default is byte-identical to `b3391c92e8` with no flag; gate below):
  - Astral raw NN seed S: `astral-gbmfix-20261123/runs/rawNN_sS` (RESULT_ASTRAL_LAPTOP / GBMFIX G0);
  - HeLa D r ≈ 7 seed S: `entlib-gbt-20261122/runs/D_gbt_sS`;
  - HeLa H r ≈ 7 seed S: `entlib-gbt-20261122/runs/H_gbt_sS`.

## Runs (fixed before any result; every seed counts, none replaced)

- **Container:** `opendia-blend-20261125.sif`, built from the pushed, clean `b3391c92e8`. The local, Docker and SIF OpenDIA sha256 must be equal.
- **Work dir:** `/nfs/research/juan/singj/opendia/blend-20261125/` (new). Nothing outside it is written.
- **Astral:** raws 01, 02 (dev) and **03 (held out)**, seeds 0-4. The settings are as in `gbmfix_astral.sbatch`:
  - 16 threads, 32 GB;
  - `-lib gap-20261114/home/lib/library.pqp`, `-scoring_lib rtwin2-20261116/lib/library12_full.pqp`.
- **HeLa D r ≈ 7 (held out):** 2737.d with `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`, seeds 0-4, 8 threads, 160 GB, as in PLAN_GBT.
- **HeLa H r ≈ 7 (held out):** 22474.d, same library, seeds 0-4, 8 threads, 160 GB, as in PLAN_GBT.
- **Gate runs** (the same container with no flag must equal G0):
  - Astral raw01_s0;
  - HeLa D r ≈ 7 s0;
  - HeLa H r ≈ 7 s0.
- **Total:** 15 Astral + 5 D + 5 H blend jobs, plus 3 gate jobs, so 28 jobs.

## Evaluation (unchanged evaluator)

- **Evaluator:** `pdec/eval_pdec.py` (sha256 `bb03387d…`), run on a layout `<set>_<arm>_s<seed>` per set, with arms G0 and B5.
- **Libraries:**
  - Astral: the arm-E library (r = 1.0018);
  - HeLa: the r ≈ 7 library (r = 7.0539).
- **Metrics:** per arm and set, at the tool's own q ≤ 0.01:
  - the pooled lower bound E/N;
  - the pooled combined bound E(1 + 1/r)/N;
  - the worst single-run combined bound;
  - the median N @ FDP_c 1 %.

## Gate

The three no-flag runs must be byte-identical (`precursors.tsv` and the per-run tables) to their G0 references.
- **If all three match:** the existing G0 outputs are the reference.
- **If any differs:** the default has drifted since those references, so they are not used. Fresh G0 runs (no flag, the same container) are submitted for all 25 cells, and they become the reference. This is fixed now and is not a choice made after the fact. The difference is reported.

## Decision rule (written before any held-out result)

B5 **passes** if all of these hold:

1. **Astral held out (raw 03):**
   - the B5 pooled combined bound is ≤ 1.00 %;
   - the B5 median N @ FDP_c 1 % is ≥ 0.98 × G0.
2. **Astral all three raws:**
   - the B5 pooled combined bound over 15 runs is ≤ 1.00 %;
   - the B5 median N @ FDP_c 1 % on each raw is ≥ 0.98 × G0.
3. **HeLa D and H (held out), each:**
   - the B5 median N @ FDP_c 1 % is ≥ 0.98 × G0;
   - the B5 pooled lower bound is ≤ G0 lower + 0.10 pp. The blend must not loosen HeLa.

**Reported, not ruled on:**
- the pooled combined bound on HeLa (r ≈ 7, so it sits close to the lower bound);
- the worst single run per set;
- own-q N;
- the pair fairness P;
- wall time and RSS (descriptive; different nodes).

**Consequences:**
- **Pass:** a frozen plan to make B5 the default. It would be byte-checked both ways, as with PLAN_GBT.
- **Fail:** B5 stays opt-in, with every set reported. Nothing is retuned: no other w and no other combination is tried on these runs. Any next variant needs a new held-out dataset.
- **HeLa-only failure:** if B5 passes on Astral but fails on HeLa, the result is reported as instrument-dependent. It is not "fixed" by a rule chosen after seeing it.

**Honesty note:** criterion 2 includes the dev raws. A pass there is expected and is not evidence by itself. Criteria 1 and 3 carry the test.

## Cost

28 jobs:
- Astral: about 7-15 min each at 16 threads.
- HeLa D/H: about 20-40 min each at 8 threads.
- The blend itself adds well under a second per run (one calibrateTail and one TDC).
