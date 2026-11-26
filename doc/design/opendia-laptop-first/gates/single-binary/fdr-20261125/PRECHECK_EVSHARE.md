# PRECHECK_EVSHARE: does `-evidence_share` (A) keep or gain IDs without loosening FDR on held-out data? (frozen 26 Nov 2026)

**Why.** DEV_AB screened A and B on the 10 Astral dev runs (raw 01/02 x seeds 0-4):
- A raised pooled N @ FDP_c 1 % from 20,871 to 21,603 (1.035x) at an unchanged combined bound (1.13 to 1.12 %).
- A did **not** bring the bound under 1 %. So this precheck does not test an FDR fix. It tests a narrower claim: on data not used to build A, A gains IDs or keeps them, and it does not loosen the entrapment bound.
- A was never run on diaPASEF.

The dev runs were used to design and compare A, so they are not reused here.

**Written** before any remote run with A, before any held-out result, and while the new raws are still downloading (never opened).
- **Seen:** the dev numbers (DEV_AB), and the G0 outcomes on raw 03, HeLa D and HeLa H (RESULT_ASTRAL_LAPTOP, RESULT_GBT_D/H, RESULT_BLEND).
- **Not seen:** any A output on raw 03, D, H, N1 or N2. No output of any kind on N1 or N2.

## Sets

| Set | Data | Status |
|---|---|---|
| raw03 | Astral `LFQ_Astral_DIA_15min_50ng_Human_03.raw` | held out from A (G0 seen before) |
| D | HeLa diaPASEF 2737.d (200SPD, 200 ng) | held out from A (G0 seen before) |
| H | HeLa diaPASEF 22474.d (PXD027359, 200 ng, 21 min) | held out from A (G0 seen before) |
| **N1** | PXD027359 `20210510_TIMS03_EVO03_PaSk_SA_HeLa_50ng_5_6min_DIA_high_speed_S1-B2_1_25186.d` (50 ng, 5.6 min) | **new**, never run |
| **N2** | PXD027359 `20210507_TIMS01_LC01_PaSk_SA_HeLa_400ng_diaPASEF_high_speed_A1_1_7546.d` (400 ng, other instrument and LC) | **new**, never run |

- N1 and N2 are downloaded from the PRIDE FTP into `evshare-20261126/data` and checked against the PRIDE SHA-1 (`pride.sha1`) before every submission.
- All four HeLa sets use the same library, `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`, with the same homolog-free foreign entrapments (r about 7). N1 and N2 are scored against that library as it is. Nothing is rebuilt for them.

## Arms (fixed now)

| Arm | Flags |
|---|---|
| G0 | none (current default) |
| **A** | `-evidence_share` |

No other arm is run or reported in this precheck. AB2 and the B variants stay dev-only.

## Runs (fixed before any result; every seed counts, none replaced)

- **Container:** `opendia-evshare-20261126.sif`, built from the pushed, clean `4c6ca367fd` by `evshare/build_container.sh`. The local, Docker and SIF OpenDIA sha256 must be equal, and the SIF `--helphelp` must list `-evidence_share`.
  - Built before this freeze (the build produces no result): SIF sha256 `c43320a6e0d350efa29d80db2cafe2470bd026edd5a4fbafe6eba6f5792f191a`; OpenDIA sha256 `733516ac5cbb222834fbebd6bad160146ddb51f55a25df493f0b720620798083` locally, in Docker and in the SIF. This is the same binary that produced the DEV_AB runs (`devAB/binary.sha256`). Evidence in `evshare/container-evidence/`.
  - The first build attempt stopped at its own help check (it grepped `--help`, where advanced flags are hidden). The check was changed to `--helphelp`; nothing else changed.
- **Work dir:** `/nfs/research/juan/singj/opendia/evshare-20261126/` (new). Nothing outside it is written.
- **Settings:** as in PRECHECK_BLEND.
  - Astral: 16 threads, `-lib gap-20261114/home/lib/library.pqp`, `-scoring_lib rtwin2-20261116/lib/library12_full.pqp` (`evshare/evs_astral.sbatch`).
  - HeLa: 8 threads, the r5 library (`evshare/evs_hela.sbatch`).
- **Seeds:** 0-4 for every set and arm.
- **G0 references:**
  - raw03 seed S: `astral-gbmfix-20261123/runs/raw03_sS`;
  - D and H seed S: `entlib-gbt-20261122/runs/{D,H}_gbt_sS`;
  - N1, N2: fresh G0 runs from this container (no reference exists).
- **Jobs** (`evshare/evs_submit.sh`):
  - gate: 3 no-flag cells, raw03_s0, D_s0 and H_s0;
  - main: 25 A cells (5 sets x 5 seeds) and 10 fresh G0 cells (N1, N2 x 5 seeds);
  - 38 jobs in all, plus 12 fresh G0 cells only if the gate rule below requires them.
- **Resource failures:** a job killed by the scheduler for memory or time is resubmitted once for the same cell, with twice the memory or time. This keeps the cell, it does not replace a seed. Any other failure is kept and counts as a failed run (N = 0). Every resubmission is recorded.

## Gate

The three no-flag runs must be byte-identical (all files) to their G0 references.
- **If all three match:** the existing G0 outputs are the reference for raw03, D and H.
- **If any differs:** the default has drifted. The existing references are not used. The 12 remaining G0 cells on raw03, D and H are run fresh from this container (`evs_submit.sh freshg0`) and become the reference. This is fixed now. The difference is reported.

## Evaluation

- **Evaluation script:** `evshare/evs_eval.sh`, written with this precheck.
- **Evaluator:** `stepA/eval_pdec_relabel.py` (sha256 `759ebe69…`), the PLAN_PDEC evaluator with the instr label fix. It runs with the corrected labels:
  - Astral: `RELABEL=stepA/label_audit_astral.tsv` and the arm-E library (`armE_raw01/library.pqp`, sha256 `d9e3fbe0…`);
  - HeLa: `RELABEL=stepA/label_audit_hela.tsv` and the r5 library (sha256 `76f2a6c2…`).
- **Secondary:** the same runs evaluated without RELABEL (`stepA/eval_pdec_instr.py`, sha256 `f85096fa…`) are reported but not ruled on.
- **Layout:** `<set>_<arm>_s<seed>` with arms G0 and A.
- **Metrics:** per set and arm, at the tool's own q <= 0.01:
  - pooled lower bound E/N and combined bound E(1 + 1/r)/N over 5 seeds;
  - median N @ FDP_c 1 % over 5 seeds (the ID metric);
  - pooled N @ FDP_c 1 % (sum over seeds), reported.

## Decision rule (written before any held-out result)

For each of the five sets raw03, D, H, N1 and N2:

1. **No ID loss:** A median N @ FDP_c 1 % >= 1.00 x G0.
2. **No FDR loosening:** A pooled combined bound <= G0 pooled combined bound + 0.10 pp.

**Outcomes:**
- **A passes** if both criteria hold on all five sets. Then a frozen plan to make A the default would follow, byte-checked both ways as with PLAN_GBT. It would be filed as an ID gain at unchanged FDR control, **not** as an FDR fix.
- **A fails** if any criterion fails on any set. A then stays opt-in, with every set reported. Nothing is retuned: no other claimant rule, threshold or feature set is tried on these runs. Any next A variant needs a new held-out set.
- **Instrument split:** if A passes on raw03 but fails on a HeLa set, or the reverse, the result is reported as instrument-dependent. It is not fixed by a rule chosen after seeing it.
- **Unusable new set:** if all five G0 runs of N1 or N2 fail (non-zero exit after the one resource resubmission), that set is reported as unusable and the rule is applied to the other four. If only some G0 or A runs fail, they count as N = 0, as above.

**Reported, not ruled on:**
- the lower bound;
- the worst single-run combined bound per set;
- own-q N;
- the pair fairness P;
- the secondary no-relabel evaluation;
- wall time and peak RSS (descriptive; different nodes).

**Honesty notes:**
- With r about 7 on HeLa, the combined bound sits near the lower bound, so criterion 2 is close to a check on E/N there.
- On Astral, 5 runs give about 60 entrapments in pool. The 0.10 pp margin is smaller than the Poisson noise of the bound. Criterion 2 can fail by chance, and a pass is weak evidence of no loosening.
- G0 seed-to-seed spread on Astral is large (dev single-run N @ FDP_c 1 % 1,728 to 2,270). The median over 5 seeds is used for that reason.

## Cost

- Astral: about 5-15 min per run at 16 threads; A adds one resolver pass. Dev peak RSS was 13.3-14.4 GB; the jobs ask for 48 GB.
- HeLa D and H: about 15-25 min per run at 8 threads, with G0 peak RSS up to 32 GB. The jobs ask for 200 GB, because A keeps the resolver inputs alive and N2 is a larger raw (15.5 GB zipped) never run before.
