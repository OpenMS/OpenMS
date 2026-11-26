# RESULT_EVSHARE: `-evidence_share` (A) on held-out data. FAIL under the frozen rule (26 Nov 2026)

- **Precheck:** PRECHECK_EVSHARE, frozen at docs `c8ba28183d` before any job. Submission record: SUBMITTED_EVSHARE.
- **Container:** `opendia-evshare-20261126.sif` (sha256 `c43320a6…`; OpenDIA `733516ac…`, the DEV_AB binary).
- **Gate: passed.** No-flag raw03_s0, D_s0 and H_s0 are byte-identical to the G0 references, 9/9 files (`evshare/gate.txt`). So G0 for raw03, D and H is the existing default output, as frozen.
- **Evaluation:** `evshare/evs_eval.sh` with the frozen evaluators.
  - Primary: relabelled, `eval_pdec_relabel.py`. Secondary: instr only, `eval_pdec_instr.py`. Run tables were fetched with a sha256 check (`evshare/remote.sha256`).
  - **Procedural note:** the script aborted at set N1, because N1 has no run tables (the evaluator raises on an empty layout, `evs_eval_run.log`). The remaining calls (N2 primary, all secondary) were then run by hand with the same evaluators, libraries and layouts. No rule, metric or input changed. Outputs: `evshare/EVAL_EVSHARE_stage1.txt` and `evshare/EVAL_EVSHARE_stage2.txt`.

## N1 is unusable

- All 10 N1 runs (5 G0, 5 A) exited with code 8 after about 32 s: `Calibration found only 22 whole-run stage-1 anchors at q <= 0.05 of 40000 panel targets`.
- N1 is the 5.6 min, 50 ng run (3,176 frames). The default fails on it as well, so this is a calibration limitation of OpenDIA on very short gradients, not an effect of A.
- Under the frozen rule (all five G0 runs failed), N1 is reported as unusable and the rule is applied to the other four sets.

## Numbers (primary, relabelled; own q <= 0.01; 5 seeds per set)

| Set | G0 combined | A combined | A - G0 | G0 median N@FDP_c 1 % | A median N@FDP_c 1 % | A / G0 | Criterion 1 (>= 1.00x) | Criterion 2 (<= +0.10 pp) |
|---|---:|---:|---:|---:|---:|---:|---|---|
| Astral raw03 | 1.21 % | 1.19 % | -0.02 pp | 2,246 | 2,258 | 1.005 | pass | pass |
| HeLa D | 1.07 % | 1.04 % | -0.03 pp | 21,805 | 22,270 | 1.021 | pass | pass |
| HeLa H | 1.28 % | 1.26 % | -0.02 pp | 37,897 | 38,031 | 1.004 | pass | pass |
| **HeLa N2 (new)** | 1.27 % | 1.27 % | +0.00 pp | 45,771 | 45,755 | **0.9997** | **fail** | pass |
| HeLa N1 (new) | n/a | n/a | | n/a | n/a | | unusable | unusable |

**Verdict: FAIL.** A misses criterion 1 on N2 by 16 precursors out of 45,771 (0.035 %). As frozen, A stays opt-in and nothing is retuned.

The miss is well inside the seed-to-seed noise: N2 single runs span 45,657-45,765 for A and 45,694-45,846 for G0. So the honest reading is "no measurable ID effect on N2", not "A loses IDs on N2". The rule was written as >= 1.00x and is applied as written.

## Reported, not ruled on

- **Lower bound:**

  | | raw03 | D | H | N2 |
  |---|---|---|---|---|
  | G0 | 0.60 % | 0.94 % | 1.12 % | 1.11 % |
  | A | 0.59 % | 0.91 % | 1.10 % | 1.11 % |

  The HeLa 90 % intervals overlap between arms on every set.
- **Worst single-run combined bound:**
  - raw03: A 1.79 % (s4), G0 1.46 %;
  - D: A 1.11 %, G0 1.12 %;
  - N2: A 1.32 %, G0 1.31 %.
- **raw03_A_s4 is an outlier:** 1,456 IDs at own q (G0 seeds give 2,261-2,343) and 692 at FDP_c 1 %. Its log shows a third GBM round that was kept, while most other runs stopped at round 2. This is a GBM round-selection instability made visible by A, not a crash. Because the median is used, it does not decide criterion 1, but it is a real per-run risk of A on Astral.
- **Own-q N:**
  - raw03: A 10,314 vs G0 11,447 (with s4);
  - D: 111,675 vs 109,910;
  - H: 192,712 vs 192,239;
  - N2: 231,036 vs 230,993.
- **Pair fairness P:** unchanged within ±0.006 on every set.
- **Secondary (instr, no relabel):** same direction everywhere.

  | Median N@FDP_c 1 % | raw03 | D | H | N2 |
  |---|---:|---:|---:|---:|
  | G0 | 2,246 | 21,495 | 37,376 | 45,304 |
  | A | 2,258 | 21,847 | 37,648 | 45,276 |
  | A / G0 | 1.005 | 1.016 | 1.007 | 0.999 |

  Combined bounds move by -0.04 to +0.01 pp.
- **Cost (descriptive; different nodes; `evshare/COST_EVSHARE.txt`):**
  - A wall was 5.5-8.7 min on raw03, 15.5-16.7 min on D, 19.5-22.6 min on H and 30.5-52.9 min on N2. G0 was 30.5-32.7 min on N2.
  - Peak RSS was about equal between arms: raw03 12.4-13.1 GB, D 31 GB, H 23 GB, N2 60.4 GB.

## Reading

- **A does not fix FDR control on any held-out set.** The combined bound moves by -0.03 to 0.00 pp. HeLa H and N2 stay at about 1.27 %, and the H and N2 lower bounds stay above 1 %. This matches the dev screen.
- **ID gains seen on dev do not generalise as a robust gain.** They are +2.1 % on D, +0.4-0.5 % on raw03 and H, and zero on N2.
- **New sets:** a short-gradient raw (N1) is outside OpenDIA's calibration range, and a long 400 ng TIMS01 raw (N2) behaves like H, with a combined bound of about 1.27 % for both arms.
- **Next:** this result, together with CALSEARCH, leaves no candidate that meets the goal on both instruments. The remaining lead is a label-free feature or competition step aimed at the 4+ shared-fragment class on diaPASEF (DIAG_STEPA, CALSEARCH). Any variant needs a new held-out set: N2 is now used, and N1 cannot be used.
