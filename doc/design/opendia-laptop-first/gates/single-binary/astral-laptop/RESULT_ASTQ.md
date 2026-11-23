# RESULT_ASTQ: the GBM rescore is what makes own q ≤ 0.01 loose on Astral (23 Nov 2026, diagnostic)

The plan is `PRECHECK_ASTQ.md`. It was frozen and pushed in docs `8d5868da49` before any job ran (`PRECHECK_ASTQ.sha256`).

- **Runs:** codon `astq-20261123/runs`, 52 jobs. All exited 0.
- **Container:** `opendia-pd10-gbmdef.sif` (OpenMS `bcb156e036`, sha256 `7dd8a1d6…`).
- **Evaluator:** `eval_pdec.py` (sha256 `bb03387d…`), unchanged, r = 1.0018. Full output: `EVAL_ASTQ_full.txt`. Cost: `COST_ASTQ.txt`.

**Anchor gate: passed.** A0 (`gbmdef`, no flags) is byte-identical to L0 (current default, `88e5da3023`) on raw01_s0 and raw02_s2. So the ladder measures the current code.

## Numbers (raws 01 + 02, 5 seeds each; pooled E and N at own q ≤ 0.01)

| Arm | Stages removed | E / N | Lower bound | **Combined bound** | Step | Pair fairness P | Median N @ FDP_c 1 % (raw01 / raw02) | Median wall |
|---|---|---:|---:|---:|---:|---:|---|---:|
| L0 | none (current default) | 118 / 20,121 | 0.59 % | **1.17 %** | | 0.473 | 2136 / 2146 | (batch 2) |
| L1 | GBM | 68 / 18,546 | 0.37 % | **0.73 %** | **-0.44 pp** | 0.496 | 1943 / 2285 | 6.3 min |
| L2 | + sibling | 87 / 18,669 | 0.47 % | 0.93 % | +0.20 pp | 0.512 | 1802 / 1947 | 5.5 min |
| L3 | + RT refinement | 82 / 18,755 | 0.44 % | 0.87 % | -0.06 pp | 0.531 | 1735 / 2006 | 4.5 min |
| L4 | + matched filter | 70 / 18,568 | 0.38 % | 0.75 % | -0.12 pp | 0.510 | 1848 / 1957 | 4.2 min |
| L5 | + 16 Nov trace, 1/K0 and evidence settings | 50 / 17,471 | 0.29 % | 0.57 % | -0.18 pp | 0.506 | 1804 / 1909 | 3.6 min |

The steps are the change in the combined bound when a stage is removed. Negative means removing it tightens the bound. Peak RSS is 17 GB in every arm.

Per raw, the combined bound with and without the GBM (L0 to L1) is 1.22 to 0.96 % on raw 01 and 1.11 to 0.52 % on raw 02. Without the GBM, 1 of 10 runs is above 1 % (raw01_s4, 1.10 %; its q ≤ 0.01 set is the same with and without the GBM). With the GBM, 6 of 10 runs are above 1 %. L5 reproduces the 16 Nov medians exactly (raw01 1804, raw02 1909).

## Attribution (frozen rule)

| Condition | Bar | GBM step (L0 to L1) | |
|---|---|---:|---|
| Step size | ≥ 0.20 pp | 0.44 pp | yes |
| Share of the total drop L0 to L5 (0.60 pp) | ≥ half | 0.73 | yes |
| Combined bound after the step | ≤ 1.10 % | 0.73 % | yes |

**Verdict: the GBM rescore is named the cause.** No other single step meets the size bar. Removing sibling rescoring loosens the bound instead (+0.20 pp).

**Fairness supporting evidence: not met.** The plan's supporting sign was "decoys easier than false targets after that stage", which means P rising with the stage. With the GBM, P is lower (0.473 against 0.496 without it). So the GBM does not make entrapment targets beat their own decoys more often. The looseness is not a pair asymmetry. It is consistent with the GBM ranking false targets above the decoy level as a group (more entrapments pass at the same decoy count).

## What the GBM gives and costs on Astral (reported, not ruled on)

- **IDs at true FDP 1 %:** raw 01 +9.9 % (1943 to 2136), but raw 02 -6.1 % (2285 to 2146). Without the GBM, raw 02 has more IDs at true FDP 1 %.
- **Seed instability is the GBM's too:** raw02_s4 has 1322 with the GBM and 2290 without it.
- **On HeLa** the GBM was +8.6 % (D) and +7.1 % (H) with bounds +0.01 / +0.02 pp on D (PLAN_GBT), but H's lower bound went from 1.15 to 1.29 % (RESULT_GBT_H). So the looseness may be the same effect, larger on Astral.

## Next (needs a decision; not run)

A fix to the GBM must keep the HeLa gain and bring Astral's own q back to ≤ 1 %. Candidate directions for a new precheck (frozen before any run):
- a single GBM round instead of three, since the logs show degenerate later rounds on Astral;
- a stricter positive set (for example q ≤ 0.005 instead of 0.01);
- dropping the GBM's q-values and keeping its ranking with the network's calibrated tail.

None of these changes the default without a frozen plan on D, then H, and an Astral check.
