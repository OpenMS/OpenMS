# RESULT_GBMFIX: neither GBM variant fixes Astral's own q (23 Nov 2026, FAIL)

The plan is `PRECHECK_GBMFIX.md`. It was frozen and pushed in docs `f1a926003e` before any build or run.

- **Code:** OpenMS `17430ba8eb` (`-gbm_rounds`, `-gbm_positive_q`, opt-in).
- **Container:** `opendia-pd10-gbmvar.sif`, sha256 `a39e1f62…`.
- **Runs:** codon `gbmfix2-20261123/runs`, 27 jobs. All exited 0.
- **Evaluator:** `pdec/eval_pdec.py` (sha256 `bb03387d…`), unchanged. Output is in `EVAL_GBMFIX_astral.txt` and `EVAL_GBMFIX_helaD.txt`.

**Gate: passed.** With no flags, the new container is byte-identical to the current default on Astral raw01_s0 and raw02_s2, and on HeLa D r ≈ 1.56 (3/3 files; `GATE_GBMFIX.txt`). The defaults are unchanged.

## Astral (raws 01 + 02, 5 seeds each)

| Arm | Pooled combined bound | Worst run | Median N @ FDP_c 1 % raw01 / raw02 | vs G0 | Crit 1 (≤ 1.00 %, no run > 1.50 %) | Crit 2 (≥ 0.97 × G0 each) |
|---|---:|---:|---|---|---|---|
| G0 (current) | 1.17 % | 1.77 % | 2136 / 2146 | | | |
| G1 (`-gbm_rounds 1`) | 1.12 % | 1.77 % | 2121 / 2146 | 0.993 / 1.000 | **fail** | pass |
| GQ (`-gbm_positive_q 0.005`) | 0.81 % | 2.94 % | 1913 / 2174 | 0.896 / 1.013 | **fail** | **fail** |
| NG (no GBM, reference) | 0.73 % | 1.10 % | 1943 / 2285 | 0.910 / 1.065 | | |

## HeLa D r ≈ 7 (seeds 0, 1)

| Arm | Mean N @ FDP_c 1 % | vs G0 | Combined bound | Crit 3 (≥ 0.98 ×) | Crit 4 (≤ G0 + 0.10 pp) | Wall / RSS |
|---|---:|---:|---:|---|---|---|
| G0 | 21,242 | | 1.31 % | | | (PLAN_GBT) |
| G1 | 21,216 | 0.999 | 1.16 % | pass | pass | 17:15-17:32 / 32.6 GB |
| GQ | 20,460 | 0.963 | 1.32 % | **fail** | pass | 17:23-17:54 / 32.7 GB |

**Verdict: FAIL. Neither variant passes all four criteria, so both stay opt-in and the default is unchanged.**

## Reading

- **G1 (one round)** is harmless on HeLa. It keeps all the IDs and even tightens the combined bound (1.31 to 1.16 %). On Astral it changes almost nothing:
  - 6 of 10 Astral runs are byte-identical in N to G0, because in those runs rounds 2 and 3 were already discarded.
  - raw01_s0 stays at 1.77 %.
  - It does fix the low seed raw02_s4 (1322 to 2217).
  - So **the Astral looseness is already present after the first GBM round**. It is not caused by self-reinforcing later rounds.
- **GQ (cleaner positives)** breaks on Astral:
  - with only about 370-700 positives in round 1, the GBM collapses on some seeds (raw01_s1 420, raw02_s2 145 at FDP_c 1 %);
  - the worst run reaches 2.94 %.
  - On HeLa it costs 3.7 % of IDs.
- **Interpretation.** On Astral the first-round GBM itself ranks some entrapment (false) targets above the decoy level. With about 2,000 IDs and 17,000 hard negatives per run, the model has too few positives to generalise across folds. The HeLa runs have 18,000-38,000 positives and do not show this.

## Next (needs a decision; not run)

The candidates all address a GBM trained on too few positives:
1. **Gate the GBM on the number of positives.** Skip it, or fall back to the network score, when round-1 positives are below a floor. A floor of about 5,000 would keep every HeLa run (≥ 18,000) and skip every Astral run (≤ 2,300). On Astral this gives NG (0.73 %, raw 02 +6.5 %, raw 01 -9 %) and leaves HeLa unchanged. It is the smallest change, and it is byte-identical on HeLa by construction.
2. Score the GBM output through the network's calibrated tail instead of its own (a bigger change).

Option 1 would get its own precheck: Astral pooled bound ≤ 1.00 %, plus HeLa byte-identity.
