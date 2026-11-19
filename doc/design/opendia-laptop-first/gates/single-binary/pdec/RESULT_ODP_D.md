# RESULT_ODP_D: OpenDIA's built-in predicted decoys calibrate like DIALibGen's at r ≈ 7. PLAN_ODP fails criteria 1 and 2

- **Plan:** PLAN_ODP, frozen c3e71db6a9 at 2026-11-18T23:45:29Z.
- **Runs:** D seeds 0-4, CUR (`-decoys predict`) and NEW (`-decoys predict -pg_scores`). Codon `entlib-odp5-20261118`,
  container `opendia-pd2-01d5ded04a.sif`.
- **Library:** every run wrote the same predicted library (`hela_r5.predicted_decoys.pqp`, sha256 76f2a6c2…, 10/10
  identical).
- **Evaluation:** `eval_pdec.py` (sha256 bb03387d…) on the tool's own q. Output is in `EVAL_ODP_D.txt`.
- **Cost:** 1:45-2:28 wall and 89-92 GB peak per run.

## Criteria (D, pooled over 5 seeds, q ≤ 0.01)
| criterion | limit | CUR | NEW | verdict |
|---|---|---|---|---|
| 1 / 2 lower bound | ≤ 0.010 | **0.0106** (90 % 0.0098-0.0114) | **0.0106** (0.0099-0.0114) | **fail / fail** |
| 1 / 2 combined bound | ≤ 0.012 | **0.0121** (0.0112-0.0131) | **0.0121** (0.0113-0.0130) | **fail / fail** |
| 3 CUR median N@FDP_c 1 % | ≥ 6,235 | **8,233** (7,968-8,286) | (NEW 10,100) | **pass** |
| 4 coverage | ≥ 99 % | **99.85 %** (1,388,026 / 1,390,163) | same library | **pass** |

**Verdict: FAIL** (criteria 1 and 2, narrowly). Nothing is retuned.

**Reported, not criteria:**
- `-pg_scores` gives +22.7 % N@FDP_c 1 % (10,100 against 8,233).
- P(entrapment beats its own decoy | q ≤ 0.05) is 0.534 / 0.539.

## Comparison (D, 5 seeds, r = 7.05, same targets)
| q ≤ 0.01, CUR / NEW | pseudo-reverse (BOUNDS_R5, s0) | DIALibGen-predicted (PLAN_PDEC) | **OpenDIA `-decoys predict`** |
|---|---|---|---|
| lower bound | 0.0158 / 0.0161 | 0.0109 / 0.0111 | **0.0106 / 0.0106** |
| combined bound | 0.0180 / 0.0184 | 0.0124 / 0.0127 | **0.0121 / 0.0121** |
| N@FDP_c 1 % (median) | 6,563 / 9,880 | 8,318 / 9,629 | 8,233 / **10,100** |
| pair win (q ≤ 0.05) | n/a | 0.542 / 0.547 | 0.534 / 0.539 |
| coverage | n/a | 99.93 % | 99.85 % |

## Reading
1. **OpenDIA alone matches the DIALibGen route.** Both bounds are slightly lower, which is within noise. N@FDP_c 1 % is the
   same for CUR and +5 % for NEW. So the single-binary route is no worse than the two-tool route. Both remove about 85 %
   of the excess that pseudo-reverse decoys left.
2. **Both fail by the same small margin.** The true FDP at q ≤ 0.01 is about 1.06-1.21 %, against limits of 1.0 / 1.2 %.
   This residual is common to both decoy sources, so it is not the predictor's doing.
3. **The PAIRDIAG prediction held only partly** (`PAIRDIAG_ODP.txt`).
   - OpenMS never gives a decoy more fragments than its target, so the 0.648 class from PLAN_PDEC is gone. The overall
     pair win fell 0.545 → 0.537.
   - Equal-count pairs still sit at **0.542** (n 10,643), and the residual lives there.
   - Splits that move: internal K/R (0.553 without, 0.505 with). Decoys kept short win less often (0.455, n 745).
   - The RT-shift quartiles are non-monotone again. The cause is not located.
4. **What remains is common to both routes and is not about fragment counts.** Candidates for a separate precheck:
   - the scorer's self-training on an entrapment-rich target set at r ≈ 7;
   - k/r-ending effects in the pseudo-reverse sequence;
   - near-homolog structure inside the foreign proteomes.
   - None of these is tested here.

**H** (held out) is submitted after this verdict commit, with amendment 1 below. It is reported either way.

## Amendment 1 (mechanical, before H): trace tensor cap
- **Why:** H gives 13.6 M candidate windows on the r ≈ 7 library, more than the 64 GiB trace tensor cap allows (about
  11.2 M). See RESULT_PDEC_D amendment 1.
- **What changes:** H runs with OpenMS bc4e4f8d9b (`opendia-pd3-bc4e4f8d9b.sif`), `-trace_tensor_gb 96` and 192 GB.
  - bc4e4f8d9b only adds the option. Its default is unchanged, and `-decoys predict` is identical to 01d5ded04a.
- **What stays the same:** criteria, library recipe and evaluation.
H submitted 2026-11-19T02:17:41Z (jobs 35153852..61), after verdict commit 480991c28a (02:17:25Z).
