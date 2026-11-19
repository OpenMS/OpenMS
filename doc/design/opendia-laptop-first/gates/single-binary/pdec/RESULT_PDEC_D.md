# RESULT_PDEC_D: predicted decoys (DIALibGen) cut the diaPASEF excess from 1.6 % to 1.1 % at r ≈ 7, but PLAN_PDEC fails criteria 1 and 2

- **Plan:** PLAN_PDEC, frozen 9feed18a6d at 2026-11-18T21:49:48Z.
- **Runs:** D seeds 0-4, CUR and NEW. Codon `entlib-pdec5-20261118`, library `hela_pdec5.pqp` (sha256 fe85c22b…, 99.93 %
  of targets paired).
- **Evaluation:** `eval_pdec.py`, frozen sha256 bb03387d…, on the tool's own q as written. Output is in `EVAL_PDEC_D.txt`.
- **Cost:** about 2 h wall and about 92 GB peak per run.

## Criteria (D, pooled over 5 seeds, q ≤ 0.01)
| criterion | limit | CUR | NEW | verdict |
|---|---|---|---|---|
| 1 / 2 lower bound N_E/N | ≤ 0.010 | **0.0109** (90 % 0.0101-0.0117) | **0.0111** (0.0103-0.0119) | **fail / fail** |
| 1 / 2 combined bound | ≤ 0.012 | **0.0124** (0.0115-0.0134) | **0.0127** (0.0118-0.0136) | **fail / fail** |
| 3 CUR median N@FDP_c 1 % | ≥ 0.95 × 6,563 = 6,235 | **8,318** (8,186-8,467) | (NEW 9,629) | **pass** |

**Verdict: FAIL** (criteria 1 and 2). Per the plan, nothing is retuned.

**Reported, not criteria:**
- `-pg_scores` gives +15.8 % N@FDP_c 1 % (9,629 against 8,318).
- P(entrapment beats its own decoy | best q ≤ 0.05) is **0.542 / 0.547** (n 5,470 / 6,244, ±0.011). Fair would be 0.50.

## Reading
1. **Large improvement, but not calibrated.** At q ≤ 0.01 the true FDP is 1.09-1.24 % (CUR). With pseudo-reverse decoys on
   the same targets (BOUNDS_R5, s0) it was 1.58-1.80 %.
   - About 85 % of the excess over 1 % is gone. The remainder is small but clear: the 90 % intervals exclude 1.0 %.
   - The excess holds at every cut. Lower / combined bounds (CUR) are 0.0071 / 0.0082 at q ≤ 0.005 and 0.0205 / 0.0234 at
     q ≤ 0.02. So q is about 1.1-1.4× optimistic, against 1.6-1.9× before.
2. **Sensitivity at true FDP rose sharply.** N@FDP_c 1 % went from 6,563 to 8,318 (+27 %, CUR). With honest decoys, the
   scorer's ranking keeps more true targets above the real 1 % line.
3. **Pair fairness at r ≈ 7 is 0.54, not 0.50.** At r = 1.56 it was 0.50 (PRECHECK_PDEC). A residual asymmetry therefore
   remains on the harder library: entrapments still beat their own predicted decoy slightly more often.
   - Candidates to test, not done here:
     - DIALibGen decoys predicted with `no cleavage` may differ subtly from targets predicted with trypsin rules, for
       example in the charge or fragment handling of internal K/R.
     - With 7× more entrapments, near-homolog structure among the foreign peptides may matter.
     - The scorer's self-training may learn target-like features from the large entrapment-rich target set.
4. **PLAN_ODP is running in parallel**, with OpenDIA's own decoys on the same targets. It tests whether OpenMS's
   predictor behaves differently.

**H** (held out) is submitted now, after this verdict commit, as the plan says. It is reported either way.
H submitted 2026-11-19T00:45:54Z (jobs 34992382..34992392), after verdict commit 1f393c2184 (00:45:45Z).

## Amendment 1 (19 Nov 01:48): H could not run with the frozen binary (mechanical, no criteria change)
- **Failure:** all 10 H runs (34992382..92) stopped with exit 6 after the search: "Trace sample count exceeds the raw
  tensor byte budget".
  - The r ≈ 7 library on the 21-min H run gives 13.6 M candidate windows. The hard-coded trace tensor cap is 64 GiB,
    about 11.2 M windows. D (9.9 M) fits.
  - The cap was never reached before, because BOUNDS_R5 did not run H.
- **Fix:** OpenMS bc4e4f8d9b adds `-trace_tensor_gb` (default 64, so output is unchanged). H is rerun with
  `-trace_tensor_gb 96` and 192 GB of memory.
- **Container:** `entlib-pdsif-20261118/opendia-pd3-bc4e4f8d9b.sif` (sha256 b635b2ee…). Compared with a2aa8fc56c it adds:
  - 2982bd647d, the gene-join fix. It does not affect this library: 0 of 1,493,505 peptides lack a gene mapping.
  - fa75f8f46d and 01d5ded04a, which are off by default.
  - bc4e4f8d9b, the cap option.
- **What stays the same:** the frozen library, D verdict and criteria. H is still reported only, as the plan says.
