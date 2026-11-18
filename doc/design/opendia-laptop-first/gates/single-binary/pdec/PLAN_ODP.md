# PLAN_ODP: OpenDIA's built-in predicted decoys (`-decoys predict`) make diaPASEF q hold at r ≈ 7

Written 18 November 2026, 23:46, before any PLAN_PDEC result was seen. No r ≈ 7 search with predicted decoys has been
evaluated. This plan does not depend on PLAN_PDEC's verdict.

**What is already known (so D is not blind):** DESIGN_OPENDIA_DECOYS precheck 2 ran on the r = 1.56 library, D seed 0.
- With `-decoys predict`, q ≤ 0.01 gave a lower bound of 0.62 / 0.53 % and a combined bound of 1.02 / 0.86 % (CUR / NEW).
- 99.9 % of targets were paired.
- H has never been searched with `-decoys predict`.

## Input
- **Library:** `entlib-r5-20261118/hela_r5.pqp` (sha256 9f2f0622…), the BOUNDS_R5 library.
  - Targets: 58,998 HeLa library peptides plus 294,990 homolog-free foreign entrapments, r = 7.05.
  - Decoys: pseudo-reverse, which OpenDIA drops and replaces.
- **Binary:** OpenMS 01d5ded04a, container `entlib-pdsif-20261118/opendia-pd2-01d5ded04a.sif` (sha256 3cb281b9…), with the
  PeptDeep models in place.
- **Comparison with PLAN_PDEC:** same targets. The only intended difference is who makes the decoys: OpenDIA at run time
  here, DIALibGen there. The binaries also differ, but 01d5ded04a changes nothing on the default path. The
  `predict_keep_short` and gene-join fixes are inside it.

## Arms and runs
- **CUR:** `-decoys predict`. **NEW:** `-decoys predict -pg_scores`.
- **Resources:** codon, 8 threads, 128 GB, 4 h.
- **D** (design): seeds 0-4, both arms.
- **H** (held out): seeds 0-4, both arms. Submitted only after the D verdict is committed.
- **Evaluation:** `pdec/eval_pdec.py` (sha256 bb03387d…), on the library each run wrote
  (`<out_dir>/hela_r5.predicted_decoys.pqp`), using the tool's own q as written.
  - The written library must be byte-identical across runs, which is checked by sha256.
  - Entrapment labels are read only in the evaluation.

## Criteria (on D, at the tool's q ≤ 0.01; pooled over 5 seeds)
1. **CUR calibrated:** lower bound ≤ 0.010 and combined bound ≤ 0.012.
2. **NEW calibrated:** lower bound ≤ 0.010 and combined bound ≤ 0.012.
3. **Sensitivity:** CUR median N@FDP_c 1 % ≥ 0.95 × 6,563. 6,563 is BOUNDS_R5 D s0 CUR with pseudo-reverse decoys on the
   same targets, a single-seed comparator.
4. **Coverage:** at least 99 % of target precursors have a predicted decoy (from the `Decoy prediction` log line).

## Outcome
- **Criteria 1, 3 and 4 pass:**
  - `-decoys predict` is proposed as the default for predicted libraries. Making it the default is a separate commit with
    ctests.
  - The switch to tell predicted from empirical libraries is still open. Its design item 3 needs its own test.
- **Criterion 2 also passes:** `-pg_scores` is proposed as the default if NEW's median N@FDP_c 1 % ≥ 1.05 × CUR's.
- **Any failure:** reported as is, with nothing retuned.
- **H** is reported either way.

**Diagnostic, not a criterion:** P(entrapment beats its own decoy | best q ≤ 0.05), expected 0.50.

**Caveats:**
- r ≈ 7 is harder than a user's library.
- Near-homolog false matches are outside this benchmark.
