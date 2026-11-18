# PLAN_PDEC: predicted decoys (DIALibGen) make OpenDIA's diaPASEF q-value hold at r ≈ 7

Written 18 November 2026, after PRECHECK_PDEC (including its "Full coverage" section) and before this plan's library is
built or searched.

**What is already known (so D is not blind):**
- PRECHECK_PDEC used the PGENT2 library (r = 1.56), D seed 0.
- There, the tool's own q ≤ 0.01 gave a lower bound of 0.69 % / 0.74 % and a combined bound of 1.13 % / 1.21 %
  (CUR / NEW).
- The r ≈ 7 library with predicted decoys has never been searched, and H has never been searched with predicted
  decoys.

**Decision (user, 18 Nov):** predicted decoys come from DIALibGen for now. Generating them inside OpenDIA is the
downstream goal and a separate work package.

## Library (`entlib-pdec5-20261118` on codon)
- **Targets and entrapments:** exactly `entlib-r5-20261118/targets.pqp`, the same one as BOUNDS_R5.
  - 58,998 HeLa library peptides plus 294,990 homolog-free foreign peptides, giving 1.39 M target precursors, r = 7.05.
  - The timsTOF prediction settings are those of BOUNDS_R5.
- **Decoys:** for every unique target peptide in targets.pqp, the pseudo-reversed sequence (last residue fixed) is
  predicted by DIALibGen 0.13.1 like a target.
  - The prediction uses the same settings as the targets, with `-generation:enzyme 'no cleavage'` and
    `-generation:decoys none`, via `make_decoy_fasta2.py`.
  - Each target precursor is paired with the decoy precursor of the same charge and the exactly reversed modified
    sequence (each modification moves with its residue), via `merge_predicted_decoys.py`.
  - **Pairing v2 was fixed before this plan.** On the PRECHECK_PDEC library it pairs 441,694 / 441,910 targets
    (99.95 %). Targets without a decoy stay in the library.
  - Every decoy peptide gets a gene mapping, which works around the PQP-reader INNER JOIN.
- **Hashes:** script hashes are in `inputs.sha256`. The built library's sha256 is recorded in `merge.log` on codon
  before any search.

## Arms and runs
- **CUR:** a2aa8fc56c without `-pg_scores`. **NEW:** the same binary with `-pg_scores`.
- **Setup:** container `opendia-pg-a2aa8fc56c.sif`, the same binary as BOUNDS_R5. No code change in this plan. Codon,
  8 threads, 128 GB.
- **D** (design): seeds 0-4, both arms.
- **H** (held out): seeds 0-4, both arms. Submitted only after the D verdict is committed.
- **Evaluation:** `eval_pdec.py <library> runs D|H`, run on the tool's output as written. Entrapment labels are read
  only there. Pooling is over 5 seeds.

## Criteria (on D, at the tool's q ≤ 0.01)
1. **CUR calibrated:** pooled lower bound ≤ 0.010 **and** pooled combined bound ≤ 0.012.
2. **NEW calibrated:** pooled lower bound ≤ 0.010 **and** pooled combined bound ≤ 0.012.
3. **No sensitivity cost at true FDP:** CUR median N@FDP_c 1 % ≥ 0.95 × 6,563.
   - 6,563 is the pseudo-reverse decoy result on the same targets (BOUNDS_R5, D s0, CUR).
   - Caveat: the comparator is a single seed.

## Outcome
- **Criteria 1 and 3 pass:**
  - Predicted decoys become the documented way to search predicted libraries with OpenDIA. The DIALibGen recipe:
    `make_decoy_fasta2.py`, predict, then `merge_predicted_decoys.py`.
  - The OpenDIA-side decoy predictor becomes the next work package.
- **Criterion 2 also passes:** `-pg_scores` is judged on the same bound.
  - If NEW's median N@FDP_c 1 % ≥ 1.05 × CUR's on D, `-pg_scores` is proposed as the default. Astral is restated from
    PGENT2: pooled FDP_c 0.0073, with gains on raws 01-03.
  - Making it the default is a separate commit with ctests.
- **Any failure:** reported as is, with nothing retuned. The bounds and the pair-fairness number from `eval_pdec.py`
  show which part failed.
- **H** is reported either way, with the same numbers.

**Diagnostic, not a criterion:** P(entrapment beats its own decoy | best q ≤ 0.05) per arm, expected 0.50.

**Caveats:**
- r ≈ 7 makes the scorer's job harder than a user's library, as noted in BOUNDS_R5.
- Near-homolog false matches are outside this benchmark (sibling competition is still an open item).
