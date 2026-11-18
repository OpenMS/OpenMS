# BOUNDS_R5: with r ≈ 7 the bounds are tight, and OpenDIA's diaPASEF q ≤ 0.01 fails (about 1.6 % true FDP)

Written 18 November 2026. HeLa D, seed 0, CUR and NEW (`-pg_scores`), container `opendia-pg-a2aa8fc56c.sif`, codon.

## Library (`entlib-r5-20261118`, `hela_r5.pqp`, sha256 9f2f0622…)
- **Targets:** the 58,998 HeLa library peptides.
- **Entrapments:** 294,990 homolog-free foreign peptides, five per target length slot.
  - Sources: A. thaliana, O. sativa, Z. mays and S. cerevisiae reviewed proteomes.
  - Excluded: any human tryptic peptide, any substring of a human protein, and anything within Hamming ≤ 2 of a
    same-length human or HeLa-library peptide.
  - Built by `make_fasta5_env.py`.
- **Prediction:** the PGENT timsTOF settings. Decoys: pseudo-reverse.
- **Size:** 1.39 M target precursors, r = 7.05.
- **Bounds:** with r ≈ 7 the combined (upper) and lower bounds differ by only 1/r ≈ 14 %.

| arm | q cut | N | entrapments | lower N_E/N (90 % Poisson) | combined (upper) | N@FDP_c 1 % |
|---|---|---|---|---|---|---|
| CUR | 0.005 | | | 0.0091 | 0.0103 | |
| CUR | 0.01 | 9,430 | 149 | **0.0158** (0.0137-0.0181) | 0.0180 | 6,563 |
| CUR | 0.02 | | | 0.0266 | 0.0304 | |
| NEW | 0.005 | | | 0.0084 | 0.0096 | |
| NEW | 0.01 | 11,375 | 183 | **0.0161** (0.0142-0.0182) | 0.0184 | 9,880 |
| NEW | 0.02 | | | 0.0294 | 0.0335 | |

## Reading
1. **FDR control fails on diaPASEF with a predicted library.**
   - At q ≤ 0.01 the true FDP is 1.6-1.8 %: the lower bound is 1.58 % and its 90 % interval starts at 1.37 %.
   - At q ≤ 0.005 the true FDP is 0.9-1.0 %, and at q ≤ 0.02 it is 2.7-3.4 %.
   - So the reported q is optimistic by about a factor of 1.6-1.9 across the range.
   - This settles the open question from BOUNDS.md and PRECHECK_FOREIGN. Option (b) is the right work.
2. **The excess is the same with and without `-pg_scores`** (1.58 % against 1.61 %). It is a property of the pipeline,
   not of the new scores.
   - `-pg_scores` gives **+51 %** IDs at 1 % true FDP on this harder library: 6,563 → 9,880.
3. **The larger entrapment space costs IDs.** At q ≤ 0.01, N falls to 9,430 from 11,313 on the r = 1.56 library. The
   larger library also makes the scorer's job harder, and that may change the calibration. This is the standard cost of
   high-r entrapment and is noted as a caveat: the r ≈ 7 setting is harder than a user's library.
4. **Cost:** each run took 1 h 45 min to 1 h 56 min wall, 88-90 GB peak, and 9.6 M candidate windows. That compares with
   8-10 min and 24 GB on the plain library. For a routine benchmark, r ≈ 2-3 is the affordable choice.

## Next (b3)
Locate where the asymmetry enters. Entrapments are scored as targets, so the question is why a false target outranks a
decoy at the same rank about 1.6 times as often as TDC assumes. Candidates, all testable on these outputs or on one
`-features_out` run:
1. **Decoy generation for predicted libraries.**
   - Decoys reuse the target's predicted RT and IM, and permute the fragment m/z with the target's intensities.
   - A predicted false target has its own predicted RT, IM and intensities, all self-consistent.
   - Test: build decoys by predicting decoy sequences with DIALibGen, the DIA-NN-style approach, and measure rho.
2. **Self-training.** Positives at q 0.15 / 0.05 contain false targets. That teaches the scorer target-like features
   that false targets also have. Test: measure rho by scorer round.
3. **The tail calibration** (`calibrateTail`). It maps scores by the decoy survival function per fold, which is
   monotone, so it should not change TDC ranks. That needs checking.
