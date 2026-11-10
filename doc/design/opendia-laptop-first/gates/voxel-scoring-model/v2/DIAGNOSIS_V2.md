# Voxel v2 diagnosis: where the IDs go (after the frozen result, the verdict does not change)

Scripts: `diagnose_v2.py` (raw 01), `diagnose_v2b.py`, `diagnose_v2c.py`. Outputs: `DIAG_V2_raw01.txt`, `DIAG_V2b.txt`, `DIAG_V2c.txt`.

1. **The lost hits do not have low scores.** On raw 01, for every arm, the base hits that the arm loses score above the 99th percentile of the arm's decoy winners (≥ 98 % of them).
2. **The arm puts false targets at the top, and these set the 1 % cut.** Entrapment winners in the top 1000 targets: base 5 / 1 (raw 01 / 02), V0 39 / 24, S0 24 / 19, F0 14 / 10, FS 6 / 3. Decoy winners above target #1000: base 8 / 3, V0 65 / 33, FS 7 / 11.
3. **These false targets look poor on the heuristics** (V0 top entrapment against V0 top non-entrapment, raw 01): coelution 0.18 against 0.78, 4 against 6 transitions, library RT error 49 s against 19 s, b/y ions 2/3 against 4/6. The base ranks them low (median base rank 2428).
4. **The voxel logit hardly uses coelution.** Spearman with the logit, over all decoy winners: V0 coelution +0.04, library correlation +0.18. FS is higher (+0.14 / +0.23) because the head gets the heuristics.
5. **Ranking is the main cause, and window choice is second.** With F0 / FS, 99 % of lost hits have the same winner window as the base. With V0 / S0, 18–25 % move (median 11.6 s). Swap test, raw 01: V0 scores on the base windows give 375, and base scores on the V0 windows give 817 (base 1024).
6. **No decoy tell in the library input.** Library intensities alone: target against decoy AUC 0.500.
7. **The rounds make the error stronger.** S0 takes the round-0 target winners as positives, and these include the high false targets. S0 (24 entrapment hits in the top 1000) is between V0 and F0. This agrees with the earlier TOPAZ `ss_iterations 2` result.

Possible cause (not proven): round 0 has 504 anchors only. These are bright and clean peaks with a small RT error. The voxel input has no library RT, and the CNN gets no explicit rule for "all 6 fragments must co-elute". So the model learns "a strong peak shape in some fragments". A window where 2–4 fragments hit real signal from a different peptide gets a high score. Entrapment peptides look like real sequences, so it is possible that their fragments hit real signal more often than decoy fragments do.

## Addendum: why the entrapment hits rank high (scripts `diagnose_v2d/e/f.py`, outputs `DIAG_V2d/e/f.txt`)

**Correction:** the model does not prefer entrapment over decoys. It lets in false candidates of both kinds near the top.

1. **Entrapment and decoys enter the top at almost the same rate.** V0, per 1000 winners of each kind, above the 1500th target: entrapment 2.56 against decoys 2.03 (raw 01), 1.50 against 1.22 (raw 02). V0's entrapment AUC (0.512 / 0.518) is *below* the base (0.531 / 0.535). The entrapment count looks large only because the library has about as many entrapment as human targets.
2. **The top false hits of the two kinds look the same.** On raw 01, the top entrapment hits against the decoys above the cut: MS2 sum 3.6e5 against 3.5e5, fragments > 5 % 4 against 4, library cosine 0.89 against 0.88.
3. **No sign of borrowed signal or related sequences.** On raw 01, a co-eluting base hit at ≤ 2 Th is present for 0.21 of the top entrapment hits and 0.23 of the top decoys. Shared fragments at ≤ 20 ppm: 0.05 against 0.07. A same-composition human parent: 0.09 against 0.03 (0.05 for all entrapment). The "real sequence" cause in the possible-cause paragraph above is not supported.
4. **The top false hits are bright windows where a few fragments carry the signal.** Against base human hits (raw 01): MS2 sum 3.6e5 against 1.7e5, top fragment share 0.48 against 0.34, fragments > 5 % 4 against 6, library cosine 0.89 against 0.99.
5. **Round 0 allows an intensity shortcut.** Between the anchors and the round-0 negatives, the MS2 sum alone gives AUC 0.936 / 0.920, and the count of fragments > 5 % alone gives 0.922 / 0.926. The negatives are mostly dim random decoy windows. Bright windows are 1 % of the decoy windows but 21 % / 17 % of the decoy logit top 0.1 % (enriched about 20×).

Possible cause (more likely now, still not proven): the round-0 negatives contain few bright interference windows. So the model learns "bright, peak-shaped, several fragments" and does not learn "all 6 fragments co-elute in library proportions". Any bright interference window gets a high score, whether the library entry is a decoy or entrapment.
