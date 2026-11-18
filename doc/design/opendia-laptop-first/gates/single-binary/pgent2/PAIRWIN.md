# PAIRWIN: where target-decoy competition breaks (b3, diagnostic)

Written 18 November 2026. Data: the PGENT2 library (homolog-free Arabidopsis entrapment), D, seeds 0-4. Script:
`pairwin.py`, output in `PAIRWIN.txt`.

OpenDIA's q comes from paired TDC: each target precursor competes with its own pseudo-reverse decoy. That is valid only if,
for a **false** target, the target and its decoy are equally likely to win. Entrapment targets are false by construction,
so this can be measured directly.

| pairs with best run q ≤ | entrapment pairs | P(entrapment beats its own decoy), CUR | P, NEW |
|---|---|---|---|
| 1.0 (all) | 1,175,895 | **0.500** ± 0.001 | 0.501 |
| 0.2 | 29,115 / 32,324 | 0.551 ± 0.005 | 0.547 |
| 0.05 | 5,094 / 5,954 | 0.579 ± 0.011 | 0.581 |
| 0.01 | 773 / 997 | **0.630** ± 0.029 | **0.644** |

## Reading
- **Over all pairs, the competition is exactly fair (0.500).** The decoy construction is not biased in general:
  - not by m/z, length or charge, which matches the AUC 0.49-0.51 found earlier;
  - not by the decoy's reuse of the target's RT and IM, which would show up here too.
- **The bias appears only among the high-scoring pairs.** Among pairs good enough to be accepted, the false target beats
  its decoy 63-64 % of the time against an expected 50 %. That is exactly the size of the FDP excess:
  - win share 0.63 ⇒ false targets per decoy ≈ 0.63 / 0.37 ≈ 1.7;
  - rho measured as 1.5-1.8.
- **So the asymmetry is in which pairs the scorer pushes to the top, not in the decoys as such.** Candidates:
  1. **Self-training.** The scorer learns from target winners at q ≤ 0.15 and then 0.05. False targets among those
     positives teach it the peculiarities of real target sequences that decoys lack. The pseudo-reverse decoy of a
     tryptic peptide has K/R at the C-terminus but a reversed internal sequence, and its *predicted* RT/IM is the
     target's, not the decoy's own. Among high scores, the model can pick out target-like features shared with false
     targets.
  2. **Library-coupling features.** `abs_library_rt_error`, `f12_libcos`, library correlation and IM error compare the
     peak with the library entry. A false target's library entry is a self-consistent prediction for its own sequence,
     while the decoy's entry is a mismatch: target RT and IM, target intensities on shifted fragments. A feature like
     "the predicted intensities agree with the observed ones when some peak is there" favours false targets over
     decoys only when there is a peak, which is exactly the high-score subset.
- **Next tests** (cheap, on the next `-features_out` run):
  - among high-scoring entrapment-versus-own-decoy pairs, which feature's difference predicts the winner;
  - and the predicted-decoy library already running, which gives decoys their own self-consistent predictions and so
    removes candidate 2.
