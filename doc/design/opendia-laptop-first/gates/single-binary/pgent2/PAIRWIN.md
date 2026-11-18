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

## Addendum: what decides the high-score pairs (`PAIRFEAT*.txt`, D s0 CUR `-features_out`)
- **Winners.** Among the 1,016 entrapment-versus-own-decoy pairs with best q ≤ 0.05, the entrapment wins 58 %. The
  winner is decided almost entirely by fragment co-elution evidence:

  | feature (entrapment minus decoy) | AUC for "entrapment wins" |
  |---|---|
  | `f12_coelution` | 0.97 |
  | `f12_support` | 0.97 |
  | `boundary_coelution` | 0.93 |

  This is expected, since co-elution is the signal.
- **Over all winners, false targets get co-eluting fragment evidence more often than their decoys, and the ratio grows
  with the evidence threshold** (`PAIRFEAT2.txt`):

  | threshold | entrapment / decoy ratio |
  |---|---|
  | `f12_support` ≥ 2 | 1.10 |
  | `f12_support` ≥ 4 | 1.32 |
  | `f12_support` ≥ 6 | 1.57 |
  | support ≥ 3 and co-elution > 0.7 | 1.69 |
  | support ≥ 3 and co-elution > 0.9 | 2.9 |

  This is the mechanism: the strong-evidence tail of false targets is heavier than that of decoys.
- **What it is not:**
  - **Library RT/IM.** Decoys copy the target's library RT and IM exactly, as `rtwin.py` and `rtcheck.py` show. They
    search the same RT window as their own target, so where they search is not the cause. Entrapments are more often at
    early RT, where accepted human density is lowest, which works *against* them.
  - **Fragment m/z overlap with human fragments.** Entrapment 0.984 and decoy 0.982 per fragment, so this is equal.
  - **Ion types.** The top-6 ion types and ordinals are close: the y share is 0.81 against 0.79, and y1-y3 0.16 against
    0.14.
  - **Residue context before the C-terminus.** Total variation distance to human is 0.07 for both.
- **Likely cause (not yet tested): sequence realism of the fragment ladder.**
  - A natural false peptide shares longer sub-sequences, and so runs of consecutive b/y ions, with real peptides more
    often than a reversed sequence does. The usual cause is common motifs at K/R-proximal positions.
  - Consecutive co-eluting fragments from one real interferer are what f12_support and co-elution reward.
  - **Test:** for each accepted entrapment and its decoy, count library fragments (± ppm) that match fragments of a
    single co-eluting human precursor (same RT window). If entrapments share k ≥ 3 fragments with one real interferer
    more often, the fix is interference-aware evidence: fragments claimed by a co-eluting identified precursor should not
    count. OpenDIA already has the `i_*` interference features, but they work at score level, not in TDC.
- **Running in parallel:** the predicted-decoy test (decoys with their own predicted RT, IM and intensities).

## Addendum 2: the shared-ladder test (`interferer.py`, `INTERFERER.txt`)
The question: does a winner share ≥ 3 of its 12 library fragments (± 15 ppm) with a single accepted human precursor in
the same isolation window (± 12.5 Th) that elutes within 10 s?

| group | entrapment | decoy |
|---|---|---|
| high band (score ≥ the q 0.05 cut) | 3.9 % (23 / 593) | 1.8 % (14 / 793) |
| random winners | 0.3 % | 0.1 % |

- **Shared ladders are about twice as common for false targets as for decoys,** but they explain only a small share of
  the excess.
- In the high band, entrapments exceed decoys (scaled by r ≈ 1.56) by about 160 ≈ 593 − 793 / 1.56 × … in this
  q ≤ 0.05 band. Shared-ladder cases account for about 10-20 of them.
- **So one real interferer is not the main mechanism.** Most strong false-target evidence is fragments spread over
  several co-eluting sources, or partial matches. Natural sequences hit those slightly more often than reversed ones
  across many fragments.
- **The pragmatic consequence.** The asymmetry is a property of natural versus reversed sequences on dense diaPASEF
  data. The usual remedies are:
  - decoys that look natural: predicted decoys with their own self-consistent RT, IM and intensities (being tested now);
  - or an empirical correction of q estimated from the run.
  - A run-derived correction cannot use entrapment labels in the run, so it would need a decoy model that already
    matches false targets. The predicted-decoy test decides whether that is achievable.
