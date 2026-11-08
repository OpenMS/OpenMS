# Post-mortem on the spent raw 02 run (descriptive, not a gate)

Raw 02 was used once, for `astral-entrapment-raw02-20261108/GATE.md` (FAIL). This analysis reads that run's outputs together with its entrapment labels. It is used only to choose the next design question. Nothing here is a new test.

Notation: z = per-fold decoy-standardised logit (the `-pooled_tdc` scale). There is one precursor winner per precursor. Classes are human target, entrapment target, and decoy.

## 1. Can the CNN tell targets from decoys? Partly, and only in the extreme tail

Per-capita tail rate of precursor winners (per 10^4, count in brackets):

| z > | human | entrapment | decoy of entrapment | decoy of human |
|---|---|---|---|---|
| 4 | 594 (8199) | 480 (6522) | 400 (5291) | 393 (5274) |
| 6 | 185 (2555) | 79 (1072) | 64 (846) | 61 (824) |
| 8 | 77 (1058) | 8.6 (117) | 7.2 (95) | 6.3 (85) |
| 9 | 51 (697) | 1.7 (23) | 1.8 (24) | 1.1 (15) |
| 10 | 31 (421) | 0.5 (7) | 0.3 (4) | 0.3 (4) |

- **The CNN does separate real signal:** human targets are enriched about 10 to 100 times over decoys above z 6.
- **Entrapment targets, which are false by construction, are scored above decoys at every depth:**
  - about 1.2 times at z 4 to 8;
  - about 1.5 to 1.7 times at z 10, though there are only 7 against 4.
- **Decoys of entrapment and decoys of human behave alike,** so the decoy generator treats the two sides equally. The excess is target-specific.
- **Consequence:** a decoy count underestimates the number of false targets by about 20 to 50% in the tail. That fits the 2 times excess FDP seen at q ≤ 0.01, given the small counts.

## 2. FDP_c along the q cut

| q cut | accepted | N_e | FDP_c | FDP_c / q |
|---|---|---|---|---|
| 0.005 | 340 | 5 | 0.030 | 5.9 |
| 0.01 | 360 | 5 | 0.028 | 2.8 |
| 0.02 | 423 | 7 | 0.033 | 1.7 |
| 0.03 | 509 | 8 | 0.032 | 1.1 |
| 0.05 | 635 | 14 | 0.044 | 0.89 |
| 0.1 | 886 | 52 | 0.118 | 1.2 |

**Correction to RESULT.md:** at q ≤ 0.05, FDP_c (0.044) is below q. The q-values are anti-conservative only in the shallow tail, q ≤ 0.03.
- **Why the shallow tail fails:** the smallest reachable q is 1/340, because only 2 decoys score above the accepted cut (z 10.39). So the q ≤ 0.01 set is effectively "everything above the 3rd decoy".
- **What that means:** at this depth the error estimate rests on 2 to 3 decoys. The real problem is the low power of the run, roughly 360 IDs where OpenSWATH-scale runs give thousands, not a miscalibrated q at depth.

## 3. The 5 accepted entrapment hits

- All have z 10.6 to 11.3, and their paired decoys have z 1.4 to 3.6, so each wins its pair clearly.
- **Chimeric or isobaric human signal explains only 1 of 5:**
  - AALEEAALPESR/2 has the same composition as the human EALAEAALESPR/2. Both apexes are at 837.0 s, and 4 of 6 fragments are shared.
  - The other 4 have no human winner within 20 ppm and 15 s that shares fragments.
- **Isobaric share in the wider entrapment set:**
  - 4.2% of entrapment precursors have a same-composition human precursor (I = L);
  - 6.3% of entrapment winners above z 6 do;
  - so isobaric interference is a minor contributor, not the main one.
- **Other 4 hits, not tested here:** most likely co-eluting unrelated signal that matches 1 to 3 fragments, with the CNN rewarding trace shape. The `.xipm` export (below) exists to check this directly.

## 4. What this says about the model

- **The ranking is not the main failure.** It is real but weak, and on this run there are about 360 precursors above the decoy floor.
- **Shuffled decoys are a mildly easier null than false targets,** by about 1.2 to 1.7 times in the tail, so decoy counts are optimistic in the extreme tail.
  - Wider entrapment-style checks are needed before trusting q ≤ 0.01 on weak runs.
- **Candidate directions,** each needing its own design run and gate:
  1. **Model:** add fragment-level evidence the CNN lacks, such as MS1 precursor trace, isotope pattern, fragment m/z error and library intensity correlation. These are the scores OpenSWATH uses to reject partial matches.
  2. **Decoys:** a harder null, for example decoys scored at target-like RTs, or shuffled windows taken from target-dense regions.
  3. **Error control:** report q with a minimum decoy support, or use a mixture-model or PEP estimate that does not depend on 2 to 3 decoys.

The `.xipm` review below checks direction 1 on extraction data.

## 5. Extraction review with `OpenDIA -out_xipm` (commit `e781cf6064`)

- **Export run:**
  - The frozen raw 02 search was re-run with `-out_xipm` on 290 selected precursors plus their paired partners (`run_xipm.sh`, `xipm_selection.tsv`).
  - Search output is byte-identical to the evaluation run (`cmp`).
  - The export produced 2610 peak maps and 260482 raw points, 2.1 MB, in 2:30 and 4.4 GB.
  - `review.py` decodes the arrays (`numpress.py`, which reproduces the logged point count exactly, with every point inside its m/z and RT window).
- **Measured per precursor, within ±4 s of its CNN winner apex:**
  - share of library fragments with signal;
  - mean pairwise fragment-trace correlation (co-elution);
  - correlation of apex fragment intensities with library intensities (sqrt scale);
  - fragment and MS1 ppm error;
  - MS1 monoisotope and isotope presence.

| group | n | fragments with signal | co-elution | library corr | abs fragment ppm | MS1 present | abs MS1 ppm |
|---|---|---|---|---|---|---|---|
| accepted human | 40 | 1.00 | **0.91** | **0.94** | 0.3 | 1.00 | 1.2 |
| accepted entrapment | 5 | 1.00 | **0.29** | **0.29** | 1.0 | 1.00 | 7.1 |
| top decoys | 20 | 0.92 | 0.20 | 0.28 | 1.6 | 0.85 | 5.5 |
| human, z 6 to 10.4 (missed at 1%) | 40 | 1.00 | 0.35 | 0.63 | 0.6 | 0.90 | 2.7 |
| entrapment, z 6 to 10.4 | 20 | 0.83 | 0.25 | 0.23 | 1.2 | 0.95 | 5.9 |
| human, z < 1 | 20 | 0.50 | −0.22 | 0.26 | 4.1 | 0.70 | 7.4 |

- **Simple rule (co-elution > 0.6 and library correlation > 0.6):**
  - keeps 34 of 40 accepted human IDs;
  - removes all 5 accepted entrapment hits, all 20 top decoys and all 20 entrapment winners above z 6.
- **Missed humans:** among humans missed at 1% (z 6 to 10.4), 12 of 40 pass the rule.
- **This is a hand-made post-hoc rule on the spent run.** It is not a result. It shows that the information that rejects the false hits is in the extraction but is not used by the CNN score.
- **What the five entrapment hits look like:** each has signal on 4 to 6 fragments and an MS1 peak. But the fragments do not co-elute (0.29) and do not match the library intensity pattern (0.29). Their MS1 m/z is 4 to 20 ppm off, against 1.2 ppm for accepted human IDs.
  - This looks like unrelated signal that happens to fall on each fragment m/z. It is not a peptide eluting.

## 6. Answer

- **The CNN can tell targets from decoys,** with seed AUC near 1 and human targets enriched 10 to 100 times. **It cannot tell true targets from random co-located signal well enough,** because it scores per-trace shape and intensity.
  - It does not see what OpenSWATH-style scores check: whether fragments co-elute with each other and with MS1, whether relative intensities match the library, and the MS1 mass error.
- **Next step:** add those as explicit inputs or scores, and test on a new design run under a new gate.
  - A first form would be a second-stage score: the CNN logit plus co-elution, library correlation and MS1 ppm, combined by a linear or semi-supervised discriminant as in PyProphet.
  - Applying the hand rule above to raw 02 would be tuning on evaluation data, so it is not reported as an ID count.
