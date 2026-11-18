# PRECHECK_SHUF: do shuffle decoys fix diaPASEF q calibration? (diagnostic, D seed 0 only; H untouched)

Written 18 November 2026. The library is the same predicted targets, with decoys from OpenSwathDecoyGenerator
`-method shuffle` (sha256 dd893ecb…) instead of pseudo-reverse. Files: `EVAL_D_shuf.txt`, `RHO_shuf.txt`.

| D seed 0 | decoys | N at q ≤ 0.01 | FDP_c at q ≤ 0.01 | N@FDP_c 1 % |
|---|---|---|---|---|
| CUR | pseudo-reverse | 10,932 | 0.0163 | 9,596 |
| CUR | shuffle | 10,202 | 0.0132 | 9,766 |
| NEW | pseudo-reverse | 12,804 | 0.0161 | 12,033 |
| NEW | shuffle | 12,849 | 0.0185 | 11,660 |

rho (false targets from entrapments per decoy) near the cut, q 0-0.02:

| arm | pseudo-reverse | shuffle |
|---|---|---|
| CUR | 1.0-2.2 | 1.2-1.5 |
| NEW | 1.4-1.8 | 1.3-1.9 |

## Reading
- **Shuffle decoys do not bring calibration to nominal.** CUR moves from 1.63 % to 1.32 %, and NEW from 1.61 % to 1.85 %.
  - On one seed, both changes are within run-level noise: about 70-130 entrapments, ±25 %.
  - rho stays at 1.1-1.9 with either decoy type.
- **The excess is not caused by pseudo-reverse versus shuffle.** Even decoys built exactly like the entrapments (a
  shuffle) score below the entrapments. What separates them is that the entrapments are in the FASTA, so they go
  through the same digestion, prediction and RT/IM windows as the targets.
  - Decoys are derived from targets after prediction, and they keep the target's predicted RT and IM.
  - That is identical for both decoy methods, so it is not the difference.
- **Remaining candidates:**
  - **(i) Fragment intensities.** Entrapments get their own predicted spectra and RT. Decoys reuse the target's library
    intensities on permuted fragment m/z. So a decoy's spectrum shape is a real peptide's, while an entrapment's is a
    predicted one for a peptide that is not there.
  - **(ii) The scorer's self-training on a predicted library.**
- **The shared excess is present in the default as much as with `-pg_scores`.** It is a property of OpenDIA's
  target-decoy competition (TDC) with OpenSWATH-style decoys on predicted libraries. It is not specific to the new
  scores.

## What this means for q9
- A decoy-method change is **not** a ready fix, so a plan for (b) with shuffle decoys is not supported.
- Fixing calibration needs a deeper look at how decoys are formed for predicted libraries. One option is predicting
  decoy spectra like targets, which DIA-NN does.
- That is a separate work package and does not block `-pg_scores`.
- The case for (a), a relative plan, is now stronger. The excess is a pipeline-wide property that `-pg_scores` does not
  worsen on D and on Astral, and it worsens on H only within noise (Q7).
