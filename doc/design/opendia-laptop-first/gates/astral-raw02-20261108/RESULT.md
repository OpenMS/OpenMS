# Result: Astral raw 02 single evaluation (frozen Option B design) - FAIL

Run: 2026-11-08 02:05:46Z to 02:27Z, once. `GATE.md`, `run.sh` and `evaluate.py` pass `sha256sum -c GATE.md.sha256`, so they were unchanged since 02:05:30Z. No bring-up fixes were needed.

## Gate quantity

| pooled paired q | accepted targets | N_h | N_e | FDP_c | FDP_lb |
|---|---|---|---|---|---|
| ≤ 0.001 | 0 | 0 | 0 | n/a | n/a |
| **≤ 0.01** | **360** | 355 | 5 | **0.0280** | 0.0139 |
| ≤ 0.05 | 635 | 621 | 14 | 0.0444 | 0.0220 |

**GATE: FAIL** (FDP_c 0.0280 > 0.015). r = 148983 / 151240 = 0.9851.

- FDP_c 95% interval (Clopper-Pearson on N_e / 360, times 1 + 1/r): 0.009 to 0.065.
  - If the true FDP_c were exactly 0.015, there is a 13% chance of 5 or more entrapment hits in 360.
  - So the excess is not conclusive. The pass rule was fixed in advance, however, and it fails.
- At q ≤ 0.05 the estimate (0.044, 14 hits) is also near 1.9 times the nominal level in the same direction. The experimental q underestimates the entrapment FDP at both cuts.

## Other recorded quantities (not part of the rule)

- **Calibration:**
  - Stage 1 gave 0 anchors and the windowed second pass gave **604 anchors**.
  - RT σ 20.7 s, MS2 offset −2.44 ppm.
  - 4:27, 16.5 GB.
  - Anchor entrapment share 0.0099 (6 of 604), so the purity lower bound is 0.980.
- **Search:** 2:52, 4.3 GB.
- **CNN crossfit:** 13:31, 22.1 GB.
  - Held-out seed-decoy AUC is 0.9995 to 0.99997 per fold, so the seed task is close to trivial.
- **Accepted targets per fold:** 197 / 8 / 155, with entrapment 2 / 0 / 3.
  - All 5 entrapment hits have q = 0.0029, the smallest q reached: 1/340 at this decoy count, since nothing reaches q ≤ 0.001.
  - Their logits (13.0 to 15.5) are among the highest. They are not borderline calls.

## Reading (descriptive, after the fact, not a new design)

- The calibration fix transfers to raw 02: 604 anchors, against 512 on raw 01.
- Seed-decoy AUC near 1 together with confident entrapment hits suggests the following. The CNN separates seeds from the shuffled-decoy windows on features that entrapment targets share, for example trace shape of any co-eluting signal. Decoys then do not model false targets well enough at this depth.
- This is a hypothesis only. Testing it needs a new design run (raw 01 or another non-evaluation run) and a new gate. Raw 02 is now spent, and raw 03 is still held back.
