# Voxel v3 phase 1b: loss or flexibility? (written before the run, after RESULT_P1.md)

One arm, as N11 with one change: the head is linear (11 -> 1, no hidden layer, no dropout). Same BCE loss, optimizer, sampled negatives, schedule and early stopping as N11.

- G11 within 5 % of L11 at N@1 % on both raws, and entrapment AUC within 0.005 of L11: the gap comes from the flexible head, and the BCE loss is not the cause. Phase 2 then uses a linear final layer over the learned features.
- G11 near N11 (gap to L11 > 5 %): the BCE loss / schedule is the cause. Phase 2 then changes the loss, not the head.
