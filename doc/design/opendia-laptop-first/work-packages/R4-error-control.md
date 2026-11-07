# R-WP4: Error control

[Master plan](../README.md)

**Objective.** Report precursor q-values with independent evidence that they control error.

**Delivered:**
- Family outer-fold cross-fitting with fold-wise q-values (`c355109de2`).
- Opt-in paired target-decoy competition, where `DECOY_<X>` competes with `<X>` (`ae4687e508`).
- Recorded failure: the log(n)-window correction failed its gate (`0c757c3d2d`).

**Open gates:**
1. An entrapment or two-species check on HeLa-like data.
2. An independent-run result without a same-run teacher.
3. A best-of-n-windows decoy bias analysis. More windows per precursor raises the best decoy score.
