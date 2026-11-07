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

## Cross-run design (proposal, not implemented)

Per-run calibration and score iteration (R-WP2, R-WP3) stay the first stage. The cross-run step comes after, in two separable parts.

**A. Global score model.**
- Pool round-0 seeds from every run's calibration anchors.
- Train one CNN with outer folds keyed by peptide family across runs, so a family's windows in every run share one fold. This keeps cross-fitting honest across runs.
- Score every run with the fold model that did not see that family.
- Score iterations run on the pooled population.
- This puts scores on one scale and lets weak runs borrow from strong ones.
- Cost: one training on pooled seeds, which can be capped per run, and scoring is unchanged.

**B. Global error control** (PyProphet global context).
- For each precursor, take the best score over runs within each fold, then run one experiment-wide paired target-decoy competition on those maxima.
- Report run-specific q-values (per run, as now) and experiment-wide q-values, with run-level IDs filtered to experiment-wide q ≤ 0.01.

A without B still needs per-run q-values. B can use per-run models, but then it compares scores from different models, so it needs A or per-run score normalisation.

Gates:
- Run-specific counts with A must not fall on any run, compared with per-run models.
- B must be validated by entrapment on a multi-run set. The ProteoBench Astral entrapment runs 01–03 are a candidate, once the rewrite reads Thermo `.raw` with no 1/K0.
