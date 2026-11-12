# V18 look: charge-state siblings as an entrapment-neutral signal

Status: descriptive pre-plan look, no training, no gate changed, no IDs recovered. 12 November 2026.

## Why this signal

Every earlier evidence arm (MS1, mass error, native counts, unique fragments, profile correction, deconvolution, cross-run) scored a precursor from its own isolation-window data. Weak true peaks and entrapment matched to interference look alike in that data, so each arm rewarded entrapment about as much as weak human peaks.

A different charge state of the same peptide is measured independently: different precursor m/z, usually a different isolation window, and different fragment m/z values. A real peptide elutes all its charge states at the same RT. A random match does not, because the second charge state is an independent draw. Decoys also have charge siblings (89% of decoy winners have one in the library, the same rate as targets), so the signal is symmetric for target-decoy competition. Entrapment peptides are not shuffled human sequences (4.6% share a composition), so they get no shared-sequence bonus.

Feature (label-free): for each precursor winner, the best per-fold decoy-normalised D11F score among candidate windows of its sibling charge states whose apex lies within 3 s (`sib_z_near`).

## Result

| | raw01 | raw02 |
|---|---:|---:|
| Winners with a sibling in the population, D / all | 47% / 41% | 42% / 41% |
| ENT vs DEC AUC of `sib_z_near` | 0.508 | 0.504 |
| All human vs DEC AUC | 0.510 | 0.506 |
| D vs DEC AUC | 0.632 | 0.656 |
| E vs DEC AUC | 0.701 | 0.690 |
| Own score at decoy p99 to p99.99 **and** `sib_z_near` > 2: human / entrapment / decoy | **78 / 0 / 1** | **99 / 2 / 3** |
| of which D | 54 | 71 |

The raw02 criteria and thresholds (2 and 3 s) were fixed from raw01 and hashed before raw02 was read (`PLAN_LOOK_RAW02.*`). Raw02 met all of them: entrapment ≤ 3, decoy ≤ 3, human ≥ 40, ENT-vs-DEC AUC ≤ 0.52.

## Limits

- **Modest yield.** About 80 to 100 precursors per raw fall in that band. Most of the deep D group (below decoy p90) has no strong sibling, so this signal does not reach it.
- **Partly redundant at peptide level.** On raw01, 38 of the 54 band D already have a sibling charge state above p99.99. These are new precursor IDs of already-found peptides, not new peptides.
- **Only about 40% of winners have a sibling in the population**, so the feature is missing for most rows. It must be coded as missing, not as weak.
- The band and threshold come from a look, not a frozen rule. A model must use the feature as an input, scored the same way for targets and decoys, not a hand rule on this band.
- Raw02 is replication after many experiments, not a pristine holdout.

## Proposed V18 (not run)

Add `sib_z_near` and a has-sibling indicator to D11F. Compute them from a sibling-free first pass to avoid circularity. Keep the same training, W, T and the unchanged pass rule. Freeze before training. No hard-negative changes.
