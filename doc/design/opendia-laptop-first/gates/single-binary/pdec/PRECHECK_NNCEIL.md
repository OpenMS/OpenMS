# PRECHECK_NNCEIL: does a learned trace model add information for faint peaks? (20 Nov 2026, diagnostic)

This is a diagnostic, not a plan. DIA-NN labels are used for evaluation only, and the DIA-NN-labelled arms are
ceilings, not methods. It follows the user's request to revisit the neural-network scorer.

**History, in brief.**
- R-WP3 and the voxel line tried trace CNNs and TOPAZ-style models as the *whole* scorer.
- On timsTOF HeLa, TOPAZ beat grid2d by +14 %. On Astral it learned a decoy artefact (entrapment vs decoy AUC 0.52) and was stopped.
- The voxel pilot lost to CNN + LDA, and V18 (charge siblings) promoted decoys.
- The heuristic-feature network (`CrossFitScorer`) won and is the current scorer.
- What had not been tested: a trace model as **one extra input** next to the hand features, on the faint peaks of the timsTOF run.

**Data.**
- The PRECHECK_FAINT sample: 6,000 DIA-NN targets and 30,000 decoys, r ≈ 1.56 D s0, with search valley 0.8 / max 8. That gives 266,207 windows, 7,668 of them right windows.
- **Input per window:** the search's own trace matrix, 6 fragments × 33 native scans centred on the apex. Three channels: log1p and linear intensity, each scaled to the window maximum, plus the library intensities.
- **Model:** a small 1-D CNN (three convolutions, max, mean and centre pooling), about 10k parameters. It is trained per workflow fold, so it is cross-fitted, and takes about 30 s per fold on CPU.
- **Label-free training** ("lf"): positives are the tool's own q ≤ 0.01 target winners (2,861 windows), negatives are decoy windows. This is how the C++ scorer is seeded.
- Scripts: `precheck_nnceil.py`, `precheck_nnceil2.py`.

## Results

TPR of each DIA-NN target's right window against the decoy precursors' best windows:

| Arm | TPR @ FPR 1e-3 | 3e-3 | 1e-2 | Q1 @ 3e-3 |
|---|---|---|---|---|
| Tool score | 0.519 | 0.558 | 0.590 | 0.253 |
| CNN alone, DIA-NN-labelled | 0.354 | 0.431 | 0.516 | 0.195 |
| CNN alone, label-free | 0.332 | 0.406 | 0.498 | 0.161 |
| GBM ceiling, 36 features | 0.509 / 0.508 | 0.585 / 0.589 | 0.658 / 0.661 | 0.276 / 0.284 |
| + 7 faint features (PRECHECK_FAINT) | 0.535 / 0.545 | 0.614 / 0.618 | 0.693 / 0.693 | 0.304 / 0.313 |
| **+ label-free CNN logit** | **0.536 / 0.537** | **0.616 / 0.617** | **0.701 / 0.705** | **0.327 / 0.323** |
| + both | 0.554 / 0.543 | 0.622 / 0.620 | 0.708 / 0.707 | 0.331 / 0.327 |

Pairs are two GBM seeds.

## Reading

- **On its own, a trace CNN is much weaker than the hand features.** That matches the earlier history and is why replacing the scorer with a CNN failed.
- **As one extra input, the label-free CNN adds as much as the 7 hand-made faint features.** That is about +3-4 TPR points (+5-7 % relative). The gain is larger in Q1 (+0.04-0.05).
- **The two overlap almost entirely.** Both together add only about +1 point more. So the CNN learns roughly what `ps_cos_apex` and `local_sn` already measure, and there is no large hidden signal in the 1-D pooled traces.
- **The real limit is probably the input.** These traces are mobility-pooled 6-fragment XICs. Two things a richer net could use are not in them: the mobility profile, which PRECHECK_IMNARROW shows matters, and MS1. A test with mobility-resolved patches (fragments × scans × 1/K0 bins) is the next ceiling question for a net, not a bigger 1-D model.

## Decision

- No C++ net from this precheck alone. Its gain equals the cheaper hand features, which still need their own plan.
- If the net line continues, the next precheck is the same stacking test with **mobility-resolved** input, after PLAN_IMNARROW settles the window width.
