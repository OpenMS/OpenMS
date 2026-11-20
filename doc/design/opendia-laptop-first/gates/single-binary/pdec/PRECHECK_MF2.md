# PRECHECK_MF2: does raw-point information add anything on top of the new default? (20 Nov 2026, diagnostic)

**Question.** PRECHECK_MF found that raw-point matched filtering (`mf`) beat the tensor version (`mt`) against the old 36-feature baseline. Now that `mt` is the default, is there still extra signal in raw points? If there is, which ingredient carries it (mass-error weighting or own-mobility centring), and can the extraction pass keep it without storing a point cloud?

## Setup

- **Sample:** faint2, the same 6,000 DIA-NN targets and 30,000 decoys on r ≈ 1.56 HeLa D as PRECHECK_MF. Raw points come from `maps.xipm`.
- **Baseline B45:** the new default, which is the 36 features plus the 9 C++ `mt_` features, from `mtcheck-20261120/runs/D_mt_s0`. That run is byte-identical to the new default (OpenMS 01c5891ce0). "Missed" means DIA-NN sample IDs the new default does not accept at q ≤ 0.01 (2,269).
- **Script:** `precheck_mf2.py` (sha256 `313e0e6c…`), codon job 39208795. The gate was written into the script header before the job was submitted.
- **Variants:** each computes the same 9 matched-filter statistics on a different trace.
  - Diagnostic, from raw points on the native frame grid:
    - `mf`: ppm weighting plus a Gaussian around the window's own mobility apex.
    - `p`: ppm weighting only.
    - `m`: own-mobility Gaussian only.
  - Implementable in the extraction pass, on the workflow's 0.94 s RT grid:
    - `tp`: today's 8-bin tensor with ppm-weighted accumulation.
    - `t16`: 16 IM bins over ±0.06, which doubles tensor memory.
    - `t16p`: `t16` with ppm weighting.
- **Gate:** an implementable variant had to add at least +0.030 TPR at decoy FPR 3e-3 in the label-free tool arm on missed IDs, and at least +0.015 in the ceiling arm over all IDs, on both seeds.

## Result (TPR at decoy FPR 3e-3, delta over B45, seed 0 / seed 1)

| Variant | Ceiling, all | Ceiling, missed | Tool, all | Tool, missed |
|---|---|---|---|---|
| B45 (absolute) | 0.783 / 0.788 | 0.296 / 0.314 | 0.771 / 0.764 | 0.266 / 0.245 |
| mf (raw, both) | +0.024 / +0.015 | +0.089 / +0.061 | +0.002 / +0.009 | +0.018 / +0.034 |
| p (raw, ppm only) | +0.011 / +0.004 | +0.038 / +0.015 | −0.004 / −0.002 | −0.003 / −0.001 |
| m (raw, own mobility only) | +0.018 / +0.004 | +0.065 / +0.030 | +0.004 / +0.001 | +0.017 / +0.015 |
| tp | +0.001 / +0.003 | +0.009 / +0.012 | −0.002 / −0.007 | −0.003 / −0.015 |
| t16 | +0.007 / +0.001 | +0.025 / +0.006 | +0.008 / −0.003 | +0.027 / +0.000 |
| t16p | +0.010 / +0.001 | +0.039 / +0.010 | +0.005 / −0.003 | +0.020 / +0.001 |

**Gate: all three implementable variants fail.** No C++ change follows.

## Reading

- **Most of the raw-point advantage was already captured once `mt` became the default.** PRECHECK_MF compared `mf` with the old 36 features, and that comparison overstated what raw points add on top of `mt`. Over B45, even the full raw-point `mf` adds only +0.018 / +0.034 on missed IDs in the label-free proxy, about a fifth of what `mt` itself added.
- **Mass-error weighting adds nothing** (`p`: −0.003 / −0.001). The 10 ppm extraction window already removes most off-mass points.
- **Own-mobility centring carries what little is left** (`m`: +0.017 / +0.015), mostly in the DIA-NN-labelled ceiling. A wider and finer mobility grid in the tensor (`t16`) does not keep it: one seed gains and the other does not.
- **The label-free tool proxy is now the binding constraint.** The ceiling arm (DIA-NN labels) still sees +0.06 to +0.09 on missed IDs for `mf`, but a scorer trained only on targets and decoys realises little of it. Further gains on faint peaks are more likely to come from how the scorer learns, for example seeds and hard negatives, than from more extraction detail.

## Decision

No point-cloud or extended-tensor implementation. `-matched_filter` stays as the default as it is. The faint-peak (B1) route via extraction features is closed for now. The next candidates are the RT-calibration weakness on H (95th-percentile residual 52 s against 8 s on D) and wrong-peak or no-candidate cases, which this route does not touch.
