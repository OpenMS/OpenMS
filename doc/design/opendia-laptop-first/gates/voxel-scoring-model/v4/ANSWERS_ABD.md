# A, B and D: follow-up diagnostics (descriptive; entrapment and DIA-NN used only to measure)

## A+B: C++ calibration
- **Code:** the constant end value is in `RTMap::operator()` (`OpenDIARTMap.cpp`).
- **Users of that map:**
  - the OpenDIA search, through `readCalibrationRT`;
  - OpenDIACalibrate: the stage-1 map, the windowed second pass and the final anchor residuals;
  - `OpenDIACNNRTWindow`.
- **Knot range is narrower than the anchors.** The knots are bin medians, so they stop half a bin inside the anchor range: 0.098–0.752 against anchors 0.061–0.847 on raw 01.
- **Default change:** a byte-identical default needs an opt-in option, for example in OpenDIACalibrate, that writes the extended map to cal.json.

## A+B: end shape (`RT_END_COMPARE.txt`)
- **Method:** each map is fitted on the middle 70 % of the anchors and scored on the held-out 15 % at each end. This is label-free.

| Map | Held-out end anchors: median \|res\| (s), raw 01 / 02 | Within tolerance, raw 01 / 02 |
|---|---|---|
| Constant (current) | 67 / 69 | 0.66 / 0.60 |
| Linear (L3, LF, LOWESS ends) | 17–22 | 0.99–1.00 |
| Quadratic | 30 / 23 | 0.79 / 0.99 |
| PCHIP | 38 / 23 | 0.79 / 0.80 |
| Extreme-anchor end knots | 25 / 20 | 0.81 / 0.99 |

- **Distance-dependent tolerance:** the label-free fit gives a slope of 0, so the data do not support a wider window.
- **Limit:** beyond the anchor range, about half of the DIA-NN reference stays outside the tolerance with every map (median 67–124 s). Curve shape does not fix this. Only anchors at the gradient ends can show what happens there.

## D: why cross-run evidence failed (`DIAG_ASYM.txt`, `DIAG_OWNERSHIP.txt`, `DIAG_MBR.txt`)
- **Asymmetric for false targets.** At matched own score, entrapment winners have an other-run window more often than decoys:
  - raw 01: 0.52 against 0.46 (q99–q99.9 bin) and 0.74 against 0.68 (above q99.9);
  - raw 02: 0.48 against 0.43.
  - So decoys cannot calibrate cross-run evidence on this benchmark.
- **Not shared signal.** A higher-scoring winner with shared fragments explains only 1–4 % of the top winners, about the same for decoys (2–4 %) and entrapment (1.3–1.6 %).
- **Not same-mass partners.** These are as common for entrapment as for decoys (about 0.7–0.85).
- **Mechanism not known.**
- **Pair-level MBR selection, symmetric by construction, also fails.** N@2 % falls (1716 → 1616 and 1736 → 1626) and FDP_c at q ≤ 0.01 is erratic (0.011–0.019), because the decoy counts become small.

## D: symmetric levers (`DIAG_D_LEVERS.txt`, `DIAG_RT_COMPOSITION.txt`)
- **Peptide-specific features look the same for entrapment and decoys but differ for human.** Medians above the decoy q99.9 score (human / decoy / entrapment):

| Feature | raw 01 | raw 02 |
|---|---|---|
| \|RT error\| (s) | 15 / 22 / 20 | 15 / 20 / 17 |
| 12-fragment library cosine | 0.95 / 0.83 / 0.84 | 0.95 / 0.83 / 0.86 |

- **Label-free RT refinement.** A ridge fit of the anchor residuals on amino-acid composition, length and charge, applied to non-anchor reference precursors, cuts the robust SD of the RT residual from 23.7 to 20.6 s (raw 01) and from 24.2 to 20.1 s (raw 02).
