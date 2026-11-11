# Stage-D deep dive (D11FWT, old population): findings

All descriptive. DIA-NN's report is used only as a reference. Entrapment labels define the groups but never reach training.
Groups:
- **D**: right window, below the cut (866 on raw 01, 929 on raw 02).
- **E**: accepted (1413 / 1362).
- **DEC**: decoy winners.
- **ENT**: entrapment winners.

Scripts: `dd_evidence.py`, `dd_separation.py`, `dd_diann.py`, `dd_band.py`, plus inline checks recorded in the `.txt` outputs.

## 1. D is the right peak (DD_DIANN.txt)
- Our winning apex agrees with DIA-NN's RT for both groups. The median difference is 0 s. The 90th percentile is 3.3 s for D and 1.6 s for E.
- DIA-NN's FWHM is 3.2 to 3.3 s. So D is not a peak-picking or RT problem.

## 2. D is low-abundance signal that DIA-NN is also less sure of (DD_DIANN.txt, DD_DIANN_EXTRA.txt)

| median, raw 01 (raw 02 similar) | D | E |
|---|---|---|
| DIA-NN log10 Precursor.Quantity | 5.30 | 6.37 |
| DIA-NN log10 Ms1.Area | 5.83 | 6.62 |
| DIA-NN Evidence | 2.76 | 4.72 |
| DIA-NN Ms1.Profile.Corr | 0.53 | 0.82 |
| DIA-NN PEP | 0.011 | 0.000 |
| DIA-NN q ≤ 0.001 share | 43 % (46 %) | 85 % (91 %) |
| our lib_cos (6 fragments) | 0.93 | 0.99 |
| our corr_med | 0.12 | 0.33 |

- The D share follows abundance and DIA-NN's own confidence closely (the shares below are for raw 01 / raw 02):

| DIA-NN quantity tertile | D share, DIA-NN q ≤ 0.001 | D share, DIA-NN q > 0.001 |
|---|---|---|
| high | 4–5 % | 34–36 % |
| mid | 23–25 % | 63–86 % |
| low | 55–56 % | 88–93 % |

- Charge 3+ is over-represented in D: 29 / 27 % against 20 % in E.
- At matched DIA-NN quantity decile, our extracted area is about the same for D and E. Our extraction therefore sees the signal that DIA-NN quantifies.
- At that matched quantity, D still has lower lib_cos, corr_med and DIA-NN Evidence than E. So D are the noisier precursors at a given abundance, also by DIA-NN's own measures.
- DIA-NN's best fragment is our top library fragment for 33 % of D and 20 % of E. It is outside our 6 library fragments for 9–11 % of D and 17 % of E.

## 3. Nothing measured carries information the score leaves unused (DD_BAND.txt, DD_BAND_MODEL.txt, DD_BAND_ALLHUMAN.txt, DD_LIBCOS.txt)

**D against decoys at matched score.** We compared D with decoys at the same score: 10 strata at the D score quantiles, over 18 raw-evidence measures and the 17 model inputs.
- lib_cos separates D from decoys with AUC 0.80 / 0.79. Entrapment against decoys is 0.57 / 0.54.
- `top1_share`, `n_sig` and `n_peak` also separate them (0.6–0.7).
- `log_area` is lower for D (AUC 0.23 / 0.30, i.e. smaller than score-matched decoys).

**Selection control.** D is chosen by DIA-NN, which itself rewards library agreement. We repeated the comparison on **all** human target winners against decoys at matched score.
- The lib_cos gap shrinks to 0.62 against 0.49 on raw 01 and 0.58 against 0.52 on raw 02.
- Every other measure is within ±0.08.

**Direct size test** (diagnostic, not a proposal). We ranked by score + k × z(lib_cos) for k = 0.5 to 4.
- N@FDP_c 1 % changes by +14 / +7 at most.
- N@2 % changes by −18 to +17.
- FDP_c at q ≤ 0.01 changes erratically.

So the D-against-decoy separation on lib_cos is mostly DIA-NN's selection, not usable information.

## Conclusion (softened)

On this population, stage D looks like a **sensitivity limit at low abundance**, not a defect in a single step.
- Our window is DIA-NN's peak.
- Our extraction sees comparable signal.
- At the D score, none of 35 measures of fragment evidence, MS1, co-elution, mass error, interference or library agreement separates true IDs from decoys beyond what the selection effect explains.

DIA-NN accepts these precursors mostly at its marginal q values (0.001–0.01), with its own Evidence and PEP close to its threshold.

The plausible remaining difference is evidence we do not measure. Examples are DIA-NN's use of more fragments, its own learned per-fragment quality, and its MS1 model. **This is not shown here.**

No feature is proposed from this dive. The next useful step is the rtlin rebuild (PLAN_V8). It changes the population, so the D set is recomputed there.
