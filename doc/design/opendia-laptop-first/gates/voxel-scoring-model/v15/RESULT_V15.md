# RESULT_V15: D11 heuristic amplitude audit

**Scope:** raw01 source audit and descriptive feature comparisons; no model was refit, no counterfactual features were scored, and no labels were used in training. The code-level conclusions come from the current `feature/opendia-clean-rewrite` implementation. The DIA-NN D/E annotations below were used only for evaluation summaries.

## Finding

**None of the 11 D11 heuristic inputs directly encodes absolute intensity magnitude.** Under a uniform positive rescaling of the measured peaks, the features are unchanged, provided peaks remain positive and candidate identities are held fixed:

- `transition_support` counts matched transitions. `log_matched_scans` counts matched scans/points. Neither uses an intensity magnitude threshold. `log_peak_width` uses the candidate time span; local-peak splitting, when enabled, uses a relative peak-height fraction. `abs_library_rt_error` is determined by the apex RT and library RT map.
- `boundary_log_fwhm` and `boundary_log_width` are measured relative to peak height/fractional boundaries. `boundary_coelution` is Pearson correlation. `boundary_library_correlation` is Pearson correlation of square-root fragment areas and square-root library intensities. `boundary_apex_agreement` uses fragment argmax positions. These are invariant under uniform positive scale.
- `byion_b` and `byion_y` count theoretical b/y ions with a matched peak whose intensity exceeds `byion_intensity_min_`. The default floor is 0.0, so they encode detected-ion presence, not how far above zero the ion is.

Thus the existing scorer can learn evidence about fragment presence, chromatographic shape, alignment, and agreement with the library. It cannot infer that an otherwise unchanged measured peak is, for example, four times weaker, because that amplitude is not present in these inputs. Fragment counts can still reflect signal loss when a peak disappears or crosses zero, but they are not continuous amplitude measurements.

## DIA-NN Stage-D versus accepted-E feature summary

Raw01 rows: 866 Stage D and 1,413 accepted E. The value `P(D>E)` is the pairwise probability that a randomly selected D row has a larger feature value than a randomly selected E row. This is descriptive, not a causal or error-control measure.

| Feature | D median | E median | P(D > E) |
|---|---:|---:|---:|
| transition_support | 4 | 6 | 0.1303 |
| log_matched_scans | 2.890 | 2.708 | 0.5716 |
| log_peak_width | 3.979 | 3.812 | 0.5635 |
| abs_library_rt_error | 19.11 | 14.97 | 0.5716 |
| boundary_log_fwhm | 1.820 | 1.689 | 0.6704 |
| boundary_log_width | 2.579 | 2.648 | 0.3917 |
| boundary_coelution | 0.2855 | 0.7665 | 0.0879 |
| boundary_library_correlation | 0.5062 | 0.8485 | 0.2625 |
| boundary_apex_agreement | 0.8333 | 1.0000 | 0.2033 |
| byion_b | 2 | 4 | 0.1772 |
| byion_y | 4 | 7 | 0.0961 |

D rows have fewer supported transitions and b/y ions, weaker co-elution, lower library-intensity agreement, and lower apex agreement than accepted E rows. Those are quality/presence differences, not absolute amplitude values. They also do not show whether the features separate true D targets from decoys.

## Stage-D versus baseline-score-matched decoys

As a descriptive control, each D row was matched to the closest D11FWT decoy winner from the same fold, without reusing a decoy, requiring a score gap <=0.5. This produced 617 of 866 possible matches; median score gap was 0.00018 and the 90th percentile was 0.0666. Matching is conditional on the existing model score, so it is not a prospective classifier test.

| Feature | Median D - matched decoy | Paired fraction D > decoy |
|---|---:|---:|
| transition_support | +1 | 0.507 |
| log_matched_scans | -0.118 | 0.447 |
| log_peak_width | +0.041 | 0.520 |
| abs_library_rt_error | -7.85 | 0.397 |
| boundary_log_fwhm | +0.011 | 0.517 |
| boundary_log_width | +0.080 | 0.532 |
| boundary_coelution | -0.012 | 0.481 |
| boundary_library_correlation | +0.129 | 0.592 |
| boundary_apex_agreement | 0 | 0.452 |
| byion_b | 0 | 0.298 |
| byion_y | 0 | 0.384 |

Most single features are close to chance in this score-conditioned comparison. Library correlation has a modest positive separation, while b/y counts remain lower for D than for matched decoys. This is consistent with the earlier finding that native fragment evidence can favor entrapment/decoys as well as weak human targets. It does not justify training an augmentation model.

## Decision

Do not train a model on intensity-scaled profiles while leaving these heuristic inputs unchanged. That augmentation would not teach absolute abundance to the model. The preceding V14b diagnostic also showed the F12 path contributes almost no median score change under attenuation. The observed D/E feature differences are real, but they describe missing-fragment and shape/coherence quality and do not cleanly separate D from similarly scored decoys.

A defensible next experiment, if continued, should be a frozen no-training counterfactual that recomputes the five boundary descriptors and F12 from the same local-baseline-attenuated six- and twelve-fragment profiles, leaving count-based search and b/y inputs fixed for their proven scale invariance. It should first validate recomputed boundary features against the stored values and stop if baseline reproduction is not close. This would test whether baseline-preserving attenuation changes peak-shape inputs, not whether the model can infer absolute intensity. Existing V9-V11 native-evidence failures remain a hard caution against promoting any new evidence without an entrapment gate.
