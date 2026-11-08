# Anchor peak-boundary labels and boundary model prototype (raw 01, diagnostic)

These are design-data diagnostics, not a gate. They spend no unused run.
- Labels come only from our own calibration anchors (targets).
- No decoys, search candidates or external tools are used as training signals.
- OpenSWATH boundaries are compared only to check scale.

## Step A: labels from a simple picker (`OpenDIAPeakPicker`, commit on `feature/opendia-clean-rewrite`)

**The picker:**
- sums the fragment traces into a consensus trace;
- smooths it with one [1 2 1]/4 pass;
- takes the apex within ±2 points of the hint;
- walks outwards while the trace falls, stopping at a valley or at 5% of the apex;
- interpolates the half-height crossings;
- counts the fragments whose maximum is within 1 point of the apex.

It is simpler and faster than MRMTransitionGroupPicker or PeakPickerChromatogram (238 ns per 6×32 window, so the full raw 01 set of 2.8M windows takes under 1 s).

**Run:** `OpenDIACalibrate -out_anchor_peaks -out_anchor_traces` on the anchors' raw-scan traces (`run.sh`). The default outputs `cal.json`, `anchors.tsv` and `diag.tsv` are byte-identical to `../window300k` (`CMP.txt`).

**Label quality** (`labels.py` → `labels.txt`):

| | value |
|---|---|
| anchors / clean labels | 512 / 512 (apex within 1 scan of calibration apex, ≥4 agreeing fragments, no edge) |
| boundary reason, left / right | floor 471 + valley 41 / floor 358 + valley 154 |
| FWHM median (p10 to p90) | 4.4 s (4.1 to 5.0); stage-1 FWHM in `cal.json` is 3.3 s (integer spectra count) |
| base width median (p10 to p90) | 14.7 s (9.9 to 18.0), 10 scans |
| tailing (right / left half-width) | 2.0 |
| same peak in our search output | median width 37.7 s, about 2.6× wider |
| same peak in OSW output (495) | median width 13.1 s; boundary differences: left 1.6 s, right 3.3 s (median) |

**Reading:** these labels are close to OpenSWATH's boundaries, which we did not use. The search's 30 s gap rule is the cause of the wide widths.

## Step B: augmentation and a small segmentation model (prototype in Python)

**Augmentation** (`src/topp/OpenDIA/augment_peak_labels.py`): 64 windows per clean anchor, 32768 in total, on the CNN tensor geometry (6 × 32 bins, ±15 s).
- RT shift of ±6 s;
- peak scale 10^U(−3, 0) on real off-peak background from the same fragments;
- noise proportional to √I;
- drop 0 to 2 fragments;
- 30% of windows get interference on one fragment;
- 20% get a neighbouring peak;
- 20% are background-only negatives.

Each anchor is assigned to one of 3 folds by a hash of its precursor ID. The peak signal is the anchor's whole raw slice. An earlier version cut it at the label boundaries, which left a step the model could learn; that version was rerun and its numbers are not used.

**Model** (`segmenter.py` → `segmenter.txt`):
- a per-fragment convolution shared across fragments, with mean and max pooling over fragments;
- 4 dilated RT convolutions;
- per-bin outputs for inside-peak, start and end;
- 7.6K parameters, inputs normalised to their window maximum then log-scaled;
- AdamW, one-cycle LR schedule, gradient clipping at 5, 12 epochs.

Training takes 59 s for all 3 folds and inference 44 µs per window (Python, 4 threads). It is scored on anchors from the held-out fold only.

| held-out, log10 SNR | n | IoU model | IoU rule | start within 1 bin, model / rule | end within 1 bin, model / rule |
|---|---|---|---|---|---|
| < 0 | 1483 | **0.66** | 0.34 | 45% / 18% | 26% / 12% |
| 0 to 0.5 | 1313 | **0.75** | 0.35 | 65% / 21% | 26% / 12% |
| 0.5 to 1 | 1373 | **0.80** | 0.40 | 73% / 27% | 33% / 12% |
| ≥ 1 | 21983 | **0.84** | 0.52 | 80% / 42% | 37% / 15% |

- **Peak detection:** the model finds 97% of held-out positives and calls 66% of background-only windows "no peak". The rule picker always returns a peak.

## Limits
- **The labels are a rule's output.** The model learns to reproduce the rule on raw scans from noisier, binned inputs. That is why it beats the rule applied to the same bins, but it is not ground truth.
- **The augmented windows are synthetic.** On real candidate windows the gap may be smaller.
- **The end boundary is the weak point:** tailing makes the right edge diffuse, and only 37% of ends are within 1 bin.
- **Negative rejection is modest** (66%).
- **Not yet shown:** that better boundaries add IDs. That needs native integration (boundaries feed width, area, co-elution and library correlation over [left, right]), an opt-in stage-2 feature group, and a predeclared gate on an unused design run.
- **Folds:** the hash split stands in for the CNN's outer family folds. It must be the same split before any boundary-derived feature reaches stage 2.
