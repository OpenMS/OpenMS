# PRECHECK_JOINT: does a joint (presence-weighted) interference fit have room to help? (16 November 2026)

Diagnostic, written before any plan for the joint fit of `DESIGN_INTERFERENCE_GRAPH.md`.
- **Inputs:** arm-E arrays (`gap-20261114/armE_raw0?`) and the C++ block I (`cpp/featI_*`), with the C++ CrossFitScorer
  at seed 0 (`cppscore/raw0?_s0.npy`).
- **Grouping:** DIA-NN 2.2 (same targets, own decoys) is used only to pick the DIA-NN-supported right-peak winners that
  we reject ("R"). They are compared with the decoy winners nearest in score ("D", 3 per R) and with entrapment winners
  in the same band ("E", with their own matched decoys, as a specificity control).
- **Scripts:** `claimant_split.py`, plus the inline scripts behind `RESIDUAL_SIGNAL.txt` and `CEILING_current.txt`.

## 1. Who explains the shared signal in block I? (`CLAIMANT_SPLIT_raw0?.txt`)
Explained share of a window's central area, split by the claimant window's own score.
- "present" means the claimant scores above today's 1 % cut. That is about 1,850 windows per run.
- "absent" means everything else.

| group (median explained share) | all claimants | present | absent, its precursor's winner | absent, not its winner |
|---|---|---|---|---|
| R, raw01 | 0.223 | 0.000 | 0.037 | 0.185 |
| D, raw01 | 0.270 | 0.000 | 0.049 | 0.217 |
| R, raw02 | 0.214 | 0.000 | 0.044 | 0.176 |
| D, raw02 | 0.238 | 0.000 | 0.050 | 0.191 |

AUC(R > D) of −explained is 0.538 / 0.514 with all claimants, and 0.504 / 0.505 with present claimants only.

- Almost all of block I's explaining-away comes from claimants that are themselves probably absent. These are other
  precursors' non-winning candidate peaks, which sit on the same fragment ions.
- A joint fit that weights each claimant by its presence would remove most of that explanation. At matched score it
  barely separates R from D. Present claimants almost never overlap R or D.
- So block I mostly works as a **"this fragment ion is crowded" signal, not as real assignment**. A presence-weighted
  joint fit would discard that signal and gains nothing in its place.

## 2. What still separates R from D at matched score? (`RESIDUAL_SIGNAL.txt`)
Inputs with |AUC − 0.5| > 0.1 for R against D, while entrapment-in-band against its decoys stays near 0.5 (specific to
real peptides). The AUCs are raw01 / raw02.

| input | AUC R > D | AUC E > D |
|---|---|---|
| transition_support | 0.82 / 0.82 | 0.49 / 0.51 |
| f12_libcos | 0.80 / 0.77 | 0.51 / 0.47 |
| i_libcos_resid | 0.77 / 0.76 | 0.50 / 0.49 |
| abs_library_rt_error (lower = R) | 0.33 / 0.34 | 0.50 / 0.55 |
| d_max, d6_max (lower = R) | 0.28-0.33 | 0.47-0.53 |
| boundary_library_correlation | 0.68 / 0.61 | 0.53 / 0.53 |

- The rejected real peptides differ from their score-matched decoys in **library pattern agreement, fragment support
  and RT agreement**. The scorer already has these inputs, but it does not combine them well enough for this weak
  class.
- This is a scorer-capacity question rather than a missing-evidence one.

## 3. Labelled ceiling on the current 24 inputs (`CEILING_current.txt`)
This is an upper-bound probe, not a method: HistGradientBoosting, 3-fold by the family fold.
- Positives are the DIA-NN-supported right-peak target winners; negatives are all decoy winners.
- The table gives the TPR among the currently rejected R.

| | at today's decoy rate | at 4× that rate |
|---|---|---|
| raw01 | 0.197 (109 of 553) | 0.461 |
| raw02 | 0.203 (100 of 493) | 0.552 |

With labels and the same inputs, about 100-110 more per raw are separable at today's decoy rate. That is roughly the
size of the gain to expect from a better use of the inputs. It is about 6 % of the current N, and a label-free method
will get less.

## Reading
- **Not worth a plan:** the full presence-weighted joint fit. The evidence it would use (present claimants) almost
  never touches the rejected IDs.
- **Worth a small plan:** a score-band or local feature normalisation.
  - The inputs that separate R from D are the ones whose meaning changes with signal level: transition_support,
    library cosine, RT error.
- **Worth a larger plan:** new evidence, which the ceiling says is needed for most of the 500.
  - Fragments the instrument actually sees, from run-specific fragment selection.
  - MS1 isotope evidence, which DIA-NN uses and our MS1 block barely does.
  - Ion-mobility agreement, not available on this Astral data.
