# BOUNDS: what the entrapment numbers prove (correction to PGENT, PGENT2, PRECHECK_* readings)

Written 18 November 2026. Scripts: `bounds.py` (shuffled 1:1 library v1, D and H, 5 seeds each) and the PGENT2 numbers.

**The FDP_c used in every PGENT/PGENT2 table so far is the "combined" estimator, N_E (1 + 1/r) / N.** Wen et al.
(Nat. Methods 22, 1454-1463, 2025) show it is an **upper bound** on the FDP. It can only show that FDR control is
**valid**. Their lower bound, N_E / N, can only show that control **failed**. When the two bounds fall on either side of
the nominal level, the analysis is inconclusive. Their paired estimator gives a tighter upper bound at peptide level when
targets and entrapments are paired one-to-one.

## Precursor level at q ≤ 0.01 (mean over 5 seeds)
| library | run | arm | lower N_E/N | combined (upper) | paired upper (peptide level) |
|---|---|---|---|---|---|
| shuffled 1:1 (v1), r = 1.12 | D | CUR | 0.0093 | 0.0176 | 0.0180 |
| shuffled 1:1 (v1), r = 1.12 | D | NEW | 0.0088 | 0.0166 | 0.0174 |
| shuffled 1:1 (v1), r = 1.12 | H | CUR | 0.0088 | 0.0167 | 0.0178 |
| shuffled 1:1 (v1), r = 1.12 | H | NEW | 0.0099 | 0.0187 | 0.0200 |
| homolog-free Arabidopsis (PGENT2), r = 1.56 | D | CUR | 0.0087 | 0.0143 | |
| homolog-free Arabidopsis (PGENT2), r = 1.56 | D | NEW | 0.0097 | 0.0159 | |

## Correct reading
- **In no run does the lower bound exceed 1 %.** So none of the data shows that OpenDIA's q ≤ 0.01 fails on diaPASEF.
  The upper bounds exceed 1 %, so the data does not show that it holds either.
- **The true FDP lies between about 0.9 % and about 1.6-1.9 %. It is inconclusive.**
- "q is about 1.5 times optimistic" (PGENT2 verdict) and "about 1.7 % true FDR" (PGENT result, README, board) overstate
  it. Each upper bound was reported as if it were the value.
- Both PGENT and PGENT2 failed their frozen criterion 1 (combined FDP ≤ 0.012). That stands. But criterion 1 was an
  upper-bound test, so failing it does not mean FDR control failed.
- `-pg_scores` and CUR are indistinguishable on every bound: the lower bounds are 0.88-0.99 % against 0.87-0.93 %.
- **The decoy-versus-false-target "rho" tables use the same combined estimate,** so rho > 1 there is also only an upper
  bound on the asymmetry.
- The gap between the two bounds is about 1/r. Wen et al. recommend r ≥ 5, for a gap of ≤ 20 %. Every library here
  has r ≈ 1.1-1.6, too small to decide.

## Next (diagnostic, no plan needed)
- Build a homolog-free Arabidopsis entrapment library with r ≈ 5. The bounds then sit within about 20 % of each other,
  which settles whether diaPASEF q is calibrated.
- Then decide:
  - either `-pg_scores` can be judged against the bound;
  - or there is a calibration defect to fix after all.
