# RT window from the density crossing: rule undefined on these runs (pre-plan check)

Done 16 November 2026, before any plan was frozen or any code written. It uses no labels: only the stage-4 mixture fit
written in each cal.json (`rt.spread_s`, `rt.spread_signal_share`) and the stage-3 half window.

**Proposed rule.** Set the RT tolerance where the fitted true-peak density pi N(x; 0, s) equals the chance-peak density
(1 − pi) / (2 hw), the way the fragment mass tolerance is set.

| calibration | pi | spread s | hw | true-peak / chance density at the centre |
|---|---|---|---|---|
| refined, arm E raw01 | 0.119 | 24.8 | 110.7 | 0.48 |
| refined, arm E raw02 | 0.111 | 24.9 | 97.6 | 0.39 |
| refined, PAIRLIB raw01 | 0.121 | 22.6 | 85.1 | 0.41 |
| refined, PAIRLIB raw02 | 0.111 | 25.2 | 94.9 | 0.37 |
| current code, arm E raw01 | 0.218 | 37.7 | 110.7 | 0.65 |
| current code, arm E raw02 | 0.214 | 30.9 | 97.6 | 0.68 |

The true-peak density stays below the chance density everywhere, so the rule has no crossing and is undefined. Counting
both sides of the window does not help either: there is still no crossing.

**Why it differs from the mass tolerance.**
- The mass tolerance is fitted over centroids of real fragment matches, where signal dominates near the centre.
- The RT mixture is over each panel member's best stage-3 peak. Most panel targets are not detectable in the run, so
  only about 11 % of the winners are true peaks (pi).
- A precursor's chance competitors inside the window do not enter this mixture at all.
- The density comparison therefore does not describe the trade-off the window controls. That trade-off is a true peak
  outside the window against more chance peaks competing with it inside.

**A well-posed version (not run, needs approval).** Choose w to maximise the expected share of detectable precursors
whose true peak is the best candidate in their window:

  P_win(w) = P(true peak within ± w) × P(true score > every chance score in a 2w window).

- The first factor comes from the spread fit.
- The second comes from the stage-3 score distributions: decoy winners' scores per unit RT give the chance-score rate,
  and the anchors' scores give the true-score distribution.
- It uses no labels and no fixed multiple, and it is computed per run. It would be tested as proposed: five seeds on
  raws 01 and 02, and raw 03 held out.

## The P_win version, computed offline (label-free, 6 calibration runs with diag.tsv)

For each run, w is chosen to maximise P(true peak within ± w) × E[F_chance(true score)^(w / hw)]:
- the true-score distribution comes from the target winners minus (1 − pi) × the decoy winners;
- the chance-score CDF over a 2 hw window is the decoy-winner CDF.

| run | spread | best w | in spreads | P_win at best w / at 4 spreads |
|---|---|---|---|---|
| arm E raw01 (current spread) | 37.7 | 71.8 | 1.90 | 0.659 / 0.559 |
| arm E raw01 (refined spread) | 24.8 | 50.3 | 2.03 | 0.728 / 0.638 |
| PAIRLIB raw01 / raw02 (refined) | 22.6 / 25.2 | 47.0 / 52.4 | 2.08 / 2.08 | 0.741 / 0.735 |
| arm D raw01 / raw02 | 25.5 / 26.7 | 53.7 / 54.9 | 2.11 / 2.06 | 0.732 / 0.712 |
| arm C raw02 | 25.7 | 54.1 | 2.11 | 0.737 |

**The criterion points the wrong way.** It asks for a window of about 2 spreads, which would cover about 90 % of true
peaks instead of about 96 %. The evidence runs the other way: narrowing raw02 from 124 to 100 s already cost IDs over
five seeds (PLAN_CALREFINE2).

**Reason.** P_win assumes the window's best stage-3 score decides. The final scorer does not work that way: it weighs
the RT error as a feature and uses 24 features. So a wider window costs less than the stage-3 competition implies.

**Conclusion.** Neither label-free rule is sound on this data. The window width is a trade-off that only the final ID
counts measure. No plan was frozen and no code was changed.
