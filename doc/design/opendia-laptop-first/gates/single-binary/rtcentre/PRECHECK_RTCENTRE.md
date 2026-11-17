# PRECHECK_RTCENTRE: stage-3 window centre on HeLa diaPASEF (diagnostic, 2026-11-17)

Data: HeLa 200 ng, `OpenDIACalibrate -panel_targets 40000 -out_diagnostics` at 0e9e4fc3a5 (`hela_diag.tsv`); DIA-NN 2.2
q ≤ 0.01 IDs (28,568 in the library) used for evaluation only.

The single calibration pass works in three stages.
- Stage 1 fits a map (map1) to the anchors that win their target-decoy pair without any RT information.
- Stage 2 runs an EM over all target candidates and gives a smooth 11-knot map (em_map).
- Stage 3 searches ±max(4σ1, 30 s) around **em_map** and fits the final map to the anchors it finds.

The windowed second pass (Astral) centres stage 3 on **map1**.

Map error against DIA-NN's RT for DIA-NN's IDs:

| map | abs err p50 | abs err p95 | within 3σ1 (6.8 s) | early (library RT < 12) median err |
|---|---|---|---|---|
| stage-1 map1 | 1.46 s | 4.88 s | 0.990 | −0.1 s |
| stage-2 em_map (current stage-3 centre) | 6.16 s | 19.56 s | 0.571 | −19.7 s |
| final isotonic map | 2.26 s | 19.31 s | 0.837 | −19.4 s |

- em_map is 16-29 s late below library RT 12, and 5-12 s late up to 50. It is 7 s early at the end.
  - It has 11 knots, a fixed σ (7 % of the run span, 23.6 s) and a prior/curvature penalty. This is too stiff for the
    steep start of a 5.6-min Evosep gradient.
- The stage-1 anchors (4,625 at q ≤ 0.01) cover library RT from −4.6. Their robust spread is 2.55 s around map1 and
  7.85 s around em_map.
- So stage 3 searches the early precursors in a ±30 s window that is off-centre by about 20-29 s. On top of that, the
  stage-3 RT likelihood (`0e4a234762`, 14 Nov) penalises anything far from that centre with σ1 = 2.3 s.
  - Result: below library RT 12 only 22 targets pass q ≤ 0.01 at stage 3, against 363 at stage 1.
  - The final map has no knot before 110.7 s.
- This happens with every one-pass calibration. The Astral runs all use the windowed pass (2 passes on raws 01-03),
  which already centres on map1, so they are not affected.

Proposed change (PLAN_RTCENTRE): in the one-pass path, centre stage 3 on map1 and treat library RTs outside the
stage-1 anchor support as the windowed pass does (search to the run edge, no RT term there). em_map stays as a
diagnostic. The windowed path is unchanged, so Astral output is byte-identical.
