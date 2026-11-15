# PLAN_CALSPREAD result: fails offline; no C++ change

The frozen design was a window-frame mixture with a member-shifted Gaussian and the stage-3 anchors excluded
(`calspread_offline.py`, output in `OFFLINE.txt`). It fails criteria 1 and 2 on 9 of 10 runs:
- sigma ranges from 4.6 to 20.8 s, against DIA-NN's robust spread of 24.4-27.1 s, a ratio of 0.18-0.81;
- coverage of DIA-NN 1 % IDs within 4 sigma is 0.52-0.95.

It passes only on arm E raw01, with sigma 25.8 s (DIA-NN 27.2 s) and coverage 0.952. Per the plan, the C++ change is not
made.

**Reading:** with the anchors removed, the remaining stage-3 target winners carry little peak signal (pi 0.01-0.16). The
Gaussian then fits whatever local excess is left near the map, which is narrower than the true RT error. The fit is
unstable at this signal level (arm C raw01: pi 0.009, sigma 4.6 s). Removing the anchors removes most of the signal that
measures the spread.

**Where this leaves the self-calibrating workflow:**
- The current spread fit stays as it is. It matched DIA-NN within 0.9-1.1 on the earlier runs, but gave 1.39 (raw01) and
  1.22 (raw02) on the arm-E library with the current calibration code (RESULT_WORKFLOW).
- None of the three offline variants tried is better on all runs.
- A different approach would need its own plan. Candidates:
  - fit the spread on held-out anchors, such as the score-only anchors of sigma_check.py;
  - iterate stage 3 once around the final map, so the window centre and the final map coincide and the existing fit's
    frames match.
