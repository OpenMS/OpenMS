# HeLa diaPASEF final benchmark (2026-11-17)

Run: 200 ng HeLa, Evosep 200SPD (`20200505_..._2737.d`), library `200623_hela_evosep_fractions_library_ptypic_decoy.pqp`
(73,922 targets, 6 fragments each, library 1/K0). Same laptop, 8 threads, run one after the other
(`bench_hela.sh`, `eval_hela.py`, `HELA_EVAL.txt`).
- OpenDIAWorkflow at OpenMS eaa2c95cd0, default settings, the library's own decoys. There is no 12-fragment scoring
  library for HeLa, so it scores the 6 library fragments.
- DIA-NN 2.2 on the library's detecting target transitions, with its own decoys (`--gen-decoys`).
- OpenSWATH + PyProphet SVM and XGBoost: the existing published-library references (`openswath-hela-published`).

There is no entrapment here. Every count is the tool's own q ≤ 0.01, a decoy estimate, not a true FDR.

| | Precursors at q ≤ 0.01 | Also in DIA-NN | Wall | Peak memory |
|---|---|---|---|---|
| OpenDIAWorkflow (this branch) | **11,653** | 93.3 % | 3:08 | 4.5 GB |
| DIA-NN 2.2 | 28,950 | n/a | 2:46 | 2.3 GB |
| OpenSWATH + PyProphet SVM | 19,132 | 88.6 % | not timed here | |
| OpenSWATH + PyProphet XGBoost | 19,725 | 87.4 % | | |
| Earlier rewrite, frozen teacher-free design (2026-11-07) | 13,253 | not computed | | |

- **OpenDIA reaches 40 % of DIA-NN and 61 % of OpenSWATH + PyProphet on HeLa.** That is far below the 81-84 % seen on
  Astral at equal true FDR.
- **It is also below the earlier HeLa rewrite result (13,253).** The Astral-driven changes since then (windowed
  calibration, D11 scorer, block I, 24 features) did not carry over to diaPASEF. This is a regression on the original
  design run.
- What OpenDIA does report is consistent: 93 % of it is in DIA-NN's list and 86 % in the SVM list. Only 301 IDs are in
  neither.
- 7,466 precursors that DIA-NN and the SVM both find are missing from OpenDIA.
- 382 of DIA-NN's 28,950 keys do not match the library's modified-sequence form. This is a matching detail and does
  not change the picture.
- Caveat: DIA-NN's q-values were anti-conservative on Astral (FDP_c 1.1-1.4 %). Without entrapment, its HeLa count is
  not a true-FDR comparison.
