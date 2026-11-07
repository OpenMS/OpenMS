# R-WP5: Close the gap to OpenSWATH + PyProphet

[Master plan](../README.md)

**Objective.** Explain and close the gap between 11830 rewrite precursors and 19132 SVM / 19725 XGBoost precursors at 1% on HeLa.

**Method:** a loss funnel over the teacher's q ≤ 0.01 IDs. At each step, count how many teacher IDs survive:
1. Library and decoy pairing eligibility.
2. A calibrated-search candidate exists.
3. A candidate apex lies inside the teacher peak bounds.
4. The precursor is in the CNN scoring manifest.
5. The CNN-selected window is the teacher peak.
6. The precursor wins its paired competition and passes q ≤ 0.01.

For each loss, record the size, an RT/IM/intensity/support stratification, and whether the loss is in candidates, choice or separation.
Fix the cause with a predeclared HeLa gate, then reapply unchanged to an independent run.

**State:** active. The findings go in this file and in the branch README.
