# Label-use map (source inspection, no new outcomes)

OpenMS revision 3b1d90ed67. No production modifications.

- OpenDIACalibration.cpp: paired target-decoy anchors determine run calibration and RT mapping. The panel and later evidence calibration use decoy labels. This is upstream of the network.
- OpenDIAWorkflow.cpp: target-only anchor positives, paired and sampled decoy negatives initialize the network. RT refinement selects target q<=0.01 results. Target/decoy siblings are grouped separately by unmodified sequence.
- OpenDIAScorer.cpp: fold-specific network fitting, self-training positives chosen by paired TDC, all decoy precursor winners negative. Median/MAD and tail calibration use decoy winners. Final paired competition uses (D+1)/T and reverse minima.
- OpenDIAGBM.cpp: target q positives, hardest 2% decoy precursor negatives, folds for scoring, decoys for calibration and final TDC.

Consequences: masking only final GBM labels or splitting final q-estimation decoys is not automatically an end-to-end RESET implementation or proof. Earlier calibration, representative selection, feature generation, and normalization must be included in any such design. Cross-fitting does not by itself establish all null exchangeability conditions. Candidate counts and sibling collection sizes can affect maxima, so audit the actual population before score filtering when available.
