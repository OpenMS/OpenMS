# Local diagnostic rules and provenance

26 November 2026. Selection tests were specified in executable source before running and are recorded here after their outputs. The known-null stress below is specified before its first run. Used data only. No remote submission, no held-out raw processing, no candidate claimed validated.

1. Production GBM training selection is challenged by changing scores only in the excluded fold. Compare selected training row sets and fixed probe predictions. The local-selection control excludes that fold BEFORE TDC/negative mining. This tests a direct dependence only, not fully nested stacking.
2. Apply the same selection comparison to D/H NET seed 0 full precursor-winner tables. Export uses no protein annotations. Negative counts are precursor winners, not all negative windows used by production. Frozen perturbation is +0.5 to held-out decoys and -0.5 to held-out targets, in the supplied score units. Perturbed data are scratch-only.
3. Known-null stress: seeds 0..9, 12,000 numeric pairs, 3 paired folds, one candidate each, 6 features. For 80% of pairs the target's first feature is shifted by +3; all other target and decoy features are iid standard normal. Seed score is feature 0. Run production GBM defaults, report every seed and pooled N/E and realized FDP, retaining any failures. Ground truth is evaluation-only. This is a synthetic sanity check, not a proteomics or cross-species validation.
4. Parallel pair calibration screen uses two fixed label-free rules documented separately before evaluation. No retuning on its results.

Success for a real calibration candidate remains original corrected-label combined bound <= 1% with no loss in median N@FDP_c 1%, on fixed complete seed sets. No human/organism/entrapment classification may influence training, groups, thresholds or selection.
