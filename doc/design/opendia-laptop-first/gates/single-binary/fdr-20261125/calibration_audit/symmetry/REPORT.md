# Target-decoy search and competition diagnostic

Run `bash /home/singjc/.jcode/scratch/symmetry-20261126/run.sh`.

This compiles this checkout's actual OpenDIASearch.cpp, OpenDIAPeakPicker.cpp and OpenDIAScorer.cpp, linking the preexisting wp03-release libOpenMS solely for core spectrum types. It writes only under this scratch directory. Synthetic data uses neutral A/B identities and numeric pair IDs, no species, proteins, or entrapment labels. No jobs are submitted, source changed, or refs updated.

## Executed checks

- Same two precursors with three identical transitions, 21 Gaussian-intensity MS2 scans, calibrated RT=10 s, RT half-window 20 s. Swapping only target/decoy flags preserves both candidate counts and every SearchCandidate field exactly.
- Identical target/decoy evidence gets identical area and support.
- Twelve pairs in three folds: strict higher scores win independently of target/decoy flags; reverse row iteration leaves ranks and q-values unchanged.
- All-target winners give q=1/12, confirming +1 correction and reverse minimum. Label-swapped all-decoy winners give q=1, intentionally asymmetric FDR reporting.
- Exact within-pair ties go to the decoy even after labels swap.
- A target with its decoy absent remains a valid competitor. This is not a proof of bias: search can legitimately have no candidate for one member, but pair-level population gating does not create a missing candidate.
- Pair-max support gate replica keeps both existing rows for support [4,1], drops both for [2,2] at threshold 3. This check is explicitly a replica, not the private workflow function.

## Source observations, not full workflow execution

All paths below relative to src/topp/OpenDIA/source/.

- OpenDIASearch.cpp:40-85 validates the same settings for all queries and stores/checks is_decoy consistency. The label has no role in matching or peak picking (the only occurrences are lines 77,80). Search identity is (precursor_id, precursor_mz), not target flag.
- OpenDIAWorkflow.cpp:107-169 makeQueries applies calibration and RT correction to every precursor, with no target-specific branch. :702-716 applies pair maximum transition support, preserving every existing candidate of a passing pair.
- OpenDIAWorkflow.cpp:383-403 strips DECOY_ for pair identity and hashes the target unmodified peptide family for a decoy if its target exists. fold=hash%3, rt_fold=hash%5, shared pair_seed from stripped name. :417-419 explicitly copies target fold to all paired rows. Orphan decoys fall back to their own sequence family. These folds were inspected, not exercised via library loading.
- OpenDIAScorer.cpp:352-363 picks one row per precursor, equal scores choose lower row. :365-386 standardizes by each fold's decoy-winner median/MAD. :388-415 pair competition favors decoys on exact ties and cross-pair equal-score sorting also places decoys first, then q=min(1,(D+1)/T). Missing decoy fold or zero spread throws rather than silently proceeding.
- OpenDIAScorer.cpp:628-648 feature mean/SD is computed across all training-fold rows without a label branch and then applied to all rows. Distribution is window-weighted, not precursor- or pair-balanced. Decoy-specific calibration is intentional, not label-symmetric preprocessing.
- OpenDIAWorkflow.cpp:724-805 feature extraction takes all surviving precursors through identical top-fragment and grid setup. OpenDIAStage2.cpp, OpenDIAInterference.cpp, OpenDIAFragmentGrid.cpp have no decoy-label branches. OpenDIATraces.cpp:93-94 legacy seed path rejects positive decoys; :131 copies a decoy label to metadata.

## Material search asymmetry to review

OpenDIAWorkflow.cpp:598 admits only target winners q<=0.01 to RT-refinement training. :637 calls fitRTRefine, :640-649 overwrites correction for the entire library using refine.fit.predict (full fit), including the training targets. :658 uses out_of_fold only to estimate residuals/narrowed tolerance. :609 runs the second search with full-fit corrections. Thus selected targets have in-sample RT centers while decoys are never training observations. This is a concrete code-path difference, but its magnitude and FDR impact were NOT measured here. It is distinct from species-dependent selection, which is absent. Decoy and target sequences can also yield different predicted corrections as intended.

The diagnostic demonstrates label independence conditional on fixed query coordinates and evidence, not end-to-end null exchangeability, decoy-library quality, or validity of GBM training/scoring. Those need separate controls.
