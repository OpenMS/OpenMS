# NASE reproducibility and score ties

Candidate hits with equal scores at the `report:top_hits` boundary are retained. Before ambiguity resolution, ties are ordered by candidate sequence, source oligo sequence, precursor charge/isotope/adduct and mass error. Exact ties in `IdentificationData::getBestMatchPerObservation` are resolved by molecule type, sequence/compound identifier, charge and adduct values. Neither allocation order nor target/decoy status is used as a tie-breaker.

When FDR uses one best match per observation, the chosen representative receives the q-value. The other equally scoring candidates remain available in identification outputs but may not appear in BED. A lexical representative is a reproducibility convention, not evidence for unique modification localization. Future runs may therefore differ from arbitrary representatives chosen in historical exports. Near-equal floating-point scores are not rounded into ties.

BED mapping CSVs may use mxA/mxC/mxG/mxU as aliases for the engine's mA?/mC?/mG?/mU? codes. An explicit engine-code entry takes precedence. IDs are taken from the supplied CSV, not hard-coded.

## Historical preprocessing

The base add_bedrmod branch changed the hard-coded WindowMower movement from `jump` to `slide`. This changes the experimental peaks used for scoring even when all previously exposed INI settings match. `preprocessing:window_mower:movetype` now exposes that choice: the default remains `slide` for compatibility with the branch, while `jump` reproduces the historical native/Docker filtering. The Becker regression explicitly selects `jump`. New precursor-removal and theoretical-ion blacklist options remain disabled/empty for the historical comparison.

## Validation

Run `IdentificationData_test`, `BedRModFile_test`, `FalseDiscoveryRate_test`, and the NucleicAcidSearchEngine TOPP tests after building this branch. `TOPP_NucleicAcidSearchEngine_reproducibility` uses 16 Becker rRNA spectra with tied localizations, a fixed target/decoy FASTA, and the original 11/6 ppm, three-missed-cleavage settings. It requires both known tied candidates for scan 43600 to survive and byte-identical BED output across thread counts 1, 4, 4. The fixture is a subset of the Human RNOME Becker experiment, with source-file paths replaced by test metadata. The FASTA was generated once with the original decoy settings so this test isolates search/post-processing rather than decoy regeneration.

The changed C++ code must be compiled and this regression run before claiming that a rebuilt container has reproducible output. Same-build repeatability is distinct from byte identity across compiler/library versions, which can introduce tiny raw-score changes.
