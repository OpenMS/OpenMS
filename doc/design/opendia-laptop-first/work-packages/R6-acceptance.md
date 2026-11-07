# R-WP6: Acceptance and merge

[Master plan](../README.md)

**Gates before any PR to `develop`:**
1. Independent acquisitions (500 pg and at least one new held-out run), with every setting frozen beforehand.
2. Entrapment-validated error control (R-WP4).
3. Native Thermo `.raw` and Bruker `.d`.
4. Packaging: LibTorch as an optional dependency, defaults unchanged, CTD/INI docs.
5. Resource budget recorded: wall time, peak RSS and disk, for the whole pipeline.
6. Split the branch into reviewable PRs (search, calibration, CNN, error control).

Inherited from the legacy plan: no conversion requirement, native vendor input, and no claim of parity with DIA-NN or Spectronaut without matched scopes.
