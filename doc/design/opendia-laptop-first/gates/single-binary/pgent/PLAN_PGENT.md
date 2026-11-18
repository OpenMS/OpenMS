# PLAN_PGENT: `-pg_scores` as the default, judged on true FDR (Astral entrapment + a new HeLa entrapment library)

Written 18 November 2026. The plan was written:
- after the library was built;
- before any search with it, so no ID count, FDP or score on it has been seen;
- with the Astral numbers below already known from PLAN_PGSCORES (they are restated, not new);
- after reading Q7_NEXT.md, which changed the Astral criterion from per-raw to pooled.

The pooled Astral criterion is therefore partly post hoc on Astral data. The new evidence in this plan is HeLa.

## Library (HeLa entrapment, `entlib-hela-20261118/` on codon)
- **FASTA:** for each of the 58,998 unique unmodified target peptides of the HeLa experimental library, one entry
  `<pep>_target` followed by `<shuf>_p_target`.
  - The shuffle is the same peptide with every K and R and the C-terminal residue kept in place, seeded by sha256 of
    the peptide.
  - It is rejected if it equals a target or an earlier shuffle. All 58,998 got one; none were dropped.
  - Built by `make_fasta.py` (sha256 in `fasta.sha256`).
- **Prediction:** DIALibGen 0.13.1 on codon, CPU (the CUDA provider did not load). Settings:
  - timsTOF instrument defaults: NCE 30, and CCS converted to 1/K0 with the instrument's mobility correction;
  - precursor m/z 400-1000, charges 2-4, lengths 7-44, 2 missed cleavages;
  - Oxidation (M) at most once, Carbamidomethyl (C) fixed, no N-terminal Met excision;
  - fragments 1-6 with fragment charge ≤ 2, no decoys (`-generation:decoys none`).
  - Targets and entrapments are predicted the same way. The experimental library's empirical RT, IM and intensities are
    **not** used, so parents and entrapments are on equal footing.
- **Decoys:** OpenSwathDecoyGenerator `-method pseudo-reverse`, the same as the reference HeLa library's decoys.
- **Entrapment label:** a precursor is an entrapment if every protein it maps to ends in `_p_target`.
  - r = entrapment / human target precursors, read from the PQP.
  - FDP_c = N_e (1 + 1/r) / N, the combined estimator of Wen et al. 2025, as on Astral.

**Caveat, stated before any result:** a predicted library is not the experimental one. The ID counts are not comparable
with the 16,116 / 24,000 numbers of the experimental library. CUR and NEW are compared on the same predicted library.

## Arms and runs
- **CUR:** OpenDIAWorkflow at a2aa8fc56c without `-pg_scores`.
- **NEW:** the same binary with `-pg_scores`.
- **Container:** `opendia-pg-a2aa8fc56c.sif`. OpenDIAWorkflow inside it is byte-identical to the PLAN_PGSCORES test
  binary.
- **Runs:** HeLa D (design) and H (held out), self-calibrating, codon, 8 threads.
  - Seeds 0-4 on D, both arms.
  - H is run only after the design verdict is committed.
- **Astral:** reuse the PLAN_PGSCORES / PLAN_MS1 runs on raws 01-03, seeds 0-4. Nothing new is run.

## Criteria (all must hold)
1. **HeLa D, true FDR:** NEW FDP_c at q ≤ 0.01, pooled over the 5 seeds' entrapment counts and IDs, ≤ 0.012.
   - It must also be ≤ CUR's pooled FDP_c + 0.003, so NEW is not materially worse calibrated than today.
   - Seeds of one run are not independent, so this is a per-run point estimate, judged on its value. It is not a
     significance test.
2. **HeLa D, IDs at equal true FDR:** NEW median N@FDP_c 1 % ≥ 1.05 × CUR median.
3. **Astral calibration:** NEW FDP_c at q ≤ 0.01 pooled over raws 01-03 (seed 0) ≤ 0.010. The value is already known:
   0.0073.
4. **Astral IDs:** NEW median N@FDP_c 1 % ≥ CUR median on each of raws 01-03. The values are already known: +11.8 %,
   +2.4 % and +6.2 %.

**Outcome:**
- Pass makes `-pg_scores` the default. The flag is removed, and the old path is deleted.
- Fail keeps the flag opt-in, and nothing is retuned.
- H is reported either way: FDP_c and N@FDP_c 1 % for both arms.

Criteria 3 and 4 are restated known results and carry no new evidence. The decision rests on 1 and 2.
