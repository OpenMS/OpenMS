# #10693, lever 2: Biosaur2 seeding instead of FeatureFinderMultiplex

**Setup.**
- **Build:** develop `5edf8394`, Release, system packages; `diag_link_develop.patch` is applied for the PIP-ECHO arms with the
  spot dump. It only writes the `PIPECHO_DIAG` file and leaves the BSA output byte-identical.
- **Recipe:** #10585 (`plfq.sh`, `-psmFDR 0.01 -proteinFDR 0.01 -picked_proteinFDR true -protein_inference aggregation`),
  with the IDs of `../../pr10637` (Comet/Percolator).
- **Arms:**
  - `-Seeding:algorithm multiplex` (default);
  - `biosaur2` with ProteomicsLFQ's Biosaur2 defaults (`mini` 500, `minlh` 3; Biosaur2 charge 1–6 and m/z 350–1500,
    `Seeding:charge` only applies to multiplex);
  - one identified-only arm (`-targeted_only true`, no seeding).
- **Checkpoints:** each arm uses its own `-feat_dir`. `Seeding:algorithm` is part of the checkpoint fingerprint, checked on
  BSA: combining biosaur2 from multiplex checkpoints exits 11 with "Seeding:algorithm: 'biosaur2' vs 'multiplex'".
- **Run:** `master_seed.sh`. `memwatch.sh` logs the RSS.

Biosaur2 detection was OOM-killed at about 12 GB RSS on this 15 GB machine with 4 and with 2 threads (memory grows across
runs within one process). It ran with 1 thread, resuming from the checkpoints already written. Runtimes are therefore not
directly comparable: multiplex took 3 h 57 min with 4 threads (13.2 GB peak); one run alone with biosaur2 took 31–33 min
with 1 thread (10.8–10.9 GB). Combine time per PIP-ECHO arm was 59 min (multiplex) against 1 h 31 min (biosaur2).

Tables: `results/proteobench_hye.tsv` (raw = `out_msstats.csv`, rt_local = `norm.csv`), `results/spots_*.tsv`
(`../candidates.py` on the PIPECHO_DIAG dumps), `results/transfer_stats_entrap.tsv` (`transfer_stats.py` of
`../../pr10633`, frozen entrapment DB `fb77e3df`).
