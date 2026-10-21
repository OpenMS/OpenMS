# PR #10640: open research questions

PR #10640 (variant knockoffs) was closed after this benchmark: the method is calibrated on PXD022787, but more
research is needed before it goes into OpenMS. The branch `claude/andes-variant-knockoffs` (head 0398d8d9) keeps
the code.

1. **Sensitivity.** At q <= 0.05 one tier keeps 169 of ~1,840 donor-carried variant PSMs that global TDC accepts.
   53,814 of 83,473 variants get no knockoff position and are never accepted. How much can knockoff placement
   (looser context matching, other enzymes' peptides, fewer distance constraints) recover without breaking the null?
2. **Realistic priors.** The tier gain (3.8x at 0.15% FDP) uses genotype-derived tiers, an upper bound. Repeat with
   priors available in practice: population allele frequency, sample RNA-seq/exome calls with realistic error.
3. **Other data.** One TMT dataset with ion-trap MS2. Needed: high-resolution MS2 (Orbitrap HCD, Astral), label-free
   DDA, DIA, and a second genotyped cohort, ideally with all runs.
4. **Small counts.** FDP rests on 1-22 absent PSMs per arm. Bootstrap over runs or larger cohorts for intervals.
5. **Database size.** VariantDatabase writes full-length protein copies (247 M residues here); ProSE ran out of
   memory at 15 GB. A windowed output (as `scripts/trim_variant_entries.py`, 27.8 M residues) is needed.
6. **Comparison.** Against alternatives beyond #10636: mass-shift-aware open search filtering, unrestricted
   modification search of variant spectra, group FDR with mimic-aware decoys.
