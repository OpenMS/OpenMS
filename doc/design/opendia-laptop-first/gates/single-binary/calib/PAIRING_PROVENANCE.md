# Pairing provenance check, 23 November 2026

No paired FDP statistic was computed on real runs.

## What was checked

- Local search PQP: `/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01/library.pqp`.
- 151,377 entrapment precursors, 151,102 original target precursors. Composition excluding the C-terminal residue plus that residue, charge, and modification multiset give 136,509 unique one-to-one candidate groups. There are 6,847 groups with multiplicity and 135 entrapment-only groups. Sequence composition is not an original peptide_pair_index.
- Searched local benchmark tree for pair/entrapment text and TSV maps. No original FDRBench peptide_pair_index file was found. `samecomp_pairs.py` is a descriptive composition diagnostic, not a generator map.
- Searched the remote OpenDIA tree to depth five. The source artifacts are under `proteobench-astral-entrapment/fasta/`: full FASTA, 100k randomly sampled FASTA, and ZIP. `unzip -l` shows the ZIP contains only the full FASTA, not a pairing TSV.

## Outcome

The prior statement that matching composition establishes applicability of the published paired estimator was too strong and is withdrawn. Complete generator pairing and the sampling/precursor-expansion history are not authenticated. The requested paired estimator is blocked on that prerequisite, rather than filled in with arbitrary matches or a selected subset. Lower/combined curves proceed as the frozen fallback. No claim is made that this bounded search proves the original map does not exist elsewhere.

If the original pairing table becomes available, validate sequence identities, modification treatment, charge expansion, multiplicities, and retained/missing partners before freezing a separate paired-analysis plan. Do not change this plan in response to the observed curves.
