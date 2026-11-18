# PLAN_PGENT2: `-pg_scores` as the default, judged on a homolog-free foreign entrapment (HeLa diaPASEF)

Written 18 November 2026. The plan was written:
- after PRECHECK_FOREIGN, which fixed the homolog rule on D seed 0 of a different library build;
- before this plan's library is predicted and before any search with it.

PRECHECK_FOREIGN saw D seed 0 numbers on a related library (same targets; different entrapment sample and seed). The D
design here is therefore not blind: its expected result is known roughly. The held-out run H has never been searched with
any Arabidopsis library.

## Library (`entlib-ath2-20261119` on codon)
- **Targets:** the 58,998 HeLa library peptides.
- **Entrapments:** Arabidopsis UP000006548 (reviewed) tryptic peptides, 0-2 missed cleavages, length 7-44 (I=L).
  - Excluded: any human tryptic peptide, any substring of a human protein, and any peptide within Hamming ≤ 2 of a
    same-length human tryptic or HeLa library peptide.
  - Length-matched to the targets, seed 20261119.
  - Built by `make_fasta2.py` (hashes in `inputs2.sha256`).
- **Prediction:** DIALibGen 0.13.1 with the PGENT timsTOF settings. Decoys: OpenSwathDecoyGenerator pseudo-reverse.
- **Estimator:** FDP_c = N_e (1 + 1/r) / N, with r read from the PQP.

## Arms and runs
- **CUR:** a2aa8fc56c without `-pg_scores`. **NEW:** the same binary with `-pg_scores`.
- **Container:** `opendia-pg-a2aa8fc56c.sif`. Codon, 8 threads, 96 GB.
- **D:** seeds 0-4, design.
- **H:** seeds 0-4, held out, run only after the design verdict is committed.

## Criteria (on D, all must hold)
1. NEW FDP_c at q ≤ 0.01, pooled over 5 seeds, ≤ 0.012.
2. NEW median N@FDP_c 1 % ≥ 1.05 × CUR median.
3. **Astral (restated, already known):**
   - pooled FDP_c at q ≤ 0.01 ≤ 0.010 (0.0073);
   - N@FDP_c 1 % ≥ CUR on raws 01-03 (+11.8 / +2.4 / +6.2 %).

**Outcome:**
- Pass makes `-pg_scores` the default. The flag and the old path are removed, the commit runs ctests and the Astral
  byte-check is replaced by a NEW-default regression run.
- Fail keeps it opt-in, and nothing is retuned.
- H is reported either way, with the same two numbers and CUR's.

**Caveat:** this benchmark measures false targets that are unrelated to the sample. Near-homolog false matches are
outside it, and are a separate open item (sibling competition).
