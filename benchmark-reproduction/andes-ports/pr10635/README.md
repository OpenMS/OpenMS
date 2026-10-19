# PR #10635 — ProSE `annotate:isotope_mimic` (port of bigbio/andes#102 and #152)

Question: does what andes reported for its isotope-mimic rule also show up in OpenMS?

andes's claim (bigbio/andes#101/#102, Astral HYE, Percolator 3.7.1, `--refine`): deamidation (+0.984016 Da) is 19 mDa
from one 13C step (+1.003355 Da); "8,042 deamidated PSMs pass at 1%: 14% have a precursor residual near 0, 61% sit at
+12 to +26 mDa, where an M+1 mis-pick lands". andes then **drops** Pass-2 PSMs whose isotope-like modifications an
isotope-step precursor error explains better: accepted deamidations at the mis-pick position 4,869 → 43, Pass-2
(modified) PSMs 9,879 → 3,940, total 50,013 → 44,804 (`--refine` gain over no-refine +21% → +8%). bigbio/andes#152 later
added the same flag as PIN column `IsotopeMimic`; no numbers for it. No entrapment numbers were given.

The PR ports the flag (`isotope_mimic` meta value + Percolator feature, option off by default), not the drop.

## Builds

| build | commit | |
|---|---|---|
| base | `74048ef5` | `git merge-base origin/develop` of the PR |
| pr | `6ceb7553` | PR head |

Release, Ubuntu 24.04 system packages, `-DWITH_THERMO_RAW=OFF`, TOPP tool `ProSE`. `ProSEAlgorithm_test` and the
ProSE TOPP tests that do not need FuzzyDiff pass on the PR head.
andes: `3dded7b` (bigbio/andes#102, rule) and its parent `72ee0a5` (no rule), Rust 1.87.0, `cargo build --release`.

## Data and protocol

The #10364 suite (frozen protocol of the `prose-andes-reproduction` bundle, sha256 `e24cbe74…`), groups that have
deamidation-relevant high-resolution precursors: Velos CID (3), HF-X HCD (3), Astral HCD (3), Lumos HCD LFQ (3);
8,000 MS2 per file prepared by the bundle's `reproduce.py inputs` (hashes in `results/inputs.sha256`, all equal to
`inputs.json`). HYE FASTA `d9ac434d…`; `hy`/`hye` target-decoy databases from `prepare_databases.py`.

- **Search:** frozen `native_dedup:<file>:default` parameters (Trypsin/P, 2 MC, length 7–40, z 2–5, fixed CAM,
  ±20 ppm, isotope errors −2..1, fragment 0.5 Da Velos else 20 ppm, 50 candidates, top 10, deduplicate) with
  variable `Oxidation (M)`, `Deamidated (N)`, `Deamidated (Q)`, at most 2 per peptide. Arms:
  - `off`: `annotate:isotope_mimic=false` (default), `on`: `true`;
  - `nodeam`: frozen modifications (≤1 Oxidation (M)), context only;
  - `post`: `off`, flagged PSMs removed after rescoring;
  - `drop`: the `off` PIN without the rows the `on` PIN flags, rescored again — the andes#102 behaviour (the
    spectrum's next candidate, usually the unmodified peptide, competes instead).
- **Rescoring:** ProSE `-out_pin`, external Percolator from OpenMS/THIRDPARTY (`Linux/x86_64/Percolator`, reports
  3.06.0; the frozen protocol used 3.09.0), `-Y -U --seed 1|42|137`, target PSMs at q ≤ 0.01, mean over seeds per
  file, then over files per group. Per file z = diff / sqrt(var_arm/3 + var_off/3), up z ≥ 2, down z ≤ −2.
- **Entrapment:** every search also on `<db>_entrap_td.fasta` from `make_entrapment.py` (one shuffled, mass-matched
  protein per target, K/R/P kept, r = 1, decoys reversed), combined FDP = 2·N_E/N. Also for Velos (the frozen protocol
  uses natural entrapment there). Per-class FDP and decoy rate (decoys at q ≤ 0.01 / targets at q ≤ 0.01) for
  flagged / deamidated unflagged / no deamidation.
- **andes like-for-like:** frozen #10364 andes command plus `--refine`, on the three Astral files; changed: no mmap
  candidate index (incompatible with `--refine`) and andes-generated reversed decoys on the target-only HYE FASTA
  (`--refine` refuses `--decoy-strategy none`). Same Percolator and seeds.

## Results

### PSMs at q ≤ 0.01 (suite databases, mean of 3 seeds per file, group mean; All = sum of 12 files)

| group | no deamidation | off | on (PR) | post (exclude flagged) | drop (andes) |
|---|---:|---:|---:|---:|---:|
| Velos CID | 2977.6 | 2913.9 | 2925.6 (+0.40%, 0/3/0) | 2460.8 (−15.6%, 0/0/3) | 2954.4 (+1.39%, 1/2/0) |
| HF-X HCD | 4597.6 | 4452.6 | 4445.6 (−0.16%, 0/3/0) | 4027.3 (−9.6%, 0/0/3) | 4525.4 (+1.64%, 2/1/0) |
| Astral HCD | 3049.8 | 2610.8 | 2568.0 (−1.64%, 0/3/0) | 2415.7 (−7.5%, 0/0/3) | 2740.8 (+4.98%, 3/0/0) |
| Lumos HCD LFQ | 5553.6 | 5583.8 | 5578.9 (−0.09%, 0/3/0) | 4570.9 (−18.1%, 0/0/3) | 5593.8 (+0.18%, 2/0/1) |
| **All 12** | 48535.3 | 46683.0 | 46554.0 (**−0.28%**, 0/12/0) | 40424.0 (−13.4%, 0/0/12) | 47443.3 (**+1.63%**, 8/3/1) |

(up/flat/down per file vs `off`.)

### Entrapment (entrapment databases, all 12 files)

| arm | PSMs | combined FDP | est. true PSMs | flagged PSMs | FDP / decoy rate flagged | deamidated unflagged | no deamidation |
|---|---:|---:|---:|---:|---|---|---|
| no deamidation | 46283.0 | 1.43% | 45619 | – | – | – | 1.43% / 0.96% |
| off | 44206.0 | 1.33% | 43616 | 5597.7 | 1.86% / 1.83% | 7.58% / 7.51% | 0.95% / 0.51% |
| on (PR) | 44135.3 | 1.35% | 43541 | 8416.0 | 2.29% / 2.01% | 7.24% / 7.16% | 0.82% / 0.40% |
| post | 38608.3 | 1.26% | 38122 | 0 | – | 7.58% / 7.51% | 0.95% / 0.51% |
| drop | 44846.7 | 1.38% | 44226 | 0 | – | 11.94% / 10.05% | 0.89% / 0.54% |

Per group: `results/per_group.tsv`; per file and seed: `results/per_file.tsv`.

### What the flagged PSMs are (`off`, suite databases)

| group | accepted deamidated / file | residual near 0 (\|err\| ≤ 5 mDa) | at mis-pick (−26..−12 mDa) | flagged | isotope offset −1 / 0 / +1 |
|---|---:|---:|---:|---:|---|
| Velos CID | 664 | 19% | 48% | 68% | 69 / 30 / 1% |
| HF-X HCD | 586 | 13% | 58% | 73% | 74 / 22 / 3% |
| Astral HCD | 285 | 15% | 63% | 68% | 64 / 25 / 7% |
| Lumos HCD LFQ | 1189 | 8% | 79% | 85% | 75 / 19 / 5% |

Residual = peptide − spectrum at the PSM's isotope offset (andes reports the opposite sign). The unmodified form of the
same peptide is among the spectrum's 10 candidates for 85–97% of accepted flagged PSMs (`results/competitor_off.tsv`).
Most flagged PSMs sit at isotope offset −1 (observed one 13C step below the deamidated candidate): the precursor is the
unmodified peptide's monoisotopic peak, matched to the deamidated candidate through the −1 isotope-error window; the
rest are M+1 picks. Their decoy rate (1.6–1.8%) and entrapment FDP (1.9%) are close to the background: the peptide
backbone is right and only the modification state is wrong, which neither decoys nor entrapment can see.

### andes like-for-like (3 Astral files, `results/andes_astral.tsv`)

| andes | PSMs | Pass-2 PSMs | deamidated | at mis-pick |
|---|---:|---:|---:|---:|
| no `--refine` | 3136.5 | – | – | – |
| `--refine`, no rule (`72ee0a5`) | 2468.7 | 66.1 | 43.8 | 47% |
| `--refine`, rule (`3dded7b`) | 2475.9 (+0.3%) | 39.3 (−41%) | 19.3 | 0% |

andes dropped 5,297–5,699 Pass-2 winners per file (targets and decoys). On these 8,000-spectrum files `--refine` loses
PSMs against the no-refine control in both variants.

### Option off and runtime

- `off` on the PR build vs the merge-base build: PIN and Percolator results byte-identical, idXML identical apart from
  date and `search_engine_version`, on all 12 files (`results/isolation_base_vs_pr_off.tsv`).
- The flag in the `on` PIN equals the rule recomputed from the PIN residual on every row; `on` PIN = `off` PIN plus
  the column (`results/checks.tsv`).
- Search time, median of 3 interleaved runs, 4 threads: Astral B1 20.84 s off / 20.05 s on, HF-X A2 15.02 / 14.35 s;
  peak RSS equal (4.18 / 4.55 GB) (`results/timing.tsv`).

## Reproduce

```bash
# bundle: prose-andes-reproduction.zip (sha256 e24cbe74…), then inside it:
python3 reproduce.py inputs astral_A2 astral_B1 astral_B3 hfx_A2 hfx_B1 hfx_B3 lumos_lfq_5192 lumos_lfq_5194 \
  lumos_lfq_5199 velos_125_R1 velos_5000_R2 velos_25000_R3       # needs pyopenms for prepare_databases.py
export S=<work dir> W=<bundle dir> PERC=<THIRDPARTY>/Linux/x86_64/Percolator/percolator
mkdir -p $S/db && cp $W/broad_benchmark/data/{hy,hye}_td.fasta $S/db/
for d in hy hye; do python3 scripts/make_entrapment.py $S/db/${d}_td.fasta $S/db/${d}_entrap_td.fasta; done
# builds installed to /opt/openms-pr and /opt/openms-base (bin, lib, share)
for db in entrap plain; do for arm in off on nodeam; do bash scripts/run_all.sh pr $arm $db; done; done
for db in entrap plain; do for ds in <12 files>; do python3 scripts/make_drop.py $S/runs/pr/$db $ds $PERC; done; done
bash scripts/run_all.sh base off plain
python3 scripts/compare_base.py $S/runs results/isolation_base_vs_pr_off.tsv
python3 scripts/evaluate.py $S/runs results
# andes: binaries of 3dded7b and 72ee0a5 as $S/ext/bin/andes-rule / andes-norule; target-only HYE FASTA
for ds in astral_A2 astral_B1 astral_B3; do bash scripts/run_andes.sh norule $ds; bash scripts/run_andes.sh rule $ds
  REFINE= bash scripts/run_andes.sh rule $ds; done
python3 scripts/evaluate_andes.py $S/runs/andes results/andes_astral.tsv
bash scripts/timing.sh
```
