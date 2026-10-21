# PR #10636 — FalseDiscoveryRate `-PSM_class` (port of bigbio/andes#164)

Question: andes reported that one global 1% PSM FDR lets small classes through at a much higher error rate (on an Astral
file, isoform-specific PSMs at 5.5% class FDR, variants at 5.1%) and that a class FDR keeps them at the nominal
level. Does the OpenMS port show this, measured with an entrapment database, and what does it cost in IDs?

## Builds and arms

| arm | build | what |
|---|---|---|
| baseline | `77cccda7` (merge-base with `develop`) | `FalseDiscoveryRate` without the PR |
| global | `349328c1` (PR head) | `FalseDiscoveryRate`, `-PSM_class` empty (default) |
| class | `349328c1` | the same + `-PSM_class modification` or `-PSM_class <meta value>` |

Release builds with Ubuntu 24.04 system packages (`-DWITH_THERMO_RAW=OFF`), 4 cores, 15 GB RAM. The PR build also had
`claude/andes-entrapment` merged locally (not pushed) for `EntrapmentEvaluator`; that merge touches neither ProSE nor
FalseDiscoveryRate. The search and Percolator are run once per file and seed and shared by all arms; only
`FalseDiscoveryRate` differs.

## Pipeline (`scripts/run.py`)

1. `ProSE` with the frozen #10364 `prose_auto` parameters from `jobs.json` (`-Search:rescore false`, everything else as
   frozen: Trypsin/P, 2 MC, 7–40, z 2–5, CAM fixed, ≤1 Ox(M), ±20 ppm (TMTpro ±50), fragment 0.5 Da Velos/Lumos CID
   else 20 ppm, 50 candidates, top 10).
2. `PercolatorAdapter -seed {1,42,137} -score_type svm -post_processing_tdc` (in-process Percolator, the backend ProSE's
   default rescoring uses since #10578; external Percolator 3.09 is not downloadable here).
3. `annotate.py` adds class meta values to every peptide hit: `bench_mc` (missed cleavages 0/1/2), `bench_z` (2/3/4+),
   `bench_iso` (canonical / isoform-specific: no protein of the hit is a canonical entry), `bench_term` (specific /
   semi-tryptic), `bench_len`. Decoys and entrapment hits carry their target's class (prefixes stripped).
4. `FalseDiscoveryRate -PSM true -protein false -FDR:PSM 1 -algorithm:add_decoy_peptides true` (global), and the same
   with `-PSM_class modification|bench_mc|bench_z|bench_iso|bench_term|bench_len` (class). So class FDR runs on the
   Percolator SVM score, and the global arm is FalseDiscoveryRate's TDC q-value on the same score.
5. `extract.py` / `evaluate.py`: best hit per spectrum, target (incl. entrapment), q ≤ 0.01. Entrapment FDP = 2·N_E/N
   (combined estimator, r = 1; a PSM is entrapment when all its proteins are `ENTRAP_`). Per file the mean over the
   three seeds, per group N and N_E summed over files. A file is up/down when |z| ≥ 2 (z = Δ / sqrt(var_a/3 + var_b/3)).
   `decoyFDR` = accepted decoys / accepted targets within the class (andes's metric); only meaningful for the global
   arm, because ClassFDR raises target q-values only.

Cross-check: `EntrapmentEvaluator -class_meta_value bench_mc` gives the same N, N_E and FDP as `evaluate.py` on
`hfx_A2` seed 1 (`results/crosscheck_entrapmentevaluator_hfx_A2_s1.tsv`; with `-database` T/E = 0.952 instead of 1).

## Data and settings

- Spectra: the frozen #10364 suite, 20 files × 8,000 MS2 (`reproduce.py inputs` via `scripts/inputs_driver.py`, all
  selected-mzML hashes verified against `inputs.json`). FASTA `d9ac434d…` → `hy`/`hye`/`human` (`prepare_db_nopyopenms.py`,
  same sha256 as `databases.json`).
- Entrapment databases (`make_entrapment.py`, the #10364 supplement script): targets + one shuffled, mass-matched
  `ENTRAP_` protein per target + reversed `DECOY_` of both. `hye`/`human` are byte-identical to the frozen ones
  (`fb77e3df…`, `76d25071…`); `hy` `9f118291…`.
- `std`: frozen parameters, all 20 files. Classes: modification, missed cleavages, charge.
- `hardmods`: + Deamidated (N, Q) and Phospho (S, T, Y), ≤2 variable mods (`database:chunk_size 10000` for memory),
  HF-X, Astral, TMTpro, Velos (12 files). Unenriched samples, so phospho PSMs are a small, mostly false class.
- `iso`: HYE + 22,140 Swiss-Prot human isoforms (UniProt 2026_03, `UP000005640_9606_additional.fasta`, `sp|` only),
  entrapment by `make_entrapment_iso.py` (segments shuffled with a seed from their own sequence, so peptides shared by
  a canonical protein and its isoforms stay shared in the entrapment). Astral, HF-X (6 files). Class canonical vs
  isoform-specific: the andes setting.
- `semi`: `peptide:enzyme_specificity semi`, HF-X, Astral, Velos (9 files). Class specific vs semi-tryptic.
- `immuno` (second dataset, see below): id_bruker MHC class I, unspecific, 1 file.
- Follow-up `ext`: `std` with external Percolator (OpenMS THIRDPARTY 3.06, `PercolatorAdapter -use_subprocess true`;
  3.09 is not downloadable here), all 20 files. This route is above nominal overall, like the frozen protocol.
- Follow-up `isofull`: `iso` on whole runs instead of the 8,000-spectrum subsets (`fetch_full.sh`, archive `.mzML.gz`
  sha256 = `inputs.json`): `hfx_A2` (111,413 MS2), `astral_A2` (123,621), `eclipse_tmtpro_10855` (33,013, against
  human + isoforms, `humaniso_entrap_td.fasta`).

## Results

Target PSMs at q ≤ 0.01 (seed means summed over files), entrapment hits, entrapment FDP; global → class.
`results/per_group.tsv` has every group and class, `results/per_file.tsv` every file.

| setting / class | global PSMs, N_E, FDP | class PSMs, N_E, FDP | Δ PSMs |
|---|---|---|---|
| **std, 20 files, all PSMs** (`-PSM_class modification`) | 59,718.7, 183.3, 0.61% | 59,464.3, 173.7, 0.58% | −0.43% (0/19/1 files up/flat/down) |
| std, Oxidation | 4,202.0, 11.0, 0.52% | 3,994.0, 3.3, 0.17% | −4.9% |
| std, Oxidation, Astral | 168.7, 3.0, 3.56% (decoys 2.57%) | 118.0, 0, 0% | −30% |
| std, Oxidation, Velos | 72.7, 1.0, 2.75% | 13.7, 0, 0% | −81% |
| std, all PSMs (`-PSM_class bench_mc`) | 59,718.7, 183.3, 0.61% | 59,056.0, 164.3, 0.56% | −1.11% (0/13/7) |
| std, 2 missed cleavages | 632.7, 10.0, 3.16% | 419.0, 2.7, 1.27% | −34% |
| std, 1 missed cleavage | 7,034.7, 36.3, 1.03% | 6,585.7, 24.7, 0.75% | −6.4% |
| std, all PSMs (`-PSM_class bench_z`) | 59,718.7, 183.3, 0.61% | 59,519.3, 176.0, 0.59% | −0.33% (0/19/1) |
| std, charge 4+ | 2,655.7, 9.3, 0.70% | 2,580.0, 6.3, 0.49% | −2.9% |
| **hardmods, 12 files, all PSMs** | 33,150.3, 144.3, 0.87% | 32,368.7, 80.7, 0.50% | −2.36% (0/2/10) |
| hardmods, Phospho | 209.0, 16.3, **15.6%** (decoys 8.0%) | 47.3, 0, 0% | −77% |
| hardmods, Deamidated | 4,511.3, 47.3, **2.10%** | 4,026.7, 14.7, 0.73% | −11% |
| hardmods, Deamidated+Phospho | 43.7, 12.3, 56% | 0 | −100% |
| hardmods, unmodified | 26,336.3, 64.7, 0.49% | 26,325.3, 64.3, 0.49% | −0.04% |
| **iso, 6 files, all PSMs** | 20,286.7, 85.3, 0.84% | 20,205.3, 81.3, 0.81% | −0.40% (0/5/1) |
| iso, isoform-specific | 53.3, 1.3, **5.0%** (decoys 3.75%) | 32.3, 0.3, 2.1% | −39% |
| iso, isoform-specific, Astral | 25.7, 1.0, 7.8% (decoys 7.8%) | 15.7, 0.3, 4.3% | −39% |
| **semi, 9 files, all PSMs** | 28,730.7, 108.0, 0.75% | 27,980.0, 101.3, 0.72% | −2.61% (0/3/6) |
| semi, semi-tryptic | 1,892.7, 7.0, 0.74% | 1,142.0, 0.3, 0.06% | −40% |
| semi, semi-tryptic, Velos | 105.7, 2.3, 4.4% | 0 | −100% |
| **ext (external Percolator), 20 files, all** (`modification`) | 63,875.3, 447.3, **1.40%** | 63,424.3, 405.3, 1.28% | −0.71% (0/15/5) |
| ext, Oxidation | 4,542.0, 44.7, 1.97% | 4,211.7, 11.0, 0.52% | −7.3% |
| ext, unmodified | 59,333.3, 402.7, 1.36% | 59,212.7, 394.3, 1.33% | −0.20% |
| ext, all (`bench_mc`) | 63,875.3, 447.3, 1.40% | 63,079.7, 395.0, 1.25% | −1.25% (0/11/9) |
| ext, 2 missed cleavages | 691.0, 17.3, 5.0% | 552.3, 3.7, 1.3% | −20% |
| ext, 1 missed cleavage | 7,629.0, 80.7, 2.1% | 7,033.3, 44.0, 1.25% | −7.8% |
| ext, all (`bench_z`) | 63,875.3, 447.3, 1.40% | 62,951.0, 385.7, 1.23% | −1.45% (0/13/7) |
| ext, charge 4+ | 2,867.3, 42.7, 3.0% | 2,470.0, 10.7, 0.86% | −14% |
| **isofull, 3 whole runs, all** | 112,602.7, 525.0, 0.93% | 112,540.3, 518.7, 0.92% | −0.06% (0/3/0) |
| isofull, isoform-specific | 310.3, 6.3, **4.1%** (decoys 3.0%) | 248.0, 0, 0% | −20% |
| isofull, isoform-specific, Astral | 100.0, 3.7, 7.3% | 69.0, 0, 0% | −31% |
| isofull, isoform-specific, HF-X | 140.3, 1.7, 2.4% | 117.0, 0, 0% | −17% |
| isofull, isoform-specific, TMTpro | 70.0, 1.0, 2.9% | 62.0, 0, 0% | −11% |

With external Percolator (`ext`) the global FDR is above nominal (1.40%), and 63,875 PSMs match the frozen develop
reference (63,754, no entrapment). Class FDR fixes the small classes (Oxidation 1.97 → 0.52%, 2 MC 5.0 → 1.3%, charge
4+ 3.0 → 0.86%) but not the large one: unmodified PSMs stay at 1.33%, so the overall FDP only falls to 1.23–1.28%.
That excess is in Percolator / the global estimate, which a class FDR cannot correct.

The in-process route (in-process Percolator, FalseDiscoveryRate TDC on the SVM score, entrapment database) is calibrated below
nominal overall (0.61% FDP on `std`; the #10629 session measured 1.37–1.41% on Astral with external Percolator 3.09 on
the frozen PINs), and gives 59,719 PSMs on the 20 files against 63,754 for the frozen develop reference (no entrapment,
external Percolator). The comparison here is between the two FDR arms on identical Percolator output.

Per-class entrapment counts are small (N_E of 0.3–16 per class); the FDPs of small classes carry large uncertainty.

### Option off = baseline

FalseDiscoveryRate of `77cccda7` and of the PR without `-PSM_class` give byte-identical idXML on all 20 files (seed 1,
`results/identity_check.tsv`).

## andes reference

- andes#164 (PR description; one Astral file, UniProt isoforms, variants and chains, Percolator; class FDR = class
  decoys above the lowest accepted score / accepted class targets, from andes's decoy output; no entrapment):

  | accepted at 1% | total | canonical | isoform-specific | chain | variant |
  |---|---:|---:|---:|---:|---:|
  | global q-value only | 46,927 | 46,561 | 183 (5.5%) | 85 | 99 (5.1%) |
  | class FDR (#164) | 47,108 | 46,884 (0.97%) | 92 (0%) | 72 (0%) | 60 (1.7%) |

  (The rendered table has one column fewer in its header than in its rows; read as above, totals add up.) The text
  also gives isoform-specific 5.5–6.7% and variant 2.9–5.1% at global 1%; the commit message `fabb9e0` says 4–6%.
  With target-decoy competition inside a class, 85 chain PSMs with no decoy above them still got q ≈ 1.2%, hence the
  pooled-null Storey estimator. On that file the guard sent the canonical class to target-decoy inside the class.
  The total rose by 0.4% (canonical +323), which q = max(global, class) cannot do by itself: the two rows must come
  from runs with different global q-values (e.g. separate Percolator runs), so the total change is not attributable
  to the class step.
- OpenMS, `iso`: isoform-specific PSMs at global 1% have 3.75% decoys / 5.0% entrapment FDP over 6 files (Astral 7.8%),
  class FDR brings them to 2.1% at −39% of that class (andes: 183 → 92, −50%, 5.5% → 0%) and −0.4% overall (andes
  +0.4%, see above). Class calibration: same direction, similar size (reproduced, on far fewer PSMs). Total: OpenMS
  loses PSMs, as the max rule implies; andes's +0.4% is not a class-FDR effect.
- OpenMS, `isofull` (3 whole runs, 6.3 entrapment hits in the class): isoform-specific PSMs at global 1% have 4.1%
  entrapment FDP and 3.0% decoys (Astral 7.3%); with `-PSM_class` 248 of 310 remain with no entrapment hit (−20%;
  andes −50%, 5.5% → 0%), −0.06% overall. This is the closest match to andes's experiment, and the clearest
  reproduction.
- andes itself on the same spectra (`scripts/run_andes.sh`, andes `fa0195a`, UniProt human reference proteome XML with
  `--uniprot-isoforms`, Percolator 3.06, `results/andes_reference.txt`): 9–18 isoform-specific PSMs per file at
  global 1% with no decoys among them, on 3 Astral (HYE mix, human-only database) and 3 TMTpro files; andes keeps 0–13
  of them after its class FDR. The effect andes reported is not visible on these files with andes either: too few
  isoform-specific PSMs.

## Second dataset: immunopeptidomics

- Data: archive.openms.de `testfiles/id_bruker/classI/high/DN12_Liver_classI_techRep1.mzML` (PXD058436, MHCquant2,
  timsTOF Pro), 8,000 MS2 with selected-ion m/z chosen like `prepare.py` (smallest SHA256 of dataset + native ID,
  `scripts/subsample_mzml.py`, streamed; charges 1–5 kept).
- Search: `human_entrap_td.fasta`, `enzyme 'unspecific cleavage'`, `enzyme_specificity none`, `low_memory true`, length
  8–12, z 1–4, ±20 ppm, fragment 30 ppm, no fixed modification, ≤1 Ox(M) (`IMMUNO` in `run.py`). Classes: peptide
  length, charge, modification. For unspecific digestion the entrapment keeps protein composition but does not pair
  peptide masses.
- Result: 506.3 PSMs at q ≤ 0.01 with 1.7 entrapment hits (0.66%). `-PSM_class bench_len` removes a third of the
  8-mers (45.0 → 30.3) and nothing else (−2.9% overall); charge and modification classes change nothing. With under
  two entrapment hits per run this dataset cannot show a miscalibrated class or its correction: inconclusive.
  The #10630 session saw a similar yield on this data (about 720 PSMs per file, charge 2–5).

## Verdict

- Calibration inside hard classes: **reproduced**. Classes that a global 1% lets through at 2–56% FDP (phospho and
  deamidation in a search with more modifications, 2 missed cleavages, charge 4+, oxidised PSMs on Astral/Velos,
  semi-tryptic PSMs on Velos, isoform-specific PSMs) come down to ≤ 1.3% with `-PSM_class`; overall FDP drops too
  (0.87 → 0.50% in `hardmods`, 1.40 → 1.23–1.28% with external Percolator). Isoform-specific PSMs on whole runs:
  4.1% → 0% at −20% of the class, matching andes (5.5% → 0%).
- A class FDR does not fix a miscalibrated large class: with external Percolator the unmodified PSMs stay at 1.33%.
- Cost: −0.3 to −1.1% PSMs on the standard suite (−0.7 to −1.45% with external Percolator), −2.4% with more variable modifications, −2.6% semi-specific.
- The estimator is conservative where a class is already calibrated: semi-tryptic PSMs at 0.74% FDP lose 40% (to
  0.06%). In classes whose targets are mostly false (Velos Oxidation: ~510 targets against ~450 decoys per file), the
  Storey null fraction is clipped at 1, the class q-values become Benjamini–Hochberg over all class targets, and few
  pass (Velos Oxidation −81%, Astral Oxidation −30%). These classes all used the pooled-null estimator; the
  target-decoy fallback was used only for the fully tryptic class of the semi search on the HF-X files (its decoys are
  stronger than the pool, as the guard intends) and for one hardmods class with fewer than 20 decoys.

## Reproduce

```
python3 -I scripts/prepare_db_nopyopenms.py <bundle>/broad_benchmark          # after fetching hye_targets.fasta
python3 -I scripts/inputs_driver.py <bundle>                                   # the 20 frozen mzML
python3 -I scripts/make_entrapment.py <bundle>/broad_benchmark/data/hye_td.fasta db/hye_entrap_td.fasta   # hy, human alike
# iso: cat hye targets + sp| isoforms > hyeiso_targets.fasta; make_entrapment_iso.py hyeiso_targets.fasta db/hyeiso_entrap_td.fasta
OPENMS_BIN=... PROTO=<bundle> DBDIR=db WORK=work python3 -I scripts/run.py std     # hardmods/iso/semi/immuno: see chain.sh
python3 -I scripts/evaluate.py work results
```
