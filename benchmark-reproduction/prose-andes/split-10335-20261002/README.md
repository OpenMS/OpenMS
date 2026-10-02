# #10335 split into #10394, #10397 and #10399 (2026-10-02)

Questions:
- Which parts of #10335 carry its gain on develop with the deisotoping fix (#10391)?
- What do its parts give on their own: peptide deduplication (#10394), multiply charged fragments with
  HyperScore (#10397), raw-spectrum retrieval and local fragment evidence (#10399)?
- Can #10335 be closed once these are in?

All numbers: frozen 20-file suite of #10364, 8,000 MS2 per file, Percolator 3.09.0 `-Y -U` with seeds
1/42/137, target PSMs at q ≤ 0.01; mean over seeds, then over files. A file is up or down when
z = Δ / SE ≥ 2 in absolute value, with SE = sqrt(var_arm/3 + var_base/3) over the seeds.
`split_summary.txt` holds every comparison with per-file deltas and z; `split_per_file.tsv` the per-file rows.

## Arms

Each build compiles ProSE from a validation tree, `f5ea2d04` plus the ProSE sources named below
(`validation-tree-*.patch`), against pyOpenMS nightly `3.7.0.dev20261001`; develop itself cannot be linked
against the nightly (#10223). Class tests pass for every build (`build_records.json`).

| Arm | Build: ProSE sources | Setting | Files |
| --- | --- | --- | --- |
| `d_base` | `dev_fix`: develop `3a47278` | frozen parameters | 20 |
| `d_dedup` | `dd_val`: develop `3a47278` + #10394 (= develop `9517361`) | frozen | 20 |
| `d_dedup_off`, `d_dedup_E` | `dd_val` | `peptide:deduplicate=false`; shuffled entrapment databases | 3; 9 |
| `d_fc`, `d_fc_hr` | `fc_val`: develop `9517361` + #10397 (= develop `2f30b40`) | frozen | 6 low-res + 2; 12 high-res |
| `d_35` | `p35_val`: #10335 `136dde3` merged with `3a47278` (not pushed) | frozen | 20 |
| `a35_hs_single`, `a35_hs_multi`, `a35_cal_single` | `p35_val` | `scoring:method` × `scoring:fragment_charges` | 6 low-res |
| `a35_raw`, `a35_raw_local`, `a35_local`, `a35_local_hr` | `p35_val` | #10335's opt-in options | 14 high-res or 6 low-res |
| `a35_local_E`, `a35_raw_local_E` | `p35_val` | the same on shuffled entrapment databases | 9 |
| `re_off`, `re_r`, `re_rl` | `re_val`: develop `2f30b40` + #10399 `0ab4e2b` (both options off by default) | off; raw; raw + local | 3; 6 low-res; 20 |
| `re_auto` | `re_on`: develop `7ae6b88` + #10399 `8475b25` (defaults: `auto`, local on) | frozen | 20 |

`uc_val` (develop `7ae6b88`, with #10398) is the TOPP-replay baseline for #10399.

## Results

### Where #10335's gain came from (low-resolution files, against develop `9517361`)

| | HyperScore, charge 1 | HyperScore, multiple charges (#10397) | calibrated, charge 1 | calibrated, multiple (#10335 default) |
| --- | --- | --- | --- | --- |
| Velos CID | identical | +8.84% | −9.29% | +9.40% |
| Lumos CID TMT | identical | +1.22% | −4.19% | +1.28% |
| Velos entrapment FDP | 1.03% | 1.10% | 1.08% | 1.15% |

On high-resolution files #10335 equals develop with #10394 byte for byte (14/14 files): its only
high-resolution effect by default was deduplication.

### Merged and proposed parts, each against the develop it was added to

| Group | #10394 dedup | #10397 fragment charges | #10399 defaults |
| --- | --- | --- | --- |
| Velos CID | +2.09% | +8.84% | +1.76% |
| HF-X HCD | +0.17% | identical | +1.12% |
| Astral HCD | +1.84% | identical | +14.33% |
| Lumos HCD LFQ | +0.03% | identical | +0.90% |
| Lumos CID TMT | +1.74% | +1.22% | +0.76% |
| Exploris 480 TMTpro | +0.48% | identical | +0.73% |
| timsTOF HT | +1.67% | identical | −0.14% |
| Files up / flat / down | 5 / 15 / 0 | 4 / 16 / 0 | 8 / 12 / 0 |

- **#10399 on low resolution:** raw retrieval loses 39.5% (Velos) and 64.5% (Lumos CID TMT), with every file
  down (`re_r`), so its default `auto` uses it only for high-resolution tolerances. Local fragment evidence
  alone gives the low-resolution gain.
- **#10399 on high resolution:** on Astral, raw retrieval alone gives +3.9% and local evidence alone +9.2%;
  together +14.3%. Raw retrieval gives 1600–1800 more spectra a candidate per Astral file.
- **History over all 20 files** (`history_per_group.tsv`, `history_table.py`): develop `f5ea2d04` 59745.7,
  `3a47278` 60902.3, `9517361` 61409.3, `2f30b40` 62200.3; #10399 63754.0. Sage 62715.3, ANDES 66025.6.

### Error rate

Combined FDP at q ≤ 0.01. Velos: natural UPS1/yeast entrapment, r = 3.38, 3 files × 3 seeds
(`split_entrapment_velos.tsv`). Others: shuffled paired-mass entrapment, r = 1, seed means summed over files.

| Group | develop `3a47278` | + #10394 | + #10397 | + #10399 |
| --- | --- | --- | --- | --- |
| Velos | 0.95% | 1.03% | 1.10% | 1.03% |
| HF-X | 1.35% | 1.49% | (unchanged) | 1.35% |
| Astral | 1.72% | 1.94% | (unchanged) | 1.40% |
| TMTpro | 1.06% | 1.11% | (unchanged) | 1.07% |

- **#10394:** the error rate rises on all four groups. The entrapment count rises in 18 of 27 shuffled and
  5 of 9 Velos file × seed pairs; the cause is not yet understood.
- **#10399 on the shuffled designs:** measured with #10335's build (`a35_raw_local_E`). #10399 equals that
  build byte for byte on all 14 high-resolution files of the normal databases.
- **Local fragment evidence alone** (`a35_local_E`): HF-X 1.48%, Astral 1.71%, TMTpro 1.17%.
- **Raw retrieval on low-resolution data** (`re_r`, `re_rl`): the Velos FDP falls to 0.84% and 0.66%, but
  with 40% fewer PSMs.

### Identity checks (`split_summary.txt`, all pass)
- `d_dedup_off` = `d_base` (3 files); `d_fc` / `d_fc_hr` = `d_dedup` on all 14 high-resolution files;
  `d_fc` = `a35_hs_multi` on all 6 low-resolution files; `d_35` = `d_dedup` on all 14 high-resolution files.
- `re_off` = develop (3 files); `re_rl` = `a35_raw_local` (14 high-resolution files);
  `re_auto` = `re_rl` (14 high-resolution files).

### TOPP and Bruker tests of #10399
- **ProSE TOPP tests** replayed through `ProSEAlgorithm::search()` (`topp_emulate.py uc_val re_on`): no hit or
  score changes. The references were updated with `update_refs_full.py`, which applies only the replay's own
  line changes and allows only the new parameters and features.
- **`TOPP_ProSE_DDA*`** (Bruker HeLa, `ENABLE_OPENTIMS_TESTS`, off in CI): emulated with `run_hela.py`
  (`hela_summary.txt`). Native 1% FDR PSMs change by −0.33%, +0.22% and −0.79%. This emulation counts develop
  itself at 4606 on the calibrated test, below that test's floor of 4650.

## Files
- `run.py` (arm table), `chain_*.sh`: the searches, inside the harness of
  [pr-untangled-20261001](../pr-untangled-20261001).
- `evaluate_split.py` → `split_summary.txt`, `split_per_file.tsv`, `split_entrapment_velos.tsv`;
  `history_table.py` → `history_per_group.tsv`, `history_per_file.tsv`.
- `topp_emulate.py`, `update_refs_full.py`: TOPP replay and reference update; `run_hela.py`: Bruker DDA emulation.
- `split-results.zip`: per job `summary.json`, parameters, resolved settings, logs and per-seed Percolator
  results. `accepted_entrapment_psms.tsv.gz`: accepted PSMs (q ≤ 0.01) of the shuffled entrapment arms.
