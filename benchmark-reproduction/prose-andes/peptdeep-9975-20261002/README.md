# #9975 (PeptDeep rescoring features) on develop (2026-10-02)

Questions:
- Does #9975, rebased onto develop `f6c680f`, still work, and what do its five Percolator features
  (`ms2_cosine`, `ms2_spectral_angle`, `ms2_pearson`, `ms2_frac_pred_found`, `rt_abs_error`) give on the
  20-file suite of #10364?
- What do they cost, and what do they do to the error rate?

All numbers: frozen 20-file suite of #10364, 8,000 MS2 per file, `report:top_hits` 10, Percolator 3.09.0
`-Y -U` with seeds 1/42/137, target PSMs at q ≤ 0.01; mean over seeds, then over files. A file is up or down
when z = Δ / SE ≥ 2 in absolute value, SE = sqrt(var_arm/3 + var_base/3) over the seeds. `pd_summary.txt`
holds every comparison with per-file deltas and z; `pd_per_file.tsv` the per-file rows.

## Build

`f5ea2d04` + the ProSE and PeptDeep sources of each arm (`validation-tree-*.patch`), compiled against pyOpenMS
nightly `3.7.0.dev20261001` as before, but with the ONNX variant of the harness (`harness_onnx/`): ONNX Runtime
1.23.2 (develop's vcpkg pin), `WITH_ONNX=1`, and the PeptDeep inference sources compiled in, because the
nightly is built without ONNX. The models are the ones the OpenMS build downloads from
`https://archive.openms.de/openms/models/` (`models.sha256`), found the way an OpenMS build tree finds them
(`OPENMS_BINARY_PATH/share/OpenMS/models`, `harness_onnx/openms_data_path.h`).

Class tests pass for every build (`build_records.json`): HyperScore, FragmentIndex, ProSEAlgorithm,
PeptDeepRescoring (with the end-to-end sections, models present) and PeptDeepInference, whose parity checks
against Python ONNX Runtime pass (largest differences: RT 6e-8, CCS 1.2e-4, MS2 intensity 9.8e-7).

| Arm | Build | Setting | Files |
| --- | --- | --- | --- |
| `re_auto`, `re_auto_E` | `re_on`: develop `f6c680f` ProSE sources | frozen parameters | 20; 9 shuffled-entrapment |
| `pd_off` | `pdr_val`: #9975 rebased (`dcba9c8`) | frozen (PeptDeep off, its default) | 4 |
| `pd_on`, `pd_on_E` | `pdr_val` | `peptdeep:enable=true` (instrument QE, NCE auto) | 20; 9 |
| `pd_inst` | `pdr_val` | `peptdeep:enable=true`, peptdeep's own instrument group | 14 (where it is not QE) |
| `pd_cal`, `pd_cal_E` | `pdr_cal`: + calibration on each spectrum's best target hit (= #9975 `9c5b1b0`) | `peptdeep:enable=true` | 20; 9 |
| `<arm>-nort`, `<arm>-noms2` | none: Percolator rerun on the arm's PIN (`ablate_pd.py`) | without `rt_abs_error`; without the four MS2 features | as the arm |

peptdeep's instrument groups (its `default_settings.yaml`): Velos, Fusion Lumos and Astral → `Lumos`,
Q Exactive HF-X and Exploris 480 → `QE`, timsTOF → `timsTOF`.

## Results

### Yield

| Instrument / acquisition | n | develop `f6c680f` | #9975 rebased | **#9975 + calibration fix** | Sage 0.14.7 | ANDES auto | fix vs Sage | fix vs ANDES |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LTQ Orbitrap Velos, CID | 3 | 2940.0 | 3067.6 (+4.3%) | **3128.2** (+6.4%) | 2919.3 | 3000.6 | +7.2% | +4.3% |
| Q Exactive HF-X, HCD | 3 | 4445.6 | 4685.8 (+5.4%) | **4828.4** (+8.6%) | 4480.9 | 4703.9 | +7.8% | +2.6% |
| Orbitrap Astral, HCD | 3 | 2655.3 | 3049.4 (+14.8%) | **3178.7** (+19.7%) | 2204.4 | 3073.8 | +44.2% | +3.4% |
| Fusion Lumos, HCD, LFQ | 3 | 5552.6 | 5617.8 (+1.2%) | **5629.0** (+1.4%) | 5604.3 | 5616.2 | +0.4% | +0.2% |
| Fusion Lumos, CID, TMT6plex | 3 | 2414.9 | 2411.3 (−0.1%) | **2410.6** (−0.2%) | 2399.8 | 2468.9 | +0.5% | −2.4% |
| Exploris 480, HCD, TMTpro | 3 | 2548.9 | 2566.2 (+0.7%) | **2574.0** (+1.0%) | 2596.3 | 2502.6 | −0.9% | +2.9% |
| timsTOF HT, plasma | 2 | 1041.2 | 1043.7 (+0.2%) | **1057.7** (+1.6%) | 1050.0 | 964.0 | +0.7% | +9.7% |
| **All 20 files (sum)** | 20 | 63754.0 | 66281.7 (+4.0%) | **67362.0** (+5.7%) | 62715.3 | 66025.6 | **+7.4%** | **+2.0%** |

- Files up / flat / down against develop: #9975 rebased 12 / 8 / 0; with the calibration fix 13 / 7 / 0.
  The fix against #9975 rebased: 7 / 13 / 0 (+1.6%).
- Sage and ANDES are the frozen 2026-09-30 measurements, without any predicted-spectrum or predicted-RT
  features of their own. Prediction-based rescoring of their results would also gain.
- TMT groups are flat. Their retention-time calibration stays poor (median residual 730–1930 s) and their MS2
  agreement low (median cosine at the selected NCE 0.46 on Lumos CID TMT, 0.64 on TMTpro; 0.91–0.98 on the
label-free HCD groups, 0.63 on Velos CID).

### Which features carry the gain (Percolator reruns on the same PINs, against develop)

| | all five | without `rt_abs_error` | without the four MS2 features |
| --- | --- | --- | --- |
| #9975 rebased | +3.96% | +4.04% | +0.04% |
| + calibration fix | +5.66% | +4.05% | +2.12% |

As rebased, `rt_abs_error` adds nothing. The calibration fix makes it work.

### Calibration set (`diagnostics/calib_check.py`)

`PeptDeepRescoring` fitted the RT calibration and chose the NCE on the better-scoring half of *all* hits.
With 10 hits per spectrum, only 12.3% (Velos) and 13.9% (Astral) of that set were best-ranked target hits;
43–44% were decoys. The fix (#9975 `9c5b1b0`) takes only each spectrum's best hit, and only if it is not a decoy.
With `report:top_hits` 1 (ProSE's default) it changes only by leaving out the decoys. Median RT residual of
the calibration set before → after (median over the group's files): HF-X 820 → 147 s, Lumos LFQ 2030 → 282 s,
Velos 1207 → 582 s, Astral 139 → 43 s. A new class test checks that adding lower-ranked hits and high-scoring decoys does not
move the calibration; it fails on the unfixed code (12 failed checks) and passes with the fix.

### Error rate

Combined FDP at q ≤ 0.01. Velos: natural UPS1/yeast entrapment, r = 3.38, 3 files × 3 seeds
(`pd_entrapment_velos.tsv`). Others: shuffled paired-mass entrapment, r = 1, seed means summed over files;
estimated true PSMs = N − 2 N_E.

| Design | develop | #9975 rebased | + fix | + fix, MS2 only | + fix, RT only |
| --- | --- | --- | --- | --- | --- |
| Velos (natural) | 1.03% | 1.16% | 1.10% | 1.12% | 1.22% |
| HF-X (shuffled) | 1.35% | 1.47% | 1.73% | 1.45% | 1.93% |
| Astral (shuffled) | 1.40% | 1.82% | 1.55% | 1.61% | 1.44% |
| TMTpro (shuffled) | 1.07% | 1.21% | 1.29% | 1.19% | 1.24% |
| Estimated true PSMs, HF-X / Astral / TMTpro | 12649 / 7031 / 7381 | 13459 / 8133 / 7410 | 13949 / 8430 / 7421 | | |

- The features raise the estimated error rate on every design, by 0.07–0.4 percentage points. The gain in
  estimated true PSMs is larger than the rise in false ones everywhere: HF-X +10.3%, Astral +19.9% with the fix.
- **MS2 features:** decoys and entrapment hits have the same feature medians in every native-score quintile
  (`diagnostics/decoy_vs_entrap.py`), so the reversed decoys are not crudely worse represented by the model.
  The cause of the rise is not understood.
- **RT feature on the shuffled designs:** a shuffled entrapment peptide has the composition of its target, so its
  predicted retention time is that of the true peptide. About 40% of the accepted HF-X entrapment PSMs have the
  unshuffled target among the spectrum's candidates (`diagnostics/sibling_check.py`). Rank-1 entrapment hits
  have a small RT error more often than decoys on HF-X (best-scoring third, entrapment / decoy: 13.7% / 7.7% and
  11.0% / 7.4%), not on Astral (11.7% / 12.2%) and not on the natural Velos design (15.9% / 18.9%,
  17.8% / 15.5%; `diagnostics/rt_tail.py`).
  The HF-X rise with the fix is therefore largely the design. On the natural design the fix lowers the
  error rate (1.16% → 1.10%).
- **peptdeep's instrument groups** (`pd_inst`) raise the Velos estimate to 1.33%.

### Instrument

Against instrument QE (`pd_inst` vs `pd_on`): Velos (`Lumos`) +4.4%, timsTOF (`timsTOF`) +3.0% (per-file
z 0.9 and 1.6), Astral (`Lumos`) −2.5%, Lumos files unchanged. peptdeep's own mapping is not uniformly
better, so `QE` stays the default.

### Cost

Search only (Percolator excluded), OMP_NUM_THREADS=4, two searches at a time on 4 cores (`pd_cost.tsv`).
Median CPU time per file: develop 114 s, #9975 467 s (×3.6, range 2.7–4.4), with the fix 442 s (×3.4).
Peak memory is unchanged. With 10 hits per spectrum PeptDeep predicts every distinct candidate peptide of
every spectrum: 180–450 CPU seconds more per 8,000-spectrum file (median 340–360). develop was measured in an earlier container;
`pd_off` (the same computation, this container) agrees with it within 6% of CPU time on its 4 files. `pd_inst` ran
while the Percolator ablations used the same cores, so its cost is not comparable.

### Identity and TOPP checks
- `pd_off` = develop byte for byte (native TSV and PIN) on 4 files, one per acquisition type (`pd_summary.txt`).
- ProSE TOPP tests replayed through `ProSEAlgorithm::search()` (`topp_emulate.py re_on pdr_val` in
  [split-10335-20261002](../split-10335-20261002)): all 8 idXML outputs identical. PeptDeep is off by
  default and none of its parameters is written to the search parameters.

## Files
- `run.py` (arm table), `chain_pd*.sh`, `build.py --harness harness_onnx`, `harness_onnx/`: the searches and
  the ONNX build, inside the harness of [pr-untangled-20261001](../pr-untangled-20261001).
- `evaluate_pd.py` → `pd_summary.txt`, `pd_per_file.tsv`, `pd_cost.tsv`, `pd_entrapment_velos.tsv`;
  `ablate_pd.py`: Percolator reruns without feature groups.
- `diagnostics/`: calibration set, decoy and entrapment feature distributions, shuffled-sibling and RT-tail
  checks, with their output in `diagnostics.txt`.
- `pd-results.zip`: per job `summary.json`, parameters, resolved settings, logs and per-seed Percolator
  results. `accepted_entrapment_psms.tsv.gz`: accepted PSMs of the shuffled entrapment arms.
