# Local raw trace CNN pilot

This is an offline, CPU-only supervised experiment, not a production identification or FDR pipeline. Training and frozen scoring emit **logits**, not calibrated correct-ID probabilities. Native cross-fitting additionally provides experimental fold-wise decoy-count q-value estimates, not independently validated FDR control. Native C++ training and inference use optional CPU LibTorch, with no Python or NumPy at analysis time. The original Python implementation remains a developer reference. No model or data uploads are performed.

## Native C++ training and scoring

`OpenDIACNN` implements the same 6141-parameter network and preprocessing using LibTorch's native convolution, autograd and Adam. It streams the acquisition through the existing C++ `TraceCollector`, then passes its raw tensor directly to the native model. No `.npy` export, Python subprocess or Python environment is needed for training or scoring. LibTorch is a native runtime dependency, not a dependency-free handwritten neural framework.

Build this experimental companion explicitly using a compatible CPU LibTorch distribution. It is **off by default**. Other OpenDIA executables do not link LibTorch:

```bash
cmake -S . -B ../OpenMS-build -DOPENDIA_WITH_TORCH=ON -DTorch_DIR=/path/to/libtorch/share/cmake/Torch
cmake --build ../OpenMS-build --target OpenDIACNN OpenDIACNN_test -j1
ctest --test-dir ../OpenMS-build -R '^OpenDIACNN_test$' --output-on-failure

../OpenMS-build/bin/OpenDIACNN -mode train -in run.mzML -lib library.pqp -windows seed_windows.tsv -model native.pt -out seed_native_logits.tsv -max_samples 6000 -threads 5
../OpenMS-build/bin/OpenDIACNN -mode score -in another_run.mzML -lib library.pqp -windows independent_windows.tsv -model native.pt -out independent_logits.tsv -max_samples 300000 -threads 5
# Diagnostic same-run full scoring must explicitly preserve train-group provenance:
../OpenMS-build/bin/OpenDIACNN -mode score -in run.mzML -lib library.pqp -windows all_windows.tsv -model native.pt -out all_native_logits.tsv -max_samples 300000 -threads 5 -allow_training_overlap
```

The window TSV is an explicit input manifest with `sample_id`, `precursor_id`, `group_id`, `label` and `rt_apex`. Use an existing labelled manifest or prepare one separately. Automatic teacher-label generation is not part of this native command. Do not label every target positive. `-mode train` splits related families into train/validation/heldout, fits only labelled training rows, stops using validation loss, restores the best model, and writes a native checkpoint and logits. `-mode score` reads that frozen native checkpoint. Geometry, preprocessing, LibTorch version and family provenance are checked. Checkpoints also record the seed, thread count, epoch cap, patience and batch size. Geometry options default to the export contract below and must match between training and scoring.

Native family assignment uses a versioned, deterministic FNV-1a hash rather than Python's SHA-256 partitioner. Both keep families disjoint, but their particular train/validation/heldout memberships differ. Native and Python training are therefore not expected to produce identical fitted weights. Their checkpoints are separate formats. The developer-only identical-weight parity check loads native parameters into the Python reference:

```bash
uv run --offline --with torch --with numpy --python 3.12 python src/topp/OpenDIA/test/check_native_trace_parity.py --native-test ../OpenMS-build/src/tests/class_tests/bin/OpenDIACNN_test
```

Training normalization and inference allocate only tensor batches, not a second full normalized tensor. The collector still stores the raw tensor, subject to a hard cap set by `-max_tensor_gib` (default **2 GiB**, 6 KiB per window), plus O(N) metadata and logits. Runs that do not fit in RAM need acquisition tiling, not a larger cap. No calibrated probabilities, end-to-end peak boundary learning or unbiased precursor identification counts are claimed.

### Native outer-family cross-fitting

Use `crossfit` to give every scoring window a prediction from a model that excluded
its entire peptide/related-decoy family from **both fitting and early stopping**:

```bash
../OpenMS-build/bin/OpenDIACNN -mode crossfit -in run.mzML -lib library.pqp \
  -windows all_windows.tsv -seeds seed_windows.tsv -max_samples 300000 \
  -model new_crossfit_model -out crossfit.tsv -folds 3 -seed 42 -threads 5 \
  -target_decoy_ratio 1.0
```

Use the searched PQP's target/decoy precursor count ratio, not the observed
candidate or accepted-ID ratio, for `target_decoy_ratio`. For a balanced library it
is 1. All labels in `-windows` are ignored. Only the explicitly sampled `-seeds`
labels enter fitting or early stopping. Each seed must occur in `-windows` with
exactly matching sample ID, precursor ID, family and RT coordinate. Do not label
all target windows positive or change seeds/hyperparameters after inspecting
outer-fold results. The caller still owns correct peptide/decoy family mapping.

The acquisition is streamed once, then three models are fitted sequentially with
the same bounded raw tensor and batch preprocessing. Outer folds use FNV1a64 over
decimal seed, NUL, `outer`, tab and family ID, modulo `folds`. Within each remaining
family set, the original hash's buckets below 6 are training and all other buckets
are validation. The native checkpoints save the outer-fold count/index and every
family's train/validation/heldout role. The command writes
`new_crossfit_model.fold0.pt`, `.fold1.pt`, etc. Existing fold checkpoints are
refused. Full-window inference currently visits all rows per model, but only that
model's outer-heldout predictions enter the final table.

Output contains `sample_id`, `precursor_id`, `group_id`, `is_decoy`, `fold`,
`logit`, `precursor_best`, `q_value`. The maximum-logit window wins independently
for each target or decoy precursor, with lexical sample-ID tie breaking. Within
each fold, exact score ties are processed together and the cumulative estimate is
`min(1, (decoys + 1) * target_decoy_ratio / targets)`. Reverse cumulative minima
give q-values; zero target counts receive 1 and nonwinning windows receive `NA`.
Both target and decoy precursor winners are required in each fold. **Raw logits
from separately trained models are never pooled**, so independently changing a
fold's score scale cannot change its q-values. This is separate-assay decoy-count
estimation, not paired target-versus-decoy competition.

`-paired_tdc` (crossfit only, off by default) switches to paired competition.
It requires `target_decoy_ratio` 1 because pairing balances the classes. A decoy
precursor whose PQP `TRAML_ID` is `DECOY_<X>` pairs with the target whose
`TRAML_ID` is `<X>`. After per-precursor best-window selection, only the stronger
member of each pair competes (exact ties go to the decoy). Unpaired precursors
compete alone, and within each fold `q = min(1, (D + 1) / T)` with the same tie
grouping and reverse minima. Pair losers keep `precursor_best = 1` but get
`q_value = 1`. A pair split across folds or a key shared by two targets is
rejected. Paired members share a family, so they share a fold. This is also an
experimental estimate, with the same caveats as below.

`-rt_window_calibration` (crossfit only, off by default) adds a signal-derived RT
window before precursor selection. For each fold, the 500 target precursors with
the highest best-window logit are anchors at their best window's apex. Fold k's
map is fitted only on the other folds' anchors: 20 equal-count library-RT bins,
median per bin, piecewise linear, constant outside the knots, up to 3 refits that
drop anchors beyond max(3 sigma, 3 s), sigma = 1.4826 x median absolute residual.
Windows with `|apex - map(library RT)| > max(4 sigma, 30 s)` are excluded from
winner selection and q-values and are written with `precursor_best = 0` and an
extra `in_rt_window = 0` column. The log reports each fold's sigma and window.
No teacher identifications enter the map, but the anchors come from the same
teacher-seeded CNN scores.

Cross-fitting removes CNN fitting/early-stopping family overlap. It does **not**
remove dependence introduced by a teacher trained on the whole acquisition,
teacher-selected positive windows, candidate generation, correlated peptides or
unequal decoy detectability. Combining fold-wise discovery sets is not a proof of
global FDR control. These estimates differ from PyProphet's statistical model.
Report them as experimental same-run teacher-seeded results, alongside baseline
overlap and candidate coverage, not as independently validated 1% FDR IDs. A fully
nested teacher or independent-run/entrapment test is still needed. The CNN still
uses fragment traces, without the baseline's MS1 feature evidence.

## Export contract and labels

`<prefix>.traces.npy` must be little-endian float32, C contiguous, shape `[N,6,32,8]`, finite and nonnegative raw matched intensities. Six fragments are sorted by decreasing library intensity. The RT axis covers a configured fixed extent around the candidate apex (default ±15 seconds), not a variable peak width. The IM axis covers ±configured tolerance around library IM. With missing mobility, the exporter puts signal in central IM bin 4 and sets `has_im=0`.

`<prefix>.samples.tsv` has these columns:

```
sample_id precursor_id group_id is_decoy label rt_apex has_im library_intensities fragment_mz
```

Use real tab delimiters. The last two fields contain six comma-separated floats. Missing fragment slots may have m/z zero only when their library intensity and entire raw trace channel are zero. Negative m/z values are rejected. Sample IDs must be unique. **The caller must assign group_id to peptide/related-decoy families consistently across all datasets**, not candidate IDs or independent target/decoy groups. The trainer cannot infer family relationships.

Labels are explicitly supplied seed annotations: `1` positive, `0` negative, `-1` unknown. Positive decoy seeds are rejected. Target status is never converted to a positive label. Training requires both positive and negative labeled training rows and at least one labeled validation row. Unknown rows receive scores but never contribute to fitting, early stopping, or AUC. Seed label quality remains the caller's responsibility.

`<prefix>.geometry.json` is required alongside every tensor/TSV pair. The exporter writes this numeric sampling contract, for example:

```json
{"version":1,"fragments":6,"rt_bins":32,"im_bins":8,"rt_half_width_seconds":15,"im_half_width_vssc":0.03,"mz_tolerance_da":0.02}
```

Version and dimensions must match, and all tolerances must be finite positive JSON numbers (not strings). The complete geometry is saved in the checkpoint. Scoring requires exact geometry equality, including numeric RT, IM and m/z tolerances. Do not edit a sidecar to conceal incompatible extraction settings.

## Python reference workflow

From the repository root, first prepare explicit teacher-derived seed windows using the **identical PQP** used by OpenDIA and the scored teacher OSW. `prepare_trace_windows.py` pairs `DECOY_<target TRAML_ID>` with its parent unmodified peptide, grouping charge states/modifications/decoys together. It selects seed-positive targets only inside a rank-one teacher identification's RT bounds at the requested seed threshold. Unconfirmed targets remain unknown. Teacher scores and candidate area can select seed windows, but are never CNN inputs.

```bash
python3 src/topp/OpenDIA/prepare_trace_windows.py --library library.pqp --teacher teacher.osw --candidates candidates.tsv --output seed_windows.tsv --all-windows all_windows.tsv --max-samples 6000 --seed-q 0.001
# Substitute your local build directory for ../OpenMS-build.
../OpenMS-build/bin/OpenDIATraces -in run.mzML -lib library.pqp -windows seed_windows.tsv -out_prefix seed_traces -max_samples 6000 -rt_half_width 15 -im_half_width 0.03 -mz_tolerance 0.02
../OpenMS-build/bin/OpenDIATraces -in run.mzML -lib library.pqp -windows all_windows.tsv -out_prefix all_traces -max_samples 271792 -rt_half_width 15 -im_half_width 0.03 -mz_tolerance 0.02
```

The full export's `max_samples` must cover the actual window count. Export currently allocates its raw tensors under the `-max_tensor_gib` cap (default 2 GiB), independently of the Python scorer's mmap behavior. Local `OpenDIATraces` is an experimental companion tool, not a stable TOPP registry entry. Build it with `cmake --build ../OpenMS-build --target OpenDIATraces`. Use distinct output prefixes.

Train/score using only cached dependencies:

```bash
uv run --offline --with torch --with numpy --python 3.12 python src/topp/OpenDIA/train_trace_model.py train --prefix /path/seeds --checkpoint /path/model.pt --output /path/seed_logits.tsv --seed 42 --epochs 40
uv run --offline --with torch --with numpy --python 3.12 python src/topp/OpenDIA/train_trace_model.py score --prefix /path/new_data --checkpoint /path/model.pt --output /path/new_logits.tsv
uv run --offline --with torch --with numpy --python 3.12 python -m unittest discover -s src/topp/OpenDIA/test -p test_trace_model.py -v
```

The network uses two 2D convolutions over raw trace images. Preprocessing is `log1p`, divided by each sample's maximum log intensity (zero windows stay zero). A seven-value auxiliary branch contains six sum-normalized library intensities and the `has_im` mask. No IDs, labels, decoy flags, RT coordinates/errors, fragment m/z, peak summaries, or handcrafted scores enter the network.

Group hashes with the supplied seed assign approximately 60%/20%/20% to train/validation/heldout. Allocation is deterministic and independent of row order, not stratified. Very small datasets can lack usable partitions and fail explicitly. Validation loss selects the best checkpoint and stops after six nonimproving epochs by default. The epoch cap is 100, default 40. CPU intra-op threads are fixed to five. The checkpoint records seed, runtime versions, geometry, preprocessing version, sampling assumptions, and all source group partitions. Reproducibility is tested on the same local CPU/runtime, not promised bit-for-bit across architectures or torch versions.

Training writes a checkpoint, `<checkpoint>.summary.json` and TSV `sample_id,group_id,logit,split`. Summary includes per-partition total/labeled/positive sample counts, unknown count, validation loss, epochs, seed, and heldout AUC (`null` if either class is absent). Heldout rows are not used for fitting or model selection. Train and validation predictions are explicitly labeled, never pooled into heldout metrics.

Scoring rejects any training-group overlap by default. This includes the full-window export of the same run used for seed training. Either export an independently grouped scoring dataset, or use `--allow-training-overlap` solely for explicitly marked diagnostic results: overlapping rows still say `train`, not `heldout`. Example for diagnostic scoring of the full same-run export:

```bash
uv run --offline --with torch --with numpy --python 3.12 python src/topp/OpenDIA/train_trace_model.py train --prefix seed_traces --checkpoint model.pt --output seed_logits.tsv --seed 42
uv run --offline --with torch --with numpy --python 3.12 python src/topp/OpenDIA/train_trace_model.py score --prefix all_traces --checkpoint model.pt --output all_logits.tsv --allow-training-overlap
```

Source validation/heldout group membership is preserved, genuinely new groups say `external`. External is not a claim of independent validation: prevent cross-run family leakage and repeated tuning yourself. Never rename groups merely to bypass the overlap guard. Load only trusted local checkpoints.

## Limits

Raw tensors are opened read-only with numpy memory mapping. Validation scans at most 128 windows at once. Normalization and inference use batches of at most 128 windows, not a whole-dataset `log1p` copy. Training and validation loss are also batched. Metadata, seven-value auxiliary arrays, group provenance and one logit per sample remain in memory with O(N) size, so this is bounded tensor working memory, not constant total memory. OS-managed mapped pages may remain cached. This supports scoring large exports such as 271792 windows without allocating a full normalized tensor.

Heldout AUC measures discrimination of the supplied **seed labels only**. Teacher-selected positives and decoy negatives are a selected subset, not ground truth for all candidate windows. Full same-run scores, including diagnostic overlapping groups, cannot be treated as an unbiased identification benchmark. The Python reference does not cross-fit or estimate q-values. Neither implementation mines labels, tunes hyperparameters automatically, calibrates probabilities, validates real identification accuracy, or trains on a GPU. The checkpoint binds numeric extraction geometry but does not encode instrument identity or exporter software provenance. Synthetic tests establish shape learning with identical total intensities and auxiliary inputs, not transfer to experimental DIA data. Candidate apices still come from the existing window generator. This model ranks windows, it does not yet learn peak boundaries or replace candidate generation.

## Local HeLa pilot, November 7, 2026

The published Evosep diaPASEF HeLa run completed the preparation, raw export, training and full frozen-model scoring workflow locally. These are single-run timings, not controlled performance medians:

| Stage | Observation |
|---|---|
| Sample | 3000 teacher-seeded targets and 3000 decoys, 5688 parent-peptide families |
| Training | 28.04 seconds, five CPU threads, 6141 parameters, validation early stop after 14 epochs |
| Group-disjoint seed evaluation | 1183 heldout windows, seed-versus-decoy AUC 0.9885 |
| Raw-zero intervention | Heldout AUC 0.4873 when raw inputs are zeroed with library/mask inputs unchanged |
| Full tensor export | 271792 windows, 81.75 seconds, peak RSS 2443476 KiB |
| Frozen full scoring | 36.29 seconds, peak RSS 2427128 KiB |

All 6000 sampled tensors match the full export exactly. Sample/full frozen logits agree within 2.9e-6. Full scores preserve 20390 training-group, 6664 validation-group, 6679 heldout-group and 238059 previously unseen-family window labels. Do not pool these partitions into an identification count. Raw-zero intervention demonstrates input dependence, not a separately trained auxiliary-only baseline.

Positive seeds came from the existing OpenSWATH/PyProphet SVM baseline at q<=0.001, with candidate apices inside teacher peak boundaries. The reported AUC is agreement with these selected seeds versus decoys, not proof of calibrated 1% precursor FDR or superiority to OpenSWATH. The next scientific gate is independent-run or cross-fitted, family-grouped target-decoy competition with an entrapment/truth check before reporting IDs. No such result is claimed here.

### Native C++ HeLa validation

The same 6000-window manifest and raw acquisition completed native training and full scoring with CPU LibTorch 2.12.1 and five threads, without invoking Python. Single-run observations, not a controlled speed comparison:

| Stage | Observation |
|---|---|
| Native training | 8.53 seconds model work, 23.82 seconds including acquisition, 20 epochs, peak RSS 843 MiB |
| Group-disjoint seed evaluation | 1224 heldout windows, seed-versus-decoy AUC 0.9936 |
| Raw-zero intervention | Heldout AUC 0.4806 with raw inputs zeroed and library/mask unchanged |
| Native full scoring | 271792 windows, 13.98 seconds model/output work, 44.21 seconds including acquisition, peak RSS 2548980 KiB |
| Identical-weight parity | Maximum native/Python logit difference 1.92e-6 across all 6000 sampled real tensors |
| Sample/full parity | Maximum logit difference 2.8e-6, preserving all sample IDs and family partitions |

Full native scores retain 20283 training-family, 6658 validation-family, 6792 heldout-family and 238059 external-family rows. Native family assignments differ from the Python reference, so the AUC difference is not evidence of a better model. All seven focused native/OpenDIA CTests and eleven Python reference tests passed. A separate Torch-disabled OpenDIA build completed the streaming smoke test with byte-identical output. The synthetic identical-weight parity fixture produced exactly matching logits. These checks establish the native port and input dependence, not precursor FDR or identification performance.

### Native crossfit HeLa evaluation

The frozen protocol used the existing 6000 seed windows, seed 42, three outer
family folds, five threads, at most 40 epochs, patience 6 and batch size 32. No
hyperparameters or seed membership were changed after inspecting outer results.
The PQP contained 73922 target and 72866 decoy precursors, giving the explicitly
supplied ratio 1.014492355831252. All 271792 windows received exactly one
outer-heldout prediction, yielding 138542 precursor winners, including 70003
targets. Independent post-run checks matched every family to its PQP parent and
verified that its scoring checkpoint marked it heldout, never train/validation.
All q-values and winner selections were independently recomputed.

| Method on the same HeLa acquisition/library | Target precursors at q <= 0.01 |
|---|---:|
| OpenSWATH + PyProphet SVM, MS1+MS2 | 19132 |
| OpenSWATH + PyProphet XGBoost, MS1+MS2 | 19725 |
| Native fragment-trace CNN, experimental fold-wise decoy counts | 4535 |

The CNN count is **not equivalent evidence of calibrated FDR** to the baseline
q-value model. Fold-specific accepted counts were 1656, 1566 and 1313. CNN overlap
was 3959 precursors with SVM and 3942 with XGBoost; the 576/593 CNN-only precursors
are unvalidated, not demonstrated additional correct IDs. The complete native
crossfit command took 115.63 seconds with peak RSS 2549020 KiB, excluding candidate
generation and teacher/seed preparation. This is one observation, not a timing
median or an equivalent end-to-end baseline speed comparison.

There is a concrete proposal-coverage limitation before scoring: although every
SVM baseline ID had some candidate window, only **6308 of 19132** had a proposed
apex inside that baseline identification's RT boundaries. The CNN selected an
in-bound apex for 6240 of these baseline precursors before the q-value filter.
The corresponding XGBoost numbers were 6270 proposed and 6199 selected. Even the
baseline apex falls inside any candidate tensor's RT extent for only 7850 (SVM)
and 7891 (XGBoost) precursors. These
are agreement-with-baseline diagnostics, not truth labels. The next design gate
is improving fixed-window/apex coverage, then reevaluating raw-trace learning on
an independent or nested-heldout experiment, rather than enlarging the network
and tuning repeatedly on these outer-fold results.

Eight focused CTests now cover outer-family raw/label perturbation invariance,
checkpoint roundtrips, fold-specific score-scale invariance, exact ties, invalid
IDs/classes and CLI input/checkpoint protection. Identical-weight native/Python
synthetic prediction parity remains exact. The real result table, three native
checkpoints and developer-only audit are retained in the local
`raw-trace-cnn-hela/native-crossfit-20261107` benchmark directory. The analysis
command itself ran directly in C++ without Python.

### Multi-apex candidate proposals

Diagnosis of the `-local_peak_split` HeLa candidates showed the coverage loss was
not missing signal. Candidate groups are long (median span 125.6 s, against a
mean baseline peak width of 7.8 s), and 11269 of 19132 SVM IDs had a candidate
group spanning the baseline apex while the single reported global apex lay
elsewhere. Only 13 IDs had no group covering the baseline apex.

`OpenDIA -max_apexes_per_group <n> -min_apex_separation <s>` now reports up to n
separated local maxima of the summed fragment trace per group, strongest first.
Each row shares the group's boundaries, area and support. Separation uses
streaming non-maximum suppression with O(n) state per precursor. The default
`n = 1` is byte-identical to the previous peak and feature tables on this run.

| `-local_peak_split`, separation 10 s | n = 1 | n = 3 | n = 5 |
|---|---:|---:|---:|
| Candidate rows | 271792 | 677810 | 988075 |
| SVM IDs with a candidate apex inside baseline bounds | 6308 | 11280 | 13508 |
| Wall time (s), single observation | 104.9 | 104.7 | 107.2 |
| Peak RSS (KiB) | 811856 | 818632 | 827504 |

These are agreement-with-baseline counts on the same acquisition used to choose
n, not truth labels or an identification result. More rows per precursor also
gives the downstream best-window selection more chances to find a decoy max, so
the CNN must be re-evaluated (new seeds, traces and crossfit) before any change
in accepted precursors is claimed.

#### Re-evaluation on n = 5 candidates

The unchanged seed recipe, which reproduces the original seed and scoring
manifests byte-exactly on `n = 1`, was applied to the `n = 5` candidates: 3000
positive seeds, 3000 decoys, 988075 scoring windows. The crossfit protocol was
frozen (seed 42, 3 folds, 40 epochs, patience 6, batch 32, same ratio), with only
`-max_tensor_gib 8` added. The independent audit (family roles, PQP parents,
winner selection, exact q-values, OSW identities) passed.

| Native CNN crossfit | n = 1 | n = 5 |
|---|---:|---:|
| Target precursors at experimental q <= 0.01 | 4535 | 6262 |
| Per fold | 1656 / 1566 / 1313 | 2101 / 2196 / 1965 |
| Overlap with SVM / XGBoost IDs | 3959 / 3942 | 5611 / 5587 |
| CNN-only vs SVM / XGBoost (unvalidated) | 576 / 593 | 651 / 675 |
| SVM IDs whose CNN-selected apex is in baseline bounds (any q) | 6240 | 12060 |
| Wall time (s) / peak RSS (KiB) | 115.6 / 2549020 | 377.9 / 7090904 |

Accepted targets went up 38%, but this is still well below the 19132 (SVM) and
19725 (XGBoost) baseline IDs. The bottleneck has moved: the CNN now selects an
in-bound apex for 12060 baseline IDs but accepts only 5611 of them at 1%.
Separation between targets and the best-of-more-windows decoys is now the
limiting factor, not proposal coverage. n was chosen on this acquisition, so
this is still not independent evidence.

### Independent acquisition (500 pg HeLa diaPASEF)

Nothing was tuned on this run: n = 5, separation 10 s, the seed recipe and the
crossfit settings are the HeLa ones. Inputs were the
`HeLa02DDM_500pg_5x3_PyDIA` mzML and the IP2 Bruker Human library, subsampled
by a rule fixed before any result (keep peptide family F iff
`fnv1a64("7\0" + F) % 5 == 0`, decoys follow parents): 74261 targets, 58571
decoys. The original 7176-target sampled library was dropped as too small
(SVM 341 / XGBoost 287 IDs at 1%, none at q <= 0.001 for seeds). This mzML stores
1/K0 in a unitless `Ion Mobility` array, so all OpenDIA tools ran with
`-im_array_as_vssc`. Without it, they silently fall back to m/z-only matching.
OpenSwath reads the array either way. The teacher is still same-run.

| Result | Value |
|---|---:|
| OpenSwath + PyProphet SVM / XGBoost, targets at q <= 0.01 | 3008 / 3514 |
| OpenDIA n = 5 candidate rows, mean windows per precursor | 2447410, 23.9 (HeLa: 7.1) |
| SVM IDs with a candidate apex in bounds / CNN-selected apex in bounds | 2820 / 2587 |
| Seeds (SVM q <= 0.001), heldout seed-vs-decoy AUC per fold | 1525 + 1525, 0.997 / 0.992 / 0.993 |
| **Native CNN targets at experimental fold-wise q <= 0.01** | **0** |
| Same at q <= 0.02 / 0.05 / 0.1 | 383 / 656 / 1136 |
| OpenSwath / CNN crossfit wall time, peak RSS | 19.6 min, 19.9 GB / 9.6 min, 17.0 GB |

Proposal coverage carried over (94% of SVM IDs have an in-bound apex), and the
CNN still separates seeds from decoys. It fails the 1% gate. Fold 0's top 25
precursors contain no decoys (24 are SVM IDs), but with fold sizes of about
19000 targets and ratio 1.27, the conservative `(D + 1) * ratio / T` estimate
needs about 127 decoy-free targets before any q reaches 0.01. Decoys enter by
rank 100. With about 24 windows per precursor, best-of-windows decoy maxima are
much higher than on HeLa. The estimator also differs from PyProphet's, so 0 vs
3008 is an upper bound on the gap, not a calibrated comparison. Fixing this
without tuning on the evaluation run would mean making decoy and target
best-of-window selection comparable (paired target-decoy competition or
per-precursor window normalization) and validating that on HeLa first.

### Paired target-decoy competition

Design gate: HeLa only. On the frozen n = 5 HeLa crossfit (same seeds, folds and
hyperparameters, so logits are byte-identical to the separate run), `-paired_tdc`
was compared against the separate `(D + 1) * 1.0145 / T` estimator. 2674 of the
precursor winners had no partner winner and competed alone.

| HeLa n = 5, targets at experimental fold-wise q <= 0.01 | Separate | Paired |
|---|---:|---:|
| Fold 0 / 1 / 2 | 2101 / 2196 / 1965 | 2103 / 2309 / 2119 |
| Total | 6262 | 6531 |
| Overlap with SVM 19132 / XGBoost 19725 | 5611 / 5587 | 5846 / 5816 |
| Native CLI run (4.7 min wall, 7.0 GB peak) | | 6531 |

Paired was at least as high in every fold, so the rule was then applied once,
unchanged, to the existing independent 500 pg crossfit logits:

| 500 pg, targets at experimental fold-wise q <= 0.01 | Separate | Paired |
|---|---:|---:|
| Fold 0 / 1 / 2 | 0 / 0 / 0 | 0 / 110 / 120 |
| Total (SVM 3008 / XGBoost 3514 baseline) | 0 | 230 |
| Overlap with SVM / XGBoost | 0 / 0 | 211 / 213 |

Pairing removes the library-ratio penalty (1.27 on 500 pg) and cancels decoys
that lose to their own target. It does not fix the main 500 pg problem: with
about 24 windows per precursor, decoy best-of-window maxima still compete with
real targets, and fold 0 stays at 0. Paired TDC is a different experimental
estimator, not evidence of calibration. Both numbers carry the same-run teacher
caveats above, and baseline agreement is not ground truth. The next HeLa-gated
step is to make best-of-window selection comparable across precursors with
different window counts.

A window-count correction was then tested offline on HeLa only, with the rule
and gate fixed in advance: precursor score `best_logit - log(n_windows)` (the
max-of-n tail shift for a logistic null), with paired TDC. It had to beat
uncorrected paired TDC in every fold. Folds went 2103 / 2309 / 2119 to
2241 / 2290 / 2220 (total 6531 to 6751, SVM overlap 5846 to 6034). Fold 1
dropped, so the gate failed, and the correction was **not** applied to the
500 pg run or implemented. HeLa has almost equal target and decoy window counts
(7.12 / 7.15 per precursor), so this run cannot strongly test the correction.

### Signal-derived RT window (offline feasibility gate)

Diagnosis: the candidate search ran without an RT prior (`-rt_window 0`), and the
CNN sees only +-15 s of traces, so high-scoring decoys sat far from their
expected RT (HeLa median 43-55 s from a library-RT fit vs 8 s for targets; 500 pg
267-286 s vs about 100 s).

Rule fixed before running, after DIATraceExtractor's anchor calibration and
without teacher IDs: anchors are the top 500 target precursor winners by CNN
logit in each fold; the map for fold k is fitted only on the other folds'
anchors (20 equal-count library-RT bins, median per bin, piecewise linear,
3 rejection rounds at max(3 sigma, 3 s), sigma = 1.4826 MAD). Each precursor
keeps its best window inside `|apex - map(library RT)| <= max(4 sigma, 30 s)`,
then paired TDC. Script: `rt_window_gate.py` (benchmark directory).

| Targets at experimental fold-wise q <= 0.01 | Paired | Paired + RT window |
|---|---:|---:|
| HeLa folds (sigma 2.6-2.7 s, half window 30 s) | 2103 / 2309 / 2119 | 2946 / 2978 / 2630 |
| HeLa total, overlap SVM / XGBoost | 6531, 5846 / 5816 | 8554, 7502 / 7434 |
| 500 pg folds (sigma 50-58 s, half window 198-232 s), applied once | 0 / 110 / 120 | 144 / 173 / 162 |
| 500 pg total, overlap SVM / XGBoost | 230, 211 / 213 | 479, 440 / 447 |

The HeLa gate passed in every fold. The native `-paired_tdc
-rt_window_calibration` crossfit reproduces 8554 exactly, with identical per-fold
in-window counts (4.2 min wall, 7.1 GB peak). On 500 pg the window helps but stays wide:
the same piecewise fit on the teacher's own SVM q <= 0.001 IDs gives sigma
46.8 s (HeLa: 2.4 s), so this library's RTs predict this run poorly and the
CNN-anchor map is already close to the achievable precision. The CNN anchors
come from the same out-of-fold scores, so this is an experimental estimate
with the same teacher caveats. The ID-free v0.13 latent-mixture calibration was
not used: its own handoff marks it NO-GO for narrowing windows, and on the
Astral run01 retrospective proxy the anchor ("auto") calibration was closer for
1018 of 1462 precursors.

### Identification-free run calibration (`OpenDIACalibrate`, experimental)

`OpenDIACalibrate` maps library RT, fragment m/z and 1/K0 to one run and estimates extraction tolerances.
Modern runs often have no iRT peptides, so it uses no IDs, iRT kits or external RTs.

```
OpenDIACalibrate -in run.d -lib lib.pqp -out calibration.json -out_anchors anchors.tsv [-im_array_as_vssc]
```

One streaming pass extracts a panel of 10000 targets (stratified over library RT) and their `DECOY_` partners:
top six fragments, summed within ±15 ppm and ±0.05 1/K0 of the library mobility.

1. **Stage 1, after DIATraceExtractor.** Candidate peaks come from the whole run.
   - Score: mean of the four best fragment-vs-rest correlations × sqrt library cosine. No RT term.
   - Anchors: paired target-decoy competition at q ≤ 0.05, an RT-mode filter, then a binned RT map.
2. **Stage 2, v0.13 latent peak/no-peak mixture.** It runs over every target candidate.
   - Fixed π 0.3, β 0.8, σ = 0.07 × run duration, 11 equally spaced knots.
   - Penalties: the stage-1 map is the prior, plus curvature. The map is monotone with a 0.002 minimum step.
   - **Difference from v0.13:** the source uses SLSQP. This port solves the normal equations, then makes the knots monotone with an isotonic fit, so numerical parity with v0.13 is not shown.
   - The v0.13 prior is a fixed linear map. Here the prior is the stage-1 map instead.
3. **Stage 3.** Targets and paired decoys are rescored within max(4 σ_stage1, 30 s) of the stage-2 map. Anchors are kept at paired q ≤ 0.01.
4. **Stage 4.** On the anchors it fits:
   - the RT map, with tolerance 4 σ;
   - the median MS2 ppm offset, with tolerance 5 robust σ (per-fragment ppm is intensity-weighted over mobility scans);
   - a robust 1/K0 line, with tolerance 4 σ.

The design was set on HeLa diaPASEF only and frozen in `f5ffc1f578`, before the 500 pg run.

HeLa-only design steps:
- **Search windows.** ±50 ppm / ±0.1 1/K0 left only 11 strict anchors. Wide windows admit interference.
- **Correlation score.** All-pair Pearson failed under interference.
- **Stage 1 selection.** Fixed correlation and cosine thresholds were replaced by decoy competition.
- **Panel size.** 2000 to 20000 was checked. 10000 runs in 15 s and uses 1.5 GB.

Teacher files were read only after the freeze, for the post-hoc audit (`calibration-frozen-f5ffc1f578/audit.py`).
The audit compares against same-run OpenSWATH+PyProphet SVM IDs at q ≤ 0.01. Those IDs are a proxy, not truth.

| HeLa 200 ng (design run) | Value |
|---|---|
| Stage 1 anchors (decoys above cut) | 1526 (75) |
| Final anchors (q ≤ 0.01) | 1278; 1130 are teacher IDs, 1115 of those within 10 s of the teacher apex |
| RT map vs 19132 teacher EXP_RT | \|error\| p50 1.8 s, p95 7.3 s; 96.3% within the ±8.8 s tolerance |
| Stage-2 EM map alone | \|error\| p50 6.5 s, p95 19.1 s |
| 1/K0 offset | +0.0216 (teacher median EXP_IM − library: +0.0203); tolerance ±0.033 |
| MS2 | +0.35 ppm, tolerance ±10.7 ppm (no teacher reference) |

**Independent 500 pg run, applied once and unchanged: no calibration.** The run failed closed with 0 anchors at q ≤ 0.01.
- Stage 1 kept only 36 anchors (σ 11.3 s), and stage 3 then had 26 decoys among the top 100.
- This matches the earlier finding: that run's library RTs predict poorly (teacher fit σ 46.8 s), and its signal is weak.
- The tool writes nothing and leaves downstream behaviour unchanged. No setting was changed after this result.

**Limits.**
- Anchor q-values only control the anchor set. They are not an identification FDR.
- The tolerances are fixed multiples of the anchor spread, and on HeLa that spread was measured on the design run itself.
- Downstream m/z and 1/K0 tolerances (OpenDIA, OpenDIATraces, OpenDIACNN) are still the fixed defaults, and the calibrated m/z and 1/K0 offsets are not applied yet. Only the RT map is consumed, by `OpenDIACNN -rt_calibration` below.

#### Downstream RT window from the ID-free calibration (`OpenDIACNN -rt_calibration`)

`-mode crossfit -paired_tdc -rt_calibration calibration.json` keeps only windows whose apex is within the calibrated RT tolerance of the calibrated map.
The window applies before precursor winners and q-values. It cannot be combined with `-rt_window_calibration`.
The calibration uses no CNN score, so the window needs no fold separation. Default output is unchanged.

The offline rule was fixed before running (`coverage-20261107/idfree_window_gate.py`). On HeLa n = 5 with paired TDC, at experimental q ≤ 0.01:

| RT restriction | Targets (folds) | SVM-baseline overlap |
|---|---|---|
| none | 6531 | 5846 |
| CNN-anchor out-of-fold window, max(4σ, 30 s) | 8554 | 7502 |
| ID-free map, calibrated ±8.76 s | **10340** (3355 / 3536 / 3449) | **8875** |
| ID-free map, max(tolerance, 30 s) | 8991 | 7900 |

The native `OpenDIACNN` crossfit run reproduces this exactly: 10340 targets and 73084 / 138542 precursors in window.
This is a decoy-count estimate on the run used to design the calibration, not a calibrated FDR.
There is no independent-run result: the 500 pg calibration failed closed, so there is no map to apply there.

#### Calibrated m/z and 1/K0 coordinates (`-apply_calibration_coordinates`)

With `-rt_calibration`, this flag shifts library coordinates before trace extraction. Tolerances do not change.
- Fragment m/z is multiplied by (1 + calibrated ppm offset × 1e-6).
- Precursor 1/K0 is mapped through the calibrated line, so the ±0.03 mobility window centres on the run's 1/K0.

The HeLa gate rule was fixed before running: accept only if the paired-TDC total exceeds 10340 and no fold drops below 3355 / 3536 / 3449.

The rule passed: total 10679, folds 3358 / 3706 / 3615. Run: 0.35 ppm, 1/K0 −0.0327 + 1.0543 × library.

| HeLa n = 5, paired TDC, experimental q ≤ 0.01 | Targets | SVM-baseline overlap |
|---|---|---|
| ID-free RT window only | 10340 | 8875 |
| plus calibrated m/z and 1/K0 coordinates | **10679** | **9146** |

Both use the design run's own calibration. There is no independent-run evidence: the 500 pg calibration failed closed.

#### Calibrate first, then extract the full library (`OpenDIA -calibration`)

The order above was wrong. The full library was extracted over the whole gradient with fixed tolerances, and the calibration only filtered the candidates afterwards.
OpenSwath calibrates first and then extracts the full library inside the calibrated windows.

`OpenDIA -calibration calibration.json` does that. Before extraction it maps library RT through the calibrated map, fragment m/z by the ppm offset, and precursor 1/K0 through the calibrated line.
It also sets the extraction tolerances from the calibration:
- RT: ±8.76 s around the predicted run RT.
- Fragments: ±10.67 ppm, using the new `SearchOptions::fragment_tolerance_ppm`.
- 1/K0: ±0.033.

These replace `-rt_window`, `-mz_tolerance` and `-im_tolerance`.
`apex_rt_error_seconds` is then relative to the predicted run RT. An invalid or partial calibration file is rejected and no output is written.
Without `-calibration` the output is byte-identical: the default k5 run reproduces `k5.tsv`, `k5.features.tsv` and the seed manifests exactly.

The HeLa gate was declared before running (`raw-trace-cnn-hela/calibrated-extraction-20261107/GATE.md`): accept if the total exceeds 10679.
The pipeline was otherwise unchanged: same seed recipe, crossfit settings and CNN trace geometry, plus `-rt_calibration -apply_calibration_coordinates`.

| HeLa n = 5, paired TDC, experimental q ≤ 0.01 | Calibrate after extraction | Calibrate first |
|---|---|---|
| Candidate rows | 988075 | 198915 |
| OpenDIA extraction wall time | 107 s | 40 s |
| CNN crossfit wall time / peak RSS | 4.1 min / 7.1 GB | 1.2 min / 2.1 GB |
| Targets (folds) | 10679 (3358 / 3706 / 3615) | **11830** (4104 / 4343 / 3383) |
| SVM / XGBoost baseline overlap | 9146 / 9114 | **10158 / 10121** |

The gate passed. Fold 2 fell from 3615 to 3383, although the rule only required a higher total, so this is a mixed result per fold.
The seeds were re-derived from the calibrated candidates with the same recipe, so the CNN training data changed as well as the candidates.
This is still the design run, with the same-run teacher. It is a decoy-count estimate, not a calibrated FDR, and there is no independent-run result: the 500 pg calibration failed closed.

### Teacher-free seeds: calibration anchors plus score iteration

`-ss_iterations <n>` (crossfit only, default 0, off) adds PyProphet-style score iterations after the `-seeds` round.
Each iteration relabels from the previous round's held-out scores and retrains every fold:
- positives: target precursor winners at paired q ≤ `-ss_fdr` (0.01);
- negatives: each decoy precursor's best window;
- capped at `-ss_max_seeds` (6000) per class, in deterministic FNV-1a order.

Outputs and checkpoints come from the last round. With the default of 0, output is byte-identical: `cmp` matched the earlier run with the same seeds.

The HeLa gates were declared before each run (`raw-trace-cnn-hela/anchor-seeds-20261107/GATE.md`), with a pass mark of ≥ 10647, which is 90% of the teacher-seeded count.
Round-0 seeds are the OpenDIACalibrate anchors at q ≤ 0.01 with a containing candidate window (1186), plus equal numbers of decoy windows.
No OpenSWATH output enters training or scoring. The scoring manifest `cal.all_windows.tsv` lists all 198914 candidates; it was written by `prepare_trace_windows.py --teacher`, so its label column holds teacher labels, but crossfit resets every manifest label to -1 and uses only `-seeds`.
`prepare_trace_windows.py --anchors <OpenDIACalibrate -out_anchors TSV>` (with `--anchor-q`, default 0.01) now writes both seeds and manifest without a teacher.
- The positive is the target window whose bounds contain the anchor run RT, preferring the apex closest to it.
- Decoys follow the existing rule.
- The panel-40000 seeds are reproduced byte-for-byte, and the manifest identities match the teacher-written one.
- The `--teacher` output is unchanged (`cmp`).

A teacher-free pipeline is therefore: `OpenDIACalibrate -out_anchors` → `OpenDIA -calibration` → `prepare_trace_windows.py --anchors` → `OpenDIACNN -mode crossfit -ss_iterations 1`.

| Seeds | Targets at paired q ≤ 0.01 | Folds |
|---|---|---|
| OpenSWATH SVM q ≤ 0.001 (teacher) | 11830 | 4104 / 4343 / 3383 |
| Calibration anchors only | 9627 (fail) | 3748 / 3146 / 2733 |
| Script round 2: round-1 q ≤ 0.001 targets vs largest-area decoy windows | 7093 (fail) | 2676 / 2291 / 2126 |
| Script round 2b: round-1 q ≤ 0.01 targets vs best decoy windows | 11332 (pass) | 3629 / 3998 / 3705 |
| **Native `-ss_iterations 1`** | **11118 (pass)** | 3221 / 4164 / 3733 |

The native run (2:42, 2.0 GB) uses the same rule as script round 2b but a different deterministic subset order. The difference from 11332 (−214) is its sampling variability.
The model is trained on the HeLa design run. Like PyProphet's semi-supervised loop, an iteration's labels for fold k come from models that were trained on fold k.
It is a decoy-count estimate, not a calibrated FDR, and there is no independent run yet.

### More calibration anchors: panel scaling and opt-in escalation

On HeLa, the anchor count grows linearly with the calibration panel, at about 0.13 anchors per panel target, and the map stays accurate (`calibration-panel-scaling-20261107/GATE.md`).

| Panel targets | Anchors | Teacher IDs within RT tolerance | Wall time / RSS |
|---|---|---|---|
| 500 / 1000 | fail (0) | n/a | 16 s / 0.8 GB |
| 2500 | 316 | 95.8% | 19 s / 0.9 GB |
| 10000 (default) | 1278 | 96.3% | 27 s / 1.5 GB |
| 40000 | 5479 | 96.9% | 57 s / 4.0 GB |

Small panels fail on the (D+1)/T floor: q ≤ 0.01 needs at least 100 winning targets.
`OpenDIACalibrate -max_panel_targets <n>` (default 0, off) retries a failed calibration with a 4× panel, up to `n`. Each retry is one more raw pass, and memory grows by about 75 KB per pair.
With the default of 0, output is byte-identical. On HeLa, `-panel_targets 1000 -max_panel_targets 40000` escalated to 4000 and calibrated (477 anchors, 97.6% of teacher IDs within tolerance).
This addresses runs that fail for lack of anchors. It does not address the 500 pg failure mode (26 decoys in the top 100), which has not been retried, by design.

Seeding the CNN from the panel-40000 anchors (3000 pairs) followed by `-ss_iterations 1` gave **13611** targets at paired q ≤ 0.01 (folds 4546/4729/4336), with no OpenSWATH input to training or scoring (manifest labels are ignored, as above).
For comparison, panel-10000 anchors gave 11118 and teacher seeds gave 11830 (`raw-trace-cnn-hela/anchor-seeds-p40000-20261107/GATE.md`).
84.7% of these are in the SVM baseline, against 85.9% for the teacher-seeded run, and the 2445 gained are 73.4% in SVM.
Candidates and the RT window still come from the panel-10000 calibration. This is the design run, a decoy-count estimate, with no independent run or entrapment check.
