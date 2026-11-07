# Local raw trace CNN pilot

This is an offline, CPU-only supervised experiment, not a production identification or FDR pipeline. It emits **logits**, not q-values or calibrated correct-ID probabilities. Only numpy and torch are required. No model or data uploads are performed.

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

## Run locally

From the repository root, first prepare explicit teacher-derived seed windows using the **identical PQP** used by OpenDIA and the scored teacher OSW. `prepare_trace_windows.py` pairs `DECOY_<target TRAML_ID>` with its parent unmodified peptide, grouping charge states/modifications/decoys together. It selects seed-positive targets only inside a rank-one teacher identification's RT bounds at the requested seed threshold. Unconfirmed targets remain unknown. Teacher scores and candidate area can select seed windows, but are never CNN inputs.

```bash
python3 src/topp/OpenDIA/prepare_trace_windows.py --library library.pqp --teacher teacher.osw --candidates candidates.tsv --output seed_windows.tsv --all-windows all_windows.tsv --max-samples 6000 --seed-q 0.001
# Substitute your local build directory for ../OpenMS-build.
../OpenMS-build/bin/OpenDIATraces -in run.mzML -lib library.pqp -windows seed_windows.tsv -out_prefix seed_traces -max_samples 6000 -rt_half_width 15 -im_half_width 0.03 -mz_tolerance 0.02
../OpenMS-build/bin/OpenDIATraces -in run.mzML -lib library.pqp -windows all_windows.tsv -out_prefix all_traces -max_samples 271792 -rt_half_width 15 -im_half_width 0.03 -mz_tolerance 0.02
```

The full export's `max_samples` must cover the actual window count. Export currently allocates its raw tensors with a hard 2 GiB tensor cap, independently of the Python scorer's mmap behavior. Larger exports need tiling before they can use this pilot. Local `OpenDIATraces` is an experimental companion tool, not a stable TOPP registry entry. Build it with `cmake --build ../OpenMS-build --target OpenDIATraces`. Use distinct output prefixes.

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

Heldout AUC measures discrimination of the supplied **seed labels only**. Teacher-selected positives and decoy negatives are a selected subset, not ground truth for all candidate windows. Full same-run scores, including diagnostic overlapping groups, cannot be treated as an unbiased identification benchmark. This pilot does not mine labels, cross-fit, tune hyperparameters, estimate FDR, calibrate probabilities, validate real identification accuracy, or train on a GPU. The checkpoint binds numeric extraction geometry but does not encode instrument identity or exporter software provenance. Synthetic tests establish shape learning with identical total intensities and auxiliary inputs, not transfer to experimental DIA data. Candidate apices still come from the existing window generator. This model ranks windows, it does not yet learn peak boundaries or replace candidate generation.

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
