# Native stage-2 with search peak features: raw 01 design run (written before the run)

Target/decoy labels only. The run is fixed by the offline ablation (`../scoring-ablation/ABLATION.md`, `ablate4.txt`), which was done before any native code existed.

- **What runs:** the frozen stage-2 command from `../stage2-design/run.sh`, plus `-stage2_feature_set search -stage2_search_features ../pipeline/search.features.tsv`. That gives the CNN logit plus transition support, log matched scans, log peak width and |library RT error|.
  - The `search` set is used rather than `extraction_search`: offline, `search` gave 717 / 1014 / 1412 and `extraction_search` gave 640 / 1008 / 1349.
- **Expected:** the native run reproduces the offline replica, 717 / 1014 / 1412 at q ≤ 0.01 / 0.02 / 0.05, within ±1% for rounding. The CNN logits must be byte-identical to `../stage2-design` (columns 1–6 and 9).
- **Decision rule:**
  - If the run reproduces the replica and every fold gives finite weights with positive `cnn_z`, the search set is frozen as the stage-2 design candidate.
  - No parameter is tuned on this run.
  - The evaluation on raw 03 needs a new hashed gate and the user's approval.
