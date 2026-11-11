# Voxel v3-seq result: ReDeeM-style sequence encoder added to the single decoy-whitened network

Plan: `PLAN_V3_SEQ.md` (sha256 `d46262125a…`, written 2026-11-10T17:06:02Z, sha check OK after the run). Run notes and deviations: `PLAN_V3_SEQ.notes` (N1–N6, I1). Data: ProteoBench Astral raws 01 and 02. Entrapment labels are used for evaluation only.

## Verdict: no arm passes

| Arm | Rule 1 (PLAN_V2 1–4 against BASE) | Rule 2 (null-winner AUC ≤ control + 0.005) | Rule 3 (N@2 % above control, both raws) | Result |
|---|---|---|---|---|
| DS11 (11 features + sequence) | FAIL: ΔN@2 % on raw 02 is −39..+247; FDP_c at q ≤ 0.01 on raw 02 is 0.018 | PASS (0.5016 / 0.5006 against 0.5015 / 0.5022) | FAIL (1544 / 1535 against 1572 / 1572) | FAIL |
| DVS (voxel + 11 features + sequence) | FAIL: ΔN@1 % on raw 02 is −968..−236; ΔN@2 % includes 0 on both raws; FDP_c 0.016 / 0.020 | PASS (0.5011 / 0.4997 against 0.4996 / 0.4994) | FAIL (1028 / 1055 against 1517 / 1438) | FAIL |

The sequence input does not add IDs at FDP_c 2 % on either raw. D11 (PLAN_V3_P2b, 11 features, no sequence) stays the best single network.

## Numbers (precursor target winners after pooled paired TDC; FDP_c = N_e (1 + 1/r) / N; decoy-count estimates)

| Raw | Arm | N@0.5 % | N@1 % | N@1.5 % | N@2 % | q ≤ 0.01 N (FDP_c) | Entrapment AUC | Null-winner AUC |
|---|---|---|---|---|---|---|---|---|
| 01 | BASE | 900 | 1024 | 1267 | 1314 | 1081 (0.013) | 0.531 | 0.518 |
| 01 | D11 | 1256 | 1348 | 1496 | 1572 | 1504 (0.015) | 0.511 | 0.5015 |
| 01 | DS11 | 1318 | 1449 | 1513 | 1544 | 1423 (0.010) | 0.509 | 0.5016 |
| 01 | DV | 1254 | 1388 | 1447 | 1517 | 1497 (0.019) | 0.500 | 0.4996 |
| 01 | DVS | 779 | 904 | 1003 | 1028 | 1009 (0.016) | 0.502 | 0.5011 |
| 02 | BASE | 1215 | 1288 | 1333 | 1376 | 1241 (0.008) | 0.535 | 0.523 |
| 02 | D11 | 980 | 1334 | 1417 | 1572 | 1353 (0.012) | 0.514 | 0.5022 |
| 02 | DS11 | 1252 | 1337 | 1394 | 1535 | 1422 (0.018) | 0.503 | 0.5006 |
| 02 | DV | 1027 | 1291 | 1394 | 1438 | 1131 (0.005) | 0.498 | 0.4994 |
| 02 | DVS | 331 | 443 | 941 | 1055 | 1025 (0.020) | 0.504 | 0.4997 |

Paired bootstrap against the control (1,000 precursor resamples, seed 5, 90 % interval; `VS_CONTROL.txt`):

| Raw | Comparison | ΔN@1 % | ΔN@2 % |
|---|---|---|---|
| 01 | DS11 − D11 | +101 (−115..+200) | −28 (−120..+88) |
| 02 | DS11 − D11 | +3 (−133..+373) | −37 (−195..+75) |
| 01 | DVS − DV | −484 (−611..−307) | −489 (−574..+2) |
| 02 | DVS − DV | −848 (−975..−258) | −383 (−847..+289) |

Against BASE (`EVALUATION_SEQ.txt`): DS11 ΔN@1 % +145..+556 / −58..+147, ΔN@2 % +174..+393 / −39..+247.

## Observations (not part of the verdict)

1. **No sequence tell.** On both raws the null-winner AUC of both arms stays within 0.002 of 0.5 and within 0.0015 of its control. The encoder did not learn "target-like sequence" over null windows, even though the library sequence alone separates targets from decoys (AUC 0.537, `SEQUENCE_TELL.txt`).
2. **DS11 is noise-level against D11.** All four bootstrap intervals include 0. The +101 at 1 % on raw 01 does not repeat on raw 02 (+3). The entrapment AUC goes down a little (0.509 / 0.503 against 0.511 / 0.514).
3. **DVS collapses when the folds are pooled, not inside the folds.** Per-fold N@2 % for DVS (531 / 454 / 516 on raw 01, 556 / 576 / 512 on raw 02) is close to DV (548 / 540 / 506, 609 / 556 / 522). But the pooled curve drops to 1028 / 1055. After the per-fold decoy median/MAD normalisation, the decoy-winner tails differ between folds: the q99.99 is 185 to 1598 on raw 01 and 456 to 1933 on raw 02. So one fold's high-scoring decoys set the cut for the other folds. DS11 also has wider tails than D11 (q99.99 up to 2117 against 670). A probable cause is the larger, more variable embedding when the sequence code is added, but this is not proven. Median/MAD sets only the scale near the decoy bulk, not the extreme tail. This is a pooling weakness that the 11-feature arms hid. A future plan should calibrate per fold on the decoy tail (for example, per-fold decoy-quantile or per-fold TDC), and must freeze that before any run.
4. **Compute.** Round 0 to 2 per fold: DS11 about 40 min, DVS about 25 min (the voxel arm uses 50 steps per epoch). DS11 about 2 h per raw, DVS 1.2–1.4 h per raw (one process, 10 threads), with a peak RSS of 2.1 / 4.8 GB. Encoding each distinct precursor once per batch, with length buckets, gave the same values (max difference 2.7e-7) at about 15× less cost (N4).

## Deviations and run notes (from `PLAN_V3_SEQ.notes`)

- N1: storage only. The dense per-position atom tensors (15 GB per raw, disk 96 % full) were replaced by a table of the 62 distinct residue graphs plus per-position ids. The encoder computes the same values.
- N2: the ported ReDeeM Trp template has 10 C (real Trp has 11). It was kept as in ReDeeM, and it applies equally to targets and decoys.
- N3: the laptop suspended from 17:09 to 21:47 during the first run (the process was paused and then stopped in N4). No scores came from that run.
- N4–N6: the exact speed-ups. Running folds in parallel thrashed (50× slower) and was stopped before any score was written. The final run was sequential, with the seeds unchanged (100 + k, 200 + 10k + r). The plan said "about 10 min" for DS11. It took 2 h per raw.
- I1: interim raw-01 observation before raw 02 finished. The chain and the plan were not changed.

## Files

`PLAN_V3_SEQ.md` (+ `.sha256`, `.written_at`, `.notes`), `seq_featurize.py`, `train_seq.py`, `chain_seq.sh`, `train_{DS11,DVS}_raw0{1,2}.log`, `EVALUATION_SEQ.txt`, `NULL_SEQ.txt`, `VS_CONTROL.txt`, `vs_control.py`, `null_auc.py`, `NULL_BASELINE.txt`, `sequence_tell.py`, `SEQUENCE_TELL.txt`.

## Next step: needs a user decision

Options:
1. Promote D11 (11 features, decoy-whitened single network, no voxel, no sequence) to C++ `OpenDIAScoringModel` under a new frozen plan. It is the only arm that meets all PLAN_V2 numbers on both raws (it was not a named candidate in PLAN_V3_P2b).
2. Add a frozen per-fold tail calibration (observation 3), and then re-test DS11 / DVS / DV. The current pooling may hide a part of the per-fold signal of the larger arms.
3. Stop the sequence encoder line. It gives no measurable gain at 2 % on these two raws.
