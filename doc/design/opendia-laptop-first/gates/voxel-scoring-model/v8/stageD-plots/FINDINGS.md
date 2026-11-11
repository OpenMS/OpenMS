# Stage-D misses in the raw data (raw 01, D11FWT, old population): are they real weak peaks or false peaks?

PDF: `stageD_examples_raw01.pdf`, made with `plot_examples.py`. It shows 17 precursors chosen by fixed rules (seed 11).
- **Data sources:**
  - fragment and MS1 XICs from the OpenDIA `.xipm` export, read with pyopenms `XIPMParquetFile`;
  - the apex MS2 scan from `raw01.mzML` (FileConverter), read with `OnDiscMSExperiment`;
  - plots drawn with pyopenms-viz 1.2.0.
- **Native-sampling check:** `apex_scan_evidence.py` → `APEX_SCAN_EVIDENCE.txt` and `APEX_SCAN_MATCHED.txt`.
  - It uses the raw scans whose isolation window holds the precursor. Astral 2 Th windows revisit every 1.65 s.
  - It counts library fragments within 10 ppm in the apex scan and its neighbours.

## What "low abundance" means here

- **Scale:** DIA-NN log10 Precursor.Quantity in this run has median 5.72 and 10th–90th percentile 3.74–6.92.
- **Definition:** "low" = the lowest tertile of the DIA-NN reference IDs we evaluate, log10 quantity < 5.57.
- **D medians:** 5.30 for D against 6.37 for accepted (E).
- **XIC scale:** in the examples, the fragment XIC apex is about 10³–10⁴ counts for D, against 10⁵ for strong E.

## They are real peaks, not false peaks (`APEX_SCAN_EVIDENCE.txt`)

| Group | n | Score median | Library fragments (of 6) in apex scan, median | Share with ≥ 4 in apex scan | Share with ≥ 2 scans (apex ± 3 s) holding ≥ 3 fragments |
|---|---|---|---|---|---|
| E accepted (sample) | 600 | 42.2 | 6 | 1.00 | 0.98 |
| **D missed** | 866 | 7.0 | **4** (best of 3 scans: 5) | **0.75** | **0.80** |
| Decoy winners, score in the D interquartile range | 9266 | 3.2 | 2 | 0.10 | 0.15 |
| Entrapment winners, same range | 4769 | 3.3 | 2 | 0.13 | 0.18 |

- The pages show this directly for D:
  - co-eluting fragment XICs at DIA-NN's RT, with our apex on it;
  - an MS1 M+0/M+1 peak at the same RT;
  - 4–6 library fragments in the apex scan.
- **Decoy and entrapment pages** show something different:
  - scattered, non-co-eluting traces;
  - 1–2 fragments in the apex scan;
  - no MS1 peak.

## But at the same score, decoys look almost as good (`APEX_SCAN_MATCHED.txt`)

- **Score-matched comparison** (strata at D score deciles):
  - fragments in the apex scan (best of 3): AUC 0.68 for D against decoys, 0.51 for entrapment against decoys;
  - points across the peak: 0.57 against 0.47;
  - TIC share: 0.44 against 0.49.
- **The decoys that reach D's score have the same raw evidence:** 4–5 fragments and 2–3 points.
- **The cause is the score cut.** The 1 % cut (22.0) sits at the top 0.1 % of about 31,000 decoy winners per fold. A real but weak peak with 4 of 6 fragments and 2–3 points looks like the luckiest decoys. Decoys get those 4–5 fragments by chance coincidence in a dense 2 Th Astral window.

## About a third of D score below most decoys (`APEX_SCAN_LOWD.txt`, `LOWD_MODEL.txt`)

- 287 of 866 D precursors score below the decoy q90.
- They still hold 4 fragments in the apex scan (decoys in the same score range hold 2), with 63 % having ≥ 2 good scans against 11 %.
- **Their model inputs are target-like:**

| Model input (median) | Low D | Decoys at same score |
|---|---|---|
| f12_libcos | 0.83 | 0.57 |
| boundary_library_correlation | 0.40 | −0.06 |
| \|RT error\| (s) | 22 | 37 |

- The W rule keeps the window with the best linear part. On that kept window, the **full** (nonlinear) D11F score is at the decoy 6th percentile, but the linear part is at the decoy 94th percentile. So the network's nonlinear part pulls these windows down.
- **Linear-only diagnostic** (`LOWD_LINEAR.txt`, not a proposal): ranking by the linear part alone recovers 40 / 15 old-D precursors. Overall N drops from 1496 to 1181 at 1 %, so the full model is better overall.

## Did the linear RT extension work?

- **Partly.** It did what it was meant to do: A+B misses fell from 351 / 331 to 230 / 203, and the new end region is decoy-clean.
- **It did not raise accepted IDs at 1 %.** Most recovered precursors are weak and land in D.
- **PLAN_V8 failed on FDP_c** (0.016 / 0.019). The excess entrapment sits inside the old RT range, not in the new ends.

## How these could be recovered (ideas; none tested)

1. **More fragments and better per-fragment evidence.** D has 4 of 6 library fragments. Extra library fragments (the 12-fragment set) and fragment-level scoring at native sampling (each scan, not 32 bins over ±15 s) would raise the count that separates real from chance. The binning smooths a 3 s peak into 1–2 bins.
2. **Fewer chance matches for decoys.** Decoys reach 4–5 fragments by chance. A per-fragment null rate from local spectrum density would make "4 fragments" mean more for a sparse scan than for a dense one.
3. **The nonlinear part pulls low D down.** About a third of D get a much lower full score than their linear part. Training on more decoy-like negatives may be over-penalising weak true windows. A calibration or regularisation check of the network on weak targets would test this.
4. **MS1 on the same peak.** Most D pages show an M+0/M+1 trace at the same RT. MS1 features did not help before (V4), but they were averaged over 32 bins.
