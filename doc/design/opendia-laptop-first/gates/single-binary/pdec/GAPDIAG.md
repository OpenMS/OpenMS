# GAPDIAG: where the DIA-NN IDs we miss go, and why (20 Nov 2026)

**Setup:**
- r≈7 HeLa D, seed 0, current default (OpenMS 0a49d8aa8e, codon `entlib-mem4-20261120/runs/D_new`).
- DIA-NN 2.2 report on the same library (`entlib-diann5-20261119/D`).
- Script: `gap_diag.py` and `ppm.py`, codon `gapdiag-20261120` (jobs 38606833, 38609344).
- DIA-NN labels are used for this diagnosis only, never in the tool.

## 1. Where DIA-NN's IDs go

There are 25,881 DIA-NN non-entrapment IDs at q ≤ 1 %. Classes:

| Class | Meaning | n | Share |
|---|---|---|---|
| A | we also identify it | 11,815 | 45.7 % |
| **B1** | **our best window is the right peak (within 6 s of DIA-NN), but it scores too low** | **8,787** | **34.0 %** |
| B2 | a right window exists, but we picked another peak | 2,111 | 8.2 % |
| C | no window within 6 s | 3,168 | 12.2 % |

- In class A, 99.5 % of our IDs sit on DIA-NN's peak.
- **B1 is the gap.** These are right peaks that we score low: 5,346 of the 8,787 have q > 0.5, and only 907 sit between 1 and 2 %. A better threshold or q correction will not bring them in. Their scores have to rise.

## 2. What DIA-NN sees in B1

- **B1 is fainter for DIA-NN too,** but only partly. The AUC of A over B1 is 0.78 on DIA-NN quantity, 0.82 on Evidence and 0.71 on MS1 profile correlation. Even so, in DIA-NN's middle quantity decile 40 % of its IDs are B1 and 36 % are A. So B1 is not just the noise floor.
- **RT and IM are not the cause.**
  - DIA-NN's RT against our RT map has the same median error in A and B1 (5.5 and 5.8 s), and almost none is outside our 41 s window.
  - DIA-NN's 1/K0 falls outside our narrowed window for 1.7 % of A and 3.9 % of B1. That is small.
- **The mass error is larger for the missed classes.** DIA-NN's best-fragment error, after our 0.65 ppm offset, is outside our 8.3 ppm calibrated half width for 2.8 % of A but 8-9 % of B1, B2 and C. The median error is 1.8 ppm for A against 2.8-3.1 ppm for B and C.

## 3. Our view of the same right windows

These are AUCs against the decoys' best windows. On B1 our identity features fall to near-decoy values:

| Feature | A median | B1 median | Decoy median | AUC A | AUC B1 |
|---|---|---|---|---|---|
| x_fragment_coelution | 0.69 | 0.24 | 0.09 | 0.955 | 0.703 |
| f12_support | 5 | 1 | 0 | 0.957 | 0.663 |
| f12_coelution | 0.69 | 0.32 | 0.17 | 0.966 | 0.728 |
| x_im_coelution | 0.76 | 0.38 | 0.20 | 0.935 | 0.701 |
| x_ms1_coelution | 0.79 | 0.39 | 0.16 | 0.843 | 0.609 |
| x_ms2_abs_ppm | 1.16 | 2.63 | 3.90 | | |
| score | 41.2 | 1.41 | 0.69 | 0.995 | 0.590 |

**Our extraction tracks DIA-NN much worse on B1.**

| Spearman, ours vs DIA-NN | A | B1 |
|---|---|---|
| Our MS2 area vs DIA-NN quantity | 0.78 | 0.54 |
| Our MS1 area vs DIA-NN MS1 area | 0.71 | 0.53 |
| Our fragment co-elution vs DIA-NN Evidence | 0.46 | 0.22 |

So on faint peptides our traces carry more foreign signal than DIA-NN's. DIA-NN still sees co-eluting fragments there, and we mostly do not.

## 4. Lead: the trace tensor m/z window

- The 6-fragment trace tensor and its co-elution and evidence features sum all points within a **fixed ±0.02 Da** of each fragment (`TraceCollector(..., 0.02, ...)` in OpenDIAWorkflow.cpp).
- At this library's fragment m/z (10th/50th/90th percentile 344/602/894) that is **±58/33/22 ppm**. The run's calibrated MS2 half width is **8.3 ppm**, and our own 12-fragment grid already uses 17.5 ppm.
- So the trace window is about 4× wider than the data needs. On faint peaks the extra width collects interfering ions and dilutes co-elution. That matches the B1 picture: lower co-elution, a worse area match and a larger mean ppm error.
- **Test:** an opt-in `-trace_tolerance_ppm` (OpenMS e19261618f). At 0 it keeps the old 0.02 Da window, so the default is byte-identical to `newdef` (25/25 ctests).

r≈1.56 D results:

| Arm | Seed | N @ FDP_c 0.5 / 1 / 2 % | Lower / combined bound at own q ≤ 1 % |
|---|---|---|---|
| default | 0 | 11,576 / **13,388** / 14,536 | 0.57 / 1.08 % |
| default | 1 | 13,545 at 1 % | 0.64 / 1.22 % |
| 8 ppm | 0 | 14,087 at 1 % | 0.87 / 1.65 % |
| **10 ppm** | 0 | 13,126 / **14,948** / 16,603 | 0.82 / 1.55 % |
| **10 ppm** | 1 | **14,867** at 1 % | 0.83 / 1.57 % |
| 20 ppm | 0 | 12,695 / 14,253 / 15,812 | 0.81 / 1.53 % |

- **This is the first lever in this session that moves IDs at true FDR:** +11.7 % (s0) and +9.8 % (s1) at FDP_c 1 %, with gains at every FDP level.
- **The tool's own q is looser,** though: the combined bound rises 1.15 → 1.56 %.
- No change in cost: wall time 5:12-5:43 against 5:44, peak about 10 GB.
- Frozen as **PLAN_TRPPM** (r≈7 D, 5 seeds), submitted after the freeze.

## 5. Other findings

- **B2 (8 %)** comes from the W rule picking a wrong peak about 22 s away, with a larger library RT error (20 against 6 s). It is a smaller, separate lever.
- **C (12 %)** is the candidate-picker miss already characterised in PRECHECK_NOCAND (valley merges).
