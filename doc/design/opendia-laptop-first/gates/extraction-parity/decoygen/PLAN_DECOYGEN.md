# PLAN_DECOYGEN: do decoys built like real peptides (own predicted spectrum and RT) remove the entrapment excess?

Written 2026-11-13, before any search with these libraries. Approved by the user ("proceed with your plan; if it works
we implement it in C++"). Root: `decoygen-20261113/` (local), `/nfs/research/juan/singj/opendia/decoygen-20261113` (codon).

## Question
- `pairlib/quirk/RESULT_QUIRK.md` showed two things:
  - In null windows, the current decoys hit random signal 15–25 % less often than any target class.
  - At a matched score, decoy and entrapment winners cannot be told apart.
- The current decoys:
  - are pseudo-reversed with the C-terminal K/R swapped, which shifts precursor mass by ±28 Da;
  - copy their target's predicted fragment intensities by ordinal, and its RT.
- The literature suggests "shuffle-and-predict" for predicted libraries: the decoy sequence gets its own predicted
  spectrum (Chan et al., J Proteome Res 2025). DIA-NN mutates the residues next to the termini.
- **Question:** which decoy construction makes decoys sample false targets correctly (rho about 1), and does it bring rule
  1c into range?

## Libraries (the only change: the decoys; targets byte-identical to PAIRLIB)
Targets are the PAIRLIB target rows, copied verbatim: 302,480 precursors, 151,245 entrapment. One decoy per target
precursor, named `DECOY_<target TransitionGroupId>` as before, so pairing, folds and every downstream script are unchanged.

| arm | decoy sequence (`scripts/make_decoys.py`) | precursor mass |
|---|---|---|
| **A** shuffle | residues 1..n−1 permuted (FNV-seeded), C-term kept. Of up to 100 permutations, the one with the fewest b/y ions within 0.01 Da of the target's (y1 and b(n−1) excluded); all of the first 20,000 reach 0 | same as target |
| **B** mutate | DIA-NN's rule: residues 2 and n−1 mutated GAVLIFMPWSCTYHKRQEND→LLLVVLLLLTSSSSLLNDQE (a modified residue is passed over as DIA-NN does) | changed |
| **C** pseudo-reverse, no K/R swap | residues 1..n−1 reversed, C-term kept | same as target |
| **D** pseudo-reverse with K/R swap | the current rule (sequences equal the PAIRLIB decoys for 98.8 %; the inspected differences are peptides not ending in K/R, where the library decoy has a different C-terminal residue) | ±28 Da |

**Every arm** gets the decoy's own DIALibGen 0.13.1 prediction for its sequence at the target's charge. That covers
precursor m/z, RT, the top 6 fragments (the 6-fragment library) and the top 12 (the 12-fragment export). Settings:
- **6 fragments:** NCE 23, fragment m/z 100–1800. Probe (`probe/`): these settings reproduce the PAIRLIB target
  transitions best, with median cosine 0.999 and the identical 6 fragments for 72 % of targets. The PAIRLIB targets'
  NormalizedRetentionTime is DIALibGen's RT exactly.
- **12 fragments:** the `pred12.parquet` settings: NCE 25, 100–2000, 1–12 fragments.
- **Compute:** the codon A100 (`lib/predict_gpu.sbatch`), 35 s per prediction.

So D differs from PAIRLIB only in its predicted decoy intensities and RT, and A, B and C differ from D only in the sequence
rule. A target loses its decoy when the decoy equals a target sequence (I = L) or another decoy, has fewer than 6
predicted fragments, or falls outside 400–1000 m/z. Counts are in `arm*/home/lib/build_lib.txt` (A: 99.46 % of targets
have a decoy). Library hashes: `arm*/home/lib/lib.sha256`.

## Recipe (PLAN_PAIRLIB, unchanged except the library and one script)
- **Unchanged:**
  - the binaries and the container (a8be5bed3c, the PAIRLIB SIF);
  - `chain_raw.sh` and `chain_score.sh` and every PAIRLIB script (hashes in `arm*/home/FROZEN.sha256`), except one;
  - calibration, search, seeds, BASE crossfit, population, voxels, F12, featD, `train_v4.py <raw> D11F D11FD`, make_w,
    calib_tail, evaluate_v2, null_auc, diagnose and parent.
- **The one replaced script:** `build_pqp12.py`. It now gives each decoy its own 12-fragment prediction instead of its
  target's annotations and intensities. The PAIRLIB version would contradict the arm's decoys.
- **Compute:** each arm's run root `arm<X>/home` is bound onto the PAIRLIB path, giving 8 data jobs and 4 scoring jobs
  on codon (`run_raw.sbatch`, `run_score.sbatch` with the arm as the first argument).

## Pass rule (PLAN_PAIRLIB / PLAN_XIC_IDS rules, unchanged, D11FDWT, per arm, both raws)
1. Against the arm's own BASE:
   - the N@2 % interval lies above 0;
   - the N@1 % interval is not entirely below 0;
   - **FDP_c at q ≤ 0.01 is ≤ 0.015**;
   - AUC(entrapment > decoy winner) is ≤ BASE.
2. N@FDP_c 2 % is > 1716 (raw01) and > 1736 (raw02).
3. The null-window winner AUC is ≤ 0.4983 (raw01) and ≤ 0.4986 (raw02).

## Hypothesis (ruled per arm)
- **H_X:** arm X's decoys sample false targets. Pooled over both raws, for D11FDWT, rho at the q ≤ 0.05 and q ≤ 0.10 cuts
  (E and D summed over the two `DIAG_D11FDWT.txt` files by `scripts/pooled_rho.py`, interval as `diagnose.py`) must:
  - have a point estimate ≤ 1.15;
  - have a 90 % interval that includes 1.
- **Supported** if both cuts meet this. **Rejected** if the q ≤ 0.10 point estimate is > 1.3. PAIRLIB reference: 1.59 and
  1.42.
- **Prediction (written now):**
  - A and C are supported;
  - B is not (its mass shift moves decoys off their target's isolation window, as D's does);
  - D sits between PAIRLIB and A.

## Label-free pre-check (recorded, part of the choice)
- `scripts/null_rates.py`, with no labels: the decoy/target ratio, in null windows (|RT error| > 3σ), of:
  - the share with transition_support ≥ 4;
  - the share with f12_libcos > 0.8.
- PAIRLIB (this script, measured): 0.854 and 0.803 on raw01, 0.830 and 0.771 on raw02. An arm is **eligible** for the
  C++ implementation only if both ratios are within 0.92–1.08 on both raws. This is the deployable check: it needs no
  entrapment.

## Choice for C++ (decided by this rule, written now)
1. Among the arms that are eligible, pass the pass rule on both raws and have H supported: take the one with pooled
   q ≤ 0.05 rho closest to 1. Ties within 0.05 go to the higher summed N@FDP_c 1 %.
2. If none pass the pass rule but some are supported and eligible: take the closest-rho arm, for the decoy generator
   only, and report the rule-1c failure as a separate scorer problem.
3. If no arm is supported, nothing goes to C++ and decoy construction is not the cause.

## Recorded, not ruled on
- N@FDP_c 0.5/1/1.5/2 %.
- rho per cut and stratum.
- The null-rate table split by H/E/D(H)/D(E) (`null_rates.py --labels`, evaluation only).
- Wall time and memory per step.

## Caveats written in advance
- Each arm has a different decoy set, so the BASE crossfit and the scorers are retrained per arm. N differences between
  arms include that retraining; their size is not measured here (the PAIRLIB D11FWT rerun reproduced exactly, so the
  pipeline is deterministic for a fixed library).
- Arms A and C keep the target's precursor m/z. Their decoys read the same MS1 trace and isolation window as their
  target, which is the realistic case for false targets. A decoy that shares fragments with its own target would borrow
  its signal; A is built to share none, and C shares ≥ 3 in 8.6 % of the first 20,000. The C-term y1 and b(n−1) ions are
  always shared, as in FDRBench entrapments (which share 1.06 b/y on average, ≥ 3 in 12 %, first 20,000 pairs).
- One run per arm. No retuning after a failure.

## Amendment 1 (2026-11-13T20:50Z, before any search output exists; no score, population or label seen)
- **What happened:** the first submission (jobs 27897373–27897392) re-ran OpenDIACalibrate per arm. All of it was
  cancelled at the search step.
  - Calibration failed outright for B raw01, C raw01 and C raw02 (exit 8, "0 anchors at the anchor FDR").
  - The others only passed on the 300k-panel retry. Their stage-3 anchors were 225/205 (A), 122 (B raw02) and 515/371
    (D), against 607/661 in PAIRLIB.
- **Why this changes nothing the test measures:** the calibration panel scores targets against the library's own decoys,
  so it varies with the arm. But the calibration is a property of the run and the target library (RT map, m/z and
  peak-boundary models), and those are unchanged. Re-running it per arm would add a confound (different RT maps per arm)
  rather than test the decoys.
- **Change:** every arm uses the PAIRLIB calibration of each raw.
  - `scripts/reuse_cal.py` copies `cal.json` unchanged. It holds no precursor IDs.
  - It remaps the target-only `anchors.tsv` to the arm's precursor IDs by TRAML_ID, and every anchor must map.
  - `chain_raw.sh` skips OpenDIACalibrate. Its only other change is this, so it is the PAIRLIB chain minus the cal step.
  - This is the k11 practice (`astral-commatched-decoys-20261109`, remap_anchors.py).
- **Pass rule, hypothesis, arms and libraries:** unchanged.
- **Recorded as a finding (not ruled on):** with the current calibration code, decoy construction affects calibration
  robustness. The calibrator should not depend on the decoy set this much; noted for the C++ work.
