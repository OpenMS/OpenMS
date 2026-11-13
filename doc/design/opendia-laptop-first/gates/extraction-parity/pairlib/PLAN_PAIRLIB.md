# PLAN_PAIRLIB: does a library that contains every entrapment's parent fix the entrapment/decoy excess?

Written 2026-11-13, before any search, population, model or score with this library exists.
- Root: `pairlib-20261113/` (local), mirrored at `/nfs/research/juan/singj/opendia/pairlib-20261113` on codon.
- Approved by the user: "start with 1".

## Question
In every arm so far, entrapment target winners exceed what the paired decoys predict by about 1.5× (rho, see
`RESULT_QCAL.md` and `RESULT_INTERF.md`). DIA-NN on our library shows the same excess.
- The diagnosis (`RESULT_QCAL.md`):
  - The ProteoBench 100K FASTA is a random sample of the full FASTA, so 95 % of entrapments lack their human
    parent (the same composition and mass) in the library.
  - Parents that elute give entrapments real, partly matching signal.
  - The mass-shifted decoys never see that signal.
- ProteoBench's own DIA-NN reference searched the full FASTA, where every parent is present, and is calibrated
  (FDP_c 0.0073 at q ≤ 0.01).
- Question: with the same size and recipe, does putting each entrapment's parent in the library bring rho to about 1
  and rule 1c into range?

## Library (the only change of population; label-free construction)
- **Human side:** the 50,317 human entries of the ProteoBench 100K FASTA, unchanged.
  - Their library rows are byte-identical to the current library: 907,440 / 907,440 human transitions, with max
    |Δ| 0 in m/z and intensity and 1.5e-8 in NRT.
  - 36 human transitions are added. They are shared peptides that the current library carried only under an
    entrapment accession.
- **Entrapment side:** each of those entries' own FDRBench shuffle (`_p_target`) from the full ProteoBench FASTA.
  The full FASTA comes in strict target / p_target pairs with identical composition. This was checked for all
  1,420,868 pairs.
- **Rows:** copied verbatim from the full predicted library `full_fasta/...entrapment.decoys.tsv.zst`, the same
  prediction and decoy pipeline as the current 100K library.
  - Protein lists are cut to the selected accessions whose decoy status matches the row.
  - Script: `scripts/subset_pairlib.py`.
- **Size:** 302,480 target precursors (151,245 entrapment), against 300,223 now, and 302,000 decoy precursors.
- **Decoys:** unchanged scheme (pseudo-reverse with K/R swap, mass-shifted), so the comparison isolates parent
  presence.
- **Files and hashes:** `lib/pairlib.sha256` and `lib/lib.sha256` (PQP via TargetedFileConverter; 12-fragment
  prediction via DIALibGen 0.13.1 with the astral-lib12 settings on `pairlib.fasta`).
- **r** = entrapment / human target precursors of this library, read from the PQP by the evaluator.

## Recipe (unchanged from PLAN_XIC_IDS + PLAN_INTERF except option names and compute location)
- **Binaries:** a8be5bed3c, the xic-only OpenDIA; binary hashes are in `container/payload` and `sif.sha256`.
  - `-peak_picking xic -xic_max_peaks 5` becomes `-max_peaks 5` (renamed options, same behaviour: verified
    byte-identical on raw01 for a8be5bed3c).
  - Calibration, seeds, BASE crossfit, population rule (pairs with max transition_support ≥ 3), the 12-fragment
    arrays, F12 and featD are unchanged.
- **Arms:** D11FWT (control recipe) and D11FDWT (+ deficit features), via `train_v4.py <raw> D11F D11FD`, then
  make_w and calib_tail.
- **Compute:** codon, as 2 per-raw Slurm jobs and 1 scoring job.
  - They run in a Singularity image that holds the same binaries, libraries and Python environment at the same
    absolute paths.
  - The codon run root and data are bind-mounted onto the local paths.

## Pass rule: PLAN_XIC_IDS rules, unchanged, applied to D11FDWT on both raws
1. Against this build's BASE:
   - The bootstrap 90 % interval of N@2 % − BASE lies above 0.
   - The interval of N@1 % − BASE is not entirely below 0.
   - **FDP_c at q ≤ 0.01 is ≤ 0.015.**
   - AUC(entrapment > decoy winner) is ≤ the BASE AUC.
2. N@FDP_c 2 % is > 1716 (raw01) and > 1736 (raw02).
3. The null-window winner AUC is ≤ 0.4983 (raw01) and ≤ 0.4986 (raw02).

D11FWT is reported against the same rules but not ruled on.

## Hypothesis test (written prediction, ruled on)
- **H:** the parent-presence diagnosis is right. Then, pooled over both raws and for D11FDWT, rho at the q ≤ 0.05 and
  q ≤ 0.10 cuts falls inside or below its 90 % interval around 1, i.e. the upper bound of the interval is ≥ 1 and the
  point estimate is ≤ 1.2.
- Pooled rho at q ≤ 0.05 above 1.3 rejects H: the excess then has a cause other than missing parents.
- Reference: about 1.6–1.8 in every earlier arm.

## Recorded, not ruled on
- N@FDP_c 0.5/1/1.5/2 %.
- E, D and rho per q cut and stratum (`diagnose.py`).
- The parent-explained share among entrapment winners (`parent.py`).
- Counts of entrapment winners whose parent is also a winner in the same run.
- Wall time and peak memory per step.

## Caveats written in advance
- The N@FDP_c and FDP_c numbers are not comparable one-to-one with the 100K results. The entrapment set differs:
  every entrapment now has a library-resident parent.
- Parents in the library now compete for the same signal. OpenDIA has no cross-precursor interference step yet, so an
  entrapment can still win on its parent's apex if it scores higher. That is a recorded limitation, and the target of
  the next work item (interference graph).
- One run. No retuning after a failure.
