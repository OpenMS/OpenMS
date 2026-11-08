# Gate: frozen Option B design, single evaluation on ProteoBench Astral raw 02 (written before raw 02 is read)

## Frozen design

Code is `feature/opendia-clean-rewrite` at `46be19510d`, built in `OpenMS-r1-full-ninja-build`. The design was set on raw 01 (`astral-design-raw01-20261107/PROTOCOL.md`, `DIAGNOSIS.md`, `pipeline/RESULT.md`) without entrapment labels.
The steps are exactly those in `run.sh`:
1. `OpenDIACalibrate -panel_targets 40000 -window_panel_targets 300000 -out_anchors -out_diagnostics`.
2. `OpenDIA -calibration` with `-min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10`.
3. `prepare_trace_windows.py --anchors --max-samples 6000`.
4. `OpenDIACNN -mode crossfit -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc -rt_calibration -apply_calibration_coordinates -ss_iterations 0 -max_tensor_gib 24`.
   - `-max_tensor_gib` is a resource cap only. Raw 01 needed 16.1 GiB.

Input: `raw/LFQ_Astral_DIA_15min_50ng_Human_02.raw`, read natively. Library: `astral-entrapment-20261107/lib/library.pqp`, unchanged.
The run is done once. `run.sh` refuses a second start (the `STARTED` marker). Raw 03 stays untouched.

## Entrapment labels (library only, `evaluate.py`)

- A target precursor is entrapment when every protein mapped to its peptide has an accession containing `_p_target`. Otherwise it is human, so a shared precursor counts as human.
- Counts from the PQP: 148983 entrapment and **151240** human, out of 300223 targets.
  - The raw 01 gate wrote 151239. The library-only check, made before this run, gives 151240, so r = 148983 / 151240 = 0.9851. The difference is below 1e-5 in r.

## Measured quantities and pass rule (unchanged from the raw 01 gate)

Accepted set: target precursor winners at the tool's pooled paired q ≤ 0.01, with N_h human and N_e entrapment.
- Combined entrapment FDP estimate: FDP_c = N_e (1 + 1/r) / (N_h + N_e). Lower-bound estimate: FDP_lb = N_e / (N_h + N_e).
- **Pass if FDP_c ≤ 0.015.** If there are no accepted targets, or any step fails, the gate is reported as FAIL, as it stands.
- Recorded but not part of the rule:
  - N_h, N_e and the accepted counts at q ≤ 0.001 and q ≤ 0.05;
  - anchor count and anchor entrapment share;
  - runtime and RSS per step.
- If FDP_c is computed on fewer than 100 accepted targets, that is stated next to the result. Its sampling error is large.

## Allowed changes during the run

None to science settings. A crash, I/O or unit failure ends the evaluation as an operational FAIL. It is recorded, and raw 03 then becomes the next candidate under a new gate.

## Limits stated in advance

- Raws 01 to 03 are technical replicates on one instrument, one gradient and one sample. A pass shows transfer across replicates, not across instruments.
- The weak-run design (`-window_panel_targets`, `-pooled_tdc`, `-ss_iterations 0`) was chosen on raw 01. It has not been checked by entrapment before this run.
