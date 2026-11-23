# Label-rule amendment after attempt-2 preparation failure, 23 November 2026

Written and pushed before any search runs. Amends only structural label queries and job wiring. Generation, decoy
method, fragment selection, OpenDIA settings, seeds, q grid, estimators and endpoints are unchanged.

## Attempt 2, retained failure

Jobs 41518555 (GPU prediction), 41518556 (CPU preparation), 41518557 (15-search array), 41518558 (evaluation).

- Prediction COMPLETED: 38:02 wall, 108.6 GiB MaxRSS, CUDA RT/MS2 providers, 8,588,443 target precursors after
  fragment selection, `targets20.pqp` 15.45 GB. Inference took about 2 min; the single-threaded PQP write took most of
  the rest.
- OpenSwathDecoyGenerator inside preparation exited 0: 45:44 wall, 157.5 GiB MaxRSS, 8,573,370 predicted decoys
  (5 dropped, 2,491,338 kept with fewer fragments, 0 model rejections), `predicted20.pqp` 27.7 GB.
- `prepare_library.py` then FAILED after 7:49 with `missing or mixed target entrapment protein mapping`, before writing
  any view. The array was dependency-cancelled, and evaluation failed with no library. No search ran and no ID,
  calibration or quantification outcome exists.

## Diagnosis (read-only queries on the completed `predicted20.pqp`, jobs 41590757, 41599386, 41601408, 41602157)

- The FASTA has 2,841,736 records: 1,420,868 `_p_target` (entrapment) and 1,420,868 `_target`, none `DECOY_`.
  `targets20.pqp` has no `DECOY_` protein.
- In `predicted20.pqp`, 518 peptides are shared by a target and a decoy precursor. A pseudo-reversed decoy sequence
  equals a target sequence (for example decoy of `LQEETER` is `ETEEQLK`, itself a target). The decoy generator reuses
  the PEPTIDE row, so it carries a `DECOY_` protein and `PEPTIDE.DECOY=1`. That affects 859 target and 520 decoy
  precursors, 0.01 % of each side.
- Counting all proteins, 367 targets look mixed (a `_p_target` protein plus a `DECOY_..._target` protein, or the
  reverse). If `DECOY_` proteins are ignored for targets, 0 targets have a missing or mixed label: 4,294,504
  entrapment, 4,293,939 non-entrapment.
- The 859 `inconsistent peptide decoy flags` are exactly those shared peptides. Outside them there are 0.
- The flags OpenDIA actually reads agree everywhere: TRANSITION.DECOY equals PRECURSOR.DECOY for every transition, and
  PRECURSOR.DECOY=1 iff TRAML_ID starts with `DECOY_`. No duplicate TRAML_ID, one peptide per precursor, no shared
  transitions, max 20 fragments.

The frozen validator's intent (each target has one unambiguous entrapment label, flags agree) holds. It failed on
target/decoy sequence collisions, a structure the frozen synthetic fixture did not include.

## Amended rules

1. Target entrapment labels (in `prepare_library.py` validation and audit, `eval_full.py`, and `eval_quant.py`
   population grading) ignore proteins whose accession starts with `DECOY_`. Every target must still map to at least
   one non-`DECOY_` protein and be all `_p_target` or all non-`_p_target`, otherwise preparation fails.
2. The peptide decoy flag check skips only peptides shared by a target and a decoy precursor. Those are counted in the
   audit as `target_decoy_shared_peptides` (peptides, target and decoy precursors). They are reported and never
   filtered, relabelled or removed. Any other mismatch still fails.
3. Collision targets stay in the library and in every count. Their RT and fragments are their own; only the shared
   PEPTIDE row is common. That a target and a decoy have the same sequence (518 peptides) is a known limitation of
   pseudo-reverse decoys, disclosed, not corrected.

New synthetic test `test_shared_target_decoy_peptide_from_decoy_generator` reproduces the structure. It checks pass
with shared rows, unchanged entrapment ratio, and that real flag mismatches and real mixed labels still fail. All 10
synthetic tests pass.

## Execution of the retry (attempt 3)

- Reuse the completed, unmodified `library/targets20.pqp` and `library/predicted20.pqp`. Their SHA256 is recorded
  before the retry in `diagnostics/attempt2_library.sha256`. Prediction and decoy generation are not rerun: same
  inputs, same frozen tools, both exited 0. Rerunning would only add cost and could change CUDA numerics.
- `prepare_views.sbatch` runs only the frozen view construction step of `prepare.sbatch` (same command, same 512 GiB
  memory request) after verifying that hash. `submit_retry.sh` chains views, then the unchanged 15-run array, then the
  unchanged evaluation. Attempt-2 logs, accounting and failed `structure.json`/`views.time.txt` are archived under
  `attempt2/` before resubmission.
- Cost reporting adds the attempt-2 preparation job (53:33 wall, its decoy step is the decoy cost of record) and
  attempt-1 to the totals.
- No further retry is authorized by this amendment. A new failure is reported as is.

## Pre-submission structural dry run (job 41602703, read-only, no search)

The amended `validate` and `audit` on the unmodified `predicted20.pqp` (SHA256 `04112671...a905a912`) pass in 9:35 wall,
2.0 GiB MaxRSS: max 20 fragments; 8,588,443 targets (6-20 fragments), 8,573,370 decoys (1-20 fragments), none without
fragments; 8,573,370 own pairs, 15,073 targets without own decoy, 0 decoys without own target; 2,491,338 pairs with
fewer decoy than target fragments; 527 sibling-count imbalances; shared peptides 518 (859 target, 520 decoy
precursors); entrapment 4,294,504 vs 4,293,939 non-entrapment targets, r = 1.000132. These are structural library
facts, not outcomes. The unpaired targets and fragment imbalance come from the frozen decoy recipe and stay as is.
