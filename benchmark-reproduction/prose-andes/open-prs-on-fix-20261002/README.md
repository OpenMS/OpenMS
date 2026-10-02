# Open ProSE PRs on develop with the deisotoping fix (2026-10-02)

Question: which gains of the open ProSE PRs still hold once develop carries #10391, the deisotoping fix?

## What was compared

Every arm is compared with develop `3a47278`, which carries #10391; the earlier comparison (2026-10-01) is
each PR against develop `f5ea2d04` without the fix.

| Arm | Source | Setting | Files |
| --- | --- | --- | --- |
| `d_base` | develop `3a47278` ProSE sources | frozen parameters | 20 |
| `d_priors` | #10378 rebased onto `3a47278` (`6bef29f`) | `annotate:self_trained_ion_priors=true` | 20 |
| `d_79` | #10379 rebased onto `3a47278` (`2b7c8f2`) | `calibration:enabled=auto` (its default) | 14 high-res + 2 low-res |
| `d_79_mass` | same | `auto` + `scoring:method=mass_accuracy` | 14 high-res |
| `d_35` | #10335 `136dde3` merged locally with `3a47278` (not pushed) | frozen parameters | 20 |
| `d_priors_E`, `d_79_E` | as above | shuffled paired-mass entrapment databases (r = 1) | HF-X, Astral, TMTpro |

Protocol as before: frozen 20-file suite of #10364, 8,000 MS2 per file, Percolator 3.09.0 `-Y -U` with seeds
1/42/137, target PSMs at q ≤ 0.01; mean over seeds, then over files. A file is up or down when
z = Δ / SE ≥ 2 in absolute value, with SE = sqrt(var_arm/3 + var_base/3) over the seeds.

**Harness.** develop cannot be linked against the pyOpenMS nightly (#10223), so each build compiles ProSE
from a validation tree: `f5ea2d04` plus the arm's ProSE sources (`validation-tree-*.patch`; `dev_fix` is
develop's own ProSE sources, identical to `3a47278` on every compiled file). Class tests pass for all four
builds (`build_records.json`). The #10378 and #10379 validation trees also replay the ProSE TOPP tests to
their merged references exactly (`cmp_refs.py`).

**Identity checks** (`checks_d.txt`, all 42 pass):
- `d_base` gives byte-identical native TSV and PIN to `ds1` (the fix on `f5ea2d04`) on all 14 high-resolution
  files, and to `base` on the 6 low-resolution files.
- `d_79` equals `d_base` on `velos_125_R1` and `lumos_tmt_5058` (`auto` does not calibrate low-res data).
- `d_priors` has the same native TSV as `d_base`, and the same PIN once the `ion_prior_*` columns are removed.

## Results (`holds_summary.txt`, `holds_per_file.tsv`)

Group change in Percolator PSMs at q ≤ 0.01; files up / flat / down of the files searched.

| Group | #10378 old → new | #10379 default old → new | #10379 + mass accuracy old → new | #10335 old → new |
| --- | --- | --- | --- | --- |
| Velos CID | +5.5 → +5.5% | 0 → 0 | — | +11.7 → +11.7% |
| HF-X HCD | +0.1 → +0.4% | −0.9 → −0.5% | −0.9 → −0.8% | +0.2 → +0.2% |
| Astral HCD | +3.5 → +6.2% | +5.5 → +4.9% | +6.5 → +4.2% | −0.2 → +1.8% |
| Lumos HCD LFQ | −0.2 → −0.1% | −0.4 → −0.3% | −0.4 → −0.4% | +0.1 → +0.0% |
| Lumos CID TMT | +3.1 → +3.1% | 0 → 0 | — | +3.0 → +3.0% |
| Exploris 480 TMTpro | +0.8 → +0.4% | +1.9 → +0.8% | +2.2 → +0.6% | +0.3 → +0.5% |
| timsTOF HT | +0.7 → +0.5% | 0 → 0 | −0.6 → +0.2% | −1.3 → +1.7% |
| Files up / flat / down | 8/11/1 → **9/11/0** | 2/15/3 → 3/15/2 | 5/6/3 → 2/8/4 | 4/16/0 → **5/15/0** |

- Low-resolution files are not deisotoped, so every Velos and Lumos CID TMT number is unchanged.
- #10379's Astral gain comes mostly from `astral_B3` (+10.5%); `astral_A2` and `astral_B1` gain +1.1 and +2.8%.

## Entrapment (`entrapment_d.tsv`, `entrapment_d_seeds.txt`)

Combined FDP 2·N_E/N at q ≤ 0.01, where N and N_E are seed means summed over the group's three files.
Baseline: develop with the fix (`hds1_E` / `tds1_E`; `d_base_E` on `hfx_A2` reproduces `hds1_E` exactly).

| Group | develop + fix | #10378 priors | #10379 default |
| --- | --- | --- | --- |
| HF-X HCD | 12571.0 PSMs, 85.0 entrapment, **1.35%** | 12574.0, 81.7, **1.30%** | 12577.3, 88.3, **1.40%** |
| Astral HCD | 6463.3, 55.7, **1.72%** | 6731.0, 68.7, **2.04%** | 6734.0, 57.7, **1.71%** |
| Exploris 480 TMTpro | 7346.7, 39.0, **1.06%** | 7457.3, 51.0, **1.37%** | 7415.7, 40.0, **1.08%** |

- **#10378 on Astral and TMTpro:** the extra PSMs bring proportionally many entrapment hits. On Astral, +268
  PSMs bring +13 entrapment hits; on TMTpro, +111 PSMs bring +12. The entrapment count rises in 8 of 9
  file × seed pairs on Astral and 6 of 9 on TMTpro (one falls in each). Seeds are technical repeats of the
  same data, so these counts overstate the evidence; the direction is nevertheless consistent.
- **#10379's default** keeps the FDP on all three groups: on Astral, +271 PSMs bring +2 entrapment hits.
- **Earlier Velos check** (natural UPS1/yeast entrapment, 2026-10-01): #10378 kept 0.91% against 0.95%.
  That check covers only low-resolution CID data.
- **Untested explanation:** #10378's noise model scores reversed sequences, and the search decoys are
  protein-reversed. Features that separate reversed sequences from targets can reject decoys better than
  other false matches, which would make target-decoy FDR estimates optimistic.

## Files

- `run.py` (arms `d_*`), `chain_d.sh`, `chain_dE.sh`: the searches. They run inside the harness of
  [pr-untangled-20261001](../pr-untangled-20261001) (`evaluate.py`, `build.py`, `harness/`).
- `evaluate_d.py`, `checks_d.py`, `evaluate_entrap_d.py`, `entrap_seeds_d.py`, `cmp_refs.py`: evaluation.
- `open-prs-results.zip`: per job `summary.json`, parameters, resolved settings, logs and per-seed Percolator
  results. `accepted_entrapment_psms.tsv.gz`: accepted PSMs (q ≤ 0.01) of the entrapment arms and their
  baselines, with an entrapment flag.
