# Sage vs. ProSE on the 20-file suite: where ProSE loses (2026-10-02)

Supplement to OpenMS issue #10364. Sage v0.14.7 (official binary, SHA256 `d3820543…`) was rerun on the frozen
20-file suite with the matched search space of the 2026-09-30 Sage supplement (`run_sage.py`). All 20 files reproduce
that supplement's per-seed Percolator counts and normalized-PIN hashes exactly. ProSE arms use the harness, PIN export
and Percolator protocol of `pr-untangled-20261001` (`run.py`); `develop` is `f5ea2d04`.

## Per-spectrum comparison (`compare.py`)

A spectrum counts as accepted when Percolator accepts it (q ≤ 0.01) in at least two of three seeds. For spectra only
Sage accepts, ProSE's top-10 candidate list gives the reason: `no_candidate`, `absent_top10`, `below_rank1`,
`rank1_not_accepted`. Tables: `per_spectrum_base_pr10335.tsv` (develop, #10335) and `per_spectrum_ds1.tsv` (deisotoping
fix). `breakdown.py` splits them by precursor charge and peptide length, and `features.py` compares PSM features.

## Findings, in order of impact

1. **Deisotoping removes real fragment ions** (`probe/trace.py`, `probe/cluster.py`, `ionloss.py`). ProSE called
   `Deisotoper::deisotopeAndSingleCharge` with `start_intensity_check=2`, which lets an envelope's second peak exceed
   the first. A small peak 1.00335 Da below a fragment ion became the monoisotopic peak and the ion was removed.
   - **Labelled fragments:** every TMT/TMTpro-labelled fragment has such a peak, from reagent isotope impurity.
   - **Share of the identified peptide's matchable b/y peaks removed:**

     | Data | Removed now | With `start_intensity_check=1` |
     | --- | --- | --- |
     | TMTpro | 18–21% | 0.7–0.8% |
     | timsTOF | 6.4% | about 1% |
     | Astral | 5.4% | about 1% |
     | HF-X | 1.5% | about 0.5% |
     | Lumos LFQ | 0.8% | about 0.5% |
   - **Arm `ds1`:** high-resolution files 6 up, 8 flat, 0 down. TMTpro +12.3% (native TDC +25.7%), Astral +4.3%,
     timsTOF +1.0%. This closes 78% of the TMTpro gap to Sage and overtakes Sage on Astral.
   - **TMTpro entrapment** (`entrapment_tmtpro.tsv`, targets + shuffled paired-mass entrapment, r = 1): develop
     6324.7 PSMs at combined FDP 1.60%, against 7346.7 (+16.2%) at 1.06% with the fix.
   - **PR patch:** `pr-deisotoping-fix.patch` gives byte-identical PIN and native TSV to `ds1` on TMTpro and Astral,
     and to develop on Velos (arm `dfix`).
2. **Retrieval gate of 5 matched fragment ions** (`fragment:min_matched_ions`; Sage uses 4). Most spectra ProSE
   gives no candidate (HF-X 550, Lumos LFQ 391, Astral 273) are weak spectra that Sage accepts with 4–7 matched peaks.
   - **Probe on 128 HF-X spectra:** at 5, none has a candidate. At 4, 103 do, 85 with Sage's peptide at rank 1.
   - **`ds1_mm4` vs develop:** 12 up, 2 flat, 0 down of 14. HF-X +3.2%, Astral +5.6%, Lumos LFQ +2.2%,
     TMTpro +13.2%, timsTOF +1.9%. ProSE then matches or exceeds Sage on HF-X, Astral and Lumos LFQ.
   - **`ds1_mm3`:** weaker. An entrapment check of the gate is pending.
3. **Only singly charged fragments are scored.** Retrieval uses fragment charges up to 2, but HyperScore scores
   charge 1 only.
   - **`fz2`** (scoring up to charge 2, capped at precursor charge − 1): Velos +7.3% and TMTpro +3.3%, but Astral
     −3.7%, and native TDC drops 10–37% in every group.
   - **Low resolution:** #10335's calibrated multi-charge scoring is the better route (Velos +11.7%).
   - Not proposed as is.

`prototype-switches.patch` holds the benchmark switches: `scoring:max_fragment_charge` and
`fragment:deisotope_start_check`. `group_summary.txt` and `arms_vs_develop.txt` give group means and per-file z-scores.
`sage-comparison-results.zip` holds `summary.json`, parameters, logs and per-seed results for every ProSE arm and
Sage file.
