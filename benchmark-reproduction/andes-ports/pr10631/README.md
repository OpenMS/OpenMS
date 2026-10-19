# PR #10631 — isobaric precursor purity and reporter scan windows (port of bigbio/andes#139)

Question: does the improvement andes reported for this fix (bigbio/andes@0a458a2, bigbio/andes@40f3ee1) also show up in
OpenMS?

## Arms

| arm | commit | what |
|---|---|---|
| base | `77cccda7` | merge-base of the PR with `develop` |
| pr | `59428021` | PR head as benchmarked first |
| fix | `6842d7d6` | PR head after the fix this benchmark found (only peaks inside an off-centre window count; an empty follow-up survey scan is left out of the interpolation) |

The same Release build (Ubuntu 24.04 system packages, `-DWITH_THERMO_RAW=ON`) of each commit; identical inputs, parameters
and identifications in all arms. Only the quantification differs.

## Data

| set | files | why |
|---|---|---|
| PXD017803 TKOpro16 (andes's benchmark) | `Fig2_m23145_SPSms3`, `Fig2_m23146_SPSms3` (SPS-MS3), `Fig2_m23142_hrMS2`, `Fig2_m23141_hrMS2` (hrMS2), .raw from PRIDE, converted with OpenMS `FileConverter` (native Thermo reader: FAIMS CV, scan windows, isolation target, spectrumRef) | FAIMS (3 CVs), SPS-MS3, hrMS2 with automatic first mass above 126 Th; known design: MET6 knocked out in 126–128N, HIS4 in 128C–130N, URA2 in 130C–132N |
| PXD007683 (archive.openms.de `tmt/OrbitrapFusionLumos`) | `a05058.mzML` | SPS-MS3 TMT11, ±0.2 Th MS2 windows whose target is the 2nd isotope for 11,742 of 37,330 MS2 (selected ion outside the window); yeast spiked into human |
| PXD054559 (archive.openms.de `tmt/OrbitrapEclipse`) | `xb10855.mzML` | FAIMS (3 CVs) hrMS2 TMTpro, ±0.25 Th windows, selected ion outside the window for 8,968 of 33,013 MS2 |

`results/acquisition.txt` lists these properties per file; `results/inputs.sha256` the input hashes (HYE FASTA
`d9ac434d…`, as required by the ProSE protocol). In every PXD017803 run the last MS1 before an MS2 is at another FAIMS
CV (14,377 of 14,377 MS2 in m23142).

## Pipeline (`scripts/run_all.sh`)

1. Identification once, shared by all arms: ProSE (in-process Percolator, the default since #10578), Trypsin/P, 2 missed
   cleavages, length 7–40, z 2–5, fixed Carbamidomethyl (C) and TMTpro/TMT6plex (K, N-term), ≤1 Oxidation (M),
   precursor ±20 ppm (andes's setting for this data; the ProSE protocol's ±50 ppm for TMTpro was not used), isotope
   error −2..1, fragment 0.5 Da (ion-trap MS2) or 20 ppm, 50 candidates, top 10, deduplicate, 1% PSM FDR. Databases:
   UniProt reviewed yeast (6,733, downloaded 2026-10-19) for PXD017803, ProteoBench HYE for the archive files.
   PSMs (q ≤ 0.01): m23145 4,408, m23146 4,402, m23142 5,557, m23141 5,708, a05058 11,685, xb10855 11,009 (andes on
   m23145: 4,570).
2. Quantification per arm, default parameters: `IsobaricAnalyzer -type tmt16plex` (PSMs joined by scan) and
   `IsobaricWorkflow -in <mzML> -in_id <idXML> -type tmt16plex -calculate_id_purity` (tmt11plex for a05058,
   tmt18plex for xb10855).
3. Metrics (`evaluate.py`): PSMs quantified, precursor purity (absent = 0), values changed vs base, and the andes
   knockout interference: per PSM of MET6/HIS4/URA2, mean of its knockout channels / mean of the other 12 channels,
   median over PSMs, overall and at purity ≥ 0.9 (andes `docs/benchmarks`). The two replicates of each acquisition
   are pooled. `ref_purity.py` is a Python transcription of andes's `precursor_purity` (window centred on the
   isolation target, only peaks inside it); `compare_ref.py` compares OpenMS's single-scan purities with it.

## Results

### SPS-MS3 (PXD017803 m23145 + m23146, pooled)

| | base IA | PR IA | base IW | PR IW | andes (m23145) |
|---|---:|---:|---:|---:|---:|
| PSMs quantified | 8,739 | 8,739 | 8,739 | 8,739 | 4,570 (one run) |
| purity: median / Q1 | 1.000 / 1.000 | 0.970 / 0.894 | 0.961 / 0.702 | 0.970 / 0.894 | 0 → 0.98 (median) |
| PSMs with purity ≥ 0.9 | 7,698 (88%) | 6,460 (74%) | 5,458 (62%) | 6,460 (74%) | 3,487 (76%) |
| KO interference, all PSMs | 0.081 (80) | 0.081 (80) | 0.081 (80) | 0.081 (80) | 0.127 |
| KO interference, purity ≥ 0.9 | 0.067 (75) | **0.036** (51) | 0.067 (59) | **0.036** (51) | **0.040**; reference filter 0.139 |

PR = PR head with the fix (`fix`); the `pr` arm gives the same numbers here (`results/quant_summary.tsv`). m23145
alone (andes's file): IA 0.081 (36 PSMs) → 0.019 (23).

- IsobaricAnalyzer base computed the purity of the SPS fragment in the MS1 scan, 1.0 for 88% of PSMs, so a purity
  filter removed almost nothing. The PR computes the MS2 precursor's purity in the survey scan of its spectrumRef
  (same FAIMS CV); IsobaricAnalyzer and IsobaricWorkflow now agree exactly.
- andes's "before" (purity missing for 2,308 of 4,570) has no counterpart in OpenMS: OpenMS did not leave purity
  missing but reported wrong values (SPS-fragment purity in IsobaricAnalyzer; in hrMS2 the survey scan of another CV,
  see below).

### hrMS2 (PXD017803 m23142 + m23141, pooled)

| | base IA | PR IA | base IW | PR IW | andes (m23142) |
|---|---:|---:|---:|---:|---:|
| identified PSMs whose scan window starts above 126 Th | 8,321 of 11,264 | | | | 4,345 of 6,361 |
| of these with no reporter m/z in the window (dropped in all arms) | 6,807 | | | | |
| PSMs quantified | 4,457 | 2,943 | 4,457 | 2,943 | |
| … of them missing reporter channels (zeros) | 1,514 | 0 | 1,514 | 0 | |
| KO interference, all PSMs | 0.160 (40) | 0.178 (30) | 0.160 (40) | 0.178 (30) | 0.210 |
| KO interference, purity ≥ 0.9 | 0.238 (14) | 0.152 (20) | 0.238 (14) | 0.152 (20) | |

The 1,514 partially covered scans had their low channels (126–128N, MET6's knockout channels) written as 0; that
pulled the base interference down (0.160) without any real improvement. On covered scans base and PR agree (0.178).
IsobaricAnalyzer base took the last MS1 in time, which is at another FAIMS CV for every MS2 of these runs.

### Window centre (selected ion outside the isolation window)

Single survey scan (`-extraction:purity_interpolation false`), IsobaricWorkflow, PSMs whose selected ion is outside
the window; agreement with the andes-style reference purity (`results/purity_vs_reference.tsv`):

| run | PSMs | within ±0.05 of reference: base / pr / fix | > 0.2 above reference: base / pr / fix |
|---|---:|---:|---:|
| a05058 | 2,962 | 1,661 / 2,121 / **2,935** | 238 / 42 / **2** |
| xb10855 | 2,609 | 1,902 / 2,050 / **2,446** | 216 / 170 / 112 |
| m23145 (SPS-MS3) | 181 | 126 / 155 / **181** | 7 / 4 / **0** |
| m23142 (hrMS2) | 28 | 14 / 14 / **28** | 5 / 3 / **0** |

On PSMs whose selected ion is inside the window all arms agree with the reference (a05058: 8,591 of 8,661 within
±0.05). The remaining xb10855 differences occur in all arms and on centred windows too (187 PSMs), so they are not
from this PR.

### Spike-in ratios (a05058, yeast 3 : 1.5 : 1 over 126–127C / 128N–129C / 130N–131C, human 1 : 1)

Channel groups and levels inferred from the data and the stated 3-, 2- and 1.5-fold changes
(`results/spikein_ratios_a05058.tsv`). Default quantification does not change; with a purity filter:

| purity ≥ 0.9 | yeast PSMs | 3× | 2× | 1.5× |
|---|---:|---:|---:|---:|
| no filter (all arms) | 1,009 | 3.135 | 2.149 | 1.433 |
| base IA (filter has no effect, purity 1.0) | 1,009 | 3.135 | 2.149 | 1.433 |
| base IW | 722 | 3.296 | 2.239 | 1.463 |
| PR (fix) IA = IW | 701 | 3.341 | 2.252 | 1.466 |

The filter expands all ratios; the 1.5× pair moves towards 1.5, the 3× and 2× pairs were already above their nominal
value, so this design does not separate accuracy from expansion.

## Verdict

| andes claim | OpenMS | verdict |
|---|---|---|
| SPS-MS3: purity present for all PSMs, median 0.98, 76% at ≥ 0.9 | present for all, median 0.97, 74% at ≥ 0.9 (base: 1.0 for 88% in IsobaricAnalyzer) | reproduced (OpenMS's baseline defect differs: wrong values, not missing ones) |
| SPS-MS3: KO interference 0.127 → 0.040 at purity ≥ 0.9; reference filter 0.139 | 0.081 → 0.036 at ≥ 0.9; base filter 0.067 (IsobaricAnalyzer: no effect) | reproduced in direction and size; few knockout PSMs (51–80), read as a direction |
| hrMS2: scans whose first mass is above 126 Th wrote zeros | 1,514 such PSMs quantified with zeros in base, skipped in PR; interference 0.160 (artificially low) → 0.178 | reproduced |
| window centred on the isolation target | correct after commit `6842d7d6`; before it, the out-of-window selected ion was still counted | reproduced with the fix |

## Not done

- andes itself was not re-run; its numbers come from `docs/benchmarks/README.md` on the same PXD017803 files (other
  database version, impurity sheet and extraction).
- Timing was not measured (not what the PR changes); IsobaricWorkflow took 5–12 s per TKO run in every arm.
