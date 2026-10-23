# ProSE publication: benchmark datasets

Status: draft, 2026-10-23. Every download link below was checked (HTTP 200, size as listed).

## 1. Main benchmark: 20 DDA runs, seven acquisition groups, full runs

This is the suite of issue [#10364](https://github.com/OpenMS/OpenMS/issues/10364), with one change: every run is
searched **in full** instead of the 8,000 MS2 spectra sampled per file there. All numbers published in #10364 and
the PRs that cite it are on the 8,000-spectrum subsets and have to be re-measured on the full runs.

The mzML files are the exact inputs frozen on 2026-09-30 (size and SHA256 from the reproduction package's
`inputs.json`, [README](https://github.com/OpenMS/OpenMS/blob/a5705ec8e1cb19d40ebfa7bb16f6b5651aa1dc65/benchmark-reproduction/prose-andes/README.md)).
Spectrum counts are from the package's per-file provenance (`broad_benchmark/provenance/<id>.json`).

**Totals:** 1,509,392 MS2 spectra; 1,497,401 of them with a recorded precursor charge of 2–5.

| Group | ID | Run | MS2 | Charge 2–5 | MS1 | mzML (as frozen) | Vendor raw |
|---|---|---|---:|---:|---:|---|---|
| LTQ Orbitrap Velos, ion-trap CID | `velos_125_R1` | PXD001819 UPS1_125amol_R1 | 36,802 | 36,768 | 7,277 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/LTQOrbitrapVelos/PXD001819/UPS1_125amol_R1.mzML.gz) (398 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/12/PXD001819/UPS1_125amol_R1.raw) (1.7 GB) |
| | `velos_5000_R2` | PXD001819 UPS1_5000amol_R2 | 38,046 | 38,005 | 7,195 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/LTQOrbitrapVelos/PXD001819/UPS1_5000amol_R2.mzML.gz) (405 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/12/PXD001819/UPS1_5000amol_R2.raw) (1.7 GB) |
| | `velos_25000_R3` | PXD001819 UPS1_25000amol_R3 | 40,439 | 40,398 | 6,973 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/LTQOrbitrapVelos/PXD001819/UPS1_25000amol_R3.mzML.gz) (418 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/12/PXD001819/UPS1_25000amol_R3.raw) (1.7 GB) |
| Q Exactive HF-X, HCD | `hfx_A2` | ProteoBench Module 2, Condition_A_Sample_Alpha_02 | 111,413 | 111,413 | 26,090 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/QExactiveHF/ProteoBench_Module_2/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_02.mzML.gz) (501 MB) | [raw](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/QExactiveHF/ProteoBench_Module_2/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_02.raw) (3.6 GB) |
| | `hfx_B1` | ProteoBench Module 2, Condition_B_Sample_Alpha_01 | 111,088 | 111,088 | 25,931 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/QExactiveHF/ProteoBench_Module_2/LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01.mzML.gz) (485 MB) | [raw](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/QExactiveHF/ProteoBench_Module_2/LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01.raw) (3.5 GB) |
| | `hfx_B3` | ProteoBench Module 2, Condition_B_Sample_Alpha_03 | 118,174 | 118,174 | 24,689 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/QExactiveHF/ProteoBench_Module_2/LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_03.mzML.gz) (514 MB) | [raw](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/QExactiveHF/ProteoBench_Module_2/LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_03.raw) (3.7 GB) |
| Orbitrap Astral, HCD | `astral_A2` | ProteoBench Module 8, 15min_50ng_Condition_A_REP2 | 123,621 | 123,607 | 3,680 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/OrbitrapAstral/ProteoBench_Module_8/LFQ_Astral_DDA_15min_50ng_Condition_A_REP2.mzML.gz) (1.0 GB) | [raw](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/OrbitrapAstral/ProteoBench_Module_8/LFQ_Astral_DDA_15min_50ng_Condition_A_REP2.raw) (2.6 GB) |
| | `astral_B1` | ProteoBench Module 8, 15min_50ng_Condition_B_REP1 | 122,686 | 122,659 | 3,702 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/OrbitrapAstral/ProteoBench_Module_8/LFQ_Astral_DDA_15min_50ng_Condition_B_REP1.mzML.gz) (1.0 GB) | [raw](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/OrbitrapAstral/ProteoBench_Module_8/LFQ_Astral_DDA_15min_50ng_Condition_B_REP1.raw) (2.6 GB) |
| | `astral_B3` | ProteoBench Module 8, 15min_50ng_Condition_B_REP3 | 122,647 | 122,625 | 3,697 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/OrbitrapAstral/ProteoBench_Module_8/LFQ_Astral_DDA_15min_50ng_Condition_B_REP3.mzML.gz) (1.0 GB) | [raw](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/OrbitrapAstral/ProteoBench_Module_8/LFQ_Astral_DDA_15min_50ng_Condition_B_REP3.raw) (2.6 GB) |
| Orbitrap Fusion Lumos, HCD, LFQ | `lumos_lfq_5192` | PXD007683 a05192 | 73,549 | 73,549 | 13,613 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/OrbitrapFusionLumos/PXD007683/a05192.mzML.gz) (294 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD007683/a05192.raw) (663 MB) |
| | `lumos_lfq_5194` | PXD007683 a05194 | 74,198 | 74,198 | 13,450 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/OrbitrapFusionLumos/PXD007683/a05194.mzML.gz) (295 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD007683/a05194.raw) (666 MB) |
| | `lumos_lfq_5199` | PXD007683 a05199 | 73,666 | 73,666 | 13,443 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/OrbitrapFusionLumos/PXD007683/a05199.mzML.gz) (292 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD007683/a05199.raw) (660 MB) |
| Orbitrap Fusion Lumos, ion-trap CID, TMT6plex | `lumos_tmt_5058` | PXD007683 a05058 | 37,330 | 37,330 | 6,870 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/tmt/OrbitrapFusionLumos/PXD007683/a05058.mzML.gz) (376 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD007683/a05058.raw) (537 MB) |
| | `lumos_tmt_5059` | PXD007683 a05059 | 37,536 | 37,536 | 6,471 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/tmt/OrbitrapFusionLumos/PXD007683/a05059.mzML.gz) (380 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD007683/a05059.raw) (534 MB) |
| | `lumos_tmt_5066` | PXD007683 a05066 | 41,078 | 41,078 | 6,051 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/tmt/OrbitrapFusionLumos/PXD007683/a05066.mzML.gz) (425 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD007683/a05066.raw) (593 MB) |
| Orbitrap Exploris 480, HCD, TMTpro | `eclipse_tmtpro_10855` | PXD054559 xb10855 | 33,013 | 33,013 | 8,764 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/tmt/OrbitrapEclipse/PXD054559/xb10855.mzML.gz) (128 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/10/PXD054559/xb10855.raw) (361 MB) |
| | `eclipse_tmtpro_10858` | PXD054559 xb10858 | 32,637 | 32,637 | 8,689 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/tmt/OrbitrapEclipse/PXD054559/xb10858.mzML.gz) (117 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/10/PXD054559/xb10858.raw) (333 MB) |
| | `eclipse_tmtpro_10863` | PXD054559 xb10863 | 33,517 | 33,517 | 8,667 | [mzML.gz](https://archive.openms.de/openms/benchmarks/pride-benchmarks/tmt/OrbitrapEclipse/PXD054559/xb10863.mzML.gz) (121 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/10/PXD054559/xb10863.raw) (344 MB) |
| timsTOF HT, DDA-PASEF, plasma | `tims_plasma_30min` | PXD053748 DDA_30min_200ng_QCPlasma_Slot1-36_1_4363 | 115,047 | 109,511 | 2,533 | [mzml](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD053748/DDA_30min_200ng_QCPlasma_Slot1-36_1_4363.mzml) (2.2 GB) | [d.zip](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD053748/DDA_30min_200ng_QCPlasma_Slot1-36_1_4363.d.zip) (1.6 GB) |
| | `tims_plasma_50min` | PXD053748 DDA_50min_200ng_QCPlasma_Slot1-36_1_4367 | 132,905 | 126,629 | 3,440 | [mzml](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD053748/DDA_50min_200ng_QCPlasma_Slot1-36_1_4367.mzml) (2.7 GB) | [d.zip](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD053748/DDA_50min_200ng_QCPlasma_Slot1-36_1_4367.d.zip) (2.0 GB) |

Notes on the files:

- **Project pages:**
  - [PXD001819](https://www.ebi.ac.uk/pride/archive/projects/PXD001819): UPS1 spiked into yeast.
  - [PXD007683](https://www.ebi.ac.uk/pride/archive/projects/PXD007683)
  - [PXD054559](https://www.ebi.ac.uk/pride/archive/projects/PXD054559)
  - [PXD053748](https://www.ebi.ac.uk/pride/archive/projects/PXD053748)
  - [ProteoBench](https://proteobench.readthedocs.io): Module 2 is the HF-X HYE data of PXD028735; Module 8 is the Astral DDA HYE data.
- **ProteoBench:** the HF-X and Astral files are the ProteoBench copies on archive.openms.de. The raw files next to them were not hash-checked against ProteoBench.
- **Folder names:** the archive folder `OrbitrapEclipse` holds Exploris 480 data. The suite IDs keep the `eclipse_` prefix.
- **TMT6plex runs:** each file also holds one MS3 scan per MS2 scan (37–41k). Only the CID MS2 spectra are searched.
- **timsTOF files:** the mzML is the exporter output deposited on PRIDE; it has 1+ precursors (5,536 and 6,276), which the 8,000-spectrum selection excluded. ProSE also reads the `.d.zip` natively (#10680). For the publication, decide on one of the two and state it.
- **Checksums:** the SHA256 of each frozen mzML is in `inputs.json` of the reproduction package (`prose-andes-reproduction.zip` at the README's commit).

### Database and search settings (as frozen in #10364)

- **FASTA:** [ProteoBenchFASTA_MixedSpecies_HYE.fasta](https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/QExactiveHF/ProteoBench_Module_2/ProteoBenchFASTA_MixedSpecies_HYE.fasta)
  (16,670,994 bytes, 31,889 targets, SHA256 `d9ac434d88492c10c8e9a587ee7dbc9480fa0995fa07a6ba35a7da8abf39aa25`).
  It is split by accession suffix (contaminants `Cont_` kept in all three), with one full-sequence reversed `DECOY_` per protein:

  | Database | Used for | Targets |
  |---|---|---:|
  | `hy` (human + yeast) | Velos, Lumos LFQ, Lumos TMT | 27,488 |
  | `hye` (human + yeast + E. coli) | HF-X, Astral | 31,889 |
  | `human` | TMTpro, timsTOF | 20,767 |

- **Digestion and modifications:** Trypsin/P, 2 missed cleavages, peptide length 7–40, fixed Carbamidomethyl (C), at most one Oxidation (M); TMT6plex or TMTpro fixed on the labelled groups.
- **Tolerances:** precursor 20 ppm (TMTpro 50 ppm); fragment 20 ppm, or 0.5 Da for the ion-trap CID groups (Velos, Lumos TMT).
- **Counting:** target PSMs at q ≤ 0.01 after Percolator, mean over seeds 1/42/137. #10364 used Percolator 3.09 `-Y -U`; ProSE now rescores in-process by default (#10578), so the publication should say which rescoring is counted.
- **Entrapment, the only estimate of the true error rate here:**
  - **Velos:** natural entrapment, UPS1/yeast against all human proteins, r = 3.376.
  - **HF-X, Astral, TMTpro:** one shuffled, mass-matched entrapment protein per target, r = 1.
  - **timsTOF, Lumos:** added later the same way (#10639).
- **Comparison engines measured on the suite so far:** Sage 0.14.7, ANDES (b7eaece, auto models), MSFragger (only in #10403).

### Caveat for the paper

ProSE's defaults were tuned on these 20 files (the PRs listed in #10364). ANDES was developed on the same Velos, Lumos and Astral studies.
A held-out set is therefore needed. The datasets in section 2 were not used for tuning the defaults and can serve as that set.

## 2. Further datasets (candidates)

### 2.1 Electron-based fragmentation (ETD, AI-ETD, EThcD, ETciD)

All numbers here are from full runs, measured in September 2026 (PR [#10261](https://github.com/OpenMS/OpenMS/pull/10261) comments
and PR [#10268](https://github.com/OpenMS/OpenMS/pull/10268)). Since then the scoring changed (#10391, #10397, #10399, in-process Percolator #10578), so they have to be re-measured.

**[PXD011695](https://www.ebi.ac.uk/pride/archive/projects/PXD011695):** one mouse brain digest on an Orbitrap Elite with ion-trap MS2, one run per activation method.
Searched against Swiss-Prot mouse (17,277), 0.5 Da fragments, 20 ppm precursor.
The reference is the authors' OMSSA/COMPASS search, which used a larger database, so it is context only.

| Run | Raw | PSMs at 1% FDR, ProSE + Percolator | Authors | Result |
|---|---|---:|---:|---|
| HCD | [nriley_Mouse_HCD.raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/05/PXD011695/nriley_Mouse_HCD.raw) (1.6 GB) | 65,519 | 63,686 | +3% |
| EThcD | [nriley_Mouse_EThcD.raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/05/PXD011695/nriley_Mouse_EThcD.raw) (618 MB) | 12,445 | 12,292 | +1% |
| ETciD | [nriley_Mouse_ETciD.raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/05/PXD011695/nriley_Mouse_ETciD.raw) (628 MB) | 11,373 | 10,436 | +9% |
| ETD (22,056 MS2) | [nriley_Mouse_ETD.raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/05/PXD011695/nriley_Mouse_ETD.raw) (614 MB) | 7,819 (#10268 defaults); 5,696 before (#10261) | 6,873 | +14% |
| AI-ETD | [nriley_Mouse_AI-ETD.raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/05/PXD011695/nriley_Mouse_AI-ETD.raw) (651 MB) | 10,890 | 12,921 | −16%, not re-measured with #10268 |

- **ETD:** 7,819 PSMs with the defaults merged in #10268 (b/y, plus c and z+1 by activation).
- **EThcD, ETciD, AI-ETD:** these counts are from #10261, with b/y + c/z+1 (EThcD, ETciD) or c/z+1 only (AI-ETD), before #10268 kept b/y for every activation.
- **Native FDR on ETD 2+ precursors is weak:** ProSE's top hit is the authors' peptide for 2,427 of their 3,482 2+ PSMs. The native FDR, one threshold for all charges, keeps only 159 of them; Percolator keeps 2,740.
- **The project holds more runs:** `nriley_Mouse_AI-ETDplus.raw` and yeast runs (`dbrademan_Yeast_HCD.raw`, `nriley_Yeast_AI-NETD.raw`).

**[PXD065289](https://www.ebi.ac.uk/pride/archive/projects/PXD065289):** HEK293, ETciD with Orbitrap MS2. The file names (`ASND01_…`) suggest an Orbitrap Ascend; PRIDE lists an Exploris 480.

| Run | PSMs at 1% FDR (native / Percolator) |
|---|---|
| ETciD run, z+1 ions (#10261 local switch, b/y + c/z+1) | 7,466 / 9,667 |
| Same run with z instead of z+1 | 3,901 / 6,154 |

- **Which file:** #10261 does not name the file. The project has 500+ raw files (fractionated trypsin, LysN, … runs with ETciD, UVPD and others), e.g. the 20 fractions `ASND01_NGL_YD1_241129_Expi_trypsin_ETciD_frac*.raw`.
- **No external reference:** no comparison to the authors' FragPipe/MSBooster results has been made yet. Their search results are deposited as `fragpipe_msbooster.zip`.
- **For the publication:** pick and name the runs, e.g. a few trypsin ETciD fractions, and compare against `fragpipe_msbooster.zip`.
