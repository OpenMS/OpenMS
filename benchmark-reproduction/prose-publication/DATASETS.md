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
  - [PXD028735](https://www.ebi.ac.uk/pride/archive/projects/PXD028735) (Van Puyvelde et al., Sci Data 2022): source of ProteoBench Module 2, HF-X, HYE mix, "Alpha" series.
  - [PXD070049](https://www.ebi.ac.uk/pride/archive/projects/PXD070049) ("LFQ Benchmark Dataset – Generation Beta", public since 2026-02-03): source of the ProteoBench Astral DDA module (Module 8 here). Mix A = H65/Y30/E5, B = H65/Y15/E20.
  - [ProteoBench](https://proteobench.readthedocs.io) documents both modules.
- **ProteoBench:** the HF-X and Astral files are the ProteoBench copies on archive.openms.de. Each of the six raw files there has the same byte size as the file of the same name on PRIDE (PXD028735, PXD070049); the hashes were not compared. For the paper, cite the PRIDE accessions.
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
  - **timsTOF:** added later with the same `human` database (#10639).
  - **Lumos:** only a whole-protein shuffle so far (#10639).
  - For the paper the entrapment databases are rebuilt; see section 3.
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

### 2.2 Orbitrap Astral

**Recommendation: the ProteoBench Astral DDA files from [PXD070049](https://www.ebi.ac.uk/pride/archive/projects/PXD070049) as the Astral dataset of the paper.** They are already in section 1 (`astral_A2`, `astral_B1`, `astral_B3`). Reasons:

- **Ground truth:** the HYE mix gives a known species composition for entrapment and known fold changes for quantification.
- **Same protocol as the rest of the suite:** the same species-split FASTA and the same entrapment design.
- **Comparable to other tools:** the public ProteoBench submissions of other engines use these files, and ANDES benchmarks `Condition_A_REP1` of the same set.
- **Results so far (all on these files):**
  - #10364 (8,000 MS2 per file): ProSE +20.5% PSMs over Sage, −13.6% against ANDES, which uses an Astral-specific model.
  - #10638 (chimeric search, whole `Condition_A_REP1` run, 123,539 MS2): PSMs at q ≤ 1% 42,469 → 91,875, peptides 26,349 → 31,703. The rank-2 FDP alone is 2.23%, about twice nominal.

**Optional extensions from the same project (not used so far):**
- The other replicates (REP4–6), Condition C (H65/Y3/E32), and DDA at 5 min and at 250 pg, for robustness across gradient length and load.
- The single-species runs (Human/Yeast/Ecoli_01–03) as natural entrapment, e.g. human-only runs searched against HYE: every yeast or E. coli hit is false.

**Secondary: [PXD046453](https://www.ebi.ac.uk/pride/archive/projects/PXD046453)** (Guzman et al., Nat Biotechnol 2024).
- An independent lab: 200 ng HeLa, about 14 min, 12 DDA raw files (`20230108_AST_Neo1_DDA_UHG_HeLa_200ng_2th2p5ms_{Cycle,top50,top75,top100}_0[1-3]*.raw`, 1.6–1.8 GB each).
- It has no ground truth, so it needs a shuffled entrapment.
- #10261 measured one of its runs (which one is not stated): 26,916 PSMs natively, 35,057 with Percolator, on develop of 2026-09.

### 2.3 Low resolution (ion-trap MS2)

**iPRG2015 and iPRG2016 are not low-res.** Both are Q Exactive HCD with Orbitrap MS2:
- **iPRG2015** ([PXD010981](https://www.ebi.ac.uk/pride/archive/projects/PXD010981), MassIVE MSV000079843): MS2 at 17,500 resolution; 200 ng yeast with spiked proteins.
- **iPRG2016** (PrEST standard, [PXD008425](https://www.ebi.ac.uk/pride/archive/projects/PXD008425)): MS2 at 30,000 resolution; PrEST pools in an E. coli background. That "B1" in #10261 is `mixtureBrep1.raw` is inferred.

Both would fit a high-resolution block. External references exist: MS-GF+ on iPRG2015 and Comet on iPRG2016 in the OpenMS archive.

**Genuinely low-resolution candidates:**

| Dataset | Instrument / fragmentation | Sample, ground truth | Status | ProSE so far (develop 2026-09, native / Percolator, full runs) |
|---|---|---|---|---|
| [PXD001819](https://www.ebi.ac.uk/pride/archive/projects/PXD001819) | LTQ Orbitrap Velos, ion-trap CID | UPS1 at 9 levels in yeast; UPS1 against human proteins as natural entrapment | in section 1 (3 runs); 24 more runs in the project | UPS1 2500 amol R1: 10,827 / 13,739; 50000 amol R1: 11,673 / 14,512 (#10261) |
| [PXD007683](https://www.ebi.ac.uk/pride/archive/projects/PXD007683) TMT arm | Fusion Lumos, ion-trap CID, TMT6plex | human/yeast | in section 1 (3 runs) | see #10364 |
| Percolator paper sets (Käll et al. 2007), [data](https://noble.gs.washington.edu/proj/percolator/data/) | LTQ, ion-trap CID, precursor charge unknown (±3 Da) | yeast-01…04 (about 35k spectra each), worm-01/02, yeast chymotrypsin; no ground truth; MS2 format only, no raw files | used in #10261 | yeast-01: 7,388 / 9,954; worm-01: 8,604 / 9,903; chymotrypsin: 1,254 / 2,036 |
| [PXD042233](https://www.ebi.ac.uk/pride/archive/projects/PXD042233) (Webel et al. 2024, HeLa QC) | Fusion Lumos, ion-trap CID (and many other instruments) | HeLa; no ground truth | one run in #10261; which one is not stated | 21,501 / 28,173 |
| [PXD011695](https://www.ebi.ac.uk/pride/archive/projects/PXD011695) HCD | Orbitrap Elite, HCD with ion-trap MS2 | mouse brain; authors' PSMs as reference | in section 2.1 | 56,375 / 65,519 |
| ProteomeTools ([PXD004732](https://www.ebi.ac.uk/pride/archive/projects/PXD004732), [PXD010595](https://www.ebi.ac.uk/pride/archive/projects/PXD010595)) | Fusion Lumos, CID/HCD/ETD with ion-trap and Orbitrap readout | synthetic human peptides, so the sequence is ground truth | not used so far | — |

**Suggested low-res block:**
- **Core:** PXD001819 and the Lumos CID TMT runs, both already in section 1.
- **Historical comparability:** the Percolator yeast and worm sets, which Percolator, Crux/Tide and Comet were benchmarked on.
- **Optional:** ProteomeTools ion-trap CID pools. The same peptides were also measured with Orbitrap readout, which gives a resolution-matched comparison.

### 2.4 Immunopeptidomics (HLA class I and II, unspecific cleavage)

The two datasets come from the same donors (IDs `DN02`…`DN17`; that they are the same people is inferred from the IDs, not confirmed)
on two instruments, so they can be compared directly.

**[PXD058436](https://www.ebi.ac.uk/pride/archive/projects/PXD058436): MHCquant2 benign tissues** (Scheid et al., Genome Biol 2025)
- **Data:** timsTOF Pro 2, DDA-PASEF. 633 runs (312 class I, 321 class II), 12 donors, about 30 tissues, usually 3 technical replicates.
- **Reference results on PRIDE:** the deposit contains result tables of MHCquant2, FragPipe and PEAKS.
- **OpenMS copy:** a 30-sample subset lies under `https://archive.openms.de/openms/testfiles/id_bruker/<classI|classII>/<low|middle|high>/`. Its mzML files match the PRIDE mzML in byte size.

**ProSE on PXD058436 so far** (no comparison with another engine):

| PR | Data | Result |
|---|---|---|
| #10630, #10636 | 8,000-MS2 subsets of DN02_Spleen, DN13_Corpus and DN12_Liver class I | A/B comparisons only |
| #10639 | DN08_Brain_classI_techRep1, full run (67,944 MS2) | 2,011 PSMs at 1% FDR |
| #10639 | DN09_Pancreas_classII_techRep1, full run (120,442 MS2) | 8,399 PSMs at 1% FDR |
| #9203 | DN17_Liver_classI_techRep2 | speed only |

- **Memory:** an unspecific human index does not fit into 12.5 GiB. The search needs `low_memory` or `database:mass_slices`, which is in PR #10639, still open.
- **FDR control:** entrapment hits on MHC data were too few to check it (2–7 per run in #10630/#10636). The paper uses no entrapment for immunopeptidomics (section 3).

**[PXD019643](https://www.ebi.ac.uk/pride/archive/projects/PXD019643): HLA Ligand Atlas** (Marcu et al., J Immunother Cancer 2021)
- **Data:** Orbitrap Fusion Lumos, CID. 1,470 raw files: 751 class I (W6/32, 400–650 m/z) and 718 class II (Tü39 + L243, 400–1000 m/z). Usually 3 DDA replicates per donor, tissue and class.
- **Reference results on PRIDE:** MHCquant v1 (Comet + Percolator).

**ProSE on PXD019643 so far:** only in PR #9975 (PeptDeep rescoring), on six class I files: prostate, liver, muscle, pancreas, bladder, aorta.
- 97,656 MS2 in total. PSMs at 1% FDR after Percolator: 9,764 with the 14 default features, 10,893 with the five PeptDeep features (+11.6%).
- The exact file names are not recorded in the PR.

**Suggested subset:** three donors, both classes, high- and middle-yield samples, the same tissues on both instruments.

PXD058436 (techRep1 of each; mzML from archive.openms.de, `.d.tar.gz` from PRIDE):

| Run | Class, yield | mzML | Bruker .d |
|---|---|---|---|
| DN12_Liver_classI_techRep1 | I, high | [mzML](https://archive.openms.de/openms/testfiles/id_bruker/classI/high/DN12_Liver_classI_techRep1.mzML) (2.2 GB) | [d.tar.gz](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058436/DN12_Liver_classI_techRep1.d.tar.gz) (3.1 GB) |
| DN02_Spleen_classI_techRep1 | I, high | [mzML](https://archive.openms.de/openms/testfiles/id_bruker/classI/high/DN02_Spleen_classI_techRep1.mzML) (2.0 GB) | [d.tar.gz](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058436/DN02_Spleen_classI_techRep1.d.tar.gz) (2.7 GB) |
| DN13_Corpus_classI_techRep1 | I, high | [mzML](https://archive.openms.de/openms/testfiles/id_bruker/classI/high/DN13_Corpus_classI_techRep1.mzML) (2.1 GB) | [d.tar.gz](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058436/DN13_Corpus_classI_techRep1.d.tar.gz) (2.6 GB) |
| DN09_Pancreas_classI_techRep1 | I, middle | [mzML](https://archive.openms.de/openms/testfiles/id_bruker/classI/middle/DN09_Pancreas_classI_techRep1.mzML) (1.6 GB) | [d.tar.gz](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058436/DN09_Pancreas_classI_techRep1.d.tar.gz) (1.8 GB) |
| DN08_Brain_classI_techRep1 | I, middle | [mzML](https://archive.openms.de/openms/testfiles/id_bruker/classI/middle/DN08_Brain_classI_techRep1.mzML) (1.3 GB) | [d.tar.gz](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058436/DN08_Brain_classI_techRep1.d.tar.gz) (1.4 GB) |
| DN12_Liver_classII_techRep1 | II, middle | [mzML](https://archive.openms.de/openms/testfiles/id_bruker/classII/middle/DN12_Liver_classII_techRep1.mzML) (2.1 GB) | [d.tar.gz](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058436/DN12_Liver_classII_techRep1.d.tar.gz) (2.7 GB) |
| DN09_Pancreas_classII_techRep1 | II, middle | [mzML](https://archive.openms.de/openms/testfiles/id_bruker/classII/middle/DN09_Pancreas_classII_techRep1.mzML) (2.0 GB) | [d.tar.gz](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058436/DN09_Pancreas_classII_techRep1.d.tar.gz) (2.4 GB) |
| DN08_Spleen_classII_techRep1 | II, high | [mzML](https://archive.openms.de/openms/testfiles/id_bruker/classII/high/DN08_Spleen_classII_techRep1.mzML) (1.9 GB) | [d.tar.gz](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058436/DN08_Spleen_classII_techRep1.d.tar.gz) (2.3 GB) |
| DN13_Corpus_classII_techRep1 | II, high | [mzML](https://archive.openms.de/openms/testfiles/id_bruker/classII/high/DN13_Corpus_classII_techRep1.mzML) (2.1 GB) | [d.tar.gz](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058436/DN13_Corpus_classII_techRep1.d.tar.gz) (2.6 GB) |

PXD019643 (DDA_1 of each; the DDA_2 and DDA_3 replicates lie next to them on PRIDE):

| Run | Class | mzML | raw |
|---|---|---|---|
| 170413_AM_AUT01-DN12_Liver_W6-32_10_DDA_1_400-650mz_msms17_standard | I | [mzML](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170413_AM_AUT01-DN12_Liver_W6-32_10_DDA_1_400-650mz_msms17_standard.mzML) (248 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170413_AM_AUT01-DN12_Liver_W6-32_10_DDA_1_400-650mz_msms17_standard.raw) (478 MB) |
| 170413_AM_AUT01-DN12_Liver_Tue39L243_10_DDA_1_400-1000mz_msms24_standard | II | [mzML](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170413_AM_AUT01-DN12_Liver_Tue39L243_10_DDA_1_400-1000mz_msms24_standard.mzML) (436 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170413_AM_AUT01-DN12_Liver_Tue39L243_10_DDA_1_400-1000mz_msms24_standard.raw) (765 MB) |
| 180921_AM_AUT01-DN09_Pancreas_W6-32_10_DDA_1_400-650mz_msms10_standard | I | [mzML](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/180921_AM_AUT01-DN09_Pancreas_W6-32_10_DDA_1_400-650mz_msms10_standard.mzML) (217 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/180921_AM_AUT01-DN09_Pancreas_W6-32_10_DDA_1_400-650mz_msms10_standard.raw) (493 MB) |
| 180921_AM_AUT01-DN09_Pancreas_Tue39L243_10_DDA_1_400-1000mz_msms17_standard | II | [mzML](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/180921_AM_AUT01-DN09_Pancreas_Tue39L243_10_DDA_1_400-1000mz_msms17_standard.mzML) (296 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/180921_AM_AUT01-DN09_Pancreas_Tue39L243_10_DDA_1_400-1000mz_msms17_standard.raw) (592 MB) |
| 160325_DK_AUT01-DN08_Spleen_W6-32_20_DDA_1_400-650mz_msms1_SpeedVac | I | [mzML](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/160325_DK_AUT01-DN08_Spleen_W6-32_20_DDA_1_400-650mz_msms1_SpeedVac.mzML) (320 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/160325_DK_AUT01-DN08_Spleen_W6-32_20_DDA_1_400-650mz_msms1_SpeedVac.raw) (603 MB) |
| 160325_DK_AUT01-DN08_Spleen_Tue39L243_20_DDA_1_400-1000mz_msms11_SpeedVac | II | [mzML](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/160325_DK_AUT01-DN08_Spleen_Tue39L243_20_DDA_1_400-1000mz_msms11_SpeedVac.mzML) (479 MB) | [raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/160325_DK_AUT01-DN08_Spleen_Tue39L243_20_DDA_1_400-1000mz_msms11_SpeedVac.raw) (826 MB) |

**For a fair comparison with the published references:**
- **Reference numbers:** MHCquant2, FragPipe and PEAKS (PXD058436) and MHCquant v1 (PXD019643) report unique peptides at 1% peptide-level FDR, with all technical replicates of a sample pooled. For example, MHCquant2 reports 18,472 peptides for DN12_Liver_classI; ProSE has so far counted PSMs on single replicates.
- **What ProSE would need:** all replicates of a sample, peptide-level FDR, and the same database as MHCquant2 (Swiss-Prot human; 20 ppm; z 1–4; 8–14 aa for class I).
- **Possible ground truth:** PXD058436 also has six synthetic peptide mix runs (`*_Synpepmix*_classI_*.raw`, Thermo), which have not been examined yet.

## 3. Databases and entrapment

Status 2026-10-23: **every dataset has a FASTA. Entrapment databases exist only for the main benchmark**, and only as scripts plus SHA256, not as files.

### 3.1 Tools

Since 2026-10-23, PR [#10629](https://github.com/OpenMS/OpenMS/pull/10629) is merged into develop. It adds two tools:

- **`DecoyDatabase -entrapment_ratio N`:** writes N shuffled entrapment copies of each target protein (prefix `ENTRAP_`), each with its own decoys.
  - Each copy is shuffled peptide by peptide. The C-terminal residue and K, R and P stay in place, and residues are permuted, never replaced, so entrapment peptides are mass-matched.
  - Shuffles that equal a target peptide or exceed the identity threshold are rejected.
  - The output is deterministic for a given `-seed`.
- **`EntrapmentEvaluator`:** estimates the FDP at each q-value threshold. It uses the combined estimator `N_E / N · (1 + T/E)`, with T/E measured from the database (`-database`, I = L).

**How the design was validated (#10629):** on HF-X, Astral and TMTpro it was compared with a foreign-proteome entrapment (reviewed *Arabidopsis*, UP000006548).
- The estimates were within 0.15–0.45 pp of the foreign reference.
- They were close to the #10364 `make_entrapment.py` design.
- The ANDES whole-protein shuffle underestimated the FDP by 20–33%. ANDES has since switched to the same design (bigbio/andes#237).

**Plan for the paper:** rebuild all shuffled entrapment databases with `DecoyDatabase -entrapment_ratio 1` from develop, instead of the per-PR scripts.
- **Several realisations:** build **3 realisations** (`-seed 1/2/3`) per database. A single realisation scatters by 0.1–0.27 pp (SD) on 3 files × 8,000 spectra.
- **Reference:** add the *Arabidopsis* foreign entrapment as a reference for 2–3 groups.
- **Archive:** freeze every FASTA (with UniProt release and download date) and every entrapment database on archive.openms.de, with SHA256 in this file.

### 3.2 Main benchmark

| Group | Target database | Entrapment so far | For the paper |
|---|---|---|---|
| Velos | `hy` | natural: UPS1/yeast against all human proteins except UPS1, r = 3.376 (`make_velos_natural.py` adds `ENTRAP_` to those accessions) | keep |
| HF-X, Astral | `hye` | shuffled r = 1, `make_entrapment.py` (#10364 supplement `sage-comparison-20261002`), SHA256 `fb77e3df…` | DecoyDatabase, 3 seeds |
| TMTpro, timsTOF | `human` | same design, SHA256 `76d25071…` | DecoyDatabase, 3 seeds |
| Lumos LFQ, Lumos TMT | `hy` | only a whole-protein 1:1 shuffle (#10639, SHA256 `1462fb72…`), the design that underestimates | DecoyDatabase, 3 seeds |

- **Target databases:** `hy`, `hye` and `human` are regenerated from the HYE FASTA (`d9ac434d…`) by `prepare_databases.py` (reproduction package), with SHA256 in its `databases.json`.
- **Decoys, open decision:** the frozen suite uses full-protein reversal (`DECOY_`). DecoyDatabase reverses peptide by peptide. In #10629, decoy co-mapping changed the estimate by 0.05–0.1 pp, so one of the two decoy designs has to be chosen for all arms.
- **Scripts:** the current scripts are on branch `claude/andes-benchmarks` (`benchmark-reproduction/andes-ports/pr10629/scripts/`: `prepare_databases_nopyopenms.py`, `frozen_make_entrapment.py`, `make_velos_natural.py`, `make_foreign.py`).

### 3.3 Further datasets

| Dataset | FASTA used so far | FASTA from the authors (PRIDE) | Entrapment for the paper |
|---|---|---|---|
| PXD011695 (ETD, mouse) | Swiss-Prot mouse, 17,277 entries; file and release not recorded | [010090_MusMusculus_Uniprot_ISOFORMS_CONCAT_20161109.fasta](https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/05/PXD011695/010090_MusMusculus_Uniprot_ISOFORMS_CONCAT_20161109.fasta) (95 MB, with isoforms and concatenated decoys) | DecoyDatabase shuffle |
| PXD065289 (ETciD, HEK293) | Swiss-Prot human, 20,416 entries; not recorded | [fasta.zip](https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/03/PXD065289/fasta.zip) (34 MB) | DecoyDatabase shuffle |
| PXD070049 (Astral HYE) | ProteoBench HYE FASTA (section 1) | [uniprotkb_proteome_HYE_UniversalContaminants.fasta](https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/02/PXD070049/uniprotkb_proteome_HYE_UniversalContaminants.fasta) (20 MB) | as HF-X/Astral; optional natural entrapment from the single-species runs |
| PXD046453 (Astral HeLa) | Swiss-Prot human, 20,416 entries; not recorded | none deposited | DecoyDatabase shuffle |
| Percolator LTQ sets | [yeast.fasta.gz](https://noble.gs.washington.edu/proj/percolator/data/yeast.fasta.gz), [worm.fasta.gz](https://noble.gs.washington.edu/proj/percolator/data/worm.fasta.gz) | same page (also `decoy0–9.fasta.gz`) | DecoyDatabase shuffle |
| iPRG2015 (PXD010981) | archive reference searches | [iPRG2015_decoy.fasta](https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/09/PXD010981/iPRG2015_decoy.fasta) (8 MB) | DecoyDatabase shuffle |
| iPRG2016 (PXD008425) | archive reference searches | [prest_pool_a.fasta](https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/12/PXD008425/prest_pool_a.fasta), [prest_pool_b.fasta](https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/12/PXD008425/prest_pool_b.fasta), [prest_1000_random.fasta](https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/12/PXD008425/prest_1000_random.fasta) | possibly natural (?): if the 1,000 random PrESTs are absent from the samples, as the file names suggest; to be checked against the study design |
| PXD042233 (Lumos CID HeLa) | Swiss-Prot human, 20,416 entries; not recorded | [UP000005640_9606.fasta](https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/12/PXD042233/UP000005640_9606.fasta), [contaminants.fasta](https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/12/PXD042233/contaminants.fasta) | DecoyDatabase shuffle |
| ProteomeTools | — | peptide lists of the pools | not needed: synthetic peptides are ground truth |
| PXD058436 (immuno, timsTOF) | [human_sp.fasta](https://archive.openms.de/openms/testfiles/human_sp.fasta) on archive.openms.de (20,431 entries, used in #10634, #10639); MHCquant2 used Swiss-Prot of 2020-10-14 | none deposited | **none**: target-decoy FDR only |
| PXD019643 (immuno, Lumos) | — | [sp_21_04_2020_decoy.fasta](https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/sp_21_04_2020_decoy.fasta) (27 MB, with decoys) | **none**: target-decoy FDR only |

**Immunopeptidomics without entrapment:** two reasons.
- With unspecific cleavage, entrapment counts per run are too small to estimate an FDP (2–7 hits per run so far).
- DecoyDatabase's peptide-wise shuffle is designed for enzymatic peptides and has not been checked on unspecific searches.

There, FDR control rests on target-decoy only. The comparison is with the published MHCquant2, FragPipe and PEAKS results at the same FDR level.

## 4. Open decisions

- **Rescoring:** which rescoring is counted, Percolator 3.09 as an external step or ProSE's in-process Percolator (#10578).
- **timsTOF input:** the PRIDE mzML or ProSE's native `.d`/`.d.zip` reader.
- **PXD065289:** which runs.
- **PXD046453:** whether this secondary Astral set goes in, and which runs.
- **ProteomeTools:** whether the ion-trap CID pools go in.
- **Immunopeptidomics:** whether to search all technical replicates per sample (needed for comparison with MHCquant2, FragPipe and PEAKS) or one replicate.
- **Comparison engines and versions:** Sage, ANDES, MSFragger, Comet, … on the full runs.
- **Decoys:** full-protein reversal, as in #10364, or DecoyDatabase's peptide-wise reversal, for all arms.
- **FASTAs of section 2:** the authors' FASTA from PRIDE or a frozen current Swiss-Prot release.
