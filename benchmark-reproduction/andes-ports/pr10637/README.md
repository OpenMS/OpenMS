# PR #10637 benchmark: PIP-ECHO `shape_features`

Does andes' gain from peak-shape transfer features (bigbio/andes#167: +23 % precursors at transfer q 0.01) show up in
OpenMS ProteomicsLFQ / PIP-ECHO? Paired A/B on the ProteoBench DDA quant ion module (QExactive HF-X, PXD028735,
3 + 3 runs), following the #10585 recipe exactly. The only setting that changes between arms is `-PipEcho:shape_features`
(and `-PipEcho:fdr` for the operating curve). All arms of one database use the same IDs and the same `-feat_dir` checkpoints.

## Builds

| build | commit | how |
|---|---|---|
| PR | `ce6eb11b` (= feature `3f421d61` + checkpoint fix `ce6eb11b`) | system packages, Release, `-DWITH_GUI=OFF -DWITH_THERMO_RAW=OFF` |
| baseline | `74048ef5` = `git merge-base origin/develop` PR head | same tree; the 9 files the PR changes in `src/openms` and `src/topp/ProteomicsLFQ.cpp` checked out at the merge-base and rebuilt without reconfiguring (the build stamp keeps the PR's revision, so the baseline can combine the same checkpoints; detection code differs only by the shape snapshot, which is dropped at combine time without the option) |

## Inputs (sha256 in `results/inputs_sha256.txt`)

- Spectra: ProteoBench's msconvert mzML (vendor peak picking), `https://proteobench.cubimed.rub.de/raws/DDA/mzml/LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}.mzML`,
  as in #10585. (The archive.openms.de mzML are converted from the same .raw files, same SHA-1, but with ThermoRawFileParser, so not identical.)
- FASTA: ProteoBench HYE, sha256 `d9ac434d…`.
- Plain HYE database (ProteoBench arms, = #10585): `DecoyDatabase -decoy_string DECOY_ -decoy_string_position prefix -method reverse -enzyme Trypsin`, sha256 `63410dbe…`.
- Entrapment database (transfer FDP): the frozen #10364 HYE target-decoy FASTA (`hye_td.fasta`, sha256 `efbf5ac2…`, targets + `DECOY_` full reversals)
  → `frozen_make_entrapment.py` (from `../pr10629/scripts/`): one shuffled entrapment protein per target (K, R, P fixed, the rest permuted
  within each segment, so every tryptic peptide is mass-matched), r = 1, plus `DECOY_` full reversals of targets and entrapment;
  sha256 `fb77e3df…` (= the frozen DB; same DB as the #10633 benchmark). `check_entrapment.py` confirms composition and cleavage sites.

## Pipeline

The scripts carry the paths of the benchmark machine (`/home/user/data`, `/home/user/runs`, `/opt/openms-{fix,base}`); adapt them.

```bash
# identification per run (#10585): Comet 2025.01 rev 1 -> PSMFeatureExtractor -> Percolator 3.06 (THIRDPARTY) -> q <= 0.10 -> PEP
search_all.sh /opt/openms-fix/bin <db> <ids_dir>            # calls search_one.sh, all parameters there
# ProteomicsLFQ arm (#10585 plfq.sh) + ConsensusMapNormalizer -algorithm_type rt_local + MSstatsConverter
plfq.sh <prefix> <db> <ids> <feat_dir> <out> -targeted_only false -pip_echo true -detect_only          # once per db
plfq.sh <prefix> <db> <ids> <feat_dir> <out> -targeted_only false -pip_echo true -PipEcho:shape_features {false,true} -PipEcho:fdr {0.01,0.02,0.03,0.05}
#   common: -psmFDR 0.01 -proteinFDR 0.01 -picked_proteinFDR true -protein_inference aggregation -threads 4
master.sh / master2.sh                                       # the whole chain as run
score_proteobench.py <msstats.csv>                           # proteobench 0.18.6, "quantms" format, min 3 runs
transfers.py <consensusXML> 0.01 0.05                        # transfers (mbr_transfer_qvalues) and entrapment transfers
compare_cxml.sh a.consensusXML.gz b.consensusXML.gz          # differences beyond process-random ids / paths
collect.sh <runs> results/                                   # results/*.tsv
```

`score_proteobench.py` reproduces the public quantms submission `d01e87b9…` exactly (33,040 / 0.2088 / 0.3866 / 0.2236).
A transfer counts as entrapment when every protein its peptide refers to starts with `ENTRAP_`; transfer FDP = 2·N_E / N (r = 1, combined estimator).
Entrapment counts only false-donor transfers, the error andes' entrapment numbers count; a transfer of the right peptide to the wrong peak is invisible to it.

## Results

`results/proteobench.tsv` (all arms, raw and rt_local), `results/transfers.tsv`, `results/output_sha256.txt`.
Option off is identical to the baseline: MSstats CSV and rt_local CSV byte-identical at q 0.05 and 0.01; consensusXML
identical up to process-random ids, timestamps and the parameter lines (`compare_cxml.sh` → 0 lines).

Run times (4 cores, 15 GB): Comet search 9–10 min/run (HYE), 13–14 min/run (entrapment DB); PIP-ECHO detection 3 h 12 min (HYE) /
2 h 23 min (entrapment), 13.2 GB peak; one combine arm 35–70 min single core, ~2 GB; shape features do not change the combine time.

### ProteoBench (plain HYE DB, min 3 runs; raw / rt_local)

| arm | transfer q | transfers | precursors | median \|eps\| | E. coli \|eps\| | CV |
|---|---|---|---|---|---|---|
| off (= baseline) | 0.01 | 13,782 | 43,586 | 0.2205 / 0.1513 | 0.4051 / 0.3105 | 0.2327 / 0.0905 |
| **on** | 0.01 | 22,187 | **44,507** | 0.2193 / 0.1519 | 0.4009 / 0.3166 | 0.2334 / 0.0917 |
| off | 0.02 | 29,894 | 45,390 | 0.2199 / 0.1540 | 0.3986 / 0.3278 | 0.2344 / 0.0930 |
| on | 0.02 | 33,845 | 45,771 | 0.2199 / 0.1549 | 0.3960 / 0.3252 | 0.2347 / 0.0937 |
| off | 0.03 | 36,627 | 46,044 | 0.2205 / 0.1560 | 0.4009 / 0.3258 | 0.2351 / 0.0942 |
| on | 0.03 | 36,906 | 46,046 | 0.2205 / 0.1560 | 0.4009 / 0.3277 | 0.2351 / 0.0941 |
| off (= baseline) | 0.05 | 41,549 | 46,401 | 0.2214 / 0.1572 | 0.4061 / 0.3359 | 0.2358 / 0.0951 |
| on | 0.05 | 39,349 | 46,235 | 0.2211 / 0.1567 | 0.4088 / 0.3395 | 0.2354 / 0.0947 |

#10585 (PIP-ECHO, fdr 0.05): 46,375 / 0.2214 / 0.4090 / 0.2356, rt_local 0.157 / 0.338 / 0.095.

### Entrapment DB (transfer FDP)

| transfer q | arm | transfers | entrapment transfers | FDP |
|---|---|---|---|---|
| 0.01 | off | 14,485 | 21 | 0.29 % |
| 0.01 | on | 24,242 | 31 | 0.26 % |
| 0.05 | off | 43,515 | 48 | 0.22 % |
| 0.05 | on | 41,435 | 46 | 0.22 % |

Identified precursors (>= 3 runs, no transfers): 80 entrapment of 48.5k, FDP 0.33 %.

### andes (#167 / #172, 6-run public LFQ DDA benchmark, normalised, MBR, min 3 runs)

| transfer q | baseline → #167 precursors | median \|eps\| | E. coli | CV | entrapment transfers |
|---|---|---|---|---|---|
| 0.01 | 54,932 → 67,586 (+23 %) | 0.155 → 0.164 | 0.212 → 0.302 | 0.110 → 0.120 | 48 → 728 |
| 0.05 | 67,035 → 73,203 (+9 %) | 0.164 → 0.177 | 0.317 → 0.471 | 0.124 → 0.134 | 537 → 1,199 |
