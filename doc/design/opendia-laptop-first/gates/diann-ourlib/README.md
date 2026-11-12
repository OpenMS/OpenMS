# DIA-NN on our exact library (diagnostic)

Question: how many precursors does DIA-NN itself accept on raw01 with OUR library (same 300,223 target precursors,
6 fragments, library RT and intensities), single run, no MBR, no second pass? This splits the gap to the ProteoBench
DIA-NN reference (8.65M predicted library, 3 runs, MBR/reanalyse) into a library/setup part and an algorithm part.

- DIA-NN 2.2.0 (local docker image `diann:2.2.0`; ProteoBench used 2.3.1).
- Input: `dd-plots-20261111/raw01.mzML` (FileConverter of the same .raw).
- Library: target rows of the exact library TSV (`lib_targets.tsv`); DIA-NN generates its own decoys.
- Settings mirror ProteoBench where applicable: auto mass accuracy, --unimod4, charges 1-4 as in the library,
  --qvalue 1 so we can filter ourselves; no --reanalyse, no --rt-profiling, no MBR (single file).
- Evaluation: DIA-NN's own Q.Value <= 0.01 precursors, entrapment FDP_c with the same r as our evaluator,
  overlap with the ProteoBench reference and with D11FWT accepted IDs.
