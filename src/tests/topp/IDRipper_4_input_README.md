# IDRipper_4 fixture — Fix-#2a end-to-end test

## Purpose

Verify that loader-side `synthesizeRunIdentifiers` (mirroring `IdXMLFile.cpp:530`)
prevents `IDMerger` from collapsing two ripped `.idparquet` files that share a
stored `run_identifier` into a single `IdentificationRun`.

## Fixture: `IDRipper_4_input.idparquet`

A minimal merged `.idparquet` directory containing:

- **search_params.parquet** — 1 row, `run_identifier = "TopPerc_test"`
- **proteins.parquet** — 2 rows, both with `run_identifier = "TopPerc_test"`
- **protein_groups.parquet** — 0 rows
- **psms.parquet** — 4 rows in 2 PeptideIdentifications:
  - PSMs 0–1: `peptide_identification_index = 0`, `spectrum_metavalues`
    contains `file_origin = "fileA.idXML"`
  - PSMs 2–3: `peptide_identification_index = 1`, `spectrum_metavalues`
    contains `file_origin = "fileB.idXML"`
  - All 4 rows share `run_identifier = "TopPerc_test"`

The stored identifier collision is the post-Fix-#2b unreachable state — once
Fix #2b is in, code can never produce this on disk. The fixture is committed
as a binary so the test runs without depending on pre-fix code being present.

## Regeneration

Run the script `IDRipper_4_input_generate.py` in this directory. Requires only
`pyarrow`; no OpenMS build needed.

## Test flow

1. `TOPP_IDRipper_4_split` — `IDRipper -in IDRipper_4_input.idparquet
   -out tmp_dir/` produces `tmp_dir/fileA.idparquet` and
   `tmp_dir/fileB.idparquet`, each carrying stored `run_identifier =
   "TopPerc_test"`.
2. `TOPP_IDRipper_4_remerge` — `IDMerger -in tmp_dir/fileA.idparquet
   tmp_dir/fileB.idparquet -out tmp.idXML`. Each load synthesizes a fresh
   identifier (Fix #2a), so the two ProtIDs end up with distinct in-memory
   identifiers and IDMerger keeps them as separate IdentificationRuns.
3. `TOPP_IDRipper_4_assert` — `grep -c '<IdentificationRun' tmp.idXML` must
   equal 2. Without Fix #2a, both ProtIDs would have the same stored
   identifier in memory after load and IDMerger would collapse them to 1
   IdentificationRun.

The synthesized identifier suffix is non-deterministic, so we count
IdentificationRun elements rather than byte-DIFFing the idXML output.
