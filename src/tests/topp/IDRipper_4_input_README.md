# IDRipper_4 fixture

## Purpose

Verify that two ripped `.idparquet` bundles that share a run identifier are not
collapsed into a single `IdentificationRun` when `IDMerger` merges them again.

## Fixture: `IDRipper_4_input.idparquet`

A native identification bundle generated from `IDRipper_4_input.idXML`:

```
IDFileConverter -in IDRipper_4_input.idXML -out IDRipper_4_input.idparquet
```

It holds one run ("Comet", q-value) with two proteins and two identifications
of two candidates each. The identifications carry the meta value `file_origin`
("fileA.idXML" and "fileB.idXML"), so `IDRipper` splits them into two files.

## Test flow

1. `TOPP_IDRipper_4_split`: `IDRipper` writes `fileA.idparquet` and
   `fileB.idparquet`. Both are made from the same protein run and so carry
   the same run identifier.
2. `TOPP_IDRipper_4_remerge`: `IDMerger` merges the two native bundles. The
   runs keep their UUIDs and identifications; a repeated run identifier gets a
   numeric suffix, so the merged bundle has two runs.
3. `TOPP_IDRipper_4_assert`: the merged bundle, converted to idXML, must have
   two `IdentificationRun`s.

Run identifiers of legacy loads are synthesized (like `IdXMLFile` does), so the
test counts `IdentificationRun` elements rather than comparing the idXML text.
