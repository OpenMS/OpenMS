# OpenMS XML schemas (XSD)

This directory bundles **only the XML schemas that OpenMS loads at runtime**.

A schema here is read from disk when a file's `isValid()` method is called
(e.g. by the `XMLValidator` TOPP tool or by unit tests), via
`File::find("/SCHEMAS/<name>.xsd")`. Regular reading/writing of OpenMS files
does **not** validate against these schemas, so a schema only needs to be
shipped if something loads it — usually the *current* version of each format,
plus any older version still under a runtime or tooling contract (for example
`Param_1_7_0.xsd`, emitted for CTD/CWL export, and the `mzIdentML` versions
selected by auto-detection). Every entry in the table below is kept for one of
those reasons.

Because the whole `share/OpenMS/` tree is installed wholesale, every file in
this folder ends up in every OpenMS package. To keep that payload minimal, we
no longer bundle superseded schema versions or schemas for formats OpenMS does
not validate.

## Schemas bundled here and who loads them

| Schema(s) | Loaded by |
|---|---|
| `ConsensusXML_1_7.xsd` | `ConsensusXMLFile` |
| `FeatureXML_1_9.xsd` | `FeatureXMLFile` |
| `IdXML_1_5.xsd` | `IdXMLFile` |
| `TraML1.0.0.xsd` | `TraMLFile` |
| `TrafoXML_1_1.xsd` | `TransformationXMLFile` |
| `mzData_1_05.xsd` | `MzDataFile` |
| `mzML_1_10.xsd`, `mzML_idx_1_10.xsd` | `MzMLFile` / `ImzMLFile` |
| `mzXML_idx_3.1.xsd` → `mzXML_3.1_mod.xsd` → `separation_technique_1.0.xsd`, `general_types_1.0.xsd` | `MzXMLFile` (via `<xs:include>`) |
| `pepXML_v114.xsd` | `PepXMLFile` |
| `protXML_v6.xsd` | `ProtXMLFile` |
| `xQuest_1_0.xsd` | `XQuestResultXMLFile` |
| `Param_1_8_0.xsd` | `ParamXMLFile` |
| `Param_1_7_0.xsd` | tool_description_lib (`convertToCTD` / `convertToCWL`) |
| `ToolDescriptor_1_0.xsd` | `ToolDescriptionFile` |
| `mzIdentML1.3.0.xsd` (default) + `1.2.0` / `1.1.0` / `1.0.0` → `FuGElightv1.0.0.xsd` | `MzIdentMLFile` (version auto-detected) |

## Removed / archived schemas

Older schema versions and schemas for formats OpenMS does not validate at
runtime were removed from `develop`. They remain permanently available from the
last release that shipped them, **`release/3.5.0`**:

- Browse: <https://github.com/OpenMS/OpenMS/tree/release/3.5.0/share/OpenMS/SCHEMAS>
- Raw: `https://raw.githubusercontent.com/OpenMS/OpenMS/release/3.5.0/share/OpenMS/SCHEMAS/<file>`

Removed: superseded versions of ConsensusXML (1.0–1.6), FeatureXML (1.0–1.8),
IdXML (1.0–1.4), Param (1.0–1.4, 1.6.2), TrafoXML (1.0), mzML (1.00 + idx),
mzXML (2.1, 3.1), TraML (0.9.3), pepXML (v122); and the unused
`CTD_0_3.xsd`, `CvMapping.xsd` (a PSI schema, also at
<http://www.psidev.info/sites/default/files/CvMapping.xsd>), `qcML_0.0.7.xsd`,
`mzQCML_0_0_5.xsd`, and `mzQuantML_1_0_0-rc2.xsd`.

The historical `xsi:noNamespaceSchemaLocation` URLs embedded in older data files
(pointing at `.../develop/share/OpenMS/SCHEMAS/...`) are informational only and
were never fetched by OpenMS; use the `release/3.5.0` URLs above to retrieve
those schemas.

## Adding a new schema version

When you bump a format's schema, add the new `<name>_<version>.xsd` here, point
the corresponding `*File` class at it, and remove the superseded version **only
if nothing still references it** — check for version auto-detection (e.g.
`MzIdentMLFile`), an `<xs:include>` from another bundled schema, and tool
contracts (e.g. `Param_1_7_0.xsd`). Anything you remove stays reachable through
the release tag that last shipped it.
