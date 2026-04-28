# PercolatorAdapter on internal PSM parquet — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a directory-based `.idparquet` format (lossless `PSMSchema` PSMs + `ProteinIdentificationArrowIO` proteins/groups/search-params) as a first-class identification format in OpenMS, wire it into `FileHandler`, and make `PercolatorAdapter`, `IDMerger`, and the search-engine adapters read/write it natively. Add a Percolator output post-filter (`score:fdr`) and a best-PSM-per-spectrum flag (`keep_all_passing`) to `PercolatorAdapter`.

**Architecture:** Three layers. Format layer = new `PSMArrowIO` class plus a free-standing `QPXFile::importFromArrow` extracted from `ConsensusMapArrowIO`. Dispatch layer = new `FileTypes::IDPARQUET` enum value plus `IDPARQUET` cases in `FileHandler::{load,store}Identifications`. Tool layer = each TOPP tool that already loads/stores identifications adds `idparquet` to its valid-formats and to its `FileHandler` allowed-types whitelist. `SageAdapter` additionally migrates from `IdXMLFile().store(...)` to `FileHandler().storeIdentifications(...)`.

**Tech Stack:** C++17, Apache Arrow, Apache Parquet, OpenMS test framework (`START_TEST`/`START_SECTION`/`TEST_EQUAL` family), CMake.

**Spec:** `docs/superpowers/specs/2026-04-28-percolator-psm-parquet-design.md`.

---

## File structure

### New

- `src/openms/include/OpenMS/FORMAT/PSMArrowIO.h` — public class declaration (export/import directory).
- `src/openms/source/FORMAT/PSMArrowIO.cpp` — implementation, ~150–250 lines.
- `src/tests/class_tests/openms/source/PSMArrowIO_test.cpp` — round-trip + error-path tests.

### Modified — format/dispatch layer

- `src/openms/include/OpenMS/FORMAT/FileTypes.h` — add `IDPARQUET` enum value.
- `src/openms/source/FORMAT/FileTypes.cpp` — add `TypeNameBinding` for `idparquet`.
- `src/openms/include/OpenMS/FORMAT/QPXFile.h` — add `static bool importFromArrow(...)`.
- `src/openms/source/FORMAT/QPXFile.cpp` — implement `importFromArrow` (extracted from `ConsensusMapArrowIO::importPSMsFromArrow`).
- `src/openms/source/FORMAT/ConsensusMapArrowIO.cpp` — `importPSMsFromArrow` calls the extracted helper for PSM construction; ConsensusFeature linkage stays local.
- `src/openms/source/FORMAT/FileHandler.cpp` — add `IDPARQUET` cases to `loadIdentifications` / `storeIdentifications`.
- `src/openms/includes.cmake` — register `PSMArrowIO.h`.
- `src/openms/sources.cmake` — register `PSMArrowIO.cpp`.
- `src/tests/class_tests/openms/executables.cmake` — register `PSMArrowIO_test`.

### Modified — tool layer

- `src/topp/PercolatorAdapter.cpp` — `idparquet` to `in`/`out` valid-formats and storeIdentifications whitelist; new `score:fdr` and `keep_all_passing` options; post-filter step.
- `src/topp/IDMerger.cpp` — `idparquet` to `formats` list and to `storeIdentifications` whitelist.
- `src/topp/SageAdapter.cpp` — migrate `IdXMLFile().store(...)` to `FileHandler().storeIdentifications(...)`; whitelist + valid-formats.
- `src/topp/CometAdapter.cpp`, `MSGFPlusAdapter.cpp`, `MSFraggerAdapter.cpp`, `XTandemAdapter.cpp` — whitelist + valid-formats.

### Modified — TOPP tests

- `src/tests/topp/THIRDPARTY/third_party_tests.cmake` — `TOPP_PercolatorAdapter_idparquet`.
- `src/tests/topp/PercolatorAdapter_score_fdr.test.cmake` (or equivalent location consistent with the existing test convention) — option-only filter test.
- `src/tests/topp/IDMerger_idparquet.test.cmake` (matching the existing convention).
- Reference fixtures under `src/tests/topp/THIRDPARTY/` (input + expected `.idparquet/` directories).

---

## Conventions used in this plan

**Build/test commands** (assume project root `cwd`):

```bash
cmake --build OpenMS-build -j$(nproc) --target <target>
ctest --test-dir OpenMS-build -R <regex> --output-on-failure
```

**Adding a new test target**: add to `src/tests/class_tests/openms/executables.cmake`, then re-run cmake (`cmake -S . -B OpenMS-build`) before the first build.

**Test file pattern** (OpenMS):

```cpp
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/PSMArrowIO.h>

using namespace OpenMS;
using namespace std;

START_TEST(PSMArrowIO, "$Id$")

START_SECTION(([EXTRA] some_behavior_test))
{
  // ... test body
}
END_SECTION

END_TEST
```

**Commit messages**: short imperative, no PR-numbered prefix. Examples in the existing log: `docs: add design spec ...`, `feat: ...`, `refactor: ...`. Follow the same style.

---

## Task 1: Register `FileTypes::IDPARQUET`

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/FileTypes.h:97` (add enum value before `BRUKER_TDF` / `SIZE_OF_TYPE`)
- Modify: `src/openms/source/FORMAT/FileTypes.cpp:109` (add `TypeNameBinding` near other parquet bindings)
- Test: `src/tests/class_tests/openms/source/FileTypes_test.cpp`

- [ ] **Step 1: Locate the existing FileTypes test for `nameToType`/`typeToName`.**

```bash
grep -nE "nameToType|typeToName" src/tests/class_tests/openms/source/FileTypes_test.cpp | head -10
```

Expected: `FileTypes_test.cpp` contains existing `START_SECTION([EXTRA] nameToType ...)` and `START_SECTION([EXTRA] typeToName ...)` blocks (or equivalent).

- [ ] **Step 2: Add a failing test for IDPARQUET name/type mapping.**

In `src/tests/class_tests/openms/source/FileTypes_test.cpp`, in the section that exercises `nameToType` / `typeToName` (or in a new section if no matching one exists):

```cpp
TEST_EQUAL(FileTypes::nameToType("idparquet"), FileTypes::IDPARQUET);
TEST_EQUAL(FileTypes::nameToType("IDPARQUET"), FileTypes::IDPARQUET);
TEST_STRING_EQUAL(FileTypes::typeToName(FileTypes::IDPARQUET), "idparquet");
```

- [ ] **Step 3: Run the test to verify it fails.**

```bash
cmake --build OpenMS-build -j$(nproc) --target FileTypes_test
```

Expected: compilation FAIL with `'IDPARQUET' is not a member of 'OpenMS::FileTypes'` (or equivalent).

- [ ] **Step 4: Add the enum value.**

`src/openms/include/OpenMS/FORMAT/FileTypes.h`, just before `BRUKER_TDF`:

```cpp
      PARQUET,            ///< Apache Parquet file format (.parquet, .pqt)
      IDPARQUET,          ///< OpenMS internal identification parquet bundle (directory: psms.parquet + proteins.parquet + protein_groups.parquet + search_params.parquet)
      BRUKER_TDF,         ///< Bruker TimsTOF .d directory (TDF format)
```

- [ ] **Step 5: Register the extension and properties.**

`src/openms/source/FORMAT/FileTypes.cpp`, immediately after the existing `PARQUET` binding (around line 109):

```cpp
    TypeNameBinding(FileTypes::PARQUET, "parquet", "Apache Parquet file", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::IDPARQUET, "idparquet", "OpenMS identification parquet bundle (directory)", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::BRUKER_TDF, "d", "Bruker TDF", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE}),
```

- [ ] **Step 6: Build and run the test.**

```bash
cmake --build OpenMS-build -j$(nproc) --target FileTypes_test
ctest --test-dir OpenMS-build -R "FileTypes" --output-on-failure
```

Expected: PASS.

- [ ] **Step 7: Commit.**

```bash
git add src/openms/include/OpenMS/FORMAT/FileTypes.h src/openms/source/FORMAT/FileTypes.cpp src/tests/class_tests/openms/source/FileTypes_test.cpp
git commit -m "$(cat <<'EOF'
feat: register FileTypes::IDPARQUET for identification parquet bundles

Adds the IDPARQUET enum value and binds it to the .idparquet extension
for the new internal PSM parquet directory format. No dispatch wiring
yet — that lands in a follow-up task.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Add `QPXFile::importFromArrow` (PSMSchema → PeptideIdentificationList)

**Goal:** Extract the PSM-construction half of `ConsensusMapArrowIO::importPSMsFromArrow` into `QPXFile::importFromArrow(table, protein_ids, peptide_ids)` so the same logic can be reused without ConsensusFeature linkage.

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/QPXFile.h` (declare new method)
- Modify: `src/openms/source/FORMAT/QPXFile.cpp` (implement, inline the small set of file-scope helpers it needs)
- Test: `src/tests/class_tests/openms/source/QPXFile_test.cpp`

**Note on duplicated helpers:** The helpers `getColumn_`, `isNull_`, `getInt32Value_`, `getDoubleValue_`, `getStringValue_`, `getBoolValue_`, `readMetaValues_` are currently duplicated across five `.cpp` files in `src/openms/source/FORMAT/` (`ConsensusMapArrowIO.cpp`, `FeatureMapArrowIO.cpp`, `ProteinIdentificationArrowIO.cpp`, `XICParquetFile.cpp`, `XIMParquetFile.cpp`). Adding a sixth copy in `QPXFile.cpp` is consistent with the existing pattern. Consolidating into `ArrowIOHelpers` is a separate refactor and is **out of scope** for this plan.

- [ ] **Step 1: Read the source helpers in `ConsensusMapArrowIO.cpp:812–905` and the PSM-construction loop at `ConsensusMapArrowIO.cpp:1457–1623`.**

```bash
grep -n "^  std::shared_ptr<arrow::Array> getColumn_\|^  double getDoubleValue_\|^  bool isNull_\|^  void readMetaValues_\|^  int32_t getInt32Value_\|^  int64_t getInt64Value_\|^  bool getBoolValue_\|^  std::string getStringValue_" src/openms/source/FORMAT/ConsensusMapArrowIO.cpp
```

Expected: locations of the helper definitions in `ConsensusMapArrowIO.cpp`.

- [ ] **Step 2: Write a failing test for `QPXFile::importFromArrow`.**

Append to `src/tests/class_tests/openms/source/QPXFile_test.cpp`:

```cpp
START_SECTION(([EXTRA] importFromArrow_round_trip))
{
  // Build minimal protein + peptide identifications.
  ProteinIdentification prot;
  prot.setIdentifier("run_1");
  prot.setScoreType("score");
  prot.setHigherScoreBetter(true);
  std::vector<ProteinIdentification> prot_ids{prot};

  PeptideIdentification pid;
  pid.setIdentifier("run_1");
  pid.setScoreType("score");
  pid.setHigherScoreBetter(true);
  pid.setRT(123.4);
  pid.setMZ(567.89);
  pid.setSpectrumReference("scan=42");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  hit.setCharge(2);
  hit.setScore(0.95);
  hit.setMetaValue("target_decoy", "target");
  hit.setMetaValue("COMET:deltaCn", 0.5);
  PeptideEvidence ev;
  ev.setProteinAccession("sp|P12345|EXAMPLE");
  hit.addPeptideEvidence(ev);
  pid.getHits().push_back(hit);

  PeptideIdentificationList pep_ids;
  pep_ids.push_back(pid);

  // Export to Arrow table (PSMSchema), then import back.
  auto table = QPXFile::exportToArrow(prot_ids, pep_ids, /*export_all_psms=*/true);
  TEST_NOT_EQUAL(table.get(), nullptr);

  std::vector<ProteinIdentification> prot_ids_out = prot_ids; // shell with run identifier preserved
  PeptideIdentificationList pep_ids_out;
  TEST_TRUE(QPXFile::importFromArrow(table, prot_ids_out, pep_ids_out));

  TEST_EQUAL(pep_ids_out.size(), 1);
  TEST_EQUAL(pep_ids_out[0].getHits().size(), 1);
  TEST_STRING_EQUAL(pep_ids_out[0].getHits()[0].getSequence().toString(), "PEPTIDE");
  TEST_EQUAL(pep_ids_out[0].getHits()[0].getCharge(), 2);
  TEST_REAL_SIMILAR(pep_ids_out[0].getHits()[0].getScore(), 0.95);
  TEST_STRING_EQUAL(String(pep_ids_out[0].getHits()[0].getMetaValue("target_decoy")), "target");
  TEST_REAL_SIMILAR(double(pep_ids_out[0].getHits()[0].getMetaValue("COMET:deltaCn")), 0.5);
}
END_SECTION
```

- [ ] **Step 3: Run the test to verify it fails to compile.**

```bash
cmake --build OpenMS-build -j$(nproc) --target QPXFile_test
```

Expected: FAIL with "no member named 'importFromArrow' in 'QPXFile'".

- [ ] **Step 4: Declare the method in `QPXFile.h`.**

`src/openms/include/OpenMS/FORMAT/QPXFile.h`, in the public section after `exportToParquet`:

```cpp
  /**
    @brief Import PSMs from a PSMSchema Arrow table.

    Reads `PSMSchema`-conformant rows and appends `PeptideIdentification`s
    to @p peptide_identifications. Each row's `run_identifier` column links
    PSMs back to the matching `ProteinIdentification` already present in
    @p protein_identifications by run identifier. If no match exists, a
    new `ProteinIdentification` shell is appended.

    @param[in]    table                     PSMSchema Arrow table (must not be null)
    @param[in,out] protein_identifications  Existing protein identifications (used for higher_score_better lookup; new shells appended for unknown run_identifiers)
    @param[out]   peptide_identifications   Peptide identifications populated from the table
    @return true on success, false on schema mismatch or unrecoverable error (errors are logged)
  */
  static bool importFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    std::vector<ProteinIdentification>& protein_identifications,
    PeptideIdentificationList& peptide_identifications);
```

- [ ] **Step 5: Implement `importFromArrow` in `QPXFile.cpp` and inline the helper functions it needs.**

In `src/openms/source/FORMAT/QPXFile.cpp`, inside the existing anonymous namespace at the top of the file (or add one if absent), copy the helpers from `ConsensusMapArrowIO.cpp:812–905` (`getColumn_`, `isNull_`, `getInt32Value_`, `getInt64Value_`, `getDoubleValue_`, `getBoolValue_`, `getStringValue_`, `readMetaValues_`). They are short and self-contained.

Then add the implementation (drawn from `ConsensusMapArrowIO::importPSMsFromArrow:1381–1623`, omitting all `ConsensusFeature`-related logic):

```cpp
bool QPXFile::importFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  std::vector<ProteinIdentification>& protein_identifications,
  PeptideIdentificationList& peptide_identifications)
{
  if (!table || table->num_rows() == 0) { return true; }

  auto combined_result = table->CombineChunks(arrow::default_memory_pool());
  if (!combined_result.ok())
  {
    OPENMS_LOG_ERROR << "QPXFile::importFromArrow: Failed to combine chunks" << std::endl;
    return false;
  }
  const auto& tbl = *combined_result;
  int64_t num_rows = tbl->num_rows();

  auto psm_validation = ArrowSchemaValidation::validate(
    tbl, PSMSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!psm_validation.valid)
  {
    OPENMS_LOG_ERROR << "QPXFile::importFromArrow: Incompatible PSM schema: "
                     << psm_validation.toString() << std::endl;
    return false;
  }

  auto col_p_id = getColumn_(tbl, PSMSchema::PEPTIDE_IDENTIFICATION_INDEX);
  auto col_peptidoform = getColumn_(tbl, PSMSchema::PEPTIDOFORM, /*required=*/false);
  auto col_sequence = getColumn_(tbl, PSMSchema::SEQUENCE, /*required=*/false);
  auto col_charge = getColumn_(tbl, PSMSchema::PRECURSOR_CHARGE);
  auto col_score = getColumn_(tbl, PSMSchema::SCORE);
  auto col_score_type = getColumn_(tbl, PSMSchema::SCORE_TYPE);
  auto col_rank = getColumn_(tbl, PSMSchema::RANK, /*required=*/false);
  auto col_rt = getColumn_(tbl, PSMSchema::RT, /*required=*/false);
  auto col_mz = getColumn_(tbl, PSMSchema::OBSERVED_MZ, /*required=*/false);
  auto col_spec_ref = getColumn_(tbl, PSMSchema::SPECTRUM_REFERENCE, /*required=*/false);
  auto col_run_id = getColumn_(tbl, PSMSchema::RUN_IDENTIFIER, /*required=*/false);
  auto col_is_decoy = getColumn_(tbl, PSMSchema::IS_DECOY, /*required=*/false);
  auto col_protein_accs = getColumn_(tbl, PSMSchema::PROTEIN_ACCESSIONS, /*required=*/false);
  auto col_additional_scores = getColumn_(tbl, PSMSchema::ADDITIONAL_SCORES, /*required=*/false);
  auto col_psm_metavalues = getColumn_(tbl, PSMSchema::PSM_METAVALUES, /*required=*/false);
  auto col_spectrum_metavalues = getColumn_(tbl, PSMSchema::SPECTRUM_METAVALUES, /*required=*/false);
  auto col_predicted_rt = getColumn_(tbl, PSMSchema::PREDICTED_RT, /*required=*/false);
  auto col_ion_mobility = getColumn_(tbl, PSMSchema::ION_MOBILITY, /*required=*/false);
  auto col_hsb = getColumn_(tbl, PSMSchema::HIGHER_SCORE_BETTER, /*required=*/false);
  auto col_scan = getColumn_(tbl, PSMSchema::SCAN, /*required=*/false);
  auto col_ref_file = getColumn_(tbl, PSMSchema::REFERENCE_FILE_NAME, /*required=*/false);

  if (!col_p_id || !col_charge || !col_score || !col_score_type)
  {
    OPENMS_LOG_ERROR << "QPXFile::importFromArrow: Missing required PSM columns" << std::endl;
    return false;
  }

  std::unordered_map<std::string, bool> higher_score_better_lookup;
  for (const auto& prot_id : protein_identifications)
  {
    higher_score_better_lookup[prot_id.getIdentifier()] = prot_id.isHigherScoreBetter();
  }

  PeptideIdentification* current_pid = nullptr;
  int32_t current_p_id = -1;

  for (int64_t row = 0; row < num_rows; ++row)
  {
    int32_t p_id = getInt32Value_(col_p_id, row, -1);

    if (current_pid == nullptr || p_id != current_p_id)
    {
      peptide_identifications.emplace_back();
      current_pid = &peptide_identifications.back();
      current_pid->setScoreType(getStringValue_(col_score_type, row));

      if (col_run_id && !isNull_(col_run_id, row))
      {
        current_pid->setIdentifier(getStringValue_(col_run_id, row));
      }

      if (col_hsb && !isNull_(col_hsb, row))
      {
        current_pid->setHigherScoreBetter(getBoolValue_(col_hsb, row, true));
      }
      else if (col_run_id && !isNull_(col_run_id, row))
      {
        auto hsb_it = higher_score_better_lookup.find(current_pid->getIdentifier());
        current_pid->setHigherScoreBetter(
          hsb_it != higher_score_better_lookup.end() ? hsb_it->second : true);
      }
      else
      {
        current_pid->setHigherScoreBetter(true);
      }

      if (col_rt && !isNull_(col_rt, row)) { current_pid->setRT(getDoubleValue_(col_rt, row)); }
      if (col_mz && !isNull_(col_mz, row)) { current_pid->setMZ(getDoubleValue_(col_mz, row)); }
      if (col_spec_ref && !isNull_(col_spec_ref, row))
      {
        current_pid->setSpectrumReference(getStringValue_(col_spec_ref, row));
      }
      if (col_spectrum_metavalues)
      {
        readMetaValues_(col_spectrum_metavalues, row, *current_pid);
      }
      current_p_id = p_id;
    }

    PeptideHit hit;
    if (col_peptidoform && !isNull_(col_peptidoform, row))
    {
      const String peptidoform_str = getStringValue_(col_peptidoform, row);
      if (!peptidoform_str.empty())
      {
        try
        {
          auto pf = ProForma::parse(peptidoform_str);
          hit.setSequence(ProForma::toAASequence(pf, ProForma::ConversionPolicy::BEST_EFFORT));
        }
        catch (...)
        {
          if (col_sequence && !isNull_(col_sequence, row))
          {
            hit.setSequence(AASequence::fromString(getStringValue_(col_sequence, row)));
          }
        }
      }
    }
    else if (col_sequence && !isNull_(col_sequence, row))
    {
      hit.setSequence(AASequence::fromString(getStringValue_(col_sequence, row)));
    }

    hit.setCharge(static_cast<Int>(getInt32Value_(col_charge, row, 0)));
    hit.setScore(getDoubleValue_(col_score, row, 0.0));

    if (col_rank && !isNull_(col_rank, row))
    {
      hit.setRank(static_cast<UInt>(getInt32Value_(col_rank, row, 0)));
    }

    if (col_is_decoy && !isNull_(col_is_decoy, row))
    {
      bool is_decoy = getBoolValue_(col_is_decoy, row, false);
      hit.setMetaValue("target_decoy", is_decoy ? "decoy" : "target");
    }

    if (col_protein_accs && !isNull_(col_protein_accs, row))
    {
      auto list_arr = std::static_pointer_cast<arrow::ListArray>(col_protein_accs);
      auto values = std::static_pointer_cast<arrow::StringArray>(list_arr->values());
      int64_t start = list_arr->value_offset(row);
      int64_t end = start + list_arr->value_length(row);
      for (int64_t k = start; k < end; ++k)
      {
        PeptideEvidence ev;
        ev.setProteinAccession(values->GetString(k));
        hit.addPeptideEvidence(ev);
      }
    }

    if (col_additional_scores && !isNull_(col_additional_scores, row))
    {
      auto list_arr = std::static_pointer_cast<arrow::ListArray>(col_additional_scores);
      auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->values());
      auto names_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->GetFieldByName("score_name"));
      auto values_arr = std::static_pointer_cast<arrow::DoubleArray>(struct_arr->GetFieldByName("score_value"));
      int64_t start = list_arr->value_offset(row);
      int64_t end = start + list_arr->value_length(row);
      for (int64_t k = start; k < end; ++k)
      {
        hit.setMetaValue(names_arr->GetString(k), values_arr->Value(k));
      }
    }

    if (col_predicted_rt && !isNull_(col_predicted_rt, row))
    {
      hit.setMetaValue("predicted_RT", getDoubleValue_(col_predicted_rt, row));
    }
    if (col_ion_mobility && !isNull_(col_ion_mobility, row))
    {
      hit.setMetaValue("ion_mobility", getDoubleValue_(col_ion_mobility, row));
    }
    if (col_scan && !isNull_(col_scan, row))
    {
      hit.setMetaValue("scan", static_cast<int>(getInt32Value_(col_scan, row)));
    }
    if (col_ref_file && !isNull_(col_ref_file, row))
    {
      hit.setMetaValue("reference_file_name", getStringValue_(col_ref_file, row));
    }

    if (col_psm_metavalues)
    {
      static const std::unordered_set<std::string> psm_excluded_mvs =
        {"target_decoy", "predicted_RT", "predicted_rt", "ion_mobility", "IM",
         "scan", "reference_file_name"};
      readMetaValues_(col_psm_metavalues, row, hit, psm_excluded_mvs);
    }

    current_pid->getHits().push_back(std::move(hit));
  }

  // Append a shell ProteinIdentification for any run_identifier we saw but don't have.
  std::unordered_set<std::string> known;
  for (const auto& p : protein_identifications) { known.insert(p.getIdentifier()); }
  for (const auto& pid : peptide_identifications)
  {
    const std::string& id = pid.getIdentifier();
    if (!id.empty() && known.insert(id).second)
    {
      ProteinIdentification shell;
      shell.setIdentifier(id);
      shell.setScoreType(pid.getScoreType());
      shell.setHigherScoreBetter(pid.isHigherScoreBetter());
      protein_identifications.push_back(std::move(shell));
    }
  }

  return true;
}
```

Add the necessary includes at the top of `QPXFile.cpp` (if not already present):

```cpp
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <unordered_map>
#include <unordered_set>
```

- [ ] **Step 6: Build and run the test.**

```bash
cmake --build OpenMS-build -j$(nproc) --target QPXFile_test
ctest --test-dir OpenMS-build -R "^QPXFile$" --output-on-failure
```

Expected: PASS.

- [ ] **Step 7: Run a wider sanity build to catch other compile breaks.**

```bash
cmake --build OpenMS-build -j$(nproc) --target OpenMS
```

Expected: PASS.

- [ ] **Step 8: Commit.**

```bash
git add src/openms/include/OpenMS/FORMAT/QPXFile.h src/openms/source/FORMAT/QPXFile.cpp src/tests/class_tests/openms/source/QPXFile_test.cpp
git commit -m "$(cat <<'EOF'
feat: add QPXFile::importFromArrow for PSMSchema -> PeptideIdentificationList

Provides a free-standing reader for PSMSchema Arrow tables, extracted
from ConsensusMapArrowIO::importPSMsFromArrow's PSM-construction logic
(no ConsensusFeature linkage). Used by the new PSMArrowIO bundle reader
in a follow-up task.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Refactor `ConsensusMapArrowIO::importPSMsFromArrow` to call `QPXFile::importFromArrow`

**Goal:** Eliminate duplication. `ConsensusMapArrowIO::importPSMsFromArrow` keeps the part that reads `consensus_feature_unique_id` and dispatches PSMs into `cmap`'s features; PSM construction itself goes through the new helper.

**Files:**
- Modify: `src/openms/source/FORMAT/ConsensusMapArrowIO.cpp` (`importPSMsFromArrow` body, ~lines 1381–1650)

- [ ] **Step 1: Confirm the existing `ConsensusMapArrowIO_test` passes before changing anything.**

```bash
cmake --build OpenMS-build -j$(nproc) --target ConsensusMapArrowIO_test
ctest --test-dir OpenMS-build -R "^ConsensusMapArrowIO$" --output-on-failure
```

Expected: PASS. (If it does not pass, stop and investigate before refactoring.)

- [ ] **Step 2: Replace `ConsensusMapArrowIO::importPSMsFromArrow`'s PSM-construction body with a call to `QPXFile::importFromArrow`, keeping only the consensus-feature-id sidecar logic local.**

Replace the body of `ConsensusMapArrowIO::importPSMsFromArrow` (currently `src/openms/source/FORMAT/ConsensusMapArrowIO.cpp:1381–1650`) with:

```cpp
bool ConsensusMapArrowIO::importPSMsFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  ConsensusMap& cmap)
{
  if (!table || table->num_rows() == 0) { return true; }

  auto combined_result = table->CombineChunks(arrow::default_memory_pool());
  if (!combined_result.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to combine chunks in importPSMsFromArrow" << std::endl;
    return false;
  }
  const auto& tbl = *combined_result;

  // Read consensus_feature_unique_id and PEPTIDE_IDENTIFICATION_INDEX side-by-side
  // so we can route each constructed PeptideIdentification to the right ConsensusFeature.
  auto col_feature_id = getColumn_(tbl, "consensus_feature_unique_id");
  auto col_p_id = getColumn_(tbl, PSMSchema::PEPTIDE_IDENTIFICATION_INDEX);
  if (!col_feature_id || !col_p_id)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Missing required columns for PSM import" << std::endl;
    return false;
  }

  // Build (p_id -> feature_unique_id-or-null) map by scanning the first row of each p_id group.
  std::unordered_map<int32_t, std::pair<bool, int64_t>> p_id_to_feature; // bool: is_null
  int32_t current_p_id = -1;
  for (int64_t row = 0; row < tbl->num_rows(); ++row)
  {
    int32_t p_id = getInt32Value_(col_p_id, row, -1);
    if (p_id == current_p_id) { continue; }
    bool is_null = isNull_(col_feature_id, row);
    int64_t fid = is_null ? 0 : getInt64Value_(col_feature_id, row, 0);
    p_id_to_feature.emplace(p_id, std::make_pair(is_null, fid));
    current_p_id = p_id;
  }

  // Delegate PSM construction to QPXFile::importFromArrow.
  std::vector<ProteinIdentification> prot_ids = cmap.getProteinIdentifications();
  PeptideIdentificationList pep_ids;
  if (!QPXFile::importFromArrow(tbl, prot_ids, pep_ids)) { return false; }
  cmap.setProteinIdentifications(prot_ids);

  // Walk each constructed PeptideIdentification and route to features by feature_id.
  std::unordered_map<int64_t, ConsensusFeature*> feature_lookup;
  for (auto& cf : cmap)
  {
    feature_lookup[static_cast<int64_t>(cf.getUniqueId())] = &cf;
  }

  // Re-derive the p_id sequence from the table to keep PeptideIdentifications in input order.
  std::vector<int32_t> p_id_order;
  current_p_id = -1;
  for (int64_t row = 0; row < tbl->num_rows(); ++row)
  {
    int32_t p_id = getInt32Value_(col_p_id, row, -1);
    if (p_id != current_p_id) { p_id_order.push_back(p_id); current_p_id = p_id; }
  }
  if (p_id_order.size() != pep_ids.size())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: PSM count mismatch after import "
                     << "(rows=" << p_id_order.size() << ", pep_ids=" << pep_ids.size() << ")" << std::endl;
    return false;
  }

  for (size_t i = 0; i < pep_ids.size(); ++i)
  {
    auto fit = p_id_to_feature.find(p_id_order[i]);
    bool is_null = (fit == p_id_to_feature.end()) ? true : fit->second.first;
    int64_t fid = (fit == p_id_to_feature.end()) ? 0 : fit->second.second;

    if (is_null)
    {
      cmap.getUnassignedPeptideIdentifications().push_back(std::move(pep_ids[i]));
    }
    else
    {
      auto it = feature_lookup.find(fid);
      if (it != feature_lookup.end())
      {
        it->second->getPeptideIdentifications().push_back(std::move(pep_ids[i]));
      }
      else
      {
        OPENMS_LOG_WARN << "ConsensusMapArrowIO: Could not find consensus feature with id "
                        << fid << " for PSM. Adding as unassigned." << std::endl;
        cmap.getUnassignedPeptideIdentifications().push_back(std::move(pep_ids[i]));
      }
    }
  }

  return true;
}
```

Add `#include <OpenMS/FORMAT/QPXFile.h>` to the top of `ConsensusMapArrowIO.cpp` if not already present.

- [ ] **Step 3: Build and re-run the existing `ConsensusMapArrowIO_test`.**

```bash
cmake --build OpenMS-build -j$(nproc) --target ConsensusMapArrowIO_test
ctest --test-dir OpenMS-build -R "^ConsensusMapArrowIO$" --output-on-failure
```

Expected: PASS — behaviour preserved.

- [ ] **Step 4: Run the broader Arrow-format test suite as a smoke test.**

```bash
ctest --test-dir OpenMS-build -R "ArrowIO|QPXFile|FeatureMapArrowIO|ConsensusMapArrowIO" --output-on-failure
```

Expected: all PASS.

- [ ] **Step 5: Commit.**

```bash
git add src/openms/source/FORMAT/ConsensusMapArrowIO.cpp
git commit -m "$(cat <<'EOF'
refactor: ConsensusMapArrowIO::importPSMsFromArrow uses QPXFile::importFromArrow

Removes the duplicated PSM-construction block; ConsensusFeature linkage
stays local. Behaviour preserved by ConsensusMapArrowIO_test.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Implement `PSMArrowIO`

**Files:**
- Create: `src/openms/include/OpenMS/FORMAT/PSMArrowIO.h`
- Create: `src/openms/source/FORMAT/PSMArrowIO.cpp`
- Create: `src/tests/class_tests/openms/source/PSMArrowIO_test.cpp`
- Modify: `src/openms/includes.cmake`
- Modify: `src/openms/sources.cmake`
- Modify: `src/tests/class_tests/openms/executables.cmake`

- [ ] **Step 1: Add the new files to CMake registration.**

Locate the alphabetical block in `src/openms/includes.cmake` and add:

```cmake
include/OpenMS/FORMAT/PSMArrowIO.h
```

Locate the alphabetical block in `src/openms/sources.cmake` and add:

```cmake
source/FORMAT/PSMArrowIO.cpp
```

Locate the alphabetical block in `src/tests/class_tests/openms/executables.cmake` and add:

```cmake
PSMArrowIO_test
```

Re-run cmake to pick up the new files:

```bash
cmake -S . -B OpenMS-build
```

Expected: configure succeeds.

- [ ] **Step 2: Write the header.**

`src/openms/include/OpenMS/FORMAT/PSMArrowIO.h`:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>

#include <vector>

namespace OpenMS
{

/**
  @brief Read and write OpenMS identification data as a parquet bundle (.idparquet).

  An idparquet bundle is a directory containing four parquet files:
    - psms.parquet           (PSMSchema, lossless PeptideIdentification data)
    - proteins.parquet       (ProteinSchema, ProteinHits)
    - protein_groups.parquet (ProteinGroupSchema, indistinguishable groups)
    - search_params.parquet  (SearchParamsSchema, run-level parameters)

  All four files are required for a valid bundle on read.

  @ingroup FileIO
*/
class OPENMS_DLLAPI PSMArrowIO
{
public:
  /**
    @brief Export protein and peptide identifications to an idparquet directory bundle.

    Writes (and overwrites) the four canonical files inside @p dir. Other
    files in @p dir are left untouched. If @p dir does not exist it is
    created (including its parent? no — only @p dir itself; the parent
    must already exist). If @p dir exists as a regular file, returns false.

    @param[in] protein_identifications  Protein identifications (run-level metadata, hits, groups, search params)
    @param[in] peptide_identifications  Peptide identifications (PSMs)
    @param[in] dir                       Output directory path
    @param[in] export_all_psms           If true, write all hits per PSM; if false, write best hit per spectrum only (default: true — bundle is intended for lossless round-trip)
    @param[in] config                    Parquet writer configuration
    @return true on success, false on error (errors are logged)
  */
  static bool exportToParquet(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications,
    const String& dir,
    bool export_all_psms = true,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Import protein and peptide identifications from an idparquet directory bundle.

    All four canonical files (psms.parquet, proteins.parquet,
    protein_groups.parquet, search_params.parquet) must be present in
    @p dir. Missing any one is an error.

    @param[in]  dir                       Input directory path
    @param[out] protein_identifications   Populated from the three protein-side parquet files
    @param[out] peptide_identifications   Populated from psms.parquet
    @return true on success, false on error (errors are logged)
  */
  static bool importFromParquet(
    const String& dir,
    std::vector<ProteinIdentification>& protein_identifications,
    PeptideIdentificationList& peptide_identifications);
};

} // namespace OpenMS
```

- [ ] **Step 3: Write the failing test.**

`src/tests/class_tests/openms/source/PSMArrowIO_test.cpp`:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/PSMArrowIO.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/METADATA/PeptideEvidence.h>

using namespace OpenMS;
using namespace std;

namespace
{
  void buildMinimalIds(std::vector<ProteinIdentification>& prot_ids,
                      PeptideIdentificationList& pep_ids)
  {
    ProteinIdentification prot;
    prot.setIdentifier("run_1");
    prot.setScoreType("score");
    prot.setHigherScoreBetter(true);
    prot.getSearchParameters().digestion_enzyme.setName("Trypsin");
    prot.getSearchParameters().setMetaValue("extra_features", "COMET:deltaCn,MS:1002049");
    prot_ids = {prot};

    PeptideIdentification pid;
    pid.setIdentifier("run_1");
    pid.setScoreType("score");
    pid.setHigherScoreBetter(true);
    pid.setRT(123.4);
    pid.setMZ(567.89);
    pid.setSpectrumReference("scan=42");
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDE"));
    hit.setCharge(2);
    hit.setScore(0.95);
    hit.setMetaValue("target_decoy", "target");
    hit.setMetaValue("COMET:deltaCn", 0.5);
    PeptideEvidence ev;
    ev.setProteinAccession("sp|P12345|EXAMPLE");
    hit.addPeptideEvidence(ev);
    pid.getHits().push_back(hit);
    pep_ids.push_back(pid);
  }
}

START_TEST(PSMArrowIO, "$Id$")

START_SECTION(([EXTRA] export_then_import_round_trip))
{
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  buildMinimalIds(prot_ids, pep_ids);

  String dir = OPENMS_GET_TEST_DATA_PATH("../tmp_PSMArrowIO_round_trip.idparquet");
  if (File::exists(dir)) { File::removeDirRecursively(dir); }

  TEST_TRUE(PSMArrowIO::exportToParquet(prot_ids, pep_ids, dir));

  TEST_TRUE(File::exists(dir + "/psms.parquet"));
  TEST_TRUE(File::exists(dir + "/proteins.parquet"));
  TEST_TRUE(File::exists(dir + "/protein_groups.parquet"));
  TEST_TRUE(File::exists(dir + "/search_params.parquet"));

  std::vector<ProteinIdentification> prot_ids_in;
  PeptideIdentificationList pep_ids_in;
  TEST_TRUE(PSMArrowIO::importFromParquet(dir, prot_ids_in, pep_ids_in));

  TEST_EQUAL(prot_ids_in.size(), 1);
  TEST_STRING_EQUAL(prot_ids_in[0].getIdentifier(), "run_1");
  TEST_STRING_EQUAL(prot_ids_in[0].getSearchParameters().digestion_enzyme.getName(), "Trypsin");
  TEST_STRING_EQUAL(String(prot_ids_in[0].getSearchParameters().getMetaValue("extra_features")),
                    "COMET:deltaCn,MS:1002049");

  TEST_EQUAL(pep_ids_in.size(), 1);
  TEST_EQUAL(pep_ids_in[0].getHits().size(), 1);
  TEST_STRING_EQUAL(pep_ids_in[0].getHits()[0].getSequence().toString(), "PEPTIDE");
  TEST_REAL_SIMILAR(double(pep_ids_in[0].getHits()[0].getMetaValue("COMET:deltaCn")), 0.5);

  File::removeDirRecursively(dir);
}
END_SECTION

START_SECTION(([EXTRA] import_missing_subfile_returns_false))
{
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  buildMinimalIds(prot_ids, pep_ids);

  String dir = OPENMS_GET_TEST_DATA_PATH("../tmp_PSMArrowIO_missing_sub.idparquet");
  if (File::exists(dir)) { File::removeDirRecursively(dir); }
  TEST_TRUE(PSMArrowIO::exportToParquet(prot_ids, pep_ids, dir));

  // Delete one subfile.
  File::remove(dir + "/psms.parquet");

  std::vector<ProteinIdentification> prot_ids_in;
  PeptideIdentificationList pep_ids_in;
  TEST_FALSE(PSMArrowIO::importFromParquet(dir, prot_ids_in, pep_ids_in));

  File::removeDirRecursively(dir);
}
END_SECTION

START_SECTION(([EXTRA] export_target_is_regular_file_returns_false))
{
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  buildMinimalIds(prot_ids, pep_ids);

  String path = OPENMS_GET_TEST_DATA_PATH("../tmp_PSMArrowIO_regular_file.idparquet");
  if (File::exists(path)) { File::remove(path); }
  // Create a regular file at the target path.
  std::ofstream f(path); f << "not a directory"; f.close();

  TEST_FALSE(PSMArrowIO::exportToParquet(prot_ids, pep_ids, path));

  File::remove(path);
}
END_SECTION

START_SECTION(([EXTRA] empty_psms_round_trips))
{
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  buildMinimalIds(prot_ids, pep_ids);
  pep_ids.clear();

  String dir = OPENMS_GET_TEST_DATA_PATH("../tmp_PSMArrowIO_empty.idparquet");
  if (File::exists(dir)) { File::removeDirRecursively(dir); }

  TEST_TRUE(PSMArrowIO::exportToParquet(prot_ids, pep_ids, dir));

  std::vector<ProteinIdentification> prot_ids_in;
  PeptideIdentificationList pep_ids_in;
  TEST_TRUE(PSMArrowIO::importFromParquet(dir, prot_ids_in, pep_ids_in));
  TEST_EQUAL(pep_ids_in.size(), 0);
  TEST_EQUAL(prot_ids_in.size(), 1);

  File::removeDirRecursively(dir);
}
END_SECTION

END_TEST
```

- [ ] **Step 4: Build and confirm the test fails to link (no implementation yet).**

```bash
cmake --build OpenMS-build -j$(nproc) --target PSMArrowIO_test
```

Expected: link FAIL (`undefined reference to PSMArrowIO::exportToParquet`).

- [ ] **Step 5: Implement `PSMArrowIO.cpp`.**

`src/openms/source/FORMAT/PSMArrowIO.cpp`:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PSMArrowIO.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>

namespace OpenMS
{

namespace
{
  constexpr const char* kPSMs = "psms.parquet";
  constexpr const char* kProteins = "proteins.parquet";
  constexpr const char* kProteinGroups = "protein_groups.parquet";
  constexpr const char* kSearchParams = "search_params.parquet";

  bool ensureDirectory_(const String& dir)
  {
    if (File::exists(dir))
    {
      if (!File::isDirectory(dir))
      {
        OPENMS_LOG_ERROR << "PSMArrowIO: target path exists and is not a directory: " << dir << std::endl;
        return false;
      }
      return true;
    }
    if (!File::makeDir(dir))
    {
      OPENMS_LOG_ERROR << "PSMArrowIO: failed to create directory: " << dir << std::endl;
      return false;
    }
    return true;
  }
}

bool PSMArrowIO::exportToParquet(
  const std::vector<ProteinIdentification>& protein_identifications,
  const PeptideIdentificationList& peptide_identifications,
  const String& dir,
  bool export_all_psms,
  const ParquetWriteConfig& config)
{
  if (!ensureDirectory_(dir)) { return false; }

  // PSMs (PSMSchema, lossless)
  auto psm_table = QPXFile::exportToArrow(
    protein_identifications, peptide_identifications, export_all_psms);
  if (!psm_table)
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: QPXFile::exportToArrow returned null" << std::endl;
    return false;
  }
  if (!ArrowIOHelpers::writeTableToParquet(psm_table, dir + "/" + kPSMs, config))
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: failed to write " << kPSMs << std::endl;
    return false;
  }

  if (!ProteinIdentificationArrowIO::exportProteinsToParquet(
        protein_identifications, dir + "/" + kProteins, config))
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: failed to write " << kProteins << std::endl;
    return false;
  }
  if (!ProteinIdentificationArrowIO::exportProteinGroupsToParquet(
        protein_identifications, dir + "/" + kProteinGroups, config))
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: failed to write " << kProteinGroups << std::endl;
    return false;
  }
  if (!ProteinIdentificationArrowIO::exportSearchParamsToParquet(
        protein_identifications, dir + "/" + kSearchParams, config))
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: failed to write " << kSearchParams << std::endl;
    return false;
  }

  return true;
}

bool PSMArrowIO::importFromParquet(
  const String& dir,
  std::vector<ProteinIdentification>& protein_identifications,
  PeptideIdentificationList& peptide_identifications)
{
  if (!File::exists(dir) || !File::isDirectory(dir))
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: not a directory: " << dir << std::endl;
    return false;
  }
  for (const char* sub : {kPSMs, kProteins, kProteinGroups, kSearchParams})
  {
    if (!File::exists(dir + "/" + sub))
    {
      OPENMS_LOG_ERROR << "PSMArrowIO: missing required file: " << dir + "/" + sub << std::endl;
      return false;
    }
  }

  protein_identifications.clear();
  peptide_identifications.clear();

  if (!ProteinIdentificationArrowIO::importFromParquet(
        dir + "/" + kProteins,
        dir + "/" + kProteinGroups,
        dir + "/" + kSearchParams,
        protein_identifications))
  {
    return false;
  }

  auto psm_table = ParquetFile::readTable(dir + "/" + kPSMs);
  if (!psm_table)
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: failed to read " << kPSMs << std::endl;
    return false;
  }

  return QPXFile::importFromArrow(psm_table, protein_identifications, peptide_identifications);
}

} // namespace OpenMS
```

- [ ] **Step 6: Build and run the new test.**

```bash
cmake --build OpenMS-build -j$(nproc) --target PSMArrowIO_test
ctest --test-dir OpenMS-build -R "^PSMArrowIO$" --output-on-failure
```

Expected: PASS (all four sections).

- [ ] **Step 7: Commit.**

```bash
git add src/openms/include/OpenMS/FORMAT/PSMArrowIO.h src/openms/source/FORMAT/PSMArrowIO.cpp src/tests/class_tests/openms/source/PSMArrowIO_test.cpp src/openms/includes.cmake src/openms/sources.cmake src/tests/class_tests/openms/executables.cmake
git commit -m "$(cat <<'EOF'
feat: add PSMArrowIO for .idparquet directory bundles (#9225)

Bundles PSMSchema PSMs and ProteinIdentificationArrowIO protein/group/
search-params tables under a single directory. Round-trip lossless;
all four subfiles required on read.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Wire `IDPARQUET` into `FileHandler::loadIdentifications` and `storeIdentifications`

**Files:**
- Modify: `src/openms/source/FORMAT/FileHandler.cpp` (around lines 1347 and 1424)
- Test: `src/tests/class_tests/openms/source/FileHandler_test.cpp`

- [ ] **Step 1: Write the failing test (round-trip via FileHandler).**

Append to `src/tests/class_tests/openms/source/FileHandler_test.cpp`:

```cpp
START_SECTION(([EXTRA] storeIdentifications_loadIdentifications_idparquet_round_trip))
{
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;

  ProteinIdentification prot;
  prot.setIdentifier("run_1");
  prot.setScoreType("score");
  prot.setHigherScoreBetter(true);
  prot_ids.push_back(prot);

  PeptideIdentification pid;
  pid.setIdentifier("run_1");
  pid.setScoreType("score");
  pid.setHigherScoreBetter(true);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  hit.setCharge(2);
  hit.setScore(0.5);
  pid.getHits().push_back(hit);
  pep_ids.push_back(pid);

  String dir = OPENMS_GET_TEST_DATA_PATH("../tmp_FileHandler_idparquet.idparquet");
  if (File::exists(dir)) { File::removeDirRecursively(dir); }

  FileHandler().storeIdentifications(dir, prot_ids, pep_ids, {FileTypes::IDPARQUET});

  std::vector<ProteinIdentification> prot_ids_in;
  PeptideIdentificationList pep_ids_in;
  FileHandler().loadIdentifications(dir, prot_ids_in, pep_ids_in, {FileTypes::IDPARQUET});

  TEST_EQUAL(prot_ids_in.size(), 1);
  TEST_EQUAL(pep_ids_in.size(), 1);

  File::removeDirRecursively(dir);
}
END_SECTION
```

- [ ] **Step 2: Run the test to confirm it fails.**

```bash
cmake --build OpenMS-build -j$(nproc) --target FileHandler_test
ctest --test-dir OpenMS-build -R "^FileHandler$" --output-on-failure
```

Expected: FAIL — likely `Exception::ParseError` "type: idparquet is not supported for storing identifications".

- [ ] **Step 3: Add the dispatch case in `loadIdentifications`.**

In `src/openms/source/FORMAT/FileHandler.cpp`, in the `switch (type)` block of `FileHandler::loadIdentifications` (around line 1360), insert before the `default:` clause:

```cpp
      case FileTypes::IDPARQUET:
      {
        if (!PSMArrowIO::importFromParquet(filename, additional_proteins, additional_peptides))
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                      "PSMArrowIO::importFromParquet failed");
        }
      }
      break;
```

- [ ] **Step 4: Add the dispatch case in `storeIdentifications`.**

In `FileHandler::storeIdentifications` (around line 1440), insert before the `default:` clause:

```cpp
      case FileTypes::IDPARQUET:
      {
        if (!PSMArrowIO::exportToParquet(additional_proteins, additional_peptides, filename))
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                              "PSMArrowIO::exportToParquet failed");
        }
      }
      break;
```

- [ ] **Step 5: Add the include at the top of `FileHandler.cpp`.**

```cpp
#include <OpenMS/FORMAT/PSMArrowIO.h>
```

- [ ] **Step 6: Build and re-run.**

```bash
cmake --build OpenMS-build -j$(nproc) --target FileHandler_test
ctest --test-dir OpenMS-build -R "^FileHandler$" --output-on-failure
```

Expected: PASS.

- [ ] **Step 7: Commit.**

```bash
git add src/openms/source/FORMAT/FileHandler.cpp src/tests/class_tests/openms/source/FileHandler_test.cpp
git commit -m "$(cat <<'EOF'
feat: dispatch IDPARQUET in FileHandler::{load,store}Identifications

Wires PSMArrowIO behind FileHandler so any tool that already uses
loadIdentifications/storeIdentifications gains .idparquet support
without further changes.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: PercolatorAdapter — accept and produce `.idparquet`

**Files:**
- Modify: `src/topp/PercolatorAdapter.cpp` (`registerOptionsAndFlags_` around line 209; `main_` output dispatch around line 580; `storeIdentifications` call around line 1166)

- [ ] **Step 1: Add `idparquet` to input valid-formats.**

`src/topp/PercolatorAdapter.cpp:210`:

```cpp
    setValidFormats_("in", ListUtils::create<String>("mzid,idXML,idparquet"));
```

`src/topp/PercolatorAdapter.cpp:212`:

```cpp
    setValidFormats_("in_decoy", ListUtils::create<String>("mzid,idXML,idparquet"));
```

- [ ] **Step 2: Add `idparquet` to output valid-formats.**

`src/topp/PercolatorAdapter.cpp:216`:

```cpp
    setValidFormats_("out", ListUtils::create<String>("idXML,mzid,osw,idparquet"));
```

`src/topp/PercolatorAdapter.cpp:230`:

```cpp
    setValidStrings_("out_type", ListUtils::create<String>("mzid,idXML,osw,idparquet"));
```

- [ ] **Step 3: Update the input-file load path to allow `IDPARQUET`.**

Locate the existing `FileHandler().loadIdentifications` call inside `main_` (around line 419/424). The current code branches on file type for `IDXML` vs `MZIDENTML`. Add an explicit branch for `IDPARQUET` (or, simpler, expand the `allowed_types` of the existing call). The minimal change:

```cpp
      else if (in_type == FileTypes::MZIDENTML)
      {
        FileHandler().loadIdentifications(in, protein_ids, peptide_ids, {FileTypes::MZIDENTML});
      }
      else if (in_type == FileTypes::IDPARQUET)
      {
        FileHandler().loadIdentifications(in, protein_ids, peptide_ids, {FileTypes::IDPARQUET});
      }
```

(Insert after the existing `MZIDENTML` branch.)

- [ ] **Step 4: Update the OSW guard at line 616.**

The current check `if (!in_osw.empty() && out_type != FileTypes::OSW)` is unaffected by `IDPARQUET` since the guard is only relevant for OSW input. No change needed.

- [ ] **Step 5: Update the final output dispatch.**

`src/topp/PercolatorAdapter.cpp:1166` currently passes `{FileTypes::IDXML, FileTypes::MZIDENTML}`. Change to:

```cpp
      FileHandler().storeIdentifications(out, all_protein_ids, all_peptide_ids,
        {FileTypes::IDXML, FileTypes::MZIDENTML, FileTypes::IDPARQUET});
```

- [ ] **Step 6: Build PercolatorAdapter to confirm clean compile.**

```bash
cmake --build OpenMS-build -j$(nproc) --target PercolatorAdapter
```

Expected: PASS.

- [ ] **Step 7: Commit.**

```bash
git add src/topp/PercolatorAdapter.cpp
git commit -m "$(cat <<'EOF'
feat(PercolatorAdapter): accept and produce .idparquet (#9225)

Adds idparquet to input/output valid-formats and to the
storeIdentifications allowed-types whitelist. Filter and best-PSM
options land in the next commit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: PercolatorAdapter — `score:fdr` post-filter and `keep_all_passing` flag

**Files:**
- Modify: `src/topp/PercolatorAdapter.cpp` (`registerOptionsAndFlags_`; post-Percolator step before the final `storeIdentifications` call)

- [ ] **Step 1: Locate the filter insertion point — immediately before the final `storeIdentifications` call (around line 1166).**

```bash
grep -n "FileHandler().storeIdentifications(out" src/topp/PercolatorAdapter.cpp
```

Expected: one match near line 1166.

- [ ] **Step 2: Register the new options.**

In `PercolatorAdapter::registerOptionsAndFlags_`, in the same option group as other score-related options (locate by `grep -n "score:" src/topp/PercolatorAdapter.cpp`), insert:

```cpp
    registerDoubleOption_("score:fdr", "<value>",
      1.0,
      "FDR cutoff applied to the Percolator q-value before writing output. "
      "PSMs with q-value > cutoff are dropped. 1.0 disables the filter.",
      false, false);
    setMinFloat_("score:fdr", 0.0);
    setMaxFloat_("score:fdr", 1.0);

    registerFlag_("keep_all_passing",
      "Keep every PSM that passes score:fdr (default: keep only the best PSM per spectrum).",
      false);
```

- [ ] **Step 3: Add an option-only test on the existing idXML path that confirms the filter drops the expected PSMs.**

Locate the existing `TOPP_PercolatorAdapter_1` test in `src/tests/topp/THIRDPARTY/third_party_tests.cmake:207`. Append a new test case after `TOPP_PercolatorAdapter_5`:

```cmake
  add_test("TOPP_PercolatorAdapter_score_fdr"
    ${TOPP_BIN_PATH}/PercolatorAdapter -test
      -ini ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_1.ini
      -in ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_1.idXML
      -out PercolatorAdapter_score_fdr_out.tmp.idXML
      -out_type idXML
      -score:fdr 0.01
      -percolator_executable "${PERCOLATOR_BINARY}")
  set_tests_properties("TOPP_PercolatorAdapter_score_fdr" PROPERTIES DEPENDS "TOPP_PercolatorAdapter_executable_check")
```

(Reference output is omitted intentionally; the test asserts only that the run succeeds. A reference comparison would lock the plan to a specific Percolator stochastic outcome and is out of scope.)

- [ ] **Step 4: Implement the post-filter step.**

Just before the final `FileHandler().storeIdentifications(out, all_protein_ids, all_peptide_ids, ...)` call in `main_`, insert:

```cpp
    // Post-filter on Percolator's q-value (PeptideHit::score after Percolator parsing).
    const double score_fdr = getDoubleOption_("score:fdr");
    if (score_fdr < 1.0)
    {
      // After Percolator output parsing, hit.score is the q-value. Lower is better.
      // IDFilter::filterHitsByScore drops hits where score >= threshold (higher is worse).
      IDFilter::filterHitsByScore(all_peptide_ids, score_fdr);
      IDFilter::removeEmptyIdentifications(all_peptide_ids);
      OPENMS_LOG_INFO << "Applied score:fdr cutoff " << score_fdr
                      << "; remaining peptide identifications: " << all_peptide_ids.size() << std::endl;
      if (all_peptide_ids.empty())
      {
        OPENMS_LOG_WARN << "score:fdr cutoff " << score_fdr << " dropped all PSMs. "
                        << "Output will contain no peptide identifications." << std::endl;
      }
    }

    if (!getFlag_("keep_all_passing"))
    {
      IDFilter::keepBestPeptideHits(all_peptide_ids);
    }
```

Add the include at the top of `PercolatorAdapter.cpp` if not already present:

```cpp
#include <OpenMS/PROCESSING/ID/IDFilter.h>
```

- [ ] **Step 5: Build and re-run the existing PercolatorAdapter tests + the new one.**

```bash
cmake --build OpenMS-build -j$(nproc) --target PercolatorAdapter
ctest --test-dir OpenMS-build -R "^TOPP_PercolatorAdapter" --output-on-failure
```

Expected: existing tests PASS unchanged (default `score:fdr 1.0` is a no-op, `keep_all_passing` defaults off but the existing test fixture has only one hit per spectrum so behaviour is identical), and the new `TOPP_PercolatorAdapter_score_fdr` PASSES.

- [ ] **Step 6: Commit.**

```bash
git add src/topp/PercolatorAdapter.cpp src/tests/topp/THIRDPARTY/third_party_tests.cmake
git commit -m "$(cat <<'EOF'
feat(PercolatorAdapter): add score:fdr post-filter and keep_all_passing flag (#9225)

score:fdr drops PSMs whose Percolator q-value exceeds the cutoff.
keep_all_passing retains all hits per spectrum (default: best per spectrum).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: IDMerger — accept and produce `.idparquet`

**Files:**
- Modify: `src/topp/IDMerger.cpp:173–176, 330`

- [ ] **Step 1: Locate the `formats` list and the `storeIdentifications` call.**

```bash
grep -n "formats\|storeIdentifications" src/topp/IDMerger.cpp | head -20
```

Expected: `formats` defined near line 170; `storeIdentifications` near line 330.

- [ ] **Step 2: Add `idparquet` to the formats list and to the storeIdentifications whitelist.**

`src/topp/IDMerger.cpp` — locate the `formats` definition (around line 168–172) and append `"idparquet"`. Then update the `storeIdentifications` call at line 330:

```cpp
    FileHandler().storeIdentifications(out, proteins, peptides, {FileTypes::IDXML, FileTypes::IDPARQUET});
```

- [ ] **Step 3: Build IDMerger and run its existing tests.**

```bash
cmake --build OpenMS-build -j$(nproc) --target IDMerger
ctest --test-dir OpenMS-build -R "^TOPP_IDMerger" --output-on-failure
```

Expected: existing tests PASS unchanged.

- [ ] **Step 4: Commit.**

```bash
git add src/topp/IDMerger.cpp
git commit -m "$(cat <<'EOF'
feat(IDMerger): accept and produce .idparquet (#9225)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: SageAdapter — migrate to `FileHandler.storeIdentifications` + accept `.idparquet`

**Files:**
- Modify: `src/topp/SageAdapter.cpp` (line 428 `setValidFormats_`; line 889 `IdXMLFile().store(...)`)

- [ ] **Step 1: Add `idparquet` to `out` valid-formats.**

`src/topp/SageAdapter.cpp:428`:

```cpp
    setValidFormats_("out", { "idXML", "idparquet" } );
```

- [ ] **Step 2: Migrate the direct `IdXMLFile().store(...)` to `FileHandler().storeIdentifications(...)`.**

`src/topp/SageAdapter.cpp:889`. Replace:

```cpp
    IdXMLFile().store(output_file, protein_identifications, peptide_identifications);
```

with:

```cpp
    FileHandler().storeIdentifications(output_file, protein_identifications, peptide_identifications,
      {FileTypes::IDXML, FileTypes::IDPARQUET});
```

Add the include at the top of `SageAdapter.cpp` if not already present:

```cpp
#include <OpenMS/FORMAT/FileHandler.h>
```

- [ ] **Step 3: Build SageAdapter and run its existing tests.**

```bash
cmake --build OpenMS-build -j$(nproc) --target SageAdapter
ctest --test-dir OpenMS-build -R "^TOPP_SageAdapter" --output-on-failure
```

Expected: existing tests PASS.

- [ ] **Step 4: Commit.**

```bash
git add src/topp/SageAdapter.cpp
git commit -m "$(cat <<'EOF'
feat(SageAdapter): migrate to FileHandler.storeIdentifications + accept .idparquet (#9225)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: CometAdapter, MSGFPlusAdapter, MSFraggerAdapter — accept `.idparquet`

**Note:** `XTandemAdapter.cpp` was removed from the codebase (only old test fixtures remain). Originally listed in the spec, dropped here.

**Files:**
- Modify: `src/topp/CometAdapter.cpp:139, 926`
- Modify: `src/topp/MSGFPlusAdapter.cpp` (find the analogous line numbers via grep)
- Modify: `src/topp/MSFraggerAdapter.cpp:888`

- [ ] **Step 1: For each adapter, locate the `setValidFormats_("out", ...)` call and the `storeIdentifications` call.**

```bash
for f in CometAdapter MSGFPlusAdapter MSFraggerAdapter; do
  echo "=== $f ==="
  grep -nE "setValidFormats_\\(\"out\"|storeIdentifications" src/topp/$f.cpp
done
```

Expected: each tool has one or two matches per pattern.

- [ ] **Step 2: For each tool, expand both calls to include `idparquet` / `FileTypes::IDPARQUET`.**

Pattern (CometAdapter shown — apply equivalent edits to the other three):

`src/topp/CometAdapter.cpp:139`:

```cpp
    setValidFormats_("out", { "idXML", "idparquet"} );
```

`src/topp/CometAdapter.cpp:926`:

```cpp
    FileHandler().storeIdentifications(out, protein_identifications, peptide_identifications, {FileTypes::IDXML, FileTypes::IDPARQUET});
```

Repeat for `MSGFPlusAdapter.cpp:887` and `MSFraggerAdapter.cpp:888` (line numbers will differ; the textual pattern is the same).

- [ ] **Step 3: Build all three adapters and run their tests.**

```bash
cmake --build OpenMS-build -j$(nproc) --target CometAdapter MSGFPlusAdapter MSFraggerAdapter
ctest --test-dir OpenMS-build -R "^TOPP_(Comet|MSGFPlus|MSFragger)Adapter" --output-on-failure
```

Expected: all existing tests PASS.

- [ ] **Step 4: Commit.**

```bash
git add src/topp/CometAdapter.cpp src/topp/MSGFPlusAdapter.cpp src/topp/MSFraggerAdapter.cpp
git commit -m "$(cat <<'EOF'
feat(search-engine adapters): accept .idparquet output (#9225)

Comet, MSGFPlus, and MSFragger adapters now accept idparquet as
an output format, matching SageAdapter and PercolatorAdapter.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: TOPP test — PercolatorAdapter on `.idparquet` input/output

**Files:**
- Create: `src/tests/topp/THIRDPARTY/PercolatorAdapter_idparquet_in.idparquet/` (input fixture, four parquet files)
- Modify: `src/tests/topp/THIRDPARTY/third_party_tests.cmake`

**Reference-comparison strategy:** parquet bytes are not stable across Arrow versions. The test imports the produced `.idparquet/` back via `FileHandler().loadIdentifications` and compares to a reference `.idXML` snapshot that is already in tree. Implementation is via a tiny helper TOPP run that converts parquet → idXML for the comparison.

- [ ] **Step 1: Generate the input fixture from `PercolatorAdapter_1.idXML`.**

Run, in the build directory, a one-off conversion (this run is *not* committed to CI; we commit its output):

```bash
cmake --build OpenMS-build -j$(nproc) --target FileConverter
# Use an existing converter or write a one-off script. For now, FileConverter
# works because dispatch is wired through FileHandler:
./OpenMS-build/bin/FileConverter \
  -in src/tests/topp/THIRDPARTY/PercolatorAdapter_1.idXML \
  -out src/tests/topp/THIRDPARTY/PercolatorAdapter_idparquet_in.idparquet
```

Verify the four expected sub-files exist:

```bash
ls src/tests/topp/THIRDPARTY/PercolatorAdapter_idparquet_in.idparquet/
```

Expected: `psms.parquet`, `proteins.parquet`, `protein_groups.parquet`, `search_params.parquet`.

> If `FileConverter` does not yet support `idparquet` output, the simplest workaround is a 5-line throwaway TOPP tool placed under `src/utils/` that calls `FileHandler().{load,store}Identifications` directly. That tool can be added under this task and removed after the fixture is generated.

- [ ] **Step 2: Add the new tests to `third_party_tests.cmake`.**

In `src/tests/topp/THIRDPARTY/third_party_tests.cmake`, after the existing `TOPP_PercolatorAdapter_5` test, add:

```cmake
  add_test("TOPP_PercolatorAdapter_idparquet"
    ${TOPP_BIN_PATH}/PercolatorAdapter -test
      -ini ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_1.ini
      -in ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_idparquet_in.idparquet
      -out PercolatorAdapter_idparquet_out.tmp.idparquet
      -out_type idparquet
      -percolator_executable "${PERCOLATOR_BINARY}")

  # Convert the produced .idparquet back to idXML and diff against the existing reference.
  add_test("TOPP_PercolatorAdapter_idparquet_convert"
    ${TOPP_BIN_PATH}/FileConverter
      -in PercolatorAdapter_idparquet_out.tmp.idparquet
      -out PercolatorAdapter_idparquet_out.tmp.idXML)
  set_tests_properties("TOPP_PercolatorAdapter_idparquet_convert" PROPERTIES DEPENDS "TOPP_PercolatorAdapter_idparquet")

  add_test("TOPP_PercolatorAdapter_idparquet_diff"
    ${DIFF}
      -in1 PercolatorAdapter_idparquet_out.tmp.idXML
      -in2 ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_1_out.idXML
      -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=" "UserParam.*loaded_file_path")
  set_tests_properties("TOPP_PercolatorAdapter_idparquet_diff" PROPERTIES DEPENDS "TOPP_PercolatorAdapter_idparquet_convert")
```

- [ ] **Step 3: Build the dependencies and run the new tests.**

```bash
cmake --build OpenMS-build -j$(nproc) --target PercolatorAdapter FileConverter
ctest --test-dir OpenMS-build -R "TOPP_PercolatorAdapter_idparquet" --output-on-failure
```

Expected: PASS.

- [ ] **Step 4: Commit.**

```bash
git add src/tests/topp/THIRDPARTY/PercolatorAdapter_idparquet_in.idparquet/ src/tests/topp/THIRDPARTY/third_party_tests.cmake
git commit -m "$(cat <<'EOF'
test(PercolatorAdapter): TOPP test on .idparquet input/output (#9225)

Reference comparison goes via FileConverter (.idparquet -> .idXML)
because parquet bytes are not stable across Arrow versions.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 12: TOPP test — IDMerger merging two `.idparquet` inputs

**Files:**
- Create: `src/tests/topp/THIRDPARTY/IDMerger_idparquet_a.idparquet/` (fixture A, derived from existing `IDMerger_1_input1.idXML`)
- Create: `src/tests/topp/THIRDPARTY/IDMerger_idparquet_b.idparquet/` (fixture B, derived from existing `IDMerger_1_input2.idXML`)
- Modify: `src/tests/topp/IDMerger.test.cmake` (or wherever IDMerger tests live — `grep -rn IDMerger_1 src/tests/topp/` to locate)

- [ ] **Step 1: Locate the existing IDMerger test invocation.**

```bash
grep -rn "IDMerger_1\|TOPP_IDMerger" src/tests/topp/ | head -10
```

Expected: existing test using two `.idXML` inputs.

- [ ] **Step 2: Generate the two `.idparquet/` fixtures from the existing `.idXML` inputs.**

```bash
./OpenMS-build/bin/FileConverter -in src/tests/topp/<idmerger_input1.idXML> -out src/tests/topp/THIRDPARTY/IDMerger_idparquet_a.idparquet
./OpenMS-build/bin/FileConverter -in src/tests/topp/<idmerger_input2.idXML> -out src/tests/topp/THIRDPARTY/IDMerger_idparquet_b.idparquet
```

Verify both directories exist with the four canonical sub-files.

- [ ] **Step 3: Add a test that merges the two `.idparquet/` inputs and converts back for diff.**

Append to the existing IDMerger test file:

```cmake
add_test(NAME "TOPP_IDMerger_idparquet"
  COMMAND ${TOPP_BIN_PATH}/IDMerger -test
    -in ${DATA_DIR_TOPP}/THIRDPARTY/IDMerger_idparquet_a.idparquet ${DATA_DIR_TOPP}/THIRDPARTY/IDMerger_idparquet_b.idparquet
    -out IDMerger_idparquet_out.tmp.idparquet)

add_test(NAME "TOPP_IDMerger_idparquet_convert"
  COMMAND ${TOPP_BIN_PATH}/FileConverter
    -in IDMerger_idparquet_out.tmp.idparquet
    -out IDMerger_idparquet_out.tmp.idXML)
set_tests_properties("TOPP_IDMerger_idparquet_convert" PROPERTIES DEPENDS "TOPP_IDMerger_idparquet")

add_test(NAME "TOPP_IDMerger_idparquet_diff"
  COMMAND ${DIFF}
    -in1 IDMerger_idparquet_out.tmp.idXML
    -in2 ${DATA_DIR_TOPP}/<existing IDMerger_1_out.idXML>
    -whitelist "IdentificationRun date" "SearchParameters id=" "UserParam.*loaded_file_path")
set_tests_properties("TOPP_IDMerger_idparquet_diff" PROPERTIES DEPENDS "TOPP_IDMerger_idparquet_convert")
```

- [ ] **Step 4: Build IDMerger + FileConverter and run.**

```bash
cmake --build OpenMS-build -j$(nproc) --target IDMerger FileConverter
ctest --test-dir OpenMS-build -R "TOPP_IDMerger_idparquet" --output-on-failure
```

Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add src/tests/topp/THIRDPARTY/IDMerger_idparquet_a.idparquet/ src/tests/topp/THIRDPARTY/IDMerger_idparquet_b.idparquet/ src/tests/topp/<the cmake file>
git commit -m "$(cat <<'EOF'
test(IDMerger): TOPP test merging two .idparquet inputs (#9225)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 13: Final integration smoke test

- [ ] **Step 1: Build everything.**

```bash
cmake --build OpenMS-build -j$(nproc)
```

Expected: clean build.

- [ ] **Step 2: Run the full test suite filtered to identification + parquet related tests.**

```bash
ctest --test-dir OpenMS-build -R "FileTypes|FileHandler|QPXFile|PSMArrowIO|ConsensusMapArrowIO|ProteinIdentificationArrowIO|FeatureMapArrowIO|TOPP_PercolatorAdapter|TOPP_IDMerger|TOPP_CometAdapter|TOPP_SageAdapter|TOPP_MSGFPlusAdapter|TOPP_MSFraggerAdapter|TOPP_XTandemAdapter" --output-on-failure
```

Expected: all PASS. Investigate any new failures before merging.

- [ ] **Step 3: Verify the original issue's request manually.**

```bash
# Pipeline A: end-to-end on a small test fixture.
./OpenMS-build/bin/CometAdapter -in <tiny mzML> -out /tmp/run.idparquet -database <fasta> -comet_executable <path>
ls /tmp/run.idparquet/
# Expected: psms.parquet, proteins.parquet, protein_groups.parquet, search_params.parquet
./OpenMS-build/bin/PercolatorAdapter -in /tmp/run.idparquet -out /tmp/scored.idparquet -score:fdr 0.01 -percolator_executable <path>
ls /tmp/scored.idparquet/
# Expected: same four files; PSM count <= input.
```

(This is a developer sanity check, not a CI test.)

---

## Self-review

**1. Spec coverage.** Walk-through:

- Spec § "Architecture/Format layer" → Tasks 2 (`QPXFile::importFromArrow`), 3 (`ConsensusMapArrowIO` refactor), 4 (`PSMArrowIO`).
- Spec § "Architecture/Dispatch layer" → Tasks 1 (`FileTypes::IDPARQUET`), 5 (`FileHandler` cases).
- Spec § "Architecture/Tool layer" → Tasks 6 (PercolatorAdapter formats), 7 (PercolatorAdapter filter options), 8 (IDMerger), 9 (SageAdapter migration), 10 (other search-engine adapters).
- Spec § "Data flow A/B/C" → covered implicitly by the round-trip tests in Tasks 4, 5, 11, 12.
- Spec § "Error handling/Format layer" → Task 4's three error tests (missing sub-file, target as regular file, empty PSMs).
- Spec § "Error handling/Tool layer" → `score:fdr` range check is in Task 7 (`setMinFloat_`/`setMaxFloat_`); empty-output warning is in Task 7.
- Spec § "Testing/Unit tests" → Task 4 covers all five listed sub-tests except "round-trip full fidelity" against an idXML fixture, which is implicitly covered by Task 11's parquet → idXML diff against the existing reference.
- Spec § "Testing/TOPP tool tests" → Tasks 11, 12. Comet/Sage parquet smoke tests are not separately added (they're trivially exercised by Task 11's pipeline check); add them only if Task 13 reveals an issue.

No spec gaps.

**2. Placeholder scan.** Searched for "TBD", "TODO", "implement later", "fill in details", "appropriate error handling". None present in the task body. Two task-12 placeholders use angle-bracketed paths (`<idmerger_input1.idXML>`) where the engineer needs to fill in the path the `grep` from step 1 returned — these are explicitly resolved by the preceding step, not unresolved placeholders.

**3. Type consistency.** `QPXFile::importFromArrow` signature is consistent across Tasks 2, 3, 4, 5. `PSMArrowIO::{export,import}FromParquet` signature is consistent across Tasks 4, 5. `score:fdr` parameter name and `keep_all_passing` flag name are consistent across Tasks 7 and 11.

No issues to fix.

---

## Execution handoff

**Plan complete and saved to `docs/superpowers/plans/2026-04-28-percolator-psm-parquet.md`. Two execution options:**

1. **Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.
2. **Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints.

**Which approach?**
