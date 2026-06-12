# OpenMS POLS Audit — Batch 3: FORMAT (I/O layer)

**Confirmed findings:** 150 · 14 high · 60 medium · 76 low.  
**Method:** 22 header-cluster finders → adversarial per-finding verification against source (retry-enabled). 1 verify dropped to API-overload.

> Post-verification severity/category/ABI. Recommendations favor ABI-safe fixes.

### [FORM-94] XICParquetFile::getChromatograms — Six adjacent same-typed Int64 filter parameters (all defaulting to -1) are trivially swappable at call sites
`high` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/XICParquetFile.h · _format-arrow-parquet_

```cpp
void getChromatograms(std::vector<XICChromatogram>& output, Int64 precursor_id = -1, Int64 transition_id = -1, const std::string& modified_sequence = "", Int64 precursor_charge = -1, Int64 product_charge = -1, Int64 ms_level = -1, Int64 run_id = -1, const std::string& filter = "") const
```
- **Expectation:** Filter arguments should be hard to confuse; a developer calling getChromatograms(out, 1318, 7) expects an obvious mapping of values to columns.
- **Actual:** getChromatograms has precursor_id, transition_id, precursor_charge, product_charge, ms_level, run_id as six positional Int64 params all defaulting to -1, interleaved with two string params. Passing e.g. a run_id where a transition_id is expected, or swapping precursor_id/transition_id, compiles cleanly and silently filters on the wrong column, returning a wrong (often empty) result set with no error.
- **Evidence:** XICParquetFile.h:198-206 declares the eight positional params; XIMParquetFile.h:207-218 is even longer (adds mobilogram_type/feature_id/feature_rt). All Int64 ID/charge/level params default to -1.
- **Fix:** Promote the typed ParquetFilter / ParquetFilterBuilder overloads (already present) as the primary API and deprecate or de-emphasize the positional overload; or introduce a small XICQuery struct of named fields. Additive (new overload/struct) keeps ABI; the existing overload can stay for compatibility.
- **Verified:** Verified directly in source. XICParquetFile.h:198-206 declares getChromatograms with six adjacent same-typed Int64 filter params (precursor_id, transition_id, precursor_charge, product_charge, ms_level, run_id), all defaulting to -1, interleaved with two std::string params (modified_sequence, filter). XICParquetFile.cpp confirms each positional Int64 maps to a DISTINCT column: buildFilterFromArgs_

### [FORM-129] MSChromatogramParquetConsumer::~MSChromatogramParquetConsumer / finalize — Destructor writes the file but swallows all write errors; only finalize() surfaces them
`high` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h · _format-dataaccess-misc_

```cpp
~MSChromatogramParquetConsumer() override; void finalize();
```
- **Expectation:** Constructing the consumer and letting it go out of scope (the RAII-natural usage) would be expected to either write the file successfully or make a failure visible.
- **Actual:** The destructor calls finalize() and catches every exception, logging to OPENMS_LOG_ERROR and continuing (cpp:191-212). If the user never calls finalize() explicitly, a failed/partial Parquet write is silently downgraded to a log line. The header does flag this ("Call this explicitly to surface write errors during normal control flow"), so the surprise is partially mitigated by docs.
- **Evidence:** cpp:191-212 destructor body catches BaseException/std::exception/... and only logs; header line 87-89: "Call this explicitly to surface write errors during normal control flow."
- **Fix:** Keep the swallow-in-destructor (throwing from a destructor is worse), but the header could state more prominently that relying on destructor-only finalization hides I/O errors. Documentation-only; no ABI impact.
- **Verifier correction:** The destructor genuinely swallows all Parquet write/close errors (cpp:191-212), and write_() performs terminal I/O including the metadata-writing Close() calls (cpp:506,519). The claim's recommendation ('Documentation-only; no ABI impact', severity ~medium) understates the issue: finalize() is neither part of the IMSDataConsumer interface nor virtual, so the OpenSwath workflow — which holds the consumer as IMSDataConsumer* (OpenSwathWorkflow.cpp:1329) and destroys it via 'delete chromatogramConsumer;' (line 1407) after prepareChromOutput new's it (OpenSwathBase.cpp:362) — has no way to invoke finalize() explicitly and is structurally forced onto the swallow path. The header note ('call finalize() explicitly to surface errors') does not help this in-tree caller. A real fix would add finalize()/error-state to the IMSDataConsumer interface (or have OpenSwathWorkflow dynamic_cast and finalize, as it already does a dynamic_cast for MSDataSqlConsumer at line 1367), not just documentation. Severity: high (silent loss/corruption of scientific .xic output for the main workflow). ABI of the current code: none.
- **Verified:** Verified against the actual code. The Impl destructor (cpp:191-212) calls finalize() and catches BaseException/std::exception/... logging only to OPENMS_LOG_ERROR — exactly as quoted. finalize()->write_() (cpp:495-528) does the real I/O: it opens the file, WriteTable (cpp:650), and critically calls parquet_writer_->Close() which writes Parquet metadata (cpp:506) and outfile_->Close() (cpp:519); an

### [FORM-9] Base64::decodeIntegersUncompressed_ (reached via decodeIntegers, zlib_compression=false) — Integer decode silently returns an empty/partial result for inputs shorter than 4 chars instead of signalling, and (unlike the float path) does not validate length % 4
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/Base64.h · _format-encoding_

```cpp
template <typename ToType> static void decodeIntegers(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out, bool zlib_compression = false)
```
- **Expectation:** Symmetric with decode(): malformed Base64 (length not a multiple of 4) should be reported the same way for integers as for floats.
- **Actual:** decodeUncompressed_ (float path) throws ConversionError when `in.size() % 4 != 0` (Base64.h lines 320-323). But decodeIntegersUncompressed_ (integer path) has NO such check — for `in.size() < 4` it just returns leaving `out` empty (lines 513-516) and otherwise proceeds, so a malformed integer buffer is decoded into garbage or a short result with no error. The two sibling decoders are inconsistent: one validates length, the other does not.
- **Evidence:** Float path Base64.h lines 320-323: `if (in.size() % 4 != 0) { throw Exception::ConversionError(...); }`. Integer path lines 513-516: `if (in.size() < 4) { return; }` with no `% 4` validation anywhere in decodeIntegersUncompressed_ (lines 507-645).
- **Fix:** Add the same `in.size() % 4 != 0` length check (and matching ConversionError) to decodeIntegersUncompressed_ so the integer and float decoders fail consistently. Source-compatible additive validation; only affects previously-silently-accepted malformed inputs.
- **Verifier correction:** The asymmetry is real (float decoder validates length % 4 and throws; integer decoder does not). However the consequence is worse than "garbage or short result with no error": because the loop reads in[i+2] and in[i+3] before its bounds guards, a malformed integer buffer whose length mod 4 is 1 or 2 causes an out-of-bounds read of `in` past size() (undefined behavior) and a secondary out-of-bounds index into the 80-element decoder_[] table. This is a memory-safety hazard on untrusted input, not just a wrong-value silent failure. Fix as recommended: add `if (in.size() % 4 != 0) throw Exception::ConversionError(...)` to decodeIntegersUncompressed_, matching decodeUncompressed_; this both restores symmetry and prevents the OOB read. The change is internal to a private template (signatures unchanged) so it is source-compatible, but it does alter runtime behavior for previously-silently-accepted malformed inputs (now they throw).
- **Verified:** Evidence confirmed verbatim. decodeUncompressed_ (float path) throws ConversionError when in.size() % 4 != 0 (Base64.h:320-323), but decodeIntegersUncompressed_ (integer path, Base64.h:507-645) has only an in.size() < 4 early-return (513-516) and NO % 4 validation. Both are private siblings of the symmetric public APIs decode()/decodeIntegers(), whose doc comments (lines 68-92) describe the same o

### [FORM-38] XMassFile::load — load() reads fid intensities with no EOF/short-file check; truncated fid yields zero-padded peaks instead of an error
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/XMassFile.h · _format-feature-fasta_

```cpp
void load(const std::string& filename, MSSpectrum& spectrum)
```
- **Expectation:** load() of a truncated or corrupt fid file (fewer intensity values than acqus declares) should signal a parse/IO error.
- **Actual:** The loop runs `while (spectrum.size() < acqus.getSize())` and pushes a peak for every index up to the acqus-declared size, calling FidHandler::getIntensity(). getIntensity() does `read(...,4)` with no success check and returns 0 on a failed/EOF read (`return (result > 0) ? result : 0`). So a fid shorter than acqus.getSize() silently produces a full-length spectrum padded with zero-intensity peaks; no exception is thrown. Only a completely unopenable fid triggers an exception (the early `if (!fid)` block).
- **Evidence:** XMassFile.h:85-91 loop `while (spectrum.size() < acqus.getSize()) { ... p.setIntensity(... fid.getIntensity()); spectrum.push_back(p); }`. FidHandler.cpp:47-57 getIntensity(): `read((char *) &result, 4); ... index_++; return (result > 0) ? result : 0;` (no stream-state check).
- **Fix:** Have getIntensity() check the stream state after read() and have load() stop / throw Exception::ParseError when the fid stream reaches EOF before acqus.getSize() peaks are read, rather than zero-padding. The fix is internal to the .cpp/handler and ABI-neutral at the XMassFile API level.
- **Verified:** Independently verified against the actual code. XMassFile.h:85-91: load() loops `while (spectrum.size() < acqus.getSize())`, pushing one peak per index up to the acqus-declared count, with intensity from FidHandler::getIntensity(). AcqusHandler::getSize() (AcqusHandler.cpp:82-85) returns td_, parsed from the $TD field of the separate `acqus` file (line 72) — entirely independent of the actual fid 

### [FORM-113] GNPSQuantificationFile::store — store() indexes consensus_map[0] unconditionally, throwing/UB on an empty ConsensusMap
`high` · `surprising-throw` · ABI: `none` · src/openms/source/FORMAT/GNPSQuantificationFile.cpp · _format-gnps-deconv-qc_

```cpp
static void store(const ConsensusMap& consensus_map, const std::string& output_file)
```
- **Expectation:** Storing an empty ConsensusMap writes a header-only (or empty) table, not a crash.
- **Actual:** The first statement dereferences consensus_map[0] to probe for the IIMN_ROW_ID meta value before any size check, so calling store() with an empty map indexes out of bounds.
- **Evidence:** GNPSQuantificationFile.cpp:25-26 'bool iimn = false; if (consensus_map[0].metaValueExists(Constants::UserParam::IIMN_ROW_ID)) iimn = true;' with no 'if (consensus_map.empty())' guard anywhere above.
- **Fix:** Guard with 'if (!consensus_map.empty() && consensus_map[0].metaValueExists(...))'. Additive, no signature change.
- **Verifier correction:** store() dereferences consensus_map[0] before any size check to probe for IIMN_ROW_ID. ConsensusMap::operator[] is noexcept and forwards to std::vector::operator[], so on an empty map this is undefined behavior (out-of-bounds read), not a thrown exception. Fix: guard with `if (!consensus_map.empty() && consensus_map[0].metaValueExists(...))` so an empty map writes the header-only table as documented.
- **Verified:** Confirmed by reading the source. GNPSQuantificationFile::store (src/openms/source/FORMAT/GNPSQuantificationFile.cpp:25-26) unconditionally evaluates `consensus_map[0].metaValueExists(Constants::UserParam::IIMN_ROW_ID)` as its first action, before any size check; a full read of the function shows no `if (consensus_map.empty())` guard anywhere. ConsensusMap::operator[] (via EXPOSED_VECTOR_INTERFACE 

### [FORM-60] InspectOutfile::load — load() appends to the output identification containers without clearing them
`high` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/FORMAT/InspectOutfile.h · _format-id-search-out_

```cpp
std::vector<Size> load(const std::string & result_filename, PeptideIdentificationList & peptide_identifications, ProteinIdentification & protein_identification, const double p_value_threshold, const std::string & database_filename = "")
```
- **Expectation:** A load() that produces peptide_identifications / protein_identification is expected to populate fresh results; reusing the same objects across two load() calls should not merge two files' hits.
- **Actual:** load() only push_back/insertHit onto the passed containers; it never clears them. InspectOutfile.cpp:183 'protein_identification.insertHit(protein_hit);', :213 and :269 'peptide_identifications.push_back(peptide_identification);'. A second load() into the same containers accumulates hits from both files. This is the same inconsistency as OMSSACSVFile but differs from PercolatorOutfile which clears.
- **Evidence:** InspectOutfile.cpp:183  protein_identification.insertHit(protein_hit);\nInspectOutfile.cpp:213  peptide_identifications.push_back(peptide_identification);\nInspectOutfile.cpp:269  peptide_identifications.push_back(peptide_identification);\n(no clear of the output containers)
- **Fix:** Document the append behavior or clear the output containers at entry to match the clearing loaders. If callers rely on append, at minimum state it in the header. Clearing is a source-compatible behavior change; abi none.
- **Verified:** Independently verified against the source. InspectOutfile::load (src/openms/source/FORMAT/InspectOutfile.cpp:55) appends to its output containers and never clears them: protein_identification.insertHit at :183, peptide_identifications.push_back at :213 and :269 — exactly as quoted. The only clear() calls in the function (:97, :141, :263) act on the result_file ifstream, and :293-296 clear local te

### [FORM-42] MzIdentMLFile::load — MzIdentMLFile::load APPENDS to the output containers, while sibling *File::load methods CLEAR them first
`high` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzIdentMLFile.h · _format-id-xml_

```cpp
void load(const std::string& filename, std::vector<ProteinIdentification>& poid, PeptideIdentificationList& peid)
```
- **Expectation:** A competent caller expects load() to populate the output vectors with the file's contents, i.e. to clear/replace whatever was passed in. This is reinforced by the class doc: 'The file adapter interface is kept the same as idXML file adapter for downward capability reasons.' IdXMLFile::load (the adapter it claims to mirror) explicitly clears: 'protein_ids.clear(); peptide_ids.clear();'. PepXMLFile::load and PepXMLFileMascot::load and PTMXMLFile::load all clear their outputs too.
- **Actual:** MzIdentMLFile::load does NOT clear poid/peid. It forwards them to MzIdentMLDOMHandler, which only ever appends via pro_id_->emplace_back() / pep_id_->push_back() (e.g. MzIdentMLDOMHandler.cpp lines 884, 1368, 1811, 2074). readMzIdentMLFile contains no .clear(). Re-using the same vectors (or passing a non-empty vector) silently accumulates duplicate identifications.
- **Evidence:** MzIdentMLFile.cpp: 'void MzIdentMLFile::load(...) { Internal::MzIdentMLDOMHandler handler(poid, peid, schema_version_, *this); handler.readMzIdentMLFile(filename); }' -- no clear. Contrast IdXMLFile.cpp: 'protein_ids.clear(); peptide_ids.clear();'. MzIdentMLDOMHandler.cpp: 'pro_id_->emplace_back();' (l.884), 'pep_id_->push_back(PeptideIdentification());' (l.1368).
- **Fix:** Additive, ABI-safe fix: clear poid and peid at the top of MzIdentMLFile::load before constructing the handler (matching IdXMLFile/PepXMLFile). This is a behavior change but the correct/consistent one; document it. No signature/ABI change.
- **Verified:** Independently verified in source. MzIdentMLFile::load (src/openms/source/FORMAT/MzIdentMLFile.cpp:32-36) does NOT clear poid/peid: it constructs MzIdentMLDOMHandler (constructor at MzIdentMLDOMHandler.cpp:87-92 only binds pointers, no clear) and calls readMzIdentMLFile (l.163), which performs no .clear() on pro_id_/pep_id_ and only ever appends: pro_id_->emplace_back() at l.884, pep_id_->push_back

### [FORM-22] IndexedMzMLFileLoader::load — load() returns bool success that is easy to ignore; failure produces no exception
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/IndexedMzMLFileLoader.h · _format-mzml_

```cpp
bool load(const std::string& filename, OnDiscPeakMap& exp)
```
- **Expectation:** Sibling *File classes (MzMLFile, MzDataFile, MzXMLFile) all declare load() as 'void' and signal a missing/garbage file by throwing FileNotFound/ParseError. A caller used to those classes would write `loader.load(fn, exp);` and assume an exception on failure.
- **Actual:** load() returns a bool (forwarded from exp.openFile(filename)) and throws nothing on failure. The doc even says 'Tries to parse the file, success needs to be checked with the return value.' A caller who ignores the return value silently proceeds with an empty/invalid OnDiscPeakMap.
- **Evidence:** Header: `bool load(const std::string& filename, OnDiscPeakMap& exp);` with doc '@return Indicates whether parsing was successful (if it is false, the file most likely was not an mzML or not indexed).' Impl: `bool IndexedMzMLFileLoader::load(...) { return exp.openFile(filename); }`
- **Fix:** Keep the bool (ABI), but document the divergence from the other *File classes prominently and mark the return value [[nodiscard]] (source-compatible additive change) so callers cannot silently drop the failure indication.
- **Verifier correction:** Severity confirmed high: a caller migrating from the sibling void-load classes writes `loader.load(fn, exp)` expecting an exception on a missing/garbage file, but instead receives `false` (commonly ignored) plus an empty OnDiscPeakMap and silently proceeds with no data — silently wrong results / data loss for reasonable use, evidenced by the in-tree TICCalculator caller that ignores the return value. Recommendation (mark the bool [[nodiscard]] and document the divergence) is correct and source-compatible: [[nodiscard]] is an attribute that only emits a warning at discarding call sites and does not change the signature or ABI.
- **Verified:** All quoted evidence verified verbatim. Header line 57: `bool load(const std::string& filename, OnDiscPeakMap& exp);` with doc "success needs to be checked with the return value" (line 50) and "@return Indicates whether parsing was successful (if it is false, the file most likely was not an mzML or not indexed)" (line 55). Impl IndexedMzMLFileLoader.cpp:39: `return exp.openFile(filename);`. Sibling

### [FORM-66] MzTabBoolean::setNull — MzTabBoolean::setNull has inverted polarity vs every sibling setNull
`high` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/MzTabBase.h · _format-mztab_

```cpp
void MzTabBoolean::setNull(bool b)
```
- **Expectation:** setNull(true) marks the cell as null; setNull(false) clears the null state. This is the contract of MzTabString/MzTabInteger/MzTabDouble/MzTabModification::setNull in the same module.
- **Actual:** Implementation is `if (!b) value_ = -1; else value_ = 0;`. So setNull(false) MAKES the value null (-1) and setNull(true) makes it NOT null (and sets it to boolean false/0). The boolean argument is inverted relative to its name and to all sibling classes.
- **Evidence:** MzTabBase.cpp:503-509: `void MzTabBoolean::setNull(bool b){ if (!b) value_ = -1; else value_ = 0; }` vs MzTabInteger.cpp:691-694 `state_ = b ? MZTAB_CELLSTATE_NULL : MZTAB_CELLSTATE_DEFAULT;` and MzTabDouble:721-724 identical, and MzTabString:448-454 `if (b) value_.clear();`.
- **Fix:** Fix the body to `value_ = b ? -1 : 0;` so it matches the sibling convention. This is a source-compatible behavior fix (the signature is unchanged); callers currently calling setNull(true) to nullify are themselves buggy. Guard with a test. If existing callers depend on the inverted behavior, add a deprecation note but correct the polarity.
- **Verified:** Verified directly in source. MzTabBase.cpp:503-509 reads `void MzTabBoolean::setNull(bool b){ if (!b) value_ = -1; else value_ = 0; }`, while isNull() (line 498-501) is `value_ < 0`. So setNull(false) sets value_=-1 → isNull()==true (nullifies), and setNull(true) sets value_=0 → isNull()==false (NOT null, literal boolean false). This is inverted relative to its name and to every sibling: MzTabInte

### [FORM-69] MzTabFile::load — load/store round-trip silently drops nucleic-acid, oligonucleotide and OSM sections
`high` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzTabFile.h · _format-mztab_

```cpp
void load(const std::string& filename, MzTab& mz_tab)
```
- **Expectation:** store() writes NUC/OLI/OSM sections (the header and store() fully support MzTabNucleicAcid/Oligonucleotide/OSM rows), so a load() of a file written by store() should restore them - a symmetric round-trip.
- **Actual:** load() only assigns metadata, protein, peptide, PSM and small-molecule sections back into the MzTab object. It never calls setNucleicAcidSectionRows / setOligonucleotideSectionRows / setOSMSectionRows, so those sections are silently discarded on read even though store() emits them. store()->load() is lossy with no error.
- **Evidence:** MzTabFile.cpp:1548-1554 assigns only setMetaData/setProteinSectionRows/setPeptideSectionRows/setPSMSectionRows/setSmallMoleculeSectionRows/setEmptyRows/setCommentRows. store() at MzTabFile.cpp:3222-3259 writes nucleic_acid_section/oligonucleotide_section/osm_section. No NUC/OLI/OSM parsing exists in load.
- **Fix:** Either implement parsing for the NUC/OLI/OSM sections in load(), or at minimum emit a warning/throw when such sections are encountered so the data loss is not silent. Additive to load() impl; no ABI change.
- **Verified:** Independently verified in source. MzTabFile::load (MzTabFile.cpp:96-1555) dispatches on a 3-char section prefix through an if-chain handling only COM/MTD/PRH/PRT/PEH/PEP/PSH/PSM/SMH/SML. NUH/NUC, OLH/OLI, OSH/OSM lines match no branch and fall through the loop body, silently dropped (not even captured as comment/empty rows). The final commit at lines 1548-1554 calls setMetaData/setProteinSectionRo

### [FORM-33] MSPFile::load(const std::string&, PeptideIdentificationList&, PeakMap&) — load() resets the output PeakMap but appends to the output ids list without clearing it
`high` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MSPFile.h · _format-peakfiles_

```cpp
void load(const std::string & filename, PeptideIdentificationList & ids, PeakMap & exp)
```
- **Expectation:** A load() that fills two output containers should treat them consistently — both cleared before being populated, so reusing the same objects across two load() calls yields only the second file's content.
- **Actual:** exp is cleared via `exp.reset()` (MSPFile.cpp:72), but ids is never cleared; the loop only does `ids.push_back(id)` (MSPFile.cpp:119). Reusing the same PeptideIdentificationList across two load() calls silently accumulates identifications from both files, desynchronizing ids from exp.
- **Evidence:** MSPFile.cpp:72 `exp.reset();` vs MSPFile.cpp:119 `ids.push_back(id);` with no `ids.clear()` anywhere (grep for ids.clear returns nothing).
- **Fix:** Add `ids.clear();` next to `exp.reset();` at the top of load(). Behavior-fix in the .cpp, no signature/ABI change. If preserving append semantics is intentional, document it explicitly in the header.
- **Verifier correction:** No correction needed; claim is accurate. Severity confirmed high because object reuse across files is a reasonable usage pattern and the failure mode is silent index-level desynchronization of ids vs exp (no throw, no log), i.e. silently-wrong results. The recommended fix (add ids.clear() next to exp.reset() at the top of load()) is correct and a behavior-only .cpp change with no ABI impact. Note: the AnnotatedMSRun overload (line 327) is immune because it builds fresh local containers each call, but the public PeptideIdentificationList/PeakMap overload is directly exposed and affected.
- **Verified:** Independently verified against the actual code. MSPFile::load(const std::string&, PeptideIdentificationList&, PeakMap&) is defined in src/openms/source/FORMAT/MSPFile.cpp lines 54-325. Line 72 calls exp.reset(); MSExperiment::reset() (MSExperiment.cpp:858-864) fully wipes spectra_, chromatograms_, ranges, and meta info. By contrast ids is only appended via ids.push_back(id) (line 119) inside the p

### [FORM-12] CsvFile (class doc) vs CsvFile::load / constructor — CsvFile doc claims no comment support, but '#'-prefixed lines are silently dropped
`high` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/FORMAT/CsvFile.h · _format-text-streams_

```cpp
@brief This class handles csv files. ... Does NOT support comment lines!
```
- **Expectation:** Given the class documentation 'Does NOT support comment lines!', a caller expects every line of the CSV to be loaded verbatim, including any line that happens to start with '#'.
- **Actual:** Both the constructor and load() pass `"#"` as the comment_symbol to TextFile::load, so every line beginning with '#' is silently skipped and never appears in the buffer or in rowCount(). A field-less CSV whose first column legitimately starts with '#' loses rows with no error.
- **Evidence:** CsvFile.cpp: `TextFile::load(filename, false, first_n, false, "#");` (line 27) and `TextFile::load(filename, true, first_n, false, "#");` (line 34); class doc in CsvFile.h line 17: 'Does NOT support comment lines!'
- **Fix:** Either honor the documentation by passing an empty comment_symbol, or update the documentation to state that '#'-prefixed lines are treated as comments and skipped. Prefer making the comment behavior an explicit, documented parameter. ABI-safe (constant/string change or additive param).
- **Verifier correction:** The category is not "misleading-name" (no symbol is misnamed); it is a documentation/contract contradiction. CsvFile's class doc states "Does NOT support comment lines!" while its constructor and load() unconditionally pass "#" as the comment_symbol to TextFile::load (CsvFile.cpp:27,34), and TextFile::load (TextFile.cpp:66-69) silently skips any line prefixed with "#" — not added to buffer_, not counted in rowCount(). Result: rows whose first field legitimately begins with '#' are silently dropped with no error or diagnostic. Minimal ABI-safe fix: either pass an empty comment_symbol to honor the doc, or correct the doc to state that '#'-prefixed lines are treated as comments and skipped (both are constant/doc-only changes, abi_impact none). A fuller fix exposing a defaulted comment_symbol parameter on CsvFile::load/ctor would be source-compatible/additive.
- **Verified:** Verified line-by-line. CsvFile.h:17 documents "@brief This class handles csv files. ... Does NOT support comment lines!" yet CsvFile.cpp:27 (constructor) and CsvFile.cpp:34 (load) both pass "#" as TextFile::load's comment_symbol argument. TextFile.cpp:66-69 confirms: when comment_symbol is non-empty and StringUtils::hasPrefix(str, comment_symbol) matches, the line is `continue`d — never pushed to 

### [FORM-108] RationalScan2ImConverter::convert / getCalibration — Unknown/out-of-range frame_id silently falls back to the FIRST calibration instead of failing
`high` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/RationalScan2ImConverter.h · _format-vendor_

```cpp
void convert(uint32_t frame_id, double* inv_ion_mobilities, const double* scans, uint32_t size) override
```
- **Expectation:** Converting with a frame_id that has no calibration (out of range, or absent from the map) would be a hard error, since using a wrong frame's calibration yields silently wrong ion-mobility values.
- **Actual:** getCalibration() returns calibrations_.begin()->second (the first/arbitrary calibration) for any frame_id not found, logging only a WARN. convert() then produces plausible-looking but physically wrong 1/K0 values. The header documents Coefficients/strategy but the public convert() override gives no hint that it tolerates bad frame_ids. Additionally, if the converter was constructed (its ctor is public) with an empty calibrations_ map, calibrations_.begin() dereferences end() — undefined behavior.
- **Evidence:** RationalScan2ImConverter.cpp:38-44 `// Fallback: use first calibration entry ... return calibrations_.begin()->second;` reached whenever frame_id is out of range or frame_to_cal_[frame_id] is not in calibrations_; header :58-60 public ctor accepts arbitrary maps.
- **Fix:** Throw on unknown frame_id (or at least guard calibrations_.begin() against an empty map and document the silent-fallback policy). If fallback must stay for robustness, state it explicitly in the convert() doc.
- **Verified:** Verified against the actual code. getCalibration() (RationalScan2ImConverter.cpp:26-45) returns calibrations_.begin()->second for ANY out-of-range frame_id or any frame_to_cal_[frame_id] not present in calibrations_, logging only OPENMS_LOG_WARN (and not even that when frame_id==0). convert() (cpp:69-77) then runs applyFormula() with this wrong frame's coefficients, yielding plausible but physical

### [FORM-109] BrukerTimsImagingFile::load — load() silently drops MALDI pixels whose frame id is absent from the loaded spectra (warn-only)
`high` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/BrukerTimsImagingFile.h · _format-vendor_

```cpp
void load(const std::string& path, MSImagingExperiment& exp, const Config& config)
```
- **Expectation:** A load() that builds an imaging geometry from MaldiFrameInfo would either map every declared pixel or fail when pixels cannot be matched to spectra, so the resulting image is complete.
- **Actual:** Pixels in MaldiFrameInfo whose frame_id has no corresponding loaded spectrum are skipped and only counted into a one-shot WARN; load() still returns 'successfully' with an image that is missing those pixels. With a frame_id_min/frame_id_max filter set in inner_config this is the common case, and the geometry silently omits the out-of-range pixels. The header's @throws list does not mention any partial-result condition.
- **Evidence:** BrukerTimsImagingFile.cpp:226-236 `if (it == frame_to_index.end()) { ++dropped; continue; } ... if (dropped > 0) { OPENMS_LOG_WARN << "... dropped " << dropped << " pixel(s) ..."; }`
- **Fix:** Document the warn-and-drop partial-image behavior in load()'s header contract (it currently lists only throwing failure modes), and consider returning the dropped count or offering a strict mode that throws when any declared pixel is unmatched.
- **Verified:** Verified against the actual code. BrukerTimsImagingFile.cpp:223-236 contains exactly the quoted warn-and-drop logic: rows from MaldiFrameInfo whose frame_id is absent from frame_to_index are counted into `dropped` and `continue`d, emitting only a one-shot OPENMS_LOG_WARN, then load() finishes normally. The mechanism the claim describes is real and reachable in normal use: readMaldiFrameInfo (SQL a

### [FORM-93] ConsensusMapArrowIO::importFeaturesFromArrow / importPSMsFromArrow — importFeaturesFromArrow APPENDS to the ConsensusMap instead of clearing it (unlike importFromParquet)
`medium` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/ConsensusMapArrowIO.h · _format-arrow-parquet_

```cpp
static bool importFeaturesFromArrow(const std::shared_ptr<arrow::Table>& table, ConsensusMap& cmap)
```
- **Expectation:** A method documented as '@param[out] cmap ConsensusMap to populate' is expected to produce a ConsensusMap that contains exactly the table's rows, i.e. to replace/clear prior contents (matching the sibling importFromParquet, which does `cmap = ConsensusMap{}`).
- **Actual:** importFeaturesFromArrow does `cmap.push_back(std::move(cf))` per row without ever clearing cmap, so calling it on a non-empty map silently concatenates. The directory-level importFromParquet (same class) instead resets via `cmap = ConsensusMap{};`. The two import entry points have opposite clear/append semantics.
- **Evidence:** ConsensusMapArrowIO.cpp:1233 `cmap.push_back(std::move(cf));` (no clear in importFeaturesFromArrow body, lines 1159-1237) vs ConsensusMapArrowIO.cpp:1357 `cmap = ConsensusMap{};` in importFromParquet. Header marks param as `@param[out] cmap ConsensusMap to populate` with no mention of append.
- **Fix:** Either clear cmap at the top of importFeaturesFromArrow (and importPSMsFromArrow operates on the just-cleared map) to match importFromParquet, or document the append semantics explicitly in the header. Behavior change is source-compatible (no signature change); if append is intentional for chunked reads, at minimum amend the Doxygen to say 'appends to cmap'.
- **Verifier correction:** importFeaturesFromArrow (public) appends via push_back without clearing cmap, while the sibling public importFromParquet resets via `cmap = ConsensusMap{}` — a real undocumented clear/append asymmetry, and the `@param[out] cmap` doc contradicts the append behavior. Correction to the claim: importPSMsFromArrow is NOT a candidate for clearing — it depends on features already present in cmap (builds feature_lookup at cpp:1302-1306 to route PSMs), so it is genuinely [in,out]; the recommendation should target only importFeaturesFromArrow (clear at top) or amend the Doxygen of both helpers to reflect their true in/out semantics. Severity is medium, not high: the primary documented entry point importFromParquet is correct, the class is @experimental, and the silent concatenation requires directly reusing a non-empty map with the low-level helper.
- **Verified:** Evidence verified exactly against source. importFeaturesFromArrow (public, ConsensusMapArrowIO.h:101, doc'd `@param[out] cmap ConsensusMap to populate`) loops `cmap.push_back(std::move(cf))` (ConsensusMapArrowIO.cpp:1207-1234) with NO clear at the top, so calling it on a non-empty map silently concatenates features with no error/warning. The sibling public entry point importFromParquet does `cmap 

### [FORM-95] ProteinGroupArrowExport::exportToArrow — Two overloads of exportToArrow return opposite things for the 'no groups' case (nullptr vs empty-non-null table)
`medium` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/ProteinGroupArrowExport.h · _format-arrow-parquet_

```cpp
static std::shared_ptr<arrow::Table> exportToArrow(const ConsensusMap& cmap)  vs  static std::shared_ptr<arrow::Table> exportToArrow(const std::vector<ProteinIdentification>&, const PeptideIdentificationList&)
```
- **Expectation:** Two overloads of the same function name should agree on what 'no protein groups' yields, so a caller can uniformly do `auto t = exportToArrow(...); writeTable(t);` without a null check that depends on which overload was picked.
- **Actual:** The ConsensusMap overload returns nullptr when there are no protein identifications or no indistinguishable groups (an empty-input, not an error). The std::vector<ProteinIdentification> overload is explicitly documented to return 'an empty table if no groups, never nullptr'. So identical 'no groups' input yields a null pointer on one overload and a valid empty table on the other.
- **Evidence:** ProteinGroupArrowExport.cpp:37 `return nullptr;` and :46 `return nullptr;` (after OPENMS_LOG_WARN 'No indistinguishable protein groups found'); header line 86 documents the other overload: 'Shared pointer to Arrow Table (empty table if no groups, never nullptr)'.
- **Fix:** Make the ConsensusMap overload also return an empty-but-valid table for the no-groups case (matching the ident overload), reserving nullptr for genuine build errors; or document the divergence prominently. Source-compatible (return type unchanged).
- **Verified:** Verified against the actual source. Two same-name overloads of ProteinGroupArrowExport::exportToArrow disagree on what the identical "no protein groups" (a benign empty input, not an error) yields. ConsensusMap overload (ProteinGroupArrowExport.cpp): line 37 `return nullptr;` when getProteinIdentifications() is empty, and line 46 `return nullptr;` after OPENMS_LOG_WARN "No indistinguishable protei

### [FORM-96] ParquetFile::getInt64 / getDouble / getBool vs getString / getStringList — String value accessors silently return empty on null while the numeric accessors in the same family force the caller to opt into null handling
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/ParquetFile.h · _format-arrow-parquet_

```cpp
static int64_t getInt64(const std::shared_ptr<arrow::Array>&, int64_t row, int64_t default_value, bool allow_null) ; static std::string getString(const std::shared_ptr<arrow::Array>&, int64_t row)
```
- **Expectation:** A coherent type-safe accessor family should treat null uniformly; if getInt64 makes the caller pass allow_null and a default to avoid a throw, getString should too (or at least signal null).
- **Actual:** getInt64/getDouble/getBool require explicit (default_value, allow_null) and throw MissingInformation on null when allow_null is false. getString/getStringList have no such params and silently return "" / empty vector for null OR out-of-range rows, so a missing string is indistinguishable from a present empty string and never raises.
- **Evidence:** ParquetFile.h:178-215 (numeric accessors with default_value+allow_null, '@throws MissingInformation if value is null and allow_null is false') vs ParquetFile.h:229/243 (getString/getStringList, no null params); ParquetFile.cpp:317 `return "";` and :319 `if (array->IsNull(row)) return "";`.
- **Fix:** Document the divergence explicitly on getString/getStringList and/or add an allow_null/has-value variant for strings to bring the family in line. Additive overloads keep ABI.
- **Verifier correction:** getInt64/getDouble/getBool force the caller to pass (default_value, allow_null) and throw MissingInformation on null when allow_null=false, whereas getString/getStringList provide no null-handling parameter and document-and-deliver a silent null->""/empty result, with no overload to mark a string as required. The divergence is the absence of any opt-in to throw on a missing string (so null and present-empty are indistinguishable), NOT a difference in out-of-range handling — both families index via array->IsNull(row) without bounds checks and are equally unsafe on out-of-range rows.
- **Verified:** Verified against ParquetFile.h and ParquetFile.cpp. The core asymmetry is real and provable: getInt64/getDouble/getBool (h:178-215) require (default_value, allow_null) and throw Exception::MissingInformation on null when allow_null is false (cpp:208-213, 254-258, 288-292), while getString/getStringList (h:229/243) take no null parameters and silently coerce null to ""/empty (cpp:319, 336). There i

### [FORM-119] MSDataSqlConsumer::consumeSpectrum / consumeChromatogram — consumeSpectrum/consumeChromatogram silently wipe the caller's spectrum/chromatogram data (undocumented)
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h · _format-dataaccess-consumers_

```cpp
void consumeSpectrum(SpectrumType & s) override; void consumeChromatogram(ChromatogramType & c) override;
```
- **Expectation:** The header doc reads only "@brief Write a spectrum to the output file" / "@brief Write a chromatogram to the output file". A caller passing a non-const reference for writing would not expect its own object to be emptied as a side effect; the sibling class MSDataCachedConsumer explicitly documents "@note May delete data from spectrum (if clearData is set)" for exactly this behavior.
- **Actual:** The implementation unconditionally clears the input after buffering it: consumeSpectrum does `spectra_.push_back(s); s.clear(false);` and consumeChromatogram does `chromatograms_.push_back(c); c.clear(false);`. There is no flag to disable this and no mention of it in the header.
- **Evidence:** MSDataSqlConsumer.cpp:77-80 `void MSDataSqlConsumer::consumeSpectrum(SpectrumType & s) { spectra_.push_back(s); s.clear(false);` and :91-94 `void MSDataSqlConsumer::consumeChromatogram(ChromatogramType & c) { chromatograms_.push_back(c); c.clear(false);`
- **Fix:** Document the mutation in the header (mirror MSDataCachedConsumer's `@note May delete data...`), e.g. add `@param[in,out] s ... @note Clears the data arrays of the passed spectrum after buffering`. Pure doc fix; ABI-neutral. Ideal (breaking) fix would be to operate on a copy like MSDataWritingConsumer, but that changes observable behavior.
- **Verifier correction:** consumeSpectrum/consumeChromatogram unconditionally clear the passed object's data arrays (via s.clear(false)/c.clear(false), keeping meta-data) after buffering a copy, with no flag to disable and no mention in the header — unlike the sibling MSDataCachedConsumer which both documents this ("@note May delete data...") and gates it behind a clearData constructor flag. Severity is medium rather than high: the consumer's own output is correct (it copies before clearing) and the dominant fire-and-forget streaming idiom is unaffected; silent peak-data loss only occurs if a caller reuses the spectrum/chromatogram after consuming it. Fix is a pure header-doc addition (mirror MSDataCachedConsumer's @note), ABI-neutral.
- **Verified:** Independently verified. Header (MSDataSqlConsumer.h:77-85) documents only "@brief Write a spectrum/chromatogram to the output file" with no mention of mutation and no disabling flag. Implementation (MSDataSqlConsumer.cpp:77-103) confirms the exact quoted evidence: consumeSpectrum does `spectra_.push_back(s); s.clear(false);` and consumeChromatogram does `chromatograms_.push_back(c); c.clear(false)

### [FORM-123] SiriusFragmentAnnotation::extractAnnotationsFromSiriusFile — Doc-comment describes the `decoy` flag's true/false branches exactly backwards from what the code does
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h · _format-dataaccess-misc_

```cpp
static std::vector<MSSpectrum> extractAnnotationsFromSiriusFile(const std::string& path_to_sirius_workspace, Size max_rank = 1, bool decoy = false, bool use_exact_mass = false)
```
- **Expectation:** The Doxygen block (lines 56-57) states: "If @p decoy is true, uses fragment annotation (./spectra/1_sumformula.tsv) from SIRIUS output ... else uses fragment annotation (./decoy/1_sumformula.tsv) from SIRIUS/PASSATUTTO output." A caller reading this would pass decoy=true to read TARGET spectra from ./spectra and decoy=false to read DECOYS from ./decoy.
- **Actual:** The implementation does the opposite: `std::string subfolder = decoy ? "/decoys/" : "/spectra/";` (SiriusFragmentAnnotation.cpp:268). decoy=true reads ./decoys/ (decoys), decoy=false reads ./spectra/ (targets). The documented mapping of the boolean to the folder is inverted, and the folder name in the doc ("./decoy/") doesn't even match the real folder ("/decoys/").
- **Evidence:** Header lines 56-57: "If @p decoy is true, uses fragment annotation (./spectra/1_sumformula.tsv) from SIRIUS output (per compound) else uses fragment annotation (./decoy/1_sumformula.tsv) from SIRIUS/PASSATUTTO output (per compound)." vs cpp:268 `std::string subfolder = decoy ? "/decoys/" : "/spectra/";`
- **Fix:** Fix the doc-comment to match behavior: decoy=true -> ./decoys/ (PASSATUTTO decoys), decoy=false -> ./spectra/ (targets). Documentation-only change; no ABI impact.
- **Verifier correction:** The Doxygen prose at lines 56-57 inverts the decoy flag's branches and misnames the folder. Correct mapping: decoy==true -> reads ./decoys/ (PASSATUTTO decoy annotations); decoy==false -> reads ./spectra/ (SIRIUS target annotations). This contradicts the function's own @param doc at line 79 ("Extract annotations for decoys? Or else targets."), which is correct. Documentation-only fix; no behavior or ABI change.
- **Verified:** Verified directly. Header lines 56-57 state "If @p decoy is true, uses fragment annotation (./spectra/1_sumformula.tsv) from SIRIUS output ... else uses fragment annotation (./decoy/1_sumformula.tsv) from SIRIUS/PASSATUTTO output." The implementation at SiriusFragmentAnnotation.cpp:268 is `std::string subfolder = decoy ? "/decoys/" : "/spectra/";` — the exact inverse: decoy==true reads /decoys/ (P

### [FORM-128] MSChromatogramParquetConsumer (class) — Chromatogram data is always written lossy (MSNumpress); there is no lossless option exposed
`medium` · `surprising-default` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h · _format-dataaccess-misc_

```cpp
class MSChromatogramParquetConsumer : public Interfaces::IMSDataConsumer ... lossy compression hardcoded
```
- **Expectation:** The header documents three RT compression schemes (0 none, 1 zlib, 5 MSNumpress lossy) and three intensity schemes including lossless ones, implying the writer can produce lossless output. A user wanting exact RT/intensity round-trip would expect a way to select scheme 1 (zlib raw doubles).
- **Actual:** The implementation hardcodes `const bool use_lossy_compression = true;` (cpp:432), so RT is always written as scheme 5 (MSNumpress linear, lossy, 0.05 Da fixed mass accuracy) and intensity as scheme 6 (slof, lossy). The non-lossy code path (cpp:468-475) is dead code; no constructor parameter or setter selects it.
- **Evidence:** cpp:432-434: `const bool use_lossy_compression = true; const int64_t rt_compression = use_lossy_compression ? 5 : 1; const int64_t intensity_compression = use_lossy_compression ? 6 : 1;` Header lists schemes 0/1 as if available.
- **Fix:** Either expose a compression-mode constructor argument (additive, ABI-safe via new overload) or document that output is always lossy MSNumpress and the schemes 0/1 are not produced by this writer. Prefer documenting now; additive overload is the ideal fix.
- **Verified:** Independently verified against the actual source. The class MSChromatogramParquetConsumer (header lines 62-107) exposes only a constructor (filename, run_id, source_file, transition_exp), destructor, consume*, finalize, setExpectedSize, setExperimentalSettings — no compression-mode parameter or setter anywhere in the public API; the constructor (cpp:699-705) just forwards to the impl. In encodeChr

### [FORM-130] MobilogramParquetConsumer::consumeMobilogram — Six defaulted, mostly same-typed (Int64/string) parameters are easily transposed at the call site
`medium` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h · _format-dataaccess-misc_

```cpp
void consumeMobilogram(const Mobilogram& m, const std::string& mobilogram_type = "", Int64 ms_level = -1, Int64 transition_id = -1, const std::string& transition_native_id = "", double feature_rt = quiet_NaN(), Int64 feature_id = -1)
```
- **Expectation:** With all-defaulted trailing parameters of overlapping types, a caller expects positional safety; e.g. they might pass transition_id where ms_level is expected, or swap the two string parameters (mobilogram_type vs transition_native_id).
- **Actual:** The signature has, in order: string mobilogram_type, Int64 ms_level, Int64 transition_id, string transition_native_id, double feature_rt, Int64 feature_id — three Int64s and two strings, all defaulted. Transposing transition_id and ms_level, or the two string args, compiles silently and writes wrong columns. Sentinels (-1, "", NaN) mean a wrong value is also indistinguishable from "not set".
- **Evidence:** Header lines 57-63: `void consumeMobilogram(const Mobilogram& m, const std::string& mobilogram_type = "", Int64 ms_level = -1, Int64 transition_id = -1, const std::string& transition_native_id = "", double feature_rt = std::numeric_limits<double>::quiet_NaN(), Int64 feature_id = -1);`
- **Fix:** Consider an options/struct parameter (e.g. MobilogramRecordInfo) to make call sites self-documenting; additive overload keeps ABI. At minimum the existing per-param docs help. Mark as source-compatible if an overload is added.
- **Verified:** Signature verified exactly as quoted (header lines 57-63; mirrored in .cpp 662-668): consumeMobilogram(const Mobilogram&, string mobilogram_type="", Int64 ms_level=-1, Int64 transition_id=-1, string transition_native_id="", double feature_rt=NaN, Int64 feature_id=-1). Real production call sites in IonMobilityScoring.cpp (lines 343, 352, 440, 562, 739, 748) pass all of these positionally, e.g. cons

### [FORM-1] Base64::encode — encode() mutates its input vector in place (endian-swaps it), and the doc note claiming the input is 'empty' afterwards is false
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/Base64.h · _format-encoding_

```cpp
template <typename FromType> static void encode(std::vector<FromType> & in, ByteOrder to_byte_order, std::string & out, bool zlib_compression = false)
```
- **Expectation:** A caller reading `encode(in, ...)` expects the source vector to be read-only input, or at worst (per the @note) to be cleared/emptied so it can't be reused.
- **Actual:** The non-const `in` is endian-swapped in place when host/target byte order differ (lines 206-229: `in[i] = tmp.f;`), corrupting the caller's data, but `in` is never cleared. The header @note 'in will be empty after this method' (line 63) is simply wrong: there is no `in.clear()` anywhere in the body. So a caller who trusts the note and reuses `in` finds it neither empty nor unchanged, but silently byte-swapped garbage.
- **Evidence:** Header line 63: `@note @p in will be empty after this method`. Body lines 210-216 swap in place: `Reinterpreter32_ tmp; tmp.f = in[i]; tmp.i = endianize32(tmp.i); in[i] = tmp.f;` — and no `in.clear()` exists in encode() (lines 194-244).
- **Fix:** ABI-safe: fix the doc note to state the truth (input is byte-order-mutated on big-endian or BIGENDIAN target, otherwise unchanged; it is NOT emptied). Ideally add a `const std::vector<FromType>&` overload that copies before swapping (source-compatible additive change). Long term, take input by const-ref and do the swap on an internal copy; that is a breaking signature change (abi_impact breaking) but removes the surprise entirely.
- **Verifier correction:** encode() (and encodeIntegers()) take a non-const std::vector<FromType>& and, when the host byte order differs from the requested ByteOrder, endian-swap the input vector IN PLACE (Base64.h lines 213-215/223-225: `in[i] = tmp.f;`). The vector is never cleared — the header @note "@p in will be empty after this method" (lines 63 and 81) is simply false. Net effect: on the default little-endian host + LITTLEENDIAN target the input is left unchanged (note still wrong, but harmless); on a BIGENDIAN target (e.g. mzXML "network" order) the input is silently mutated to byte-swapped garbage and a caller who reuses it gets corrupt data with no error. Minimal fix (abi none): correct the @note to state the truth — input is byte-order-mutated when host/target endianness differ, otherwise unchanged, and is never emptied. Better fix (source-compatible): add a const std::vector<FromType>& overload that swaps an internal copy. Severity medium, not high, because the corrupting branch is non-default and requires input reuse; no in-tree caller is affected.
- **Verified:** Verified against src/openms/include/OpenMS/FORMAT/Base64.h. The claim is factually accurate on both counts: (1) The header @note (line 63, and identically line 81 for encodeIntegers) states "@p in will be empty after this method" — this is FALSE. The body (lines 196-244) calls out.clear() but contains no in.clear(); the only .clear() calls in the entire header are out.clear() (lines 196,278,312,35

### [FORM-2] Base64::encodeIntegers — encodeIntegers() endian-swaps the caller's input vector in place; @note again claims it becomes 'empty' but it never does
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/Base64.h · _format-encoding_

```cpp
template <typename FromType> static void encodeIntegers(std::vector<FromType> & in, ByteOrder to_byte_order, std::string & out, bool zlib_compression = false)
```
- **Expectation:** Same as encode(): a competent caller expects the input vector to be unchanged, or (per the @note) emptied.
- **Actual:** When byte orders differ, `in` is mutated element by element (lines 366-368, 375-377: `in[i] = tmp;`) and left byte-swapped, not emptied. The @note 'in will be empty after this method' (line 81) is false — no clear() is performed (lines 349-394). Reusing `in` after the call yields corrupted integers on big-endian hosts or BIGENDIAN output.
- **Evidence:** Header line 81: `@note @p in will be empty after this method`. Body lines 364-369: `UInt32 tmp = in[i]; tmp = endianize32(tmp); in[i] = tmp;` with no subsequent `in.clear()`.
- **Fix:** Same as encode(): correct the doc note and, ideally, add a const-ref overload (source-compatible). A const-ref-only signature is the clean fix but is ABI-breaking.
- **Verifier correction:** The @note 'in will be empty after this method' (Base64.h:81) is false: encodeIntegers never calls in.clear() (lines 348-394). Separately, when and ONLY when the host byte order differs from to_byte_order (the branch at line 360 — i.e. BIGENDIAN output on a little-endian host, or any output on a big-endian host), the caller's vector is mutated in place to byte-swapped values (lines 366-368, 375-377) and left that way. On the default little-endian-host + LITTLEENDIAN (mzML) path the input is left unchanged. Net: the documented side effect (emptying) never happens, and an undocumented conditional in-place corruption does. Fix the note and document/eliminate the mutation; a const-ref overload is source-compatible, a const-ref-only signature would be ABI-breaking.
- **Verified:** Evidence verified against the actual code. Header line 81 carries `@note @p in will be empty after this method` for encodeIntegers, and the body (lines 348-394) contains NO in.clear() — only out.clear() at line 351 — so the note is false in every case. Lines 366-368 / 375-377 do mutate in place (in[i] = endianize...(in[i])), leaving the caller's vector byte-swapped, but ONLY inside the conditional

### [FORM-3] MSNumpressCoder::encodeNP / encodeNPRaw — Encoding failures (exceptions, error-tolerance violations) are swallowed; caller gets an empty/unmodified string with no error signal
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/MSNumpressCoder.h · _format-encoding_

```cpp
void encodeNP(const std::vector<double>& in, std::string& result, bool zlib_compression, const NumpressConfig& config)
```
- **Expectation:** An encode entry point that can fail (overflow, accuracy loss beyond tolerance) should signal failure via return value or exception, so the caller does not write a silently-empty binary data array into an mzML file.
- **Actual:** On any thrown exception the implementation only writes to std::cerr and returns, leaving `result` empty (MSNumpressCoder.cpp lines 270-281). When the post-encode error-tolerance check fails it likewise prints to std::cerr and leaves `result` empty (lines 257-261). The header documents this only as a terse '@note In case of error, result string is empty' (line 143) and '@note In case of error, result is given back unmodified' (line 188). There is no boolean return and no throw, so a caller who does not inspect emptiness of `result` silently produces a spectrum with no encoded data.
- **Evidence:** MSNumpressCoder.cpp lines 257-261: `if (n >= 0) { ... std::cerr << "Error occurred at position n = " ...; }` (result left empty). Lines 270-281: catch blocks only `std::cerr << ... ` and fall through. Header line 143: `@note In case of error, result string is empty`.
- **Fix:** ABI-safe: at minimum route the std::cerr diagnostics through OPENMS_LOG_ERROR so failures are visible in normal logging. Ideally add an overload returning bool (or throwing) on encode failure (source-compatible additive). Changing the existing void signatures to throw is ABI source-compatible but a behavioral break for callers relying on the empty-string contract.
- **Verifier correction:** The surprise is real but mis-stated on severity and impact. The primary mzML writer is NOT affected: it checks result.empty() and falls back to base64 (no data loss). The actual residual risk is a swallowed C++ exception in encodeNPRaw (lines 270-281) producing an empty data array in unchecked callers — MzMLSqliteHandler (sqMass) and the Parquet/Mobilogram consumers — all of which use encodeNPRaw's output without an emptiness check. Those callers disable the error-tolerance check (numpressErrorTolerance = -1.0), so the lines-257-261 path is dead for them; only the exception path can silently corrupt their output. Severity is medium (loud-and-recoverable for the main path; latent data loss only for binary-container writers on a rare exception), not high. Minimal fix (route std::cerr through OPENMS_LOG_ERROR) is ABI: none; adding a bool/throwing overload is source-compatible; the unchecked encodeNPRaw callers should additionally verify result non-emptiness before storing.
- **Verified:** Evidence verified line-for-line. MSNumpressCoder::encodeNPRaw (src/openms/source/FORMAT/MSNumpressCoder.cpp) does swallow all encode failures: the error-tolerance check at lines 257-261 prints to std::cerr and leaves `result` empty, and the catch blocks at lines 270-281 (catch int / char const* / ...) only write to std::cerr and fall through, leaving `result` empty. The void signatures never retur

### [FORM-5] ZlibCompression::compressString — compressString takes its INPUT by non-const reference; the header itself admits 'the reference is read-only despite the non-const type'
`medium` · `param-order-or-bool` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/ZlibCompression.h · _format-encoding_

```cpp
static void compressString(std::string& raw_data, std::string& compressed_data)
```
- **Expectation:** A function named compressString(input, output) should accept any std::string input, including const strings and temporaries, since compression does not modify the source.
- **Actual:** The first parameter `raw_data` is a non-const `std::string&`, so `compressString(my_const_str, out)` or `compressString("literal", out)` fails to compile even though the data is only read. The header doc even acknowledges this: 'treated as a raw byte buffer; the reference is read-only despite the non-const type' (line 40). The non-const ref is only there because the impl takes `&str[0]` (ZlibCompression.cpp line 23), not because the input is mutated.
- **Evidence:** Header line 40: `@param[in]  raw_data ... the reference is read-only despite the non-const type`. ZlibCompression.cpp line 23: `compressData(reinterpret_cast<Bytef*>(&str[0]), str.size(), compressed);` — read-only use.
- **Fix:** Change the parameter to `const std::string&` and use `str.data()` (legal and non-mutating since C++11). This is source-compatible (callers passing lvalues still compile) but technically an ABI/mangling change for this symbol (`Ss` vs `RKSs`), so it is a breaking ABI change for the exported function and should be staged accordingly. The non-const overload's existence is the real surprise.
- **Verifier correction:** The surprise is real and the evidence is exact, but the severity should be medium, not high: passing a const string, temporary, or literal produces a loud compile-time error (or forces callers to keep a mutable copy), not silently-wrong output, data loss, or a crash. The recommended fix (const std::string& + str.data()) is source-compatible for all 22 in-tree callers (all pass mutable lvalues) but changes the exported symbol's mangling (R -> RK), making it an ABI-breaking change for this OPENMS_DLLAPI function that should be staged with other ABI breaks.
- **Verified:** All quoted evidence is real and verified independently. Header line 40 verbatim: "@param[in]  raw_data ... (treated as a raw byte buffer; the reference is read-only despite the non-const type)." ZlibCompression.cpp line 23 verbatim: compressData(reinterpret_cast<Bytef*>(&str[0]), str.size(), compressed); — the input str is only read, never mutated. compressString(input, output) takes its read-only

### [FORM-6] ms::numpress::MSNumpress::decodeLinear / decodePic / decodeSlof / encode* / decode* — Decode/encode routines throw a bare `const char*` (not a std::exception) on corrupt input
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MSNUMPRESS/MSNumpress.h · _format-encoding_

```cpp
size_t decodeLinear(const unsigned char *data, const size_t dataSize, double *result)
```
- **Expectation:** A C++ decode helper, when it fails on corrupt input, throws a std::exception-derived type (or returns an error code) so a normal `catch (const std::exception&)` handler catches it.
- **Actual:** These functions throw a raw `const char*` string literal on corruption, as documented: 'this method may throw a const char* if it deems the input data to be corrupt' (header lines 157-159, 174-176, 258-260, 330). A caller using the idiomatic `catch (const std::exception&)` will NOT catch this and will hit std::terminate via an uncaught exception. OpenMS's own wrapper has to add a `catch (...)` to survive (MSNumpressCoder.cpp lines 356-359).
- **Evidence:** Header lines 157-159: `Note that this method may throw a const char* if it deems the input data to be corrupt`. MSNumpressCoder.cpp line 356: `catch (...)` (must use catch-all because the type is not a std::exception).
- **Fix:** This is vendored third-party (pwiz/ms-numpress) code; keep the upstream signature for sync but the surprise is real. Document prominently in the OpenMS wrapper that direct callers must use catch(...) / catch(const char*). The MSNumpressCoder wrapper already converts to Exception::ConversionError — steer public callers to that wrapper. Doc-only / no ABI change.
- **Verifier correction:** Severity is medium, not high. The functions fail loudly (they throw rather than returning silently-wrong data), and the const char* throw is prominently documented at every declaration. The danger is that the non-std::exception type defeats the idiomatic `catch(const std::exception&)`, causing std::terminate for a direct caller on corrupt input — recoverable only via catch(...) or catch(const char*). It is mitigated because the public OpenMS path (MSNumpressCoder) already wraps these in catch(...) and rethrows Exception::ConversionError. Doc-only fix; keep the vendored upstream signature for sync.
- **Verified:** Independently verified. Header lines 157-159, 174-176, 258-260, 275-277, 330, 345 document that decodeLinear/decodePic/decodeSlof/decodeSafe "may throw a const char*" on corrupt input. The .cpp confirms bare string-literal throws (e.g. MSNumpress.cpp:213 `throw "[MSNumpress::decodeInt] Corrupt input data! ";`, plus :429/:436/:451/:585/:629/:822 for decode and :349/:358/:658/:799 for encode overflo

### [FORM-36] ChromeleonFile::removeCommasAndParseDouble — const parse helper takes its argument by non-const reference and mutates it in place
`medium` · `hidden-side-effect` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/ChromeleonFile.h · _format-feature-fasta_

```cpp
double removeCommasAndParseDouble(std::string& number) const;
```
- **Expectation:** A method named 'removeCommasAndParseDouble' that is declared const and 'parses' a number reads the string and returns a double. The Doxygen says '@param[in] number', strongly implying the input is not modified.
- **Actual:** The implementation calls StringUtils::remove(number, ',') on the caller's string, permanently stripping the commas from the argument before parsing. The method is const (it doesn't modify the ChromeleonFile object) yet it mutates a caller-owned string passed as a documented [in] parameter.
- **Evidence:** Header: `double removeCommasAndParseDouble(std::string& number) const;` with `@param[in] number`. Impl (ChromeleonFile.cpp:125): `StringUtils::remove(number, ','); return StringUtils::toDouble(number);`
- **Fix:** Take the argument by value (`double removeCommasAndParseDouble(std::string number) const`) or by const-ref and copy internally, so the caller's string is not clobbered. If ABI stability forbids changing the signature, at minimum change the Doxygen tag to '@param[in,out]' to warn callers. The by-value fix is source-compatible for typical call sites passing temporaries but is an ABI break (mangled signature change).
- **Verifier correction:** The surprise is real: a const-qualified helper documented as taking `@param[in] number` mutates the caller's string in place (strips all commas via StringUtils::remove, which erases on the reference). But severity is medium, not high — the only call site passes throwaway split substrings, so no incorrect data or loss results in practice, and the non-const-ref parameter prevents binding temporaries/literals (the worst misuse would fail to compile rather than silently corrupt). The minimal correct fix is documentation (change tag to @param[in,out]); changing the signature to by-value or const-ref to truly honor [in] is an ABI break (mangled-name change).
- **Verified:** Code confirms the claim. Header (ChromeleonFile.h:50) declares `double removeCommasAndParseDouble(std::string& number) const;` with Doxygen `@param[in] number`. Impl (ChromeleonFile.cpp:125-129) calls `StringUtils::remove(number, ',')` then `StringUtils::toDouble(number)`. Verified StringUtils::remove (StringUtils.h:665) does `s.erase(std::remove(...))` on the by-reference argument, so the caller'

### [FORM-18] FileHandler::computeFileHash — computeFileHash returns empty string on a missing/unreadable file instead of signaling failure
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/FileHandler.h · _format-filehandling_

```cpp
static std::string computeFileHash(const std::string& filename)
```
- **Expectation:** A function documented as '@return The SHA-1 hash of the given file' should either return a valid hash or signal failure (throw FileNotFound, like the rest of FileHandler's file ops). A caller storing/comparing the checksum has no reason to expect a sentinel.
- **Actual:** If the file cannot be opened, the function returns "" (FileHandler.cpp:754-757). The empty string is a silent sentinel: it is not a valid SHA-1 (which is 40 hex chars; even the hash of empty input is da39a3ee...), so a caller that stores it into SourceFile::setChecksum() records a meaningless empty checksum with no error. loadExperiment(..., compute_hash=true) calls this and stores the result unconditionally (FileHandler.cpp:1032).
- **Evidence:** FileHandler.cpp:751 'std::string FileHandler::computeFileHash(const std::string& filename) { std::ifstream file{...}; if (!file.is_open()) { return ""; } ...'. Header lines 405-409 document only '@return The SHA-1 hash of the given file' with no failure semantics.
- **Fix:** Throw Exception::FileNotFound when the file cannot be opened (consistent with the rest of FileHandler). If ABI/behavior compatibility with existing callers tolerating "" is a concern, at minimum document the empty-string-on-failure contract in the header. The throw is source-compatible (callers currently cannot meaningfully use "").
- **Verifier correction:** computeFileHash returns "" on an unopenable file, an undocumented out-of-band sentinel (no valid SHA-1 is empty). The header promises only a SHA-1 with no failure contract, and the project's own unit test treats "" as the failure marker, confirming the behavior is intentional-but-undocumented. Real POLS gap. However the primary internal caller (loadExperiment) only hashes a file the loader has already opened successfully, so it does not routinely store meaningless checksums; the realistic exposure is a direct external/pyOpenMS caller passing a missing path, where the empty result is silent but trivially detectable rather than a plausible wrong value. Fix: throw Exception::FileNotFound (source-compatible) or at minimum document the empty-string-on-failure contract in the header. Severity medium, not high.
- **Verified:** Verified against source. FileHandler.cpp:754-757 returns "" when the file cannot be opened; the header (FileHandler.h:404-409) documents only "@return The SHA-1 hash of the given file" with no failure semantics. "" is a genuine out-of-band sentinel: a valid SHA-1 is always 40 hex chars (even empty input -> da39a3ee...), so it can never collide with a real result. The method is public, static, and 

### [FORM-19] ControlledVocabulary::CVTerm::isHigherBetterScore — isHigherBetterScore is a static method taking a CVTerm by value, and returns true for ANY (even non-score) term by default
`medium` · `misleading-name` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/ControlledVocabulary.h · _format-filehandling_

```cpp
static bool isHigherBetterScore(ControlledVocabulary::CVTerm term)
```
- **Expectation:** Named like a predicate on 'this' CVTerm ('is this score higher-better?'), a caller expects to call term.isHigherBetterScore() and to get a meaningful answer only for score terms (the inline comment even says 'if it is a score type, lookup has_order'). One would expect non-score terms to be rejected or to return false.
- **Actual:** It is a static method that takes the term BY VALUE (a copy) rather than operating on *this (CVTerm.cpp:93). It never checks whether the term is a score type at all; it scans unparsed lines only for the 'has_order MS:1002109' (lower-better) relationship and returns true otherwise. So for any term lacking that exact line — including non-score terms and unknown terms — it returns true ('higher is better'). The 'is...Score' name and the documented intent ('if it is a score type') are both contradicted.
- **Evidence:** Header line 69: 'static bool isHigherBetterScore(ControlledVocabulary::CVTerm term); ///if it is a score type, lookup has_order'. CVTerm.cpp:93-106: 'bool ControlledVocabulary::CVTerm::isHigherBetterScore(ControlledVocabulary::CVTerm term) { ... for (... term.unparsed ...) if (hasPrefix(*unp, "relationship: has_order MS:1002109")) return false; return true; }'.
- **Fix:** Rename to something honest like hasLowerBetterOrder() / defaultsToHigherBetter(), and make it a const member taking no argument (or accept by const-ref to avoid the by-value copy). At minimum, change the signature to 'bool isHigherBetterScore() const' (uses *this) and pass-by-const-ref if kept static. Renaming is breaking; an additive const overload plus a [[deprecated]] alias preserves ABI.
- **Verified:** All quoted evidence is real and accurate. Header (ControlledVocabulary.h:69): `static bool isHigherBetterScore(ControlledVocabulary::CVTerm term); ///if it is a score type, lookup has_order`. Implementation (ControlledVocabulary.cpp:93-106): the live code is static, takes the term BY VALUE (a copy), and does NOT operate on *this — the original *this-based body is commented out at lines 95-99. It p

### [FORM-112] GNPSMGFFile::store — store() never checks the output stream and silently produces nothing on an unwritable path
`medium` · `silent-failure` · ABI: `none` · src/openms/source/FORMAT/GNPSMGFFile.cpp · _format-gnps-deconv-qc_

```cpp
void store(const std::string& consensus_file_path, const StringList& mzml_file_paths, const std::string& out) const
```
- **Expectation:** A *File::store() that cannot create/open the destination signals failure (throw UnableToCreateFile), as the sibling MzQCFile::store does.
- **Actual:** It opens 'ofstream output_file(out);' with no success check; writeMSMSBlockHeader_/writeMSMSBlock_ are guarded by 'if (output_file.is_open())' and silently no-op when the file is not open. A bad/locked/permission-denied 'out' yields a successful-looking call that wrote nothing.
- **Evidence:** GNPSMGFFile.cpp:248 'ofstream output_file(out);' with no check; writeMSMSBlockHeader_ (line 149) and writeMSMSBlock_ (line 173) both begin 'if (output_file.is_open())'. Compare MzQCFile.cpp:46-50 which throws Exception::UnableToCreateFile when '!os'.
- **Fix:** Add 'if (!output_file) throw Exception::UnableToCreateFile(...);' right after opening, matching MzQCFile. Additive, no signature change.
- **Verified:** Evidence confirmed verbatim. GNPSMGFFile.cpp:248 'ofstream output_file(out);' has no stream-state check, and writeMSMSBlockHeader_ (line 149) and writeMSMSBlock_ (line 173) both guard with 'if (output_file.is_open())', so on an unwritable/locked/permission-denied 'out' the ConsensusMap still loads, the loop still runs with progress logged, every write silently no-ops, and store() returns normally 

### [FORM-114] GNPSMetaValueFile::store — GNPS TSV writer does not verify the output stream opened (inconsistent with MzQCFile)
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/source/FORMAT/GNPSMetaValueFile.cpp · _format-gnps-deconv-qc_

```cpp
static void store(const ConsensusMap& consensus_map, const std::string& output_file)
```
- **Expectation:** Across the FORMAT *File writers, store() to an unopenable path fails loudly.
- **Actual:** It opens 'std::ofstream outstr(output_file.c_str());' and writes via SVOutStream without ever testing the stream state, so an unwritable destination silently yields a truncated/empty file. The neighboring MzQCFile::store throws UnableToCreateFile in the same situation, making the cluster inconsistent.
- **Evidence:** GNPSMetaValueFile.cpp:26-27 'std::ofstream outstr(output_file.c_str()); SVOutStream out(outstr, ...);' with no '!outstr' check. Same pattern in GNPSQuantificationFile.cpp:35. Contrast MzQCFile.cpp:46-50.
- **Fix:** Add a '!outstr' check throwing Exception::UnableToCreateFile after opening, in both GNPS TSV/TXT writers, to match MzQCFile. Additive.
- **Verifier correction:** Severity is medium rather than high: the harm (silently truncated/empty output) is real but triggered by an environmental edge case (unwritable path / full disk), not normal-path use, and downstream GNPS upload typically surfaces the bad file. The claim is otherwise fully accurate. Worth noting the codebase already provides a safe alternative — the SVOutStream(const std::string& file_out, ...) constructor throws Exception::FileNotWritable on open failure — so the cleanest fix is either adding the `if (!outstr) throw Exception::UnableToCreateFile(...)` check (as recommended, additive) or switching to that file-path SVOutStream constructor in both GNPS writers.
- **Verified:** Independently verified all quoted evidence. GNPSMetaValueFile.cpp:26-27 opens `std::ofstream outstr(output_file.c_str())` and passes it to the `SVOutStream(ostream&, ...)` constructor (SVOutStream.cpp:40-48), which performs NO stream-state check; the function writes and returns with no `!outstr` test. GNPSQuantificationFile.cpp:35 is identical. MzQCFile.cpp:46-50, same FORMAT dir and same maintain

### [FORM-116] QcMLFile::collectSetParameter — collectSetParameter appends to 'ret' without clearing and silently inserts empty map entries for unknown set/run names
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/QcMLFile.h · _format-gnps-deconv-qc_

```cpp
void collectSetParameter(const std::string& setname, const std::string& qp, std::vector<std::string>& ret)
```
- **Expectation:** A 'collect...into ret' accessor fills ret with the requested values (clearing prior content, like the other QcMLFile out-param getters) and does not mutate the file's stored data.
- **Actual:** It does not clear 'ret' (it only push_back's, so prior caller content survives) -- unlike existsRunQualityParameter/getRunIDs/getRunNames which all clear() first. It is also non-const and uses operator[] on setQualityQPs_members_[setname] and runQualityQPs_[*it], so querying an unknown set or member name silently default-inserts empty entries into the maps as a side effect of a logically read-only 'collect'.
- **Evidence:** QcMLFile.cpp:750-762: loops 'setQualityQPs_members_[setname]' and 'runQualityQPs_[*it]' (both operator[] inserts), 'ret.push_back(jt.value)' with no 'ret.clear()'. Header line 126 even carries a vestigial 'void/* std::vector<std::string>& */ collectSetParameter' signature comment.
- **Fix:** Clear 'ret' at entry to match sibling getters, and replace operator[] with find()/at() to avoid mutating the maps on lookup (then the method can be made const). Behavior-only fix; source/ABI compatible.
- **Verifier correction:** collectSetParameter has two confirmed hidden side effects: it appends to 'ret' without clearing (unlike all sibling getters, which clear()), and it uses std::map::operator[] on setQualityQPs_members_ and runQualityQPs_, silently default-inserting empty map entries for unknown set/run names and thereby mutating the in-memory model from a read-only 'collect'. The phantom inserted entries can later surface as spurious empty runs/sets in getRunIDs/getSetNames iteration. Severity is medium rather than high: it does not corrupt persisted files, the main pyOpenMS caller passes a fresh empty vector, and the condition requires either a reused output vector or an unknown-key query — recoverable rather than a guaranteed crash/persistent data loss. Fix: clear() at entry and replace operator[] with find()/at() (allowing the method to be made const, which is source-compatible for callers but changes the symbol's mangling).
- **Verified:** Independently verified against QcMLFile.cpp:750-762 and QcMLFile.h:126,185-189. All quoted evidence is accurate: (1) collectSetParameter does NOT clear 'ret' (only push_back), whereas every sibling out-param getter in the same class clears first — getRunNames (line 313), getRunIDs (line 319), existsRunQualityParameter (line 362), existsSetQualityParameter (line 386) all call ids.clear() at entry. 

### [FORM-117] MzQCFile::store — store() leads with two adjacent same-typed path strings (input_file, output_file) that are trivially swappable, and the class doc promises load+store but only store exists
`medium` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/MzQCFile.h · _format-gnps-deconv-qc_

```cpp
void store(const std::string& input_file, const std::string& output_file, ... ) const
```
- **Expectation:** A store() either takes the destination first (consistent with other *File::store(out, data)) or makes the in/out paths hard to confuse; a class documented for 'load and store' offers both.
- **Actual:** The first two parameters are 'input_file' (mzML being summarized) then 'output_file' (mzQC destination), both std::string with no type-level protection; swapping them compiles and would read the destination path as the mzML and write QC into the source path. Additionally the header brief says 'used to load and store mzQC files' but the class exposes only store() (no load), an asymmetric pair.
- **Evidence:** MzQCFile.h:45-54 'void store(const std::string& input_file, const std::string& output_file, ...)'. MzQCFile.h:20 brief: 'File adapter for mzQC files used to load and store mzQC files'. No load() member is declared.
- **Fix:** Cannot reorder without breaking ABI; at minimum fix the brief to say 'store' only (or add a load()). Consider a future typed wrapper for in/out paths. Doc fix is non-breaking.
- **Verifier correction:** store() takes the destination path SECOND (output_file), breaking the OpenMS convention where every other *File::store leads with the single destination 'filename'. The two adjacent same-typed std::string paths are swappable; a swap does NOT cause the destination to be read as mzML (mzML data comes from the separate `exp` MSExperiment arg) — instead input_file is used only for provenance (path/basename/file-hash), so a swap silently writes QC JSON to the intended-source path (potential overwrite) and records the destination as provenance. Separately, the class brief claims "load and store" but only store() exists (no load()), a misleading asymmetric doc.
- **Verified:** The core surprise is real and provable, but the claim mis-states the swap mechanism. CONFIRMED: MzQCFile::store leads with two adjacent same-typed path strings (input_file, output_file), and uniquely among OpenMS *File classes the FIRST string is NOT the write destination. Every other store() (MzMLFile, MzXMLFile, FeatureXMLFile, IdXMLFile, MzTabFile) takes a single destination 'filename' first; M

### [FORM-118] UnimodXMLFile::load — load()'s 'modifications' is documented '@param[in]' but is the OUTPUT and yields caller-owned raw pointers
`medium` · `ownership-lifetime` · ABI: `none` · src/openms/include/OpenMS/FORMAT/UnimodXMLFile.h · _format-gnps-deconv-qc_

```cpp
void load(const std::string& filename, std::vector<ResidueModification*>& modifications)
```
- **Expectation:** An '@param[in]' vector is read by the function; load() output should be marked [out], and the ownership of any returned heap objects should be stated.
- **Actual:** The Doxygen says '@param[in] modifications the modifications which are read from the file', but load() populates this vector with newly-allocated ResidueModification* (the handler 'new's them). It is an output parameter, not an input, and the caller becomes responsible for deleting the raw pointers, which the [in] annotation and signature do not convey.
- **Evidence:** Header UnimodXMLFile.h:35,42 '@param[in] modifications the modifications which are read from the file' / 'void load(const std::string & filename, std::vector<ResidueModification *> & modifications);'. Impl UnimodXMLFile.cpp:31 'Internal::UnimodXMLHandler handler(modifications, file);' fills the vector with handler-allocated pointers.
- **Fix:** Change the doc to '@param[out]' and document ownership (caller must delete, or prefer unique_ptr in a new overload). Doc-only fix is non-breaking; switching to smart pointers would be source-breaking and should be additive.
- **Verifier correction:** load()'s `modifications` parameter is mislabeled `@param[in]` but is actually an output parameter into which load() places newly heap-allocated `ResidueModification*` objects whose ownership transfers to the caller (caller must delete or wrap them). The recommended fix is doc-only: change to `@param[out]` and state that the caller owns and must delete (or wrap) the returned pointers; this is non-breaking (abi_impact none). Severity is medium (leak hazard, recoverable) rather than high.
- **Verified:** Verified against source. Header UnimodXMLFile.h:36 documents `@param[in] modifications the modifications which are read from the file`, while load() (UnimodXMLFile.cpp:27-33) passes the vector to UnimodXMLHandler, which `new`s ResidueModification objects (UnimodXMLHandler.cpp:40 and :167) and pushes raw owning pointers into it. So `modifications` is an OUTPUT, not input, and the caller becomes res

### [FORM-133] XMLHandler::cvParamToValue — cvParamToValue signals conversion/lookup failure only via a sentinel DataValue::EMPTY
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/HANDLERS/XMLHandler.h · _format-handlers-core_

```cpp
DataValue cvParamToValue(const ControlledVocabulary& cv, const std::string& parent_tag, const std::string& accession, const std::string& name, const std::string& value, const std::string& unit_accession) const
```
- **Expectation:** A converter named cvParamToValue that '...converts the value ... to the DataValue with the correct type' would, on a hard failure (e.g. an integer-typed CV term with a non-integer value, or an entirely unknown accession), either throw or otherwise force the caller to notice.
- **Actual:** On a type-conversion failure it logs a warning and `return DataValue::EMPTY;` (lines 254, 267, 281, 294, 311). On an unknown accession in a non-'sample' tag it likewise warns and returns DataValue::EMPTY (line 320). The doc does note 'DataValue::EMPTY if a conversion error occured', but EMPTY is also a legitimate value (a cvParam with no value yields a string DataValue), so a caller that does not explicitly compare against DataValue::EMPTY silently drops malformed values. Warnings are additionally suppressed to DEBUG level in release builds (XMLHandler.cpp:145-150), so the failure is effectively invisible.
- **Evidence:** XMLHandler.cpp:251-255 `catch (Exception::ConversionError&) { warning(LOAD, ...); return DataValue::EMPTY; }`; XMLHandler.cpp:314-321 catch InvalidValue -> warning + `return DataValue::EMPTY;`. Header doc lines 465/477: 'DataValue::EMPTY if a conversion error occured ... or the DataValue upon success'.
- **Fix:** Keep the EMPTY-sentinel for ABI stability but make the contract explicit in the docs that callers MUST test the result against DataValue::EMPTY, and consider an additive overload that takes a `bool& ok` out-param or throws on hard conversion errors. Documenting/adding overload is source-compatible.
- **Verifier correction:** cvParamToValue signals all hard conversion/lookup failures by returning the sentinel DataValue::EMPTY after only a warning() that is downgraded to DEBUG in release builds. The result is documented (header 465/477) but the doc does not state that callers MUST test for EMPTY; consequently one of the two in-tree callers (MzIdentMLDOMHandler.cpp:2185-2186) does not check and silently stores an empty meta value for malformed input, while the other (MzMLHandler.cpp:1426) does check. The genuine semantic-overlap hazard is that DataValue::EMPTY is also the standard 'meta value not found' return (MetaInfoInterface.cpp:122,140), so dropped values become indistinguishable from never-set ones — NOT, as the claim states, that a valueless cvParam yields EMPTY (a valueless cvParam actually yields a string "" which compares unequal to EMPTY). Impact is limited to malformed/recoverable input files and is loud in debug builds, hence medium rather than high. The recommended fix (clarify the doc that callers must test against EMPTY, plus an additive bool& ok / throwing overload) is source-compatible; the existing signature/ABI is unchanged.
- **Verified:** Code confirmed. XMLHandler::cvParamToValue (src/openms/source/FORMAT/HANDLERS/XMLHandler.cpp:201-343) returns DataValue::EMPTY on every hard failure: integer cast fail (254), double cast fail (267), date parse fail (281), bool fail (294), missing-numerical-value (311), and unknown accession in a non-'sample' tag (320). Each is preceded only by warning(LOAD,...). warning() (XMLHandler.cpp:130-152) 

### [FORM-134] XMLHandler::asInt_ / asUInt_ / asDouble_ / asFloat_ / asBool_ — as*_ conversion helpers swallow ConversionError and return 0/false on bad input
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/HANDLERS/XMLHandler.h · _format-handlers-core_

```cpp
Int asInt_(const std::string & in) const; UInt asUInt_(const std::string & in) const; double asDouble_(const std::string & in) const; float asFloat_(const std::string & in) const; bool asBool_(const std::string & in) const
```
- **Expectation:** A `const` conversion helper `asInt_("abc")` would be expected to throw (like StringUtils::toInt32 does) or otherwise propagate the failure, so a malformed attribute does not become a valid-looking 0.
- **Actual:** Each helper catches ConversionError, calls the non-fatal `error()` (which only logs — and only at DEBUG level in release builds), and returns the default-initialized value (`res = 0` / `0.0` / `false`). So `asInt_("not-a-number")` returns 0 and parsing continues as if the value were legitimately zero. `asBool_` returns `false` for unparseable input after merely logging.
- **Evidence:** XMLHandler.h:536-548 `Int asInt_... { Int res = 0; try { res = StringUtils::toInt32(in); } catch (Exception::ConversionError&) { error(LOAD, ...); } return res; }`; identical shape for asUInt_/asDouble_/asFloat_/asBool_ (lines 557-628). `error(LOAD,...)` at XMLHandler.cpp:111-127 only logs, never throws.
- **Fix:** These are protected helpers, but they participate in load correctness. Document that they return a default on failure and rely on error() routing; for stricter formats prefer the attributeAs*_ paths (which fatalError on missing). If behavior change is acceptable, route through fatalError for required numeric fields. Doc clarification is non-breaking; behavior change is source-compatible at most.
- **Verifier correction:** The as*_ helpers do swallow ConversionError and return default 0/0.0/false on bad input, but error(LOAD,...) is NOT logged 'only at DEBUG level in release builds' — it uses OPENMS_LOG_ERROR unconditionally (always ERROR level; LogStream.h:577, XMLHandler.cpp:126). The DEBUG-in-release downgrade applies to warning() (XMLHandler.cpp:145-150), not error(). Thus malformed input is always reported loudly at error level, but the call still returns 0/false and parsing continues, so the value becomes a valid-looking zero. The surprise is the silent return-of-default, not silent logging.
- **Verified:** Verified independently. The std::string overloads asInt_/asUInt_/asDouble_/asFloat_/asBool_ (XMLHandler.h:536-628) catch Exception::ConversionError (asBool_ handles the unrecognized-value else branch), call the non-fatal error(LOAD,...) which only logs and never throws (XMLHandler.cpp:111-127), and return the default-initialized 0/0.0/false. The underlying StringUtils::toInt32/toDouble/toFloat gen

### [FORM-135] XMLHandler::optionalAttributeAsString_ — Two optionalAttributeAsString_ overloads disagree on what the bool return means
`medium` · `return-value` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/HANDLERS/XMLHandler.h · _format-handlers-core_

```cpp
bool optionalAttributeAsString_(std::string& value, const Attributes& a, const char* name) const  vs  bool optionalAttributeAsString_(std::string& value, const Attributes& a, const XMLCh* name) const
```
- **Expectation:** All `optionalAttributeAs*_` overloads share the documented contract '@return if the attribute was present'. The two string overloads (char* vs XMLCh* name) should return the same thing for the same input.
- **Actual:** The `const char*` overload returns `true` whenever the attribute is present, regardless of value (line 715-724: sets value, `return true;`). The `const XMLCh*` overload returns `!value.empty()` (line 877-886: `value = sm_.convert(val); return !value.empty();`). So for a present-but-empty attribute, the two overloads return different booleans (true vs false), even though both are documented identically as 'if the attribute was present'.
- **Evidence:** XMLHandler.h:715-724 (char* overload) `value = sm_.convert(val); return true;` versus XMLHandler.h:877-886 (XMLCh* overload) `value = sm_.convert(val); return !value.empty();`.
- **Fix:** Make both overloads consistent (return `val != nullptr`, i.e. presence only, matching the documented contract) and treat empty-value handling separately. Changing the XMLCh* overload to return presence is source-compatible (same signature); callers relying on the empty-string distinction should be audited.
- **Verifier correction:** Severity is medium, not high. The divergence only manifests for the edge case of a present-but-EMPTY attribute value; for any non-empty or absent value both overloads agree. When it does trigger on the XMLCh* path, it fails closed (field left unset, same as the "absent" branch) rather than producing wrong numeric results, data corruption of valid data, or a crash. It is a silent, recoverable inconsistency that violates the documented presence contract and invites misuse — hence medium. The recommended fix (return `val != nullptr` to match the rest of the family) preserves the signature and is source-compatible; only runtime behavior on the empty-string edge case changes, so any caller depending on the present-but-empty=>false distinction (e.g. the if-guarded XMLCh* sites in MzMLHandler) should be audited.
- **Verified:** Independently verified in XMLHandler.h. The char* overload (lines 715-724) does `value = sm_.convert(val); return true;` (presence only), while the XMLCh* overload (lines 877-886) does `value = sm_.convert(val); return !value.empty();`. Both carry the identical doc comment "@return if the attribute was present". Critically, EVERY other overload in the family (optionalAttributeAsInt_, UInt, Double,

### [FORM-138] XMLHandler::warning — warning() is silently downgraded to DEBUG-level logging in release builds
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/HANDLERS/XMLHandler.h · _format-handlers-core_

```cpp
void warning(ActionMode mode, const std::string & msg, UInt line = 0, UInt column = 0) const
```
- **Expectation:** A method named warning() that takes a 'Warning handler' role would emit a visible warning to the user/log when called.
- **Actual:** In release builds (no OPENMS_ASSERTIONS) the message is routed to OPENMS_LOG_DEBUG rather than OPENMS_LOG_WARN, with the in-source justification '// warn only in Debug mode but suppress warnings in release mode (more happy users)'. Combined with cvParamToValue/as*_ helpers that report malformed data only via warning(), this means data-quality problems during load are effectively invisible to end users.
- **Evidence:** XMLHandler.cpp:145-150 `#ifdef OPENMS_ASSERTIONS\n OPENMS_LOG_WARN << error_message ...\n #else\n OPENMS_LOG_DEBUG << error_message ...\n #endif`.
- **Fix:** Document on the header declaration that warning() is suppressed to DEBUG in release builds so callers do not assume their warnings are user-visible; or gate the suppression behind an option rather than build type. Doc/behavioral; non-breaking to document.
- **Verifier correction:** Minor wording fix: the silent-warning-only path is cvParamToValue_ (malformed CV term values -> warning() + return DataValue::EMPTY). The numeric/date string conversion helpers asInt_/asDouble_/asFloat_/asBool_/asDateTime_ use error() (OPENMS_LOG_ERROR, always visible), not warning(). Severity is medium rather than high: malformed CV data silently becomes DataValue::EMPTY with no user-visible diagnostic in default release builds, which is silent data-quality loss on reasonable input, but the affected path returns a defined sentinel rather than corrupting otherwise-valid results, and it is recoverable by raising the debug level. Recommendation stands: document the release-build DEBUG-downgrade on the warning() declaration, or gate the suppression behind a runtime option rather than build type. Doc/behavioral; non-breaking.
- **Verified:** Independently verified in source. XMLHandler.cpp:129-152 defines warning(ActionMode,...) and at lines 145-150 routes the message to OPENMS_LOG_WARN only #ifdef OPENMS_ASSERTIONS, else OPENMS_LOG_DEBUG, with the exact in-source justification '// warn only in Debug mode but suppress warnings in release mode (more happy users)'. OPENMS_ASSERTIONS is gated Debug-only in src/openms/CMakeLists.txt:27-31

### [FORM-139] AcqusHandler::getParam — getParam() is the only non-const accessor and silently inserts an empty entry for unknown keys
`medium` · `hidden-side-effect` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/HANDLERS/AcqusHandler.h · _format-handlers-misc_

```cpp
std::string getParam(const std::string & param);
```
- **Expectation:** getParam() on a 'Read-only acqus File handler' should be a const lookup that returns the stored value (or empty) without mutating the handler; siblings getPosition()/getSize() are const.
- **Actual:** getParam() is non-const and returns params_[param]. Because params_ is a std::map, operator[] default-inserts an empty string for any key not present, mutating the internal map on a pure lookup. This is also why it cannot be const. It silently grows state and 'succeeds' (returns "") for misspelled/absent params instead of signaling not-found.
- **Evidence:** Header comment line 19 "Read-only acqus File handler"; getPosition()/getSize() declared const but "std::string getParam(const std::string & param);" is not. Impl AcqusHandler.cpp: "return params_[param];"
- **Fix:** Change the body to use params_.find()/at() and mark the method const (ABI-breaking signature change because const-qualification mangles into the symbol). If ABI must be preserved, at minimum switch the lookup to find() to stop the silent insertion (source-compatible, fixes the const-correctness violation in spirit) and document that unknown keys return "".
- **Verifier correction:** getParam() is non-const and returns params_[param], so std::map::operator[] default-inserts an empty string for absent keys, mutating the internal map on what looks like a pure read of a "Read-only" handler (this is why it can't be const, unlike the const getPosition()/getSize()). Severity is medium, not high: the only callers (XMassFile.h) intentionally use it for optional metadata and treat the empty return as "not present", so unknown/misspelled keys yield empty metadata strings rather than wrong scientific results, data loss, or crashes. Recommended fix (use find()/at() and mark const) is ABI-breaking due to const name mangling; a find()-only change preserving the signature is source-compatible.
- **Verified:** All quoted evidence is accurate. Header (src/openms/include/OpenMS/FORMAT/HANDLERS/AcqusHandler.h) documents the class as a "Read-only acqus File handler" (line 19); getPosition() and getSize() are declared const (lines 45, 51) but getParam() is not (line 48). The impl (src/openms/source/FORMAT/HANDLERS/AcqusHandler.cpp line 116) is `return params_[param];`, and std::map::operator[] default-insert

### [FORM-140] CachedMzMLHandler::writeMetadata — writeMetadata() takes the (large) MSExperiment BY VALUE while its sibling writeMetadata_x() takes it by const-ref
`medium` · `pass-by-value-vs-const-ref / overload-naming` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h · _format-handlers-misc_

```cpp
void writeMetadata(MapType exp, const std::string& out_meta, bool addCacheMetaValue=false);
```
- **Expectation:** A write/store method named writeMetadata(exp, out) should not copy a whole PeakMap; const& is expected, and the two near-identical overloads writeMetadata vs writeMetadata_x should share a calling convention.
- **Actual:** writeMetadata takes MapType exp by value (a full deep copy of the experiment) purely so it can clear the spectra in place, whereas writeMetadata_x(const MapType& exp, ...) takes the same argument by const reference. The .cpp even carries a '// TODO : remove copy' comment. A caller picking writeMetadata pays a silent full-map copy; the _x suffix gives no hint that it is the const-ref (cheaper) variant.
- **Evidence:** Header: "void writeMetadata(MapType exp, ...)" vs "void writeMetadata_x(const MapType& exp, ...)". Impl writeMetadata: "// TODO : remove copy" and mutates exp[i].clear(false) on the by-value copy.
- **Fix:** Keep the by-value overload for ABI but document the copy cost; recommend writeMetadata_x for hot paths and consider eventually deprecating the by-value overload. A true fix (take const& and build the stripped copy internally) changes the signature/mangling and is ABI-breaking, so do it additively or at a major version.
- **Verifier correction:** writeMetadata(MapType exp, ...) takes the full MSExperiment/PeakMap BY VALUE (deep-copying all peak data) only to immediately strip it via exp[i].clear(false), as flagged by the in-source '// TODO : remove copy'. Its public sibling writeMetadata_x(const MapType& exp, ...) is the const-ref variant that avoids the peak-data copy by building the stripped map internally, but the '_x' suffix gives no hint it is the cheaper/preferred form. Almost all callers use the expensive by-value overload. This is a silent performance cliff (correct output, but an avoidable full-experiment copy) selected by the more natural-looking name; not a param-order or bool-trap issue. A true fix (take const& and strip internally) changes the mangled signature and is therefore ABI-breaking; it must be done additively or at a major version.
- **Verified:** Evidence verified verbatim. Header (src/openms/include/OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h) declares two PUBLIC overloads: line 77 `void writeMetadata(MapType exp, const std::string& out_meta, bool addCacheMetaValue=false)` (by value) and line 80 `void writeMetadata_x(const MapType& exp, ...)` (const-ref). MapType is PeakMap/MSExperiment (line 47), so by-value is a full deep copy of all pea

### [FORM-141] MzMLSqliteHandler::readSpectra — readSpectra/readChromatograms APPEND to the output vector and then throw if it was non-empty, instead of clearing it
`medium` · `misleading-diagnostic` · ABI: `none` · src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h · _format-handlers-misc_

```cpp
void readSpectra(std::vector<MSSpectrum> & exp, const std::vector<int> & indices, bool meta_only = false) const;
```
- **Expectation:** A read*(out, ...) method is expected to fill the output vector with exactly the requested items, clearing any prior contents (load-clears-then-fills convention used elsewhere).
- **Actual:** readSpectra does not clear exp; prepareSpectra_ appends via spectra.resize(spectra.size()+1) per row. readSpectra then enforces indices.size() == exp.size(), so passing a pre-filled vector causes a (correct-but-surprising) Exception::IllegalArgument 'Illegal spectral indices detected'. The error message points at the indices, not at the real cause (caller passed a non-empty vector). A reused output vector silently turns a valid call into a throw.
- **Evidence:** MzMLSqliteHandler.cpp readSpectra: "prepareSpectra_(conn.getDB(), exp, indices); if (indices.size() != exp.size()) throw ...\"Illegal spectral indices detected\"". prepareSpectra_ body: "spectra.resize(spectra.size() + 1);" (no exp.clear()).
- **Fix:** Clear the output vector at the start of readSpectra/readChromatograms (and have prepareSpectra_/prepareChroms_ not assume an empty input), or document that the vector must be empty on entry. Adding exp.clear() is source-compatible and ABI-safe; it only changes behavior for the already-broken reuse case.
- **Verifier correction:** readSpectra/readChromatograms do not clear their output vector and rely on prepareSpectra_/prepareChroms_ appending (resize(size()+1)). Combined with the post-condition `indices.size() != exp.size()`, passing a non-empty output vector throws Exception::IllegalArgument whose message ("Illegal spectral indices detected ...") misattributes the failure to the indices rather than the reused (non-empty) vector. This is a LOUD failure (exception), not silent — re-categorize from "silent-failure" to a misleading-diagnostic / non-clearing-output-contract issue. No wrong results, data loss, or crash; the reuse case is recoverable. All current OpenMS callers pass fresh empty vectors, so the surprise is only reachable by external/pyOpenMS code. Recommended fix (exp.clear() at function start) is correct, source-compatible, and ABI-impact none (a .cpp body change, no signature/header change). The header @param[in]/@param[out] tagging for `exp` is also inconsistent (readSpectra tags it @param[in]).
- **Verified:** Code facts verified and accurate. In src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp, readSpectra (lines 391-412) and readChromatograms (414-437) never clear the output vector. prepareSpectra_ appends via `spectra.resize(spectra.size() + 1)` (line 767) and prepareChroms_ via `chromatograms.resize(chromatograms.size()+1)` (line 646). The guard `if (indices.size() != exp.size()) throw Except

### [FORM-48] SequestInfile::setEnzyme — "setEnzyme" setter returns an inverted error code (0 = success, non-zero = enzyme NOT found)
`medium` · `return-value` · ABI: `none` · src/openms/include/OpenMS/FORMAT/SequestInfile.h · _format-id-search-in_

```cpp
Size setEnzyme(const std::string& enzyme_name)
```
- **Expectation:** A method documented "sets the enzyme used for cleavage" and named setEnzyme reads as a void setter; if it returns anything, a caller assumes a truthy value means success.
- **Actual:** It returns Size. The implementation returns `(einfo_i == enzyme_info_.end()) ? enzyme_info_.size() : 0;` — i.e. 0 when the enzyme IS found, and enzyme_info_.size() (a non-zero count) when it is NOT found. The return is an inverted success flag with no documentation, and the enzyme silently defaults to number 0 on an unknown name.
- **Evidence:** Header: `/// sets the enzyme used for cleavage (by means of the number from a list of enzymes)\n    Size setEnzyme(const std::string& enzyme_name);`  Impl (SequestInfile.cpp:542-554): `enzyme_number_ = 0; ... for (...) { if (einfo_i->first == enzyme_name) break; } return (einfo_i == enzyme_info_.end()) ? enzyme_info_.size() : 0;`
- **Fix:** Document the return value explicitly (0 == success / set; non-zero == unknown enzyme), and ideally invert it to a `bool` success or throw on an unknown enzyme. ABI-safe interim fix: add a Doxygen note clarifying the inverted semantics. A cleaner additive fix is a new `bool trySetEnzyme(const std::string&)` returning true on success; deprecate the surprising return.
- **Verifier correction:** setEnzyme returns Size as an UNDOCUMENTED, INVERTED status: 0 == enzyme found/set (success), non-zero (== enzyme_info_.size()) == unknown enzyme (failure) — opposite of the usual truthy==success reading. Correction to the original claim: on an unknown name the enzyme does NOT silently default to number 0; instead enzyme_number_ is left equal to enzyme_info_.size(), an out-of-range index (the loop runs to end() incrementing the counter). Severity is medium (invites backwards `if(setEnzyme(...))` misuse and leaves an out-of-range enzyme index, but is recoverable and the status is at least returned). ABI impact of a documentation/Doxygen clarification is none; the recommended bool-inversion or throw-on-unknown would be source-breaking and should be done additively (e.g. a new bool trySetEnzyme).
- **Verified:** Verified against source. Header (SequestInfile.h:117-118) documents only "sets the enzyme used for cleavage" yet returns Size; the impl (SequestInfile.cpp:542-554) returns `(einfo_i == enzyme_info_.end()) ? enzyme_info_.size() : 0` — i.e. 0 on SUCCESS (enzyme found) and a non-zero count on FAILURE (unknown enzyme). This is an undocumented INVERTED error flag, the opposite of the conventional C++ r

### [FORM-49] SequestInfile::addEnzymeInfo — "addEnzymeInfo" mutates and shrinks the vector the caller passes in
`medium` · `hidden-side-effect` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/SequestInfile.h · _format-id-search-in_

```cpp
void addEnzymeInfo(std::vector<std::string>& enzyme_info)
```
- **Expectation:** A method named `addEnzymeInfo` taking a vector reads as consuming the vector as input to register an enzyme; the caller does not expect their vector to be rewritten and to lose its first element.
- **Actual:** The argument is taken by non-const reference and modified in place: element [2] (the cut-after residue list) is de-duplicated and rewritten, and `enzyme_info.erase(enzyme_info.begin())` removes the first element (the name) before the rest is stored. After the call the caller's vector is altered/shortened.
- **Evidence:** Header: `void addEnzymeInfo(std::vector<std::string> & enzyme_info);`  Impl (SequestInfile.cpp:379-398): `... if (enzyme_info[2].length() != aas.size()) { enzyme_info[2].clear(); ... } std::string enzyme_name = enzyme_info[0]; enzyme_info.erase(enzyme_info.begin());`
- **Fix:** Take the parameter by value or `const&` and work on a local copy so the caller's vector is not clobbered. ABI-break if the signature changes; ABI-safe interim fix is to document that the argument is consumed/modified, or add an overload `addEnzymeInfo(const std::vector<std::string>&)`.
- **Verified:** Verified directly in source. Header SequestInfile.h:199 declares `void addEnzymeInfo(std::vector<std::string>& enzyme_info)` (non-const ref). SequestInfile.cpp:378-406 confirms the hidden side effect: enzyme_info[2] is de-duplicated and rewritten in place when duplicates exist (lines 387-395), then `enzyme_info.erase(enzyme_info.begin())` (line 398) removes the first element (the name) before the 

### [FORM-52] MascotXMLFile::load — Output parameters that load() clears and overwrites are documented as @param[in]
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MascotXMLFile.h · _format-id-search-in_

```cpp
void load(const std::string& filename, ProteinIdentification& protein_identification, PeptideIdentificationList& id_data, const SpectrumMetaDataLookup& lookup)
```
- **Expectation:** A `@param[in]` tag tells the caller the argument is read-only input. A developer might pre-populate `protein_identification`/`id_data` expecting load() to append to or merge with them.
- **Actual:** Both arguments are outputs that are reset/cleared on entry: `protein_identification = ProteinIdentification();` and `id_data.clear();`. Any pre-existing content is silently discarded. The Doxygen marks them `@param[in]`, contradicting the actual data direction (should be `@param[out]`).
- **Evidence:** Header lines 44-46: `@param[in] protein_identification ...` `@param[in] id_data ...`. Impl (MascotXMLFile.cpp:39-40): `protein_identification = ProteinIdentification();  id_data.clear();`
- **Fix:** Change the Doxygen direction tags to `@param[out]` and state that load() clears the outputs first (so callers do not pass data expecting an append/merge). Documentation-only fix; fully ABI-safe.
- **Verified:** Verified directly. Header MascotXMLFile.h lines 44-45 document `@param[in] protein_identification` and `@param[in] id_data`, both passed as non-const references. The 4-arg load() overload (MascotXMLFile.cpp:22-30) forwards to the 5-arg overload, which on entry executes `protein_identification = ProteinIdentification();` (line 39) and `id_data.clear();` (line 40), fully resetting/clearing both. So 

### [FORM-54] InspectInfile::getBlind / getMulticharge — Boolean-sounding "blind"/"multicharge" flags are tri-state UInt with an undocumented 2 = "not set" sentinel
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/FORMAT/InspectInfile.h · _format-id-search-in_

```cpp
UInt getBlind() const; void setBlind(UInt blind); UInt getMulticharge() const; void setMulticharge(UInt)
```
- **Expectation:** The Doxygen describes these as on/off switches ("run Inspect in a blind mode", "If set to true, ..."), so a caller expects a bool with the usual 0/1 meaning.
- **Actual:** They are `UInt` with a third hidden state: the member docs say `0 - false, 1 - true, 2 - not set`. A caller treating the return as a bool (`if (getBlind())`) misreads the "not set" sentinel (2) as true, and there is no enum or named constant exposing this.
- **Evidence:** Header line 88-89: `UInt getBlind() const; void setBlind(UInt blind);` with prose "run Inspect in a blind mode ... If true". Member comments (lines 140-141, 149-150): `/// 0 - false, 1 - true, 2 - not set`.
- **Fix:** Document the tri-state on the public getters/setters (or expose named constants), since callers cannot see the private member comments. Documentation-only fix is ABI-safe; a cleaner additive option is an `isBlindSet()`-style query.
- **Verifier correction:** Claim is accurate; only severity is recalibrated from high to medium. The boolean-sounding UInt flags getBlind()/getMulticharge() are tri-state with an undocumented 2 = "not set" sentinel that is also the default-constructed value, so `if (getBlind())` reads the unset default as truthy ("blind on"). The sentinel is documented solely in private member comments (header lines 141, 150) and confirmed in InspectInfile.cpp (ctor sets 2; store() emits only when != 2). The fix (document tri-state on public getters/setters, or add an isBlindSet()-style query and/or migrate to std::optional<bool>) is ABI-safe/additive, so abi_impact = none rather than the claim's framing of a silently-wrong-results high-severity issue — OpenMS's own output handling is correct; only external callers treating the value as a plain bool are misled.
- **Verified:** Independently confirmed in both header and .cpp. The public getters return UInt getBlind()/getMulticharge() (header lines 88, 116) with Doxygen prose that describes pure on/off switches: "run Inspect in a blind mode ... If true" and "If set to true, attempt to guess the precursor charge ...". The tri-state semantics (0=false, 1=true, 2=not set) are documented ONLY in the private member comments at

### [FORM-58] OMSSACSVFile::load — protein_identification out-parameter is documented as populated but is entirely ignored
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/OMSSACSVFile.h · _format-id-search-out_

```cpp
void load(const std::string & filename, ProteinIdentification & protein_identification, PeptideIdentificationList & id_data) const
```
- **Expectation:** Header documents '@param[in] protein_identification the protein ProteinIdentification data' and the parameter is a non-const reference, so a caller expects load() to fill protein identification data into it.
- **Actual:** The parameter name is commented out in the definition and the body never touches it. OMSSACSVFile.cpp:24: 'void OMSSACSVFile::load(const std::string & filename, ProteinIdentification & /* protein_identification */, PeptideIdentificationList & id_data) const'. Nothing is ever written to it, so callers silently receive whatever they passed in.
- **Evidence:** OMSSACSVFile.cpp:24  void OMSSACSVFile::load(const std::string & filename, ProteinIdentification & /* protein_identification */, PeptideIdentificationList & id_data) const
- **Fix:** Either populate protein_identification (search engine/score type/identifier, as the sibling InspectOutfile/SequestOutfile parsers do) or document explicitly that the parameter is currently unused and reserved. Behavior fix is source-compatible; a doc-only fix is abi none. Marking the param [[maybe_unused]] and updating the doc is the minimal honest fix.
- **Verifier correction:** The non-const reference out-parameter `protein_identification` is silently ignored: its name is commented out in the definition (OMSSACSVFile.cpp:24) and the body never writes to it, so callers (including the pyOpenMS load() wrapper) receive an empty/default ProteinIdentification with no search-engine, score-type, or identifier set — unlike sibling parsers. The claim's framing that it is "documented as populated" is imprecise: the header doc actually tags it `@param[in]` (an input marker), not `@param[out]`; that doc is itself wrong/misleading for a non-const ref but does not literally promise population. Minimal honest fix: mark the param [[maybe_unused]] and correct the doc to state it is currently unused/reserved (doc-only, ABI none), or populate it as the recommendation suggests (source-compatible).
- **Verified:** Code facts confirmed. Header (OMSSACSVFile.h:47) declares a non-const out/in-out reference `ProteinIdentification & protein_identification`; the definition (OMSSACSVFile.cpp:24) comments the name out `/* protein_identification */` and the body (lines 26-99) never writes to it — it only populates id_data. So callers receive whatever ProteinIdentification they passed in (typically empty/default). Th

### [FORM-62] XQuestResultXMLFile::getNumberOfHits / getMinScore / getMaxScore — 'Returns ... in the file' getters actually return cached state from the last load() and return sentinels if load() was never called
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/XQuestResultXMLFile.h · _format-id-search-out_

```cpp
int getNumberOfHits() const; double getMinScore() const; double getMaxScore() const;
```
- **Expectation:** The header docs say 'Returns the total number of hits in the file', 'Returns minimum score among the hits in the file' — phrasing implies these inspect the file/result and are usable standalone.
- **Actual:** They return private members populated only as a side effect of load(). Without a prior load(), getNumberOfHits() returns the constructor sentinel -1 (XQuestResultXMLFile.cpp:24 'n_hits_(-1)'), and getMinScore()/getMaxScore() return uninitialized doubles (min_score_/max_score_ are never initialized in the ctor). Values reflect the most recently loaded file, not 'the file' in any standalone sense.
- **Evidence:** XQuestResultXMLFile.cpp:24  n_hits_(-1)   (min_score_/max_score_ not initialized)\nXQuestResultXMLFile.cpp:37-39  this->n_hits_ = handler.getNumberOfHits(); this->min_score_ = ...; this->max_score_ = ...;\nXQuestResultXMLFile.cpp:56/61/66  return this->n_hits_ / min_score_ / max_score_;
- **Fix:** Document that these return statistics of the LAST loaded file and are only valid after load(); initialize min_score_/max_score_ in the constructor to a defined sentinel (e.g. NaN) so a pre-load call is detectable instead of returning indeterminate memory. Initializing members is abi none.
- **Verifier correction:** The getters return cached state from the most recent load(), not properties of "the file" as the docs imply. n_hits_ returns the loud sentinel -1 before any load (XQuestResultXMLFile.cpp:24), but getMinScore()/getMaxScore() return INDETERMINATE values: min_score_/max_score_ are declared without initializers (header lines 125-126) and never set by the constructor, so a pre-load call is undefined behavior returning garbage rather than a recognizable sentinel. Fix: initialize min_score_/max_score_ in the constructor (e.g. to NaN) and document that all three accessors are valid only after a successful load() and describe the last-loaded file.
- **Verified:** Confirmed against the source. The three getters (XQuestResultXMLFile.cpp:54-67) merely return cached members. The constructor (lines 22-26) initializes ONLY n_hits_(-1); min_score_/max_score_ have no in-class initializer (header lines 125-126) and no ctor initialization, so they are populated solely as a side effect of load() (lines 37-39). The header docs (lines 59,65,71) say "Returns the total n

### [FORM-26] CachedmzML::CachedmzML(const std::string&) — Single-argument constructor is non-explicit, allowing implicit string-to-CachedmzML conversion (with heavy disk I/O)
`medium` · `implicit-conversion` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/CachedMzML.h · _format-mzml_

```cpp
CachedmzML(const std::string& filename)
```
- **Expectation:** A type that performs file loading in its constructor should not implicitly convert from a string; a stray `std::string` passed where a CachedmzML is expected should not silently open and index files.
- **Actual:** The constructor is not marked explicit, so a std::string implicitly converts to a CachedmzML, triggering load_() (metadata parse + side-car index build + opening an ifstream). This can silently fire expensive I/O at unexpected call sites.
- **Evidence:** Header: `CachedmzML(const std::string& filename);` (no `explicit`). Impl: `CachedmzML::CachedmzML(const std::string& filename) { load_(filename); }`.
- **Fix:** Mark the constructor `explicit` (source-compatible for all well-formed call sites that intentionally construct; only breaks accidental implicit conversions, which is the point).
- **Verifier correction:** The constructor is genuinely non-explicit and does perform heavy I/O, confirming the core claim. However, the severity is medium, not high: no current OpenMS API accepts a CachedmzML by value or const-reference, so the implicit string-to-CachedmzML conversion cannot actually fire at any existing call site — the footgun is latent. The strongest argument for the fix is internal inconsistency: all sibling FORMAT file constructors and CachedmzML's own subclass SpectrumAccessOpenMSCached already mark the identical std::string constructor `explicit`. Marking CachedmzML's constructor `explicit` is source-compatible (all existing call sites use direct/base initialization) and ABI-neutral.
- **Verified:** Evidence is verbatim-accurate. Header (src/openms/include/OpenMS/FORMAT/CachedMzML.h:55) declares CachedmzML(const std::string& filename) with no `explicit`; the impl (CachedMzML.cpp:25-28) delegates to load_(), which (lines 44-60) builds the side-car index via CachedMzMLHandler::createMemdumpIndex, opens an ifstream, and parses the full mzML metadata via FileHandler().loadExperiment — genuine hea

### [FORM-28] MzMLFile::transform(filename, consumer, map, ...) — transform() has two adjacent same-typed bool flags whose meaning is invisible at the call site
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzMLFile.h · _format-mzml_

```cpp
void transform(const std::string& filename_in, Interfaces::IMSDataConsumer * consumer, PeakMap& map, bool skip_full_count = false, bool skip_first_pass = false)
```
- **Expectation:** At a call like `f.transform(fn, &c, map, true, true)` a reader cannot tell which bool is skip_full_count and which is skip_first_pass; the defaults make it easy to set the wrong one.
- **Actual:** Both trailing parameters are `bool` with default false and adjacent; swapping them compiles silently and changes behavior (skipping the metadata first pass vs skipping the count pass).
- **Evidence:** Header: `void transform(..., bool skip_full_count = false, bool skip_first_pass = false);` Used internally as `transform(filename, &c, true, true);` in getCentroidInfo — illustrating two positional trues with no naming.
- **Fix:** Consider a strongly-typed flags/enum parameter in a new additive overload; at minimum keep the doc's @param ordering exact. ABI-safe additive change.
- **Verifier correction:** The named symbol (map overload, MzMLFile.h:136) is correctly described, but the quoted internal evidence `transform(filename, &c, true, true)` in getCentroidInfo (MzMLFile.cpp:268) calls the 3-arg overload (line 119), not the map overload. The two-adjacent-defaulted-bool surprise is identical in both overloads; real call sites of the map overload include MzMLFile_test.cpp:1213 and OpenSwathMzMLFileCacher.cpp:207. The boolean-trap concern is strongly corroborated by in-codebase workarounds: inline /*skip_full_count=*/ /*skip_first_pass=*/ annotations at FileConverter.cpp:1032 and a clarifying comment at SwathFile.cpp:189. Severity is medium (boolean trap, recoverable), not high.
- **Verified:** VERIFIED against source. The claimed symbol exists at MzMLFile.h:136 with exactly the signature given: `void transform(const std::string& filename_in, Interfaces::IMSDataConsumer* consumer, PeakMap& map, bool skip_full_count = false, bool skip_first_pass = false)`. Both trailing params are bool, both default false, and adjacent — the classic boolean-trap POLS surprise. The two flags have genuinely

### [FORM-67] MzTabBoolean::get — Boolean getter returns Int and encodes null as -1, so `if (x.get())` treats null as true
`medium` · `return-value` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzTabBase.h · _format-mztab_

```cpp
Int MzTabBoolean::get() const
```
- **Expectation:** A getter on a class named MzTabBoolean returns a bool (or at least a 0/1 value), so `if (row.unique.get())` is safe.
- **Actual:** get() returns Int and simply returns the raw internal `value_`, which is -1 when the cell is null (isNull() == value_ < 0). A caller doing `if (row.unique.get())` gets `true` for a NULL cell because -1 is truthy, and `if (row.unique.get() == 1)` is the only safe form. The boolean-ness is silently lost.
- **Evidence:** MzTabBase.cpp:493-496 `Int MzTabBoolean::get() const { return value_; }`; MzTabBase.cpp:498-501 `bool isNull() const { return value_ < 0; }`. Header MzTabBase.h:168 declares `Int get() const;` inside `class MzTabBoolean`.
- **Fix:** Document loudly that get() may return -1 (null) and callers must check isNull() first, or add a `bool getBool() const` / `std::optional<bool> value() const` accessor. Changing the existing return type would be ABI-breaking; the additive accessor is the safe fix.
- **Verifier correction:** The code facts are exactly as quoted. The behavior (get() returns Int, -1 for null, truthy) is real and undocumented. However, the claimed concrete harm ("a caller doing if (row.unique.get()) gets true for NULL") is hypothetical: no caller anywhere in OpenMS invokes MzTabBoolean::get() in a boolean context, and serialization goes through toCellString() which handles null correctly. The finding is a latent API foot-gun affecting external consumers, mitigated by an adjacent isNull(), and consistent with the whole MzTab* tri-state family — hence medium severity, not high. abi_impact of the existing-behavior finding is none; the recommended additive bool getBool()/std::optional<bool> accessor would be source-compatible, while changing get()'s return type would be breaking.
- **Verified:** Code confirmed exactly as quoted: MzTabBase.h:168 declares `Int get() const;` in class MzTabBoolean; MzTabBase.cpp:493-496 `Int MzTabBoolean::get() const { return value_; }`; MzTabBase.cpp:498-501 `isNull() { return value_ < 0; }`; default ctor and setNull(false) set value_ = -1. So a getter on a class named MzTabBoolean returns Int and yields -1 (truthy) for a null cell, making `if (x.get())` wro

### [FORM-68] MzTabSpectraRef::setSpecRefFile — setSpecRefFile is a byte-for-byte duplicate of setSpecRef (does not set any file)
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzTabBase.h · _format-mztab_

```cpp
void MzTabSpectraRef::setSpecRefFile(const std::string& spec_ref)
```
- **Expectation:** Given the sibling setMSFile(Size) sets the ms_run index, a method named setSpecRefFile sounds like it sets the file/run associated with the spectrum reference (i.e. something distinct from setSpecRef, which sets the textual native spectrum reference).
- **Actual:** setSpecRefFile writes to the exact same member `spec_ref_` with the exact same body as setSpecRef; the only difference is setSpecRef also logs a warning on empty input. There is no 'file' semantics at all - the name promises a different field than it touches.
- **Evidence:** MzTabBase.cpp:190-201 setSpecRef sets `spec_ref_ = spec_ref;`; MzTabBase.cpp:215-222 setSpecRefFile sets the identical `spec_ref_ = spec_ref;`. The 'File' in the name maps to no separate member (the file/run is `ms_run_`, set only by setMSFile).
- **Fix:** Deprecate setSpecRefFile (it is a confusing alias) and point callers to setSpecRef; or, if a file-setting overload was intended, implement it to set ms_run_/a path. Removal is source-breaking; mark `[[deprecated]]` first.
- **Verifier correction:** setSpecRefFile is a never-called, undocumented public alias of setSpecRef: both set spec_ref_ (the textual native spectrum reference), not any file/run. The file/run is ms_run_, set only by setMSFile(Size). The "File" in the name maps to no separate member. It is not currently corrupting data (no callers), but it invites a wrong-field bug. Recommend marking [[deprecated]] (source-compatible) pointing to setSpecRef, or removing it; full removal would be source-breaking.
- **Verified:** Verified directly. setSpecRefFile (MzTabBase.cpp:215-222) has a body identical to setSpecRef (MzTabBase.cpp:190-201): both write `spec_ref_ = spec_ref;` to the textual native spectrum reference. The only difference is setSpecRef logs a warning on empty input. The parameter is even named `spec_ref`, and nothing file/run-related (ms_run_, set only by setMSFile(Size)) is touched. The header (line 308

### [FORM-74] MzTabModification::getModOrSubstIdentifier — Getter asserts !isNull(), so it can abort (debug) / mislead (release) on a null modification
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzTab.h · _format-mztab_

```cpp
MzTabString MzTabModification::getModOrSubstIdentifier() const
```
- **Expectation:** A plain const getter named get...Identifier() returns the stored identifier (possibly a null MzTabString) without precondition surprises.
- **Actual:** It begins with `assert(!isNull())`. On a null modification this aborts the process in debug builds; in release the assert is stripped and it just returns the (null) identifier. The getter therefore has an undocumented precondition and inconsistent debug/release behavior. The setter pair is also asymmetric: setModificationIdentifier exists but the getter is named getModOrSubstIdentifier.
- **Evidence:** MzTab.cpp:65-69 `MzTabString getModOrSubstIdentifier() const { assert(!isNull()); return mod_identifier_; }`; setter at MzTab.cpp:60-63 is named setModificationIdentifier. Same assert-in-getter pattern in MzTabSpectraRef::getSpecRef (MzTabBase.cpp:203-207) and getMSFile (209-213).
- **Fix:** Drop the assert from these getters (let callers check isNull()) or document the precondition; align the get/set name pair (getModificationIdentifier vs setModificationIdentifier). Removing the assert is a behavior fix with no ABI change; renaming would need a deprecated alias.
- **Verifier correction:** The getter does not throw; it calls assert(!isNull()), which abort()s only in debug builds and is compiled out in release. So the surprise is: an undocumented precondition with inconsistent debug-abort vs release-silent-null behavior, not a thrown exception. Severity is medium (loud in debug, recoverable null in release; no silent data corruption), not high. Fix: drop or document the assert (body-only change, ABI none); separately, the get/set name pair is misaligned (getModOrSubstIdentifier vs setModificationIdentifier) — renaming for symmetry would require a deprecated alias, which is source-compatible.
- **Verified:** Evidence verified verbatim. MzTab.cpp:65-69 defines `MzTabString MzTabModification::getModOrSubstIdentifier() const { assert(!isNull()); return mod_identifier_; }`. The setter at MzTab.cpp:60-63 is named setModificationIdentifier (asymmetric get/set naming). The identical assert-in-getter pattern exists in MzTabSpectraRef::getSpecRef (MzTabBase.cpp:203-207) and getMSFile (209-213). The header decl

### [FORM-87] OMSFile::load — load() appends/merges into the passed object instead of constructing/replacing it
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/OMSFile.h · _format-oms-sqlite_

```cpp
void load(const std::string& filename, IdentificationData& id_data); void load(const std::string& filename, FeatureMap& features); void load(const std::string& filename, ConsensusMap& consensus);
```
- **Expectation:** The Doxygen says 'Read in an OMS file and construct an IdentificationData object / construct a feature map'. A caller reading 'construct' reasonably expects the output object to be replaced (cleared then filled), so calling load() twice, or on a non-empty object, yields exactly the file's contents.
- **Actual:** load() never clears the target. OMSFileLoad::load(FeatureMap&) calls features.push_back(feature) for every feature (OMSFileLoad.cpp:994) and load(IdentificationData&) merges into the existing IdentificationData (OMSFileLoad.cpp:758-786, no clear). Pre-existing features/IDs in the passed object are kept and the loaded data is appended/merged on top.
- **Evidence:** OMSFileLoad.cpp:994 `features.push_back(feature);`; OMSFileLoad.cpp:1035 `consensus.push_back(feature);`; OMSFileLoad.cpp:758 load(IdentificationData&) just calls loadInputFiles_/loadScoreTypes_/... with no clear; header doc OMSFile.h:55-60 'Read in an OMS file and construct an IdentificationData object'.
- **Fix:** Document the append/merge semantics explicitly in the header (replace 'construct' with 'load into / append to'), or clear the target at the top of each load() for replace semantics. If round-trip-into-fresh-object is the only tested path, an additive doc fix is ABI-safe; clearing would change behavior for callers who currently rely on merge.
- **Verifier correction:** Severity adjusted from high to medium. The append/merge-without-clear behavior and "construct" doc mismatch are real and silent, but: (1) every shipped call site passes a freshly-constructed object (IDFileConverter.cpp:682, NucleicAcidSearchEngine.cpp:1066, MapAlignerIdentification.cpp:275, and all three OMSFile_test.cpp load sections), so there is no live bug today; (2) the only path to wrong results is the atypical case of loading into a non-empty/reused object, which yields duplicated/merged data but is additive and recoverable, not data-loss or crash. It "invites misuse and produces silently wrong results" (medium) rather than "silently wrong for typical/reasonable use" (high). The recommendation stands: either fix the doc to say "load into / append to" (additive, ABI-safe) or clear the target at the top of each load() to match FeatureXMLFile's replace semantics (behavior change, not ABI). abi_impact is none: signatures are unchanged and the only currently-tested path (round-trip into a fresh object) is unaffected either way.
- **Verified:** Independently verified the code. OMSFile.h:55-74 documents load() as "construct an IdentificationData object / construct a feature map / construct a consensus map" with @param[out]-style intent. The implementations in OMSFileLoad.cpp do NOT clear the target: load(IdentificationData&) (lines 758-786) only calls loadInputFiles_/loadScoreTypes_/... which register* into existing collections (e.g. load

### [FORM-89] OSWFile::getRunID — const getter opens a second connection in READWRITE_OR_CREATE mode, which can create a file
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/OSWFile.h · _format-oms-sqlite_

```cpp
UInt64 getRunID() const;
```
- **Expectation:** getRunID() const on an OSWFile whose primary connection was opened READ_ONLY should be a pure read and never create or write a file.
- **Actual:** getRunID() constructs a fresh SqliteConnector(filename_) using the default SqlOpenMode::READWRITE_OR_CREATE (SqliteConnector.h:50). Per openDatabase_, that adds SQLITE_OPEN_CREATE, so on a path that doesn't exist a 0-byte SQLite DB is created instead of failing fast, and an existing file is opened for writing rather than read-only.
- **Evidence:** OSWFile.cpp:899 `SqliteConnector conn(filename_);` (no mode -> default); SqliteConnector.h:50 default `SqlOpenMode mode = SqlOpenMode::READWRITE_OR_CREATE`; SqliteConnector.cpp:45-46 `case READWRITE_OR_CREATE: flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE;`. Contrast the OSWFile member conn_ which is opened READ_ONLY (OSWFile.cpp:527).
- **Fix:** Pass SqlOpenMode::READ_ONLY explicitly in getRunID() (and other read-only helpers) so a missing file errors instead of being silently created. Implementation-only change, ABI-safe.
- **Verifier correction:** getRunID() const opens a second SqliteConnector in the default READWRITE_OR_CREATE mode rather than READ_ONLY. Because the OSWFile constructor already opens conn_ READ_ONLY and requires the file to exist, the "silently creates a missing file" outcome is effectively unreachable outside a TOCTOU race; the real, reachable surprise is that a pure const read requests a write-capable open (SQLITE_OPEN_READWRITE), so the operation fails or is blocked on read-only files/media that the READ_ONLY constructor accepted. Fix: pass SqlOpenMode::READ_ONLY in getRunID() (and the other read-only helpers such as readRunBasenames() at line 868 which share the same bare default). Implementation-only; signature unchanged, ABI-safe.
- **Verified:** Evidence is accurate and verified independently. OSWFile.cpp:899 constructs `SqliteConnector conn(filename_);` with no mode; SqliteConnector.h:50 defaults that to SqlOpenMode::READWRITE_OR_CREATE, which SqliteConnector.cpp:45-46 maps to SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE. By contrast the class member conn_ is opened READ_ONLY (OSWFile.cpp:527). So a const getter documented as "extract the 

### [FORM-98] ParamXMLFile::load — load() merges into the existing Param instead of clearing it first
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/ParamXMLFile.h · _format-param_

```cpp
void load(const std::string& filename, Param& param);
```
- **Expectation:** A caller passing an existing Param to load() reasonably expects load to replace its contents with the file's parameters (the usual *File::load contract), so that reusing the same Param object across two load() calls yields just the second file's parameters.
- **Actual:** ParamXMLHandler is constructed with `param_(param)` (a reference to the caller's object) and only ever calls param_.setValue()/param_.addSection(); it never clears param_. So load() merges the file into whatever the Param already contained. Loading file B into a Param that previously held file A leaves A's non-overlapping entries in place, with values silently overwritten only where keys collide.
- **Evidence:** ParamXMLFile::load: `Internal::ParamXMLHandler handler(param, filename, schema_version_); parse_(filename, &handler);` — handler ctor: `XMLHandler(filename, version), param_(param) {}` (no clear). All mutations are `param_.setValue(...)` / `param_.addSection(...)`; grep shows no `param_.clear()` anywhere in ParamXMLHandler.cpp.
- **Fix:** Document explicitly in the header that load() merges into (does not reset) the passed Param, or add `param.clear();` at the start of load() to match the conventional clear-then-fill semantics. If changing behavior is too risky for ABI, at minimum fix the doc; ideally add an overload or a bool clear-first parameter (defaulted to preserve current behavior).
- **Verified:** Evidence verified independently. ParamXMLFile::load (ParamXMLFile.cpp:357-361) constructs Internal::ParamXMLHandler(param, ...) which stores a reference param_(param) (ParamXMLHandler.cpp:19-23) and only ever calls param_.setValue()/addSection()/setValidStrings()/setMin*/setMax* — there is no param_.clear() anywhere in the handler. So load() merges the file into whatever the Param already held; re

### [FORM-99] ParamJSONFile::load — load() documented to return false on an unknown parameter, but actually throws ParseError (the false path is dead)
`medium` · `contract-mismatch` · ABI: `none` · src/openms/include/OpenMS/FORMAT/ParamJSONFile.h · _format-param_

```cpp
static bool load(const std::string& filename, Param& param);
```
- **Expectation:** The header states: '@return returns true if file was successfully loaded; false if an unknown (non-default) parameter name was encountered in the JSON file.' A caller therefore writes `if (!ParamJSONFile::load(f, p)) { handle unknown-param }` expecting a false return, not an exception.
- **Actual:** When an unknown key is found, the implementation throws `Exception::ParseError` (ParamJSONFile.cpp line 89), it does not return false. The only `return false;` (line 165) sits immediately after a `throw` inside the json::exception catch block and is unreachable. The function in practice returns only `true` or throws — the documented false-return contract can never occur.
- **Evidence:** Header: '@return ... false if an unknown (non-default) parameter name was encountered'. Impl: `if (!param.exists(key)) { ... throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", msg); }` and later `throw Exception::ParseError(...); return false;` (return after throw = dead code).
- **Fix:** Reconcile the contract: either return false (no throw) on an unknown parameter as documented, or change the doc to state that an unknown parameter throws ParseError and that the function returns true on success / throws otherwise. Pure doc+impl-return fix; no signature change, so ABI-neutral.
- **Verifier correction:** The header @return documents a false return for an unknown parameter, but the implementation throws Exception::ParseError (ParamJSONFile.cpp:86-90) and never returns false; the lone `return false;` (line 165) is unreachable dead code after a throw (line 164). The function returns only `true` (line 167) or throws. This is a documented-return-that-never-occurs / dead-code contract mismatch, not a "silent" failure (the rejection is loud — it throws). It is genuinely observable: the two real callers in TOPPBase.cpp (lines 274, 354) use `if (!ParamJSONFile::load(...)) { return ILLEGAL_PARAMETERS; }`, so their false-branches are dead; the thrown ParseError is instead caught at TOPPBase.cpp:524 and returned as INPUT_FILE_CORRUPT rather than the intended ILLEGAL_PARAMETERS. Severity medium (loud and recoverable, wrong exit code / dead branches, but no silent wrong results or data loss). ABI impact none (doc + impl-return fix, no signature change).
- **Verified:** Independently verified against the actual code. Header ParamJSONFile.h line 71 documents: "@return returns true if file was successfully loaded; false if an unknown (non-default) parameter name was encountered in the JSON file". The implementation (ParamJSONFile.cpp) contradicts this: on an unknown key, lines 86-90 do `if (!param.exists(key)) { ... throw Exception::ParseError(...); }` — it throws,

### [FORM-103] ParamJSONFile::load — load() dereferences param.begin().getTrace().front() before any file/empty check — undefined behavior on an empty Param
`medium` · `surprising-undefined-behavior (crash before documented exception)` · ABI: `none` · src/openms/source/FORMAT/ParamJSONFile.cpp · _format-param_

```cpp
static bool load(const std::string& filename, Param& param);
```
- **Expectation:** load(filename, param) should first validate the file (FileNotFound etc., which the header documents) and operate safely regardless of whether the passed Param already has entries. A caller passing a freshly default-constructed (empty) Param expects a clean error, not a crash.
- **Actual:** The very first statements call `param.begin().getTrace()` and then `traces.front().name` to derive the tool namespace, before the file-existence checks. If param is empty, `traces` can be empty and `.front()` is undefined behavior (crash), and this happens even before the documented FileNotFound check that the header advertises.
- **Evidence:** ParamJSONFile.cpp lines 42-43: `auto traces = param.begin().getTrace(); std::string toolNamespace = traces.front().name + ":1:";` — executed before the `if (!ifs.good())` / File::exists checks at lines 45-60.
- **Fix:** Guard against an empty trace (throw a clear InvalidValue/Precondition exception) and move the file-existence validation ahead of namespace derivation so the documented FileNotFound path is honored. Implementation-only fix, ABI-neutral.
- **Verifier correction:** load() does not 'throw' on an empty Param — it invokes std::vector::front() on an empty trace vector (UB / crash) at lines 42-43, before any file-existence check. This violates the header's documented @exception FileNotFound ordering when the Param is empty. However, the header documents `param` as having 'pre-filled defaults', and all real callers pass a populated Param, so the documented/reasonable path is unaffected; the defect is silent UB on a precondition violation (empty/default-constructed Param) rather than a clean error. Fix is to validate the file first and guard against an empty trace; implementation-only, ABI-neutral.
- **Verified:** Verified against the actual code. ParamJSONFile.cpp lines 42-43 are exactly as quoted: `auto traces = param.begin().getTrace();` then `traces.front().name + ":1:"`, executed BEFORE the ifstream/File::exists checks at lines 45-60. I traced the empty-Param path: Param::begin() -> ParamIterator(root_); in Param.cpp the ctor (line 1473) detects an empty root (no entries, no nodes), sets root_=nullptr 

### [FORM-32] DTAFile::store / DTAFile::load — load/store use different proton masses, so a load->store->load round-trip is lossy
`medium` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/DTAFile.h · _format-peakfiles_

```cpp
void store(const std::string&, const SpectrumType&) const  /  void load(const std::string&, SpectrumType&)
```
- **Expectation:** load() and store() are an inverse pair: storing what load() produced and reloading should recover the same precursor MH+/charge values.
- **Actual:** load (line 120) computes m/z with the accurate constant `Constants::PROTON_MASS_U` (1.0072764667710), but store (line 204) reconstructs MH+ with a hardcoded `1.0`: `os << ((precursor.getMZ() - 1.0) * precursor.getCharge() + 1.0)`. The two are not inverse, so the stored MH+ is off and a re-load shifts the precursor m/z by ~0.0073*charge per round-trip.
- **Evidence:** Load: `precursor.setMZ((mh_mass - Constants::PROTON_MASS_U) / charge + Constants::PROTON_MASS_U);` (DTAFile.h:120). Store: `os << ((precursor.getMZ() - 1.0) * precursor.getCharge() + 1.0);` (DTAFile.h:204).
- **Fix:** Source/ABI-compatible: change the two literal `1.0` values in store() to `Constants::PROTON_MASS_U` to match load(). Header-only inline template change; recompiles dependents but no API signature change.
- **Verifier correction:** store() uses a hardcoded proton mass of 1.0 while load() uses Constants::PROTON_MASS_U (1.0072764667710), making them non-inverse for multiply-charged precursors. The round-trip m/z error is (PROTON_MASS_U - 1.0)*(charge-1)/charge = 0.0072765*(z-1)/z Th, independent of m/z: 0 Th at z=1, ~0.0049 Th at z=3, ~0.0064 Th at z=8 — NOT '~0.0073*charge per round-trip' as stated. Fix: replace both literal 1.0 values in store() (line 204) with Constants::PROTON_MASS_U to match load(). Severity is medium rather than high: the error is small (sub-0.007 Th), silent (no crash/warning), zero for singly-charged spectra, confined to a legacy single-spectrum format, and recoverable; but it does silently corrupt multiply-charged precursor m/z on round-trip and currently slips past the existing unit test only because the test value's relative error falls just under the 1e-5 TEST_REAL_SIMILAR tolerance.
- **Verified:** Independently confirmed in DTAFile.h. load() (line 120) reconstructs precursor m/z using the accurate constant Constants::PROTON_MASS_U (1.0072764667710): `precursor.setMZ((mh_mass - Constants::PROTON_MASS_U) / charge + Constants::PROTON_MASS_U);`. store() (line 204) reconstructs the written MH+ using a hardcoded integer proton mass 1.0: `os << ((precursor.getMZ() - 1.0) * precursor.getCharge() + 

### [FORM-75] MRMFile::loadMzML — Documented 'tmp' (temporary directory for cached data) parameter is silently ignored
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MRMFile.h · _format-quant_

```cpp
static std::vector<::OpenSwath::SwathMap> loadMzML(const std::string& file, const std::string& tmp, std::shared_ptr<ExperimentalSettings>& exp_meta)
```
- **Expectation:** The header documents '@param[in] tmp Temporary directory (for cached data)'. A caller passing a tmp dir expects caching/spill behavior to honor it (e.g. for large SRM files), as is the convention for cached-load APIs elsewhere in OpenMS (SwathFile).
- **Actual:** The implementation drops the argument entirely: the parameter is commented out as 'const std::string& /*tmp*/' and never referenced. loadMzML always loads the full experiment into memory via FileHandler regardless of what is passed.
- **Evidence:** MRMFile.cpp:23-25: 'std::vector<::OpenSwath::SwathMap> MRMFile::loadMzML(const std::string& file, const std::string& /*tmp*/, std::shared_ptr<ExperimentalSettings>& exp_meta)'. Body has no use of tmp.
- **Fix:** Either honor the parameter (cache to tmp like SwathFile does) or update the doc to state it is currently unused/reserved. If it is genuinely dead, an additive deprecation comment is the ABI-safe step; removing the parameter is source-breaking, so keep the signature but document the no-op.
- **Verifier correction:** The parameter is genuinely ignored and the doc is misleading, but the impact is memory/scalability (always loads in-memory) rather than incorrect output, data loss, or crash — hence medium, not high. The ABI-safe fix is to keep the signature and correct the Doxygen comment to state the tmp directory is currently unused/reserved (or, better, implement caching to match SwathFile). Removing the parameter would be source-breaking; documenting the no-op is ABI-neutral (none).
- **Verified:** Evidence is verbatim accurate. MRMFile.h:27 documents '@param[in] tmp Temporary directory (for cached data).', while MRMFile.cpp:23-25 comments the parameter out as 'const std::string& /*tmp*/' and never references it; the body unconditionally loads the full experiment into memory via FileHandler::loadExperiment (line 31). This is a genuine surprise, not a domain convention or idiom: the conventio

### [FORM-76] MRMFeatureQCFile::store — store() silently writes no file at all when the selected QC vector is empty
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MRMFeatureQCFile.h · _format-quant_

```cpp
void store(const std::string& filename, const MRMFeatureQC& mrmfqc, const bool is_component_group)
```
- **Expectation:** A store(filename, ...) call is expected to produce an output file at 'filename' (at minimum a header-only CSV). Callers commonly rely on the file existing afterwards (e.g. to feed a downstream tool).
- **Actual:** When the relevant vector is empty, store() returns immediately without creating or truncating any file and without raising the documented Exception::UnableToCreateFile. A stale file at that path is left untouched; a missing file stays missing. No signal to the caller.
- **Evidence:** MRMFeatureQCFile.cpp:192 'if (mrmfqc.component_group_qcs.empty()) return;' and :250 'if (mrmfqc.component_qcs.empty()) return;' before any CsvFile::store(filename) call.
- **Fix:** Write at least the header row (current code already builds headers from the first element, so guard only that part), or document the early-return contract explicitly. ABI-safe: behavior/doc change only.
- **Verifier correction:** store() returns silently without writing or truncating any file (and without the documented UnableToCreateFile) when the selected QC vector is empty; the header documents only the UnableToCreateFile exception, so this no-op is undocumented. Real but medium severity (missing/stale file, surfaced loudly on downstream load — not silently-wrong data). Note the empty-guards also prevent the subsequent .at(0) from throwing on empty input, so a correct fix must still build and write the header row rather than simply remove the guard.
- **Verified:** Evidence is verified verbatim: MRMFeatureQCFile.cpp:192 `if (mrmfqc.component_group_qcs.empty()) return;` and :250 `if (mrmfqc.component_qcs.empty()) return;` both sit before the only `CsvFile::store(filename)` calls (lines 246, 285). CsvFile::store -> TextFile::store writes unconditionally, so the non-empty path always creates a file, but the empty path returns first: no file written, no stale fi

### [FORM-77] MRMFeatureQCFile::load — Bare bool is_component_group selects which half of MRMFeatureQC is (re)loaded; a single load leaves the other half untouched
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MRMFeatureQCFile.h · _format-quant_

```cpp
void load(const std::string& filename, MRMFeatureQC& mrmfqc, const bool is_component_group) const
```
- **Expectation:** load(filename, mrmfqc) reads 'the file' into the object. A reader of load(file, qc, true) cannot tell what 'true' means at the call site, and would expect one call to fully populate mrmfqc.
- **Actual:** is_component_group toggles between two entirely different parse paths. load() clears and fills ONLY component_qcs (false) or ONLY component_group_qcs (true); the other vector is left as-is. Fully populating an MRMFeatureQC requires two load() calls on two different files with opposite flags - non-obvious from the signature.
- **Evidence:** MRMFeatureQCFile.cpp:31-48: branch on is_component_group clears+fills only one of mrmfqc.component_qcs / mrmfqc.component_group_qcs.
- **Fix:** Prefer named methods loadComponentQCs()/loadComponentGroupQCs() (additive overloads, ABI-safe) or an enum instead of bool. At minimum document that each call populates only one half. Splitting is additive; removing the bool overload is source-breaking.
- **Verifier correction:** is_component_group toggles between two of THREE parse targets, not two halves: false clears+fills only mrmfqc.component_qcs; true clears+fills only mrmfqc.component_group_qcs; a third vector, component_group_pair_qcs, is populated by neither call. A single load() therefore populates exactly one vector and leaves the rest untouched; fully populating component- and group-level QCs requires two load() calls on two separate files with opposite flags, which is non-obvious from the bare-bool signature and undocumented in the method's Doxygen.
- **Verified:** Code confirmed at MRMFeatureQCFile.cpp:31-48. load() branches on is_component_group: false clears+fills ONLY mrmfqc.component_qcs (lines 33-38); true clears+fills ONLY mrmfqc.component_group_qcs (lines 42-47). The non-selected vector is left untouched. MRMFeatureQCFile_test.cpp:105-106 and 258-259 confirm real usage requires TWO load() calls on TWO different files (MRMFeatureQCFile_1.csv with fals

### [FORM-78] AbsoluteQuantitationMethodFile::load — load() only WARNS (does not throw) when required columns are missing, then silently emits zero/empty-filled records
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h · _format-quant_

```cpp
void load(const std::string & filename, std::vector<AbsoluteQuantitationMethod> & aqm_list)
```
- **Expectation:** The header documents '@exception Exception::ParseError is thrown if an error occurs during parsing'. A caller expects a malformed/incomplete file (e.g. missing llod/uloq/transformation_model columns) to be flagged.
- **Actual:** Missing required columns produce only an OPENMS_LOG_WARN and parsing continues; parseLine_ substitutes 0 / empty string for every absent field. The returned aqm_list looks fully populated but silently carries default values, with no ParseError raised.
- **Evidence:** AbsoluteQuantitationMethodFile.cpp:31-57 emits only OPENMS_LOG_WARN on missing headers; parseLine_ (lines 79-93) defaults missing fields to ""/0 instead of erroring.
- **Fix:** Throw Exception::ParseError (as documented) when required columns are absent, or downgrade the header's exception contract to match the warn-and-default behavior. Behavior/doc change, ABI-safe.
- **Verifier correction:** The behavior is correctly described, but it is not fully silent: missing required columns trigger an OPENMS_LOG_WARN (cpp:45-56) before parsing continues with 0/empty defaults. The documented Exception::ParseError contract is genuinely violated for the missing-column case, but because a warning is logged the failure is loud-and-recoverable rather than silent — severity medium, not high.
- **Verified:** Verified against source. Header (AbsoluteQuantitationMethodFile.h:36-37) documents `@exception Exception::ParseError is thrown if an error occurs during parsing`. In the .cpp, load() (lines 31-57) only emits OPENMS_LOG_WARN when any of the 11 required columns is missing and then continues; it never throws ParseError (grep confirms the file contains no ParseError throw, and CsvFile::load only raise

### [FORM-79] AbsoluteQuantitationMethodFile::parseLine_ (via load) — Numeric parsing uses std::stod/std::stoi, so a non-empty garbage cell throws std::invalid_argument, not the documented Exception::ParseError
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h · _format-quant_

```cpp
void load(const std::string & filename, std::vector<AbsoluteQuantitationMethod> & aqm_list)
```
- **Expectation:** The header promises Exception::ParseError on parse failure. Callers catching OpenMS exceptions expect a malformed numeric cell (e.g. llod="abc") to surface as ParseError.
- **Actual:** parseLine_ calls std::stod/std::stoi directly on non-empty cells with no try/catch, so a malformed numeric value escapes as a raw std::invalid_argument/std::out_of_range, bypassing the documented OpenMS exception type.
- **Evidence:** AbsoluteQuantitationMethodFile.cpp:82-92 'std::stod(tl[headers.at("llod")])' etc.; :87 'std::stoi(...)'. No conversion-error handling. setCastValue_ (152-170) likewise uses raw std::stod/std::stoi.
- **Fix:** Wrap conversions and rethrow as Exception::ParseError to honor the documented contract (matches OpenMS convention). ABI-safe (implementation-only).
- **Verifier correction:** Severity is medium rather than high: the failure is loud (an exception still propagates and aborts the load — no silent wrong results or data loss). The real defect is the exception TYPE: std::invalid_argument/std::out_of_range escape instead of the documented Exception::ParseError and are not caught by catch(Exception::BaseException&), so a malformed method CSV can turn into an uncaught exception (std::terminate/crash) for callers that only guard against OpenMS exceptions. Fix is implementation-only (route conversions through String::toDouble/toInt or wrap stod/stoi and rethrow as Exception::ParseError); signatures unchanged, so ABI impact is none.
- **Verified:** Verified against source. Header (AbsoluteQuantitationMethodFile.h:36-42) documents that load() throws Exception::ParseError "if an error occurs during parsing". The implementation (AbsoluteQuantitationMethodFile.cpp) parses numeric cells with raw std::stod (lines 82-85, 91) and std::stoi (line 87), and setCastValue_ does the same (std::stod line 160, std::stoi line 164). There is no try/catch anyw

### [FORM-80] MSstatsFile::storeLFQ — Bare bool 'is_isotope_label_type' is unreadable at call sites and is a boolean despite its noun-like name
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MSstatsFile.h · _format-quant_

```cpp
void storeLFQ(const std::string& filename, const ConsensusMap&, const ExperimentalDesign&, const StringList& reannotate_filenames, const bool is_isotope_label_type, const std::string& bioreplicate, const std::string& condition, const std::string& retention_time_summarization_method, const bool remove_shared_peptides = true)
```
- **Expectation:** A parameter named is_isotope_label_type reads like it asks 'does this have an isotope label type?'. At a call site, storeLFQ(..., true, "rep", "cond", ...) gives no hint of meaning.
- **Actual:** Per the doc it actually encodes a binary choice of label string: 'If true, IsotopeLabelType is H, else L'. It is a value selector, not a presence predicate, and sits adjacent to several same-typed std::string args raising swap risk.
- **Evidence:** MSstatsFile.h:45 '@param[in] is_isotope_label_type If true, IsotopeLabelType is H, else L' and :55 'const bool is_isotope_label_type'.
- **Fix:** Prefer an enum {Light, Heavy} (additive overload to keep ABI) or rename to 'use_heavy_label'. The default-less bool wedged between string params is the swap hazard; an overload is ABI-safe, signature change is source-breaking.
- **Verified:** Evidence verified in source. Header MSstatsFile.h:45 documents '@param[in] is_isotope_label_type If true, IsotopeLabelType is 'H', else 'L'' and :55 declares 'const bool is_isotope_label_type'. The .cpp (MSstatsFile.cpp:307-313) confirms the bool is a VALUE SELECTOR, not a presence predicate: it picks the literal label string "H" vs "L" (default "L"), and even carries a '@todo remove? ... I think 

### [FORM-81] MRMFeaturePickerFile (public CsvFile inheritance) — MRMFeaturePickerFile publicly inherits CsvFile (exposing load/store/addRow), inconsistent with sibling File adapters that inherit privately
`medium` · `asymmetric-api` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/MRMFeaturePickerFile.h · _format-quant_

```cpp
class MRMFeaturePickerFile : public CsvFile
```
- **Expectation:** Sibling adapters in the same cluster (AbsoluteQuantitationMethodFile, MRMFeatureQCFile) inherit CsvFile privately so only their own load()/store() are visible. A reader expects MRMFeaturePickerFile to behave the same.
- **Actual:** MRMFeaturePickerFile inherits CsvFile publicly, so a MRMFeaturePickerFile object also exposes the entire CsvFile API (CsvFile::load(file, sep, ...), store(), addRow(), getRow()). A caller can invoke the raw CsvFile::load with a different signature than the class's own load(), and mutate the internal row buffer - surprising for a typed parameter-set loader.
- **Evidence:** MRMFeaturePickerFile.h:46-48 'class OPENMS_DLLAPI MRMFeaturePickerFile : public CsvFile'. Contrast AbsoluteQuantitationMethodFile.h:24-25 'private CsvFile' and MRMFeatureQCFile.h:25-26 'private CsvFile'.
- **Fix:** Change to 'private CsvFile' for consistency and encapsulation. This is source-breaking for any code relying on the exposed CsvFile API, so flag as the ideal fix while acknowledging ABI/source impact; an additive doc note is the conservative interim step.
- **Verifier correction:** MRMFeaturePickerFile publicly inherits CsvFile (MRMFeaturePickerFile.h:46-47), leaking the full CsvFile public API (load/store/addRow/getRow/clear/rowCount), whereas its sibling FORMAT adapters AbsoluteQuantitationMethodFile and MRMFeatureQCFile inherit CsvFile privately. This is an asymmetric-encapsulation surprise. The practical risk is misuse (calling the leaked base load/store/addRow yields empty/raw-buffer results instead of typed parameter sets), but it is loud/recoverable rather than silent corruption, because the class's own load(filename, cp_list, cgp_list) is a distinct 3-arg overload that is not shadowed; hence medium, not high. Recommended fix (switch to `private CsvFile`) is source-breaking, including for the pyOpenMS binding which currently declares CsvFile as the public base (bind_misc.cpp:4240).
- **Verified:** Evidence verified against source. MRMFeaturePickerFile.h:46-47 declares `class OPENMS_DLLAPI MRMFeaturePickerFile : public CsvFile`, while the two sibling FORMAT adapters by the same authors inherit privately: AbsoluteQuantitationMethodFile.h:24-25 `private CsvFile` and MRMFeatureQCFile.h:25-26 `private CsvFile`. CsvFile (CsvFile.h:23-96) publicly exposes load(filename,sep,ie,first_n), store, addR

### [FORM-82] SwathFile::loadSqMass — loadSqMass takes an exp_meta out-parameter but never writes to it
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/SwathFile.h · _format-targeted_

```cpp
std::vector<OpenSwath::SwathMap> loadSqMass(const std::string& file, std::shared_ptr<ExperimentalSettings>& /* exp_meta */)
```
- **Expectation:** Like every sibling SwathFile loader (loadMzML, loadMzXML, loadSplit, loadFromMSExperiment, loadBrukerTdf), the caller-supplied exp_meta shared_ptr is populated with the file's ExperimentalSettings so downstream code (e.g. SwathQC) has metadata.
- **Actual:** The parameter name is commented out both in the header declaration and in the definition (line 310 of SwathFile.cpp), and the body never assigns exp_meta. The caller's shared_ptr is left exactly as passed in (typically null or stale), with no diagnostic.
- **Evidence:** Header: `std::vector<OpenSwath::SwathMap> loadSqMass(const std::string& file, std::shared_ptr<ExperimentalSettings>& /* exp_meta */);`. SwathFile.cpp line 310: `SwathFile::loadSqMass(const std::string& file, std::shared_ptr<ExperimentalSettings>& /* exp_meta */)` — the body (lines 312-337) returns swath_maps and never touches exp_meta. Compare loadMzML line 140 `exp_meta = exp_stripped;`.
- **Fix:** Either populate exp_meta from the sqMass RUN/metadata (additive, source- and ABI-compatible — same signature) so behavior matches siblings, or, if metadata genuinely cannot be recovered, drop the parameter from the public signature (ABI-breaking) and document the absence. At minimum, set exp_meta to a non-null empty ExperimentalSettings so callers do not dereference null.
- **Verifier correction:** The signature and the never-written exp_meta are confirmed exactly as quoted. The functional impact is correctly identified (out-param silently never populated, unlike all siblings), but the claim's recommended "set exp_meta to non-null so callers do not dereference null" is partly misdirected: the actual in-tree caller (OpenSwathWorkflow.cpp:1079) already default-constructs a non-null ExperimentalSettings, so there is no null-dereference risk — the real defect is that for sqMass inputs the experimental metadata (source files, instrument, sample) is silently empty in downstream output (prepareChromOutput, OpenSwathBase.cpp:372), inconsistent with mzML/mzXML/Bruker/Thermo paths. Severity downgraded high→medium (silent provenance/metadata loss, recoverable, no quant-result corruption or crash). The recommended additive fix (populate exp_meta from the sqMass RUN/metadata, keeping the same signature) is source- and ABI-compatible: abi_impact none.
- **Verified:** Independently verified against the actual code. The evidence is exact: header line 100 and definition line 310 both comment out the parameter name (`std::shared_ptr<ExperimentalSettings>& /* exp_meta */`), and the body (SwathFile.cpp lines 312-337) returns swath_maps without ever assigning exp_meta. This is a genuine POLS violation: loadSqMass shares the identical out-parameter signature and role 

### [FORM-84] SqMassFile::SqMassConfig::write_full_meta — Config field named/documented 'write_full_meta' also changes load() (read) behavior
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/SqMassFile.h · _format-targeted_

```cpp
bool write_full_meta{true}; ///< Write full meta data (verbose, lossless).
```
- **Expectation:** A flag named write_full_meta and documented as "Write full meta data" affects only store(). When reading, the file's actual contents determine what is loaded.
- **Actual:** SqMassFile::load() applies the same config (setConfig(config_.write_full_meta, ...)) and the read path consults write_full_meta_: if a caller set write_full_meta=false on the config and then calls load(), readExperiment skips reading the embedded full-meta mzML blob and infers structure from raw SQL tables, yielding a metadata-poorer MSExperiment than the file actually contains.
- **Evidence:** SqMassFile.cpp load() lines 24-26: `sql_mass.setConfig(config_.write_full_meta, config_.use_lossy_numpress, config_.linear_fp_mass_acc); sql_mass.readExperiment(map);`. MzMLSqliteHandler.cpp readExperiment uses it at line 286 `if (write_full_meta_)` and line 342 `if (!write_full_meta_ || nr_results == 0 || exp_empty)`.
- **Fix:** Rename/clarify the field (e.g. document it as 'controls full-meta on both write and read') or, ABI-safely, split read vs. write semantics so load() ignores write_full_meta and always reads whatever the file provides. Doc-only fix is source/ABI-compatible; behavioral split is source-compatible.
- **Verifier correction:** Graded medium rather than high: the wrong/poorer read is silent, but it is NOT the default. The default value is true (reads the full blob correctly); triggering the surprise requires the caller to affirmatively set write_full_meta=false and reuse the same SqMassFile/config for load(). On that path the result is silently a metadata-poorer experiment (data-quality loss on read), recoverable only by re-reading with the flag set true. ABI: the recommended doc-only clarification is ABI/source-compatible (none); the meaningful behavioral fix (make load() ignore write_full_meta and always read whatever the file provides) is source-compatible — no signature change — hence abi_impact = source-compatible.
- **Verified:** Independently verified in code. SqMassFile.h line 44 declares `bool write_full_meta{true}; ///< Write full meta data (verbose, lossless).` — write-only language. But SqMassFile::load() (SqMassFile.cpp:24-26) applies this same config to the handler before calling readExperiment. In MzMLSqliteHandler::readExperiment (MzMLSqliteHandler.cpp), line 286 `if (write_full_meta_)` gates whether the embedded

### [FORM-86] TransformationXMLFile::load — load() by default fits a model as a side effect of reading the file
`medium` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/FORMAT/TransformationXMLFile.h · _format-targeted_

```cpp
void load(const std::string& filename, TransformationDescription& transformation, bool fit_model=true)
```
- **Expectation:** A File::load reads stored data points (and any stored model type/params) into the object. Constructing/fitting a numerical model is a separate, potentially expensive and failure-prone step the caller would opt into.
- **Actual:** With the default fit_model=true, load() reads the data points and then calls transformation.fitModel(model_type_, params_) — i.e. running a regression/interpolation fit as a hidden side effect of 'loading'. Fitting can fail or be costly; a caller who only wanted the stored pairs gets more than load implies.
- **Evidence:** Header line 48: `void load(const std::string& filename, TransformationDescription& transformation, bool fit_model=true);`. TransformationXMLFile.cpp lines 36-41: `transformation.setDataPoints(data_); if (fit_model) { transformation.fitModel(model_type_, params_); }`.
- **Fix:** Keep the parameter (ABI) but strengthen the header doc to state that the default triggers a model fit, so the side effect is discoverable at the declaration. Optionally flip guidance toward fit_model=false for pure round-trip reads. Doc-only fix is source/ABI-compatible.
- **Verifier correction:** load() with the default fit_model=true does set the stored data points AND then re-fits the file's stored model type via TransformationDescription::fitModel (a genuine regression/interpolation that can be costly and can throw IllegalArgument for unsupported model types). This is a deliberate, tested round-trip design (a TransformationDescription's apply() requires a fitted model; store() always writes a model type), not a silent data hazard: data points are always populated and failures are loud. The actionable flaw is that fit_model and its fitting side effect are entirely absent from the header doc block, so the default is undiscoverable at the declaration. Fix is doc-only (source/ABI-compatible): document that the default triggers a model fit and note fit_model=false for pure data-point reads.
- **Verified:** Evidence verified exactly. Header line 48: `void load(const std::string& filename, TransformationDescription& transformation, bool fit_model=true);` and TransformationXMLFile.cpp:36-41 do `transformation.setDataPoints(data_); if (fit_model) { transformation.fitModel(model_type_, params_); }`. fitModel (TransformationDescription.cpp:68-103) performs real numerical fitting (constructs Transformation

### [FORM-11] CsvFile::load / CsvFile::CsvFile(const std::string&, char, bool, Int) — CsvFile constructor and load() disagree on line trimming (one trims, the other does not)
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/source/FORMAT/CsvFile.cpp · _format-text-streams_

```cpp
void load(const std::string& filename, char is = ',', bool ie = false, Int first_n = -1); CsvFile(const std::string& filename, char is = ',', bool ie = false, Int first_n = -1);
```
- **Expectation:** A caller expects `CsvFile f(name, is, ie, n);` and `CsvFile f; f.load(name, is, ie, n);` to produce identical buffers, since they take the same arguments and are documented identically.
- **Actual:** The constructor calls `TextFile::load(filename, false, first_n, false, "#")` (trim_lines = false), while `load()` calls `TextFile::load(filename, true, first_n, false, "#")` (trim_lines = true). The two entry points silently differ in whether leading/trailing whitespace is stripped from every line, which changes the parsed field values for files with surrounding spaces.
- **Evidence:** Constructor: `TextFile::load(filename, false, first_n, false, "#");`  vs  load(): `TextFile::load(filename, true, first_n, false, "#");` (CsvFile.cpp lines 27 and 34). The header docs for both are word-for-word identical and mention nothing about trimming.
- **Fix:** Make both call sites pass the same `trim_lines` value (and ideally expose it as a documented parameter). Pick one behavior, document it, and route the constructor through `load()` so they cannot drift. ABI-safe: changing the constant in the constructor body is purely an implementation change.
- **Verifier correction:** Severity is medium rather than high: the divergence silently produces different field values (no error/warning), but it only affects whitespace at the very beginning/end of a line — not whitespace around interior fields or around the delimiter — and only for inputs that actually carry such surrounding whitespace. Well-formed CSVs are unaffected, and the discrepancy is recoverable. The claim's evidence, line numbers, documentation observation, and recommendation (route the constructor through load() so the two cannot drift; changing the constant is ABI-safe) are all accurate.
- **Verified:** Independently verified in src/openms/source/FORMAT/CsvFile.cpp. Constructor (line 27) calls TextFile::load(filename, false, first_n, false, "#") while load() (line 34) calls TextFile::load(filename, true, first_n, false, "#"). The 2nd argument is TextFile's documented `trim_lines` parameter ("Whether or not the lines are trimmed when reading them from file"), and StringUtils::trim (StringUtils.h:5

### [FORM-13] CsvFile::getRow — getRow blindly strips first+last char of every field when itemenclosed_ is set, corrupting unquoted/empty fields
`medium` · `silent-failure` · ABI: `none` · src/openms/source/FORMAT/CsvFile.cpp · _format-text-streams_

```cpp
bool getRow(Size row, StringList& list) const
```
- **Expectation:** With enclosure enabled, a caller expects getRow to remove the enclosing quote characters from quoted fields and leave other fields intact (or to report an error on a malformed field).
- **Actual:** For every field it executes `list[i] = list[i].substr(1, list[i].size() - 2)` unconditionally when itemenclosed_ is true, without checking that the field is actually quoted. An unquoted field has its first and last characters chopped off; an empty field (size 0) produces a `substr(1, (size_t)-2)` call with surprising/UB-adjacent arguments. The function still returns true, so the caller sees corrupted data as success.
- **Evidence:** CsvFile.cpp lines 74-80: `if (itemenclosed_) { list[i] = list[i].substr(1, list[i].size() - 2); }` with no check that list[i] starts/ends with the enclosing character and no guard against size < 2.
- **Fix:** Guard the strip: only remove enclosing characters when the field actually begins and ends with the configured quote char and has length >= 2; otherwise leave it or signal an error. ABI-safe (implementation-only fix inside the .cpp).
- **Verified:** The quoted evidence is real and verified at CsvFile.cpp lines 74-80: when itemenclosed_ is true, getRow runs `list[i] = list[i].substr(1, list[i].size() - 2)` for EVERY field with no check that the field actually starts/ends with a quote char and no size guard, then returns true. StringList is std::vector<std::string> (ListUtils.h:44), so these are plain std::string substr calls. I compiled the ex

### [FORM-14] SVOutStream::operator<<(std::string) — Stream insertion operator<< throws on an embedded newline instead of setting the stream fail state
`medium` · `surprising-throw` · ABI: `source-compatible` · src/openms/source/FORMAT/SVOutStream.cpp · _format-text-streams_

```cpp
SVOutStream& operator<<(std::string str)
```
- **Expectation:** A `std::ostream`-derived `operator<<` is expected to be no-throw on data content and to signal problems via the stream's fail/bad bits; callers writing `svout << somefield;` do not wrap inserts in try/catch.
- **Actual:** operator<<(std::string) throws Exception::IllegalArgument whenever the argument contains '\n'. Since SVOutStream is publicly derived from std::ostream and used in normal `<<` chains (often inside loops over user/derived data), a single data value containing a newline aborts the whole write with an exception from an insertion operator.
- **Evidence:** SVOutStream.cpp lines 60-65: `if (str.contains('\n')) { throw Exception::IllegalArgument(..."argument must not contain newline characters"); }`. The class is `class OPENMS_DLLAPI SVOutStream : public std::ostream`.
- **Fix:** This is a deliberate contract, but it surprises ostream users. Document the throw prominently at the class level (currently only on the operators), and consider an opt-in mode that substitutes/escapes newlines (like the separator replacement) instead of throwing. ABI-safe to add such a mode additively.
- **Verified:** Confirmed by direct read of SVOutStream.cpp lines 60-65 and SVOutStream.h line 30-31; corroborated by uncaught caller usage in src/topp/TextExporter.cpp.

### [FORM-146] PeakFileOptions::setRTRange — setRTRange flips has_rt_range_ on emptiness; sibling range setters always set true
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/source/FORMAT/OPTIONS/PeakFileOptions.cpp · _format-validators-options_

```cpp
void setRTRange(const DRange<1>& range)
```
- **Expectation:** After calling setRTRange(range), hasRTRange() returns true -- just like setMZRange/setIntensityRange/setPrecursorMZRange make their corresponding has*Range() return true. The four range setters are documented as a parallel family ('restricts the range of ... values').
- **Actual:** setRTRange does `has_rt_range_ = !rt_range_.isEmpty();` (conditional), while setMZRange, setIntensityRange and setPrecursorMZRange all unconditionally do `has_*_range_ = true;`. DIntervalBase::isEmpty() returns true for a default-constructed/cleared DRange (min=+inf,max=-inf). So setRTRange(DRange<1>()) leaves hasRTRange()==false, whereas setMZRange(DRange<1>()) makes hasMZRange()==true. The RT setter is silently special-cased relative to its siblings.
- **Evidence:** setRTRange: `rt_range_ = range; has_rt_range_ = !rt_range_.isEmpty();` (PeakFileOptions.cpp:63-64) vs setMZRange: `mz_range_ = range; has_mz_range_ = true;` (lines 79-80), setIntensityRange: `has_intensity_range_ = true;` (line 96), setPrecursorMZRange: `has_precursor_mz_range_ = true;` (line 112). FeatureFileOptions::setRTRange uses the unconditional form (`has_rt_range_ = true;`, FeatureFileOptions.cpp:61), so PeakFileOptions is also inconsistent with its sibling class.
- **Fix:** Pick one convention across all four setters. The least-surprising, ABI-safe fix is to make setMZRange/setIntensityRange/setPrecursorMZRange match setRTRange's empty-check (so an empty range never registers as a filter) OR make setRTRange unconditional like the others; either is a source/binary-compatible body-only change. Document the chosen semantics on all four setters.
- **Verifier correction:** The asymmetry and quoted lines are exactly correct. Refinement: hasFilters() (PeakFileOptions.cpp:301) gates on has_rt_range_ || hasMSLevels() || has_precursor_mz_range_, and reader handlers (MzMLHandler/MzXMLHandler/MzDataHandler/FeatureXMLHandler/ConsensusXMLHandler) all test `hasXRange() && !getXRange().encloses(...)`. An empty range encloses nothing, so on the unconditional setters an empty/default range silently filters out everything (data loss at load), while the RT setter safely no-ops. The surprise is real but the worst consequence requires the unusual act of setting an empty range; with valid ranges all four behave consistently, so severity is medium, not high. ABI: every proposed fix (make all four match either convention) is a body-only change to a non-virtual member defined in the .cpp; no signature/layout/vtable change -> binary- and source-compatible (ABI impact: none).
- **Verified:** Independently verified in source. PeakFileOptions::setRTRange (PeakFileOptions.cpp:63-64) does `has_rt_range_ = !rt_range_.isEmpty();` while setMZRange (80), setIntensityRange (96) and setPrecursorMZRange (112) all unconditionally do `has_*_range_ = true;`. DIntervalBase::isEmpty() (DIntervalBase.h:216) returns true for a default-constructed/cleared DRange, so setRTRange(DRange<1>()) leaves hasRTR

### [FORM-147] XMLValidator::isValid — isValid swallows all parse exceptions and can return true on a parse that never completed; leaks parser on throw
`medium` · `silent-failure` · ABI: `none` · src/openms/source/FORMAT/VALIDATORS/XMLValidator.cpp · _format-validators-options_

```cpp
bool isValid(const std::string& filename, const std::string& schema, std::ostream& os = std::cerr)
```
- **Expectation:** A validator named isValid() that documents throwing FileNotFound and ParseError should surface failures: if the underlying parse blows up, the caller should see a thrown error or a false result, not a silent 'valid'.
- **Actual:** The whole parse is wrapped in `try { parser->parse(source); delete(parser); } catch (...) { /* nothing to do here */ }`. Any exception thrown by parse() is silently discarded and the function then returns valid_. valid_ is initialized to true and is only set false by the SAX error callbacks; if parse() throws before those callbacks fire (e.g. a low-level Xerces failure), isValid() returns true for a file that was never validated. Additionally, because `delete(parser)` is inside the try before the catch, a throw from parse() leaks the SAX2XMLReader.
- **Evidence:** `try { parser->parse(source); delete(parser); } catch (...) { /// nothing to do here }` then `return valid_;` (XMLValidator.cpp:71-81). valid_ starts true (`valid_(true)`, line 23).
- **Fix:** Move `delete(parser)` out of the try (or use unique_ptr) so it is not leaked on throw, and on catch either set valid_=false or rethrow as ParseError so a failed parse cannot masquerade as a valid file. Body-only change; ABI-safe.
- **Verifier correction:** isValid does not swallow ALL invalid-file cases: genuine schema violations and malformed XML are reported correctly via the fatalError/error SAX callbacks (valid_ set false). The real defect is narrower but still genuine: when parser->parse(source) throws a Xerces-level exception (e.g. low-level I/O failure mid-read, encoding/memory/internal error) that bypasses the SAX ErrorHandler callbacks, catch(...) discards it and the function returns the still-true valid_, so a parse that never completed reports the file as valid. Secondly, delete(parser) sits inside the try after parse(), so a throw leaks the SAX2XMLReader (a single, short-lived object — minor). Note loadGrammar (line 65) is outside the try entirely, so a schema-load failure there propagates normally. Fix is body-only/ABI-safe: own the parser with unique_ptr (or move delete out of the try), and on catch set valid_=false or rethrow as ParseError to honor the documented contract.
- **Verified:** Evidence verified verbatim in src/openms/source/FORMAT/VALIDATORS/XMLValidator.cpp. valid_ is initialized true (line 23) and is set false ONLY by the SAX warning/error/fatalError callbacks (lines 89/98/107). Lines 71-79 wrap the parse in `try { parser->parse(source); delete(parser); } catch (...) { /// nothing to do here }`, then line 81 `return valid_;`. Both claimed defects are real: (1) Silent-

### [FORM-105] HDF5Connector::HDF5Connector — Default constructor opens an existing HDF5 file READ-WRITE (must exist), not read-only
`medium` · `surprising-default` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/HDF5Connector.h · _format-vendor_

```cpp
HDF5Connector(const std::string& filename, bool createNewFile = false)
```
- **Expectation:** A 'Connector' to an existing file with the default flag (createNewFile=false) would open the file for reading, or at worst fail gracefully on a read-only file.
- **Actual:** HDF5Connector.cpp sets 'unsigned int openFlag = H5F_ACC_RDWR;' and only switches to H5F_ACC_TRUNC when createNewFile is true. There is no read-only path at all. Opening any existing HDF5 file with the default ctor acquires a read/write handle, fails on a read-only file/filesystem, and the ctor's own comment documents H5F_ACC_RDONLY as a flag that is never used.
- **Evidence:** HDF5Connector.cpp:42 `unsigned int openFlag = H5F_ACC_RDWR;` / :43-46 `if (createNewFile) { openFlag = H5F_ACC_TRUNC; }` ; header :35 `HDF5Connector(const std::string& filename, bool createNewFile = false);`
- **Fix:** Additive/ABI-safe: add a third mode (e.g. an enum OpenMode{READ_ONLY, READ_WRITE, CREATE} parameter with a new overload, or a bool read_only) so callers can open existing files read-only. At minimum document loudly in the header that the default opens RDWR and requires write access. Do not silently change the existing default (would break write callers).
- **Verifier correction:** The default constructor (createNewFile=false) opens an existing HDF5 file READ-WRITE (H5F_ACC_RDWR), not read-only; there is no read-only code path and the only RDONLY mention is a dead .cpp comment. This violates least-surprise for a class named 'Connector' whose default is meant to connect to (read) an existing file. Severity is medium, not high: opening a read-only file/filesystem fails LOUDLY via a thrown HDF5 exception rather than silently producing wrong data or losing data (RDWR does not truncate). The recommended fix (new OpenMode enum overload, or a new defaulted read_only bool / overload) is additive and source-compatible — do not change the existing default, which write callers rely on; at minimum, document in the header that the default opens RDWR and requires write access.
- **Verified:** Evidence verified exactly against source. Header (HDF5Connector.h:35) declares HDF5Connector(const std::string& filename, bool createNewFile = false) as a public OPENMS_DLLAPI, Doxygen-documented FileIO API. HDF5Connector.cpp:42 sets `unsigned int openFlag = H5F_ACC_RDWR;` and :43-46 only switches to H5F_ACC_TRUNC when createNewFile is true. There is no read-only branch; H5F_ACC_RDONLY appears sol

### [FORM-97] ParquetFile (class doc) — Class doc advertises 'Zip / unzip of Parquet directory bundles' as a capability, but no zip/unzip methods exist in the public API
`low` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/FORMAT/ParquetFile.h · _format-arrow-parquet_

```cpp
class OPENMS_DLLAPI ParquetFile { ... }
```
- **Expectation:** A capabilities list in the class brief should correspond to callable public methods.
- **Actual:** The brief lists 'Zip / unzip of Parquet directory bundles (e.g. .oswpq archives)' and discusses store-only compression, but the public interface exposes only Arrow builder helpers, writeTable/readTable, column/value accessors, jsonEscape, and rowCount — no zip or unzip entry points. A developer scanning the header for the advertised archive support finds nothing to call.
- **Evidence:** ParquetFile.h:38 '- Zip / unzip of Parquet directory bundles (e.g. .oswpq archives)' and :41 'All Parquet zip archives use store-only compression'; the only public members are appendOrThrow/finishArray/writeTable/readTable/getColumn/getOptionalColumn/getInt64/getDouble/getBool/getString/getStringList/jsonEscape/rowCount (lines 62-264). grep for zip/unzip in the header finds only the doc lines.
- **Fix:** Either remove the zip/unzip bullet from the brief or expose the documented archive helpers as public static methods. Doc-only change is ABI-neutral.
- **Verifier correction:** The class doc brief in ParquetFile.h (lines 38, 41) advertises zip/unzip of .oswpq archive bundles and store-only compression, but ParquetFile exposes no such methods. The capability actually lives in the separate ZipArchiveFile class (zipDirectory/unzipDirectory/etc.). The ParquetFile brief is a stale doc remnant from before the archive logic was extracted; the .cpp also retains a dead `#include <zip.h>`/`OPENMS_HAVE_LIBZIP` and a misleading "Parquet archive utilities" section header. Fix: remove the zip/unzip bullet and store-only-compression paragraph from the brief (and the dead include/macro), optionally adding a @see ZipArchiveFile cross-reference. Doc-only, ABI-neutral.
- **Verified:** Evidence verified against the actual code. ParquetFile.h:38 and :41 do advertise "Zip / unzip of Parquet directory bundles (e.g. .oswpq archives)" and "All Parquet zip archives use store-only compression (-0)" in the class capabilities brief, yet the public interface (lines 62-264) exposes only appendOrThrow/finishArray/writeTable/readTable/getColumn/getOptionalColumn/getInt64/getDouble/getBool/ge

### [FORM-120] MSDataWritingConsumer::getNrSpectraWritten / getNrChromatogramsWritten — Pure read-only counter getters are non-const
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h · _format-dataaccess-consumers_

```cpp
virtual Size getNrSpectraWritten(); virtual Size getNrChromatogramsWritten();
```
- **Expectation:** A getter named getNrSpectraWritten() that just returns a member should be callable on a `const MSDataWritingConsumer&` and be marked const.
- **Actual:** Both are declared non-const, so they cannot be called through a const reference/pointer, even though the implementation only returns a member (`return spectra_written_;`).
- **Evidence:** MSDataWritingConsumer.h:138 `virtual Size getNrSpectraWritten();` / :143 `virtual Size getNrChromatogramsWritten();`; MSDataWritingConsumer.cpp:147-149 `Size MSDataWritingConsumer::getNrSpectraWritten() {return spectra_written_;}`
- **Fix:** Add `const` to both declarations and definitions. This is technically an ABI change to the virtual signature (mangled name changes), so for strict ABI stability add const-qualified overloads or schedule for the next ABI-breaking release; mark honestly as source-compatible for most callers, ABI-breaking for the vtable.
- **Verifier correction:** The two getters are genuinely non-const and only return a member, so they cannot be invoked on a `const MSDataWritingConsumer&`. This is a legitimate but minor const-correctness wart, severity LOW (no wrong results/data loss/crash possible; only an inability to call through a const handle, and all current callers use non-const objects). The fix (adding `const` to the virtual declarations and definitions) changes the virtual functions' mangled names and vtable layout, so it is ABI-breaking while remaining source-compatible for the existing non-const call sites; schedule for an ABI-breaking release or add const-qualified overloads.
- **Verified:** Evidence verified verbatim. MSDataWritingConsumer.h:138 `virtual Size getNrSpectraWritten();` and :143 `virtual Size getNrChromatogramsWritten();` are both non-const; MSDataWritingConsumer.cpp:147/149 define them as `{return spectra_written_;}` / `{return chromatograms_written_;}` — pure read-only member returns. They are NOT part of the IMSDataConsumer base interface (added in this class), so con

### [FORM-121] MSDataWritingConsumer::consumeSpectrum / consumeChromatogram — consume methods document @param[out] s but operate on a copy, never writing back to the caller's object
`low` · `documentation` · ABI: `none` · src/openms/include/OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h · _format-dataaccess-consumers_

```cpp
void consumeSpectrum(SpectrumType & s) override; (doc: @param[out] s)
```
- **Expectation:** The header annotates the argument `@param[out] s The spectrum to be written to mzML`, and the IMSDataConsumer interface declares `@param[in,out] s ... possibly modified`. A caller reading `[out]` would expect processSpectrum_ transformations (or dataprocessing) to be reflected in their own spectrum after the call.
- **Actual:** The implementation copies the input first (`SpectrumType scpy = s;`), processes and writes the copy, and never assigns back to `s`. The caller's object is left unchanged, contradicting the `@param[out]` annotation. This is also inconsistent with sibling consumers (Cached/Sql) which DO mutate the caller's object.
- **Evidence:** MSDataWritingConsumer.h:110 `@param[out] s The spectrum to be written to mzML`; MSDataWritingConsumer.cpp:62-63 `SpectrumType scpy = s; processSpectrum_(scpy);` (and analogously :109 for chromatograms) — no write-back to `s`.
- **Fix:** Change the doc annotation from `@param[out] s` to `@param[in] s` (it is effectively read-only). Pure doc fix; ABI-neutral. Note the cross-class inconsistency (mutates vs copies) for callers writing generic IMSDataConsumer code.
- **Verifier correction:** consumeSpectrum/consumeChromatogram treat their reference argument as read-only (process a local copy and never write back), so the Doxygen `@param[out] s`/`@param[out] c` annotations in MSDataWritingConsumer.h:110/120 are incorrect and should be `@param[in]`. Documentation defect only (ABI-neutral, no runtime impact: written mzML is correct). The behavior also differs from sibling consumers (MSDataCachedConsumer, MSDataSqlConsumer) which DO mutate the caller's object via clear(), but the base IMSDataConsumer interface documents `@param[in,out] ... possibly modified`, which permits both, so generic callers are already warned not to rely on mutation.
- **Verified:** All quoted evidence verified independently. Header MSDataWritingConsumer.h:110 documents `@param[out] s The spectrum to be written to mzML` and :120 `@param[out] c`. The implementation (MSDataWritingConsumer.cpp:62-63) does `SpectrumType scpy = s; processSpectrum_(scpy);` and :109-110 `ChromatogramType ccpy = c; processChromatogram_(ccpy);` — it processes/writes the COPY and never assigns back to 

### [FORM-122] MSDataSqlConsumer::MSDataSqlConsumer (buffer_size) — buffer_size is a signed int stored/compared as unsigned size_t (flush_after_), so a negative buffer flushes on every element
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h · _format-dataaccess-consumers_

```cpp
MSDataSqlConsumer(const std::string& sql_filename, UInt64 run_id, int buffer_size = 500, bool full_meta = true, bool lossy_compression=false, double linear_mass_acc=1e-4);
```
- **Expectation:** buffer_size is a count; passing a small or zero/negative value should be rejected or behave sanely. The parameter is also renamed to flush_after in the .cpp, obscuring intent.
- **Actual:** `int buffer_size` is stored into `size_t flush_after_` and used as `spectra_.size() >= flush_after_`. A negative buffer_size wraps to a huge size_t (data never auto-flushes until destructor); buffer_size==0 makes `reserve(0)` and flush after every spectrum. No validation. Header name (buffer_size) and impl name (flush_after) diverge.
- **Evidence:** MSDataSqlConsumer.h:55 `int buffer_size = 500`; MSDataSqlConsumer.cpp:16-19 constructor param named `int flush_after` assigned to `flush_after_(flush_after)`; :85 `if (spectra_.size() >= flush_after_)`. flush_after_ is `size_t` (header:96 `size_t flush_after_;`).
- **Fix:** Use an unsigned type (Size/size_t) for the parameter or validate `buffer_size > 0`, and align the .cpp parameter name with the header (buffer_size). Changing the parameter type is ABI-breaking; a guard/clamp inside the body is ABI-neutral and prevents the wrap-around surprise.
- **Verifier correction:** Corrected claim: buffer_size is declared `int` (header) but stored in the unsigned `size_t flush_after_` member and used only in `>=` comparisons; there is no validation. A negative value does NOT "flush on every element" — it wraps to a huge size_t and throws std::length_error at `spectra_.reserve(flush_after_)` in the constructor (lines 22-23), i.e. a loud immediate crash, not silent buffering. The actual flush-every-element behavior occurs for buffer_size==0 (size>=0 is always true), which is functionally correct output but defeats batching (performance only). The header/impl parameter-name divergence (buffer_size vs flush_after) is cosmetic and has zero effect on callers or ABI. Net: a low-severity unvalidated-parameter readability/robustness nit; fixing via a body-side guard/clamp is ABI-neutral. (Separately and out of scope: the call site src/topp/OpenSwathMzMLFileCacher.cpp:162 passes only 5 args and shifts bool/double arguments into the wrong parameters — a distinct, more serious bug not covered by this claim.)
- **Verified:** The structural facts are verified: header (line 55) declares `int buffer_size = 500`; the .cpp constructor (line 16) renames the param to `int flush_after` and assigns it into the `size_t flush_after_` member (line 19), which is compared as `spectra_.size() >= flush_after_` (lines 85, 99). The signed-int-into-unsigned-size_t conversion with no validation, and the header/impl name divergence, are r

### [FORM-124] SiriusFragmentAnnotation::extractAnnotationsFromSiriusFile — Input-only `decoy` flag is documented as @param[out]
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h · _format-dataaccess-misc_

```cpp
static std::vector<MSSpectrum> extractAnnotationsFromSiriusFile(const std::string& path_to_sirius_workspace, Size max_rank = 1, bool decoy = false, bool use_exact_mass = false)
```
- **Expectation:** A reader of `@param[out] decoy` (line 79) expects `decoy` to be an output the function writes to, i.e. the function tells you whether decoys were found.
- **Actual:** `decoy` is a plain `bool` passed by value and is read-only input that selects which subfolder to parse (cpp:264-268). It is never written. The @param[out] tag is wrong (same mistake on `decoy_generation` in the sibling method).
- **Evidence:** Header line 79: "@param[out] decoy Extract annotations for decoys? Or else targets. Run twice if you want both" for a `bool decoy = false` by-value parameter.
- **Fix:** Change tag to @param[in]. Documentation-only; no ABI impact. (The fact that it is a bare bool that also flips use_exact_mass off internally — cpp:266 `if (decoy) use_exact_mass = false;` — is a separate hidden side effect worth a note in the doc.)
- **Verifier correction:** The @param[out] decoy tag (header line 79) is wrong; decoy is a read-only by-value input selecting the /decoys/ vs /spectra/ subfolder (cpp line 268) and is never written, so it should be @param[in]. Severity is low (not high/medium): the by-value bool type signature directly contradicts the [out] tag, so the doc error is obvious and self-correcting at the call site rather than silently producing wrong results. The same mistaken tag exists on decoy_generation in the sibling extractAndResolveSiriusAnnotations (header line 46). A genuinely more useful doc note would flag the hidden side effect at cpp line 266 (`if (decoy) use_exact_mass = false;`), which silently ignores the caller's use_exact_mass when decoy is true.
- **Verified:** All quoted evidence verified independently. Header line 79 reads "@param[out] decoy Extract annotations for decoys? Or else targets. Run twice if you want both" for the declaration at line 82: `bool decoy = false` passed by value. In the cpp (extractAnnotationsFromSiriusFile, line 264), `decoy` is a by-value bool that is only ever READ: line 268 `std::string subfolder = decoy ? "/decoys/" : "/spec

### [FORM-125] SiriusFragmentAnnotation::extractAndResolveSiriusAnnotations — Input-only `decoy_generation` flag documented as @param[out]
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h · _format-dataaccess-misc_

```cpp
static std::vector<SiriusTargetDecoySpectra> extractAndResolveSiriusAnnotations(const std::vector<std::string>& sirius_workspace_subdirs, double score_threshold, bool use_exact_mass, bool decoy_generation)
```
- **Expectation:** `@param[out] decoy_generation` (line 46) implies the function returns/writes this flag.
- **Actual:** `decoy_generation` is a by-value `bool` input that merely gates whether decoy spectra are also extracted (cpp:41-48). It is never an output.
- **Evidence:** Header line 46: "@param[out] decoy_generation Extract decoy spectra from SIRIUS subdirectories." for parameter `bool decoy_generation`.
- **Fix:** Retag as @param[in]. Documentation-only; no ABI impact.
- **Verifier correction:** The claim is correct as stated. Minor refinement: severity is low because `decoy_generation` is pass-by-value, so the type alone makes the `@param[out]` tag self-evidently a typo (a by-value bool cannot be an output); a reader cannot be functionally misled. Same documentation error also appears on line 79 (`@param[out] decoy` for by-value bool), confirming it is a recurring doc-tag typo. Fix: retag both as @param[in].
- **Verified:** Verified against source. Header line 46 documents `@param[out] decoy_generation Extract decoy spectra from SIRIUS subdirectories.` for a by-value `bool decoy_generation` parameter (header line 49). The .cpp (SiriusFragmentAnnotation.cpp:21,41) confirms it is a pass-by-value bool that is only READ — used solely in `if (decoy_generation)` (line 41) to gate whether decoys are extracted (extractAnnota

### [FORM-126] MSChromatogramParquetConsumer::consumeChromatogram — consumeChromatogram clears the caller's chromatogram, but the doc only says it appends
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h · _format-dataaccess-misc_

```cpp
void consumeChromatogram(ChromatogramType& c) override
```
- **Expectation:** The header comment "Consume a chromatogram and append it to the parquet output" reads as a write-only export; a caller would not expect their chromatogram object to be emptied. The sibling MobilogramParquetConsumer goes out of its way to document `const Mobilogram& m` as "not modified", setting an expectation of non-mutation across this family.
- **Actual:** consumeChromatogram encodes then calls `c.clear(false)` (cpp:283), wiping the peak data from the caller's object. consumeSpectrum is documented "no-op; spectra are ignored" but also clears: `s.clear(false)` (cpp:216-217). The base IMSDataConsumer does document [in,out]/"possibly modified", but this consumer's own header comments omit the clear and one even says "no-op".
- **Evidence:** Header line 81: "Consume a spectrum (no-op; spectra are ignored for chromatogram export)." vs cpp:214-217 `consumeSpectrum(...) { s.clear(false); }`. Header line 83: "append it to the parquet output" vs cpp:283 `c.clear(false);`.
- **Fix:** Update the header comments to state that both the spectrum and chromatogram arguments are cleared after consumption (matching the IMSDataConsumer [in,out] contract). Documentation-only; no ABI impact.
- **Verified:** Evidence verified line-for-line. cpp:283 `c.clear(false)` wipes the caller's chromatogram peak data after encoding, and cpp:214-217 `consumeSpectrum(){ s.clear(false); }` mutates the spectrum despite the header (line 81) calling it a "no-op". The header for consumeChromatogram (line 83) says only "append it to the parquet output" with no mention of the clear. So the surprise is real at the DERIVED

### [FORM-127] MSChromatogramParquetConsumer::consumeChromatogram — consumeChromatogram throws if a chromatogram's native ID has no matching transition/precursor metadata
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h · _format-dataaccess-misc_

```cpp
void consumeChromatogram(ChromatogramType& c) override
```
- **Expectation:** A streaming consumer's consumeChromatogram is typically a best-effort sink; a caller feeding chromatograms through MzMLFile().transform(...) would not expect a single unrecognized native ID to abort the whole transform with an exception mid-stream.
- **Actual:** If the native ID looks like a precursor but no compound metadata entry is found, or a transition native ID has no metadata entry, consumeChromatogram throws Exception::InvalidValue (cpp:244-246 and cpp:277-279). The header documents none of this throwing behavior on consumeChromatogram.
- **Evidence:** cpp:277-279: `throw Exception::InvalidValue(..., "Chromatogram native ID '" + native_id + "' does not have a matching transition metadata entry...")`; header line 83-84 only says "append it to the parquet output".
- **Fix:** Document in the header that consumeChromatogram throws Exception::InvalidValue when the chromatogram's native ID cannot be matched to the transition experiment. Documentation-only; no ABI impact.
- **Verified:** Independently verified. The two throw sites exist verbatim: MSChromatogramParquetConsumer.cpp:244-246 (precursor native ID with no matching compound metadata) and cpp:277-279 (transition native ID with no matching metadata), both throwing Exception::InvalidValue with the quoted messages. The public consumeChromatogram (cpp:714-717) delegates straight to impl_->consumeChromatogram(c) with no try/ca

### [FORM-4] MSNumpressCoder::encodeNPRaw — encodeNPRaw with np_compression==NONE returns leaving `result` empty/unmodified, indistinguishable from a real failure
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MSNumpressCoder.h · _format-encoding_

```cpp
void encodeNPRaw(const std::vector<double>& in, std::string& result, const NumpressConfig& config)
```
- **Expectation:** Calling an 'encode' function with a non-empty input and a valid config should produce some output; if NONE means 'pass through', a caller might expect the raw bytes back.
- **Actual:** When `config.np_compression == NONE`, encodeNPRaw returns immediately without touching `result` (MSNumpressCoder.cpp lines 83-86), so a non-empty input yields an empty result. Because the documented failure contract is also 'result is given back unmodified' (header line 188), the caller cannot distinguish 'NONE means no-op' from 'encoding failed'. Both produce the same empty/unmodified output.
- **Evidence:** MSNumpressCoder.cpp lines 83-86: `if (config.np_compression == NONE) { return; }` (before any write to result). Header line 188: `@note In case of error, "result" is given back unmodified`.
- **Fix:** Document explicitly that NONE is a no-op that leaves result untouched, distinct from failure, so callers don't treat empty output as an error (or vice versa). ABI-safe doc-only fix; no signature change needed.
- **Verifier correction:** When config.np_compression == NONE, encodeNPRaw returns at lines 83-86 leaving `result` UNMODIFIED (it does not clear it, unlike the wrapper encodeNP), which is indistinguishable from the documented error contract (header line 188). The header documents neither that NONE is a no-op nor what state `result` is left in, and the precise behavior is "unchanged/stale", not "empty" as the claim states. Severity is low, not high: all in-tree callers either set a concrete LINEAR/SLOF scheme or guard `!= NONE` before calling, so no real path suffers wrong results or data loss; it is an advanced raw primitive with an under-documented NONE/error overlap. Recommended doc-only clarification (NONE is a no-op leaving result untouched, distinct from failure) is ABI-safe.
- **Verified:** Evidence verified verbatim: MSNumpressCoder.cpp lines 83-86 contain `if (config.np_compression == NONE) { return; }` placed before any write to `result`, and header line 188 documents the error contract as "result is given back unmodified". So NONE and a real failure produce the identical outcome, and the enum is documented only as "No compression is applied" (line 49) without specifying that the 

### [FORM-7] ms::numpress::MSNumpress::decodeLinear — decodeLinear returns size_t but documents returning -1 on bad size; -1 becomes a gigantic unsigned count
`low` · `documentation` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MSNUMPRESS/MSNumpress.h · _format-encoding_

```cpp
size_t decodeLinear(const unsigned char *data, const size_t dataSize, double *result)
```
- **Expectation:** A function returning the 'number of decoded doubles' as size_t returns a small non-negative count; an error is signalled out of band.
- **Actual:** The doc states the return is 'the number of decoded doubles, or -1 if dataSize < 4 or 4 < dataSize < 8' (header line 164). But the return type is `size_t` (unsigned), so the documented -1 is actually SIZE_MAX (~1.8e19). A caller writing `for (size_t i=0;i<decodeLinear(...);++i)` or `out.resize(decodeLinear(...))` on a short buffer triggers a catastrophic loop/allocation rather than seeing an error. The OpenMS wrapper avoids this only because it controls sizes, but the public-header contract is internally contradictory.
- **Evidence:** Header lines 166-169 declare `size_t decodeLinear(...)`; line 164 doc: `@return the number of decoded doubles, or -1 if dataSize < 4 or 4 < dataSize < 8`.
- **Fix:** Vendored code, so prefer not to change the signature for upstream-sync reasons; at minimum correct the doc to say the sentinel is `(size_t)-1` == SIZE_MAX (not a negative number) and that callers MUST guard dataSize >= 8 before use. Doc-only, no ABI impact.
- **Verifier correction:** The header doc (MSNumpress.h:164) is wrong/contradictory, but not in the way claimed. decodeLinear never returns -1/SIZE_MAX. Per MSNumpress.cpp:404-490 it returns a valid count (0, 1, or ri) and THROWS a const char* on insufficient input (dataSize<8, <12, <16). The documented "-1 if dataSize < 4 or 4 < dataSize < 8" sentinel does not exist in the implementation, and the size threshold is 8 (throw), not 4. So there is no real path to a gigantic unsigned count; a caller trusting the doc would get a thrown exception, not a catastrophic resize/loop. The OpenMS wrapper already guards via try/catch (MSNumpressCoder.cpp:304/317-318). Doc-only defect, low severity. Fix: change line 164 to document that the function returns the number of decoded doubles and throws const char* when there are too few bytes — do not claim a (size_t)-1 sentinel, which is also nonexistent.
- **Verified:** Header evidence is verbatim-confirmed: src/openms/include/OpenMS/FORMAT/MSNUMPRESS/MSNumpress.h line 164 documents "@return the number of decoded doubles, or -1 if dataSize < 4 or 4 < dataSize < 8" for a function declared (lines 166-169) as size_t decodeLinear(...). So the doc promises a negative sentinel from an unsigned return type — internally contradictory, and a real documentation defect. HOW

### [FORM-8] OpenMS::compress (MSVC linker export pragma) — Base64.h force-exports a symbol literally named `compress` from the DLL on MSVC, colliding with zlib's global compress()
`low` · `other` · ABI: `none` · src/openms/include/OpenMS/FORMAT/Base64.h · _format-encoding_

```cpp
#pragma comment(linker, "/export:compress")
```
- **Expectation:** A format header should not silently add an exported public symbol with a generic name like `compress` to the shared library.
- **Actual:** Line 33 of Base64.h emits `#pragma comment(linker, "/export:compress")`, which re-exports the C symbol `compress` (zlib's compression function) from the OpenMS DLL. Any downstream consumer linking OpenMS on MSVC gets a `compress` symbol in OpenMS's export table — a generic name highly prone to collision and very surprising for anyone who merely included Base64.h for the Base64 class.
- **Evidence:** Base64.h lines 32-34: `#ifdef OPENMS_COMPILER_MSVC` / `#pragma comment(linker, "/export:compress")` / `#endif`.
- **Fix:** Document why this export exists (it works around zlib `compress` not being re-exported for pyOpenMS/plugins) directly at the pragma, and if possible scope it to the single TU that needs it rather than a public header so every includer doesn't inherit the export. Build-system change; no source ABI of OpenMS classes affected.
- **Verifier correction:** The exported symbol is the global C function `compress` from zlib (invoked at ZlibCompression.cpp:34), not an `OpenMS::compress`; the claim's symbol field is mislabeled but its prose is accurate. The pragma adds `compress` to the OpenMS DLL's export table (at OpenMS link time) on MSVC only; it does not add `compress` to a non-DLL consumer's export table merely by inclusion. No OpenMS C++ class/function ABI is affected. Severity is low, not high: the realistic bad outcome is a loud MSVC linker collision in a downstream DLL that also includes Base64.h and links zlib, not silent miscompilation or data loss.
- **Verified:** Verified directly in src/openms/include/OpenMS/FORMAT/Base64.h: lines 32-34 are exactly `#ifdef OPENMS_COMPILER_MSVC` / `#pragma comment(linker, "/export:compress")` / `#endif`, in a public header (12 real source TUs include it, plus pyOpenMS bindings). The `compress` referenced is the global C function from zlib, actually called at ZlibCompression.cpp:34 (`int zlib_error = compress(...)`); there 

### [FORM-10] Base64::decodeSingleString — decodeSingleString omits the zlib_compression default that every other decode* method provides, breaking the call-convention pattern of the class
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/Base64.h · _format-encoding_

```cpp
static void decodeSingleString(const std::string& in, std::string& out, bool zlib_compression)
```
- **Expectation:** Within one class, sibling decode methods should share the same defaulting convention; decode/decodeIntegers/decodeStrings all default zlib_compression=false.
- **Actual:** decode (line 74), decodeIntegers (line 92), and decodeStrings (line 117) all declare `bool zlib_compression = false`. decodeSingleString (line 126) is the lone decoder that requires the bool explicitly (no default). A caller used to the rest of the class will be surprised that `decodeSingleString(in, out)` fails to compile, and there is no encodeSingleString counterpart (the encode side has encodeStrings) — an asymmetric, inconsistent member of the family.
- **Evidence:** Base64.h line 126: `static void decodeSingleString(const std::string& in, std::string& out, bool zlib_compression);` (no default), vs line 117 `decodeStrings(..., bool zlib_compression = false)` and line 74 `decode(..., bool zlib_compression = false)`.
- **Fix:** Add `= false` to decodeSingleString's zlib_compression parameter to match the rest of the class. Adding a default argument is source-compatible and not an ABI change.
- **Verified:** Evidence verified by reading Base64.h directly. Line 74 decode, line 92 decodeIntegers, and line 117 decodeStrings all declare `bool zlib_compression = false`; line 126 `decodeSingleString(const std::string& in, std::string& out, bool zlib_compression)` is the lone decoder without a default — exactly as claimed. The encode side (encode/encodeIntegers/encodeStrings, all with `= false`) reinforces t

### [FORM-37] MsInspectFile::store — store() second parameter named/typed 'SpectrumType spectrum' while load() uses 'FeatureMapType feature_map'
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/MsInspectFile.h · _format-feature-fasta_

```cpp
void store(const std::string& filename, const SpectrumType& spectrum) const
```
- **Expectation:** In a load/store pair on the same File class, the two halves operate on the same logical object, so the parameter name and template type should match (here a feature map).
- **Actual:** load() declares `template <typename FeatureMapType> ... FeatureMapType& feature_map` but store() declares `template <typename SpectrumType> ... const SpectrumType& spectrum`. The store side mislabels the payload as a spectrum, contradicting the class purpose (MsInspect 'features') and the load() naming.
- **Evidence:** load: `template <typename FeatureMapType> void load(const std::string& filename, FeatureMapType& feature_map)`. store: `template <typename SpectrumType> void store(const std::string& filename, const SpectrumType& spectrum) const`.
- **Fix:** Rename store()'s template parameter and argument to FeatureMapType/feature_map to match load(). Header-only template, so this is source-compatible and ABI-neutral.
- **Verifier correction:** The load/store template-parameter and argument naming mismatch is real (load: FeatureMapType/feature_map; store: SpectrumType/spectrum), but store() is an explicitly documented "NOT IMPLEMENTED" stub that throws Exception::NotImplemented unconditionally, so it causes no silent wrong results. The same stub signature is a repeated codebase convention (identical in SpecArrayFile.h and KroenikFile.h). Severity downgraded from medium to low: cosmetic naming nit on dead, loud-failing code, not a misuse trap.
- **Verified:** The quoted evidence is verbatim accurate: load() declares `template <typename FeatureMapType> void load(const std::string&, FeatureMapType& feature_map)` and store() declares `template <typename SpectrumType> void store(const std::string&, const SpectrumType& spectrum) const` (MsInspectFile.h lines 49-50, 139-140). So the load/store naming-and-type inconsistency physically exists. However, the cla

### [FORM-39] FASTAFile::setPosition / writeEnd — setPosition() returns a bool whose meaning is undocumented and not obviously a success flag
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/FORMAT/FASTAFile.h · _format-feature-fasta_

```cpp
bool setPosition(const std::streampos& pos); void writeEnd();
```
- **Expectation:** A bool-returning setPosition() returns whether the seek succeeded, so callers can branch on it.
- **Actual:** There is no documentation of what the bool means; the impl clears the stream and seekg()s, then returns a value. Callers cannot tell from the header whether false means 'past EOF', 'seek failed', or 'previous EOF state'. The sibling stream-mutating call writeEnd() returns void, so the bool here looks meaningful but is unspecified.
- **Evidence:** FASTAFile.h:130 `bool setPosition(const std::streampos& pos);` with no @return doc. Impl (FASTAFile.cpp:227-233): `infile_.clear(); infile_.seekg(pos); ...`
- **Fix:** Document the @return contract (true = seek succeeded / stream good) in the header. Doc-only, ABI-neutral.
- **Verifier correction:** setPosition()'s bool is an in-bounds check, not a seek-success flag: it returns true iff pos <= fileSize_ (after clear()+seekg()), and false only when the requested offset exceeds the file size. It does not report seekg() failure and is independent of EOF state. The real surprise is narrower than claimed and amounts to a missing `@return` doc on a 4-line function whose sole caller already handles false correctly. Recommendation reduces to: add `@return true if pos is within the file (<= file size); false otherwise` to the header. Doc-only, ABI-neutral, low severity.
- **Verified:** Evidence is real: FASTAFile.h:130 declares `bool setPosition(const std::streampos& pos);` with only `/// seek stream to @p pos` and no @return; the impl (FASTAFile.cpp:227-236) is as quoted. BUT the claim mischaracterizes the bool's meaning. The implementation is a trivial 4-line bounds check: `if (pos <= fileSize_) { infile_.clear(); infile_.seekg(pos); return true; } return false;`. So the bool 

### [FORM-40] FeatureXMLFile::loadSize — loadSize() neither loads data nor returns a byte size; it parses to count features
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/FORMAT/FeatureXMLFile.h · _format-feature-fasta_

```cpp
Size loadSize(const std::string& filename);
```
- **Expectation:** A method named loadSize on a File class reads something about size — a competent caller might expect it to return the file's byte size, or to 'load' and report the loaded element count. The name and the fact that it sits next to load() suggest data is loaded.
- **Actual:** loadSize() performs a full XML parse in size-only mode and returns the NUMBER OF FEATURES the file contains, discarding the data into a throwaway dummy FeatureMap. It does not return bytes, and the FeatureMap is not retained. The expensive full-parse side effect and the 'count of features' semantics are invisible from the bare declaration (also undocumented — no Doxygen block).
- **Evidence:** FeatureXMLFile.h:56 `Size loadSize(const std::string& filename);` (no doc). Impl (FeatureXMLFile.cpp:33-43): `FeatureMap dummy; ... handler.setSizeOnly(true); parse_(filename, &handler); return handler.getSize();`
- **Fix:** Add a Doxygen comment clarifying it returns the feature count (not bytes) and that it parses the whole file; consider renaming to countFeatures()/getNumberOfFeatures() in a future ABI-breaking cleanup with the old name kept as a [[deprecated]] alias. Doc addition is ABI-neutral; rename is breaking.
- **Verifier correction:** loadSize() does NOT perform a full or expensive XML parse and does not count features by iteration. It parses only the file header until the `<featureList count="N">` opening tag, reads the declared `count` attribute into expected_size_, and throws EndParsingSoftly (a clean early-termination caught by parse_). It returns that declared feature count (not bytes) using a throwaway dummy FeatureMap. The genuine surprise is the misleading name next to load() plus the absence of any Doxygen comment; the "expensive full parse" framing in the claim is incorrect.
- **Verified:** Evidence quoted is accurate: FeatureXMLFile.h:56 declares `Size loadSize(const std::string& filename);` with NO Doxygen block (the doc at lines 48-53 belongs to load()), and FeatureXMLFile.cpp:33-43 uses a throwaway `FeatureMap dummy`, sets `handler.setSizeOnly(true)`, parses, and returns `handler.getSize()`. The name "loadSize" sitting directly beside `load()` is genuinely misleading: it returns 

### [FORM-41] PEFFFile::getSequence / getModifiedSequence (PEFFEntry) — PEFFEntry getter getSequence()/getModifiedSequence() can throw on malformed sequence content
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/FORMAT/PEFFFile.h · _format-feature-fasta_

```cpp
AASequence getSequence() const; AASequence getModifiedSequence() const;
```
- **Expectation:** A const getter named getSequence() returns the stored sequence; getModifiedSequence() returns it with mods applied. A caller would not expect a getter to throw.
- **Actual:** getSequence() constructs an AASequence from the raw stored 'sequence' string (which the class doc explicitly says may contain 'unusual symbols' kept verbatim). AASequence::fromString throws on unparseable residues, so this getter can throw Exception on entries that loaded fine. The header documents no @exception for these getters.
- **Evidence:** PEFFFile.h:243-247 `AASequence getSequence() const;` doc says only 'Get the base AASequence'. FASTAFile.h class doc (shared FASTA-derived semantics) states unusual symbols 'will be kept'. (Confirmed AASequence construction is the only way to produce an AASequence from a string.)
- **Fix:** Document the throwing behavior (@exception) on getSequence()/getModifiedSequence(), or provide a non-throwing variant returning an optional/empty AASequence. Doc addition is ABI-neutral.
- **Verifier correction:** getSequence()/getModifiedSequence() do call AASequence::fromString(sequence) and can propagate Exception::ParseError, and the header omits an @exception note. But the surprise is narrower than claimed: fromString uses permissive=true, which converts the translation-end '*' (and '#','+') to 'X' rather than throwing, and X/B/Z/J/U/O are valid residues. The throw only occurs for genuinely non-residue characters (e.g. '-', digits, stray punctuation). The fix is an ABI-neutral @exception doc note (optionally a noexcept/optional-returning variant).
- **Verified:** Code confirmed: PEFFEntry::getSequence() is `return AASequence::fromString(sequence);` (PEFFFile.cpp:502-504) and getModifiedSequence() begins the same (510). readEntry_ stores every non-whitespace sequence char verbatim with no alphabet check (PEFFFile.cpp:1643), so an entry that loads fine can hold characters outside the residue alphabet. fromString defaults to permissive=true and throws Excepti

### [FORM-17] FileHandler::getTypeByFileName — getTypeByFileName documents @exception FileNotFound but never touches the filesystem
`low` · `misleading-documentation` · ABI: `none` · src/openms/include/OpenMS/FORMAT/FileHandler.h · _format-filehandling_

```cpp
static FileTypes::Type getTypeByFileName(const std::string& filename)
```
- **Expectation:** Header says '@exception Exception::FileNotFound is thrown if the file is not present'. A caller would expect this function to check that the file exists and to throw on a missing file.
- **Actual:** The implementation (FileHandler.cpp:223-265) is pure string manipulation on the basename/extension. It never calls File::exists and never throws FileNotFound. A non-existent path with a known extension is happily classified (e.g. 'does_not_exist.mzML' -> MZML). The documented exception can never fire.
- **Evidence:** Header line 69: '@exception Exception::FileNotFound is thrown if the file is not present'. Implementation FileHandler.cpp:223 'FileHandler::getTypeByFileName(const std::string& filename) { std::string basename = File::basename(filename) ... return FileTypes::nameToType(tmp); }' — no File::exists / no throw anywhere in the function.
- **Fix:** Fix the doc: remove the '@exception FileNotFound' line from getTypeByFileName (it is a copy-paste from getTypeByContent). This is a doc-only, source/ABI-compatible change. Do NOT add an existence check — callers rely on the pure-string behavior (e.g. for output filenames that do not yet exist).
- **Verifier correction:** getTypeByFileName never throws FileNotFound (the category is more precisely 'misleading/incorrect documentation', not a surprising throw, since no throw exists). The @exception line is a copy-paste from getTypeByContent. Fix is doc-only and source/ABI-compatible: delete the '@exception Exception::FileNotFound' line at FileHandler.h:69. Do not add an existence check — the pure-string behavior is intentional and relied upon for output filenames. Severity is low, not high: zero behavioral impact, no silent wrong results.
- **Verified:** Verified independently. FileHandler.h line 69 documents "@exception Exception::FileNotFound is thrown if the file is not present" on getTypeByFileName, but the implementation (FileHandler.cpp:223-265) is pure string manipulation: File::basename, StringUtils suffix/prefix checks, toUpper, and FileTypes::nameToType. There is no File::exists call and no throw of any kind, so the documented exception 

### [FORM-20] FileHandler::storeSpectrum — storeSpectrum takes a non-const MSSpectrum& even though it only writes (does not modify) the spectrum
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/FileHandler.h · _format-filehandling_

```cpp
void storeSpectrum(const std::string& filename, MSSpectrum& spec, const std::vector<FileTypes::Type> allowed_types = {})
```
- **Expectation:** A store/serialize method takes its payload by const-ref (as storeExperiment takes 'const PeakMap&' and storeFeatures takes 'const FeatureMap&' in this same class). A non-const reference signals to the caller that the spectrum will be modified, and forbids passing a const spectrum or a temporary.
- **Actual:** storeSpectrum takes 'MSSpectrum& spec' (non-const) but the body only forwards it to DTAFile().store(filename, spec) / XMassFile().store(filename, spec), which do not need a mutable spectrum; nothing in storeSpectrum mutates spec (FileHandler.cpp:808-842). The non-const parameter is gratuitous and inconsistent with the other store* methods.
- **Evidence:** Header line 214: 'void storeSpectrum(const std::string& filename, MSSpectrum& spec, ...)'. Compare header line 191 'storeExperiment(..., const PeakMap& exp, ...)' and line 248 'storeFeatures(..., const FeatureMap& map, ...)'. Implementation FileHandler.cpp:808 forwards spec to .store() with no mutation.
- **Fix:** Change the parameter to 'const MSSpectrum& spec' to match the other store* methods and allow storing const/temporary spectra. This is source-compatible for almost all callers (a const-ref is more permissive) but technically an ABI/signature change; do it at the next ABI break or add a const-ref overload now.
- **Verifier correction:** The finding is genuine but severity should be LOW, not a hazard: the non-const MSSpectrum& never causes silently wrong results, data loss, or crashes — it only mildly surprises and over-restricts (cannot store a const or temporary spectrum), and any misuse fails loudly at compile time. Recommend changing the parameter to `const MSSpectrum& spec` to match every sibling store* method. This is source-compatible (all current callers — only the pyOpenMS bindings — keep compiling), though it changes the mangled symbol, so land it at the next ABI break or add a const-ref overload now.
- **Verified:** Independently verified. Header FileHandler.h:214 declares `void storeSpectrum(const std::string&, MSSpectrum& spec, ...)` with a non-const MSSpectrum&. Its own Doxygen at line 209 marks it `@param[in] spec` ("the spectrum to store the data from") — read-only intent. The implementation (FileHandler.cpp:808-842) only forwards spec to DTAFile().store(filename, spec) and XMassFile().store(filename, sp

### [FORM-21] FileHandler::loadIdentifications — loadIdentifications neither clears nor uniformly appends; the OMSSAXML path indexes additional_proteins[0] after a single push_back, corrupting prior contents on a non-empty vector
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/FileHandler.h · _format-filehandling_

```cpp
void loadIdentifications(const std::string& filename, std::vector<ProteinIdentification>& additional_proteins, PeptideIdentificationList& additional_peptides, ...)
```
- **Expectation:** The 'additional_' prefix suggests the method APPENDS to whatever the caller already has. A caller would expect each loaded file's proteins/peptides to be added on top of existing entries, consistently across formats.
- **Actual:** Behavior is format-dependent and not append-consistent. The OMSSAXML branch does additional_proteins.push_back(ProteinIdentification()) and then loads into additional_proteins[0] (FileHandler.cpp:1542-1544) — index 0, NOT the just-added back element. If the caller passes a non-empty vector (the documented 'additional' use case), the OMSSA result is written into the pre-existing first protein ID, silently overwriting it, while peptides are appended elsewhere. PROTXML correctly uses back(). So the same call can clobber or append depending on file type, with no error.
- **Evidence:** FileHandler.cpp:1542 'additional_proteins.push_back(ProteinIdentification()); OMSSAXMLFile().load(filename, additional_proteins[0], additional_peptides, true, true);' vs FileHandler.cpp:1556-1558 (PROTXML) 'additional_proteins.push_back(ProteinIdentification()); ... ProtXMLFile().load(filename, additional_proteins.back(), additional_peptides.back());'.
- **Fix:** Change the OMSSAXML branch to load into additional_proteins.back() (matching PROTXML), so 'additional' semantics hold for all formats. This is an internal bug fix, source/ABI-compatible. Also document in the header whether loadIdentifications appends or replaces (currently unstated).
- **Verifier correction:** The OMSSAXML/PROTXML branch inconsistency (additional_proteins[0] vs .back()) is real, but the claimed corruption mechanism is wrong: OMSSAXMLFile::load (OMSSAXMLFile.cpp:31-33) clears BOTH inputs (protein_identification = ProteinIdentification(); peptide_identifications.clear()), so on a non-empty caller vector it would WIPE all pre-existing peptides (not append them) and reset additional_proteins[0] while leaving the pushed back-element dead — not the "overwrite first protein, append peptides" described. More importantly, no caller passes a non-empty vector with an OMSSA file (all TOPP tools and pyOpenMS bindings pass freshly-empty vectors), and on an empty vector [0]==back() so the branch is correct. The header never documents append semantics. Recommend the cosmetic fix to use additional_proteins.back() for uniformity (source/ABI compatible), but this is a low-severity latent smell, not a high-severity data-corruption-for-reasonable-use bug.
- **Verified:** The quoted evidence is real and accurate: FileHandler.cpp:1542-1544 (OMSSAXML) does additional_proteins.push_back(...) then loads into additional_proteins[0], while FileHandler.cpp:1556-1558 (PROTXML) pushes then loads into additional_proteins.back(). The index-[0]-vs-back() inconsistency genuinely exists. HOWEVER the claim's failure mechanism is materially wrong and the impact is overstated. (1) 

### [FORM-111] QcMLFile::existsRunQualityParameter / existsSetQualityParameter — exists*QualityParameter() is named like a bool predicate but returns void via an out-parameter
`low` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/QcMLFile.h · _format-gnps-deconv-qc_

```cpp
void existsRunQualityParameter(const std::string& filename, const std::string& qpname, std::vector<std::string>& ids) const
```
- **Expectation:** A method named exists... returns bool answering whether the quality parameter exists.
- **Actual:** It returns void and writes the matching parameter IDs into the out-parameter 'ids' (cleared first). A caller must inspect ids.empty() to learn existence. The sibling existsRun()/existsSet() in the very same class DO return bool, so the naming is self-contradictory within one class.
- **Evidence:** Header: 'void existsRunQualityParameter(const std::string& filename, const std::string& qpname, std::vector<std::string>& ids) const'. Impl QcMLFile.cpp:360-382 does 'ids.clear(); ... ids.push_back(qit->id);' and returns void. Contrast 'bool existsRun(const std::string& filename, bool checkname = false) const' (QcMLFile.cpp:323).
- **Fix:** Additive fix: add a bool-returning overload/companion (e.g. hasRunQualityParameter(...)) and keep the void out-param version for ID retrieval, or rename the retrieval method to getRunQualityParameterIDs(...). Pure rename is source-breaking; prefer adding the bool helper. No ABI removal needed.
- **Verifier correction:** The naming inconsistency is real (exists* prefix on a void out-param method, alongside bool-returning existsRun/existsSet in the same class), but it is not a silent-failure hazard: the void return type makes any bool-predicate misuse a compile error, and the Doxygen at lines 141/143 already describes the ID-returning behavior. Both existing callers use it correctly. Therefore severity is low (mild compiler-enforced surprise), not medium. The remedy is additive and source-compatible.
- **Verified:** Evidence verified independently. Header QcMLFile.h:142-144 declares void existsRunQualityParameter(...)/existsSetQualityParameter(...) with a std::vector<std::string>& ids out-param. Impl QcMLFile.cpp:360-382 and 384-408 confirm: ids.clear() then ids.push_back(qit->id), returning void; existence must be inferred from ids.empty(). The sibling predicates existsRun (line 323) and existsSet (line 342)

### [FORM-115] ExperimentalDesignFile::load(const TextFile&, const bool, std::string) — Third 'filename' argument has a default value defined only in the .cpp, invisible/unusable from the header
`low` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/FORMAT/ExperimentalDesignFile.h · _format-gnps-deconv-qc_

```cpp
static ExperimentalDesign load(const TextFile& text_file, const bool require_spectra_file, std::string filename)
```
- **Expectation:** If load() can be called with two arguments, the default for the third parameter is declared in the header so every caller can use it.
- **Actual:** The header declares 'std::string filename' with NO default, but the definition supplies 'std::string filename = "--no design file provided--"'. Callers that include only the header cannot omit the argument; the default is effectively private to the defining translation unit, and defaulting in the definition rather than the declaration is fragile/surprising.
- **Evidence:** Header ExperimentalDesignFile.h:35 'static ExperimentalDesign load(const TextFile& text_file, const bool require_spectra_file, std::string filename);' (no default). Definition ExperimentalDesignFile.cpp:451 'ExperimentalDesignFile::load(const TextFile &text_file, const bool require_spectra_file, std::string filename = "--no design file provided--")'.
- **Fix:** Move the default to the header declaration (and remove it from the definition). Source-compatible and makes the two-arg call usable by all TUs; no ABI change.
- **Verifier correction:** The facts are accurate, but the claim overstates impact. The default in the .cpp is effectively dead: no caller invokes the TextFile overload with only two arguments (all 2-arg `load(file, false)` calls bind to the separate `load(const std::string&, bool)` overload), and even the in-TU caller passes the third argument explicitly. So the "default invisible to callers" condition is real, but it never causes wrong behavior and any misuse would be a loud compile error, not a silent fault. Severity LOW. Fix as recommended: move `= "--no design file provided--"` to the header declaration and remove it from the definition — source-compatible, no ABI change.
- **Verified:** Evidence verified verbatim. Header ExperimentalDesignFile.h:35 declares `static ExperimentalDesign load(const TextFile& text_file, const bool require_spectra_file, std::string filename);` with NO default; the definition at ExperimentalDesignFile.cpp:451 supplies `std::string filename = "--no design file provided--"`. This is legal but smelly C++: a default argument placed in the definition is usab

### [FORM-131] MzDataHandler::MzDataHandler — MzDataHandler constructor doc-comments label read vs write backwards
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/FORMAT/HANDLERS/MzDataHandler.h · _format-handlers-core_

```cpp
MzDataHandler(MapType & exp, const std::string & filename, const std::string & version, ProgressLogger & logger) /*doc: write-only*/; MzDataHandler(const MapType & exp, ...) /*doc: read-only*/
```
- **Expectation:** The non-const `MapType& exp` constructor is the read-into-experiment (loading) ctor; the `const MapType& exp` ctor is the write-out (storing) ctor. The doc comments should match.
- **Actual:** Header comments are swapped: the `MapType & exp` overload is documented '/// Constructor for a write-only handler' and the `const MapType & exp` overload '/// Constructor for a read-only handler'. The .cpp confirms the opposite: the `MapType&` ctor sets `exp_(&exp), cexp_(nullptr)` (read target) and the `const MapType&` ctor sets `exp_(nullptr), cexp_(&exp)` (write source).
- **Evidence:** Header lines 43-47: '/// Constructor for a write-only handler\n MzDataHandler(MapType & exp, ...)' and '/// Constructor for a read-only handler\n MzDataHandler(const MapType & exp, ...)'. MzDataHandler.cpp:18-21 `MzDataHandler(MapType & exp, ...) : ... exp_(&exp), cexp_(nullptr)` and :30-33 `MzDataHandler(const MapType & exp, ...) : ... exp_(nullptr), cexp_(&exp)`.
- **Fix:** Swap the two doc comments so the non-const overload reads 'read-only handler' and the const overload reads 'write-only handler', matching MzMLHandler.h:100-104 and MzXMLHandler.h:46-50. Comment-only change, no ABI impact.
- **Verifier correction:** The claim is accurate as stated. Severity is low rather than higher: the two overloads differ by const-qualification of MapType&, so const-correctness and overload resolution catch any mismatch loudly at compile time and the actual runtime behavior is correct — the swapped comments cause mild surprise/potential confusion but no silently wrong results, data loss, or crash. Fix is a comment-only swap (header lines 43 and 47) to match MzMLHandler.h/MzXMLHandler.h; no ABI impact.
- **Verified:** Independently verified. Header MzDataHandler.h:43-47 labels the non-const `MapType& exp` ctor '/// Constructor for a write-only handler' and the `const MapType& exp` ctor '/// Constructor for a read-only handler'. The .cpp confirms the opposite semantics: the `MapType&` ctor sets exp_(&exp), cexp_(nullptr) and exp_ (declared `MapType*` at h:86) is the mutated parse target (exp_->addSpectrum(spec_)

### [FORM-132] MzMLSpectrumDecoder::setSkipXMLChecks — setSkipXMLChecks parameter named 'only' contradicts its meaning
`low` · `naming-inconsistency` · ABI: `none` · src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h · _format-handlers-core_

```cpp
void setSkipXMLChecks(bool only)
```
- **Expectation:** A setter named setSkipXMLChecks(bool) takes a flag that, when true, skips XML checks. The parameter name should read like 'skip' so `setSkipXMLChecks(true)` clearly means 'skip the checks'.
- **Actual:** The header declares the parameter as `bool only`, which reads at the call site as if it controls some 'only' mode. The implementation in MzMLSpectrumDecoder.cpp:616 names it `bool skip` and just does `skip_xml_checks_ = skip;`. The member doc on line 43/162 says 'Whether to skip some XML checks ... and be fast instead', confirming the boolean is a plain skip flag, not an 'only' selector.
- **Evidence:** Header line 163: `void setSkipXMLChecks(bool only);` vs MzMLSpectrumDecoder.cpp:616 `void MzMLSpectrumDecoder::setSkipXMLChecks(bool skip) { skip_xml_checks_ = skip; }`.
- **Fix:** Rename the header parameter from `only` to `skip` to match the implementation and the member semantics. Parameter-name-only change; no ABI impact.
- **Verifier correction:** The declared parameter name `only` (header line 163) is a stray/leftover token; the definition uses `skip` and assigns straight to skip_xml_checks_, so the boolean's polarity is correct (true = skip checks). This is a pure naming inconsistency, not an inverted/contradictory flag. The name never appears at call sites and is fully disambiguated by the function name and the inline doc comment. Recommend renaming `only` -> `skip` for clarity; no behavioral or ABI impact (parameter names in declarations are non-binding in C++).
- **Verified:** Evidence is literally accurate: header line 163 declares `void setSkipXMLChecks(bool only);` while the definition (MzMLSpectrumDecoder.cpp:616-618) names it `bool skip` and does `skip_xml_checks_ = skip;`. The member doc (lines 43/162) and the constructor (line 98, `bool skip_xml_checks = false`) all confirm true = skip checks. So the parameter name `only` is a genuine inconsistency/leftover. BUT 

### [FORM-136] XMLHandler::reset — reset() documented to 'release internal memory' is an empty no-op in the base class
`low` · `misleading-documentation` · ABI: `none` · src/openms/include/OpenMS/FORMAT/HANDLERS/XMLHandler.h · _format-handlers-core_

```cpp
void reset()
```
- **Expectation:** A method documented 'Release internal memory used for parsing' would actually free or clear parser state when called.
- **Actual:** The base implementation is completely empty: `void XMLHandler::reset() {}`. A caller invoking reset() to reclaim memory or re-initialize parse state on the base handler gets nothing, with no indication. (The doc comment is also truncated: '/// Release internal memory used for parsing (call').
- **Evidence:** Header line 351-352: '/// Release internal memory used for parsing (call\n void reset();'. XMLHandler.cpp:62-64 `void XMLHandler::reset() {}`.
- **Fix:** Either remove the dead method or implement it to clear `open_tags_`, `sm_`, `cv_terms_`, etc.; at minimum finish the doc comment and state it is a no-op unless overridden. Removing would be breaking; documenting/implementing is non-breaking.
- **Verifier correction:** reset() is a confirmed empty no-op with a truncated doc comment, but it is not dead — XMLFile.cpp's XMLCleaner_ RAII guard calls it through a base XMLHandler* during (exception) cleanup. Crucially reset() is non-virtual, so even a hypothetical derived override would not be dispatched via that base-pointer call. Practical impact is benign (memory is reclaimed by the handler destructor moments later), so this is a misleading-doc/dead-design smell at low severity, not a correctness or memory-safety failure. Recommendation stands: finish the doc comment to state it is a no-op unless overridden, and either make it virtual + implement clearing of open_tags_/cv_terms_/sm_ or remove it; making it non-no-op/virtual is the only change that could affect ABI.
- **Verified:** Evidence confirmed verbatim. Header XMLHandler.h:351-352 has the truncated doc comment "/// Release internal memory used for parsing (call" above "void reset();", and XMLHandler.cpp:62-64 defines it as an empty body. The surprise is real and actually sharper than stated: reset() is NOT virtual (unlike the neighboring writeTo/getLoadDetail which are), and it IS wired into a live RAII guard — XMLFil

### [FORM-137] MzMLHandler::getOptions — getOptions() returns a non-const reference to internal options with no const counterpart
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLHandler.h · _format-handlers-core_

```cpp
PeakFileOptions& getOptions()
```
- **Expectation:** A getter named getOptions() on a handler would, by sibling convention, offer a const read accessor (FeatureXMLHandler and ConsensusXMLHandler both provide `const ... getOptions() const`).
- **Actual:** MzMLHandler exposes only `PeakFileOptions& getOptions();` — a single non-const overload returning a mutable reference into the object's internals (`options_`). There is no `const PeakFileOptions& getOptions() const;`, so the options cannot be read from a const MzMLHandler, and any caller can mutate internal parse/write options through the returned reference unexpectedly. This is inconsistent with FeatureXMLHandler.h:63-66 and ConsensusXMLHandler.h:50-55 which both provide const + non-const pairs.
- **Evidence:** MzMLHandler.h:142-145 `void setOptions(...); PeakFileOptions& getOptions();` (no const overload). Contrast FeatureXMLHandler.h:63-66 `FeatureFileOptions& getOptions(); const FeatureFileOptions& getOptions() const;`.
- **Fix:** Add an additive `const PeakFileOptions& getOptions() const;` overload to MzMLHandler (and MzXMLHandler/MzDataHandler which expose no getter at all) to match the FeatureXML/ConsensusXML convention. Purely additive; no ABI break.
- **Verifier correction:** Corrected claim: MzMLHandler (an Internal, "do not use directly" handler) provides only a non-const `PeakFileOptions& getOptions();` and lacks the `const PeakFileOptions& getOptions() const;` overload that its sibling handlers FeatureXMLHandler and ConsensusXMLHandler provide — a minor const-correctness inconsistency. The mutable getter itself is NOT a defect: it is the deliberate, documented configuration idiom shared by every public *File class and used as `getOptions().addMSLevel(1)`; the "any caller can mutate internal options unexpectedly" framing is incorrect. The only genuine surprise is the missing const overload, which has near-zero practical impact because MzMLHandler is never used as a const object and all real callers go through MzMLFile/IndexedMzMLFileLoader, which already expose both overloads. Severity: low (mild API ergonomics surprise, not wrong results/data loss). Recommended fix (purely additive `const` overload) is fine but optional.
- **Verified:** Facts verified independently. MzMLHandler.h:142-145 declares only `void setOptions(const PeakFileOptions&); PeakFileOptions& getOptions();` with no `const` overload. Sibling handlers FeatureXMLHandler.h:63-66 and ConsensusXMLHandler.h:50-55 both provide const+non-const pairs (def confirmed in their .cpp), so MzML is genuinely the odd one out among handlers that offer a getter (MzXML/MzData handler

### [FORM-142] MzMLSqliteHandler::readSpectra / readChromatograms — Doxygen @param directions are wrong: the output vector is tagged [in] and the input 'indices' is tagged [out]
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h · _format-handlers-misc_

```cpp
void readSpectra(std::vector<MSSpectrum> & exp, const std::vector<int> & indices, bool meta_only = false) const;
```
- **Expectation:** @param[in]/@param[out] should match data flow; the result vector should be [out] and the index selection should be [in].
- **Actual:** readSpectra documents the result vector as "@param[in] exp The result" (it is actually the output) and readChromatograms documents both the result and the selection wrongly: "@param[out] exp The result" but "@param[out] indices A list of indices restricting the resulting chromatograms" — indices is an input filter, not an output. A caller trusting the docs could believe indices is populated by the call.
- **Evidence:** Header readSpectra: "@param[in] exp The result". Header readChromatograms: "@param[out] indices A list of indices restricting the resulting chromatograms only to those specified here".
- **Fix:** Doc-only fix: change to @param[out] exp and @param[in] indices for both methods. No ABI or source impact.
- **Verifier correction:** Confirmed doc bug but lower impact than stated. readSpectra: "@param[in] exp" should be "@param[out] exp" (exp is the written result). readChromatograms: "@param[out] indices" should be "@param[in] indices" (indices is a const input filter). The claim's assertion that a caller could believe indices is populated by the call overstates risk: indices is declared const std::vector<int>&, so the compiler prevents any mutation regardless of the doc tag. Doc-only fix; no source or ABI impact. Severity is low, not high/medium.
- **Verified:** Evidence verified verbatim in src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h. readSpectra (line 87) tags the output result vector as "@param[in] exp The result" — wrong; the .cpp (lines 391-411) writes exp via prepareSpectra_/populateSpectraWithData_, so it is an output (the sibling readExperiment at line 75 and readChromatograms at line 96 correctly use @param[out] exp). readChroma

### [FORM-143] MzMLSqliteHandler::writeExperiment / writeRunLevelInformation — Write methods tag a const-ref INPUT experiment as @param[out]
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h · _format-handlers-misc_

```cpp
void writeExperiment(const MSExperiment & exp); / void writeRunLevelInformation(const MSExperiment& exp, bool write_full_meta);
```
- **Expectation:** On a write* method, the source experiment is the input; @param[in] is expected. @param[out] signals the function fills it.
- **Actual:** writeExperiment documents "@param[out] exp The data to write" and writeRunLevelInformation documents "@param[out] exp The result data structure", but both take const MSExperiment& and only read from it. The [out] tag is contradicted by the const& signature and by 'The data to write'.
- **Evidence:** Header writeExperiment: "@param[out] exp The data to write ... void writeExperiment(const MSExperiment & exp);". writeRunLevelInformation: "@param[out] exp The result data structure ... void writeRunLevelInformation(const MSExperiment& exp, bool write_full_meta);"
- **Fix:** Doc-only fix: change @param[out] to @param[in] for both. No ABI/source impact.
- **Verified:** Verified directly in the header: writeExperiment (line 177) documents "@param[out] exp The data to write" with signature void writeExperiment(const MSExperiment & exp) (line 179), and writeRunLevelInformation (line 211) documents "@param[out] exp The result data structure" with signature void writeRunLevelInformation(const MSExperiment& exp, bool write_full_meta) (line 214). The @param[out] tag is

### [FORM-144] MzMLSqliteSwathHandler::readSwathWindows / readMS1Spectra / readSpectraForWindow — Class is documented 'Read-only accessor' but all three read methods are non-const
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteSwathHandler.h · _format-handlers-misc_

```cpp
std::vector<OpenSwath::SwathMap> readSwathWindows(); (etc.)
```
- **Expectation:** A class whose brief is 'Read-only accessor for SWATH/DIA spectrum indices' should offer const read methods, allowing use through a const handle.
- **Actual:** readSwathWindows(), readMS1Spectra() and readSpectraForWindow() are all declared non-const, so the handler cannot be used via a const reference despite the 'read-only' contract. (The class also carries spec_id_/chrom_id_ 'append' members that are never read by any public method, hinting at copy-paste from the writer.)
- **Evidence:** Header brief: "Read-only accessor for SWATH/DIA spectrum indices stored in an sqMass file." Methods: "std::vector<OpenSwath::SwathMap> readSwathWindows();", "std::vector<int> readMS1Spectra();", "std::vector<int> readSpectraForWindow(const OpenSwath::SwathMap & swath_map);" — none const.
- **Fix:** Add const to the three read methods (each opens its own connection on demand per the constructor doc, so there is no logical mutation). const-qualification changes the mangled symbol, so this is technically ABI-breaking on these symbols; safe to do for an Internal class but should be flagged. Source-compatible for callers.
- **Verifier correction:** The three read methods are indeed non-const despite the "Read-only accessor" brief, and the .cpp proves they perform no member mutation (each opens its own connection on filename_), so const-qualifying them is sound. But this is a low-severity const-correctness papercut, not a functional defect: it causes only a compile-time inconvenience for const handles, the lone caller uses a non-const local, and no behavior/results are affected. Adding const would change the mangled symbols, so the recommended fix is ABI-breaking (acceptable for this Internal class, source-compatible for callers). Secondary: spec_id_/chrom_id_ are dead copy-paste leftovers from the writer sibling MzMLSqliteHandler, never read and left uninitialized — harmless but worth removing.
- **Verified:** Every quoted fact checks out against the actual source. Header line 26 reads "Read-only accessor for SWATH/DIA spectrum indices stored in an sqMass file." yet readSwathWindows() (l.72), readMS1Spectra() (l.85) and readSpectraForWindow(...) (l.108) are all non-const. The .cpp confirms none mutate any member: each opens its own SqliteConnector(filename_), reads, finalizes, returns; filename_ is only

### [FORM-145] FidHandler (class) — FidHandler publicly inherits std::ifstream, exposing the raw stream API (seekg/read/close) alongside its index-tracking methods
`low` · `ownership-lifetime` · ABI: `none` · src/openms/include/OpenMS/FORMAT/HANDLERS/FidHandler.h · _format-handlers-misc_

```cpp
class FidHandler : public std::ifstream
```
- **Expectation:** A 'fid File handler' that tracks a position index would be expected to encapsulate the stream so its index_ stays consistent with the file cursor.
- **Actual:** Public inheritance from std::ifstream means callers can call seekg()/read()/close() directly on a FidHandler, moving the underlying file position without updating index_. getIndex() then reports a position out of sync with the actual cursor, and getIntensity() reads from wherever the raw stream was left. The stream's own API silently desynchronizes the handler's tracked index.
- **Evidence:** Header: "class OPENMS_DLLAPI FidHandler : public std::ifstream" with private member "Size index_;" updated only inside getIntensity().
- **Fix:** Prefer composition (hold an ifstream member) or protected/private inheritance so the raw cursor cannot be moved out from under index_. This changes the type's public base and is ABI/source-breaking, so for an Internal class flag it as the ideal fix; short term, document that callers must not use the inherited ifstream positioning API.
- **Verifier correction:** FidHandler does publicly inherit std::ifstream and exposes the raw stream API alongside its getIndex()/getIntensity() index tracking, which is a real public-inheritance-of-std::ifstream anti-pattern that can desync index_ from the file cursor IF a caller uses the inherited positioning API. But this is a low-severity code-quality smell, not a correctness hazard: the class is namespace Internal, explicitly documented '@note Do not use this class directly', read-only, and has a single caller (XMassFile::load) that never misuses the inherited API, so no reasonable use silently misbehaves. The class as-is causes no ABI/source break; only the proposed fix would. Severity downgraded from the implied high/medium to low.
- **Verified:** Code is read correctly: FidHandler.h line 25-26 declares `class OPENMS_DLLAPI FidHandler : public std::ifstream`, and `index_` (line 52) is incremented only inside getIntensity() (FidHandler.cpp line 55). Public inheritance from std::ifstream genuinely exposes seekg/read/close, and index_ tracks a logical point counter that the raw stream API can desync from the cursor. That structural observation

### [FORM-50] XTandemInfile::setPrecursorErrorType / getPrecursorErrorType — "PrecursorErrorType" getter/setter actually control the precursor mass TYPE (mono/average), not an error/tolerance type
`low` · `misleading-name` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/XTandemInfile.h · _format-id-search-in_

```cpp
void setPrecursorErrorType(MassType mono_isotopic); MassType getPrecursorErrorType() const
```
- **Expectation:** In this very class, names ending in "ErrorType"/"ErrorUnit" denote the mass error model: `setPrecursorMassErrorUnit(ErrorUnit)` takes Da/ppm. So `setPrecursorErrorType` reads as another error-related setting, and its parameter (named `mono_isotopic`) compounds the confusion.
- **Actual:** It is the precursor mass TYPE selector (monoisotopic vs average). It takes a `MassType` and writes `precursor_mass_type_`. The sibling for fragments is correctly named `setFragmentMassType`-style via the enum but here the word "Error" is wrong, and the parameter is named `mono_isotopic` in the header yet `mass_type` in the impl.
- **Evidence:** Header: `/// sets the precursor mass type\n    void setPrecursorErrorType(MassType mono_isotopic);` and `/// returns the precursor mass type\n    MassType getPrecursorErrorType() const;`. Impl (XTandemInfile.cpp:649-657): `MassType getPrecursorErrorType() const { return precursor_mass_type_; }  void setPrecursorErrorType(const MassType mass_type) { precursor_mass_type_ = mass_type; }`. Compare sibling `setPrecursorMassErrorUnit(ErrorUnit)`.
- **Fix:** Rename to `setPrecursorMassType`/`getPrecursorMassType` (matching the member `precursor_mass_type_` and the doc) for clarity. ABI-break to rename; ABI-safe path is to add the correctly-named methods as inline forwarders and deprecate the misnamed pair. At minimum, fix the header parameter name to `mass_type`.
- **Verifier correction:** Claim is accurate. Refinement: the misuse risk is type-guarded — setPrecursorErrorType takes MassType (MONOISOTOPIC/AVERAGE) while the error-unit sibling takes ErrorUnit (DALTONS/PPM), so passing the wrong kind of value is a compile error, not a silent miscompute. Hence severity is low (clarity-only), not medium. abi_impact of the recommended rename is breaking; a source-compatible path exists by adding inline forwarders getPrecursorMassType/setPrecursorMassType and deprecating the misnamed pair. Minimal non-ABI fix: change the header parameter name from `mono_isotopic` to `mass_type` to match the impl.
- **Verified:** Independently verified in the source. Header (XTandemInfile.h:73-77) declares `void setPrecursorErrorType(MassType mono_isotopic)` and `MassType getPrecursorErrorType() const`, both documented as "sets/returns the precursor mass type". The impl (XTandemInfile.cpp:649-657) reads/writes `precursor_mass_type_` and names the parameter `mass_type` (mismatching the header's `mono_isotopic`). The `MassTy

### [FORM-51] XTandemInfile::setSemiCleavage / setAllowIsotopeError — Write-only bool options with no matching getter, unlike sibling option pairs in the same class
`low` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/XTandemInfile.h · _format-id-search-in_

```cpp
void setSemiCleavage(const bool); void setAllowIsotopeError(const bool)
```
- **Expectation:** Every other boolean/scalar option on this config class is a symmetric set/get pair (e.g. `setNoiseSuppression`/`getNoiseSuppression`, `setSemi...` would have `getSemi...`). A caller expects to be able to read back a value they set.
- **Actual:** `setSemiCleavage` and `setAllowIsotopeError` exist with no corresponding `getSemiCleavage`/`getAllowIsotopeError`. The values (`semi_cleavage_`, `allow_isotope_error_` members) cannot be retrieved through the public API, breaking the round-trip a caller assumes from the rest of the class.
- **Evidence:** Header lines 158-167: `void setSemiCleavage(const bool semi_cleavage);` and `void setAllowIsotopeError(const bool allow_isotope_error);` have no getters, while immediately below `bool getNoiseSuppression() const;` / `void setNoiseSuppression(const bool ...)` form a pair. grep confirms no `getSemiCleavage`/`getAllowIsotopeError` in header or cpp.
- **Fix:** Add the missing const getters `getSemiCleavage()` and `getAllowIsotopeError()`. This is purely additive and ABI-safe.
- **Verifier correction:** setSemiCleavage and setAllowIsotopeError are indeed write-only (no getters; members read only internally in writeTo_). The asymmetry is real but the cited symmetric foil is imperfect: getNoiseSuppression itself has no backing member/setter, so the class already has incidental asymmetries. Impact is a mild read-back surprise on a one-way file-writer config object, not silently wrong results — hence low severity, not medium/high. Recommendation to add const getSemiCleavage()/getAllowIsotopeError() is valid and ABI-safe (source-compatible, additive).
- **Verified:** Evidence verified. Header lines 158/161 declare void setSemiCleavage(const bool) and void setAllowIsotopeError(const bool); cpp lines 696-704 implement them by assigning the private members semi_cleavage_ (line 255) and allow_isotope_error_ (line 257). grep confirms no getSemiCleavage/getAllowIsotopeError anywhere in header, cpp, or pyOpenMS bindings. The members are read only internally inside wr

### [FORM-53] OMSSAXMLFile::load — Output identification containers documented @param[in]; load() silently clears them
`low` · `documentation-param-direction` · ABI: `none` · src/openms/include/OpenMS/FORMAT/OMSSAXMLFile.h · _format-id-search-in_

```cpp
void load(const std::string& filename, ProteinIdentification& protein_identification, PeptideIdentificationList& id_data, bool load_proteins = true, bool load_empty_hits = true)
```
- **Expectation:** `@param[in]` signals read-only input; a caller may reasonably reuse the same containers across multiple load() calls expecting accumulation.
- **Actual:** load() resets the outputs on entry: `protein_identification = ProteinIdentification();` and `peptide_identifications.clear();`. Pre-existing content is discarded, yet the parameters are tagged `@param[in]` and the `id_data` Doxygen name does not even match the actual parameter name `peptide_identifications`.
- **Evidence:** Header lines 47-49: `@param[in] protein_identification ...` `@param[in] id_data The identifications with m/z and RT`. Impl (OMSSAXMLFile.cpp:32-33): `protein_identification = ProteinIdentification(); peptide_identifications.clear();`
- **Fix:** Mark these `@param[out]`, document that load() clears them first, and fix the Doxygen name mismatch (`id_data` -> `peptide_identifications`). Documentation-only, ABI-safe.
- **Verifier correction:** The output containers in OMSSAXMLFile::load are documented `@param[in]` but are reset on entry (protein_identification reassigned, peptide_identifications cleared), so they should be `@param[out]` per the established OpenMS convention (cf. PepXMLFile.h/MzMLFile.h). There is also a Doxygen/declaration name `id_data` that mismatches the .cpp parameter name `peptide_identifications`. This is a low-severity documentation inconsistency, not a misuse trap: clearing-on-entry is the universal OpenMS loader pattern (PepXMLFile::load, ProtXMLFile::load do the same) and is explicitly commented in the source, so expecting accumulation is not a reasonable use. Fix: relabel `@param[out]`, note that load() clears inputs first, and reconcile the `id_data`/`peptide_identifications` name. Doc-only, ABI-safe.
- **Verified:** All literal facts verified. Header OMSSAXMLFile.h lines 48-49 tag the output containers `@param[in] protein_identification` and `@param[in] id_data`, and OMSSAXMLFile.cpp lines 32-33 reset them on entry (`protein_identification = ProteinIdentification(); peptide_identifications.clear();`). The `id_data` (header/Doxygen) vs `peptide_identifications` (.cpp definition) name mismatch is also real. The

### [FORM-55] PepNovoInfile::getModifications — "getModifications" returns via an out-parameter while every sibling Infile getter returns by const reference
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/PepNovoInfile.h · _format-id-search-in_

```cpp
void getModifications(std::map<std::string, std::string>& modification_key_map) const
```
- **Expectation:** Sibling classes expose modifications as `const std::map<...>& getModifications() const` (InspectInfile line 128, SequestInfile line 202). A caller of PepNovoInfile would expect the same return-by-const-ref shape.
- **Actual:** PepNovoInfile::getModifications takes a non-const out-parameter map and fills it, returning void — an inconsistent convention versus the other *Infile classes in the same module, and it also requires the caller to pre-declare a map.
- **Evidence:** Header (PepNovoInfile.h:65): `void getModifications(std::map<std::string, std::string> & modification_key_map) const;` vs InspectInfile.h:128 `const std::map<std::string, std::vector<std::string> >& getModifications() const;` and SequestInfile.h:202 same pattern.
- **Fix:** Add an overload (or new method) returning `const std::map<std::string,std::string>&` to match the sibling convention; keep the out-param version for ABI stability. Purely additive, ABI-safe.
- **Verifier correction:** The inconsistency is real but narrower than stated: PepNovoInfile::getModifications uses a void out-parameter while sibling InspectInfile/SequestInfile getters return by const ref. The return TYPE cannot be made identical (PepNovo's map<string,string> carries different semantics than the siblings' map<string,vector<string>>); only the out-param-vs-return calling convention is the genuine inconsistency. It is a low-severity convention nit (loud compile-time mismatch, no silent runtime impact). The recommended additive overload returning `const std::map<std::string,std::string>&` would be source-compatible (purely additive, existing symbol unchanged).
- **Verified:** Verified directly. PepNovoInfile.h:65 declares `void getModifications(std::map<std::string,std::string>& modification_key_map) const`, and the .cpp (line 163) just does `modification_key_map = mods_and_keys_;` — an out-parameter + void return. Its two sibling Infile getters in the same FORMAT module return by const ref: InspectInfile.h:128 and SequestInfile.h:202 both expose `const std::map<std::s

### [FORM-56] XTandemInfile::write — Serialization method is named write() here but store() in every sibling *Infile class
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/FORMAT/XTandemInfile.h · _format-id-search-in_

```cpp
void write(const std::string& filename, bool ignore_member_parameters = false, bool force_default_mods = false)
```
- **Expectation:** The other input-file adapters in this cluster serialize via `void store(const std::string&)` (InspectInfile, SequestInfile, PepNovoInfile). A developer moving between these classes expects `store()` on XTandemInfile too.
- **Actual:** XTandemInfile uses `write(...)` with two trailing bare-bool flags, breaking the naming convention of the sibling classes and presenting two same-typed bool defaults that are easy to transpose at the call site.
- **Evidence:** XTandemInfile.h:186 `void write(const std::string& filename, bool ignore_member_parameters = false, bool force_default_mods = false);` vs InspectInfile.h:47 `void store(const std::string& filename);`, SequestInfile.h:47, PepNovoInfile.h:50 (all `store`).
- **Fix:** Provide a `store(const std::string&)` alias for cross-class consistency (ABI-safe additive), and/or document the two bool flags clearly. Renaming write()->store() would be an ABI break, so prefer the alias.
- **Verifier correction:** The naming inconsistency (write vs store) is confirmed and genuine, but the title overstates impact. The method is fully documented at its declaration (XTandemInfile.h:175-187, including @param for both bools), and the practical call surface is now minimal: the former XTandemAdapter TOPP tool that called this has been removed, so the only remaining callers are the pyOpenMS binding (bind_misc.cpp:5075, using named keyword args) and the unit test (XTandemInfile_test.cpp:195). The transposition-of-two-bools concern is therefore a low-severity ergonomic wart, not a data-loss/silent-wrong-result hazard. Recommendation stands: add an additive `void store(const std::string&)` (ABI-safe) for cross-class consistency rather than renaming write() (which would be an ABI break).
- **Verified:** Evidence is accurate and verified independently. XTandemInfile.h:186-187 declares `void write(const std::string& filename, bool ignore_member_parameters = false, bool force_default_mods = false)`, while the sibling input-file adapters all serialize via `store()`: InspectInfile.h:47, SequestInfile.h:47, PepNovoInfile.h:50, and PercolatorInfile.h:29 (static `store`). So the naming divergence (write 

### [FORM-57] InspectOutfile::getWantedRecords — Doc says 'loads only results which exceeds a given P-value threshold' but it keeps records BELOW/at the threshold
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/FORMAT/InspectOutfile.h · _format-id-search-out_

```cpp
std::vector<Size> getWantedRecords(const std::string & result_filename, double p_value_threshold)
```
- **Expectation:** From the header brief ('loads only results which exceeds a given P-value threshold' / '@param p_value_threshold Only identifications exceeding this threshold are read') a caller expects records with p-value GREATER than the threshold to be kept.
- **Actual:** The implementation keeps records whose p-value is <= threshold and skips those that exceed it. InspectOutfile.cpp:1172-1175: comment '// take only those peptides whose p-value is less or equal the given threshold' followed by 'if (StringUtils::toFloat(substrings[p_value_column]) > p_value_threshold) { continue; }'. The same inverted-from-doc logic is in load() at lines 164-166.
- **Evidence:** InspectOutfile.cpp:1172  // take only those peptides whose p-value is less or equal the given threshold\nInspectOutfile.cpp:1173  if (StringUtils::toFloat(substrings[p_value_column]) > p_value_threshold)\nInspectOutfile.cpp:1175    continue;   (vs header brief 'loads only results which exceeds a given P-value threshold')
- **Fix:** Fix the Doxygen brief and the @param text to say records with p-value <= threshold are retained (lower p-value = more significant). Behavior is correct for p-values; only the documentation contradicts it. Documentation-only change: abi_impact none.
- **Verifier correction:** The header documentation for getWantedRecords (and the analogous load) is inverted relative to the code. Brief "loads only results which exceeds a given P-value threshold" and @param "Only identifications exceeding this threshold are read" imply p > threshold is kept, but the implementation keeps records with p-value <= threshold (skips those exceeding it; comment at InspectOutfile.cpp:1172 confirms "less or equal"). The behavior is domain-correct (lower p-value = more significant); only the documentation contradicts it. Documentation-only fix; abi_impact none.
- **Verified:** Independently verified against the actual code. Header (InspectOutfile.h:63,66): brief = "loads only results which exceeds a given P-value threshold"; @param = "Only identifications exceeding this threshold are read". Plain reading of "exceeds/exceeding a threshold" means p > threshold is kept. The implementation does the OPPOSITE: getWantedRecords (InspectOutfile.cpp:1172-1176) has the comment "t

### [FORM-59] OMSSACSVFile::load — load() appends to id_data instead of clearing it, unlike sibling Percolator load()
`low` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/FORMAT/OMSSACSVFile.h · _format-id-search-out_

```cpp
void load(const std::string & filename, ProteinIdentification & protein_identification, PeptideIdentificationList & id_data) const
```
- **Expectation:** A FORMAT load(file, data) that fills an output container is normally expected to replace the container's contents (PercolatorOutfile::load does 'peptides.clear();' first). A caller reusing the same id_data vector across calls expects fresh contents.
- **Actual:** OMSSACSVFile::load never clears id_data; it only emplace_back/insertHit, so contents accumulate across calls. OMSSACSVFile.cpp:91 'id_data.emplace_back();' and :96 'id_data.back().insertHit(p);' with no prior clear. Contrast PercolatorOutfile.cpp:205 'peptides.clear();'.
- **Evidence:** OMSSACSVFile.cpp:91  id_data.emplace_back();\nOMSSACSVFile.cpp:96  id_data.back().insertHit(p);\n(no id_data.clear() anywhere; cf. PercolatorOutfile.cpp:205  peptides.clear();)
- **Fix:** Make the clear-vs-append contract consistent across the *Outfile/*File loaders. Either clear id_data at the start of load() (matching PercolatorOutfile) or document the append semantics in the header. Adding a clear() is source-compatible behavior change; abi none.
- **Verifier correction:** OMSSACSVFile::load (OMSSACSVFile.cpp:91/96) appends to id_data without clearing it first, whereas its sibling PercolatorOutfile::load (PercolatorOutfile.cpp:205), taking the identical PeptideIdentificationList& output parameter, calls peptides.clear() before populating. The asymmetry is real and undocumented. Severity is low, not medium: the only caller (the pyOpenMS binding) and the typical OpenMS usage pass a fresh container, so the accumulation only manifests if a C++ caller explicitly reuses the same vector across load() calls — a recoverable, non-crashing, non-silent-corruption misuse with no current in-tree occurrence. Adding id_data.clear() at the start of load() is a source-compatible behavior change; ABI impact none.
- **Verified:** Evidence verified in the actual code. OMSSACSVFile::load (src/openms/source/FORMAT/OMSSACSVFile.cpp:24-99) never clears id_data: line 91 `id_data.emplace_back();` (guarded by a new-spectrum-number check) and line 96 `id_data.back().insertHit(p);`, with no clear() anywhere in the function. The cited sibling, PercolatorOutfile::load (PercolatorOutfile.cpp:170-205), takes the IDENTICAL output type `P

### [FORM-61] InspectOutfile::load — Doxygen tags the input filename as [out] and the actual output containers' direction is also mislabeled
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/InspectOutfile.h · _format-id-search-out_

```cpp
std::vector<Size> load(const std::string & result_filename, PeptideIdentificationList & peptide_identifications, ProteinIdentification & protein_identification, const double p_value_threshold, const std::string & database_filename = "")
```
- **Expectation:** The @param direction tags should tell callers which arguments are read and which are written, so a caller knows result_filename is input and peptide_identifications/protein_identification are outputs.
- **Actual:** The header marks the input filename as [out]: '@param[out] result_filename Input parameter which is the file name of the input file' (self-contradictory: tagged [out] but text says 'Input parameter'). The two genuine output containers are also tagged [out] but described as if input. The direction annotations are unreliable across the whole signature.
- **Evidence:** InspectOutfile.h:52  @param[out] result_filename Input parameter which is the file name of the input file\nInspectOutfile.h:53  @param[out] peptide_identifications Output parameter ...
- **Fix:** Correct the @param[in]/@param[out] tags: result_filename, p_value_threshold, database_filename are [in]; peptide_identifications, protein_identification are [out]. Documentation-only; abi none.
- **Verifier correction:** Only one tag is incorrect: `result_filename` is documented `@param[out]` but is an input (opened read-only via ifstream; the text on the same line even reads "Input parameter"). The fix is to change that one tag to `@param[in]`. The claim's assertion that the output containers (peptide_identifications, protein_identification) are also mislabeled is wrong — they are correctly tagged [out] and correctly described as outputs in the .cpp. The self-contradiction "[out] ... Input parameter" is visible in the same line, so it loudly flags itself rather than silently misleading a caller; hence low severity, not high.
- **Verified:** Verified against InspectOutfile.h:50-61 and InspectOutfile.cpp:55-79. The quoted evidence is real: line 52 tags the input filename as `@param[out] result_filename Input parameter which is the file name of the input file` — self-contradictory (tagged [out] but text says "Input parameter"). The .cpp confirms result_filename is opened via ifstream and only read (line 64), so it is genuinely an input.

### [FORM-63] PercolatorOutfile::load — load() mutates the caller's SpectrumMetaDataLookup by registering reference-format regexes
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/PercolatorOutfile.h · _format-id-search-out_

```cpp
void load(const std::string& filename, ProteinIdentification& proteins, PeptideIdentificationList& peptides, SpectrumMetaDataLookup& lookup, ScoreType output_score = ScoreType::QVALUE)
```
- **Expectation:** A load() taking a SpectrumMetaDataLookup& would be read as 'use this lookup to resolve spectra'; a caller does not expect load() to permanently add regex reference formats to their lookup object.
- **Actual:** If lookup.reference_formats is empty, load() injects three hard-coded regexes into the caller's lookup (MS-GF+, Mascot, X! Tandem formats) that persist after the call. PercolatorOutfile.cpp:181-190: 'if (lookup.reference_formats.empty()) { lookup.addReferenceFormat(...); ... }'. The header gives no hint that lookup is modified.
- **Evidence:** PercolatorOutfile.cpp:181  if (lookup.reference_formats.empty())\nPercolatorOutfile.cpp:184  lookup.addReferenceFormat(R"(_SII_...)");  (plus two more)
- **Fix:** Document in the header that load() seeds default reference formats into lookup when it has none, so callers reusing the lookup across formats understand the persistent mutation. Documentation-only; abi none.
- **Verifier correction:** load() seeds default reference-format regexes into the passed SpectrumMetaDataLookup only when it has none (reference_formats.empty()), so the formats persist after the call. This is signaled by the non-const reference parameter and uses the documented addReferenceFormat API; it is skipped entirely if the caller already registered any format. The real surprise is limited to reusing one lookup object across multiple loaders; document the seeding in the header. Not a high/medium-impact issue — additive, conditional, recoverable, and functionally required for spectrum-reference resolution.
- **Verified:** Evidence verified verbatim: PercolatorOutfile.cpp:181-190 conditionally injects three hard-coded reference-format regexes (MS-GF+, Mascot, X! Tandem) into the caller's lookup via addReferenceFormat when lookup.reference_formats is empty. reference_formats is a public member (SpectrumLookup.h:64) and addReferenceFormat (decl SpectrumLookup.h:187) permanently appends, so the mutation persists and is

### [FORM-64] PercolatorInfile::load — extra_scores is documented as @param[out] but is read-only input (const&), misleading the direction contract
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/PercolatorInfile.h · _format-id-search-out_

```cpp
static PeptideIdentificationList load(const std::string& pin_file, bool higher_score_better, const std::string& score_name, const StringList& extra_scores, StringList& filenames, std::string decoy_prefix = "", double threshold = 0.01, bool SageAnnotation = false)
```
- **Expectation:** An @param[out] tag tells callers the function writes results into that argument; a caller may pass an empty StringList expecting load() to fill the list of extra score names.
- **Actual:** extra_scores is 'const StringList& extra_scores' and is only iterated as input (PercolatorInfile.cpp:175 'for (const std::string& s : extra_scores)'). It is the NAMES the caller wants extracted, never written. Only 'filenames' is a true out-parameter (push_back at :249). The [out] tag on extra_scores is wrong.
- **Evidence:** PercolatorInfile.h:51  @param[out] extra_scores A list of additional score names that should be extracted ...\nPercolatorInfile.cpp:78  const StringList& extra_scores,\nPercolatorInfile.cpp:175  for (const std::string& s : extra_scores)
- **Fix:** Change the Doxygen tag for extra_scores to @param[in]. Documentation-only; the const& already protects the caller. abi none.
- **Verifier correction:** extra_scores is the input list of additional score NAMES the caller wants extracted (read-only, const&, iterated at cpp:175); it is never written. It is mis-tagged @param[out] in the header (h:51) and should be @param[in]. The genuine out-parameter is `filenames` (cpp:249). Doc-only fix; const& already prevents any caller-side corruption, so the surprise is mild rather than dangerous.
- **Verified:** Verified directly in source. PercolatorInfile.h:51 tags `@param[out] extra_scores`, but the signature (h:67, cpp:78) declares it `const StringList& extra_scores`, and cpp:175 only reads it: `for (const std::string& s : extra_scores)` — these are the caller-supplied score NAMES to look up in the pin header (cpp:177), with matches copied into the local `found_extra_scores` set (cpp:179). It is never

### [FORM-65] InspectOutfile::getExperiment — getExperiment is a heavyweight I/O + file-type-detection routine, not a cheap accessor as 'get' implies, and resets its out-params
`low` · `naming-vs-behavior` · ABI: `none` · src/openms/include/OpenMS/FORMAT/InspectOutfile.h · _format-id-search-out_

```cpp
void getExperiment(PeakMap & exp, std::string &type, const std::string & in_filename)
```
- **Expectation:** A method named getExperiment(exp, type, in_filename) reads like a lightweight retrieval; a caller might not expect content-sniffing I/O and might assume passing pre-filled exp/type is harmless.
- **Actual:** The inline body clears both out-parameters (type.clear(); exp.reset();), then sniffs file type from content (FileHandler::getTypeByContent), throws ParseError on UNKNOWN, and performs a full experiment load (fh.loadExperiment). It is a parse-and-load operation, not a getter, and silently discards whatever was in exp/type. See InspectOutfile.h:116-129.
- **Evidence:** InspectOutfile.h:118  type.clear();\nInspectOutfile.h:119  exp.reset();\nInspectOutfile.h:122  FileTypes::Type in_type = fh.getTypeByContent(in_filename);\nInspectOutfile.h:128  fh.loadExperiment(in_filename, exp, {in_type});
- **Fix:** Consider renaming to loadExperiment()/readExperiment() to signal I/O, or at minimum document that it performs content-based type detection and a full file load and resets the out-parameters. A rename is breaking; doc + a deprecated alias is the additive path. Documentation-only fix abi none.
- **Verifier correction:** getExperiment(exp, type, in_filename) does perform content-based file-type detection (getTypeByContent) and a full file load (loadExperiment), throwing ParseError on an undeterminable type — so the 'get' prefix understates that it is I/O. But it is documented as reading 'from a file' with a documented ParseError, the out-param reset (type.clear/exp.reset) is the normal output-parameter idiom rather than data loss, and the method fails loudly rather than silently. The defensible surprise is purely naming/altitude (get vs load/read), low severity. Additive fix: clarify the doc comment that it sniffs type and loads the whole file; optionally add a loadExperiment()/readExperiment() alias (rename would be breaking). Doc/alias is ABI-neutral.
- **Verified:** Evidence is accurate: InspectOutfile.h:116-129 inline body does type.clear() (118), exp.reset() (119), content-sniff via FileHandler::getTypeByContent (122), throw ParseError on UNKNOWN (123-126), and a full fh.loadExperiment (128). So it is genuinely a parse-and-load, not a cheap accessor — the 'get' prefix is misleading. However the claim is overstated on two points. (1) 'Hidden': the doc commen

### [FORM-43] ProtXMLFile::store — ProtXMLFile::store is a fully declared public method whose only behavior is to throw NotImplemented
`low` · `documentation-mismatch` · ABI: `none` · src/openms/include/OpenMS/FORMAT/ProtXMLFile.h · _format-id-xml_

```cpp
void store(const std::string& filename, const ProteinIdentification& protein_ids, const PeptideIdentification& peptide_ids, const std::string& document_id = "")
```
- **Expectation:** The Doxygen says '@brief [not implemented yet!] Stores the data in an ProtXML file' and '@exception Exception::UnableToCreateFile is thrown if the file could not be created'. A caller scanning the header sees a normal store() with an UnableToCreateFile contract and would reasonably attempt to use it (it is part of the standard load/store pair).
- **Actual:** The implementation unconditionally throws Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION) for ANY input, regardless of filename validity. The documented UnableToCreateFile exception is never thrown; the real exception type (NotImplemented) is not documented in the @exception list.
- **Evidence:** ProtXMLFile.cpp: 'void ProtXMLFile::store(...) { throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }'. Header @exception only lists UnableToCreateFile.
- **Fix:** Keep throwing (the feature is unimplemented), but fix the doc: replace the @exception UnableToCreateFile line with '@exception Exception::NotImplemented is always thrown (writing protXML is not supported)'. Ideally mark the method [[deprecated("protXML writing not implemented")]] or remove it until implemented; removal is breaking, so the doc fix is the ABI-safe path.
- **Verifier correction:** ProtXMLFile::store() is an unimplemented stub that unconditionally throws Exception::NotImplemented for any input (loud, fail-fast — not silent). Its unimplemented status is clearly and repeatedly documented in the same header (@brief '[not implemented yet!]' x2, class doc 'storing not supported, yet', @todo). The genuine (minor) defect is only that the @exception line documents Exception::UnableToCreateFile while the code actually throws Exception::NotImplemented. Fix: replace the @exception line with '@exception Exception::NotImplemented is always thrown (writing protXML is not supported)'. ABI-safe, doc-only.
- **Verified:** Code confirmed: ProtXMLFile::store (header line 72) is a fully-declared public method, and the impl (cpp lines 45-49) unconditionally `throw Exception::NotImplemented(...)` for any input. The header @exception (line 70) documents `UnableToCreateFile`, which is never thrown — a genuine doc/contract mismatch (NotImplemented and UnableToCreateFile are confirmed distinct classes in Exception.h). Howev

### [FORM-44] PepXMLFile::store — PepXMLFile::store takes protein_ids/peptide_ids by NON-const reference but never modifies them; inconsistent with sibling store() methods
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/PepXMLFile.h · _format-id-xml_

```cpp
void store(const std::string& filename, std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids, const std::string& mz_file = "", const std::string& mz_name = "", bool peptideprophet_analyzed = false, double rt_tolerance = 0.01)
```
- **Expectation:** store() writes data to a file and must not modify its inputs; callers expect to be able to store a const collection. Every sibling does: IdXMLFile::store, MzIdentMLFile::store, ProtXMLFile::store all take 'const std::vector<ProteinIdentification>&' / 'const ...&'. A non-const ref here signals 'I may mutate your data' and prevents passing const/temporary collections.
- **Actual:** The whole store() body only reads protein_ids/peptide_ids: 'protein_ids.begin()->getSearchParameters()', const_iterator loops over peptide_ids, and 'for (const PeptideIdentification& pep : peptide_ids)'. Nothing is mutated. The non-const reference is gratuitous and breaks the const-store convention of the module.
- **Evidence:** PepXMLFile.cpp store body uses only reads: 'search_params = protein_ids.begin()->getSearchParameters();', 'for (PeptideIdentificationList::const_iterator it = peptide_ids.begin(); ...)', 'for (const PeptideIdentification& pep : peptide_ids)'. No assignment/sort/push on either parameter across the full body (lines 351-922).
- **Fix:** Change the two parameters to const references to match IdXMLFile/MzIdentMLFile/ProtXMLFile. This is source-compatible for callers passing non-const lvalues but is technically an ABI/signature change (mangled name changes). If strict ABI must be preserved, add a new const overload and forward; deprecate the non-const one.
- **Verifier correction:** The finding is correct that PepXMLFile::store takes protein_ids/peptide_ids by non-const reference while never modifying them and while every sibling (IdXMLFile, MzIdentMLFile, ProtXMLFile) uses const references. The severity should be 'low' rather than implied higher: there is no incorrect output, data loss, or crash — only an interface-quality inconsistency that prevents passing const/temporary collections. ABI impact is 'breaking' (mangled name changes) though source-compatible for all existing callers, which pass non-const lvalues. The safe minimal fix is to change both parameters to const references; if strict ABI must be preserved, add a const overload and forward.
- **Verified:** Independently verified. Header (PepXMLFile.h:84-86) declares store() with std::vector<ProteinIdentification>& and PeptideIdentificationList& by NON-const reference. Read the full .cpp body (lines 351-922): no mutation of either parameter occurs — only reads via .empty()/.size(), 'search_params = protein_ids.begin()->getSearchParameters();' (copied into a local value), const_iterator loops over pep

### [FORM-45] PTMXMLFile::store — PTMXMLFile::store is a const method but takes ptm_informations by NON-const reference, falsely implying it may mutate the map
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/PTMXMLFile.h · _format-id-xml_

```cpp
void store(const std::string& filename, std::map<std::string, std::pair<std::string, std::string> >& ptm_informations) const
```
- **Expectation:** A const store() that 'Stores the data in an PTMXML file' should not modify the data it stores, and a caller should be able to store a const map. The non-const reference parameter signals possible mutation and blocks const/temporary arguments. Sibling store() methods (IdXMLFile, MzIdentMLFile) take their payload by const reference.
- **Actual:** store() forwards the map to PTMXMLHandler, whose writeTo() only reads ptm_informations_ to emit XML; nothing in the store path mutates the map. The non-const reference is gratuitous and inconsistent with the load() counterpart's intent (load needs non-const to fill; store does not).
- **Evidence:** PTMXMLFile.cpp: 'void PTMXMLFile::store(const std::string& filename, map<...>& ptm_informations) const { Internal::PTMXMLHandler handler(ptm_informations, filename); save_(filename, &handler); }'. Header doc says only '@throw UnableToCreateFile' and 'The data is read in and stored'.
- **Fix:** Change the parameter to 'const std::map<...>&'. Note PTMXMLHandler's ctor also takes a non-const ref (shared load/store ctor), so this requires either a const ctor overload or a const_cast at the boundary. Source-compatible for existing callers but changes the mangled signature (ABI). If ABI-frozen, add a const-map overload.
- **Verified:** Verified directly. PTMXMLFile.h:51 declares `void store(const std::string&, std::map<std::string, std::pair<std::string,std::string>>&) const` — a const method taking the payload by NON-const reference. PTMXMLFile.cpp:27-31 forwards the map to PTMXMLHandler's ctor and calls save_(); PTMXMLHandler.cpp writeTo() iterates ptm_informations_ via a const_iterator and only emits XML — nothing in the stor

### [FORM-46] PepXMLFile::store — PepXMLFile::store doc says 'Stores idXML as PepXML file' but the method writes pepXML; @exception text mentions wrong condition
`low` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/FORMAT/PepXMLFile.h · _format-id-xml_

```cpp
void store(const std::string& filename, ... double rt_tolerance = 0.01)
```
- **Expectation:** The brief should describe what store() does. A reader expects 'Stores identifications as a pepXML file'.
- **Actual:** The @brief reads 'Stores idXML as PepXML file' (it does not consume an idXML file; it serializes in-memory ProteinIdentification/PeptideIdentification to pepXML). The @exception says 'UnableToCreateFile is thrown if the file could not be opened for writing', which is accurate, but the brief conflates the in-memory data model with the idXML on-disk format and can mislead a caller into thinking a path/conversion of an existing idXML file is involved.
- **Evidence:** PepXMLFile.h: '@brief Stores idXML as PepXML file'. Implementation serializes the passed protein_ids/peptide_ids vectors directly to a pepXML stream; no idXML file is read.
- **Fix:** Doc-only, ABI-safe: reword the brief to 'Stores the given identifications as a pepXML file'. No signature change.
- **Verifier correction:** The '@brief Stores idXML as PepXML file' on PepXMLFile::store (PepXMLFile.h:80) is genuinely misleading: store() serializes the in-memory ProteinIdentification/PeptideIdentification containers passed by reference to a pepXML output stream (PepXMLFile.cpp:351-357) and never reads an idXML file. Because OpenMS does have a real idXML->pepXML file-conversion concept (IDFileConverter), the wording invites a wrong reading. But severity is low, not the implied higher confusion: the parameter types (filename as output + in-memory ID vectors as input) make the intent unambiguous on inspection and prevent any actual misuse. Doc-only, ABI-safe fix: reword to 'Stores the given identifications as a pepXML file', matching the correct phrasing already used in MzIdentMLFile.h:60.
- **Verified:** Evidence is verbatim accurate. PepXMLFile.h line 80 reads '@brief Stores idXML as PepXML file'. The implementation (PepXMLFile.cpp:351-357) opens an ofstream on the given filename and serializes the in-memory protein_ids/peptide_ids vectors to pepXML — it never reads an idXML file. The brief is genuinely anomalous, not a domain convention: sibling classes phrase it correctly (IdXMLFile.h:77 'Store

### [FORM-47] ProtXMLFile::matchModification_ — matchModification_ doc tags 'mass' as @param[in,out] and the description param is documented [out] in the wrong slot
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/FORMAT/ProtXMLFile.h · _format-id-xml_

```cpp
void matchModification_(const double mass, const std::string& origin, std::string& modification_description)
```
- **Expectation:** Doxygen @param[in,out] mass implies the function reads and writes mass. The third parameter is the output (the modification description).
- **Actual:** The signature declares 'const double mass' (read-only, cannot be modified), so '@param[in,out] mass' is wrong; mass is purely [in]. The doc also annotates the wrong parameter: it lists '@param[in] modification_description [out] Name of the modification' -- modification_description is the actual output (taken by non-const std::string&), but is tagged @param[in]. The implementation only assigns to modification_description.
- **Evidence:** ProtXMLFile.h doc: '@param[in,out] mass Modified AA's mass' and '@param[in] modification_description [out] Name of the modification'. ProtXMLFile.cpp: 'void ProtXMLFile::matchModification_(const double mass, ...) { ... modification_description = mods[0]; }' (mass never assigned).
- **Fix:** Doc-only, ABI-safe: change to '@param[in] mass' and '@param[out] modification_description'. This is a protected helper, so impact is limited to maintainers, but the inverted in/out tags are actively misleading.
- **Verified:** Verified directly. Header (ProtXMLFile.h:96-100) doc reads '@param[in,out] mass Modified AA's mass' and '@param[in] modification_description [out] Name of the modification'. The signature is 'void matchModification_(const double mass, const std::string& origin, std::string& modification_description)'. 'mass' is const double — provably read-only — so '@param[in,out]' is self-contradictory; the impl

### [FORM-23] IndexedMzMLFileLoader::store — store() takes its source map by non-const reference (cannot store a const map)
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/IndexedMzMLFileLoader.h · _format-mzml_

```cpp
void store(const std::string& filename, OnDiscPeakMap& exp); void store(const std::string& filename, PeakMap& exp);
```
- **Expectation:** store() is a write-out operation; every other *File class declares it const-taking: `store(const std::string&, const PeakMap&) const`. A caller expects to be able to store a `const PeakMap&`.
- **Actual:** Both store() overloads take a non-const `OnDiscPeakMap& exp` / `PeakMap& exp` (and the methods themselves are non-const because they mutate options_). The non-const PeakMap& in particular is gratuitous — the impl only reads from it via f.store(). This breaks the const-correct call sites that work for MzMLFile.
- **Evidence:** Header: `void store(const std::string& filename, OnDiscPeakMap& exp);` and `void store(const std::string& filename, PeakMap& exp);` (note doc tag `@param[out] exp MS data to be stored` — labeling a store *source* as [out] is itself wrong). Impl `store(..., PeakMap& exp)` only reads exp: `f.store(filename, exp);`
- **Fix:** Add a `store(const std::string&, const PeakMap&)` overload (additive, source-compatible) so const maps can be stored consistently with the other File classes; fix the `@param[out]` doc tags to `@param[in]`.
- **Verifier correction:** Narrowed claim: Only the in-memory overload's parameter is gratuitously non-const — store(const std::string&, PeakMap& exp) reads exp solely to forward to MzMLFile::store(..., const PeakMap&), so it should take const PeakMap&. This means a const PeakMap cannot be stored, unlike MzMLFile/MzXMLFile/MzDataFile which all take const MapType&. The OnDiscPeakMap& overload is NOT gratuitous (it calls non-const OnDiscMSExperiment::getSpectrum/getChromatogram), and the methods cannot simply be made const because they mutate options_ (setWriteIndex). Recommended fix is additive and source-compatible: change/overload the PeakMap parameter to const, and fix the copy-pasted @param[out] doc tags (the exp source should be @param[in]). Severity is low: a const map fails loudly at compile time (no silent data corruption or data loss), and the misuse is fully recoverable.
- **Verified:** Evidence is confirmed and the core surprise is real, but the claim overreaches and must be narrowed. CONFIRMED: header lines 65/73 declare both store() overloads with non-const refs; the PeakMap& overload (cpp lines 62-68) only reads exp (forwards to MzMLFile::store which takes const PeakMap&); sibling File classes (MzMLFile.h:90, MzXMLFile.h:63, MzDataFile.h:66) all declare store(const std::strin

### [FORM-24] IndexedMzMLFileLoader::store / setOptions — store() silently overwrites the user-configured WriteIndex option on the object
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/IndexedMzMLFileLoader.h · _format-mzml_

```cpp
void store(const std::string& filename, PeakMap& exp); void setOptions(const PeakFileOptions&)
```
- **Expectation:** A caller who does setOptions(opts) expects the options object to reflect exactly what they set; store() should not mutate the persistent options_.
- **Actual:** Both store() overloads call `options_.setWriteIndex(true)` on the member options before storing, permanently flipping the user's WriteIndex flag as a side effect of the store call (it is never restored). A subsequent getOptions() returns a different value than what was setOptions().
- **Evidence:** Impl: `void IndexedMzMLFileLoader::store(...) { ... options_.setWriteIndex(true);  // ensure that we write the index ... }` (in both overloads), with no save/restore.
- **Fix:** Apply WriteIndex(true) to a local copy of the options passed to the consumer/MzMLFile rather than mutating the persistent member options_ (source- and ABI-compatible internal change).
- **Verifier correction:** store() does mutate the persistent options_ member by forcing setWriteIndex(true) and never restoring it, so getOptions() no longer round-trips a prior setOptions() — confirmed. But the practical impact is minimal, not high: WriteIndex already defaults to true, the override is semantically correct for an indexed-mzML writer (a false value is contradictory and rightly ignored), no output is ever wrong, and no caller reads the flag back. Recommended fix (apply WriteIndex(true) to a local copy passed to the consumer/MzMLFile instead of mutating options_) is correct and ABI-/source-compatible.
- **Verified:** Evidence confirmed verbatim: both IndexedMzMLFileLoader::store() overloads (IndexedMzMLFileLoader.cpp lines 48 and 65) call options_.setWriteIndex(true) on the persistent member options_ with no save/restore. So getOptions().getWriteIndex() after a store() returns true regardless of a prior setOptions() that set it false — the round-trip violation is factually real and not a documented/idiomatic/i

### [FORM-25] CachedmzML::CachedmzML(const CachedmzML&) — Copy constructor silently drops the stream position and does not copy chrom_index/filename_cached state
`low` · `ownership-lifetime` · ABI: `none` · src/openms/include/OpenMS/FORMAT/CachedMzML.h · _format-mzml_

```cpp
CachedmzML(const CachedmzML & rhs)
```
- **Expectation:** Copying a fully-loaded CachedmzML should yield a copy in which getSpectrum/getChromatogram work identically to the original.
- **Actual:** The copy ctor opens a fresh ifstream from `rhs.filename_cached_` but its initializer list copies `filename_` (not `filename_cached_`) and `spectra_index_`/`chrom_index_` — yet it never copies `filename_cached_` itself. So the copy's `filename_cached_` is left empty even though it reads rhs.filename_cached_ to open the stream; any later error message / re-open via filename_cached_ on the copy is wrong/empty. Stream open errors are also not checked.
- **Evidence:** Impl: `CachedmzML::CachedmzML(const CachedmzML & rhs) : meta_ms_experiment_(rhs.meta_ms_experiment_), ifs_(rhs.filename_cached_.c_str(), std::ios::binary), filename_(rhs.filename_), spectra_index_(rhs.spectra_index_), chrom_index_(rhs.chrom_index_) { }` — `filename_cached_` is not in the member-init list, so the copy's filename_cached_ stays default-empty.
- **Fix:** Add `filename_cached_(rhs.filename_cached_)` to the initializer list and verify the freshly-opened stream's state (additive/internal; source- and ABI-compatible).
- **Verifier correction:** The copy constructor omits filename_cached_ from its member-init list, so the copy's filename_cached_ is default-empty even though ifs_ is correctly opened from rhs.filename_cached_, and spectra_index_/chrom_index_/filename_/meta_ms_experiment_ are all correctly copied. Consequence: getSpectrum/getChromatogram work identically on the copy in the normal path (the stated 'breaks identical behavior' expectation is NOT violated for success cases). The only defect is that, if seekg later fails (rare, e.g. >2GB cache on 32-bit systems), the thrown Exception::ParseError on a copied instance carries an empty filename field instead of the cached path (lines 71/88), degrading the diagnostic. The claim's 'drops stream position', 'does not copy chrom_index_/filename', and 'later re-open via filename_cached_' assertions are incorrect. Fix is still valid: add filename_cached_(rhs.filename_cached_) to the init list (additive, ABI-safe).
- **Verified:** Code confirmed at src/openms/source/FORMAT/CachedMzML.cpp:35-42. The copy ctor's init list genuinely omits filename_cached_: it opens ifs_ from rhs.filename_cached_ (line 37) and copies meta_ms_experiment_, filename_, spectra_index_, chrom_index_, but filename_cached_ (declared at header line 180) is left default-empty in the copy. So the bug — a member silently not copied — is real and provable. 

### [FORM-27] MzMLFile::transform vs MzXMLFile::transform — transform() signature drifts between MzMLFile and MzXMLFile (extra skip_first_pass bool)
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzMLFile.h · _format-mzml_

```cpp
MzMLFile::transform(const std::string&, IMSDataConsumer*, bool, bool); MzXMLFile::transform(const std::string&, IMSDataConsumer*, bool)
```
- **Expectation:** Sibling File classes presenting the 'same interface' (as IndexedMzMLFileLoader's doc claims for the family) should expose matching transform() signatures so call sites are portable.
- **Actual:** MzMLFile::transform has a 4th trailing parameter `bool skip_first_pass = false` that MzXMLFile::transform lacks. Code written against one will not compile against the other, and the two trailing bools (skip_full_count, skip_first_pass) are same-typed and easily swapped at call sites.
- **Evidence:** MzMLFile.h: `void transform(const std::string& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count = false, bool skip_first_pass = false);` vs MzXMLFile.h: `void transform(const std::string& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count = false);`
- **Fix:** Document the divergence; if parity is desired, add the skip_first_pass overload to MzXMLFile (additive, source-compatible). The two adjacent bools are inherently swap-prone — consider an enum/flags type in any new overload.
- **Verifier correction:** transform() signature drifts between MzMLFile and MzXMLFile: MzMLFile adds a 4th trailing `bool skip_first_pass = false` (also on its map-taking overload) that MzXMLFile lacks. The two classes share no virtual/base transform() contract — parity is only an informal, documented family convention (IndexedMzMLFileLoader.h:22). Consequence is a COMPILE-TIME portability break (a 4-arg MzMLFile call won't bind to MzXMLFile), not silently wrong results; severity is low. Secondary, milder hazard: skip_full_count and skip_first_pass are adjacent same-typed bools, swap-prone at MzMLFile call sites. Fix is additive/source-compatible (add the skip_first_pass overload to MzXMLFile, ideally via an enum/flags type).
- **Verified:** Evidence verified verbatim. MzMLFile.h:119 declares `void transform(const std::string& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count = false, bool skip_first_pass = false);` while MzXMLFile.h:83 declares only `void transform(const std::string& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count = false);` — a real 4th-param divergence (mirrored in 

### [FORM-29] Internal::XMLFile::isValid — isValid() throws NotImplemented on a default-constructed instance instead of returning a bool
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/FORMAT/XMLFile.h · _format-mzml_

```cpp
bool isValid(const std::string& filename, std::ostream& os)
```
- **Expectation:** A predicate named isValid returning bool is expected to answer true/false (or at most throw FileNotFound for a missing file), not throw because the *object itself* has no schema bound.
- **Actual:** If the XMLFile was default-constructed (schema_location_ empty), isValid throws Exception::NotImplemented rather than returning false or throwing a file-related error. The throw condition depends on object construction state, not the file argument, which is surprising for a bool query.
- **Evidence:** Header doc: '@throws OpenMS::Exception::NotImplemented if no schema is bound (default-constructed instance — schema_location_ is empty).' and constructor doc: 'so isValid cannot be used until derived-class logic initializes schema_location_'.
- **Fix:** Document is sufficient for ABI stability, but the throw-from-a-bool-predicate based on construction state is the surprise; consider asserting/documenting that derived classes must bind a schema. ABI-safe (doc-only).
- **Verifier correction:** isValid() does throw Exception::NotImplemented on a default-constructed (schema-less) instance based on construction state rather than the file argument — verified in src/openms/source/FORMAT/XMLFile.cpp:391-393 and locked in by src/tests/class_tests/openms/source/XMLFile_test.cpp:46. But this is a low-severity surprise, not high: XMLFile lives in the OpenMS::Internal namespace and is never used directly for validation; every real caller uses a schema-binding derived class, the error is loud/immediate (a thrown exception, no silent corruption), and NotImplemented is OpenMS's standard sentinel for an unconfigured base path. The behavior is also already documented at the throw site (XMLFile.h:67/42/127).
- **Verified:** Code confirmed: XMLFile.cpp:389-397 — isValid() throws Exception::NotImplemented when schema_location_ is empty, i.e. on a default-constructed instance, before ever touching the file argument. The throw depends on object construction state, not on the file, which is the stated surprise and is accurate. The header doc (XMLFile.h:67, also :42 and :127) documents this exact throw, and XMLFile_test.cp

### [FORM-30] MzMLFile::getCentroidInfo — getCentroidInfo() (a 'get' method) opens/parses a file and temporarily mutates the object's options
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzMLFile.h · _format-mzml_

```cpp
std::map<UInt, SpecInfo> getCentroidInfo(const std::string& filename, const Size first_n_spectra_only = 10)
```
- **Expectation:** A `get...Info` method named like an accessor is expected to be cheap/read-only and not perform file I/O or change the object's configuration.
- **Actual:** getCentroidInfo parses the file (via transform) and toggles `options_.setFillData(true)` during the call, restoring it afterward. It is also non-const. The 'get' name hides a full (partial) file parse plus transient mutation of the persistent PeakFileOptions member.
- **Evidence:** Impl: `bool oldoption = options_.getFillData(); options_.setFillData(true); ... transform(filename, &c, true, true); ... options_.setFillData(oldoption);` — mutates and restores the member options; performs a parse.
- **Fix:** Doc-only / naming: the I/O cost is implied by the filename arg, but the transient options_ mutation is non-obvious — note it is not thread-safe w.r.t. concurrent option reads. ABI-safe.
- **Verifier correction:** getCentroidInfo's file I/O and non-const-ness are NOT surprising — they are the OpenMS *File-class convention (load/store/transform are likewise non-const I/O), the filename argument signposts the parse, and the header docstring documents it. The defensible surprise is narrower: the method transiently mutates the persistent `options_.setFillData(true)` member and restores it only on the success path, with no RAII guard. This makes it (1) not thread-safe versus a concurrent getOptions() reader, and (2) not exception-safe — if transform/safeParse_ throws (FileNotFound or ParseError, both reachable) the restore is skipped and options_ is left permanently flipped. On the normal path the effect is fully invisible (restored). Recommendation: wrap the option in an RAII/scope-guard restore so the state survives exceptions; document non-thread-safety. Doc/impl-only, ABI-safe.
- **Verified:** Verified against the actual code. Header line 189: `std::map<UInt, SpecInfo> getCentroidInfo(const std::string& filename, const Size first_n_spectra_only = 10)` — non-const, confirmed. Source MzMLFile.cpp lines 233-274 match the quoted evidence exactly: `oldoption = options_.getFillData(); options_.setFillData(true); ... transform(filename, &c, true, true); ... options_.setFillData(oldoption);`. `

### [FORM-70] MzTabFile::load — load() overwrites/replaces the passed MzTab rather than the more cautious clear-then-fill being explicit; pre-populated sections are clobbered or stale
`low` · `incomplete-parse` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzTabFile.h · _format-mztab_

```cpp
void load(const std::string& filename, MzTab& mz_tab)
```
- **Expectation:** Conventionally a load(file, out) clears `out` and fills it from the file. A caller might reasonably pre-set some sections (e.g. NUC rows) and expect load to leave unparsed sections untouched or to fully reset the object.
- **Actual:** load() unconditionally overwrites exactly five sections via set...() but leaves any pre-existing nucleic-acid/oligo/OSM rows in `mz_tab` in place (because it never touches them). So if the same MzTab object is reused, the result is an inconsistent mix of freshly-loaded sections and stale leftover sections - neither a clean replace nor an append.
- **Evidence:** MzTabFile.cpp:1548-1554 (the only writes back into mz_tab) cover 5 sections + empty/comment rows; the NUC/OLI/OSM members of mz_tab are never reset, so whatever was there before load() survives.
- **Fix:** Document that load() requires a fresh MzTab, or have load() reset all sections (assign empty NUC/OLI/OSM) before filling. Behavior fix in the .cpp; no ABI change.
- **Verifier correction:** load() does not clear NUC/OLI/OSM sections, but the practically relevant defect is that the loader never parses nucleic-acid/oligonucleotide/OSM (NUH/NUC, OLH/OLI, OSH/OSM) section lines at all, so such data in an input file is silently dropped even for a fresh MzTab. The claimed "stale leftover on object reuse" outcome is technically true but unreachable in the codebase: all callers pass a freshly default-constructed (empty) MzTab, and the 5 parsed sections are cleanly replaced via set...(). Recommendation reduces to: document that load() is a partial parser ignoring oligonucleotide sections (or implement those sections), rather than a data-loss-on-reuse hazard. No ABI change.
- **Verified:** Evidence at MzTabFile.cpp:1548-1554 is accurate: load() writes back exactly 5 sections (MetaData, Protein, Peptide, PSM, SmallMolecule) plus empty/comment rows via plain-assignment set...() calls, and never touches the NUC/OLI/OSM members of MzTab. I confirmed by reading the whole load body (lines 96-1555): the section dispatch only handles COM, MTD, PRH/PRT, PEH/PEP, PSH/PSM, SMH/SML (greps at 23

### [FORM-71] MzTabMFile — MzTabMFile provides store() but no load() - a write-only *File adapter
`low` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzTabMFile.h · _format-mztab_

```cpp
class MzTabMFile { void store(const std::string&, const MzTabM&) const; }
```
- **Expectation:** By analogy with MzTabFile (which has both load and store) and the general OpenMS *File convention, a class named MzTabMFile that has store() would also expose load() to read MzTab-M files back.
- **Actual:** MzTabMFile declares only store() (and protected header/section generators). There is no load() anywhere in the header or .cpp. The adapter is write-only, surprising for anyone expecting load/store symmetry from the sibling MzTabFile.
- **Evidence:** MzTabMFile.h:33 `void store(const std::string& filename, const MzTabM& mztab_m) const;` is the only public I/O method; grep for 'load' in MzTabMFile.cpp finds none. Contrast MzTabFile.h:75 `void load(const std::string& filename, MzTab& mz_tab);`.
- **Fix:** Document the write-only nature in the class brief, or add a load() to match MzTabFile. Purely additive; no ABI impact.
- **Verifier correction:** MzTabMFile is intentionally a write-only/export-only adapter: MzTab-M in OpenMS is generated (MzTabM::exportFeatureMapToMzTabM) and stored, never parsed back, so the absence of load() is by design rather than an oversight. The genuine, low-severity surprise is the undocumented asymmetry versus the identically-named/briefed sibling MzTabFile, which does provide load(). Recommend documenting the write-only nature in the class brief; adding load() is optional. Severity is low because the gap is compile-time obvious (no runtime hazard, no data corruption).
- **Verified:** Code claim is factually correct. MzTabMFile.h exposes only public store() (line 33); no load() exists in the header or in src/openms/source/FORMAT/MzTabMFile.cpp (grep returns zero 'load' hits). The sibling MzTabFile has both store() and load() (MzTabFile.h:75, implemented at MzTabFile.cpp:96), and shares the near-identical brief "File adapter for MzTab[-M] files", so the asymmetry is a genuine, u

### [FORM-72] MzTabSpectraRef::setMSFile — setMSFile silently ignores index 0 instead of signaling the invalid 1-based value
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzTabBase.h · _format-mztab_

```cpp
void MzTabSpectraRef::setMSFile(Size index)
```
- **Expectation:** Setting the MS-run index to a value should either store it or report an error. The parameter is a 1-based mzTab ms_run index (toCellString emits `ms_run[index]`).
- **Actual:** For index < 1 (i.e. 0) the method does nothing at all: it neither stores nor warns nor throws (only a debug-build assert). The MzTabSpectraRef silently keeps its previous/default ms_run_, so a caller passing 0 gets a wrong-but-silent reference. The 1-based requirement is also undocumented in the header.
- **Evidence:** MzTabBase.cpp:181-188 `void setMSFile(Size index){ assert(index >= 1); if (index >= 1){ ms_run_ = index; } }` - no else branch, no log, no throw (asserts are compiled out in release).
- **Fix:** Document the 1-based contract in the header and, on index==0, throw Exception::InvalidValue or at least OPENMS_LOG_WARN instead of silently no-op'ing. Behavior fix; no ABI change.
- **Verifier correction:** setMSFile(0) is a silent no-op in release builds (only a debug assert), and the 1-based mzTab contract is undocumented in the header. But because ms_run_ defaults to 0 and isNull() treats ms_run_<1 as the null/unset state, passing 0 to a fresh object yields a correct "null" cell, not a wrong reference. The only genuine silent failure is overwriting an already-set valid index with 0 (a stale value then persists), which is uncommon and recoverable; no in-tree caller does this. Severity is low (mild surprise / invites misuse), not high. Documenting the 1-based contract and warning/throwing on 0 is a worthwhile hardening with no ABI change.
- **Verified:** Evidence is accurate: MzTabBase.cpp:181-188 shows `setMSFile(Size index){ assert(index>=1); if(index>=1){ms_run_=index;} }` with no else/log/throw, and the header (MzTabBase.h:300) does not document the 1-based contract. The no-op on index 0 in release builds is real (assert is compiled out). However, the severity framing ("a caller passing 0 gets a wrong-but-silent reference") is overstated. ms_r

### [FORM-73] MzTabPSMSectionRow::addPepEvidenceToRows — addPepEvidenceToRows operates on the single *this row and overwrites fields (name says add-to-Rows)
`low` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/MzTab.h · _format-mztab_

```cpp
void addPepEvidenceToRows(const std::vector<PeptideEvidence>& peptide_evidences)
```
- **Expectation:** A member named addPepEvidenceToRows (plural Rows, verb 'add') sounds like it appends evidence into a collection of rows, or accumulates onto existing data.
- **Actual:** It is a member on a single MzTabPSMSectionRow and unconditionally (re)writes this row's pre/post/start/end/accession fields - including resetting them to empty MzTabString() when the evidence vector is empty. It neither touches a collection of rows nor 'adds' to existing values; it overwrites them.
- **Evidence:** MzTab.cpp:510-519 sets `pre = MzTabString(); post = MzTabString(); start = ...; end = ...;` on empty input and otherwise rebuilds those members from the evidence; the method is declared on the row struct (MzTab.h:287).
- **Fix:** Rename to something like `setPepEvidenceFields` / `fillFromPepEvidence` (keep the old name as a deprecated inline forwarder for ABI). Doc-only fix is the minimal step; rename is source-compatible if a forwarder is kept.
- **Verifier correction:** The name misleads on two counts: "add" (it overwrites/replaces this row's pre/post/start/end/accession rather than appending) and "Rows" (it mutates a single *this MzTabPSMSectionRow, not a collection). On empty input it clears pre/post/start/end (but NOT accession) and returns; on non-empty input it rebuilds those five fields. Minimal fix is doc-only (clarify it sets/replaces this row's fields, abi_impact none); a clean rename such as setPepEvidenceFields/fillFromPepEvidence with a deprecated inline forwarder is source-compatible.
- **Verified:** Independently verified. addPepEvidenceToRows is a member of struct MzTabPSMSectionRow (declared MzTab.h:287, defined MzTab.cpp:510) and operates solely on *this, never on a collection of rows — so the plural "Rows" in the name is unjustified (the single call site at MzTab.cpp:1286 passes one freshly-built row, with a comment that also wrongly says "to Rows"). The verb "add" is also misleading: on 

### [FORM-88] OMSFile::load — load()'s filename is documented as @param[out]
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/OMSFile.h · _format-oms-sqlite_

```cpp
void load(const std::string& filename, IdentificationData& id_data);
```
- **Expectation:** filename is the input path the data is read FROM; it should be documented as @param[in]. The data object is the out-parameter.
- **Actual:** Every load() overload tags the input path as @param[out] filename and the populated data object as @param[in], inverting the actual data-flow direction. store() correctly documents filename as @param[in].
- **Evidence:** OMSFile.h:57-58 `@param[out] filename The input file` / `@param[in] id_data The IdentificationData object`; contrast OMSFile.h:36-37 store() `@param[in] filename The output file`.
- **Fix:** Swap the direction tags: filename -> @param[in], data object -> @param[out]. Doc-only, ABI-safe.
- **Verifier correction:** All three load() overloads (OMSFile.h:57, 65, 71) document the input path as @param[out] filename and the populated data object as @param[in], inverting the data-flow tags. Fix: filename -> @param[in], data object (id_data/features/consensus) -> @param[out]. Doc-only, ABI-safe. Real but low severity: filename is a const std::string& and the comment text itself says "The input file," so the tag is self-evidently wrong on inspection and cannot cause silently wrong behavior.
- **Verified:** Verified against OMSFile.h:57-58 and the OMSFile.cpp implementation. All three load() overloads tag the input path as @param[out] filename and the populated data object as @param[in], inverting the real data flow. The signature `void load(const std::string& filename, IdentificationData& id_data)` proves it: filename is const (cannot be an out-param) and is fed to OMSFileLoad helper(filename,...) a

### [FORM-90] OSWFile::getRunID — getRunID throws on zero runs too, but doc only mentions 'more than one run'
`low` · `surprising-throw` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/OSWFile.h · _format-oms-sqlite_

```cpp
UInt64 getRunID() const; // @throws Exception::SqlOperationFailed more than on run exists
```
- **Expectation:** From the header note '@throws Exception::SqlOperationFailed more than on run exists', a caller expects an exception only when there are 2+ runs, and presumably a normal return otherwise.
- **Actual:** The implementation throws whenever nr_results != 1, i.e. also when the RUN table is empty (zero runs). Additionally, if zero rows are returned, the returned UInt64 'id' would be uninitialized were the throw ever bypassed; the value is only assigned inside the row loop.
- **Evidence:** OSWFile.cpp:916-919 `if (nr_results != 1) { throw Exception::SqlOperationFailed(... "contains more than one run..."); }`; `UInt64 id;` declared uninitialized at OSWFile.cpp:907 and only set at :911 inside the while loop.
- **Fix:** Update the doc/throw message to cover the 'no run' case (e.g. 'requires exactly one run'), and initialize id to avoid relying on the throw. Doc + impl change, ABI-safe.
- **Verifier correction:** getRunID() throws Exception::SqlOperationFailed whenever the RUN table does not contain exactly one row — i.e. on zero runs as well as on 2+ runs — but both the header doc and the exception message only mention the "more than one run" case. The failure is loud and recoverable (an exception, not silent wrong data). The uninitialized `UInt64 id;` is never actually returned because `return id;` is reachable only when exactly one row was read (and thus `id` was assigned), so there is no real uninitialized-read defect; the recommendation to initialize `id` is defensive/cosmetic only. Suggested fix: correct the doc and message to "requires exactly one run" and optionally zero-initialize `id`.
- **Verified:** Code confirmed at OSWFile.cpp:897-921. The header doc (OSWFile.h:156) says "@throws Exception::SqlOperationFailed more than on run exists", but the impl throws on `nr_results != 1` (line 916), which also fires when the RUN table is EMPTY (zero runs), and the throw message (line 918) only mentions "more than one run". So the doc and the exception message both fail to cover the zero-run case — a gen

### [FORM-91] SqliteConnector::tableExists / columnExists / countTableRows / executeStatement / prepareStatement / getDB — Read-only query helpers are non-const, so they can't be called on a const SqliteConnector
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/SqliteConnector.h · _format-oms-sqlite_

```cpp
bool tableExists(const std::string& tablename); bool columnExists(const std::string&, const std::string&); Size countTableRows(const std::string&); sqlite3* getDB();
```
- **Expectation:** tableExists/columnExists/countTableRows are logical reads and should be usable on a const SqliteConnector& (e.g. when a class holds the connection and exposes const accessors). getDB() returning a mutable handle is expected, but the pure-query helpers being non-const is surprising.
- **Actual:** None of these member overloads are const-qualified, even tableExists/columnExists/countTableRows which only SELECT. They forward to the static overloads taking the raw db_ pointer, so const-ness is purely a missing qualifier, not a technical limitation.
- **Evidence:** SqliteConnector.h:76 `bool tableExists(const std::string& tablename) { return tableExists(db_, tablename); }`, :83 `Size countTableRows(...)`, :93 `bool columnExists(...)`; all non-const; underlying statics SqliteConnector.h:157,168 are free of any mutation contract.
- **Fix:** Add const overloads (or make tableExists/columnExists/countTableRows const) for the read-only query helpers. Adding const member functions is source-compatible but changes the mangled symbols / vtable-independent ABI of those inline methods; since they are inline header methods the practical ABI risk is low.
- **Verifier correction:** Severity is low, not higher: the limitation is enforced at compile time (fails loudly, no silent wrong results / data loss / crash) and is trivially worked around by holding the connection non-const. The recommendation's ABI note is slightly off: tableExists and columnExists are inline header methods (no exported symbol, so adding const is effectively ABI-neutral), but countTableRows is an out-of-line exported symbol in SqliteConnector.cpp whose mangled name WOULD change if made const — that single symbol is an ABI break even though all source callers recompile unchanged. Cleanest fix is to make the three read-only helpers const (member db_ stays usable as sqlite3* const), mirroring SQLiteCpp::Database::tableExists() const already shipping in the repo.
- **Verified:** Evidence verified in source. SqliteConnector.h:76 tableExists, :93 columnExists are non-const inline members forwarding to static sqlite3*-taking overloads; :83 countTableRows is a non-const member defined out-of-line in SqliteConnector.cpp:92. The .cpp confirms all three are pure reads: tableExists runs "SELECT 1 FROM sqlite_master..." (line 82), columnExists runs "PRAGMA table_info(...)" (line 6

### [FORM-92] OSWFile::readProtein — @param[out] index on a by-value Size index parameter
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/OSWFile.h · _format-oms-sqlite_

```cpp
void readProtein(OSWData& swath_result, const Size index);
```
- **Expectation:** index is a by-value input selecting which protein to populate; it is read, not written. swath_result is the in/out data being populated.
- **Actual:** The header tags index as @param[out] although it is a const-by-value Size used only to index into swath_result.getProteins(); the populated container swath_result is tagged @param[in]. The direction tags are inverted, mirroring the OMSFile::load defect.
- **Evidence:** OSWFile.h:74-75 `@param[in] swath_result OSWData ...` / `@param[out] index Index into swath_result.getProteins()[index]`; signature OSWFile.h:78 `void readProtein(OSWData& swath_result, const Size index);`.
- **Fix:** Fix the direction tags: swath_result -> @param[in,out], index -> @param[in]. Doc-only, ABI-safe.
- **Verifier correction:** The header direction tags are inverted but harmlessly. `index` is a `const Size` (size_t) passed by value (OSWFile.h:78), used only to index into swath_result.getProteins()[index] (OSWFile.cpp:575-582); it is a read-only input yet tagged @param[out]. `swath_result` is the OSWData& populated with peptides (getFullProteins_, OSWFile.cpp:579) but tagged @param[in]. Correct tags: `@param[in] index`, `@param[in,out] swath_result`. Doc-only, ABI-safe. Severity is low, not high: the const-by-value signature plus surrounding prose make the out-tag obviously a typo, so it cannot mislead into incorrect or silently-wrong usage.
- **Verified:** The quoted evidence is accurate. OSWFile.h:74-75 tags `@param[in] swath_result` and `@param[out] index`, while the signature (line 78) is `void readProtein(OSWData& swath_result, const Size index)`. OSWFile.cpp:573-584 confirms the real data flow: `index` is `const Size` (Size = typedef size_t, a by-value unsigned scalar, provably read-only here, used only as swath_result.getProteins()[index]), an

### [FORM-100] ParamJSONFile::load — load() is static while store() is a non-static member, and the flatHierarchy member only affects store
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/FORMAT/ParamJSONFile.h · _format-param_

```cpp
static bool load(const std::string& filename, Param& param);
```
- **Expectation:** Sibling load/store pairs in the *File classes are both non-static instance methods (cf. ParamXMLFile, ToolDescriptionFile), so a caller expects `ParamJSONFile f; f.load(...); f.store(...);` to work uniformly and for the object's `flatHierarchy` flag to influence both directions.
- **Actual:** load() is `static` (called as `ParamJSONFile::load(...)`), store()/writeToStream() are non-static and read the `flatHierarchy` member. So `flatHierarchy` silently has no effect on load(); a caller who sets `f.flatHierarchy=true` then calls `f.load(...)` gets the flag ignored, and the call style differs between the two methods of the same class.
- **Evidence:** Header: `static bool load(...)` vs `void store(...) const;` and member `bool flatHierarchy{};`. Impl: `bool ParamJSONFile::load(...)` does not reference flatHierarchy; only writeToStream/store use `if (!flatHierarchy)` (line 295).
- **Fix:** Make load() a non-static member for symmetry with store() and the other *File classes, or document that load is static and unaffected by flatHierarchy. Converting static→non-static changes the mangled name (ABI break for the static symbol) — if ABI must hold, keep static but clearly document the asymmetry; the clean fix is a non-static load.
- **Verifier correction:** The asymmetry is real and the evidence is correct, but the practical impact is minor. ParamJSONFile::load is static (called as ParamJSONFile::load(...)) while store()/writeToStream() are non-static const and read the flatHierarchy member, unlike the sibling ParamXMLFile which has both load and store as non-static instance methods. However, flatHierarchy is documented (header lines 60-63) as a write-only setting and load() is intentionally format-agnostic (line 74 maps `__`->`:`, accepting both flat and nested key styles), so the flag has no defined meaning for load and its absence there is documented/by-design rather than a silent bug. The actual surprise is limited to the inconsistent call style (static vs instance method) versus the other *File classes — a loud, cosmetic convention mismatch, not a source of wrong results, data loss, or crashes. Recommendation: for symmetry, either make load() a non-static member (note: static->non-static changes the mangled name = ABI break for that symbol) or simply document the static asymmetry; either way severity is low. abi_impact of the current code is none; only the proposed non-static fix would be ABI-breaking.
- **Verified:** The code facts are all independently confirmed. Header (ParamJSONFile.h): `static bool load(const std::string&, Param&)` (line 76) vs non-static `void store(...) const` and `void writeToStream(...) const` (lines 85-86), with public member `bool flatHierarchy{}` (line 64). Impl (ParamJSONFile.cpp): load() (lines 38-168) never references flatHierarchy; only writeToStream() reads it via `if (!flatHie

### [FORM-101] ParamJSONFile::store — Doxygen tags store()'s param as @param[in,out] but it is const Param& (read-only)
`low` · `documentation-defect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/ParamJSONFile.h · _format-param_

```cpp
void store(const std::string& filename, const Param& param, const ToolInfo& tool_info) const;
```
- **Expectation:** An `@param[in,out]` annotation tells the caller their Param may be modified by store(). A developer might pass a non-const working copy expecting store to update it.
- **Actual:** param is declared `const Param&`; store() cannot and does not modify it. The `[in,out]` tag is simply wrong.
- **Evidence:** Header: '@param[in,out] param The param data structure that should be stored.' against signature `const Param& param`.
- **Fix:** Change the Doxygen direction to `@param[in] param`. Doc-only fix, ABI-neutral.
- **Verifier correction:** The header's Doxygen tag `@param[in,out] param` contradicts the `const Param& param` signature and should be `@param[in]`. This is a trivial, ABI-neutral documentation typo, not a behavioral/naming surprise: the const-qualified reference makes mutation impossible and the signature (right next to the comment) is authoritative, so no competent C++ dev would be misled into expecting store() to update their Param. Fix: change [in,out] to [in] (and ideally fix load()'s filename tag from [out] to [in] for consistency). Recommend reclassifying from misleading-name/high to documentation-defect/low.
- **Verified:** The quoted evidence is factually accurate: ParamJSONFile.h line 81 tags param as `@param[in,out]` while line 85 declares it `const Param& param`, and the .cpp (store→writeToStream) only ever reads it, so the [in,out] direction is wrong and should be [in]. However, this is NOT a genuine POLS surprise to a competent C++ dev. The signature `const Param&` is on the very next line and is compiler-enfor

### [FORM-102] ParamXMLFile::store / writeXMLToStream — Doxygen marks store()/writeXMLToStream() param (and the output stream) as @param[out] when they are read-only inputs
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/FORMAT/ParamXMLFile.h · _format-param_

```cpp
void store(const std::string& filename, const Param& param) const; void writeXMLToStream(std::ostream* os_ptr, const Param& param) const;
```
- **Expectation:** `@param[out] param` signals to the caller that store() fills/overwrites their Param. A reader of the header would think store() is a load-like operation populating param.
- **Actual:** param is `const Param&` (input that gets serialized, never written). The `[out]` tag inverts the actual data direction. (writeXMLToStream additionally tags BOTH os_ptr and param as [out], but only os_ptr is genuinely output.)
- **Evidence:** Header: store(): '@param[out] param The Param class that should be stored in the file.' with `const Param& param`. writeXMLToStream(): '@param[out] os_ptr ...' and '@param[out] param ...' with `const Param& param`.
- **Fix:** Tag param as `@param[in]` in both store() and writeXMLToStream(); keep os_ptr as `[out]`. Same wrong-direction tags appear in ParamCTDFile and ParamCWLFile store/write methods (`@param[out] param`/`@param[out] tool_info` against const&) and should be fixed identically. Doc-only, ABI-neutral.
- **Verified:** Evidence verified verbatim: ParamXMLFile.h line 34 tags `@param[out] param` and line 44 tags both `@param[out] os_ptr` AND `@param[out] param` for `void store(const std::string&, const Param&)` and `void writeXMLToStream(std::ostream*, const Param&)`. The .cpp (ParamXMLFile.cpp:30,55) confirms param is purely iterated/serialized and never mutated; only os_ptr is genuine output. So the [out] tags o

### [FORM-104] ParamJSONFile::flatHierarchy — flatHierarchy member doc is garbled and applies only to writing, not loading
`low` · `documentation-typo` · ABI: `none` · src/openms/include/OpenMS/FORMAT/ParamJSONFile.h · _format-param_

```cpp
bool flatHierarchy{};
```
- **Expectation:** A public bool named flatHierarchy with a doc comment should clearly say what toggling it does and in which direction (load vs store).
- **Actual:** The doc reads 'If set to true, all parameters will be listed on when writing the JSON file.' — a broken sentence ('listed on when'), and it silently applies only to store/writeToStream; load() ignores it entirely. ParamCWLFile's identically-named member has different (correct) wording 'listed without nesting', so the sibling docs disagree.
- **Evidence:** ParamJSONFile.h: '\brief If set to true, all parameters will be listed on when writing the JSON file.' vs ParamCWLFile.h: '\brief If set to true, all parameters will be listed without nesting when writing the CWL File.' load() impl never reads flatHierarchy.
- **Fix:** Fix the wording to match ParamCWLFile ('listed without nesting (flattened, names include the hierarchy)') and state explicitly it only affects store()/writeToStream(), not load(). Doc-only, ABI-neutral.
- **Verifier correction:** The doc comment for ParamJSONFile::flatHierarchy (ParamJSONFile.h:61) contains a grammatical typo ("listed on when writing the JSON file") and uses wording inconsistent with the identically-named ParamCWLFile::flatHierarchy ("listed without nesting"). This is a cosmetic documentation defect, not a misleading name: the member name is accurate and the second sentence of the same comment ("The names will be expanded to include the nesting hierarchy.") correctly describes the effect, and both comments already state it applies "when writing", so the write-only direction is documented. Recommended fix is to align the JSON wording with the CWL wording and optionally note explicitly that it affects store()/writeToStream() only. Doc-only, ABI-neutral.
- **Verified:** The factual core checks out: ParamJSONFile.h:61 does read "all parameters will be listed on when writing the JSON file." which is grammatically broken ("listed on when"), and the sibling ParamCWLFile.h:24 uses different wording ("listed without nesting"). Confirmed in the .cpp that flatHierarchy is read only in writeToStream() (ParamJSONFile.cpp:295) and load() (lines 38-168) never references it. 

### [FORM-31] EDTAFile::store(const std::string&, const ConsensusMap&) const — store(ConsensusMap) documented "NOT IMPLEMENTED" but is fully implemented
`low` · `stale-documentation` · ABI: `none` · src/openms/include/OpenMS/FORMAT/EDTAFile.h · _format-peakfiles_

```cpp
void store(const std::string & filename, const ConsensusMap & map) const
```
- **Expectation:** The Doxygen comment says "NOT IMPLEMENTED", so a caller reading the header expects the call to throw NotImplemented (as the sibling SpecArrayFile/KroenikFile stores do) or be a no-op.
- **Actual:** EDTAFile.cpp:272 contains a complete, working implementation that writes a consensus EDTA file (header with RT/m/z/intensity/charge quadruplets, sub-feature columns, NA padding) and calls tf.store(filename).
- **Evidence:** Header line 90-95: "@brief Stores a ConsensusMap as an enhanced DTA file.\n\n NOT IMPLEMENTED" above `void store(const std::string & filename, const ConsensusMap & map) const;`. Implementation EDTAFile.cpp:272 `void EDTAFile::store(...) { ... tf.store(filename); }` is real.
- **Fix:** Source-compatible doc-only fix: remove the "NOT IMPLEMENTED" line from the Doxygen block and document the actual output format. No ABI change. (The inverse risk — callers who coded around a believed no-op — is why the stale doc is dangerous.)
- **Verifier correction:** store(const std::string&, const ConsensusMap&) const is fully implemented and round-trip-tested (EDTAFile.cpp:272-317; EDTAFile_test.cpp:74-99). The Doxygen "NOT IMPLEMENTED" line in EDTAFile.h:91 is stale/contradictory documentation, not a silent failure: the call succeeds and writes a valid consensus EDTA file. The danger is doc-driven (a caller trusting the doc may wrongly assume the call throws/is a no-op, mirroring the genuine SpecArrayFile/KroenikFile convention where the same phrase pairs with a real throw). Fix is comment-only: delete the "NOT IMPLEMENTED" line and document the RT/m/z/intensity/charge + sub-feature/NA-padded output format. Re-categorize as stale-documentation, severity low (no wrong results/crash possible), ABI impact none.
- **Verified:** Verified against source. EDTAFile.h:88-95 documents store(const std::string&, const ConsensusMap&) with the Doxygen line "NOT IMPLEMENTED", but EDTAFile.cpp:272-317 is a complete, working writer: it validates the .edta extension, computes max sub-features, emits a header of RT/m/z/intensity/charge quadruplets with per-sub-feature suffixes, writes consensus rows plus sub-feature columns with NA pad

### [FORM-34] SpecArrayFile::store(const std::string&, const SpectrumType&) const — store() takes a parameter named 'spectrum' although load() reads/writes a FeatureMap; it also unconditionally throws
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/SpecArrayFile.h · _format-peakfiles_

```cpp
template <typename SpectrumType> void store(const std::string& filename, const SpectrumType& spectrum) const
```
- **Expectation:** store() should be the inverse of load(), which fills a FeatureMap; the symmetric store signature should accept a FeatureMap (parameter conventionally named feature_map), matching the documented "Stores a featureXML as a SpecArray file".
- **Actual:** store() is templated on SpectrumType, names its data parameter `spectrum`, prints to std::cerr, and then always throws NotImplemented regardless of input. The name/signature suggest a per-spectrum store, contradicting the FeatureMap-oriented load() and the "featureXML" doc.
- **Evidence:** SpecArrayFile.h:101-106 `template <typename SpectrumType> void store(...const SpectrumType& spectrum...) { std::cerr << "Store() for SpecArrayFile not implemented..."; throw Exception::NotImplemented(...); }`, vs load() taking `FeatureMapType& feature_map` (line 50).
- **Fix:** Source-compatible: keep the throw but rename/retype the unused parameter to `const FeatureMap&` (or document it) so the not-implemented stub matches the format's data model; or drop the cerr print since the thrown exception already signals it. No external ABI surface (templated, always throws).
- **Verifier correction:** The inconsistency is real and confirmed, but the "unconditionally throws" aspect is an intentional documented NOT-IMPLEMENTED stub, not a surprise. The genuine, surprising part is the convention mismatch: store() is templated on `SpectrumType` with parameter `spectrum` (and dereferences spectrum.size()), while load() uses `FeatureMapType feature_map`, the class/store docs describe a featureXML/feature format, and the sole real caller (FileHandler::storeFeatures, line 1349) passes a FeatureMap. Severity is low (loud, always-throwing, never silent), not high/medium. Recommended fix is source-compatible: rename/retype the parameter to a FeatureMap-oriented name to match load() and the doc.
- **Verified:** Evidence verified verbatim at SpecArrayFile.h:101-106. store() is templated on SpectrumType, names its parameter `spectrum`, prints to std::cerr (including spectrum.size()), and unconditionally throws NotImplemented. load() (line 49-50) is templated on FeatureMapType / feature_map and builds Feature objects; class doc says the format is feature/.pepList-oriented and store()'s own doc says "Stores 

### [FORM-35] MSPFile::load(const std::string&, PeptideIdentificationList&, PeakMap&) — Doxygen marks the output PeakMap exp as @param[in] and lists params in an order that does not match the signature
`low` · `documentation` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MSPFile.h · _format-peakfiles_

```cpp
void load(const std::string & filename, PeptideIdentificationList & ids, PeakMap & exp)
```
- **Expectation:** The Doxygen @param entries should match the signature order (filename, ids, exp) and mark the filled PeakMap as output.
- **Actual:** The doc block lists `@param[in] exp PeakMap which contains the spectra after reading` (exp is an output, and load() does `exp.reset()` then fills it) before `@param[in] filename` and `@param[out] ids`. So exp is mis-marked [in] and the listing order (exp, filename, ids) is scrambled relative to the (filename, ids, exp) signature.
- **Evidence:** MSPFile.h:54-62 `@param[in] exp PeakMap which contains the spectra after reading` / `@param[in] filename ...` / `@param[out] ids ...` above `void load(const std::string & filename, PeptideIdentificationList & ids, PeakMap & exp);`.
- **Fix:** Doc-only fix: reorder to match signature and mark `@param[out] exp`. No ABI impact.
- **Verifier correction:** The Doxygen block for MSPFile::load(const std::string&, PeptideIdentificationList&, PeakMap&) (MSPFile.h:54-56) lists @param entries in order (exp, filename, ids) which does not match the signature order (filename, ids, exp), and marks the filled output PeakMap exp as @param[in] although load() calls exp.reset() and populates it (MSPFile.cpp:72ff). Correct fix is doc-only: reorder to filename, ids, exp and change exp to @param[out]. Impact is a low-severity documentation inconsistency; the descriptive text already states exp is filled, so semantics are not silently misrepresented. No ABI or source-compatibility impact.
- **Verified:** Evidence is accurate. MSPFile.h:54-62 documents the three-arg load() with @param entries in the order (exp, filename, ids) while the signature order is (filename, ids, exp), and it marks `exp` as @param[in] even though MSPFile.cpp:54-72 shows load() does exp.reset() and then fills exp — so exp is genuinely an output. By contrast, ids is correctly tagged @param[out] in the same block, and OpenMS us

### [FORM-83] SwathFile::loadSqMass — loadSqMass omits the tmp/readoptions parameters every other SwathFile loader has
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/FORMAT/SwathFile.h · _format-targeted_

```cpp
std::vector<OpenSwath::SwathMap> loadSqMass(const std::string& file, std::shared_ptr<ExperimentalSettings>& exp_meta)
```
- **Expectation:** All SwathFile::load* members share the shape (file, tmp, exp_meta, readoptions): loadMzML, loadMzXML, loadSplit, loadFromMSExperiment, loadBrukerTdf. A caller iterating over loaders or templating over them expects loadSqMass to follow the same convention (e.g. honoring readoptions == "cache").
- **Actual:** loadSqMass has only (file, exp_meta). It silently ignores caching/readoptions entirely and always builds SpectrumAccessSqMass directly. A caller who selected readoptions="cache" for the run gets a different access strategy for sqMass inputs with no signal.
- **Evidence:** Header line 100: `std::vector<OpenSwath::SwathMap> loadSqMass(const std::string& file, std::shared_ptr<ExperimentalSettings>& /* exp_meta */);` vs line 69-73 loadMzML which has `tmp`, `readoptions`, `plugin_consumer`.
- **Fix:** Document explicitly in the header that loadSqMass ignores readoptions/tmp by design (it streams from the SQLite DB), so the asymmetry is intentional and visible. A future-additive overload accepting (file, tmp, exp_meta, readoptions) for signature parity would be source- and ABI-compatible.
- **Verifier correction:** loadSqMass really does omit (tmp, readoptions) that all other SwathFile loaders take, and additionally ignores its own exp_meta out-param (commented out / unused). This is a real inconsistent-convention surprise. But the "silently different access strategy / no signal" harm is largely benign: readoptions=="cache" exists to convert streamed/in-memory data into an on-disk lazy-access format, which sqMass already is via SpectrumAccessSqMass. Hence honoring "cache" would be a near no-op for sqMass; ignoring it does not produce wrong results, data loss, or a crash. Severity downgraded high->low. The recommended fix (document the intentional asymmetry in the header and/or add an additive (file, tmp, exp_meta, readoptions) overload) is sound and source-compatible; current ABI impact is none.
- **Verified:** The factual core is verified. Header line 100 and SwathFile.cpp line 310 confirm loadSqMass(const std::string& file, std::shared_ptr<ExperimentalSettings>& exp_meta) — it is the only SwathFile loader lacking the (tmp, readoptions) parameters that loadMzML/loadMzXML/loadSplit/loadFromMSExperiment/loadBrukerTdf all carry. The implementation always builds SpectrumAccessSqMass directly and never consu

### [FORM-85] SqMassFile::transform — transform exposes two dead bool parameters that have no effect
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/SqMassFile.h · _format-targeted_

```cpp
void transform(const std::string& filename_in, Interfaces::IMSDataConsumer* consumer, bool skip_full_count = false, bool skip_first_pass = false) const
```
- **Expectation:** A boolean parameter named skip_first_pass=true actually skips the first pass; skip_full_count=true skips the full count. Callers tune performance with them.
- **Actual:** Both parameters are entirely ignored — the implementation comments them out (`bool /* skip_full_count */, bool /* skip_first_pass */`) and the first-pass call is dead-commented. Setting either to true changes nothing. This is honestly disclosed in the header @note, but the bare-bool defaults at call sites remain misleading.
- **Evidence:** SqMassFile.cpp line 37: `void SqMassFile::transform(const std::string& filename_in, Interfaces::IMSDataConsumer* consumer, bool /* skip_full_count */, bool /* skip_first_pass */) const`; header @note lines 126-128: "@p skip_full_count and @p skip_first_pass parameters are currently unused and have no effect".
- **Fix:** Keep them for ABI but consider [[maybe_unused]] and a clearer name, or deprecate the overload and add a 2-arg transform. Removing the params would be ABI-breaking; the documented @note already mitigates, so this is low severity.
- **Verifier correction:** The pyOpenMS evidence is slightly off: bind_misc.cpp:4669 and :4710 bind MzMLFile::transform and MzXMLFile::transform (where the flags are functional), not SqMassFile::transform. The dead-flag finding is correct for SqMassFile specifically. Real SqMassFile::transform callers are FileConverter.cpp:1032 (passes false/false) and OpenSwathMzMLFileCacher.cpp:167/175 (passes true/true, expecting a speedup that never happens). The flags affect only performance, never output correctness — output data is identical regardless of flag value — hence low severity. ABI impact of the finding itself is none; only the recommended fix (removing params) would be breaking, so keeping them and applying [[maybe_unused]]/deprecation as suggested is appropriate.
- **Verified:** Evidence verified verbatim in the actual code. SqMassFile.cpp:37 declares the params commented out: `void SqMassFile::transform(const std::string& filename_in, Interfaces::IMSDataConsumer* consumer, bool /* skip_full_count */, bool /* skip_first_pass */) const`. Line 43 dead-comments the first-pass call: `// if (!skip_first_pass) transformFirstPass_(filename_in, consumer, skip_full_count);`. The t

### [FORM-15] GzipIfstream::read (doc) and class @brief — GzipIfstream docs reference bzip2 and the wrong file extension (*.gzip)
`low` · `incorrect-documentation` · ABI: `none` · src/openms/include/OpenMS/FORMAT/GzipIfstream.h · _format-text-streams_

```cpp
size_t read(char * s, size_t n)
```
- **Expectation:** Documentation for a gzip decompressor should describe gzip and the conventional '.gz' extension.
- **Actual:** The class @brief says 'Decompresses files which are compressed in the gzip format (*.gzip)' (gzip files use '.gz'), and read()'s doc is a copy-paste from Bzip2Ifstream: 'Reads n bytes from the bzip2 compressed file into buffer s'. A maintainer/caller could be misled about which format/handle is in play.
- **Evidence:** GzipIfstream.h line 18 '(*.gzip)' and line 33 '@brief Reads n bytes from the bzip2 compressed file into buffer s'.
- **Fix:** Fix the doc to say gzip and '.gz'. Pure documentation change, no ABI impact.
- **Verifier correction:** GzipIfstream.h contains two documentation errors (not a misleading name): the class @brief (line 18) says the gzip extension is '*.gzip' when gzip files conventionally use '.gz', and the read() @brief (line 33) is a copy-paste from Bzip2Ifstream reading 'Reads n bytes from the bzip2 compressed file into buffer s' despite this being the zlib/gzip class (uses <zlib.h> and gzFile). Pure documentation fix; no ABI or behavioral impact. Severity is low because the class name, includes, and member types make the correct format obvious, and no extension checking exists.
- **Verified:** Both quoted strings exist verbatim. Line 18: '@brief Decompresses files which are compressed in the gzip format (*.gzip)' — gzip uses '.gz', not '.gzip'. Line 33: '@brief Reads n bytes from the bzip2 compressed file into buffer s' — a verbatim copy-paste from the bzip2 class, even though GzipIfstream includes <zlib.h>, owns a 'gzFile gzfile_', and is unambiguously the gzip/zlib decompressor. So th

### [FORM-16] TextFile::store — store() is non-const although it only reads the internal buffer
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FORMAT/TextFile.h · _format-text-streams_

```cpp
void store(const std::string& filename)
```
- **Expectation:** A store/serialize method that does not mutate the object is conventionally const, so it can be called on a const TextFile and signals read-only intent.
- **Actual:** store() is declared and defined non-const even though it only iterates buffer_ to write the file and never mutates the object. This prevents storing a const TextFile and is inconsistent with begin()/end() which provide const overloads.
- **Evidence:** TextFile.h line 80 `void store(const std::string& filename);`; TextFile.cpp lines 79-109 only read `buffer_`. Same pattern in CsvFile::store.
- **Fix:** Mark store() const (and CsvFile::store). This is source-compatible for callers but changes the mangled signature, so it is technically an ABI break on the symbol; if strict ABI must be preserved, add a const overload or defer to the next ABI-breaking release.
- **Verifier correction:** store() is non-const and only reads buffer_, which is verified. But this is the dominant OpenMS convention for serializer classes that hold their own buffer (5 of 6 such single-arg store(filename) methods are non-const), not a genuine surprise. The only effect of the missing const is a loud, immediate compile error if one tries to store a const TextFile — never silent wrong behavior. Marking store() const is a reasonable cleanup but is a low-severity style improvement, and it is an ABI break on the mangled symbol (store -> name mangles differently with the const qualifier), so it should be deferred to an ABI-breaking release or done via a const overload.
- **Verified:** Code facts are accurate: TextFile.h:80 declares `void store(const std::string& filename)` non-const, and TextFile.cpp:79-109 only iterates buffer_ read-only (`for (const std::string& it : buffer_)`) and never mutates the object; begin()/end() do have const overloads; CsvFile::store (CsvFile.cpp:37-40) is likewise non-const and just forwards. So the claim is not wrong about the code. However, the f

### [FORM-148] SemanticValidator::setNameAttribute / setValueAttribute / setUnitAccessionAttribute — Doc comments for the attribute-name setters are copy-pasted from setAccessionAttribute and describe the wrong attribute
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/FORMAT/VALIDATORS/SemanticValidator.h · _format-validators-options_

```cpp
void setNameAttribute(const std::string& name); void setValueAttribute(const std::string& value)
```
- **Expectation:** setNameAttribute should be documented as setting the 'name' attribute and setValueAttribute as setting the 'value' attribute, consistent with what each setter actually writes.
- **Actual:** Both carry the copy-pasted comment 'Sets the name of the attribute for accessions in the CV parameter tag name'. setNameAttribute actually sets name_att_ and setValueAttribute sets value_att_ (SemanticValidator.cpp:57-65), not the accession attribute. A reader trusting the doc would think three different setters all configure the accession attribute. setUnitAccessionAttribute's '(default: unitAccession)' note is fine, but the name/value setters' bodies contradict their doc.
- **Evidence:** Header: line 80 '/// Sets the name of the attribute for accessions ... (default: 'name')' above setNameAttribute; line 83 '/// Sets the name of the attribute for accessions ... (default: 'value')' above setValueAttribute. Impl: `name_att_ = name;` and `value_att_ = value;` (SemanticValidator.cpp:59,64).
- **Fix:** Fix the doc comments to describe the name and value attributes respectively. Comment-only change; no ABI/source impact.
- **Verifier correction:** The two doc comments are not byte-for-byte identical as the claim states; their '(default: 'name')' / '(default: 'value')' suffixes are correct. The actual defect is the copy-pasted leading clause 'Sets the name of the attribute for accessions in the CV parameter tag name' appearing above setNameAttribute (sets name_att_) and setValueAttribute (sets value_att_), where it should describe the name and value attributes respectively. setUnitAccessionAttribute (line 105) and setUnitNameAttribute (line 108) are documented correctly.
- **Verified:** Verified against actual source. Header (SemanticValidator.h:80-84): setNameAttribute is documented '/// Sets the name of the attribute for accessions in the CV parameter tag name (default: 'name')' and setValueAttribute '/// Sets the name of the attribute for accessions ... (default: 'value')'. Both reuse the 'for accessions' phrasing copied from setAccessionAttribute (line 77). The implementation

### [FORM-149] PeakFileOptions::getNumpressConfigurationMassTime / getNumpressConfigurationIntensity / getNumpressConfigurationFloatDataArray — Set/Get doc labels are swapped on all three numpress configuration accessor pairs
`low` · `incorrect-or-swapped-doc-comment` · ABI: `none` · src/openms/include/OpenMS/FORMAT/OPTIONS/PeakFileOptions.h · _format-validators-options_

```cpp
NumpressConfig getNumpressConfigurationMassTime() const; void setNumpressConfigurationMassTime(NumpressConfig)
```
- **Expectation:** The doc above a `get...() const` accessor should say 'Get', and the doc above a `set...()` mutator should say 'Set'.
- **Actual:** Each of the three pairs has the labels reversed: the comment 'Set numpress configuration options...' sits above the getter and 'Get numpress configuration options...' sits above the setter. A reader scanning comments would call the wrong member. (The setter for MassTime also has a hidden side effect -- it prints a warning to std::cerr for SLOF/PIC compression -- which the header does not mention.)
- **Evidence:** Header lines 181-192: '/// Set numpress configuration options for m/z or rt dimension' precedes `getNumpressConfigurationMassTime() const`; '/// Get numpress configuration options for m/z or rt dimension' precedes `setNumpressConfigurationMassTime(...)`; same swap for Intensity (185-188) and FloatDataArray (189-192). Impl side effect: setNumpressConfigurationMassTime emits 'Warning, compression of m/z or time dimension with pic or slof algorithms can lead to data loss' (PeakFileOptions.cpp:252-255).
- **Fix:** Swap the Set/Get words so each label matches its method, and add a note to setNumpressConfigurationMassTime that it warns on std::cerr for lossy SLOF/PIC configs. Comment-only change; no ABI/source impact.
- **Verifier correction:** All three numpress accessor pairs in PeakFileOptions.h (lines 181-192) have their Doxygen Set/Get labels swapped: '/// Set ...' sits above each `get...() const` getter and '/// Get ...' above each `void set...()` setter. This is a real comment-only defect, but low severity: the method names, get/set prefixes, const-ness, return types, and parameters all correctly indicate each member's role on the same line, so actual mis-use is unlikely. Additionally, setNumpressConfigurationMassTime (cpp:250-256) emits an undocumented std::cerr warning for SLOF/PIC compression (loud and recoverable, appropriate behavior). Recommend swapping the Set/Get words and documenting the warning. No ABI/source impact.
- **Verified:** Both pieces of evidence are confirmed verbatim. Header lines 181-192: all three getter declarations (getNumpressConfigurationMassTime/Intensity/FloatDataArray, each `const` returning NumpressConfig) are preceded by a '/// Set ...' comment, and all three setters (`void set...(config)`) by a '/// Get ...' comment — the Set/Get words are swapped on every pair. cpp:250-256 confirms setNumpressConfigur

### [FORM-150] PeakFileOptions::setFillData / setAlwaysAppendData / setSkipXMLChecks — Boolean setters declare their parameter as 'only', unrelated to what the flag controls
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/OPTIONS/PeakFileOptions.h · _format-validators-options_

```cpp
void setFillData(bool only); void setAlwaysAppendData(bool only); void setSkipXMLChecks(bool only)
```
- **Expectation:** A setter's bool parameter name should describe the flag (e.g. setFillData(bool fill), setSkipXMLChecks(bool skip)) so call sites and IDE hints read sensibly.
- **Actual:** In the header these three setters declare the parameter as `bool only`, copied from setMetadataOnly(bool only). setFillData(bool only) controls whether to fill data, setAlwaysAppendData(bool only) controls always-append, setSkipXMLChecks(bool only) controls skipping checks -- none of them is an 'only' flag. The .cpp definitions even use sensible names (fill_data, always_append_data, skip), so the public header is the misleading surface.
- **Evidence:** Header: `void setFillData(bool only);` (line 139), `void setAlwaysAppendData(bool only);` (line 135), `void setSkipXMLChecks(bool only);` (line 143). Impl uses different names: `void PeakFileOptions::setFillData(bool fill_data)` (PeakFileOptions.cpp:210), `setAlwaysAppendData(bool always_append_data)` (line 170), `setSkipXMLChecks(bool skip)` (line 180).
- **Fix:** Rename the header parameters to match the implementation (fill/always_append/skip). Parameter-name changes are source- and binary-compatible.
- **Verifier correction:** Evidence confirmed verbatim. Severity adjusted from any elevated level down to low: the misleading `only` parameter name is purely a readability/IDE-hint issue. Each setter carries a correct descriptive doc-comment on the line directly above it (header lines 134, 138, 142), and the parameter is a single bool that maps directly to the documented flag, so there is no path to silently wrong results, data loss, crashes, or argument-order misuse. Renaming the header parameters to fill_data/always_append_data/skip (matching the .cpp) is both source- and binary-compatible (parameter names are not part of the ABI signature).
- **Verified:** Verified directly. Header PeakFileOptions.h declares all three setters with parameter `only`: setAlwaysAppendData(bool only) (line 135), setFillData(bool only) (line 139), setSkipXMLChecks(bool only) (line 143). The .cpp uses sensible names: setAlwaysAppendData(bool always_append_data) (line 170), setSkipXMLChecks(bool skip) (line 180), setFillData(bool fill_data) (line 210). The `only` name is cl

### [FORM-106] MascotRemoteQuery::run — run() returns void and swallows every failure into an out-of-band error flag
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MascotRemoteQuery.h · _format-vendor_

```cpp
void run()
```
- **Expectation:** A method documented as 'Execute the full Mascot workflow synchronously (login -> query -> export results)' would either return a success/failure status or throw on failure (network down, bad login, server error, empty reply).
- **Actual:** run() returns void. On any failure (curl init failure, login error, HTTP >=400, empty server reply, unparseable response) it merely assigns error_message_ and returns normally. A caller who does not separately remember to call hasError() afterwards will treat a completely failed run as success and then read empty getMascotXMLResponse().
- **Evidence:** MascotRemoteQuery.cpp:81-151 run() — e.g. :88-91 `if (!curl) { error_message_ = "Failed to initialize libcurl"; return; }` and :135-139 `if (hasError()) { curl_easy_cleanup(curl); return; }`. Header :72 `void run();` with no return value or @throws.
- **Fix:** ABI-safe minimum: document in the header that run() never throws and that callers MUST check hasError()/getErrorMessage() afterward. Better (source-compatible): add a `bool runChecked()` returning success, or have run() throw on hard failures. Keep the existing void run() for ABI but redirect it.
- **Verifier correction:** run() does return void and reports all failures via the hasError()/getErrorMessage() flag rather than a return code or exception, as claimed. But this is a conventional OpenMS error-flag idiom, the sole in-tree caller (MascotAdapterOnline.cpp:329-336) already checks hasError() correctly, and the error channel is documented on the adjacent getErrorMessage()/hasError() members. The real, narrowly-scoped defect is that run()'s OWN doc comment omits the "never throws; callers MUST check hasError()" contract and the API offers no enforcement (no [[nodiscard]], no status return), making it a low-severity documentation/clarity footgun for a careless new caller — not a high-severity silent-wrong-results bug. Recommended fix (abi-safe, abi_impact none): add the contract to run()'s header doc comment. Optional source-compatible hardening: add a bool runChecked() or mark related accessors appropriately.
- **Verified:** The mechanical claim is accurate and verified. MascotRemoteQuery.h:71-72 declares `void run()` documented as "Execute the full Mascot workflow synchronously (login -> query -> export results)" with no @throws and no status return. MascotRemoteQuery.cpp:81-151 confirms run() returns void on every failure path: curl init failure (:88-91 `error_message_ = "Failed to initialize libcurl"; return;`), lo

### [FORM-107] MascotRemoteQuery::run — run() silently clears previously fetched responses/identifier via updateMembers_()
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MascotRemoteQuery.h · _format-vendor_

```cpp
void run()
```
- **Expectation:** Calling run() executes the configured workflow; results retained from a prior successful run (mascot_xml_, search_identifier_) would persist until explicitly overwritten by a new successful run.
- **Actual:** run() calls updateMembers_() as its first statement, which clears mascot_xml_, mascot_decoy_xml_, error_message_ and search_identifier_ before doing anything. If a second run() fails early (e.g. login error), the previously obtained XML/search-identifier are already wiped, so getMascotXMLResponse()/getSearchIdentifier() return empty even though valid data existed a moment earlier.
- **Evidence:** MascotRemoteQuery.cpp:83 `updateMembers_();` (first line of run()); updateMembers_ :521-524 `mascot_xml_.clear(); mascot_decoy_xml_.clear(); error_message_.clear(); search_identifier_.clear();`
- **Fix:** Document that run() resets all prior results up front. If unintended, move the result-clearing out of updateMembers_() (which is meant for param sync) into run()'s success path only, or clear lazily.
- **Verifier correction:** The mechanism is real but the scenario is mis-stated. The actual (mild) surprise is not "a second run() wipes results" — no caller reuses the object (single-use, copy/assign =deleted, MascotAdapterOnline creates a fresh instance per batch). The defensible anomaly is that result-clearing is placed inside updateMembers_(), a DefaultParamHandler override conventionally meant only to sync config from param_ and auto-invoked on setParameters(). Consequently calling setParameters() alone — without run() — silently discards prior results. This is recoverable (re-run regenerates) and never triggered in practice, so it is a low-severity design nit, not high-severity data loss. Recommendation stands in spirit: move result-clearing out of updateMembers_() into run()'s start so param syncing and result lifecycle are decoupled.
- **Verified:** Evidence verified exactly: run() (MascotRemoteQuery.cpp:83) calls updateMembers_() as its first statement, and updateMembers_() (lines 521-524) clears mascot_xml_, mascot_decoy_xml_, error_message_, and search_identifier_. So the described mechanism is real. However the claim's framing ("a second run() fails early and silently wipes previously valid results") is overstated and contrived: the only 

### [FORM-110] BrukerTimsFile::transform — transform() (streaming-consumer API) actually loads the entire file into memory first
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h · _format-vendor_

```cpp
void transform(const std::string& path, Interfaces::IMSDataConsumer* consumer)
```
- **Expectation:** A transform(path, consumer) modeled on OpenMS streaming consumers implies constant-memory, frame-by-frame feeding to the consumer — the whole point of a consumer interface over load().
- **Actual:** Per the header's own @note, transform() currently materializes the full dataset in memory before streaming to the consumer, so it has the same peak memory cost as load() despite the streaming-consumer signature. A caller choosing transform() specifically to bound memory will be surprised by the full-file allocation.
- **Evidence:** BrukerTimsFile.h:226-228 `/// @note Currently loads the full dataset into memory before streaming to the consumer. /// A future optimization should iterate frame-by-frame for true constant-memory operation.`
- **Fix:** The @note already documents this honestly (good). Keep the note prominent until true streaming exists; consider naming/marking it so callers don't assume constant-memory. No ABI change needed.
- **Verifier correction:** The behavior is accurately described, but the implied severity (a streaming API that silently allocates the full file) overstates the hazard. The non-streaming behavior is documented in three places at the point of use: the header @note (lines 226-228), a duplicate NOTE in the implementation (BrukerTimsFile.cpp:2239-2241), and a runtime OPENMS_LOG_INFO message ("loading full dataset (streaming optimization pending)", line 2242). Because it is loud, recoverable, and produces no incorrect data, this is a low-severity (mild) surprise rather than a high-severity silent hazard. The recommendation (keep the note prominent / mark it so callers don't assume constant memory; no ABI change) is sound. Optionally consider naming the runtime log a one-shot WARNING rather than INFO so memory-bounded callers notice it.
- **Verified:** Evidence is genuine and verified in both files. Header BrukerTimsFile.h:226-228 carries the quoted @note verbatim, and the implementation in BrukerTimsFile.cpp:2239-2256 confirms it: transform() builds a full MSExperiment via loadFrames_/loadDIA_/loadDDA_, sorts it, then iterates feeding the consumer — same peak memory as load(). The surprise is real against an established OpenMS convention: IMSDa
