# OpenMS POLS Audit — All HIGH-severity findings (62)

High = silently produces wrong results / data loss / crashes for reasonable use. Each verified against source.


## B1 (3)

### [DATA-24] DistanceMatrix::setValue
`hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h

**setValue only updates the cached minimum when neither index equals the current min coordinate, risking a stale minimum**

- Expectation: A setValue that advertises 'keep min_element_ up-to-date' should keep the cached minimum correct after any single write.
- Actual: The fast path `if (i != min_element_.first && j != min_element_.second)` updates min only when the new value is smaller. But when you overwrite a non-min cell that shares ONLY i (or ONLY j) with the min coordinate, the branch is still taken and a full updateMinElement() is never triggered even if you raise the previous-min via the other branch; the asymmetric guard makes the min-tracking guarantee non-obvious and easy to violate.
- Fix: Document precisely under which writes min_element_ stays valid and when callers must call updateMinElement(); or simplify the guard to compare the (i,j) pair against the min coordinate pair rather than each index independently. ABI: doc-only fix is non-breaking.
- Verifier note: The defect is not merely a non-obvious/doc-only guarantee: setValue silently leaves the cached minimum stale, violating the documented invariant ("setValue alone keeps min_element_ up to date"). Specifically, the else branch `if (value <= matrix_[min_element_.first][min_element_.second]) { matrix_[i][j] = value; }` writes a value into a cell (i,j) that shares exactly one index with the min coordin

### [DATA-25] DistanceMatrix::operator==
`surprising-throw` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h

**Equality operator throws (in debug) instead of returning false for differently-sized matrices**

- Expectation: operator== is expected to be total: comparing two matrices of different sizes should return false, never throw.
- Actual: It asserts size equality via OPENMS_PRECONDITION; in a debug build comparing different-sized matrices aborts/throws, and in release it then iterates rhs's dimensions reading this->matrix_ out of bounds.
- Fix: Return false when sizes differ before the loop: `if (dimensionsize_ != rhs.dimensionsize_) return false;`. ABI: additive/behavioral fix, source-compatible.
- Verifier note: operator== is non-total: in DEBUG it throws Exception::Precondition on size mismatch (this is documented at the call site, so half-expected), but in RELEASE the precondition macro is compiled out and the loop iterates over rhs.dimensionsize() while indexing this->matrix_[i][j]; when this is smaller than rhs this reads this->matrix_[i] beyond the allocated dimensionsize_ pointers -> heap out-of-bou

### [KERN-29] MassTrace::getIntensity
`misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/MassTrace.h

**getIntensity() returns a quantification value (FWHM area/median/height), not an intensity, and silently returns 0 if FWHM was never estimated**

- Expectation: A method named getIntensity(bool smoothed) returns the (raw or smoothed) intensity of the mass trace, e.g. the apex or summed intensity.
- Actual: It dispatches on the internal quant_method_: MT_QUANT_AREA returns computeFwhmArea()/computeFwhmAreaSmooth(), MT_QUANT_MEDIAN returns the median peak intensity, MT_QUANT_HEIGHT returns getMaxIntensity(). For AREA, computeFwhmArea() returns 0 unless estimateFWHM() was called first (it early-returns when fwhm_start_idx_==fwhm_end_idx_==0). The 'smoothed' flag is ignored entirely for MT_QUANT_MEDIAN.
- Fix: Rename to getQuantitatedIntensity()/computeQuantity() with a [[deprecated]] alias for getIntensity; document that AREA mode requires a prior estimateFWHM() call and that 'smoothed' is a no-op for MEDIAN. ABI: source-compatible if alias retained.
- Verifier note: All claim components verified. Severity is high: the silent-0 failure is reachable in-tree (MassTraceExtractor no-EPD path computes feature intensity 0 because getIntensity(false) runs before estimateFWHM under the default MT_QUANT_AREA), producing wrong quantitation with no error. Recommendation stands: rename to getQuantity()/computeQuantity() with a [[deprecated]] getIntensity alias (source-com


## B2 (4)

### [META-9] SpectrumMetaDataLookup::addMissingSpectrumReferences
`param-order-or-bool` · ABI: `breaking` · src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h

**Documented [in,out] 'proteins' is taken by value, so spectra_data updates are silently discarded**

- Expectation: The Doxygen marks proteins as '@param[in,out] proteins Protein IDs corresponding to the Peptide IDs' and says with override_spectra_data the 'ProteinIdentifications should be updated with new spectra_data values'. A caller passing its protein vector expects the prot.setMetaValue("spectra_data", ...) writes to be visible after the call.
- Actual: The parameter is declared by value (std::vector<ProteinIdentification> proteins), not by reference. In SpectrumMetaDataLookup.cpp the function loops 'for (auto& prot : proteins) prot.setMetaValue("spectra_data", spectra_data);' but mutates only the local copy, which is destroyed on return. The caller's protein IDs are never updated.
- Fix: Change the parameter to std::vector<ProteinIdentification>& (true in/out), or if ABI must be preserved add an overload taking a reference and deprecate the by-value one. At minimum, fix the doc to stop claiming [in,out]. The current state is a silent bug at real call sites (IDFileConverter).

### [META-12] SpectrumMetaDataLookup::getSpectrumMetaData
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h

**Silently returns leaving 'meta' untouched when no reference format matches, unlike findByReference which throws ParseError**

- Expectation: Given the class's other parse path SpectrumLookup::findByReference throws Exception::ParseError when 'no reference format matched', and this function's doc only mentions '@throw ElementNotFound if a spectrum look-up was necessary', a caller passing an unparseable reference expects either an exception or a documented signal of failure.
- Actual: If none of reference_formats matches the spectrum_ref, the for-loop completes without entering the body and the function returns normally, leaving the output 'meta' at its default-constructed values (rt/precursor_rt/precursor_mz = NaN, charge=0, ms_level=0, scan_number=-1). No exception, no return value to inspect.
- Fix: Throw Exception::ParseError on no-match (mirroring findByReference) or document the silent no-op and that callers must check meta for sentinel/NaN values. Behavioral fix; flag abi as source-compatible since it adds a throw.
- Verifier note: Field is `precursor_charge` (not `charge`); it defaults to 0. All other claimed sentinels are exact. Severity kept high: on an unparseable reference the function silently returns leaving `meta` fully at sentinels (NaN RT/mz, ms_level=0 which is an invalid MS level, scan_number=-1), with no exception and no return value — a caller that does not defensively test for NaN/sentinels will silently propa

### [CHEM-24] IMSIsotopeDistribution::size
`misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h

**size() caps the reported peak count at the global static SIZE rather than reporting the actual peak count**

- Expectation: size() on a container-like class returns the number of stored elements, consistent with empty() returning peaks_.empty().
- Actual: size() returns min(peaks_.size(), SIZE) where SIZE is a mutable static class member shared across all instances. An instance holding 12 peaks reports size()==SIZE (e.g. 10), so iterating [0,size()) skips real peaks; the value also changes globally if any code reassigns the static SIZE. empty() (peaks_.empty()) and size() can therefore disagree about how many peaks are accessible.
- Fix: ABI-safe: fix the doc comment (it currently claims the wrong inequality) and clearly state size() is clamped to the static SIZE truncation length, distinct from peaks_.size(). Consider adding a separate peakCount()/rawSize() accessor (additive) for the true count.
- Verifier note: size() does not merely cap at "e.g. 10"; the static SIZE is default-initialized to 0 and is never assigned anywhere in the repo, so size() returns 0 for any distribution constructed from peaks (until a fold operation grows peaks_ to SIZE or a client mutates the static). This yields silently empty results from getMasses()/getAbundances()/operator<</pyOpenMS len() despite empty()==false. The doc com

### [META-31] ExperimentalDesign::getSample
`silent-failure` · ABI: `none` · src/openms/source/METADATA/ExperimentalDesign.cpp

**getSample dereferences find_if result without checking end() (UB / crash on missing combination)**

- Expectation: A lookup that takes a fraction_group and label and 'returns the sample index' would either return the index, throw a clear ElementNotFound, or return a sentinel when no matching row exists.
- Actual: The body is `return std::find_if(msfile_section_.begin(), msfile_section_.end(), ...)->sample;` with no comparison to end(). If no entry matches (fraction_group, label), the end iterator is dereferenced, which is undefined behavior (out-of-bounds read or crash), not a diagnosable error.
- Fix: Capture the iterator, compare to end(), and throw Exception::ElementNotFound (mirroring the rest of the class) when not found. Additive/behavioral fix; does not change the signature so ABI is preserved.


## B3 (14)

### [FORM-9] Base64::decodeIntegersUncompressed_ (reached via decodeIntegers, zlib_compression=false)
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/Base64.h

**Integer decode silently returns an empty/partial result for inputs shorter than 4 chars instead of signalling, and (unlike the float path) does not validate length % 4**

- Expectation: Symmetric with decode(): malformed Base64 (length not a multiple of 4) should be reported the same way for integers as for floats.
- Actual: decodeUncompressed_ (float path) throws ConversionError when `in.size() % 4 != 0` (Base64.h lines 320-323). But decodeIntegersUncompressed_ (integer path) has NO such check — for `in.size() < 4` it just returns leaving `out` empty (lines 513-516) and otherwise proceeds, so a malformed integer buffer is decoded into garbage or a short result with no error. The two sibling decoders are inconsistent: one validates length, the other does not.
- Fix: Add the same `in.size() % 4 != 0` length check (and matching ConversionError) to decodeIntegersUncompressed_ so the integer and float decoders fail consistently. Source-compatible additive validation; only affects previously-silently-accepted malformed inputs.
- Verifier note: The asymmetry is real (float decoder validates length % 4 and throws; integer decoder does not). However the consequence is worse than "garbage or short result with no error": because the loop reads in[i+2] and in[i+3] before its bounds guards, a malformed integer buffer whose length mod 4 is 1 or 2 causes an out-of-bounds read of `in` past size() (undefined behavior) and a secondary out-of-bounds

### [FORM-12] CsvFile (class doc) vs CsvFile::load / constructor
`misleading-doc` · ABI: `none` · src/openms/include/OpenMS/FORMAT/CsvFile.h

**CsvFile doc claims no comment support, but '#'-prefixed lines are silently dropped**

- Expectation: Given the class documentation 'Does NOT support comment lines!', a caller expects every line of the CSV to be loaded verbatim, including any line that happens to start with '#'.
- Actual: Both the constructor and load() pass `"#"` as the comment_symbol to TextFile::load, so every line beginning with '#' is silently skipped and never appears in the buffer or in rowCount(). A field-less CSV whose first column legitimately starts with '#' loses rows with no error.
- Fix: Either honor the documentation by passing an empty comment_symbol, or update the documentation to state that '#'-prefixed lines are treated as comments and skipped. Prefer making the comment behavior an explicit, documented parameter. ABI-safe (constant/string change or additive param).
- Verifier note: The category is not "misleading-name" (no symbol is misnamed); it is a documentation/contract contradiction. CsvFile's class doc states "Does NOT support comment lines!" while its constructor and load() unconditionally pass "#" as the comment_symbol to TextFile::load (CsvFile.cpp:27,34), and TextFile::load (TextFile.cpp:66-69) silently skips any line prefixed with "#" — not added to buffer_, not c

### [FORM-22] IndexedMzMLFileLoader::load
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/IndexedMzMLFileLoader.h

**load() returns bool success that is easy to ignore; failure produces no exception**

- Expectation: Sibling *File classes (MzMLFile, MzDataFile, MzXMLFile) all declare load() as 'void' and signal a missing/garbage file by throwing FileNotFound/ParseError. A caller used to those classes would write `loader.load(fn, exp);` and assume an exception on failure.
- Actual: load() returns a bool (forwarded from exp.openFile(filename)) and throws nothing on failure. The doc even says 'Tries to parse the file, success needs to be checked with the return value.' A caller who ignores the return value silently proceeds with an empty/invalid OnDiscPeakMap.
- Fix: Keep the bool (ABI), but document the divergence from the other *File classes prominently and mark the return value [[nodiscard]] (source-compatible additive change) so callers cannot silently drop the failure indication.
- Verifier note: Severity confirmed high: a caller migrating from the sibling void-load classes writes `loader.load(fn, exp)` expecting an exception on a missing/garbage file, but instead receives `false` (commonly ignored) plus an empty OnDiscPeakMap and silently proceeds with no data — silently wrong results / data loss for reasonable use, evidenced by the in-tree TICCalculator caller that ignores the return val

### [FORM-33] MSPFile::load(const std::string&, PeptideIdentificationList&, PeakMap&)
`asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MSPFile.h

**load() resets the output PeakMap but appends to the output ids list without clearing it**

- Expectation: A load() that fills two output containers should treat them consistently — both cleared before being populated, so reusing the same objects across two load() calls yields only the second file's content.
- Actual: exp is cleared via `exp.reset()` (MSPFile.cpp:72), but ids is never cleared; the loop only does `ids.push_back(id)` (MSPFile.cpp:119). Reusing the same PeptideIdentificationList across two load() calls silently accumulates identifications from both files, desynchronizing ids from exp.
- Fix: Add `ids.clear();` next to `exp.reset();` at the top of load(). Behavior-fix in the .cpp, no signature/ABI change. If preserving append semantics is intentional, document it explicitly in the header.
- Verifier note: No correction needed; claim is accurate. Severity confirmed high because object reuse across files is a reasonable usage pattern and the failure mode is silent index-level desynchronization of ids vs exp (no throw, no log), i.e. silently-wrong results. The recommended fix (add ids.clear() next to exp.reset() at the top of load()) is correct and a behavior-only .cpp change with no ABI impact. Note:

### [FORM-38] XMassFile::load
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/XMassFile.h

**load() reads fid intensities with no EOF/short-file check; truncated fid yields zero-padded peaks instead of an error**

- Expectation: load() of a truncated or corrupt fid file (fewer intensity values than acqus declares) should signal a parse/IO error.
- Actual: The loop runs `while (spectrum.size() < acqus.getSize())` and pushes a peak for every index up to the acqus-declared size, calling FidHandler::getIntensity(). getIntensity() does `read(...,4)` with no success check and returns 0 on a failed/EOF read (`return (result > 0) ? result : 0`). So a fid shorter than acqus.getSize() silently produces a full-length spectrum padded with zero-intensity peaks; no exception is thrown. Only a completely unopenable fid triggers an exception (the early `if (!fid)` block).
- Fix: Have getIntensity() check the stream state after read() and have load() stop / throw Exception::ParseError when the fid stream reaches EOF before acqus.getSize() peaks are read, rather than zero-padding. The fix is internal to the .cpp/handler and ABI-neutral at the XMassFile API level.

### [FORM-42] MzIdentMLFile::load
`asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzIdentMLFile.h

**MzIdentMLFile::load APPENDS to the output containers, while sibling *File::load methods CLEAR them first**

- Expectation: A competent caller expects load() to populate the output vectors with the file's contents, i.e. to clear/replace whatever was passed in. This is reinforced by the class doc: 'The file adapter interface is kept the same as idXML file adapter for downward capability reasons.' IdXMLFile::load (the adapter it claims to mirror) explicitly clears: 'protein_ids.clear(); peptide_ids.clear();'. PepXMLFile::load and PepXMLFileMascot::load and PTMXMLFile::load all clear their outputs too.
- Actual: MzIdentMLFile::load does NOT clear poid/peid. It forwards them to MzIdentMLDOMHandler, which only ever appends via pro_id_->emplace_back() / pep_id_->push_back() (e.g. MzIdentMLDOMHandler.cpp lines 884, 1368, 1811, 2074). readMzIdentMLFile contains no .clear(). Re-using the same vectors (or passing a non-empty vector) silently accumulates duplicate identifications.
- Fix: Additive, ABI-safe fix: clear poid and peid at the top of MzIdentMLFile::load before constructing the handler (matching IdXMLFile/PepXMLFile). This is a behavior change but the correct/consistent one; document it. No signature/ABI change.

### [FORM-60] InspectOutfile::load
`asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/FORMAT/InspectOutfile.h

**load() appends to the output identification containers without clearing them**

- Expectation: A load() that produces peptide_identifications / protein_identification is expected to populate fresh results; reusing the same objects across two load() calls should not merge two files' hits.
- Actual: load() only push_back/insertHit onto the passed containers; it never clears them. InspectOutfile.cpp:183 'protein_identification.insertHit(protein_hit);', :213 and :269 'peptide_identifications.push_back(peptide_identification);'. A second load() into the same containers accumulates hits from both files. This is the same inconsistency as OMSSACSVFile but differs from PercolatorOutfile which clears.
- Fix: Document the append behavior or clear the output containers at entry to match the clearing loaders. If callers rely on append, at minimum state it in the header. Clearing is a source-compatible behavior change; abi none.

### [FORM-66] MzTabBoolean::setNull
`misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/FORMAT/MzTabBase.h

**MzTabBoolean::setNull has inverted polarity vs every sibling setNull**

- Expectation: setNull(true) marks the cell as null; setNull(false) clears the null state. This is the contract of MzTabString/MzTabInteger/MzTabDouble/MzTabModification::setNull in the same module.
- Actual: Implementation is `if (!b) value_ = -1; else value_ = 0;`. So setNull(false) MAKES the value null (-1) and setNull(true) makes it NOT null (and sets it to boolean false/0). The boolean argument is inverted relative to its name and to all sibling classes.
- Fix: Fix the body to `value_ = b ? -1 : 0;` so it matches the sibling convention. This is a source-compatible behavior fix (the signature is unchanged); callers currently calling setNull(true) to nullify are themselves buggy. Guard with a test. If existing callers depend on the inverted behavior, add a deprecation note but correct the polarity.

### [FORM-69] MzTabFile::load
`asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/FORMAT/MzTabFile.h

**load/store round-trip silently drops nucleic-acid, oligonucleotide and OSM sections**

- Expectation: store() writes NUC/OLI/OSM sections (the header and store() fully support MzTabNucleicAcid/Oligonucleotide/OSM rows), so a load() of a file written by store() should restore them - a symmetric round-trip.
- Actual: load() only assigns metadata, protein, peptide, PSM and small-molecule sections back into the MzTab object. It never calls setNucleicAcidSectionRows / setOligonucleotideSectionRows / setOSMSectionRows, so those sections are silently discarded on read even though store() emits them. store()->load() is lossy with no error.
- Fix: Either implement parsing for the NUC/OLI/OSM sections in load(), or at minimum emit a warning/throw when such sections are encountered so the data loss is not silent. Additive to load() impl; no ABI change.

### [FORM-94] XICParquetFile::getChromatograms
`param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FORMAT/XICParquetFile.h

**Six adjacent same-typed Int64 filter parameters (all defaulting to -1) are trivially swappable at call sites**

- Expectation: Filter arguments should be hard to confuse; a developer calling getChromatograms(out, 1318, 7) expects an obvious mapping of values to columns.
- Actual: getChromatograms has precursor_id, transition_id, precursor_charge, product_charge, ms_level, run_id as six positional Int64 params all defaulting to -1, interleaved with two string params. Passing e.g. a run_id where a transition_id is expected, or swapping precursor_id/transition_id, compiles cleanly and silently filters on the wrong column, returning a wrong (often empty) result set with no error.
- Fix: Promote the typed ParquetFilter / ParquetFilterBuilder overloads (already present) as the primary API and deprecate or de-emphasize the positional overload; or introduce a small XICQuery struct of named fields. Additive (new overload/struct) keeps ABI; the existing overload can stay for compatibility.

### [FORM-108] RationalScan2ImConverter::convert / getCalibration
`silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/RationalScan2ImConverter.h

**Unknown/out-of-range frame_id silently falls back to the FIRST calibration instead of failing**

- Expectation: Converting with a frame_id that has no calibration (out of range, or absent from the map) would be a hard error, since using a wrong frame's calibration yields silently wrong ion-mobility values.
- Actual: getCalibration() returns calibrations_.begin()->second (the first/arbitrary calibration) for any frame_id not found, logging only a WARN. convert() then produces plausible-looking but physically wrong 1/K0 values. The header documents Coefficients/strategy but the public convert() override gives no hint that it tolerates bad frame_ids. Additionally, if the converter was constructed (its ctor is public) with an empty calibrations_ map, calibrations_.begin() dereferences end() — undefined behavior.
- Fix: Throw on unknown frame_id (or at least guard calibrations_.begin() against an empty map and document the silent-fallback policy). If fallback must stay for robustness, state it explicitly in the convert() doc.

### [FORM-109] BrukerTimsImagingFile::load
`silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/BrukerTimsImagingFile.h

**load() silently drops MALDI pixels whose frame id is absent from the loaded spectra (warn-only)**

- Expectation: A load() that builds an imaging geometry from MaldiFrameInfo would either map every declared pixel or fail when pixels cannot be matched to spectra, so the resulting image is complete.
- Actual: Pixels in MaldiFrameInfo whose frame_id has no corresponding loaded spectrum are skipped and only counted into a one-shot WARN; load() still returns 'successfully' with an image that is missing those pixels. With a frame_id_min/frame_id_max filter set in inner_config this is the common case, and the geometry silently omits the out-of-range pixels. The header's @throws list does not mention any partial-result condition.
- Fix: Document the warn-and-drop partial-image behavior in load()'s header contract (it currently lists only throwing failure modes), and consider returning the dropped count or offering a strict mode that throws when any declared pixel is unmatched.

### [FORM-113] GNPSQuantificationFile::store
`surprising-throw` · ABI: `none` · src/openms/source/FORMAT/GNPSQuantificationFile.cpp

**store() indexes consensus_map[0] unconditionally, throwing/UB on an empty ConsensusMap**

- Expectation: Storing an empty ConsensusMap writes a header-only (or empty) table, not a crash.
- Actual: The first statement dereferences consensus_map[0] to probe for the IIMN_ROW_ID meta value before any size check, so calling store() with an empty map indexes out of bounds.
- Fix: Guard with 'if (!consensus_map.empty() && consensus_map[0].metaValueExists(...))'. Additive, no signature change.
- Verifier note: store() dereferences consensus_map[0] before any size check to probe for IIMN_ROW_ID. ConsensusMap::operator[] is noexcept and forwards to std::vector::operator[], so on an empty map this is undefined behavior (out-of-bounds read), not a thrown exception. Fix: guard with `if (!consensus_map.empty() && consensus_map[0].metaValueExists(...))` so an empty map writes the header-only table as documente

### [FORM-129] MSChromatogramParquetConsumer::~MSChromatogramParquetConsumer / finalize
`silent-failure` · ABI: `none` · src/openms/include/OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h

**Destructor writes the file but swallows all write errors; only finalize() surfaces them**

- Expectation: Constructing the consumer and letting it go out of scope (the RAII-natural usage) would be expected to either write the file successfully or make a failure visible.
- Actual: The destructor calls finalize() and catches every exception, logging to OPENMS_LOG_ERROR and continuing (cpp:191-212). If the user never calls finalize() explicitly, a failed/partial Parquet write is silently downgraded to a log line. The header does flag this ("Call this explicitly to surface write errors during normal control flow"), so the surprise is partially mitigated by docs.
- Fix: Keep the swallow-in-destructor (throwing from a destructor is worse), but the header could state more prominently that relying on destructor-only finalization hides I/O errors. Documentation-only; no ABI impact.
- Verifier note: The destructor genuinely swallows all Parquet write/close errors (cpp:191-212), and write_() performs terminal I/O including the metadata-writing Close() calls (cpp:506,519). The claim's recommendation ('Documentation-only; no ABI impact', severity ~medium) understates the issue: finalize() is neither part of the IMSDataConsumer interface nor virtual, so the OpenSwath workflow — which holds the co


## B4a (5)

### [ANID-10] IDBoostGraph::getProteinScores_ / GetScoreTgTVisitor
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/IDBoostGraph.h

**Missing target_decoy meta value silently classifies a protein/peptide as decoy instead of signalling 'unknown'**

- Expectation: A hit whose target_decoy status is unknown/undefined should be reported as unknown (the GetScoreTgTVisitor docstring even promises '(-1.0, false)' for that case), or an error should be raised so the caller knows the FDR input is incomplete.
- Actual: Both getProteinScores_ and the GetScoreTgTVisitor read `getMetaValue("target_decoy").toString()[0] == 't'`. getMetaValue returns DataValue::EMPTY when the key is absent (MetaInfoInterface.cpp:118-125), DataValue::toString() of EMPTY is the empty string (DataValue.cpp:702-708), and indexing [0] of an empty std::string yields '\0' which is != 't'. The hit is therefore silently counted as a *decoy* (label 0) and still fed into FDR, rather than being flagged. The documented '(-1.0, false)' fallback only fires for non-hit node types, never for a real ProteinHit/PeptideHit that simply lacks the annotation.
- Fix: Use getMetaValue("target_decoy", DataValue::EMPTY) and explicitly test for EMPTY, then either throw Exception::MissingInformation or surface a tri-state; at minimum stop scoring an un-annotated hit as a confident decoy. ABI-safe (implementation-only change).

### [ANID-13] AA::AA(const char)
`silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h

**AA(char) reads out of bounds for bytes >= 128 instead of yielding an invalid AA**

- Expectation: The constructor doc states 'All other chars produce an invalid AA (?)', so any char value should map to a valid or the invalid('?') AA without undefined behavior.
- Actual: CharToAA is declared with exactly 128 entries (`constexpr char const CharToAA[128]`), but the constructor indexes it with `(unsigned char)c`, which ranges 0..255. Any non-7-bit byte (extended ASCII, UTF-8 continuation bytes, or a signed char that was >127) indexes past the end of the array — an out-of-bounds read rather than the promised invalid AA. The table comment itself says it only covers '7-bit ASCII'.
- Fix: Either size the table to 256 (filling 128..255 with the invalid value 27) or mask with `c & 0x7F` only after confirming high bit clear / branch to invalid for c<0 or >=128. The 256-entry table is the additive, ABI-safe fix.

### [ANID-34] NeighborSeq::NeighborSeq
`ownership-lifetime` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h

**rvalue-ref constructor stores a const-reference (no move); passing a temporary dangles immediately**

- Expectation: A constructor taking std::vector<AASequence>&& signals it takes ownership by moving the vector in, so NeighborSeq ns(makeDigest()); is the idiomatic call and the object owns its data for its lifetime.
- Actual: The member is `const std::vector<AASequence>& digested_relevant_peptides_;` (a reference) and the ctor does `digested_relevant_peptides_(std::move(digested_relevant_peptides))`. Binding a const& member to the rvalue-ref parameter does NOT extend lifetime and does NOT move — it just aliases the caller's object. Pass a temporary (exactly what `&&` invites) and the reference dangles the instant the constructor returns; every later isNeighborPeptide()/getNeighborStats() call is UB.
- Fix: Either actually own the data (member `std::vector<AASequence> digested_relevant_peptides_;` initialized from std::move — ABI/source change but matches the `&&` contract and removes the dangling trap), or change the parameter to `const std::vector<AASequence>&` so the signature stops promising a move. The first is the least-surprising fix. abi_impact breaking due to member-layout/signature change; do it on the next ABI break, and in the meantime the @note is the only guard.

### [ANID-36] SimpleSearchEngineAlgorithm::search
`asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h

**Doc claims outputs 'are not cleared', but prot_ids is overwritten while pep_ids is appended**

- Expectation: The header explicitly states: 'Existing contents of @p prot_ids and @p pep_ids are not cleared by this call.' A caller would expect to accumulate results from several search() calls into the same two vectors.
- Actual: postProcessHits_ does `protein_ids = vector<ProteinIdentification>(1);` — it destroys/replaces any pre-existing prot_ids — while pep_ids are only appended (push_back). So the documented contract is true for pep_ids but FALSE for prot_ids, and the two outputs behave asymmetrically. Worse, it then re-stamps the identifier on ALL pep_ids (`for (auto & pid : peptide_ids) pid.setIdentifier(identifier);`), rewriting identifiers of any pre-existing PSMs.
- Fix: Fix the documentation to state that prot_ids is replaced (and pre-existing PSM identifiers are re-stamped), OR make the behavior symmetric (append a new ProteinIdentification run instead of overwriting). Documentation-only fix is non-breaking; making it symmetric is a behavior change. Flagging the ideal as: append a run so both outputs accumulate.

### [ANID-41] NeighborSeq::NeighborStats::unfindable / noNB / oneNB / multiNB
`surprising-throw` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h

**Percentage formatters divide by total() and crash (integer divide-by-zero) on an empty/default-constructed NeighborStats**

- Expectation: Trivial const string formatters on a small POD struct should be safe to call on any instance, including a default-constructed/empty one.
- Actual: Each formatter computes `count * 100 / total()`; total() is the sum of the four counters, which is 0 for a default-constructed NeighborStats (or any case with no registered peptides). Integer division by zero is undefined behavior (SIGFPE on typical hardware). The header even @warning-documents it but ships it anyway.
- Fix: Guard total()==0 (return '0 (0%)' or 'n/a') inside each formatter. Purely additive, non-breaking, and removes a crash on a documented-but-live code path.
- Verifier note: Minor refinement, not a rejection: the formatters do not 'throw' in the C++ exception sense; they trigger SIGFPE via integer division by zero (UB) when total()==0. The recommended fix (guard total()==0 inside each inline formatter, returning e.g. "0 (0%)" or "n/a") is purely additive to inline bodies and changes no signatures, struct layout, or exported symbols, hence abi_impact=none. The reachabl


## B4b (7)

### [ANSW-1] SwathWindowLoader::readSwathWindows
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h

**readSwathWindows silently returns empty vectors when the file is missing or unreadable**

- Expectation: A static reader named readSwathWindows, taking a filename, throws (e.g. FileNotFound) when the file does not exist or cannot be opened, so a caller that mistyped a path or pointed at a missing window file learns about it.
- Actual: The implementation does `std::ifstream data(filename.c_str());` and immediately `std::getline(data, line);` with no open()/good() check. A missing or unreadable file produces a failed stream; the getline fails, the do/while body never runs (or runs once on garbage), and the function returns leaving both output vectors empty with no error. The header documents this: 'A missing or unreadable file is not reported as an exception; the output vectors are left empty.'
- Fix: Additive/source-compatible fix: after constructing the ifstream, check `if (!data) throw Exception::FileNotFound(...)`. The sibling annotateSwathMapsFromFile would then surface a clear error for a bad path instead of failing the downstream count check with a confusing message. If preserving silent behavior is mandatory for ABI/back-compat, at minimum emit a WARN; the doc already flags this as surprising.
- Verifier note: readSwathWindows does lack any open/good check on the ifstream (real silent-failure defect), but the precise consequence claimed is incorrect. For a missing/unreadable file, getline fails leaving line empty; StringUtils::split on the empty string returns early leaving headerSubstrings empty (StringUtils.h:761); the immediately following StringUtils::toDouble(headerSubstrings[0]) at SwathWindowLoad

### [ANSW-5] SwathQC::getSpectraProcessingFunc
`hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SwathQC.h

**Uniform subsampling index counts only sampled spectra, not all spectra seen, biasing which spectra get analyzed**

- Expectation: The member-function path is documented as the speed-efficient way to 'sample the spectra as they are loaded'; uniform subsampling implies the decision uses the running index of ALL MS1 spectra streamed in.
- Actual: The returned lambda calls isSubsampledSpectrum_(nr_ms1_spectra_, cd_spectra_, ms1_spectra_seen_) and only increments ms1_spectra_seen_ AFTER a spectrum passes the filter (it is incremented inside the accepted branch, line 69, never on rejection). So ms1_spectra_seen_ is the count of ACCEPTED spectra, but it is used as the global spectrum index for uniform spacing. Because rejected spectra never advance the index, the index never reaches the spacing thresholds for later spectra the way the algorithm assumes — sampling is front-loaded/non-uniform. The private member is even documented as 'keeps track of number of spectra passed to getSpectraProcessingFunc()', which is not what it counts. The static getChargeDistribution avoids the bug entirely by doing its own indexing with the true loop variable i (SwathQC.cpp:131) and forcing sample-all mode.
- Fix: Increment a true 'seen' counter for every MS1 spectrum entering the lambda and pass THAT as idx to isSubsampledSpectrum_ (separate from the accepted count). This is an internal-impl fix with no signature/ABI change. Note this only manifests when nr_ms1_spectra_>0 (external getExpSettingsFunc/setNrMS1Spectra path); the in-tree static caller is unaffected, which is likely why it went unnoticed.
- Verifier note: Minor correction to the claim's 'went unnoticed' explanation: it is not strictly true that 'the in-tree caller is unaffected.' A member-path unit test exists (SwathQC_test.cpp storeJSON section) and exercises nr_ms1_spectra_>0, but it is masked because the test data has only 3 MS1 spectra against a subsample target of 10 (subsample>=total => isSubsampledSpectrum_ always returns true, so ms1_spectr

### [ANSW-26] SpectrumAccessSqMass::SpectrumAccessSqMass(handler, indices)
`surprising-default` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h

**Passing an empty index vector silently means 'all spectra', not 'no spectra'**

- Expectation: Constructing a subset view with an explicit list of indices and passing an empty list yields an empty view (zero visible spectra) — the natural reading of 'expose only the spectra at the given indices'.
- Actual: An empty `indices` vector is treated as the sentinel for 'all spectra'. getNrSpectra() then returns the full file count, not 0. A caller that computed an empty index set (e.g. a SWATH window that matched nothing) silently gets the entire file instead of an empty selection.
- Fix: The behavior is documented, but the overload signature does not telegraph it. Prefer a named factory or an explicit 'all spectra' default-constructed marker so an empty selection cannot be confused with 'all'. Additive: add a factory method; do not change existing ctor semantics (ABI/behavior depended on by ChromatogramExtractor). Flag honestly: the ideal fix (empty => empty) would be breaking.
- Verifier note: Two minor refinements to the claim, neither undermining the surprise: (1) The (parent, indices) ctor's empty case (cpp 30-33) inherits the PARENT's subset unchanged (sidx_ = sp.sidx_), not literally 'all file spectra' — though if the parent itself had an empty subset it transitively means all-file, so the empty==don't-narrow spirit is the same. (2) The claim that the behavior is 'depended on by Ch

### [ANSW-41] MRMFeatureFilter::calculateIonRatio
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h

**calculateIonRatio returns 0.0 sentinel on missing data and divides without guarding against zero denominator**

- Expectation: A function named calculateIonRatio returns the ratio of two transition intensities; a caller would expect either a valid ratio or a clear error/NaN when it cannot be computed.
- Actual: When neither component has the requested key/native_id the function silently returns the sentinel 0.0 (indistinguishable from a real ratio of 0). When only component_1 has the value it returns component_1's raw value (NOT a ratio). When both values exist it computes `feature_1 / feature_2` with no check that feature_2 != 0, so a zero denominator yields +inf/NaN with no signal.
- Fix: Document the 0.0/single-value fallbacks explicitly, and guard the division (return NaN or throw on feature_2==0). Additive/source-compatible fix: keep signature, add a denominator==0 guard returning std::numeric_limits<double>::quiet_NaN() and document the missing-key behavior; ABI unchanged.
- Verifier note: Minor refinement: for the "intensity" branch the missing-data gate tests the "native_id" metavalue, not the intensity value (intensity is always available via getIntensity()); for other feature_names it tests the actual key. The recommended additive fix (denominator==0 guard returning quiet_NaN + documenting fallbacks) keeps the signature, so ABI is unchanged but it is a behavior change, hence sou

### [ANSW-46] ReactionMonitoringTransition::getPrediction
`silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/MRM/ReactionMonitoringTransition.cpp

**getPrediction()'s guard checks the WRONG predicate (hasPrecursorCVTerms instead of hasPrediction), and in release builds dereferences a possibly-null pointer**

- Expectation: A caller who follows the header's instruction ('You first need to check whether the object is accessible using hasPrediction()') expects getPrediction() to be safe whenever hasPrediction() is true, and to be guarded against the null case.
- Actual: The precondition is `OPENMS_PRECONDITION(hasPrecursorCVTerms(), "...has no Prediction object, check first with hasPrediction()")` -- it checks hasPrecursorCVTerms(), an unrelated member, due to copy-paste from getPrecursorCVTermList(). The function then `return *prediction_;`. OPENMS_PRECONDITION compiles to nothing unless OPENMS_ASSERTIONS is set, so in a normal release build getPrediction() unconditionally dereferences prediction_, which is nullptr by default. A transition with no prediction (the default) gives an immediate null-deref crash, and even a debug build checks the wrong flag.
- Fix: Fix the precondition to `OPENMS_PRECONDITION(hasPrediction(), ...)`. ABI-safe (function body change only). Ideally also harden against null in release (e.g. throw Exception::InvalidValue when prediction_ == nullptr) for parity with the documented contract; that is also source/ABI compatible.
- Verifier note: Claim is accurate as stated. Minor clarification: even in a DEBUG (OPENMS_ASSERTIONS) build the bug is harmful — the wrong predicate can throw spuriously when a prediction exists but precursor CV terms do not, and can fail to throw (then null-deref) when a prediction is absent but precursor CV terms are present. In RELEASE the precondition compiles to nothing, so any default transition (prediction

### [ANSW-47] TargetedExperiment::getProteinByRef / getPeptideByRef / getCompoundByRef
`silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/TARGETED/TargetedExperiment.cpp

**By-ref getters return a reference into a map populated with operator[], silently inserting/returning a null entry for an unknown ref in release builds**

- Expectation: getProteinByRef(ref) for an unknown ref should fail loudly (throw / assert), as the message implies; at worst it should not corrupt internal state.
- Actual: The only guard is `OPENMS_PRECONDITION(protein_reference_map_.contains(ref), "Could not find protein in map")`, which is compiled out in release. The next line is `return *(protein_reference_map_[ref]);`. For an unknown ref, std::unordered_map::operator[] INSERTS a default-constructed entry (a null `const Protein*`) and then the code dereferences that nullptr -> crash/UB. So an unknown ref is not a clean lookup failure but a null-deref plus a mutation of the cache map. getPeptideByRef/getCompoundByRef have the identical pattern.
- Fix: Replace operator[] with `.at(ref)` or an iterator find + explicit throw (e.g. Exception::ElementNotFound) so an unknown ref is a deterministic error in all build types and does not mutate the cache. Body-only change; ABI-safe. Callers already must check has*() first, so behavior for valid refs is unchanged.

### [ANSW-77] TargetedSpectraExtractor::annotateSpectra
`unit-or-index` · ABI: `none` · src/openms/source/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.cpp

**ppm 'mz_tolerance' is divided by 1e6 but never multiplied by m/z, yielding a tolerance ~1e6x too small in ppm mode**

- Expectation: With the parameter in ppm mode (mz_unit_is_Da_ == false), a tolerance of e.g. 20 ppm at m/z 500 should produce an absolute window of about +/-0.01 Th (ppm * mz / 1e6).
- Actual: The code computes 'mz_tolerance = mz_tolerance_ / 1e6' and then uses it directly as an absolute window around spectrum_mz, omitting the multiplication by the m/z. A 20 ppm setting becomes +/-2e-5 Th regardless of m/z, so essentially no spectra match. The same pattern appears for fwhm_threshold ('fwhm_threshold_ / 1e6' at line 462).
- Fix: In ppm mode compute 'mz_tolerance = mz_tolerance_ * spectrum_mz / 1e6' (ppm is relative to the measured m/z). Behavioral fix, ABI-neutral. Add a regression test asserting the ppm window scales with m/z.
- Verifier note: In ppm mode (mz_unit_is_Da_ == false), line 370 should compute the absolute tolerance as mz_tolerance_ * spectrum_mz / 1e6 (i.e. Math::ppmToMass(mz_tolerance_, spectrum_mz)), not mz_tolerance_ / 1e6. As written, the matching window is the m/z-independent value mz_tolerance_/1e6 Th, which is too small by a factor equal to the precursor m/z (roughly 100-2000x in practice, NOT ~1e6x as the title stat


## B4c (8)

### [ANAL-1] FeatureGroupingAlgorithmUnlabeled::getResultMap
`asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h

**Incremental path (setReference/addToGroup/getResultMap) yields a different, un-postprocessed result than group()**

- Expectation: The header documents two equivalent ways to run: (a) group(), (b) setReference + addToGroup* + getResultMap. A caller expects getResultMap() to return the same finished consensus map that group() produces.
- Actual: group() ends with postprocess_(maps, out) which transfers protein IDs, re-indexes unassigned peptide IDs, and applies canonical sorting (sortByQuality/sortByMaps/sortBySize), and also fixes up column headers. The incremental path (addToGroup in the .cpp) never calls postprocess_ and never populates ColumnHeaders, so getResultMap() returns raw, unsorted pairfinder output with no protein IDs and empty/partial column headers. The TOPP tool FeatureLinkerUnlabeled.cpp has to manually patch this (out_map.setColumnHeaders(dummy.getColumnHeaders()); transferSubelements(...)). A caller following the header's 'two ways' contract gets a silently inferior result.
- Fix: Document that getResultMap() returns the raw incremental result and that the caller must perform postprocessing/header setup themselves (as FeatureLinkerUnlabeled does), or add a finalizeResult() helper. Do not silently change getResultMap to run postprocess_ (would change existing TOPP output ordering). Doc fix is source-compatible.
- Verifier note: The incremental path (setReference/addToGroup*/getResultMap) performs the SAME StablePairFinder linking as group(), so the consensus features themselves are equivalent. However, getResultMap() returns the raw pairfinder result WITHOUT the postprocessing that group() applies via postprocess_(): it lacks transferred protein IDs, lacks re-indexed unassigned peptide IDs (map_index), is not canonically

### [ANAL-9] MapAlignmentEvaluationAlgorithmRecall::evaluate
`silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.cpp

**evaluate() divides by ground-truth size (and by cons_map_tool.size()) with no guard, returning NaN/inf instead of signaling on empty input**

- Expectation: A caller expects a recall value in [0,1], or a clear error when the ground truth or the tool's consensus map is empty.
- Actual: Same pattern as Precision: `out = (1.0 / double(cons_map_gt.size())) * sum;` divides by zero when no GT consensus feature has size >= 2. Additionally, line 85 `gt.push_back(gt_i / cons_map_tool.size());` divides by cons_map_tool.size(), which is zero (UB / crash or trap) when the input tool map is empty. Neither case is checked or signaled.
- Fix: Guard against empty cons_map_gt and empty cons_map_tool before dividing (throw or return a documented sentinel). Additive early checks, ABI-safe.

### [ANAL-29] ItraqEightPlexQuantitationMethod::getReferenceChannel / updateMembers_
`silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.cpp

**Setting reference_channel=120 (an allowed param value) silently leaves the reference channel at a stale index**

- Expectation: The 'reference_channel' param is range-validated to [113,121] (setMinInt(113)/setMaxInt(121)), and the description says '113-121'. A caller setting reference_channel=120 expects either an error or a defined reference channel index back from getReferenceChannel().
- Actual: 120 passes range validation and reaches updateMembers_, where the branch 'else if (ref_ch == 120) { OPENMS_LOG_WARN << "Invalid channel selection." }' only logs a warning and does NOT assign reference_channel_. The member keeps whatever value it had before (0 after construction, or a previously-set index), so getReferenceChannel() returns a silently wrong/stale index instead of signaling failure.
- Fix: Make the invalid 120 selection a hard error: throw Exception::InvalidParameter in updateMembers_() for ref_ch==120 instead of only logging, or tighten validation so 120 cannot be set (e.g. setValidStrings/an explicit allowed-values check). Additive and ABI-safe (no signature change).

### [ANAL-33] AbsoluteQuantitation::calculateRatio
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h

**calculateRatio returns 0.0 instead of signaling a missing feature/value**

- Expectation: The doc advertises '@exception Exception::UnableToFit' and the name implies a real ratio. A caller expects either a valid ratio or a thrown error when the requested feature_name/IS is absent.
- Actual: When neither component has the requested metavalue/native_id, the function silently returns the initial 'ratio = 0.0' with only an OPENMS_LOG_DEBUG message; it never throws. A 0.0 ratio is indistinguishable from a legitimately measured zero and propagates into fitCalibration/applyCalibration as a real data point.
- Fix: Either throw Exception::UnableToFit (as documented) when the feature/IS value is missing, or return std::optional / a NaN sentinel that callers must check. Behavior change; gate behind a flag if existing pipelines rely on the 0.0 fallback. Doc-only fix (drop the bogus @exception) is the minimal source-compatible step.
- Verifier note: Claim is accurate. Minor refinement: the silent-0.0 return is not merely incidental — it is explicitly enshrined by an existing unit test (AbsoluteQuantitation_test.cpp:201 asserts the missing-feature case == 0.0 via TEST_REAL_SIMILAR). This both strengthens the finding (the silent-0.0 path is exercised and intended) and means the 'throw Exception::UnableToFit' fix is more invasive than the recomm

### [ANAL-42] NuXLFDR::splitIntoPeptidesAndXLs
`silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFDR.h

**splitIntoPeptidesAndXLs keeps only the single first hit overall, not the first hit of each class as documented**

- Expectation: The header states an input identification 'contributes to pep_pi or xl_pi (or both, if it had hits of both classes)' and keeps 'the first hit encountered for each class'. A caller therefore expects a PSM with a plain-peptide hit and a cross-link hit to appear in BOTH output lists.
- Actual: The loop guards every push with `if (pep_ph.empty() && xl_ph.empty())` for BOTH branches (NuXLFDR.cpp lines 47-54). Once ANY first hit is added (to either class), this condition is false for all remaining hits, so a second hit of the OTHER class is never captured. An identification can therefore never contribute to both lists; only its single best hit is retained. This silently drops cross-link evidence whenever a plain-peptide hit ranks above it (and vice versa).
- Fix: If both-class capture is intended (as the header claims), gate each branch on its own vector being empty (`if (pep_ph.empty())` / `if (xl_ph.empty())`). If single-best-hit is intended, fix the header doc to say so. Behavior change is the FDR-correct fix; doc fix is the conservative one. ABI: either is source-compatible.
- Verifier note: The claim is accurate. Precise statement: with report_top_hits >= 2 (use_all_hits), a PeptideIdentification carrying hits of both classes is split into only one output list (whichever class's hit appears first in the hit order), and only that single hit is retained — silently discarding the documented per-class first hit of the other class before FDR is computed. Fix per recommendation: gate each 

### [ANAL-43] NuXLFDR::calculatePeptideAndXLQValueAndFilterAtPSMLevel
`silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFDR.h

**Tie-breaker divides by score range and silently produces inf/NaN scores when all XL hits share one score**

- Expectation: A deterministic 'add a tiny tie-breaker to q-value scores' step should leave scores well-defined for any input, including the degenerate case where all XL hits have equal score.
- Actual: The tie-breaker computes `score_range = max_score - min_score` and then `(1.0 - (score - min_score) / score_range) * 1e-5` (NuXLFDR.cpp lines 176-190). When every XL hit shares the same svm_score/NuXL:score, score_range == 0 and the division yields inf/NaN, corrupting every XL hit's stored score with no error. The header even flags this as a @warning but the code neither guards nor throws.
- Fix: Guard `if (score_range <= 0) skip tie-break` (or add an epsilon). Additive and non-breaking. ABI: source-compatible.
- Verifier note: In the degenerate case (all XL hits share one score) score_range == 0 and the numerator (score - min_score) is also 0, so the division is 0/0 = NaN — the corrupted scores are NaN, not 'inf'. NaN is then added to every XL hit's score via p.setScore(...), with no guard or throw. The header @warning at NuXLFDR.h:131-133 documents the hazard but the code neither guards (e.g. `if (score_range <= 0) ski

### [ANAL-51] NuXLPresets::getPresets
`hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLPresets.h

**getPresets mutates the global ResidueDB singleton (adds methionine loss) for DEB/NM presets**

- Expectation: A function named getPresets that fills output StringList/string parameters reads configuration; it should not mutate global chemistry state shared by the whole process.
- Actual: When the preset name contains the substring 'DEB' or 'NM', getPresets reaches into the global ResidueDB singleton, casts away const on the Methionine residue, and permanently adds a CH4S1 loss formula to it. This silently changes the chemistry of residue 'M' for every other consumer of ResidueDB in the process, and is not reversible.
- Fix: Document the ResidueDB mutation prominently in the header doc (@note that calling this with a DEB/NM preset alters global Methionine chemistry). Ideally move the residue-loss registration out of the preset getter into an explicit, separately-named setup call (e.g. applyPresetChemistry(...)) so a pure 'get' has no global effect — additive and source-compatible if the old behavior is preserved behind a deprecated path.
- Verifier note: Claim is accurate as stated. Strengthening note: the mutation not only changes 'M' chemistry globally but is also cumulative — because addLossFormula does an unconditional push_back, invoking getPresets more than once with a DEB/NM preset in the same process appends the CH4S1 loss multiple times, compounding the corruption. It also directly contradicts the immutability invariant explicitly documen

### [ANAL-58] DeconvolvedSpectrum::toSpectrum
`surprising-crash` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h

**toSpectrum() dereferences peak_groups_[0] with no empty check, crashing on an empty DeconvolvedSpectrum**

- Expectation: Converting a (possibly empty) DeconvolvedSpectrum to an MSSpectrum should return an empty/header-only spectrum, not crash. Default-constructed DeconvolvedSpectrum objects exist and are common.
- Actual: Early in toSpectrum the code unconditionally reads peak_groups_[0].isPositive() to build the metadata string. If the spectrum has zero peak groups, this is out-of-bounds access (UB / crash) before any peaks are even iterated.
- Fix: Add an early `if (peak_groups_.empty()) return out_spec;` (or branch the chargemass metadata) before indexing peak_groups_[0]. Pure implementation fix, no ABI/signature change.
- Verifier note: toSpectrum does not throw a catchable exception; line 41's peak_groups_[0] uses the unchecked std::vector::operator[], so calling toSpectrum on an empty/default-constructed DeconvolvedSpectrum is out-of-bounds access (undefined behavior / crash, not a thrown OpenMS exception). The defect and recommended one-line guard fix are otherwise exactly as claimed; this is an undocumented non-empty precondi


## B5 (9)

### [FEAT-28] SeedListGenerator::generateSeedList
`silent-failure` · ABI: `none` · src/openms/source/FEATUREFINDER/SeedListGenerator.cpp

**generateSeedList dereferences a possibly-end precursor iterator and indexes precursors[0] without bounds checks**

- Expectation: A 'generate a seed list from MS2 precursors' helper should robustly skip MS2 spectra that lack a usable precursor or a preceding MS1 spectrum, not crash on them.
- Actual: For every MS2 spectrum it calls experiment.getPrecursorSpectrum(exp_it) and immediately dereferences prec_it->getRT(), and reads precursors[0].getMZ(). getPrecursorSpectrum can legitimately return experiment.end() (e.g. an MS2 at the start of the run, or no preceding MS1), and getPrecursors() can be empty for an MS2 with no recorded precursor. Both cases are undefined behavior / out-of-bounds access.
- Fix: Guard the loop body: skip the spectrum when precursors.empty() or prec_it == experiment.end(). Additive, no ABI change (implementation-only fix).
- Verifier note: Category "silent-failure" is partly imprecise: Path 1 (prec_it->getRT() on an end() iterator for an MS2 at run start / no preceding MS1) is more accurately undefined-behavior / likely hard crash than a silent failure, while Path 2 (precursors[0] on an empty vector) is the silent out-of-bounds read that can produce garbage seed positions without crashing. Both are real and reachable. The recommenda

### [FEAT-35] Biosaur2Algorithm::run
`hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/FEATUREFINDER/Biosaur2Algorithm.h

**run() destructively mutates the stored MS data set via setMSData(); getMSData() no longer returns what was set**

- Expectation: After setMSData(exp) and run(...), getMSData() returns the experiment that was set, and run() can be called again on the same stored data.
- Actual: run() erases all non-MS1 spectra from the internal ms_data_, and may centroid (PeakPickerHiRes) and apply TOF filtering in place. The stored experiment is consumed/altered; a second run() or a subsequent getMSData() sees the reduced/centroided data, not the original.
- Fix: Operate on a local copy of ms_data_ inside run() so the stored data is preserved and run() is idempotent (behavior change, but source-compatible). Minimally, strengthen the header note to say the stored data is irreversibly modified and run() is not re-entrant on the same instance.
- Verifier note: Accurate but understated: in addition to MS1 erasure, profile centroiding, and TOF filtering on ms_data_, the FAIMS path (Biosaur2Algorithm.cpp:284) std::move's ms_data_ into IMDataConverter::splitByFAIMSCV, leaving the stored experiment empty (moved-from). Thus after a FAIMS run() getMSData() returns an emptied experiment, and after any run() a second run() sees reduced/centroided (or empty) data

### [FEAT-38] ElutionPeakDetection::filterByPeakWidth
`hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/ElutionPeakDetection.h

**filterByPeakWidth unconditionally prints to std::cout and indexes the result vector without an empty check**

- Expectation: A library filter quietly writes the surviving traces to the output vector and returns; it does not print to stdout and does not crash on edge inputs.
- Actual: It always writes a 'pw low/pw high' line to std::cout, and then indexes filt_mtraces[0] and filt_mtraces[filt_mtraces.size()-1]. If the input is empty (or all traces are filtered out), filt_mtraces is empty and these accesses are out-of-bounds (UB/crash).
- Fix: Remove the std::cout diagnostic (or gate it behind a debug macro) and guard the indexing with `if (!filt_mtraces.empty())`. Both are additive and ABI-safe.

### [PROC-15] NLargest::filterSpectrum
`hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/FILTERING/NLargest.h

**Filtering reorders the spectrum to descending-intensity order, not just removing peaks**

- Expectation: A caller of an 'N largest peaks' filter expects the surviving peaks to remain in the spectrum's natural (m/z / position) order; only count changes.
- Actual: filterSpectrum calls spectrum.sortByIntensity(true) and then spectrum.select(indices) with indices 0..n-1. MSSpectrum::select preserves the order of the given indices (confirmed in MSSpectrum.cpp), so the spectrum is left sorted by DESCENDING intensity, not by m/z. The peaks are never re-sorted back by position.
- Fix: Re-sort the spectrum by position after select() (e.g. spectrum.sortByPosition()), or clearly document in the header that the spectrum is left in descending-intensity order. ABI-safe: it is an implementation-body change in an inline template. If callers depend on the current order, add a bool keep_order param with a default preserving today's behavior.

### [PROC-16] RankScaler::filterSpectrum
`hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/SCALING/RankScaler.h

**RankScaler leaves the spectrum sorted by intensity (not m/z) after scaling**

- Expectation: A 'scaler' rewrites intensities in place; a caller expects peak ordering (by m/z) to be unchanged.
- Actual: filterSpectrum calls spectrum.sortByIntensity() and never restores positional order, so on return the spectrum is sorted by intensity. Any downstream code that assumes m/z-sorted peaks (e.g. for alignment/scoring) will silently operate on a reordered spectrum.
- Fix: Call spectrum.sortByPosition() at the end, or document that the spectrum is returned in intensity-sorted order. ABI-safe inline-body change.

### [PROC-23] estimateNoiseFromRandomScans
`silent-failure` · ABI: `none` · src/openms/source/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.cpp

**Random-scan picker indexes exp[] directly with a random number, not via the filtered scan-index list, so it can sample wrong-MS-level / empty spectra (and overrun exp)**

- Expectation: A caller reading 'estimateNoiseFromRandomScans(exp, ms_level, ...)' expects it to draw n_scans random spectra OF THE REQUESTED ms_level (it built spec_indices precisely for that) and average their percentile intensity.
- Actual: It builds spec_indices of valid same-level scans, then computes `UInt scan = distribution * (spec_indices.size()-1)` and dereferences `exp[scan]` -- the random value indexes the FULL experiment, NOT spec_indices. So it samples arbitrary spectra (wrong MS level, empty, or chromatograms) and `scan` is bounded by spec_indices.size(), unrelated to exp.size(). If a sampled spectrum is empty, `tmp[idx]` reads tmp[0] of an empty vector (UB). The whole ms_level filtering is silently discarded.
- Fix: Index through the filter list: `UInt scan = spec_indices[(Size)(distribution(generator) * (spec_indices.size()-1))];` then `exp[scan]`. This is a pure source-compatible bug fix (no signature change). Also clamp idx for tiny spectra.
- Verifier note: The bug is real but one detail is wrong: exp[scan] does NOT overrun exp, because scan is bounded by spec_indices.size()-1 and spec_indices.size() <= exp.size(). The actual out-of-bounds/UB is tmp[idx] when a sampled spectrum is empty (empty-vector read). The core defect stands: exp[scan] must be exp[spec_indices[scan]]; the ms_level and non-empty filtering is silently bypassed, yielding wrong nois

### [PROC-31] FeatureOverlapFilter::mergeFAIMSFeatures
`hidden-side-effect` · ABI: `none` · src/openms/source/PROCESSING/FEATURE/FeatureOverlapFilter.cpp

**mergeFAIMSFeatures silently wipes protein IDs, peptide IDs, data processing and meta info via feature_map.clear()**

- Expectation: A function documented as 'Merge FAIMS features ... (modified in place)' rebuilds the feature list but otherwise preserves the FeatureMap's attached metadata (protein/peptide identifications, document/unique IDs, ranges, data processing, MetaInfo).
- Actual: Line 515 calls feature_map.clear() with the default argument. FeatureMap::clear(bool clear_meta_data = true) (FeatureMap.cpp:462) clears not only data_ but also clearMetaInfo(), clearRanges(), the DocumentIdentifier, clearUniqueId(), protein_identifications_, unassigned_peptide_identifications_, data_processing_ and id_data_. So all of this metadata is silently destroyed whenever any FAIMS features are present, while the same call is a no-op (returns early) when no FAIMS features exist - giving inconsistent, data-loss behavior.
- Fix: Call feature_map.clear(false) to keep attached metadata, or repopulate data in place via swap of the feature container only. This is a bug fix that is source- and ABI-compatible (no signature change). Also document that unassigned peptide IDs are not carried through the merge.

### [PROC-32] IDFilter::filterHitsByScore(std::vector<IdentificationType>&, double, IDScoreSwitcherAlgorithm::ScoreType)
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/ID/IDFilter.h

**Doc and warning claim hits are removed when score_type is missing, but hits are silently kept**

- Expectation: The Doxygen note states 'Removes a hit if the @p score_type is not found at all' and the emitted warning says 'No hit with the given score_type found. All hits removed.' A caller therefore expects that IDs lacking the requested score type have all their hits removed.
- Actual: When neither the main score nor any secondary score of an ID matches score_type, the else-branch finds an empty score_name and does nothing - keepMatchingItems is never called, so every hit of that ID is retained. The 'All hits removed' warning only fires if NOT A SINGLE id in the whole vector matched, and even then no hits are actually removed. The actual behavior is the opposite of both the note and the warning.
- Fix: Decide on one behavior and make code, note and warning agree: either clear hits of IDs whose score_type is absent (matching the documented contract) or change the note/warning to say such hits are left untouched. Body-only change; source- and ABI-compatible.

### [PROC-33] MorphologicalFilter::filterRange / applyErosion_ / applyDilation_
`hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h

**Public filter uses function-local static buffers, making filter()/filterExperiment() non-reentrant and not thread-safe**

- Expectation: A non-static member filter method taking explicit input/output ranges is expected to be reentrant and safe to call concurrently on independent MorphologicalFilter instances (a common pattern when parallelizing per-spectrum filtering).
- Actual: filterRange uses `static std::vector<...> buffer;` (line 171) and applyErosion_/applyDilation_ each use their own `static std::vector<ValueType> buffer;` (lines 330, 440). These statics are shared across all instances and threads, so concurrent calls corrupt each other's intermediate results. filter() and filterExperiment() route through these methods, so the whole public API is silently unsafe to run in parallel despite nothing in the signature suggesting shared state.
- Fix: Make the scratch buffers per-instance members (or stack-local) instead of static. Note this changes templated inline code that is recompiled by clients, so it is source-compatible but technically a header/inline-definition change; gate behind a minor version. Document thread-safety expectations in the class docs regardless.


## B6 (8)

### [MATH-7] Math::absdev
`misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/StatisticFunctions.h

**absdev ("absolute deviation") sums SIGNED deviations, not absolute ones - always ~0 for the default mean**

- Expectation: A function named 'absdev' / documented as 'absolute deviation' returns the mean absolute deviation, i.e. mean(|x_i - mean|), a positive spread measure (like MeanAbsoluteDeviation a few lines above which DOES use fabs).
- Actual: The loop does `sum_value += *iter - mean;` with NO fabs. When `mean` is the actual mean (the default path computes it internally), the signed deviations cancel and the result is ~0 (up to floating point error). It is mathematically a near-zero number, never the absolute deviation.
- Fix: This is a latent bug, not just a naming issue. Add fabs: `sum_value += std::fabs(*iter - mean);` so it actually computes mean absolute deviation. Since the current result is essentially meaningless (~0), fixing it is source/ABI-compatible at the signature level; only the (broken) numeric output changes. Alternatively deprecate in favor of the already-correct MeanAbsoluteDeviation.

### [MATH-14] Math::MultipleTesting::computeModelFDR
`silent-failure` · ABI: `none` · src/openms/include/OpenMS/MATH/STATISTICS/MultipleTesting.h

**computeModelFDR returns an all-NaN vector (no throw) when any input PEP is NaN**

- Expectation: A method named computeModelFDR taking posterior error probabilities would, on a malformed input (a single NaN PEP), either skip it or signal an error.
- Actual: If ANY element is NaN, the function returns immediately with a vector of all-NaN of the same length and no diagnostic - one bad value silently wipes every FDR estimate. A caller that does not check for NaN gets silently corrupted results for the whole set.
- Fix: Document the all-or-nothing NaN propagation in the brief, or filter NaNs like qValue() does instead of poisoning the whole output. Doc-only fix is ABI-safe; filtering changes output values (source-compatible).

### [ML-10] ClusteringGrid::getClusters
`surprising-throw` · ABI: `source-compatible` · src/openms/include/OpenMS/ML/CLUSTERING/ClusteringGrid.h

**getClusters dereferences map end() (UB / crash) for an empty cell, undocumented**

- Expectation: A const getter named getClusters returning a std::list<int> for a cell should, for a cell with no clusters, return an empty list (the natural reading of '@return list of cluster indices ... centred in this cell').
- Actual: Implementation is `return cells_.find(cell_index)->second;` with no existence check. For a cell that is not present in `cells_`, `find` returns `end()` and the `->second` dereference is undefined behavior (typically a crash). The only reason internal callers are safe is that every call site first guards with isNonEmptyCell(); a public caller has no way to know this from the signature.
- Fix: Guard the lookup: return an empty list when the cell is absent (additive, ABI-safe, no signature change). At minimum document the precondition that the cell must be non-empty.
- Verifier note: The behavior is UB/crash (dereferencing map end()), not strictly a thrown exception; the category label "surprising-throw" is the closest available bucket but the actual failure mode is undefined behavior rather than a documented throw. Severity high because a public const getter crashes for a reasonable input (an empty/absent cell) with no loud signal. The recommended fix (guard and return an emp

### [ML-21] ROCCurve::cutoffNeg
`misleading-name` · ABI: `source-compatible` · src/openms/source/ML/ROCCURVE/ROCCurve.cpp

**cutoffNeg() iterates only over positive samples but divides by the negative count**

- Expectation: cutoffNeg(fraction) -- the symmetric partner of cutoffPos -- should compute a cutoff based on the negative class (true negatives among negatives), as the name and 'trueNeg' variable suggest.
- Actual: The loop only enters its body when 'cit->second' is true (i.e. the sample is a POSITIVE), yet increments a counter named trueNeg and divides it by neg_ (the negative count). It therefore mixes positive samples with the negative population size; the returned score is computed from the wrong class. cutoffPos uses the identical 'if (cit->second)' guard with truePos/pos_, which is correct there, making cutoffNeg's copy-paste error clear.
- Fix: The negative-class loop should test '!cit->second'. Fix the predicate to '!cit->second' (and verify against the rocN convention). Source-compatible behavioral fix; document the threshold semantics. The class is already marked '[buggy and usage is discouraged]', which corroborates this.

### [COMP-1] BinnedSharedPeakCount::operator()(const BinnedSpectrum&, const BinnedSpectrum&)
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/COMPARISON/BinnedSharedPeakCount.h

**Documented @throw IncompatibleBinning never happens; mismatched binning is silently unchecked in release builds**

- Expectation: Per the header's @throw documentation, passing two spectra with different binning throws IncompatibleBinning, so a caller can rely on the call to validate compatibility and fail loudly on incompatible input.
- Actual: The implementation only checks binning via OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(...), ...). OPENMS_PRECONDITION is compiled out unless OPENMS_ASSERTIONS is defined (debug-only), so in a normal release build there is NO check at all. Even when active, it throws Exception::Precondition, not IncompatibleBinning. With incompatible bins the code proceeds to cwiseProduct on mismatched sparse vectors and returns a meaningless score.
- Fix: Either (a) restore the documented contract by throwing Exception::IncompatibleBinning (or a real exception) unconditionally when !isCompatible(...), independent of OPENMS_ASSERTIONS, or (b) correct the header to drop the @throw IncompatibleBinning line and state that compatibility is only checked in assertion-enabled builds. Option (a) is source-compatible (adds a throw on a currently-UB path); (b) is doc-only.
- Verifier note: Refinement only (not a contradiction): the claim says OPENMS_PRECONDITION is "debug-only." More precisely, on Linux/non-MSVC it is disabled even in Debug builds because config.h gates OPENMS_ASSERTIONS behind `#if (0)`; it is active only in MSVC Debug. Net effect for the documented contract is the same or worse: the @throw is never honored in any standard Linux build, release or debug. Recommended

### [COMP-3] BinnedSumAgreeingIntensities::operator()(const BinnedSpectrum&, const BinnedSpectrum&)
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/COMPARISON/BinnedSumAgreeingIntensities.h

**Documented @throw IncompatibleBinning never thrown; incompatible bins silently produce a garbage score in release builds**

- Expectation: Header says passing spectra with different binning throws IncompatibleBinning, so callers can pass arbitrary pairs and trust the call to validate.
- Actual: Only OPENMS_PRECONDITION guards compatibility, which is a no-op in release builds (no OPENMS_ASSERTIONS) and throws Exception::Precondition (not IncompatibleBinning) in debug. With mismatched bins the elementwise add/subtract on differently-sized sparse vectors yields an undefined/meaningless score with no error.
- Fix: Throw a real exception (Exception::IncompatibleBinning / IllegalArgument) unconditionally on !isCompatible, or fix the header to remove the false @throw guarantee. Adding an unconditional throw on a currently-UB path is source-compatible.

### [COMP-4] BinnedSpectralContrastAngle::operator()(const BinnedSpectrum&, const BinnedSpectrum&)
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/COMPARISON/BinnedSpectralContrastAngle.h

**Documented @throw IncompatibleBinning never thrown; division by sqrt of empty-spectrum norms can return NaN with no error**

- Expectation: Header promises IncompatibleBinning on mismatched binning, and a similarity functor would be expected to either validate input or return a defined score for degenerate input.
- Actual: Compatibility is only checked via the debug-only OPENMS_PRECONDITION (no-op in release; throws Exception::Precondition not IncompatibleBinning in debug). Additionally, for an empty spectrum (all-zero bins) sum1 or sum2 is 0, so score = numerator/sqrt(0) = 0/0 = NaN is returned silently instead of a defined value or an error.
- Fix: Throw unconditionally on incompatible binning, and guard the normalization (return 0.0 when either norm is 0) so the score never silently becomes NaN. Both are additive/source-compatible. Also correct the header @throw line.
- Verifier note: Minor refinement: the documented IncompatibleBinning is not merely "thrown as a different exception type" — that named exception class does not exist anywhere in OpenMS (it appears only in three @throw doc comments). So the contract is unfulfillable as written. The compatibility check is debug-only (OPENMS_PRECONDITION, no-op in release; throws Exception::Precondition, not IncompatibleBinning, in 

### [COMP-6] SpectraSTSimilarityScore::dot_bias
`surprising-default` · ABI: `none` · src/openms/include/OpenMS/COMPARISON/SpectraSTSimilarityScore.h

**dot_bias default -1 does NOT trigger recomputation; only 0 does**

- Expectation: The Doxygen says '@param[in] dot_product if -1 this value will be calculated as well.' and the default argument is -1, so a caller passing the default (or explicitly -1) expects the dot product to be computed internally.
- Actual: The implementation branches on `if (dot_product != 0)` and only recomputes via `(*this)(bin1, bin2)` when dot_product == 0. With the documented/default sentinel -1, it takes the `dot_product != 0` branch and divides by -1, returning a negative, meaningless dot_bias. The recompute path is only reached for dot_product == 0, the value the doc does not mention.
- Fix: Fix the implementation to honor the documented sentinel: branch on `if (dot_product >= 0)` (or `!= -1`) and recompute otherwise; or change the default and doc to 0. ABI-neutral (implementation/doc fix). If changing the default value, that is source-compatible only.


## B7 (4)

### [QC-3] FragmentMassError::compute
`inconsistent-convention` · ABI: `none` · src/openms/source/QC/FragmentMassError.cpp

**Two compute() overloads compute the average/variance differently, yielding different results**

- Expectation: Both compute() overloads carry identical doc ('Stores average FME over all spectra and its variance') and should produce the same statistics for the same data.
- Actual: The FeatureMap overload accumulates all PPM errors, computes result.average_ppm ONCE after the full first pass (cpp:238), then runs a SECOND pass to compute variance against that final average (cpp:241). The PeptideIdentificationList overload instead computes result.average_ppm and calls calculateVariance_ INSIDE the per-pep_id loop (cpp:303-305): the average is recomputed on every iteration and variance is accumulated against a moving/partial average. calculateVariance_ also divides each contribution by counter_ppm (cpp:166), which keeps changing across iterations. The two overloads therefore produce different variance (and a wrong one in the list overload), despite identical contracts.
- Fix: Restructure the list overload to mirror the FeatureMap overload: full accumulation pass, compute the final average once, then a second pass for variance. This is an implementation fix only; no signature/ABI change.

### [QC-8] PeptideMass::compute
`misleading-name` · ABI: `none` · src/openms/include/OpenMS/QC/PeptideMass.h

**PeptideMass annotates the OBSERVED mass from precursor m/z, not the theoretical mass the name/docs promise**

- Expectation: The class is documented as 'QC metric calculating theoretical mass of a peptide sequence' and the method @brief says 'Sets the mass metavalue to all PeptideHits by computing the theoretical mass'. A caller expects the theoretical neutral mass of the peptide SEQUENCE (e.g. AASequence::getMonoWeight()), independent of the measured precursor.
- Actual: The implementation computes the observed/experimental neutral mass from the precursor m/z: hit.setMetaValue("mass", (pi.getMZ() - Constants::PROTON_MASS_U) * hit.getCharge()). It never touches the peptide sequence at all, so the 'mass' value depends on measured m/z and charge, not on the theoretical sequence mass. For a mis-assigned/charge-wrong PSM the two differ; a caller reading 'mass' as theoretical mass gets a silent bug.
- Fix: Either rename/redocument the metavalue and method to reflect that it is the observed neutral mass derived from precursor m/z (source-compatible doc fix; metavalue rename is breaking), or change the implementation to compute the sequence theoretical mass via AASequence (behavior change). At minimum fix the doc to stop calling it 'theoretical'.

### [QC-11] QCBase::names_of_requires / QCBase::Requires
`silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/QC/QCBase.h

**names_of_requires[] is missing the entry for Requires::ID, causing an out-of-bounds read in isRunnable()**

- Expectation: names_of_requires[] is documented as 'strings corresponding to enum Requires' and is indexed by the enum value (QCBase::names_of_requires[i] for i in [0, SIZE_OF_REQUIRES)). A caller/maintainer expects one string per enum member up to (but excluding) SIZE_OF_REQUIRES.
- Actual: The enum has 7 real members (NOTHING=0 .. ID=6), but names_of_requires only contains 6 strings: {"fail","raw.mzML","postFDR.featureXML","preFDR.featureXML","contaminants.fasta","trafoAlign.trafoXML"}. QCBase::isRunnable() loops 'for (i=0; i < (UInt64)SIZE_OF_REQUIRES; ++i)' and reads names_of_requires[i]; when a metric requires ID (i=6) the warning path reads names_of_requires[6] — past the end of the 6-element array (undefined behavior). The ID requirement also has no human-readable name.
- Fix: Add the missing 'ID' string (and an entry for the NOTHING/SIZE bookkeeping if intended) so the array length matches the enum; keep it append-only to preserve indices. Source-compatible additive fix.
- Verifier note: Confirmed as stated. Minor clarification: the OOB read occurs only on the diagnostic warning path inside isRunnable() (the function still correctly returns false), and it is reachable via IdentificationSummary (the only metric requiring Requires::ID, used by MzQCFile) when the ID input is missing. It reads a std::string at past-the-end memory and streams it to the log (UB: garbage output, potentia

### [SYST-6] File::getPathLocations / File::executableExtensions_
`surprising-crash` · ABI: `source-compatible` · src/openms/include/OpenMS/SYSTEM/File.h

**Default argument std::getenv("PATH") constructs a std::string from a possibly-null pointer (UB/crash) when PATH is unset**

- Expectation: Calling getPathLocations() with no argument reads the PATH environment variable and returns its components; if PATH is unset it should yield an empty list, not crash.
- Actual: The default argument is 'std::getenv("PATH")', whose result is implicitly converted to std::string. If the environment variable is not set, getenv returns nullptr and constructing std::string(nullptr) is undefined behavior (typically a crash). The same pattern is used for the Windows %PATHEXT% default in executableExtensions_(const std::string& ext = std::getenv("PATHEXT")).
- Fix: Change the default to a const char* / nullptr-safe wrapper (e.g. default to "" and resolve getenv internally with a null check), so an unset PATH/PATHEXT yields an empty result instead of UB. The public getPathLocations is the ABI-relevant one; a null-guard inside the body plus a safe default is source-compatible.
- Verifier note: The category should be 'surprising-crash'/UB rather than 'surprising-throw': constructing std::string from the null pointer returned by std::getenv when PATH/PATHEXT is unset is undefined behavior (a crash), not a thrown C++ exception. Everything else in the claim holds. Reachable via File::findExecutable -> getPathLocations()/executableExtensions_() (File.cpp:863,869), called from TOPPBase.cpp:14
