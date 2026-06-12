# OpenMS POLS Audit — Batch 2: METADATA · CHEMISTRY

**Confirmed findings:** 124 (METADATA 55, CHEMISTRY 69).  
**Severity:** 4 high · 69 medium · 51 low.  
**Method:** 18 header-cluster finders → adversarial per-finding verification against source → (run in two passes; a transient network burst killed 9 finders in pass 1, re-run with retry in pass 2).

> Post-verification severity/category/ABI shown. Recommendations favor ABI-safe fixes.


---

# METADATA

**Counts:** 3 high · 33 medium · 19 low

### [META-9] SpectrumMetaDataLookup::addMissingSpectrumReferences — Documented [in,out] 'proteins' is taken by value, so spectra_data updates are silently discarded
`high` · `param-order-or-bool` · ABI: `breaking` · src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h

```cpp
static bool addMissingSpectrumReferences(PeptideIdentificationList& peptides, const std::string& filename, bool stop_on_error = false, bool override_spectra_data = false, bool override_spectra_references = false, std::vector<ProteinIdentification> proteins = std::vector<ProteinIdentification>())
```
- **Expectation:** The Doxygen marks proteins as '@param[in,out] proteins Protein IDs corresponding to the Peptide IDs' and says with override_spectra_data the 'ProteinIdentifications should be updated with new spectra_data values'. A caller passing its protein vector expects the prot.setMetaValue("spectra_data", ...) writes to be visible after the call.
- **Actual:** The parameter is declared by value (std::vector<ProteinIdentification> proteins), not by reference. In SpectrumMetaDataLookup.cpp the function loops 'for (auto& prot : proteins) prot.setMetaValue("spectra_data", spectra_data);' but mutates only the local copy, which is destroyed on return. The caller's protein IDs are never updated.
- **Evidence:** Header: 'std::vector<ProteinIdentification> proteins = std::vector<ProteinIdentification>()' with doc '@param[in,out] proteins'. Cpp lines 324-332 set spectra_data on the by-value copy. IDFileConverter.cpp:420-426 passes protein_identifications with override_spectra_data=true, expecting the update to take effect.
- **Fix:** Change the parameter to std::vector<ProteinIdentification>& (true in/out), or if ABI must be preserved add an overload taking a reference and deprecate the by-value one. At minimum, fix the doc to stop claiming [in,out]. The current state is a silent bug at real call sites (IDFileConverter).
- **Verified:** Independently verified in the actual code. Header line 321 declares the parameter by value: `std::vector<ProteinIdentification> proteins = std::vector<ProteinIdentification>()`. The Doxygen directly above it (lines 308, 310) explicitly promises in/out semantics: `@param[in,out] proteins Protein IDs...` and `override_spectra_data ... ProteinIdentifications should be updated with new "spectra_data" values`. In SpectrumMetaDataLookup.cpp lines 324-3

### [META-12] SpectrumMetaDataLookup::getSpectrumMetaData — Silently returns leaving 'meta' untouched when no reference format matches, unlike findByReference which throws ParseError
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h

```cpp
void getSpectrumMetaData(const std::string& spectrum_ref, SpectrumMetaData& meta, MetaDataFlags flags = MDF_ALL) const
```
- **Expectation:** Given the class's other parse path SpectrumLookup::findByReference throws Exception::ParseError when 'no reference format matched', and this function's doc only mentions '@throw ElementNotFound if a spectrum look-up was necessary', a caller passing an unparseable reference expects either an exception or a documented signal of failure.
- **Actual:** If none of reference_formats matches the spectrum_ref, the for-loop completes without entering the body and the function returns normally, leaving the output 'meta' at its default-constructed values (rt/precursor_rt/precursor_mz = NaN, charge=0, ms_level=0, scan_number=-1). No exception, no return value to inspect.
- **Evidence:** Cpp lines 71-153: the work and the only 'return' (line 150) are inside 'if (found)'. There is no else/fallthrough error path, in contrast to SpectrumLookup::findByReference (SpectrumLookup.cpp:199-201) which throws Exception::ParseError.
- **Fix:** Throw Exception::ParseError on no-match (mirroring findByReference) or document the silent no-op and that callers must check meta for sentinel/NaN values. Behavioral fix; flag abi as source-compatible since it adds a throw.
- **Verifier correction:** Field is `precursor_charge` (not `charge`); it defaults to 0. All other claimed sentinels are exact. Severity kept high: on an unparseable reference the function silently returns leaving `meta` fully at sentinels (NaN RT/mz, ms_level=0 which is an invalid MS level, scan_number=-1), with no exception and no return value — a caller that does not defensively test for NaN/sentinels will silently propagate invalid data, while the documented-as-equivalent findByReference would throw. Recommendation stands: either throw Exception::ParseError on no-format-match to mirror findByReference, or document the silent no-op and require callers to check meta for sentinel/NaN values. ABI: adding a throw leaves the signature, vtable, and exported symbol unchanged, so it is source-compatible (recompiled callers are unaffected at the API level; only runtime behavior changes).
- **Verified:** Verified against the actual source. In SpectrumMetaDataLookup.cpp:67-153, getSpectrumMetaData(const std::string&, ...) iterates reference_formats; ALL work and the ONLY return (line 150) are inside `if (found)` (line 76). If no regex matches, the loop falls through and the function returns normally, leaving `meta` at its default-constructed sentinels — confirmed via the SpectrumMetaData ctor (header lines 148-153): rt/precursor_rt/precursor_mz=Na

### [META-31] ExperimentalDesign::getSample — getSample dereferences find_if result without checking end() (UB / crash on missing combination)
`high` · `silent-failure` · ABI: `none` · src/openms/source/METADATA/ExperimentalDesign.cpp

```cpp
unsigned getSample(unsigned fraction_group, unsigned label = 1)
```
- **Expectation:** A lookup that takes a fraction_group and label and 'returns the sample index' would either return the index, throw a clear ElementNotFound, or return a sentinel when no matching row exists.
- **Actual:** The body is `return std::find_if(msfile_section_.begin(), msfile_section_.end(), ...)->sample;` with no comparison to end(). If no entry matches (fraction_group, label), the end iterator is dereferenced, which is undefined behavior (out-of-bounds read or crash), not a diagnosable error.
- **Evidence:** ExperimentalDesign.cpp:599-606: `return std::find_if(msfile_section_.begin(), msfile_section_.end(), [&fraction_group, &label](const MSFileSectionEntry& r){ return r.fraction_group == fraction_group && r.label == label; })->sample;`
- **Fix:** Capture the iterator, compare to end(), and throw Exception::ElementNotFound (mirroring the rest of the class) when not found. Additive/behavioral fix; does not change the signature so ABI is preserved.
- **Verified:** Verified verbatim at src/openms/source/METADATA/ExperimentalDesign.cpp:599-606: getSample(unsigned fraction_group, unsigned label) does `return std::find_if(msfile_section_.begin(), msfile_section_.end(), [...](const MSFileSectionEntry& r){ return r.fraction_group == fraction_group && r.label == label; })->sample;` with no comparison against end(). If no row matches the (fraction_group, label) pair, find_if returns end() and `->sample` dereferenc

### [META-1] SourceFile::getFileSize / SourceFile::setFileSize — File size accessors use float while the backing member is double (silent precision narrowing)
`medium` · `return-value` · ABI: `breaking` · src/openms/include/OpenMS/METADATA/SourceFile.h

```cpp
float getFileSize() const; void setFileSize(float file_size);
```
- **Expectation:** A getter/setter pair for a stored value uses the same type as the member, so round-tripping set then get is lossless and matches what other numeric accessors in the module use (Sample::getMass etc. all use double).
- **Actual:** The member is `double file_size_;` (SourceFile.h:99) but the public API is `float getFileSize()` / `void setFileSize(float)`. setFileSize even casts: `file_size_ = static_cast<double>(file_size);` (SourceFile.cpp:77). A caller passing a double file size silently narrows to float, and getFileSize() truncates the stored double to float on return. The widened member buys nothing because the public surface is float on both ends.
- **Evidence:** SourceFile.h:70-72: `float getFileSize() const; ... void setFileSize(float file_size);` vs member `double file_size_;` (line 99). SourceFile.cpp:75-78 casts float arg up to double.
- **Fix:** Add `double getFileSizeInBytes()`/double-typed overloads or change the accessors to `double` to match the member (and the rest of the METADATA numeric accessors). If ABI must be preserved, add a new double-typed overload and deprecate the float one; at minimum document that values are stored/compared as double but truncated through the float API.
- **Verifier correction:** The float/double mismatch is real and the double member is genuinely pointless, but the claimed impact is overstated: the value is "file size in MB" and float comfortably represents realistic file sizes, so there is no practical loss of scientific accuracy or data — the surprise is API inconsistency and misleading precision (the static_cast suggests double is preserved when it never is), not silent result corruption. Severity is medium, not high.
- **Verified:** Evidence verified exactly. SourceFile.h:99 declares `double file_size_;` while the public API is `float getFileSize()` (h:70) and `void setFileSize(float)` (h:72). SourceFile.cpp:77 confirms the tell-tale `file_size_ = static_cast<double>(file_size);` — a no-op widening of an already-narrowed float, so the double member and the cast buy nothing: every value entering or leaving the object goes through float. operator== (cpp:37) and the hash helper

### [META-3] ContactPerson::setName — setName silently parses, reorders, and can drop name tokens instead of storing the name verbatim
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/METADATA/ContactPerson.h

```cpp
void setName(const std::string & name);
```
- **Expectation:** Among `setFirstName`/`setLastName`, a `setName(full)` convenience setter would store the full name (or split it predictably). The header only says 'gets split into first and last name internally'.
- **Actual:** setName applies surprising heuristics: if the string contains a comma it is treated as 'Last, First' and the order is swapped (first_name_ = tmp[1], last_name_ = tmp[0]); otherwise it splits on the FIRST space only via the two-element pattern and silently discards everything after the second token. So setName("Jane Q Public") yields first="Jane", last="Q" and drops "Public"; setName("Doe, John") swaps to first="John", last="Doe". A caller round-tripping getFirstName()+getLastName() will not recover the input.
- **Evidence:** ContactPerson.cpp:53-73: comma branch does `first_name_=trim(tmp[1]); last_name_=trim(tmp[0]);` (order swap); space branch does `first_name_=tmp[0]; last_name_=tmp[1];` using only the first two tokens (tmp[2+] ignored).
- **Fix:** Document the comma=Last,First convention and the lossy multi-token behavior explicitly in the header, or rename to `setNameFromCombinedString`/add a clearly-documented parsing contract. Ideally store the trailing tokens into last_name_ instead of discarding them. Doc/behavior fix is ABI-safe.
- **Verified:** Independently confirmed in ContactPerson.cpp:53-73. setName does apply undocumented heuristics: (1) comma branch swaps order, first_name_=trim(tmp[1]), last_name_=trim(tmp[0]) implementing an unstated "Last, First" convention; (2) space branch sets first_name_=tmp[0], last_name_=tmp[1] using ONLY the first two tokens. StringUtils::split (StringUtils.h:756-822) returns ALL tokens, so tmp[2+] are genuinely discarded. The named examples are correct:

### [META-4] Gradient::clearEluents / Gradient::clearTimepoints — clearEluents/clearTimepoints leave the percentages matrix desynchronized, breaking the class invariant
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/METADATA/Gradient.h

```cpp
void clearEluents(); void clearTimepoints();
```
- **Expectation:** addEluent/addTimepoint keep the eluents/timepoints vectors and the percentages matrix in lockstep (addEluent pushes a new percentages row, addTimepoint pushes a new column), so the matching clear* methods would also keep them consistent.
- **Actual:** clearEluents() only does `eluents_.clear()` and clearTimepoints() only does `times_.clear()`; neither touches `percentages_`. After clearEluents(), percentages_ still holds rows for the now-removed eluents, so getPercentages() reports rows with no matching eluent, and isValid() iterates `percentages_[i][j]` for i up to the OLD eluent count via stale data. A subsequent addTimepoint() pushes a 0 onto only `eluents_.size()` (=0) rows, permanently corrupting matrix shape. getPercentage()/setPercentage() then index a matrix whose dimensions no longer match the lookup vectors.
- **Evidence:** Gradient.cpp:46-49 clearEluents clears only eluents_; lines 71-74 clearTimepoints clears only times_. Contrast addEluent (lines 41-43) which appends `percentages_.emplace_back(times_.size(),0)` and addTimepoint (lines 62-68) which appends a column to each row. clearPercentages() (180-185) is the only method that reshapes percentages_.
- **Fix:** Have clearEluents()/clearTimepoints() also reset/reshape percentages_ to stay consistent (or call clearPercentages() afterwards). This is an implementation fix with no signature change, hence ABI-safe.
- **Verified:** Verified against src/openms/source/METADATA/Gradient.cpp. clearEluents() (lines 46-49) does only eluents_.clear() and clearTimepoints() (71-74) does only times_.clear(); neither touches percentages_. By contrast addEluent (41-43) appends percentages_.emplace_back(times_.size(),0) and addTimepoint (62-68) pushes a column onto each of eluents_.size() rows, so the add* methods deliberately keep eluents_/times_/percentages_ in lockstep (invariant doc

### [META-6] IdentificationData::registerProcessingSoftware — registerProcessingSoftware silently discards assigned_scores of a same-name/version software (no merge), unlike sibling register* functions
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/source/METADATA/ID/IdentificationData.cpp

```cpp
ProcessingSoftwareRef registerProcessingSoftware(const ProcessingSoftware& software)
```
- **Expectation:** Per the class contract documented in IdentificationData.h ('If items with an existing key are registered subsequently, attempts are made to merge new information ... into the existing entry'), registering a ProcessingSoftware whose name/version already exists but with a different/longer assigned_scores list should at least merge the new score assignments into the stored entry, or signal that the data was dropped.
- **Actual:** The body is just `return processing_softwares_.insert(software).first;`. Because ProcessingSoftwares is `std::set<ProcessingSoftware>` ordered by `Software::operator<` (which compares only name_ and version_ - confirmed in Software.cpp: `tie(name_, version_) < tie(...)`), the 'key' ignores assigned_scores. A second registration with the same name/version returns the pre-existing iterator and the new assigned_scores are silently thrown away. No merge() call is made, in contrast to registerInputFile/registerObservation which call existing.merge(file/obs).
- **Evidence:** registerProcessingSoftware: `return processing_softwares_.insert(software).first;` (IdentificationData.cpp:242). ProcessingSoftware.h:40 `// ordering is done using "operator<" inherited from "Software"`. Software.cpp:37-40 `return tie(name_, version_) < tie(rhs.name_, rhs.version_);`. Contrast registerInputFile (IdentificationData.cpp:213-220) which calls `existing.merge(file)`.
- **Fix:** Mirror the sibling register* functions: on `!result.second`, modify the stored element to merge in the new assigned_scores (e.g. append score refs not already present). This is an additive, source/ABI-compatible behavior fix. At minimum, document in ProcessingSoftware.h / the registerProcessingSoftware doxygen that assigned_scores is NOT part of the uniqueness key and will be ignored on re-registration.
- **Verifier correction:** Claim is accurate as stated; only the severity is adjusted from high to medium. The silently discarded data is the `assigned_scores` list (a vector of ScoreTypeRef = which score types the software assigns), not the actual numeric scores or identification results. On re-registration of a same name/version ProcessingSoftware, the pre-existing entry is kept and the new assigned_scores are dropped with no merge, no modify, and no warning — reachable both directly and via the public IdentificationData::merge() (IdentificationData.cpp:1064). Notably, ProcessingSoftware lacks any merge() method (unlike its siblings), so the documented contract in IdentificationData.h:61-62 ("attempts are made to merge new information (e.g. additional scores)") is not honored for this type. This is genuine silent data loss but limited to auxiliary score-type associations and recoverable, hence medium. Fix is additive/source/ABI-compatible: add the `!result.second` + `.modify(...)` pattern (and a ProcessingSoftware::merge that unions assigned_scores), or at minimum document that assigned_scores is not part of the uniqueness key.
- **Verified:** Independently verified every piece of evidence and it all holds. IdentificationData.cpp:242 is a bare `return processing_softwares_.insert(software).first;` with no `result.second` check and no merge/modify call. Software::operator< (Software.cpp:37-40) orders only by `tie(name_, version_)`, so assigned_scores is not part of the set key — a re-registration with same name/version returns the existing iterator and the incoming assigned_scores are s

### [META-8] IdentificationDataInternal::ScoredProcessingResult::getScore — getScore returns NaN with a success flag; ignoring the bool yields a NaN that silently poisons arithmetic
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/METADATA/ID/ScoredProcessingResult.h

```cpp
std::pair<double, bool> getScore(ScoreTypeRef score_ref) const
```
- **Expectation:** A function named getScore returning a pair would commonly be read as 'returns the score'. A caller who writes `double s = result.getScore(ref).first;` expects a usable number.
- **Actual:** On a missing score it returns `{quiet_NaN(), false}`. The first element is the documented sentinel NaN; correctness depends entirely on the caller checking `.second`. The header documents this ('@return A pair: score (or NaN), success indicator'), so it is a documented contract, but the pair<double,bool> ordering (value first, success second) makes it easy to consume `.first` while ignoring `.second`, silently propagating NaN.
- **Evidence:** getScore -> getScoreAndStep; on not-found `return std::make_tuple(std::numeric_limits<double>::quiet_NaN(), std::nullopt, false);` (ScoredProcessingResult.h:150-151). getScore(ScoreTypeRef, optional) directly returns `make_pair(quiet_NaN(), false)` (ScoredProcessingResult.h:127).
- **Fix:** This is documented, so keep it; consider additionally offering `std::optional<double> findScore(ScoreTypeRef) const` so callers cannot accidentally consume the NaN. Additive, ABI-safe. At minimum keep the '(or NaN)' note prominent.
- **Verified:** Code matches the claim exactly. ScoredProcessingResult::getScore(ScoreTypeRef) (h:101-106) delegates to getScoreAndStep (h:137-152), which on not-found returns make_tuple(quiet_NaN(), nullopt, false) at h:150-151; the getScore(ScoreTypeRef, optional) overload returns make_pair(quiet_NaN(), false) at h:127. Doxygen documents '@return A pair: score (or NaN), success indicator' (h:99,111), a documented contract. The pair<double,bool> ordering (value

### [META-10] SpectrumMetaDataLookup::addMissingIMToPeptideIDs — Documented bool-returning lookup actually throws ElementNotFound on an unmatched native ID, unlike its siblings
`medium` · `contract-violation/undocumented-throw` · ABI: `none` · src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h

```cpp
static bool addMissingIMToPeptideIDs(PeptideIdentificationList& peptides, const MSExperiment& exp)
```
- **Expectation:** The doc states '@return True if all missing IM information was successfully added ... false otherwise.' A caller reading this contract (and matching the sibling functions addMissingRTsToPeptideIDs / addMissingFAIMSToPeptideIDs, which catch ElementNotFound and return false) expects a non-throwing bool result.
- **Actual:** The implementation calls 'Size index = lookup.findByNativeID(native_id);' with no try/catch (SpectrumMetaDataLookup.cpp:211). If a peptide's spectrum_reference is not found among the spectra, findByNativeID throws Exception::ElementNotFound, which escapes the function instead of yielding 'false'.
- **Evidence:** Cpp lines 208-223: 'std::string native_id = pep.getSpectrumReference(); Size index = lookup.findByNativeID(native_id);' — no exception handling, contrasted with addMissingRTsToPeptideIDs (lines 176-188) and addMissingFAIMSToPeptideIDs (lines 282-296) which both wrap the lookup in try/catch.
- **Fix:** Wrap the findByNativeID call in try/catch like the sibling functions, setting all_ids_have_im=false on ElementNotFound; or document that this function may throw. Additive doc/impl fix, no signature change.
- **Verifier correction:** The category 'silent-failure' is inaccurate: the function does not fail silently — it throws Exception::ElementNotFound loudly. The real surprise is that a function whose doc (header line 279) promises a bool 'false otherwise' result, and whose two sibling functions (addMissingRTsToPeptideIDs, addMissingFAIMSToPeptideIDs) catch ElementNotFound, instead lets the exception propagate on an unmatched native ID. Recommended fix is additive (wrap findByNativeID in try/catch setting all_ids_have_im=false, or add an @throw clause to the doc); no signature/ABI change.
- **Verified:** Verified against the actual source. In SpectrumMetaDataLookup.cpp, addMissingIMToPeptideIDs (lines 195-225) calls `Size index = lookup.findByNativeID(native_id);` at line 211 with NO try/catch. findByNativeID (SpectrumLookup.cpp:65-75) throws Exception::ElementNotFound when the native ID is absent from its map (and its header explicitly documents @throw). The two sibling functions follow the identical lookup pattern but DO guard it: addMissingRTs

### [META-11] SpectrumMetaDataLookup::addMissingIMToPeptideIDs — 'addMissing...' unconditionally overwrites IM on all IDs, including those that already have it
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h

```cpp
static bool addMissingIMToPeptideIDs(PeptideIdentificationList& peptides, const MSExperiment& exp)
```
- **Expectation:** The name 'addMissingIMToPeptideIDs' and the doc 'adds missing ion mobility (IM) information' imply the function only fills in IDs lacking an IM value, leaving already-annotated IDs untouched.
- **Actual:** The loop processes every peptide ID unconditionally and calls 'pep.setMetaValue(Constants::UserParam::IM, spec.getDriftTime())' whenever the spectrum has IM, with no check for a pre-existing IM value. Already-set IM values are silently overwritten (unlike addMissingFAIMSToPeptideIDs, which skips IDs where the meta value already exists).
- **Evidence:** Cpp lines 208-223 have no 'metaValueExists' guard, whereas addMissingFAIMSToPeptideIDs lines 274-279 do: 'if (pep.metaValueExists(Constants::UserParam::FAIMS_CV)) { ++annotated_count; continue; }'.
- **Fix:** Either skip IDs that already carry a Constants::UserParam::IM meta value (to match the 'missing' semantics and the FAIMS sibling), or rename/redocument to make the overwrite behavior explicit. Impl-only fix; no ABI change.
- **Verified:** Verified directly against the source. In SpectrumMetaDataLookup.cpp the loop (lines 208-223) iterates every PeptideIdentification and calls pep.setMetaValue(Constants::UserParam::IM, spec.getDriftTime()) whenever determineIMFormat(spec)==IM_SPECTRUM, with NO metaValueExists(Constants::UserParam::IM) guard. So an ID that already carries an IM meta value is silently overwritten with the spectrum-level drift time. The sibling addMissingFAIMSToPeptid

### [META-13] USI::toString — toString() returns an empty string for an invalid USI instead of a representation of its fields
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/METADATA/USI.h

```cpp
std::string toString() const
```
- **Expectation:** A toString() on an object is expected to render the object's contents; for a partially-filled USI a caller would expect to see whatever fields are set (the doc says '@return Complete USI string'), or at least a non-silent indication of the problem.
- **Actual:** If isValid() is false (any of collection/ms_run/index empty), toString() returns "" with no other signal. The streaming operator<< (which just calls toString()) therefore prints nothing for an invalid USI, which can silently emit empty fields into output files.
- **Evidence:** USI.cpp lines 168-185: 'if (!isValid()) { return ""; }' before building the string. operator<< (lines 387-391) forwards to toString().
- **Fix:** Document the empty-string-on-invalid contract prominently on toString() (and operator<<), or render the partial fields. Doc-only fix preferred for ABI stability.
- **Verifier correction:** toString() returns "" (and operator<< prints nothing) when isValid() is false — i.e. when collection, ms_run, or index is empty — and this empty-on-invalid contract is undocumented; the header doc only states "@return Complete USI string". A real misuse path exists in TextExporter.cpp:590, where buildUSI(...).toString() is written directly into a TSV column and buildUSI explicitly returns an invalid USI() when no spectrum reference/MS-run mapping is available, so an invalid USI silently emits an empty column. Severity is medium (silent empty field, recoverable) rather than high. Fix: document the empty-string-on-invalid behavior on toString()/operator<< (or render the partial fields / signal the problem); a doc-only fix is ABI-neutral.
- **Verified:** Verified against source. USI.cpp:168-173 contains exactly `std::string USI::toString() const { if (!isValid()) { return ""; } ... }`, and isValid() (USI.cpp:89) returns `!collection_.empty() && !ms_run_.empty() && !index_.empty()`, so any missing required field yields an empty string. operator<< (USI.cpp:387-391) just forwards to toString(), so it prints nothing for an invalid USI. The header doc (USI.h:206-211) documents only `@return Complete U

### [META-14] SpectrumNativeIDParser::extractScanNumber — Same function name returns a real scan number, an index+1, or a synthetic cycle*1000+experiment value depending on accession
`medium` · `contract-violation` · ABI: `source-compatible` · src/openms/include/OpenMS/METADATA/SpectrumNativeIDParser.h

```cpp
static Int extractScanNumber(const std::string& native_id, const std::string& native_id_type_accession)
```
- **Expectation:** A function named extractScanNumber returning Int is expected to return THE scan number found in the native ID. A caller passing different accessions expects a consistent 'scan number' meaning.
- **Actual:** For WIFF (MS:1000770) it returns 'cycle * 1000 + experiment' — a synthetic composite, not any number literally in the ID; for index-based IDs (MS:1000774) it returns index + 1; for others it returns the raw scan integer. The same return slot thus carries three different semantics, and the WIFF path can even throw Exception::InvalidValue when experiment >= 1000.
- **Evidence:** Cpp lines 187-195 (WIFF: 'cycle * 1000 + experiment', throws InvalidValue if experiment >= 1000) and lines 164-167 ('index + 1 ... especially for pepXML'). The header @note documents both, but the unified name/return type still hides the unit shift.
- **Fix:** The header notes already disclose this; keep them and consider surfacing the WIFF composite via a clearly-named helper. The InvalidValue throw from a function otherwise documented to '@return ... -1 on failure' is the more actionable surprise — catch it and return -1, or document the throw. Source-compatible doc/impl fix.
- **Verifier correction:** The "three semantics in one Int" overload behavior is real but is explicitly documented at the declaration (header @note lines 90-91 and the format table lines 30-43) and reflects an accepted pepXML domain convention, so it is a low-severity disclosed quirk rather than a hidden trap. The genuine, actionable surprise is the throw of Exception::InvalidValue when experiment >= 1000 on the WIFF path (cpp lines 192-195): it bypasses the surrounding catch(ConversionError&), propagates through the SpectrumLookup wrapper, and contradicts the documented "@return ... -1 on failure" contract (header line 88, SpectrumLookup.h line 213) — the function carries no @throw clause and several callers only guard for the -1 sentinel, so they would crash on an uncaught exception. Fix is source-compatible: either catch and return -1, or add a @throw note. Re-graded category from unit-or-index to contract-violation and severity to medium because the failure is loud (an exception, recoverable) on a rare edge case rather than a silent wrong result.
- **Verified:** Evidence verified in src/openms/source/METADATA/SpectrumNativeIDParser.cpp. WIFF (MS:1000770) returns cycle*1000+experiment (lines 187-190); index-based MS:1000774 returns index+1 (lines 164-166); others return the raw scan integer (lines 168-171). So the single Int return slot does carry three semantics. HOWEVER, the headline "unit-or-index" surprise is largely defused by exclusion (b): both special cases are explicitly documented at the point o

### [META-16] MetaInfo::getValue — getValue returns a const reference to the caller-supplied default; dangles for temporaries
`medium` · `ownership-lifetime` · ABI: `source-compatible` · src/openms/include/OpenMS/METADATA/MetaInfo.h

```cpp
const DataValue& getValue(const std::string& name, const DataValue& default_value = DataValue::EMPTY) const
```
- **Expectation:** A getter with a default value can be safely called with an inline temporary, e.g. `const auto& v = mi.getValue("x", DataValue(0));` and the returned reference stays valid as long as it is read in the same statement / bound to a const ref.
- **Actual:** On miss the function returns the `default_value` parameter by reference (MetaInfo.cpp:87 `return default_value;`). If the caller passed a temporary, the returned reference dangles the moment the full-expression ends. The sibling MetaInfoInterface::getMetaValue(name, default) was deliberately changed to return BY VALUE for exactly this reason (header comment line 71: "Note: return needs to be by value to prevent life-time issues at caller site (e.g. if he passes a temporary to default-value)"), but MetaInfo::getValue still returns by reference, so the underlying class retains the trap.
- **Evidence:** MetaInfo.h:90 `const DataValue& getValue(const std::string& name, const DataValue& default_value = DataValue::EMPTY) const;`  MetaInfo.cpp:80-88 returns `default_value` by reference on miss. Compare MetaInfoInterface.h:71 comment explicitly documenting the temporary-default lifetime hazard.
- **Fix:** Either match MetaInfoInterface and add a by-value overload for the (name, default) / (index, default) forms, or at minimum document in the header that `default_value` must outlive the returned reference. ABI-safe additive fix: keep the existing reference overloads, add value-returning overloads or improve docs; the ideal fix (changing the return type to by-value) would be source/ABI breaking.
- **Verified:** Verified against actual code. MetaInfo.h:90 declares `const DataValue& getValue(const std::string& name, const DataValue& default_value = DataValue::EMPTY) const;` and MetaInfo.cpp:87 returns `default_value` by reference on a miss (same for the UInt overload at cpp:97). So if the caller passes a temporary as default_value and the lookup misses, the returned const& binds to a parameter that itself bound to a now-destroyed temporary — a genuine dan

### [META-17] MetaInfo::begin / MetaInfo::end (non-const) — Non-const begin()/end() hand out mutable iterators into a sorted flat_map, letting callers silently break the key-ordering invariant
`medium` · `const-correctness` · ABI: `none` · src/openms/include/OpenMS/METADATA/MetaInfo.h

```cpp
iterator begin(); iterator end();
```
- **Expectation:** Mutating iterators over an associative/sorted container should not let the caller corrupt the container's sort invariant; either keys are immutable (as in std::map, where key_type is const) or no mutable key access is given.
- **Actual:** begin()/end() return `boost::container::flat_map<UInt,DataValue>::iterator` whose `first` (the registry index / key) is assignable. The header itself admits the hazard (line 158: "@note Modifying the key (first element) may invalidate the sorted order invariant"). A caller can write `it->first = N;` and silently destroy the flat_map's required sorted order, after which find()/getValue()/exists() give wrong results with no error.
- **Evidence:** MetaInfo.h:160-166 mutable begin()/end(); MetaInfo.h:47 MapType is a sorted `boost::container::flat_map`; the @note at MetaInfo.h:158 acknowledges the invariant can be broken.
- **Fix:** Prefer not exposing mutable key access at all. Lowest-risk additive options: document the invariant hazard (already partially done) more strongly, or provide a value-only mutable accessor. Removing the mutable begin()/end() would be source-breaking, so flag the ideal fix but keep it documented.
- **Verifier correction:** The claim is accurate in substance and evidence. The only correction is severity: high → medium. Silent wrong-result corruption is real but is not triggered by normal use of the mutable iterators (mutating .second is safe and is the intended use); it requires a caller to deliberately assign to it->first, an atypical key write. The boost-vs-std::map divergence (flat_map value_type is std::pair<Key,T>, key not const — confirmed at flat_map.hpp:88-89,181) makes this a genuine, non-idiomatic footgun rather than a high-likelihood data-loss path.
- **Verified:** Independently verified against the actual code. MetaInfo.h line 47 defines MapType = boost::container::flat_map<UInt, DataValue>; lines 160 and 166 expose non-const begin()/end() returning MapType::iterator; the @note at line 158 verbatim acknowledges "Modifying the key (first element) may invalidate the sorted order invariant" — all quoted evidence is exact. The technical crux is real and provable from the boost header itself (/usr/include/boost

### [META-18] MetaInfoDescription::operator< — operator< is not a strict weak ordering consistent with operator==, breaking ordered-container use
`medium` · `other` · ABI: `source-compatible` · src/openms/include/OpenMS/METADATA/MetaInfoDescription.h

```cpp
bool operator<(const MetaInfoDescription& rhs) const
```
- **Expectation:** For a type that defines both operator< and operator== and the full relational set, the ordering must be a strict weak ordering whose equivalence (`!(a<b) && !(b<a)`) matches operator==, otherwise std::set/std::map/std::sort behave incorrectly.
- **Actual:** operator< compares DataProcessing only by vector size and per-element `getProcessingActions().size()` (MetaInfoDescription.cpp:31-88), never by content, while operator== compares DataProcessing pointees for value equality (cmpPtrSafe) and meta-value equality. Two objects with identical sizes but different DataProcessing/meta CONTENT are unequal under operator== yet ordering-equivalent (`!(a<b) && !(b<a)`). Inserting them into a std::set<MetaInfoDescription> would treat distinct objects as duplicates. operator>= is `!(a<b)` and operator<= is `a<b || a==b`, which are mutually inconsistent given this equivalence/equality mismatch.
- **Evidence:** MetaInfoDescription.cpp:31-88 (size-only DP comparison), :72-85 only compares `getProcessingActions().size()`; operator== at :20-29 compares pointee content via Helpers::cmpPtrSafe. operator<= at :90-93 vs operator>= at :100-103.
- **Fix:** Make operator< a true total order consistent with operator== (compare the same fields == compares, in a fixed order), ideally via a defaulted `operator<=>` over the same members; <compare> is already included (MetaInfoDescription.h:11). Body-only change, ABI-compatible.
- **Verifier correction:** The mismatch and the operator<=/operator>= inconsistency are confirmed exactly as described. However, the audited type is not currently used in any ordered container (no std::set<MetaInfoDescription>, std::map keyed on it, or std::sort over it exists in the codebase; only an unordered vector<pair<string,MetaInfoDescription>> in MzDataHandler). Therefore the impact is latent (invites future misuse, ships a broken strict-weak-ordering and self-inconsistent relational set) rather than producing silently wrong results today, so severity is medium, not high. The fix recommendation (make operator< consistent with operator==, e.g. a defaulted operator<=> over the same members) is sound and body-level/ABI-friendly; note the equality compares pointee contents, so a defaulted <=> over the raw DataProcessingPtr vector would compare pointers, not pointees — a content-aware comparison is required to truly match operator==.
- **Verified:** Independently confirmed in src/openms/source/METADATA/MetaInfoDescription.cpp. operator< (lines 31-88) compares the DataProcessing vector first by size (34-35), then per element ONLY by getProcessingActions().size() (76-79) plus a null-vs-nonnull tiebreak (81-84) — never by actual content (no software_, no actions-set contents, no completion_time_, no per-element meta). operator== (20-29) compares pointee value equality via Helpers::cmpPtrSafe, w

### [META-19] CVTermList::setCVTerms — setCVTerms appends instead of replacing existing terms
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/METADATA/CVTermList.h

```cpp
void setCVTerms(const std::vector<CVTerm>& terms)
```
- **Expectation:** A 'set' accessor (paired with getCVTerms) is expected to overwrite/replace the container's contents, leaving the object holding exactly the supplied terms.
- **Actual:** The implementation loops calling addCVTerm(tr), which does cv_terms_[accession].push_back(...). Existing terms are NOT cleared, so terms accumulate; calling setCVTerms twice, or calling it on a non-empty object, leaves stale terms in the map. Even within one call, two input terms sharing an accession are both kept (append), not deduplicated/replaced.
- **Evidence:** CVTermList.cpp: void CVTermList::setCVTerms(const vector<CVTerm>& cv_terms){ for (const CVTerm& tr : cv_terms){ addCVTerm(tr); } } and addCVTerm: cv_terms_[cv_term.getAccession()].push_back(cv_term);
- **Fix:** Make setCVTerms clear cv_terms_ before adding (matching 'set' semantics), or if accumulation is intended, deprecate-and-add a clearly named overload (e.g. appendCVTerms) and document setCVTerms as additive. A pure clear-then-add fix is behavior-changing but source/ABI compatible; choose based on existing callers.
- **Verifier correction:** setCVTerms is additive, not replacing: it appends via addCVTerm without clearing cv_terms_ first. The surprise is real and amplified by sibling methods replaceCVTerm/replaceCVTerms (which assign/overwrite) and consumeCVTerms (documented additive) — setCVTerms duplicates consume semantics under a 'set' name. Severity downgraded from high to medium: the only path that yields wrong results is calling it on a reused/pre-populated object (uncommon; no non-test in-tree callers do this) and the outcome is duplicate/stale terms that are loud-ish and recoverable, not silent cross-state corruption or a crash. Recommended fix (clear cv_terms_ at the start of setCVTerms) is source- and ABI-compatible — same signature/symbol, behavior-only change — so abi_impact is none.
- **Verified:** Evidence confirmed verbatim. CVTermList.cpp:32-39 setCVTerms loops over addCVTerm, and addCVTerm (line 29) does cv_terms_[accession].push_back(...). There is no clear() before adding, so the method is additive: calling it on a non-empty/reused object accumulates stale terms, and two input terms sharing an accession are both retained rather than replaced. This is a genuine misleading-name surprise, sharpened by the fact that the SAME class provide

### [META-24] CVTermList::addCVTerm — addCVTerm silently keys terms by accession, storing accession-less terms under an empty key
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/METADATA/CVTermList.h

```cpp
void addCVTerm(const CVTerm& term)
```
- **Expectation:** Adding a CV term that lacks an accession (the map key) should either be rejected or signaled, since CV terms are identified by accession; a caller would not expect a malformed term to be silently filed under the empty-string key.
- **Actual:** addCVTerm does cv_terms_[cv_term.getAccession()].push_back(cv_term) with no validation. A term with an empty accession is silently stored under key "", where it will never be found by a meaningful hasCVTerm() lookup. The source itself flags this: '// TODO exception if empty'.
- **Evidence:** CVTermList.cpp line 26-30: 'void CVTermList::addCVTerm(const CVTerm& cv_term){ // TODO exception if empty\n    cv_terms_[cv_term.getAccession()].push_back(cv_term); }'
- **Fix:** Either throw/log on empty accession (behavior change) or document explicitly that addCVTerm requires a non-empty accession and that the accession is used as the map key. A throwing change is source-compatible at the type level; prefer documenting plus an opt-in validating overload.
- **Verifier correction:** addCVTerm does silently store an accession-less CVTerm under the empty-string map key (confirmed, including the in-source "// TODO exception if empty"). However the impact is over-stated as data loss: the term is still retrievable via getCVTerms() under key "", no other entries are corrupted, there is no crash, and serialization/equality round-trips. It is a latent misfiling/lookup-invisibility trap rather than guaranteed silently-wrong results in normal use, so severity is medium, not high. Recommendation stands: document the non-empty-accession precondition and/or add a validating overload or throw/log on empty accession.
- **Verified:** Evidence is verbatim-correct. CVTermList.cpp:26-30 is exactly `void CVTermList::addCVTerm(const CVTerm& cv_term){ // TODO exception if empty\n  cv_terms_[cv_term.getAccession()].push_back(cv_term); }`. The map is keyed by accession (cv_terms_ is map<string, vector<CVTerm>>), hasCVTerm() is a pure key lookup (cv_terms_.contains(accession)), and CVTerm freely allows an empty accession (default-constructible std::string accession_, setAccession does

### [META-26] Precursor::getDriftTimeWindowLowerOffset — Drift-time lower offset is ADDED to get the window start, but the m/z lower offset is SUBTRACTED
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/METADATA/Precursor.h

```cpp
double getDriftTimeWindowLowerOffset() const
```
- **Expectation:** By analogy with the m/z isolation window in the same class, a 'lower offset' is the distance below the target: window start = target - lowerOffset. The two offset families are named identically (getIsolationWindowLowerOffset / getDriftTimeWindowLowerOffset) so a caller naturally assumes the same sign convention.
- **Actual:** The m/z doc says the start is 'p.getMZ() - p.getIsolationWindowLowerOffset()' (line 131), but the drift-time doc says the start is 'p.getDriftTime() + p.getDriftTimeWindowLowerOffset()' (line 181). The lower offset is ADDED, not subtracted. Both the lower and upper drift offsets use '+' in their docs (lines 181 and 195), so the lower offset cannot be a signed magnitude below the target the way the m/z lower offset is.
- **Evidence:** Precursor.h line 131: 'p.getMZ() - p.getIsolationWindowLowerOffset()'; line 181: 'p.getDriftTime() + p.getDriftTimeWindowLowerOffset()'; line 195: 'p.getDriftTime() + p.getDriftTimeWindowUpperOffset()'. setDriftTimeWindowLowerOffset enforces bound >= 0 (Precursor.cpp line 273), confirming a non-negative magnitude, so '+' for the lower bound is genuinely inconsistent with the m/z lower-offset '-'.
- **Fix:** Treat this as a documentation/contract bug: the drift-time lower-offset doc should compute the window start as 'getDriftTime() - getDriftTimeWindowLowerOffset()' to match the m/z convention (or, if the asymmetry is intentional, prominently document why the sign differs). No signature change needed; this is a source-compatible doc fix. Verify against existing callers in the codebase before changing the documented arithmetic.
- **Verified:** Verified directly in the source. Precursor.h:131 documents the m/z isolation window start as `p.getMZ() - p.getIsolationWindowLowerOffset()` (minus), while Precursor.h:181 documents the ion-mobility window start as `p.getDriftTime() + p.getDriftTimeWindowLowerOffset()` (plus); line 195 also uses plus for the upper drift offset. Both lower offsets are constrained to be non-negative magnitudes — `setDriftTimeWindowLowerOffset` has OPENMS_PRECONDITI

### [META-27] Precursor::getUnchargedMass — getUnchargedMass() silently fabricates charge 2 when charge is unknown (0)
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/METADATA/Precursor.h

```cpp
inline double getUnchargedMass() const
```
- **Expectation:** A getter named getUnchargedMass() on a precursor with unknown charge should signal the missing information (return a sentinel, throw, or at least be named to convey the guess), not silently return a plausible-looking mass computed from an invented charge.
- **Actual:** When charge_ == 0 the method substitutes c = 2 and returns getMZ()*2 - 2*PROTON_MASS_U. The returned value looks like a real measured mass but is a guess; a caller that does not read the doc comment will treat unknown-charge precursors as confidently doubly charged and silently get wrong masses.
- **Evidence:** Precursor.h lines 216-221: 'int c = charge_; (c == 0) ? c = 2 : c = charge_; return getMZ() * c - c * Constants::PROTON_MASS_U;' with comment 'if charge is unknown, i.e. 0, our best guess is doubly charged'. (The ternary is also dead-code-confusing: the false branch assigns c = charge_ which it already equals.)
- **Fix:** Keep the existing method for ABI stability but make the guess explicit at the call site: add an overload/parameter like getUnchargedMass(int default_charge) or a separate hasKnownCharge() check, and document the silent z=2 fallback in the brief (not only mid-sentence). At minimum simplify the ternary to 'int c = (charge_ == 0) ? 2 : charge_;'. All additive/source-compatible.
- **Verifier correction:** getUnchargedMass() on Precursor silently substitutes charge 2 when charge_==0, returning getMZ()*2 - 2*PROTON_MASS_U — a plausible-looking but fabricated mass. The behavior is real and surprising (notably the sibling LogMzPeak::getUnchargedMass() returns 0.0 as a sentinel for unknown charge, the opposite contract for an identically-named method). It is, however, documented inline at the declaration and in the pyOpenMS docstring, and the Precursor variant has no production caller exercising the charge-0 path, so the practical blast radius is limited: severity medium (invites misuse, documented/recoverable) rather than high. Fix is additive/source-compatible: add an explicit overload getUnchargedMass(int default_charge) and/or a hasKnownCharge()/getCharge()==0 guard, surface the z=2 fallback in the brief, and simplify the confusing ternary to `int c = (charge_ == 0) ? 2 : charge_;`.
- **Verified:** Evidence is verbatim-correct. Precursor.h:215-221 contains exactly: comment "if charge is unknown, i.e. 0, our best guess is doubly charged", then `int c = charge_; (c == 0) ? c = 2 : c = charge_; return getMZ() * c - c * Constants::PROTON_MASS_U;`. So a default-constructed Precursor (charge_{} == 0) returns a confident-looking mass computed from an invented z=2 rather than signaling missing info. This is a genuine least-surprise violation, reinf

### [META-29] SpectrumSettings::getDataProcessing — const getDataProcessing() returns a NEW vector by value (a copy), unlike its sibling const getters
`medium` · `return-value` · ABI: `breaking` · src/openms/include/OpenMS/METADATA/SpectrumSettings.h

```cpp
const std::vector<std::shared_ptr<const DataProcessing>> getDataProcessing() const
```
- **Expectation:** Every other const accessor in this class (getPrecursors, getProducts, getInstrumentSettings, getSourceFile, getAcquisitionInfo) returns 'const T&' to the internal member. By the same naming pattern, const getDataProcessing() looks like a cheap reference accessor, so a caller may bind it to 'const auto&' assuming it aliases the member.
- **Actual:** The const overload returns 'std::vector<std::shared_ptr<const DataProcessing>>' BY VALUE, constructed fresh each call via Helpers::constifyPointerVector (a per-call allocation/copy of the pointer vector). Its element type also differs from the member (shared_ptr<const DataProcessing> vs DataProcessingPtr). Binding the result to 'const auto&' extends a temporary's lifetime but never aliases the member, and the call is silently O(n)-allocating, not O(1).
- **Evidence:** SpectrumSettings.h lines 163-166: non-const overload returns 'std::vector<DataProcessingPtr>&' but const overload returns 'const std::vector<std::shared_ptr<const DataProcessing>>' (no '&'); SpectrumSettings.cpp lines 212-215 'return OpenMS::Helpers::constifyPointerVector(data_processing_);'. Identical pattern in ChromatogramSettings.h lines 134-137 / ChromatogramSettings.cpp 181-189.
- **Fix:** Cannot change return type without ABI break; document clearly that the const overload returns a freshly-built copy (not a reference) and is not a cheap accessor, so callers in hot loops cache the result. The const-correctness of constifying is reasonable; the surprise is the silent copy and the type asymmetry — call those out in the doc.
- **Verifier correction:** The per-call copy is NOT performed by Helpers::constifyPointerVector — that helper returns a `const std::vector<std::shared_ptr<const T>>&` reference via reinterpret_cast and is O(1) zero-copy. The copy is instead materialized at the return statement of `SpectrumSettings::getDataProcessing() const`, because the member function's declared return type is by-value (`const std::vector<std::shared_ptr<const DataProcessing>>`, no `&`) while the helper yields a reference; the by-value return type forces a fresh vector copy (O(n) shared_ptr copies, each an atomic refcount increment). Everything else in the claim holds: by-value return vs sibling `const T&` getters, `const auto&` binds a temporary (safe, lifetime-extended, never aliases the member), and the element-type asymmetry (shared_ptr<const DataProcessing> vs DataProcessingPtr=shared_ptr<DataProcessing>). Severity is medium, not high: binding to `const auto&` does not dangle and does not yield wrong values — the only harms are a silent O(n) copy in hot loops (perf footgun, fully recoverable) and potential surprise in generic code from the differing element type. Recommendation stands: document that this const overload returns a freshly-built copy (not a member reference) and that callers in hot loops should cache it; do not claim the cost lives inside constifyPointerVector.
- **Verified:** Verified against source. The const overload (SpectrumSettings.h:166, .cpp:212-215) genuinely returns `const std::vector<std::shared_ptr<const DataProcessing>>` BY VALUE (no `&`), while every sibling const getter (getPrecursors/getProducts/getInstrumentSettings/getSourceFile/getAcquisitionInfo, .h:125-153) returns `const T&`. So a caller who writes `const auto&` gets a lifetime-extended temporary that never aliases data_processing_, and pays a per

### [META-30] Acquisition::getIdentifier / setIdentifier — 'identifier' is documented as an index/number but typed and accessed as a free-form string
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/METADATA/Acquisition.h

```cpp
const std::string& getIdentifier() const; void setIdentifier(const std::string& identifier)
```
- **Expectation:** The accessor names and docs describe a numeric scan index/number ('return the identifier/index/number of the acquisition', 'sets the index/number of the scan'), so a caller may expect to set/read it as an integer index and to be able to use it for ordering or array lookup.
- **Actual:** The value is a std::string with no numeric validation or guarantee; it is whatever the acquisition software wrote. The getter/setter names disagree with each other too: getIdentifier 'return the identifier/index/number' vs setIdentifier 'sets the index/number'. A caller treating it as a numeric index can silently get a non-numeric native id.
- **Evidence:** Acquisition.h lines 46-49: '/// return the identifier/index/number of the acquisition' const std::string& getIdentifier() const; '/// sets the index/number of the scan' void setIdentifier(const std::string& identifier); backing field 'std::string identifier_;' (line 52).
- **Fix:** Source-compatible doc fix: describe the field consistently as a free-text/native acquisition identifier (string), and drop the 'index/number' wording or state that it is not guaranteed numeric. No signature change.
- **Verifier correction:** The field is a free-form native acquisition/scan identifier (std::string), not a guaranteed numeric index. On mzML load it is populated from the native externalSpectrumID (often non-numeric); only mzData populates it from the numeric acqNumber. Numeric misuse does not fail silently: StringUtils::toInt32(getIdentifier()) throws (or, as in the mzData writer, is caught and substituted with 0 plus a warning), so the failure mode is loud/recoverable rather than silently-wrong — hence medium, not high. Recommended fix is doc-only: describe getIdentifier/setIdentifier consistently as a native/free-text identifier and drop or qualify the "index/number" wording. ABI: none.
- **Verified:** Verified against source. Header (Acquisition.h:46-52) matches the quoted evidence exactly: getIdentifier doc says "return the identifier/index/number of the acquisition", setIdentifier doc says "sets the index/number of the scan", backing field is plain std::string identifier_. The .cpp (Acquisition.cpp:26-34) is a trivial string passthrough with zero numeric validation. The misleading-name surprise is real and provable from OpenMS itself: on mzM

### [META-34] DocumentIdentifier::operator== — operator== only compares the id_ field, silently ignoring loaded file path and type
`medium` · `silent-failure` · ABI: `none` · src/openms/source/METADATA/DocumentIdentifier.cpp

```cpp
bool operator==(const DocumentIdentifier& rhs) const
```
- **Expectation:** Equality on a value type holding {identifier, loaded file path, loaded file type} would compare all observable members, so two DocumentIdentifiers with different file paths/types are unequal.
- **Actual:** The operator returns `id_ == rhs.id_` only; two objects with identical (often empty) ids but different file_path_ / file_type_ compare equal. Callers using == to detect 'same document including provenance' get false positives.
- **Evidence:** DocumentIdentifier.cpp:82-85 `bool DocumentIdentifier::operator==(const DocumentIdentifier & rhs) const { return id_ == rhs.id_; }`
- **Fix:** Either document explicitly that equality is id-only, or extend it to compare file_path_ and file_type_. Behavioral change; widening the comparison risks breaking callers that rely on id-only semantics, so prefer documenting the current behavior plus an explicit equalsFull() helper.
- **Verifier correction:** operator== returns id_ == rhs.id_ only, ignoring file_path_ and file_type_, and this id-only semantics is nowhere documented (header line 54 just says "Equality operator"). The surprise is real but its practical blast radius is narrower than the claim implies: the operator is consumed by FeatureMap/ConsensusMap operator== chains (FeatureMap.cpp:149, ConsensusMap.cpp:568), which also compare full data content, so the effect is a silent provenance-only blind spot rather than broad data-level false positives. No current caller was found relying on == alone to assert "same document including provenance." Recommendation stands: document the id-only contract and/or add an explicit equalsFull() helper; widening == in place would be a source-compatible behavioral change that could affect the map equality chains.
- **Verified:** Evidence confirmed verbatim at DocumentIdentifier.cpp:82-85: `bool DocumentIdentifier::operator==(const DocumentIdentifier& rhs) const { return id_ == rhs.id_; }`. The class (header lines 89-93) holds three value-significant members (id_, file_path_, file_type_); copy ctor, assignment, and swap (lines 75-80) all preserve all three, but == compares only id_. The header comment (line 54) says merely "Equality operator" with no note that equality is

### [META-35] DocumentIdentifier::setLoadedFileType — setLoadedFileType opens and reads the file from disk instead of recording a passed-in type
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/METADATA/DocumentIdentifier.h

```cpp
void setLoadedFileType(const std::string& file_name)
```
- **Expectation:** A setter named setLoadedFileType, sibling to setLoadedFilePath, would simply store the type the caller provides (and the comment 'set the file_type according to the type of the file loaded from' reinforces that it just records it).
- **Actual:** The implementation calls FileHandler::getTypeByContent(file_name), which opens the file on disk and sniffs its content to detect the type. A trivial 'setter' thus performs file I/O, can be slow, and can fail/throw on a missing or unreadable path; the param is named file_name though it is the path that gets read.
- **Evidence:** DocumentIdentifier.cpp:65-68 `void DocumentIdentifier::setLoadedFileType(const std::string & file_name){ file_type_ = FileHandler::getTypeByContent(file_name); }`
- **Fix:** Document the I/O side effect in the header, and add an overload taking a FileTypes::Type directly (e.g. setLoadedFileType(FileTypes::Type)) for callers that already know the type. Additive overload; no ABI break.
- **Verified:** Evidence is exact: DocumentIdentifier.cpp:65-68 is `void DocumentIdentifier::setLoadedFileType(const std::string & file_name){ file_type_ = FileHandler::getTypeByContent(file_name); }`. Independently verified that FileHandler::getTypeByContent (FileHandler.cpp:340+) opens the file via ifstream, reads bytes, transparently handles bzip2/gzip/zip decompression, and sniffs content to detect the type — genuine, potentially slow disk I/O. FileHandler.h

### [META-37] IdentifierMSRunMapper::getIdentifier / getMSRunPaths — Reverse-direction lookups handle 'not found' inconsistently: one throws, the sibling returns empty
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/METADATA/IdentifierMSRunMapper.h

```cpp
const std::string& getIdentifier(const StringList&) const; const StringList& getMSRunPaths(const std::string&) const
```
- **Expectation:** The two symmetric reverse lookups in the same class (path->identifier and identifier->paths) should treat a missing key the same way, so callers can reason about error handling uniformly.
- **Actual:** getIdentifier(ms_run_paths) throws Exception::ElementNotFound when the key is absent, while the mirror getMSRunPaths(identifier) silently returns an empty StringList. A caller who learned one convention will mis-handle the other.
- **Evidence:** IdentifierMSRunMapper.cpp:87-95 getIdentifier throws ElementNotFound; IdentifierMSRunMapper.cpp:107-115 getMSRunPaths returns empty_stringlist_. Header lines 67-68 vs 76-77 document the divergence.
- **Fix:** Keep the behaviors (both are individually documented) but make the asymmetry explicit in adjacent docs, or add a throwing getMSRunPaths variant / non-throwing tryGetIdentifier for parity (the latter already exists). Documentation/additive; no ABI break.
- **Verified:** Code matches the claim exactly. In the same class (a documented "two-way mapping"), the two symmetric reverse lookups handle a missing key oppositely: getIdentifier(const StringList&) at cpp:87-95 throws Exception::ElementNotFound, while the mirror getMSRunPaths(const std::string&) at cpp:107-115 silently returns the static empty_stringlist_. Both are named get* and neither name signals the behavioral difference (unlike std::map::at vs operator[]

### [META-38] AnnotatedMSRun::begin / cbegin (const) — const begin()/cbegin() throw if spectra and peptide-id counts disagree
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/METADATA/AnnotatedMSRun.h

```cpp
auto begin() const; auto cbegin() const
```
- **Expectation:** begin()/cbegin(), especially const ones used in range-for, are expected to be cheap, noexcept-ish accessors that never throw.
- **Actual:** Every begin()/cbegin()/begin() const calls checkPeptideIdSize_, which throws Exception::InvalidValue when data.getSpectra().size() != peptide_ids_.size(). Since the setters (setMSExperiment, setPeptideIdentifications) accept mismatched sizes without validating, a perfectly normal construction sequence can make a subsequent `for (auto x : run)` throw from begin().
- **Evidence:** AnnotatedMSRun.h:170-194 each begin/cbegin calls `checkPeptideIdSize_(OPENMS_PRETTY_FUNCTION);`; AnnotatedMSRun.cpp:53-62 checkPeptideIdSize_ throws InvalidValue on size mismatch; setters at AnnotatedMSRun.cpp:23-51 perform no size check.
- **Fix:** Document that iteration requires matching sizes and that begin() validates the invariant; consider validating in the setters (or providing a resize/sync helper) so the throw surfaces at mutation time rather than at iteration time. Doc + optional behavioral change; no ABI break.
- **Verifier correction:** Claim is accurate as written; only the severity is adjusted from (implied) high to medium. The behavior is real: const begin()/cbegin() (and non-const begin()) call checkPeptideIdSize_, which throws Exception::InvalidValue on a spectra/peptide-id size mismatch, while the setters and the MSExperiment&& constructor accept/produce mismatched sizes without validation, so a normal construction sequence makes a later range-for throw from begin(). The mitigation is documentation plus optionally validating in setters or providing a sync/resize helper; this is doc + optional behavioral change with no ABI break. Severity is medium (not high) because the invariant violation surfaces as a loud thrown exception rather than silently wrong results, data loss, or a crash.
- **Verified:** Independently confirmed against the actual code. Header lines 170-194: cbegin() const, begin(), and begin() const each call checkPeptideIdSize_(OPENMS_PRETTY_FUNCTION). CPP lines 53-62: checkPeptideIdSize_ throws Exception::InvalidValue when data.getSpectra().size() != peptide_ids_.size(). CPP lines 23-51 setters (setMSExperiment / setPeptideIdentifications) perform NO size validation, and the header constructor at line 54 (explicit AnnotatedMSRu

### [META-39] ExperimentalDesign::setMSFileSection — setMSFileSection silently reorders the input rows via sort_()
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/METADATA/ExperimentalDesign.h

```cpp
void setMSFileSection(const MSFileSection& msfile_section)
```
- **Expectation:** A plain setter storing the provided MSFileSection would keep the rows in the caller-supplied order, so a later getMSFileSection() round-trips the same sequence.
- **Actual:** setMSFileSection stores the vector and then calls sort_(), reordering entries by (fraction_group, fraction, label, sample, path). Code that sets rows and relies on positional indices, or expects get to return exactly what was set, is silently broken; nothing in the declaration hints at reordering.
- **Evidence:** ExperimentalDesign.cpp:505-509 `void ExperimentalDesign::setMSFileSection(const MSFileSection& msfile_section){ msfile_section_ = msfile_section; sort_(); }`; sort_ at ExperimentalDesign.cpp:742-750 sorts by the tuple key.
- **Fix:** Document in the header that the section is sorted into canonical order on set. Documentation-only; no ABI impact.
- **Verified:** Verified directly in source. ExperimentalDesign.cpp:505-509 shows setMSFileSection assigns msfile_section_ = msfile_section; then calls sort_(). sort_() at lines 742-750 sorts the vector by std::tie(fraction_group, fraction, label, sample, path) — exactly the key claimed. getMSFileSection() (lines 500-503) returns the stored (now sorted) vector, so set-then-get does NOT round-trip the caller's row order. The public declaration in the header (line

### [META-40] ExperimentalDesign::setSampleSection — Two-arg constructor validates the design (isValid_) but the equivalent set* sequence skips validation
`medium` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/METADATA/ExperimentalDesign.h

```cpp
void setSampleSection(const SampleSection&); void setMSFileSection(const MSFileSection&)
```
- **Expectation:** Building a design via the (msfile_section, sample_section) constructor and building it via setMSFileSection + setSampleSection should yield equally-validated objects.
- **Actual:** The constructor calls sort_() then isValid_(), which enforces consecutive fraction groups starting at 1, unique (path,label), etc. The setters call only sort_() (setMSFileSection) or nothing (setSampleSection), so an object assembled through setters can hold an invalid design that the constructor would have rejected, surprising callers who expect uniform invariants.
- **Evidence:** Ctor ExperimentalDesign.cpp:41-49 `sort_(); isValid_();`; setMSFileSection ExperimentalDesign.cpp:505-509 calls only sort_(); setSampleSection ExperimentalDesign.cpp:511-514 just assigns.
- **Fix:** Document that setters do not validate (or call isValid_() from setMSFileSection for parity). Prefer documentation plus an optional validate() entry point to avoid breaking existing incremental-build callers. No ABI break for the doc route.
- **Verifier correction:** Precise statement: setMSFileSection DOES call sort_() (so it is partially symmetric with the ctor) but skips isValid_(); setSampleSection skips both sort_() and isValid_(). The asymmetry is specifically that the validation step (isValid_()) performed by the two-arg ctor is absent from both setters, and isValid_() is private with no public validation entry point, so a setter-assembled design cannot be validated to ctor parity. Severity medium (not high): the invalid state requires a caller to actively bypass the obvious validating constructor; the primary file-ingest path goes through the validated ctor and the failure mode at set time is silent rather than corrupting (no crash), but downstream quant code (see comment at ExperimentalDesign.cpp:681 noting order is 'assumed') may then produce wrong mappings. Recommendation stands: document the non-validating behavior on the setters, optionally call isValid_() from the setters or expose a public validate() for parity; doc route is ABI-none.
- **Verified:** Evidence verified verbatim in src/openms/source/METADATA/ExperimentalDesign.cpp. Ctor (lines 41-49) runs sort_() then isValid_(); setMSFileSection (505-509) runs only sort_(); setSampleSection (511-514) is a bare assignment with no sort_() and no isValid_(). isValid_() (637-720) enforces real invariants and actively throws: fraction groups must be consecutive starting at 1, (FractionGroup,Fraction,Label) unique, (Path,Label) unique, one sample pe

### [META-44] ProteinHit::operator=(const MetaInfoInterface&) — Assignment from MetaInfoInterface overwrites only meta values, silently leaving score/accession/sequence/coverage intact
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/METADATA/ProteinHit.h

```cpp
ProteinHit& operator=(const MetaInfoInterface& source)
```
- **Expectation:** operator= normally produces an object equal to the source; a developer writing `proteinHit = something;` would expect a full state assignment (or, knowing the param is a MetaInfoInterface, at least a clearly-documented partial copy).
- **Actual:** This public assignment operator copies only the MetaInfoInterface sub-object and leaves all ProteinHit-specific fields (score_, rank_, accession_, sequence_, coverage_, modifications_) untouched. The header comment is just "/// Assignment for MetaInfo" with no warning about the partial nature.
- **Evidence:** Header L120; ProteinHit.cpp L42-46 `MetaInfoInterface::operator=(source); return *this;`.
- **Fix:** Expand the header doc to state it assigns ONLY meta values and leaves hit fields unchanged; consider renaming to assignMetaInfo() in a future major version. abi_impact none for doc fix.
- **Verified:** Evidence verified exactly. ProteinHit.h L120 declares `ProteinHit& operator=(const MetaInfoInterface& source)` with only the comment `/// Assignment for MetaInfo`. ProteinHit.cpp L42-46 implements it as `MetaInfoInterface::operator=(source); return *this;` — copying ONLY the meta sub-object and leaving score_/rank_/accession_/sequence_/coverage_/modifications_ untouched. The unit test (ProteinHit_test.cpp L97-110) confirms this is intentional: af

### [META-47] PeptideHit::getTargetDecoyType / isDecoy — getTargetDecoyType() (and isDecoy()) throws on an unrecognized 'target_decoy' meta value string
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/METADATA/PeptideHit.h

```cpp
TargetDecoyType getTargetDecoyType() const
```
- **Expectation:** A getter that is documented to "Return UNKNOWN if the meta value does not exist" would be expected to also return UNKNOWN (or some non-throwing result) for an unparseable/foreign value — getters generally don't throw.
- **Actual:** If the 'target_decoy' meta value exists but is not one of decoy/target/target+decoy (case-insensitive), getTargetDecoyType() throws Exception::InvalidValue. Because isDecoy() calls getTargetDecoyType(), a simple `hit.isDecoy()` can throw on data written by other tools using a different convention (e.g. "1"/"0", "reverse"). The header doc for getTargetDecoyType()/isDecoy() does not mention this throw.
- **Evidence:** PeptideHit.cpp L407-420 throws Exception::InvalidValue for unrecognized strings; isDecoy() L383-386 delegates to it. Header L335-347 doc mentions only the UNKNOWN-when-missing case.
- **Fix:** Document the throw on getTargetDecoyType()/isDecoy(), or make isDecoy()/getTargetDecoyType() return UNKNOWN for unrecognized values instead of throwing. abi_impact none for doc; behavior change is source-compatible.
- **Verifier correction:** Severity is medium, not high: the failure is loud (a thrown Exception::InvalidValue, recoverable) rather than silently wrong results or data corruption. It can still abort a TOPP tool mid-run on otherwise-valid input from another tool, with no documentation warning, and the 36 unguarded isDecoy() call sites assume a non-throwing bool. abi_impact of the documentation fix is none; the recommended behavior change (return UNKNOWN instead of throwing) is source-compatible (broadens runtime behavior, signature unchanged).
- **Verified:** Independently verified in the actual code. PeptideHit.cpp L407-419: getTargetDecoyType() returns UNKNOWN only when the meta value is absent (L409-411); if the value EXISTS but lowercased is not "decoy"/"target+decoy"/"target", it throws Exception::InvalidValue (L419). isDecoy() (L383-386) delegates directly to getTargetDecoyType(), so it inherits the throw. The header documents only the benign cases: getTargetDecoyType() doc (L335-347) lists four

### [META-48] ProteinModificationSummary::operator== — Equality compares modification keys by raw pointer identity, not by modification value
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/METADATA/ProteinModificationSummary.h

```cpp
bool operator==(const ProteinModificationSummary& rhs) const
```
- **Expectation:** operator== on a value-like struct should compare two summaries as equal when they describe the same modifications and statistics at the same positions.
- **Actual:** AALevelSummary is a std::map keyed on `const ResidueModification*`; map equality compares the pointer values. Two summaries describing the identical modification but holding pointers to different ResidueModification instances (e.g. copied/reconstructed objects) compare unequal, while ordering also depends on pointer addresses (non-deterministic across runs). The header does not warn about pointer-identity semantics for comparison.
- **Evidence:** Header L40 `using ModificationsToStatistics = std::map<const ResidueModification*, Statistics>;`; ProteinModificationSummary.cpp operator== compares `AALevelSummary == rhs.AALevelSummary` (pointer-keyed map equality).
- **Fix:** Document that comparison/ordering relies on ResidueModification pointer identity (valid only when mods come from the shared ModificationsDB singleton), or compare via getFullId() like the ProteinHit hash does. abi_impact none for doc fix.
- **Verifier correction:** operator== reduces to pointer-identity comparison of ResidueModification keys: summaries with chemically identical mods but distinct ResidueModification* instances compare unequal. This is correct only under the (unstated) convention that all mod pointers originate from the ModificationsDB singleton; off that path (reconstructed/copied/deserialized mods) it silently returns false. The 'non-deterministic ordering across runs' aspect affects map iteration order but not the equality result, so it does not itself make operator== nondeterministic. Recommendation stands: either document the pointer-identity precondition or compare via ResidueModification::getFullId() as ProteinHit's hashing does. Severity medium (not high) because the canonical singleton path works and the struct currently has no callers; abi_impact none for a doc-only fix.
- **Verified:** Evidence is verbatim-confirmed. Header L40: `using ModificationsToStatistics = std::map<const ResidueModification*, Statistics>;`, outer map L41 keyed by size_t position; .cpp L15 `return AALevelSummary == rhs.AALevelSummary;`. The inner map is keyed by a raw `const ResidueModification*`, so std::map equality compares keys by pointer-address equivalence (default std::less<const ResidueModification*>). Two ProteinModificationSummary objects descri

### [META-49] PeptideIdentification::setExperimentLabel — setExperimentLabel("") silently no-ops, so it cannot clear a previously set label (unlike sibling setBaseName)
`medium` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/METADATA/PeptideIdentification.h

```cpp
void setExperimentLabel(const std::string& type)
```
- **Expectation:** Consistent with the sibling setBaseName(), passing an empty string should remove the stored value, so set("") clears the field.
- **Actual:** setExperimentLabel() stores only non-empty labels and, unlike setBaseName() which calls removeMetaValue() for empty input, does NOT remove an existing 'experiment_label' meta value. So calling setExperimentLabel("") to reset leaves the old label in place, while getExperimentLabel() participates in operator==, producing surprising inequality after an apparent 'clear'.
- **Evidence:** PeptideIdentification.cpp setExperimentLabel L198-205 (only sets when non-empty; no else/remove), versus setBaseName L178-189 which removeMetaValue on empty. operator== uses getExperimentLabel() (L43).
- **Fix:** Make setExperimentLabel("") remove the meta value to mirror setBaseName, or document the asymmetry. abi_impact source-compatible (behavior change) / none for doc.
- **Verifier correction:** setExperimentLabel("") does not clear a previously set 'experiment_label' meta value (no else/removeMetaValue), unlike the sibling setBaseName("") which calls removeMetaValue. The result is a stale value surviving an apparent clear and surprising operator== inequality. Severity downgraded from high to medium: recoverable, loud-ish, and confined to a niche pepXML field rather than silently corrupting common results.
- **Verified:** Verified directly against source. PeptideIdentification.cpp setExperimentLabel (L198-205) stores only when !label.empty() and has NO else branch, so setExperimentLabel("") is a silent no-op that does not remove an existing 'experiment_label' meta value. The sibling setBaseName (L178-189) has an explicit else calling removeMetaValue() on empty input, so setBaseName("") clears the field. Because getExperimentLabel() (L191-196) is getMetaValue("expe

### [META-50] IdentifiedSequence::allParentsAreDecoys — Boolean predicate allParentsAreDecoys() throws on empty parent set instead of returning a bool
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/METADATA/ID/IdentifiedSequence.h

```cpp
bool allParentsAreDecoys() const
```
- **Expectation:** A const bool predicate named allParentsAreDecoys() returns true/false; for vacuous cases callers expect a defined boolean (conventionally 'all of empty' == true), not an exception.
- **Actual:** When parent_matches is empty it throws Exception::MissingInformation. A caller writing `if (seq.allParentsAreDecoys())` to filter decoys will get an unexpected exception for any sequence with no recorded parents.
- **Evidence:** IdentifiedSequence.h lines 64-71: `if (parent_matches.empty()) { std::string msg = "no parent found for identified molecule"; throw Exception::MissingInformation(...); }`
- **Fix:** Document the throw-on-empty contract prominently in the header, or add a defaulted parameter / overload that returns a chosen value for the empty case (e.g. `allParentsAreDecoys(bool value_if_empty)`). Changing the existing behavior to return a bool would be source-compatible but semantically breaking for current callers relying on the throw.
- **Verifier correction:** The throw is real and undocumented, but the impact is a loud, recoverable exception rather than silent corruption, so severity is medium, not high. allParentsAreDecoys() (header-defined template method, signature unchanged regardless of fix, so abi_impact = none) throws Exception::MissingInformation when parent_matches is empty -- a reachable state (de novo / pre-inference IDs pass checkParentMatches_ and live in IdentificationData), and the only caller (FalseDiscoveryRate) does not guard it.
- **Verified:** Evidence verified verbatim: IdentifiedSequence.h lines 64-77 show allParentsAreDecoys() throws Exception::MissingInformation when parent_matches.empty(), instead of returning a bool. The method has zero doxygen documentation, so the throw-on-empty contract is undocumented. The empty-parent state is genuinely reachable and normal: IdentificationData::checkParentMatches_ (IdentificationData.cpp:183-201) only validates the entries that exist and doe

### [META-52] IteratorWrapper::operator uintptr_t — IteratorWrapper implicitly converts to an integer pointer value
`medium` · `implicit-conversion` · ABI: `breaking` · src/openms/include/OpenMS/METADATA/ID/MetaData.h

```cpp
operator uintptr_t() const
```
- **Expectation:** An iterator-like wrapper should not silently convert to an integer; comparing or passing it where an integer is expected should not compile.
- **Actual:** The non-explicit conversion operator lets any IteratorWrapper (and thus every *Ref type built on it, e.g. ParentSequenceRef, ObservationRef) implicitly convert to uintptr_t. Expressions like `ref1 == 0`, `ref - 1`, or passing a ref to an integer parameter compile and silently use the element address, which is a surprising and non-portable value.
- **Evidence:** MetaData.h lines 32-35: `/// Conversion to pointer type for hashing\n operator uintptr_t() const { return uintptr_t(&(**this)); }` (not marked explicit).
- **Fix:** Mark the conversion `explicit` (source-compatible for legitimate hashing call sites that can add a cast; flags accidental integer uses). If a hash is the only need, provide a named hash functor / std::hash specialization instead of an implicit conversion.
- **Verifier correction:** The non-explicit operator uintptr_t() on IteratorWrapper (and thus all *Ref types) is a real implicit-conversion POLS surprise. However, contrary to the claim's recommendation, marking it `explicit` is NOT source-compatible: the conversion is depended upon by isValidHashedReference_ (IdentificationData.cpp:35, `lookup.count(ref)` over std::unordered_set<uintptr_t>) and ~8 call sites; those would need explicit uintptr_t() casts added. The cleanest fix is to make the operator explicit AND update isValidHashedReference_/related sites to cast explicitly, or replace the implicit conversion with a named accessor (e.g. asAddress()) or a std::hash specialization. Severity is medium (invites silent misuse, no current bug), and the recommended change is source-breaking within the module.
- **Verified:** Evidence verified exactly: MetaData.h lines 31-35 define a non-explicit `operator uintptr_t() const { return uintptr_t(&(**this)); }` on IteratorWrapper. Every *Ref type (ParentSequenceRef, ObservationRef, IdentifiedPeptideRef, ObservationMatchRef, etc.) is a typedef of IteratorWrapper<...>, so all inherit this implicit iterator->integer conversion. This is a genuine POLS violation: iterators do not implicitly convert to integers in standard C++,

### [META-53] Observation::merge — Observation::merge unconditionally overwrites rt/mz, unlike sibling merge() methods that conflict-check
`medium` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/METADATA/ID/Observation.h

```cpp
Observation& merge(const Observation& other)
```
- **Expectation:** Given the sibling merge() methods in this cluster (ParentSequence, InputFile, ObservationMatch) which preserve existing scalar fields and throw on conflicting values, a caller expects Observation::merge to merge rt/mz consistently (preserve/validate), not silently replace them.
- **Actual:** Observation::merge always assigns rt = other.rt and mz = other.mz, discarding the current object's position even when other holds NaN or a conflicting value.
- **Evidence:** Observation.h lines 48-55: `Observation& merge(...) { addMetaValues(other); rt = other.rt; mz = other.mz; return *this; }`. Contrast ParentSequence.h lines 60-66 which throw Exception::InvalidValue on conflicting non-empty values, and ObservationMatch.h lines 83-89 which throw on conflicting charge.
- **Fix:** Document that Observation::merge overwrites position from `other`, or align it with the sibling pattern (keep existing non-NaN value, take other's only when unset, throw on genuine conflict). Behavior change would be source-compatible but observably different.
- **Verifier correction:** Observation::merge (Observation.h:48-55) unconditionally overwrites rt and mz from `other` (last-writer-wins, no NaN guard, no conflict check), diverging from the consistent convention of its three sibling merge() methods (ParentSequence::merge, InputFile::merge, ObservationMatch::merge) which preserve already-set scalar fields and throw Exception::InvalidValue on conflicting values. Two qualifications to the original framing: the symbol is in the internal IdentificationDataInternal namespace (not a public stable API), and merge only triggers on re-registering an existing (input_file, data_id) observation, so a genuine rt/mz conflict is uncommon. When it does occur (e.g. a valid position overwritten by NaN depending on registration order), the position is silently lost with no diagnostic. Aligning it with the sibling pattern would be source-compatible but observably different at runtime.
- **Verified:** Verified against the actual source. Observation.h lines 48-55 confirm Observation::merge unconditionally does `addMetaValues(other); rt = other.rt; mz = other.mz;` with no NaN guard and no conflict check. The three sibling merge() methods in the same internal cluster all follow the opposite, consistent convention: keep the existing scalar if set, take other's only when unset/empty/zero, and throw Exception::InvalidValue on a genuine conflict — Pa

### [META-54] IdentifiedMolecule::getFormula — Untyped Size fragment_type is reinterpreted as two different enums and is silently ignored for compounds
`medium` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/METADATA/ID/IdentifiedMolecule.h

```cpp
EmpiricalFormula getFormula(Size fragment_type = 0, Int charge = 0) const
```
- **Expectation:** A parameter named fragment_type with a single integer type would be expected to mean the same thing regardless of the molecule, and to be honored.
- **Actual:** fragment_type is static_cast to Residue::ResidueType for peptides but to NASequence::NASFragmentType for oligos (two unrelated enums with different meanings per value), so the same integer means different fragment kinds depending on molecule type. For compounds both fragment_type and charge are silently ignored (the stored formula is returned unchanged). A caller passing e.g. 5 (BIon for peptides) to a compound or to an oligo gets surprising results.
- **Evidence:** IdentifiedMolecule.cpp lines 80-103: PROTEIN branch `static_cast<Residue::ResidueType>(fragment_type)`, RNA branch `static_cast<NASequence::NASFragmentType>(fragment_type)`, COMPOUND branch `// @TODO: what about fragment type and charge?` returns `getIdentifiedCompoundRef()->formula` ignoring both args.
- **Fix:** Document in the header that fragment_type is interpreted per-molecule-type (Residue::ResidueType vs NASequence::NASFragmentType) and that it (and charge) are ignored for compounds. Ideal additive fix: typed overloads or an enum wrapper; ABI-safe if added alongside.
- **Verifier correction:** fragment_type IS reinterpreted as two unrelated enums (Residue::ResidueType for peptides, NASequence::NASFragmentType for oligos) and is silently ignored together with charge for compounds, as the evidence states. However the claim's chosen example value 5 (BIon) actually maps to BIon in BOTH enums (the enums coincide at indices 0,1 and 4-16), so it is NOT the strongest illustration. The real per-value semantic divergence is at index 2 (NTerminal vs FivePrime) and index 3 (CTerminal vs ThreePrime), plus NAS-only values WIon/AminusB/DIon at the tail. The default value 0 (Full) is consistent across all three molecule types, so callers relying on defaults are unaffected; the surprise is for callers passing non-default peptide-style values to oligos or any non-default value to compounds. Severity is medium (wrong mass returned silently for non-default cross-molecule misuse, but recoverable and not hit on the common default path) rather than high.
- **Verified:** Verified against the actual code. IdentifiedMolecule.cpp lines 80-103 match the quoted evidence exactly: the PROTEIN branch does static_cast<Residue::ResidueType>(fragment_type), the RNA branch does static_cast<NASequence::NASFragmentType>(fragment_type), and the COMPOUND branch carries the literal `// @TODO: what about fragment type and charge?` comment and returns getIdentifiedCompoundRef()->formula, ignoring both args. The header (line 48) dec

### [META-55] ScoredProcessingResult::getScore — getScore returns NaN paired with a success flag the caller must check
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/METADATA/ID/ScoredProcessingResult.h

```cpp
std::pair<double, bool> getScore(ScoreTypeRef score_ref) const
```
- **Expectation:** A getter named getScore returning a double would be expected to yield the score; a caller may use the returned double directly.
- **Actual:** On lookup failure it returns {NaN, false}; the double is a NaN sentinel and the real failure signal is the second bool. A caller writing `double s = obj.getScore(t).first;` silently propagates NaN if they ignore the .second flag.
- **Evidence:** ScoredProcessingResult.h lines 99-106 (doc: "@return A pair: score (or NaN), success indicator") and getScore(score_ref, step) lines 126-127 `return std::make_pair(std::numeric_limits<double>::quiet_NaN(), false);`
- **Fix:** This is documented, so keep as-is for ABI; ensure header doc emphasizes that the bool must be checked and the double is NaN on miss. Consider an std::optional<double>-returning overload as an additive, clearer API.
- **Verifier correction:** The code and evidence are exactly as claimed. Adjustment is to severity only: this is medium, not high. The failure is documented (the doc string names NaN and the success indicator), recoverable (the .second bool is available), and NaN is a detectable poison value downstream; the silent-NaN propagation only occurs on an actual lookup miss rather than on normal/typical use. A concrete production caller that ignores the flag exists (src/openms/source/FORMAT/MzTabM.cpp:609), confirming the misuse is real but bounded. ABI impact is none — the suggested fix (doc clarification plus an additive std::optional<double> overload) does not change the existing signature.
- **Verified:** Verified against the actual header (src/openms/include/OpenMS/METADATA/ID/ScoredProcessingResult.h). All quoted evidence is exact: line 101 `std::pair<double, bool> getScore(ScoreTypeRef score_ref) const`, the doc "@return A pair: score (or NaN), success indicator" (line 99), and the NaN sentinel `return std::make_pair(std::numeric_limits<double>::quiet_NaN(), false);` (line 127). The surprise is genuine and not reject-worthy: it is not a mass-sp

### [META-2] HPLC::getFlux / HPLC::setFlux / HPLC::getPressure — Pressure/flux accessors are UInt but the backing members are signed Int (negative values round-trip as huge UInt)
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/METADATA/HPLC.h

```cpp
UInt getFlux() const; void setFlux(UInt flux);
```
- **Expectation:** `UInt getFlux()` returning an unsigned value implies the stored quantity is non-negative and that get(set(x))==x for the unsigned domain.
- **Actual:** The members are signed: `Int pressure_;` and `Int flux_;` (HPLC.h:87-88), but the API exposes `UInt getPressure()/setPressure(UInt)` and `UInt getFlux()/setFlux(UInt)`. setPressure/setFlux store a UInt into a signed Int, and getPressure/getFlux read a signed Int back as UInt. Any value above INT_MAX silently wraps to a negative member and re-emerges differently; nothing guarantees the unsigned invariant the signature advertises. getTemperature is Int (allows negative °C), so the type choice here is also inconsistent within the same class.
- **Evidence:** HPLC.h:62-69 declare UInt pressure/flux accessors; members at HPLC.h:86-88 are `Int pressure_; Int flux_;`. HPLC.cpp:77-95 assign UInt<->Int without range checks.
- **Fix:** Make member types match the public unsigned API (store UInt) or make the API signed to match the members; either way pick one. ABI-safe path: keep signatures, change member declarations to `UInt`, which is a non-ABI-affecting internal change and removes the silent sign reinterpretation.
- **Verifier correction:** The signed-member/unsigned-API mismatch is real, but the round-trip is lossless for all values in [0, INT_MAX]; the unsigned invariant only fails for inputs > INT_MAX (~2.1e9), which are non-physical for pressure (bar) and flux (uL/sec). Hence this is a low-severity type/consistency wart (and inconsistent with the deliberately-signed getTemperature), not a high-severity silent-corruption bug. Fixing by changing members from Int to UInt is non-ABI-affecting (same 4-byte size, unchanged public signatures).
- **Verified:** Evidence is accurate. HPLC.h:62-69 declare UInt getPressure/setPressure(UInt) and UInt getFlux/setFlux(UInt); members at HPLC.h:87-88 are signed `Int pressure_; Int flux_;` (Int=int, UInt=unsigned int per Types.h:64,72). HPLC.cpp:82-94 store UInt into signed Int and read signed Int back as UInt with no range checks, while getTemperature (line 57) is genuinely Int — so the within-class type inconsistency is real. The type mismatch is a genuine sur

### [META-5] Sample::getOrganism / Sample::setOrganism / Sample::getNumber — Doc-comments for organism/number/state accessors are copy-paste wrong ('returns the sample name')
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/METADATA/Sample.h

```cpp
const std::string & getOrganism() const; void setOrganism(const std::string & organism);
```
- **Expectation:** Each accessor's doc describes that accessor's field; getOrganism documents the organism, setOrganism documents setting the organism.
- **Actual:** getOrganism is documented `/// returns the sample name (default: "")` and setOrganism is documented `/// sets the sample name`. The doc-comments were copied from getName/setName and never updated, so generated docs/tooltips mislabel the organism field as the name. (getNumber/setNumber are likewise mismatched: setNumber doc says 'sample ID' while getNumber doc reuses name wording.)
- **Evidence:** Sample.h:70-73: `/// returns the sample name (default: "")` above `getOrganism()` and `/// sets the sample name` above `setOrganism()`. The actual impl (Sample.cpp:94-102) operates on organism_.
- **Fix:** Fix the doc-comments to describe organism/number. Doc-only change, fully ABI-safe.
- **Verifier correction:** Only getOrganism/setOrganism have wrong doc-comments (copied from getName/setName): Sample.h line 70 `/// returns the sample name (default: "")` and line 72 `/// sets the sample name` should describe the organism. The claim's assertion that getNumber/setNumber (and state) are also mismatched is incorrect — Sample.h:75-78 correctly document the sample number, and 'sample ID' there is a deliberate clarifying note, not an error.
- **Verified:** The core surprise is real and verified: in Sample.h, line 70 documents getOrganism() as `/// returns the sample name (default: "")` and line 72 documents setOrganism() as `/// sets the sample name`. These are verbatim copies of the getName/setName comments (lines 65/67), but the implementations (Sample.cpp:94-102) operate on organism_. Doxygen output and IDE tooltips would mislabel the organism accessors as the name field. However, the claim over

### [META-7] IdentificationData::getCurrentProcessingStep — getCurrentProcessingStep is a non-const getter, so it cannot be called on a const IdentificationData
`low` · `const-correctness` · ABI: `none` · src/openms/include/OpenMS/METADATA/ID/IdentificationData.h

```cpp
ProcessingStepRef getCurrentProcessingStep()
```
- **Expectation:** A getter named getCurrentProcessingStep() that just reports the currently-set processing step (a pure read) should be const, callable on a `const IdentificationData&`. All the other accessors in this class (getInputFiles, getScoreTypes, getBestMatchPerObservation, findScoreType, etc.) are const.
- **Actual:** It is declared non-const and the body is simply `return current_step_ref_;` (IdentificationData.cpp:595-599). It performs no mutation but is not marked const, so it is unusable on const instances and breaks the read-only convention of the surrounding getters.
- **Evidence:** Header: `ProcessingStepRef getCurrentProcessingStep();` (IdentificationData.h:449). Impl: `IdentificationData::getCurrentProcessingStep() { return current_step_ref_; }` (IdentificationData.cpp:596-599). Note ProcessingStepRef is an iterator into a std::set, so a const overload could still return it.
- **Fix:** Add a `const` overload `ProcessingStepRef getCurrentProcessingStep() const;` (additive). The ProcessingStepRef is an IteratorWrapper over a set::iterator that can be obtained from a const member without exposing mutability of the set's elements (set elements are const anyway), so this is safe and ABI-additive.
- **Verifier correction:** The factual claim is fully correct. Re-grade severity from any implied higher level to LOW: the missing const is a compile-time-loud convention inconsistency (the getter is unusable on a `const IdentificationData&` while all sibling getters are const), not a source of silent wrong results, data loss, or crashes. Recommended fix — adding `ProcessingStepRef getCurrentProcessingStep() const;` — is safe and source-compatible/additive (the surprise itself carries no ABI implication; abi_impact = none).
- **Verified:** All quoted evidence verified against the actual source. Header line 449 declares `ProcessingStepRef getCurrentProcessingStep();` with no const. Impl (IdentificationData.cpp:595-599) is exactly `return current_step_ref_;` with no mutation. Every sibling accessor IS const: getInputFiles (342), getProcessingSteps (354), getScoreTypes (372), getObservations (378), findScoreType (509), getBestMatchPerObservation (460), getMatchesForObservation (465). 

### [META-15] SpectrumNativeIDParser::isNativeID — isNativeID recognizes the 'frame=' prefix that getRegExFromNativeID handles, but the class doc table omits 'frame=' from both
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/METADATA/SpectrumNativeIDParser.h

```cpp
static bool isNativeID(const std::string& id)
```
- **Expectation:** isNativeID's doc enumerates the 'Recognized prefixes: scan=, scanId=, scanID=, controllerType=, function=, sample=, index=, spectrum=, file='. A caller relying on that list would not expect 'frame=' to be treated as a native ID.
- **Actual:** The implementation also accepts 'frame=' (Bruker TDF), and getRegExFromNativeID maps 'frame=' to the scan regex — but the documented prefix list in isNativeID (and the CV table in the class doc) does not mention 'frame='. The doc and code disagree on the recognized set.
- **Evidence:** Cpp isNativeID lines 22-26 include 'StringUtils::hasPrefix(id, "frame=")'; header doc line lists prefixes without 'frame='. getRegExFromNativeID cpp line 41 also handles 'frame='.
- **Fix:** Add 'frame=' to the documented recognized-prefix list (and the class CV table) so the doc matches behavior. Doc-only fix; no ABI impact.
- **Verifier correction:** The mismatch is genuine but its impact is benign: isNativeID accepting `frame=` is correct behavior for Bruker TDF native IDs (MS:1002818, pattern `frame=<int> scan=<int>`), and getRegExFromNativeID correctly maps it to the scan regex. The defect is purely a documentation-completeness gap — the isNativeID prefix list (header line 122) and the class CV table (lines 30-43) omit `frame=` even though both code paths handle it. No functional/ABI impact; doc-only fix to add `frame=` to the documented prefix list and CV table.
- **Verified:** Evidence verified against actual source. isNativeID (SpectrumNativeIDParser.cpp line 25) includes `StringUtils::hasPrefix(id, "frame=")`, and getRegExFromNativeID (line 41) maps `frame=` to the scan regex. The Doxygen doc for isNativeID (header line 122) enumerates "Recognized prefixes: scan=, scanId=, scanID=, controllerType=, function=, sample=, index=, spectrum=, file=" and omits `frame=`; the class CV table (header lines 30-43) likewise shows

### [META-20] CVTermListInterface::~CVTermListInterface — Non-virtual destructor on a class explicitly documented as a base class
`low` · `ownership-lifetime` · ABI: `breaking` · src/openms/include/OpenMS/METADATA/CVTermListInterface.h

```cpp
~CVTermListInterface(); // non-virtual
```
- **Expectation:** A class whose own documentation says 'It can be used to inherit from instead of CVTermList' and which is publicly derived from (RetentionTime, Peptide, Protein, etc. in TargetedExperimentHelper.h) should have a virtual destructor so deleting a derived object via a base pointer is well-defined.
- **Actual:** The destructor is non-virtual (header explicitly comments '// Destructor (non virtual)'). Derived classes such as RetentionTime add their own data members; deleting any of them through a CVTermListInterface* is undefined behavior (sliced/partial destruction). Note the sibling CVTermList does declare 'virtual ~CVTermList()'.
- **Evidence:** CVTermListInterface.h line 47-48: '// Destructor (non virtual)\n    ~CVTermListInterface();' combined with TargetedExperimentHelper.h: 'class OPENMS_DLLAPI RetentionTime : public CVTermListInterface'
- **Fix:** Make the destructor virtual (additive, but ABI-breaking because it changes the vtable layout) OR document explicitly that derived objects must never be deleted through a CVTermListInterface* and consider protecting the destructor. At minimum the misleading 'use to inherit from' doc and non-virtual dtor should be reconciled.
- **Verifier correction:** CVTermListInterface has a non-virtual destructor that owns/deletes cvt_ptr_, while being documented as a class to 'inherit from' and being publicly derived by RetentionTime/Interpretation/TraMLProduct/Peptide::Modification. Deleting a derived object through a CVTermListInterface* would be UB (the base dtor is statically bound). This is a real API-contract inconsistency, but in practice harmless: the codebase never heap-allocates these types nor deletes them via a base pointer (they are value types in vectors). The claim also overlooks that derived RetentionTime already declares 'virtual ~RetentionTime() = default', making the hierarchy inconsistent rather than purely non-polymorphic, and the base MetaInfoInterface likewise has a non-virtual destructor. Severity is low (latent, never triggered) not high. Making the destructor virtual is ABI-breaking (adds a vtable to a class that currently has none of its own).
- **Verified:** Evidence is fully accurate: CVTermListInterface.h lines 47-48 declare a non-virtual destructor with the explicit comment '// Destructor (non virtual)', the class-level doc (lines 26-27) says it 'can be used to inherit from instead of CVTermList', it is publicly derived by RetentionTime, Interpretation, TraMLProduct, and Peptide::Modification in TargetedExperimentHelper.h, and the sibling CVTermList does declare 'virtual ~CVTermList()' (line 46). 

### [META-21] CVTermListInterface::replaceCVTerms — replaceCVTerms map overload takes a non-const lvalue reference unlike every sibling
`low` · `inconsistent-convention` · ABI: `breaking` · src/openms/include/OpenMS/METADATA/CVTermListInterface.h

```cpp
void replaceCVTerms(std::map<std::string, std::vector<CVTerm> >& cv_terms)
```
- **Expectation:** Given the const-correct siblings replaceCVTerms(const std::vector<CVTerm>&, const std::string&) and replaceCVTerms(const std::map<...>&), this overload should also take its map argument by const reference; it is a read-only input.
- **Actual:** This overload takes std::map<...>& by non-const reference even though the body only reads it (forwards to cvt_ptr_->replaceCVTerms(cv_terms)). It cannot be called with a const map or a temporary, and the non-const signature falsely implies the argument may be mutated. It also collides confusingly with the const-ref map overload added just below it.
- **Evidence:** CVTermListInterface.h line 62: 'void replaceCVTerms(std::map<std::string, std::vector<CVTerm> > & cv_terms);' vs line 74: 'void replaceCVTerms(const std::map<std::string, std::vector<CVTerm> >& cv_term_map);'. Impl line 94-98 only reads cv_terms.
- **Fix:** Deprecate the non-const overload and rely on the const-ref map overload (they are redundant). Removing it is ABI-breaking; making it const is also a signature change. Safest additive step: mark the non-const one [[deprecated]] and route callers to the const overload.
- **Verifier correction:** The non-const-ref map overload (CVTermListInterface.h:62) is read-only — it forwards to CVTermList::replaceCVTerms which takes const std::map<...>& — so its non-const reference parameter is gratuitous and inconsistent with the const-correct sibling at line 74 and all other map-taking methods. It needlessly prevents passing a const map or temporary and misleadingly implies mutation. However, contrary to the original wording: both overloads have existed since the class was created (not 'added just below it'), and the two do NOT cause a compile-time ambiguity — overload resolution is well-defined (non-const lvalue selects line 62, const/temporary selects line 74). Severity is low: behavior is identical and read-only, misuse is caught loudly at compile time, no wrong results/data loss/crash. Safest fix is additive: leave the const overload as canonical and either route the non-const through it or mark it [[deprecated]]; const-ifying or removing the non-const overload would be ABI-breaking.
- **Verified:** Core claim verified against source. CVTermListInterface.h line 62 declares void replaceCVTerms(std::map<std::string, std::vector<CVTerm> >& cv_terms) by non-const lvalue ref, while line 74 declares the same with const std::map<...>&. The impl (CVTermListInterface.cpp:94-98) only reads cv_terms and forwards to cvt_ptr_->replaceCVTerms(cv_terms); the underlying CVTermList only has replaceCVTerms(const std::map<...>&) (CVTermList.h:67), so the argum

### [META-22] CVTermList::hasCVTerm — hasCVTerm documented as 'checks whether the term has a value' but checks accession presence
`low` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/METADATA/CVTermList.h

```cpp
bool hasCVTerm(const std::string& accession) const
```
- **Expectation:** The doc comment 'checks whether the term has a value' would lead a caller to think this tests CVTerm::hasValue() semantics; given the parameter is an accession, a reader needs the doc to clarify it tests existence.
- **Actual:** The method returns whether any term with the given accession exists in the map (cv_terms_.contains(accession)); it has nothing to do with values. The same wrong doc comment is copy-pasted in CVTermListInterface.h ('checks whether the term has a value'). The name 'hasCVTerm' is fine; the documentation actively misleads.
- **Evidence:** CVTermList.h line 94-95: '/// checks whether the term has a value\n    bool hasCVTerm(const std::string& accession) const;' and CVTermList.cpp: 'return cv_terms_.contains(accession);'. Same comment at CVTermListInterface.h line 85-86.
- **Fix:** Fix the doc comment to 'checks whether a CV term with the given accession is present' in both headers. Doc-only, fully ABI compatible.
- **Verifier correction:** The defect is a wrong/copy-pasted doc comment, not a misleading name. The comment '/// checks whether the term has a value' appears at CVTermList.h:94 and CVTermListInterface.h:85 above 'bool hasCVTerm(const std::string& accession) const', but the implementation only does 'cv_terms_.contains(accession)' — an existence check by accession with no value semantics. Fix the comment in both headers to e.g. 'checks whether a CV term with the given accession is present'. Doc-only, fully ABI-compatible. Severity is low because the function name and the 'accession' parameter make the correct behavior obvious from the signature regardless of the erroneous comment.
- **Verified:** Evidence is confirmed verbatim. CVTermList.h:94-95 has '/// checks whether the term has a value' above 'bool hasCVTerm(const std::string& accession) const;', the implementation in CVTermList.cpp:71-74 is just 'return cv_terms_.contains(accession);' (pure accession-existence check, no value semantics), and the identical wrong comment is copy-pasted at CVTermListInterface.h:85-86. The doc is objectively wrong: it claims a value-check while the meth

### [META-23] CVTermList::getCVTerms — getCVTerms documented as 'returns the accession string of the term'
`low` · `incorrect-doc-comment` · ABI: `none` · src/openms/include/OpenMS/METADATA/CVTermList.h

```cpp
const std::map<std::string, std::vector<CVTerm> >& getCVTerms() const
```
- **Expectation:** Doc should describe that this returns the full map of accession -> vector<CVTerm>.
- **Actual:** The doc comment is a copy-paste leftover from CVTerm::getAccession: '/// returns the accession string of the term'. It returns the whole CV-term map, not a single accession string. The identical wrong comment appears in CVTermListInterface.h.
- **Evidence:** CVTermList.h line 72-73: '/// returns the accession string of the term\n    const std::map<std::string, std::vector<CVTerm> >& getCVTerms() const;' and CVTermListInterface.h line 79-80 (identical).
- **Fix:** Correct the doc comment in both headers to describe the returned map. Doc-only, ABI compatible.
- **Verifier correction:** The defect is a stale, incorrect Doxygen comment, not a misleading method name. The comment '/// returns the accession string of the term' on CVTermList::getCVTerms (CVTermList.h:72) and CVTermListInterface::getCVTerms (CVTermListInterface.h:79) is a copy-paste leftover from CVTerm::getAccession (CVTerm.h:105). Both methods return 'const std::map<std::string, std::vector<CVTerm> >&' (the full accession->vector<CVTerm> map), which the signature makes plain. Correct the comment in both headers to e.g. '/// returns the map of CV-term accessions to their CVTerm vectors'. Doc-only, ABI- and source-compatible. Severity is low because the correct method name and the visible return type prevent any silent misuse.
- **Verified:** Evidence verified exactly. CVTermList.h:72-73 and CVTermListInterface.h:79-80 both carry the comment '/// returns the accession string of the term' directly above 'const std::map<std::string, std::vector<CVTerm> >& getCVTerms() const'. The comment is a confirmed copy-paste leftover from CVTerm.h:105-106, where it correctly describes getAccession() returning 'const std::string&'. So the doc is objectively wrong: it claims a single accession string

### [META-25] CVTerm::Unit::~Unit — Plain value struct CVTerm::Unit has a virtual destructor, giving it a surprising vtable
`low` · `other` · ABI: `breaking` · src/openms/include/OpenMS/METADATA/CVTerm.h

```cpp
virtual ~Unit()
```
- **Expectation:** A small aggregate-like value type holding three std::strings (accession, name, cv_ref), used as a by-value member and compared with operator==, is expected to be a trivial value type with no vtable/polymorphism.
- **Actual:** Unit declares 'virtual ~Unit()'. This adds a vtable pointer, makes the type polymorphic and larger, and signals an inheritance intent that does not exist (Unit is a leaf data holder). It also makes the otherwise-trivial value semantics non-obvious and is inconsistent with the rest of the struct being default-able.
- **Evidence:** CVTerm.h line 50-53: '/// Destructor\n      virtual ~Unit()\n      {\n      }'
- **Fix:** Make the destructor non-virtual (ideally = default and non-virtual) since Unit is not a base class. Changing virtual->non-virtual is ABI-breaking (vtable removal), so guard behind a deliberate ABI-break window; otherwise document that the virtual is vestigial.
- **Verifier correction:** CVTerm::Unit (CVTerm.h:51) declares a `virtual ~Unit()` despite being a leaf, aggregate-like value type with no subclasses and no polymorphic use anywhere in the codebase. The virtual is vestigial: it adds a vtable pointer (larger sizeof, non-trivial/non-standard-layout) and misleadingly implies inheritance intent. It does not cause incorrect behavior — the type still works correctly by value — so this is a low-severity design/space wart, not a correctness or misuse hazard. Removing the virtual (ideally `~Unit() = default;` non-virtual) is ABI-breaking due to vtable removal and should be deferred to an ABI-break window; otherwise the virtual should be annotated as intentional/vestigial.
- **Verified:** Evidence confirmed verbatim: CVTerm.h lines 50-53 declare `virtual ~Unit() {}`. CVTerm::Unit is a leaf value struct holding three std::strings (accession, name, cv_ref); it is used by-value as the CVTerm member `unit_`, returned by value/ref, compared via operator==, and hashed via std::hash specialization. A repo-wide grep finds NO class deriving from CVTerm::Unit and no polymorphic use (no Unit* base pointers, no dynamic_cast). So the `virtual`

### [META-28] AcquisitionInfo — Container exposes add/insert/resize but no element removal (no clear/erase/pop_back)
`low` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/METADATA/AcquisitionInfo.h

```cpp
class AcquisitionInfo : private std::vector<Acquisition>, public MetaInfoInterface
```
- **Expectation:** A vector-like class that publishes push_back, insert, resize, operator[], begin/end and back is expected to also expose the matching removal operations (clear, erase, pop_back) so callers can manage the collection symmetrically.
- **Actual:** Only a curated subset of std::vector is re-exported via 'using ContainerType::...': operator[], begin, end, size, push_back, empty, back, insert, resize. There is no clear, erase, pop_back, front, reserve, or at. A caller can add Acquisitions but has no in-place way to remove a single one (resize() can only truncate from the end), and the base is privately inherited so the missing members are inaccessible.
- **Evidence:** AcquisitionInfo.h lines 60-71: 'using ContainerType::operator[]; using ContainerType::begin; using ContainerType::end; using ContainerType::size; using ContainerType::push_back; using ContainerType::empty; using ContainerType::back; using ContainerType::insert; using ContainerType::resize;' — no clear/erase/pop_back exported.
- **Fix:** Additively export the missing symmetric members ('using ContainerType::clear; using ContainerType::erase; using ContainerType::pop_back;' and likely front/at/reserve) to complete the container surface. Purely additive, source- and ABI-compatible.
- **Verifier correction:** The asymmetry and missing-members evidence are accurate and confirmed. Adjustments: (1) Severity is low, not medium/high — the missing operations cause a compile-time error for any caller attempting removal (loud and immediate), with no path to silent wrong results, data loss, or crashes, and no current caller needs removal. (2) abi_impact is none (not merely source-compatible): adding 'using ContainerType::clear/erase/pop_back/...' declarations is purely additive and exports no new symbols from the shared library. (3) Strengthened justification that this is a genuine surprise rather than a convention: the sibling classes MSSpectrum (MSSpectrum.h:140,146) and MSChromatogram (MSChromatogram.h:86,92) use the identical idiom and DO export pop_back/erase/front/reserve, so AcquisitionInfo's omission is an inconsistency within the same codebase, not an intentional curated design.
- **Verified:** Evidence verified exactly: AcquisitionInfo.h lines 60-68 re-export only operator[], begin, end, size, push_back, empty, back, insert, resize from the privately-inherited std::vector<Acquisition>; clear, erase, pop_back, front, reserve, at are all absent (grep confirmed exit 1). Because the base is privately inherited, the omitted members are genuinely inaccessible and there is no std::vector& conversion to fall back on; resize() can only truncate

### [META-32] ExperimentalDesign::getSample — getSample is a read-only lookup but is declared non-const
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/METADATA/ExperimentalDesign.h

```cpp
unsigned getSample(unsigned fraction_group, unsigned label = 1)
```
- **Expectation:** A pure lookup named getSample that only reads msfile_section_ should be const and callable on a const ExperimentalDesign, like every other getter in the class.
- **Actual:** The method is declared non-const even though its body only performs a std::find_if read and returns a member value; it cannot be called on a const ExperimentalDesign, unlike the sibling getNumberOf*/getPathLabelTo* accessors which are all const.
- **Evidence:** ExperimentalDesign.h:270 `unsigned getSample(unsigned fraction_group, unsigned label = 1);` (no trailing const) vs. ExperimentalDesign.cpp:599 body only reads msfile_section_.
- **Fix:** Add const qualifier. Since this changes the mangled name, ideally add a const overload (or change to const if no callers rely on the non-const binding) and deprecate the non-const form.
- **Verifier correction:** The claim is factually accurate: getSample (ExperimentalDesign.h:270) is a pure read-only lookup (cpp:599-606 only does std::find_if and returns ->sample) yet is non-const, unlike all sibling getters. The only correction is grading: the impact is a loud compile error on const instances, never a silent miscalculation, so severity is low (mild surprise / invites mild friction), not high. The single direct caller (pyOpenMS bind_metadata.cpp:186) does not depend on the non-const binding, so simply qualifying it const is source-compatible for all callers, but it changes the mangled name and therefore breaks binary ABI; a deprecated non-const overload would preserve ABI.
- **Verified:** Verified directly. ExperimentalDesign.h:270 declares `unsigned getSample(unsigned fraction_group, unsigned label = 1);` with no trailing const, and the definition (ExperimentalDesign.cpp:599-606) only performs a std::find_if read over msfile_section_ and returns ->sample by value — it mutates nothing. Every sibling read-only accessor in the class is const: getNumberOfSamples/Fractions/Labels/MSFiles/FractionGroups (h:254-267), all getPathLabelTo*

### [META-33] ExperimentalDesign::getNumberOfFractions — getNumberOfFractions documented as 'highest fraction index' but returns count of distinct fraction ids
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/METADATA/ExperimentalDesign.h

```cpp
unsigned getNumberOfFractions() const
```
- **Expectation:** Doc '(= highest fraction index)' implies the return is the maximum fraction value seen.
- **Actual:** The implementation collects fraction ids into a std::set and returns the set's size (count of distinct fractions). For a design using fraction ids {1,3} the doc predicts 3 but the function returns 2.
- **Evidence:** Header ExperimentalDesign.h:256 `// @return the number of fractions (= highest fraction index)`; impl ExperimentalDesign.cpp:542-547 `auto fs = set<Size>(); ... fs.insert(row.fraction); ... return fs.size();`
- **Fix:** Update the doc to 'number of distinct fraction identifiers'. Documentation-only; no ABI impact.
- **Verifier correction:** The function returns the count of distinct fraction identifiers (set size), which equals the "highest fraction index" only when fraction ids are contiguous and 1-based. The doc comment at ExperimentalDesign.h:256 should read "number of distinct fraction identifiers". Severity is low (not the implied high): divergence requires non-contiguous fraction ids, an edge case the implementation itself flags as ill-defined in its TODO, and it produces no crash or silent data corruption in normal contiguous designs.
- **Verified:** Evidence is confirmed verbatim. Header ExperimentalDesign.h:256 documents the return as "the number of fractions (= highest fraction index)", but the impl ExperimentalDesign.cpp:542-547 builds a set<Size> of distinct row.fraction values and returns fs.size() — a count of distinct fraction ids, not the maximum index. The commented-out prior implementation (lines 534-540) used std::max_element(...)->fraction (matching the doc), so the doc went stal

### [META-36] IdentifierMSRunMapper::getPrimaryMSRunPath — getPrimaryMSRunPath returns an empty string for unknown identifier or out-of-range merge index
`low` · `silent-failure` · ABI: `none` · src/openms/source/METADATA/IdentifierMSRunMapper.cpp

```cpp
std::string getPrimaryMSRunPath(const PeptideIdentification& pepid) const
```
- **Expectation:** A getter named getPrimaryMSRunPath should return the run path, or clearly signal that none is available (throw / optional), so the caller cannot accidentally proceed with an empty path.
- **Actual:** Three distinct failure modes (identifier not in map, empty path list, merge_index >= size) all silently return an empty std::string(), which is indistinguishable from a legitimately empty path and easy to forward into a USI or file lookup unnoticed.
- **Evidence:** IdentifierMSRunMapper.cpp:55-58, 61-64, 74-77 each `return std::string();` for not-found / empty / invalid-index cases.
- **Fix:** Document the empty-string-means-not-found contract in the header, or add a tryGetPrimaryMSRunPath(...) bool overload (consistent with the existing tryGetIdentifier pattern) so callers can detect failure. Additive; no ABI break.
- **Verified:** Evidence confirmed verbatim. IdentifierMSRunMapper.cpp lines 55-58, 61-64, and 74-77 each `return std::string();` for the three distinct failure modes (identifier not in map, empty path list, merge_index out of range), all indistinguishable from a legitimately empty path. The header (line 62) documents the method only as "Get the primary MS run path ... (using id_merge_index metadata)" and says NOTHING about an empty-string-means-not-found contra

### [META-41] PeptideEvidence::setStart — setStart() documented as setting the position of the *last* AA (copy-paste error from setEnd)
`low` · `incorrect-documentation` · ABI: `none` · src/openms/include/OpenMS/METADATA/PeptideEvidence.h

```cpp
void setStart(const Int a)
```
- **Expectation:** setStart() sets the position of the FIRST amino acid of the peptide in the protein; the doc should describe the start/N-terminal-side position.
- **Actual:** The header comment reads: "set the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). ... N-terminal positions must be marked with N_TERMINAL_AA" — identical 'last AA' wording to setEnd. The implementation simply stores `start_ = a;`, i.e. it really is the start position, so the doc actively contradicts the field it sets.
- **Evidence:** Header L79: "/// set the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set to UNKNOWN_POSITION. N-terminal positions must be marked with N_TERMINAL_AA"; setEnd L85 has the same "last AA" wording. PeptideEvidence.cpp L98-101: setStart just assigns start_.
- **Fix:** Fix the doc comment to say "position of the FIRST AA of the peptide" (start), purely a comment change. abi_impact none.
- **Verifier correction:** The defect is real but is a documentation bug, not a misleading name. The setStart() symbol name and implementation (start_ = a) are correct; only the Doxygen comment on PeptideEvidence.h L79 is wrong, having been copy-pasted from setEnd() ("set the position of the last AA..."). It should read "set the position of the FIRST AA of the peptide in protein coordinates (starting at 0 for the N-terminus)...". Category corrected to incorrect-documentation and severity to low, since the correct function name, getStart() getter, constructor `start` parameter, and N_TERMINAL_AA marker all steer a competent developer to the right meaning, making silent misuse unlikely. Comment-only change; ABI impact none.
- **Verified:** Verified against the actual source. Header L79 documents setStart() as "set the position of the last AA of the peptide in protein coordinates..." — verbatim the same "last AA" wording as setEnd() on L85. The implementation (PeptideEvidence.cpp L98-101) does `start_ = a;`, i.e. it stores the START position. So the doc comment genuinely contradicts the field/operation: it should say FIRST AA, not last AA. This is confirmed as a copy-paste error fro

### [META-42] ProteinHit::setSequence / setAccession — setSequence()/setAccession() silently trim() the input, so the stored value differs from what was passed
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/METADATA/ProteinHit.h

```cpp
void setSequence(const std::string& sequence); void setAccession(const std::string& accession)
```
- **Expectation:** A plain setter stores the value verbatim; getSequence()/getAccession() should return exactly what was set.
- **Actual:** setSequence() and setAccession() (and the values constructor) call StringUtils::trim() on the stored string. Passing " PEPTIDE " stores "PEPTIDE"; leading/trailing whitespace is silently removed, which can surprise round-trip or equality logic. This trimming is undocumented in the header.
- **Evidence:** ProteinHit.cpp setSequence L115-119 `sequence_ = sequence; StringUtils::trim(sequence_);`, setAccession L136-140, ctor L31-39 trims accession and sequence. Header L156-161 documents none of this.
- **Fix:** Document the trimming in the header comments for setSequence/setAccession (and the constructor). abi_impact none (doc only).
- **Verifier correction:** Claim is accurate as stated. Two refinements: (1) the trimming is long-standing pre-existing behavior (old String::trim()), not introduced by the std::string migration; (2) the rvalue overload setSequence(std::string&&) at L122-126 also trims, so BOTH setSequence overloads (not just the const-ref one) are affected. Severity is low rather than medium: trimming stray whitespace from accessions/sequences is benign normalization, only matters for whitespace-padded inputs, and is observable/recoverable; the only defect is the undocumented header. Recommendation (document trimming on setSequence/setAccession and the values constructor) is appropriate; abi_impact none.
- **Verified:** Independently verified in the actual code. ProteinHit.cpp confirms: setSequence(const std::string&) L115-119 does `sequence_ = sequence; StringUtils::trim(sequence_);`; the rvalue overload L122-126 does the same after move; setAccession L136-140 trims accession_; and the values constructor L31-39 initializes accession_(StringUtils::trim(accession)) and sequence_(StringUtils::trim(sequence)). StringUtils::trim (StringUtils.h L555-568) mutates in p

### [META-43] ProteinHit::setModifications — setModifications takes a non-const lvalue reference, so const sets and temporaries cannot be passed
`low` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/METADATA/ProteinHit.h

```cpp
void setModifications(std::set<std::pair<Size, ResidueModification> >& mods)
```
- **Expectation:** A setter that only copies its argument should take a const reference (or by value + move), allowing `hit.setModifications(myConstMods)` or `hit.setModifications({...})`.
- **Actual:** The parameter is a mutable lvalue reference `std::set<...>& mods` even though the body just copies (`modifications_ = mods;`) and never mutates the argument. Callers holding a const set or a temporary fail to compile; readers may wrongly assume the input set is modified (in/out parameter).
- **Evidence:** Header L173; ProteinHit.cpp L153-156 `modifications_ = mods;` (no mutation of mods).
- **Fix:** Add an overload taking `const std::set<...>&` (additive, source-compatible); ideally deprecate the non-const-ref signature. abi_impact source-compatible.
- **Verifier correction:** setModifications takes a non-const lvalue reference (ProteinHit.h L173) although the body only copies (ProteinHit.cpp L155 `modifications_ = mods;`). This prevents passing const sets or temporaries/braced-init-lists and falsely suggests in/out semantics. It is inconsistent with the codebase's own ModificationDefinitionsSet::setModifications(const std::set<...>&) and all sibling ProteinHit const-ref setters. The defect is real but its impact is a loud compile-time failure plus a readability surprise (LOW severity), not silently wrong behavior. Fix: add an overload taking `const std::set<std::pair<Size, ResidueModification> >&` (source-compatible, additive).
- **Verified:** Evidence verified independently. Header L173 declares `void setModifications(std::set<std::pair<Size, ResidueModification> >& mods)` — a non-const lvalue reference — and ProteinHit.cpp L153-156 has a pure-copy body `modifications_ = mods;` that never mutates the argument. Because the parameter is a mutable lvalue reference, a const set cannot be passed and `hit.setModifications({...})` (temporary/braced-init) will not compile, even though copy-on

### [META-45] PeptideIdentification::empty — empty() ignores RT, MZ, spectrum_reference and hits' annotations — a spectrum-anchored ID with metadata reports empty()==true
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/METADATA/PeptideIdentification.h

```cpp
bool empty() const
```
- **Expectation:** For an object that anchors a spectrum match, empty() should report whether the object carries no meaningful content; a PeptideIdentification with an RT/MZ/spectrum_reference set (i.e. a real spectrum it points at) is not 'empty'.
- **Actual:** empty() only checks id_, hits_, significance threshold, score_type_, higher_score_better_==true (default) and base name. It does NOT consider mz_, rt_, or the spectrum_reference meta value. A PeptideIdentification that has set RT/MZ and a spectrum reference but no hits returns empty()==true, which can cause callers to silently discard spectrum-anchored records.
- **Evidence:** PeptideIdentification.cpp L225-233 (no mz_/rt_/spectrum_reference checks). setRT/setMZ store mz_/rt_ (L55-83).
- **Fix:** Document precisely what empty() considers, or add the RT/MZ/spectrum_reference checks behind a clearly named method (e.g. hasContent()). abi_impact none for doc; source-compatible if a new method is added.
- **Verifier correction:** empty() accurately omits mz_, rt_, and spectrum_reference from its check, but this reflects an intentional (though under-documented) convention that empty() means "no identification result content" — RT/MZ/spectrum_reference are spectrum-anchor metadata, not content. The real caller (FeatureFinderIdentificationAlgorithm) relies on exactly this behavior to drop hit-less records, so there is no demonstrated "silent discard of spectrum-anchored records" harm. The genuine, mild surprise is the asymmetry with operator== (which compares mz_/rt_) and the imprecise one-line doc. Remedy: document precisely what empty() considers (and that RT/MZ/spectrum_reference are deliberately ignored). abi_impact none for a doc fix; adding a clearly-named hasContent()/hasAnchor() method would be source-compatible.
- **Verified:** Code claim is factually accurate. PeptideIdentification::empty() (PeptideIdentification.cpp L225-233) checks only id_, hits_, getSignificanceThreshold()==0.0, score_type_, higher_score_better_==true, and getBaseName(). It does NOT check mz_, rt_, or the spectrum_reference meta value. setRT/setMZ store mz_/rt_ (L60-78) and setSpectrumReference stores a meta value (L168-171), confirmed. So an object with RT/MZ/spectrum_reference set but no hits/id/

### [META-46] PeptideIdentification::hasMZ — hasMZ() is documented as "shortcut for isnan(getRT())" (copy-paste from hasRT)
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/METADATA/PeptideIdentification.h

```cpp
bool hasMZ() const
```
- **Expectation:** hasMZ()'s doc should describe the MZ check (isnan(getMZ())).
- **Actual:** The header comment for hasMZ() reads "/// shortcut for isnan(getRT())" — it references RT, not MZ. The implementation correctly returns `!std::isnan(mz_)`, so only the doc is wrong, but a reader trusting the comment could believe hasMZ() reflects RT presence.
- **Evidence:** Header L104-105 (both hasRT and hasMZ carry the comment "shortcut for isnan(getRT())"). PeptideIdentification.cpp L80-83 hasMZ() checks mz_.
- **Fix:** Fix the comment to "shortcut for !isnan(getMZ())" (and note both comments say isnan but the methods return the negation). abi_impact none.
- **Verifier correction:** The hasMZ() Doxygen brief at PeptideIdentification.h L104 reads "/// shortcut for isnan(getRT())", a copy-paste from hasRT() (L97) that references RT instead of MZ; the implementation (PeptideIdentification.cpp L80-83) correctly checks mz_. Fix the comment to reference getMZ() (e.g. "/// shortcut for !isnan(getMZ())"). The separate claim that the comment is wrong because the method returns the negation of isnan is not a real defect — it is consistent informal shorthand used for both hasRT() and hasMZ().
- **Verified:** Evidence confirmed verbatim. Header L104 documents hasMZ() with "/// shortcut for isnan(getRT())" — a copy-paste of the hasRT() brief on L97, referencing RT instead of MZ. The implementation (cpp L80-83) correctly returns !std::isnan(mz_), so only the doc comment is wrong. This is a genuine, provable copy-paste documentation defect, not a code bug. However, severity is LOW, not anything higher: the method name hasMZ() is self-describing and the i

### [META-51] IdentifiedMolecule::getIdentifiedPeptideRef / getIdentifiedCompoundRef / getIdentifiedOligoRef — get*Ref() accessors throw when the variant holds a different molecule type
`low` · `surprising-throw` · ABI: `none` · src/openms/source/METADATA/ID/IdentifiedMolecule.cpp

```cpp
IdentifiedPeptideRef getIdentifiedPeptideRef() const
```
- **Expectation:** A method named get...Ref() reads a stored reference; a caller may treat it like a plain getter and not wrap it in try/catch.
- **Actual:** Each accessor throws Exception::IllegalArgument when the active variant alternative is not the requested type. Because getMoleculeType() must be checked first, the natural call site `mol.getIdentifiedPeptideRef()` on a compound/oligo throws.
- **Evidence:** IdentifiedMolecule.cpp lines 43-78: each accessor ends with `throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg);` ("matched molecule is not a peptide" etc.).
- **Fix:** Document the throwing precondition (must match getMoleculeType()) in the header, and consider adding non-throwing `tryGetIdentified*Ref()` variants returning std::optional. Additive, ABI-safe.
- **Verified:** Evidence verified exactly: IdentifiedMolecule.cpp lines 43-78 — each of getIdentifiedPeptideRef/Compound/Oligo uses std::get_if and, on a mismatched active alternative, throws Exception::IllegalArgument ("matched molecule is not a peptide/compound/oligonucleotide"). The header (lines 40-44) documents none of this. All call sites checked (IdentificationDataConverter.cpp:382-397, IdentificationData.cpp:882-888/1331-1378, FalseDiscoveryRate.cpp:693-


---

# CHEMISTRY

**Counts:** 1 high · 36 medium · 32 low

### [CHEM-24] IMSIsotopeDistribution::size — size() caps the reported peak count at the global static SIZE rather than reporting the actual peak count
`high` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h

```cpp
size_type size() const { return std::min(peaks_.size(), SIZE); }
```
- **Expectation:** size() on a container-like class returns the number of stored elements, consistent with empty() returning peaks_.empty().
- **Actual:** size() returns min(peaks_.size(), SIZE) where SIZE is a mutable static class member shared across all instances. An instance holding 12 peaks reports size()==SIZE (e.g. 10), so iterating [0,size()) skips real peaks; the value also changes globally if any code reassigns the static SIZE. empty() (peaks_.empty()) and size() can therefore disagree about how many peaks are accessible.
- **Evidence:** Header line 162: `size_type size() const { return std::min(peaks_.size(), SIZE); }` with `static size_type SIZE;` (line 125) and `bool empty() const { return peaks_.empty(); }` (line 288). Header note at 156-159 even admits 'Size is not smaller than predefined SIZE' which is the opposite (it is capped *at most* SIZE).
- **Fix:** ABI-safe: fix the doc comment (it currently claims the wrong inequality) and clearly state size() is clamped to the static SIZE truncation length, distinct from peaks_.size(). Consider adding a separate peakCount()/rawSize() accessor (additive) for the true count.
- **Verifier correction:** size() does not merely cap at "e.g. 10"; the static SIZE is default-initialized to 0 and is never assigned anywhere in the repo, so size() returns 0 for any distribution constructed from peaks (until a fold operation grows peaks_ to SIZE or a client mutates the static). This yields silently empty results from getMasses()/getAbundances()/operator<</pyOpenMS len() despite empty()==false. The doc comment (lines 156-159) is also literally backwards — size() is at most SIZE, not "not smaller than" SIZE. Recommended fix is additive/source-compatible: correct the doc and add a peakCount()/rawSize() accessor returning peaks_.size().
- **Verified:** Independently verified against the source. Header line 162 is exactly `size_type size() const { return std::min(peaks_.size(), SIZE); }`, line 125 declares `static size_type SIZE;` (mutable, non-const), and line 288 `empty()` returns `peaks_.empty()`. So size() is clamped to a mutable static shared across ALL instances and is decoupled from the actual stored-peak count, while empty() reflects the true count — they can and do disagree. The doc at 

### [CHEM-1] EmpiricalFormula::getMonoWeight / getAverageWeight / getLightestIsotopeWeight — Weight getters silently add/subtract PROTON mass per charge ('includes proton charges'), surprising for anions and for charge-carrier choice
`medium` · `unit-or-index` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/EmpiricalFormula.h

```cpp
double getMonoWeight() const; double getAverageWeight() const; double getLightestIsotopeWeight() const;
```
- **Expectation:** getMonoWeight() of an EmpiricalFormula returns the (neutral) monoisotopic mass of the summed elements; if a charge is involved a caller would expect either no charge correction or an electron-mass correction (the physically correct charge carrier for de/protonation is a proton mass minus an electron).
- **Actual:** The implementation computes 'weight = Constants::PROTON_MASS_U * charge_' and adds element masses (EmpiricalFormula.cpp:49,71,59). It uses the full PROTON mass (not proton-minus-electron) per unit charge, and for a NEGATIVE charge it SUBTRACTS proton masses. So a formula carrying charge -1 returns a mass ~1.007 Da LOWER than the neutral composition, which is the opposite of the physical [M-H]- intuition and is not the neutral mass either. The behavior is so unsuitable that AdductInfo refuses any charged EmpiricalFormula and tracks charge itself.
- **Evidence:** EmpiricalFormula.cpp:49 `double weight = Constants::PROTON_MASS_U * charge_;` ... ; AdductInfo.cpp:28-31 throws InvalidParameter: "EmpiricalFormula must not have a charge (...), since the internal weight computation of EF is unsuitable for adducts."
- **Fix:** Keep the ABI (this convention is relied upon), but sharpen the doc: state explicitly that the returned mass is shifted by +PROTON_MASS_U*charge (proton mass, not electron-corrected), including being LOWERED for negative charges, and is NOT the neutral mass nor an ion m/z. Optionally add a clearly named additive helper (e.g. getNeutralMonoWeight()) that ignores charge_. abi_impact: doc/additive only.
- **Verifier correction:** The behavior is real and as described (full proton mass per charge, subtracted for negative charge), but it is partially documented in the header as "(includes proton charges)" — the issue is that this note is incomplete/misleading, not absent. The surprise is confined to charged EmpiricalFormula objects (charge_ != 0); the common neutral case (charge_=0) correctly returns the neutral monoisotopic mass. Hence severity is medium (silent only on an edge-case construction, with a partial doc warning and a loudly-throwing primary consumer in AdductInfo), not high. Fix is doc-only (ABI none) plus an optional additive getNeutralMonoWeight() helper (source-compatible).
- **Verified:** Code reading is accurate. EmpiricalFormula.cpp:49,59,71 all init `weight = Constants::PROTON_MASS_U * charge_` then add element masses; PROTON_MASS_U=1.0072764667710 (full proton, not proton-minus-electron, and ELECTRON_MASS_U exists separately so the omission is real). For negative charge_ the term is negative, so a charge -1 formula returns ~1.007 Da BELOW the neutral composition — neither the neutral mass nor the [M-H]- intuition. AdductInfo.c

### [CHEM-4] EmpiricalFormula::operator< — operator< orders primarily by NUMBER OF DISTINCT ELEMENTS, not by mass/composition — surprising for a chemical-formula 'less than'
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/EmpiricalFormula.h

```cpp
bool operator<(const EmpiricalFormula& rhs) const;
```
- **Expectation:** For a mass-spec formula type, '<' is most naturally read as a mass ordering or at least a lexicographic composition ordering; callers using std::set<EmpiricalFormula>/std::map keys expect a meaningful, composition-based order.
- **Actual:** operator< first compares formula_.size() (count of distinct elements), then charge_, then the map itself (EmpiricalFormula.cpp:660-673). So C100 (1 distinct element) sorts BEFORE H2O (2 distinct elements) regardless of mass, and the final tiebreak compares (const Element*) pointer-keyed maps, whose order depends on allocation addresses and is not reproducible across runs.
- **Evidence:** EmpiricalFormula.cpp:662-672 `if (formula_.size() != rhs.formula_.size()) return formula_.size() < rhs.formula_.size(); ... return formula_ < rhs.formula_;` (map keyed on const Element*).
- **Fix:** Document explicitly that operator< is an arbitrary strict-weak ordering for container use (NOT a mass or canonical-composition order) and that the tiebreak is pointer-dependent / non-portable. abi_impact: doc-only; a value-based reordering would be behavior-breaking.
- **Verifier correction:** operator< is confirmed to order primarily by number of distinct elements (formula_.size()), then charge_, then a pointer-keyed map comparison whose order is non-reproducible across runs (heap-allocated Element* + std::less<const Element*>). This is a genuine misleading-name surprise and the tiebreak is genuinely non-portable. However, severity is medium, not high: every current use in the codebase (TheoreticalSpectrumGenerator's std::set<EmpiricalFormula>/std::map<EmpiricalFormula,...>) uses the container only for deduplication/lookup, where the ordering value does not influence output correctness — so no silently-wrong results today. The defect is a latent reproducibility/misuse hazard. Recommended fix is documentation-only (clarify it is an arbitrary strict-weak ordering for container use, NOT mass/composition order, and that the tiebreak is pointer-dependent/non-deterministic across runs); abi_impact for that fix is none. A value-based reordering would be behavior-breaking but is not required to address the surprise.
- **Verified:** Code matches the claim verbatim. EmpiricalFormula.cpp:660-673: operator< compares formula_.size() (count of DISTINCT elements) first, then charge_, then formula_ < rhs.formula_. So C100 (1 distinct element) sorts before H2O (2 distinct elements) regardless of mass — genuinely surprising for a mass-spec formula '<', which a competent dev would read as mass or composition ordering. Not a chemistry convention (neither Hill notation nor mass-ordering

### [CHEM-5] EmpiricalFormula::toString / toMap / operator<< — toString()/toMap() silently drop the charge, while operator<< prints it — round-trip via toString loses charge with no warning
`medium` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/EmpiricalFormula.h

```cpp
std::string toString() const; std::map<std::string,int> toMap() const;
```
- **Expectation:** Given the string constructor parses a trailing charge (e.g. 'H4C-1-'), a caller would expect toString() to be its inverse and preserve charge, especially since operator<< DOES emit the charge.
- **Actual:** toString() and toMap() omit the charge entirely (documented '(charges are not included)'), so EmpiricalFormula(ef.toString()) silently loses a non-zero charge, whereas streaming the same object via operator<< appends '+'/'-N' (EmpiricalFormula.cpp:410-437). The two serializers disagree, and the lossy one looks like the natural round-trip.
- **Evidence:** EmpiricalFormula.h:215 `/// returns the formula as a string (charges are not included)`; EmpiricalFormula.cpp:248-257 toString() never reads charge_; vs operator<< EmpiricalFormula.cpp:415-436 prints charge.
- **Fix:** Document that toString() is NOT a faithful round-trip when isCharged(), and point callers to operator<< (or add an optional bool include_charge param / a toStringWithCharge()). abi_impact: doc or additive overload.
- **Verifier correction:** toString()/toMap() drop the charge (documented inline at EmpiricalFormula.h:215,218) while operator<< (EmpiricalFormula.cpp:410-436) prints it, so EmpiricalFormula(ef.toString()) silently loses a non-zero charge. However, operator<< is NOT a clean inverse either: it omits element counts of 1 and mishandles negative counts (if (it.second > 1), cpp:408), so toString() — which always emits the count (cpp:254) — is actually the faithful ELEMENT serializer; the only field on which the two serializers disagree is the charge. Fix via doc clarification of the non-round-trip behavior or an additive include_charge/toStringWithCharge() overload (source-compatible).
- **Verified:** Code confirms the asymmetry. toString() (EmpiricalFormula.cpp:248-257) delegates to toMap() (259-267); neither reads charge_. operator<< (397-438) prints the charge as '+'/'-N' (410-436). The string ctor (cpp:36 -> parseFormula_) parses a trailing charge. So EmpiricalFormula(ef.toString()) silently drops a non-zero charge while streaming preserves it — a genuine asymmetric-serializer surprise, not a domain convention or C++ idiom. Two corrections

### [CHEM-7] EnzymaticDigestion::digestUnmodified — Pair overload documented as (start, end) but actually returns (start, length)
`medium` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/EnzymaticDigestion.h

```cpp
Size digestUnmodified(const std::string_view& sequence, std::vector<std::pair<Size, Size>>& output, Size min_length = 1, Size max_length = 0) const
```
- **Expectation:** The doc says "@param[out] output Digestion products as vector of pairs of start and end positions", so a caller expects output[i].second to be the past-the-end index and would compute the peptide length as second - first or extract via substr(first, second - first).
- **Actual:** Every code path emits (start, length), not (start, end). In the unspecific path it does `output.emplace_back(i, j - i)`; the normal path delegates to digestAfterTokenize_ which does `output.emplace_back(fragment_positions[i-1], l)` with l being a length; the semi path does `output.emplace_back(p.first, p.second - p.first)`. A caller reading second as an end index gets a number that is actually the length, silently producing wrong sub-sequences.
- **Evidence:** EnzymaticDigestion.cpp lines 538 `output.emplace_back(i, j - i);`, 341 `output.emplace_back(fragment_positions[i - 1], l);`, 559 `output.emplace_back(p.first, p.second - p.first);`. The header doc says "pairs of start and end positions" (EnzymaticDigestion.h line 115). EnzymaticDigestion_test.cpp line 334 confirms it: for product "BCDEFG" starting at index 2 it asserts `p.first == 2 && p.second == 6` (length 6), not end index 8.
- **Fix:** Do NOT change the emitted convention (callers/tests depend on (start, length)). Instead fix the Doxygen on the header to say "pairs of (start_index, length)" so it matches behavior, and ideally rename/clarify by adding a typedef or a comment cross-referencing that this differs from ProteaseDigestion::digest. Source/ABI compatible (doc-only).
- **Verifier correction:** The problem is a misleading @param[out] doc, not a misleading symbol name. EnzymaticDigestion.h:115 documents the pair overload's output as "pairs of start and end positions", but every code path emits (start, length): EnzymaticDigestion.cpp:538, :341 (via digestAfterTokenize_), and :559. This is especially trap-prone because the derived ProteaseDigestion::digest pair overload (ProteaseDigestion.cpp:119) genuinely emits (start, past-the-end) and is documented as such (ProteaseDigestion.h:62) — so the two same-shaped APIs use opposite conventions. Severity is medium (not high): the sole existing in-repo consumer of the pair overload, FragmentIndex.cpp, already uses the (start, length) interpretation correctly, so no current data is wrong; the risk is a future caller trusting the doc. Fix is Doxygen-only (change to "(start_index, length)" and note the divergence from ProteaseDigestion::digest); source/ABI compatible.
- **Verified:** Independently verified against the source. The pair overload EnzymaticDigestion::digestUnmodified (EnzymaticDigestion.cpp:505) emits (start, length) on every path: unspecific path line 538 `output.emplace_back(i, j - i)`; normal path delegates to digestAfterTokenize_ which does line 341 `output.emplace_back(fragment_positions[i-1], l)` with l a length; semi path line 559 `output.emplace_back(p.first, p.second - p.first)` (converting semiSpecificD

### [CHEM-8] ProteaseDigestion::digest / EnzymaticDigestion::digestUnmodified — Two sibling pair-of-Size digestion APIs use opposite pair conventions (start,end) vs (start,length)
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ProteaseDigestion.h

```cpp
Size ProteaseDigestion::digest(const AASequence& protein, std::vector<std::pair<size_t,size_t>>& output, ...) const  vs  Size EnzymaticDigestion::digestUnmodified(..., std::vector<std::pair<Size, Size>>& output, ...) const
```
- **Expectation:** Two digestion functions in the same class hierarchy that both fill a std::vector<std::pair<Size,Size>> of digestion-product positions should use the same convention for the pair members.
- **Actual:** ProteaseDigestion::digest emits (start, past-the-end) pairs (`output.emplace_back(begin, pep_positions[i])`, doc: "start and past-the-end indices"), while the inherited EnzymaticDigestion::digestUnmodified pair overload emits (start, length). A developer who learned one convention and switches to the other writes a silent off-by-(end-start) bug. The internal semiSpecificDigestion_ helper is shared by both and emits (start,end), which ProteaseDigestion uses directly but digestUnmodified post-converts to (start,length) — a footgun for anyone reading the code.
- **Evidence:** ProteaseDigestion.cpp lines 119 `output.emplace_back(begin, pep_positions[i]);` and 140; ProteaseDigestion.h lines 62 doc "start and past-the-end indices of peptides". Contrast EnzymaticDigestion.cpp lines 538/341/559 (length) and EnzymaticDigestion.h line 115 doc.
- **Fix:** Align the documentation explicitly in both headers and add a one-line cross-reference warning ("NOTE: unlike EnzymaticDigestion::digestUnmodified which returns (start,length), this returns (start,end)"). If a future ABI break is acceptable, standardize on (start,length) everywhere. Doc-only fix is source/ABI safe.
- **Verifier correction:** The two sibling pair-of-Size digestion APIs genuinely use opposite conventions (ProteaseDigestion::digest = (start, past-the-end); EnzymaticDigestion::digestUnmodified = (start, length)), confirmed in source. However the claim mis-states the sibling's documentation: EnzymaticDigestion.h:120's doc (line 115) reads "Digestion products as vector of pairs of start and end positions" — it does NOT document "length", and is in fact factually wrong because the implementation emits (start, length). So the public header doc actively misleads; the (start,length) reality is only disclosed in non-public .cpp comments. This makes a doc-only fix necessary on BOTH headers (correct EnzymaticDigestion.h to say (start,length) and add the cross-reference warning to ProteaseDigestion.h). Severity is medium, not high: each function is self-consistent in isolation and misuse requires cross-applying outputs; the resulting off-by-(end-start) error is silent (wrong substring extents) but not triggered on the default single-API path. ABI impact of the recommended doc fix is none; only the optional "standardize on (start,length) everywhere" would be ABI-breaking.
- **Verified:** Verified against source. ProteaseDigestion::digest (pair overload) emits (start, past-the-end): ProteaseDigestion.cpp:119/140 `output.emplace_back(begin, pep_positions[i])`, header ProteaseDigestion.h:62 doc "start and past-the-end indices". The inherited sibling EnzymaticDigestion::digestUnmodified (pair overload) emits (start, length): EnzymaticDigestion.cpp:538 `emplace_back(i, j - i)`, :341 `emplace_back(fragment_positions[i-1], l)` with l=le

### [CHEM-13] DigestionEnzymeDB::getEnzyme — Name lookup is case-sensitive except for the canonical name; synonyms only match exact case
`medium` · `documentation-mismatch` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/DigestionEnzymeDB.h

```cpp
const DigestionEnzymeType* getEnzyme(const std::string& name) const
```
- **Expectation:** Given the @note 'enzymes are registered in regular and in toLowercase() style, if unsure use toLowercase', a caller may assume lookups are case-insensitive and pass e.g. a lowercased synonym.
- **Actual:** Only the canonical name is indexed in both original and lowercased form; synonyms are indexed only in their original case (addEnzyme_ inserts `enzyme_names_[*it] = enzyme` for synonyms with no lowercase variant). So getEnzyme(toLower(synonym)) throws ElementNotFound even though getEnzyme(synonym) succeeds. The 'use toLowercase' advice silently fails for synonyms.
- **Evidence:** DigestionEnzymeDB.h lines 181-186: canonical name added as `enzyme_names_[name]` and `enzyme_names_[StringUtils::toLower(name)]`, but the synonym loop adds only `enzyme_names_[*it] = enzyme` (no lowercased synonym). Lookup (line 75) is a plain map find with no case folding. @note at lines 71-72 advertises toLowercase usage.
- **Fix:** Either also index lowercased synonyms in addEnzyme_, or correct the @note to state that the toLowercase guarantee applies only to canonical names. Implementation fix is source/ABI safe (internal); doc fix is also safe.
- **Verifier correction:** getEnzyme's @note promises lookups work after toLowercase, but addEnzyme_ only lowercases the canonical name (lines 181-182), not synonyms (lines 183-186). Thus getEnzyme(toLower(synonym)) throws Exception::ElementNotFound for any mixed-case synonym (e.g. "Clostripain", "Glu-C"), while getEnzyme(synonym) with original case succeeds. The failure is a thrown exception (loud/recoverable), not silently wrong results — so this is a documentation-vs-implementation mismatch of medium severity, not a silent-failure of high severity. Fix: either index StringUtils::toLower(*it) for synonyms in addEnzyme_, or narrow the @note to state the toLowercase guarantee applies only to the canonical name. Both fixes are internal to a header-only template (ABI: none).
- **Verified:** Code confirms the claim exactly. In DigestionEnzymeDB.h addEnzyme_ (lines 181-186), the canonical name is indexed both verbatim (`enzyme_names_[name]`) and lowercased (`enzyme_names_[StringUtils::toLower(name)]`), but the synonym loop inserts only `enzyme_names_[*it] = enzyme` with no lowercased variant. getEnzyme (line 75) is a plain `enzyme_names_.find(name)` with zero case folding, and the @note at lines 71-72 advertises a general "if unsure u

### [CHEM-14] IsotopeDistribution::IsotopeDistribution() — Default-constructed IsotopeDistribution is not empty: it contains one peak (0,1), so size()==1
`medium` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h

```cpp
IsotopeDistribution()
```
- **Expectation:** A default-constructed container is empty; size() returns 0 and iterating yields nothing, like std::vector or after clear().
- **Actual:** The default ctor pushes a single dummy peak. In IsotopeDistribution.cpp the body is `distribution_.emplace_back(0, 1);`, so a fresh `IsotopeDistribution d;` has `d.size() == 1`, `d.begin() != d.end()`, `d.getMostAbundant()` == Peak1D(0,1), and `d.getMax()/getMin()` return 0. This differs from the state after `clear()` (which is empty).
- **Evidence:** IsotopeDistribution.cpp:33-36: `IsotopeDistribution::IsotopeDistribution() { distribution_.emplace_back(0, 1); }`. Header doc (IsotopeDistribution.h:63-65) only says "Default constructor" with no mention of the seeded peak.
- **Fix:** Document the seeded (mass=0, prob=1) peak in the constructor doxygen so callers do not treat a default-constructed instance as empty. Changing the ctor to start empty would be behavior-breaking; keep ABI and just fix the docs (or add a static factory making the intent explicit).
- **Verifier correction:** Claim is factually correct in all specifics (ctor body, size()==1, getMostAbundant()==Peak1D(0,1), getMax()/getMin()==0, undocumented, and the clear()-vs-default asymmetry). Only severity is adjusted from the implied medium-to-high down to a solid medium: the wrong value is a non-physical mass=0 peak that is somewhat self-evident on inspection and, in real OpenMS usage, the distribution is virtually always repopulated by a CoarseIsotopePatternGenerator/FineIsotopePatternGenerator before consumption, so silent-wrong-result exposure is real but limited. Recommendation stands: document the seeded (mass=0, prob=1) peak in the constructor doxygen (ABI-safe, abi_impact=none); changing the ctor to start empty would be behavior-breaking but also ABI-none.
- **Verified:** Evidence verified line-by-line. IsotopeDistribution.cpp:33-36 is exactly `IsotopeDistribution::IsotopeDistribution() { distribution_.emplace_back(0, 1); }`, so a default-constructed instance has size()==1, a non-empty iteration range, getMostAbundant()==Peak1D(0,1), and getMax()/getMin()==0 (confirmed via cpp lines 62-92). The header (IsotopeDistribution.h:63-65) documents only "Default constructor" with no mention of the seeded peak; the only "s

### [CHEM-17] IsotopeDistribution::merge — merge() silently re-sorts by mass and trims both tails (and intensities) by min_prob, beyond "merging into bins"
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h

```cpp
void merge(double resolution, double min_prob)
```
- **Expectation:** From the name and doc "creates a new container and assigns each isotope to the nearest bin" with a `resolution`, a caller expects only rebinning at the given resolution; min_prob would just be the per-bin probability floor.
- **Actual:** merge() also (1) calls sortByMass() reordering the distribution, (2) calls trimLeft(min_prob) and trimRight(min_prob) deleting low-probability isotopes from BOTH ends before binning, and (3) calls trimIntensities(min_prob) afterward. So data outside [first>=min_prob, last>=min_prob] is dropped, not merged. The mass span used for binning is taken AFTER trimming, so results depend on min_prob in a non-obvious way.
- **Evidence:** IsotopeDistribution.cpp:240-273: `sortByMass(); trimLeft(min_prob); trimRight(min_prob); ... trimIntensities(min_prob);`. Header doc (IsotopeDistribution.h:123-130) mentions only nearest-bin assignment and downsampling.
- **Fix:** Document that merge() sorts by mass and discards isotopes with probability < min_prob at both tails plus per-bin, so callers do not assume probability is conserved. Doc-only fix is ABI-safe.
- **Verifier correction:** merge() does sort by mass and trim both tails plus per-bin entries by min_prob, exactly as evidenced. Refinement to severity: the dropped data is limited to isotopes below the caller's own min_prob threshold and the operation is documented as a downsample, so this is a documentation gap / non-conservation surprise (medium) rather than high-severity uncontrolled data loss. Recommendation stands: document that merge() sorts by mass, discards isotopes with probability < min_prob from both tails BEFORE binning (which also shifts the mass range used to size the output), and applies a per-bin probability floor, so probability is not conserved.
- **Verified:** Verified against IsotopeDistribution.cpp:240-273 and the header doc at IsotopeDistribution.h:123-130. The implementation does exactly what the claim states: merge() calls sortByMass() (line 243), then trimLeft(min_prob) (244) and trimRight(min_prob) (245) which delete entries from BOTH tails where intensity < min_prob (confirmed via trimLeft 210-220 / trimRight 194-208), then computes mass_range from the ALREADY-TRIMMED container (line 248) so ou

### [CHEM-19] IsotopeDistribution::insert — insert(mass, intensity) always appends (push_back), never inserts in mass order
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h

```cpp
void insert(const Peak1D::CoordinateType& mass, const Peak1D::IntensityType& intensity)
```
- **Expectation:** For a container that is conceptually mass-ordered and offers sortByMass(), `insert(mass, intensity)` suggests placing the new isotope at the correct position (sorted insertion), like std::set::insert / std::map::insert.
- **Actual:** It unconditionally appends to the end via push_back, which can leave the distribution unsorted. Callers relying on subsequent positional methods (getMax/getMin via MZLess are fine, but trimLeft/trimRight/merge assume order) may silently get wrong results.
- **Evidence:** IsotopeDistribution.h:189-192: `inline void insert(...) { distribution_.push_back(Peak1D(mass, intensity)); }`.
- **Fix:** Rename intent in docs to clarify it appends (or provide an append() alias), and note that the caller must sortByMass() afterward if order matters. Adding a clearly-named append() is additive; renaming would break ABI, so prefer doc + alias.
- **Verifier correction:** `insert(mass, intensity)` is an undocumented append (push_back) that leaves the distribution in caller-determined order. The name misleadingly suggests ordered/sorted insertion for a mass-ordered container that offers sortByMass(). The footgun is real but narrower than claimed: merge() is safe because it sortByMass()es internally; getMax/getMin/getMostAbundant/averageMass/renormalize/trimIntensities are order-independent; only trimLeft/trimRight assume mass order. No in-tree caller is broken (all insert in ascending mass order). Severity is medium (recoverable footgun, no silent corruption of existing use), not high. Recommended remedy: add a doc comment stating it appends and that callers must sortByMass() if order matters, and optionally add an additive append() alias (source-compatible). Renaming would be ABI-breaking and is not necessary.
- **Verified:** Evidence confirmed: IsotopeDistribution.h:189-192 defines `insert(mass, intensity)` as an unconditional `distribution_.push_back(...)`, with NO doc comment. For a container that is conceptually mass-ordered, exposes sortByMass(), and has order-sensitive trim operations, the name `insert` genuinely misleads — in C++ `insert` connotes ordered/keyed insertion (std::set/std::map) or positional insertion (vector::insert(pos,val)), not a bare append. T

### [CHEM-20] CoarseIsotopePatternGenerator::approximateFromPeptideWeight — approximateFromPeptideWeight scales m/z spacing by charge but the input `mass` is taken as-is (not divided by charge)
`medium` · `unit-or-index` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h

```cpp
static IsotopeDistribution approximateFromPeptideWeight(double mass, UInt num_peaks = 20, UInt charge = 1)
```
- **Expectation:** For charge > 1, a self-consistent m/z output would divide the monoisotopic mass by charge as well as the neutron spacing; or the API would clearly document that `mass` is a charge-1 (neutral+proton) value placed at the monoisotopic position unchanged.
- **Actual:** The mono peak is placed at the raw `mass` (`result[0] = Peak1D(mass, 1.0f)`), while subsequent peaks are spaced by NEUTRON_MASS_U/charge. So for charge>1 the base m/z is not divided by charge but the isotope spacing is, producing an inconsistent m/z axis unless the caller already pre-divided `mass`. The doc says `mass` is the "m/z of monoisotopic peak (with charge = 1)", which is easy to misread when also passing charge>1.
- **Evidence:** CoarseIsotopePatternGenerator.cpp:201-219: `result[0] = Peak1D(mass, 1.0f); ... result[k] = Peak1D(mass + (k * NEUTRON_MASS_U / charge), ...)`. Header doc (CoarseIsotopePatternGenerator.h:186-188) describes `mass` as charge-1 m/z but couples it with a `charge` param.
- **Fix:** Clarify in the doxygen that `mass` is used verbatim as the monoisotopic m/z and the caller is responsible for any charge division of the base mass; only the inter-peak spacing is divided by charge. Doc-only, ABI-safe.
- **Verifier correction:** For charge>1, approximateFromPeptideWeight produces an internally inconsistent m/z axis: the monoisotopic peak is placed at the raw input mass (a charge-1 m/z per the doc) while subsequent peaks are spaced by NEUTRON_MASS_U/charge. Thus the base m/z is not charge-scaled but the isotope spacing is. The header doc describes mass as a charge-1 m/z but does not warn that, for charge>1, the base position is left unscaled — yielding wrong absolute m/z positions for multiply-charged species. Intensities are unaffected (the function's tested purpose). Recommended fix: either document that the caller must pre-divide the base mass for charge>1 (only spacing is divided), or divide the base mass appropriately for self-consistency; doc-only fix is ABI-safe.
- **Verified:** Code confirmed at CoarseIsotopePatternGenerator.cpp:210,218: result[0]=Peak1D(mass,1.0f) places the mono peak at the raw input mass, while result[k]=Peak1D(mass + k*NEUTRON_MASS_U/charge, ...) divides ONLY the inter-peak spacing by charge. So for charge>1 the base m/z is left at the (documented) charge-1 position while spacing is scaled to the charge-z value — an internally inconsistent m/z axis. The header doc (line 186) calls mass the "m/z of m

### [CHEM-22] IsoSpecWrapper::run / IsoSpec*Wrapper::run — run() invalidates the wrapper object as a side effect; a second call is undefined
`medium` · `lifetime-side-effect (one-shot consume; documented)` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsoSpecWrapper.h

```cpp
virtual IsotopeDistribution run() = 0;
```
- **Expectation:** A method named run() returning a result is typically idempotent or at least re-callable; nothing in the name suggests it consumes/destroys the object.
- **Actual:** Calling run() invalidates the object: it must not be called again, and nothing other than destruction is valid afterwards. This is a strong, easily-missed lifetime constraint behind a plain verb name.
- **Evidence:** IsoSpecWrapper.h:150 "@note Calling this method invalidates the object! In future versions this limitation might be removed." and 142-149 describing the one-shot usage pattern.
- **Fix:** The constraint is documented (good), but the name is still surprising; consider rvalue-ref-qualifying run() (`run() &&`) in a future ABI break so the compiler enforces single-use, or rename to runAndConsume(). For now keep the prominent doc note. Doc is ABI-safe; ref-qualifier is breaking.
- **Verifier correction:** run() is a one-shot consuming operation on IsoSpecWrapper subclasses, and this is prominently DOCUMENTED at the declaration (IsoSpecWrapper.h:150), so it is not a "hidden" side effect. The concrete risk is real but narrower than claimed: only IsoSpecTotalProbWrapper::run() truly consumes the object -- it advances the internal IsoLayeredGenerator without reset, so a second call silently yields a truncated/empty IsotopeDistribution (no crash/exception). IsoSpecThresholdWrapper::run() calls ITG->reset() first and is effectively re-callable, contradicting the blanket "a second call is undefined." Severity is medium, not high: a second call can produce silently-wrong results, but only under a usage pattern the prominent @note warns against and the canonical temporary-object idiom (`X(...).run()`) plus deleted copy ctors actively steer away from; the Threshold subclass even self-heals. Recommendation stands: the current doc note is ABI-safe (abi_impact none for status quo). Enforcing single-use via an rvalue-ref qualifier (`run() &&`) or renaming to runAndConsume() would be a source-breaking/ABI-breaking change and should be deferred to a planned API break; a cheaper non-breaking improvement is to make the consumed state explicit (e.g. throw on second call) and to make TotalProb's run() reset like Threshold's so the documented one-shot restriction is actually consistent.
- **Verified:** Evidence verified in the actual source. Header IsoSpecWrapper.h:153 declares `virtual IsotopeDistribution run() = 0;` and lines 141-150 contain exactly the quoted note ("Calling this method invalidates the object! ... the method should not be called again, nor anything other than destroying the object should be done with it") plus the temporary-object usage idiom. The behavior is real, not merely a warning: IsoSpecTotalProbWrapper::run() (cpp 251

### [CHEM-25] IMSIsotopeDistribution::SIZE / ABUNDANCES_SUM_ERROR — Behavior-controlling parameters are public mutable static (global) members shared by all instances
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h

```cpp
static size_type SIZE; static abundance_type ABUNDANCES_SUM_ERROR;
```
- **Expectation:** Truncation length and normalization tolerance look like per-distribution settings; a caller would not expect one instance's configuration to affect every other distribution in the process.
- **Actual:** SIZE (how many peaks size()/folding expose) and ABUNDANCES_SUM_ERROR (normalize() tolerance) are non-const public statics. Writing IMSIsotopeDistribution::SIZE = n from anywhere changes size(), operator*= folding width, and getMasses()/getAbundances() for every distribution globally, with no instance-level indication.
- **Evidence:** Header lines 122-125: `static abundance_type ABUNDANCES_SUM_ERROR;` and `static size_type SIZE;`. operator*= uses `peaks_container dest(SIZE);` and setMinimumSize_() resizes to SIZE; normalize() compares against ABUNDANCES_SUM_ERROR.
- **Fix:** ABI-safe: document loudly that these are process-global, must be initialized before first use, and are not thread-safe to mutate. Longer term (breaking): make them per-instance members with accessors. Flag abi as source-compatible for any move to instance members.
- **Verifier correction:** SIZE and ABUNDANCES_SUM_ERROR are public, non-const statics defined with NO initializer (IMSIsotopeDistribution.cpp:21-22), so they zero-initialize to 0. They are genuinely process-global behavior switches: SIZE controls size() truncation, operator*= folding width / setMinimumSize_() resize, and getMasses()/getAbundances() output length; ABUNDANCES_SUM_ERROR is normalize()'s tolerance. Mutating either from anywhere affects all instances globally and is not thread-safe. The doxygen/comment wording suggests per-distribution settings and is misleading. The class is unusable (empty output, zero-length fold buffers) until SIZE is set, but the codebase contains no mutation sites; it is an obscure IMS-internal class (only IMSElement uses it). ABI-safe recommendation: document loudly that these are process-global, must be initialized before first use, and are not thread-safe to mutate. Longer term: make them per-instance members; this is source-compatible only if read accessors replace direct IMSIsotopeDistribution::SIZE access (direct readers would break, but there are none in-tree besides the read-only pyOpenMS bindings).
- **Verified:** Evidence verified line-for-line. Header (src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h:122,125) declares `static abundance_type ABUNDANCES_SUM_ERROR;` and `static size_type SIZE;` — both public, non-const, mutable statics. They are defined with NO initializer in the .cpp (IMSIsotopeDistribution.cpp:21-22), so they are zero-initialized: SIZE=0, ABUNDANCES_SUM_ERROR=0. SIZE drives size() (h:162 std::min(peaks_.s

### [CHEM-27] MassDecompositionAlgorithm::getDecompositions — getDecompositions takes its output vector first and the input mass second, opposite to the sibling RealMassDecomposer API
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h

```cpp
void getDecompositions(std::vector<MassDecomposition> & decomps, double weight)
```
- **Expectation:** Across the same module, decomposition queries take the mass as the leading input and return the result; sibling RealMassDecomposer::getDecompositions(double mass, double error) returns the decompositions by value.
- **Actual:** This getDecompositions is void and takes the result container as the FIRST parameter and the mass as the SECOND, inverting both the return convention and the argument order of the sibling class. A caller used to RealMassDecomposer must remember decomps comes first here, and the method does not clear `decomps`, so pre-existing entries are appended to silently.
- **Evidence:** Header line 63: `void getDecompositions(std::vector<MassDecomposition> & decomps, double weight);` vs RealMassDecomposer.h line 79: `decompositions_type getDecompositions(double mass, double error);`. The pyOpenMS wrapper had to construct a fresh vector to paper over this (bind_misc.cpp:1713).
- **Fix:** ABI-safe: document the out-param-first order and whether decomps is cleared (impl should clear it). Additive: offer a return-by-value overload `std::vector<MassDecomposition> getDecompositions(double weight)` matching the sibling convention; deprecate the out-param form. abi none for the additive overload.
- **Verifier correction:** getDecompositions uses output-first / mass-second / void with an undocumented non-clearing append (confirmed by the impl push_back-only loop and by the test's manual decomps.clear() before reuse), inverting the return-by-value, mass-first convention of the delegated backend RealMassDecomposer::getDecompositions(double,double). RealMassDecomposer is vendored `ims`-namespace code (the backend this wrapper calls), not a strict peer sibling, and the practical hazard is the undocumented no-clear append rather than the argument order alone. No production C++ callers exist outside the unit test and the pyOpenMS shim. ABI-safe additive fix: document that decomps is appended-to (or have the impl clear it), and optionally add a return-by-value overload std::vector<MassDecomposition> getDecompositions(double weight) matching the backend convention. abi_impact: none.
- **Verified:** Evidence verified against the actual code. Header line 63 is exactly `void getDecompositions(std::vector<MassDecomposition> & decomps, double weight);` — output container first, mass second, void return. The impl (MassDecompositionAlgorithm.cpp:52-73) only push_backs and never clears `decomps`; the unit test itself must manually call `decomps.clear()` (MassDecompositionAlgorithm_test.cpp:56) before reusing the vector, proving the silent-append be

### [CHEM-30] TheoreticalSpectrumGenerator::generateSpectrum — generateSpectrum documents 'c- and z-ions' for ECD/ETD but emits c, z+1 and z+2 ions (no z-ions)
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h

```cpp
static MSSpectrum generateSpectrum(const Precursor::ActivationMethod& fm, const AASequence& seq, int precursor_charge)
```
- **Expectation:** Per the header doc 'Activation method ECD or ETD will generate only c- and z-ions', a caller expects classic z-ion m/z values in the result.
- **Actual:** The implementation explicitly sets add_z_ions=false and enables add_zp1_ions/add_zp2_ions, so it produces c-ions plus z+1 and z+2 radical cations (different m/z than z-ions) and NO z-ions. The z-ion peaks the doc promises never appear.
- **Evidence:** Header line 83: 'Activation method 'ECD' or 'ETD' will generate only c- and z-ions.' Impl TheoreticalSpectrumGenerator.cpp lines 281-286: theo_gen_settings.setValue("add_c_ions", "true"); setValue("add_z_ions", "false"); setValue("add_zp1_ions", "true"); setValue("add_zp2_ions", "true");
- **Fix:** Fix the doc comment to say 'c-ions and z+1/z+2 radical ions' (the chemically correct ExD products). No ABI change. If true z-ions are also intended, additionally set add_z_ions=true; that is a behavior change and should be opt-in.
- **Verifier correction:** For ECD/ETD, generateSpectrum emits c-ions plus z+1 (z., radical) and z+2 (z') ions and explicitly disables classic z-ions (add_z_ions=false). The header's "only c- and z-ions" is inaccurate within OpenMS's own ion vocabulary (where ZIon/add_z_ions is a distinct, separately-controlled ion type). Because add_metainfo defaults to false here, the shifted m/z (~+1/+2 Da relative to ZIon) appear unlabeled. Recommended fix: change the doc to "c-ions and z+1/z+2 radical ions" (chemically correct ExD products). Severity is medium, not high: the output is the chemically appropriate ExD radical series, so this is a documentation/naming mismatch that can mislead callers building expected z-ion m/z arrays, not a correctness bug producing nonsensical data. ABI: none (doc-only).
- **Verified:** Evidence verified exactly. Header line 83 says "Activation method 'ECD' or 'ETD' will generate only c- and z-ions." The implementation (TheoreticalSpectrumGenerator.cpp:279-287) for ECD/ETD explicitly sets add_c_ions=true, add_z_ions=FALSE, add_zp1_ions=true, add_zp2_ions=true, add_b_ions=false, add_y_ions=false. So classic z-ions (OpenMS ZIon, the type controlled by add_z_ions) are deliberately NOT emitted; instead c-ions plus z+1 (z., radical) 

### [CHEM-31] SpectrumAnnotator::addPeakAnnotationsToPeptideHit — 'addPeakAnnotationsToPeptideHit' overwrites existing PeakAnnotations rather than adding to them
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/SpectrumAnnotator.h

```cpp
void addPeakAnnotationsToPeptideHit(PeptideHit& ph, const PeakSpectrum& spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa, bool include_unmatched_peaks = false) const
```
- **Expectation:** A method named 'add...ToPeptideHit' whose doc says it 'directly fills the PeakAnnotations vector' suggests it appends matches to the hit's existing annotations.
- **Actual:** It builds a fresh local vector and calls ph.setPeakAnnotations(std::move(peak_annotations)), discarding any annotations the PeptideHit already had. Calling it twice (e.g. to merge two TSG configurations) silently throws away the first result.
- **Evidence:** SpectrumAnnotator.cpp line 429 'std::vector<PeptideHit::PeakAnnotation> peak_annotations;' built fresh, then line 481 'ph.setPeakAnnotations(std::move(peak_annotations));'. No read of ph.getPeakAnnotations().
- **Fix:** Either rename intent in docs to 'sets/replaces the PeakAnnotations', or make it truly additive by appending to ph.getPeakAnnotations() before setting. Safest additive fix: update the doc to say it replaces existing annotations; consider an overload/flag to append.
- **Verifier correction:** Severity is medium, not unbounded-high. The replacement is silent data loss, but the dominant usage (annotating a freshly-loaded/just-built PeptideHit, as both existing test calls and the typical workflow do) has nothing to lose, so it does not produce wrong scientific values or crashes for the common path. The risk materializes only when the hit already carries PeakAnnotations (e.g. populated by a search-engine adapter or a prior call to merge two TSG configs), where the prior set is silently discarded — recoverable by re-running, but loud nothing warns the caller. Recommended fix: rename the doc to state it sets/replaces the existing PeakAnnotations, or make it additive by appending to ph.getPeakAnnotations() (optionally via a flag), to match the additive connotation of "add...ToPeptideHit".
- **Verified:** Independently verified in src/openms/source/CHEMISTRY/SpectrumAnnotator.cpp. addPeakAnnotationsToPeptideHit builds a fresh local vector at line 429 (std::vector<PeptideHit::PeakAnnotation> peak_annotations;) and at line 481 calls ph.setPeakAnnotations(std::move(peak_annotations)). It never reads ph.getPeakAnnotations(). I confirmed PeptideHit::setPeakAnnotations (PeptideHit.cpp:378-381) is a plain replacing setter: fragment_annotations_ = std::mo

### [CHEM-33] Tagger::getTag — getTag appends to 'tags' and then sorts/dedups the caller's pre-existing entries too
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/Tagger.h

```cpp
void getTag(const std::vector<double>& mzs, std::vector<std::string>& tags) const
```
- **Expectation:** The doc '@param[out] tags the vector of tags, that is filled with this function' implies the output vector is populated with the tags for this spectrum (i.e. cleared/owned).
- **Actual:** It never clears 'tags'; it inserts new tags at tags.end() and then sort+unique the WHOLE vector. Any pre-existing content of 'tags' is silently retained, reordered, and de-duplicated against the new results. A caller reusing one vector across spectra accumulates and reorders prior tags without warning.
- **Evidence:** Tagger.cpp line 168 'tags.insert(tags.end(), tags_local.begin(), tags_local.end());' then lines 172-177 sort/unique/erase over all of 'tags'. No tags.clear() anywhere in getTag.
- **Fix:** Document the append+dedup behavior explicitly, or change the doc/contract to clear first. If preserving ABI, at minimum amend the doc to '@param[in,out] tags: new tags are appended; the combined list is sorted and de-duplicated'.
- **Verifier correction:** getTag does not clear 'tags' and appends new results before sorting/de-duplicating the combined vector (Tagger.cpp:168,172-177), contradicting the header's '@param[out]' contract (Tagger.h:54). However, no existing caller is affected: OpenNuXL and all tests/bindings pass a freshly-constructed empty container per spectrum. The real risk is latent — a future caller reusing one vector across spectra would silently accumulate and reorder prior tags. Recommend amending the doc to [in,out] (or adding tags.clear()); doc-only fix has no ABI impact. Severity downgraded high->medium since there is no current silent corruption.
- **Verified:** Code and evidence confirmed verbatim. Tagger.cpp getTag (lines 150-178) never calls tags.clear(); at line 168 it does tags.insert(tags.end(), tags_local.begin(), tags_local.end()) inside an omp critical, then at lines 172-177 sort/unique/erase over the ENTIRE tags vector. The header (Tagger.h:54) documents the parameter as '@param[out] tags the vector of tags, that is filled with this function'. The [out] annotation is meaningful here because the

### [CHEM-34] NucleicAcidSpectrumGenerator::getMultipleSpectra — getMultipleSpectra infers polarity from only the smallest charge and never validates the documented 'all positive or all negative' precondition
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h

```cpp
void getMultipleSpectra(std::map<Int, MSSpectrum>& spectra, const NASequence& oligo, const std::set<Int>& charges, Int base_charge = 1) const
```
- **Expectation:** Doc says 'All values in charges must be either positive or negative.' The sibling getSpectrum() throws IllegalArgument on mixed signs, so a caller expects the same enforcement here.
- **Actual:** Polarity is taken from a single element: 'bool negative_mode = *charges.begin() < 0;'. A mixed-sign set is not rejected; it silently processes everything in the polarity of the smallest element, producing physically meaningless spectra for the opposite-sign charges. Inconsistent with getSpectrum which throws.
- **Evidence:** NucleicAcidSpectrumGenerator.cpp line 391 'bool negative_mode = *charges.begin() < 0;' with no sign-consistency check. Contrast getSpectrum lines 349-353 which throw IllegalArgument when 'max_charge * min_charge < 0'.
- **Fix:** Add the same sign-consistency check getSpectrum uses and throw IllegalArgument on mixed signs. Body-only change; no ABI impact.
- **Verifier correction:** getMultipleSpectra infers polarity solely from the smallest set element (*charges.begin(), the most-negative due to ascending set ordering) and never validates the documented all-same-sign precondition, unlike getSpectrum which throws IllegalArgument on mixed signs. With a mixed-sign set the opposite-sign charges are silently SKIPPED by the polarity loop, yielding empty spectra (when add_metainfo_ pre-creates map entries at lines 397-404) or missing map keys — i.e. silently incomplete output, not wrong-polarity peaks. Fix: add the same sign-consistency check getSpectrum uses and throw Exception::IllegalArgument on mixed signs (body-only, no ABI impact).
- **Verified:** Evidence is verified in-tree. NucleicAcidSpectrumGenerator.cpp:391 derives polarity from a single element ('bool negative_mode = *charges.begin() < 0;') and getMultipleSpectra performs NO sign-consistency check, while the sibling getSpectrum (lines 349-353) throws Exception::IllegalArgument on mixed signs (max_charge*min_charge < 0). The header (line 64) documents the precondition "All values in charges must be either positive or negative" but do

### [CHEM-35] NucleicAcidSpectrumGenerator::getSpectrum — getSpectrum silently swaps min_charge/max_charge when |max|<|min| and appends to (does not clear) the output spectrum
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h

```cpp
void getSpectrum(MSSpectrum& spectrum, const NASequence& oligo, Int min_charge, Int max_charge) const
```
- **Expectation:** A 'getSpectrum' taking an explicit (min_charge, max_charge) pair is expected to honor the arguments as given and to fill a fresh output spectrum.
- **Actual:** If abs(max_charge) < abs(min_charge) it silently does 'swap(max_charge, min_charge)', so passing the arguments in the wrong order is masked instead of flagged. It also never clears 'spectrum'; peaks are appended to whatever the caller passed and the whole vector is then re-sorted. Additionally the loop guard 'z < oligo.size()' silently caps the max charge at the oligo length with no warning.
- **Evidence:** NucleicAcidSpectrumGenerator.cpp lines 354-357 'if (abs(max_charge) < abs(min_charge)){ swap(max_charge, min_charge); }', line 376 loop 'z <= abs(max_charge) && z < oligo.size()', no spectrum.clear() before the loop, line 382 'spectrum.sortByPosition();'.
- **Fix:** Document the append-not-clear contract (consistent with the AA TSG note) and the silent min/max swap and charge cap; or throw on inverted/out-of-range charges. Prefer a doc clarification to avoid ABI/behavior breaks.
- **Verifier correction:** getSpectrum's three surprises are real but are undocumented intentional behaviors rather than bugs: (1) it appends to (never clears) the output spectrum — a deliberate OpenMS generator convention that the sibling TheoreticalSpectrumGenerator documents but this class does not; the test itself must call spectrum.clear(true) between calls. (2) It silently swaps min/max when abs(max)<abs(min) instead of validating, despite already throwing on sign mismatch. (3) The loop guard `z < oligo.size()` silently caps the max charge at the oligo length (visible in the test's "last one is missing" assertions). Recommendation: document all three in the header (mirroring the AA TSG note) and/or throw on inverted/out-of-range charges; this is a doc/contract clarification with no ABI impact.
- **Verified:** All three quoted facts are verified verbatim in src/openms/source/CHEMISTRY/NucleicAcidSpectrumGenerator.cpp. (1) Silent swap at lines 354-357: `if (abs(max_charge) < abs(min_charge)) { swap(max_charge, min_charge); }` — confirmed. Notably the function DOES throw on sign mismatch (lines 349-352) but silently reorders inverted magnitudes, an internal inconsistency, so it is not a defensible convention. (2) No spectrum.clear(): confirmed — getSpect

### [CHEM-37] TheoreticalSpectrumGenerator::getPrefixAndSuffixIonsMZ — getPrefixAndSuffixIonsMZ appends to the output vector (no clear) despite '@param[out] spectrum'
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h

```cpp
void getPrefixAndSuffixIonsMZ(std::vector<float>& spectrum, const AASequence& peptide, int charge) const
```
- **Expectation:** An '@param[out] spectrum' that 'Generates a simple tandem MS Spectrum' implies the vector is populated with this peptide's ions (fresh).
- **Actual:** The function never clears 'spectrum'; each add* call push_backs onto it, then std::sort sorts the entire vector including any pre-existing entries. Reusing a vector across calls silently merges peptides and re-sorts prior content.
- **Evidence:** TheoreticalSpectrumGenerator.cpp lines 1172-1207: loop calls addPrefixAndSuffixIons_ (which push_back) with no spectrum.clear(); line 1203 'std::sort(spectrum.begin(), spectrum.end());'.
- **Fix:** Document the append-and-sort behavior in the header, or clear() at entry. Doc clarification is ABI-safe; clearing would be a behavior change.
- **Verified:** see above

### [CHEM-38] SpectrumAnnotator::addIonMatchStatistics — addIonMatchStatistics also mutates the input spectrum (adds DataArrays and meta values), not just the PeptideIdentification
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/SpectrumAnnotator.h

```cpp
void addIonMatchStatistics(PeptideIdentification& pi, MSSpectrum &spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const
```
- **Expectation:** Doc says it 'Adds ion match statistics to pi PeptideIdentification' and marks '@param[in] spec'. A caller expects spec to be read-only input and only pi/its hits to be modified.
- **Actual:** It calls annotateMatches(spec, ...) for each hit, which writes IonName/IonMatchError DataArrays and sets fragment_mass_tolerance meta values on 'spec'. The non-const 'MSSpectrum& spec' is mutated as an undocumented side effect.
- **Evidence:** SpectrumAnnotator.cpp line 142 'annotateMatches(spec, *ph, tg, sa);' inside addIonMatchStatistics; annotateMatches lines 129-130 'spec.setMetaValue(...)' plus DataArray writes. Header line 76 marks '@param[in] spec'.
- **Fix:** Document that 'spec' is modified (in/out) by this call, or operate on a copy if read-only input is intended. Doc/in-out annotation fix is ABI-safe.
- **Verifier correction:** addIonMatchStatistics mutates the input spectrum well beyond the claim: besides the DataArray/meta writes in annotateMatches (spec.setMetaValue for fragment_mass_tolerance and fragment_mass_tolerance_ppm at lines 129-130, and replacement of Float/String/Integer DataArrays at lines 131-133), addIonMatchStatistics itself reorders spec via spec.sortByIntensity() (line 144) and spec.sortByPosition() (line 370). Thus a caller's spectrum is left reordered with overwritten DataArrays and new meta values, despite the @param[in] doc. Fix: annotate spec as @param[in,out] (or document the mutation), or operate on a copy if read-only input is intended — a comment-only, ABI-safe change.
- **Verified:** Verified against actual code. Header line 83 declares `void addIonMatchStatistics(PeptideIdentification& pi, MSSpectrum &spec, ...)`; the Doxygen marks `@param[in] spec` (line 76) and states outputs go only to `pi`/PeptideHit (lines 73, 80). The .cpp confirms the side effect: line 142 calls `annotateMatches(spec, ...)` per hit, which at lines 129-133 calls `spec.setMetaValue("fragment_mass_tolerance"/"fragment_mass_tolerance_ppm")` and replaces t

### [CHEM-39] ProForma::Modification::resolved_mod — Public raw pointer cached in the AST is silently invalidated by JSON round-trip and by value copies that don't re-resolve
`medium` · `ownership-lifetime` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ProForma.h

```cpp
const ResidueModification* resolved_mod = nullptr;
```
- **Expectation:** A public field named resolved_mod on a value-type AST struct (Modification) suggests a stable handle to the resolved modification that travels with the data.
- **Actual:** resolved_mod is a non-owning raw pointer into ModificationsDB that is only populated by resolveModifications() and is explicitly reset to nullptr on JSON deserialization (ProFormaDataJson.h:469 `mod.resolved_mod = nullptr;  // Reset to avoid stale pointers`). toJSON() drops it entirely, so a serialize/deserialize round-trip silently loses resolution; and copying a Peptidoform copies the pointer without re-resolving. A caller reading resolved_mod after deserialization or after editing the AST can get a stale or null pointer with no signal.
- **Evidence:** ProForma.h:343 `const ResidueModification* resolved_mod = nullptr;`; ProFormaDataJson.h:469 comment 'Reset to avoid stale pointers'. The Modification to_json (ProFormaDataJson.h:450-463) never serializes resolved_mod.
- **Fix:** Document explicitly that resolved_mod is a transient, non-owning cache that is null until resolveModifications() is called and is dropped by (de)serialization, and that callers must call resolveModifications() after any deserialize/copy/edit. abi_impact none for doc; making it private with an accessor would be source-breaking.
- **Verifier correction:** resolved_mod is a transient, non-owning cache pointer into ModificationsDB (a stable process-wide singleton). It is null until resolveModifications() runs, and is silently dropped by JSON serialization: toJSON() never writes it and from_json() explicitly resets it to nullptr. Therefore a toJSON()->peptidoformFromJSON() round-trip silently loses resolution (all resolved_mod become nullptr) despite the toJSON doc claiming a complete AST. Internal mass/AASequence APIs self-heal by re-resolving on a copy, but callers reading resolved_mod directly after deserialization get a silent nullptr with no signal. The pointer is NOT dangling after a plain value-copy (the singleton keeps it valid), so the original 'stale pointer on copy/edit' framing is overstated; the real defect is silent, lossy resolution on the JSON path. Recommendation: document that resolved_mod is transient/non-owning and is nulled by (de)serialization, and that callers must call resolveModifications() after deserialize; this is a doc-only change (abi none). Making the field private with an accessor would be source-breaking.
- **Verified:** All quoted evidence is verbatim-accurate. ProForma.h:343 declares `const ResidueModification* resolved_mod = nullptr;` as a PUBLIC field on a value-type struct (Modification) with default member-wise copy. ProFormaDataJson.h:469 confirms from_json sets it to nullptr ("Reset to avoid stale pointers"), and to_json (450-463) serializes only tag/label, dropping resolution entirely. So a toJSON/peptidoformFromJSON round-trip silently loses resolution:

### [CHEM-42] MzPAF::fromPeakAnnotation — fromPeakAnnotation() returns an empty result both for 'no annotations' and for 'not valid mzPAF', conflating the two
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/MzPAF.h

```cpp
static MzPAFPeakAnnotations fromPeakAnnotation(const PeptideHit::PeakAnnotation& peak_annotation)
```
- **Expectation:** A converter from a PeakAnnotation would let the caller distinguish 'this annotation string is not mzPAF' from 'it is mzPAF but parsed to nothing'.
- **Actual:** It calls tryParseMultiple() and returns result.value_or(MzPAFPeakAnnotations{}); on any parse failure the empty default is returned. The caller cannot tell a malformed/non-mzPAF string from a legitimately empty annotation, and there is no boolean/optional signal. Documented as 'or empty if not valid mzPAF', but empty is also a normal value, so the failure is swallowed.
- **Evidence:** MzPAF.cpp:1027-1032: `auto result = tryParseMultiple(peak_annotation.annotation); return result.value_or(MzPAFPeakAnnotations{});`.
- **Fix:** Consider returning std::optional<MzPAFPeakAnnotations> (or pairing with isMzPAFFormat) so callers can detect non-mzPAF input. Additive overload keeps ABI. abi_impact source-compatible if adding an overload.
- **Verifier correction:** The claim is correct that fromPeakAnnotation() conflates 'not mzPAF' with 'empty annotation' by collapsing tryParseMultiple()'s std::optional via value_or(MzPAFPeakAnnotations{}), and the header doc 'or empty if not valid mzPAF' confirms rather than excuses it. Adjustment: severity is medium, not high — the function returns an empty (not incorrect) result with no crash or data corruption, and callers can already distinguish non-mzPAF input via the public isMzPAFFormat() or by calling tryParseMultiple() directly, so the failure is recoverable rather than silently producing wrong results. The recommended fix (add a std::optional<MzPAFPeakAnnotations>-returning overload, or pair with isMzPAFFormat) is sound and additive; abi_impact is source-compatible only if added as an overload (changing the existing return type would be breaking).
- **Verified:** Evidence is accurate. MzPAF.cpp:1030-1031 reads exactly `auto result = tryParseMultiple(peak_annotation.annotation); return result.value_or(MzPAFPeakAnnotations{});`. tryParseMultiple() (846-856) catches all exceptions and returns std::nullopt on any parse failure (malformed/non-mzPAF, empty-input throw), and fromPeakAnnotation collapses that nullopt into the same empty MzPAFPeakAnnotations{} that a legitimately-empty parse yields. Empty is genui

### [CHEM-44] IsoelectricPoint::computePI — computePI returns an in-range pH boundary (0 or 14) when no zero-crossing exists, only logging a warning
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/IsoelectricPoint.h

```cpp
static double computePI(const AASequence& seq, ProteomicsPkaScale scale = ProteomicsPkaScale::LEHNINGER, double tolerance = 1e-4)
```
- **Expectation:** A pI getter returning a double would be expected to return the actual isoelectric point, or signal failure when the net charge never crosses zero in the searched interval (e.g., a strongly basic peptide whose pI is above 14).
- **Actual:** When the net charge stays the same sign across [0,14], computePI returns whichever boundary (0.0 or 14.0) has the smaller |charge| and merely emits OPENMS_LOG_WARN. The returned value is a normal-looking pH that is not the true pI; a caller that does not scrape logs gets a silently wrong number. This is documented in the class doc but the sentinel is indistinguishable from a real result.
- **Evidence:** IsoelectricPoint.cpp:201-207: same-sign branch logs a warning and `return (std::abs(charge_low) <= std::abs(charge_high)) ? pH_low : pH_high;`.
- **Fix:** Documented behavior; consider an overload returning std::optional<double> or a status flag so non-converged cases are detectable programmatically. abi_impact source-compatible if adding an overload.
- **Verified:** Evidence confirmed verbatim. IsoelectricPoint.cpp:201-207: when computeCharge has the same sign at both pH=0 and pH=14 (true pI lies outside [0,14]), the function emits OPENMS_LOG_WARN and returns `(std::abs(charge_low) <= std::abs(charge_high)) ? pH_low : pH_high` — i.e. 0.0 or 14.0. This is a genuine silent-failure surprise: the returned double is an in-band value indistinguishable from a real pI (a peptide could legitimately have pI ≈ 14.0), w

### [CHEM-46] ResidueDB::getResidue — Two getResidue() overloads have opposite failure modes: by-name throws, by-one-letter-code silently returns nullptr
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ResidueDB.h

```cpp
const Residue* getResidue(const std::string& name) const;  /  const Residue* getResidue(const unsigned char& one_letter_code) const;
```
- **Expectation:** A caller using either overload of getResidue() to look up a residue would expect both to signal an unknown/unfound residue the same way (both throw, or both return nullptr that the caller checks).
- **Actual:** getResidue(const std::string&) throws Exception::InvalidValue when the residue is not found, while getResidue(const unsigned char&) just indexes residue_by_one_letter_code_[code] and returns nullptr for any code with no residue (e.g. lowercase or punctuation). A caller who learned 'getResidue throws' from the string overload will not null-check the char overload and dereference a nullptr.
- **Evidence:** ResidueDB.cpp:57-61 `if (r == nullptr) { throw Exception::InvalidValue(..."Residue not found: ", name); }` for the string overload, versus ResidueDB.cpp:64-69 `const Residue* ResidueDB::getResidue(const unsigned char& one_letter_code) const { //TODO why does this not throw but the std::string version does?? ... return residue_by_one_letter_code_[one_letter_code]; }` — the in-code TODO explicitly flags this inconsistency.
- **Fix:** Document the divergence prominently in the header (string throws; char returns nullptr). Ideally make the char overload's contract explicit by renaming/aliasing (e.g. add getResidueOrNull(char) and a throwing getResidue(char)). Do not silently change the throwing behavior of the char overload, as callers rely on the nullptr-return — but at minimum the header must warn that this overload returns nullptr on miss.
- **Verified:** Independently verified. ResidueDB.cpp:40-62 shows getResidue(const std::string&) throws Exception::InvalidValue when the name is not found (lines 57-60). ResidueDB.cpp:64-69 shows getResidue(const unsigned char&) just returns residue_by_one_letter_code_[one_letter_code], which is a 256-entry std::array<const Residue*,256> default-initialized to nullptr (header line 182). So any code with no registered residue (lowercase letters, punctuation, cont

### [CHEM-47] AASequence::getMonoWeight / AASequence::getFormula — getMonoWeight()/getFormula() accept any ResidueType but silently return the Internal mass/formula for several legal enum values
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/AASequence.h

```cpp
double getMonoWeight(Residue::ResidueType type = Residue::Full, Int charge = 0) const;  /  EmpiricalFormula getFormula(Residue::ResidueType type = Residue::Full, Int charge = 0) const;
```
- **Expectation:** Since the parameter type is Residue::ResidueType, a caller expects every enumerator (Precursor, BIonMinusH20, YIonMinusH20, BIonMinusNH3, YIonMinusNH3, Zp1Ion, Zp2Ion, NonIdentified, Unannotated) to either be computed correctly or to throw.
- **Actual:** The switch only handles Full/Internal/N-/C-terminal and a/b/c/x/y/z ions. All other (legal) enum values fall into `default`, which only does `OPENMS_LOG_ERROR << ... unknown ResidueType` and then returns `mono_weight` / `ef` — i.e. the Internal mass/formula with terminal additions omitted. A caller passing Residue::Precursor or Residue::Zp1Ion gets a silently wrong number, not an exception.
- **Evidence:** AASequence.cpp:586-590 `default: OPENMS_LOG_ERROR << "AASequence::getMonoWeight: unknown ResidueType\n"; } return mono_weight;` and AASequence.cpp:465-469 `default: OPENMS_LOG_ERROR << "AASequence::getFormula: unknown ResidueType\n"; } return ef;`. Residue.h:152-174 shows Precursor, BIonMinusH20, Zp1Ion etc. are valid enumerators.
- **Fix:** Throw Exception::InvalidValue (or NotImplemented) in the default branch instead of returning the Internal value, or implement the missing types. At minimum document in the header which ResidueType values are actually supported by these methods.
- **Verifier correction:** The default branch is not fully silent: it emits OPENMS_LOG_ERROR before returning the (wrong) Internal-equivalent value, so the failure is logged though the returned number is still incorrect and propagates into numeric pipelines unchecked. Severity is medium rather than high because (1) the error is logged and (2) no in-repository caller passes these unhandled ResidueType values to getMonoWeight/getFormula/getMZ — fragment generation (TheoreticalSpectrumGenerator) computes Zp1/Zp2/etc. offsets independently. The recommended fix (throw Exception::NotImplemented/InvalidValue in default, or implement Zp1Ion/Zp2Ion using the already-existing getInternalToZp1Ion()/getInternalToZp2Ion() helpers, and document the supported subset in the header) is source-compatible and ABI-neutral.
- **Verified:** Code confirmed verbatim. In AASequence.cpp, getMonoWeight (switch at lines 544-588) and getFormula (switch at lines 423-467) handle only Full/Internal/NTerminal/CTerminal and a/b/c/x/y/z ions. The 9 other valid ResidueType enumerators (Zp1Ion, Zp2Ion, Precursor, BIonMinusH20, YIonMinusH20, BIonMinusNH3, YIonMinusNH3, NonIdentified, Unannotated — Residue.h:164-172) fall into default, which returns the value computed so far. Because the terminal-mo

### [CHEM-48] Residue::getMonoWeight / Residue::getAverageWeight / Residue::getFormula — Residue weight/formula getters silently return the Full value for legal ResidueType values that hit the default branch
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/Residue.h

```cpp
double getMonoWeight(ResidueType res_type = Full) const;  /  double getAverageWeight(ResidueType res_type = Full) const;  /  EmpiricalFormula getFormula(ResidueType res_type = Full) const;
```
- **Expectation:** A caller passing any Residue::ResidueType enumerator expects the corresponding type's weight/formula, or an exception for an unsupported type.
- **Actual:** getMonoWeight/getAverageWeight/getFormula handle only a subset of the enum. Legal values such as Precursor, BIonMinusH20, YIonMinusH20, BIonMinusNH3, YIonMinusNH3, NonIdentified, Unannotated (and, for getAverageWeight, also Zp1Ion/Zp2Ion) fall through to `default`, which prints to cerr and returns the Full weight/formula. e.g. getMonoWeight(Residue::Precursor) returns the Full residue mass with only a cerr message.
- **Evidence:** Residue.cpp:372-374 `default: cerr << "Residue::getMonoWeight: unknown ResidueType" << endl; return mono_weight_;`, Residue.cpp:320-322 (getAverageWeight default returns average_weight_), Residue.cpp:273-275 (getFormula default returns formula_).
- **Fix:** Throw on unsupported ResidueType in the default branch (and remove the cerr writes), or implement the missing types. If behavior cannot change for ABI reasons, document the supported-type subset in the header.
- **Verifier correction:** The claim is factually accurate but overstates the practical danger. (1) The truly nonsensical enumerators (Precursor, BIonMinusH20/NH3, YIonMinusH20/NH3, NonIdentified, Unannotated) are CV-term/annotation markers used only in transition I/O (TraMLHandler, TransitionTSVFile) and have no meaningful single-residue weight; they are not plausibly passed to a per-Residue weight getter. The one user-facing path that takes an arbitrary type, TOPP MassCalculator, explicitly restricts fragment_type to full/internal/N-/C-terminal/a..z-ion (MassCalculator.cpp:94), so the maintainers already treat the supported subset as the contract. (2) The genuinely defensible silent-wrong-result is narrower than stated: getAverageWeight does NOT handle Zp1Ion or Zp2Ion (it stops at ZIon), so for those two legitimate per-residue ion types it silently returns the Full average weight even though getMonoWeight and getFormula compute them correctly. This getMonoWeight/getFormula-vs-getAverageWeight asymmetry is the strongest evidence of a real latent bug. Recommended fix: at minimum add the missing Zp1Ion/Zp2Ion cases to getAverageWeight; ideally replace the cerr+return-Full default branches with a thrown exception (e.g. IllegalArgument) or document the supported-type subset in the header. Throwing is source-compatible (no signature/ABI change) but is a runtime behavior change for any caller relying on the current silent fallback.
- **Verified:** Code confirmed exactly as claimed. Residue.cpp getFormula (273-275), getAverageWeight (320-322), getMonoWeight (372-374) all have a `default:` that writes "unknown ResidueType" to cerr and returns the Full member value (formula_/average_weight_/mono_weight_) instead of throwing. The ResidueType enum (Residue.h 152-173) has Full..Zp2Ion plus Precursor, BIonMinusH20/NH3, YIonMinusH20/NH3, NonIdentified, Unannotated; getMonoWeight/getFormula handle 

### [CHEM-49] ResidueDB::getResidues — getResidues() writes to std::cout when the requested residue set is empty/unknown
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ResidueDB.h

```cpp
const std::set<const Residue*> getResidues(const std::string& residue_set = "All") const;
```
- **Expectation:** A const accessor returning a set of residues should be silent; on an unknown set it should just return the documented empty set.
- **Actual:** When the result set is empty (unknown set name, or genuinely empty), the method prints 'Residue set cannot be found: ...' to std::cout. A library accessor unexpectedly emits to stdout, polluting program output, and an empty return also silently masks a typo'd set name as 'no residues'.
- **Evidence:** ResidueDB.cpp:103-106 `if (s.empty()) { cout << "Residue set cannot be found: '" + residue_set + "'" << endl; } return s;`
- **Fix:** Remove the cout (or route to OPENMS_LOG_WARN). Consider documenting that an unknown set yields an empty set so callers validate against getResidueSets(). The cout removal is behavior-only, ABI-safe.
- **Verifier correction:** getResidues() does write to std::cout ('Residue set cannot be found: ...') whenever the returned set is empty — an undocumented stdout side-effect from a const accessor, reachable via the user-configurable residue_set parameter in MassDecompositionAlgorithm. Note that the empty-set return itself IS documented in the header (line 111: 'returns an empty set if the specified residue set is not defined'); only the cout is undocumented. Fix: drop the cout or route to OPENMS_LOG_WARN. Behavior-only change, ABI-safe.
- **Verified:** Evidence confirmed verbatim. ResidueDB.cpp:103-106: `if (s.empty()) { cout << "Residue set cannot be found: '" + residue_set + "'" << endl; } return s;`. With `using namespace std;` (line 18) and `#include <iostream>` (line 16), this is genuinely std::cout/stdout, not a redefined log stream — grep shows no OPENMS_LOG/LogStream anywhere in the file. So getResidues(), a const singleton accessor (the header docstring even stresses every accessor is 

### [CHEM-52] Residue::getModificationName — getModificationName() returns the empty string for user-defined modifications even when one is set
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/Residue.h

```cpp
const std::string& getModificationName() const;
```
- **Expectation:** getModificationName() on a modified residue should return a non-empty identifier for the modification that is present (consistent with isModified() returning true).
- **Actual:** It returns modification_->getId(), and getId() is empty by definition for user-defined modifications (e.g. mass-tag mods like M[+12321] created via setModificationByDiffMonoMass). So a residue with isModified()==true can return an empty getModificationName(), and callers using emptiness as 'unmodified' get a wrong answer. ResidueModification's own toString()/getFullId() carry the user-defined name, making getId() the inconsistent choice here.
- **Evidence:** Residue.cpp:467-472 `const std::string& Residue::getModificationName() const { ... if (!isModified()) return EMPTY; return modification_->getId(); }`. ResidueModification.h:314-315 documents isUserDefined() as 'user-defined modification (empty id)'.
- **Fix:** Document that this returns the short Id (empty for user-defined mods) and that callers should use getModification()->getFullId()/toString() to get a name for user-defined mods; or add a getModificationFullId()/getModificationToString() helper. Doc/additive fix, ABI-safe.
- **Verified:** Independently verified in source. Residue::getModificationName() (src/openms/source/CHEMISTRY/Residue.cpp:467-472) returns EMPTY when !isModified(), otherwise modification_->getId(). The key premise is provably true: ResidueModification::isUserDefined() is literally `id_.empty() && !full_id_.empty()` (ResidueModification.cpp:544-547), and createUnknownFromMassString sets FullId/FullName but deliberately NOT Id (Residue.cpp:458 path; ResidueModifi

### [CHEM-53] AASequence::getPrefix / AASequence::getSuffix — The 'index' parameter of getPrefix/getSuffix is actually a length/count, not a position
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/AASequence.h

```cpp
AASequence getPrefix(Size index) const;  /  AASequence getSuffix(Size index) const;
```
- **Expectation:** A parameter named 'index' in a sequence class (where getResidue(index)/operator[](index) are 0-based positions) strongly implies a position, so getSuffix(index) might be read as 'the suffix starting at position index'.
- **Actual:** getPrefix(index) returns the FIRST index residues and getSuffix(index) returns the LAST index residues — index is a count. getSuffix(2) of 'PEPTIDE' is 'DE', not 'PTIDE'. Using the same name 'index' as the positional accessors invites an off-by-everything bug, especially for getSuffix where position vs count give different substrings.
- **Evidence:** AASequence.h:507-511 docs 'returns a peptide sequence of the first index residues' / 'of the last index residues'. AASequence.cpp:701-717 getPrefix inserts [begin, begin+index); AASequence.cpp:719-736 getSuffix inserts [begin+(size()-index), end).
- **Fix:** Rename the parameter to 'length'/'num' in the header signature and docs (doc/param-name change, ABI-safe). Optionally add clearer aliases. Do not change semantics.
- **Verifier correction:** The parameter name 'index' in getPrefix/getSuffix denotes a residue COUNT/length, contradicting the same name's positional meaning in operator[] and getSubsequence within the same class. Recommendation (rename param to 'length'/'num' in header + docs) is correct and ABI-safe: C++ parameter names are not part of the mangled symbol, and C++ has no named-argument calls, so renaming is purely a documentation/source-readability change — abi_impact none, fully source-compatible.
- **Verified:** Independently verified against the actual code. Header AASequence.h:507-511 declares getPrefix(Size index)/getSuffix(Size index) with docs "a peptide sequence of the first index residues" / "of the last index residues". AASequence.cpp:701-717 implements getPrefix as insert[begin, begin+index) (first index residues) and 719-736 getSuffix as insert[begin+(size()-index), end) (last index residues) — so 'index' is unambiguously a COUNT/length, not a 

### [CHEM-57] NASequence::getPrefix — getPrefix(length)/getSuffix(length) throw IndexOverflow when length == size(), so you cannot request the whole sequence
`medium` · `inconsistent-edge-case` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/NASequence.h

```cpp
NASequence getPrefix(Size length) const
```
- **Expectation:** For a length-based prefix/suffix API (the docs explicitly say 'given length (not end index!)'), requesting length == size() should return the entire sequence. A prefix of N elements from an N-element sequence is well-defined and is the common edge case.
- **Actual:** getPrefix throws when `length >= seq_.size()` (NASequence.cpp:86-89), and getSuffix likewise (NASequence.cpp:95-97). So `seq.getPrefix(seq.size())` throws IndexOverflow instead of returning a copy of the whole sequence. The test even encodes this: `TEST_EXCEPTION(Exception::IndexOverflow, seq.getPrefix(10))` for a length-6 sequence, but length==size() is also rejected by the same `>=`.
- **Evidence:** if (length >= seq_.size())\n{\n  throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, seq_.size() - 1);\n}
- **Fix:** Change the guard to `length > seq_.size()` so a length equal to the sequence size returns the full sequence (matching length semantics and AASequence::getPrefix conventions). Document the off-by-one behavior if kept. This is a behavior change, so guard with care/tests.
- **Verifier correction:** getPrefix/getSuffix throw IndexOverflow (a LOUD exception, not a silent failure) when length == size(), so requesting the whole sequence by length fails — contradicting the documented length semantics ("given length, not end index"). The directly analogous AASequence::getPrefix/getSuffix in the same module use `index > size()` and return `*this` when index==size(), establishing that length==size() should yield the full sequence; NASequence diverges. Note the cited test sequence is 7 residues, not 6. Recommended fix: change guard to `length > seq_.size()` (and special-case length==size to return the full sequence), matching AASequence. ABI: none (definition-only change).
- **Verified:** Evidence is confirmed verbatim. NASequence.cpp:86 `if (length >= seq_.size()) throw IndexOverflow(...)` and the identical guard at line 95 for getSuffix. The header (NASequence.h:449,458) documents pure length semantics ("given length (not end index!)") and does NOT document the length==size() rejection. So `seq.getPrefix(seq.size())` throws instead of returning the whole sequence. This is NOT a domain convention: the sibling/analog class AASeque

### [CHEM-61] ModificationsDB::getModification — Doc says "the first one is returned" on multiple matches, but the LAST match is actually returned
`medium` · `misleading-documentation` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ModificationsDB.h

```cpp
const ResidueModification* getModification(const std::string& mod_name, const std::string& residue = "", ResidueModification::TermSpecificity term_spec = NUMBER_OF_TERM_SPECIFICITY) const
```
- **Expectation:** When several modifications match, the documented behavior (header line 142: "If more than one matching modification is found, the first one is returned with a warning") leads a caller to expect the first matching DB entry.
- **Actual:** getModification() delegates to searchModificationsFast(), whose own doc (header lines 119-122) admits it returns "the one occurrence ... (the last occurrence)". The impl loop assigns `mod = it;` for every match without breaking, so the value returned is the LAST match in iteration order, not the first. The two adjacent public methods therefore document opposite tie-break behavior.
- **Evidence:** searchModificationsFast loop: `for (const auto& it : modifications->second) { if (...) { mod = it; nr_mods++; } }` (ModificationsDB.cpp:154-163) keeps the last; getModification doc: "the first one is returned with a warning" (ModificationsDB.h:142).
- **Fix:** Fix the doc to say "an arbitrary (currently the last) match is returned", or make the tie-break deterministic and consistent across getModification()/searchModificationsFast(). Pure doc/comment fix is ABI-safe; making it deterministically "first" is source-compatible.
- **Verifier correction:** The doc at ModificationsDB.h:142 and the warning at ModificationsDB.cpp:276 both claim the "first" match is returned, but getModification() actually returns the LAST element in std::set iteration order (sorted by pointer address) because searchModificationsFast()'s loop (cpp:154-163) assigns `mod = it` for every match without breaking. The two public methods (getModification "first" vs searchModificationsFast "last occurrence") document contradictory tie-break behavior. Note the order is pointer-address-based, so neither "first" nor "last" corresponds to any meaningful DB/insertion order — the pick is effectively arbitrary. A warning is always emitted on multiple matches, so the result is loud, not silent. Fix: correct the doc/warning to "an arbitrary (currently the last in pointer-sorted order) match is returned", or make the tie-break deterministic and consistent across both methods. Pure doc/comment fix is ABI-safe.
- **Verified:** Independently verified in the source. ModificationsDB.h:142 documents "If more than one matching modification is found, the first one is returned with a warning", and the runtime warning at ModificationsDB.cpp:276 also says "picking the first one only." But getModification() (cpp:255-284) delegates to searchModificationsFast(), whose loop (cpp:154-163) does `mod = it; nr_mods++;` for every match with NO break, so it keeps the LAST match in iterat

### [CHEM-62] ModificationsDB::findModificationIndex — Self-contradictory contract: documents returning numeric_limits<Size>::max() AND throwing on the same condition; actually always throws
`medium` · `documentation-contract-contradiction` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ModificationsDB.h

```cpp
Size findModificationIndex(const std::string& mod_name) const
```
- **Expectation:** The header gives two mutually exclusive promises: "return numeric_limits<Size>::max() if not exactly one matching modification was found or no matching residue or modification were found" and "@throw Exception::ElementNotFound if not exactly one matching modification was found." A caller cannot tell whether to check for a sentinel or catch an exception.
- **Actual:** The implementation throws Exception::ElementNotFound for not-found and for >1 match; it only ever returns numeric_limits<Size>::max() in an unreachable internal-consistency branch (which also throws right after). The documented sentinel-return path is effectively dead, so a caller who follows the doc and checks `== max()` will never see it and will instead get an uncaught exception.
- **Evidence:** ModificationsDB.cpp:299-335: throws if `!has(mod_name)`, throws if `size() > 1`, and the `index == max()` case at line 331 throws again. Header lines 175-178 promise both a max() return and a throw.
- **Fix:** Remove the contradictory "return numeric_limits<Size>::max()" sentence from the doc (the method throws). ABI-safe doc fix.
- **Verifier correction:** The contract is self-contradictory but the runtime behavior is not a silent failure: ModificationsDB::findModificationIndex ALWAYS throws Exception::ElementNotFound on failure (not-found, >1 match, or the unreachable internal-inconsistency branch) and never returns numeric_limits<Size>::max(). The documented sentinel-return path is dead. The correct category is a documentation/contract contradiction (loud exception, not silent), severity medium. Fix: remove the "return numeric_limits<Size>::max() ..." sentence from ModificationsDB.h so the doc states only the throwing behavior; ABI-safe doc-only change.
- **Verified:** Independently confirmed against the source. Header ModificationsDB.h:172-180 documents two mutually exclusive outcomes for the SAME condition: "return numeric_limits<Size>::max() if not exactly one matching modification was found or no matching residue or modification were found" AND "@throw Exception::ElementNotFound if not exactly one matching modification was found." The implementation (ModificationsDB.cpp:297-336) throws on !has() (line 301),

### [CHEM-64] ModificationDefinitionsSet::isCompatible — isCompatible() ignores the set's own max_mods_per_peptide limit
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/ModificationDefinitionsSet.h

```cpp
bool isCompatible(const AASequence& peptide) const
```
- **Expectation:** A ModificationDefinitionsSet stores both the allowed modifications and a max-modifications-per-peptide cap (setMaxModifications). A method named isCompatible(peptide) is expected to test the peptide against the full constraint set, including that the number of variable mods does not exceed the configured maximum.
- **Actual:** isCompatible only checks (a) that required fixed mods are present and (b) that no mod outside the variable/fixed name sets is present. It never consults max_mods_per_peptide_, so a peptide carrying more variable modifications than the set permits is reported compatible.
- **Evidence:** ModificationDefinitionsSet.cpp:195-264 uses only getVariableModificationNames()/getFixedModificationNames(); max_mods_per_peptide_ is never referenced in the function. (The related inferFromPeptides even carries a TODO at cpp:336 noting the max is unhandled.)
- **Fix:** Document that isCompatible ignores the max-per-peptide cap, or add the count check against max_mods_per_peptide_ (behavior change; gate behind a new overload/flag to stay source-compatible).
- **Verifier correction:** isCompatible() does ignore the set's max_mods_per_peptide_ cap (confirmed at cpp:195-264; never referenced). However, the cap defaults to 0 ("unlimited"/unset) and is only meaningful after an explicit setMaxModifications() call, so the silent wrong-true result requires that explicit opt-in and is recoverable — severity is medium, not high. The recommended fix (new overload/flag, or documenting the omission) is source-compatible; an in-place behavior change would be ABI-none but a behavior change.
- **Verified:** Verified against source. ModificationDefinitionsSet::isCompatible (src/openms/source/CHEMISTRY/ModificationDefinitionsSet.cpp:195-264) only (a) checks required fixed mods are present (lines 204-227) and (b) rejects mods whose getFullId() is outside the variable/fixed name sets (lines 229-261, plus N/C-terminal checks). max_mods_per_peptide_ is genuinely never referenced in the function, despite the class storing it via setMaxModifications/getMaxM

### [CHEM-66] ModificationDefinition::operator== / operator< / std::hash — operator== compares by mod pointer identity while operator< (and the set) order by modification name — equality and ordering disagree
`medium` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/ModificationDefinition.h

```cpp
bool operator==(const ModificationDefinition& rhs) const; bool operator<(const ModificationDefinition&) const
```
- **Expectation:** For a value type stored in std::set<ModificationDefinition>, callers expect == and < to be consistent: !(a<b) && !(b<a) should mean a==b. Two definitions naming the same modification should compare equal.
- **Actual:** operator== compares the raw ResidueModification* (`mod_ == rhs.mod_`) plus the flags, while operator< compares only getModificationName() (the FullId string). Two ModificationDefinitions naming the same modification but holding different ResidueModification* (e.g. one resolved via the DB, one constructed from a stack copy or a re-interned mod) are order-equivalent (neither <) yet operator== returns false. The std::hash specialization likewise hashes pointer identity, so equal-by-name definitions hash differently. This silently breaks std::set semantics (set dedups by name via <) versus explicit == / unordered containers (which use pointer identity).
- **Evidence:** ModificationDefinition.cpp:53-58 (`mod_ == rhs.mod_ && ...`) vs cpp:67-70 (`getModificationName() < rhs.getModificationName()`); std::hash uses `reinterpret_cast<uintptr_t>(&md.getModification())` (header lines 133-141, with note "pointer identity, not content").
- **Fix:** Make operator== consistent with operator< (compare by FullId, or compare the modification by value), and update std::hash to hash the FullId rather than the pointer. This is a behavioral change; prefer doing it deliberately with tests, since current callers may rely on pointer identity. Source-compatible signature-wise.
- **Verifier correction:** operator== / std::hash and operator< are inconsistent, but the decisive divergence is that operator< compares only getModificationName() (FullId) and ignores fixed_modification_ and max_occurrences_, whereas operator== and std::hash include all three fields. Thus two ModificationDefinitions with the same modification name but different fixed/max flags are order-equivalent (neither < the other) yet compare unequal and hash differently. In std::set<ModificationDefinition> (used by ModificationDefinitionsSet for fixed_mods_ and variable_mods_) such definitions silently dedup to one, contradicting ==/unordered_set semantics. The pointer-vs-name aspect the claim leads with is largely masked by the DB interning identical-named mods to the same ResidueModification* (so == and < usually agree on the name); the robust bug is the missing flag fields in operator<. Fix: make operator< compare the same tuple (name/FullId, fixed, max_occurrences) as operator==, or make == compare by FullId+flags by value, and keep std::hash consistent with the chosen equality. Behavioral change; add tests since some callers may rely on current pointer/name behavior.
- **Verified:** Independently confirmed in ModificationDefinition.cpp and .h. operator== (cpp:53-58) compares mod_ pointer AND fixed_modification_ AND max_occurrences_; operator< (cpp:67-70) compares ONLY getModificationName() (the FullId). std::hash (header:133-141) hashes pointer identity + the two flags. So == and < genuinely disagree: the equivalence !(a<b)&&!(b<a) is strictly coarser than ==. The strongest, fully provable break is NOT the pointer-vs-name an

### [CHEM-67] ModificationDefinition::ModificationDefinition(const ResidueModification&, bool, UInt) — Constructor stores the address of the passed-by-reference ResidueModification (borrowed pointer), creating a dangling risk
`medium` · `ownership-lifetime` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ModificationDefinition.h

```cpp
explicit ModificationDefinition(const ResidueModification& mod, bool fixed = true, UInt max_occur = 0)
```
- **Expectation:** A constructor taking `const ResidueModification&` by reference reads like it copies/registers the modification by value; a caller may pass a temporary or a local and expect the ModificationDefinition to remain valid afterward.
- **Actual:** The constructor stores `mod_ = &mod;` — it retains a non-owning pointer to the caller's object. If the referenced ResidueModification is a temporary or goes out of scope (and is not a stable ModificationsDB-interned entry), the ModificationDefinition dangles. getModification()/getModificationName() then read freed memory. The string-name constructor, by contrast, resolves to a stable DB pointer.
- **Evidence:** ModificationDefinition.cpp:35-40: `mod_(&mod)`. Member is `const ResidueModification* mod_;` (header line 105). No copy of the referenced object is made.
- **Fix:** Document loudly that the referenced ResidueModification must outlive the definition (ideally must be a ModificationsDB-owned pointer), or change to resolve/intern via ModificationsDB::searchModification() to a stable pointer. Doc fix is ABI-safe; interning is a behavior change.
- **Verifier correction:** The constructor stores a borrowed, non-owning pointer (mod_ = &mod) to the caller's ResidueModification rather than copying or interning it; the referent must outlive the ModificationDefinition. This is a real but currently LATENT hazard: it is undocumented and modeled by the unit test (which passes a same-scope stack local), but has no in-tree callers that pass a true temporary, so it is an API lifetime trap rather than an active crash. Recommended fix: document the lifetime requirement at the declaration, or intern via ModificationsDB::searchModification() to a stable DB pointer (the latter is a behavior change but ABI-safe; the doc-only fix has no ABI impact).
- **Verified:** Independently verified. Header line 105 declares `const ResidueModification* mod_;` (non-owning raw pointer). ModificationDefinition.cpp:35-40 confirms the constructor body `mod_(&mod)` stores the ADDRESS of the by-reference argument with no copy. getModification() (lines 87-95) and getModificationName() (lines 97-104) dereference mod_, so a destroyed referent yields use-after-free. The contrast in the claim is also correct: the string ctor (27-3

### [CHEM-68] ModificationDefinitionsSet::findMatches — Three adjacent same-typed bool flags whose doc order (consider_variable before consider_fixed) is the reverse of the signature order
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ModificationDefinitionsSet.h

```cpp
void findMatches(std::multimap<double, ModificationDefinition>& matches, double mass, const std::string& residue = "", ResidueModification::TermSpecificity term_spec = NUMBER_OF_TERM_SPECIFICITY, bool consider_fixed = true, bool consider_variable = true, bool is_delta = true, double tolerance = 0.01) const
```
- **Expectation:** The Doxygen block lists @param consider_variable before @param consider_fixed (header lines 137-138); a reader cross-referencing the docs to the call site naturally maps the first positional bool to consider_variable.
- **Actual:** The actual signature order is `consider_fixed, consider_variable, is_delta` (header line 144). A caller who follows the doc ordering will silently swap consider_fixed and consider_variable. With three bare adjacent bools (consider_fixed, consider_variable, is_delta), an accidental swap compiles cleanly and changes which mod set / mass mode is used with no diagnostic.
- **Evidence:** Doc order: ModificationDefinitionsSet.h:137 (consider_variable) then :138 (consider_fixed); signature order: :144 (consider_fixed first, then consider_variable, then is_delta).
- **Fix:** Reorder the @param lines to match the signature, and consider an enum/struct-of-options overload to make call sites self-documenting. Doc reorder is ABI-safe; an options overload is source-compatible/additive.
- **Verifier correction:** The doc-order vs signature-order mismatch is real (lines 137-138 list consider_variable before consider_fixed; line 144 declares consider_fixed first). The "silently swap" hazard is genuine but currently latent: every call site uses the defaults (consider_fixed=true, consider_variable=true), so no live bug exists today. The most damaging misuse requires a caller to explicitly pass false for exactly one flag following the (wrong) doc order. Fix: reorder the @param lines to match the signature (ABI-safe). An enum/struct-of-options overload would be source-compatible/additive but is not required to resolve the documentation defect.
- **Verified:** Verified against the actual code. Header lines 137-138 document @param consider_variable BEFORE @param consider_fixed, while the signature on line 144 declares them in the reverse order (consider_fixed, consider_variable, is_delta). The .cpp (lines 318-334) confirms the two bools select genuinely different data (fixed_mods_ vs variable_mods_), so swapping them silently changes which modification set is matched — wrong-but-plausible results, no di

### [CHEM-2] EmpiricalFormula::estimateFromWeightAndComp / estimateFromMonoWeightAndComp / estimateFromWeightAndCompAndS — 'estimate' return-flag semantics are documented inconsistently (true vs '1' vs false) and on failure the object is left in a partially-built state
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/EmpiricalFormula.h

```cpp
bool estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P);
```
- **Expectation:** A bool returning method has an unambiguous meaning, and on a documented failure the caller can rely on the object being unchanged or fully valid.
- **Actual:** The header doc for estimateFromWeightAndComp/estimateFromMonoWeightAndComp says: "@return bool flag ... true = no problems, 1 = negative hydrogens requested" (EmpiricalFormula.h:144,159) — i.e. it claims the FAILURE case returns '1' (which is true), contradicting itself, while estimateFromWeightAndCompAndS correctly documents "false = negative hydrogens requested" (EmpiricalFormula.h:175). In code, failure returns false and EXITS EARLY after clearing and inserting C/N/O/S/P but BEFORE inserting H (EmpiricalFormula.cpp:146-150), so on 'false' the formula is left without hydrogens — a partially-populated, misleading composition rather than untouched or complete.
- **Evidence:** EmpiricalFormula.h:144 `@return ... true = no problems, 1 = negative hydrogens requested.`; EmpiricalFormula.cpp:146-150 `if (adjusted_H < 0) { ... return false; }` (returns before inserting H, after formula_.clear()+C/N/O/S/P insert at 134-140).
- **Fix:** Fix the doc typo so all three methods read 'true = ok, false = negative hydrogens requested', and document that on false the formula_ is left without H (or restore prior state). abi_impact: doc-only.
- **Verified:** Both evidence claims verified against the actual code. EmpiricalFormula.h:144 and :159 document @return as "true = no problems, 1 = negative hydrogens requested." Because 1 == true in C++, this is self-contradictory: it maps BOTH success and failure to true, while the implementation actually returns false on failure (EmpiricalFormula.cpp:149 for average, :186 for mono). The third overload estimateFromWeightAndCompAndS correctly documents "false =

### [CHEM-3] EmpiricalFormula::contains — contains() doc says elements of ef must be 'LESS abundant' than this, but the predicate is actually '>=' (less-or-EQUAL) and the direction wording is reversed-sounding
`low` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/EmpiricalFormula.h

```cpp
bool contains(const EmpiricalFormula& ef) const;
```
- **Expectation:** this->contains(ef) reads as 'this is a superset of ef' and should be true when this has at least as many of every element as ef (i.e. counts are >= ef's, equality allowed).
- **Actual:** The code returns false only when this->getNumberOf(e) < ef count, i.e. returns true when this's counts are >= ef's (equality allowed) — correct behavior, but the header doc states "returns true if all elements from @p ef are LESS abundant (negative allowed) than the corresponding elements of this". 'LESS abundant ... than' wrongly excludes equality and phrases the comparison so it reads as the opposite subset relation, inviting an off-by-one/equality bug.
- **Evidence:** EmpiricalFormula.h:261 `/// returns true if all elements from @p ef are LESS abundant (negative allowed) than the corresponding elements of this EmpiricalFormula`; EmpiricalFormula.cpp:379 `if (this->getNumberOf(it.first) < it.second) return false;` (so equality returns true).
- **Fix:** Reword doc to: 'returns true if, for every element in ef, this has at least as many atoms (>=, equality allowed); i.e. this contains ef as a sub-formula'. abi_impact: doc-only.
- **Verifier correction:** The code is correct: contains(ef) returns true iff, for every element in ef, this has at least as many atoms (this-count >= ef-count, equality allowed) — i.e. this contains ef as a sub-formula. The header doc is wrong on the equality boundary: "LESS abundant ... than" implies strict '<' and excludes equality, but the implementation allows equality ('<='). The directional wording ("ef ... LESS abundant than this") is actually correct, not reversed. Recommended doc: "returns true if, for every element in @p ef, this EmpiricalFormula has at least as many atoms (>=, equality allowed); i.e. this contains @p ef as a sub-formula (element counts may be negative)." Doc-only fix; no ABI/behavioral change.
- **Verified:** Evidence confirmed verbatim. Header EmpiricalFormula.h:261 documents contains() as "returns true if all elements from @p ef are LESS abundant (negative allowed) than the corresponding elements of this", while the implementation (EmpiricalFormula.cpp:375-385) iterates ef's elements and returns false only when this->getNumberOf(e) < ef-count, i.e. returns true iff every element of ef satisfies ef-count <= this-count (equality allowed; getNumberOf r

### [CHEM-6] Element setters (setMonoWeight/setAverageWeight/setSymbol/setName/setAtomicNumber/setIsotopeDistribution) — Element exposes full mutators, but every Element a caller can obtain is owned by the singleton ElementDB and handed out as const Element* (immutable by design)
`low` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/Element.h

```cpp
void setMonoWeight(double weight); void setSymbol(const std::string&); ...
```
- **Expectation:** If a class has public setters, a caller expects to be able to mutate instances they hold and have that be meaningful.
- **Actual:** ElementDB owns all Elements via std::unique_ptr<const Element> and returns const Element* from getElement/getNames/getSymbols (ElementDB.h:53-71, members at 115-119); the DB is documented 'Singleton that stores elements' and 'immutable'. EmpiricalFormula keys its map on const Element* into these shared singleton objects. The public non-const setters on Element therefore cannot be applied to any DB element and, if they could, would corrupt the shared, pointer-identity-based lookup for the whole process. The mutable API surface is misleading about how Elements are actually used.
- **Evidence:** ElementDB.h:51-53 `returns a pointer to the (immutable) singleton instance`; ElementDB.h:106 `std::unique_ptr<const Element> e`; Element.h:69,75,93 public setters; EmpiricalFormula.h:66 `typedef std::map<const Element*, SignedSize> MapType_;`.
- **Fix:** Document on Element that instances vended by ElementDB are immutable/shared and must not be mutated; consider clarifying that setters are only for constructing custom Elements not registered in the DB. abi_impact: doc-only.
- **Verifier correction:** Element exposes full public mutators, yet every Element the library vends is const Element* owned by the ElementDB singleton and keyed by pointer identity in EmpiricalFormula, so the setters are never applicable to a library-obtained Element. The asymmetry is real and the header (Element.h) gives no immutability hint — the warning lives only in ElementDB.h. But Element is a genuine copyable value type, so setters DO work on caller-constructed local instances, and the const Element* vending path means the compiler already prevents mutating shared DB state without an explicit const_cast. Risk is therefore a misleading/unused mutable API surface (low), not silent corruption under reasonable use; const_cast-then-mutate would corrupt process-global pointer-identity lookups but is a loud, deliberate misuse. The pyOpenMS bindings already treat Element as read-only. Recommendation (doc-only on Element.h) has no ABI impact.
- **Verified:** Evidence verified against the actual code. Element.h:63-93 declares public non-const setters (setAtomicNumber/setAverageWeight/setMonoWeight/setIsotopeDistribution/setName/setSymbol) with no caveat. ElementDB owns every Element via make_unique<const Element>(...) (ElementDB.cpp:586,651, buildElement_) and vends ONLY const Element* from getElement/getNames/getSymbols/getAtomicNumbers (ElementDB.h:53-71, members 115-119, doc'd "immutable singleton 

### [CHEM-9] DigestionEnzyme::operator< / operator== — operator< orders by name only, but operator== compares name+synonyms+regex+description
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/DigestionEnzyme.h

```cpp
bool operator<(const DigestionEnzyme& enzyme) const  vs  bool operator==(const DigestionEnzyme& enzyme) const
```
- **Expectation:** For a type with both operator< and operator==, callers (and std::set/std::sort/std::map) expect consistency: !(a<b) && !(b<a) should imply a==b.
- **Actual:** operator< compares only getName(), while operator== additionally compares synonyms_, cleavage_regex_, and regex_description_. Two enzymes with the same name but different regexes are 'equivalent' under operator< (neither orders before the other) yet are operator!= unequal. Code that relies on the ordering to deduplicate or to treat ordering-equivalence as equality (e.g. std::set<DigestionEnzyme> semantics vs std::find using ==) silently disagrees.
- **Evidence:** DigestionEnzyme.cpp lines 161-164 `return this->getName() < enzyme.getName();` vs lines 138-144 comparing name_, synonyms_, cleavage_regex_, regex_description_.
- **Fix:** Document that operator< is a name-only ordering and is not consistent with operator== (or make operator< tie-break on the same fields == uses). Behavior change to operator< could affect sorted containers, so prefer documenting the contract; field-aligning the comparison is source-compatible but a subtle behavior change.
- **Verifier correction:** The inconsistency is real and correctly described, but no current OpenMS code relies on the ordering for dedup/equality: DigestionEnzymeDB uses std::set<const DigestionEnzymeType*> (pointer ordering) and name-keyed std::map with enforced unique names, and operator< has no value-typed sorted-container callers. So the "silently disagree" failure mode is latent/hypothetical, not active. Severity is low (mild surprise / invites future misuse), not a silently-wrong-results bug. Recommendation to document operator< as a name-only ordering (or tie-break on the same fields == uses) is appropriate; documenting is ABI-none and field-aligning operator< is source-compatible (subtle behavior change only).
- **Verified:** Evidence verified verbatim. DigestionEnzyme.cpp:161-164 `operator<` returns `this->getName() < enzyme.getName()` (name-only), while `operator==` at :138-144 compares name_, synonyms_, cleavage_regex_, and regex_description_. So ordering-equivalence (!(a<b)&&!(b<a)) does NOT imply operator== equality — a genuine, real C++ inconsistency between a type's strict-weak-ordering equivalence and its equality, and the header comment only says "order opera

### [CHEM-10] ProteaseDigestion::peptideCount — peptideCount is non-const although it only reads state, while the heavier digest() is const
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/CHEMISTRY/ProteaseDigestion.h

```cpp
Size peptideCount(const AASequence& protein)
```
- **Expectation:** A pure counting query named peptideCount should be const (and is expected to be, since the much heavier digest() overloads on the same class are const).
- **Actual:** peptideCount is declared non-const, so it cannot be called on a const ProteaseDigestion& even though it performs no mutation (it only calls the const tokenize_ and reads enzyme_/missed_cleavages_). This surprises callers holding a const reference and is inconsistent with digest() being const.
- **Evidence:** ProteaseDigestion.h line 72 `Size peptideCount(const AASequence& protein);` (no const) vs lines 50/69 `... digest(...) const`. ProteaseDigestion.cpp lines 48-69: body only reads members and calls const helpers.
- **Fix:** Add `const` to peptideCount. This is source-compatible for callers but changes the mangled name (ABI). Prefer adding a const overload or marking const in the next ABI-breaking window; documenting is the conservative interim fix.
- **Verifier correction:** peptideCount() is indeed non-const while it only reads state and the heavier digest() overloads are const — a real const-correctness asymmetry. However the practical impact is low, not medium: misuse is a compile-time error (loud/recoverable, never silent wrong data or crash), and the sole in-tree caller already holds a non-const reference so nothing is currently blocked. Marking it const is an ABI break (the mangled name gains the const-this qualifier), though source-compatible for callers; a const overload would avoid the ABI break.
- **Verified:** Code verified independently. ProteaseDigestion.h line 72 declares `Size peptideCount(const AASequence& protein);` with no trailing const, while both digest() overloads at lines 50 and 69 are const. The .cpp body (lines 48-69) only reads enzyme_ and missed_cleavages_, and calls the const helper tokenize_ (declared const at EnzymaticDigestion.h:180) and protein.toUnmodifiedString() on the parameter, not on *this — it performs no mutation of the obj

### [CHEM-11] DigestionEnzymeRNA::getCutsAfterRegEx / getCutsBeforeRegEx / getThreePrimeGain / getFivePrimeGain — RNA-enzyme string getters return by value while base-class string getters return const ref
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/DigestionEnzymeRNA.h

```cpp
std::string getCutsAfterRegEx() const
```
- **Expectation:** Sibling string accessors in the same DigestionEnzyme hierarchy should have a consistent return type; the base class returns const std::string& for getName/getRegEx/getRegExDescription, so callers may expect the same for the RNA-specific getters and bind a const ref without copying.
- **Actual:** DigestionEnzymeRNA's getCutsAfterRegEx/getCutsBeforeRegEx/getThreePrimeGain/getFivePrimeGain all return std::string by value, silently copying a member that exists as a plain std::string. Same inconsistency in DigestionEnzymeProtein (getPSIID/getXTandemID return by value). Not a correctness bug, but a surprising per-call allocation and an inconsistency with the base.
- **Evidence:** DigestionEnzymeRNA.h lines 38/44/50/56 return `std::string` by value (members are std::string, DigestionEnzymeRNA.h lines 66-69). Base getters return `const std::string&` (DigestionEnzyme.h lines 76, 91, 97).
- **Fix:** For consistency and to avoid copies, change to `const std::string&` returns (source-compatible for most callers, but changes ABI/mangling). At minimum document the by-value choice. Conservative: leave as-is and note the inconsistency.
- **Verifier correction:** The code facts are accurate, but the framing overstates impact. The current by-value getters are not a bug and carry no ABI/correctness problem (abi_impact of the existing code = none; only the *proposed fix* to const-ref would be ABI-breaking). By-value string returns are a standard, safe idiom; binding a const ref to them is harmless (lifetime extension), unlike the reverse. The genuine observation is a low-severity stylistic inconsistency in return-type category between base (const std::string&) and the RNA/Protein-specific string getters (std::string by value). Recommend documenting or leaving as-is rather than an ABI-breaking unification.
- **Verified:** Evidence is exactly correct. DigestionEnzymeRNA.h:38,44,50,56 declare getCutsAfterRegEx/getCutsBeforeRegEx/getThreePrimeGain/getFivePrimeGain returning std::string by value; the backing members (lines 66-69) are plain std::string, and the .cpp (DigestionEnzymeRNA.cpp:51-84) confirms each getter copies the member. The base DigestionEnzyme.h returns const std::string& for getName/getRegEx/getRegExDescription (lines 76,91,97). DigestionEnzymeProtein

### [CHEM-12] DigestionEnzymeDB::getEnzymeByRegEx vs getEnzyme — Sibling lookup-by-key methods throw different exception types for the not-found case
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/DigestionEnzymeDB.h

```cpp
const DigestionEnzymeType* getEnzymeByRegEx(const std::string& cleavage_regex) const
```
- **Expectation:** getEnzyme(name) and getEnzymeByRegEx(regex) are parallel lookup methods; a caller writing a single catch handler for 'enzyme not found' expects them to signal the same way.
- **Actual:** getEnzyme throws Exception::ElementNotFound when the name is unknown, while getEnzymeByRegEx throws Exception::IllegalArgument when the regex is unknown. A caller catching ElementNotFound around both will miss the regex case. The source even carries a TODO acknowledging this.
- **Evidence:** DigestionEnzymeDB.h line 78 `throw Exception::ElementNotFound(...)` (getEnzyme) vs lines 89-91 `// @TODO: why does this use a different exception than "getEnzyme"?` then `throw Exception::IllegalArgument(...)` (getEnzymeByRegEx).
- **Fix:** Standardize both on Exception::ElementNotFound (the regex case is genuinely 'element not found'). This is source-compatible at the call site if callers catch the common base; document if a change is deferred.
- **Verifier correction:** The inconsistency is real and confirmed (including the in-source TODO), but its severity is low rather than what the title implies: the divergence cannot silently produce wrong results — an unmatched regex case throws IllegalArgument that propagates loudly to a top-level handler. Recommendation stands: standardize getEnzymeByRegEx on Exception::ElementNotFound. ABI note: as an inline template header there is no binary-ABI break, but it is only fully source-compatible for callers catching the common BaseException; callers catching the specific IllegalArgument type would observe a behavior change. No current in-tree caller does so.
- **Verified:** Evidence verified verbatim in DigestionEnzymeDB.h: getEnzyme (line 78) throws Exception::ElementNotFound on unknown name; getEnzymeByRegEx (lines 89-91) carries the literal "@TODO: why does this use a different exception than getEnzyme?" comment and throws Exception::IllegalArgument on unknown regex. The two are documented parallel lookup methods (with matching hasEnzyme/hasRegEx predicates), and in Exception.h both ElementNotFound and IllegalArg

### [CHEM-15] IsotopeDistribution::getMax / IsotopeDistribution::getMin — getMax()/getMin() document "returns the isotope" but return only a coordinate (m/z double)
`low` · `documentation` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h

```cpp
Peak1D::CoordinateType getMax() const; Peak1D::CoordinateType getMin() const;
```
- **Expectation:** By the doc text "returns the isotope with the largest/smallest m/z" and by analogy with sibling getMostAbundant() (which returns a Peak1D), a caller expects the isotope peak (mass+probability).
- **Actual:** getMax()/getMin() return `Peak1D::CoordinateType` (a bare double m/z), not the Peak1D. The probability/intensity of that isotope is discarded. Within the same class getMostAbundant() returns a full Peak1D, making the trio inconsistent.
- **Evidence:** IsotopeDistribution.h:89-96 — getMax/getMin doc'd "returns the isotope with the largest/smallest m/z" but typed `Peak1D::CoordinateType`; getMostAbundant returns `Peak1D`. Impl IsotopeDistribution.cpp:62-78 returns `...->getMZ()`.
- **Fix:** Reword the doxygen to "returns the m/z of the isotope with the largest/smallest m/z" (the value, not the isotope). If a peak is genuinely wanted, add getMaxPeak()/getMinPeak() returning Peak1D for symmetry with getMostAbundant(). Doc fix is ABI-safe; adding overloads is additive.
- **Verifier correction:** The code behaves exactly as its return type declares: getMax()/getMin() return a `double` m/z value (Peak1D::CoordinateType), not a Peak1D. The genuine, real surprise is purely in the doxygen wording ("returns the isotope ...") plus the asymmetry with getMostAbundant() (which returns a Peak1D) — and the same loose wording leaks into the pyOpenMS docstrings. This is a low-severity documentation/API-consistency issue, not a return-value/data-loss bug: the type system makes any misread loud and immediate in C++. Recommended fix is the doc reword ("returns the m/z of the isotope with the largest/smallest m/z"), which is ABI-safe; optionally add additive getMaxPeak()/getMinPeak() returning Peak1D for symmetry.
- **Verified:** Evidence verified independently. Header (IsotopeDistribution.h:89-93) doc-comments getMax()/getMin() as "returns the isotope with the largest/smallest m/z" yet they are typed Peak1D::CoordinateType, which Peak1D.h:42 defines as a bare `double`. The .cpp (IsotopeDistribution.cpp:62-78) returns `std::max/min_element(..., MZLess())->getMZ()` — only the m/z, discarding probability. The sibling getMostAbundant() (line 80-87) returns a full Peak1D, so 

### [CHEM-16] CoarseIsotopePatternGenerator::getRoundMasses — getRoundMasses() doc contradicts setRoundMasses()/class semantics about what true vs false means
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h

```cpp
bool getRoundMasses() const;
```
- **Expectation:** Getter and setter of the same flag describe the same meaning: true => masses rounded to integers, false => accurate masses (consistent with the class doc and the round_masses ctor default of false meaning accurate masses).
- **Actual:** The getter's doc inverts and renames the meaning: "returns the current value of the flag to return expected masses (true) or atomic numbers (false)". The setter says "round masses to integer values (true) or return accurate masses (false)". So per the getter, true==expected(accurate) masses, but per the setter true==rounded integers — directly contradictory. The class header (lines 44-45) confirms false==accurate masses, so the getter doc is wrong.
- **Evidence:** CoarseIsotopePatternGenerator.h:98-99 (setter) vs 104-105 (getter): "...to return expected masses (true) or atomic numbers (false)." Impl just returns round_masses_ (CoarseIsotopePatternGenerator.cpp:54-57).
- **Fix:** Fix the getRoundMasses() doxygen to match the setter: true => integer-rounded masses, false => accurate masses. Pure doc fix, no ABI impact.
- **Verifier correction:** getRoundMasses() doc at CoarseIsotopePatternGenerator.h:104 is incorrect: it should read like the setter, i.e. "returns whether masses are rounded to integer values (true) or returned as accurate masses (false)". The current wording ("expected masses (true) or atomic numbers (false)") both inverts and renames the flag's meaning relative to the setter (line 98), the class doc (lines 42-43, accurate-mass default), and the implementation (cpp:539-541, true => round(mass)). The claim's cited class-header lines 44-45 are slightly off (correct lines are 42-43). Pure doxygen fix, no ABI/source impact.
- **Verified:** Verified directly. Getter doc (CoarseIsotopePatternGenerator.h:104): "returns the current value of the flag to return expected masses (true) or atomic numbers (false)." Setter doc (line 98): "round masses to integer values (true) or return accurate masses (false)." Class doc (lines 42-43): "use setRoundMasses accordingly. The default is to return accurate masses" with ctor default round_masses=false (line 84). The implementation (cpp:539-541) con

### [CHEM-18] IsotopeDistribution::trimRight / IsotopeDistribution::trimLeft — trimLeft/trimRight trim by container position (not by mass), and silently assume the distribution is mass-sorted
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h

```cpp
void trimRight(double cutoff); void trimLeft(double cutoff);
```
- **Expectation:** Names "Left"/"Right" plus docs ("trims the right/left side of the isotope distribution") suggest trimming by mass position (low-mass / high-mass side).
- **Actual:** Both operate purely on the current container order: trimRight resizes off the trailing elements whose intensity < cutoff; trimLeft erases leading elements with intensity < cutoff. Neither sorts. If the distribution is intensity-sorted (e.g. after sortByIntensity()) or otherwise unsorted, "Left/Right" no longer correspond to mass sides and the result is not what the name implies. The mass-sorted precondition is undocumented.
- **Evidence:** IsotopeDistribution.cpp:194-220 — trimRight loops over `distribution_.rbegin()..rend()` then `resize`; trimLeft iterates from `begin()` and `erase(begin(), iter)`; no sort call. Header doc (IsotopeDistribution.h:132-150) speaks of "right/left side".
- **Fix:** Document that these operate on the current (assumed mass-sorted) ordering and call sortByMass() first if unsure; or have them sort internally. Doc fix is ABI-safe.
- **Verifier correction:** trimRight/trimLeft do operate on current container order with no internal sort, and the mass-sorted precondition is undocumented (accurate). But "Left/Right" are NOT misleading for the canonical mass-sorted distribution that every generator produces — there they correctly mean low-mass/high-mass tails, matching the name and domain convention. The only failure mode is calling sortByIntensity() before trimming, which no in-tree caller does. This is a low-severity documentation omission (add a @note that the distribution is assumed mass-sorted, or call sortByMass() internally), not a high-impact name/behavior mismatch.
- **Verified:** Evidence is accurate. IsotopeDistribution.cpp:194-208 (trimRight) scans distribution_ from rbegin() to the first element with intensity >= cutoff and resize()s off the trailing low-intensity entries; trimLeft (210-220) scans from begin() and erase()s the leading low-intensity entries. Neither sorts. Both are purely tail/head container operations. The header docs (IsotopeDistribution.h:132-150) say "Trims the right/left side" and do NOT document a

### [CHEM-21] CoarseIsotopePatternGenerator::estimateFromPeptideWeight (and sibling estimate* non-static methods) — estimateFrom*/estimateForFragment* are non-const despite not mutating the generator, breaking const-correctness symmetry with run()
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h

```cpp
IsotopeDistribution estimateFromPeptideWeight(double average_weight)
```
- **Expectation:** Compute-only methods that derive an IsotopeDistribution from weights (and only read max_isotope_/round_masses_/overrides) should be const, like run() which is const. A caller holding a `const CoarseIsotopePatternGenerator&` expects to be able to call them.
- **Actual:** Several estimate* methods are non-const (e.g. estimateFromPeptideWeight, estimateFromWeightAndComp, estimateFromRNAWeight, estimateForFragmentFromPeptideWeight) even though their bodies only construct an EmpiricalFormula and call ef.getIsotopeDistribution(*this) / delegate, with no member mutation. Meanwhile near-identical siblings (estimateForFragmentFromPeptideWeightAndS, estimateForFragmentFromWeightAndComp) ARE const, so the const-ness is inconsistent and surprising.
- **Evidence:** CoarseIsotopePatternGenerator.h:147,150,156,168,210,218,226,233,240,257,273,310,326 declared non-const; lines 294,349 declared const for analogous logic. Impl CoarseIsotopePatternGenerator.cpp:114-175 shows no member mutation.
- **Fix:** Mark the read-only estimate* methods const to match run() and the const fragment overloads. Adding const can be source-compatible but is technically an ABI/overload change for a non-virtual member; prefer doing it in a coordinated minor-version bump, or at minimum document the inconsistency.
- **Verifier correction:** The asymmetry is real and the methods can be made const (EmpiricalFormula::getIsotopeDistribution takes the generator by const ref; the const siblings prove it compiles). But the practical consequence of the current non-const-ness is a compile-time error for a `const&` caller, not silently wrong data — hence low, not high/medium. ABI: these are non-virtual members with no const/non-const overload twins, so adding const is source-compatible for essentially all callers (non-const objects bind to const methods) but changes the mangled symbol name → it is technically an ABI break for anyone linking the precompiled library, and should be batched into a minor-version bump as the recommendation states.
- **Verified:** Independently verified. The header (lines 147,156,168,210,218,226,233,240,257,273,310,326) declares estimateFrom*/estimateForFragmentFrom* non-const, while two near-identical siblings estimateForFragmentFromPeptideWeightAndS (line 294) and estimateForFragmentFromWeightAndComp (line 349) ARE const. The .cpp bodies confirm zero member mutation: each method only constructs a local EmpiricalFormula (and optionally a local `solver` generator), reads i

### [CHEM-23] MassDecomposition::compatible — const predicate compatible() writes diagnostic text to std::cerr
`low` · `hidden-side-effect` · ABI: `none` · src/openms/source/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.cpp

```cpp
bool MassDecomposition::compatible(const MassDecomposition& deco) const
```
- **Expectation:** A const boolean predicate named compatible() should only inspect the two objects and return true/false with no observable side effect.
- **Actual:** On every negative result the method prints the offending amino acid and its count to std::cerr (`cerr << it->first << " " << it->second << endl;`) before returning false. Any caller that probes many candidate decompositions for compatibility (e.g. de-novo tagging loops) silently spams stderr.
- **Evidence:** Lines 163-175:   bool MassDecomposition::compatible(const MassDecomposition& deco) const   {     for (...) {       if (it2 == decomp_.end() || decomp_.find(it->first)->second < it->second)       {         cerr << it->first << " " << it->second << endl;         return false;       }     }     return true;   }
- **Fix:** Remove the leftover debug `cerr` statement (it looks like forgotten debugging output). This is a pure source/ABI-compatible change of an inline-less .cpp body; no signature change. If logging is intended, route it through OPENMS_LOG_DEBUG instead of unconditional std::cerr.
- **Verifier correction:** The cerr statement (line 170) genuinely exists and is an undocumented side effect of a const predicate, but there are currently NO production callers of compatible() in the OpenMS source tree (only the pyOpenMS binding, the unit test, and the header). The described "de-novo tagging loops silently spam stderr" is therefore hypothetical, not an active defect. Severity is low (mild, loud, recoverable code-smell / leftover debug output), not high. The fix recommendation stands: remove the stray cerr (matching the behavior of the analogous containsTag()), or route through OPENMS_LOG_DEBUG if logging is ever intended. This is a source- and ABI-compatible .cpp body change.
- **Verified:** Evidence confirmed verbatim: MassDecomposition.cpp lines 163-175 show the const predicate compatible() executing `cerr << it->first << " " << it->second << endl;` on every false return (line 170). This is a genuine, undocumented hidden side effect in a const boolean predicate. The header doc (MassDecomposition.h line 81) mentions no output, and the sibling method containsTag() (lines 128-161) performs the identical member-count check WITHOUT any 

### [CHEM-26] IMSIsotopeDistribution::operator*= — operator*= mutates its const-reference argument via const_cast
`low` · `const-correctness` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h

```cpp
IMSIsotopeDistribution & operator*=(const IMSIsotopeDistribution & distribution)
```
- **Expectation:** A binary operator taking `const IMSIsotopeDistribution&` does not modify the right-hand operand.
- **Actual:** The implementation const_casts the parameter and calls setMinimumSize_() on it, which can resize the argument's internal peaks_ container (padding to the static SIZE). The const-qualified parameter is observably mutated, so a caller's distribution passed as the RHS can grow in peak count as a side effect of folding.
- **Evidence:** IMSIsotopeDistribution.cpp:75-77: `IMSIsotopeDistribution & non_const_distribution = const_cast<IMSIsotopeDistribution &>(distribution); non_const_distribution.setMinimumSize_();` where setMinimumSize_() does `peaks_.resize(SIZE)`.
- **Fix:** ABI-safe: operate on a local copy of the RHS instead of const_casting the caller's object, so the const contract is honoured. No signature change.
- **Verifier correction:** operator*=(const IMSIsotopeDistribution&) does const_cast its RHS and call setMinimumSize_(), which resizes (zero-pads) the argument's peaks_ to the static SIZE — confirmed verbatim. The mutation is only observable when the public static SIZE is set >0 by a client AND the RHS had fewer peaks than SIZE; it then zero-pads the RHS, changing its size()/getMasses()/getAbundances()/operator== but leaving real-peak abundances and getAverageMass() numerically intact. With the default SIZE==0 it is a no-op. Severity is therefore low (mild, recoverable surprise; no wrong results for real peaks, no crash, no data loss), not high. Fix is .cpp-local (fold over a copy of the RHS); no ABI/signature impact.
- **Verified:** Evidence is factually exact. IMSIsotopeDistribution.cpp:75-77 const_casts the `const IMSIsotopeDistribution&` RHS and calls setMinimumSize_() on it; setMinimumSize_() (lines 225-231) does `if (peaks_.size() < SIZE) peaks_.resize(SIZE)`, so the const argument's internal peaks_ container is mutated (padded). This is a genuine const-correctness violation, not a standard idiom or domain convention: a binary fold operator silently grows its const RHS.

### [CHEM-28] MassDecomposition::operator== — operator== compares against a string and parses it, but there is no operator== against another MassDecomposition
`low` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h

```cpp
bool operator==(const std::string& deco) const
```
- **Expectation:** A class with operator< and a string equality would also offer == between two instances of the same type; equality against a string would be a convenience addition, not the only equality.
- **Actual:** The only equality operator takes a std::string and internally constructs a temporary MassDecomposition(deco) to compare. So `md == "C2 M4"` works (and silently re-parses the string each call) but `md1 == md2` does not compile. operator< exists for the type but operator== for the type does not, an asymmetric set of relational ops.
- **Evidence:** Header lines 73-76: `bool operator<(const MassDecomposition& rhs) const;` and `bool operator==(const std::string& deco) const;` (no `operator==(const MassDecomposition&)`). Impl (MassDecomposition.cpp:101-106) builds `MassDecomposition md(deco)` to compare.
- **Fix:** ABI-safe additive: add `bool operator==(const MassDecomposition&) const` (and !=) so type-to-type equality is available and matches operator<. Keep the string overload but note it parses on each call.
- **Verifier correction:** Claim is accurate but over-rated. `operator<(const MassDecomposition&)` exists while the only equality is `operator==(const std::string&)` (which constructs a temporary and re-parses each call); there is no homogeneous `operator==`/`operator!=`, and because the string ctor is `explicit`, `md1 == md2` does not compile. The project's own test must round-trip via `.toString()` to compare two instances, confirming the gap. Correct fix is the ABI-safe additive one: add `bool operator==(const MassDecomposition&) const` and `operator!=`, keeping the string overload. Severity is low (compile-time/loud failure plus a minor re-parse inefficiency), not a silent-wrong-result or crash issue. Adding the operators is purely additive — new symbols only, no signature changes, no ambiguity (string ctor is explicit) — so abi_impact is none (source-compatible).
- **Verified:** Independently verified. Header (MassDecomposition.h:73,76) declares `bool operator<(const MassDecomposition&) const` and, as the ONLY equality, `bool operator==(const std::string&) const`. No `operator==(const MassDecomposition&)` and no `operator!=` exist anywhere (grep confirms). Impl (MassDecomposition.cpp:101-106) builds a temporary `MassDecomposition md(deco)` from the string and compares internals, i.e. it re-parses on every call. The strin

### [CHEM-29] IMSAlphabet::getElement / getMass / getName (index overloads) — Index-based accessors do no bounds checking and read out of range on a bad index (UB), while the name-based siblings throw a clean exception
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h

```cpp
const element_type & getElement(size_type index) const; mass_type getMass(size_type index) const; const name_type & getName(size_type index) const
```
- **Expectation:** Within one class the two lookup styles should fail consistently; the name-based getElement(name)/getMass(name) are documented to throw Exception::InvalidValue on a miss, so the index-based ones would be expected to validate too (or this asymmetry should be obvious).
- **Actual:** getElement(index) returns `elements_[index]` with operator[] (no bounds check); getMass(index)/getName(index) forward to it. An out-of-range index is undefined behavior (silent bad read / crash), whereas the same conceptual error via a name throws a catchable exception. The header gives no @throws on the index overloads.
- **Evidence:** Header line 109-112: `const element_type & getElement(size_type index) const { return elements_[index]; }`. Contrast getElement(name) impl (IMSAlphabet.cpp:43-53) which throws Exception::InvalidValue when not found.
- **Fix:** ABI-safe: document that the index overloads require 0 <= index < size() and are unchecked (UB otherwise), making the asymmetry with the throwing name overloads explicit. Optionally add a checked accessor. No signature change needed.
- **Verifier correction:** The asymmetry is genuine and as quoted, but it is a low-severity documentation gap, not a high-impact trap. The index overloads (getElement/getMass/getName(size_type)) use the standard, well-understood unchecked std::vector::operator[] (UB only on a caller programming error), whereas the name overloads throw a documented Exception::InvalidValue. Recommended fix: document that the index overloads require 0 <= index < size() and are unchecked, making the contrast with the throwing name overloads explicit; optionally add a checked accessor. Doc-only, ABI-safe.
- **Verified:** Evidence is accurate. Header IMSAlphabet.h:109-112 defines `getElement(size_type index) const { return elements_[index]; }` using std::vector::operator[] with no bounds check; getName(size_type) (IMSAlphabet.cpp:21-24) and getMass(size_type) (cpp:26-29) both forward to it, inheriting the unchecked access. The name-based siblings genuinely differ: getElement(const name_type&) (cpp:43-53) throws Exception::InvalidValue on a miss and is documented t

### [CHEM-32] SpectrumAnnotator::annotateMatches — annotateMatches Doxygen in/out param directions are swapped (the const input is marked [out], the mutated output marked [in])
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/SpectrumAnnotator.h

```cpp
void annotateMatches(PeakSpectrum &spec, const PeptideHit& ph, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const
```
- **Expectation:** @param[out] should mark the parameter that is written; @param[in] should mark read-only inputs. A caller relying on the docs would think 'ph' receives the annotation result and 'spec' is read-only.
- **Actual:** The doc marks '@param[out] ph' but ph is 'const PeptideHit&' and cannot be written. It marks '@param[in] spec' but spec is the non-const PeakSpectrum that actually receives the IonName/IonMatchError DataArrays and meta values. The directions are reversed.
- **Evidence:** Header line 62 '@param[in] spec ...', line 63 '@param[out] ph A spectrum identifications ...' while signature is 'PeakSpectrum &spec, const PeptideHit& ph'. Impl SpectrumAnnotator.cpp line 129-130 writes spec.setMetaValue(...); ph is consumed read-only.
- **Fix:** Swap the Doxygen directions: '@param[in,out] spec' (mutated) and '@param[in] ph' (read-only). Pure doc fix, no ABI impact.
- **Verifier correction:** The Doxygen in/out directions for annotateMatches are swapped: header line 62 '@param[in] spec' should be '@param[in,out] spec' (it is the non-const PeakSpectrum that is mutated — sorted, given setMetaValue and Float/String/Integer DataArrays in SpectrumAnnotator.cpp:105,129-133), and line 63 '@param[out] ph' should be '@param[in] ph' (it is 'const PeptideHit&', read-only, used only for getSequence()/getCharge()). Genuine but low severity: pure documentation defect, compiler-enforced constness prevents any silent data error, no ABI impact.
- **Verified:** Evidence confirmed verbatim. Header SpectrumAnnotator.h line 62 marks '@param[in] spec' and line 63 marks '@param[out] ph', while the signature (line 70) is 'void annotateMatches(PeakSpectrum &spec, const PeptideHit& ph, ...)'. The directions are objectively swapped: 'ph' is 'const PeptideHit&' and physically cannot be an output, yet is marked [out]; 'spec' is the non-const PeakSpectrum that IS mutated, yet is marked [in]. The .cpp (SpectrumAnnot

### [CHEM-36] NucleicAcidSpectrumGenerator::getMultipleSpectra — Doxygen marks the read-only input oligo as [out]
`low` · `other` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h

```cpp
void getMultipleSpectra(std::map<Int, MSSpectrum>& spectra, const NASequence& oligo, const std::set<Int>& charges, Int base_charge = 1) const
```
- **Expectation:** @param[out] denotes a parameter the function writes to; 'oligo' is the const input sequence to fragment.
- **Actual:** The header documents '@param[out] oligo Target oligonucleotide sequence' although the signature is 'const NASequence& oligo' (read-only input). Misleading direction annotation.
- **Evidence:** Header lines 57-58: '@param[out] spectra Output spectra' then '@param[out] oligo Target oligonucleotide sequence'; signature line 68 has 'const NASequence& oligo'.
- **Fix:** Change '@param[out] oligo' to '@param[in] oligo'. Pure doc fix, no ABI impact.
- **Verifier correction:** Doxygen on line 58 annotates the read-only input 'oligo' (const NASequence&) as @param[out]; it should be @param[in]. The 'spectra' parameter's @param[out] is correct. Pure documentation fix, no ABI/source impact. Low severity because the const in the signature authoritatively and visibly contradicts the annotation, preventing actual misuse.
- **Verified:** Verified against the actual code. Header line 58 reads '@param[out] oligo Target oligonucleotide sequence' while the signature on line 68 is 'void getMultipleSpectra(std::map<Int, MSSpectrum>& spectra, const NASequence& oligo, ...)' — oligo is a const reference, i.e. read-only input. The implementation (NucleicAcidSpectrumGenerator.cpp:386) only reads oligo via getUnchargedSpectrum_(oligo) and cannot mutate it (const); meanwhile 'spectra' is genu

### [CHEM-40] ProForma::toAASequence — Doc promises a throw on a 'STRICT' policy that does not exist in the ConversionPolicy enum
`low` · `documentation-inaccuracy` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ProForma.h

```cpp
static AASequence toAASequence(const Peptidoform& pf, ConversionPolicy policy = ConversionPolicy::FAIL_ON_LOSS)
```
- **Expectation:** The @throws contract should reference a policy value the caller can actually pass. The enum values are FAIL_ON_LOSS, DROP_UNLOCALISED, BEST_EFFORT.
- **Actual:** The doxygen says '@throws Exception::ConversionError if STRICT policy and conversion not possible', but there is no STRICT enumerator anywhere in the header or implementation. A caller searching for STRICT to opt into throwing behavior will not find it and may not realize FAIL_ON_LOSS is the throwing mode.
- **Evidence:** ProForma.h:820 `@throws Exception::ConversionError if STRICT policy and conversion not possible`; enum ConversionPolicy (ProForma.h:63-68) has only FAIL_ON_LOSS/DROP_UNLOCALISED/BEST_EFFORT; grep for STRICT finds only this comment.
- **Fix:** Fix the doc to reference FAIL_ON_LOSS (the throwing policy) and the actual exception type thrown. abi_impact none.
- **Verifier correction:** The @throws documentation at ProForma.h:820 references a nonexistent STRICT policy. The actual throwing policy is ConversionPolicy::FAIL_ON_LOSS (which is also the default arg). The exception type (Exception::ConversionError) is correctly documented. This is a documentation-accuracy bug, not a silent failure: toAASequence's runtime behavior is correct and throws loudly. Fix: change "if STRICT policy" to "if policy == FAIL_ON_LOSS". Severity low (mildly misleading doc, recoverable, default already throws); abi_impact none.
- **Verified:** Independently confirmed in the actual source. Header ProForma.h:820 says "@throws Exception::ConversionError if STRICT policy and conversion not possible", but the ConversionPolicy enum (ProForma.h:63-68) has only FAIL_ON_LOSS, DROP_UNLOCALISED, BEST_EFFORT. A repo-wide grep for STRICT under src/openms/source/CHEMISTRY and in the header finds it ONLY in that one comment. The implementation (ProForma.cpp:2004-2042) confirms the throwing policy is 

### [CHEM-41] ProForma::getMonoWeight(const Peptidoform&) — getMonoWeight(Peptidoform) silently ignores the peptidoform's own per-chain charge field
`low` · `return-value` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/ProForma.h

```cpp
static double getMonoWeight(const Peptidoform& pf)
```
- **Expectation:** Given that Peptidoform carries an optional charge (ProForma.h:524) and there is a getMZ(Peptidoform, charge) overload, a caller might expect getMonoWeight/getMZ to honor pf.charge if set.
- **Actual:** getMonoWeight(Peptidoform) computes only the neutral mass and never inspects pf.charge. There is no getMZ(const Peptidoform&) single-arg overload (only getMZ(pf, charge)), so the per-chain charge stored on a Peptidoform is never automatically used for m/z — unlike getMZ(const PeptidoformIon&) which does read pfi.charge. This asymmetry can lead a caller to assume pf.charge participates in m/z when it silently does not.
- **Evidence:** ProForma.cpp:2232-2243 getMonoWeight(Peptidoform) never references pf.charge; getMZ(const Peptidoform&, int charge) requires an explicit charge (ProForma.cpp:2290-2298); getMZ(const PeptidoformIon&) reads pfi.charge (ProForma.cpp:2269-2288). Peptidoform::charge exists (ProForma.h:524).
- **Fix:** Either add a getMZ(const Peptidoform&) overload that uses pf.charge (symmetry with PeptidoformIon), or document that pf.charge is ignored by the mass/m-z helpers and must be passed explicitly. abi_impact source-compatible if adding an overload.
- **Verifier correction:** The surprise is not in getMonoWeight (correctly documented as neutral mass; charge-independence is expected). The genuine, low-severity surprise is the API asymmetry: ProForma provides getMZ(const PeptidoformIon&) which automatically reads pfi.charge, but provides no getMZ(const Peptidoform&) overload that reads pf.charge — even though Peptidoform::charge (ProForma.h:524) is populated by the parser exactly in chimeric contexts (where whole-ion getMonoWeight/getMZ throw, making per-chain m/z the intended use). Callers must manually unpack pf.charge (int vs vector<AdductIon>) and call getMZ(pf, charge). No silently-wrong results occur (the existing getMZ forces an explicit charge; the convenient overload is simply absent), so this is a mild discoverability/symmetry gap, not a correctness bug. Fix: either add getMZ(const Peptidoform&) mirroring the PeptidoformIon overload (additive, source-compatible), or document on getMZ and on the Peptidoform::charge field that pf.charge is not consumed by the mass/m-z helpers.
- **Verified:** Code evidence verified accurate. getMonoWeight(const Peptidoform&) (ProForma.cpp:2232-2243) never references pf.charge; getMZ(const Peptidoform&, int) (2290-2298) requires an explicit charge and ignores pf.charge; getMZ(const PeptidoformIon&) (2269-2288) reads pfi.charge; there is NO single-arg getMZ(const Peptidoform&) overload. Peptidoform::charge (ProForma.h:524) is real and IS populated by the parser — but only in chimeric contexts (parser li

### [CHEM-43] HydrophobicityProfile::computeHydrophobicMoment — computeHydrophobicMoment hardcodes the Eisenberg scale while sibling profile methods take a configurable scale parameter
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/HydrophobicityProfile.h

```cpp
std::vector<double> computeHydrophobicMoment(const AASequence& seq, const Size window_size = 11, const double angle = 100.0) const
```
- **Expectation:** In a class whose other public methods (computeProfile, computeWindowedProfile) accept a HydrophobicityScaleMethod, a caller would expect computeHydrophobicMoment to either accept the same scale parameter or to clearly state which scale it uses.
- **Actual:** computeHydrophobicMoment has no scale parameter and silently uses HydrophobicityScaleMethod::EISENBERG internally; the header documents the angle and window but not that the hydrophobicity values are fixed to Eisenberg. A caller who set KYTE_DOOLITTLE elsewhere may assume consistency.
- **Evidence:** HydrophobicityProfile.cpp:122 `double hydrophobicity = residue.getHydrophobicity(HydrophobicityScaleMethod::EISENBERG);` with no scale parameter in the signature (header lines 87-91).
- **Fix:** Document in the header that the hydrophobic moment is computed on the Eisenberg consensus scale (per Eisenberg 1984), and optionally add a scale parameter overload. abi_impact none for doc; source-compatible if adding an overload.
- **Verifier correction:** The code is exactly as quoted, but the severity is low, not a notable API surprise. Hardcoding the Eisenberg scale for the hydrophobic moment is a correct domain convention (Eisenberg 1984 defines the moment on the Eisenberg consensus scale, matching the cited R.Peptides hmoment.R reference implementation), and the class-level documentation (header lines 35-41) already cites Eisenberg specifically for hydrophobic-moment computation. There is no silent-wrong-result risk because the method exposes no scale parameter to be ignored — the signature difference is visible at every call site. The legitimate, narrow fix is to add one sentence to the per-method Doxygen block (lines 79-86) stating the moment uses the Eisenberg scale; an optional scale-parameter overload would be source-compatible but is not required and arguably less faithful to the standard definition.
- **Verified:** The quoted evidence is accurate: HydrophobicityProfile.cpp:122 hardcodes `getHydrophobicity(HydrophobicityScaleMethod::EISENBERG)`, and the header signature (lines 87-91) has no scale parameter, while the sibling methods computeProfile (lines 60-63) and computeWindowedProfile (lines 73-77) do accept a HydrophobicityScaleMethod. So the surface-level inconsistency is real. However, the claim overstates it for two reasons. (1) Domain convention: the

### [CHEM-45] MzPAFParseError — MzPAFParseError stores context_before_/context_after_ but exposes no getters, unlike the sibling ProForma::ParseError
`low` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/MzPAF.h

```cpp
class MzPAFParseError : public Exception::ParseError { ... private: std::string context_before_, context_after_; }
```
- **Expectation:** Given ProForma::ParseError exposes getContextBefore()/getContextAfter()/getExpected()/getFound(), a developer using the parallel mzPAF parser would expect the same structured context accessors for programmatic error rendering.
- **Actual:** MzPAFParseError keeps context_before_ and context_after_ as private members populated by extractContext_, but provides only getErrorCode(), getPosition(), and getFormattedMessage() — no accessors to read the context fields. The structured context is captured but unreachable, an inconsistency with ProForma::ParseError in the same cluster.
- **Evidence:** MzPAF.h:222-232 declares getErrorCode/getPosition/getFormattedMessage and private context_before_/context_after_ with no getters; compare ProForma.h:627-645 getContextBefore/getContextAfter/getExpected/getFound.
- **Fix:** Add getContextBefore()/getContextAfter() accessors to MzPAFParseError for parity. abi_impact source-compatible (additive).
- **Verifier correction:** MzPAFParseError does store context_before_/context_after_ privately with no individual getters, unlike the sibling ProForma::ParseError which exposes getContextBefore()/getContextAfter()/getExpected()/getFound(). But the context is NOT unreachable: getFormattedMessage() (MzPAF.cpp:85-92) renders both fields into its message string. The surprise is the missing structured accessors for parity, not lost/inaccessible data — severity is low (loud, recoverable, no silent harm).
- **Verified:** Independently verified against source. MzPAF.h:222-233 confirms MzPAFParseError exposes only getErrorCode()/getPosition()/getFormattedMessage() and keeps context_before_/context_after_ private, populated by extractContext_ (MzPAF.cpp:67-83). The sibling ProForma::ParseError in the same module, built from an identical constructor signature, exposes getContextBefore()/getContextAfter()/getExpected()/getFound() (ProForma.h:627-648). So the asymmetri

### [CHEM-50] Residue::isInResidueSet — Read-only predicate isInResidueSet() is not const
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/CHEMISTRY/Residue.h

```cpp
bool isInResidueSet(const std::string& residue_set);
```
- **Expectation:** A boolean predicate named 'isInResidueSet' that only queries membership should be const, callable on a const Residue& (and the rest of the predicate block — operator==, isModified, hasNeutralLoss — are all const).
- **Actual:** isInResidueSet is declared non-const even though its body only does `residue_sets_.contains(residue_set)`. It cannot be called on a const Residue, which is the common case since ResidueDB hands out const Residue*. This forces awkward copies or const_casts at call sites.
- **Evidence:** Residue.h:428 `bool isInResidueSet(const std::string& residue_set);` (note: no trailing const, unlike every sibling predicate) and Residue.cpp:586-589 `bool Residue::isInResidueSet(const std::string& residue_set){ return residue_sets_.contains(residue_set); }` — pure read.
- **Fix:** Add `const` to the declaration and definition. Adding const is source-compatible for callers but changes the mangled symbol, so technically an ABI break for that one symbol; provide it as the canonical fix and rebuild dependents.
- **Verifier correction:** The defect is real exactly as quoted, but its impact is overstated. Calling isInResidueSet on a const Residite& fails at COMPILE TIME (a loud error), it does not silently produce wrong results — so this is low severity, not a silent-correctness hazard. The fix (add const to both declaration and definition) is correct and source-compatible for callers, but changes the mangled symbol name and is therefore an ABI break for that one symbol; dependents must be rebuilt.
- **Verified:** Independently verified against the actual code. Residue.h:428 declares `bool isInResidueSet(const std::string& residue_set);` with no trailing const — and it is the ONLY predicate in the "Predicates" block missing const; siblings hasNeutralLoss(), hasNTermNeutralLosses(), operator==, operator!=, and isModified() are all const. Residue.cpp:586-589 confirms the body is a pure read: `return residue_sets_.contains(residue_set);`. The member residue_s

### [CHEM-51] AASequence::getResidue / AASequence::operator[] — Doc comments say 'returns a pointer' but both return a const reference
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/AASequence.h

```cpp
const Residue& getResidue(Size index) const;  /  const Residue& operator[](Size index) const;
```
- **Expectation:** Reading 'returns a pointer to the residue at position index', a caller may treat the result as nullable / pointer-like.
- **Actual:** Both functions actually return `const Residue&` (a reference, never null; out-of-range throws IndexOverflow). The 'pointer' wording in the doc is wrong and could mislead a caller into pointer-style null checks or `&`/`*` mismatches.
- **Evidence:** AASequence.h:467-468 `/// returns a pointer to the residue at position @p index` over `const Residue& getResidue(Size index) const;` and AASequence.h:489-490 `/// returns a pointer to the residue at given position` over `const Residue& operator[](Size index) const;`. Impl AASequence.cpp:130-137 and 646-653 return `*peptide_[index]` (a reference) and throw on overflow.
- **Fix:** Fix the doc comments to 'returns a (const) reference to the residue ... throws Exception::IndexOverflow if out of range'. Doc-only, ABI-safe.
- **Verifier correction:** The doc comments for both AASequence::getResidue(Size) const and AASequence::operator[](Size) const incorrectly say "returns a pointer". Both return `const Residue&` (a reference, never null) and throw Exception::IndexOverflow when index is out of range. Doc-only fix; the misleading wording is immediately contradicted by the visible return type, so any pointer-style misuse fails to compile rather than producing wrong results.
- **Verified:** Evidence verified exactly. AASequence.h:467-468 and 489-490 do carry `/// returns a pointer to the residue ...` directly above `const Residue& getResidue(Size index) const;` and `const Residue& operator[](Size index) const;`. The implementations (AASequence.cpp:130-137 and 646-653) return `*peptide_[index]` (a reference) and throw Exception::IndexOverflow on out-of-range, never returning null. So the "pointer" wording is genuinely wrong and the d

### [CHEM-54] Residue::getLossNames — Doc says getLossNames() returns 'an empty string' but it returns a vector of strings
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/Residue.h

```cpp
const std::vector<std::string>& getLossNames() const;
```
- **Expectation:** The doc 'gets neutral loss name (if there is one, else returns an empty string)' implies a single string return.
- **Actual:** The function returns `const std::vector<std::string>&` (a list of loss names), and on 'no loss' returns an empty vector, not an empty string. Minor but the doc describes the wrong shape and cardinality of the return.
- **Evidence:** Residue.h:304-305 `/// gets neutral loss name (if there is one, else returns an empty string)` over `const std::vector<std::string>& getLossNames() const;`.
- **Fix:** Fix the doc to 'returns the neutral loss names (empty vector if none)'. Doc-only, ABI-safe.
- **Verifier correction:** Doc comment at Residue.h:304 is inaccurate: getLossNames() returns `const std::vector<std::string>&` (a list of neutral loss names), and returns an empty vector — not an empty string — when there are none. Fix doc to e.g. "returns the neutral loss names (empty vector if none)". The mismatch is cosmetic: the return type on the adjacent line is self-documenting and the compiler prevents any string-vs-vector misuse, so severity is low.
- **Verified:** Evidence verified. Residue.h:304-305 doc reads "gets neutral loss name (if there is one, else returns an empty string)" over `const std::vector<std::string>& getLossNames() const;`. The .cpp (Residue.cpp:190-193) returns the `loss_names_` member, a vector<string>; on "no loss" it is an empty vector, not an empty string. So the doc mis-describes both cardinality (singular "name" vs a list) and the empty-case shape ("empty string" vs empty vector).

### [CHEM-55] AASequence::getMonoWeight — getMonoWeight(type, charge) returns a proton-added neutral-plus-charge mass, not an m/z, despite a sibling getMZ existing
`low` · `unit-or-index` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/AASequence.h

```cpp
double getMonoWeight(Residue::ResidueType type = Residue::Full, Int charge = 0) const;
```
- **Expectation:** getMonoWeight is paired with getMZ(charge); a caller might assume getMonoWeight always yields the neutral monoisotopic weight regardless of the charge argument.
- **Actual:** When charge != 0, getMonoWeight adds PROTON_MASS_U * charge to the mass (the mass of the charged species, summed protons) but does NOT divide by charge. So getMonoWeight(Full, 2) is the [M+2H] mass, not the neutral mass and not the m/z. The charge parameter silently shifts the returned 'weight' by added protons.
- **Evidence:** AASequence.cpp:507-511 `double mono_weight(Constants::PROTON_MASS_U * charge);` then sums residues; getMZ at AASequence.cpp:498-504 divides getMonoWeight(type,charge)/charge.
- **Fix:** Header already says 'in the given ionic form'; strengthen the doc to state explicitly that a non-zero charge adds charge*proton mass (charged-species mass, not m/z, not neutral). Doc-only, ABI-safe.
- **Verifier correction:** getMonoWeight(type, charge) returns the [M+nH] charged-species monoisotopic mass (sum of neutral mass plus charge proton masses), not an m/z — confirmed by AASequence.cpp:511 (mono_weight init = PROTON_MASS_U*charge, never divided) and getMZ dividing this by charge. The behavior is already documented at AASequence.h:476 ("mono isotopic weight ... in the given ionic form"). The realistic neutral-mass use is the default charge=0; passing a non-zero charge to obtain a neutral mass is not a plausible caller assumption (all in-tree callers intentionally use the ionic [M+nH] value). The mild surprise is purely the "Weight"-vs-"MZ" naming next to the sibling getMZ. Recommendation is doc-only/ABI-safe: optionally state explicitly that a non-zero charge adds charge*proton mass (charged-species mass, not m/z, not neutral).
- **Verified:** Code confirms the mechanics: AASequence.cpp:511 initializes mono_weight = PROTON_MASS_U * charge (proton mass = 1.0072764667710 u, not neutron), sums residue/terminal masses, and returns the total [M+nH] charged-species mass WITHOUT dividing by charge; getMZ (line 504) divides getMonoWeight(type,charge)/charge. So the evidence is accurate. However the claim mis-states the surprise. (1) The header at line 476 already documents this accurately: "re

### [CHEM-56] NASequence::get — get(index) is a non-const getter and cannot be called on a const NASequence
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/CHEMISTRY/NASequence.h

```cpp
const Ribonucleotide* get(size_t index)
```
- **Expectation:** A pure read accessor named get() that only returns seq_[index] should be const, matching getSequence() which provides both const and non-const overloads. Callers expect to read elements from a const NASequence.
- **Actual:** get() is declared only as a non-const member (line 357-360, body `return seq_[index];`). It performs no mutation yet is not const, so `const NASequence& s; s.get(0);` fails to compile. The companion set() exists, but unlike getSequence() there is no const get overload.
- **Evidence:** const Ribonucleotide* get(size_t index)\n{\n  return seq_[index];\n}  // no `const` qualifier, body only reads
- **Fix:** Add `const` to the existing method (or add a const overload to preserve ABI): `const Ribonucleotide* get(size_t index) const`. The implementation already only reads, so marking it const is source-compatible for existing non-const callers.
- **Verifier correction:** get(size_t) is a non-const-only read accessor (NASequence.h:357), inconsistent with getSequence() and operator[] which both expose const overloads, so it cannot be called on a const NASequence. Real but low severity: a const read path already exists via `operator[](size_t) const` (line 368, returns the same seq_[index]), so const access is not blocked, and the get() limitation manifests as a hard compile error rather than any silent misbehavior. Recommended fix: add a const overload `const Ribonucleotide* get(size_t index) const` (preferred — keeps the existing symbol for binary compat) rather than mutating the existing method's signature.
- **Verified:** Verified in src/openms/include/OpenMS/CHEMISTRY/NASequence.h lines 357-360: `const Ribonucleotide* get(size_t index) { return seq_[index]; }` is non-const with a read-only body, and it is the ONLY get() overload — confirmed no const companion exists in the header or NASequence.cpp. This is a genuine const-correctness inconsistency: the sibling element accessors in the same class DO provide const overloads — getSequence() has both const/non-const 

### [CHEM-58] NASequence::getPrefix/getSuffix vs getSubsequence — length handling is inconsistent: getSubsequence clamps an over-long length, getPrefix/getSuffix throw
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/NASequence.h

```cpp
NASequence getPrefix(Size) const / NASequence getSubsequence(Size start=0, Size length=Size(-1)) const
```
- **Expectation:** Sibling length-based slicing methods in the same class should treat an out-of-range length the same way.
- **Actual:** getSubsequence silently clamps: `if (length > size() - start) length = size() - start;` (NASequence.cpp:115-116), so an over-long length is accepted. getPrefix/getSuffix instead throw IndexOverflow for length >= size() (NASequence.cpp:86,95). A caller who learned the clamping behavior from getSubsequence will be surprised by an exception from getPrefix/getSuffix.
- **Evidence:** getSubsequence: if (length > size() - start) length = size() - start;\ngetPrefix: if (length >= seq_.size()) throw Exception::IndexOverflow(...)
- **Fix:** Pick one convention (preferably clamp, consistent with getSubsequence) and apply it across all three, or document the divergence prominently. Aligning getPrefix/getSuffix to clamp is source-compatible for valid inputs.
- **Verifier correction:** Claim is accurate but understates one detail and overstates severity. Detail: getPrefix/getSuffix use `>=` (line 86, 95), so they throw even for length == size() (the full sequence), making them stricter than getSubsequence in two respects, not one. Severity: the inconsistency is real but the throwing path is LOUD and recoverable (Exception::IndexOverflow), not a silent wrong-result/data-loss hazard, so it is low, not medium. The clamping in getSubsequence is intentional/load-bearing (the Size(-1) default depends on it). Recommendation stands: pick one convention (clamp is the friendlier, source-compatible choice) and apply to all three, or document the divergence in the header.
- **Verified:** The code is exactly as quoted. getSubsequence (NASequence.cpp:115-116) silently clamps an over-long length: `if (length > size() - start) length = size() - start;` — and this clamping is load-bearing, since the header default `getSubsequence(start=0, length=Size(-1))` relies on it to return the whole sequence. By contrast getPrefix (line 86) and getSuffix (line 95) throw Exception::IndexOverflow on `length >= seq_.size()`. The divergence is genui

### [CHEM-59] NASequence::getSubsequence — getSubsequence (a const accessor) prints to std::cout
`low` · `hidden-side-effect` · ABI: `none` · src/openms/source/CHEMISTRY/NASequence.cpp

```cpp
NASequence getSubsequence(Size start = 0, Size length = Size(-1)) const
```
- **Expectation:** A const accessor that returns a subsequence should have no observable side effects beyond returning a value; it must not write to standard output.
- **Actual:** When the preceding residue is a phosphorothioate, the method executes `cout << seq_[start - 1]->getCode();` (NASequence.cpp:123) — an apparent leftover debug print that writes to stdout on every such call. This pollutes program output unexpectedly for a getter.
- **Evidence:** if (start > 0 && seq_[start - 1]->getCode().back() == '*' )\n{\n  cout << seq_[start - 1]->getCode();\n  static const RibonucleotideDB* rdb = RibonucleotideDB::getInstance();
- **Fix:** Remove the stray `cout` line. Pure cleanup, no API/ABI change.
- **Verified:** Verified independently. NASequence.cpp:109 defines getSubsequence as const (header line 474: `NASequence getSubsequence(Size start = 0, Size length = Size(-1)) const`). The file has `using namespace std;` (line 19), so line 123 `cout << seq_[start - 1]->getCode();` writes the preceding residue's code to std::cout. The branch fires whenever the residue at start-1 ends in '*' (phosphorothioate), which is a normal case during W/X fragment-ion genera

### [CHEM-60] Ribonucleotide::getAvgMass / setAvgMass — Doc comments for getAvgMass/setAvgMass are swapped (getter documented as 'Set', setter as 'Get')
`low` · `swapped-doc-comment` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/Ribonucleotide.h

```cpp
double getAvgMass() const / void setAvgMass(double avg_mass)
```
- **Expectation:** The Doxygen comment above an accessor should describe what that accessor does: 'Get ...' above the getter, 'Set ...' above the setter.
- **Actual:** The header swaps them: the comment `/// Set the average mass of the ribonucleotide` sits above `double getAvgMass() const` (returns avg_mass_), and `/// Get the average mass of the ribonucleotide` sits above `void setAvgMass(double avg_mass)` (assigns avg_mass_). The implementation (Ribonucleotide.cpp:88-96) confirms getAvgMass reads and setAvgMass writes, so only the docs are wrong — but they actively mislead a reader scanning the header.
- **Evidence:** /// Set the average mass of the ribonucleotide\ndouble getAvgMass() const;\n\n/// Set the average mass of the ribonucleotide  (line 110 above setAvgMass actually reads 'Get')\nvoid setAvgMass(double avg_mass);
- **Fix:** Swap the two doc comments so each describes its own method. Documentation-only fix, no ABI impact.
- **Verifier correction:** The doc comments above getAvgMass/setAvgMass are swapped: `/// Set the average mass...` (line 106) sits above the getter `double getAvgMass() const` and `/// Get the average mass...` (line 109) sits above the setter `void setAvgMass(double)`. The methods themselves are named and implemented correctly (cpp 88-96). This is a documentation-only swap (not a misleading symbol name); fix by exchanging the two comments. Low severity, no ABI impact.
- **Verified:** Verified in the header: line 106 `/// Set the average mass of the ribonucleotide` sits above the getter `double getAvgMass() const;` (line 107), and line 109 `/// Get the average mass of the ribonucleotide` sits above the setter `void setAvgMass(double avg_mass);` (line 110). The cpp confirms getAvgMass returns avg_mass_ (lines 88-90) and setAvgMass assigns it (lines 93-95), so the methods behave correctly and only the two doc comments are transp

### [CHEM-63] ModificationsDB::getNumberOfModifications — Doc claims count is "read from the unimod.xml file" but it counts all providers plus runtime-interned mods
`low` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ModificationsDB.h

```cpp
Size getNumberOfModifications() const
```
- **Expectation:** Header line 84: "Returns the number of modifications read from the unimod.xml file" implies the count reflects only the UniMod entries.
- **Actual:** It returns mods_.size(), which is populated from every provider (UniMod, custom_mods.xml, PSI-MOD.obo, XLMOD.obo) and additionally grows whenever addModification()/addNewModification_() interns a previously unseen modification at runtime. The number is therefore neither stable nor limited to unimod.xml.
- **Evidence:** ModificationsDB.cpp:107-115 returns `mods_.size()`. loadFromProviders_ (lines 529-584) pushes from all providers; addModification (lines 482, 507) and addNewModification_ (line 523) push more into mods_.
- **Fix:** Reword doc to "Returns the total number of modifications currently in the database (all providers plus any runtime-added)." ABI-safe doc fix.
- **Verifier correction:** The doc comment, not the symbol name, is inaccurate. Line 84 should read approximately: "Returns the total number of modifications currently in the database (aggregated across all providers: unimod.xml, custom_mods.xml, PSI-MOD.obo, XLMOD.obo, plus any runtime-added via addModification())." Re-categorized from misleading-name to misleading-doc (the name getNumberOfModifications is accurate); severity lowered to low (mild prose-only surprise, no incorrect results/crashes, API internally consistent).
- **Verified:** Evidence verified verbatim. Header ModificationsDB.h:84 reads "Returns the number of modifications read from the unimod.xml file". The implementation (ModificationsDB.cpp:107-115) returns mods_.size(). mods_ is populated by loadFromProviders_ (529-584) from ALL configured providers, and initializeModificationsDB (63-86) registers four by default: UniMod (unimod.xml), a second UnimodXMLDataProvider for custom_mods.xml, PSI-MOD.obo, and XLMOD.obo. 

### [CHEM-65] ModifiedPeptideGenerator::applyVariableModifications — Parameter named keep_original in header, keep_unmodified in impl; doc claims original is the "first entry" but output is appended to an existing vector
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h

```cpp
static void applyVariableModifications(const MapToResidueType& var_mods, const AASequence& peptide, Size max_variable_mods_per_peptide, std::vector<AASequence>& all_modified_peptides, bool keep_original=true)
```
- **Expectation:** Header declares the flag as keep_original and documents "also emit @p peptide unchanged as the first entry." A caller reads this as: when true, the unmodified peptide is the first element of all_modified_peptides.
- **Actual:** The implementation parameter is keep_unmodified (silent renaming risk for anyone using designated/aggregate reading or matching by name) and, crucially, all_modified_peptides is appended to ("existing contents are preserved" per the same header). If the vector is non-empty on entry, the unmodified peptide is NOT the first entry of the vector. "first entry" only holds for an initially-empty output.
- **Evidence:** Header lines 114-119 use keep_original and say "emit @p peptide unchanged as the first entry"; impl signature uses keep_unmodified (ModifiedPeptideGenerator.cpp:138) and the original is push_back/inserted into the possibly non-empty all_modified_peptides (cpp:144, 256-265).
- **Fix:** Make the header and impl parameter name agree, and change the doc to "...as the first appended entry". ABI-safe (parameter names are not part of ABI; doc-only).
- **Verifier correction:** The header declares the flag as keep_original (ModifiedPeptideGenerator.h:119) but both implementation signatures name it keep_unmodified (ModifiedPeptideGenerator.cpp:138 and :273) — a header/impl naming inconsistency. Separately, the header doc (h:112) says the unmodified peptide is emitted "as the first entry," which is imprecise given the output vector is explicitly appended-to with existing contents preserved (h:96-97, h:111; cpp:144, :256-265): it is only the first *appended* entry, not necessarily vector index 0. This is doc/naming-only; there is no behavioral bug or ABI impact. Fix: rename to agree and reword to "...as the first appended entry."
- **Verified:** Both factual evidence items check out. (1) Name mismatch is real: header (ModifiedPeptideGenerator.h:119, and also :112/:145/:150/:156) names the flag keep_original, while the impl signatures use keep_unmodified (ModifiedPeptideGenerator.cpp:138 for applyVariableModifications and :273 for applyAtMostOneVariableModification_). (2) The append semantics are real and explicitly documented: header lines 96-97 say "all_modified_peptides is appended to 

### [CHEM-69] CrossLinksDB::getAllSearchModifications — Override silently changes the 'search-eligible' criterion and drops thread-safety vs the base class with the same name
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/CHEMISTRY/CrossLinksDB.h

```cpp
void getAllSearchModifications(std::vector<std::string>& modifications) const
```
- **Expectation:** CrossLinksDB is documented to "expose the same lookup, iteration, and matching interface as the base class" (header lines 22-25). A caller of the inherited-looking getAllSearchModifications expects the same filter semantics and the same thread-safety as ModificationsDB::getAllSearchModifications.
- **Actual:** The base class filters by `getUniModRecordId() > 0` and sorts case-INsensitively under the OpenMS_ModificationsDB critical section. The CrossLinksDB override instead filters by `!getPSIMODAccession().empty()`, sorts case-SENSITIVELY, and iterates mods_ with NO omp critical lock. Same name, different inclusion set, different sort order, and a dropped lock that the base class deliberately holds while touching the mutable mods_.
- **Evidence:** Base: ModificationsDB.cpp:586-617 (UniModRecordId>0, omp critical, case-insensitive sort). Override: CrossLinksDB.cpp:32-44 (PSIMODAccession non-empty, plain sort(), no critical section).
- **Fix:** Document the differing criterion in the header explicitly, and add the `#pragma omp critical(OpenMS_ModificationsDB)` around the mods_ iteration to match the base class's thread-safety contract. Lock addition and doc are ABI-safe.
- **Verifier correction:** Narrow the finding to: CrossLinksDB::getAllSearchModifications (CrossLinksDB.cpp:32-44) iterates the mutable mods_ vector without the #pragma omp critical(OpenMS_ModificationsDB) guard that every other mods_ access in the hierarchy uses (e.g. base ModificationsDB.cpp:586-617). It is name hiding, not an override — the base method (ModificationsDB.h:226) is non-virtual, so there is no polymorphic-dispatch surprise. The differing filter (PSI-MOD accession vs UniMod record id) and case-sensitive sort are documented (header lines 77-88) and domain-required (XLMOD entries lack UniMod record IDs), so they are not surprises. Recommendation reduced to: add the omp critical section to match the base class's thread-safety convention (ABI-safe, source-compatible).
- **Verified:** Verified the code directly. Three sub-claims, only one survives. (1) FILTER CRITERION: The override filters by !getPSIMODAccession().empty() vs the base's getUniModRecordId()>0 — TRUE, but this is explicitly documented in the header (lines 77-78: "usable for identification searches (i.e. carry a PSI-MOD accession)") AND domain-required: CrossLinksDB loads CHEMISTRY/XLMOD.obo, whose cross-linkers carry PSI-MOD/XLMOD accessions, not UniMod record I
