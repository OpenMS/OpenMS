# OpenMS POLS Audit — Batch 7: QC + SYSTEM + APPLICATIONS + IONMOBILITY/INTERFACES/IMAGING

**Confirmed findings:** 43 · 4 high · 15 medium · 24 low.

### [QC-3] FragmentMassError::compute — Two compute() overloads compute the average/variance differently, yielding different results
`high` · `inconsistent-convention` · ABI: `none` · src/openms/source/QC/FragmentMassError.cpp · _qc-metrics-a_

```cpp
void compute(PeptideIdentificationList& pep_ids, const ProteinIdentification::SearchParameters& search_params, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit, double)
```
- **Expectation:** Both compute() overloads carry identical doc ('Stores average FME over all spectra and its variance') and should produce the same statistics for the same data.
- **Actual:** The FeatureMap overload accumulates all PPM errors, computes result.average_ppm ONCE after the full first pass (cpp:238), then runs a SECOND pass to compute variance against that final average (cpp:241). The PeptideIdentificationList overload instead computes result.average_ppm and calls calculateVariance_ INSIDE the per-pep_id loop (cpp:303-305): the average is recomputed on every iteration and variance is accumulated against a moving/partial average. calculateVariance_ also divides each contribution by counter_ppm (cpp:166), which keeps changing across iterations. The two overloads therefore produce different variance (and a wrong one in the list overload), despite identical contracts.
- **Evidence:** FeatureMap path FragmentMassError.cpp:237-241 'result.average_ppm = accumulator_ppm / counter_ppm; ... fmap.applyFunctionOnPeptideIDs(fVar);'. List path FragmentMassError.cpp:302-305 inside 'for (auto& pep_id : pep_ids)': 'result.average_ppm = accumulator_ppm / counter_ppm; calculateVariance_(result, pep_id, counter_ppm);'
- **Fix:** Restructure the list overload to mirror the FeatureMap overload: full accumulation pass, compute the final average once, then a second pass for variance. This is an implementation fix only; no signature/ABI change.

### [QC-8] PeptideMass::compute — PeptideMass annotates the OBSERVED mass from precursor m/z, not the theoretical mass the name/docs promise
`high` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/QC/PeptideMass.h · _qc-metrics-b_

```cpp
void compute(FeatureMap& features)
```
- **Expectation:** The class is documented as 'QC metric calculating theoretical mass of a peptide sequence' and the method @brief says 'Sets the mass metavalue to all PeptideHits by computing the theoretical mass'. A caller expects the theoretical neutral mass of the peptide SEQUENCE (e.g. AASequence::getMonoWeight()), independent of the measured precursor.
- **Actual:** The implementation computes the observed/experimental neutral mass from the precursor m/z: hit.setMetaValue("mass", (pi.getMZ() - Constants::PROTON_MASS_U) * hit.getCharge()). It never touches the peptide sequence at all, so the 'mass' value depends on measured m/z and charge, not on the theoretical sequence mass. For a mis-assigned/charge-wrong PSM the two differ; a caller reading 'mass' as theoretical mass gets a silent bug.
- **Evidence:** PeptideMass.cpp:23: hit.setMetaValue("mass", (pi.getMZ() - Constants::PROTON_MASS_U) * hit.getCharge()); vs header: '@brief QC metric calculating theoretical mass of a peptide sequence' / 'computing the theoretical mass'.
- **Fix:** Either rename/redocument the metavalue and method to reflect that it is the observed neutral mass derived from precursor m/z (source-compatible doc fix; metavalue rename is breaking), or change the implementation to compute the sequence theoretical mass via AASequence (behavior change). At minimum fix the doc to stop calling it 'theoretical'.

### [QC-11] QCBase::names_of_requires / QCBase::Requires — names_of_requires[] is missing the entry for Requires::ID, causing an out-of-bounds read in isRunnable()
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/QC/QCBase.h · _qc-metrics-b_

```cpp
static const std::string names_of_requires[]; enum class Requires : UInt64 { NOTHING, RAWMZML, POSTFDRFEAT, PREFDRFEAT, CONTAMINANTS, TRAFOALIGN, ID, SIZE_OF_REQUIRES }
```
- **Expectation:** names_of_requires[] is documented as 'strings corresponding to enum Requires' and is indexed by the enum value (QCBase::names_of_requires[i] for i in [0, SIZE_OF_REQUIRES)). A caller/maintainer expects one string per enum member up to (but excluding) SIZE_OF_REQUIRES.
- **Actual:** The enum has 7 real members (NOTHING=0 .. ID=6), but names_of_requires only contains 6 strings: {"fail","raw.mzML","postFDR.featureXML","preFDR.featureXML","contaminants.fasta","trafoAlign.trafoXML"}. QCBase::isRunnable() loops 'for (i=0; i < (UInt64)SIZE_OF_REQUIRES; ++i)' and reads names_of_requires[i]; when a metric requires ID (i=6) the warning path reads names_of_requires[6] — past the end of the 6-element array (undefined behavior). The ID requirement also has no human-readable name.
- **Evidence:** QCBase.cpp:16 array has 6 elements; QCBase.h enum lists NOTHING..ID (7 members) before SIZE_OF_REQUIRES; QCBase.cpp:69-73 loop indexes names_of_requires[i] with i up to SIZE_OF_REQUIRES-1 (=6).
- **Fix:** Add the missing 'ID' string (and an entry for the NOTHING/SIZE bookkeeping if intended) so the array length matches the enum; keep it append-only to preserve indices. Source-compatible additive fix.
- **Verifier correction:** Confirmed as stated. Minor clarification: the OOB read occurs only on the diagnostic warning path inside isRunnable() (the function still correctly returns false), and it is reachable via IdentificationSummary (the only metric requiring Requires::ID, used by MzQCFile) when the ID input is missing. It reads a std::string at past-the-end memory and streams it to the log (UB: garbage output, potential crash, or info leak). Fix: append one string (e.g. "id.idXML") to names_of_requires[] so the array length matches SIZE_OF_REQUIRES; append-only preserves existing indices and is source/ABI-compatible.

### [SYST-6] File::getPathLocations / File::executableExtensions_ — Default argument std::getenv("PATH") constructs a std::string from a possibly-null pointer (UB/crash) when PATH is unset
`high` · `surprising-crash` · ABI: `source-compatible` · src/openms/include/OpenMS/SYSTEM/File.h · _sys-core_

```cpp
static StringList getPathLocations(const std::string& path = std::getenv("PATH"));
```
- **Expectation:** Calling getPathLocations() with no argument reads the PATH environment variable and returns its components; if PATH is unset it should yield an empty list, not crash.
- **Actual:** The default argument is 'std::getenv("PATH")', whose result is implicitly converted to std::string. If the environment variable is not set, getenv returns nullptr and constructing std::string(nullptr) is undefined behavior (typically a crash). The same pattern is used for the Windows %PATHEXT% default in executableExtensions_(const std::string& ext = std::getenv("PATHEXT")).
- **Evidence:** File.h:267 'static StringList getPathLocations(const std::string& path = std::getenv("PATH"));' and File.h:357 'static StringList executableExtensions_(const std::string& ext = std::getenv("PATHEXT"));'.
- **Fix:** Change the default to a const char* / nullptr-safe wrapper (e.g. default to "" and resolve getenv internally with a null check), so an unset PATH/PATHEXT yields an empty result instead of UB. The public getPathLocations is the ABI-relevant one; a null-guard inside the body plus a safe default is source-compatible.
- **Verifier correction:** The category should be 'surprising-crash'/UB rather than 'surprising-throw': constructing std::string from the null pointer returned by std::getenv when PATH/PATHEXT is unset is undefined behavior (a crash), not a thrown C++ exception. Everything else in the claim holds. Reachable via File::findExecutable -> getPathLocations()/executableExtensions_() (File.cpp:863,869), called from TOPPBase.cpp:1453 and PythonInfo.cpp:48. Fix: change defaults to "" (or a null-safe wrapper) and resolve getenv with a nullptr check inside the body; this is source-compatible for the public getPathLocations (only the default-argument expression changes, the mangled symbol and parameter type are unchanged).

### [APPL-1] TOPPBase::getOutputDirOption — getOutputDirOption() creates the directory on disk as a side effect of a getter
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/APPLICATIONS/TOPPBase.h · _app-topp_

```cpp
std::string getOutputDirOption(const std::string& name) const
```
- **Expectation:** A method named get...Option that returns a std::string and whose Doxygen lists only the four parameter-resolution exceptions reads as a pure accessor: it resolves and returns the value of a registered output-dir parameter.
- **Actual:** Besides resolving the value it calls File::makeDir(tmp), creating the directory tree on the filesystem. The header documentation (lines 722-730) does not mention any directory creation at all; only UnregisteredParameter/RequiredParameterNotGiven/WrongParameterType/InvalidParameter are documented.
- **Evidence:** TOPPBase.cpp:1349-1355: `std::string tmp = getParamAsString_(name, p.default_value.toString()); writeDebug_(...); // create directory if it does not exist\n    File::makeDir(tmp); return tmp;`
- **Fix:** Document the directory-creation side effect in the header (and that it may throw if creation fails), or split into a pure getOutputDirOption() plus an explicit ensureOutputDir_()/createOutputDir_() helper. Documenting is source- and ABI-compatible; splitting is additive (add new method, keep old). At minimum the surprise must be in the declared contract.
- **Verifier correction:** getOutputDirOption() does create the directory via File::makeDir() as an undocumented side effect (TOPPBase.cpp:1353), but it does NOT throw on creation failure: File::makeDir returns a bool using std::error_code and the return value is ignored, so failed directory creation is silent. The header contract (TOPPBase.h:722-730) should document the directory-creation side effect AND note that failures are currently swallowed; or split into a pure getter plus an explicit createOutputDir_() helper that checks/propagates the makeDir() result.

### [APPL-3] TOPPOpenSwathBase::loadSwathFiles — loadSwathFiles returns bool 'success' but most real failures throw instead of returning false
`medium` · `return-value` · ABI: `none` · src/openms/include/OpenMS/APPLICATIONS/OpenSwathBase.h · _app-topp_

```cpp
bool loadSwathFiles(..., const bool force, ..., Interfaces::IMSDataConsumer* plugin_consumer = nullptr)
```
- **Expectation:** @return 'Returns whether loading and sanity check was successful' implies callers can branch on the bool to detect any load/validation failure.
- **Actual:** The bool is false ONLY for the m/z window gap/overlap sanity check (and only when -force is off). All other failures throw: load errors propagate from loadSwathFiles_, and mixed-file-type detection throws Exception::IllegalArgument. A caller that only checks `if (!loadSwathFiles(...))` will miss thrown errors, and one that relies on the bool covering 'sanity check' will be surprised that mixed-type inconsistency is an exception, not a false.
- **Evidence:** OpenSwathBase.cpp:262 throws IllegalArgument for mixed types; only lines 318 and 333 `return false;` (gap/overlap with !force); otherwise `return true;` at line 341.
- **Fix:** Clarify the header @return to state explicitly that the bool reflects only the window gap/overlap check and that other errors (load failures, mixed input types) are reported via exceptions. Documentation-only; ABI-neutral.
- **Verifier correction:** The bool return of loadSwathFiles reflects only the m/z SWATH-window gap/overlap sanity check, and only when -force is off (false at OpenSwathBase.cpp:318 and :333). It does NOT signal load failures (those throw out of loadSwathFiles_) and does NOT signal the mixed-input-file-type inconsistency, which throws Exception::IllegalArgument at line 262 despite being part of the documented "sanity check." The @return "Returns whether loading and sanity check was successful" overstates the bool's coverage. Failures that bypass the bool surface as exceptions and fail loudly (non-zero TOPP exit via TOPPBase::main), so this is a documentation/contract-clarity defect (medium), not a silent-failure bug. Recommend clarifying the header @return to state the bool reflects only the window gap/overlap check and that load errors and mixed-type inputs are reported via exceptions. Documentation-only, ABI-neutral.

### [QC-4] FragmentMassError::Statistics / getResults — Zero-valued Statistics is pushed on no-data, indistinguishable from a real measurement of 0 ppm
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/QC/FragmentMassError.h · _qc-metrics-a_

```cpp
struct Statistics { double average_ppm = 0; double variance_ppm = 0; }; const std::vector<Statistics>& getResults() const
```
- **Expectation:** When compute() cannot calculate any fragment mass error (no pep IDs, or zero matched peaks), the caller should be able to tell the result is 'not computed' rather than 'measured 0 ppm error'.
- **Actual:** On the no-data paths compute() pushes a default-constructed Statistics{average_ppm=0, variance_ppm=0} onto results_ and returns silently (cpp:176-180 no pepIDs; cpp:231-235 counter_ppm==0; cpp:251-255 empty pep_ids; cpp:296-301). A consumer reading getResults() sees average_ppm==0, which is also a perfectly valid measured value (a spectrum with zero net ppm error), so genuine failure is masked as a flawless score.
- **Evidence:** FragmentMassError.cpp:231-235 'if (counter_ppm == 0) { results_.push_back(result); return; }' where 'Statistics result;' has average_ppm=0, variance_ppm=0.
- **Fix:** Initialize Statistics fields to NaN (additive, source-compatible since the struct stays the same size/layout) or add a count/valid flag so callers can distinguish 'no data' from 'zero error'. Document the no-data behavior either way.
- **Verifier correction:** The collision is real on all four no-data paths, but the practical impact is medium not high: the shipped consumer (TOPP QualityControl, FragmentMassError.cpp callers) does not read getResults() to gate decisions and the class is unbound in pyOpenMS, so the failure mode is a misleading QC value for a direct getResults() reader rather than silently wrong results in the main pipeline. Recommendation stands: initialize fields to NaN (layout-preserving, source-compatible) or add a valid/count flag, and document the no-data behavior in the header.

### [QC-6] Contaminants::ContaminantsSummary::empty_features — empty_features.first counts empty PeptideIdentifications, not features, and can exceed the total
`medium` · `misleading-name` · ABI: `none` · src/openms/source/QC/Contaminants.cpp · _qc-metrics-a_

```cpp
std::pair<Int64, Int64> empty_features; ///< (Number of features without a peptide hit, total number of features)
```
- **Expectation:** Per the field doc, empty_features.first = number of features lacking a peptide hit, and .second = total number of features, so first <= second.
- **Actual:** The counter feature_has_no_sequence is incremented once per feature with no PeptideIdentifications (cpp:80) BUT ALSO once per PeptideIdentification whose hit list is empty (cpp:87-89, inside the per-id loop). A single feature carrying several empty PeptideIdentifications therefore increments the counter multiple times, so empty_features.first counts empty-ID events, not empty features, and can exceed empty_features.second (= features.size(), cpp:103). The documented invariant first<=second is violated and the field's name/meaning is wrong.
- **Evidence:** Contaminants.cpp:83-90 'for (auto& id : f.getPeptideIdentifications()) { if (id.getHits().empty()) { ++feature_has_no_sequence; continue; } ... }' combined with cpp:78-81 '++feature_has_no_sequence' for empty feature; cpp:102-103 'final.empty_features.first = feature_has_no_sequence; final.empty_features.second = features.size();'.
- **Fix:** Either rename/redefine the field to reflect 'features or PSMs without a usable sequence' and fix the doc, or count distinct features (break out of the inner loop after the first empty id, or use a per-feature flag). Implementation/doc change only; struct layout unchanged, no ABI impact.
- **Verifier correction:** Severity adjusted to medium (claim implied higher). The affected value is a reported QC summary metric (empty_features.first), not a result that drives downstream computation, crashes, or corrupts data; the pathological over-count requires the relatively uncommon case of a feature carrying multiple empty-hit PeptideIdentifications, and in the typical case first<=second still holds. It is nonetheless a silently-wrong/mislabeled reported number whose name and documentation do not match its computation. abi_impact = none for the substantive fix (count distinct features via break/per-feature flag and/or correct the doc) — struct layout is unchanged; only a field rename (one possible remediation) would be source-breaking, but that is not required.

### [QC-7] IdentificationSummary::compute — compute() returns NaN means or throws on empty/partial input with no documented contract
`medium` · `silent-failure` · ABI: `none` · src/openms/source/QC/IdentificationSummary.cpp · _qc-metrics-a_

```cpp
Result compute(std::vector<ProteinIdentification>& prot_ids, PeptideIdentificationList& pep_ids)
```
- **Expectation:** A 'compute a summary of an idXML file' entry point should handle an empty or peptide-only input gracefully, or document that non-empty input is required.
- **Actual:** Several unguarded divisions silently yield NaN: peptide_length_mean = sum / peptides.size() when no peptides (cpp:43), missed_cleavages_mean = missed_cleavages / pep_count when pep_count==0 (cpp:57), protein_hit_scores_mean = sum / protein_hit_count when there are no protein hits (cpp:71). Additionally pep_ids.front() (cpp:76) and prot_ids.front() (cpp:80) are called unconditionally and dereference the front of a possibly-empty container (UB/throw). The header gives no precondition. A caller passing an empty or protein-less idXML gets NaN fields or a crash rather than a clear error.
- **Evidence:** IdentificationSummary.cpp:43 'result.peptide_length_mean = (double)peptide_length_sum / peptides.size();'; cpp:57 '... / pep_count;'; cpp:71 '... / protein_hit_count;'; cpp:76 'if (pep_ids.front().getScoreType() == "FDR")'; cpp:80 'if (prot_ids.front().getScoreType() == "FDR")'.
- **Fix:** Guard the divisions (return 0 or NaN deliberately and documented) and the front() calls (check empty()), and document the empty-input behavior. Pure implementation change, no ABI impact.
- **Verifier correction:** The unguarded divisions (cpp:43, 57, 71) are the real defect: they silently produce NaN, reachable even via the only legitimate caller (MzQCFile::store), because its `!prot_ids.empty() && !pep_ids.empty()` gate does not guarantee non-zero peptide-hits, missed-cleavage counts, or protein-hits — those NaNs flow into the mzQC JSON with no error or warning. The front() calls (cpp:76, 80) are undefined behavior on empty containers (noexcept forwarding to std::vector::front(), NOT a throw as the title implies), but through the in-repo caller they are unreachable on empty input due to that gate; they are only a hazard for a direct external API caller, since the public header documents no precondition. Recommended fix unchanged: guard the three divisions against zero denominators and check empty() before the two front() calls, and document the empty-input contract. Pure .cpp implementation change — no ABI impact.

### [QC-9] PeptideMass::compute — Doc says 'all PeptideHits' get annotated, but only the first hit per PeptideIdentification is touched
`medium` · `doc-contract-mismatch` · ABI: `none` · src/openms/include/OpenMS/QC/PeptideMass.h · _qc-metrics-b_

```cpp
void compute(FeatureMap& features)
```
- **Expectation:** Header: 'Each PeptideHit in the FeatureMap will be annotated with its theoretical mass as metavalue mass' and '@brief Sets the mass metavalue to all PeptideHits'. A caller expects every PeptideHit (all rank candidates) to carry a 'mass' metavalue.
- **Actual:** Only getHits()[0] is annotated: 'auto& hit = pi.getHits()[0]; hit.setMetaValue("mass", ...)'. Hits 2..n silently get no 'mass' metavalue, so a caller iterating all hits and reading 'mass' will hit missing-metavalue on lower-ranked hits.
- **Evidence:** PeptideMass.cpp:22-23: auto& hit = pi.getHits()[0]; hit.setMetaValue("mass", ...); — only index 0, despite header 'all PeptideHits'.
- **Fix:** Fix the doc to say 'first PeptideHit of each PeptideIdentification' (source-compatible), or loop over all hits if all were intended.
- **Verifier correction:** PeptideMass::compute annotates only the first PeptideHit (getHits()[0]) of each PeptideIdentification, contradicting the header docs (PeptideMass.h:20 and :33) which claim ALL PeptideHits get the 'mass' metavalue. Lower-ranked hits are left without it. A caller reading hit.getMetaValue("mass") on those hits via the single-arg overload throws (loud), while the default-value overload silently returns the default — so the effect is 'loud-or-silent surprise' rather than purely silent. Fix the documentation to say 'the first PeptideHit of each PeptideIdentification' (source-compatible), or loop over all hits if full annotation was intended.

### [QC-10] TIC::getResults — TIC::getResults() is declared but never defined, and its type is inconsistent with compute()'s actual result
`medium` · `return-value` · ABI: `source-compatible` · src/openms/include/OpenMS/QC/TIC.h · _qc-metrics-b_

```cpp
const std::vector<MSChromatogram>& getResults() const
```
- **Expectation:** A public getResults() on a QC metric, mirroring sibling metrics (Ms2IdentificationRate::getResults, PSMExplainedIonCurrent::getResults), should return the stored results and link. The return type std::vector<MSChromatogram> implies the metric retains chromatograms.
- **Actual:** There is NO definition of TIC::getResults anywhere in the source tree (grep for 'TIC::getResults' returns nothing). Calling it produces an unresolved-symbol link error. Moreover compute() returns a TIC::Result struct (intensities/relative_intensities/retention_times/area/fall/jump), and the class doc explicitly states 'results are returned to the caller and not retained by the instance' — so no vector<MSChromatogram> is ever populated. The declaration is a dead, mistyped landmine.
- **Evidence:** TIC.h:83 declares 'const std::vector<MSChromatogram>& getResults() const;'; grep 'TIC::getResults' over src/ finds zero definitions; TIC.h:28-29 doc: 'results are returned to the caller and not retained by the instance.'
- **Fix:** Remove the getResults() declaration from TIC.h (the metric is stateless by design). Removing an unimplemented member is effectively source-compatible since no caller can currently link against it.
- **Verifier correction:** Severity is medium, not high: any attempt to call TIC::getResults() fails at LINK time (unresolved symbol), which is loud and compile/link-time, not a silent wrong-result or runtime corruption. The declaration is dead, undefined, and mistyped (returns std::vector<MSChromatogram>& while compute() actually yields a TIC::Result by value; the class retains no chromatograms). Recommended fix — remove the declaration — is source-compatible: since no definition exists, no caller can currently link against it, and it is a non-virtual member so no vtable/ABI dependency exists.

### [QC-13] PSMExplainedIonCurrent::compute — Default tolerance=20 is unit-blind: passing tolerance_unit=DA without tolerance gives a 20 Dalton match window
`medium` · `surprising-default` · ABI: `source-compatible` · src/openms/include/OpenMS/QC/PSMExplainedIonCurrent.h · _qc-metrics-b_

```cpp
void compute(FeatureMap&, const MSExperiment&, const QCBase::SpectraMap&, ToleranceUnit tolerance_unit = ToleranceUnit::AUTO, double tolerance = 20)
```
- **Expectation:** The default 'tolerance = 20' reads as a sensible ppm fragment tolerance. A caller who explicitly selects ToleranceUnit::DA but leaves tolerance at its default would expect a sane Da window (or for the default to track the unit).
- **Actual:** The default value 20 is shared across units. If tolerance_unit != AUTO (e.g. DA) the AUTO auto-detection block is skipped, so 'tolerance' stays 20 and is fed directly to a Da-based MatchedIterator (DaTrait) — a 20 Dalton fragment match window, which matches essentially everything and silently inflates the explained-ion-current score. Nothing clamps or warns.
- **Evidence:** PSMExplainedIonCurrent.h:61 default 'tolerance = 20'; PSMExplainedIonCurrent.cpp:112-116 'if (tolerance_unit == ToleranceUnit::DA) { ... MIV mi(theo_spectrum, filtered_spec, tolerance); ... }' with the AUTO override only entered when tolerance_unit==AUTO (cpp:156).
- **Fix:** Make the default depend on the unit (e.g. require an explicit tolerance when DA is chosen), or document that the 20 default is ppm-only and a Da tolerance must always be passed. Doc-only change is source-compatible.
- **Verifier correction:** The mechanism and evidence are exactly as claimed. Only the implied impact level is adjusted down to medium: the unit-blind default 20 silently produces a 20 Da fragment window when ToleranceUnit::DA is selected without an explicit tolerance, but no in-tree caller hits this (all use AUTO or pass DA with an explicit tolerance), so it invites misuse on the public API rather than breaking shipped behavior. Recommended fix (require/guard an explicit tolerance when DA is chosen, or document that 20 is a ppm-only default) is source-compatible and does not change the ABI.

### [QC-16] Ms2IdentificationRate::compute — compute() appends to the results vector with no reset; repeated calls accumulate, and getResults() has no clear()
`medium` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/QC/Ms2IdentificationRate.h · _qc-metrics-b_

```cpp
void compute(const FeatureMap&, const MSExperiment&, bool assume_all_target = false)
```
- **Expectation:** A getResults() paired with compute() typically reflects the result of the most recent compute(), or there is a documented way to reset. A caller calling compute() twice (e.g. on two files) expects either replacement or a clear API.
- **Actual:** compute() -> writeResults_() always does 'rate_result_.push_back(id_rate_data)'; nothing ever clears rate_result_, and there is no public clear/reset. Two compute() calls leave getResults() with two entries; getResults()[0] is the first run, not the latest. This is intentional for multi-file batching but is unguarded and easy to misuse (stale accumulation across reuse of one instance).
- **Evidence:** Ms2IdentificationRate.cpp:147 'rate_result_.push_back(id_rate_data);' with no clear; header getResults() at Ms2IdentificationRate.h:114 returns the whole accumulated vector.
- **Fix:** Document the accumulating semantics on compute()/getResults() and/or add a public clear(). Additive (add clear()) is ABI-safe; doc clarification is source-compatible.

### [QC-20] MQEvidence::exportFeatureMap / MQMsms::exportFeatureMap — Sibling exportFeatureMap() signatures disagree on whether MSExperiment is optional
`medium` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/QC/MQEvidenceExporter.h · _qc-mq-export_

```cpp
void MQEvidence::exportFeatureMap(..., const MSExperiment& exp = {}, ...) vs void MQMsms::exportFeatureMap(..., const MSExperiment& exp, ...)
```
- **Expectation:** The two MaxQuant exporters (MQEvidence, MQMsms) are documented as parallel classes with the same exportFeatureMap entry point, so a caller would expect identical optionality of the same parameters.
- **Actual:** MQEvidence::exportFeatureMap defaults exp to {} (optional), while MQMsms::exportFeatureMap makes exp required (no default), even though prot_map has a default in both. A caller porting code from one exporter to the other gets a surprising compile error or, worse, silently passes an empty MSExperiment to MQEvidence and loses MS/MS-derived columns (msms_mz, base_peak_fraction) with no signal.
- **Evidence:** MQEvidenceExporter.h: 'const OpenMS::MSExperiment& exp= {}, const std::map<...>& prot_map = {}'. MQMsmsExporter.h: 'const OpenMS::MSExperiment& exp, const std::map<...>& prot_map = {}'.
- **Fix:** Make the two signatures consistent. Either require exp in both or default it in both. Adding a default to MQMsms is source-compatible; removing the default from MQEvidence would be a breaking source change, so prefer aligning toward the optional form.

### [IONM-4] IonImage::getIntensity — Returns 0.0 for never-written pixels, indistinguishable from a true zero intensity
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/IMAGING/IonImage.h · _small-im-interfaces-imaging_

```cpp
double getIntensity(UInt x, UInt y) const
```
- **Expectation:** Reading intensity at a pixel either returns a meaningful value or signals 'no data' for masked-out cells.
- **Actual:** An in-bounds pixel that was never set (hasPixel == false) returns 0.0, the same value a genuinely-measured zero would return. Callers that sum/average over getData() without consulting getMask()/hasPixel() silently treat absent pixels as measured zeros.
- **Evidence:** IonImage.cpp lines 46-49 return intensities_[idx] unconditionally; resize() (lines 21-28) zero-fills intensities_ and false-fills mask_. Header line 70 documents the 0.0-if-never-set behavior.
- **Fix:** Behavior is documented and ABI-stable; keep getIntensity but make the mask coupling unmissable in docs and consider an additive helper like `std::optional<double> getIntensityIfPresent(x,y)`. Additive only; no ABI break.
- **Verifier correction:** The named symbol getIntensity is itself documented ("0.0 if hasPixel is false", header lines 67/70), so a developer reading its doc is warned — the genuine, under-documented trap is bulk aggregation over getData(), whose docstring omits any mention of masked-out zeros or getMask(). Combined with the producer extractIonImage leaving non-acquired in-bounds cells unset (so absent pixels really do occur), summing/averaging getData() silently counts absent cells as measured zeros. This is real but is medium severity (recoverable via the parallel mask, partially documented), not a high silently-wrong-results defect. Module is IMAGING, not IONMOBILITY.

### [SYST-2] File::fileSize — fileSize returns -1 on error but the return type is unsigned, so the sentinel becomes the maximum UInt64
`medium` · `return-value` · ABI: `source-compatible` · src/openms/include/OpenMS/SYSTEM/File.h · _sys-core_

```cpp
static UInt64 fileSize(const std::string& file);
```
- **Expectation:** On error a caller can detect failure; a documented '-1 on error' return implies a negative/distinguishable sentinel.
- **Actual:** The function returns UInt64 and on error does 'return -1;' (File.cpp:169,173), which wraps to 0xFFFFFFFFFFFFFFFF (~1.8e19). A caller doing 'if (fileSize(f) > threshold)' or summing sizes silently treats a missing/unreadable file as an enormous file rather than an error. Only an exact '== (UInt64)-1' check works, which the header text ('-1 on error') does not make obvious for an unsigned type.
- **Evidence:** File.h:81 '/// The filesize in bytes (or -1 on error, e.g. if the file does not exist)' and 'static UInt64 fileSize(...)'; File.cpp:169 'if (!File::exists(file)) return -1;' and File.cpp:173 'if (ec) return -1;'.
- **Fix:** Either return a signed Int64 (ABI-breaking) or, additively, document the exact sentinel value 'std::numeric_limits<UInt64>::max()' and add a named constant so callers compare correctly. At minimum, change the doc to say the sentinel is the maximum UInt64, not '-1'.
- **Verifier correction:** The claim is accurate. Severity is medium rather than high: the misleading "-1 on error" sentinel on an unsigned UInt64 return is a real footgun (wraps to ~1.8e19, defeating naive `>`/sum checks), and the two existing callers (MzMLSplitter.cpp:98, MzMLHandler.cpp:1411) do not check it, but both consume already-validated input files so the wrong path is reachable yet uncommon and recoverable, not a silent-corruption/crash class. Primary recommended fix is additive and source-compatible: correct the doc to state the sentinel is `std::numeric_limits<UInt64>::max()` and expose a named constant for callers to compare against. Switching the return type to a signed Int64 would be the only ABI-breaking option and is not required.

### [SYST-3] SysInfo::MemUsage::delta / SysInfo::MemUsage::diff_str_ — delta() never reports a memory decrease: the negative sign is computed then discarded
`medium` · `silent-failure` · ABI: `none` · src/openms/source/SYSTEM/SysInfo.cpp · _sys-core_

```cpp
std::string delta(const std::string& event = "delta"); std::string diff_str_(size_t mem_before, size_t mem_after);
```
- **Expectation:** delta() reports the signed difference between the two timepoints (before -> after), so a memory drop shows as a negative value.
- **Actual:** diff_str_ first appends a '-' to s when mem_after < mem_before, but the very next statement reassigns s with '=' (not '+=') to the abs() value plus ' MB', discarding the sign. delta() therefore always reports a non-negative magnitude, so a real decrease is indistinguishable from an increase of the same magnitude.
- **Evidence:** SysInfo.cpp:326-332: 'std::string s;\n if (mem_after < mem_before) { s +=std::string("-"); }\n s =StringUtils::toStr(std::abs(((long long)mem_after - (long long)mem_before) / 1024)) + " MB";' — the conditional '-' written to s is overwritten by the subsequent assignment.
- **Fix:** Build the magnitude string first, then prepend the sign (use '+=' or compose 'sign + magnitude'). Internal helper fix; no API change.

### [SYST-8] Network::downloadFile — downloadFile returns void and the actual on-disk filename is unknowable to the caller (auto-renamed to .0/.1 on collision)
`medium` · `return-value` · ABI: `none` · src/openms/include/OpenMS/SYSTEM/Network.h · _sys-net-lang_

```cpp
static void downloadFile(const std::string& url, const std::string& download_folder)
```
- **Expectation:** A caller of `downloadFile(url, folder)` expects to know where the file landed so it can subsequently open/process it. Since the destination filename is derived from the URL and silently changed to '<name>.0', '<name>.1', ... when a file already exists, the caller cannot reliably reconstruct the path it must read next.
- **Actual:** The method returns `void`. The chosen filename is computed internally by `saveFileName_()` (strip ?query/#fragment, take basename, fall back to "download", then append .0/.1/... to avoid overwriting) and is only emitted via OPENMS_LOG_INFO ('Stored as ...'). There is no return value or out-parameter conveying the final path.
- **Evidence:** Network.cpp:50-76: `void Network::downloadFile(...)` ... `std::string filename = folder + "/" + saveFileName_(url, folder);` ... `OPENMS_LOG_INFO << "Stored as '" << filename << "'."`. saveFileName_ (lines 43-46): `while (fs::exists(dest / (basename + "." + std::to_string(i)))) ++i; return basename + "." + std::to_string(i);`
- **Fix:** Additive, ABI-safe fix: add an overload `static std::string downloadFile(const std::string& url, const std::string& download_folder)` (or a new `downloadFileTo`) that returns the final on-disk path. Keep the existing void signature for compatibility. Document that the void overload's destination filename is non-deterministic under collisions.
- **Verifier correction:** Severity is medium, not high. The collision behavior never overwrites existing files (no data loss), the final path is logged loudly, and a caller that guesses the wrong path fails loudly at file-open time (recoverable), so it invites misuse but does not produce silently-wrong results. The header already fully documents the .0/.1 renaming algorithm, so the surprise is narrowly the missing programmatic return of the path, not hidden behavior. The recommended additive fix must NOT be a return-type-only overload of downloadFile (ill-formed in C++); use a distinctly-named method (e.g. std::string downloadFileTo(url, folder)). The current sole caller is the unit test (Network_test.cpp), which only exercises the no-collision path.

### [SYST-9] JavaInfo::canRun / PythonInfo::canRun — Sibling canRun() methods (explicitly cross-referenced as analogous) have incompatible signatures and error-reporting models
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/SYSTEM/JavaInfo.h · _sys-net-lang_

```cpp
JavaInfo: static bool canRun(const std::string&, bool verbose_on_error = true)  vs  PythonInfo: static bool canRun(std::string&, std::string& error_msg)
```
- **Expectation:** Both headers explicitly say they are analogous ('Similar classes exist for other external tools, e.g. PythonInfo/JavaInfo'). A developer who learned one canRun() expects the other to follow the same convention: same parameter mutability and the same way of obtaining the error description.
- **Actual:** JavaInfo::canRun takes the executable by const-ref, does NOT mutate it, and reports errors by *logging* (controlled by a `verbose_on_error` bool). PythonInfo::canRun takes the executable by non-const ref, MUTATES it to an absolute path, and reports errors via an `error_msg` out-parameter (never logs). The two cannot be used interchangeably and have opposite mutation semantics for the same conceptual argument.
- **Evidence:** JavaInfo.h:36 `static bool canRun(const std::string& java_executable, bool verbose_on_error = true);` ; PythonInfo.h:38 `static bool canRun(std::string& python_executable, std::string& error_msg);` ; both class docs cross-reference each other (JavaInfo.h:20, PythonInfo.h:20).
- **Fix:** Document the divergence prominently in both headers. For convergence without breaking ABI, add overloads so each class offers both the logging form and the error_msg form. Avoid changing the existing signatures (ABI break).
- **Verifier correction:** The finding is correct in substance but its severity is overstated. The two signatures are incompatible enough that they cannot be accidentally interchanged silently — a `const std::string&` won't bind where a non-const `std::string&` is required, and the trailing parameters differ in type (bool verbose_on_error vs std::string& error_msg), so attempting to use one in place of the other yields a compile error, not silently wrong results or data loss. Hence medium (invites confusion / mild surprise via PythonInfo's in-place path mutation, but loud and recoverable), not high. The recommendation stands: document the divergence in both headers (especially since they cross-reference each other as analogous) and optionally add additive overloads to converge without an ABI break. abi_impact of the existing inconsistency is none.

### [APPL-2] SearchEngineBase::getRawfileName — getRawfileName() opens and scans the entire input mzML and can throw, despite a getter name
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/APPLICATIONS/SearchEngineBase.h · _app-topp_

```cpp
std::string getRawfileName(int ms_level = 2) const
```
- **Expectation:** A const getter named getRawfileName(ms_level=2) reads like it returns the configured -in path; the default ms_level argument looks incidental.
- **Actual:** It reads the -in path, then for mzML actually opens the file and runs MzMLFile::getCentroidInfo over MS level ms_level, throwing FileEmpty if no spectra of that level exist and IllegalArgument if spectra are profile / not centroided (unless -force). So a 'get the filename' call performs heavy file I/O and validation and can abort the tool.
- **Evidence:** SearchEngineBase.cpp:31-48: `std::string inputfile_name = getStringOption_("in"); ... MzMLFile mzml; ... const auto& centroid_info = mzml.getCentroidInfo(inputfile_name); ... throw Exception::FileEmpty(...)` and lines 58-59 throw IllegalArgument on profile data.
- **Fix:** Rename to something verb-like that signals work (e.g. resolveAndValidateRawFile_) or document prominently in one line that the call parses the file and may throw. The throwing behavior IS documented in the header @throws block, so this is borderline; the name is the residual surprise. Renaming is source-breaking (protected; few callers), so prefer keeping name + keeping the @throws doc; flag only the name/I/O mismatch.
- **Verifier correction:** getRawfileName(ms_level=2) does perform file I/O and can throw despite its accessor-style name, but it does NOT scan the entire mzML: MzMLFile::getCentroidInfo defaults to first_n_spectra_only=10 and aborts parsing after ~10 spectra (EndParsingSoftly). The throwing/validation behavior (FileEmpty on no spectra of the level; IllegalArgument on profile / non-centroided unless -force) is fully documented in the header @brief and @throws block. The residual issue is purely the accessor-like name vs. its resolve-and-validate work; severity is low (loud, documented, recoverable, correct domain behavior). Keep the name + existing docs; optionally add a one-line note that the call parses the file and may throw.

### [APPL-4] TOPPBase::getToolPrefix — getToolPrefix()/getIniLocation_() yield 'ToolName:-1:' before main() assigns the instance number
`low` · `lifecycle-precondition` · ABI: `none` · src/openms/include/OpenMS/APPLICATIONS/TOPPBase.h · _app-topp_

```cpp
std::string getToolPrefix() const
```
- **Expectation:** Header says the prefix is 'later found in the INI file ... f.e.: "FileConverter:1:"', implying a well-formed instance-qualified prefix.
- **Actual:** getToolPrefix() returns tool_name_ + ':' + instance_number_ + ':', and instance_number_ is constructed as -1 and only overwritten during main() startup. If called before the instance is parsed from the command line, it silently returns e.g. 'FileConverter:-1:' rather than failing or returning a documented sentinel — a plausible source of a wrong INI lookup location.
- **Evidence:** TOPPBase.cpp:92 `return tool_name_ + ":" + instance_number_ + ":";`; constructor TOPPBase.cpp:98 `instance_number_(-1)`; reassignment via const_cast at TOPPBase.cpp:191.
- **Fix:** Document that getToolPrefix()/getIniLocation_() are only valid after main()/command-line parsing has run, or assert instance_number_>=0. Documentation/assertion is ABI-neutral.
- **Verifier correction:** getToolPrefix() (public) depends on instance_number_ being assigned during main() startup; before that it returns "ToolName:-1:" using the sentinel -1. The string output described in the claim is correct (a custom operator+(std::string,int) in StringUtils.h does decimal conversion). But this is an undocumented lifecycle precondition on a public helper, not a silent-failure data bug: -1 is an unambiguous invalid-instance sentinel (real instances start at 1), so an early call produces a non-matching INI section rather than silently reading the wrong parameters, and no actual caller invokes it before main()/command-line parsing. Reasonable mild hardening: document that getToolPrefix()/getIniLocation_() are only valid after main() runs, or assert instance_number_ >= 0. Doc/assert only — ABI-neutral.

### [APPL-5] INIUpdater::getNewToolName — getNewToolName returns the unchanged name (and true) for tools that were never renamed
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/APPLICATIONS/INIUpdater.h · _app-topp_

```cpp
bool getNewToolName(const std::string& old_name, const std::string& tools_type, std::string& new_name)
```
- **Expectation:** A method named getNewToolName that 'Finds the name of the new tool' and returns true on success suggests true means 'a rename mapping was found' and new_name holds a different, updated name.
- **Actual:** When no rename mapping exists, it falls back to ToolHandler::getTOPPToolList() and, if the OLD name is a current tool, sets new_name = old_name and returns true. So a caller that interprets true as 'this tool was renamed' will wrongly treat unchanged tools as renamed; only the old==new equality reveals the no-op.
- **Evidence:** INIUpdater.cpp:87-93: `// default to ToolHandler ... if (topp.contains(old_name)) { new_name = old_name; return true; }`
- **Fix:** Document that true means 'a resolvable current tool name was produced (possibly identical to old_name)', not 'renamed'. Documentation-only; ABI-neutral. The header comment is also a C-style /* */ block (won't render in Doxygen) and lacks @param/@return for new_name.
- **Verifier correction:** getNewToolName returns true and new_name==old_name when no rename mapping exists but old_name is a current TOPP tool (INIUpdater.cpp:87-93). This is technically a "true means resolvable-current-name, possibly identical" semantic rather than "renamed". However the two in-tree callers explicitly handle this ("// might be the same name") and the test treats same-name resolution as expected success, so there is no actual misuse and no wrong-result path. The actionable point is documentation only: the header uses a non-Doxygen /* */ block and omits @param for the out-parameter new_name; it should state that true means "a current tool name was produced (possibly equal to old_name)", not "renamed".

### [QC-1] FWHM — FWHM class-level documentation describes m/z error calibration, not FWHM
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/QC/FWHM.h · _qc-metrics-a_

```cpp
class OPENMS_DLLAPI FWHM : public QCBase
```
- **Expectation:** The class Doxygen brief should describe what FWHM does: moving FWHM (peak width) metavalues from features to their PeptideIdentifications.
- **Actual:** Lines 17-27 read: "@brief QC metric calculating (un)calibrated m/z error" and "The metric sets m/z-values of the original experiment and the calculated reference m/z-values, uncalibrated m/z error (ppm) and calibrated m/z error (ppm) as metavalues of all PeptideIdentifications". This is copy-pasted from MzCalibration and is entirely unrelated. FWHM.cpp only copies the feature's 'FWHM'/'model_FWHM' metavalue onto each PeptideIdentification (FWHM.cpp lines 14-38). The compute() @brief on line 38 is correct and directly contradicts the class brief.
- **Evidence:** FWHM.h:18-27 '@brief QC metric calculating (un)calibrated m/z error ... sets m/z-values ... uncalibrated m/z error (ppm) and calibrated m/z error (ppm)'. FWHM.cpp:22 'pi.setMetaValue("FWHM", f.getMetaValue("FWHM"));'
- **Fix:** Fix the class-level Doxygen to describe FWHM propagation (peak full-width-half-maximum), matching compute()'s own brief. Documentation-only change. No ABI impact.
- **Verifier correction:** The class-level Doxygen @brief in FWHM.h (lines 17-27) is copy-pasted from MzCalibration.h and wrongly describes (un)calibrated m/z error calibration. The class actually propagates peak FWHM (full-width-half-maximum) by moving the feature's 'FWHM'/'model_FWHM' metavalue onto each PeptideIdentification, exactly as compute()'s own @brief (line 38) states. Fix is documentation-only: replace the class brief with an accurate description of FWHM propagation. No ABI impact.

### [QC-2] MissedCleavages::compute — Class doc promises 'FWHM' and 'mass' metavalue annotation that compute() never writes
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/QC/MissedCleavages.h · _qc-metrics-a_

```cpp
void compute(FeatureMap& fmap)
```
- **Expectation:** Per the class doc, PeptideHits should be augmented with 'missed_cleavages', 'FWHM', and 'mass' metavalues after compute().
- **Actual:** The class Doxygen (lines 22-31) states the PeptideHits are augmented with MetaInformation: 'missed_cleavages', 'FWHM' (from feature's 'FWHM' or 'model_FWHM'), and 'mass' (experimental mass of peptide). The implementation only ever calls setMetaValue("missed_cleavages", num_mc) (MissedCleavages.cpp:39). Neither 'FWHM' nor 'mass' is set anywhere in MissedCleavages.cpp. A caller relying on these two documented annotations gets nothing (and metaValueExists returns false / getMetaValue throws downstream).
- **Evidence:** MissedCleavages.h:27-30 "augmented with MetaInformation: - 'missed_cleavages' - 'FWHM' ... - 'mass' (experimental mass of peptide)". MissedCleavages.cpp:39 only 'pep_id.getHits()[0].setMetaValue("missed_cleavages", num_mc);'
- **Fix:** Either remove the 'FWHM' and 'mass' bullets from the class doc (documentation-only, no ABI impact) or, if the annotations are intended, add the setMetaValue calls. The doc fix is the safe additive change.
- **Verifier correction:** The class-level Doxygen of MissedCleavages (MissedCleavages.h:27-30) incorrectly states that compute() augments PeptideHits with 'FWHM' and 'mass' metavalues in addition to 'missed_cleavages'. The implementation (MissedCleavages.cpp:39) only writes 'missed_cleavages'. Crucially, the method-level doc on compute() (lines 46-54) is correct and promises only 'missed_cleavages', so a developer reading the function's own contract is not misled. The 'FWHM' and 'mass' annotations described in the class blurb are actually produced by separate QC classes (FWHM.cpp and PeptideMass.cpp), whose bullet text was apparently copied here. Impact is a class-summary documentation inaccuracy (low severity), not a behavioral defect; fix by removing the two bullets from the class doc — documentation-only, no ABI impact.

### [QC-5] Contaminants::getResults — getResults() is non-const, breaking the const-getter convention of every sibling QC metric
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/QC/Contaminants.h · _qc-metrics-a_

```cpp
const std::vector<Contaminants::ContaminantsSummary>& getResults();
```
- **Expectation:** A read-only results accessor named getResults() should be const, as it is on every other QC metric (FragmentMassError, MissedCleavages, DBSuitability, IdentificationSummary's siblings).
- **Actual:** Contaminants::getResults() is declared without const (header line 110, cpp line 149) even though it only returns a reference to results_ and mutates nothing. By contrast FragmentMassError::getResults() (FragmentMassError.h:94), MissedCleavages::getResults() (MissedCleavages.h:62), and DBSuitability::getResults() (DBSuitability.h:215) are all const. As a result one cannot read contaminant results from a const Contaminants&, surprising callers who treat the QC metrics uniformly.
- **Evidence:** Contaminants.h:110 'const std::vector<Contaminants::ContaminantsSummary>& getResults();' (no trailing const). Compare FragmentMassError.h:94 'const std::vector<Statistics>& getResults() const;'.
- **Fix:** Add const to getResults(). This is source-compatible for nearly all callers (const member function is callable on non-const objects) but is technically an ABI-breaking signature change for the symbol; safe to do at a minor version. Mark abi as source-compatible since callers recompile cleanly.
- **Verifier correction:** Contaminants::getResults() (header line 110, cpp line 149) is non-const while every other QC sibling's getResults() is const (FragmentMassError.h:94, MissedCleavages.h:62, PSMExplainedIonCurrent.h:88, Ms2IdentificationRate.h:114, DBSuitability.h:215, TIC.h:83). It returns a const reference to results_ and mutates nothing, so it should be const. Impact is limited to a compile error when called on a const Contaminants& — a loud, recoverable surprise, not silently wrong data — so this is low, not high/medium. Fix: add trailing const; source-compatible for all current callers.

### [QC-12] QCBase::Requires::NOTHING — Requires::NOTHING is not the empty/zero flag: Status(Requires::NOTHING) sets bit 0 and is non-empty
`low` · `surprising-default` · ABI: `source-compatible` · src/openms/include/OpenMS/QC/QCBase.h · _qc-metrics-b_

```cpp
enum class Requires : UInt64 { NOTHING, RAWMZML, ... }; using Status = FlagSet<Requires>
```
- **Expectation:** An enum member literally named NOTHING used with a FlagSet ('does not require anything', per the comment) should map to the empty bitset, so Status(Requires::NOTHING).empty() == true and it acts as identity under OR.
- **Actual:** FlagSet converts an enum value r to bit 2^r (getPow_ returns 'UInt64(1) << r'). NOTHING==0 therefore becomes 2^0 == 1 — the SAME bit RAWMZML would occupy if RAWMZML were value 0, and a distinct non-zero bit in general. So Status(Requires::NOTHING) has value()==1 and empty()==false; ORing NOTHING into a Status spuriously sets bit 0. The name promises 'no requirement' but the value is a real flag.
- **Evidence:** FlagSet.h:203-207 getPow_ returns 'UInt64(1) << UInt64(en)'; FlagSet.h:36-39 'explicit FlagSet(const ENUM& en) : value_(getPow_(en))'; QCBase.h:35 'NOTHING,  //< default, does not require anything'.
- **Fix:** Document clearly that NOTHING must never be used to build a Status (use a default-constructed Status() for 'no requirements'), or reserve enum value 0 as an unused sentinel and start real flags at 1. Doc/constant change is source-compatible; reordering the enum is breaking.
- **Verifier correction:** Requires::NOTHING (enum value 0) does not map to the empty FlagSet: Status(Requires::NOTHING) has value()==1 and empty()==false, because FlagSet maps enum value r to bit 1<<r (so value 0 -> bit 1). This is the documented and unit-tested FlagSet contract (the zeroth enum value is bit 1; the empty set is the default-constructed FlagSet()). The name NOTHING/"does not require anything" therefore mismatches its value. Severity is low: NOTHING is never actually used to construct or OR into a Status anywhere in the codebase (all requirements() use the empty default ctor QCBase::Status() for 'no requirements'), so no current behavior is wrong; this is a latent footgun, and misuse would fail isRunnable() loudly rather than silently corrupt results. The claim's "same bit as RAWMZML" wording is incorrect: NOTHING->bit 1, RAWMZML->bit 2, distinct bits; bit 0/value 1 is just unassigned. Recommended safe fix is a doc/comment clarification (source-compatible); renumbering the enum to reserve 0 as a sentinel would be ABI-breaking.

### [QC-14] PSMExplainedIonCurrent::compute — compute() silently pushes a zero-score Statistics{} when the FeatureMap has no PeptideIDs, indistinguishable from a real result
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/QC/PSMExplainedIonCurrent.h · _qc-metrics-b_

```cpp
void compute(FeatureMap&, const MSExperiment&, const QCBase::SpectraMap&, ToleranceUnit, double)
```
- **Expectation:** Given the documented throw behavior ('@throws MissingInformation If PSMExplainedIonCurrent couldnt be calculated for any spectrum'), a caller expects compute() to signal failure when nothing can be computed, and getResults() entries to be valid computed statistics.
- **Actual:** When the map has no peptide IDs, compute() does 'results_.push_back(result);' with a default-constructed Statistics (average_correctness=0, variance_correctness=0) and returns without warning or exception. getResults().back() then reports a fully valid-looking score of 0.0 that is indistinguishable from a genuinely computed 0. Callers like DatabaseSuitability.cpp:230 ('eic.getResults()[0]') consume this directly.
- **Evidence:** PSMExplainedIonCurrent.cpp:134-140 'bool has_pepIDs = QCBase::hasPepID(fmap); if (!has_pepIDs) { results_.push_back(result); return; }' with result default-constructed (average_correctness=0).
- **Fix:** Either skip pushing a result (so getResults() size reflects only computed metrics) or throw/flag the no-ID case, consistent with the all-spectra-failed throw later in the same function. Behavior change; gate behind doc update.
- **Verifier correction:** The behavior exists as quoted, but (a) it is an intentional, module-wide QC convention (same pattern in MissedCleavages and FragmentMassError, gated by QCBase::hasPepID): no peptide IDs => metric not applicable => push a placeholder Statistics{} to keep the results vector index-aligned; IDs present but uncomputable => throw MissingInformation. (b) The cited consumer DatabaseSuitability.cpp:229 uses the PeptideIdentificationList overload (guarded on pep_ids.empty()), not the FeatureMap overload in the claim's signature, and is fed FDR-processed IDs, so it does not consume the silent zero in normal operation. The genuine residual issue is narrow: the doc comment documents only the throw paths and omits the silent no-ID placeholder, and the two 'nothing computed' cases are handled inconsistently. This is a low-severity documentation/consistency wart, not a silently-wrong-results bug for reasonable use.

### [QC-15] RTAlignment::compute — compute() writes six undocumented feature-level meta values (rt_align_start/end, rt_raw_start/end, rt_align, rt_raw) beyond the documented PepID-only rt_raw/rt_align
`low` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/QC/RTAlignment.h · _qc-metrics-b_

```cpp
void compute(FeatureMap& fm, const TransformationDescription& trafo) const
```
- **Expectation:** The class doc states 'Sets meta values rt_raw and rt_align in PeptideIdentifications of the featureMaps PepIDs. It does not change the RT of the features.' and the method @brief promises 'sets meta values rt_raw and rt_align in all PepIDs'. A caller expects only PeptideIdentifications to gain only rt_raw and rt_align.
- **Actual:** For every Feature, compute() also sets feature-level meta values: feature.setMetaValue("rt_align", ...), "rt_raw", "rt_align_start", "rt_align_end", "rt_raw_start", "rt_raw_end". These feature meta values (and the *_start/*_end pairs) are entirely undocumented in the header; a caller relying on the doc would not expect features to be annotated at all.
- **Evidence:** RTAlignment.cpp:49-54 'feature.setMetaValue("rt_align", ...); feature.setMetaValue("rt_raw", ...); feature.setMetaValue("rt_align_start", ...); feature.setMetaValue("rt_align_end", ...); feature.setMetaValue("rt_raw_start", ...); feature.setMetaValue("rt_raw_end", ...);' vs header lines 20-26/37-44.
- **Fix:** Update the class/method documentation to list the feature-level meta values written (rt_align, rt_raw and the *_start/*_end bounding-box variants). Source-compatible doc fix; no signature change.
- **Verifier correction:** compute(FeatureMap&, const TransformationDescription&) writes six undocumented Feature-level meta values (rt_align, rt_raw, rt_align_start, rt_align_end, rt_raw_start, rt_raw_end) on every Feature, beyond the header's documented PepID-only rt_raw/rt_align. The four *_start/*_end bounding-box variants are wholly undocumented; rt_align/rt_raw are documented only for PeptideIdentifications. Severity is low: the side effect is additive, visible in featureXML output, and does not produce silently wrong results, data loss, or crashes. Fix is documentation-only (list the feature-level meta values), source-compatible, no signature change.

### [QC-17] Ms2IdentificationRate::IdentificationRateData::identification_rate — Field named 'identification_rate' / metric 'MS2 Identification Rate' is a 0..1 fraction, but mzTab export silently multiplies by 100 (percent)
`low` · `unit-or-index` · ABI: `none` · src/openms/include/OpenMS/QC/Ms2IdentificationRate.h · _qc-metrics-b_

```cpp
double identification_rate = 0.
```
- **Expectation:** Consumers of a field called 'identification_rate' need to know whether it is a fraction (0..1) or a percentage (0..100). The header doc describes the metric as '#identified PSMs divided by total number of MS2 scans' (a fraction).
- **Actual:** getResults()[i].identification_rate is the raw ratio in [0,1] (writeResults_: 'ratio = (double)pep_ids_count / ms2_spectra_count'). But addMetaDataMetricsToMzTab() emits '100 * ms2_irs[i].identification_rate', i.e. a percentage, under the same conceptual name. A caller mixing the struct value with the mzTab value (or expecting a percentage from a field literally called a 'rate') gets a 100x discrepancy.
- **Evidence:** Ms2IdentificationRate.cpp:137 'double ratio = (double)pep_ids_count / ms2_spectra_count;' stored into identification_rate; Ms2IdentificationRate.cpp:83 'ms2_ir.setValue(StringUtils::toStr(100 * ms2_irs[i].identification_rate));'.
- **Fix:** Document the field's unit (fraction 0..1) on IdentificationRateData and note the mzTab serialization uses percent; or unify the units. Source-compatible doc fix.
- **Verifier correction:** The struct field identification_rate is consistently a fraction [0,1] across all in-tree consumers (unit test asserts 1/3; DatabaseSuitability.cpp annotates "[0, 1]"). The only percent value is produced by addMetaDataMetricsToMzTab() (line 83) into a distinct mzTab custom parameter named "MS2_ID_Rate_N" — a different artifact and a different name from the struct field. No code path co-mingles the two, so the claimed "100x discrepancy for a caller mixing the values" does not occur in practice. The legitimate, low-severity improvement is to add an explicit unit doc-comment to IdentificationRateData::identification_rate (fraction 0..1) and note that the mzTab serialization presents it as a percentage. Source-compatible doc-only fix.

### [QC-18] MQExporterHelper::isValid — isValid() is a read-only-sounding predicate that creates and deletes a file on disk as a side effect
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/QC/MQExporterHelper.h · _qc-mq-export_

```cpp
static bool isValid(const std::string& filename_)
```
- **Expectation:** A static bool isValid(filename) named like a pure predicate (and documented as merely checking 'the path in the ctor was not empty and could be created') should inspect state without modifying the filesystem. Sibling isValid() methods across OpenMS (e.g. MzMLFile().isValid(in,os), FeatureXMLFile().isValid) are pure validators.
- **Actual:** isValid() forwards to File::writable(filename). In File::writable (src/openms/source/SYSTEM/File.cpp:427-452), when the file does NOT yet exist it probes writability by opening an std::ofstream on the path (creating the file) and then calling std::remove() to delete it. So calling isValid() touches the filesystem: it creates and removes the target file, and can clobber/truncate an existing-but-empty target via the ofstream open if timing differs.
- **Evidence:** File::writable: 'std::ofstream f(file.c_str()); bool ok = f.is_open() && f.good(); f.close(); if (ok) { std::remove(file.c_str()); }' ; MQExporterHelper::isValid: 'return File::writable(filename);'
- **Fix:** Rename to reflect the probe (e.g. isWritable/canWrite) or document the create+delete side effect explicitly in the header. Additive/safe fix: update the doc comment to state it probes by creating and removing the file. abi_impact none for a doc fix; a rename would be source-compatible if the old name is kept as a deprecated alias.
- **Verifier correction:** isValid(filename) is a writability probe (it forwards to File::writable), which for a non-existent path creates the file via std::ofstream and then std::remove()s it — a filesystem side effect under a predicate-style name. Two corrections to the claim: (a) the side effect is partially documented in the header ('Checks if file is writable ... could be created') and in the caller-class docs, so it is not undocumented; (b) it does NOT clobber/truncate an existing target — File::writable opens the ofstream only when fs::exists()==false, so existing files go through access(W_OK) and are never opened/truncated (a TOCTOU race aside). The probe creates/removes only the intended-output file, which the callers are about to write anyway. Recommendation: rename to isWritable/canWrite (keeping isValid as a deprecated alias, source-compatible) or sharpen the doc to state it probes by creating and removing the file. Severity is low (mild surprise, no data loss), not high.

### [QC-19] MQExporterHelper::isValid — isValid doc claims it checks 'the path in the ctor', but the parameter is a filename and call sites pass a directory, so it validates inconsistent things
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/QC/MQExporterHelper.h · _qc-mq-export_

```cpp
static bool isValid(const std::string& filename_)
```
- **Expectation:** The doc says isValid returns whether 'the path in the ctor was not empty and could be created', implying it reflects the constructed object's state.
- **Actual:** isValid takes an arbitrary string named filename_ and just calls File::writable on it. The constructor stores filename_ = path + "/evidence.txt" (a file), and exportFeatureMap calls isValid(filename_) (a file). But the public TOPP caller QualityControl.cpp:417 calls MQExporterHelper::isValid(out_txt_dir) passing a DIRECTORY. File::writable behaves differently for a directory (always exists -> access(W_OK)) vs a not-yet-existing file (create+probe+delete). The same-named call therefore means different things depending on what the caller passes, and it does not actually consult the ctor's stored path at all.
- **Evidence:** Header doc: 'Checks if file is writable (i.e. the path in the ctor was not empty and could be created)'. QualityControl.cpp:417: 'if (MQExporterHelper::isValid(out_txt_dir))'. Exporter ctor: 'filename_ = path + "/evidence.txt";'. exportFeatureMap: 'if (!MQExporterHelper::isValid(filename_))'.
- **Fix:** Clarify the parameter contract (file path, not directory) in the doc and rename the parameter from filename_ (which uses the OpenMS private-member underscore convention on a public static arg). abi_impact none.
- **Verifier correction:** isValid is a static helper that returns File::writable(arg) on any string passed; the header doc wrongly describes it as checking "the path in the ctor", which it never consults. Internal callers pass a file (dir/evidence.txt or dir/msms.txt) while the TOPP caller QualityControl.cpp:417 passes a directory, so File::writable takes different code paths. Real but low severity: the directory check is a benign proxy and any real failure is loud (FileNotWritable). Note the `filename_` underscore is only in the header declaration (line 157); the definition uses `filename`. Fix = correct the doc to state it validates the given file/dir path and is independent of ctor state, and rename the header param to drop the private-member underscore.

### [QC-21] MQEvidence::exportFeatureMap — exportFeatureMap throws mid-export after partially writing rows when a feature lacks a ConsensusFeature
`low` · `partial-output-on-documented-throw` · ABI: `none` · src/openms/include/OpenMS/QC/MQEvidenceExporter.h · _qc-mq-export_

```cpp
void exportFeatureMap(const FeatureMap&, const ConsensusMap&, const MSExperiment& = {}, const std::map<std::string,std::string>& = {})
```
- **Expectation:** An 'export' that the doc says writes 'one row per feature' would be expected to either validate inputs up front or skip/flag unmatched features, not abort partway leaving a truncated output file.
- **Actual:** exportFeatureMap iterates features and, the first time a feature's unique id is not found in the ConsensusMap mapping, throws Exception::MissingInformation. Because rows are streamed to file_ inside the loop, the on-disk evidence.txt/msms.txt already contains the header plus all rows written before the offending feature, and file_ is not flushed/cleaned on the throw path. The caller is left with a partially written, header-bearing file that looks valid.
- **Evidence:** exportFeatureMap loop: 'if (c_id != fTc.end()) { exportRowFromFeature_(...); } else { throw Exception::MissingInformation(..., "Feature in FeatureMap has no associated ConsensusFeature."); }' with the throw occurring after prior iterations have already written to file_.
- **Fix:** Document that partial output may remain on throw, and/or validate the full feature-to-CF mapping before writing any rows. Behavioral hardening is source/ABI compatible.
- **Verifier correction:** The throw is real and occurs mid-loop after partial rows have been streamed, and the partially-written header-bearing file is indeed left on disk (flushed by ~MQEvidence's file_.close()). But the throw is explicitly documented in the header (lines 97-98), so this is not a 'surprising/undocumented throw.' The actual (mild) surprise is that the export is non-transactional: on the documented throw it leaves a truncated, valid-looking evidence.txt with no flag/cleanup. The throw is loud and only triggers on input violating the documented precondition (a feature with no corresponding ConsensusFeature), so the practical hazard is limited to a caller that both swallows the exception and trusts the stale file. Recommendation (document partial-output-on-throw and/or validate the full feature-to-CF mapping before writing any rows) is source- and ABI-compatible.

### [QC-22] MQExporterHelper::hasValidPepID_, hasPeptideIdentifications_, makeFeatureUIDtoConsensusMapIndex_, proteinGroupID_ — Public API uses the OpenMS trailing-underscore private-member naming convention, signaling 'internal' for callable API
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/QC/MQExporterHelper.h · _qc-mq-export_

```cpp
static bool hasValidPepID_(...); static bool hasPeptideIdentifications_(...); static std::map<...> makeFeatureUIDtoConsensusMapIndex_(...); static Size proteinGroupID_(...)
```
- **Expectation:** In OpenMS the trailing-underscore suffix denotes private/internal members. A developer scanning the public section of a header expects underscore-suffixed names to be non-API.
- **Actual:** These four functions are declared under 'public:' yet carry the trailing-underscore convention, contradicting the visibility. This makes them look like internal helpers a caller should not use, while they are in fact the supported helper API (and proteinGroupID_/makeFeatureUIDtoConsensusMapIndex_ have public-facing doc).
- **Evidence:** All declared after 'public:' in MQExporterHelper.h with trailing underscores, e.g. 'static bool hasValidPepID_(...)', 'static OpenMS::Size proteinGroupID_(...)', 'static std::map<OpenMS::Size, OpenMS::Size> makeFeatureUIDtoConsensusMapIndex_(...)'.
- **Fix:** Drop the trailing underscore from the public names (provide deprecated aliases for source compatibility) or move genuinely-internal ones to a detail/private scope. abi_impact source-compatible with aliases.
- **Verifier correction:** The four functions do violate the documented OpenMS trailing-underscore convention (which reserves '_' suffixes for protected/private members) while being declared public, so the surprise is real. However, they are NOT a "supported public helper API": they are used only by two sibling classes in the same QC module (MQMsmsExporter.cpp, MQEvidenceExporter.cpp), are absent from pyOpenMS, and have no standalone test. The accurate framing is the inverse of the title — these are genuinely-internal helpers that were forced into 'public' merely to allow cross-class access, so the appropriate fix is to relocate them to a detail/private scope rather than to canonize them as renamed public API. Severity is low (cosmetic/convention inconsistency, no functional impact). ABI: removing/renaming symbols on this OPENMS_DLLAPI-exported class is technically breaking unless deprecated aliases are kept; with aliases (or because they have no real external consumers) the change is source-compatible.

### [IONM-1] IMDataConverter::reshapeIMFrameToMany — Documented MissingInformation throw is skipped for an empty input spectrum
`low` · `documentation-contract-mismatch` · ABI: `none` · src/openms/include/OpenMS/IONMOBILITY/IMDataConverter.h · _small-im-interfaces-imaging_

```cpp
static MSExperiment reshapeIMFrameToMany(MSSpectrum im_frame)
```
- **Expectation:** Per the header: '@throws Exception::MissingInformation if im_frame does not have IM data in floatDataArrays.' An input lacking IM data should throw.
- **Actual:** If the input spectrum is empty, the function returns an empty MSExperiment WITHOUT ever checking for or requiring an IM float-data array, so the documented exception is silently not thrown. An empty frame (which also has no IM array) is treated as success.
- **Evidence:** IMDataConverter.cpp lines 98-101: `if (im_frame.empty()) { // nothing to split (we do not even check for IM data, for robustness) return out; }` — the getIMData() call that throws is only reached for non-empty input (line 110).
- **Fix:** Either document the empty-input early-return as an explicit exception to the throw contract, or validate IM data presence before the empty short-circuit. Behavior-only change for the empty edge case; pick the documentation fix to stay ABI-safe.
- **Verifier correction:** The empty-input early return (IMDataConverter.cpp:98-101) bypasses the documented MissingInformation throw because the IM-data check (via getIMData() at line 110, which throws per MSSpectrum.cpp:794) is only reached for non-empty input. The mismatch is real but benign: an empty frame correctly yields an empty MSExperiment, so this is a documentation/contract gap rather than a silent-failure producing wrong results or data loss. Recommended fix: document the empty-input early-return as an explicit exception to the throw contract (ABI-safe, no signature change).

### [IONM-2] IMTypes::oneOverK0ToCCS / IMTypes::ccsToOneOverK0 — buffer_gas_mass is not validated; non-positive values yield NaN/garbage silently
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/IONMOBILITY/IMTypes.h · _small-im-interfaces-imaging_

```cpp
static double oneOverK0ToCCS(double one_over_k0, double mz, int charge, double buffer_gas_mass = N2_BUFFER_GAS_MASS)
```
- **Expectation:** Given the functions throw Exception::InvalidValue for every other out-of-domain argument (one_over_k0<=0, mz<=0, charge==0), a non-physical buffer_gas_mass (0 or negative) would likewise be rejected.
- **Actual:** buffer_gas_mass is used unchecked in reducedMass_; a value of 0 makes the reduced mass 0 and the result +Inf, and a negative value makes sqrt(mu) NaN, returning NaN without any exception. The asymmetry vs the other three validated arguments is surprising.
- **Evidence:** IMTypes.cpp lines 186-193 validate only one_over_k0/mz/charge; reducedMass_ (lines 177-181) divides by (ion_mass + buffer_gas_mass) and the caller takes std::sqrt(mu) with no buffer_gas_mass>0 guard. Same pattern in ccsToOneOverK0 (lines 198-205).
- **Fix:** Add `buffer_gas_mass <= 0.0` to the existing InvalidValue guard in both functions. Behavior change only for already-invalid input; ABI-safe.
- **Verifier correction:** Real but low severity. buffer_gas_mass is an explicit configuration argument with a safe default (N2_BUFFER_GAS_MASS=28.0); the common default path never triggers the bug. The silent +Inf (mass==0) / NaN (negative mass) only occurs when a caller deliberately passes a non-physical drift-gas mass, which is a programmer error rather than reasonable-use bad data — unlike the validated mz/charge/one_over_k0 which are measured per-ion inputs. Recommended fix is correct and ABI-safe (function-body comparison only): add `|| buffer_gas_mass <= 0.0` to the existing InvalidValue guards in both oneOverK0ToCCS and ccsToOneOverK0.

### [IONM-3] MSImagingGeometry::addPixel — Bounds checking is fully disabled whenever width OR height is 0
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/IMAGING/MSImagingGeometry.h · _small-im-interfaces-imaging_

```cpp
void addPixel(UInt x, UInt y, Size spectrum_index)
```
- **Expectation:** Once setDimensions has been called, addPixel rejects out-of-range coordinates (as the header states).
- **Actual:** The bounds check requires BOTH width_>0 AND height_>0. If the caller sets a 1-D-like geometry (e.g. width=5, height=0, or never calls setDimensions), out-of-range pixels at any coordinate are accepted silently. A single zero dimension turns off all coordinate validation.
- **Evidence:** MSImagingGeometry.cpp lines 56-61: `if (width_ > 0 && height_ > 0 && (x >= width_ || y >= height_)) { throw ... }`.
- **Fix:** Treat any set dimension independently: validate x against width_ when width_>0 and y against height_ when height_>0, instead of gating on both being non-zero. ABI-safe behavior tightening (only rejects coordinates that were already conceptually out of range).
- **Verifier correction:** The bounds check is genuinely disabled when EITHER width_ or height_ is 0 (the `&&` gate). This is a real but minor doc/implementation contract mismatch, not a practical data-integrity hazard. The claim's framing of "width=5,height=0" as a legitimate 1-D geometry is incorrect: a zero dimension denotes an empty grid (no pixel is ever valid), and real 1D line scans are modeled as height=1/width=1, where the check is fully active. The only production caller (BrukerTimsImagingFile) never reaches the gap with actual pixels. Recommended fix — validate each dimension independently (reject x when width_>0 && x>=width_; reject y when height_>0 && y>=height_) — is ABI-safe behavior tightening that only rejects already-out-of-range coordinates. Severity is low (degenerate-input contract gap), not a high-impact silent corruption.

### [IONM-5] Interfaces::SpectrumMeta — Default-constructed SpectrumMeta leaves RT and ms_level uninitialized
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/INTERFACES/DataStructures.h · _small-im-interfaces-imaging_

```cpp
SpectrumMeta() : index(0) {}
```
- **Expectation:** Default-constructing a small POD-like meta struct yields well-defined zero/default member values.
- **Actual:** The constructor initializes only `index`; `double RT` and `int ms_level` are left with indeterminate values. A consumer reading meta.RT or meta.ms_level right after default construction reads garbage, which is easy to do given the struct otherwise looks value-initialized.
- **Evidence:** DataStructures.h lines 147-150: `SpectrumMeta() : index(0) {}` with RT (line 142) and ms_level (line 145) not in the init list and no in-class initializers.
- **Fix:** Add in-class member initializers (`double RT = 0.0; int ms_level = 0;`) or extend the ctor init list. Source-compatible; changes only previously-indeterminate values.
- **Verifier correction:** Interfaces::SpectrumMeta's default ctor `SpectrumMeta() : index(0) {}` (DataStructures.h:147-150) initializes only `index`, leaving `double RT` (line 142) and `int ms_level` (line 145) indeterminate. The claim is accurate, but the practical risk is low because this INTERFACES struct is constructed only by the test mock (MockImplementation.cpp:46), which discards the fields; production code uses the distinct OpenSwath::SpectrumMeta struct instead. It remains a latent trap for future implementers of the ISpectrumAccess interface. Fix with in-class initializers (`double RT = 0.0; int ms_level = 0;`); this is source- and ABI-compatible.

### [IONM-6] Interfaces::ChromatogramMeta — Default-constructed ChromatogramMeta leaves precursor/product isolation targets uninitialized
`low` · `uninitialized-member` · ABI: `none` · src/openms/include/OpenMS/INTERFACES/DataStructures.h · _small-im-interfaces-imaging_

```cpp
ChromatogramMeta() : index() {}
```
- **Expectation:** Default construction yields defined member values for the m/z target fields.
- **Actual:** Only `index` is value-initialized; `double precursor_isolation_target` and `double product_isolation_target` are left indeterminate. Code that constructs a ChromatogramMeta and reads these targets before assignment reads garbage m/z values.
- **Evidence:** DataStructures.h lines 67-70: `ChromatogramMeta() : index() {}` with precursor_isolation_target/product_isolation_target (lines 64-65) absent from the init list and lacking in-class initializers.
- **Fix:** Add in-class initializers (`double precursor_isolation_target = 0.0; double product_isolation_target = 0.0;`). Source-compatible.
- **Verifier correction:** Module is INTERFACES, not IONMOBILITY; category is uninitialized-member, not return-value. The defect is real (precursor_isolation_target/product_isolation_target left indeterminate by the default constructor), but the two fields are never read or written anywhere in OpenMS, so the "reads garbage m/z values" impact is hypothetical (latent risk for external API consumers of this public interface header). Severity is low. Adding in-class initializers (double precursor_isolation_target = 0.0; double product_isolation_target = 0.0;) is the correct fix and is source-compatible with no ABI break.

### [SYST-1] SysInfo::getProcessMemoryConsumption / SysInfo::getProcessPeakMemoryConsumption — Windows doc comments for current-vs-peak memory are swapped relative to the implementation
`low` · `incorrect-documentation` · ABI: `none` · src/openms/include/OpenMS/SYSTEM/SysInfo.h · _sys-core_

```cpp
static bool getProcessMemoryConsumption(size_t& mem_virtual); static bool getProcessPeakMemoryConsumption(size_t& mem_virtual);
```
- **Expectation:** getProcessMemoryConsumption() returns the CURRENT working set and getProcessPeakMemoryConsumption() returns the PEAK working set, and the doc comments describe each accordingly.
- **Actual:** The header doc comments are inverted vs. the code. getProcessMemoryConsumption() is documented as 'On Windows, this is equivalent to "Peak Working Set (Memory)"' (SysInfo.h:30) but the implementation reads pmc.WorkingSetSize (current) (SysInfo.cpp:167). getProcessPeakMemoryConsumption() is documented as 'On Windows, this is equivalent to "Working Set (Memory)"' (SysInfo.h:38) but reads pmc.PeakWorkingSetSize (SysInfo.cpp:239). The 'current' getter is labelled 'Peak' and the 'peak' getter is labelled 'Working Set'.
- **Evidence:** SysInfo.h:30 '/// On Windows, this is equivalent to 'Peak Working Set (Memory)' in Task Manager.' on getProcessMemoryConsumption; SysInfo.cpp:167 'mem_virtual = pmc.WorkingSetSize / 1024;'. SysInfo.h:38 '/// On Windows, this is equivalent to 'Working Set (Memory)' in Task Manager.' on getProcessPeakMemoryConsumption; SysInfo.cpp:239 'mem_virtual = pmc.PeakWorkingSetSize / 1024;'.
- **Fix:** Swap the two doc comments so the names and Task-Manager descriptions agree with the code (Working Set vs Peak Working Set). Pure comment fix; no signature change.
- **Verifier correction:** The function names and first-line summaries are correct and match the implementation; only the secondary Windows Task-Manager equivalence comments are swapped. getProcessMemoryConsumption() (reads WorkingSetSize, current) should be labelled 'Working Set (Memory)', and getProcessPeakMemoryConsumption() (reads PeakWorkingSetSize, peak) should be labelled 'Peak Working Set (Memory)'. This is a documentation-comment inversion, not a misleading symbol name. Fix by swapping the two Task-Manager equivalence lines (SysInfo.h:30 and SysInfo.h:38).

### [SYST-4] SysInfo::MemUsage::delta / SysInfo::MemUsage::usage — delta()/usage() use mem_after==0 as 'not recorded yet' sentinel and silently re-sample
`low` · `hidden-side-effect` · ABI: `none` · src/openms/source/SYSTEM/SysInfo.cpp · _sys-core_

```cpp
std::string delta(const std::string& event = "delta"); std::string usage();
```
- **Expectation:** Per the header, 'When delta() or usage() are called, and the second timepoint is not recorded yet, this will be done internally' — i.e. after() is called only if the user did not already call it.
- **Actual:** The 'was after() already called?' test is 'if (mem_after == 0)'. If a legitimate after() measurement returned 0 (e.g. getProcessMemoryConsumption failed and left mem_after=0, which the API explicitly allows per SysInfo.h:34), delta()/usage() will call after() AGAIN on every invocation, re-sampling at a later time and silently producing a different/non-reproducible result than the timepoint the user intended.
- **Evidence:** SysInfo.cpp:296 'if (mem_after == 0) { after(); ... }' in delta() and SysInfo.cpp:311 'if (mem_after == 0) { after(); ... }' in usage(); contrast SysInfo.h:34 '@return True on success, false otherwise. If false is returned, then mem_virtual is set to 0.'
- **Fix:** Track a separate bool 'after_recorded_' flag instead of overloading mem_after==0 as a sentinel. Adds a private member (recompile of TU); for ABI, can be done additively in a minor release.
- **Verifier correction:** The mem_after==0 sentinel genuinely cannot distinguish "after() never called" from "after() called but the measurement failed and returned 0" (allowed per SysInfo.h:34). But the practical effect is NOT a silently-different later-timepoint resample: such failures are platform-capability (not transient), so the auto-retry almost always also returns 0, yielding a deterministic 0-based delta, and the only realistic trigger is a platform where memory querying is already broken (before() also 0). MemUsage is a diagnostic/logging-only helper, so impact is a wrong line in a human-readable memory report, not wrong data/results — hence low, not high. Note: the recommended fix (adding a bool member to the public OPENMS_DLLAPI MemUsage struct) would be ABI-BREAKING, not "additive in a minor release" as the recommendation claims; the current code as-described has no ABI impact.

### [SYST-5] File::empty — empty() returns true for a non-existent file AND for a directory, conflating 'absent', 'is-a-dir' and 'zero bytes'
`low` · `documentation` · ABI: `none` · src/openms/include/OpenMS/SYSTEM/File.h · _sys-core_

```cpp
static bool empty(const std::string& file);
```
- **Expectation:** A predicate named empty(file) reports whether an existing file has zero content; a caller commonly guards with exists() separately.
- **Actual:** empty() returns true if the path does not exist (File.cpp:146) and also returns true when file_size fails for any reason 'e.g. if it is a directory' (File.cpp:148-149). So a missing path, a directory, and a permission error all report 'empty=true', identical to a genuinely zero-byte file. A caller using empty() to decide 'safe to skip / nothing to read' cannot distinguish these and may silently treat a directory or a missing input as an empty file.
- **Evidence:** File.cpp:146 'if (!fs::exists(p, ec)) return true;' and File.cpp:147-149 'auto sz = fs::file_size(p, ec); if (ec) return true; // e.g. if it is a directory'. Header File.h:75 only states '/// Return true if the file does not exist or the file is empty'.
- **Fix:** Documented for the missing-file case but not for directories/errors. Either extend the doc to explicitly state 'returns true for directories and on any stat error' or split into exists()/size-based checks at call sites. Doc-only change is ABI-safe.
- **Verifier correction:** empty() does conflate absent/directory/permission-error/zero-byte (all return true), and the directory/stat-error case is undocumented in the public header. But this is a documentation/mild-surprise issue, not a silent-failure data hazard: the conflation collapses to the safe "nothing to read" result for all in-tree guard-style call sites, and OpenMS's own input validator (TOPPBase::inputFileReadable_) already pre-checks exists()/readable() and guards with !isDirectory(), showing the intended usage and that the directory case is handled at the framework level. Recommended fix is the documentation extension only ("returns true for directories and on any stat error"), which is ABI-safe.

### [SYST-7] StopWatch::operator< / <= / > / >= (vs operator==) — Relational operators compare only CPU time while == compares all three times plus running state, so ordering is inconsistent with equality
`low` · `misleading-documentation` · ABI: `none` · src/openms/include/OpenMS/SYSTEM/StopWatch.h · _sys-core_

```cpp
bool operator<(const StopWatch&) const; bool operator==(const StopWatch&) const;
```
- **Expectation:** For a comparable type, the relational and equality operators agree: '!(a<b) && !(b<a)' should mean a==b, and the doc says '<' is true only if the watch is lesser in ALL timings (clock, user and system).
- **Actual:** operator< is implemented as 'getCPUTime() < rhs.getCPUTime()' (StopWatch.cpp:263-266) — a single scalar (user+system) — directly contradicting the header which states '<' returns true only if 'in all timings lesser than ... (clock, user and system time)'. Meanwhile operator== compares accumulated_times_ (all components) AND last_start_ AND is_running_ (StopWatch.cpp:197-202). Two watches with equal CPU time but different wall/system split are neither < nor >, yet are also not ==, so the ordering is not a strict weak ordering consistent with ==.
- **Evidence:** StopWatch.cpp:265 'return getCPUTime() < stop_watch.getCPUTime();' vs StopWatch.h:140-145 'Return true, if the stop watch is in all timings lesser than ... (clock, user and system time).'; StopWatch.cpp:199-201 compares accumulated_times_, last_start_, is_running_.
- **Fix:** Fix the header to state that ordering is by CPU time only (doc-only, ABI-safe), or make the operators consistent with the documented all-timings semantics. Do not leave the documented contract contradicting the code.
- **Verifier correction:** The operators are NOT misnamed; the defect is that StopWatch.h:140-170 documents operator< / <= / >= / > as comparing ALL three timings (clock, user, system), but StopWatch.cpp:263-282 orders solely by getCPUTime() (= user + system, wall/clock time ignored). This makes the documented "all timings" contract false and makes the ordering inconsistent with operator== (which compares the full accumulated_times_, last_start_, and is_running_). Recommended fix is doc-only/ABI-safe: state that ordering is by CPU time only. Impact is low because this is a diagnostic timing utility whose relational operators do not drive scientific results.

### [SYST-10] RWrapper::findR — findR ignores the user-supplied `executable` name in all of its diagnostics, always reporting 'Rscript'
`low` · `misleading-diagnostics` · ABI: `none` · src/openms/include/OpenMS/SYSTEM/RWrapper.h · _sys-net-lang_

```cpp
static bool findR(const std::string& executable = "Rscript", bool verbose = true)
```
- **Expectation:** A caller passing a custom interpreter path (the parameter is documented as 'Name of the R interpreter') expects the success/failure diagnostics to refer to the interpreter they actually requested.
- **Actual:** findR runs `bp::search_path(executable)` but every status/error line hardcodes the literal string 'Rscript' regardless of the `executable` argument, so a failure for a custom interpreter prints misleading messages about 'Rscript'.
- **Evidence:** RWrapper.cpp:134 `OPENMS_LOG_INFO << "Finding R interpreter 'Rscript' ...";` ; :178 `"Error: 'Rscript' executable returned with error (command: 'Rscript " << args_str << "')"` ; :181 `"Make sure 'Rscript' is installed properly."` ; :203 `"Trying to invoke 'Rscript' ... success"` — none use the `executable` parameter.
- **Fix:** ABI-safe: substitute the `executable` parameter into the log/error strings instead of the literal 'Rscript'. Pure implementation fix, no signature change.
- **Verifier correction:** findR's actual process invocation correctly uses the user-supplied `executable` (bp::search_path at line 144), and its catch(process_error) branch at line 191 does interpolate `executable`. However four diagnostic lines (134, 178 twice, 181, 203) hardcode the literal 'Rscript' regardless of the argument, so a custom interpreter (e.g. InternalCalibration's user-facing 'rscript_executable' path param) yields inconsistent/misleading log+error text. Impact is cosmetic (verbose diagnostics only; the correct binary is still searched/run). ABI-safe fix: substitute `executable` into those literal strings.
