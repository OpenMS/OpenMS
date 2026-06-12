# OpenMS POLS Audit — Batch 4b: ANALYSIS/OPENSWATH + TARGETED

**Confirmed findings:** 82 · 7 high · 39 medium · 36 low. (assembled from solo runs after parallel-induced rate limiting was resolved)

### [ANSW-26] SpectrumAccessSqMass::SpectrumAccessSqMass(handler, indices) — Passing an empty index vector silently means 'all spectra', not 'no spectra'
`high` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h · _osw-dataaccess_

```cpp
SpectrumAccessSqMass(const OpenMS::Internal::MzMLSqliteHandler& handler, const std::vector<int>& indices)
```
- **Expectation:** Constructing a subset view with an explicit list of indices and passing an empty list yields an empty view (zero visible spectra) — the natural reading of 'expose only the spectra at the given indices'.
- **Actual:** An empty `indices` vector is treated as the sentinel for 'all spectra'. getNrSpectra() then returns the full file count, not 0. A caller that computed an empty index set (e.g. a SWATH window that matched nothing) silently gets the entire file instead of an empty selection.
- **Evidence:** Header line 93 documents `an empty vector falls back to "all spectra" semantics`; getNrSpectra() (SpectrumAccessSqMass.cpp:190-202) returns `handler_.getNrSpectra()` when `sidx_.empty()`. The (parent, indices) ctor has the same empty==inherit-all behavior (lines 30-37).
- **Fix:** The behavior is documented, but the overload signature does not telegraph it. Prefer a named factory or an explicit 'all spectra' default-constructed marker so an empty selection cannot be confused with 'all'. Additive: add a factory method; do not change existing ctor semantics (ABI/behavior depended on by ChromatogramExtractor). Flag honestly: the ideal fix (empty => empty) would be breaking.
- **Verifier correction:** Two minor refinements to the claim, neither undermining the surprise: (1) The (parent, indices) ctor's empty case (cpp 30-33) inherits the PARENT's subset unchanged (sidx_ = sp.sidx_), not literally 'all file spectra' — though if the parent itself had an empty subset it transitively means all-file, so the empty==don't-narrow spirit is the same. (2) The claim that the behavior is 'depended on by ChromatogramExtractor' is overstated: ChromatogramExtractorAlgorithm.cpp:253-257 merely calls getNrSpectra() and returns early when < 1; nothing positively relies on empty==all. Consequently the empty=>empty fix would be behavior-breaking (changes observable results, so semantically breaking) but ABI-compatible — the signature itself is unchanged, so abi_impact for the recommended additive fix (a named factory) is none. The recommendation to add an explicit factory/marker rather than mutate ctor semantics is sound.

### [ANSW-41] MRMFeatureFilter::calculateIonRatio — calculateIonRatio returns 0.0 sentinel on missing data and divides without guarding against zero denominator
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h · _osw-mrm-feature_

```cpp
double calculateIonRatio(const Feature& component_1, const Feature& component_2, const std::string& feature_name) const
```
- **Expectation:** A function named calculateIonRatio returns the ratio of two transition intensities; a caller would expect either a valid ratio or a clear error/NaN when it cannot be computed.
- **Actual:** When neither component has the requested key/native_id the function silently returns the sentinel 0.0 (indistinguishable from a real ratio of 0). When only component_1 has the value it returns component_1's raw value (NOT a ratio). When both values exist it computes `feature_1 / feature_2` with no check that feature_2 != 0, so a zero denominator yields +inf/NaN with no signal.
- **Evidence:** Impl src/openms/source/ANALYSIS/OPENSWATH/MRMFeatureFilter.cpp:911-952 — `double ratio = 0.0; ... ratio = feature_1 / feature_2; ... else if (component_1.metaValueExists(feature_name)) { ... ratio = feature_1; } ... return ratio;`. The else branches only OPENMS_LOG_DEBUG and fall through to `return 0.0` / `return feature_1`.
- **Fix:** Document the 0.0/single-value fallbacks explicitly, and guard the division (return NaN or throw on feature_2==0). Additive/source-compatible fix: keep signature, add a denominator==0 guard returning std::numeric_limits<double>::quiet_NaN() and document the missing-key behavior; ABI unchanged.
- **Verifier correction:** Minor refinement: for the "intensity" branch the missing-data gate tests the "native_id" metavalue, not the intensity value (intensity is always available via getIntensity()); for other feature_names it tests the actual key. The recommended additive fix (denominator==0 guard returning quiet_NaN + documenting fallbacks) keeps the signature, so ABI is unchanged but it is a behavior change, hence source-compatible rather than none.

### [ANSW-77] TargetedSpectraExtractor::annotateSpectra — ppm 'mz_tolerance' is divided by 1e6 but never multiplied by m/z, yielding a tolerance ~1e6x too small in ppm mode
`high` · `unit-or-index` · ABI: `none` · src/openms/source/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.cpp · _osw-transitions-helper_

```cpp
void annotateSpectra(const std::vector<MSSpectrum>&, const TargetedExperiment&, std::vector<MSSpectrum>&, FeatureMap&, bool) const
```
- **Expectation:** With the parameter in ppm mode (mz_unit_is_Da_ == false), a tolerance of e.g. 20 ppm at m/z 500 should produce an absolute window of about +/-0.01 Th (ppm * mz / 1e6).
- **Actual:** The code computes 'mz_tolerance = mz_tolerance_ / 1e6' and then uses it directly as an absolute window around spectrum_mz, omitting the multiplication by the m/z. A 20 ppm setting becomes +/-2e-5 Th regardless of m/z, so essentially no spectra match. The same pattern appears for fwhm_threshold ('fwhm_threshold_ / 1e6' at line 462).
- **Evidence:** const double mz_tolerance = mz_unit_is_Da_ ? mz_tolerance_ : mz_tolerance_ / 1e6;  // line 370\n...\nconst double mz_left_lim = spectrum_mz ? spectrum_mz - mz_tolerance : ...;  // line 373
- **Fix:** In ppm mode compute 'mz_tolerance = mz_tolerance_ * spectrum_mz / 1e6' (ppm is relative to the measured m/z). Behavioral fix, ABI-neutral. Add a regression test asserting the ppm window scales with m/z.
- **Verifier correction:** In ppm mode (mz_unit_is_Da_ == false), line 370 should compute the absolute tolerance as mz_tolerance_ * spectrum_mz / 1e6 (i.e. Math::ppmToMass(mz_tolerance_, spectrum_mz)), not mz_tolerance_ / 1e6. As written, the matching window is the m/z-independent value mz_tolerance_/1e6 Th, which is too small by a factor equal to the precursor m/z (roughly 100-2000x in practice, NOT ~1e6x as the title states), so ppm-mode annotation silently matches almost nothing. Fix is a function-body change only (ABI-neutral). The fwhm_threshold / 1e6 line (462) is a separate, weaker concern: FWHM is an absolute width threshold, so the "missing * mz" reasoning does not cleanly transfer.

### [ANSW-1] SwathWindowLoader::readSwathWindows — readSwathWindows silently returns empty vectors when the file is missing or unreadable
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h · _osw-workflow_

```cpp
static void readSwathWindows(const std::string& filename, std::vector<double>& swath_prec_lower, std::vector<double>& swath_prec_upper)
```
- **Expectation:** A static reader named readSwathWindows, taking a filename, throws (e.g. FileNotFound) when the file does not exist or cannot be opened, so a caller that mistyped a path or pointed at a missing window file learns about it.
- **Actual:** The implementation does `std::ifstream data(filename.c_str());` and immediately `std::getline(data, line);` with no open()/good() check. A missing or unreadable file produces a failed stream; the getline fails, the do/while body never runs (or runs once on garbage), and the function returns leaving both output vectors empty with no error. The header documents this: 'A missing or unreadable file is not reported as an exception; the output vectors are left empty.'
- **Evidence:** SwathWindowLoader.cpp:93-99: `std::ifstream data(filename.c_str()); ... std::getline(data, line);` with no failure check. Header doc (SwathWindowLoader.h:101-104): 'A missing or unreadable file is not reported as an exception; the output vectors are left empty.'
- **Fix:** Additive/source-compatible fix: after constructing the ifstream, check `if (!data) throw Exception::FileNotFound(...)`. The sibling annotateSwathMapsFromFile would then surface a clear error for a bad path instead of failing the downstream count check with a confusing message. If preserving silent behavior is mandatory for ABI/back-compat, at minimum emit a WARN; the doc already flags this as surprising.
- **Verifier correction:** readSwathWindows does lack any open/good check on the ifstream (real silent-failure defect), but the precise consequence claimed is incorrect. For a missing/unreadable file, getline fails leaving line empty; StringUtils::split on the empty string returns early leaving headerSubstrings empty (StringUtils.h:761); the immediately following StringUtils::toDouble(headerSubstrings[0]) at SwathWindowLoader.cpp:104 indexes [0] on an EMPTY vector — undefined behavior (crash under hardened STL, garbage otherwise), NOT a clean return with empty output vectors. The try/catch only catches Exception::ConversionError, so the UB is not converted into the 'silently empty vectors' the claim and the header doc (lines 101-104) both promise. The header documentation is therefore also wrong about its own behavior. Net: real POLS/silent-failure defect, but the outcome is worse (UB) than 'empty vectors, no error'. Fix is additive and source-compatible: after `std::ifstream data(filename.c_str());` add `if (!data) throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);` — no signature/ABI change.

### [ANSW-5] SwathQC::getSpectraProcessingFunc — Uniform subsampling index counts only sampled spectra, not all spectra seen, biasing which spectra get analyzed
`high` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SwathQC.h · _osw-workflow_

```cpp
std::function<void(const OpenMS::MSSpectrum&)> getSpectraProcessingFunc()
```
- **Expectation:** The member-function path is documented as the speed-efficient way to 'sample the spectra as they are loaded'; uniform subsampling implies the decision uses the running index of ALL MS1 spectra streamed in.
- **Actual:** The returned lambda calls isSubsampledSpectrum_(nr_ms1_spectra_, cd_spectra_, ms1_spectra_seen_) and only increments ms1_spectra_seen_ AFTER a spectrum passes the filter (it is incremented inside the accepted branch, line 69, never on rejection). So ms1_spectra_seen_ is the count of ACCEPTED spectra, but it is used as the global spectrum index for uniform spacing. Because rejected spectra never advance the index, the index never reaches the spacing thresholds for later spectra the way the algorithm assumes — sampling is front-loaded/non-uniform. The private member is even documented as 'keeps track of number of spectra passed to getSpectraProcessingFunc()', which is not what it counts. The static getChargeDistribution avoids the bug entirely by doing its own indexing with the true loop variable i (SwathQC.cpp:131) and forcing sample-all mode.
- **Evidence:** SwathQC.cpp:64-69: `if (!isSubsampledSpectrum_(nr_ms1_spectra_, cd_spectra_, ms1_spectra_seen_)) { return; } ++ms1_spectra_seen_;` (increment only on the accepted path). Contrast SwathQC.cpp:127-134 (static path) which iterates `for (size_t i...)` and calls `isSubsampledSpectrum_(nr_spec, nr_samples, i)` with the true index. Member doc SwathQC.h:139 'keeps track of number of spectra passed to getSpectraProcessingFunc()'.
- **Fix:** Increment a true 'seen' counter for every MS1 spectrum entering the lambda and pass THAT as idx to isSubsampledSpectrum_ (separate from the accepted count). This is an internal-impl fix with no signature/ABI change. Note this only manifests when nr_ms1_spectra_>0 (external getExpSettingsFunc/setNrMS1Spectra path); the in-tree static caller is unaffected, which is likely why it went unnoticed.
- **Verifier correction:** Minor correction to the claim's 'went unnoticed' explanation: it is not strictly true that 'the in-tree caller is unaffected.' A member-path unit test exists (SwathQC_test.cpp storeJSON section) and exercises nr_ms1_spectra_>0, but it is masked because the test data has only 3 MS1 spectra against a subsample target of 10 (subsample>=total => isSubsampledSpectrum_ always returns true, so ms1_spectra_seen_ advances unbroken and the bug is invisible). The bug fires only when total MS1 spectra exceed the subsample count, which is the real-world OpenSwathWorkflow/OpenSwathFileSplitter regime (cd_spectra_=30, hundreds-to-thousands of MS1 spectra). In that regime the sampling collapses from the intended 30 uniformly-spaced spectra to exactly 1. Otherwise the claim, its quoted line references, and its internal-fix recommendation (track a true 'seen' counter incremented for every MS1 spectrum and pass that as idx) are accurate.

### [ANSW-46] ReactionMonitoringTransition::getPrediction — getPrediction()'s guard checks the WRONG predicate (hasPrecursorCVTerms instead of hasPrediction), and in release builds dereferences a possibly-null pointer
`high` · `silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/MRM/ReactionMonitoringTransition.cpp · _targeted-experiment_

```cpp
const Prediction & getPrediction() const
```
- **Expectation:** A caller who follows the header's instruction ('You first need to check whether the object is accessible using hasPrediction()') expects getPrediction() to be safe whenever hasPrediction() is true, and to be guarded against the null case.
- **Actual:** The precondition is `OPENMS_PRECONDITION(hasPrecursorCVTerms(), "...has no Prediction object, check first with hasPrediction()")` -- it checks hasPrecursorCVTerms(), an unrelated member, due to copy-paste from getPrecursorCVTermList(). The function then `return *prediction_;`. OPENMS_PRECONDITION compiles to nothing unless OPENMS_ASSERTIONS is set, so in a normal release build getPrediction() unconditionally dereferences prediction_, which is nullptr by default. A transition with no prediction (the default) gives an immediate null-deref crash, and even a debug build checks the wrong flag.
- **Evidence:** Impl: `const ReactionMonitoringTransition::Prediction & ReactionMonitoringTransition::getPrediction() const { OPENMS_PRECONDITION(hasPrecursorCVTerms(), "ReactionMonitoringTransition has no Prediction object, check first with hasPrediction()") return *prediction_; }`; default `prediction_(nullptr)` in ctor.
- **Fix:** Fix the precondition to `OPENMS_PRECONDITION(hasPrediction(), ...)`. ABI-safe (function body change only). Ideally also harden against null in release (e.g. throw Exception::InvalidValue when prediction_ == nullptr) for parity with the documented contract; that is also source/ABI compatible.
- **Verifier correction:** Claim is accurate as stated. Minor clarification: even in a DEBUG (OPENMS_ASSERTIONS) build the bug is harmful — the wrong predicate can throw spuriously when a prediction exists but precursor CV terms do not, and can fail to throw (then null-deref) when a prediction is absent but precursor CV terms are present. In RELEASE the precondition compiles to nothing, so any default transition (prediction_ == nullptr) hits an unconditional null dereference at 'return *prediction_;'. Fix: change line 322 guard to hasPrediction(); optionally also throw on null in release for parity with the header contract. Both fixes are body-only and ABI-safe.

### [ANSW-47] TargetedExperiment::getProteinByRef / getPeptideByRef / getCompoundByRef — By-ref getters return a reference into a map populated with operator[], silently inserting/returning a null entry for an unknown ref in release builds
`high` · `silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/TARGETED/TargetedExperiment.cpp · _targeted-experiment_

```cpp
const Protein & getProteinByRef(const std::string & ref) const
```
- **Expectation:** getProteinByRef(ref) for an unknown ref should fail loudly (throw / assert), as the message implies; at worst it should not corrupt internal state.
- **Actual:** The only guard is `OPENMS_PRECONDITION(protein_reference_map_.contains(ref), "Could not find protein in map")`, which is compiled out in release. The next line is `return *(protein_reference_map_[ref]);`. For an unknown ref, std::unordered_map::operator[] INSERTS a default-constructed entry (a null `const Protein*`) and then the code dereferences that nullptr -> crash/UB. So an unknown ref is not a clean lookup failure but a null-deref plus a mutation of the cache map. getPeptideByRef/getCompoundByRef have the identical pattern.
- **Evidence:** `OPENMS_PRECONDITION(protein_reference_map_.contains(ref), "Could not find protein in map") return *(protein_reference_map_[ref]);` (and equivalent for peptide/compound using hasPeptide/hasCompound preconditions).
- **Fix:** Replace operator[] with `.at(ref)` or an iterator find + explicit throw (e.g. Exception::ElementNotFound) so an unknown ref is a deterministic error in all build types and does not mutate the cache. Body-only change; ABI-safe. Callers already must check has*() first, so behavior for valid refs is unchanged.

### [ANSW-17] ChromatogramExtractor::extractChromatograms — "bartlett" filter is advertised and accepted but throws NotImplemented at extraction time
`medium` · `documentation-mismatch` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h · _osw-chrom_

```cpp
void extractChromatograms(const OpenSwath::SpectrumAccessPtr input, std::vector<OpenSwath::ChromatogramPtr>& output, const std::vector<ExtractionCoordinates>& extraction_coordinates, double mz_extraction_window, bool ppm, const std::string& filter)
```
- **Expectation:** The class doc says "Multiple filter types available (Bartlett, tophat)" and the @param doc says filter is "Which filter to use (bartlett or tophat)". A caller passing filter="bartlett" expects a working Bartlett-windowed extraction.
- **Actual:** getFilterNr_("bartlett") returns 2 (passes validation), but ChromatogramExtractorAlgorithm::extractChromatograms only implements used_filter==1 (tophat) and the used_filter==2 branch does `throw Exception::NotImplemented(...)`. "bartlett" therefore passes argument validation and then fails mid-run. (The sibling ChromatogramExtractorAlgorithm.h is honest: "Which function to apply in m/z space (currently \"tophat\" only)".)
- **Evidence:** ChromatogramExtractorAlgorithm.cpp getFilterNr_: `else if (filter == "bartlett") { return 2; }`; extractChromatograms: `else if (used_filter == 2) { throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }`. ChromatogramExtractor.h class doc line 39: "Multiple filter types available (Bartlett, tophat)" and @param line 73: "Which filter to use (bartlett or tophat)".
- **Fix:** Source-compatible fix: update the ChromatogramExtractor.h class doc and the two @param descriptions to state tophat is the only implemented filter (matching ChromatogramExtractorAlgorithm.h), and make getFilterNr_ reject "bartlett" with a clear "not implemented" message at validation time so failure is immediate and explicit rather than mid-extraction. No ABI change.
- **Verifier correction:** "bartlett" is documented as a usable filter in ChromatogramExtractor.h (class doc line 39 and @param on lines 73/101) but is unimplemented: getFilterNr_ accepts it (returns 2) so validation passes, and the extraction loop's used_filter==2 branch throws Exception::NotImplemented mid-run. The failure is a LOUD exception (not silent data corruption), so the accurate category is documentation-mismatch / deferred-validation, severity medium. Fix is source-compatible (no ABI change): correct the ChromatogramExtractor.h doc to state tophat is the only implemented filter (matching ChromatogramExtractorAlgorithm.h) and make getFilterNr_ reject "bartlett" at validation time so the failure is immediate and explicit.

### [ANSW-18] ChromatogramExtractorAlgorithm::ExtractionCoordinates::ion_mobility — Default ion_mobility = 0.0 is treated as a real IM extraction target, not "no IM"
`medium` · `surprising-default` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h · _osw-chrom_

```cpp
double ion_mobility = 0.0; ///< ion mobility value around which should be extracted
```
- **Expectation:** A default-constructed ExtractionCoordinates (ion_mobility = 0.0) should mean "no ion-mobility constraint", so a caller who only fills mz/rt expects normal m/z extraction.
- **Actual:** The extraction gate is `use_im = (extraction_coordinates[k].ion_mobility >= 0.0 && has_im)` and return_chromatogram uses `if (coord.ion_mobility >= 0 && im_extraction_width > 0.0)`. The 'IM disabled' sentinel is therefore a NEGATIVE value, but the struct default is 0.0 (which is >= 0.0). So whenever an IM array is present, a default coordinate is silently extracted as a valid IM target at drift time 0.0, and return_chromatogram annotates a drift time of 0.0.
- **Evidence:** ExtractionCoordinatesAlgorithm.cpp: `const bool use_im = (extraction_coordinates[k].ion_mobility >= 0.0 && has_im);`. ChromatogramExtractor.h line 245: `if (coord.ion_mobility >= 0 && im_extraction_width > 0.0)`. Header default: `double ion_mobility = 0.0;`.
- **Fix:** Source-compatible: change the in-class default to a negative sentinel (e.g. `double ion_mobility = -1.0;`) so the documented/used 'IM off' sentinel matches the default, and document that negative means 'no IM'. Changing only the default initializer is source-compatible (no signature change); review call sites that explicitly set 0.0 expecting IM-at-zero. Lowest-risk alternative: document the sentinel convention prominently on the member.
- **Verifier correction:** The default ion_mobility = 0.0 conflicts with the module's disabled-IM sentinel of -1 (used by LightCompound::drift_time and propagated through prepare_coordinates). A default-constructed ExtractionCoordinates is therefore treated as a valid IM target at drift time 0.0, not as 'no IM'. However the misuse requires a hand-built coordinate (default 0.0 not overwritten) combined with an explicitly positive im_extraction_window, since IM extraction is also gated on im_extraction_window > 0.0; the standard prepare_coordinates path is safe because it sets ion_mobility = getDriftTime() = -1. Fix is source-compatible: change the in-class initializer to `double ion_mobility = -1.0;` to align with the documented/used 'IM off' sentinel, document that negative means 'no IM', and review any call site that deliberately sets 0.0 to extract at drift time 0. Note the now-dead `if (ion_mobility < 0) warn` branch at ChromatogramExtractorAlgorithm.cpp:335 should be re-examined.

### [ANSW-24] SpectrumAccessQuadMZTransforming::getSpectrumById — getSpectrumById mutates the m/z values of the spectrum buffer in place, corrupting cached spectra of an in-memory backing store
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h · _osw-dataaccess_

```cpp
OpenSwath::SpectrumPtr getSpectrumById(int id) override
```
- **Expectation:** A `getSpectrumById` accessor returns a (logically read-only) view/copy of the requested spectrum; calling it twice returns the same data, and it does not mutate the underlying store.
- **Actual:** The transform writes back into the SpectrumPtr returned by the wrapped accessor: `s->getMZArray()->data[i] = predict;`. When the underlying accessor is SpectrumAccessOpenMSInMemory, getSpectrumById hands back the *cached* shared SpectrumPtr (`return spectra_[id];`), so the transform permanently overwrites the cached m/z array; a second call re-applies the quadratic to already-transformed m/z, and any other holder of that shared spectrum sees corrupted m/z.
- **Evidence:** SpectrumAccessQuadMZTransforming.cpp:39-58 mutates `s->getMZArray()->data[i]`; SpectrumAccessOpenMSInMemory.cpp:59-64 `getSpectrumById` returns the shared cached `spectra_[id]` without copying.
- **Fix:** Deep-copy the m/z array before transforming (allocate a fresh BinaryDataArray, copy then transform), so the call is idempotent and side-effect-free regardless of backing store. This is an internal .cpp fix with no ABI change.
- **Verifier correction:** getSpectrumById is not strictly read-only: it transforms m/z in place on whatever SpectrumPtr the wrapped accessor returns. This is only safe when the backing accessor hands back a freshly allocated spectrum per call (e.g. SpectrumAccessOpenMS, which the unit test exercises). With a caching backing store such as SpectrumAccessOpenMSInMemory (which returns the shared cached spectra_[id] without copying), the call permanently mutates the cached m/z array seen by every holder and is non-idempotent — calling it twice for the same id applies the quadratic correction twice. Fix by deep-copying the m/z BinaryDataArray (and constructing a fresh Spectrum) before transforming, making the accessor side-effect-free and idempotent regardless of backing store. This is an internal .cpp change with no public-signature/ABI impact.

### [ANSW-28] MRMFeatureOpenMS::getMetaValue — getMetaValue returns double and throws ConversionError for non-numeric meta values
`medium` · `surprising-throw` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h · _osw-dataaccess_

```cpp
double getMetaValue(std::string name) const
```
- **Expectation:** A generic getMetaValue(name) returns the stored value; for a meta value that is a string or list, it should not throw merely because the declared return type is double.
- **Actual:** The body returns `mrmfeature_.getMetaValue(name)`, a DataValue, implicitly converted to double. If the named meta value is a String (or missing), DataValue's conversion throws Exception::ConversionError (or the key-miss path throws), which a caller of a plain 'getMetaValue' would not anticipate. The double return type is dictated by the OpenSwath::IMRMFeature interface but is not reflected/warned in this class's header.
- **Evidence:** MRMFeatureAccessOpenMS.cpp:329-332 `double MRMFeatureOpenMS::getMetaValue(std::string name) const { return mrmfeature_.getMetaValue(name); }`; MetaInfoInterface::getMetaValue returns DataValue (MetaInfoInterface.h:71).
- **Fix:** Document that name must reference a numeric meta value and that a missing/non-numeric value throws. Note this method is an override of IMRMFeature::getMetaValue but the header omits the `override` keyword (line 104), which masks the interface contract; add `override`. The `override` add is source-compatible and ABI-neutral.
- **Verifier correction:** getMetaValue throws Exception::ConversionError only on the MISSING-KEY path: MetaInfoInterface::getMetaValue returns DataValue::EMPTY for an absent name, and DataValue::operator double() throws on EMPTY_VALUE. It does NOT throw for a String/list meta value — operator double() falls through to `return data_.dou_;`, yielding undefined behavior / silent garbage from the inactive union member. So the hazards are: (1) undocumented throw on missing key, and (2) silent UB (not a throw) on non-numeric values. The header should document that name must reference an existing numeric meta value, and `override` should be added on line 104 (source-compatible, ABI-neutral).

### [ANSW-70] OpenSwathOSWParquetWriter::write — write() to an existing .oswpq archive deletes prior runs by default despite docs promising append
`medium` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.h · _osw-export_

```cpp
void write(const std::string& output_path, const OpenSwath::LightTargetedExperiment& assay_library, const FeatureMap& feature_map, UInt64 run_id, const std::string& input_filename, bool enable_uis_scoring) const
```
- **Expectation:** Per the class and method docs: "If the output directory already exists, the writer will append a new run partition ... Existing runs are not modified" / "the new run is appended under a fresh run_id= partition".
- **Actual:** Appending only happens for an existing *directory* output, or for an archive when setPreserveExisting(true) was called. preserve_existing_ defaults to false (header line 88). For an existing *archive file* with the default, write() ignores the old contents (OpenSwathOSWParquetWriter.cpp:316-329 routes to a fresh temp dir) and at the end removes the existing zip: 'if (File::exists(output_zip_abs) && !preserve_existing_) { File::remove(output_zip_abs); }' (cpp:1317-1319). All previously written runs are silently lost.
- **Evidence:** header doc: "append a new run partition under runs/run_id=<id>/ ... Existing runs are not modified."; cpp:1317 'if (File::exists(output_zip_abs) && !preserve_existing_) File::remove(output_zip_abs);'
- **Fix:** Document on write() that for archive outputs the default overwrites/replaces the archive, and that setPreserveExisting(true) is required for the append behavior the class doc describes; or flip the default to preserve+append for consistency with the documented directory behavior. ABI-safe doc change; the default flip is source-compatible but alters runtime semantics.
- **Verifier correction:** For archive (.oswpq file) outputs, write() with the default (preserve_existing_ = false) does NOT append: it ignores the existing archive's contents (cpp:314-329 routes to a fresh temp dir, bypassing the run_id-collision guard at cpp:396) and deletes the prior archive at cpp:1317-1319 before re-zipping only the new run, so all previously written runs are dropped. This directly contradicts the unconditional append/"Existing runs are not modified" promise in the write() and class Doxygen (h:47-49, 59-60). For DIRECTORY outputs the documented append behavior does hold. The real caller (OpenSwathWorkflow) mitigates via the --append_oswpq flag whose help states "(default: overwrite)", so the data loss is documented at the tool level and recoverable, but the library-level write() doc remains misleading. Fix: scope the method/class Doxygen to state that archive outputs overwrite by default and require setPreserveExisting(true) for append, OR make archive and directory outputs consistent.

### [ANSW-71] OpenSwathOSWParquetReader::fetchPeakGroupFeatures / fetchTransitionFeatures / fetchUnscoredData — fetch* methods ignore the path loaded by the constructor/load() and require the path to be re-supplied
`medium` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetReader.h · _osw-export_

```cpp
PeakGroupFeatureScoresResult fetchPeakGroupFeatures(const std::string& oswpq_dir, const std::string& level="ms2", const std::string& main_score="") const
```
- **Expectation:** The convenience constructor doc says "the object is ready to use after construction" and oswpqPath() exists "so Python-side code can call fetch methods without re-supplying the path", implying fetch* operate on the already-loaded archive.
- **Actual:** Each fetch* takes a mandatory oswpq_dir argument (no default) and reads exclusively from that argument; the stored oswpq_dir_ (set by load()/the constructor) is never consulted by fetch* (OpenSwathOSWParquetReader.cpp:143-164 etc. open from the passed-in oswpq_dir). load() also populates rows_ which fetch* do not use. A caller who built the reader with the convenience constructor must still pass the path again, and could even pass a *different* path than the one loaded.
- **Evidence:** header: "It is provided for Python ergonomics so callers can create an instance with a single argument"; oswpqPath() doc: "so Python-side code can call fetch methods without re-supplying the path"; yet fetchPeakGroupFeatures(const std::string& oswpq_dir, ...) is required and the cpp uses oswpq_dir directly.
- **Fix:** Add overloads (or default the argument to empty and fall back to oswpq_dir_ when empty) so fetch* honor the loaded path, matching the documented ergonomics. Additive zero-path overloads are source-compatible; otherwise a doc fix.
- **Verifier correction:** Claim stands. Minor precision: the mandatory argument means omitting the path is a compile error (loud), not a silent failure; the real silent hazard is that a fetch call can read a path that differs from the one the constructor/load() loaded, while the stored oswpq_dir_ and rows_ remain unused. Recommended fix (default oswpq_dir="" and fall back to oswpq_dir_ when empty) is source-compatible; an additive zero-/loaded-path overload would be fully source-compatible; a doc-only correction has no ABI impact.

### [ANSW-72] OpenSwathOSWWriter::prepareRows / prepareRowsInto / prepareLine — pep and transition parameters are documented as inputs but are completely ignored
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h · _osw-export_

```cpp
OSWData prepareRows(const OpenSwath::LightCompound& pep, const OpenSwath::LightTransition* transition, const FeatureMap& output, const std::string& id) const
```
- **Expectation:** From the docs, @p pep ("The compound used for extraction") and @p transition ("The transition used for extraction") influence the emitted rows.
- **Actual:** Both parameters are ignored in the implementation: prepareRows forwards a default-constructed LightCompound and nullptr to prepareRowsInto (OpenSwathOSWWriter.cpp:1008-1017), and prepareRowsInto's pep/transition are commented out and unused (cpp:1019-1023). The transition pointer may even be null with no effect. Callers passing a real transition expecting per-transition output silently get identical results regardless.
- **Evidence:** cpp:1015 'prepareRowsInto(rows, OpenSwath::LightCompound(), nullptr, output, id);' and signatures 'const OpenSwath::LightCompound& /* pep */, const OpenSwath::LightTransition* /* transition */'
- **Fix:** Either remove the dead parameters (source-breaking; provide deprecated forwarding overloads) or document them clearly as currently ignored / reserved. Minimal ABI-safe fix: header note that pep/transition are unused.
- **Verifier correction:** The pep and transition parameters of prepareRows/prepareRowsInto/prepareLine are genuinely unused (commented out in signatures and never referenced; prepareRows forwards LightCompound()/nullptr to prepareRowsInto), yet the public Doxygen documents them as meaningful extraction inputs — a real doc/implementation mismatch. However, severity is medium, not high: ignoring them does not produce silently-wrong or lost data (all rows are correctly derived from the FeatureMap and member run_id_), and the sole in-tree caller already passes empties with an explicit "currently unused" comment. The defect is misleading documentation that invites no-op misuse, not incorrect results. ABI impact is none (signatures unchanged); the minimal fix is a header note that pep/transition are currently ignored/reserved.

### [ANSW-75] OpenSwathMatrixExporter::buildMatrix / writeMatrix — Matrix exporter throws on empty input while its sibling exporters explicitly allow empty exports
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathMatrixExporter.h · _osw-export_

```cpp
static OpenSwathQuantMatrix buildMatrix(const std::vector<OpenSwathExportRow>& rows, const OpenSwathMatrixExportConfig& config); static void writeMatrix(...)
```
- **Expectation:** Sibling classes in the same cluster document empty handling: OpenSwathResultsExporter::write and OpenSwathParquetExporter::writeFeatureScores/writeTransitionScores state 'Empty exports are allowed and result in a header-only TSV or an empty Parquet table.' A developer would expect buildMatrix/writeMatrix to behave the same.
- **Actual:** buildMatrix throws Exception::Precondition when rows is empty (cpp:284-288) and on empty peptide/protein/gene summarization (cpp:347-351, 381-386); writeMatrix throws on an empty matrix (cpp:448-452). None of this throw-on-empty behavior is mentioned in OpenSwathMatrixExporter.h, so callers porting from the results/parquet exporters get an unexpected exception instead of an empty output.
- **Evidence:** cpp:286 'throw Exception::Precondition(... "No OpenSWATH rows were available for matrix export.")'; cpp:450 'throw ... "Cannot write an empty OpenSWATH quantification matrix."' vs OpenSwathResultsExporter.h "Empty exports are allowed".
- **Fix:** Document in OpenSwathMatrixExporter.h that empty input throws (and why), making the inconsistency with the other exporters explicit; or align behavior to emit a header-only/empty matrix. Doc-only fix is ABI-safe.
- **Verifier correction:** Severity is medium, not high: the failure is loud (a thrown Exception::Precondition, not silent wrong results or data loss) and recoverable (callers can catch or pre-check emptiness). It rises above 'low/mild surprise' because it can abort a real OpenSwathExport TOPP run on legitimate input (all rows excluded by user filters) where the Results/Parquet branches in the same tool succeed, and because the sibling cluster's explicitly documented empty-tolerance contract actively misleads a developer porting between these APIs. The recommended fix (document the throw in OpenSwathMatrixExporter.h, or align behavior to emit a header-only/empty matrix) is ABI-safe (doc-only or behavior-only; no signature/symbol/layout change).

### [ANSW-21] OpenSwathGeneInference::infer / OpenSwathPeptideInference::infer / OpenSwathProteinInference::infer — The level-specific wrappers throw if config.level does not match the class, instead of inferring at the class's level
`medium` · `surprising-throw` · ABI: `none` · src/openms/source/ANALYSIS/OPENSWATH/OpenSwathGeneInference.cpp · _osw-inference_

```cpp
std::vector<LevelContextResultRow> infer(const std::vector<LevelContextInputRow>& input, const LevelContextInferenceConfig& config) const
```
- **Expectation:** Calling OpenSwathGeneInference::infer(input, config) would perform gene-level inference; the class name already encodes the level, so config.level would be ignored or defaulted.
- **Actual:** The wrapper ignores its own identity and requires the caller to ALSO set config.level == InferenceLevel::Gene (resp. Peptide/Protein); a mismatch throws Exception::Precondition. The class adds no behavior beyond this redundant assertion plus delegation to LevelContextInference::infer. A caller who reuses one config across the three wrappers gets a throw rather than gene/peptide/protein results.
- **Evidence:** OpenSwathGeneInference.cpp:19-23 `if (config.level != InferenceLevel::Gene) { throw Exception::Precondition(..., "OpenSwathGeneInference requires config.level = InferenceLevel::Gene."); } return LevelContextInference::infer(input, config);` (identical pattern in OpenSwathPeptideInference.cpp:19 and OpenSwathProteinInference.cpp:19)
- **Fix:** Make the wrappers authoritative: copy config and force cfg.level to the class level before delegating (additive, source-compatible behavior change that removes a surprising throw), or document prominently that config.level is a required, redundant precondition. Forcing the level is the least-surprising fix and is ABI-safe.
- **Verifier correction:** Severity is medium, not high: the precondition fails loud and immediate (Exception::Precondition), never silently produces wrong results or data loss. It is an undocumented, identity-contradicting precondition that traps a default-constructed config (default level=Peptide) for the Gene and Protein wrappers and invites cross-wrapper config-reuse misuse, but every failure is recoverable and observable. The recommended fix (copy config and force cfg.level to the class level before delegating) touches only the function body, so abi_impact is none (source-compatible behavior change), not breaking.

### [ANSW-22] LevelContextInference::infer (non-finite score handling) — Rows with non-finite score silently get NaN p/q/PEP while an all-NaN distribution throws
`medium` · `silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/OPENSWATH/LevelContextInference.cpp · _osw-inference_

```cpp
static std::vector<LevelContextResultRow> infer(const std::vector<LevelContextInputRow>& input, const LevelContextInferenceConfig& config)
```
- **Expectation:** Given the function throws Precondition when the target or decoy distribution is empty or all-NaN, a caller would expect consistent error signaling for malformed scores; the header documents only normal p/q/PEP outputs and gives no hint that individual rows can come back NaN.
- **Actual:** An individual row whose score is non-finite is silently written with NaN pvalue/qvalue/pep and skipped from matching, whereas a wholly-NaN target/decoy distribution throws. The header doc ('Result rows containing the original entity and score plus derived p-value, q-value, and PEP estimates') does not mention the NaN sentinel, so a caller can propagate NaN q-values undetected.
- **Evidence:** LevelContextInference.cpp:318-324 `if (!std::isfinite(row.score)) { result.pvalue = quiet_NaN(); result.qvalue = quiet_NaN(); result.pep = quiet_NaN(); continue; }` vs buildErrorStatistics_ throwing at LevelContextInference.cpp:248-255.
- **Fix:** Document the NaN-on-non-finite-score behavior in the header @return text so callers know to check std::isfinite on returned q-values. Doc-only, ABI-safe.
- **Verifier correction:** The behavior and code reading are accurate, but severity is medium rather than high: only rows with already non-finite scores get the NaN sentinel; valid rows are unaffected, and a caller can recover by checking std::isfinite on returned q-values. The genuine defect is the undocumented, inconsistent error signaling (silent per-row NaN vs. throwing per-distribution), which can let NaN q-values propagate undetected through downstream threshold filters. Fix is documentation-only in the header @return text (ABI-safe).

### [ANSW-61] MRMAssay::detectingTransitions / restrictTransitions / uisTransitions (and *Light variants) — Doxygen tags exp as @param[in] (input) but the methods overwrite it in place
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h · _osw-mrm-assay_

```cpp
void detectingTransitions(OpenMS::TargetedExperiment& exp, int min_transitions, int max_transitions); void restrictTransitions(OpenMS::TargetedExperiment& exp, ...); void uisTransitions(OpenMS::TargetedExperiment& exp, ...)
```
- **Expectation:** A parameter documented '@param[in] exp the input, unfiltered transitions' (line 111-112 for detectingTransitions, line 97 for restrictTransitions, line 138 for uisTransitions) reads as a read-only input; the result must be returned somewhere else.
- **Actual:** Every one of these methods rebuilds a transition list and writes it back with exp.setTransitions(...), destroying the caller's original library. There is no separate output: the in-param IS the output. reannotateTransitions even tags the same parameter '@param[out]' (line 75), so the convention is internally contradictory across siblings of the same class.
- **Evidence:** MRMAssay.cpp:831 `exp.setTransitions(std::move(transitions));` (reannotateTransitions), :888 `exp.setTransitions(std::move(transitions));` (restrictTransitions), :1072 `exp.setTransitions(transitions);` (uisTransitions). Header MRMAssay.h:75 `@param[out] exp` vs :97/:111/:138 `@param[in] exp`.
- **Fix:** ABI-neutral: fix the Doxygen to '@param[in,out] exp ... modified in place; transition list is replaced' on all of detectingTransitions/restrictTransitions/uisTransitions/reannotateTransitions and their Light variants, matching the (already in-place) behavior. No signature change needed; this is a docs-only fix that removes the surprise.
- **Verifier correction:** The claim is correct on the facts but should be re-graded medium (not high): no silent wrong results, no crash, no data loss under correct use — the non-const TargetedExperiment& is itself a mutation signal and the replacement happens on the caller's own still-owned object. The fix is docs-only: relabel exp as '@param[in,out] exp ... modified in place; the transition (and peptide/protein) lists are replaced' on reannotateTransitions (H75), restrictTransitions (H97), detectingTransitions (H111), uisTransitions (H138) and their Light variants (H198, H223, H239; H253 is already correct). Additionally fix the unrelated copy-paste tags where the by-value 'int min_transitions' is wrongly tagged '@param[out]' (H112, H239) — should be '@param[in]'.

### [ANSW-63] MRMDecoy::shufflePeptide / shufflePeptideLight — Returns a decoy that may still exceed identity_threshold (or be identical to target) after max_attempts, with no signal
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h · _osw-mrm-assay_

```cpp
OpenMS::TargetedExperiment::Peptide shufflePeptide(OpenMS::TargetedExperiment::Peptide peptide, const double identity_threshold, int seed = -1, const int max_attempts = 100) const
```
- **Expectation:** The Doxygen states 'This function will shuffle the given peptide sequences ... such that the resulting relative sequence identity is below identity_threshold.' A caller reasonably assumes the returned peptide always satisfies identity < identity_threshold.
- **Actual:** The loop runs at most max_attempts times; when the attempt budget is exhausted it returns whatever 'shuffled' currently holds, which can still have identity > identity_threshold (or, for sequences with too many fixed K/R/P residues, can equal the input). No exception, no flag, no bool. The caller silently gets a poor/identity decoy that inflates FDR estimates.
- **Evidence:** MRMDecoy.cpp:390 `while (AASequenceIdentity(peptide_seq, shuffled) > identity_threshold && attempts < max_attempts)` then unconditional `return std::make_pair(shuffled, shuffled_mods);` at :490. No post-loop check of the achieved identity.
- **Fix:** Additive (ABI-safe): keep the signature but document that the identity bound is best-effort and not guaranteed when max_attempts is exhausted. Optionally add an overload/out-param returning the achieved identity or a bool 'converged' so callers can detect failure. Do not silently change the existing return.
- **Verifier correction:** The function does silently return a decoy that can still exceed identity_threshold after max_attempts, with no exception/flag/bool (evidence at MRMDecoy.cpp:390-391 and :490 is exact and correct). However, the claim's "can equal the input" / "identity decoy" sub-case is inaccurate for the actual call sites: a decoy whose sequence exactly equals its target is caught by the exact-match dedup (allPeptideSequences pre-populated with all targets at MRMDecoy.cpp:542, checked at :605 / :958) and excluded — though only silently at DEBUG level. The real, demonstrable silent failure is a high-but-not-exactly-identical decoy (e.g. 70-90% identity) that is not an exact duplicate and is therefore added with no signal that it failed the documented identity bound. Severity is medium rather than high: the result is a degraded-quality decoy affecting FDR statistics for a minority of hard-to-shuffle short peptides, recoverable and not a crash/data-loss. Recommendation stands and is ABI-safe (none): fix the Doxygen to state the bound is best-effort and capped by max_attempts, and optionally add an overload/out-param exposing achieved identity or a 'converged' bool so callers can detect and drop non-converged decoys.

### [ANSW-64] MRMIonSeries::annotateIon / getIon — On no match returns sentinel ("unannotated", -1.0) instead of signaling failure; -1 is a valid-looking m/z
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMIonSeries.h · _osw-mrm-assay_

```cpp
std::pair<std::string, double> annotateIon(const IonSeries& ionseries, const double product_mz, const double mz_threshold); std::pair<std::string, double> getIon(IonSeries& ionseries, const std::string& ionid)
```
- **Expectation:** The header says the return is 'the annotation and product m/z of the queried fragment ion'. A caller expects a real annotation/m/z, and would not expect a magic m/z of -1 that can flow into a transition's product m/z.
- **Actual:** Both methods return std::make_pair("unannotated", -1) when nothing matches. The numeric -1 is an out-of-band sentinel that callers must know to test against the string "unannotated"; nothing in the type or name hints at this. reannotateTransitions actually does setProductMZ(targetion.second) (i.e. -1) before checking the string, relying on a downstream skip.
- **Evidence:** MRMIonSeries.cpp:30 `return make_pair(std::string("unannotated"), -1);` (getIon), :89 `ion = make_pair(unannotated, -1);` and :136 `return ion;` (annotateIon). MRMAssay.cpp:804 `tr.setProductMZ(targetion.second);` runs before the :807 unannotated check.
- **Fix:** ABI-neutral docs fix: document the sentinel return ('returns {"unannotated", -1.0} if no ion is within mz_threshold; callers must check the string before using the m/z'). Longer term (source-compatible additive) consider an overload returning std::optional. Keep existing signature for ABI.
- **Verifier correction:** Both getIon (MRMIonSeries.cpp:22-32) and annotateIon (MRMIonSeries.cpp:81-137) return std::make_pair("unannotated", -1.0) when no ion matches, and the header documentation does not mention this sentinel — a real POLS/silent-failure surprise. However, the strongest framing of the MRMAssay risk overstates current impact: at MRMAssay.cpp:804 the -1 IS written via setProductMZ before the :807 string check, but the following continue (:812) drops the transition so the -1 is never emitted; the second call site (MRMAssay.cpp:1332-1354) checks "unannotated" and continues before assigning product_mz, so it is unconditionally safe. Net: no current call site ships a -1 m/z, so this is a latent foot-gun (medium), not active data corruption (high). Recommended fix is docs-only (sentinel must be documented; callers MUST check string=="unannotated" before using the m/z), with an optional additive std::optional overload — both ABI-preserving (abi_impact: none for the docs fix).

### [ANSW-65] MRMDecoy::AASequenceIdentity — Out-of-bounds read when the two sequences differ in length (length check is debug-only assert)
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h · _osw-mrm-assay_

```cpp
float AASequenceIdentity(const std::string& sequence, const std::string& decoy) const
```
- **Expectation:** A public 'compute identity between two sequences' helper, documented only as 'Compute relative identity ... between two sequences', should handle (or clearly reject) sequences of unequal length.
- **Actual:** The equal-length requirement is enforced only by OPENMS_PRECONDITION, which compiles out in release builds. The loop iterates over sequence.size() and indexes decoy_v[i]; if decoy is shorter than sequence this is an out-of-bounds vector access (undefined behavior). The header documents no length precondition.
- **Evidence:** MRMDecoy.cpp:86 `OPENMS_PRECONDITION(sequence.size() == decoy.size(), ...)` then :91-93 `for (Size i = 0; i < sequence_v.size(); i++){ if (sequence_v[i] == decoy_v[i]) ...}` and :98 divides by sequence_v.size().
- **Fix:** Document the equal-length precondition in the header, and/or replace the debug-only OPENMS_PRECONDITION with a hard check (throw IllegalArgument or clamp to min length) so release builds do not read OOB. The header signature is unchanged (ABI-safe).
- **Verifier correction:** Real defect, but downgrade from high to medium. The equal-length requirement of public MRMDecoy::AASequenceIdentity is enforced only by OPENMS_PRECONDITION (debug-only; empty in release), so unequal-length inputs cause an unchecked decoy_v[i] out-of-bounds read (UB) in release builds, and the header documents no length precondition. However the sole internal caller (shufflePeptideLight) only ever passes equal-length strings and the method is not bound in pyOpenMS, so there is no currently-reachable crash from normal OpenMS usage; the exposure is to external/future C++ callers following the silent contract. Fix: document the equal-length precondition in the header and/or replace the debug-only OPENMS_PRECONDITION with a hard runtime check (throw IllegalArgument) or clamp to min length. Signature unchanged, ABI-safe.

### [ANSW-13] MRMFeatureFilter::calculateIonRatio — calculateIonRatio returns a non-ratio (numerator alone, or 0.0) instead of signaling a missing denominator
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h · _osw-mrm-feature_

```cpp
double calculateIonRatio(const Feature& component_1, const Feature& component_2, const std::string& feature_name) const
```
- **Expectation:** A method named calculateIonRatio returning the 'ratio between two transitions' should return component_1/component_2, or otherwise clearly signal that a ratio could not be computed (NaN/throw) when the second component lacks the requested value.
- **Actual:** When only component_1 has the value, the method returns component_1's raw value as the 'ratio' (ratio = feature_1; lines 926-927 and 942-943). When neither has the value, or feature_name=='intensity' and component_1 lacks native_id, it returns the initial sentinel ratio = 0.0. In all these cases the caller receives a finite double indistinguishable from a legitimate ratio.
- **Evidence:** double ratio = 0.0; ... else if (component_1.metaValueExists(feature_name)) { ... const double feature_1 = (double)component_1.getMetaValue(feature_name); ratio = feature_1; } else { OPENMS_LOG_DEBUG << "Feature metaValue ... not found ..."; } return ratio;
- **Fix:** Document and ideally return std::numeric_limits<double>::quiet_NaN() (or throw) when the denominator value is missing, instead of returning the numerator or 0.0. Callers (e.g. FilterFeatureMap line 164) compare the result against ion_ratio_l/ion_ratio_u bounds; a silent 0.0 or numerator-only value can pass/fail QC incorrectly. ABI: changing the return value on the missing-data path is source-compatible (signature unchanged) but a behavioral change, so guard behind a documented note or a new overload if strict back-compat is needed.
- **Verifier correction:** Severity is medium rather than high: the function does not crash or corrupt data wholesale; it produces a silently-wrong QC pass/fail decision only on the specific missing-denominator path (one transition has the feature value, the other does not — e.g. a missing internal standard or undetected second transition). The numerator-only return (lines 926-927, 942-943) is the clearer defect (a bare intensity/metaValue masquerading as a ratio); the neither-found 0.0 return is a milder surprise but still collides with a legitimate ratio of 0.0. Recommendation stands: return std::numeric_limits<double>::quiet_NaN() (which checkRange would then reject, failing QC loudly) on any missing-denominator path, and document the behavior; the change is source-compatible (signature unchanged) but behavioral.

### [ANSW-16] MRMFeatureFinderScoring::prepareProteinPeptideMaps_ — Public API method carries a trailing-underscore (private) name and must be called by external users
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h · _osw-mrm-feature_

```cpp
void prepareProteinPeptideMaps_(const OpenSwath::LightTargetedExperiment& transition_exp)
```
- **Expectation:** In OpenMS, a trailing underscore denotes a non-public/internal member (e.g. updateMembers_, computeQuality_). A method documented as 'Calling this method _is_ required before calling scorePeakgroups' is part of the supported public workflow and should not look private.
- **Actual:** prepareProteinPeptideMaps_ sits in the public section with a trailing underscore, and its own doc says callers MUST invoke it before scorePeakgroups. The naming convention thus actively contradicts its public, mandatory-to-call status, surprising users who skip underscore-suffixed members as internal.
- **Evidence:** /** @brief Prepares the internal mappings of peptides and proteins.\n     * Calling this method _is_ required before calling scorePeakgroups.\n     */\n    void prepareProteinPeptideMaps_(const OpenSwath::LightTargetedExperiment& transition_exp);
- **Fix:** Add a public, non-underscored alias (e.g. prepareProteinPeptideMaps) that forwards to the existing one, and document the underscored form as deprecated. Purely additive, ABI-safe. Renaming outright would break callers.
- **Verifier correction:** The naming inconsistency is real and the method is genuinely public/mandatory (called by OpenSwathWorkflow.cpp:1266 and exposed in pyOpenMS bind_misc.cpp:5173). However, the consequence of skipping it is NOT silent wrong results: scorePeakgroups throws Exception::IllegalArgument naming the exact method to call (MRMFeatureFinderScoring.cpp:533-537). Severity should therefore be medium (invites confusion, but loud and recoverable), not high. The recommended additive alias is ABI-safe (abi_impact none); an outright rename would break the 3 in-tree callers plus the pyOpenMS binding.

### [ANSW-43] MRMFeatureFilter::checkMetaValue — checkMetaValue returns true ("in range") when the metaValue is absent
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h · _osw-mrm-feature_

```cpp
bool checkMetaValue(const Feature& component, const std::string& meta_value_key, const double& meta_value_l, const double& meta_value_u, bool& key_exists) const
```
- **Expectation:** Per the header (`@return True if the metaValue is within the bounds, and False otherwise`), a caller expects false when the value cannot be checked. A missing required QC value would intuitively be a failed check.
- **Actual:** When the key is absent, the function sets key_exists=false but leaves the local `check = true` and returns true — i.e. a feature missing the QC metaValue silently PASSES the bound check. The only way to detect this is to also inspect the separate `key_exists` out-param, which the @return doc does not mention.
- **Evidence:** Impl src/openms/source/ANALYSIS/OPENSWATH/MRMFeatureFilter.cpp:1006-1017 — `bool check = true; if (component.metaValueExists(...)) { ... check = checkRange(...); } else { key_exists = false; ... } return check;`
- **Fix:** Document in @return that an absent key returns true and that callers must consult key_exists; or (behavior change, flag honestly) return false when key_exists is false. Doc fix is source/ABI-safe; behavior change is breaking for current callers that rely on pass-through.
- **Verifier correction:** checkMetaValue returns true (pass/"in range") when the metaValue key is absent: `check` is initialized true and the else-branch leaves it unchanged, only setting the separate key_exists=false out-param. The @return doc ("True if within bounds, False otherwise", MRMFeatureFilter.h:196) does not document the absent-key case, so callers relying solely on the return value silently treat features missing a required QC metaValue as passing. Note: OpenMS's own callers mitigate this by checking key_exists to exclude missing values from the QC-score denominator (not to fail the feature), so the practical hazard is for external/future callers — hence medium, not high. Recommend documenting in @return that an absent key returns true and that callers MUST inspect key_exists (source/ABI-safe); changing the return to false on absence would be a breaking behavior change for current pass-through callers.

### [ANSW-53] OpenSwath_Scores::calculate_lda_prescore / calculate_swath_lda_prescore / calculate_lda_single_transition — LDA prescore member functions ignore *this and score only the argument
`medium` · `misleading-name` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h · _osw-scoring-core_

```cpp
double calculate_lda_prescore(const OpenSwath_Scores& scores) const;
```
- **Expectation:** A non-static const member like myScores.calculate_lda_prescore(...) would be expected to compute the LDA score of the object it is called on (this->library_corr, this->norm_rt_score, ...). The functions take a parameter named 'scores' that a caller could reasonably read as auxiliary/extra context.
- **Actual:** The implementation reads ONLY the passed-in 'scores' argument (scores.library_corr, scores.norm_rt_score, ...) and never touches any member of *this. So `a.calculate_lda_prescore(b)` silently scores b and completely ignores a. The methods are effectively static but declared as const instance methods, inviting the call `feature_scores.calculate_lda_prescore()`-style usage where the obvious operand (the receiver) is dropped.
- **Evidence:** OpenSwathScores.cpp:62-68: `return scores.library_corr * -0.34664267 + scores.library_norm_manhattan * 2.98700722 + ...` — every term is `scores.X`, none is `this->X` or a bare member; the function body uses the parameter exclusively. Declared `double calculate_lda_prescore(const OpenSwath_Scores& scores) const;`.
- **Fix:** Make these `static` (additive, source-compatible: existing `obj.calculate_lda_prescore(x)` calls still compile) so the receiver-independence is explicit at the call site, and/or rename the parameter to make clear it is the sole operand. Ideal fix: `static double calculate_lda_prescore(const OpenSwath_Scores& scores);`.
- **Verifier correction:** The three functions are genuinely receiver-independent (they read only the `scores` parameter, never *this), confirming the misleading-name defect. However: (1) severity is medium, not high — all current callers pass the receiver as the argument (scores.calculate_lda_prescore(scores)), so no results are silently wrong today; the hazard is latent misuse of the redundant API. (2) The recommended `static` fix is ABI-breaking (symbol mangling changes, implicit this removed), though source-compatible for the existing obj.f(x) call form. A safer fix is to drop the parameter and use the members (so the receiver becomes meaningful) or convert to a documented free/static function plus rename the parameter to make it the sole operand.

### [ANSW-54] MRMScoring::calcXcorrCoelutionScore / calcXcorrShapeScore / calcMIScore (and siblings) — Score functions return 0 (a 'perfect coelution' value) on degenerate matrices in release builds
`medium` · `silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/OPENSWATH/MRMScoring.cpp · _osw-scoring-core_

```cpp
double calcXcorrCoelutionScore();  // guarded by OPENSWATH_PRECONDITION(rows>1, ...)
```
- **Expectation:** The header documents `calcXcorrCoelutionScore` as 'a distance where zero indicates perfect coelution' and the implementation has OPENSWATH_PRECONDITION(...rows>1...), so a caller assumes too-small input is rejected.
- **Actual:** OPENSWATH_PRECONDITION is `assert(...)` (Macros.h), which is a no-op under NDEBUG (release). With a 0x0 or 1x1 cross-correlation matrix (e.g. a single-transition peak group), `meanStdevUpperTriangle` feeds 0 or 1 samples to `mean_and_stddev`, whose `sample_stddev()` returns 0 for count<=1 and whose `mean()` is 0 for count 0. The function then silently returns 0.0 — i.e. the value the docs define as 'perfect coelution' — on essentially no data, instead of throwing or signalling.
- **Evidence:** MRMScoring.cpp:392-396 `calcXcorrCoelutionScore(){ OPENSWATH_PRECONDITION(xcorr_matrix_max_peak_.rows() > 1, ...); return meanStdevUpperTriangle(...); }`; Macros.h:19-20 `#define OPENSWATH_PRECONDITION(condition, message) assert((condition)&&(message));`; StatsHelpers.h:166-176 `sample_variance(){ return (c_>1u)?(q_/(c_-1)):0; }` and mean starts at 0.0.
- **Fix:** Either document that callers must guarantee >=2 transitions, or replace the assert-only guard with a real runtime throw (Exception) for these public-facing score entry points, mirroring what calcSNScore already half-does. At minimum return a clearly-degenerate sentinel (NaN) rather than the 'perfect' value 0.
- **Verifier correction:** The assert-only guard (OPENSWATH_PRECONDITION == assert, no-op under NDEBUG) on the public, documented score entry points (calcXcorrCoelutionScore/Shape/MI and siblings) does cause them to silently return 0.0 — the documented "perfect coelution" value — on a 0x0 or 1x1 matrix in release builds, instead of rejecting too-small input. This is a real API footgun for direct callers of MRMScoring who trust the precondition. However, it does NOT silently corrupt OpenMS's own results: the in-tree callers guard the degenerate single-transition / <2-series case explicitly (OpenSwathScores.cpp:71-77 uses a separate LDA model that omits coelution/shape; IonMobilityScoring.cpp:613-618 short-circuits series<2 and itself assigns coelution=0 / shape=NaN; OpenSwathScoring.cpp:519 guards precursor_ids.size()>1). Severity is therefore medium (invites misuse by external API consumers; the shipped pipeline is protected), not high. Recommendation stands: document the >=2-series requirement on these public methods, and/or return NaN rather than the "perfect" 0 on degenerate input so the no-data case is loud rather than indistinguishable from a perfect score.

### [ANSW-60] DIAScoring::dia_ms1_massdiff_score — On no-signal the ppm_score out-param is set to a plausible-looking finite ppm value, not a sentinel
`medium` · `return-value` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h · _osw-scoring-core_

```cpp
bool dia_ms1_massdiff_score(double precursor_mz, const SpectrumSequence& spectrum, const RangeMobility& im_range, double& ppm_score) const
```
- **Expectation:** A failing mass-error measurement should leave the out-param at an obviously-invalid sentinel (the sibling dia_massdiff_score uses -1 for missing transitions), so a caller that ignores the bool return value does not mistake the failure value for a real, good ppm error.
- **Actual:** On failure the function returns false but writes ppm_score = getPPMAbs(precursor_mz + window_span, precursor_mz) — a finite positive ppm roughly equal to the extraction window width, which looks like a legitimate (merely large) measurement. A caller who reads ppm_score without checking the bool silently treats a no-signal case as a measured mass error. (This IS documented in the header, but the value is easy to consume by accident given the parallel dia_massdiff_score uses -1.)
- **Evidence:** DIAScoring.cpp:299-303 `if (!signalFound){ ppm_score = Math::getPPMAbs(precursor_mz + mz_range.getSpan(), precursor_mz); return false; }`; contrast dia_massdiff_score.cpp:268 `diff_ppm.push_back(-1);` sentinel. Header DIAScoring.h:144-145 documents the worst-case-bound behavior.
- **Fix:** Keep behavior but ensure all callers branch on the bool; consider [[nodiscard]] on the return to force callers to handle the false case (additive, source-compatible). Alternatively document more loudly that ppm_score on failure is a worst-case bound, not a measurement.
- **Verifier correction:** The claim is accurate and the evidence checks out, with one important nuance the claim already concedes: the behavior is explicitly and prominently documented at the declaration (DIAScoring.h:144-145, including the literal phrase 'not -1'), and it is an intentional worst-case-bound design rather than an oversight. The genuinely surprising/actionable part is not 'undocumented sentinel confusion' but that (a) it is inconsistent with the sibling dia_massdiff_score's -1 convention, and (b) both production callers in OpenSwathScoring.cpp (lines 338 and 461) ignore the bool and consume the failure value as a real score that propagates to feature metadata and output writers. Severity is medium (recoverable, loud, intentional), not high. The recommendation to add [[nodiscard]] (source-compatible) to force callers to handle the false case is sound and would directly fix the two existing offending call sites.

### [ANSW-7] MasstraceCorrelator::matchMassTraces_ — Boolean flag 'padEnds' is declared as 'double' with default 'true'
`medium` · `param-order-or-bool` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h · _osw-scoring-extra_

```cpp
void matchMassTraces_(const MasstracePointsType& hull_points1, const MasstracePointsType& hull_points2, std::vector<double>& vec1, std::vector<double>& vec2, double mindiff, double padEnds = true)
```
- **Expectation:** A flag named 'padEnds' that is documented '@param[in] padEnds Whether to pad ends with zeros' and defaulted to 'true' reads as a bool; a caller expects padEnds(true/false) to toggle behavior cleanly.
- **Actual:** The parameter type is 'double', and the default 'true' is implicitly converted to 1.0. The body uses it only as a truthiness test (if (!padEnds)), so any nonzero double (e.g. an accidentally-swapped numeric argument) silently means 'pad', and 0.0 means 'don't pad'. The double type also invites confusion with the adjacent 'double mindiff' parameter directly before it.
- **Evidence:** Header: 'double padEnds = true'. cpp: 'if (!padEnds) { ... }'. Doc comment: '@param[in] padEnds Whether to pad ends with zeros'.
- **Fix:** Change the type to 'bool padEnds = true'. This is a protected method, so the blast radius is internal; still ABI-breaking for the mangled signature. If strict ABI must be preserved on this protected symbol, at minimum fix the type in a coordinated change. Prefer the bool fix.
- **Verifier correction:** Boolean flag 'padEnds' is declared as 'double padEnds = true' (implicitly 1.0) and used only as a truthiness test (`if (!padEnds)`), while documented as "Whether to pad ends with zeros". It sits directly after `double mindiff`, creating two adjacent doubles where one is secretly a bool — a swapped argument would compile silently. The method is protected and its only internal caller omits the argument (always defaulting to pad), so this is a latent type-safety hazard, not an active wrong-result bug; severity medium. Fix: change to `bool padEnds = true` (ABI-breaking on the mangled symbol, but source-compatible for true/false callers).

### [ANSW-8] MasstraceCorrelator::scoreHullpoints — 'max_lag' parameter is silently ignored
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h · _osw-scoring-extra_

```cpp
void scoreHullpoints(..., const int max_lag, ...)
```
- **Expectation:** A parameter named 'max_lag' on a cross-correlation scoring routine reads as bounding the lag search / output; a caller would expect passing a different max_lag to change the result.
- **Actual:** The implementation completely ignores the argument (it is commented out in the parameter list). The cross-correlation is always computed with a hard-coded window. The lag bounding (lag >= -max_lag && lag <= max_lag) happens only in the caller createPseudoSpectra, not inside scoreHullpoints. The header even admits '@param[in] max_lag Currently unused' but the public signature still advertises the parameter.
- **Evidence:** cpp signature: 'const double min_corr, const int /* max_lag */, const double mindiff'. Header doc: '@param[in] max_lag Currently unused'.
- **Fix:** Either honor max_lag inside the function (bound the xcorr lag search), or remove the parameter. Since removal is ABI-breaking for a public method, the low-risk fix is to keep documenting it as unused but consider deprecating; ideal fix is to actually use it or drop it.
- **Verifier correction:** scoreHullpoints does ignore its max_lag argument (commented out in the .cpp parameter list), but this is explicitly documented at the declaration as "@param[in] max_lag Currently unused," and the lag bound is correctly enforced by the sole in-tree caller createPseudoSpectra. So it is a documented dead public parameter rather than a hidden silent failure producing wrong results; severity is medium, not high. ABI is unaffected as-is (removing the parameter would be the breaking change, not the current state).

### [ANSW-9] SpectrumAddition::addUpSpectra — Single-element input is returned as an alias to the caller's spectrum, not a fresh summed spectrum
`medium` · `ownership-lifetime` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h · _osw-scoring-extra_

```cpp
static OpenSwath::SpectrumPtr addUpSpectra(const SpectrumSequence& all_spectra, double sampling_rate, bool filter_zeros)
```
- **Expectation:** A static 'addUpSpectra' that 'adds up a list of Spectra by resampling them and then addition of intensities' should return a newly-constructed, independently-owned result spectrum that the caller can mutate freely.
- **Actual:** When all_spectra.size() == 1 the function returns all_spectra[0] directly. Since SpectrumPtr is std::shared_ptr<Spectrum>, the returned pointer aliases the caller's input spectrum (shared ownership of the same object), unlike the multi-spectrum path which builds a new resampled container and the filter_zeros path which builds a fresh filtered spectrum. A caller who mutates the 'result' (e.g. sorts or filters it) will silently mutate their input, and no resampling to sampling_rate / no zero-filtering is applied in this branch.
- **Evidence:** cpp: 'if (all_spectra.size() == 1) return all_spectra[0];' with 'typedef std::shared_ptr<Spectrum> SpectrumPtr;' and 'using SpectrumSequence = std::vector<OpenSwath::SpectrumPtr>;'.
- **Fix:** For the size()==1 case, deep-copy into a new Spectrum (and honor sampling_rate / filter_zeros) before returning, so the result is always independently owned and consistently processed. Source-only change; no ABI impact.
- **Verifier correction:** For all_spectra.size()==1, addUpSpectra returns the caller's own spectrum pointer (a shared_ptr alias), and that branch also silently skips both resampling-to-sampling_rate and zero-filtering — so even with filter_zeros=true the result may retain zero-intensity points and unresampled m/z, inconsistent with the function's documented contract and with every other branch (which build fresh, processed Spectrum objects). The shared-aliasing-mutates-the-caller hazard is real in principle but latent: the only in-tree SpectrumPtr-overload caller treats the result as read-only and already uses the same single-element aliasing idiom for concatenate, so it is best characterized as a recognized convention plus an undocumented processing-inconsistency rather than an active data-corruption bug. Fix: in the size==1 case, deep-copy into a fresh Spectrum and honor sampling_rate/filter_zeros for consistent, independently-owned output. Source-only; no ABI impact.

### [ANSW-10] ConfidenceScoring::scoreMap / getAssayRT_ — Missing-RT assays throw instead of being skipped; the -1.0 'missing value' guard in scoreMap is dead
`medium` · `dead-code-unhandled-edge-case` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/ConfidenceScoring.h · _osw-scoring-extra_

```cpp
void scoreMap(FeatureMap& features); double getAssayRT_(const TargetedExperiment::Peptide& assay)
```
- **Expectation:** scoreMap()'s loop 'if (current_rt == -1.0) continue; // indicates a missing value' implies getAssayRT_ returns -1.0 for an assay without an RT, so library assays lacking RT are silently skipped when computing the RT range.
- **Actual:** getAssayRT_ only has an OPENMS_PRECONDITION (debug-only no-op in release) and then calls TargetedExperiment::Peptide::getRetentionTime(), which THROWS Exception::IllegalArgument ('No retention time information available') when no RT is set. It never returns -1.0. So the '== -1.0' guard is dead code, and a single RT-less assay makes the public scoreMap() throw mid-loop rather than skip it.
- **Evidence:** ConfidenceScoring.cpp getAssayRT_: 'OPENMS_PRECONDITION(assay.hasRetentionTime(), ...) return assay.getRetentionTime();'. TargetedExperimentHelper.h getRetentionTime(): 'if (!hasRetentionTime()) { throw Exception::IllegalArgument(...,"No retention time information available"); }'. scoreMap(): 'double current_rt = getAssayRT_(*it); if (current_rt == -1.0) continue; // indicates a missing value'.
- **Fix:** Make getAssayRT_ actually return the documented missing sentinel: 'if (!assay.hasRetentionTime()) return -1.0;' before reading the RT, matching the guard in scoreMap(). Source-only; both are non-public/inline so no external ABI impact.
- **Verifier correction:** File is at src/openms/include/OpenMS/ANALYSIS/OPENSWATH/ConfidenceScoring.h (not QUANTITATION). The core claim holds: getAssayRT_ never returns the -1.0 'missing value' sentinel that scoreMap()'s 'if (current_rt == -1.0) continue;' guard checks for — OPENMS_PRECONDITION is a release no-op and getRetentionTime() throws Exception::IllegalArgument when no RT is set. So the guard is dead code and one RT-less library assay makes public scoreMap() throw mid-loop instead of skipping it. The category is more accurately dead-code/unhandled-missing-value than 'silent-failure' (the failure is a loud throw); the silently-broken part is the never-firing guard. Severity medium (loud, recoverable). ABI impact none: the fix only changes the out-of-line body of protected getAssayRT_ in the .cpp; no signature/declaration changes. Recommended fix: add 'if (!assay.hasRetentionTime()) return -1.0;' in getAssayRT_ before reading the RT.

### [ANSW-76] TargetedSpectraExtractor::targetedMatching — targetedMatching silently overrides the user-configured 'top_matches_to_report' parameter (and is non-const/non-reentrant because of it)
`medium` · `hidden-side-effect` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h · _osw-transitions-helper_

```cpp
void targetedMatching(const std::vector<MSSpectrum>& spectra, const Comparator& cmp, FeatureMap& features)
```
- **Expectation:** A 'matching' method that takes (spectra, comparator, features) reads the object's configured parameters and is logically const, like every sibling scoring/selecting method in this class (all the other analysis methods are declared const).
- **Actual:** It is the only matching method declared non-const because it mutates the member 'top_matches_to_report_': it saves the value, forces it to 1, runs, then restores it (TargetedSpectraExtractor.cpp:779-803). The user-set 'top_matches_to_report' parameter is therefore silently ignored inside this call, and the temporary mutation makes the method non-reentrant / not thread-safe even though it appears to only fill 'features'.
- **Evidence:** const Size tmp = top_matches_to_report_;\n    top_matches_to_report_ = 1;\n    for (...) { matchSpectrum(...); ... }\n    top_matches_to_report_ = tmp;  // TargetedSpectraExtractor.cpp:779-803
- **Fix:** Refactor matchSpectrum to take an explicit 'max_matches' argument (defaulting to top_matches_to_report_) so targetedMatching can pass 1 locally without touching the member; then targetedMatching/untargetedMatching can be const. Additive, ABI-safe path: add a const overload of matchSpectrum(input, cmp, matches, n) and call it; keep the member untouched.

### [ANSW-78] OpenSwathHelper::selectSwathTransitions — selectSwathTransitions appends to its '@param[out]' output instead of clearing it first
`medium` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h · _osw-transitions-helper_

```cpp
static void selectSwathTransitions(const OpenMS::TargetedExperiment&, OpenMS::TargetedExperiment& selected_transitions, double, double, double); (and LightTargetedExperiment overload)
```
- **Expectation:** A function documented with '@param[out] selected_transitions Selected transitions for SWATH window' should fully populate (and therefore reset) the output, so calling it twice with the same object yields only the second window's transitions. Sibling readers/matchers in this module (annotateSpectra, matchSpectrum, TransitionParquetFile reader) explicitly clear their out-params.
- **Actual:** Both overloads only append: the TargetedExperiment version calls addTransition() in a loop (and overwrites peptides/proteins via setPeptides/setProteins, creating an inconsistent half-overwrite/half-append result); the Light version push_back()s into transitions/compounds/proteins. Pre-existing contents are silently retained and accumulate across calls.
- **Evidence:** transition_exp_used.addTransition(tr);  // OpenSwathHelper.cpp:31 (no clear of transitions)\n...\ntransition_exp_used.transitions.push_back(tr);  // OpenSwathHelper.cpp:127 (Light overload, no clear)
- **Fix:** Clear the output experiment at function entry (transitions/peptides/proteins/compounds) to match the '@param[out]' contract and the rest of the module. If existing callers rely on accumulation this is a behavior change; the safe documented fix is to state '@param[in,out] ... appended to' in the header. Prefer clearing.
- **Verifier correction:** Severity is medium rather than high. Every current in-tree caller (OpenSwathWorkflow.cpp:966 and :579 via a fresh per-context member; DIAChromHandler.cpp:69) declares a brand-new output object immediately before each call, so no existing call site reuses the object across calls and no live bug is triggered today. The hazard is for reasonable future/external (e.g. pyOpenMS) callers who, trusting the '@param[out]' contract, reuse one object across SWATH windows in a loop and silently get an accumulated (and, for the TargetedExperiment overload, internally inconsistent) experiment with no error. That is medium: invites silently-wrong results but is not currently producing them in production paths. Recommended fix (clear transitions/peptides/proteins/compounds at function entry) is source-compatible — it changes runtime behavior, not the signature or exported symbol; abi_impact is therefore source-compatible, not breaking.

### [ANSW-79] OpenSwathHelper::sampleExperiment — sampleExperiment silently drops all decoys and throws if fewer than 3 candidates remain, neither of which is documented
`medium` · `documentation-contract-gap` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h · _osw-transitions-helper_

```cpp
static OpenSwath::LightTargetedExperiment sampleExperiment(const OpenSwath::LightTargetedExperiment&, Size bins, Size peptides_per_bin, unsigned int seed = 0, bool sort_by_intensity = false, double top_fraction = 1.0, const std::unordered_set<std::string>& = {})
```
- **Expectation:** A generic 'sample a subset of peptides uniformly across the RT range' helper should sample from the provided experiment and return a (possibly small/empty) subset; the doc lists no filtering and no throwing precondition.
- **Actual:** It silently excludes every transition/compound flagged as decoy before sampling (good_ids only collects !getDecoy()), and it throws Exception::IllegalArgument when the surviving candidate count is < 3. Neither the decoy filtering nor the throw-on-small-input is mentioned in the header documentation, so a caller passing a small or decoy-heavy library gets a surprising exception instead of a subset.
- **Evidence:** if (!tr.getDecoy()) good_ids.insert(tr.getPeptideRef());  // OpenSwathHelper.cpp:191\nif (all_candidates.size() < 3) { throw Exception::IllegalArgument(...); }  // OpenSwathHelper.cpp:245-250
- **Fix:** Document both behaviors in the header (@note decoys are excluded; @throws IllegalArgument if fewer than 3 non-decoy candidates). Doc-only change, ABI-neutral. Optionally expose decoy inclusion as a parameter.
- **Verifier correction:** The function does silently exclude decoy-flagged compounds (OpenSwathHelper.cpp:191,198) and throws Exception::IllegalArgument when fewer than 3 non-decoy candidates remain (lines 245-250), and neither is documented in the header (OpenSwathHelper.h:184-204). However, the original category "silent-failure" overstates it: the <3 case is a loud, descriptive exception rather than a silent failure; only the decoy filtering is silent. Severity is medium (loud/recoverable, clear message, decoy exclusion is contextually reasonable for the iRT-calibration purpose the doc references) rather than high. Fix is doc-only and ABI-neutral.

### [ANSW-81] TransitionTSVFile::convertTSVToTargetedExperiment — convertTSVToTargetedExperiment (TargetedExperiment overload) appends transitions to a non-empty output but overwrites peptides/proteins/compounds
`medium` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h · _osw-transitions-helper_

```cpp
void convertTSVToTargetedExperiment(const char* filename, FileTypes::Type filetype, OpenMS::TargetedExperiment& targeted_exp)
```
- **Expectation:** A 'convert file into targeted_exp' reader with '@param[out] targeted_exp' should yield exactly the file's contents regardless of the object's prior state, and behave like the sibling TransitionParquetFile reader which the docs say is '**reset** to empty' before filling.
- **Actual:** TSVToTargetedExperiment_ calls exp.addTransition() in a loop (appending to any pre-existing transitions) while calling setCompounds/setPeptides/setProteins at the end (overwriting those). Reusing a TargetedExperiment across two reads thus yields duplicated/orphaned transitions plus only the last file's peptides/proteins. The Light overload streams and likewise does not clear. The inconsistency with the Parquet reader's explicit reset is itself a cross-class surprise.
- **Evidence:** exp.addTransition(rm_trans);  // TransitionTSVFile.cpp:597\nexp.setCompounds(compounds); exp.setPeptides(peptides); exp.setProteins(proteins);  // TransitionTSVFile.cpp:642-644 (no clear of transitions)
- **Fix:** Clear targeted_exp at the start of convertTSVToTargetedExperiment (both overloads) to match the Parquet reader and the '@param[out]' contract. Behavioral change for any caller relying on accumulation; document explicitly either way. ABI-neutral.
- **Verifier correction:** Both convertTSVToTargetedExperiment overloads fail to reset the output object, contrary to the '@param[out]' contract and the sibling TransitionParquetFile reader which is explicitly 'cleared before being filled'. Heavy (TargetedExperiment) overload: appends transitions (addTransition=push_back, line 597) but overwrites peptides/proteins/compounds (setPeptides/setProteins/setCompounds, lines 642-644) -> on reuse with a non-empty object, duplicated/orphaned transitions plus only the last file's peptides/proteins. NOTE: a debug-only POSTCONDITION (line 646; active only under OPENMS_ASSERTIONS) makes this LOUD (throws) in assertion builds and SILENT in release builds. Light (LightTargetedExperiment) overload: never resets and push_backs transitions+compounds+proteins -> pure accumulation of everything on reuse, with no postcondition guard at all. Recommendation stands: reset the output at the start of both overloads (e.g. assign a default-constructed object, mirroring TransitionParquetFile:209) and document the chosen semantics. ABI-neutral (function-body change, no signature change).

### [ANSW-4] SwathQC::setNrMS1Spectra — setNrMS1Spectra(0) means 'inspect ALL spectra', not 'expect zero spectra'
`medium` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SwathQC.h · _osw-workflow_

```cpp
void setNrMS1Spectra(size_t nr)
```
- **Expectation:** A setter named setNrMS1Spectra(size_t) records the count of MS1 spectra; passing 0 would naturally mean 'there are zero MS1 spectra' (so nothing is sampled).
- **Actual:** Passing 0 switches the sampler into 'sample everything' mode: isSubsampledSpectrum_ returns true for every index when total_spec_count==0. So setNrMS1Spectra(0) maximizes work instead of disabling it. This 0-as-sentinel inversion is only discoverable in the doc.
- **Evidence:** Header SwathQC.h:103 'If nr is set to 0, all spectra passed into getSpectraProcessingFunc() will be inspected'. Impl SwathQC.cpp:172 'if (total_spec_count == 0) return true;' and SwathQC.cpp:117 the static helper sets `sq.setNrMS1Spectra(0); // leave at 0, such that all incoming spectra are sampled`.
- **Fix:** Behavior is load-bearing (used by the static getChargeDistribution path), so keep it but make the sentinel explicit: document the 0 meaning in the @brief and ideally add a named helper like sampleAllSpectra() that forwards to setNrMS1Spectra(0). Source-compatible/additive.
- **Verifier correction:** setNrMS1Spectra(0) does not mean "there are zero MS1 spectra"; it is a sentinel meaning "total count unknown -> inspect EVERY spectrum fed to getSpectraProcessingFunc()" (isSubsampledSpectrum_ returns true for all indices when total_spec_count==0, SwathQC.cpp:172). This silently maximizes work (peak-pick + deisotope per MS1) and populates a charge distribution rather than disabling sampling. The behavior is load-bearing: the static getChargeDistribution() path deliberately calls setNrMS1Spectra(0) (line 117) because it does its own external subsampling. The surprise is documented in the detailed doxygen body (line 103) but NOT in the @brief (line 98), so it is missed by anyone reading the signature/brief. Severity is medium (silent semantic inversion, but loud by CPU cost and recoverable) rather than high.

### [ANSW-39] MetaboTargetedTargetDecoy::generateMissingDecoysByMassShift — mass_to_add documented as 'maximum number of transitions per assay' but is a Da mass shift
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TARGETED/MetaboTargetedTargetDecoy.h · _targeted-chrom-handlers_

```cpp
static void generateMissingDecoysByMassShift(TargetedExperiment& t_exp, std::vector<MetaboTargetDecoyMassMapping>& mappings, const double& mass_to_add)
```
- **Expectation:** A caller reads '@param[in] mass_to_add the maximum number of transitions required per assay' and would pass an integer transition count.
- **Actual:** mass_to_add is a m/z mass offset (Da) added directly to product m/z values: 'return d + mass_to_add' and 'setProductMZ(it.decoy_product_masses[i])' (with masses = target + mass_to_add). The doc text is copy-pasted/wrong; the sibling resolveOverlappingTargetDecoyMassesByDecoyMassShift correctly documents 'mass_to_add (e.g. CH2)'.
- **Evidence:** MetaboTargetedTargetDecoy.cpp: '[mass_to_add](double d) -> double { return d + mass_to_add; }' applied to target_product_masses; header doc line: '@param[in] mass_to_add the maximum number of transitions required per assay'.
- **Fix:** Fix the Doxygen for mass_to_add to 'mass (Da) added to target product m/z to form the shifted decoy fragment (e.g. CH2 = 14.0157)'. Documentation-only; no ABI impact.
- **Verifier correction:** mass_to_add is a Da mass offset added directly to each target product m/z to form the shifted decoy fragment (e.g. CH2 = 14.0157 Da), not a transition count. Fix the Doxygen at MetaboTargetedTargetDecoy.h:65 to read e.g. "mass (Da) added to target product m/z to form the shifted decoy fragment (e.g. CH2 = 14.0157)". Documentation-only; double parameter type means an erroneously-passed integer count would compile and silently corrupt decoy m/z, but no crash/data loss and the single in-tree caller is already correct.

### [ANSW-30] ReactionMonitoringTransition::getPrediction — getPrediction() guards against the wrong member (checks hasPrecursorCVTerms instead of hasPrediction); dereferences a null prediction_ in release
`medium` · `silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/MRM/ReactionMonitoringTransition.cpp · _targeted-experiment_

```cpp
const ReactionMonitoringTransition::Prediction & getPrediction() const
```
- **Expectation:** The header documents 'You first need to check whether the object is accessible using hasPrediction()' and the precondition should verify prediction_ is non-null before returning *prediction_.
- **Actual:** The body asserts the wrong condition: OPENMS_PRECONDITION(hasPrecursorCVTerms(), "ReactionMonitoringTransition has no Prediction object, check first with hasPrediction()"). It checks precursor_cv_terms_ != nullptr, NOT prediction_ != nullptr. A caller who correctly calls hasPrediction() (false) and skips the call is fine, but if the assert ever fires it fires on the wrong member, and in release builds (where OPENMS_PRECONDITION is a no-op) any caller hitting getPrediction() on a transition with prediction_==nullptr dereferences a null pointer and gets UB.
- **Evidence:** Line 320-324: `const ReactionMonitoringTransition::Prediction & ReactionMonitoringTransition::getPrediction() const { OPENMS_PRECONDITION(hasPrecursorCVTerms(), "ReactionMonitoringTransition has no Prediction object, check first with hasPrediction()") return *prediction_; }`
- **Fix:** Change the precondition predicate to hasPrediction() (one-token fix). ABI-neutral, source-compatible. Optionally throw a real Exception instead of a precondition so a missing Prediction signals in release builds too.
- **Verifier correction:** The documented correct usage (call hasPrediction(), skip getPrediction() when false) works in both debug and release. The defect is that the precondition guards the wrong member: a caller who misuses the API (calls getPrediction without checking) is supposed to hit a loud assert in debug, but because the assert tests precursor_cv_terms_ instead of prediction_, it can silently pass when precursor CV terms exist, letting *prediction_ dereference a nullptr (UB). It can also spuriously fire when precursor CV terms are absent. In release the no-op precondition means a correct guard would not help anyway, so the unique harm is the defeated debug-time safety net. Severity adjusted to medium (defeated guard / invites silent UB only on API misuse), not high (correct documented usage is unaffected). Fix: replace hasPrecursorCVTerms() with hasPrediction() at line 322.

### [ANSW-31] TargetedExperiment::getProteinByRef — getProteinByRef/getPeptideByRef/getCompoundByRef silently insert and dereference a null pointer for an unknown ref in release builds
`medium` · `silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/TARGETED/TargetedExperiment.cpp · _targeted-experiment_

```cpp
const Protein & getProteinByRef(const std::string & ref) const
```
- **Expectation:** A const 'getXByRef' lookup for a nonexistent reference should signal failure (throw / documented precondition), not corrupt the lookup table.
- **Actual:** The lookup uses non-const std::unordered_map::operator[] on the mutable reference map: `return *(protein_reference_map_[ref]);`. For a missing ref, operator[] value-initializes a new entry to a null `const Protein*`, then the code dereferences it (null-pointer UB). It also silently mutates the cache (inserts a bogus entry). The only guard is OPENMS_PRECONDITION, which is compiled out in release. Same pattern in getPeptideByRef (line 464) and getCompoundByRef (line 474).
- **Evidence:** Lines 400-408: `if (protein_reference_map_dirty_) { createProteinReferenceMap_(); } OPENMS_PRECONDITION(protein_reference_map_.contains(ref), "Could not find protein in map") return *(protein_reference_map_[ref]);`
- **Fix:** Use find() and throw Exception::ElementNotFound (or return after checking contains()) on a missing ref instead of operator[]; mark the cache map members truly const-lookup-safe. Additive/behavioral fix, ABI-neutral (no signature change).
- **Verifier correction:** The mechanism and code are exactly as claimed. Only the severity is adjusted from high to medium: the API provides hasProtein/hasPeptide/hasCompound and the assertion messages explicitly direct callers to check first, so a recoverable correct-usage path exists, and the misuse failure mode is a null-pointer dereference crash (plus silent cache insertion) rather than silently-wrong scientific results. Confirmed locations: getProteinByRef line 407, getPeptideByRef line 464, getCompoundByRef line 474; mutable maps header lines 300/304/308; Debug-only assertions src/openms/CMakeLists.txt lines 27-32.

### [ANSW-32] TargetedExperiment::hasProtein / hasPeptide / hasCompound / getProteinByRef — const 'has*'/'get*ByRef' lookups lazily rebuild and mutate internal reference maps (non-thread-safe surprise on a const object)
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h · _targeted-experiment_

```cpp
bool hasProtein(const std::string & ref) const
```
- **Expectation:** A const predicate named hasProtein() reads state and is safe to call concurrently on a shared const TargetedExperiment.
- **Actual:** On first call after any modification these const methods run createProteinReferenceMap_()/createPeptideReferenceMap_()/createCompoundReferenceMap_(), which write the mutable reference_map_ and clear the *_dirty_ flag. Two threads calling has*/get*ByRef on the same const object concurrently race on the mutable maps. The cost (full O(N) rebuild) is also hidden behind a cheap-looking const predicate.
- **Evidence:** TargetedExperiment.cpp lines 410-417: `bool TargetedExperiment::hasProtein(...) const { if (protein_reference_map_dirty_) { createProteinReferenceMap_(); } return protein_reference_map_.contains(ref); }`; createProteinReferenceMap_ writes `protein_reference_map_[...] = ...; protein_reference_map_dirty_ = false;` (lines 678-685). Members are declared `mutable` in the header (lines 300-310).
- **Fix:** Document the lazy-cache side effect on these const accessors; for thread safety guard the rebuild with a mutex or build the maps eagerly in the setters/adders. ABI-neutral if implemented internally.
- **Verifier correction:** Severity is medium, not high. The unsynchronized mutable-cache mutation in const methods is real UB under concurrent const-read, BUT exposure in current code is latent rather than demonstrated: the main OpenSWATH hot path (ChromatogramExtractor.h:417-425, extract_id_) reaches these via a NON-const TargetedExperiment& specialization, and no current call site was found sharing one const TargetedExperiment across threads for concurrent has*/get*ByRef. So it invites misuse and can crash/corrupt if hit, but is not a demonstrated silent-wrong-result in shipping code paths -> medium. On ABI: a documentation-only fix is ABI-neutral (none), and eager map-building in the setters/adders is source-compatible (no signature change); however the claim's alternative remedy of adding a std::mutex/std::atomic member to guard the rebuild WOULD change class layout and is ABI-breaking — graded source-compatible as the representative non-breaking remediation exists.

### [ANSW-35] ReactionMonitoringTransition::getPrecursorCVTermList / getProductMZ accessors family — getPrecursorCVTermList() dereferences a possibly-null precursor_cv_terms_ guarded only by a precondition that is a no-op in release builds
`medium` · `release-only-precondition-nullptr-deref` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h · _targeted-experiment_

```cpp
const CVTermList & getPrecursorCVTermList() const
```
- **Expectation:** Given the header note 'You first need to check whether they exist using hasPrecursorCVTerms()', a caller who forgets the check expects a thrown exception, not undefined behavior.
- **Actual:** The body is `OPENMS_PRECONDITION(hasPrecursorCVTerms(), ...) return *precursor_cv_terms_;`. precursor_cv_terms_ is a raw pointer defaulting to nullptr. In release builds OPENMS_PRECONDITION compiles away, so calling getPrecursorCVTermList() without first calling hasPrecursorCVTerms() dereferences nullptr (UB) rather than signalling. Same release-build hazard applies to getProductChargeState() -> product_.getChargeState() which is precondition-guarded on charge_set_ and returns garbage/0 silently.
- **Evidence:** ReactionMonitoringTransition.cpp lines 243-247: `const CVTermList & ReactionMonitoringTransition::getPrecursorCVTermList() const { OPENMS_PRECONDITION(hasPrecursorCVTerms(), ...) return *precursor_cv_terms_; }`. Header lines 97-108 document the hasPrecursorCVTerms() contract.
- **Fix:** For these pointer-backed getters, throw a real Exception (e.g. Exception::InvalidValue / Precondition) when the backing object is absent so release builds fail loudly instead of corrupting memory. Behavioral hardening, ABI-neutral.
- **Verifier correction:** Confirmed: getPrecursorCVTermList() (and getPrediction()) dereference a possibly-null raw pointer guarded only by OPENMS_PRECONDITION, which compiles to a no-op in release builds (Macros.h:94), so violating the documented hasPrecursorCVTerms() contract causes a nullptr dereference (UB/crash) in production instead of a thrown exception. Reject the secondary getProductChargeState() example: that accessor is NOT precondition-guarded in ReactionMonitoringTransition; it returns product_.getChargeState(), whose underlying charge_ is value-initialized to 0, yielding a well-defined sentinel 0 (not garbage, not UB) when unset — and 0 is the intended 'unknown charge' value callers rely on. Severity downgraded to medium: the failure requires caller misuse (skipping the documented check), the contract is documented in the header, and the failure mode is a hard crash (loud) rather than silently corrupted results. Recommendation stands for the pointer-backed getters only: throw a real Exception (e.g. Exception::Precondition/InvalidValue) unconditionally in the body so release builds fail loudly; this is ABI-neutral (.cpp body change, no signature/layout change).

### [ANSW-36] TargetedExperiment::operator+= / operator+ — operator+= concatenates proteins/peptides/transitions without de-duplicating refs, producing a merged experiment that can fail containsInvalidReferences() (duplicate ids)
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h · _targeted-experiment_

```cpp
TargetedExperiment& operator+=(const TargetedExperiment & rhs)
```
- **Expectation:** The header says 'Proteins, peptides and transitions are merged'; merging two valid libraries is expected to yield a valid library, or to reconcile shared ids.
- **Actual:** operator+= does a raw vector concatenation of every collection (cvs_, proteins_, peptides_, transitions_, ...). It performs no duplicate-id check; the source even leaves `// todo: check for double entries`. If both operands share a protein/peptide/compound/transition id (common when merging overlapping assay libraries), the result silently contains duplicate ids, which then makes containsInvalidReferences() report the merged experiment as invalid and makes getXByRef ambiguous.
- **Evidence:** TargetedExperiment.cpp lines 153-185: a series of `proteins_.insert(proteins_.end(), rhs.proteins_.begin(), rhs.proteins_.end());` etc., ending with `// todo: check for double entries // transitions, peptides, proteins`. containsInvalidReferences() returns true on any duplicate id (lines 584-633).
- **Fix:** Document that operator+/+= performs naive concatenation and does NOT deduplicate (caller must ensure disjoint ids), or add an opt-in de-duplicating merge. Doc fix is ABI-neutral; behavioral merge change should be additive (new method).
- **Verifier correction:** operator+= concatenates all collections without de-duplicating ids (author's own `// todo: check for double entries` confirms it is unfinished), so merging two libraries that share an id produces an experiment that containsInvalidReferences() flags as invalid. The original claim that this is uniformly silent is overstated: TSV and PQP export (TransitionTSVFile.cpp:1780, TransitionPQPFile.cpp:1038) call containsInvalidReferences() and throw IllegalArgument, so those paths fail loudly. The genuinely silent damage is limited to (1) TraML export / in-memory experiments, which are never validated (e.g. FileMerger.cpp:392 stores merged TraML directly), and (2) getXByRef, whose reference map (createProteinReferenceMap_, lines 678-693) is a last-write-wins std::map so duplicate ids silently resolve to a single arbitrary entry while the underlying vector still holds both. Recommendation stands and is ABI-neutral: document the disjoint-id precondition on operator+/+= (the current "see operator+= for details" doc is misleading because operator+= documents no details), and/or add an opt-in de-duplicating merge as a new method.

### [ANSW-48] TargetedExperimentHelper::PeptideCompound::getChargeState / TraMLProduct::getChargeState — getChargeState() silently returns 0 when charge was never set (release builds), so '0' is ambiguous between 'unset' and a genuine value
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h · _targeted-experiment_

```cpp
int getChargeState() const
```
- **Expectation:** Either getChargeState() always returns a meaningful charge, or it fails when no charge was set -- as the header strongly implies with hasCharge() and the precondition message 'Cannot return charge which was never set'.
- **Actual:** getChargeState() is `OPENMS_PRECONDITION(charge_set_, "...") return charge_;` with `charge_{0}` default. The precondition is a no-op in release, so a caller that forgets hasCharge() gets 0 with no signal. Since 0 is not a valid peptide charge, callers may treat 0 as a sentinel -- but ReactionMonitoringTransition::getProductChargeState() forwards to product_.getChargeState(), and a real product charge could legitimately be reported the same way after default construction, making the unset case indistinguishable from a populated-but-zero state.
- **Evidence:** `int getChargeState() const { OPENMS_PRECONDITION(charge_set_, "Cannot return charge which was never set") return charge_; }` with `int charge_{0};`. Same pattern in TraMLProduct.
- **Fix:** Document in the header that the return value is meaningless unless hasCharge() is true (it partly does for the analyte). For a stronger fix, throw in release when charge_set_ is false; that is a behavior change but ABI-safe (inline body). Keep hasCharge() as the documented guard.
- **Verifier correction:** The claim is correct and, if anything, understated. The ambiguity is not merely hypothetical: a real in-tree caller, OpenSwathScoring.cpp:446, reads transition.getProductChargeState() WITHOUT the isProductChargeStateSet() guard and tests `!= 0`, depending on the charge_{0} default leaking through the no-op precondition to mean 'unset'. Thus the precondition message 'Cannot return charge which was never set' is contradicted by release behavior the codebase itself relies on. Severity is medium (not high): the documented guards hasCharge()/isProductChargeStateSet() exist and are honored by most callers, and 0 is non-physical so the sentinel convention is internally workable -- the failure requires ignoring the documented guard. Recommended fix (header doc clarifying the value is meaningless unless hasCharge()) is ABI/source neutral; the stronger throw-in-release variant is source-compatible (inline body, behavior change only).

### [ANSW-50] TargetedExperiment::hasProtein / hasPeptide / hasCompound / getProteinByRef (const) — const lookup methods lazily build and mutate internal reference maps (rebuilt wholesale on any add*, never invalidated on direct vector mutation)
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/source/ANALYSIS/TARGETED/TargetedExperiment.cpp · _targeted-experiment_

```cpp
bool hasProtein(const std::string & ref) const
```
- **Expectation:** A const has*/get*ByRef call is a pure read; calling it concurrently from multiple threads on a const TargetedExperiment is safe.
- **Actual:** hasProtein()/getProteinByRef() (and peptide/compound variants) check a `mutable bool *_dirty_` flag and, if dirty, call createProteinReferenceMap_() which populates the `mutable` map. Thus a const method performs a non-thread-safe lazy mutation (data race if two threads call has*() concurrently on a shared const object). Additionally the dirty flag is only set by add*/set* methods -- the caches store raw `const Protein*` into proteins_ elements, so any reallocation of the underlying vector (e.g. via a non-tracked path) would leave dangling pointers; the design assumes callers never mutate the vectors except through add*/set*.
- **Evidence:** `mutable bool protein_reference_map_dirty_;` and `if (protein_reference_map_dirty_) { createProteinReferenceMap_(); }` inside const hasProtein/getProteinByRef; createProteinReferenceMap_() stores `&getProteins()[i]`.
- **Fix:** Document that const lookups are not thread-safe due to lazy cache init, or build the map eagerly / guard with a mutex. ABI-safe (internal). At minimum add a header note next to has*/get*ByRef.

### [ANSW-51] TargetedExperiment::operator+= / operator+ — operator+= concatenates all lists without de-duplicating ids, so merging can silently create duplicate (invalid) protein/peptide/transition references
`medium` · `silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/TARGETED/TargetedExperiment.cpp · _targeted-experiment_

```cpp
TargetedExperiment& operator+=(const TargetedExperiment & rhs)
```
- **Expectation:** Merging two assay libraries with '+' yields a consistent library; at minimum, identical ids are not duplicated.
- **Actual:** operator+= does a plain vector insert/append of cvs, proteins, peptides, transitions, etc. with an explicit `// todo: check for double entries`. So combining two libraries that share a protein/peptide/transition id produces duplicate ids, which containsInvalidReferences() then reports as invalid. The header for operator+ only says 'Proteins, peptides and transitions are merged' with no warning that duplicates are not collapsed, so a caller reasonably expects set-union semantics but gets list concatenation.
- **Evidence:** Impl: `proteins_.insert(proteins_.end(), rhs.proteins_.begin(), rhs.proteins_.end()); ... // todo: check for double entries // transitions, peptides, proteins`.
- **Fix:** Document that '+' is a concatenation (multiset union), not a de-duplicating merge, and that the result may contain duplicate ids; point callers to containsInvalidReferences(). ABI-safe doc fix. A de-duplicating merge would be a behavior change and should be opt-in.
- **Verifier correction:** operator+/operator+= perform plain list concatenation (multiset union) of all member vectors with no id de-duplication, as flagged by the in-source "// todo: check for double entries". Merging two libraries that share a protein/peptide/compound/transition id therefore produces duplicate ids that the class's own containsInvalidReferences() reports as invalid (returns true with OPENMS_LOG_ERROR). The header doc ("Proteins, peptides and transitions are merged" / "Add one targeted experiment to another") gives no warning and arguably implies set-union. Severity is medium (not high): the invalid state is detectable and recoverable via containsInvalidReferences() and surfaces as loud LOG_ERROR downstream, rather than silently producing wrong numbers or crashing. Fix is documentation-only (ABI-safe); a de-duplicating merge would be a behavior change and should be opt-in.

### [ANSW-52] ReactionMonitoringTransition::getProductMZ / setProductMZ — Product m/z (and precursor m/z) default to a silent sentinel (0.0) that getProductMZ()/getPrecursorMZ() returns with no 'isSet' indicator
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h · _targeted-experiment_

```cpp
double getProductMZ() const
```
- **Expectation:** getProductMZ()/getPrecursorMZ() either return a real m/z or signal that none was set.
- **Actual:** The header documents 'The default values for precursor m/z is 0.0 which indicates that it is uninitialized.' product_ defaults mz_{0} and precursor_mz_(0.0). getProductMZ()/getPrecursorMZ() just return the stored double, so an uninitialized transition reports m/z 0.0 indistinguishable from an (invalid) real 0.0 m/z. There is no hasProductMZ()/hasPrecursorMZ(). MRMMapping relies on exactly this: it treats |mz|<1e-5 as 'not set'. Callers without that knowledge can map or compute against a 0.0 m/z silently. (Contrast: IncludeExcludeTarget uses numeric_limits<double>::max() as its uninitialized sentinel -- inconsistent sentinel convention between the two sibling transition classes.)
- **Evidence:** RMT header: 'default values for precursor m/z is 0.0 which indicates that it is uninitialized'; IncludeExcludeTarget ctor: `precursor_mz_(std::numeric_limits<double>::max()), product_mz_(std::numeric_limits<double>::max())`.
- **Fix:** Document the 0.0 sentinel on getPrecursorMZ()/getProductMZ() and note the divergent sentinel (0.0 vs DBL_MAX) versus IncludeExcludeTarget. Optionally add additive hasPrecursorMZ()/hasProductMZ() helpers. ABI-safe.
- **Verifier correction:** getProductMZ()/getPrecursorMZ() return a stored double whose uninitialized value is the silent sentinel 0.0 (product m/z 0.0 is undocumented; precursor 0.0 is noted only at class level), with no hasProductMZ()/hasPrecursorMZ() helper despite the class providing isProductChargeStateSet() for charge. The sentinel also diverges from sibling class IncludeExcludeTarget (DBL_MAX). Recommended fix is additive/documentation only (source-compatible). Downgraded from high to medium: the MRMMapping 1e-5 check is on the chromatogram, not the transition, and transition matching is tolerance-based, so a default-0.0 transition will not silently mismap real chromatograms; the real exposure is undocumented computation against a 0.0 m/z by unaware callers, which is recoverable, not guaranteed-wrong.

### [ANSW-19] PeakPickerMobilogram::pickMobilogram — 3-arg pickMobilogram takes input Mobilogram by value while the 2-arg overload and the sibling PeakPickerChromatogram take const&
`low` · `inconsistent-convention` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/PeakPickerMobilogram.h · _osw-chrom_

```cpp
void pickMobilogram(Mobilogram mobilogram, Mobilogram& picked_mobilogram, Mobilogram& smoothed_mobilogram)
```
- **Expectation:** Two same-named overloads in one class (and the parallel PeakPickerChromatogram class) should accept the input the same way; a caller expects `const Mobilogram&` like the 2-arg overload, avoiding a silent deep copy of a potentially large mobilogram.
- **Actual:** The 2-arg overload is `pickMobilogram(const Mobilogram& mobilogram, ...)` but the 3-arg overload is `pickMobilogram(Mobilogram mobilogram, ...)` (by value). The parallel PeakPickerChromatogram declares BOTH overloads as `const MSChromatogram&`. The by-value parameter forces an extra copy on every call and reads as if the input might be consumed/mutated, when in fact it is only used locally.
- **Evidence:** PeakPickerMobilogram.h: `void pickMobilogram(const Mobilogram& mobilogram, Mobilogram& picked_mobilogram);` vs `void pickMobilogram(Mobilogram mobilogram, Mobilogram& picked_mobilogram, Mobilogram& smoothed_mobilogram);`. PeakPickerChromatogram.h: both overloads use `const MSChromatogram& chromatogram`. Impl PeakPickerMobilogram.cpp line 67 confirms by-value signature.
- **Fix:** Change the 3-arg parameter to `const Mobilogram& mobilogram` to match the 2-arg overload and the PeakPickerChromatogram sibling. This is an ABI break (mangled signature changes) but source-compatible for nearly all callers; if ABI must be preserved, add a new const& overload and deprecate the by-value one.

### [ANSW-20] PeakIntegrator::estimateBackground — estimateBackground docs name a baseline type "vertical_sum" that does not exist
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h · _osw-chrom_

```cpp
PeakBackground estimateBackground(const MSChromatogram& chromatogram, const double left, const double right, const double peak_apex_pos) const
```
- **Expectation:** A caller reading estimateBackground's documentation expects the named baseline types to be the actual accepted values of the `baseline_type` parameter.
- **Actual:** All four estimateBackground doxygen blocks say "The user can choose to compute one of two background types: \"vertical_sum\" and \"base_to_base\"". But "vertical_sum" appears nowhere in the implementation or the parameter's valid strings; the real accepted values are base_to_base, vertical_division, vertical_division_min, vertical_division_max. A caller who sets baseline_type="vertical_sum" trusting the docs gets rejected by setValidStrings / hits the InvalidParameter throw in estimateBackground_.
- **Evidence:** PeakIntegrator.h lines 341/370/399/428: '...one of two background types: "vertical_sum" and "base_to_base"'. PeakIntegrator.cpp getDefaultParameters: `params.setValidStrings("baseline_type", {"base_to_base","vertical_division","vertical_division_min","vertical_division_max"});` — no "vertical_sum". estimateBackground_ branches only on BASELINE_TYPE_BASETOBASE / VERTICALDIVISION{,_MIN,_MAX}.
- **Fix:** Doc-only fix: replace "vertical_sum" with the actual valid values (vertical_division_min / vertical_division_max / base_to_base) in all four estimateBackground doc blocks. No ABI or signature change.
- **Verifier correction:** Severity adjusted from the claim to low: the misuse fails loudly and immediately (the invalid string is rejected at parameter-set time by setValidStrings, and again at runtime via the InvalidParameter throw at PeakIntegrator.h:734). There is no silent wrong result or data loss — the caller gets an immediate clear error and cannot proceed with a bad computation, so this is a mild surprise rather than a dangerous one. Recommendation stands: doc-only fix replacing "vertical_sum" with the actual accepted values (base_to_base / vertical_division / vertical_division_min / vertical_division_max) in all four estimateBackground doc blocks; no ABI or signature change.

### [ANSW-25] SpectrumAccessSqMass::getNrChromatograms — getNrChromatograms reports the full-file count while spectra are subset, and getChromatogramById always throws
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h · _osw-dataaccess_

```cpp
size_t getNrChromatograms() const override
```
- **Expectation:** Within one accessor, the count getters and the by-id getters share the same index space and honour the same subset, mirroring getNrSpectra()/getSpectrumById(). A non-zero getNrChromatograms() implies getChromatogramById(0..n-1) is callable.
- **Actual:** getNrSpectra() returns the subset size (`sidx_.size()`), but getNrChromatograms() ignores the subset and returns the underlying file's total chromatogram count, while getChromatogramById()/getChromatogramNativeID() unconditionally throw NotImplemented. A caller iterating `for (i=0; i<getNrChromatograms(); ++i) getChromatogramById(i)` (the standard idiom that works on every other ISpectrumAccess) throws on the first iteration.
- **Evidence:** SpectrumAccessSqMass.cpp:209-215 getNrChromatograms returns `handler_.getNrChromatograms()` with a `// TODO: currently chrom indices are not supported`; lines 204-207 and 217-220 throw NotImplemented for getChromatogramById/getChromatogramNativeID.
- **Fix:** Header already documents this; strengthen by making getNrChromatograms() return 0 to match the unusable by-id accessors (so the standard count-then-iterate loop is a no-op rather than a throw), or document the asymmetry at getNrChromatograms() itself. Returning 0 is a behavior change; document-only is the safe additive fix. No signature/ABI change.
- **Verifier correction:** The inconsistency is real but the claim overstates severity. The "always throws" half (getChromatogramById/getChromatogramNativeID) is explicitly documented in the header (lines 54-60, 184-200) and fails loudly via Exception::NotImplemented — not silent wrong data. The genuinely under-documented surprise is narrower: getNrChromatograms() returns the underlying file's full chromatogram count and ignores the spectrum subset (sidx_), whereas getNrSpectra() honours it — the getNrChromatograms() doc claims "Total chromatogram count" with no caveat. This is practically inert because by-id chromatogram access always throws, so the count can never be used to index anything. The canonical count-then-iterate loop does exist in-tree (SpectrumAccessOpenMSInMemory.cpp:32-36) and would throw, but that consumer guards sqMass specifically with a dynamic_cast (lines 20-24); the throw is only reachable via a wrapper such as SpectrumAccessTransforming. Recommendation stands: document the count-vs-subset asymmetry at getNrChromatograms() (safe additive fix), or return 0 to make the loop a no-op (behavior change). No signature/ABI change.

### [ANSW-27] SpectrumAccessOpenMSInMemory::SpectrumAccessOpenMSInMemory(ISpectrumAccess&) — In-memory copy silently drops all chromatograms when the source is an sqMass accessor
`low` · `documentation-overclaim` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h · _osw-dataaccess_

```cpp
explicit SpectrumAccessOpenMSInMemory(OpenSwath::ISpectrumAccess& origin)
```
- **Expectation:** A constructor documented as reading 'every spectrum and chromatogram' and exposing 'exactly the same data as the source' copies the chromatograms too, regardless of the concrete source type.
- **Actual:** There is a SqMass fast path that calls `getAllSpectra(...)` only; the chromatogram-copy loop is in the else branch and is skipped entirely for SqMass sources. The resulting in-memory view reports getNrChromatograms()==0 even if the source file contained chromatograms — but for SqMass this is masked because SqMass's own getChromatogramById throws anyway, so the data loss is invisible until a different source type is wrapped the same way.
- **Evidence:** SpectrumAccessOpenMSInMemory.cpp:19-37: the `if (dynamic_cast<SpectrumAccessSqMass*>...)` branch only fills spectra_/spectra_meta_; the chromatogram loop (`for i<getNrChromatograms() chromatograms_.push_back(...)`) is in the `else` only.
- **Fix:** Header claims 'exactly the same data as the source' and 'every spectrum and chromatogram'; either copy chromatograms in the SqMass branch too, or narrow the doc to say chromatograms are taken only from the generic path. Internal .cpp/doc fix, no ABI change.
- **Verifier correction:** The SqMass fast path skips the chromatogram-copy loop, so an in-memory view built from a SpectrumAccessSqMass reports getNrChromatograms()==0 even when the source file had chromatograms. But this is NOT recoverable silent data loss versus a working alternative: SpectrumAccessSqMass itself cannot serve chromatograms at all (getChromatogramById/getChromatogramNativeID throw NotImplemented, as documented in SpectrumAccessSqMass.h), so the generic path would throw for a SqMass source rather than copy them — the fast path avoids that crash. The genuine, minor issue is that the SpectrumAccessOpenMSInMemory header overclaims "every spectrum and chromatogram"/"exactly the same data as the source" while silently dropping the chromatogram count for SqMass sources. Recommendation: narrow the InMemory doc string to state chromatograms are copied only via the generic ISpectrumAccess path (SqMass sources carry spectra only). The claim's stated recommendation (copy chromatograms in the SqMass branch) is infeasible. Internal .cpp/doc fix, no ABI change.

### [ANSW-29] OpenSwathDataAccessHelper::convertToOpenMSChromatogramFilter — Output chromatogram is the first parameter while the sibling convert* helpers put the output last
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h · _osw-dataaccess_

```cpp
static void convertToOpenMSChromatogramFilter(MSChromatogram& chromatogram, const ChromatogramPtr& cptr, double rt_min, double rt_max)
```
- **Expectation:** Within the same helper class, the conversion functions use a consistent argument order. The non-filter overload is convertToOpenMSChromatogram(const ChromatogramPtr& cptr, MSChromatogram& chromatogram) — source first, OpenMS output last.
- **Actual:** convertToOpenMSChromatogramFilter reverses this: the OpenMS output `chromatogram` is the FIRST parameter and the source `cptr` is second, with two trailing same-typed doubles (rt_min, rt_max). Both the source and destination are easily swappable at the call site, and the order is inconsistent with convertToOpenMSChromatogram / convertToOpenMSSpectrum (output last).
- **Evidence:** DataAccessHelper.h:37 `convertToOpenMSChromatogram(const OpenSwath::ChromatogramPtr& cptr, OpenMS::MSChromatogram & chromatogram)` vs line 42-45 `convertToOpenMSChromatogramFilter(OpenMS::MSChromatogram & chromatogram, const OpenSwath::ChromatogramPtr& cptr, double rt_min, double rt_max)`.
- **Fix:** Add an overload with the consistent (cptr, chromatogram, rt_min, rt_max) order and deprecate the current ordering; or at minimum document the inversion. Additive overload keeps ABI; changing the existing signature would be source-breaking.
- **Verifier correction:** The output/source ordering of convertToOpenMSChromatogramFilter (chromatogram first) is inconsistent with its siblings convertToOpenMSChromatogram and convertToOpenMSSpectrum (output last), and is undocumented; the inconsistency is visibly confusing (adjacent opposite-order calls in MRMFeatureFinderScoring.cpp:1398/1402). But the source/dest pair is NOT silently swappable — MSChromatogram& vs const ChromatogramPtr& are distinct types, so a transposition fails to compile rather than corrupting data. Hence low severity. abi_impact for the recommended additive overload is 'none'; only altering the existing signature would be source-breaking.

### [ANSW-73] OpenSwathOSWWriter::addRun / setRunId — run_id is silently mutated (sign bit cleared) before being stored/written
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h · _osw-export_

```cpp
void addRun(const UInt64 run_id, const std::string& input_filename); void setRunId(const UInt64 run_id)
```
- **Expectation:** A UInt64 run_id passed to addRun/setRunId is stored and written verbatim as the RUN.ID / FEATURE.RUN_ID key.
- **Actual:** Both methods apply Internal::SqliteHelper::clearSignBit(run_id) before storing: addRun binds/stores 'rid = clearSignBit(run_id)' (cpp:951-983) and setRunId sets 'run_id_ = clearSignBit(run_id)' (cpp:986-989). A run_id with the high bit set is silently transformed, so the value read back from the OSW file differs from what was passed, and two ids differing only in the top bit collide.
- **Evidence:** cpp:951 'const UInt64 rid = Internal::SqliteHelper::clearSignBit(run_id);' ; cpp:988 'run_id_ = Internal::SqliteHelper::clearSignBit(run_id);'
- **Fix:** Document in the header that the top/sign bit of run_id is cleared (to keep SQLite INTEGER storage non-negative for joins). ABI-safe doc-only fix.
- **Verifier correction:** run_id IS mutated (bit 63 cleared) before storage/write, and this is undocumented in the addRun/setRunId header doc. But the title's "silently wrong" framing overstates impact: clearSignBit is an established OpenMS-wide convention (documented at SqliteConnector.h:256-258) applied identically to RUN.ID, FEATURE.RUN_ID and all feature uniqueIds to keep SQLite INTEGER joins valid; the consistency across keys means joins work, and the practical collision/divergence risk is negligible because run_ids come from UniqueIdGenerator. Correct severity is low (mild undocumented surprise), category hidden-side-effect, ABI none. Fix: document bit-clearing in the header comments for addRun/setRunId.

### [ANSW-74] OpenSwathOSWWriter::getSeparateScore — Doc @brief says 'concatenated scores' but the method (and its name) return separated per-element scores
`low` · `doc-typo` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h · _osw-export_

```cpp
std::vector<std::string> getSeparateScore(const Feature& feature, const std::string& score_name) const
```
- **Expectation:** The method name getSeparateScore and return type std::vector<std::string> imply one entry per score-list element (separated).
- **Actual:** The doc @brief contradicts the name: '@brief Prepare concatenated scores for SQLite insertion'. The implementation iterates the score list and pushes one string per element (cpp:996-1006), i.e. separated, matching the name but not the brief. The contradictory @brief will mislead readers about whether the result is one concatenated string or a per-transition vector.
- **Evidence:** header line 240 '@brief Prepare concatenated scores for SQLite insertion' vs cpp:1001-1004 loop 'separated_scores.push_back(oswValueToString(score_values.at(i)));'
- **Fix:** Fix the @brief to say 'separated' (one entry per list element). Doc-only, ABI-safe.
- **Verifier correction:** The @brief on line 240 says "concatenated" but the method returns separated per-element scores. This is a stale doc word, not a misleading name: the method name (getSeparateScore), the return type (std::vector<std::string>), and the @returns line ("A vector of strings with the queried scores", line 247) all correctly convey "separated". The word "concatenated" should be changed to "separated". Doc-only, ABI-safe, low impact.

### [ANSW-23] OpenSwathPeptidoformInference::preparePrecursorBM — Doc says it emits true/false rows per feature, but a feature with both MS1 and MS2 precursor evidence yields four rows
`low` · `api-contract-imprecision` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathPeptidoformInference.h · _osw-inference_

```cpp
static std::vector<BayesianModelRow> preparePrecursorBM(const std::vector<IPFPrecursorRow>& rows)
```
- **Expectation:** From 'For each feature this emits rows for the true (1) and false (0) precursor hypotheses', a caller expects two rows per feature (one true, one false).
- **Actual:** If both ms1_precursor_pep and ms2_precursor_pep are present, the function emits FOUR rows per feature (a true/false pair for each evidence source); applyBM later multiplies the two same-hypothesis evidence terms together. The header gives no indication that the row count depends on how many precursor-evidence sources are populated.
- **Evidence:** OpenSwathPeptidoformInference.cpp:176-189: separate `if (row.ms1_precursor_pep.has_value()) { bm_rows.push_back(...1...); bm_rows.push_back(...0...); }` and `if (row.ms2_precursor_pep.has_value()) { bm_rows.push_back(...1...); bm_rows.push_back(...0...); }` both execute when both are present.
- **Fix:** Clarify in the header that one true/false pair is emitted per available precursor-evidence source (so up to four rows when both MS1 and MS2 PEPs are set) and that applyBM combines same-hypothesis evidence multiplicatively. Doc-only, ABI-safe.
- **Verifier correction:** preparePrecursorBM emits one true(1)/false(0) row pair PER available precursor-evidence source, so up to four rows per feature when both MS1 and MS2 PEPs are present (two rows when only one source is set, and a synthetic fallback pair when neither is set, cpp:201-202). The header doc implies two rows per feature and should state the per-source cardinality. This is a pure-function output-cardinality/doc-contract imprecision, not a hidden side-effect: applyBM correctly collapses the rows by (feature_id, hypothesis) and multiplies same-hypothesis evidence (the multiplicative combination is already documented on applyBM), and precursorInference consumes only the hypothesis==1 posterior, so no wrong results arise from intended use. Severity is low; the fix is a doc clarification on preparePrecursorBM, which is ABI-safe.

### [ANSW-62] MRMAssay::detectingTransitions / detectingTransitionsLight — min_transitions documented as @param[out] though it is a plain by-value int input
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h · _osw-mrm-assay_

```cpp
void detectingTransitions(OpenMS::TargetedExperiment& exp, int min_transitions, int max_transitions)
```
- **Expectation:** A parameter marked '@param[out]' is expected to be a non-const reference/pointer the function writes a result into; a competent caller may pass an uninitialized variable expecting it to be filled.
- **Actual:** min_transitions (and max_transitions) are pass-by-value int inputs that are only read (used as the minimum/maximum transition count). The '@param[out]' tag is simply wrong and could mislead a reader into thinking the function returns a count there.
- **Evidence:** Header MRMAssay.h:112 `@param[out] min_transitions the minimum number of transitions required per assay` and :239 same for the Light version; MRMAssay.cpp:922 `if (m->second.size() >= (Size)min_transitions)` shows it is read-only.
- **Fix:** ABI-neutral docs fix: change '@param[out] min_transitions' to '@param[in] min_transitions' (and likewise for the Light variant). No code/ABI change.
- **Verifier correction:** The claim is accurate but the severity is low, not a functional hazard: because min_transitions is declared by value (int, not int& or int*), the @param[out] tag is a self-evidently inert documentation error. The compiler/type system prevents any write-back, so a caller cannot actually be led into a silent bug — at most mild reader confusion. Fix: change `@param[out] min_transitions` to `@param[in] min_transitions` at MRMAssay.h:112 and :239.

### [ANSW-66] MRMRTNormalizer::removeOutliersRANSAC / removeOutliersIterative — @return clause is logically inverted: an outlier-removal routine throws instead of returning when the fit is bad
`low` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h · _osw-mrm-assay_

```cpp
static std::vector<std::pair<double,double>> removeOutliersRANSAC(const std::vector<std::pair<double,double>>& pairs, double rsq_limit, double coverage_limit, size_t max_iterations, double max_rt_threshold, size_t sampling_size)
```
- **Expectation:** The documented return contract reads: 'A vector of pairs is returned if the R^2 limit was reached without reaching the coverage limit. If the limits are reached, an exception is thrown.' A reader parses this as 'success returns, hitting the targets throws' — backwards — and would not expect a routine named removeOutliers* to throw on the common 'fit not good enough' case.
- **Actual:** The function returns the corrected pairs on success and throws Exception::UnableToFit when the R^2 limit is NOT reached (bestrsq < rsq_limit) or when too few points survive (coverage). The doc wording ('if the R^2 limit was reached', 'If the limits are reached, an exception is thrown') contradicts the code's actual condition.
- **Evidence:** MRMRTNormalizer.cpp:47 `if (bestrsq < rsq_limit){ throw Exception::UnableToFit(... "WARNING: rsq: " ...)` and :56 `if (new_pairs.size() < d){ throw ... }`; same throw-on-failure pattern in removeOutliersIterative at :193. Header @return text at MRMRTNormalizer.h:82-87 and :115-120.
- **Fix:** ABI-neutral docs fix: rewrite the @return/@exception to: 'Returns the outlier-filtered pairs when both the R^2 and coverage limits are satisfied; throws Exception::UnableToFit if the R^2 limit cannot be met or coverage drops below coverage_limit.' No signature change.
- **Verifier correction:** The @return text in the header (lines 82-87 for RANSAC, 115-120 for Iterative) is internally contradictory and inverted relative to the code: it says "returned if the R^2 limit was reached without reaching the coverage limit. If the limits are reached, an exception is thrown," whereas the code returns when both limits ARE satisfied and throws Exception::UnableToFit when R^2 < rsq_limit or surviving points < coverage_limit*size. The adjacent @exception clause already states the correct behavior. This is a docs-only defect (the throw itself is documented and is expected by callers), not a hidden/surprising throw; severity low. Suggested fix: "@return Returns the outlier-filtered pairs when both the R^2 and coverage limits are satisfied; throws Exception::UnableToFit if the R^2 limit cannot be met or coverage drops below coverage_limit." No signature/ABI change.

### [ANSW-67] MRMRTNormalizer::removeOutliersRANSAC — Throws UnableToFit for too-few input/sample points (<30 input, <5 sampled) regardless of rsq/coverage limits
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h · _osw-mrm-assay_

```cpp
static std::vector<std::pair<double,double>> removeOutliersRANSAC(..., double max_rt_threshold, size_t sampling_size)
```
- **Expectation:** Given the documented parameters (rsq_limit, coverage_limit, sampling_size), a caller expects the routine to attempt RANSAC on whatever data is supplied and only fail when the quality thresholds are not met.
- **Actual:** Before any fitting, the method hard-throws if sampling_size < 5 or pairs.size() < 30. These hidden minimum-count gates are not mentioned in the header at all; a caller passing a small but otherwise valid RT-peptide set gets an exception that is unrelated to rsq_limit/coverage_limit.
- **Evidence:** MRMRTNormalizer.cpp:29 `if (n < 5){ throw Exception::UnableToFit(... "is below limit of 5 peptides required ...")` and :36 `if (pairs.size() < 30){ throw ... "is below limit of 30 peptides required ...")`.
- **Fix:** ABI-neutral docs fix: document the hard minimums (>=5 sampled, >=30 input pairs) in the header and the @exception clause so the throw is not a surprise. No signature/ABI change.

### [ANSW-68] MRMDecoy::switchKR / switchKRLight — Name says 'switch K<->R' but on a non-K/R C-terminus it mutates to an arbitrary (hash-chosen) amino acid
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h · _osw-mrm-assay_

```cpp
void switchKR(OpenMS::TargetedExperiment::Peptide& peptide) const; static void switchKRLight(std::string& sequence)
```
- **Expectation:** switchKR reads as 'swap the terminal K for R or R for K'. A caller would expect it to be a no-op (or a documented swap) and certainly not to replace the residue with an unrelated amino acid.
- **Actual:** If the last residue is neither K nor R, switchKRLight replaces it with an amino acid selected by an FNV-1a hash of the sequence (one of 17 residues). This silent identity-changing mutation is only mentioned in a @note, not in the name. The behavior is reachable from the public generateDecoys 'switchKR' flag path.
- **Evidence:** MRMDecoy.cpp:308-329: `if (sequence[lastAA]=='K'){...='R';} else if (...=='R'){...='K';} else { ... hash ... sequence[lastAA] = aa[res_pos][0]; }`. Header @note at MRMDecoy.h:150 and :248 documents the random replacement only in prose.
- **Fix:** ABI-neutral: the @note already exists; strengthen it (and the generateDecoys doc) to state that non-tryptic C-termini are replaced by a deterministically-chosen amino acid, so the precursor-mass-shift side effect is explicit. Renaming would break ABI and is not recommended.
- **Verifier correction:** Accurate as stated, with two refinements: (1) the mutation affects DECOY sequences only (the intended output is deliberately a different/fake sequence), not target/truth data, so it does not silently corrupt analysis results — it produces a decoy that is more different than the name implies. (2) The behavior is documented in an adjacent header @note (h:150, h:248), though NOT in the method name, signature, or the user-facing TOPP "switchKR" param description (OpenSwathDecoyGenerator.cpp:126), which mention only the K<->R swap. Recommendation (doc-strengthening of the @note and the TOPP param help, and optionally the generateDecoys docstring) is correct and ABI-neutral; renaming would break ABI and is not warranted at this severity.

### [ANSW-69] MRMIonSeries::getIon — Read-only lookup takes IonSeries by non-const reference (and the method itself is non-const)
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMIonSeries.h · _osw-mrm-assay_

```cpp
std::pair<std::string, double> getIon(IonSeries& ionseries, const std::string& ionid)
```
- **Expectation:** A 'getIon' that merely selects an existing ion from a series reads like a const query that should accept a const IonSeries& and not require a mutable object.
- **Actual:** The parameter is a non-const IonSeries& and the member function is non-const, so a caller holding a `const IonSeries` (or calling on a const MRMIonSeries) cannot use it, even though the body only reads (contains() + operator[] on an already-present key). The non-const ref also leaves open the impression that the series is modified.
- **Evidence:** MRMIonSeries.cpp:22-32: body only does `if (ionseries.contains(ionid)) return make_pair(ionid, ionseries[ionid]); else return make_pair("unannotated", -1);` — no mutation. Header declares `getIon(IonSeries& ionseries, ...)` (non-const param, non-const method).
- **Fix:** Source-compatible improvement: change the parameter to `const IonSeries&` (replace operator[] with .at()/.find()) and make the method const. This is source-compatible for callers but changes the mangled signature, so it is ABI-breaking; if ABI must be frozen, leave as-is and at minimum document that the series is not modified.
- **Verifier correction:** Minor amplification of the claim: the misleading-mutation impression is not merely "left open" by the non-const ref — the Doxygen comment explicitly tags the parameter `@param[in,out]` (line 62), falsely documenting an output/mutation that the body never performs. Severity is LOW, not higher: the code cannot produce wrong results, data loss, or crashes (operator[] on a present key is safe), so the only impact is reduced usability with const objects plus a misleading signature/doc. On ABI: the recommendation's framing is right that the SOURCE change is compatible (all current callsites — MRMDecoy.cpp:684, :1069, and the unit test — pass mutable lvalues, and adding const-qualification only relaxes requirements). I classify abi_impact as source-compatible to denote that callers need no source changes; note however that because this is a public OPENMS_DLLAPI member, altering the parameter type and adding the const method qualifier changes the mangled name and IS a binary-compatibility break, so it should be batched into a release where ABI breaks are permitted. The suggested fix (const IonSeries& param + const method + replace operator[] with .at()/.find()) is correct and is the reason the non-const ref exists today (operator[] is non-const on std::map).

### [ANSW-14] MRMFeatureFilter::FilterFeatureMap — FilterFeatureMap overwrites its '[in]' features argument in place (it is really in/out)
`low` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h · _osw-mrm-feature_

```cpp
void FilterFeatureMap(FeatureMap& features, const MRMFeatureQC& filter_criteria, const TargetedExperiment& transitions)
```
- **Expectation:** The Doxygen tags it '@param[in] features FeatureMap to flag or filter', so a caller would expect the input map to be read-only and the filtered result delivered elsewhere. There is no separate output parameter, so a reader might even assume nothing is written back.
- **Actual:** The implementation builds a local features_filtered and then assigns features = features_filtered; (MRMFeatureFilter.cpp line 300), destructively replacing the caller's input. The same in-place overwrite pattern is in FilterFeatureMapPercRSD (line 493) and FilterFeatureMapBackgroundInterference (line 615).
- **Evidence:** // from MRMFeatureFilter.cpp\n    FeatureMap features_filtered;\n    ...\n      features = features_filtered;   // line 300, overwrites the [in] argument
- **Fix:** Re-tag the parameter as @param[in,out] features (and likewise for the two sibling methods) in the header so the in-place mutation is documented. ABI-neutral. If a future major version allows, splitting into separate in/out params would be clearer but is a breaking change.
- **Verifier correction:** The `@param[in] features` tag is wrong: `features` is in/out. The mutation is broader than the claim states — it happens in BOTH modes, not only "filter": in "flag" mode FilterFeatureMap writes QC metavalues onto features and their subordinates (e.g. cpp lines 244, 253-254, 262-263, 268, 283-284, 292-293), and in "filter" mode it overwrites the whole map via `features = features_filtered;` (line 300). The two sibling methods FilterFeatureMapPercRSD (line 493) and FilterFeatureMapBackgroundInterference (line 615) repeat the overwrite pattern and also carry the misleading `@param[in]` tag. Mitigating context: all three take a non-const `FeatureMap&`, which already signals possible mutation, so this is a documentation defect rather than a silent trap. Fix is ABI-neutral: re-tag as `@param[in,out] features`.

### [ANSW-15] MRMFeatureFilter::calculateRTDifference — calculateRTDifference / calculateResolution take non-const Feature& but only read them
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h · _osw-mrm-feature_

```cpp
double calculateRTDifference(Feature& component_1, Feature& component_2) const
```
- **Expectation:** A read-only 'calculate' query on two features should take const Feature& so callers can pass const features and temporaries, and so the signature signals it does not mutate its inputs.
- **Actual:** Both take non-const Feature& even though their bodies only call getRT()/getWidth()/getMetaValue()/metaValueExists(). calculateRTDifference is just std::abs(component_1.getRT() - component_2.getRT()). This forces callers to hold mutable, lvalue Feature objects and misleadingly implies the features may be modified.
- **Evidence:** double MRMFeatureFilter::calculateRTDifference(Feature& component_1, Feature& component_2) const { return std::abs(component_1.getRT() - component_2.getRT()); }  // and calculateResolution(Feature&, Feature&) which only reads getRT()/getWidth()/getMetaValue()
- **Fix:** Change the parameters to const Feature&. This is source-compatible for all current call sites (they already pass lvalues) and only widens what callers may pass. Strictly it is an ABI change to the mangled symbol, but as a non-virtual public method on a library type the practical risk is low; do it in a minor release.

### [ANSW-42] MRMFeatureFilter::FilterFeatureMap — FilterFeatureMap mutates its 'features' argument in place but documents it as @param[in]
`low` · `documentation-mismatch` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h · _osw-mrm-feature_

```cpp
void FilterFeatureMap(FeatureMap& features, const MRMFeatureQC& filter_criteria, const TargetedExperiment& transitions)
```
- **Expectation:** The Doxygen tag `@param[in] features FeatureMap to flag or filter` marks features as input-only; a caller would expect the FeatureMap to be unchanged and look for the result elsewhere.
- **Actual:** The method modifies `features` in place: in flag mode it writes meta-values/flags onto the features, and in filter mode it replaces the entire map with the filtered copy (`features = features_filtered;`). There is no separate output parameter, so the documented `[in]` is wrong — it is really `[in,out]`.
- **Evidence:** Header line 75-76 `@param[in] features FeatureMap to flag or filter`; impl src/openms/source/ANALYSIS/OPENSWATH/MRMFeatureFilter.cpp:298-301 `if (flag_or_filter_ == "filter") { features = features_filtered; }`. Same pattern in FilterFeatureMapPercRSD / FilterFeatureMapBackgroundInterference.
- **Fix:** Change the Doxygen annotation to `@param[in,out] features` for all three Filter* methods. Pure doc fix; no signature/ABI change.
- **Verifier correction:** The `@param[in] features` Doxygen tag is wrong for all three Filter* methods (FilterFeatureMap, FilterFeatureMapPercRSD, FilterFeatureMapBackgroundInterference): the function both writes flags/meta-values onto the input features (flag mode) and reassigns the entire map (`features = features_filtered;`) in filter mode, so it is `[in,out]`. However, the surprise is low-severity, not the implied silent-misuse trap: the parameter is a non-const `FeatureMap&` (which already signals possible mutation to a competent C++ dev), and there is no separate output parameter, so a caller has nowhere else to look for the result. Fix is purely cosmetic: change `@param[in] features` to `@param[in,out] features` (and `@param[out] transitions` on the const-ref `transitions` is also wrong and should be `@param[in]`).

### [ANSW-44] MRMFeatureFinderScoring::prepareProteinPeptideMaps_ — Public method carries a trailing-underscore (private-convention) name and is a mandatory pre-call step
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h · _osw-mrm-feature_

```cpp
void prepareProteinPeptideMaps_(const OpenSwath::LightTargetedExperiment& transition_exp)
```
- **Expectation:** In OpenMS the trailing underscore marks a protected/private/internal member (e.g. updateMembers_, computeScore_). A caller scanning the public section would not expect to have to call a `_`-suffixed method, and tooling/bindings treat such names as internal.
- **Actual:** prepareProteinPeptideMaps_ is declared in the public section and its own doc states `Calling this method _is_ required before calling scorePeakgroups`, i.e. it is a required part of the public call sequence yet is named like an internal helper.
- **Evidence:** Header lines 135-142: public section, `/** ... Calling this method _is_ required before calling scorePeakgroups. */ void prepareProteinPeptideMaps_(const OpenSwath::LightTargetedExperiment& transition_exp);`
- **Fix:** Provide a public alias without the underscore (e.g. prepareProteinPeptideMaps) delegating to the existing one, and document the required ordering on scorePeakgroups. Additive alias is source/ABI-compatible; renaming in place would be source-breaking.

### [ANSW-45] MRMFeatureQC::ComponentGroupPairQCs — ComponentGroupPairQCs bounds are left uninitialized, unlike the sibling QC structs which default to {0.0 .. 1e12}
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h · _osw-mrm-feature_

```cpp
struct ComponentGroupPairQCs { double resolution_l; double resolution_u; double rt_diff_l; double rt_diff_u; ... }
```
- **Expectation:** The two sibling structs ComponentQCs and ComponentGroupQCs initialize every bound member (e.g. `retention_time_l { 0.0 }`, `retention_time_u { 1e12 }`). A developer default-constructing a ComponentGroupPairQCs reasonably expects the same safe defaults.
- **Actual:** resolution_l, resolution_u, rt_diff_l, rt_diff_u have NO default member initializers. A default-constructed ComponentGroupPairQCs holds indeterminate doubles, so a forgotten field silently feeds garbage bounds into checkRange/resolution filtering rather than a benign wide-open default.
- **Evidence:** Header lines 250-265: `struct ComponentGroupPairQCs { ... double resolution_l; double resolution_u; double rt_diff_l; double rt_diff_u; };` (no `{ ... }` initializers), vs ComponentQCs lines 115-125 `double retention_time_l { 0.0 }; double retention_time_u { 1e12 }; ...`.
- **Fix:** Add brace-initializers matching the siblings, e.g. `resolution_l{0.0}; resolution_u{1e12}; rt_diff_l{0.0}; rt_diff_u{1e12};`. Source-compatible; changes default-constructed values but those were previously UB so no sane caller relied on them.
- **Verifier correction:** ComponentGroupPairQCs indeed lacks the default member initializers (and operator==/!=) its sibling QC structs have, so a default-constructed instance holds indeterminate doubles in its four bound fields — a genuine convention inconsistency. But contrary to the claim, these members (resolution_l/u, rt_diff_l/u, resolution_pair_name) and the component_group_pair_qcs vector are never read by any C++ code (MRMFeatureFilter, MRMFeatureQCFile, TOPP tools); they are only assigned in a pyOpenMS test. Thus the indeterminate values cannot currently feed garbage bounds into any filtering/checkRange logic — the UB is latent rather than active, making this a low-severity inconsistency/latent-trap, not a high-severity silently-wrong-result bug. The recommended fix (add resolution_l{0.0}/resolution_u{1e12}/rt_diff_l{0.0}/rt_diff_u{1e12}, and arguably operator==/!= for parity) is correct and source-compatible.

### [ANSW-55] DIAScoring::dia_by_ion_score — Scoring method takes AASequence by non-const reference but never modifies it
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h · _osw-scoring-core_

```cpp
void dia_by_ion_score(const SpectrumSequence& spectrum, AASequence& sequence, int charge, const RangeMobility& im_range, double& bseries_score, double& yseries_score) const
```
- **Expectation:** A read-only b/y-ion counting score taking the peptide as `AASequence& sequence` signals that it mutates the sequence; callers cannot pass a `const AASequence` or a temporary.
- **Actual:** The body only calls `sequence.toString()` and passes `sequence` to `DIAHelpers::getBYSeries(const AASequence&, ...)`. It never mutates the sequence. The non-const reference is gratuitous and forces callers to keep a mutable lvalue.
- **Evidence:** DIAScoring.cpp:373 `const std::string sequence_key = sequence.toString();` and :380 `OpenMS::DIAHelpers::getBYSeries(sequence, ...)`; DIAHelper.h:96 `void getBYSeries(const AASequence& a, ...)` takes const. No assignment to `sequence` anywhere in the function.
- **Fix:** Change the parameter to `const AASequence& sequence`. This is source-compatible for existing callers passing lvalues and additionally enables const/temporary arguments. Minor mangled-signature change for the (rare) external caller.

### [ANSW-56] DIAScoring::dia_isotope_scores / score_with_isotopes — DIA isotope/score functions take the spectrum sequence by non-const reference but treat it read-only
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h · _osw-scoring-core_

```cpp
void dia_isotope_scores(... SpectrumSequence& spectrum, ...) const; void score_with_isotopes(SpectrumSequence& spectrum, ...) const;
```
- **Expectation:** Taking `SpectrumSequence& spectrum` (non-const) in a scoring routine signals the spectrum is modified (e.g. sorted/filtered in place); a caller must hold a mutable copy.
- **Actual:** `dia_isotope_scores` forwards `spectrum` to `diaIsotopeScoresSub_(const SpectrumSequence&, ...)` and `score_with_isotopes` forwards it to `DiaPrescore::scorePrepared(const SpectrumSequence&, ...)`; neither mutates the sequence. The non-const reference is inconsistent with the const overloads in the same class (e.g. dia_massdiff_score takes `const SpectrumSequence&`).
- **Evidence:** DIAScoring.cpp:241 `diaIsotopeScoresSub_(transitions, spectrum, ...)` with :424 `diaIsotopeScoresSub_(..., const SpectrumSequence& spectrum, ...)`; DIAScoring.cpp:418 `DiaPrescore::scorePrepared(spectrum, ...)` with DIAPrescoring.h:118 `static void scorePrepared(const SpectrumSequence& spec, ...)`.
- **Fix:** Take `const SpectrumSequence&` to match dia_massdiff_score / dia_ms1_massdiff_score in the same header and to stop implying in-place mutation. ABI-breaking for these two symbols only; could be done additively by overloading.
- **Verifier correction:** Claim is accurate. Two functions are non-const-ref but read-only and inconsistent with their const-ref siblings. Note the original expectation ("a caller must hold a mutable copy") slightly overstates the harm: SpectrumSequence is a vector of shared_ptr, so a const& would not actually grant deep immutability of the spectra, and the sole production caller already passes a mutable pooled buffer, so no copy is forced. The fix to `const SpectrumSequence&` is ABI-breaking (mangled-name change) for exactly these two symbols, as claimed, and is source-compatible for lvalue callers; can be done additively via overload.

### [ANSW-57] OpenSwathScoring::calculateLibraryScores — calculateLibraryScores is non-const while its sibling calculate*Scores methods are const
`low` · `missed-const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h · _osw-scoring-core_

```cpp
void calculateLibraryScores(OpenSwath::IMRMFeature*, const std::vector<TransitionType>&, const CompoundType&, const double, OpenSwath_Scores&);  // non-const
```
- **Expectation:** All `calculate...Scores` entry points in this class compute scores into an out-parameter without mutating the scoring object; calculateChromatographicScores and calculateChromatographicIdScores are declared const, so a developer expects calculateLibraryScores to be const too (and to be callable on a const OpenSwathScoring).
- **Actual:** calculateLibraryScores is declared non-const, yet its body only reads `su_` and `rt_normalization_factor_` and calls the static methods `MRMScoring::calcLibraryScore` and `MRMScoring::calcRTScore`. It mutates nothing, so the missing const is an unexplained asymmetry that needlessly prevents calling it on a const instance.
- **Evidence:** OpenSwathScoring.h:173-177 declares it without `const`; :127-132 and :152-156 declare the chromatographic siblings `const`. OpenSwathScoring.cpp:614-641 body reads `su_.use_library_score_`, `rt_normalization_factor_` and calls only static MRMScoring methods.
- **Fix:** Add `const` to calculateLibraryScores to match its siblings. Source-compatible for callers; technically mangles the symbol so flag as a (minor) ABI change.
- **Verifier correction:** calculateLibraryScores does mutate nothing and could be const, but the claim wrongly frames non-const as breaking a class-wide const convention. Of the six calculate*Scores entry points, only two (the chromatographic pair) are const; the three DIA variants are non-const and legitimately so (they mutate the fetch_spectrum_tmp member). Non-const is therefore the majority pattern, so a developer does NOT strongly expect calculateLibraryScores to be const. This is a minor const-correctness omission (low severity), not a medium inconsistent-convention defect. Adding const is source-compatible for the single non-const caller and only minorly mangles the exported symbol.

### [ANSW-58] OpenSwathScoring::calculateDIAIdScores (drift_target parameter) — drift_target documented as [out] but passed by value as const double (cannot return anything)
`low` · `doc-vs-signature-mismatch` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h · _osw-scoring-core_

```cpp
void calculateDIAIdScores(..., const double drift_target, ...);  // doc: @param[out] drift_target
```
- **Expectation:** A parameter annotated `@param[out] drift_target target drift value` should be writable and observable by the caller after the call.
- **Actual:** The parameter is declared `const double drift_target` (by value, const). It is physically impossible for the function to communicate a result back through it; the `[out]` Doxygen tag is wrong and will mislead callers into thinking the function returns the drift target. (The sibling calculateDIAScores documents the same-typed parameter as `[in]`.)
- **Evidence:** OpenSwathScoring.h:256 `const double drift_target,` vs the doc block at :244 `@param[out] drift_target target drift value`. Compare :207 `const double drift_target` documented at :192 as `@param[in]`.
- **Fix:** Fix the Doxygen annotation to `@param[in] drift_target`. No code/ABI change. If an actual output of the resolved drift target is intended, add a separate `double& drift_target_out` parameter (additive overload).
- **Verifier correction:** The drift_target parameter of calculateDIAIdScores is annotated @param[out] (header line 244) but declared `const double drift_target` by value (line 256), making any output through it impossible. The implementation confirms it is used only as input. The fix is to change the annotation to @param[in], matching the identical parameter in the sibling calculateDIAScores. This is a documentation-only inconsistency with low severity: the const-by-value type guarantees no caller can be silently given wrong data — at worst they observe an unchanged input. The actual drift output is returned via the OpenSwath_Scores& scores object.

### [ANSW-59] OpenSwathScoring::OpenSwathScoring (default constructor) — Default constructor leaves several scoring members uninitialized; scoring before initialize() is UB
`low` · `silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/OPENSWATH/OpenSwathScoring.cpp · _osw-scoring-core_

```cpp
OpenSwathScoring();
```
- **Expectation:** After default-constructing an OpenSwathScoring, member state is fully defined; calling a scoring/fetch method without an explicit initialize() should at worst use documented defaults.
- **Actual:** The constructor initializes only a subset of members; `merge_spectra_by_peak_width_fraction_`, `use_ms1_ion_mobility_`, and `apply_im_peak_picking_` are left uninitialized and are only assigned inside initialize(). The class exposes no indication that initialize() is mandatory before scoring, so a caller who skips it reads indeterminate values (e.g. use_ms1_ion_mobility_ branches in calculateDIAScores).
- **Evidence:** OpenSwathScoring.cpp:122-130 initializer list sets rt_normalization_factor_, spacing_for_spectra_resampling_, add_up_spectra_, spectra_addition_method_, spectra_merge_method_type_, im_drift_extra_pcnt_ only — not merge_spectra_by_peak_width_fraction_, use_ms1_ion_mobility_, apply_im_peak_picking_ (declared OpenSwathScoring.h:65,71,72). These are assigned only in initialize() (:176,178,179).
- **Fix:** Add default member initializers in the header (e.g. `bool use_ms1_ion_mobility_ = false;`) so the object is always in a valid state. Additive, source- and ABI-compatible (no layout change).
- **Verifier correction:** The constructor does leave merge_spectra_by_peak_width_fraction_, use_ms1_ion_mobility_, and apply_im_peak_picking_ uninitialized, and these are read in calculateDIAScores/calculatePrecursorDIAScores — so invoking those methods without initialize() is UB. However, this is a latent landmine rather than a live bug: every current code path that reaches the DIA scoring methods calls initialize() first (MRMFeatureFinderScoring.cpp:611), and the only caller that skips initialize() (SwathMapMassCorrection) uses only fetchSpectrumSwath, which reads exclusively constructor-initialized members. Adding default member initializers in the header (e.g. bool use_ms1_ion_mobility_ = false;) is the correct, additive, ABI-compatible fix and removes the indeterminate state. Severity is low (mild surprise / invites future misuse), not high.

### [ANSW-6] MasstraceCorrelator::scoreHullpoints — Output params 'lag' and 'lag_intensity' are left untouched when correlation is below min_corr
`low` · `api-contract` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h · _osw-scoring-extra_

```cpp
void scoreHullpoints(const MasstracePointsType& hull_points1, const MasstracePointsType& hull_points2, int& lag, double& lag_intensity, double& pearson_score, const double min_corr, const int max_lag, const double mindiff = 0.1)
```
- **Expectation:** A caller passing three out-params (lag, lag_intensity, pearson_score) expects all three to be written by the call, so they can be read afterwards regardless of the scores.
- **Actual:** When pearson_score <= min_corr the function returns early after assigning only pearson_score; 'lag' and 'lag_intensity' are never written. They retain whatever the caller passed in. In createPseudoSpectra the locals 'int lag; double lag_intensity;' (declared uninitialized and reused across loop iterations) thus carry stale values from a previous high-correlation pair into the next low-correlation pair.
- **Evidence:** MasstraceCorrelator.cpp: 'pearson_score = Math::pearsonCorrelationCoefficient(...);\n if (pearson_score <= min_corr) \n { \n return;\n }\n ... lag = pt->first; ... lag_intensity = pt->second;' — lag/lag_intensity assigned only after the early return. Header documents all three as '@param[out]'.
- **Fix:** Initialize lag = 0 and lag_intensity = 0.0 (or a documented sentinel) at the top of the function before the early-return guard, so the out-params have well-defined values on every code path. Source-only change, no ABI impact.
- **Verifier correction:** scoreHullpoints does leave the out-params 'lag' and 'lag_intensity' unwritten when pearson_score <= min_corr (confirmed: early return at MasstraceCorrelator.cpp:138-141 before the assignments at 145-146), which contradicts the header's @param[out] documentation — a real least-astonishment/API-contract issue. But it does NOT cause a silent wrong-result bug: both in-tree callers (createPseudoSpectra line 283 and ClusterMassTracesByPrecursor.cpp line 255) read 'lag' only behind a short-circuited 'pearson_score > min_correlation && ...' guard whose first term is precisely the negation of the early-return condition, so a stale 'lag' is never consumed; the primary caller documents this ordering dependency in the code (lines 280-282), and the secondary caller declares the locals fresh inside the loop. Category corrected from 'silent-failure' to 'api-contract' and severity from high to low. Recommendation to initialize lag=0/lag_intensity=0 at the top is still a reasonable hardening and is source-only / ABI-neutral.

### [ANSW-11] IonMobilityScoring::driftScoring / driftScoringMS1 / driftIdScoring — 'drift_target' documented as @param[out] but is a by-value const input
`low` · `documentation` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/IonMobilityScoring.h · _osw-scoring-extra_

```cpp
static void driftScoring(..., const double drift_target, ...)
```
- **Expectation:** A Doxygen '@param[out] drift_target Ion Mobility extraction target' tells a caller the function will compute and write back the drift target it actually used.
- **Actual:** In driftScoring, driftScoringMS1 and driftIdScoring the parameter is declared 'const double drift_target' — passed by value and const, so it is physically impossible for the function to return anything through it. It is a pure input. The '@param[out]' annotation contradicts the signature and would mislead a caller into expecting an output they will never receive.
- **Evidence:** Header: '@param[out] drift_target Ion Mobility extraction target' alongside 'static void driftScoring(..., const double drift_target, RangeMobility im_range, ...)' (and identical mismatch in driftScoringMS1 line 106/117 and driftIdScoring line 161/175).
- **Fix:** Change the Doxygen tags to '@param[in] drift_target'. Documentation-only fix; no ABI impact.
- **Verifier correction:** The factual core is correct: `drift_target` is documented `@param[out]` in driftScoring, driftScoringMS1 and driftIdScoring while declared `const double` by value (a pure input), and the implementation only reads it to compute delta scores. The fix (change to `@param[in]`) is correct and doc-only. The severity should be 'low' rather than implied higher: because the parameter is a by-value const double, the compiler makes it impossible to mislead a caller into reading back a value, and the function's real outputs (the `scores` object) are correctly annotated, so the bug is a mild documentation inconsistency with no runtime or ABI consequence. Category is documentation, not return-value (there is no value-returning channel involved).

### [ANSW-12] ConfidenceScoring::scoreAssay_ — 'scoreAssay_' normalizes the caller's feature_intensities in place (non-const ref) as a side effect of scoring
`low` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/ConfidenceScoring.h · _osw-scoring-extra_

```cpp
double scoreAssay_(const TargetedExperiment::Peptide& assay, double feature_rt, DoubleList& feature_intensities, const std::set<std::string>& transition_ids = std::set<std::string>())
```
- **Expectation:** A scoring function named scoreAssay_ taking feature intensities reads as: compute a score from these intensities and return it. The header even documents it only as 'reordered to match the assay's transition list'.
- **Actual:** Beyond reordering, scoreAssay_ calls OpenSwath::Scoring::normalize_sum(feature_intensities) which rescales the caller's DoubleList to sum to 1 in place. Because the same feature_intensities object is passed to scoreAssay_ repeatedly across the true assay and every decoy in scoreFeature_, the vector is re-normalized on each call; the header's '@param[in,out]' note mentions reordering but not the normalization mutation, which a caller reusing the vector elsewhere would not expect.
- **Evidence:** ConfidenceScoring.cpp: 'OpenSwath::Scoring::normalize_sum(feature_intensities); OpenSwath::Scoring::normalize_sum(assay_intensities); double dist_int = manhattanDist_(...);'. Header: '@param[in,out] feature_intensities Feature's transition intensities; reordered to match the assay's transition list.'
- **Fix:** Document the in-place normalization in the @param[in,out] note, or normalize a local copy to avoid mutating the caller's buffer. Both are protected/source-only; no external ABI impact.
- **Verifier correction:** scoreAssay_ does mutate the caller's feature_intensities in place via normalize_sum, but the parameter is a non-const reference explicitly annotated @param[in,out] and in-place intensity normalization is an established convention in this module (cf. NormalizedManhattanDist), so the mutation is signaled, not hidden, and idempotency means it causes no wrong results in the actual (only) caller scoreFeature_. The real issue is that the @param[in,out] description says "reordered to match the assay's transition list" whereas scoreAssay_ actually leaves feature_intensities in element order and instead normalizes it to sum 1; fix the doc text to say "normalized to sum 1 in place." Severity low (doc accuracy only); no ABI impact (protected method).

### [ANSW-80] OpenSwathHelper::computeTransitionGroupId — computeTransitionGroupId returns an empty string (not the input) when the id does not have the expected '..._Precursor_iX' shape
`low` · `undocumented-behavior` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h · _osw-transitions-helper_

```cpp
static std::string computeTransitionGroupId(const std::string& precursor_id)
```
- **Expectation:** Documented as 'reversing the operation performed by computePrecursorId()', a caller expects either the original transition group id back, or an error, for a malformed input.
- **Actual:** If the underscore-split yields fewer than 3 parts the function returns "" with no error, silently producing an empty group id that downstream code may use as a map key. computePrecursorId itself never round-trips cleanly when the group id contains the literal substring '_Precursor_i', but the failure mode here is a silent empty return rather than a signal.
- **Evidence:** if (substrings.size() == 3) return substrings[0];\nelse if (substrings.size() > 3) { ... }\nreturn "";  // OpenSwathHelper.h:61-68
- **Fix:** Document the empty-string return as the failure sentinel in the header, or assert/throw on malformed input. Doc-only is ABI-neutral; behavioral throw would be source-compatible for well-formed callers.
- **Verifier correction:** computeTransitionGroupId does return "" for inputs with fewer than 3 underscore-delimited parts, and the header's @return ("Original transition group id") does not document this empty-string sentinel. But the claim's asserted harm is not realized: both non-test callers guard the result (ChromatogramExtractor.h:201 checks for empty and skips; MSChromatogramParquetConsumer.cpp:231 only runs on ids containing "_Precursor_i" and throws on a metadata-map miss). The round-trip critique is also overstated — ids with embedded underscores reverse correctly (proven by the test "tr_gr2__test_Precursor_i0" -> "tr_gr2__test"). Correct fix is doc-only: note in the header that malformed input returns "" as a sentinel. Severity low, category undocumented-behavior, ABI-neutral.

### [ANSW-82] TargetedSpectraExtractor::untargetedMatching / targetedMatching (spectral_library_score on no-match) — A 'no match found' result is encoded as spectral_library_score = 0.0, indistinguishable from a genuine zero-similarity match
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h · _osw-transitions-helper_

```cpp
void targetedMatching(const std::vector<MSSpectrum>&, const Comparator&, FeatureMap&)
```
- **Expectation:** When no library spectrum clears the threshold, a caller iterating features needs to tell 'no match' apart from 'matched, but score happened to be 0'.
- **Actual:** On no match the method writes spectral_library_name = "" and spectral_library_score = 0.0; a real match whose contrast-angle score is 0 (orthogonal spectra, allowed when min_match_score is 0) produces the same 0.0 score, so the numeric score alone cannot signal failure. The only reliable discriminator (empty name) is documented but easy to miss given the explicit '[0-1]' score doc.
- **Evidence:** features[i].setMetaValue("spectral_library_name", "");\nfeatures[i].setMetaValue("spectral_library_score", 0.0);  // TargetedSpectraExtractor.cpp:797-799
- **Fix:** Document in the header that callers must check spectral_library_name (not score) to detect a miss, or use a sentinel like NaN / a separate boolean meta value for 'matched'. Doc-only fix is ABI-neutral.
- **Verifier correction:** The 0.0 no-match sentinel collides with a genuine 0.0 contrast-angle match ONLY when the user explicitly lowers min_match_score to 0 (default is 0.8, which makes score=0.0 unambiguously "no match"). The no-match value encoding is already documented in the header; what is missing is a warning that the score becomes an ambiguous miss-indicator when min_match_score=0, plus guidance to test spectral_library_name (empty) rather than the score. Recommended fix is doc-only (or an explicit NaN/boolean "matched" meta value) and is ABI-neutral. Severity is low, not high: default config is immune, the condition is narrow, and a documented discriminator (empty name) already exists.

### [ANSW-2] SwathMapMassCorrection::correctMZ — correctMZ Doxygen param directions are inverted: the modified swath_maps is marked [in], the const map is marked [out]
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h · _osw-workflow_

```cpp
void correctMZ(const std::map<std::string, MRMTransitionGroupType*>& transition_group_map, const OpenSwath::LightTargetedExperiment& targeted_exp, std::vector<OpenSwath::SwathMap>& swath_maps, const bool pasef)
```
- **Expectation:** Doxygen @param[out] marks the argument the function writes, @param[in] marks read-only inputs; a competent reader uses these to know which arguments are mutated.
- **Actual:** The doc annotates `transition_group_map` (a `const&` map — cannot be an out parameter) as `@param[out]`, and annotates `swath_maps` as `@param[in]` even though correctMZ REPLACES the shared_ptrs inside swath_maps with SpectrumAccessQuadMZTransforming wrappers (the prose even says 'will be modified (replaced with a corrected version)'). A reader trusting the [in]/[out] tags would assume their swath_maps survive unchanged and that the transition group map is an output, both wrong.
- **Evidence:** Header SwathMapMassCorrection.h:54 '@param[out] transition_group_map' (declared `const ... &`), :56 '@param[in] swath_maps The raw swath maps ... will be modified (replaced with a corrected version)'. Impl SwathMapMassCorrection.cpp:738-739 replaces `swath_maps[i].sptr = std::shared_ptr<...>(new SpectrumAccessQuadMZTransforming(...))`.
- **Fix:** Doc-only, source/ABI-neutral: change `swath_maps` to `@param[in,out]` and `transition_group_map` to `@param[in]`. Same correction applies to correctIM (transition_group_map marked [out] but const&). No signature change needed.
- **Verifier correction:** The Doxygen direction tags on correctMZ (and correctIM) are inverted: transition_group_map is marked @param[out] but is a const& (read-only), and swath_maps is marked @param[in] but is mutated (its sptr members are replaced with SpectrumAccessQuadMZTransforming wrappers per impl 736-741). Fix: swath_maps -> @param[in,out], transition_group_map -> @param[in]. However, the mutation of swath_maps is also stated in plain prose on the very same line ("will be modified (replaced with a corrected version)") and in the function-level description, and the [out]-on-const contradiction is self-evident, so a competent reader is not silently misled — hence low severity, not high. (For correctIM, swath_maps is already marked [in,out].)

### [ANSW-3] SwathMapMassCorrection::correctIM — correctIM marks a const-ref argument as @param[out]
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h · _osw-workflow_

```cpp
void correctIM(const std::map<std::string, MRMTransitionGroupType*>& transition_group_map, const OpenSwath::LightTargetedExperiment& targeted_exp, const std::vector<OpenSwath::SwathMap>& swath_maps, const bool pasef, TransformationDescription& im_trafo)
```
- **Expectation:** @param[out] should only be applied to arguments the function writes through (non-const ref/pointer).
- **Actual:** `transition_group_map` is passed by `const&` yet annotated `@param[out]`. The only genuine output, `im_trafo`, IS correctly marked [out], so the bogus [out] on the const map is purely misleading. A caller may believe correctIM fills the transition group map.
- **Evidence:** Header SwathMapMassCorrection.h:72 '@param[out] transition_group_map' with declaration `const std::map<...>& transition_group_map` (line 78).
- **Fix:** Doc-only fix: change to `@param[in] transition_group_map`. No ABI/source impact.
- **Verifier correction:** The claim is accurate as stated. transition_group_map is passed by const& (header line 78, cpp line 141) and is never written; the @param[out] on header line 72 should be @param[in]. The only correctly-marked output is im_trafo. Recommended fix is doc-only: change @param[out] transition_group_map to @param[in] (and, while at it, @param[in,out] swath_maps to @param[in], since swath_maps is also const& in correctIM).

### [ANSW-38] DIAChromHandler::extractAndMapChromatogramsForTransitions — "extractAndMap" performs no mapping/filtering and ignores mrm_mapping_param
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TARGETED/DIAChromHandler.h · _targeted-chrom-handlers_

```cpp
std::vector<MSChromatogram> extractAndMapChromatogramsForTransitions(const std::vector<OpenSwath::SwathMap>&, const OpenSwath::LightTargetedExperiment& transition_exp, const ChromExtractParams& cp, const Param& mrm_mapping_param) override
```
- **Expectation:** The base interface IChromatogramHandler documents this method as 'Extract (or select) chromatograms for the given transitions and return mapped & filtered chromatograms' and 'Implementations should return only chromatograms that are mapped and ready for processing'. The name 'AndMap' plus a Param mrm_mapping_param argument lead a caller to expect a mapping/filtering step against the transition list.
- **Actual:** The DIA implementation runs NO MRMMapping pass and does NO per-window or per-transition filtering. Per its own header doc: 'the entire transition_exp is handed to ChromatogramExtractor::prepare_coordinates (no per-window precursor-m/z filtering is performed here)' and 'no separate MRMMapping pass is run ... the result is the concatenation of everything extracted (no per-map filtering down to a "mapped subset" either)'. The mrm_mapping_param argument is 'Currently unused ... kept only to satisfy the IChromatogramHandler interface.'
- **Evidence:** DIAChromHandler.h doc: 'no separate MRMMapping pass is run, and the result is the concatenation of everything extracted (no per-map filtering down to a "mapped subset" either)'; and '@param[in] mrm_mapping_param Currently unused (no MRMMapping pass is run in the DIA handler)'. Contrast IChromatogramHandler.h: 'return mapped & filtered chromatograms' / 'return only chromatograms that are mapped'.
- **Fix:** Either drop 'AndMap' from the DIA-applicable name or document at the interface level that mapping is implicit (native-ID assignment) and not all implementations filter. Cheapest ABI-neutral fix: keep the signature but rename in a follow-up and relax the IChromatogramHandler contract wording so it no longer promises 'mapped & filtered'. Removing the unused mrm_mapping_param would be source-breaking, so keep it but document the inconsistency prominently.
- **Verifier correction:** "AndMap" is defensible: the DIA handler DOES map chromatograms to transitions implicitly, since ChromatogramExtractor::return_chromatogram assigns native IDs from the transition list (DIAChromHandler.cpp:234) and DIA extracts exactly one chromatogram per transition coordinate (so filtering to a "mapped subset" is a vacuous no-op, not a missing step). The accurate, provable surprise is that the mrm_mapping_param argument is silently ignored (void-cast at line 215), so MRMMapping tuning passed by a caller has no effect in the DIA path — while the IChromatogramHandler base contract still promises "mapped & filtered" without caveat. Recoverable (no wrong/lost data; the param has no meaningful DIA job), hence low severity. ABI-neutral fix: relax the interface wording and document the inert param; keep the signature.

### [ANSW-40] IChromatogramHandler::createDefault — createDefault() comment says 'SRM/MRM-based' but returns a runtime SRM/MRM-or-DIA dispatcher
`low` · `stale-documentation` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TARGETED/IChromatogramHandler.h · _targeted-chrom-handlers_

```cpp
static std::unique_ptr<IChromatogramHandler> createDefault()
```
- **Expectation:** Header comment '/// Factory: create the default handler (currently SRM/MRM-based)' tells a caller they get an SRM/MRM handler.
- **Actual:** It returns a DefaultChromHandler, which is a dispatcher that inspects the input maps on every call and routes to either the SRM/MRM or the DIA handler.
- **Evidence:** MRMChromHandler.cpp: 'return std::unique_ptr<IChromatogramHandler>(new DefaultChromHandler());' with comment 'Return the default delegating handler which will pick SRM/MRM vs DIA at runtime'. Header comment contradicts: '(currently SRM/MRM-based)'.
- **Fix:** Update the header comment to 'create the default dispatching handler (SRM/MRM or DIA chosen at runtime)'. Documentation-only; no ABI impact.
- **Verifier correction:** The symbol name createDefault() is not misleading; only the parenthetical header comment "(currently SRM/MRM-based)" is stale. createDefault() returns a DefaultChromHandler, a runtime dispatcher that inspects the swath_maps on each call and delegates to either an internal SRM/MRM (MRMChromHandler) or DIA (DIAChromHandler) handler. The mismatch is comment-only and the behavior is a superset of what the comment implies, so it cannot cause wrong results; severity is low and the fix is documentation-only.

### [ANSW-33] TargetedExperimentHelper::PeptideCompound::getDriftTime — getDriftTime() returns the -1 sentinel for 'unset' but offers no hasDriftTime(), unlike the parallel charge / retention-time accessors which all have presence checks
`low` · `api-inconsistency` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h · _targeted-experiment_

```cpp
double getDriftTime() const
```
- **Expectation:** Consistent with hasCharge()/getChargeState() and hasRetentionTime()/getRetentionTime() in the same class, a caller would expect a hasDriftTime() guard (or a clearly-documented sentinel) before relying on getDriftTime().
- **Actual:** drift_time_ is default-initialized to -1 and getDriftTime() returns it verbatim. There is no hasDriftTime(); setDriftTime() does not set a 'set' flag. A caller reading getDriftTime() on a peptide/compound that never had a drift time set silently receives -1 (a physically meaningless drift time) with no way to distinguish 'unset' from a real value short of magic-number comparison.
- **Evidence:** Lines 225-234 / 287: `void setDriftTime(double dt){ drift_time_ = dt; }  double getDriftTime() const { return drift_time_; }` and `double drift_time_{-1};`. Compare hasCharge()/getChargeState() (lines 212-222) and hasRetentionTime()/getRetentionTime() (lines 239-256).
- **Fix:** Add a bool hasDriftTime() const (drift_time_ >= 0 or a dedicated set flag) and at minimum document the -1 sentinel on getDriftTime(). Additive method, ABI-safe.
- **Verifier correction:** The factual claim holds (no hasDriftTime(), undocumented -1 sentinel, asymmetric vs hasCharge()/hasRetentionTime() in the same class), but the category/severity is overstated as 'silent-failure'. It is an api-inconsistency / discoverability surprise, severity low: the -1 sentinel is a deliberate, codebase-wide convention that every in-tree caller of getDriftTime() already guards with `>= 0.0`/`> 0` (TransitionTSVFile.cpp:1543/1577, MRMFeatureFinderScoring.cpp:572/580, FeatureFinderAlgorithmMetaboIdent.cpp:601-605). A negative drift time is physically impossible, so it self-identifies as 'unset' and fails safe downstream rather than producing silently-wrong quantitation. Recommendation stands and is good hygiene: add an additive `bool hasDriftTime() const { return drift_time_ >= 0; }` and document the -1 sentinel on get/setDriftTime(); this is source-compatible and ABI-safe (no layout/vtable change).

### [ANSW-34] IncludeExcludeTarget (class doc) vs IncludeExcludeTarget::IncludeExcludeTarget — IncludeExcludeTarget copies ReactionMonitoringTransition's class doc claiming default charge = numeric_limits<Int>::max(), but it has no charge member and its default m/z is double::max() while the sibling RMT defaults m/z to 0.0
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h · _targeted-experiment_

```cpp
IncludeExcludeTarget()
```
- **Expectation:** The doc-stated defaults should match the class. A developer reading 'Default values for precursor an product charge is set to numeric_limits<Int>::max()' expects charge accessors and that default; and a developer moving between IncludeExcludeTarget and ReactionMonitoringTransition expects consistent default m/z semantics.
- **Actual:** IncludeExcludeTarget's @brief is a verbatim copy of ReactionMonitoringTransition's ('The default values for precursor and product m/z values are set to numeric_limits<double>::max(). Default values for precursor an product charge is set to numeric_limits<Int>::max().') but IncludeExcludeTarget exposes no charge accessor/member at all. Meanwhile the m/z default differs across the two sibling classes: IncludeExcludeTarget ctor sets precursor_mz_/product_mz_ to numeric_limits<double>::max(), whereas ReactionMonitoringTransition defaults precursor_mz_ to 0.0 (per its header doc and ctor). Callers using max() as an 'unset' marker for one class and 0.0 for the other will write silent bugs.
- **Evidence:** IncludeExcludeTarget.h lines 21-27 doc; IncludeExcludeTarget.cpp lines 15-20: `precursor_mz_(std::numeric_limits<double>::max()), product_mz_(std::numeric_limits<double>::max())`. Contrast ReactionMonitoringTransition.h lines 28-30 ('default values for precursor m/z is 0.0') and ReactionMonitoringTransition.cpp line 26 `precursor_mz_(0.0)`.
- **Fix:** Fix the IncludeExcludeTarget class doc to describe its actual members/defaults (no charge; m/z default double::max()). Doc-only change, ABI-neutral. Aligning the m/z 'unset' sentinel between the two classes would be a behavioral change to weigh separately.
- **Verifier correction:** IncludeExcludeTarget's class @brief (IncludeExcludeTarget.h:24-26) wrongly documents a "precursor an product charge" default of numeric_limits<Int>::max(), but the class has no charge member or accessor — that sentence is copied boilerplate and should be removed. The m/z portion of the doc is correct for this class (ctor sets precursor_mz_/product_mz_ to numeric_limits<double>::max(), IncludeExcludeTarget.cpp:17-18). A separate, real inconsistency exists between sibling classes: ReactionMonitoringTransition uses 0.0 as its uninitialized-m/z sentinel (RMT.cpp:26) while IncludeExcludeTarget uses double::max(); harmonizing that sentinel would be a behavioral change to evaluate independently.

### [ANSW-37] ReactionMonitoringTransition (default library_intensity_) — getLibraryIntensity() returns the magic sentinel -101 when never set, with no hasLibraryIntensity() and no documentation of the sentinel
`low` · `leaky-undocumented-sentinel` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h · _targeted-experiment_

```cpp
double getLibraryIntensity() const
```
- **Expectation:** An intensity getter should return a real (non-negative) intensity, or expose a presence check; a caller would not expect a negative magic number that means 'unset'.
- **Actual:** The default constructor sets library_intensity_(-101) and getLibraryIntensity() returns it unchanged. The header doc ('Returns the library intensity (ion count or normalized ion count from a spectral library)') never mentions that -101 is an unset sentinel, and there is no hasLibraryIntensity() to distinguish unset from a real value. Code that sums or filters by library intensity silently incorporates -101 for unset transitions.
- **Evidence:** ReactionMonitoringTransition.cpp lines 22-35 ctor: `library_intensity_(-101)`; lines 349-352: `double ReactionMonitoringTransition::getLibraryIntensity() const { return library_intensity_; }`. Header lines 149-153 give no hint of the sentinel.
- **Fix:** Document the -101 'unset' sentinel on getLibraryIntensity() and/or add hasLibraryIntensity(). Additive/doc fix, ABI-safe.
- **Verifier correction:** getLibraryIntensity() returns the magic value -101 when never set. This sentinel IS a real, used convention (writers guard with `> -100` in TraMLHandler.cpp:806 and TransitionTSVFile.cpp:1683), but it is undocumented on the getter and not exposed as a named presence check (no hasLibraryIntensity()). The value can propagate unfiltered into the OpenSwath scoring path (DataAccessHelper.cpp:136 → LightTransition.library_intensity), though the main DIA normalization clamps negatives to 0 (DIAPrescoring.cpp:53). Recommend documenting the -101 'unset' sentinel on getLibraryIntensity() and adding hasLibraryIntensity(); additive and ABI-safe.

### [ANSW-49] ReactionMonitoringTransition::setNativeID / getNativeID vs setName / getName — setName/getName and setNativeID/getNativeID are silent aliases for the SAME field; setting one clobbers the other
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h · _targeted-experiment_

```cpp
void setNativeID(const std::string & name); const std::string & getNativeID() const
```
- **Expectation:** A transition has a human-readable 'name' and a separate 'native id' (as everywhere else in OpenMS where getNativeID and getName are distinct concepts). A developer setting setName(...) then reading getNativeID() expects two independent values.
- **Actual:** Both pairs read/write the single member name_. `setName{ name_ = name; }`, `setNativeID{ name_ = name; }`, and both getters return name_. So setNativeID() silently overwrites the value returned by getName() and vice-versa. There is no separate native-id storage. A caller mixing the two (common when porting code that distinguishes them) gets silent data loss.
- **Evidence:** Impl: `void ReactionMonitoringTransition::setNativeID(const std::string & name){ name_ = name; }` and `const std::string & ReactionMonitoringTransition::getNativeID() const { return name_; }`, identical to setName/getName.
- **Fix:** Document explicitly in the header that name and nativeID are the same underlying identifier (alias), so callers don't expect two fields. Changing them to distinct fields would be an ABI break and would alter TraML round-tripping, so prefer the doc fix; flag the aliasing as intentional.
- **Verifier correction:** setName/getName and setNativeID/getNativeID are intentional aliases for a transition's single required identifier (the TraML Transition `id`, member name_ commented as "id"). This is dictated by the data model — a TraML transition has no second "name" field. It is a documentation/POLS gap (two undocumented aliased accessor pairs), not silent data loss: the getters always return equal values, and no caller in the codebase sets both or relies on them being distinct (OpenSwath uses getNativeID exclusively; getName is unused for transitions). Recommended fix is a header doc note that name and nativeID are the same identifier (ABI: none). Severity low, not high.
