# OpenMS POLS Audit — Batch 5: FEATUREFINDER + PROCESSING

**Confirmed findings:** 89 (FEATUREFINDER 53, PROCESSING 36) · 9 high · 35 medium · 45 low.

### [FEAT-28] SeedListGenerator::generateSeedList — generateSeedList dereferences a possibly-end precursor iterator and indexes precursors[0] without bounds checks
`high` · `silent-failure` · ABI: `none` · src/openms/source/FEATUREFINDER/SeedListGenerator.cpp · _ff-identification_

```cpp
void generateSeedList(const PeakMap& experiment, SeedList& seeds)
```
- **Expectation:** A 'generate a seed list from MS2 precursors' helper should robustly skip MS2 spectra that lack a usable precursor or a preceding MS1 spectrum, not crash on them.
- **Actual:** For every MS2 spectrum it calls experiment.getPrecursorSpectrum(exp_it) and immediately dereferences prec_it->getRT(), and reads precursors[0].getMZ(). getPrecursorSpectrum can legitimately return experiment.end() (e.g. an MS2 at the start of the run, or no preceding MS1), and getPrecursors() can be empty for an MS2 with no recorded precursor. Both cases are undefined behavior / out-of-bounds access.
- **Evidence:** PeakMap::ConstIterator prec_it = experiment.getPrecursorSpectrum(exp_it); const vector<Precursor>& precursors = exp_it->getPrecursors(); DPosition<2> point(prec_it->getRT(), precursors[0].getMZ()); — and MSExperiment::getPrecursorSpectrum returns 'spectra_.end()' in several branches (iterator==begin, ms_level==1, none found).
- **Fix:** Guard the loop body: skip the spectrum when precursors.empty() or prec_it == experiment.end(). Additive, no ABI change (implementation-only fix).
- **Verifier correction:** Category "silent-failure" is partly imprecise: Path 1 (prec_it->getRT() on an end() iterator for an MS2 at run start / no preceding MS1) is more accurately undefined-behavior / likely hard crash than a silent failure, while Path 2 (precursors[0] on an empty vector) is the silent out-of-bounds read that can produce garbage seed positions without crashing. Both are real and reachable. The recommendation is correct: guard the loop body with `if (precursors.empty() || prec_it == experiment.end()) continue;`. This is implementation-only and additive, so ABI impact is none.

### [FEAT-35] Biosaur2Algorithm::run — run() destructively mutates the stored MS data set via setMSData(); getMSData() no longer returns what was set
`high` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/FEATUREFINDER/Biosaur2Algorithm.h · _ff-metabo_

```cpp
void run(FeatureMap& feature_map, std::vector<Hill>& hills, std::vector<PeptideFeature>& peptide_features)
```
- **Expectation:** After setMSData(exp) and run(...), getMSData() returns the experiment that was set, and run() can be called again on the same stored data.
- **Actual:** run() erases all non-MS1 spectra from the internal ms_data_, and may centroid (PeakPickerHiRes) and apply TOF filtering in place. The stored experiment is consumed/altered; a second run() or a subsequent getMSData() sees the reduced/centroided data, not the original.
- **Evidence:** Biosaur2Algorithm.cpp:204-217 erases `ms_data_` spectra with MSLevel!=1 and calls `centroidProfileSpectra_(ms_data_)` / `processTOF_(ms_data_)` directly on the member. The header notes it (line 189-190) but the destructive, non-idempotent nature is easy to miss given a getMSData() that returns a mutable reference.
- **Fix:** Operate on a local copy of ms_data_ inside run() so the stored data is preserved and run() is idempotent (behavior change, but source-compatible). Minimally, strengthen the header note to say the stored data is irreversibly modified and run() is not re-entrant on the same instance.
- **Verifier correction:** Accurate but understated: in addition to MS1 erasure, profile centroiding, and TOF filtering on ms_data_, the FAIMS path (Biosaur2Algorithm.cpp:284) std::move's ms_data_ into IMDataConverter::splitByFAIMSCV, leaving the stored experiment empty (moved-from). Thus after a FAIMS run() getMSData() returns an emptied experiment, and after any run() a second run() sees reduced/centroided (or empty) data. The header @notes (189-190) cover only MS1 filtering and profile centroiding; the move-to-empty and non-idempotency are undocumented.

### [FEAT-38] ElutionPeakDetection::filterByPeakWidth — filterByPeakWidth unconditionally prints to std::cout and indexes the result vector without an empty check
`high` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/ElutionPeakDetection.h · _ff-metabo_

```cpp
void filterByPeakWidth(std::vector<MassTrace>&, std::vector<MassTrace>&)
```
- **Expectation:** A library filter quietly writes the surviving traces to the output vector and returns; it does not print to stdout and does not crash on edge inputs.
- **Actual:** It always writes a 'pw low/pw high' line to std::cout, and then indexes filt_mtraces[0] and filt_mtraces[filt_mtraces.size()-1]. If the input is empty (or all traces are filtered out), filt_mtraces is empty and these accesses are out-of-bounds (UB/crash).
- **Evidence:** ElutionPeakDetection.cpp:375 `std::cout << "pw low: " << filt_mtraces[0].estimateFWHM(true) << ... << filt_mtraces[filt_mtraces.size() - 1].estimateFWHM(true) << '\n';` executed unconditionally after the filtering loop.
- **Fix:** Remove the std::cout diagnostic (or gate it behind a debug macro) and guard the indexing with `if (!filt_mtraces.empty())`. Both are additive and ABI-safe.

### [PROC-23] estimateNoiseFromRandomScans — Random-scan picker indexes exp[] directly with a random number, not via the filtered scan-index list, so it can sample wrong-MS-level / empty spectra (and overrun exp)
`high` · `silent-failure` · ABI: `none` · src/openms/source/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.cpp · _proc-calibration-noise_

```cpp
float estimateNoiseFromRandomScans(const MSExperiment& exp, const UInt ms_level, const UInt n_scans = 10, const double percentile = 80)
```
- **Expectation:** A caller reading 'estimateNoiseFromRandomScans(exp, ms_level, ...)' expects it to draw n_scans random spectra OF THE REQUESTED ms_level (it built spec_indices precisely for that) and average their percentile intensity.
- **Actual:** It builds spec_indices of valid same-level scans, then computes `UInt scan = distribution * (spec_indices.size()-1)` and dereferences `exp[scan]` -- the random value indexes the FULL experiment, NOT spec_indices. So it samples arbitrary spectra (wrong MS level, empty, or chromatograms) and `scan` is bounded by spec_indices.size(), unrelated to exp.size(). If a sampled spectrum is empty, `tmp[idx]` reads tmp[0] of an empty vector (UB). The whole ms_level filtering is silently discarded.
- **Evidence:** vector<Size> spec_indices; ... if (exp[i].getMSLevel()==ms_level && !exp[i].empty()) spec_indices.push_back(i); ... UInt scan = (UInt)(distribution(generator) * (spec_indices.size() - 1)); ... for (const auto& peak : exp[scan]) ...  -- note exp[scan], never spec_indices[scan].
- **Fix:** Index through the filter list: `UInt scan = spec_indices[(Size)(distribution(generator) * (spec_indices.size()-1))];` then `exp[scan]`. This is a pure source-compatible bug fix (no signature change). Also clamp idx for tiny spectra.
- **Verifier correction:** The bug is real but one detail is wrong: exp[scan] does NOT overrun exp, because scan is bounded by spec_indices.size()-1 and spec_indices.size() <= exp.size(). The actual out-of-bounds/UB is tmp[idx] when a sampled spectrum is empty (empty-vector read). The core defect stands: exp[scan] must be exp[spec_indices[scan]]; the ms_level and non-empty filtering is silently bypassed, yielding wrong noise estimates for any non-trivial (e.g. interleaved MS1/MS2) experiment. Fix is internal (function body only), so ABI impact is none.

### [PROC-15] NLargest::filterSpectrum — Filtering reorders the spectrum to descending-intensity order, not just removing peaks
`high` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/FILTERING/NLargest.h · _proc-filtering-scaling_

```cpp
template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
```
- **Expectation:** A caller of an 'N largest peaks' filter expects the surviving peaks to remain in the spectrum's natural (m/z / position) order; only count changes.
- **Actual:** filterSpectrum calls spectrum.sortByIntensity(true) and then spectrum.select(indices) with indices 0..n-1. MSSpectrum::select preserves the order of the given indices (confirmed in MSSpectrum.cpp), so the spectrum is left sorted by DESCENDING intensity, not by m/z. The peaks are never re-sorted back by position.
- **Evidence:** NLargest.h: 'spectrum.sortByIntensity(true);' then 'for (Size i = 0; i != peakcount_; ++i) { indices.push_back(i); } spectrum.select(indices);'. MSSpectrum::select builds tmp in indices order and swaps it in, so output order == intensity-sorted order.
- **Fix:** Re-sort the spectrum by position after select() (e.g. spectrum.sortByPosition()), or clearly document in the header that the spectrum is left in descending-intensity order. ABI-safe: it is an implementation-body change in an inline template. If callers depend on the current order, add a bool keep_order param with a default preserving today's behavior.

### [PROC-16] RankScaler::filterSpectrum — RankScaler leaves the spectrum sorted by intensity (not m/z) after scaling
`high` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/SCALING/RankScaler.h · _proc-filtering-scaling_

```cpp
template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
```
- **Expectation:** A 'scaler' rewrites intensities in place; a caller expects peak ordering (by m/z) to be unchanged.
- **Actual:** filterSpectrum calls spectrum.sortByIntensity() and never restores positional order, so on return the spectrum is sorted by intensity. Any downstream code that assumes m/z-sorted peaks (e.g. for alignment/scoring) will silently operate on a reordered spectrum.
- **Evidence:** RankScaler.h: 'spectrum.sortByIntensity(); ... it->setIntensity(count); ... while (it != spectrum.begin());' with no subsequent sortByPosition().
- **Fix:** Call spectrum.sortByPosition() at the end, or document that the spectrum is returned in intensity-sorted order. ABI-safe inline-body change.

### [PROC-31] FeatureOverlapFilter::mergeFAIMSFeatures — mergeFAIMSFeatures silently wipes protein IDs, peptide IDs, data processing and meta info via feature_map.clear()
`high` · `hidden-side-effect` · ABI: `none` · src/openms/source/PROCESSING/FEATURE/FeatureOverlapFilter.cpp · _proc-misc_

```cpp
static void mergeFAIMSFeatures(FeatureMap& feature_map, double max_rt_diff = 5.0, double max_mz_diff = 0.05)
```
- **Expectation:** A function documented as 'Merge FAIMS features ... (modified in place)' rebuilds the feature list but otherwise preserves the FeatureMap's attached metadata (protein/peptide identifications, document/unique IDs, ranges, data processing, MetaInfo).
- **Actual:** Line 515 calls feature_map.clear() with the default argument. FeatureMap::clear(bool clear_meta_data = true) (FeatureMap.cpp:462) clears not only data_ but also clearMetaInfo(), clearRanges(), the DocumentIdentifier, clearUniqueId(), protein_identifications_, unassigned_peptide_identifications_, data_processing_ and id_data_. So all of this metadata is silently destroyed whenever any FAIMS features are present, while the same call is a no-op (returns early) when no FAIMS features exist - giving inconsistent, data-loss behavior.
- **Evidence:** FeatureOverlapFilter.cpp:515 `feature_map.clear();` followed by repopulating only the feature list; FeatureMap.cpp:462-477 shows clear() default-clears protein_identifications_, unassigned_peptide_identifications_, data_processing_, id_data_, meta info, ranges and document/unique IDs.
- **Fix:** Call feature_map.clear(false) to keep attached metadata, or repopulate data in place via swap of the feature container only. This is a bug fix that is source- and ABI-compatible (no signature change). Also document that unassigned peptide IDs are not carried through the merge.

### [PROC-32] IDFilter::filterHitsByScore(std::vector<IdentificationType>&, double, IDScoreSwitcherAlgorithm::ScoreType) — Doc and warning claim hits are removed when score_type is missing, but hits are silently kept
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/ID/IDFilter.h · _proc-misc_

```cpp
template<class IdentificationType> static void filterHitsByScore(std::vector<IdentificationType>& ids, double threshold_score, IDScoreSwitcherAlgorithm::ScoreType score_type)
```
- **Expectation:** The Doxygen note states 'Removes a hit if the @p score_type is not found at all' and the emitted warning says 'No hit with the given score_type found. All hits removed.' A caller therefore expects that IDs lacking the requested score type have all their hits removed.
- **Actual:** When neither the main score nor any secondary score of an ID matches score_type, the else-branch finds an empty score_name and does nothing - keepMatchingItems is never called, so every hit of that ID is retained. The 'All hits removed' warning only fires if NOT A SINGLE id in the whole vector matched, and even then no hits are actually removed. The actual behavior is the opposite of both the note and the warning.
- **Evidence:** IDFilter.h:873 note 'Removes a hit if the @p score_type is not found at all'; lines 892-911 the else branch only filters when result.score_name is non-empty; line 913 warning 'All hits removed' although no removal occurs.
- **Fix:** Decide on one behavior and make code, note and warning agree: either clear hits of IDs whose score_type is absent (matching the documented contract) or change the note/warning to say such hits are left untouched. Body-only change; source- and ABI-compatible.

### [PROC-33] MorphologicalFilter::filterRange / applyErosion_ / applyDilation_ — Public filter uses function-local static buffers, making filter()/filterExperiment() non-reentrant and not thread-safe
`high` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h · _proc-misc_

```cpp
template<typename InputIterator, typename OutputIterator> void filterRange(InputIterator input_begin, InputIterator input_end, OutputIterator output_begin)
```
- **Expectation:** A non-static member filter method taking explicit input/output ranges is expected to be reentrant and safe to call concurrently on independent MorphologicalFilter instances (a common pattern when parallelizing per-spectrum filtering).
- **Actual:** filterRange uses `static std::vector<...> buffer;` (line 171) and applyErosion_/applyDilation_ each use their own `static std::vector<ValueType> buffer;` (lines 330, 440). These statics are shared across all instances and threads, so concurrent calls corrupt each other's intermediate results. filter() and filterExperiment() route through these methods, so the whole public API is silently unsafe to run in parallel despite nothing in the signature suggesting shared state.
- **Evidence:** MorphologicalFilter.h:171 `static std::vector<typename InputIterator::value_type> buffer;`; :330 and :440 `static std::vector<ValueType> buffer;`.
- **Fix:** Make the scratch buffers per-instance members (or stack-local) instead of static. Note this changes templated inline code that is recompiled by clients, so it is source-compatible but technically a header/inline-definition change; gate behind a minor version. Document thread-safety expectations in the class docs regardless.

### [FEAT-19] FeatureFinderAlgorithmMetaboIdent::getMSData — getMSData() returns an EMPTY map after a successful run(), contradicting its 'empty if run was not executed' doc
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h · _ff-algorithms_

```cpp
PeakMap& getMSData(); const PeakMap& getMSData() const;
```
- **Expectation:** The doc comment '/// @brief Retrieve chromatograms (empty if run was not executed)' tells a caller that getMSData() is empty ONLY when run() did not execute, so after a successful run() it should hold the (processed) MS data.
- **Actual:** run() moves ms_data_ out (splitByFAIMSCV(std::move(ms_data_))) and runSingleGroup_() explicitly calls ms_data_.reset() at the end ('// not needed anymore, free up the memory'). After a successful run(), getMSData() therefore returns an EMPTY map -- the exact opposite of what the doc implies. Additionally the doc text says 'Retrieve chromatograms' although this accessor returns the MS spectra (ms_data_), not the chromatograms (that is getChromatograms()).
- **Evidence:** Header: '/// @brief Retrieve chromatograms (empty if run was not executed)\n  PeakMap& getMSData() { return ms_data_; }'. Source FeatureFinderAlgorithmMetaboIdent.cpp:408: 'ms_data_.reset(); // not needed anymore, free up the memory'; run() at :173 'splitByFAIMSCV(std::move(ms_data_))'.
- **Fix:** Fix the doc to state that getMSData() is empty BOTH before run() and after run() (the data is consumed/freed). If callers need the post-processing data, do not reset ms_data_ in run(), or add an explicit note. Doc-only change is source/ABI-compatible; keeping the data would be a behavioral change.

### [FEAT-21] FeatureFinderAlgorithmMetaboIdent::run — run() has no MS-data parameter; the required spectra must be injected via setMSData() first or run() silently returns unchanged
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h · _ff-algorithms_

```cpp
void run(const std::vector<FeatureFinderMetaboIdentCompound>& metaboIdentTable, FeatureMap& features, const std::string& spectra_file = "")
```
- **Expectation:** A feature-finder run(table, features, spectra_file) that takes a 'spectra_file' string and an output FeatureMap looks self-contained -- a caller reasonably expects it to load/use the spectra it needs.
- **Actual:** spectra_file is only a fall-back primaryMSRunPath annotation, NOT an input loader. The actual LC-MS input must be supplied beforehand via setMSData(); run() consumes ms_data_ via std::move. If setMSData() was never called (or the map has no MS1 scans), run() returns 'features' unchanged with no error -- a silent no-op. The hidden ordering dependency (setMSData before run) is not visible in run()'s signature.
- **Evidence:** Header doc: 'If @p spectra_file is provided it will be used as a fall-back to setPrimaryMSRunPath ... If there are no MS1 scans in the MSData return @p features unchanged.' Source run() at :173 moves ms_data_; runSingleGroup_ at :276 'if (ms_data_.empty())' returns early.
- **Fix:** Document the mandatory setMSData()-before-run() contract on run() itself, and consider signalling the empty-input case (log already exists at runtime). Doc-only change is ABI-safe; adding an overload run(PeakMap&&, ...) would be additive.
- **Verifier correction:** run(table, features, spectra_file) does NOT load spectra_file; the LC-MS input is the member ms_data_ that must be set via setMSData() before run() (run() consumes it via std::move at :173). spectra_file is used only as a fallback primaryMSRunPath annotation (runSingleGroup_ :274). If setMSData() was not called or the map has no MS1 scans, run() emits a log-only WARN and returns 'features' unchanged (:276-280) -- recoverable and partly documented in the header doc on run() (lines 148-151) and via the WARN log, but the mandatory setMSData()-before-run() precondition is not stated in run()'s signature or doc. Fix is doc-only (ABI-safe).

### [FEAT-23] FeatureFinderAlgorithmPickedHelperStructs::MassTraces::computeIntensityProfile — computeIntensityProfile() dereferences begin() unconditionally -> UB on an empty MassTraces, while sibling queries throw a documented Precondition
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h · _ff-algorithms_

```cpp
void computeIntensityProfile(std::list<std::pair<double, double> >& intensity_profile) const
```
- **Expectation:** Sibling const queries on MassTraces document and enforce the empty case: getTheoreticalmaxPosition() and getRTBounds() both '@exception Exception::Precondition is thrown if there are no mass traces'. A caller reasonably expects computeIntensityProfile() to be equally safe (throw or no-op) on empty input.
- **Actual:** computeIntensityProfile() begins with 'trace_it = this->begin();' and immediately iterates 'trace_it->peaks' with no empty check. On an empty MassTraces this dereferences end() -> undefined behavior (crash), instead of the documented Precondition exception its siblings raise. The header gives no precondition.
- **Evidence:** Source :196-201: 'TTraceIterator trace_it = this->begin(); ... for (TTracePeakIterator trace_peak_it = trace_it->peaks.begin(); ...'. Compare getRTBounds() :162 'if (this->empty()) throw Exception::Precondition(...)'.
- **Fix:** Add an empty-guard at the top of computeIntensityProfile() (early return on empty, or throw Exception::Precondition to match siblings) and document it in the header. Behavior/doc change; ABI-safe.

### [FEAT-27] FeatureFinderAlgorithmMetaboIdent::FeatureFinderMetaboIdentCompound — Constructing a compound that violates the documented contract is silently dropped at run() time, not reported to the caller
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h · _ff-algorithms_

```cpp
FeatureFinderMetaboIdentCompound(const std::string& _name, const std::string& _formula, double _mass, const std::vector<int>& _charges, const std::vector<double>& _rts, const std::vector<double>& _rt_ranges, const std::vector<double>& _iso_distrib, const std::vector<double>& _ion_mobilities = {}, const std::string& _adduct = "")
```
- **Expectation:** Either the constructor/accessors validate the contract (mass>0 or non-empty formula, non-empty rts, non-zero charges, rt_ranges/ion_mobilities sized 1 or rts.size()), or run() reports rejected compounds so the caller knows their target was not searched.
- **Actual:** The struct is a plain container that performs no validation (per its own doc). Compounds violating the contract are 'silently dropped with an OPENMS_LOG_ERROR message at run time'. A programmatic caller building the table gets no return value, exception, or count of dropped targets from run() -- the target simply does not appear in the results, which is easy to miss.
- **Evidence:** Header doc on the struct: 'No validation is performed in the constructor or accessors ... Compounds that violate the contract are silently dropped with an @c OPENMS_LOG_ERROR message at run time.'
- **Fix:** Provide feedback to the caller (e.g. return/collect the count or list of dropped targets from run(), or add a validate() method on the compound). Adding a getter/return path is additive and ABI-safe; the current behavior should at least be cross-referenced from run()'s doc.
- **Verifier correction:** The behavior is real and accurately described, but it is not strictly "silent": each dropped target/charge emits an OPENMS_LOG_ERROR ("... - skipping this target."). The surprise is that there is no PROGRAMMATIC feedback path — run()/runSingleGroup_/addTargetToLibrary_ all return void, no exception is thrown, and no dropped-target count or list is exposed (getNShared() reports shared identifications, not drops). A library/pyOpenMS caller assembling the table in memory must scrape logs to learn a target was not searched. Severity downgraded from high to medium accordingly: data is effectively missing from results, but the failure is loud (LOG_ERROR) and the input remains intact/recoverable. The recommended fix (collect/return a dropped-target count or list from run(), or add a validate() method on the compound) is additive and source-compatible; changing run()'s return type instead would be ABI-breaking, so the additive route is preferred.

### [FEAT-8] EmgModel::getCenter — EmgModel::getCenter() returns statistics_.mean(), not the EMG peak position (which is governed by retention_)
`medium` · `misleading-name` · ABI: `none` · src/openms/source/FEATUREFINDER/EmgModel.cpp · _ff-fitters_

```cpp
EmgModel::CoordinateType EmgModel::getCenter() const
```
- **Expectation:** getCenter() on an exponentially-modified-Gaussian elution model returns the position of the peak maximum, as the base-class doc for getCenter promises ('position of the maximum'). For an EMG the apex is determined by the height/width/symmetry/retention parameters, not by the Gaussian 'mean'.
- **Actual:** It returns statistics_.mean(). But the sampled curve in setSamples() is computed entirely from height_, width_, symmetry_ and retention_ (the 'emg:*' params); statistics_.mean() ('statistics:mean') never enters the shape and is left at its default/fit value. The header doc even copy-pastes 'get the center of the Gaussian model i.e. the position of the maximum', which is wrong for an EMG. The returned value is unrelated to where the modelled peak actually peaks.
- **Evidence:** Header: '/// get the center of the Gaussian model i.e. the position of the maximum\n CoordinateType getCenter() const override;'. Cpp getCenter(): 'return statistics_.mean();'. setSamples() builds data from 'height_ * width_ / symmetry_ ... exp(... diff = pos - retention_ ...)' with no use of statistics_.mean().
- **Fix:** Either compute and return the actual EMG apex (e.g. argmax of the sampled interpolation, or the analytic mode) or, if statistics_.mean() is intentionally the reference center, fix the doc and rename intent. Returning the apex is the behaviour callers expect from getCenter on a skewed peak. ABI: changing the return value is source-compatible (same signature); a doc-only fix is none.
- **Verifier correction:** getCenter() does return statistics_.mean() and that value is disconnected from the EMG apex (governed by emg:retention/width/symmetry) — confirmed. But the misleading "position of the maximum" promise lives in EmgModel's own header (EmgModel.h:49), copy-pasted from GaussModel; the base InterpolationModel::getCenter doc only says the center's definition is model-particular. Additionally, through the EmgFitter1D production path statistics:mean is never set to a peak value (stays at the Fitter1D default 1.0), making the returned center even more clearly meaningless than "left at its fit value." abi_impact: a doc-only fix is none; changing the return to the true apex would be source-compatible (same signature, changed value).

### [FEAT-9] EmgModel::setSamples — EmgModel::setSamples() does not normalize to unit integral and ignores scaling_, unlike GaussModel/BiGaussModel
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/source/FEATUREFINDER/EmgModel.cpp · _ff-fitters_

```cpp
void EmgModel::setSamples() override
```
- **Expectation:** Sibling InterpolationModel subclasses normalize their sampled curve and apply the intensity scaling factor: GaussModel and BiGaussModel both compute 'factor = scaling_ / interpolation_step_ / sum' and multiply every sample by it, so setScalingFactor() controls the model height and the discrete integral is well-defined. A caller would expect EmgModel to behave the same.
- **Actual:** EmgModel::setSamples() pushes raw EMG values and applies NO normalization factor and NEVER reads scaling_. As a result setScalingFactor() (inherited, documented public API that writes 'intensity_scaling') has no effect on an EmgModel, and the curve magnitude is whatever height_/width_/symmetry_ produce. This silently diverges from the Gauss/BiGauss contract in the same module.
- **Evidence:** GaussModel::setSamples: 'IntensityType factor = scaling_ / interpolation_step_ / std::accumulate(...); for (auto& value : data) value *= factor;'. EmgModel::setSamples: the loop ends with 'data.push_back((part1 * sqrt_2pi * exp(...) / (1 + exp(...))));' followed directly by 'interpolation_.setScale(...); interpolation_.setOffset(min_);' with no factor and no reference to scaling_.
- **Fix:** Document explicitly that EmgModel is non-normalized and that setScalingFactor() is a no-op for it, or apply the same scaling_/normalization the sibling models use. At minimum the silent no-op of the inherited setScalingFactor() on this subclass should be called out. ABI: additive/doc fix is none; applying scaling changes numeric output (source-compatible).
- **Verifier correction:** The factual claim is correct and provable. Severity should be medium, not high: setScalingFactor() being a silent no-op on EmgModel (while functional on GaussModel/BiGaussModel) is a Principle-of-Least-Surprise violation for callers using the documented InterpolationModel base API, but it does not produce silently-wrong results in normal use — no in-tree caller sets scaling on an EmgModel, and EMG encodes its height in emg:height which the fitter (EmgFitter1D) optimizes against the raw data. abi_impact is 'none' for the recommended doc/contract fix; applying scaling_ to match the siblings would instead be source-compatible (changes numeric output only). Recommended fix: either document that EmgModel is non-normalized and setScalingFactor() has no effect on it, or read scaling_ and normalize like the sibling models.

### [FEAT-12] EmgScoring::elutionModelFit / calcElutionFitScore — Fit-failure sentinel -1 collides with a legitimate worst-case Pearson correlation score
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/source/FEATUREFINDER/EmgScoring.cpp · _ff-fitters_

```cpp
double elutionModelFit(const ConvexHull2D::PointArrayType& current_section, bool smooth_data) const
```
- **Expectation:** A scoring function returning a Pearson-correlation-style fit quality should signal 'could not fit' distinguishably from 'fit succeeded but correlated -1'. A caller cannot otherwise tell an aborted fit from a genuinely anti-correlated profile.
- **Actual:** elutionModelFit() returns -1 when there are fewer than 2 points ('if (current_section.size() < 2) return -1;'), and the underlying EmgFitter1D::fit1d() also returns -1.0 on NaN correlation. -1 is simultaneously the lowest valid Pearson value, so calcElutionFitScore() averages real and sentinel -1 values together and the caller has no way to detect failure. The code's own comment flags this ('Currently -1 is just the "lowest" pearson correlation to a fit that you can have.').
- **Evidence:** EmgScoring.h elutionModelFit: 'if (current_section.size() < 2) { return -1; }' and TODO comment '//TODO think about penalizing aborted fits even more. Currently -1 is just the "lowest" pearson correlation ...'. EmgFitter1D.cpp: 'if (std::isnan(correlation)) { correlation = -1.0; }'.
- **Fix:** Use a distinct out-of-range sentinel (e.g. NaN, or a documented value < -1) for the failure/too-few-points path, or return a status alongside the score, and document the contract. Behaviour/return-value change only; source-compatible.

### [FEAT-14] EmgFitter1D::symmetric_ — fit1d() silently skips Levenberg-Marquardt optimization when an internal symmetric_ flag is set by initial-parameter estimation
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/EmgFitter1D.h · _ff-fitters_

```cpp
if (!symmetric_) { EgmFitterFunctor functor(4, &d); optimize_(x_init, functor); }
```
- **Expectation:** fit1d() is the public entry point that fits the EMG to the data via Levenberg-Marquardt (per the class brief). A caller expects the LM optimization to run and the returned model to be an actual fit.
- **Actual:** Whether LM runs at all is gated on an internal 'symmetric_' flag that setInitialParameters_()/setInitialParametersMOM_() flip to true when the moment/symmetry estimate is inf/NaN. In that case fit1d() skips optimize_() entirely and just builds an EmgModel from the rough initial guess (symmetry_ forced to 10.0), then returns a Pearson correlation as if a fit happened. This data-dependent 'no actual fit' branch is invisible from the signature and the returned quality does not flag it.
- **Evidence:** EmgFitter1D.cpp fit1d: 'if (!symmetric_) { EgmFitterFunctor functor(4, &d); optimize_(x_init, functor); }' (no else). setInitialParameters_: 'if (std::isinf(symmetry_) || std::isnan(symmetry_)) { symmetric_ = true; symmetry_ = 10.0; }'. symmetric_ is an inherited internal member, not a documented public knob.
- **Fix:** Document this fallback in fit1d()'s contract (degenerate input -> no LM optimization, returns correlation of the initial-guess model), and consider signalling it via the returned quality or a log warning. Doc/behaviour clarification; source-compatible.
- **Verifier correction:** The mechanism is correctly described. Severity is medium, not high: the no-LM branch is reached only on degenerate input where the symmetry estimate is inf/NaN (e.g. zero denominator when the weighted/intensity median position equals the first sample position). It does not produce silently wrong results in the corrupting sense — the EmgModel built from the initial guess is valid, and the returned Pearson correlation genuinely measures the (likely poor) model-vs-data agreement, so a caller checking the quality score can detect a bad fit. The surprise is that LM is silently skipped (no else, no warning, no dedicated flag) despite fit1d's brief promising Levenberg-Marquardt optimization. Recommendation stands: document the degenerate-input fallback in fit1d's contract and consider a log warning or an explicit quality signal; change is doc/behavior-clarification only and source-compatible.

### [FEAT-30] FeatureFinderIdentificationAlgorithm::setMSData — setMSData silently discards all non-MS1 spectra from the stored data
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h · _ff-identification_

```cpp
void setMSData(const PeakMap& ms_data); void setMSData(PeakMap&& ms_data)
```
- **Expectation:** setMSData / getMSData is a simple cached-data setter/getter pair; a caller setting data and reading it back expects the same PeakMap.
- **Actual:** Both overloads erase every spectrum whose MS level != 1 from the internal copy, so getMSData() returns a strictly smaller PeakMap than was set. The const-ref overload doc only says 'copied into the internal cache'; the MS1-only filtering is not mentioned (only the move overload vaguely notes getMSData gives a 'processed/modified version').
- **Evidence:** ms_data_ = ms_data; vector<MSSpectrum>& specs = ms_data_.getSpectra(); specs.erase(std::remove_if(specs.begin(), specs.end(), [](const MSSpectrum& s){ return s.getMSLevel() != 1; }), specs.end());
- **Fix:** Document the MS1-only filtering explicitly in both setMSData overload briefs (the algorithm only operates on MS1, so the behavior is intentional but surprising). Doc-only change, no ABI impact.

### [FEAT-31] Internal::FFIDAlgoExternalIDHandler::addExternalPeptide — addExternalPeptide stores a raw pointer to the caller's peptide and truncates its hits, with no lifetime warning
`medium` · `ownership-lifetime` · ABI: `none` · src/openms/source/FEATUREFINDER/FFIDAlgoExternalIDHandler.cpp · _ff-identification_

```cpp
void addExternalPeptide(PeptideIdentification& peptide)
```
- **Expectation:** 'Add an external peptide to the handler's map' reads like the handler copies/records the peptide's data; the input would not be expected to be mutated nor to need to outlive the call.
- **Actual:** It calls peptide.sort(), then peptide.getHits().resize(1) (destroying all but the best hit of the caller's object), and stores &peptide (a raw pointer into the caller's container) in external_peptide_map_. The handler later dereferences these pointers (e.g. in fillExternalRTMap_/dummy-ID creation). If the caller's PeptideIdentificationList is reallocated or destroyed, the stored pointers dangle. The sibling private method addPeptideToMap_ in the main algorithm documents exactly this hazard ('stores a pointer ... Make sure it stays valid until destruction'), but addExternalPeptide does not.
- **Evidence:** peptide.sort(); PeptideHit& hit = peptide.getHits()[0]; peptide.getHits().resize(1); ... external_peptide_map_[hit.getSequence()][charge].emplace(rt, &peptide);
- **Fix:** Document both the in-place mutation (sort + truncate to best hit) and the pointer-retention/lifetime contract in the header, mirroring the existing addPeptideToMap_ caveat. Doc-only; no ABI change.
- **Verifier correction:** Core facts confirmed. One emphasis correction: the headline danger is not primarily a crash — the consuming code paths (fillExternalRTMap_, addDummyPeptideID_) that dereference the stored pointers are not currently called by any live code path, so the dangling-pointer deref is latent. The guaranteed, immediately-realizable surprise is the silent in-place mutation: addExternalPeptide unconditionally sorts and truncates the caller's PeptideIdentification to a single hit. The fix remains doc-only (mirror addPeptideToMap_'s CAUTION comment for both the mutation and the pointer-retention contract on addExternalPeptide and addExternalPeptideToMap_); no signature/ABI change.

### [FEAT-16] IsotopeModel::setSamples(const EmpiricalFormula&) vs InterpolationModel::setSamples() — Inherited no-arg setSamples() throws NotImplemented while the model also exposes a working setSamples(formula) overload
`medium` · `api-asymmetry` · ABI: `source-compatible` · src/openms/include/OpenMS/FEATUREFINDER/IsotopeModel.h · _ff-isotope-models_

```cpp
virtual void setSamples(const EmpiricalFormula & formula); using InterpolationModel::setSamples;
```
- **Expectation:** Given `using InterpolationModel::setSamples;` brings the no-arg setSamples() into scope, and the base class describes setSamples() as 'set sample/supporting points of interpolation wrt params', a caller would expect model.setSamples() (the parameter-driven form siblings use) to build the IsotopeModel from its configured params.
- **Actual:** IsotopeModel does NOT override the no-arg setSamples(); the inherited InterpolationModel::setSamples() body is `throw Exception::NotImplemented(...)`. The model can only be built via the formula-taking overload setSamples(formula), which is asymmetric with every other InterpolationModel subclass (GaussModel/EmgModel/BiGaussModel all override the no-arg setSamples()). A caller who follows the sibling pattern and calls setSamples() with no args gets a runtime throw instead of a built model.
- **Evidence:** IsotopeModel.h:71-74 declares `virtual void setSamples(const EmpiricalFormula & formula);` then `using InterpolationModel::setSamples;`. InterpolationModel.h:136-139 `virtual void setSamples() { throw Exception::NotImplemented(...); }`. IsotopeModel does not redeclare the no-arg form, so it stays the throwing base version. IsotopeFitter1D.cpp:104 must explicitly call the formula overload to work around this.
- **Fix:** Either override the no-arg setSamples() in IsotopeModel to call setSamples(getFormula()) (restoring the sibling contract), or document clearly on the class/overload that only setSamples(formula) builds the IsotopeModel and the inherited no-arg form throws. Adding an override is source-compatible and ABI-additive.
- **Verifier correction:** The defect is real but the category is mis-stated: it is NOT a silent failure. A caller who follows the sibling contract (configure params, then call setSamples() with no args) gets a loud runtime throw of Exception::NotImplemented from the inherited InterpolationModel::setSamples(), because IsotopeModel never overrides the no-arg form and its updateMembers_() (unlike every sibling) does not rebuild the interpolation. The model can only be built via setSamples(const EmpiricalFormula&), which IsotopeFitter1D.cpp:104 explicitly works around. Correct category is api-asymmetry / surprising-but-loud, severity medium (recoverable, throws rather than producing wrong data). Recommended fix (override no-arg setSamples() to call setSamples(getFormula()), or document the asymmetry) is source-compatible and ABI-additive — the virtual slot already exists in the base vtable, so no layout change.

### [FEAT-17] ElutionModelFitter::fitElutionModels — fitElutionModels() can clear the entire FeatureMap and overwrites each feature's intensity in place
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/source/FEATUREFINDER/ElutionModelFitter.cpp · _ff-isotope-models_

```cpp
void fitElutionModels(FeatureMap& features)
```
- **Expectation:** A method named 'fit elution models to all features' is expected to fit models and annotate features (e.g. add meta values / model parameters). The verb 'fit' does not suggest it may delete all features or replace their primary intensity values.
- **Actual:** fitElutionModels() (1) overwrites each feature's intensity with the model area (or 0.0 / an imputed value) via setIntensity (ElutionModelFitter.cpp:414, 428, 447), destroying the input raw intensity (only preserved into a meta value 'raw_intensity'); and (2) if no feature yields a valid model it calls `features.clear(); return;` (line 323), silently emptying the entire map. A caller passing a FeatureMap to be 'fitted/validated' will be surprised that the map can come back empty and that intensities are mutated.
- **Evidence:** ElutionModelFitter.cpp:322-323 '// no valid feature ... return empty features. ... if (!has_valid_models) { features.clear(); return; }'; setIntensity calls at lines 414 ('feat_it->setIntensity(0.0)'), 428 ('feat_it->setIntensity(area)'), 447 ('(*it)->setIntensity(area)').
- **Fix:** Document on the public fitElutionModels() declaration that it mutates feature intensities in place (raw value moved to meta 'raw_intensity') and that the FeatureMap is cleared when no model is valid. Doc-only change; ABI-safe.

### [FEAT-36] Biosaur2Algorithm::Hill::ion_mobility_median / PeptideFeature::ion_mobility — Field named *_median actually stores an intensity-weighted MEAN, not a median
`medium` · `misleading-name` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/Biosaur2Algorithm.h · _ff-metabo_

```cpp
double ion_mobility_median = -1.0;  /* and */ double ion_mobility = -1.0;
```
- **Expectation:** A field named ion_mobility_median holds the median ion mobility of the hill's points.
- **Actual:** It holds the intensity-weighted MEAN ion mobility; the median is only used as a fallback when total intensity is zero. The member's own doc comment even contradicts the name: 'Intensity-weighted mean ion mobility (median fallback; -1 if not available)'.
- **Evidence:** Biosaur2Algorithm.cpp:1762-1777 computes `weighted_im_sum += hill.ion_mobilities[i] * intensity` then `processed_hill.ion_mobility_median = weighted_im_sum / intensity_sum_im;` (median only in the else branch). Header line 78 doc says 'Intensity-weighted mean ion mobility (median fallback)'.
- **Fix:** Rename to ion_mobility_weighted_mean (and PeptideFeature already exposes ion_mobility, which is fine). Renaming a public struct member is source-breaking; if ABI/source stability is required, keep the field but fix the name in a future major version and meanwhile correct the doc to drop 'median' from the primary description. The drift_time_median field really is a median, so the inconsistency is doubly surprising.
- **Verifier correction:** Severity is medium, not high: the value is still a valid central-tendency estimate of ion mobility (intensity-weighted mean vs. median), so the consequence is a silently-wrong statistic for callers who rely on the name (e.g. expecting median robustness), not a crash or data loss. ABI/source impact of the recommended rename is breaking (public member of an OPENMS_DLLAPI public struct), though the field as it stands has no ABI defect — the breakage is incurred only by the corrective rename. Interim mitigation per the recommendation is sound: fix the doc to drop "median" from the primary description now, and rename to ion_mobility_weighted_mean in a future major version.

### [FEAT-37] ElutionPeakDetection::computeMassTraceSNR — computeMassTraceSNR divides by noise area with no zero guard, returning +inf/NaN for a clean trace
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/ElutionPeakDetection.h · _ff-metabo_

```cpp
double computeMassTraceSNR(const MassTrace&)
```
- **Expectation:** A signal-to-noise getter returns a finite number; if noise is zero (smoothed fit is perfect) it returns a large but finite/sentinel value, consistent with the sibling computeApexSNR which guards against zero noise.
- **Actual:** When the trace is non-empty but the noise RMSE is 0 (e.g. a perfectly-smoothed/short trace), noise_area = 0 and the method computes signal_area / 0 = +inf (or 0/0 = NaN), with no guard. Its sibling computeApexSNR explicitly guards `if (noise_level > 0.0)`.
- **Evidence:** ElutionPeakDetection.cpp:69-84 `double noise_area = computeMassTraceNoise(tr) * tr.getTraceLength(); double signal_area = tr.computePeakArea(); snr = signal_area / noise_area;` with no zero check, vs computeApexSNR (lines 86-99) which checks `if (noise_level > 0.0)`.
- **Fix:** Add the same `if (noise_area > 0.0)` guard as computeApexSNR so the two SNR methods behave consistently (additive, ABI-safe). Document the value returned when noise is zero.
- **Verifier correction:** computeMassTraceSNR can return +inf or NaN when noise_area == 0 (RMSE of raw-vs-smoothed intensities is 0 for a perfectly-smoothed or flat trace), with no zero guard, inconsistent with the sibling computeApexSNR which guards if (noise_level > 0.0). Real but narrower in impact than stated: the internal EPD pipeline uses the guarded computeApexSNR for SNR thresholds (lines 444, 525); computeMassTraceSNR's non-finite return is reachable primarily by external/pyOpenMS callers on degenerate input. Recommended fix (guard noise_area > 0.0, return 0/sentinel otherwise) is internal to the .cpp body — ABI impact none, not source-compatible-only.

### [FEAT-41] FeatureHypothesis::getFWHM / getChromatograms vs getCentroidMZ / getCentroidRT — Empty-hypothesis handling is inconsistent across getters: throw vs return 0.0 vs undefined-behavior deref
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFindingMetabo.h · _ff-metabo_

```cpp
double getFWHM() const;  std::vector<OpenMS::MSChromatogram> getChromatograms(UInt64 feature_id) const;
```
- **Expectation:** All accessors of an empty FeatureHypothesis behave the same way (all throw, or all return a documented sentinel).
- **Actual:** getCentroidMZ() and getCentroidRT() throw InvalidValue on an empty hypothesis; getFWHM() silently returns 0.0; getChromatograms() and getMaxIntensity() never check emptiness and dereference iso_pattern_[0] / iterate, so getChromatograms() on an empty hypothesis is out-of-bounds (UB).
- **Evidence:** FeatureFindingMetabo.cpp:225-233 getCentroidMZ throws on empty; :235-243 getCentroidRT throws on empty; :245-252 getFWHM returns 0.0 on empty; :102-104 getChromatograms `double mz = iso_pattern_[0]->getCentroidMZ();` with no empty guard.
- **Fix:** Pick one contract and apply it uniformly. At minimum add an empty guard to getChromatograms() to avoid UB (additive, ABI-safe), and document whether empty-hypothesis getters throw or return a sentinel.
- **Verifier correction:** getChromatograms() is the only accessor with true UB on an empty hypothesis (unguarded iso_pattern_[0] deref at line 104). getMaxIntensity() is NOT UB — its loop runs zero times on an empty vector and returns 0.0, making it a third silent-sentinel behavior alongside getFWHM(). Accurate statement: empty-hypothesis handling spans THREE behaviors — throw (getCentroidMZ/getCentroidRT/getMonoisotopicFeatureIntensity), silent 0.0 (getFWHM/getMaxIntensity), and unguarded deref/UB (getChromatograms). Severity is medium, not high: in the actual FeatureFindingMetabo pipeline hypotheses always contain >=1 mass trace, so the UB is not triggered internally and produces no silently-wrong OpenMS output for normal runs; it is only reachable by an external caller constructing an empty FeatureHypothesis via the public default ctor. Recommendation stands: add an empty guard to getChromatograms() (additive, ABI-safe) and document/uniformize the empty contract.

### [FEAT-43] MultiplexDeltaMassesGenerator::getLabelShort / getLabelLong — Non-const 'getLabel*' lookups MUTATE the internal map and silently return empty string for unknown labels
`medium` · `misleading-name` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h · _ff-multiplex-a_

```cpp
std::string getLabelShort(const std::string& label); std::string getLabelLong(const std::string& label);
```
- **Expectation:** A 'get' translation accessor on a name->name map is read-only and either returns the translated label or fails loudly when the label is unknown.
- **Actual:** Both methods are non-const and implemented as 'return label_long_short_[label];' / 'return label_short_long_[label];'. std::map::operator[] DEFAULT-CONSTRUCTS and INSERTS an empty entry when the key is absent, so calling getLabelShort("bogus") mutates the object (grows the map with key->"") and returns an empty string instead of signalling 'not found'. This is why the methods cannot be const.
- **Evidence:** Header: 'std::string getLabelShort(const std::string& label);' (line 131), 'std::string getLabelLong(const std::string& label);' (line 138). Impl (MultiplexDeltaMassesGenerator.cpp:545-553): 'std::string ...getLabelShort(const std::string& label){ return label_long_short_[label]; }' and 'getLabelLong{ return label_short_long_[label]; }'.
- **Fix:** Make both methods const and use map::at() (throwing Exception::ElementNotFound on miss) or map::find() returning a clearly-documented empty sentinel without insertion. Const-correctness fix changes the method signature's const-qualifier (mangled name), so do it additively if ABI matters: add 'getLabelShort(...) const' overloads using find(), keep the old ones [[deprecated]].
- **Verifier correction:** The facts are exactly as claimed. Adjustment is only to severity (high -> medium): unknown labels cannot silently corrupt normal output because they are rejected upstream by an explicit validation loop (MultiplexDeltaMassesGenerator.cpp:167-182, throws Exception::InvalidParameter), and the only internal callers iterate over keys that were inserted into the maps. The genuine defects remain: operator[] read silently inserts key->"" (mutation + empty sentinel), which is why the methods are non-const and not callable on a const object. Recommended fix (make const + use at()/find()) is correct; note that changing the const-qualifier changes the mangled symbol name, so it is ABI-breaking unless done additively via new const overloads with the old ones kept [[deprecated]].

### [FEAT-45] MultiplexDeltaMassesGenerator::getSamplesLabelsList — Non-const getSamplesLabelsList() returns a copy by value while the const overload returns a reference
`medium` · `asymmetric-api` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h · _ff-multiplex-a_

```cpp
std::vector<std::vector<std::string> > getSamplesLabelsList(); / const std::vector<std::vector<std::string> >& getSamplesLabelsList() const;
```
- **Expectation:** Same-named const/non-const getter pair: the non-const one returns a modifiable reference to internal state.
- **Actual:** The non-const overload returns 'std::vector<std::vector<std::string> >' by value (copy of samples_labels_); only the const overload returns a reference. Edits to the non-const return value are lost. Identical footgun to getDeltaMassesList.
- **Evidence:** Header lines 115 & 124: 'std::vector<std::vector<std::string> > getSamplesLabelsList();' vs 'const std::vector<std::vector<std::string> >& getSamplesLabelsList() const;'. Impl (MultiplexDeltaMassesGenerator.cpp:535-543): non-const path 'return samples_labels_;' returns by value.
- **Fix:** Return a non-const reference from the non-const overload, or remove it in favour of the const ref overload. Either is an ABI change to the symbol's return type; introduce a new ref-returning accessor and deprecate the value one if ABI stability is required.

### [FEAT-51] MultiplexFilteringProfile::MultiplexFilteringProfile — Constructor takes exp_profile by non-const reference and documents it @param[in,out], but only reads it
`medium` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/MultiplexFilteringProfile.h · _ff-multiplex-b_

```cpp
MultiplexFilteringProfile(MSExperiment& exp_profile, const MSExperiment& exp_centroided, ...)
```
- **Expectation:** A non-const MSExperiment& parameter documented as [in,out] tells the caller the profile experiment will be modified, so callers must pass a mutable, non-shared copy.
- **Actual:** The body only iterates `for (const auto& spectrum : exp_profile)` and feeds each spectrum to `SplineInterpolatedPeaks(const MSSpectrum&)` (const ctor). exp_profile is never written. The non-const ref and the [in,out] doc are both wrong, forcing callers to give up const-correctness for no reason.
- **Evidence:** Header line 48 doc: `@param[in,out] exp_profile`. Header line 66 signature `MSExperiment& exp_profile`. Impl MultiplexFilteringProfile.cpp:86-89 only reads: `for (const auto& spectrum : exp_profile) { exp_spline_profile_.emplace_back(spectrum); }`. SplineInterpolatedPeaks.h:48 takes `const MSSpectrum&`.
- **Fix:** Change the parameter to `const MSExperiment& exp_profile` and the doc to @param[in]. This is source-compatible for callers passing lvalues (const-ref binds more, not less), but changes the mangled constructor signature (ABI). If ABI must be frozen, at minimum fix the @param[in,out] doc to @param[in].

### [FEAT-52] MultiplexFilteredPeak::MultiplexFilteredPeak — Constructor argument order (mz, rt, mz_idx, rt_idx) is reversed relative to the class's own addSatellite/checkSatellite (rt_idx, mz_idx)
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/MultiplexFilteredPeak.h · _ff-multiplex-b_

```cpp
MultiplexFilteredPeak(double mz, float rt, size_t mz_idx, size_t rt_idx)
```
- **Expectation:** Within one class, the two interchangeable size_t index parameters (rt index of spectrum, mz index within spectrum) should appear in a consistent order so call sites are not silently swappable.
- **Actual:** The constructor orders the indices m/z-then-RT (mz_idx, rt_idx), while addSatellite(size_t rt_idx, size_t mz_idx, ...) and checkSatellite(size_t rt_idx, size_t mz_idx) order them RT-then-m/z. Both pairs are bare size_t, so swapping them compiles cleanly and silently mislabels which axis the index refers to — a classic silent-bug shape.
- **Evidence:** Header line 50: `MultiplexFilteredPeak(double mz, float rt, size_t mz_idx, size_t rt_idx);` vs line 75 `void addSatellite(size_t rt_idx, size_t mz_idx, size_t pattern_idx);` and line 89 `bool checkSatellite(size_t rt_idx, size_t mz_idx) const;`. Sibling class MultiplexSatelliteCentroided.h:35 also uses (rt_idx, mz_idx).
- **Fix:** Cannot reorder the ctor without breaking ABI; instead document the order explicitly in the header and consider a named factory or strong index types (e.g. RTIndex/MZIndex strong typedefs, which this module already includes boost strong_typedef for) so the axes cannot be transposed. Flagging the consistent order as the ideal fix.
- **Verifier correction:** The intra-class index-order inconsistency is real and verified: ctor uses (mz_idx, rt_idx) while addSatellite/checkSatellite (and sibling MultiplexSatelliteCentroided) use (rt_idx, mz_idx); even the private member declaration order differs between the two classes (mz_idx_,rt_idx_ vs rt_idx_,mz_idx_). But the original severity is downgraded from an implied high to medium: both real call sites pass arguments correctly so there is no current defect, and a transposition would most likely trigger out-of-bounds/visibly-wrong blacklist indexing rather than silently plausible results. abi_impact corrected to 'none' for the flagged issue (only the suggested reorder fix would be ABI-breaking, which is already noted in the recommendation).

### [PROC-26] SignalToNoiseEstimator::estimate_ — estimate_() divides by element count without an emptiness guard: on an empty range it divides by zero, producing NaN/inf rather than a defined result
`medium` · `surprising-throw` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h · _proc-calibration-noise_

```cpp
GaussianEstimate estimate_(const PeakIterator& scan_first_, const PeakIterator& scan_last_) const
```
- **Expectation:** A statistics helper computing mean/variance over [first,last) should handle an empty range deterministically (return 0, or signal). The S/N base class is documented to operate on raw data intervals.
- **Actual:** estimate_ counts `size` while summing, then does `m = m / size;` and `v = v / (double)size;`. If scan_first_ == scan_last_ (empty spectrum, which init()/computeSTN_ can legitimately receive), size==0 and these are 0/0 -> NaN, silently flowing into max_intensity_ estimation (AUTOMAXBYSTDEV path) and downstream S/N values.
- **Evidence:** `int size = 0; ... while(run!=scan_last_){ m += ...; ++size; ++run; } m = m / size; ... v = v / ((double)size);` with no `if (size==0)` guard.
- **Fix:** Guard the empty range (return {0,0} or throw a documented exception) before dividing. Source-compatible; ABI unchanged. computeSTN_ callers already special-case empty windows, so a defined empty-range result is the least surprising behavior.
- **Verifier correction:** estimate_() has no guard for an empty input range and divides by size==0, producing NaN mean/variance that silently propagate into max_intensity_ (the default AUTOMAXBYSTDEV path) rather than no exception being thrown. The defect is real but two supporting statements in the original claim are inaccurate: the existing 'noise_for_empty_window' handling covers sparse sliding windows, not an empty whole-spectrum input; and the main centroiding caller (PeakPickerHiRes) already guards input.size() < 5 before init(), so the NaN is triggered chiefly via the unguarded callers (PeakPickerIterative, PeakPickerChromatogram, EICExtractor, TargetedSpectraExtractor) and direct/pyOpenMS use on empty or near-empty spectra. Recommended fix: add `if (scan_first_ == scan_last_) return {0.0, 0.0};` (or throw a documented exception) at the top of estimate_(); this is an inline-body change, source-compatible, ABI unchanged.

### [PROC-27] PrecursorCorrection::correctToNearestFeature — Default rt_tolerance_s and mz_tolerance are both 0.0, so the default-argument call silently corrects (almost) nothing
`medium` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/CALIBRATION/PrecursorCorrection.h · _proc-calibration-noise_

```cpp
static std::set<Size> correctToNearestFeature(const FeatureMap& features, MSExperiment& exp, double rt_tolerance_s = 0.0, double mz_tolerance = 0.0, bool ppm = true, bool believe_charge = false, bool keep_original = false, bool all_matching_features = false, int max_trace = 2, int debug_level = 0)
```
- **Expectation:** Calling correctToNearestFeature(features, exp) with defaults would be expected to do a reasonable nearest-feature correction. The name implies an active correction step.
- **Actual:** With rt_tolerance_s = 0.0 and mz_tolerance = 0.0 (and ppm = true), the RT/m/z acceptance windows collapse to zero width, so essentially no precursor can fall inside a feature's (zero-extended) bounding box / mass trace tolerance. The function returns an empty set and corrects nothing, while appearing to have run. The zero defaults give a silent no-op rather than a useful default.
- **Evidence:** Header defaults: `double rt_tolerance_s = 0.0, double mz_tolerance = 0.0, bool ppm = true`. overlaps_() extends the bbox by rt_tolerance (0) and compatible_() uses mz_tolerance (0).
- **Fix:** Either drop the defaults (force callers to choose tolerances) or pick non-zero sensible defaults, and document that 0 means 'no match'. Removing defaults is source-breaking for default-arg callers; safest is to document the no-op behavior. ABI of the symbol itself is unaffected by default-arg changes.
- **Verifier correction:** The default-argument call correctToNearestFeature(features, exp) is a silent no-op, but specifically because the default mz_tolerance=0.0 makes compatible_()'s strict comparison `mass_error < mz_tolerance` (mass_error = fabs(...) >= 0) always false, so no feature is ever accepted and an empty set is returned. The default rt_tolerance_s=0.0 does NOT collapse the RT acceptance window: overlaps_() uses the feature's actual convex-hull bounding box (which already has real RT/mz extent) and only adds rt_tolerance as extra margin (m/z is even hardcoded-extended by 0.01). Recommendation stands: pick non-zero sensible defaults or document that mz_tolerance=0 yields no matches; note real callers always pass explicit tolerances.

### [PROC-28] SignalToNoiseEstimatorMedianRapid::NoiseEstimator::get_noise_value — get_noise_value() clamps noise to a minimum of 1.0 and its window lookups are guarded only by assert(), so out-of-range m/z is UB in release and a floor of 1.0 is silently imposed
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h · _proc-calibration-noise_

```cpp
double get_noise_value(double mz)
```
- **Expectation:** A 'get noise value at this m/z' accessor would be expected to return the estimated noise for that m/z, and to handle an m/z outside the spectrum's range safely.
- **Actual:** get_noise_value returns `std::max(1.0, (even+odd)/2.0)` -- it silently floors the noise at 1.0 (documented only in a @note, surprising if intensities are <1). More seriously, get_noise_even/odd compute a window index and only `assert(window_nr >= 0/< size)`; asserts are compiled out in release (NDEBUG), so passing an mz < mz_start or beyond the last window indexes the vector out of bounds (UB) instead of erroring.
- **Evidence:** `return std::max(1.0, (get_noise_even(mz)+get_noise_odd(mz))/2.0);` and in get_noise_even: `int window_nr = (int)((mz - mz_start)/window_length); assert(window_nr >= 0); assert(window_nr < (int)result_windows_even.size()); double noise = result_windows_even[window_nr];`
- **Fix:** Replace asserts with real bounds handling (clamp window_nr to [0, size-1] or throw) so release builds are safe; document the 1.0 floor in the function brief, not just a note. Source-compatible; these are header-inline so changing them is an ABI consideration only for inlined call sites (recompile suffices).
- **Verifier correction:** The genuine, load-bearing defect is that the public header-inline accessors get_noise_even/get_noise_odd guard their computed window index with assert() only; since OpenMS ships Release with -DNDEBUG, an out-of-range m/z on the public get_noise_value() accessor is undefined behavior (out-of-bounds vector access) in release rather than an error. The 1.0 floor is a documented (@note), deliberate divide-by-zero-avoidance convention and is NOT itself a surprise. Severity is medium, not high: the out-of-range path is not reachable by OpenMS's own callers (FIAMSDataProcessor only queries m/z values that were in the input array) and is not exposed via pyOpenMS, so no in-tree reasonable use produces wrong results or crashes; the risk is to third-party C++ consumers who query outside the spectrum range. Recommended fix stands: replace the asserts with real bounds handling (clamp window_nr to [0,size-1] or throw) and document any in-range precondition on the public brief. ABI: bodies are header-inline, so the change is source-compatible and only requires recompilation of inlined call sites (no layout/signature change).

### [PROC-3] PeakPickerIterative::pick — pick() leaves 'output' untouched (not cleared) when input has <3 points, but fully overwrites it otherwise — inconsistent output contract
`medium` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/CENTROIDING/PeakPickerIterative.h · _proc-centroiding_

```cpp
void pick(const MSSpectrum& input, MSSpectrum& output)
```
- **Expectation:** An output-parameter pick(input, output) should leave 'output' in a well-defined state on every path: either always cleared-then-filled, or documented to require an empty output. A caller reusing an 'output' object across spectra expects stale peaks to be cleared.
- **Actual:** On the normal path pick() copies meta, clears float arrays, and rebuilds 'output' from scratch. But when 'input.size() < 3' it returns immediately, leaving whatever peaks/meta 'output' already contained — so a reused output silently keeps the previous spectrum's peaks for short inputs.
- **Evidence:** PeakPickerIterative.h:276-282: 'if (input.size() < 3) return; // copy the spectrum meta data; copySpectrumMeta(input, output); output.setType(...); output.getFloatDataArrays().clear();' — the early return happens before any clearing of output.
- **Fix:** Move 'output.clear(true)' (or copySpectrumMeta) above the size check so 'output' is always reset, matching the normal path. Header-only inline change; ABI: source-compatible (behavior change only).

### [PROC-5] PeakPickerIM::pickIMCluster / pickIMElutionProfiles / pickIMTraces — All three pickIM* methods silently no-op (return without modifying the spectrum) when the input has no IM data, instead of signaling that nothing was done
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h · _proc-centroiding_

```cpp
void pickIMTraces(MSSpectrum&); void pickIMCluster(MSSpectrum&) const; void pickIMElutionProfiles(MSSpectrum&) const;
```
- **Expectation:** Given the documented in-place contract ('reduces it to a single spectrum', 'extracts IM elution profiles'), a caller passing a spectrum that turns out to lack IM data would expect either an exception (as the doc for pickIMCluster even advertises: '@throws MissingInformation if input spectrum lacks ion mobility data') or at least a return value indicating no work was done.
- **Actual:** validateIMFormatForPicking returns false for IMFormat::NONE and the method returns early, leaving the spectrum unchanged with no error and no return value. For pickIMCluster this directly contradicts its '@throws Exception::MissingInformation if input spectrum lacks ion mobility data' documentation — no such exception is thrown for the NONE case.
- **Evidence:** PeakPickerIM.cpp validateIMFormatForPicking: 'case IMFormat::NONE: return false;' (lines 737-738); callers: pickIMTraces 'if (!validateIMFormatForPicking(spectrum)) { return; }' (765-769), pickIMCluster (994-997), pickIMElutionProfiles (1236-1239). The header for pickIMCluster (PeakPickerIM.h:79) claims '@throws Exception::MissingInformation if input spectrum lacks ion mobility data'.
- **Fix:** Either make the no-IM-data case throw (matching the documented MissingInformation contract) or remove the @throws claim and document the silent-no-op behavior. All methods return void, so callers cannot detect the no-op; consider returning bool. ABI: doc fix = none; throwing or return-type change = source/behavior breaking.
- **Verifier correction:** For pickIMCluster the surprise is a documented-contract violation: its header advertises '@throws Exception::MissingInformation if input spectrum lacks ion mobility data', but when determineIMFormat yields IMFormat::NONE, validateIMFormatForPicking returns false and pickIMCluster returns void silently — no MissingInformation (or any) exception is thrown (the only nearby guard, containsIMData(), merely LOG_WARNs and returns). For pickIMTraces and pickIMElutionProfiles there is no @throws claim; their issue is the narrower one of an undocumented silent no-op (void return leaves the spectrum unchanged with no signal). Fix: either throw on the NONE case to honor the documented contract, or remove the @throws line and document the silent-skip behavior; a bool/return-value signal would be source-breaking.

### [PROC-6] PeakPickerIM::pickIMCluster — pickIMCluster silently discards all StringDataArrays and IntegerDataArrays and all-but-one FloatDataArray, beyond just averaging peaks
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h · _proc-centroiding_

```cpp
void pickIMCluster(MSSpectrum& spec) const
```
- **Expectation:** A method documented as 'reduces it to a single spectrum where peaks ... are averaged together' with intensity-weighted m/z and IM is expected to rewrite peaks and the IM array, but a caller would not expect unrelated per-peak metadata (e.g. annotated charges, labels, extra float arrays) to be wiped.
- **Actual:** The method clears all StringDataArrays and IntegerDataArrays and then strips every FloatDataArray except the IM-centroid array. Any caller-attached metadata is lost without mention in the public doc.
- **Evidence:** PeakPickerIM.cpp pickIMCluster: 'spectrum.getStringDataArrays().clear(); spectrum.getIntegerDataArrays().clear();' (1206-1207) and 'removeAllFloatDataArraysExcept(spectrum, Constants::UserParam::ION_MOBILITY_CENTROID);' (1228). pickIMElutionProfiles does the same (1300-1301, 1324). The header doc (lines 63-92) does not mention discarding data arrays.
- **Fix:** Document in the public header that all non-IM data arrays are dropped (the .cpp explains why, the header does not). Doc-only; ABI: none.
- **Verifier correction:** Severity is medium rather than high: the dropped per-peak arrays (charges/labels/extra floats) are silently lost and unrecoverable, but the operation merges peaks and thereby genuinely invalidates any array indexed by the original peak positions, so a careful caller could anticipate per-peak arrays being invalidated. The surprising part is that it is silent and undocumented in the public header (the .cpp explains the rationale; the header does not), and that even unrelated extra FloatDataArrays are stripped down to only the IM-centroid array. Recommendation stands: document in the header that pickIMCluster and pickIMElutionProfiles drop all StringDataArrays, IntegerDataArrays, and every FloatDataArray except the IM-centroid array. Doc-only fix; no signature/ABI change.

### [PROC-17] RankScaler::filterSpectrum — 'Scaler' overwrites intensities with integer ranks rather than scaling them
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/SCALING/RankScaler.h · _proc-filtering-scaling_

```cpp
template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
```
- **Expectation:** Given the class name 'RankScaler' and the SCALING module, a developer expects intensities to be multiplied/transformed while retaining proportional meaning (like SqrtScaler, Normalizer).
- **Actual:** Each peak's intensity is replaced by its (dense) rank count: the most intense peak gets intensity == number of distinct intensity levels, the least intense gets 1. Original intensity magnitudes are discarded; this is a categorical re-labeling, not a scaling. Equal intensities share a rank (dense ranking).
- **Evidence:** RankScaler.h: 'typename SpectrumType::size_type count = spectrum.size(); ++count; ... if (it->getIntensity() != last_int) { --count; } ... it->setIntensity(count);'
- **Fix:** Document precisely in the (currently sparse) class doc that intensities are replaced by descending dense ranks (1 = weakest). Renaming to RankTransformer would be clearer but is ABI-breaking; prefer the doc fix plus noting the highest-rank == n_distinct convention.
- **Verifier correction:** Each peak's intensity is overwritten with an integer counter, NOT a 1-based dense rank. After ascending sort, the algorithm starts count=n+1 and walks from the most intense peak to the least, decrementing count by 1 only when the intensity value changes. Result: the MOST intense peak always gets n (number of peaks), and the value decreases by 1 at each distinct intensity level going down; tied peaks share a value. The least intense peak gets n-(d-1) (where d = number of distinct levels), which equals 1 only when all intensities are distinct. Example: [1,2,2,2,9] (sorted) -> [3,4,4,4,5]. The claimed 'most intense == number of distinct levels' and 'least intense == 1' hold only in the no-ties case. The genuine, well-supported point stands: original magnitudes are discarded and replaced by integer rank labels despite the 'Scaler'/'Scales each peak' naming and SCALING-module placement. Fix via precise documentation (ABI-safe); a rename to RankTransformer would be ABI-breaking and is not required.

### [PROC-18] SqrtScaler::filterSpectrum — Negative intensities are silently clamped to 0 and a warning is written to std::cerr
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/SCALING/SqrtScaler.h · _proc-filtering-scaling_

```cpp
template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
```
- **Expectation:** A sqrt transform on out-of-domain (negative) input should signal an error through the library's logging/exception channel; a caller does not expect raw std::cerr output or silent data alteration.
- **Actual:** Negative intensities are set to 0 (so sqrt(0)=0) and a one-line warning is emitted directly to std::cerr, bypassing OpenMS's LogStream/OPENMS_LOG_WARN facility. The clamping is lossy and not reflected in any return value or status.
- **Evidence:** SqrtScaler.h: 'if (intens < 0) { intens = 0; warning = true; } it->setIntensity(std::sqrt(intens)); ... if (warning) { std::cerr << "Warning negative intensities were set to zero" << std::endl; }'
- **Fix:** Route the warning through OPENMS_LOG_WARN (consistent with the rest of OpenMS) and document the clamp-to-zero behavior in the header. ABI-safe inline-body change.
- **Verifier correction:** Negative intensities are clamped to 0 (yielding sqrt(0)=0) and a one-line warning is emitted to std::cerr instead of OPENMS_LOG_WARN. The behavior is loud (a warning fires, once per spectrum) rather than truly silent, but it bypasses OpenMS's configurable LogStream facility and performs an undocumented lossy data alteration with no count, context, or return status. Fix: route through OPENMS_LOG_WARN and document the clamp-to-zero behavior in the header; this is an inline template-body change (source-compatible, no exported-symbol/layout change).

### [PROC-19] WindowMower::filterPeakSpectrumForTopNInJumpingWindow — Last-window peak count is scaled by window-fill fraction and can round to zero, dropping all peaks at the spectrum tail
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/FILTERING/WindowMower.h · _proc-filtering-scaling_

```cpp
template <typename SpectrumType> void filterPeakSpectrumForTopNInJumpingWindow(SpectrumType& spectrum)
```
- **Expectation:** A 'top N in window' filter keeps up to peakcount peaks per window, including the final window.
- **Actual:** For the trailing window the retained count is round(fraction_of_window_filled * peakcount). If the last window covers a small m/z fraction, last_window_peakcount can round to 0, so ALL peaks in the final window are discarded even if they are intense. This per-window asymmetry (full peakcount for all windows except a fractional count for the last) is not evident from the name/signature.
- **Evidence:** WindowMower.h: 'double last_window_size_fraction = last_window_size / windowsize_; Size last_window_peakcount = static_cast<Size>(std::round(last_window_size_fraction * peakcount_));' then 'if (peaks_in_window.size() > last_window_peakcount) { std::partial_sort(... begin()+last_window_peakcount ...) }'.
- **Fix:** Document the fractional last-window behavior in the header, and consider std::max<Size>(1, ...) so a non-empty final window keeps at least one peak. ABI-safe inline-body change.

### [PROC-20] WindowMower::filterPeakSpectrumForTopNInJumpingWindow — Retained peaks are matched back by full-peak equality (std::find), conflating duplicate peaks
`medium` · `return-value` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/FILTERING/WindowMower.h · _proc-filtering-scaling_

```cpp
template <typename SpectrumType> void filterPeakSpectrumForTopNInJumpingWindow(SpectrumType& spectrum)
```
- **Expectation:** The window selection picks specific peaks; the final selection should keep exactly those chosen peaks.
- **Actual:** After choosing 'out' peaks, the function maps them back onto the original spectrum via 'std::find(out.begin(), out.end(), *it)', i.e. by value equality of the whole peak (m/z + intensity). Two distinct original peaks with identical m/z AND intensity are indistinguishable: both are kept (or both dropped) regardless of which the window logic intended. This is also O(n^2).
- **Evidence:** WindowMower.h: 'if (std::find(out.begin(), out.end(), *it) != out.end()) { Size index(it - spectrum.begin()); indices.push_back(index); } ... spectrum.select(indices);'
- **Fix:** Track selected indices directly through the window loop instead of value-matching, or document the duplicate-peak ambiguity. ABI-safe inline-body change.
- **Verifier correction:** After choosing the top-N peaks per window into `out`, the function recovers original indices via `std::find(out.begin(), out.end(), *it)`, i.e. by full-value equality (Peak1D::operator== is defaulted: m/z + intensity). When two or more original peaks share identical m/z AND intensity, they are indistinguishable: if the window logic selected only some of them, the value match keeps ALL of them, producing a wrong retained-peak count and mis-associating the parallel data-arrays (which select() reorders by index). It is also O(n^2). Note the peak m/z/intensity values kept are not themselves 'wrong' (the duplicates are value-identical); the observable defects are the count and data-array metadata mapping. Fix by tracking selected indices through the window loop instead of value-matching. Inline template body change in a header: signature unchanged, source-compatible (callers recompile, no symbol/ABI break).

### [PROC-21] Normalizer::filterSpectrum — 'to_one'/'to_TIC' divide by a divisor that can be zero or negative without guarding
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/SCALING/Normalizer.h · _proc-filtering-scaling_

```cpp
template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum) const
```
- **Expectation:** A normalizer either produces values in [0,1] (to_one) or summing to 1 (to_TIC), or signals when normalization is impossible.
- **Actual:** For 'to_one' the divisor is the max intensity but is seeded from the first peak's intensity 'as a safety measure'; if all intensities are <= 0 the divisor can be <= 0, producing negative or inf/NaN results. For 'to_TIC' divisor is the sum of intensities and is not checked for 0; a spectrum whose intensities sum to 0 yields division by zero (inf/NaN) silently. No exception or status is raised for these degenerate inputs.
- **Evidence:** Normalizer.h: to_one 'divisor = spectrum.begin()->getIntensity(); ... if (divisor < it->getIntensity()) divisor = it->getIntensity();'; to_TIC 'divisor += it->getIntensity();'; then unconditional 'it->setIntensity(it->getIntensity() / divisor);'.
- **Fix:** Guard divisor <= 0 (skip or throw with a clear message) before the division loop, and document behavior for all-zero/negative spectra. ABI-safe inline-body change.
- **Verifier correction:** The 'to_one' seeding from the first peak's intensity is a partial mitigation explicitly noted in the line-71 comment, but it is incomplete: (a) it does not prevent division by zero when the maximum intensity is exactly 0 (all-zero or zero-max spectra), and (b) on all-negative spectra it produces a negative divisor, yielding sign-flipped 'normalized' values outside [0,1] rather than failing loudly. 'to_TIC' has no guard at all and divides by zero when intensities sum to 0. No exception or status signals these degenerate cases; only the unknown-method path throws (and is the only documented throw). Severity is medium (not high): the failures require degenerate but plausible spectra, produce silently-wrong data (inf/NaN/negative) that propagates downstream, but do not crash and are recoverable. Recommended fix: guard divisor <= 0 (skip the spectrum or throw a clear message) before the division loop and document all-zero/negative behavior; this is an inline header-body change, source-compatible (recompile of dependents only, no signature/layout change).

### [PROC-22] SpectraMerger::average — average() reads Gaussian-named parameters even when average_type is 'tophat'
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h · _proc-filtering-scaling_

```cpp
template <typename MapType> void average(MapType& exp, const std::string& average_type, int ms_level = -1)
```
- **Expectation:** Selecting average_type=='tophat' should configure the run entirely from the average_tophat:* parameters.
- **Actual:** The defaults for ms_level and spectrum_type are first read from average_gaussian:* and only overwritten by average_tophat:* when average_type=='tophat'. More surprisingly, the Gaussian-only parameters (average_gaussian:rt_FWHM, cutoff, precursor_mass_tol, precursor_max_charge) are always fetched and the precursor-mass matching gate (precursor_mass_ppm >= 0) is applied for BOTH gaussian and tophat runs, because precursor_mass_ppm is taken unconditionally from average_gaussian:precursor_mass_tol. So a tophat average is silently influenced by a gaussian-namespaced parameter.
- **Evidence:** SpectraMerger.h: 'double precursor_mass_ppm = param_.getValue("average_gaussian:precursor_mass_tol");' (no tophat branch) followed by the shared loop 'if (precursor_mass_ppm >= 0 && ms_level >= 2 ...) { ... add = areMassesMatched(...); }' used for both average types.
- **Fix:** Read precursor_mass_tol/max_charge from the active method's namespace (or document that average_gaussian:precursor_mass_tol governs tophat too). Behavioral/source change; keep default value so existing configs are unaffected.
- **Verifier correction:** The cross-namespace coupling is real: tophat averaging reads precursor_mass_tol/precursor_max_charge only from the average_gaussian namespace (lines 371-372) and applies the precursor-mass-matching gate (lines 410, 463) to tophat runs, with no average_tophat:precursor_* parameter existing to control or disable it. However, this influence is NOT unconditional as the claim implies: with the default value 0.0, areMassesMatched (line 298) returns true immediately (tol_ppm <= 0), making the gate a no-op, so default tophat behavior is unchanged. The surprise only occurs when a user explicitly sets average_gaussian:precursor_mass_tol > 0 and runs tophat. The ms_level/spectrum_type "default-from-gaussian-then-overwrite" asymmetry (lines 351-365) is also real but harmless because those are correctly overwritten for tophat. Recommendation stands: either add average_tophat:precursor_mass_tol/precursor_max_charge parameters and read from the active method's namespace, or document that average_gaussian:precursor_mass_tol governs tophat too.

### [PROC-34] FeatureOverlapFilter::filter(FeatureMap&, FeatureComparator, FeatureOverlapCallback, bool) — Legacy filter() defaults to TRACE_LEVEL overlap, contradicting the enum's documented CONVEX_HULL default
`medium` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/FEATURE/FeatureOverlapFilter.h · _proc-misc_

```cpp
static void filter(..., bool check_overlap_at_trace_level = true)
```
- **Expectation:** FeatureOverlapMode documents CONVEX_HULL as '(default)'. A developer expects the default overlap-detection strategy across the API to be convex-hull bounding boxes.
- **Actual:** The bool overload defaults check_overlap_at_trace_level = true, and the cpp maps true -> FeatureOverlapMode::TRACE_LEVEL (line 128). So the convenient default-argument call performs trace-level overlap (which additionally throws MissingInformation if subordinate convex hulls are absent), not the CONVEX_HULL mode the enum calls the default. Two sibling APIs in the same class disagree on what 'default' overlap means.
- **Evidence:** FeatureOverlapFilter.h:19 `CONVEX_HULL, ///< ... (default)`; :71 `bool check_overlap_at_trace_level = true`; FeatureOverlapFilter.cpp:128 `mode = check_overlap_at_trace_level ? TRACE_LEVEL : CONVEX_HULL`.
- **Fix:** Align the documented default: either change the enum comment to reflect TRACE_LEVEL being the de-facto default of the legacy entry point, or flip the bool default to false (behavior-changing). Safest additive fix is to clarify the docs; if behavior should change, do it in a minor release. ABI of the signature itself is unchanged.
- **Verifier correction:** The inconsistency is real and exactly as quoted, but it is a documentation/API-design conflict rather than a silent-wrong-results bug: the bool default of true is intentional backward-compatibility (historically the legacy filter always did trace-level overlap; the only production caller passes true explicitly). The mismatch is between the newly-added enum's "(default)" annotation on CONVEX_HULL and the legacy bool overload's de-facto TRACE_LEVEL default. Severity is medium, not high — the divergent path that lacks subordinate convex hulls throws Exception::MissingInformation (loud/recoverable) rather than silently corrupting data; when hulls are present the two modes simply yield different (trace-level being stricter) results. ABI of the signature is unchanged; the safest fix is to align the enum doc comment (or flip the bool to false in a minor release, which would be behavior-changing but still source/ABI-compatible).

### [PROC-10] GaussFilter::filter — filter() silently leaves data unmodified when the Gaussian is too narrow, despite documenting that it throws
`medium` · `silent-failure` · ABI: `none` · src/openms/source/PROCESSING/SMOOTHING/GaussFilter.cpp · _proc-smoothing_

```cpp
void GaussFilter::filter(MSSpectrum & spectrum)
```
- **Expectation:** The header doc states: '@exception Exception::IllegalArgument is thrown, if the @em gaussian_width parameter is too small.' A caller therefore expects either a correctly smoothed spectrum or a thrown exception when the width is too small relative to the data spacing.
- **Actual:** When the filter produces no signal (found_signal==false, i.e. the gaussian_width is smaller than the data spacing) the method throws nothing. It at most emits a log warning, and only if the non-default param write_log_messages_ is true. Critically, in the !found_signal branch the smoothed output is NOT copied back, so the spectrum is silently returned UNCHANGED. The only throw IllegalArgument in the file (lines 99/148) is for the unrelated 'ppm tolerance on chromatograms/mobilograms' case, not the documented 'gaussian_width too small' case.
- **Evidence:** Header GaussFilter.h:60 '@exception Exception::IllegalArgument is thrown, if the @em gaussian_width parameter is too small.' ; cpp:72 'if (!found_signal && spectrum.size() >= 3) { if (write_log_messages_) { ... OPENMS_LOG_WARN ... } } else { /* copy the new data into the spectrum */ }'. The copy-back only happens in the else branch, so the failure path leaves the input untouched and no exception is thrown.
- **Fix:** Either throw Exception::IllegalArgument in the !found_signal path (matching the documented contract) or update the header to remove the false @exception promise and document that the spectrum is left unchanged on failure. Additive/ABI-safe fix: correct the documentation and optionally make the warning unconditional; behavior-changing fix (throwing) would be source-breaking for callers relying on silent no-op, so flag it but prefer doc correction plus a settable strict mode.
- **Verifier correction:** Severity is medium rather than high: the input spectrum is left unchanged (original profile data is preserved, not zeroed/corrupted, and there is no crash), so the worst case is downstream code silently processing UNsmoothed data when it believed smoothing succeeded. It is recoverable and a too-narrow width is a misconfiguration, but the documented @exception contract is outright false and the default path is a silent no-op (warning is off by default), which is the surprise. The finding as filed is a documentation/behavior mismatch; abi_impact is none for the discrepancy itself and for the recommended doc-correction fix (the alternative "throw" fix would be source-breaking, but that is a proposed remediation, not the finding).

### [PROC-11] LowessSmoothing::smoothData — Output vector is appended-to (push_back) in the normal path but assigned in the degenerate path, so it is not cleared
`medium` · `asymmetric-api` · ABI: `none` · src/openms/source/PROCESSING/SMOOTHING/LowessSmoothing.cpp · _proc-smoothing_

```cpp
void LowessSmoothing::smoothData(const DoubleVector& input_x, const DoubleVector& input_y, DoubleVector& smoothed_output)
```
- **Expectation:** An out-parameter named smoothed_output of a smoothData() call is expected to be (re)filled with exactly the smoothed series, i.e. cleared/overwritten so its final size equals the input size regardless of its prior contents.
- **Actual:** In the main path the method does smoothed_output.push_back(...) for each point WITHOUT first clearing smoothed_output. If the caller reuses a non-empty vector (e.g. calls smoothData twice with the same buffer), the new results are appended after the old ones and the size becomes wrong. In the degenerate path (input_x.size() <= 2) it instead does 'smoothed_output = input_y', which overwrites. The two paths thus disagree on whether prior contents survive.
- **Evidence:** cpp:38 'smoothed_output = input_y; return;' (overwrite path) vs cpp:77 'smoothed_output.push_back(qr.eval(rt));' inside the per-point loop with no preceding smoothed_output.clear().
- **Fix:** Add 'smoothed_output.clear();' (and optionally reserve(input_size)) at the top of smoothData so both paths consistently overwrite. ABI-safe, body-only change.
- **Verifier correction:** The asymmetry and missing clear are real and exactly as quoted (cpp:38 overwrite vs cpp:77 append-without-clear). However, no current production caller is affected: pyOpenMS passes a fresh empty buffer, the only in-tree C++ caller is commented-out, and the unit test manually clears the buffer before each reuse (LowessSmoothing_test.cpp:77/86/95) — which itself evidences the missing internal clear. Hence the impact is "silent wrong-sized output on buffer reuse" for a plausible-but-currently-unused pattern, warranting medium (invites silent misuse, recoverable, not presently triggered) rather than high. Fix is ABI-safe and body-only.

### [PROC-13] SavitzkyGolayFilter::filter — filter() is a silent no-op when the spectrum has fewer points than the frame size
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h · _proc-smoothing_

```cpp
void filter(MSSpectrum & spectrum)
```
- **Expectation:** Calling filter(spectrum) is expected to smooth the spectrum or signal an error if it cannot (e.g. throw if the frame is wider than the data).
- **Actual:** The low-level template begins with 'if (frame_size_ > n) { return; }', returning before writing anything. Because the public filter() pre-copies the spectrum into 'output' and then swaps, a too-short spectrum is silently returned unchanged with no warning or exception. The caller gets no indication that no smoothing was applied.
- **Evidence:** SavitzkyGolayFilter.h:94 'if (frame_size_ > n) { return; }'; the public wrapper (h:163-171) does 'MSSpectrum output = spectrum; filter(spectrum.begin(), spectrum.end(), output.begin()); std::swap(spectrum, output);' so an early return leaves output==spectrum (unchanged).
- **Fix:** Document that spectra shorter than frame_size are returned unchanged, and optionally emit a ProgressLogger/LOG warning. ABI-safe (doc + optional warning). A stricter throw would be source-breaking, so flag but do not force it.
- **Verifier correction:** filter() (and the equivalent MSChromatogram and Mobilogram overloads, plus filterExperiment) silently returns the input unchanged — no warning, log, or exception — whenever the container has fewer points than frame_size_, because the low-level template returns at line 94 before writing any output and the wrapper's pre-copied `output` (== input) is then swapped back. Recommended fix is doc + optional warning (ABI-safe, none).

### [FEAT-20] FeatureFinderAlgorithmMetaboIdent::setMSData — setMSData() silently drops all non-MS1 spectra from the supplied map
`low` · `incomplete-documentation` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h · _ff-algorithms_

```cpp
void setMSData(const PeakMap& m); void setMSData(PeakMap&& m);
```
- **Expectation:** A plain 'setMSData(map)' setter is expected to store the map the caller passed; a subsequent getMSData() should reflect the same spectra (a setter/getter pair).
- **Actual:** Both overloads erase every spectrum whose MS level != 1 before storing, so MS2+ spectra are silently discarded. A caller who sets a full PeakMap and reads it back (where applicable) sees a different, filtered map. The const-ref overload also copies before filtering. This filtering side effect is not mentioned in the brief doc '/// @brief Set spectra'.
- **Evidence:** Source FeatureFinderAlgorithmMetaboIdent.cpp:962-986: 'ms_data_ = m; ... specs.erase(std::remove_if(... [](const MSSpectrum & s) { return s.getMSLevel() != 1; }), specs.end());' in both overloads.
- **Fix:** Document in the header that setMSData() retains only MS1 spectra. Ideal additive fix: rename intent in docs (e.g. 'sets MS data; non-MS1 spectra are dropped'). Doc-only change is ABI-safe.
- **Verifier correction:** setMSData() retains only MS1 spectra in the algorithm's internal copy (intentional, since this is an MS1-only feature finder and matches the sibling FeatureFinderIdentificationAlgorithm). The behavior is partially documented: the rvalue overload's header comment already warns getMSData() returns a "processed/modified version." The genuine gap is only the const-ref overload's terse "/// @brief Set spectra", which should also note non-MS1 spectra are dropped. No caller-side data loss (const-ref filters a copy) and no wrong results, so this is a low-severity documentation nit, not a hidden side-effect.

### [FEAT-22] FeatureFinderAlgorithmPickedHelperStructs::MassTraces::isValid — isValid(seed_mz, trace_tolerance) is a read-only check but declared non-const, unlike the sibling overload / other queries
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h · _ff-algorithms_

```cpp
bool isValid(double seed_mz, double trace_tolerance)
```
- **Expectation:** An 'isValid(...)' predicate that only inspects the traces should be const (it returns bool and reads state). The sibling MassTrace::isValid() and other queries (getPeakCount, getRTBounds, getTheoreticalmaxPosition) are const.
- **Actual:** MassTraces::isValid(double, double) is declared and defined non-const even though it only reads this->size() and this->at(j).getAvgMZ(). This prevents calling it on a const MassTraces& and is inconsistent with the const MassTrace::isValid() in the same file.
- **Evidence:** Header: 'bool isValid(double seed_mz, double trace_tolerance);' (no const) vs 'bool isValid() const;' on MassTrace. Source :97-113 only reads members.
- **Fix:** Add const: 'bool isValid(double seed_mz, double trace_tolerance) const;'. This is a source-compatible change (adding const to a member function), though it changes the mangled name so it is technically ABI-affecting; provide as the ideal fix.

### [FEAT-24] FeatureFinderAlgorithmPickedHelperStructs::MassTrace::getAvgMZ — getAvgMZ() divides by summed intensity with no guard -> returns NaN on an empty/zero-intensity trace
`low` · `robustness` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h · _ff-algorithms_

```cpp
double getAvgMZ() const
```
- **Expectation:** A getter named getAvgMZ() returning the intensity-weighted average m/z should return a valid m/z, or signal failure, when there is nothing to average.
- **Actual:** It computes 'sum / intensities' where both accumulate over peaks. For an empty trace (or all-zero intensities) this is 0.0/0.0 = NaN, returned silently as if it were a valid m/z. Callers (e.g. MassTraces::isValid compares fabs(seed_mz - getAvgMZ()) <= tol) then compare against NaN, which is always false -- a silent logic bug.
- **Evidence:** Source :65-75: 'double sum = 0.0; double intensities = 0.0; for (...) { sum += mz*int; intensities += int; } return sum / intensities;' with no empty/zero guard.
- **Fix:** Guard against empty/zero intensity (return e.g. 0.0 or throw Precondition) and document. ABI-safe behavioral fix.
- **Verifier correction:** getAvgMZ() does compute sum/intensities with no empty/zero guard, so it returns NaN (0.0/0.0) on an empty or all-zero-intensity trace, and a NaN feeding the `<= tol` comparison in MassTraces::isValid silently evaluates false. This is provable. But the path is unreachable in the actual FeatureFinderAlgorithmPicked pipeline (traces are always seeded with ≥1 positive-intensity peak before getAvgMZ is called), so it is NOT a live high-severity silent-wrong-results bug. The real, defensible finding is an API-consistency/robustness gap: this exported, default-constructible getter returns NaN on empty input while sibling methods in the same struct (getTheoreticalmaxPosition, getRTBounds) deliberately throw Exception::Precondition on empty input and document it. Recommended fix (guard returning a sentinel or throwing Precondition, plus doc) is ABI-safe — signature unchanged.

### [FEAT-25] FeatureFinderAlgorithmPickedHelperStructs::MassTraces::getTheoreticalmaxPosition — getTheoreticalmaxPosition() returns a trace INDEX, not an m/z or RT 'position'
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h · _ff-algorithms_

```cpp
Size getTheoreticalmaxPosition() const
```
- **Expectation:** In a mass-spec codebase 'Position' overwhelmingly denotes an m/z or RT coordinate. getTheoreticalmaxPosition() reads as 'the position (m/z/RT) of the theoretical maximum'.
- **Actual:** It returns the std::vector INDEX (0-based slot) of the trace with the largest theoretical_int -- a container index, not a coordinate. Only the doc comment ('Returns the theoretical maximum trace index') reveals this; the name actively misleads.
- **Evidence:** Source :115-133 returns 'Size max = 0; ... max = i;' the index. Header doc: '@brief Returns the theoretical maximum trace index'.
- **Fix:** Add a clearer additive alias such as getMaxTheoreticalIntensityTraceIndex() and deprecate the old name, or at minimum strengthen the doc. Renaming the existing symbol is ABI-breaking; an additive alias is source-compatible.
- **Verifier correction:** getTheoreticalmaxPosition() does return a std::vector index (the slot of the trace with max theoretical_int), so the name is genuinely misleading given 'Position' usually means m/z/RT. However the claim understates the documentation: the doxygen @brief on the very declaration correctly reads 'Returns the theoretical maximum trace index', and the Size (integer) return type signals an index rather than a coordinate. Severity is therefore low (mild, documented, type-visible surprise) rather than medium. An additive alias like getMaxTheoreticalIntensityTraceIndex() would be source-compatible (abi_impact: none); renaming the existing symbol would be ABI-breaking.

### [FEAT-26] TraceFitter::getArea — getArea() is non-const while every sibling fit-result getter (getHeight/getCenter/getValue/getFWHM/getLowerRTBound/getUpperRTBound) is const
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/TraceFitter.h · _ff-algorithms_

```cpp
virtual double getArea() = 0;
```
- **Expectation:** getArea() is a pure read-only query of the fitted model (like getHeight/getFWHM, which are const) and should be callable on a const TraceFitter& / const std::shared_ptr<const TraceFitter>.
- **Actual:** getArea() is declared non-const in the base and overrides (GaussTraceFitter, EGHTraceFitter), yet the implementations only read fitted members and return a value -- GaussTraceFitter::getArea() is literally 'return 2.506628 * height_ * sigma_;'. The missing const is gratuitous and breaks const-correct call sites, inconsistent with all sibling getters in the same class.
- **Evidence:** Header TraceFitter.h:198 'virtual double getArea() = 0;' (no const) vs :137 'virtual double getHeight() const = 0;' etc. Source GaussTraceFitter.cpp:116-120 'double GaussTraceFitter::getArea() { return 2.506628 * height_ * sigma_; }'.
- **Fix:** Make getArea() const across the base and overrides to match the sibling getters. Source-compatible at call sites but changes the virtual's mangled name/vtable slot signature, so technically ABI-breaking; flag as the ideal fix.
- **Verifier correction:** Claim is accurate on the facts. Adjustment is to severity: it should be LOW, not implied-higher. The missing const is a real, gratuitous inconsistency (all six sibling getters are const; both getArea overrides only read fitted members), and a const TraceFitter* already exists in-tree (ElutionModelFitter::calculateFitQuality_), so the const-correctness break is concrete rather than theoretical. However, the defect is compile-time-loud and never yields wrong values, data loss, or crashes, so it is a mild surprise (low), not a silent-correctness or misuse-trap (medium). The recommendation to add const across base + GaussTraceFitter/EGHTraceFitter overrides is correct, and the author is right that it is source-compatible at call sites but ABI-breaking (changes the virtual's mangled name and vtable slot signature).

### [FEAT-10] TraceFitter::getArea — getArea() is non-const although it only reads/derives from fitted members; cannot be called on a const fitter
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/TraceFitter.h · _ff-fitters_

```cpp
virtual double getArea() = 0;
```
- **Expectation:** A read-only query that integrates the already-fitted model ('Area under the fitted RT model', '@return Integrated intensity') should be const, like the sibling getters getCenter()/getHeight()/getFWHM()/getValue() which are all const. A caller holding a const TraceFitter& expects to be able to ask for the area.
- **Actual:** getArea() is declared non-const in the base and in both overrides, yet GaussTraceFitter::getArea() just returns '2.506628 * height_ * sigma_' and EGHTraceFitter::getArea() only reads sigma_/tau_/height_ and local constants. Nothing is mutated. The missing const blocks const usage and is inconsistent with every other query in the same interface.
- **Evidence:** GaussTraceFitter.h: 'double getArea() override;' (no const) vs 'double getValue(double rt) const override;'. GaussTraceFitter.cpp: 'double GaussTraceFitter::getArea() { return 2.506628 * height_ * sigma_; }'. EGHTraceFitter.cpp getArea() reads only members and EPSILON_COEFS_.
- **Fix:** Add const to getArea() in TraceFitter and both subclasses. This is an ABI break (mangled name and vtable slot signature change), so if strict ABI must hold, add a new 'double getArea() const' overload and deprecate the non-const one. Mark abi_impact breaking for the ideal fix.
- **Verifier correction:** The surprise is genuine but its severity is low, not the elevated impact the framing suggests. getArea() returns correct values whenever invoked (all current callers hold non-const fitters), so there is no silent miscomputation, data loss, or crash. The only consequence is that a caller holding a `const TraceFitter&` (or const subclass) cannot call getArea() even though it can call every other query in the interface — a compile-time inconsistency, not a runtime defect. The ideal fix (add `const` to the base and both overrides) is an ABI break: it alters the mangled name and the virtual function's vtable slot signature. A source-compatible alternative exists (add a new `double getArea() const` overload and deprecate the non-const one), but the honest ABI impact of the primary recommended fix is breaking.

### [FEAT-11] GaussTraceFitter::getValue / getCenter / getArea / getSigma / getLowerRTBound — Trace-fitter query methods read uninitialized members when called before fit() (no fitted-state guard)
`low` · `silent-failure` · ABI: `source-compatible` · src/openms/source/FEATUREFINDER/GaussTraceFitter.cpp · _ff-fitters_

```cpp
double GaussTraceFitter::getValue(double rt) const
```
- **Expectation:** Calling getValue()/getCenter()/getArea()/getSigma()/getLowerRTBound() before a successful fit() either throws/asserts or returns a defined sentinel. A competent caller of a 'fitter' object that hasn't fit yet expects an error, not silent garbage.
- **Actual:** GaussTraceFitter members sigma_, x0_, height_ are not in-class initialized and the default constructor body is empty, so before fit() (or after fit() threw UnableToFit) these getters evaluate expressions over indeterminate doubles, e.g. getValue() = height_ * exp(-0.5*(rt-x0_)^2/sigma_^2) -> NaN/garbage with no signal. EGHTraceFitter has the same pattern (apex_rt_, height_, sigma_, tau_ uninitialized; '= default' ctor). The result is a silent, plausible-looking number.
- **Evidence:** GaussTraceFitter.h: 'double sigma_; double x0_; double height_;' (no initializers). GaussTraceFitter.cpp: 'GaussTraceFitter::GaussTraceFitter() { //setName(...) }' (empty). getValue: 'return height_ * exp(-0.5 * pow2(rt - x0_) / pow2(sigma_));'. EGHTraceFitter.h members likewise uninitialized; 'EGHTraceFitter::EGHTraceFitter() = default;'.
- **Fix:** In-class initialize the members (e.g. = NaN) and/or set a fitted_ flag in fit()/getOptimizedParameters_() and throw Exception::Precondition (or return NaN) from the queries when unfitted. In-class initializers and an added guard are source-compatible and ABI-safe.
- **Verifier correction:** The uninitialized-members / no-fitted-guard pattern is real and confirmed for both GaussTraceFitter and EGHTraceFitter, but it does not affect normal use: every in-tree caller constructs the fitter and calls fit() before any getter. The genuine surprise is limited to the misuse path — querying a freshly constructed (or post-UnableToFit) fitter reads indeterminate doubles and returns a plausible NaN/garbage value with no signal. Recommended fix (in-class NaN initialization and/or a fitted_ flag throwing Exception::Precondition from the queries) is source-compatible and ABI-safe.

### [FEAT-13] EmgScoring::getDefaults — EmgScoring::getDefaults() is a non-const factory call that ignores the object's configured params
`low` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/FEATUREFINDER/EmgScoring.h · _ff-fitters_

```cpp
Param getDefaults()
```
- **Expectation:** A method named getDefaults() on a configurable object is expected to be a const accessor returning the default parameter set, consistent with DefaultParamHandler::getDefaults() used throughout OpenMS.
- **Actual:** It is non-const and returns 'EmgFitter1D().getDefaults()' - it constructs a brand-new EmgFitter1D each call and returns ITS defaults, completely ignoring fitter_emg1D_params_ that the user may have installed via setFitterParam(). So getDefaults() never reflects this object's state and could be const. The pairing setFitterParam()/getDefaults() (instead of getParameters()) is also asymmetric: there is no getter for the params you set.
- **Evidence:** Header: 'Param getDefaults() { return EmgFitter1D().getDefaults(); }' and 'void setFitterParam(const Param& param) { fitter_emg1D_params_ = param; }' with member 'Param fitter_emg1D_params_;' but no getter for it.
- **Fix:** Make getDefaults() const, and add a symmetric getFitterParam()/getParameters() accessor returning fitter_emg1D_params_ so callers can read back what they set. Adding const to an inline header method and adding a getter are source-compatible.
- **Verifier correction:** EmgScoring::getDefaults() is non-const and returns the defaults of a freshly constructed EmgFitter1D, ignoring any params installed via setFitterParam(). This is correctly documented ("Get default params"; "use getDefaults to see what you can set") and matches OpenMS convention (getDefaults = defaults, getParameters = configured state), so it is not a misleading name and causes no wrong results — fitRT_() does apply fitter_emg1D_params_. The real, minor issues are: (a) the method could be const but is not, and (b) the API is asymmetric (setFitterParam setter with no getter for fitter_emg1D_params_). Recommend adding const and a symmetric getFitterParam()/getParameters() accessor; both changes are source-compatible (recompile-only, non-virtual inline method).

### [FEAT-29] SeedListGenerator::generateSeedList(PeptideIdentificationList&, SeedList&, bool) — A 'generate seed list' method silently re-sorts the caller's peptide hits when use_peptide_mass is true
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/SeedListGenerator.h · _ff-identification_

```cpp
void generateSeedList(PeptideIdentificationList& peptides, SeedList& seeds, bool use_peptide_mass = false)
```
- **Expectation:** A generateSeedList(...) that takes peptides to read positions from would be expected to leave the input identifications untouched (it is conceptually a read/convert operation).
- **Actual:** When use_peptide_mass is true the method calls pep.sort() on every PeptideIdentification, reordering the caller's hit vectors in place. The non-const reference and a one-line note in the doc hint at it, but the verb 'generate' and the surrounding read-only overloads make the mutation easy to miss.
- **Evidence:** if (!pep.getHits().empty() && use_peptide_mass) { pep.sort(); const PeptideHit& hit = pep.getHits().front(); ... }; header note: 'The peptide hits in @p peptides will be sorted if @p use_peptide_mass is true.'
- **Fix:** Take peptides by const-ref and sort a local copy of the relevant hit (or document the mutation prominently in the brief). If keeping the signature, keep the existing doc note but consider a const overload. Changing the parameter to const& is source-compatible for read-only callers but is an ABI/signature change.
- **Verifier correction:** The side effect is real but documented in the method's own doxygen brief, so it is not truly 'hidden', and the effect is a benign in-place stable_sort of hits by score (best first) — recoverable, no data loss/crash. Severity is low, not high. ABI impact of the method as-is is 'none' (no change proposed/needed); the recommendation to switch to const& would be an ABI/signature change but is not the current state, so abi_impact for the claim as filed is none.

### [FEAT-32] Internal::FFIDAlgoExternalIDHandler::classifyFeaturesWithSVM — classifyFeaturesWithSVM overwrites each feature's overallQuality with the SVM probability
`low` · `hidden-side-effect` · ABI: `none` · src/openms/source/FEATUREFINDER/FFIDAlgoExternalIDHandler.cpp · _ff-identification_

```cpp
void classifyFeaturesWithSVM(FeatureMap& features, const Param& param)
```
- **Expectation:** A 'classify features using SVM' step is expected to attach class/probability annotations (predicted_class, predicted_probability), not to clobber the pre-existing overallQuality score computed upstream by OpenSWATH.
- **Actual:** After predicting, it sets features[i].setOverallQuality(prob_positive) for every feature, silently replacing the OpenSWATH overall quality. The source even flags this as questionable ('@TODO: store previous (OpenSWATH) overall quality in a meta value?'). The header brief 'Classify features using SVM' gives no hint that the quality score domain changes (OpenSWATH score -> SVM probability in [0,1]).
- **Evidence:** features[i].setMetaValue("predicted_probability", prob_positive); // @TODO: store previous (OpenSWATH) overall quality in a meta value?\n      features[i].setOverallQuality(prob_positive);
- **Fix:** Document that overallQuality is repurposed to the SVM positive-class probability after classification; ideally preserve the prior score in a meta value as the TODO suggests. Doc change is ABI-safe.
- **Verifier correction:** classifyFeaturesWithSVM does overwrite each feature's overallQuality with the SVM positive-class probability (prob_positive in [0,1]) for every feature, discarding the prior score with no meta value preserving it, and the header brief gives no hint of this domain change — so the hidden side effect is real. But the discarded score is the FeatureFinderIdentification (FFID) internal overall quality, not an OpenSWATH quality as titled; the probability remains available in the "predicted_probability" meta value; and the overwrite is intentional, with downstream filterClassifiedFeatures depending on overallQuality == predicted_probability. The class is Internal::FFIDAlgoExternalIDHandler. Recommendation stands: document the repurposing on the public brief and ideally store the prior quality in a meta value per the TODO (doc/behavior change, ABI-safe).

### [FEAT-33] PeakWidthEstimator::getPeakWidth — getPeakWidth is a pure read but is non-const, blocking use on a const PeakWidthEstimator
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/PeakWidthEstimator.h · _ff-identification_

```cpp
double getPeakWidth(double mz)
```
- **Expectation:** A getter that only evaluates a precomputed B-spline at an m/z should be const, callable on a const PeakWidthEstimator.
- **Actual:** getPeakWidth mutates nothing — it only calls bspline_->eval(...) (which is itself const) and compares mz against the stored range — yet it is declared non-const, so it cannot be called through a const reference/pointer to the estimator.
- **Evidence:** double getPeakWidth(double mz); with body width = (*bspline_).eval(mz); where BSpline2d::eval is declared 'double eval(const double x) const;'
- **Fix:** Make getPeakWidth const (double getPeakWidth(double mz) const). This is source-compatible for existing callers but changes the mangled signature, so it is technically an ABI break for that single method.
- **Verifier correction:** getPeakWidth is genuinely a pure read declared non-const and could/should be const. But the practical impact is minimal: it only blocks const-context use, which the compiler rejects loudly at compile time, and no current caller is even const, so nothing is actually blocked. This is a mild ergonomic const-correctness improvement (low severity), not a behavior bug. Making it const is source-compatible for callers but an ABI break for that single method's mangled symbol.

### [FEAT-15] IsotopeModel::getOffset / getFormula — getOffset() and getFormula() are read-only accessors but declared non-const (sibling getCenter()/getCharge() are const)
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/IsotopeModel.h · _ff-isotope-models_

```cpp
CoordinateType getOffset(); EmpiricalFormula getFormula();
```
- **Expectation:** A getter named getOffset()/getFormula() that only reads state should be const, like the other accessors in the same class (getCharge() const, getCenter() const) and in InterpolationModel. Callers expect to call them on a `const IsotopeModel&`.
- **Actual:** getOffset() merely returns getInterpolation().getOffset() (IsotopeModel.cpp:229-232) and mutates nothing, yet is declared non-const (IsotopeModel.h:65). getFormula() only reads mean_, charge_, averagine_ (IsotopeModel.cpp:68-100) and constructs a return value, yet is also non-const (IsotopeModel.h:68). This breaks const-correctness and forces callers to hold a non-const reference just to read the offset/formula. ExtendedIsotopeModel::getOffset() has the identical defect.
- **Evidence:** IsotopeModel.h:65 `CoordinateType getOffset();` (no const) and :68 `EmpiricalFormula getFormula();`; cf. :54 `UInt getCharge() const;` and :80 `getCenter() const`. Bodies IsotopeModel.cpp:229-232 and 68-100 perform no mutation.
- **Fix:** Add `const` to getOffset() and getFormula() (and ExtendedIsotopeModel::getOffset()). Adding const to a member function changes the mangled name, so it is technically an ABI break for these symbols; do it in a major/ABI-bump release, or add const-qualified overloads alongside the existing ones.
- **Verifier correction:** Confirmed: IsotopeModel::getOffset() (h:65) and getFormula() (h:68), plus ExtendedIsotopeModel::getOffset() (h:62), are non-const read-only accessors, asymmetric with const siblings getCharge()/getCenter() in the same classes. Adding const is feasible (getInterpolation() const and LinearInterpolation::getOffset() const exist). Severity is LOW, not medium/high: the defect manifests as a loud compile error when calling on a const reference, never as silent wrong behavior. ABI: mutating the existing symbols to const changes the mangled name and is an ABI break for OPENMS_DLLAPI shared libs; a const-qualified overload added alongside would instead be source-compatible.

### [FEAT-18] ElutionModelFitter::calculateFitQuality_ — Quality calculation divides by per-point model_value and by total_weights with no zero/empty guard, can return Inf/NaN
`low` · `silent-failure` · ABI: `source-compatible` · src/openms/source/FEATUREFINDER/ElutionModelFitter.cpp · _ff-isotope-models_

```cpp
double calculateFitQuality_(const TraceFitter* fitter, const MassTraces& traces)
```
- **Expectation:** A 'calculate quality of model fit (mean relative error)' helper should return a finite error; a competent caller expects it to handle degenerate cases (no points in range, zero model value).
- **Actual:** calculateFitQuality_ accumulates `mre += diff / model_value;` (line 85) where model_value = fitter->getValue(rt) can be 0 at the tails, producing Inf; and returns `mre / total_weights` (line 90) where total_weights is 0 if no peak falls in [rt_start, rt_end], producing NaN/Inf. No guard exists. The returned value is stored as feature meta 'model_error' (line 146) and silently propagates a non-finite error.
- **Evidence:** ElutionModelFitter.cpp:85 'mre += diff / model_value;', :90 'return mre / total_weights;'; total_weights starts at 0.0 (line 68) and is only incremented when a point is in range (line 86).
- **Fix:** Guard against model_value==0 and total_weights==0 (return a sentinel such as -1.0 already used elsewhere in the class, or skip such points) and document the contract. Internal helper, so the fix is source/ABI compatible.
- **Verifier correction:** The reachable defect is total_weights==0 -> return 0.0/0.0 = NaN when no data point falls within [rt_start, rt_end] (or when summed theoretical_int is 0), not the claimed Inf via model_value==0 (that 0 only occurs in EGH's cut-off tail, which lies outside the [getLowerRTBound, getUpperRTBound] window the loop restricts to; Gauss getValue is never exactly 0). The non-finite value is stored solely as the informational, write-only meta value 'model_error' (never read by any downstream filter or the model_status validity check), so it does not silently corrupt scientific results — impact is a misleading meta value, not data loss/crash. Recommended fix: guard total_weights==0 (and optionally model_value==0) by returning the existing -1.0 sentinel; private member function, so the change is source/ABI compatible.

### [FEAT-34] FeatureFindingMetabo::run — run() silently reorders the caller's input_mtraces vector in place
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFindingMetabo.h · _ff-metabo_

```cpp
void run(std::vector<MassTrace>& input_mtraces, FeatureMap& output_featmap, std::vector<std::vector<OpenMS::MSChromatogram>>& output_chromatograms)
```
- **Expectation:** A parameter named input_mtraces passed to the 'main method' is read-only input; the caller's vector order is preserved.
- **Actual:** run() sorts the caller's vector in place by centroid m/z before processing, so the caller's container is permanently reordered as a side effect.
- **Evidence:** FeatureFindingMetabo.cpp:885 `std::sort(input_mtraces.begin(), input_mtraces.end(), CmpMassTraceByMZ());` operating directly on the non-const reference parameter `std::vector<MassTrace>& input_mtraces`.
- **Fix:** Document the reordering in the header (additive, ABI-safe). Ideally sort a local copy or restore order, but that changes behavior; at minimum the header doc must state 'input_mtraces is sorted by m/z as a side effect'. If a non-mutating contract is wanted, add a const overload that copies.
- **Verifier correction:** run() sorts the caller's input_mtraces vector in place by centroid m/z (CmpMassTraceByMZ) before processing, permanently reordering the caller's container. This is real and undocumented in the header. But because the parameter is already a non-const std::vector<MassTrace>& (signaling possible mutation), the function's own outputs are unaffected, and the only in-tree caller relies solely on .size() afterward, the practical impact is a mild documentation surprise rather than a high-severity hazard. Fix is additive/ABI-safe: state in the header doc that input_mtraces is sorted by m/z as a side effect (optionally a const overload that copies if a non-mutating contract is desired).

### [FEAT-39] FeatureHypothesis::getMonoisotopicFeatureIntensity / getSummedFeatureIntensity — smoothed flag has no default here but does on sibling getters; default value lives uselessly in the .cpp
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/FEATUREFINDER/FeatureFindingMetabo.h · _ff-metabo_

```cpp
double getMonoisotopicFeatureIntensity(bool) const;  double getSummedFeatureIntensity(bool) const;
```
- **Expectation:** Consistent with sibling intensity getters getAllIntensities(bool smoothed = false) and getMaxIntensity(bool smoothed = false), these would default smoothed to false and name the parameter.
- **Actual:** The header declares an unnamed `bool` with no default, so callers MUST pass it explicitly, unlike the sibling getters. Worse, the .cpp definition writes `bool smoothed = false`, a default-argument-on-a-definition that is invisible to every other translation unit and only misleads readers of the .cpp.
- **Evidence:** Header lines 81-82 `double getMonoisotopicFeatureIntensity(bool) const; double getSummedFeatureIntensity(bool) const;` (no default) vs lines 63/85 `getAllIntensities(bool smoothed = false)`, `getMaxIntensity(bool smoothed = false)`. FeatureFindingMetabo.cpp:34 `getMonoisotopicFeatureIntensity(bool smoothed = false)` and :44 `getSummedFeatureIntensity(bool smoothed = false)` carry the (ineffective) default in the definition.
- **Fix:** Add `bool smoothed = false` and a parameter name in the header to match the sibling getters, and remove the default from the .cpp (it belongs in the declaration). Adding a default argument is source-compatible and ABI-safe.

### [FEAT-40] MassTraceDetection::run — max_traces default of 0 silently means 'unlimited', not 'extract zero traces'
`low` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/MassTraceDetection.h · _ff-metabo_

```cpp
void run(const PeakMap&, std::vector<MassTrace>&, const Size max_traces = 0)
```
- **Expectation:** max_traces = 0 caps the output at zero traces (a literal reading), or the meaning of 0 is documented.
- **Actual:** 0 is a sentinel for 'no limit': the cap is applied only when `max_traces > 0`. A caller passing/leaving 0 gets every trace, the opposite of the literal meaning, and the header gives no hint.
- **Evidence:** MassTraceDetection.cpp:681 `if (max_traces > 0 && found_masstraces.size() == max_traces) { break; }`. Header line 73 documents nothing about the 0 sentinel.
- **Fix:** Document the 0 = unlimited convention in the header (additive, ABI-safe). Optionally provide a named constant or a clearer overload without the parameter for the unlimited case.
- **Verifier correction:** The undocumented `max_traces = 0` does mean "unlimited" (cap applied only when max_traces > 0, MassTraceDetection.cpp:681), confirmed. But this is the standard C++ idiom for a maximum-count parameter and is the sensible/expected default (the literal "zero output" reading would make the overload useless and is not what any caller assumes). The genuine issue is narrow: the header gives no hint that 0 = unlimited. Severity is low (mild documentation gap on an idiomatic sentinel), not a silently-wrong-results trap. Fix is a one-line header doc comment; ABI-safe.

### [FEAT-42] ElutionPeakDetection::computeMassTraceNoise / computeMassTraceSNR / computeApexSNR — Read-only S/N and noise compute methods are not declared const
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/ElutionPeakDetection.h · _ff-metabo_

```cpp
double computeMassTraceNoise(const MassTrace&);  double computeMassTraceSNR(const MassTrace&);  double computeApexSNR(const MassTrace&);
```
- **Expectation:** Methods named compute* that take a const MassTrace& and only read state are const member functions, so they can be called on a const ElutionPeakDetection.
- **Actual:** All three are non-const members although they neither modify *this nor the passed-in (const) MassTrace; they only read trace data. This blocks calling them through a const reference and surprises callers who treat them as pure computations.
- **Evidence:** ElutionPeakDetection.cpp:48 `double ElutionPeakDetection::computeMassTraceNoise(const MassTrace& tr)` (no const qualifier), :69 computeMassTraceSNR, :86 computeApexSNR; none modify any member.
- **Fix:** Add `const` to all three method declarations/definitions. This is source-compatible for callers but changes the mangled names, so it is technically an ABI break; batch it with the next ABI-breaking release.
- **Verifier correction:** The claim is accurate as to the code: all three methods (computeMassTraceNoise, computeMassTraceSNR, computeApexSNR) are non-const yet only read state and could be declared const. Severity is low rather than higher: it is a mild const-correctness/ergonomic surprise with no correctness, data-loss, or crash consequences, and the methods are normally invoked on a non-const object internally. Adding const is source-compatible for callers but is an ABI break (mangled-name change), so it should be batched with the next ABI-breaking release. The header path in the claim is mistyped (path segment doubled); the actual file is src/openms/include/OpenMS/FEATUREFINDER/ElutionPeakDetection.h.

### [FEAT-1] BaseModel::BaseModel() — Default constructor leaves cut_off_ uninitialized; getCutOff()/isContained() read garbage
`low` · `uninitialized-member-risk` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/BaseModel.h · _ff-models1d_

```cpp
BaseModel() : DefaultParamHandler("BaseModel")
```
- **Expectation:** After default-constructing a model, getCutOff() returns the documented default (0.0) and isContained() compares against a well-defined cutoff.
- **Actual:** The ctor only does defaults_.setValue("cutoff", 0.0, ...). It never calls defaultsToParam_() nor updateMembers_(). The base DefaultParamHandler(name) ctor does NOT copy defaults_ into param_ (confirmed in DefaultParamHandler.cpp: it only sets error_name_/flags). The member IntensityType cut_off_ has no in-class initializer. So a freshly constructed BaseModel has an indeterminate cut_off_ until some external call triggers updateMembers_(). getCutOff() then returns garbage and isContained(pos) compares getIntensity(pos) >= <garbage>.
- **Evidence:** Ctor body: 'defaults_.setValue("cutoff", 0.0, ...)' with no defaultsToParam_(); member decl 'IntensityType cut_off_;' (no initializer, line 107); updateMembers_() is the only writer of cut_off_. Contrast InterpolationModel ctor which DOES call defaultsToParam_() (InterpolationModel.h:46).
- **Fix:** Add an in-class initializer 'IntensityType cut_off_{0.0};' (ABI-neutral, fixes the read-before-init), and/or call defaultsToParam_() at the end of the BaseModel ctor for consistency with InterpolationModel. Both are source/ABI-compatible.
- **Verifier correction:** cut_off_ is indeed left without an in-class initializer and the BaseModel ctor does not initialize it. But BaseModel is abstract and every concrete subclass calls defaultsToParam_() (→ updateMembers_()) in its own ctor, setting cut_off_=0.0 before any use. So getCutOff()/isContained() do NOT return garbage for any reachable, fully-constructed object; the issue is a latent uninitialized-member / future-maintenance fragility, severity low. Fix recommendation (in-class initializer IntensityType cut_off_{0.0};, ABI-neutral) remains valid and worthwhile.

### [FEAT-2] Fitter1D::fit1d — Base fit1d() throws NotImplemented instead of being pure virtual
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/Fitter1D.h · _ff-models1d_

```cpp
virtual QualityType fit1d(const RawDataArrayType& range, std::unique_ptr<InterpolationModel>& model)
```
- **Expectation:** An abstract base 'model fitter' whose central operation is fitting would declare fit1d() pure virtual (= 0), so a non-overriding subclass fails to compile, and the doc comment 'return interpolation model' implies it returns a quality/model.
- **Actual:** fit1d is a concrete virtual that unconditionally 'throw Exception::NotImplemented(...)' (Fitter1D.cpp:61-64). A caller holding a Fitter1D& that forgot to override gets a runtime throw rather than a compile-time error. The header doc '/// return interpolation model' also mislabels it: it returns a QualityType (correlation/quality), the model is an out-param.
- **Evidence:** Header: 'virtual QualityType fit1d(const RawDataArrayType& /* range */, std::unique_ptr<InterpolationModel>& /* model */);' (not =0). Impl: 'throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);'
- **Fix:** Prefer making fit1d pure virtual (= 0) to surface the contract at compile time; if ABI/legacy subclass concerns block that, at minimum fix the doc comment to '/// fit the model and return the fit quality; resulting model returned via out-param'. Pure-virtualizing is source-breaking for any subclass that relied on inheriting the throwing default (rare).
- **Verifier correction:** Base fit1d() is a concrete virtual that throws Exception::NotImplemented rather than being pure virtual — a deliberate, tested throwing-default idiom (asserted in Fitter1D_test.cpp:95), used because intermediate classes (LevMarqFitter1D, MaxLikeliFitter1D) inherit it without overriding. It fails loudly and recoverably at runtime, not silently. The header doc `/// return interpolation model` is stale: the function returns QualityType (the fit quality), and the model is returned via the out-param. Recommended fix is purely cosmetic/clarity: correct the doc comment to e.g. `/// fit the model to the data; returns the fit quality, model returned via out-param`. Pure-virtualizing is optional and would be source-breaking for the two intermediate subclasses that rely on the inherited default. As-is there is no ABI impact.

### [FEAT-3] InterpolationModel::getScalingFactor / scaling semantics — 'scaling' is documented as area-under-model, not a multiplicative intensity scale
`low` · `documentation-naming-inconsistency` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/InterpolationModel.h · _ff-models1d_

```cpp
CoordinateType getScalingFactor() const
```
- **Expectation:** getScalingFactor()/setScalingFactor() controls a multiplicative factor on intensities (the param is even named 'intensity_scaling' = 'Scaling factor used to adjust the model distribution to the intensities of the data').
- **Actual:** The doxygen on getScalingFactor states 'A scaling factor of @p scaling means that the area under the model equals @p scaling. Default is 1.' i.e. it is an absolute area normalization, not a free multiplier. The param description ('intensity_scaling') and the getter name suggest a plain intensity multiplier, which contradicts the area-equals-scaling definition.
- **Evidence:** getScalingFactor doc: 'A scaling factor of @p scaling means that the area under the model equals @p scaling.' vs param registration: defaults_.setValue("intensity_scaling", 1.0, "Scaling factor used to adjust the model distribution to the intensities of the data").
- **Fix:** Reconcile the doc and the param description so the unit (area vs multiplier) is unambiguous. Documentation-only; ABI-safe.
- **Verifier correction:** The evidence is accurate, but the framing overstates it as a unit-or-index hazard. The authoritative behavior (area under the model equals the scaling factor) is correctly documented in-place on getScalingFactor() and is consistently implemented in all four derived setSamples() methods. The only defect is that the param name "intensity_scaling" and its registration description ("...adjust the model distribution to the intensities of the data") suggest a plain multiplicative intensity scale, which is loosely worded relative to the precise area-normalization semantics. This is a low-severity documentation/naming inconsistency (recommend aligning the param description with the area-equals-scaling definition), not a silently-wrong unit/index bug. ABI-safe / documentation-only.

### [FEAT-4] InterpolationModel::getCenter / setSamples — getCenter() and setSamples() are concrete virtuals that throw NotImplemented
`low` · `surprising-throw` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/InterpolationModel.h · _ff-models1d_

```cpp
virtual CoordinateType getCenter() const; virtual void setSamples()
```
- **Expectation:** getCenter() (a const accessor for the model center) and setSamples() should be pure virtual in an abstract base, so subclasses must implement them; a const 'get the center' accessor should never throw.
- **Actual:** Both are concrete virtuals whose body is 'throw Exception::NotImplemented(...)'. A caller treating InterpolationModel polymorphically and calling getCenter() on a subclass that forgot to override gets a runtime throw from a const getter, which is exactly the surprising-throw-from-a-getter pattern.
- **Evidence:** 'virtual CoordinateType getCenter() const { throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }' and 'virtual void setSamples() { throw Exception::NotImplemented(...); }'.
- **Fix:** Make getCenter() and setSamples() pure virtual (= 0) so the missing-override is a compile error. This is source-breaking only for subclasses that intentionally inherited the throw (none expected). If ABI must hold, at least document that these throw if not overridden.
- **Verifier correction:** getCenter() and setSamples() are concrete virtuals that throw NotImplemented instead of being pure virtual; this deviates from the parent BaseModel's pure-virtual convention, and getCenter() is an undocumented const accessor that throws. But the practical impact is low: all in-tree concrete subclasses override getCenter(), the throw is a loud recoverable exception (no silent/wrong results), so no reasonable use is affected. The recommendation must be split: getCenter() can be made pure virtual safely, but setSamples() must NOT be made `= 0` — IsotopeModel intentionally inherits the throwing no-arg setSamples() (overloading only setSamples(EmpiricalFormula) + a using-declaration), so making it pure virtual would break compilation of IsotopeModel. The safe minimal fix is to make getCenter() pure virtual and to document that the base setSamples() throws if not overridden.

### [FEAT-5] ModelDescription::getName() / getParam() non-const overloads — getName()/getParam() return mutable references to internals under a 'get' name
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/ModelDescription.h · _ff-models1d_

```cpp
std::string& getName(); Param& getParam()
```
- **Expectation:** A 'getName()'/'getParam()' accessor returns a read-only view; mutation goes through setName()/setParam().
- **Actual:** Both have a non-const overload returning a mutable reference to the internal name_/parameters_ ('/// Mutable access to the model name'). Callers can silently mutate the stored state through what reads like a getter (e.g. desc.getName() += "x";), bypassing setName/setParam and any future validation. setName/setParam already exist, so the mutable getters are redundant and surprising.
- **Evidence:** 'std::string & getName() { return name_; }' and 'Param & getParam() { return parameters_; }', each labeled 'Mutable access', alongside existing setName()/setParam().
- **Fix:** Prefer removing the non-const getName()/getParam() overloads in favor of the existing setters (source-breaking for code that mutates via the getter). If ABI/source stability is required, leave them but document the mutation channel clearly. Honest ideal fix: drop the mutable overloads.
- **Verifier correction:** getName()/getParam() do expose non-const mutable-reference overloads that read like getters and bypass setName()/setParam(), redundant with the existing setters — a real const-correctness nit. But impact is low, not high: the mutation channel is documented inline ('/// Mutable access'), there are no non-test callers anywhere in the tree, and misuse cannot produce silently-wrong results, data loss, or crashes. The recommended fix (drop the mutable overloads) is reasonable and source-breaking only for code mutating via the getter — none exists outside the unit test, whose 'std::string& getName()' / 'Param& getParam()' sections would need trivial edits.

### [FEAT-6] ModelDescription::ModelDescription(const BaseModel*) — Single-argument non-explicit constructor enables implicit BaseModel*->ModelDescription conversion
`low` · `implicit-conversion` · ABI: `source-compatible` · src/openms/include/OpenMS/FEATUREFINDER/ModelDescription.h · _ff-models1d_

```cpp
ModelDescription(const BaseModel * model)
```
- **Expectation:** A converting constructor from a raw pointer should be explicit, so a BaseModel* does not silently become a ModelDescription in overload resolution / assignment contexts.
- **Actual:** The constructor 'ModelDescription(const BaseModel * model)' is non-explicit, so any BaseModel* implicitly converts to a ModelDescription<D>. Combined with operator== taking a ModelDescription, a stray pointer can implicitly construct a temporary and compare/assign unexpectedly.
- **Evidence:** '/// constructor provided for convenience\n ModelDescription(const BaseModel * model) : name_(model->getName()), parameters_(model->getParameters()) {}' (no explicit).
- **Fix:** Mark the constructor 'explicit'. Source-compatible for normal '{}'/direct construction; breaks only code relying on the implicit conversion (rare, and that code is the bug this fix prevents).
- **Verifier correction:** The non-explicit single-argument constructor `ModelDescription(const BaseModel*)` does enable an implicit BaseModel*->ModelDescription<D> conversion, but the conversion is narrower than claimed: because BaseModel is not a template and D cannot be deduced from a BaseModel*, the implicit conversion only fires in contexts where D is already fixed (e.g. comparing a stray pointer against an existing ModelDescription<2> via operator==). The constructor dereferences the pointer immediately, so a null/dangling pointer crashes rather than silently mis-comparing. No code in the repository relies on this conversion, and ModelDescription is essentially unused legacy code. Recommendation stands: mark the constructor `explicit`. This is source-compatible (no ABI change; would only reject the currently-nonexistent code that abuses the implicit conversion).

### [FEAT-7] BaseModel::getSamples(std::ostream&) — getSamples(ostream&) is non-const while getSamples(SamplesType&) it delegates to is const
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/FEATUREFINDER/BaseModel.h · _ff-models1d_

```cpp
virtual void getSamples(std::ostream& os)
```
- **Expectation:** An overload that only reads the model and writes samples to a stream should be const, matching the sibling 'virtual void getSamples(SamplesType& cont) const'.
- **Actual:** getSamples(std::ostream&) is declared non-const even though its body only calls the const getSamples(SamplesType&) and streams the result. This forces callers to hold a non-const BaseModel just to print, and is inconsistent with its own const sibling overload.
- **Evidence:** 'virtual void getSamples(SamplesType& cont) const = 0;' (const) vs 'virtual void getSamples(std::ostream& os)' (non-const) whose body is 'SamplesType samples; getSamples(samples); for (...) os << sample;'.
- **Fix:** Add 'const' to getSamples(std::ostream&). This changes the function signature/mangled name (it is virtual), so it is technically ABI-breaking for the vtable slot; for a library that values ABI, do it at the next ABI break. Source-compatible for all callers.
- **Verifier correction:** getSamples(std::ostream&) is genuinely non-const while it only reads the model (it merely calls the const getSamples(SamplesType&) and streams the samples), making it inconsistent with its const sibling. The defect is real but low-severity: it causes no incorrect behavior, only the inability to print from a const BaseModel. Adding const is source-compatible for all callers/subclasses (no override exists), but ABI-breaking because it changes the virtual function's mangled name and vtable slot.

### [FEAT-44] MultiplexDeltaMassesGenerator::getDeltaMassesList — Non-const getDeltaMassesList() returns a COPY by value while the const overload returns a reference; mutating the result is silently discarded
`low` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h · _ff-multiplex-a_

```cpp
std::vector<MultiplexDeltaMasses> getDeltaMassesList(); / const std::vector<MultiplexDeltaMasses>& getDeltaMassesList() const;
```
- **Expectation:** When a class exposes a const ref-returning getter and a non-const sibling of the same name, the non-const one returns a non-const reference so the caller can modify internal state in place (the usual OpenMS get-by-ref idiom).
- **Actual:** The non-const overload returns 'std::vector<MultiplexDeltaMasses>' BY VALUE (a copy of delta_masses_list_), whereas the const overload returns 'const std::vector<...>&'. A caller writing 'generator.getDeltaMassesList()[0] = ...;' or '.push_back(...)' on a non-const object modifies a throwaway temporary; the generator is unchanged. Overload resolution silently picks the value-returning copy for non-const objects.
- **Evidence:** Header lines 101 & 106: 'std::vector<MultiplexDeltaMasses> getDeltaMassesList();' vs 'const std::vector<MultiplexDeltaMasses>& getDeltaMassesList() const;'. Impl (MultiplexDeltaMassesGenerator.cpp:525-533) returns 'delta_masses_list_' by value in the non-const path. The pyOpenMS binding even force-casts to the const overload: 'return self.getDeltaMassesList(); }, nb::rv_policy::reference_internal' on a const self (bind_misc.cpp:1773).
- **Fix:** Either return a non-const reference from the non-const overload (matches the const sibling and the OpenMS by-ref idiom) or drop the non-const overload entirely so all reads go through the const ref. Changing the return type changes the symbol/ABI; prefer adding a returns-reference accessor and deprecating the value-returning one.
- **Verifier correction:** The asymmetry is real (non-const overload returns by value, const overload by reference) but the 'silently discarded mutation / data loss for reasonable use' framing is not supported: no caller mutates through the getter and the doc never promises in-place mutation. The non-const overload is effectively a redundant copying accessor; the worst realistic outcomes are a wasted copy or a compile/lifetime error (e.g. binding `auto&` to the returned temporary), not silent wrong results. This is a low-severity API smell — fix by either making the non-const overload return a non-const reference (ABI-breaking) or dropping it so all reads go through the const ref.

### [FEAT-46] operator<(const MultiplexDeltaMasses&, const MultiplexDeltaMasses&) — operator< returns true when the left operand has MORE delta masses (sorts larger multiplets first), contradicting a < b ascending intuition
`low` · `misleading-name` · ABI: `none` · src/openms/source/FEATUREFINDER/MultiplexDeltaMasses.cpp · _ff-multiplex-a_

```cpp
bool operator<(const MultiplexDeltaMasses &dm1, const MultiplexDeltaMasses &dm2)
```
- **Expectation:** A free operator< induces a normal ascending ordering: a < b implies a sorts before b because a is 'smaller'. With std::sort the smallest element ends up first.
- **Actual:** When the two patterns have different sizes, the comparator returns 'dm1.size() > dm2.size()', i.e. a pattern with MORE delta masses is considered 'less than' one with fewer. So std::sort places larger multiplets (quadruplets) before singlets even though they are numerically 'bigger'. This is intentional (complete multiplets searched first) but the < symbol actively inverts size ordering, which surprises any caller using this type in an ordered container or sort and reasoning about magnitude.
- **Evidence:** MultiplexDeltaMasses.cpp:65-69: 'if (dm1.getDeltaMasses().size() != dm2.getDeltaMasses().size()) { // Search first for complete multiplets, then knock-out cases. return (dm1.getDeltaMasses().size() > dm2.getDeltaMasses().size()); }'. Declared as a public free function in the header (MultiplexDeltaMasses.h:97).
- **Fix:** Keep the algorithmic intent but rename the concept: provide a named comparator (e.g. struct LessForSearchPriority / compareBySearchPriority) instead of overloading the bare operator<, or document the inverted-size semantics prominently on the operator< declaration in the header (currently undocumented there). If the type is only ever sorted for search priority, a named functor avoids surprising any future ordered-container use. Adding a named comparator is purely additive.
- **Verifier correction:** operator< is a valid strict weak ordering, not a correctness bug. It inverts ordering only on the PRIMARY (size) key — larger multiplets sort first — while the secondary key (relative mass-shift) is normal ascending. The inverted intent is documented in the function body and at the std::sort call site (MultiplexDeltaMassesGenerator.cpp:339-340), but NOT on the header declaration (line 97). The only real surprise/risk is a future caller using this internal type in an ordered container or a magnitude-based sort; current usage is correct and intentional. Fix is purely additive (document on the declaration or add a named compareBySearchPriority functor); severity low, no ABI impact.

### [FEAT-47] MultiplexClustering::MultiplexClustering(const MSExperiment&, double, bool, double) — Bare bool 'mz_tolerance_unit' selects ppm vs Da; unreadable at the call site and undocumented direction
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/MultiplexClustering.h · _ff-multiplex-a_

```cpp
MultiplexClustering(const MSExperiment& exp, double mz_tolerance, bool mz_tolerance_unit, double rt_typical)
```
- **Expectation:** A unit selector for an m/z tolerance is an enum or is at least unambiguous; 'true' should not silently mean 'ppm'.
- **Actual:** The 4th parameter is a naked bool whose only documentation is in the param block ('unit for mz_tolerance, ppm (true), Da (false)'). At a call site 'MultiplexClustering(exp, 5.0, true, 60.0)' it is impossible to tell whether 5.0 is 5 ppm or 5 Da, and swapping the polarity silently changes the grid spacing from a relative ppm step to an absolute Da step (impl branches on it for both grid spacing and rt_scaling).
- **Evidence:** Header line 68 plus doc 'mz_tolerance_unit unit for mz_tolerance, ppm (true), Da (false)' (line 63). Impl (MultiplexClustering.cpp:106-119 and 139-146) branches: 'if (mz_tolerance_unit) { ... mz * (1 + scaling*mz_tolerance/1000000) ...} else { ... mz + scaling*mz_tolerance ...}'.
- **Fix:** Add an overload or follow-up constructor taking a small enum (e.g. MZUnit::PPM / MZUnit::DA) instead of bool; keep the bool ctor for ABI but document the polarity on the declaration. Additive enum overload is non-breaking.
- **Verifier correction:** Corrected claim: The 4th parameter is a naked bool 'mz_tolerance_unit' whose polarity (ppm=true, Da=false) IS documented in the param block on the declaration (header line 63), and the only production call site (FeatureFinderMultiplexAlgorithm.cpp:1023) makes it self-documenting via '(param_.getValue("algorithm:mz_unit") == "ppm")'. The genuine, but minor, surprise is the bool-as-mode-selector idiom: a future caller writing a literal true/false could silently invert grid spacing (relative ppm step vs absolute Da step) and rt_scaling. No call site with a bare literal exists. Recommendation (additive enum overload) remains sound and non-breaking; the existing bool ctor's direction is already documented. ABI impact of the existing code is none; the suggested enum overload would be source-compatible/non-breaking.

### [FEAT-48] MultiplexClustering::MultiplexClustering(const MSExperiment&, double, bool, double) — @throw documentation on the centroid-only constructor is copy-pasted and wrong; it actually throws on out-of-range MZ/RT ranges
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/MultiplexClustering.h · _ff-multiplex-a_

```cpp
MultiplexClustering(const MSExperiment& exp, double mz_tolerance, bool mz_tolerance_unit, double rt_typical)
```
- **Expectation:** The documented @throw condition matches what the constructor checks. This constructor takes no peak-boundaries argument, so a 'boundaries vs spectra count mismatch' throw is impossible here.
- **Actual:** The header documents '@throw Exception::IllegalArgument if centroided data and the corresponding list of peak boundaries do not contain same number of spectra' (copied from the other constructor). The implementation has no boundaries parameter and instead throws Exception::IllegalArgument when MinMZ/MaxMZ/MinRT/MaxRT fall outside sensible ranges (uninitialized ranges). A caller relying on the doc would neither expect the real throw nor guard against it.
- **Evidence:** Header lines 66 (doc) vs 68 (signature without boundaries). Impl (MultiplexClustering.cpp:88-92): 'if (!RangeMZ(0.0,1.0e12).containsMZ({mz_min,mz_max}) || !RangeRT(-1.0e12,1.0e12).containsRT({rt_min,rt_max})) { throw Exception::IllegalArgument(... "MinMZ,MaxMZ,MinRT,MaxRT values outside of sensible value ranges..."); }'.
- **Fix:** Fix the @throw doc to describe the real precondition (sensible/initialized m/z and RT ranges). Documentation-only; no ABI impact.
- **Verifier correction:** The @throw on the centroid-only constructor MultiplexClustering(const MSExperiment&, double, bool, double) (header line 66, for the declaration at line 68) is a copy-paste of the boundaries-constructor's doc (line 54) and is wrong: this overload has no peak-boundaries parameter, so a 'boundaries vs spectra count mismatch' throw cannot occur. The actual throw (MultiplexClustering.cpp:88-92) is Exception::IllegalArgument when the experiment's MS1 m/z/RT ranges fall outside RangeMZ(0,1e12)/RangeRT(-1e12,1e12) — i.e. when the MSExperiment is empty or has uninitialized ranges. The throw is loud and recoverable, so this is a low-severity documentation defect, not a silent-wrong-result hazard. Fix: replace the @throw text to describe the real precondition (non-empty / properly initialized MS1 m/z and RT ranges). Documentation-only, no ABI impact.

### [FEAT-49] MultiplexIsotopicPeakPattern::getMassShifts — getMassShifts() returns a full deep copy of the MultiplexDeltaMasses by value, not a const reference
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h · _ff-multiplex-a_

```cpp
MultiplexDeltaMasses getMassShifts() const
```
- **Expectation:** A const getter named getMassShifts on a member of class type returns a const reference to the member (cheap), as is the norm for OpenMS get-by-ref accessors.
- **Actual:** getMassShifts() returns 'MultiplexDeltaMasses' by value, deep-copying the underlying std::vector<DeltaMass> (each with a std::multiset<std::string>) on every call. A caller doing 'pattern.getMassShifts().getDeltaMasses()' in a loop silently pays repeated allocation, and 'auto& ms = pattern.getMassShifts();' binds to a temporary, not the member.
- **Evidence:** Header line 50: 'MultiplexDeltaMasses getMassShifts() const;'. Impl (MultiplexIsotopicPeakPattern.cpp:44-47): 'MultiplexDeltaMasses MultiplexIsotopicPeakPattern::getMassShifts() const { return mass_shifts_; }'.
- **Fix:** Add a 'const MultiplexDeltaMasses& getMassShifts() const' accessor (additive overload is impossible since differing only by return type, so add a new name like massShifts() or change the return type in a major-version bump). Document the copy in the interim.

### [FEAT-50] MultiplexFiltering::getBlacklist — getBlacklist() is a pure read-only computation but is declared non-const
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/FEATUREFINDER/MultiplexFiltering.h · _ff-multiplex-b_

```cpp
MSExperiment getBlacklist()
```
- **Expectation:** A getter that only reads internal state and returns a freshly built value should be const, callable on a const MultiplexFiltering.
- **Actual:** getBlacklist() builds a brand-new local MSExperiment from exp_centroided_ and blacklist_ and returns it by value, mutating nothing — yet it is non-const, so it cannot be called on a const instance and falsely signals that calling it changes the object.
- **Evidence:** Header line 91: `MSExperiment getBlacklist();` (no const). Impl MultiplexFiltering.cpp:378-402 constructs a local `exp_blacklist`, only reads members, and `return exp_blacklist;`.
- **Fix:** Add `const` to the declaration and definition: `MSExperiment getBlacklist() const;`. Pure addition of const is source-compatible for callers; it does change the mangled name (ABI) so note that.
- **Verifier correction:** getBlacklist() is genuinely a pure read-only computation declared non-const and should be `MSExperiment getBlacklist() const;`. However, severity is low, not high: the only caller uses a non-const instance so nothing breaks, and the omission causes no wrong results/crash/data loss — just a const-correctness papercut. Adding const is source-compatible for all callers; it changes the mangled symbol name (an ABI-level change requiring dependents to recompile), but no source edits are needed at call sites.

### [FEAT-53] MultiplexFilteredPeak::size — size() counts satellite peaks INCLUDING the primary peak, so it is not the count of 'satellite peaks' the name/doc imply
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/FEATUREFINDER/MultiplexFilteredPeak.h · _ff-multiplex-b_

```cpp
size_t size() const
```
- **Expectation:** Given getSatellites() and the doc 'return number of satellite peaks', size() should return the number of satellites distinct from the primary peak.
- **Actual:** The doc for satellites_ states 'The primary peak is part of the satellite peak set' (header lines 138-139), and size() returns satellites_.size() directly. So size() counts the primary peak as one of the satellites — a caller using size() as 'how many supporting satellite peaks were found' is off by the primary and per-pattern duplicates.
- **Evidence:** Header line 108 doc 'return number of satellite peaks'; lines 138-140 'Mapping from a pattern index ... The primary peak is part of the satellite peak set.' Impl MultiplexFilteredPeak.cpp:91-94 `return satellites_.size();`.
- **Fix:** Clarify in the header doc that size() includes the primary peak (and may include multiple satellites per pattern across the RT band). If a 'satellites excluding primary' count is meaningful to callers, add a separate explicitly named accessor. Doc-only fix is ABI-safe.
- **Verifier correction:** size() returns the total number of supporting centroided satellite peaks for this filtered peak — including the satellite that coincides with the primary peak, and including multiple satellites per pattern across the RT band. This is the documented design (header line 138: "The primary peak is part of the satellite peak set"; class doc lines 29-39 define satellites as the pattern-forming peaks and cross-reference @see size()). It is NOT promised to exclude the primary or to be a per-pattern count. The only genuine deficiency is that the size() accessor doc on line 108 ("return number of satellite peaks") does not cross-reference the line-138 note, so a caller reading only the public accessor doc could mistake it for a "satellites excluding primary" count. Severity is low (mild, locally documented, no silent wrong results, no real consumer affected); fix is a doc clarification only and is ABI-safe.

### [PROC-24] MZTrafoModel::predict — predict() reads coeff_[0..2] with no trained/empty check; on an untrained model it is out-of-bounds UB despite the class providing isTrained()/isValidModel()
`low` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/CALIBRATION/MZTrafoModel.h · _proc-calibration-noise_

```cpp
double predict(double mz) const
```
- **Expectation:** Given the class exposes isTrained() and isValidModel(), a caller might call predict() on a default-constructed model and expect either a throw, an identity transform, or a documented precondition error -- not silent memory corruption.
- **Actual:** predict() unconditionally evaluates `coeff_[0] + coeff_[1]*mz + coeff_[2]*mz*mz`. On a default-constructed (untrained) model coeff_ is empty, so coeff_[0] is out-of-bounds access (UB). The header only says 'Make sure the model was trained' in prose; there is no assert/throw, unlike getCoefficients() which DOES throw Precondition.
- **Evidence:** MZTrafoModel.cpp predict(): `double predict = coeff_[0] + coeff_[1]*mz + coeff_[2]*mz*mz;` with no isTrained() guard. Contrast getCoefficients(): `if (!isTrained()) throw Exception::Precondition(...)`.
- **Fix:** Add `if (!isTrained()) throw Exception::Precondition(...)` (or OPENMS_PRECONDITION) at the top of predict(), mirroring getCoefficients(). Source-compatible; ABI unchanged.
- **Verifier correction:** predict() does index coeff_[0..2] with no runtime guard, which is UB on an untrained (empty-coeff_) model, and getCoefficients() does throw Precondition while predict() does not — both verified. But the precondition is explicitly documented inline at the function (header line 143) and every in-tree production caller (InternalCalibration) gates predict() behind isValidModel(), so this is a documented-precondition inconsistency / standard unchecked-operator[] hot-path idiom rather than silent memory corruption for reasonable use. Downgrade from high to low. Adding OPENMS_PRECONDITION(isTrained(), ...) at the top of predict() to mirror getCoefficients() is a reasonable, source-compatible, ABI-neutral hardening (a release-build no-op via OPENMS_PRECONDITION, or a throw if Exception::Precondition is used).

### [PROC-25] MZTrafoModel::getCoefficients — getCoefficients() is a read-only accessor but is non-const, so it cannot be called on a const MZTrafoModel
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/PROCESSING/CALIBRATION/MZTrafoModel.h · _proc-calibration-noise_

```cpp
void getCoefficients(double& intercept, double& slope, double& power)
```
- **Expectation:** A 'get' accessor that only copies the three internal coefficients into out-params should be const, callable on a `const MZTrafoModel&` (the sibling read accessors getRT()/isTrained()/predict()/toString() are all const).
- **Actual:** getCoefficients() is declared non-const even though its body only reads coeff_ (after an isTrained() check) and writes to the caller's out-params. It cannot be invoked through a const reference, breaking const-correctness symmetry with the other getters.
- **Evidence:** Header: `void getCoefficients(double& intercept, double& slope, double& power);` (no trailing const). Body: only reads coeff_[0..2] into out refs.
- **Fix:** Add `const` to the declaration and definition. This is technically an ABI break (mangled name changes), but it is additive in spirit and source-compatible for all existing callers; if strict ABI must be preserved, add a const overload.

### [PROC-29] SignalToNoiseEstimatorMedian / SignalToNoiseEstimatorMeanIterative noise_for_empty_window — Sparse/empty windows silently receive a fixed noise of 1e20 (=> S/N ~ 0) instead of being flagged, so unreliable estimates are indistinguishable from real near-zero S/N
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h · _proc-calibration-noise_

```cpp
defaults_.setValue("noise_for_empty_window", std::pow(10.0, 20), ...)
```
- **Expectation:** When a window has too few points to estimate noise, a caller reading getSignalToNoise(i) would expect a way to know the estimate is unreliable, not a magic number silently substituted.
- **Actual:** For windows with fewer than min_required_elements points, noise = noise_for_empty_window_ (default 1e20), so stn_estimates_[i] = intensity/1e20 ~ 0. The caller gets a normal-looking double with no flag; only a stderr/log warning is emitted (and only if >20% sparse). A near-zero S/N from a sparse window is indistinguishable from a genuine low-S/N peak.
- **Evidence:** Median ctor: `defaults_.setValue("noise_for_empty_window", std::pow(10.0, 20), "noise value used for sparse windows", ...)`. computeSTN_: `if (elements_in_window < min_required_elements_) { noise = noise_for_empty_window_; ... } ... stn_estimates_[window_count] = intensity / noise;`
- **Fix:** Document this prominently in getSignalToNoise()'s contract, or expose a per-point 'was-sparse' flag. At minimum the base-class getSignalToNoise() @note should mention the sentinel behavior. Additive (new query method); ABI-safe.
- **Verifier correction:** The behavior is real but the claim overstates it as a near-undocumented magic number with corruption risk. In fact the sentinel is documented at the class level (brief, class @note line 44, and member comment lines 399-401), and the 1e20 value drives S/N toward 0, which conservatively SUPPRESSES sparse-window points (they fall below any threshold) rather than producing falsely-high/retained peaks. The genuine, narrower gap is that the base-class getSignalToNoise() contract does not surface a per-point reliability flag, so a sparse-window S/N~0 is indistinguishable per-point from a real low-S/N peak. This is a low-severity documentation/API-ergonomics gap, fully recoverable via getSparseWindowPercent() and the log warning; the fix is additive (per-point sparse flag or a getSignalToNoise() @note) and ABI-safe.

### [PROC-30] SignalToNoiseEstimatorMedian::computeSTN_ vs SignalToNoiseEstimatorMeanIterative::computeSTN_ — Sibling S/N estimators use different position accessors and window-boundary inequalities (getPos() with <= vs getMZ() with <), so the median and mean-iterative estimators define the sliding window inconsistently
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h · _proc-calibration-noise_

```cpp
void computeSTN_(const Container& c) override
```
- **Expectation:** Two estimators in the same module, both documented with 'win_len in Thomson' and a sliding m/z window, should define the window identically so results are comparable apart from the median-vs-mean statistic.
- **Actual:** SignalToNoiseEstimatorMedian uses `(*it).getPos()` and includes the right border with `<=`; SignalToNoiseEstimatorMeanIterative uses `(*it).getMZ()` and excludes the right border with `<`. The differing accessor and the half-open vs half-closed boundary mean the two estimators silently use slightly different window membership for the same data.
- **Evidence:** Median: `while ((*window_pos_borderright).getPos() <= (*window_pos_center).getPos() + window_half_size)`. MeanIterative: `while ((*window_pos_borderright).getMZ() < (*window_pos_center).getMZ() + window_half_size)`.
- **Fix:** Unify the position accessor and boundary convention across the two estimators (pick getMZ() or getPos() and one inequality). Behavioral; keep defaults so existing numeric expectations are reviewed. ABI unaffected (protected/templated).
- **Verifier correction:** getPos() and getMZ() are identical aliases (Peak1D.h: getPos is documented "Alias for getMZ()", both return position_[0]) — that part of the claim is false and has no effect. The only real difference is the right-window-boundary inequality: SignalToNoiseEstimatorMedian includes the right border (<=, closed window [c-h, c+h]) while SignalToNoiseEstimatorMeanIterative excludes it (<, half-open [c-h, c+h)). Both use inclusive left borders. The discrepancy only manifests for a data point whose m/z exactly equals center + win_len/2, affecting window membership by at most one peak — a negligible, non-silent edge case, not a meaningful "different window" for general data.

### [PROC-1] Deisotoper::deisotopeWithAveragineModel — Doxygen describes parameters 'rem_low_intensity'/'used_for_open_search' that do not exist; the real knob is the integer 'number_of_final_peaks'
`low` · `stale-documentation` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h · _proc-centroiding_

```cpp
static void deisotopeWithAveragineModel(MSSpectrum& spectrum, double fragment_tolerance, bool fragment_unit_ppm, int number_of_final_peaks = 5000, int min_charge = 1, int max_charge = 3, bool keep_only_deisotoped = false, unsigned int min_isopeaks = 2, unsigned int max_isopeaks = 10, bool make_single_charged = true, bool annotate_charge = false, bool annotate_iso_peak_count = false, bool add_up_intensity = false)
```
- **Expectation:** Reading the class-level and method-level documentation, a caller expects a boolean 'rem_low_intensity' flag plus a 'used_for_open_search' switch that selects between keeping the top 1000 (open search) or 5000 peaks.
- **Actual:** The signature has no such parameters. Low-intensity filtering is controlled solely by the integer 'number_of_final_peaks' (default 5000); the caller must pass 1000 explicitly for open search. The doc text 'If @p rem_low_intensity is true, peaks not belonging to the highest 1000/5000 peaks are removed (see @p used_for_open_search)' refers to a previous API that no longer exists.
- **Evidence:** Header doc lines 93-95: "intensity 0. If @p rem_low_intensity is true, peaks not belonging to the highest 1000/5000 peaks are removed (see @p used_for_open_search)." — but the signature parameter is 'int number_of_final_peaks = 5000' (line 119) and the .cpp uses 'if (number_of_final_peaks > 0) { NLargest nlargest_filter = NLargest(number_of_final_peaks); ... }' (Deisotoper.cpp:63-66). grep confirms 'rem_low_intensity'/'used_for_open_search' appear nowhere except the comment.
- **Fix:** Fix the @param/@p references in the doc to describe 'number_of_final_peaks' (0 = no filtering; 1000 recommended for open search, else 5000). Doc-only change, no signature change. ABI: none.
- **Verifier correction:** The header's prose (Deisotoper.h:93-95) references two parameters that do not exist anywhere in the code (rem_low_intensity, used_for_open_search) — leftover from a prior API. The actual low-intensity filtering is controlled by the integer number_of_final_peaks (default 5000; 0 disables; pass 1000 for open search), implemented in Deisotoper.cpp:63-67 via NLargest. This is a doc-only inconsistency, classified as stale-documentation (not misleading-name) and low severity: the correct @param for number_of_final_peaks already exists on line 105, the phantom names are absent from the signature, and there is no behavioral hazard. Fix: rewrite lines 93-95 to reference number_of_final_peaks. ABI: none.

### [PROC-2] PeakPickerIM::pickIMTraces — pickIMTraces is non-const while sibling pickIMCluster/pickIMElutionProfiles are const, despite all three being in-place transforms that do not mutate the picker's logical state
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h · _proc-centroiding_

```cpp
void pickIMTraces(MSSpectrum& spectrum)
```
- **Expectation:** Three sibling methods with the same '@param[in,out] spectrum' in-place contract should have the same const-qualification; a caller holding a 'const PeakPickerIM&' expects to be able to call all three.
- **Actual:** pickIMTraces is non-const; pickIMCluster and pickIMElutionProfiles are const. None of the three modify the picker's parameters; the only mutable state used (ccs_warning_shown_) is already declared 'mutable'. pickIMTraces is non-const only because it calls private helpers (sumFrame_, extractIonMobilityTraces, computeOptimalSamplingRate) that were left non-const, not for any semantic reason.
- **Evidence:** Header: 'void pickIMTraces(MSSpectrum& spectrum);' (line 56, no const) vs 'void pickIMCluster(MSSpectrum& spec) const;' (line 92) and 'void pickIMElutionProfiles(MSSpectrum& input) const;' (line 102). The private helpers it calls (e.g. 'void sumFrame_(...)' line 121, 'std::vector<MSSpectrum> extractIonMobilityTraces(...)' line 129) are non-const but only read picker members. 'mutable bool ccs_warning_shown_' (line 148) already handles the only write.
- **Fix:** Make pickIMTraces (and its private helpers) const for symmetry. This is source-compatible for callers (a non-const object can still call a const method) but changes the member function signature. ABI: breaking (mangled name changes). If strict ABI must be kept, add a const overload.
- **Verifier correction:** The inconsistency is real and correctly described, but its severity is low, not high. A caller holding a const PeakPickerIM& can call pickIMCluster/pickIMElutionProfiles but not pickIMTraces; this is enforced at compile time (loud, no silent miscompute, no data loss, no crash) and any non-const PeakPickerIM can call all three. Fixing it by making pickIMTraces (and its private helpers) const is source-compatible for call sites but ABI-breaking (mangled-name change); if strict ABI must hold, the const-overload alternative the claim mentions applies.

### [PROC-4] PeakPickerIterative (defaults) — User-facing parameter names carry a trailing underscore ('signal_to_noise_', 'sn_bin_count_', 'nr_iterations_', 'sn_win_len_'), the C++ private-member convention, breaking the established param naming
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/CENTROIDING/PeakPickerIterative.h · _proc-centroiding_

```cpp
defaults_.setValue("signal_to_noise_", ...); setValue("sn_bin_count_", ...); setValue("nr_iterations_", ...); setValue("sn_win_len_", ...)
```
- **Expectation:** TOPP/Param keys exposed to users follow lowercase_words without a trailing underscore (the trailing underscore is reserved for C++ class members). Sibling parameters in the very same constructor ('peak_width', 'spacing_difference', 'check_width_internally', 'ms1_only') follow this.
- **Actual:** Four parameters are registered with a trailing underscore, so users must write 'signal_to_noise_' / 'sn_bin_count_' / 'nr_iterations_' / 'sn_win_len_' in INI files and on the command line — inconsistent with both sibling keys here and with PeakPickerHiRes ('signal_to_noise').
- **Evidence:** PeakPickerIterative.h:92,97-100: 'defaults_.setValue("signal_to_noise_", 1.0, ...)', 'defaults_.setValue("sn_bin_count_", 30, ...)', 'defaults_.setValue("nr_iterations_", 5, ...)', 'defaults_.setValue("sn_win_len_", 20.0, ...)' vs 'defaults_.setValue("peak_width", ...)' (line 93) and 'spacing_difference' (line 96). These keys are surfaced to TOPP via getDefaults() in src/topp/PeakPickerIterative.cpp:89.
- **Fix:** Rename the keys to drop the trailing underscore. This is a user-visible/INI-breaking change, so ideally register the new names while accepting the old as deprecated aliases. ABI: none (string keys), but behavior/INI breaking — handle with alias + deprecation.
- **Verifier correction:** Severity lowered from (implied) higher to low: the four keys signal_to_noise_/sn_bin_count_/nr_iterations_/sn_win_len_ do carry an inconsistent trailing underscore (vs sibling peak_width/spacing_difference/check_width_internally/ms1_only and vs PeakPickerHiRes 'signal_to_noise'), but the picker works correctly with the actual registered names, the inconsistency is discoverable via getDefaults/CTD and surfaces as a loud unknown-parameter error rather than silent wrong results. Recommendation to rename with deprecated aliases is reasonable; ABI impact is none (string keys), though a rename would be INI-behavior-breaking.

### [PROC-7] PeakPickerHiRes::pick(const MSSpectrum&, MSSpectrum&, std::vector<PeakBoundary>&, bool) — The 'check_spacings' default differs by overload (true for MSSpectrum, false for MSChromatogram/Mobilogram) — same-named bool flag with opposite defaults depending on container type
`low` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h · _proc-centroiding_

```cpp
void pick(const MSSpectrum& input, MSSpectrum& output, std::vector<PeakBoundary>& boundaries, bool check_spacings = true) const
```
- **Expectation:** A bool flag with the same name across sibling overloads is expected to default consistently; a caller templating over container type or copying a call site between spectrum and chromatogram would assume the same default behavior.
- **Actual:** The MSSpectrum overload defaults check_spacings=true while the MSChromatogram and Mobilogram overloads default it to false. Copying a spectrum call to a chromatogram (or vice versa) and relying on the default silently changes whether spacing constraints are enforced.
- **Evidence:** PeakPickerHiRes.h line 108: 'void pick(const MSSpectrum&..., bool check_spacings = true)' vs line 119: 'void pick(const MSChromatogram&..., bool check_spacings = false)' and line 130: 'void pick(const Mobilogram&..., bool check_spacings = false)'.
- **Fix:** This is intentional and is documented inline ('yes for spectra, no for chromatograms'); keep but ensure the asymmetry is highlighted at the parameter, not only buried in the @param text, since it is an easy silent-bug source when call sites are copied. ABI: none (doc clarification only).
- **Verifier correction:** The asymmetry is genuine and behaviorally significant, but the claim overstates it as a silent/buried trap: the differing default is explicitly documented at each overload's @param ("yes for spectra, no for chromatograms"), is driven by a real domain distinction (regular m/z spacing vs. irregular chromatogram/mobilogram axes), and is partially self-correcting in pick_() which force-disables spacing checks when thresholds are infinite. It is a low-severity (mild) surprise, not medium — misuse via copy-paste is plausible but visible in review and recoverable.

### [PROC-8] Deisotoper::deisotopeAndSingleCharge — The two deisotope overloads default 'min_isopeaks' inconsistently (3 here vs 2 in deisotopeWithAveragineModel), so switching algorithms silently changes cluster acceptance
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h · _proc-centroiding_

```cpp
static void deisotopeAndSingleCharge(MSSpectrum& spectrum, double fragment_tolerance, bool fragment_unit_ppm, int min_charge = 1, int max_charge = 3, bool keep_only_deisotoped = false, unsigned int min_isopeaks = 3, unsigned int max_isopeaks = 10, bool make_single_charged = true, ...)
```
- **Expectation:** Two sibling deisotoping algorithms documented as sharing 'the same in-place spectrum contract' would be expected to share the same default for the minimum-cluster-size parameter, so swapping one for the other does not silently alter which peaks are treated as isotope clusters.
- **Actual:** deisotopeAndSingleCharge defaults min_isopeaks=3 while deisotopeWithAveragineModel defaults min_isopeaks=2. A caller switching algorithms while keeping defaults gets a different minimum isotope-peak requirement, changing results.
- **Evidence:** Header line 123: 'unsigned int min_isopeaks = 2' (deisotopeWithAveragineModel) vs line 177: 'unsigned int min_isopeaks = 3' (deisotopeAndSingleCharge). Class doc (lines 26-29) states they 'share the same in-place spectrum contract'.
- **Fix:** Document the differing default explicitly (or align them). Doc-only or default-alignment; aligning the default is observably behavior-breaking, so prefer documenting. ABI: none.
- **Verifier correction:** The two sibling deisotoping methods default min_isopeaks inconsistently (2 for deisotopeWithAveragineModel vs 3 for deisotopeAndSingleCharge) despite a class doc claiming a shared contract and identical per-parameter documentation. This is a real documentation/convention inconsistency, but lower impact than implied: the methods have differing parameter lists (so they are not drop-in positional swaps), nearly all real callers pass min_isopeaks explicitly, each method is internally consistent, and the threshold is freely tunable and recoverable. Fix is doc-only (note the differing default in both doc blocks) or default alignment (behavior-breaking, so documenting is preferred). No ABI impact.

### [PROC-9] PeakCandidate::mz — PeakCandidate stores the m/z centroid in a 'float mz' field (single precision) while the weighted m/z is computed in double, silently truncating high-resolution m/z accuracy
`low` · `unit-or-index` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/CENTROIDING/PeakPickerIterative.h · _proc-centroiding_

```cpp
struct PeakCandidate { int index; double peak_apex_intensity; double integrated_intensity; double leftWidth; double rightWidth; float mz; };
```
- **Expectation:** In a high-resolution TOF peak picker (header: 'designed for TOF-MS data', 'more accurate determination of peak center'), the peak m/z carried through the algorithm should be double precision, matching how it is computed and how Peak1D stores positions.
- **Actual:** The recomputed 'weighted_mz' is a double but is stored into 'float mz', losing precision; this same float is then compared against double leftWidth/rightWidth borders for the enclosure test, so the de-duplication threshold operates on a truncated value.
- **Evidence:** PeakPickerIterative.h struct: 'float mz;' (line 33). Assignment: 'double weighted_mz = 0; ... weighted_mz /= integrated_intensity; ... PeakCandidates[peak_it].mz = weighted_mz;' (lines 221-234). Enclosure test mixes the float against doubles: 'if (PeakCandidates[m].mz >= PeakCandidates[peak_it].leftWidth && PeakCandidates[m].mz <= PeakCandidates[peak_it].rightWidth)' (line 345).
- **Fix:** Change 'float mz' to 'double mz' to match the precision of the computation and the borders. This is an internal struct (declared in the public header but used only by the algorithm); ABI: source-compatible for users who do not use PeakCandidate directly, but the struct layout changes.
- **Verifier correction:** The float-mz narrowing is real and propagates to the output centroid (Peak1D, line 355), not just internal de-duplication. But the precision lost is ~0.01-0.03 ppm at typical m/z (400-1800), well below instrument mass accuracy and below the border spacing in the enclosure test, so it does not cause practically wrong results. Severity is low (mild, avoidable inconsistency), not high. ABI: PeakCandidate is declared in a public header but is used only internally by PeakPickerIterative and is not part of any public API signature, so changing float->double is source-compatible for users (only an internal layout change).

### [PROC-35] SplineInterpolatedPeaks::getNavigator — getNavigator() documentation tells callers to check getSplineCount(), a method that does not exist
`low` · `other` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/MISC/SplineInterpolatedPeaks.h · _proc-misc_

```cpp
SplineInterpolatedPeaks::Navigator getNavigator(double scaling = 0.7)
```
- **Expectation:** Documentation that instructs the caller to guard a throwing call ('Will throw ... Check using getSplineCount()') should reference a real, callable method.
- **Actual:** There is no getSplineCount() anywhere in the codebase; the count accessor is named size(). A caller following the docs to avoid the Exception::InvalidSize thrown when packages_ is empty would fail to compile or hunt for a nonexistent method, and may not realize size() is the guard to use.
- **Evidence:** SplineInterpolatedPeaks.h:154 `Check using getSplineCount().`; the actual accessor is `size_t size() const;` (:78); grep over src/ finds getSplineCount only in this comment.
- **Fix:** Update the doc to reference size(). Comment-only change; source- and ABI-compatible.
- **Verifier correction:** The doc comment at SplineInterpolatedPeaks.h:154 references a nonexistent method getSplineCount(); it should reference size() (which returns packages_.size()). A caller following the doc literally gets a compile error (loud/recoverable), not silently wrong behavior, so this is a low-severity documentation defect. Comment-only fix; source- and ABI-compatible (abi_impact: none).

### [PROC-36] SplineInterpolatedPeaks::SplineInterpolatedPeaks(const MSSpectrum&) / (const MSChromatogram&) / (const std::vector<double>&, const std::vector<double>&) — Single-argument constructors are non-explicit, enabling implicit MSSpectrum/MSChromatogram -> SplineInterpolatedPeaks conversions
`low` · `implicit-conversion` · ABI: `source-compatible` · src/openms/include/OpenMS/PROCESSING/MISC/SplineInterpolatedPeaks.h · _proc-misc_

```cpp
SplineInterpolatedPeaks(const MSSpectrum& raw_spectrum)
```
- **Expectation:** A heavyweight value type that builds spline fits over an entire spectrum should not be implicitly constructible from an MSSpectrum/MSChromatogram; conversions of this magnitude should be explicit to avoid accidental spline-fitting in overload resolution or function-argument passing.
- **Actual:** Both the MSSpectrum and MSChromatogram constructors (and the two-vector constructor) are declared without `explicit`, so an MSSpectrum can implicitly convert to a SplineInterpolatedPeaks wherever the latter is expected, silently running a full spline initialization.
- **Evidence:** SplineInterpolatedPeaks.h:42,48,54 declare the constructors with no `explicit` keyword.
- **Fix:** Mark the single-argument constructors `explicit`. This is source-compatible for normal direct-initialization usage but could break code relying on implicit conversion; do it in a minor release. No vtable/layout change.
- **Verifier correction:** Only the two single-argument constructors are implicit converting constructors at risk: SplineInterpolatedPeaks(const MSSpectrum&) (line 48) and SplineInterpolatedPeaks(const MSChromatogram&) (line 54). The two-vector constructor (line 42) takes two arguments and cannot participate in implicit conversion, so it should be dropped from the conversion claim. Marking the two single-argument constructors `explicit` is correct and source-compatible (no in-tree usage relies on implicit conversion; all current uses are direct-initialization). Severity is low, not high: no function/operator in the codebase takes SplineInterpolatedPeaks by value, so no silent spline-fitting occurs in overload resolution today, and any accidental conversion still yields a correct result.

### [PROC-12] FastLowessSmoothing::lowess — int return looks like a success/error status but is not one; degenerate n<2 input returns 1 and silently leaves result mostly zero
`low` · `return-value` · ABI: `none` · src/openms/source/PROCESSING/SMOOTHING/FastLowessSmoothing.cpp · _proc-smoothing_

```cpp
int OPENMS_DLLAPI lowess(const std::vector<double>& x, const std::vector<double>& y, double f, int nsteps, double delta, std::vector<double>& result)
```
- **Expectation:** A function returning int and documented only as writing its smoothing into an out-param 'result' invites callers to treat the int as an error/status code (0=ok, nonzero=error). One would expect a nonzero return to flag a problem.
- **Actual:** The return is not an error code. For valid input (n>=2) it always returns 0. For n<2 it returns 1 but that is the DEGENERATE success case where only result[0]=y[0] is written and the rest of result is left zero-initialized. So '1' (which a caller would read as an error) actually corresponds to a silently-truncated/partial result, and '0' is the normal case. The header never documents the return semantics, and the sole in-tree caller (TransformationModelLowess.cpp:165) ignores it.
- **Evidence:** Header FastLowessSmoothing.h:73 declares 'int ... lowess(...)' with no @return doc; cpp:492-496 'if (n < 2) { ys[0] = y[0]; return 1; }' and cpp:550 'return 0;' as the only two return values. result is resized(n) and only result[0] gets set in the n<2 branch.
- **Fix:** Document the return value explicitly (or change to void, since callers do not use it) and invert the convention so failure is signaled clearly, or throw on n<2 instead of partially filling result. Changing the signature to void is ABI-breaking; the ABI-safe fix is a @return doc clarification plus making the n<2 case fill result consistently (e.g. result.assign(n, y[0]) or throw).
- **Verifier correction:** The int return of FastLowessSmoothing::lowess is undocumented and resembles a 0=ok/nonzero=error status code, but actually returns 0 for all valid input and 1 only for the degenerate n<2 corner — an inverted, confusing convention. However: the value is never consulted by any in-tree code (caller and tests both ignore it), the n<2 branch is unreachable from real callers (TransformationModelLowess throws on size<2 first; OPENMS_PRECONDITION throws in assertion builds), and the claim's failure description is wrong — for n==1 result is fully written, and for n==0 the implementation hits out-of-bounds UB rather than returning 1 with a zero-filled result. This is a low-severity documentation nit, not a silent-wrong-result/data-loss bug. ABI-safe fix: add a @return doc note (and optionally throw or assign result consistently for n<2); changing the signature to void would be ABI-breaking and unnecessary.

### [PROC-14] GaussFilterAlgorithm::filter — bool return means 'any nonzero output was produced', not 'filtering succeeded' as the name suggests
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/PROCESSING/SMOOTHING/GaussFilterAlgorithm.h · _proc-smoothing_

```cpp
bool filter(ConstIterT mz_in_start, ConstIterT mz_in_end, ConstIterT int_in_start, IterT mz_out, IterT int_out)
```
- **Expectation:** A function named filter() returning bool is naturally read as 'true = filtering succeeded / applied'. A caller would write 'if (!filter(...)) handleError();'.
- **Actual:** The returned bool is 'found_signal', set true only if at least one output intensity has |value|>0. A perfectly valid all-zero (or all-flattened) region returns false even though filtering ran and wrote outputs normally. So false does not mean failure; it means 'the result was all zeros'. The GaussFilter wrapper even uses this to decide whether to copy results back, conflating 'no signal' with 'do not update data'.
- **Evidence:** GaussFilterAlgorithm.h:117 'bool found_signal = false;' ... h:137 'if (fabs(new_int) > 0) found_signal = true;' ... h:139 'return found_signal;'. The output iterators are written for every input point regardless (h:132-135).
- **Fix:** Rename the return concept in docs to 'returns true if any non-zero signal remained after filtering' (the brief currently says nothing about the return). ABI-safe doc fix; do not change the bool meaning since the GaussFilter wrapper depends on it.
- **Verifier correction:** The return bool of GaussFilterAlgorithm::filter is "found_signal": true iff at least one output intensity had |value| > 0 after convolution. It is NOT a success flag — outputs (m/z and intensity) are written to the output iterators for every input point regardless of the return value, so false never indicates that filtering failed to run or that outputs are missing. The brief @brief comment documents nothing about the return; it should state "returns true if any non-zero signal remained after filtering." The GaussFilter wrapper (GaussFilter.cpp) deliberately relies on this: when found_signal is false and the spectrum/chromatogram/mobilogram has >= 3 points, it leaves the original data untouched (does not copy the filtered all-zero result back) and optionally logs a warning; for sizes < 3 it copies regardless. Fix is documentation-only; do not change the bool's meaning since the wrapper depends on it.
