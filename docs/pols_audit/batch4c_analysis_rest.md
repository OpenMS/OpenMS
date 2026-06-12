# OpenMS POLS Audit — Batch 4c: ANALYSIS MAPMATCHING / QUANTITATION / NUXL / TOPDOWN

**Confirmed findings:** 69 · 8 high · 29 medium · 32 low. (solo re-run; clean)

> Post-verification severity/category/ABI. Recommendations favor ABI-safe fixes.

### [ANAL-9] MapAlignmentEvaluationAlgorithmRecall::evaluate — evaluate() divides by ground-truth size (and by cons_map_tool.size()) with no guard, returning NaN/inf instead of signaling on empty input
`high` · `silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.cpp · _mm-alignment_

```cpp
void evaluate(const ConsensusMap& consensus_map_in, const ConsensusMap& consensus_map_gt, const double& rt_dev, const double& mz_dev, const Peak2D::IntensityType& int_dev, const bool use_charge, double& out) override
```
- **Expectation:** A caller expects a recall value in [0,1], or a clear error when the ground truth or the tool's consensus map is empty.
- **Actual:** Same pattern as Precision: `out = (1.0 / double(cons_map_gt.size())) * sum;` divides by zero when no GT consensus feature has size >= 2. Additionally, line 85 `gt.push_back(gt_i / cons_map_tool.size());` divides by cons_map_tool.size(), which is zero (UB / crash or trap) when the input tool map is empty. Neither case is checked or signaled.
- **Evidence:** Line 85: `gt.push_back(gt_i / cons_map_tool.size());`; line 100: `double recall = (1.0 / double(cons_map_gt.size())) * sum;` — both denominators can be zero with no guard.
- **Fix:** Guard against empty cons_map_gt and empty cons_map_tool before dividing (throw or return a documented sentinel). Additive early checks, ABI-safe.

### [ANAL-1] FeatureGroupingAlgorithmUnlabeled::getResultMap — Incremental path (setReference/addToGroup/getResultMap) yields a different, un-postprocessed result than group()
`high` · `asymmetric-api` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h · _mm-grouping_

```cpp
ConsensusMap & getResultMap()
```
- **Expectation:** The header documents two equivalent ways to run: (a) group(), (b) setReference + addToGroup* + getResultMap. A caller expects getResultMap() to return the same finished consensus map that group() produces.
- **Actual:** group() ends with postprocess_(maps, out) which transfers protein IDs, re-indexes unassigned peptide IDs, and applies canonical sorting (sortByQuality/sortByMaps/sortBySize), and also fixes up column headers. The incremental path (addToGroup in the .cpp) never calls postprocess_ and never populates ColumnHeaders, so getResultMap() returns raw, unsorted pairfinder output with no protein IDs and empty/partial column headers. The TOPP tool FeatureLinkerUnlabeled.cpp has to manually patch this (out_map.setColumnHeaders(dummy.getColumnHeaders()); transferSubelements(...)). A caller following the header's 'two ways' contract gets a silently inferior result.
- **Evidence:** Header doc: 'There are two ways to run the algorithm: a) Simply call "group"... b) Call "setReference", "addToGroup" (n times), "getResultMap"'. group() in .cpp ends 'postprocess_(maps, out);'; addToGroup() ends 'pairfinder_input_[0].swap(result);' with no postprocess. FeatureLinkerUnlabeled.cpp line 262 'out_map.setColumnHeaders(dummy.getColumnHeaders());' shows the caller must fix headers itself.
- **Fix:** Document that getResultMap() returns the raw incremental result and that the caller must perform postprocessing/header setup themselves (as FeatureLinkerUnlabeled does), or add a finalizeResult() helper. Do not silently change getResultMap to run postprocess_ (would change existing TOPP output ordering). Doc fix is source-compatible.
- **Verifier correction:** The incremental path (setReference/addToGroup*/getResultMap) performs the SAME StablePairFinder linking as group(), so the consensus features themselves are equivalent. However, getResultMap() returns the raw pairfinder result WITHOUT the postprocessing that group() applies via postprocess_(): it lacks transferred protein IDs, lacks re-indexed unassigned peptide IDs (map_index), is not canonically sorted (sortByQuality/sortByMaps/sortBySize), and has column headers populated only for the reference map (setReference) — addToGroup never builds headers for the other maps. The header's "two ways to run" wording and getResultMap()'s doc imply a finished, equivalent map but deliver a metadata-incomplete, unsorted one. Fix: document that getResultMap() returns the raw incremental result requiring caller-side postprocessing (protein/peptide ID transfer, column-header setup, sorting) exactly as FeatureLinkerUnlabeled.cpp does, or add an additive finalizeResult() helper. Do not silently make getResultMap()/addToGroup call postprocess_, as that would change existing TOPP output ordering.

### [ANAL-42] NuXLFDR::splitIntoPeptidesAndXLs — splitIntoPeptidesAndXLs keeps only the single first hit overall, not the first hit of each class as documented
`high` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFDR.h · _nuxl-core_

```cpp
void splitIntoPeptidesAndXLs(const PeptideIdentificationList& peptide_ids, PeptideIdentificationList& pep_pi, PeptideIdentificationList& xl_pi) const
```
- **Expectation:** The header states an input identification 'contributes to pep_pi or xl_pi (or both, if it had hits of both classes)' and keeps 'the first hit encountered for each class'. A caller therefore expects a PSM with a plain-peptide hit and a cross-link hit to appear in BOTH output lists.
- **Actual:** The loop guards every push with `if (pep_ph.empty() && xl_ph.empty())` for BOTH branches (NuXLFDR.cpp lines 47-54). Once ANY first hit is added (to either class), this condition is false for all remaining hits, so a second hit of the OTHER class is never captured. An identification can therefore never contribute to both lists; only its single best hit is retained. This silently drops cross-link evidence whenever a plain-peptide hit ranks above it (and vice versa).
- **Evidence:** NuXLFDR.cpp: `if (static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0) { if (pep_ph.empty() && xl_ph.empty()) pep_ph.push_back(ph); } else { if (pep_ph.empty() && xl_ph.empty()) xl_ph.push_back(ph); }` -- the same combined emptiness check gates both branches.
- **Fix:** If both-class capture is intended (as the header claims), gate each branch on its own vector being empty (`if (pep_ph.empty())` / `if (xl_ph.empty())`). If single-best-hit is intended, fix the header doc to say so. Behavior change is the FDR-correct fix; doc fix is the conservative one. ABI: either is source-compatible.
- **Verifier correction:** The claim is accurate. Precise statement: with report_top_hits >= 2 (use_all_hits), a PeptideIdentification carrying hits of both classes is split into only one output list (whichever class's hit appears first in the hit order), and only that single hit is retained — silently discarding the documented per-class first hit of the other class before FDR is computed. Fix per recommendation: gate each branch on its own vector (`if (pep_ph.empty())` / `if (xl_ph.empty())`) to honor the header, or correct the header to state single-best-overall semantics. The behavioral fix is the FDR-correct one. ABI: none (function-body or comment-only change; signature and layout unchanged).

### [ANAL-43] NuXLFDR::calculatePeptideAndXLQValueAndFilterAtPSMLevel — Tie-breaker divides by score range and silently produces inf/NaN scores when all XL hits share one score
`high` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFDR.h · _nuxl-core_

```cpp
void calculatePeptideAndXLQValueAndFilterAtPSMLevel(..., std::vector<double> xl_PSM_qvalue_thresholds, std::vector<double> xl_peptidelevel_qvalue_thresholds, ...)
```
- **Expectation:** A deterministic 'add a tiny tie-breaker to q-value scores' step should leave scores well-defined for any input, including the degenerate case where all XL hits have equal score.
- **Actual:** The tie-breaker computes `score_range = max_score - min_score` and then `(1.0 - (score - min_score) / score_range) * 1e-5` (NuXLFDR.cpp lines 176-190). When every XL hit shares the same svm_score/NuXL:score, score_range == 0 and the division yields inf/NaN, corrupting every XL hit's stored score with no error. The header even flags this as a @warning but the code neither guards nor throws.
- **Evidence:** NuXLFDR.cpp: `double score_range = max_score - min_score;` then `double small_value = (1.0 - ((double)p.getMetaValue("svm_score") - min_score) / score_range) * 1e-5; p.setScore(p.getScore() + small_value);`
- **Fix:** Guard `if (score_range <= 0) skip tie-break` (or add an epsilon). Additive and non-breaking. ABI: source-compatible.
- **Verifier correction:** In the degenerate case (all XL hits share one score) score_range == 0 and the numerator (score - min_score) is also 0, so the division is 0/0 = NaN — the corrupted scores are NaN, not 'inf'. NaN is then added to every XL hit's score via p.setScore(...), with no guard or throw. The header @warning at NuXLFDR.h:131-133 documents the hazard but the code neither guards (e.g. `if (score_range <= 0) skip`) nor throws, and the NaN scores flow into IDFilter::filterHitsByScore and are written to the output _XLs.idXML. The recommended fix is a body-only guard inside the .cpp; it does not alter the function signature, so its ABI impact is 'none' rather than 'source-compatible'.

### [ANAL-51] NuXLPresets::getPresets — getPresets mutates the global ResidueDB singleton (adds methionine loss) for DEB/NM presets
`high` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLPresets.h · _nuxl-frag_

```cpp
void getPresets(const std::string& p, const std::string& custom_presets_file, StringList& nucleotides, StringList& mapping, StringList& modifications, StringList& fragment_adducts, std::string& can_cross_link);
```
- **Expectation:** A function named getPresets that fills output StringList/string parameters reads configuration; it should not mutate global chemistry state shared by the whole process.
- **Actual:** When the preset name contains the substring 'DEB' or 'NM', getPresets reaches into the global ResidueDB singleton, casts away const on the Methionine residue, and permanently adds a CH4S1 loss formula to it. This silently changes the chemistry of residue 'M' for every other consumer of ResidueDB in the process, and is not reversible.
- **Evidence:** 'if (StringUtils::hasSubstring(p, "DEB") || StringUtils::hasSubstring(p, "NM")) { auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M')); r_ptr->addLossFormula(EmpiricalFormula("CH4S1")); }' (NuXLPresets.cpp:160-165). The header documents only the six output params, no side effect.
- **Fix:** Document the ResidueDB mutation prominently in the header doc (@note that calling this with a DEB/NM preset alters global Methionine chemistry). Ideally move the residue-loss registration out of the preset getter into an explicit, separately-named setup call (e.g. applyPresetChemistry(...)) so a pure 'get' has no global effect — additive and source-compatible if the old behavior is preserved behind a deprecated path.
- **Verifier correction:** Claim is accurate as stated. Strengthening note: the mutation not only changes 'M' chemistry globally but is also cumulative — because addLossFormula does an unconditional push_back, invoking getPresets more than once with a DEB/NM preset in the same process appends the CH4S1 loss multiple times, compounding the corruption. It also directly contradicts the immutability invariant explicitly documented in ResidueDB.h (lines 31-37: "No public method changes the observable state of the singleton"), making the const_cast a deliberate circumvention rather than an oversight in an undocumented area.

### [ANAL-33] AbsoluteQuantitation::calculateRatio — calculateRatio returns 0.0 instead of signaling a missing feature/value
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h · _quant-absolute-mdv_

```cpp
double calculateRatio(const Feature & component_1, const Feature & component_2, const std::string & feature_name)
```
- **Expectation:** The doc advertises '@exception Exception::UnableToFit' and the name implies a real ratio. A caller expects either a valid ratio or a thrown error when the requested feature_name/IS is absent.
- **Actual:** When neither component has the requested metavalue/native_id, the function silently returns the initial 'ratio = 0.0' with only an OPENMS_LOG_DEBUG message; it never throws. A 0.0 ratio is indistinguishable from a legitimately measured zero and propagates into fitCalibration/applyCalibration as a real data point.
- **Evidence:** double ratio = 0.0; ... else { OPENMS_LOG_DEBUG << "Feature metaValue " << feature_name << " not found ..."; } return ratio; — no throw despite the documented '@exception Exception::UnableToFit'.
- **Fix:** Either throw Exception::UnableToFit (as documented) when the feature/IS value is missing, or return std::optional / a NaN sentinel that callers must check. Behavior change; gate behind a flag if existing pipelines rely on the 0.0 fallback. Doc-only fix (drop the bogus @exception) is the minimal source-compatible step.
- **Verifier correction:** Claim is accurate. Minor refinement: the silent-0.0 return is not merely incidental — it is explicitly enshrined by an existing unit test (AbsoluteQuantitation_test.cpp:201 asserts the missing-feature case == 0.0 via TEST_REAL_SIMILAR). This both strengthens the finding (the silent-0.0 path is exercised and intended) and means the 'throw Exception::UnableToFit' fix is more invasive than the recommendation implies: it would require updating that test. The minimal source-compatible step (dropping the bogus @exception, or returning a NaN sentinel/std::optional) remains valid; a throwing variant is a behavior change (ABI of the symbol unchanged, but runtime behavior and the test change).

### [ANAL-29] ItraqEightPlexQuantitationMethod::getReferenceChannel / updateMembers_ — Setting reference_channel=120 (an allowed param value) silently leaves the reference channel at a stale index
`high` · `silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.cpp · _quant-tmt-itraq_

```cpp
Size getReferenceChannel() const override
```
- **Expectation:** The 'reference_channel' param is range-validated to [113,121] (setMinInt(113)/setMaxInt(121)), and the description says '113-121'. A caller setting reference_channel=120 expects either an error or a defined reference channel index back from getReferenceChannel().
- **Actual:** 120 passes range validation and reaches updateMembers_, where the branch 'else if (ref_ch == 120) { OPENMS_LOG_WARN << "Invalid channel selection." }' only logs a warning and does NOT assign reference_channel_. The member keeps whatever value it had before (0 after construction, or a previously-set index), so getReferenceChannel() returns a silently wrong/stale index instead of signaling failure.
- **Evidence:** updateMembers_(): 'Int ref_ch = param_.getValue("reference_channel"); if (ref_ch == 121) { reference_channel_ = 7; } else if (ref_ch == 120) { OPENMS_LOG_WARN << "Invalid channel selection." << std::endl; } else { reference_channel_ = ref_ch - 113; }'. setDefaultParams_(): 'defaults_.setMinInt("reference_channel", 113); defaults_.setMaxInt("reference_channel", 121);' (120 is within range and the description even says 'Please note that 120 is not valid.').
- **Fix:** Make the invalid 120 selection a hard error: throw Exception::InvalidParameter in updateMembers_() for ref_ch==120 instead of only logging, or tighten validation so 120 cannot be set (e.g. setValidStrings/an explicit allowed-values check). Additive and ABI-safe (no signature change).

### [ANAL-58] DeconvolvedSpectrum::toSpectrum — toSpectrum() dereferences peak_groups_[0] with no empty check, crashing on an empty DeconvolvedSpectrum
`high` · `surprising-crash` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h · _topdown-deconv_

```cpp
MSSpectrum toSpectrum(int to_charge, double tol = 10.0, bool retain_undeconvolved = false)
```
- **Expectation:** Converting a (possibly empty) DeconvolvedSpectrum to an MSSpectrum should return an empty/header-only spectrum, not crash. Default-constructed DeconvolvedSpectrum objects exist and are common.
- **Actual:** Early in toSpectrum the code unconditionally reads peak_groups_[0].isPositive() to build the metadata string. If the spectrum has zero peak groups, this is out-of-bounds access (UB / crash) before any peaks are even iterated.
- **Evidence:** DeconvolvedSpectrum.cpp:41 `val << ... << std::to_string(FLASHHelperClasses::getChargeMass(peak_groups_[0].isPositive()));` with no preceding `if (empty()) return ...` guard.
- **Fix:** Add an early `if (peak_groups_.empty()) return out_spec;` (or branch the chargemass metadata) before indexing peak_groups_[0]. Pure implementation fix, no ABI/signature change.
- **Verifier correction:** toSpectrum does not throw a catchable exception; line 41's peak_groups_[0] uses the unchecked std::vector::operator[], so calling toSpectrum on an empty/default-constructed DeconvolvedSpectrum is out-of-bounds access (undefined behavior / crash, not a thrown OpenMS exception). The defect and recommended one-line guard fix are otherwise exactly as claimed; this is an undocumented non-empty precondition that the sole internal caller already works around with an explicit empty() check.

### [ANAL-7] MapAlignmentTransformer::transformRetentionTimes(PeakMap&, const TransformationDescription&, bool) — Chromatogram original-RT meta value uses lowercase key "original_rt" instead of the documented "original_RT"
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/source/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.cpp · _mm-alignment_

```cpp
static void transformRetentionTimes(PeakMap& msexp, const TransformationDescription& trafo, bool store_original_rt = false)
```
- **Expectation:** The header documents storeOriginalRT_: 'The original RT is written to a meta value "original_RT"', and the MapRTTransformer TOPP doc states original_RT is created 'for every major data element (spectrum, chromatogram, feature, consensus feature, peptide identification)'. A caller passing store_original_rt=true and later querying chromatograms for meta value "original_RT" expects to find the original retention times.
- **Actual:** Spectra, features, consensus features and peptide IDs all go through storeOriginalRT_ which writes meta value "original_RT" (a single scalar). But the chromatogram branch instead writes a per-point vector under a different, lowercase key "original_rt": `if (store_original_rt && !chromatogram.metaValueExists("original_rt")) { chromatogram.setMetaValue("original_rt", original_rts); }`. So chromatograms never get "original_RT", contradicting both the header doc and the TOPP doc.
- **Evidence:** Line 58-61: `if (store_original_rt && !chromatogram.metaValueExists("original_rt")) { chromatogram.setMetaValue("original_rt", original_rts); }` vs. storeOriginalRT_ (lines 25-26): `if (meta_info.metaValueExists("original_RT")) return false; meta_info.setMetaValue("original_RT", original_rt);` and MapRTTransformer.cpp doc: 'meta values with the name "original_RT" ... for every major data element (spectrum, chromatogram, ...)'.
- **Fix:** Standardize the chromatogram-level key to "original_RT" to match the documentation and the rest of the data structures (a scalar original_RT on the chromatogram, in addition to or instead of the per-point vector). Because external files may already contain the lowercase per-point key, the safe additive fix is to also write/read the documented "original_RT" key (or correct the docs to state chromatograms use a different convention). Source-compatible; only the persisted meta-value key name changes.
- **Verifier correction:** The function signature is unchanged; only a runtime/persisted string-literal meta-key differs. Severity is medium, not high: the original RT data is NOT lost (it is stored, just under the wrong/lowercase key "original_rt" and as a per-point vector rather than a scalar), and no internal code consumes it, so there is no silent corruption of computed results or broken reverse transform. The harm is that any external/downstream consumer following the header and TOPP documentation queries chromatograms for "original_RT" and silently gets an empty result — a documented contract violated without any error. Recoverable once discovered. ABI impact is none (header/signature untouched; only the literal key string changes).

### [ANAL-8] MapAlignmentEvaluationAlgorithmPrecision::evaluate — evaluate() divides by the ground-truth size with no guard, returning NaN/inf instead of signaling when no valid GT consensus features exist
`medium` · `silent-failure` · ABI: `none` · src/openms/source/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmPrecision.cpp · _mm-alignment_

```cpp
void evaluate(const ConsensusMap& consensus_map_in, const ConsensusMap& consensus_map_gt, const double& rt_dev, const double& mz_dev, const Peak2D::IntensityType& int_dev, const bool use_charge, double& out) override
```
- **Expectation:** A caller invoking the precision evaluation expects either a precision value in [0,1] or a clear error if the ground truth contains no usable consensus features (size >= 2).
- **Actual:** The GT map is first filtered to keep only consensus features of size >= 2 into cons_map_gt. If none qualify (or the GT map is empty), cons_map_gt.size() is 0 and the final line `out = (1.0 / double(cons_map_gt.size())) * sum;` performs a division by zero, writing NaN (0/0) or inf to `out` with no exception or warning. The caller silently receives a garbage score.
- **Evidence:** Lines 28-34 build cons_map_gt only from features with size >= 2; line 95: `out = (1.0 / double(cons_map_gt.size())) * sum;` with no check that cons_map_gt is non-empty.
- **Fix:** Add an explicit guard before the division: if cons_map_gt.empty(), either throw Exception::InvalidParameter/MissingInformation, or document and set out to a defined sentinel. Purely additive (an early check), fully ABI-safe.
- **Verifier correction:** evaluate() divides by cons_map_gt.size() at line 95 with no guard. cons_map_gt is the GT map filtered to consensus features of size >= 2; if it ends up empty, sum is 0 and the result is (1.0/0.0)*0.0 = inf*0 = NaN (not literally 0/0), silently written to `out`. The caller receives NaN with no exception or warning. The same unguarded division exists in the Recall sibling (line 100). Trigger requires a degenerate GT input, so severity is medium, not high.

### [ANAL-11] MapAlignmentEvaluationAlgorithm::evaluate / isSameHandle — Three consecutive same-kind numeric tolerance params (rt_dev, mz_dev, int_dev) plus a trailing out-double are trivially swappable at call sites
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h · _mm-alignment_

```cpp
virtual void evaluate(const ConsensusMap&, const ConsensusMap&, const double& rt_dev, const double& mz_dev, const Peak2D::IntensityType& int_dev, const bool use_charge, double& out)
```
- **Expectation:** A caller reading `evaluate(in, gt, a, b, c, true, score)` should be able to tell which tolerance is RT vs m/z vs intensity, and which double is the output.
- **Actual:** evaluate() and isSameHandle() take rt_dev, mz_dev (both double) and int_dev (Peak2D::IntensityType, a floating type) back-to-back, then a bare bool use_charge, then a non-const double& out used as the return slot. Two of the three tolerances are the same type (double) and the intensity type is also floating, so transposing rt_dev/mz_dev compiles silently and produces a wrong-but-plausible score. Returning the result through a trailing out-param rather than a return value compounds the surprise.
- **Evidence:** MapAlignmentEvaluationAlgorithm.h line 37: `virtual void evaluate(const ConsensusMap & conensus_map_in, const ConsensusMap & consensus_map_gt, const double & rt_dev, const double & mz_dev, const Peak2D::IntensityType & int_dev, const bool use_charge, double & out) = 0;` (note also the misspelled parameter name 'conensus_map_in').
- **Fix:** Prefer returning the score as the function result (double evaluate(...)) and/or grouping tolerances into a small struct; at minimum keep the parameter names accurate (fix 'conensus_map_in'). Changing the signature is source-breaking, so the pragmatic fix is documentation plus an additive overload returning the score; flagging the swappable-param hazard honestly.
- **Verifier correction:** The hazard is genuine: rt_dev (double), mz_dev (double) and int_dev (float = Peak2D::IntensityType) are three consecutive implicitly-convertible floating tolerances, so transposing them compiles silently and produces a wrong evaluation score (verified against isSameHandle's per-axis comparisons). Severity is medium, not high: the affected API is a precision/recall evaluation metric reached only through class tests (Precision/Recall derived classes), with no shipping TOPP tool calling evaluate, so the blast radius is a wrong benchmarking number rather than user-data loss/corruption or a crash. ABI impact of flagging is none; the recommended fix (documentation + additive `double evaluate(...)` overload / tolerance struct, plus fixing the `conensus_map_in` typo) is source-compatible, whereas mutating the existing signature in place would be breaking.

### [ANAL-3] FeatureGroupingAlgorithm::group(const std::vector<ConsensusMap>&, ConsensusMap&) — Default ConsensusMap group() silently lossily converts to FeatureMap instead of grouping consensus features
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h · _mm-grouping_

```cpp
virtual void group(const std::vector<ConsensusMap> & maps, ConsensusMap & out)
```
- **Expectation:** Calling group() with a vector<ConsensusMap> on an algorithm that does not specifically override it should either group the consensus features or fail loudly if unsupported.
- **Actual:** The base implementation logs a WARN and converts every input ConsensusMap to a FeatureMap (MapConversion::convert(maps[i], true, fm)), discarding consensus sub-structure, then forwards to the FeatureMap overload. The only signal is an OPENMS_LOG_WARN line that is easy to miss in batch runs; the method returns normally with a degraded result. A caller who passes consensus maps to e.g. FeatureGroupingAlgorithmLabeled (which has no ConsensusMap override) gets silent flattening.
- **Evidence:** FeatureGroupingAlgorithm.cpp: 'OPENMS_LOG_WARN << "FeatureGroupingAlgorithm::group() does not support ConsensusMaps directly. Converting to FeatureMaps."... MapConversion::convert(maps[i], true, fm);' then 'group(maps_f, out);'. Header doc only says 'Algorithms not supporting ConsensusMap input should simply not override this method, as the base implementation will forward the data to the FeatureMap version'.
- **Fix:** The header comment understates the loss. Document explicitly in the public doxygen that the fallback flattens consensus sub-features (lossy) rather than 'forwarding'. Behavior change (throwing) would be breaking; a doc clarification is source-compatible.
- **Verifier correction:** The fallback is not fully "silent" — it emits an OPENMS_LOG_WARN before flattening, so a caller who reads logs is warned. The genuine surprise is (a) the loss is real (ConsensusFeature sub-features are sliced away via BaseFeature::operator= during ConsensusMap->FeatureMap conversion) and (b) the header doxygen euphemistically calls it "forward the data," masking the flatten-then-regroup. The recommended fix is a doxygen clarification only (no signature/behavior change), which is ABI-neutral (none). Throwing instead would be the breaking alternative, but that is not the recommendation. Severity downgraded high->medium because the warning makes it loud/recoverable.

### [ANAL-17] ConsensusMapNormalizerAlgorithmThreshold::computeCorrelation — Silently returns all-1.0 'no normalization' vector when any single map has no surviving ratios
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h · _mm-normalizer_

```cpp
static std::vector<double> computeCorrelation(const ConsensusMap & map, const double & ratio_threshold, const std::string & acc_filter, const std::string & desc_filter)
```
- **Expectation:** If normalization factors cannot be computed for a map, the caller is signaled (exception / status) so it does not silently skip normalization for ALL maps.
- **Actual:** If the surviving-ratio list is empty for ANY one map j, the function abandons the whole computation and returns vector<double>(number_of_maps, 1.0) (a no-op factor for every map), with only a LOG_WARN. The returned 1.0-vector is indistinguishable, to a programmatic caller, from a legitimately-computed factor of 1.0.
- **Evidence:** 'if (ratios.empty()) { OPENMS_LOG_WARN << ... "Result will be unnormalized." ...; return vector<double>(number_of_maps, 1.0); }'
- **Fix:** Document the sentinel prominently (partially done via @warning) and consider an additive overload returning a success flag or an std::optional/bool out-param so callers can detect the 'unnormalized' fallback. Header doc-only change is ABI-safe; the new overload is additive.
- **Verifier correction:** If a single map j produces an empty surviving-ratio list, computeCorrelation returns vector<double>(number_of_maps, 1.0) — disabling normalization for EVERY map, not merely map j — with only a LOG_WARN. The returned unity vector is a valid factor value and thus an undetectable in-band sentinel for a programmatic caller. The behavior is documented via the header @warning, but that doc understates the all-or-nothing cross-map scope. Failure is loud and non-destructive (data left unnormalized, not corrupted), so severity is medium, not high. The proposed fix (header doc clarification + additive overload returning a success flag/optional) is ABI-safe.

### [ANAL-18] ConsensusMapNormalizerAlgorithmMedian::computeMedians — Returns 0 as a failure sentinel that collides with a valid map index, and fills medians with all-1.0 on failure
`medium` · `ambiguous-return-value` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h · _mm-normalizer_

```cpp
static Size computeMedians(const ConsensusMap & map, std::vector<double> & medians, const std::string& acc_filter, const std::string& desc_filter)
```
- **Expectation:** The return value is documented as 'index of map with largest number of features'; a caller treats it purely as an index. There is no documented way to detect failure from the return value.
- **Actual:** On the not-enough-features path the function returns 0 (which is also the legitimate index of the first map) AND sets every median to 1.0, with only a LOG_WARN. A caller that uses the returned index and the medians vector cannot distinguish 'map 0 is the largest' from 'computation failed, no normalization'.
- **Evidence:** 'medians[j] = 1.0; ... if (!enough_features_left) { OPENMS_LOG_WARN << ... "Result will be unnormalized." ...; return 0; }'. Header doc: '@return index of map with largest number of features'.
- **Fix:** Document the 0-with-all-1.0-medians fallback in the header @return, and prefer an additive overload that reports failure explicitly (bool return or out-status). Doc change is ABI-safe; new overload additive.
- **Verifier correction:** The return value 0 is an overloaded sentinel that collides with a valid map index, and the header @return documents only the success meaning — confirmed and real. But the failure is not silent: it emits OPENMS_LOG_WARN and leaves medians observably all-1.0, which makes normalization a no-op (intensities unnormalized, not corrupted). The only in-tree non-Python caller (normalizeMaps) is correct by construction because all-1.0 medians make any index value a no-op. The risk is to external API consumers (notably the pyOpenMS binding that surfaces (idx, medians)). Severity is medium (loud + recoverable), category is ambiguous/overloaded-return-value rather than silent-failure. Recommendation stands: document the 0-with-all-1.0-medians fallback in the header @return (ABI-safe) and add an additive overload reporting failure explicitly.

### [ANAL-21] ConsensusMapNormalizerAlgorithmQuantile::resample — resample dereferences data_in.front()/back() with no empty-input guard; empty source + nonzero points is undefined behavior
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h · _mm-normalizer_

```cpp
static void resample(const std::vector<double> & data_in, std::vector<double> & data_out, UInt n_resampling_points)
```
- **Expectation:** A public resample utility, documented to handle n_resampling_points == 0, would be expected to also handle an empty data_in gracefully (or document that it must be non-empty).
- **Actual:** Only n_resampling_points == 0 is guarded. With a non-empty request but empty data_in, the code calls data_in.front() and data_in.back() unconditionally -> out-of-bounds access / undefined behavior. The header documents the n==0 case but never states data_in must be non-empty. This is directly reachable from the pyOpenMS-exposed static.
- **Evidence:** 'data_out.clear(); data_out.resize(n_resampling_points); if (n_resampling_points == 0) { return; } data_out[0] = data_in.front(); data_out[n_resampling_points - 1] = data_in.back();' -- no check for data_in.empty().
- **Fix:** Add a precondition guard (throw Exception::InvalidValue or early-return for empty data_in) and document it. Body-only fix; ABI-safe.
- **Verifier correction:** resample() (src/openms/source/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.cpp:100-129) guards only n_resampling_points==0; with empty data_in and n>0 it calls data_in.front()/back() on an empty std::vector, which is undefined behavior (not a thrown exception — the original "surprising-throw" category is really "missing precondition leading to UB/crash"). The header doc explicitly enumerates the n==0 case but omits any data_in-non-empty precondition. It is reachable (a) from Python: the static is bound via nanobind at src/pyOpenMS/bindings/bind_analysis.cpp:388 as ConsensusMapNormalizerAlgorithmQuantile.resample(data_in, n_resampling_points), so resample([], 5) crashes; and (b) internally from normalizeMaps (line 45) when one map is empty while another is non-empty. Fix: early-return or throw Exception::InvalidValue when data_in.empty() && n_resampling_points>0, and document the precondition. Body-only change, ABI-safe.

### [ANAL-22] ConsensusMapNormalizerAlgorithmQuantile::resample / extractIntensityVectors / setNormalizedIntensityValues — Internal pipeline helpers are exposed as public statics with stateful lock-step preconditions that are easy to violate
`medium` · `other` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h · _mm-normalizer_

```cpp
static void resample(...); static void extractIntensityVectors(...); static void setNormalizedIntensityValues(...)
```
- **Expectation:** Public static methods of a normalizer are expected to be self-contained, order-independent operations a caller can use standalone.
- **Actual:** extractIntensityVectors and setNormalizedIntensityValues only work correctly when called in the exact same consensus-feature iteration order, and setNormalizedIntensityValues silently consumes per-map vectors positionally (feature_ints[map_idx][progress_indices[map_idx]++]). A caller mixing them, or reordering the map between calls, gets silently wrong intensities with no error. The header itself warns the caller is responsible for lock-step ordering.
- **Evidence:** setNormalizedIntensityValues comment: 'assumes the input map and feature_ints are in the same order as in the beginning'; body: 'double intensity = feature_ints[map_idx][progress_indices[map_idx]++];'. Header: 'The caller is responsible for keeping @p feature_ints in lock-step with the iteration order used by extractIntensityVectors'.
- **Fix:** Document these as internal helpers (or move to a detail namespace in a future major version). Minimum: keep header warnings (already present) and avoid promoting them as general-purpose. ABI-safe doc-only now; relocation is breaking.
- **Verifier correction:** The lock-step / positional-consumption hazard applies to the extractIntensityVectors + setNormalizedIntensityValues PAIR only, not to resample(), which is a self-contained order-independent pure function. setNormalizedIntensityValues positionally drains feature_ints[map_idx] in consensus-feature iteration order (cpp:168), so calling it standalone with reordered or independently-built vectors yields silently wrong intensities (matching-shape reorder) or out-of-bounds UB (size mismatch). The C++ header documents the precondition thoroughly (h:117-123, cpp:157-158), but the pyOpenMS def_static bindings (bind_analysis.cpp:393-394) re-expose both as general-purpose operations with one-line docstrings that omit the lock-step warning, which is the real misuse surface. Severity is medium not high: the genuine public API (normalizeMaps) always uses the helpers in correct lock-step, and misuse requires deliberate standalone use of explicitly-warned low-level helpers. Safe fix (clarify docs / mark internal) is ABI-safe doc-only (none); relocating to a detail namespace would be source-breaking.

### [ANAL-4] FeatureDistance::operator() — Distance functor operator() mutates internal state (norm_factor) and is therefore not const / not thread-safe
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h · _mm-pairfinder_

```cpp
std::pair<bool, double> operator()(const BaseFeature & left, const BaseFeature & right)
```
- **Expectation:** A distance functor's evaluation operator is a pure read-only computation: given two features it returns a distance without altering the functor. Callers (e.g. QTClusterFinder caching distances, parallel grouping) would reasonably assume operator() is reentrant and could be const, so the same FeatureDistance instance can be shared across threads.
- **Actual:** operator() is non-const and, when m/z tolerance is in ppm, writes back into the object: 'max_diff_mz *= left_mz * 1e-6; // overwrite this parameter - it will be recomputed each time anyway: params_mz_.norm_factor = 1 / max_diff_mz;' (FeatureDistance.cpp:158-160). Each call mutates params_mz_.norm_factor based on the left feature's m/z, so concurrent calls race and the functor cannot be safely shared between threads.
- **Evidence:** FeatureDistance.cpp:156-161: 'if (params_mz_.max_diff_ppm) { max_diff_mz *= left_mz * 1e-6; params_mz_.norm_factor = 1 / max_diff_mz; }'. Header declares non-const: 'std::pair<bool, double> operator()(const BaseFeature & left, const BaseFeature & right);' (FeatureDistance.h:92-93).
- **Fix:** Compute the per-call ppm norm factor as a local variable passed into distance_() instead of writing it back into params_mz_; then mark operator() const. The local-variable refactor is source-compatible, but flipping the member function to const changes its mangled name (ABI-breaking), so for strict ABI either keep non-const and just remove the mutation, or add a new const overload.
- **Verifier correction:** The mutation and non-const declaration are exactly as claimed and constitute a real hidden side effect (operator() writes params_mz_.norm_factor based on the left feature's m/z and reads it back via distance_). However, the thread-safety hazard is currently latent: no existing caller runs operator() concurrently on a shared instance, and single-threaded results are correct (norm_factor is overwritten before every read). Hence severity is medium (invites misuse / blocks safe parallelization and forced a const_cast workaround at FeatureGroupingAlgorithmKD.cpp:456), not high. ABI: the substantive fix — compute the ppm norm factor as a local and pass it into distance_() instead of writing into params_mz_ — removes the side effect and is source- and ABI-compatible. Additionally marking operator() const (or adding a const overload) would change the mangled name and is ABI-breaking, so that polish is optional/separable; the core defect can be remediated without breaking ABI.

### [ANAL-14] TransformationModel::checkDatumRange — checkDatumRange is named like a validator but silently truncates the value and is non-const
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h · _mm-transform_

```cpp
double checkDatumRange(const double& datum, const double& datum_min, const double& datum_max)
```
- **Expectation:** A method named check...() returns a bool / throws on out-of-range and does not alter data; 'check' implies validation, not mutation.
- **Actual:** checkDatumRange returns the datum clamped to [min,max] (truncates to the bound and only logs at INFO). It is the clamping primitive used inside weightData, so out-of-range coordinates are silently snapped to bounds rather than reported. It is also non-const despite reading no member state.
- **Evidence:** TransformationModel.cpp: "if (datum >= datum_max) { ... datum_checked = datum_max; } else if (datum <= datum_min) { ... datum_checked = datum_min; } return datum_checked;" Header doc admits: 'If the datum is below the min bounds, the min bound is returned.'
- **Fix:** ABI-risky to rename. Document clearly that it CLAMPS (rename ideal: clampDatumToRange). Mark const (it uses no members) — that is source-compatible. Elevate the truncation log above INFO so silent clamping is visible.
- **Verifier correction:** checkDatumRange is a clamp, not a validator: it returns the datum clamped to [min,max] (double return, not bool) and weightData writes that clamped value back into the data, so out-of-range coordinates are snapped to the user-configured bounds. It is also non-const despite using no members. Mitigating: the behavior IS documented in the header ("the min/max bound is returned"), the clamping only triggers when the user opts into x/y_datum_min/max bounds with weighting enabled, and it is logged at INFO. The surprise is real (name implies validation, not mutation/clamping) but severity is medium, not high — no silently-wrong default-path results, recoverable, logged. Recommendation: rename to clampDatumToRange (ABI-breaking, so do it on a major bump) and/or mark const (source-compatible, though it changes the symbol mangling); raise the truncation log above INFO. The abi_impact of the minimal/preferred fix (const) is source-compatible; a rename would be breaking.

### [ANAL-41] NuXLMarkerIonExtractor::extractMarkerIons — extractMarkerIons returns normalized (relative) intensities, not the spectrum's actual intensities
`medium` · `unit-or-index` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLMarkerIonExtractor.h · _nuxl-core_

```cpp
static MarkerIonsType extractMarkerIons(const PeakSpectrum& s, const double marker_tolerance)
```
- **Expectation:** A caller reading 'extract ... RNA marker ions' from a const PeakSpectrum expects the returned (m/z, intensity) pairs to carry the marker peaks' actual measured intensities found in the input spectrum.
- **Actual:** The function makes an internal copy 'PeakSpectrum spec(s)' and runs 'Normalizer normalizer; normalizer.filterSpectrum(spec);' BEFORE searching for marker peaks (NuXLMarkerIonExtractor.cpp lines 30-32). The intensity stored in the result (it->second[i].second = max_intensity) is therefore a NORMALIZED relative value, not the raw intensity. NuXLReport.cpp consequently multiplies by 100 to print a percentage (line 60: it->second[i].second * 100.0). The name/signature give no hint that intensities are rescaled.
- **Evidence:** NuXLMarkerIonExtractor.cpp: `PeakSpectrum spec(s);` `Normalizer normalizer;` `normalizer.filterSpectrum(spec);` then `it->second[i].second = max_intensity;` over the normalized `spec`.
- **Fix:** Rename to extractNormalizedMarkerIons or document in the header that returned intensities are Normalizer-relative (and against which normalization mode). Additively, expose an overload/flag that returns raw intensities. ABI: a doc/comment fix is non-breaking; adding an overload is source-compatible.
- **Verifier correction:** The factual claim is accurate and fully supported by the code. The only adjustment is severity: medium rather than high. Returned marker-ion intensities are Normalizer "to_one" relative values (peak intensity divided by the spectrum's maximum peak, range [0,1]), not the input spectrum's actual intensities. They are not silently corrupt — they are well-defined relative fractions, and the sole in-tree consumer (NuXLReport.cpp) correctly multiplies by 100 to display percentages. The surprise primarily affects a new/external caller who assumes raw intensities. Minimal fix (header doc comment stating intensities are Normalizer-to_one-relative) is ABI-non-breaking (none); the optional rename/overload would be source-compatible.

### [ANAL-45] NuXLDeisotoper::deisotopeAndSingleCharge — 16-parameter function with 9 same-typed bool/unsigned flags is unreadable and trivially mis-ordered at call sites
`medium` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLDeisotoper.h · _nuxl-core_

```cpp
static void deisotopeAndSingleCharge(MSSpectrum& spectrum, double fragment_tolerance, bool fragment_unit_ppm, int min_charge = 1, int max_charge = 3, bool keep_only_deisotoped = false, unsigned int min_isopeaks = 3, unsigned int max_isopeaks = 10, bool make_single_charged = true, ...)
```
- **Expectation:** A caller can configure deisotoping without silently swapping flags; adjacent same-typed parameters should be hard to transpose.
- **Actual:** The signature has 16 parameters including a run of interchangeable bools (keep_only_deisotoped, make_single_charged, annotate_charge, annotate_iso_peak_count, use_decreasing_model, add_up_intensity, annotate_features, preserve_high_intensity_peaks) and two adjacent unsigned ints (min_isopeaks, max_isopeaks). A call like deisotopeAndSingleCharge(s, 10, true, 1, 3, true, ...) is opaque and a swapped bool compiles silently with a completely different result (e.g. keep_only_deisotoped vs make_single_charged).
- **Evidence:** NuXLDeisotoper.h: parameter list spans `bool keep_only_deisotoped = false, unsigned int min_isopeaks = 3, unsigned int max_isopeaks = 10, bool make_single_charged = true, bool annotate_charge = false, bool annotate_iso_peak_count = false, bool use_decreasing_model = true, ... bool add_up_intensity = false, bool annotate_features = false, bool preserve_high_intensity_peaks = false`.
- **Fix:** Introduce a parameter/options struct (additive overload) so flags are named at the call site; keep the existing signature for ABI. ABI: adding a struct overload is source-compatible.
- **Verifier correction:** The function has 17 parameters (not 16; the claim omits the trailing `double preserve_low_mz_peaks_threshold = -1e10`). It contains 9 bool flags and 3 unsigned int flags. The genuinely adjacent, silently-transposable same-typed runs are: min_isopeaks/max_isopeaks (two adjacent unsigned int), annotate_charge/annotate_iso_peak_count (two adjacent bool), and add_up_intensity/annotate_features/preserve_high_intensity_peaks (three adjacent bool). The claim's headline pair (keep_only_deisotoped vs make_single_charged) are both bools in the same call and so confusable, but they are not adjacent (positions 6 and 9). Severity is medium rather than high: it is an internal NuXL-specific static helper with a single real call site (OpenNuXL.cpp:3088), and the function loudly validates tolerance and min/max isopeaks, so misuse invites wrong results but is recoverable rather than a silent data-corruption hazard across all reasonable uses. The recommended fix (an additive options-struct overload preserving the existing signature) is source-compatible, not breaking.

### [ANAL-48] NuXLProteinReport::getCrossLinkEfficiency — getCrossLinkEfficiency filters cross-links by comparing the integer meta value NuXL:isXL against the string "false"
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLReport.h · _nuxl-core_

```cpp
static std::map<char, double> getCrossLinkEfficiency(const PeptideIdentificationList& peps)
```
- **Expectation:** A skip-condition for 'is not a cross-link' should compare against the same representation the value is stored in. Elsewhere in the same module NuXL:isXL is read as an int (NuXLFDR uses static_cast<int>(...) == 0).
- **Actual:** The function skips hits with `ph.getMetaValue("NuXL:isXL") == "false"` (NuXLReport.cpp line 285) -- a string literal compare. Since NuXL:isXL is populated/consumed as an integer (0/non-0) in NuXLFDR.cpp (`static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0`), this string comparison very likely never matches, so the intended 'skip non-XL hits' filter is silently a no-op, and plain peptides may be counted as cross-links in the efficiency denominator/numerator.
- **Evidence:** NuXLReport.cpp: `if (ph.isDecoy() || ph.getMetaValue("NuXL:isXL") == "false") continue;` vs NuXLFDR.cpp: `if (static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0)`.
- **Fix:** Compare consistently with the integer form: `static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0`. Behavior fix, source-compatible. ABI: source-compatible.
- **Verifier correction:** The string compare `getMetaValue("NuXL:isXL") == "false"` is provably a dead no-op: NuXL:isXL is stored as INT_VALUE (bool promoted to DataValue(int); no bool ctor exists; test idXML shows type="int"), and DataValue::operator== short-circuits to false on the INT_VALUE-vs-STRING_VALUE type mismatch (DataValue.cpp:759,780). The intended 'skip non-XL hits' filter never executes. However, the claim's stated impact (plain peptides polluting the efficiency numerator/denominator) is mostly masked by the subsequent guard `if (best_localization >= 0)` (line 288), since non-XL hits default to best_localization_position = -1. Real bug, but practical numeric impact is conditional rather than guaranteed — hence medium. Fix: `static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0`.

### [ANAL-49] NuXLAnnotatedHit::best_localization_position — On a tie, best_localization_position stores the LAST tied residue, not the first or a documented winner
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLAnnotatedHit.h · _nuxl-core_

```cpp
int best_localization_position = -1; // UNKNOWN
```
- **Expectation:** A field named best_localization_position is expected to hold the single best-scoring cross-link site; on ties a caller would expect a defined, documented choice (commonly the first index).
- **Actual:** Per the NuXLAnnotateAndLocate.h documentation of the producing code, ties 'resolve to the last index, not the first' and the implementation even carries a caveat 'check if there are situations where multiple have the same score'. So best_localization_position is an order-dependent, last-tie value. NuXLReport.cpp and NuXLReport getProteinReportEntries then use this single position verbatim to place the cross-link in the protein (e.g. char c = aas.toUnmodifiedString()[best_localization]; xl_pos_in_protein = peptide_start_in_protein + best_localization), so the arbitrary tie resolution silently determines the reported residue.
- **Evidence:** NuXLAnnotateAndLocate.h: 'lowercases every residue that ties for the maximum score ... and stores the last-tied position in best_localization_position (... ties resolve to the last index, not the first).'
- **Fix:** Document the tie semantics on the field and ideally expose all tied positions (best_localization already lowercases all ties); downstream protein reporting should not assume a unique site. ABI: adding fields/docs is source-compatible.
- **Verifier correction:** Claim is accurate in substance but slightly overstated in impact. best_localization_position deterministically stores the LAST residue tied (within 1e-6) for the maximum site score, an order-dependent and undocumented-at-the-field choice that propagates verbatim into protein-level cross-link placement (NuXLReport.cpp lines 291, 442). It is a real silent-failure surface, but bounded: the companion best_localization string already lowercases ALL tied sites, ties at the exact float maximum are uncommon in real data, and the worst case is mislocalizing the site by a few residues within a single peptide rather than corrupting unrelated results. Severity therefore medium (silently wrong but recoverable from the sibling field), not high. Recommendation stands: document tie semantics on the field and have downstream protein reporting consume the all-ties information instead of the single scalar.

### [ANAL-52] NuXLParameterParsing::getTargetNucleotideToFragmentAdducts — getTargetNucleotideToFragmentAdducts registers N-/C-terminal modifications in the global ModificationsDB
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLParameterParsing.h · _nuxl-frag_

```cpp
static NucleotideToFragmentAdductMap getTargetNucleotideToFragmentAdducts(StringList fragment_adducts);
```
- **Expectation:** A 'get...Map' parser that turns a StringList into a map of nucleotide->fragment-adduct definitions is expected to be a pure transformation with no global side effects.
- **Actual:** For every parsed fragment adduct it also creates and inserts two ResidueModification entries (N-term and C-term) into the global ModificationsDB singleton via addModification, permanently extending the process-wide modification database. Calling the function twice or from another context silently leaves new modifications behind.
- **Evidence:** 'if (!ModificationsDB::getInstance()->has(name)) { ... ModificationsDB::getInstance()->addModification(std::move(c_term)); ... ModificationsDB::getInstance()->addModification(std::move(n_term)); }' (NuXLParameterParsing.cpp:210-227). The header comment describes only parsing: 'Parse tool parameter to create map from target nucleotide to all its fragment adducts'.
- **Fix:** Document the ModificationsDB registration in the header (@note registers each adduct as N-/C-terminal mods in the global ModificationsDB). Preferred additive fix: split out a registerFragmentAdductModifications(...) call so the getter stays pure; keep the combined behavior available for existing callers during deprecation.
- **Verifier correction:** getTargetNucleotideToFragmentAdducts is not a pure parser: for every parsed fragment adduct it registers an N-terminal and a C-terminal ResidueModification (id/name = adduct name, fullId = "<name> (N-term)"/"<name> (C-term)", diffMonoMass = adduct mass) into the global ModificationsDB singleton. This is guarded by `has(name)` and addModification is idempotent (dedupes by fullId), so the effect does not accumulate duplicates or corrupt existing reference modifications, and does not cause wrong results or crashes. The surprise is that a 'get...Map' getter mutates a process-wide singleton, undocumented in the header. Recommended additive fix: factor the registration into a separate registerFragmentAdductModifications(...) so the getter stays pure, and add a header @note documenting the global ModificationsDB registration.

### [ANAL-53] NuXLPresets::getAllPresetsNames — getAllPresetsNames returns an empty list (not an error) when the presets file is missing or unparseable
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLPresets.h · _nuxl-frag_

```cpp
StringList getAllPresetsNames(const std::string& custom_presets_file = "");
```
- **Expectation:** A function returning the set of available preset names should distinguish 'this installation has zero presets' from 'the presets file could not be found or parsed'. The header documents only '@return StringList containing all available preset names'.
- **Actual:** On a missing file, or on any JSON parse exception, the function only emits a LOG_WARN and returns an empty StringList. A caller that branches on whether a chosen preset is in the returned list cannot tell a genuinely empty preset set from a load failure, and will silently treat every preset name as 'unknown'.
- **Evidence:** 'catch (const std::exception& e) { OPENMS_LOG_WARN << "Error reading presets from " ... } ... else { OPENMS_LOG_WARN << "Presets file not found: " ... } return presets;' (NuXLPresets.cpp:55-64). presets is left empty in both failure branches.
- **Fix:** Document the failure-returns-empty behavior in the header, or add an overload/out-param signalling load success vs an empty preset set. (Note: the sibling getPresets does throw std::runtime_error on a missing/unparseable file, so the two functions handle the same error inconsistently — an additional inconsistent-convention smell worth aligning.) Additive overload keeps ABI stable.
- **Verifier correction:** The behavior is real but severity is medium, not high. The function returns an empty StringList both when the presets file is missing and when JSON parsing throws, emitting only OPENMS_LOG_WARN, and the header documents no failure contract — so callers cannot distinguish 'no presets' from 'load failure'. However this is loud (warning logged) and the single in-repo caller's downstream getPresets() throws on the resulting unknown preset, so it does not produce silently-wrong results or data loss in practice. The inconsistency with the sibling getPresets() (which throws on the same missing-file/parse-error/unknown-preset conditions) is confirmed and reinforces that empty-on-failure is an oversight, not an intended convention. Recommended fix: document the empty-on-failure behavior in the header (ABI: none) and/or align error handling with getPresets via an additive overload or out-param (ABI: source-compatible).

### [ANAL-34] AbsoluteQuantitation::applyCalibration — applyCalibration silently clamps negative back-calculated concentration to 0.0
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h · _quant-absolute-mdv_

```cpp
double applyCalibration(const Feature & component, const Feature & IS_component, const std::string & feature_name, const std::string & transformation_model, const Param & transformation_model_params)
```
- **Expectation:** applyCalibration is documented to '@return The absolute concentration' (and to throw UnableToFit on failure). A caller expects the model's back-calculated value, including the ability to detect an out-of-range / below-baseline reading.
- **Actual:** Any negative back-calculated concentration is silently replaced by 0.0, so a strongly negative (i.e. out-of-calibration / failed) reading and a true-zero reading become indistinguishable. This clamp is undocumented in the header.
- **Evidence:** if (calculated_concentration < 0.0) { calculated_concentration = 0.0; } return calculated_concentration; (AbsoluteQuantitation.cpp:271-276).
- **Fix:** Document the clamp in the header, and consider returning the raw value (or flagging below-LLOQ) so callers can detect out-of-range readings. Doc-only addition is source- and ABI-compatible; changing the clamp behavior is a behavior change to gate behind a parameter.
- **Verifier correction:** applyCalibration silently clamps any negative back-calculated concentration to 0.0 (AbsoluteQuantitation.cpp:271-274) before returning, and this clamp is undocumented in the header (AbsoluteQuantitation.h:250 only states '@return The absolute concentration'). This makes an out-of-calibration/below-baseline negative reading indistinguishable from a genuine zero. Severity is medium rather than high: clamping unphysical negative concentrations is a domain-defensible convention and the value is pinned to a boundary (information loss, recoverable) rather than corrupted into a plausibly-wrong number. Recommend documenting the clamp; optionally gate raw-value/below-LLOQ flagging behind a parameter (a behavior change, not a doc-only change).

### [ANAL-35] IsotopeLabelingMDVs::isotopicCorrection / isotopicCorrections — correction_matrix_agent documented as @param[out] but is a const input; supplied correction_matrix silently ignored when an agent is given
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h · _quant-absolute-mdv_

```cpp
void isotopicCorrection(const Feature& normalized_feature, Feature& corrected_feature, const Matrix<double>& correction_matrix, const DerivatizationAgent& correction_matrix_agent)
```
- **Expectation:** Doxygen '@param[out] correction_matrix_agent' tells a reader this is an output. And '@param[in] correction_matrix' tells the caller their matrix will be used for the correction.
- **Actual:** correction_matrix_agent is a 'const DerivatizationAgent&' input that selects an internally hard-coded matrix; it is never written. Worse, when correction_matrix_agent != NOT_SELECTED and is found (e.g. TBDMS), the caller-supplied correction_matrix is completely ignored in favor of the built-in matrix — the @param[in] correction_matrix is silently overridden.
- **Evidence:** Header: '@param[out]  correction_matrix_agent name of the derivatization agent ...'. Source: parameter is 'const DerivatizationAgent& correction_matrix_agent'; and 'if (correction_matrix_agent != DerivatizationAgent::NOT_SELECTED && correction_matrix_search != correction_matrices.end()) { ... uses built-in matrix ... } else { ... uses supplied correction_matrix ... }'.
- **Fix:** Fix the Doxygen to '@param[in]' for correction_matrix_agent, and document that a non-NOT_SELECTED agent overrides the supplied correction_matrix (or, ideally, validate/throw if both a non-identity matrix and a known agent are supplied). Doc fix is source- and ABI-compatible.
- **Verifier correction:** correction_matrix_agent is documented @param[out] (IsotopeLabelingMDVs.h:63 and :83) but is declared const DerivatizationAgent& — a read-only input, never written; the [out] direction is simply wrong. It selects an internally hard-coded matrix (only TBDMS, .cpp:60-67); when non-NOT_SELECTED and found, the caller-supplied correction_matrix is ignored (.cpp:73-95). That override is partially signaled in the param prose ('the internally stored correction matrix if the name of the agent is supplied'), so it is not fully undocumented as originally claimed; the clearly indefensible bug is the @param[out] mislabel. Fix: change to @param[in] and state explicitly that a known agent overrides correction_matrix. Comment-only fix; source- and ABI-compatible.

### [ANAL-36] KDTreeFeatureMaps::getNeighborhood vs queryRegion — Sibling range queries disagree on whether the result vector is appended or cleared
`medium` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h · _quant-absolute-mdv_

```cpp
void getNeighborhood(..., std::vector<Size>& result_indices, ...) const; void queryRegion(..., std::vector<Size>& result_indices, ...) const;
```
- **Expectation:** Two adjacent query methods on the same class that fill an out 'result_indices' vector should follow the same convention (both clear, or both append).
- **Actual:** getNeighborhood APPENDS to result_indices ('Results are appended to result_indices (the existing contents are preserved).'), while queryRegion REPLACES it ('Replaces (does not append to) result_indices.' / 'cleared first, then filled'). A caller who reuses one vector across both calls, or swaps one call for the other, will silently get wrong (duplicated or lost) results.
- **Evidence:** Header getNeighborhood: 'Results are appended to @p result_indices (the existing contents are preserved).' vs queryRegion: 'Replaces (does not append to) @p result_indices.'
- **Fix:** Make the convention uniform (clear-then-fill is the least-surprising default), or rename to encode the behavior (e.g. appendNeighborhood). At minimum, raise the discrepancy to the top of each method's brief. Renaming/behavior change affects callers; the additive safe step is to clearly document and add a clear-first variant.
- **Verifier correction:** getNeighborhood and queryRegion genuinely disagree on the out-vector convention (append vs clear-then-fill), and the inconsistency is real and reproducible in both code and docs — notably getNeighborhood wraps queryRegion yet appends where the wrappee clears. But the behavior is NOT an undocumented landmine: both conventions are explicitly stated in the Doxygen brief and reinforced via [in,out] (getNeighborhood) vs [out] (queryRegion) parameter tags, and every existing caller passes a fresh empty vector so no current code is affected. The hazard is a plausible-but-documented future misuse (reusing one vector across both calls, or assuming symmetric semantics), which yields duplicated/lost indices silently — recoverable and warned-against, hence medium not high. ABI is unaffected (signatures unchanged); the safe additive fix is documentation emphasis or a clear-first variant.

### [ANAL-38] AbsoluteQuantitation::calculateBias — calculateBias divides by actual_concentration with no zero guard
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h · _quant-absolute-mdv_

```cpp
double calculateBias(const double & actual_concentration, const double & calculated_concentration)
```
- **Expectation:** A documented '@exception Exception::UnableToFit' and a percentage return suggest a defined result or an error for invalid input.
- **Actual:** Bias is computed as fabs(actual - calculated)/actual*100 with no check for actual_concentration == 0, silently yielding inf/nan, and the documented UnableToFit exception is never thrown.
- **Evidence:** double bias = fabs(actual_concentration - calculated_concentration)/actual_concentration*100; return bias; (AbsoluteQuantitation.cpp:153-157).
- **Fix:** Guard against actual_concentration == 0 (throw the documented exception or return NaN with a documented contract) and correct the header @exception note. Adding a guard is source- and ABI-compatible.
- **Verifier correction:** The quoted code and the spurious '@exception Exception::UnableToFit' doc are both confirmed. But the title overstates impact: actual_concentration==0 is a degenerate input never produced in normal calibration (standards have non-zero known concentrations; lowest test point is 0.01), and the failure mode is a loud inf/nan, not a silently-wrong finite result. The genuinely actionable defect is the header documenting an exception the function never throws and the absence of a defined contract for zero input. Fix: correct the header @exception note (the adjacent calculateBiasAndR uses '@exception None') and either guard actual_concentration==0 or document the inf/nan contract. Such a change leaves the signature unchanged, so it is source- and ABI-compatible (abi_impact none).

### [ANAL-39] IsotopeLabelingMDVs::calculateMDV — NORM_SUM normalization of 'intensity' lacks the zero-sum guard present in every sibling branch
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h · _quant-absolute-mdv_

```cpp
void calculateMDV(const Feature& measured_feature, Feature& normalized_feature, const MassIntensityType& mass_intensity_type, const std::string& feature_name)
```
- **Expectation:** calculateMDV normalizes intensities to sum to 1. With an all-zero/empty feature the result should be well-defined (zeros) and consistent across the feature_name=='intensity' and metavalue branches.
- **Actual:** Three of the four branches guard the divisor (NORM_MAX intensity, NORM_MAX metavalue, NORM_SUM metavalue all check '!= 0.0'), but the NORM_SUM + feature_name=='intensity' branch divides by feature_peak_apex_intensity_sum unconditionally, silently producing inf/nan when the sum is zero. The inconsistency makes behavior depend on which feature_name string the caller passes.
- **Evidence:** NORM_SUM intensity branch: 'normalized_feature.getSubordinates().at(...).setIntensity((it->getIntensity() / feature_peak_apex_intensity_sum));' with no surrounding 'if (feature_peak_apex_intensity_sum != 0.0)', unlike the metavalue branch right below it which has the guard.
- **Fix:** Add the same 'if (sum != 0.0)' guard to the NORM_SUM intensity branch for consistency. Source- and ABI-compatible bug fix.
- **Verifier correction:** Severity is medium rather than high: the divide-by-zero only triggers for a degenerate all-zero/empty-intensity feature (not the normal data path), so it does not corrupt typical runs. But within that edge case it is a genuine silent failure — NaN/inf is written and propagated downstream rather than producing the well-defined zeros the three guarded sibling branches yield. Fix is to add `if (feature_peak_apex_intensity_sum != 0.0)` around the NORM_SUM intensity loop (lines 322-325), matching the metavalue branch; this is source- and ABI-compatible.

### [ANAL-24] IsobaricChannelExtractor::extractSingleSpec — channel_qc out-param doc mislabels both pair members (m/z delta vs m/z, peak count vs channel index)
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h · _quant-isobaric_

```cpp
std::vector<double> extractSingleSpec(Size spec_idx, const MSExperiment& exp, std::vector<std::pair<double, unsigned>>& channel_qc)
```
- **Expectation:** Per the header: '@param[out] channel_qc vector of pairs of m/z and channel index'. A caller would read pair.first as an absolute m/z and pair.second as the channel index.
- **Actual:** The implementation stores pair.first = mz_delta (the signed distance between expected and observed reporter m/z, NOT an absolute m/z) and pair.second = peak_count (the number of peaks found in the window, NOT a channel index). The channel index is implicit in the vector position (map_index), so storing it in .second would be redundant anyway.
- **Evidence:** IsobaricChannelExtractor.cpp:802 `channel_qc[map_index].second = peak_count;` and :807 `channel_qc[map_index].first = mz_delta;` where :805 `double mz_delta = reporter_mz - idx_nearest->getMZ();`. Header IsobaricChannelExtractor.h:85 `@param[out] channel_qc vector of pairs of m/z and channel index`.
- **Fix:** Fix the Doxygen to say '.first = signed m/z delta (expected - observed), .second = peak count within the search window'. Additively, consider a named struct (like the existing ChannelQC) instead of std::pair<double,unsigned> so the member meanings are self-documenting; the std::pair signature change would be ABI-breaking, so only the doc fix is safe now.
- **Verifier correction:** Both documented member meanings are wrong: .first is a SIGNED m/z delta (expected reporter_mz minus observed nearest-peak m/z), not an absolute m/z; .second is the peak count within the ~0.5 Da search window (>1 flags non-unique signal), not a channel index (which is the implicit vector position map_index). Severity is medium rather than high because channel_qc carries diagnostic QC statistics, not the quantitation intensities (those are the return value); a misled caller produces wrong QC reporting but does not corrupt the actual quantification, and the one existing caller already uses the correct semantics. Safe fix is doc-only (abi_impact none); the recommended named-struct refactor replacing std::pair<double,unsigned> would be ABI-breaking and should be deferred.

### [ANAL-25] IsobaricChannelExtractor::printStats — 'printStats' takes its stats by non-const ref and reorders/mutates the caller's data while only printing
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h · _quant-isobaric_

```cpp
void printStats(ChannelQCSet& stats) const; void printStatsWithMissing(std::vector<ChannelQC>& stats) const
```
- **Expectation:** A method named printStats that is declared const is expected to read the stats and emit log output without altering the caller's data.
- **Actual:** printStats(ChannelQCSet&) feeds each channel's mz_deltas vector to Math::median(begin, end, false), which std::sorts the range in place, permanently reordering the caller's stored deltas. printStatsWithMissing goes further: it erases NaN entries (std::remove_if/erase) and rewrites every delta to its absolute value (std::transform with std::abs) in the caller's vector. The const qualifier is technically honored only because the argument (not *this) is mutated.
- **Evidence:** IsobaricChannelExtractor.cpp:843 `Math::median(stats[cl_it->name].mz_deltas.begin(), ... , false)`; StatisticFunctions.h:138-141 sorts in place when sorted==false. printStatsWithMissing: cpp:888-890 erase NaN, :897 `std::transform(... std::abs ...)` rewriting cur_deltas in place. Header comment IsobaricChannelExtractor.h:110 admits '@param[in] stats the stats to print (NOT const, since we need to sort it for median calculation)'.
- **Fix:** Either copy each delta vector before computing the median (so the caller's data is untouched), or rename to make the mutation explicit (e.g. summarizeAndConsumeStats). At minimum document on printStats(ChannelQCSet&) that it sorts the passed deltas in place, as is already done for printStatsWithMissing. Taking a copy is source/ABI compatible.
- **Verifier correction:** printStats(ChannelQCSet&) const sorts the caller's mz_deltas in place (via Math::median(...,false) -> std::sort) and is entirely undocumented regarding this mutation. printStatsWithMissing(std::vector<ChannelQC>&) const additionally erases NaN entries (data loss) and rewrites all deltas to absolute values; its header comment warns only about sorting, understating the destructive erase/abs operations. The const qualifier protects *this, not the argument. Severity is medium (not high) because the sole in-tree caller discards the vector immediately after the call, so no current results are silently wrong; the risk is to future/reuse callers. Fix: copy each delta vector internally before computing the median (source- and ABI-compatible), and document or rename to make the consume-style mutation explicit.

### [ANAL-26] IsobaricChannelExtractor::extractChannels — extractChannels silently rewrites the stored 'select_activation' member, making repeat calls non-idempotent
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h · _quant-isobaric_

```cpp
void extractChannels(const PeakMap& ms_exp_data, ConsensusMap& consensus_map)
```
- **Expectation:** extractChannels reads MS/MS data and writes channels into the consensus map. A caller would not expect it to alter the object's own configuration state.
- **Actual:** When select_activation is 'auto', extractChannels overwrites the member selected_activation_ with the concrete 'HCID,HCD' string. This member now diverges from the Param value ('auto' still stored in getParameters()), and on a second call the 'auto' branch is no longer taken (the member is already the expanded list), so the behavior of the second invocation depends on whether the first ran. It also accumulates QC into channel_mz_delta and calls printStats() (log spam) without that being implied by the name.
- **Evidence:** IsobaricChannelExtractor.cpp:429-432 `if (selected_activation_ == "auto") { selected_activation_ = ...HCID + "," + ...HCD; }` mutating the member, while updateMembers_ only ever sets it from the param. printStats() called at :702.
- **Fix:** Expand 'auto' into a local variable instead of overwriting the member, so the object stays in sync with its Param and calls are idempotent. Source- and ABI-compatible (internal change only).
- **Verifier correction:** extractChannels has two undocumented side effects on the object: (1) when select_activation=="auto" it overwrites the member selected_activation_ with "HCID,HCD" (cpp:431), leaving the member out of sync with getParameters() (still "auto") and skipping the "auto" branch on later calls; (2) it accumulates QC into the member channel_mz_delta (cpp:645-646) without resetting it (the public clearStats() is never invoked by extractChannels), and prints stats via printStats() (cpp:702). Correction to the original claim: the *channel-extraction result written into consensus_map is idempotent* — "auto" and the expanded "HCID,HCD" select the same scans, so the consensus output is identical on repeat calls. The non-idempotency that actually bites is (a) config-state divergence (getParameters() vs internal member) and (b) silently accumulating/doubling QC statistics exposed via getStats() on object reuse. Severity is medium (not high): core quantitation intensities are not silently wrong, and there is log output; but QC numbers and reported config silently drift on reuse and are easy to miss. Fix is internal-only and ABI/source compatible: expand "auto" into a local string instead of the member, and reset channel_mz_delta at the start of extractChannels (or document caller's clearStats() responsibility).

### [ANAL-27] IsobaricIsotopeCorrector::correctIsotopicImpurities — Static correction requires out-map pre-sized to in-map; the only guard is a debug-only precondition
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h · _quant-isobaric_

```cpp
static IsobaricQuantifierStatistics correctIsotopicImpurities(const ConsensusMap& consensus_map_in, ConsensusMap& consensus_map_out, const IsobaricQuantitationMethod* quant_method)
```
- **Expectation:** Given the doc '@param[in,out] consensus_map_out The map where the corrected values should be stored', a caller could reasonably pass an empty or default-constructed output map and expect it to be filled, or expect a thrown exception if the sizes mismatch.
- **Actual:** The function indexes consensus_map_out[i] for every i in [0, consensus_map_in.size()), clearing and re-inserting handles, but never resizes consensus_map_out. The size match is enforced only by OPENMS_PRECONDITION, which is a no-op unless OPENMS_ASSERTIONS is defined. In a release build, passing an out-map smaller than in (e.g. empty) is out-of-bounds access / undefined behavior with no error. (The IsobaricQuantifier::quantify path is safe because it does consensus_map_out = consensus_map_in first, but the static method is public and callable directly.)
- **Evidence:** IsobaricIsotopeCorrector.cpp:44-45 `OPENMS_PRECONDITION(consensus_map_in.size() == consensus_map_out.size(), ...)`; loop :81 `for (... i < consensus_map_out.size() ...)` with :87 `consensus_map_out[i].clear()` and :90 reads `consensus_map_in[i]`. Macros.h:87-94 makes OPENMS_PRECONDITION empty without OPENMS_ASSERTIONS.
- **Fix:** Replace the debug-only precondition with an always-on size check that throws Exception::InvalidParameter (or copy/resize consensus_map_out from consensus_map_in internally as quantify() does). Document that consensus_map_out must already have the same size as consensus_map_in. The throw is source/ABI compatible.
- **Verifier correction:** The only sound caller (IsobaricQuantifier::quantify) always sets consensus_map_out = consensus_map_in before calling, so the precondition is satisfied in-tree. For the public static method called directly: passing an EMPTY out-map is a SILENT no-op (loop bounded by consensus_map_out.size()==0; returns zeroed stats, no exception, no correction) — not out-of-bounds UB as the evidence states. Out-of-bounds/UB occurs only when consensus_map_out is non-empty and LARGER than consensus_map_in, since cpp L90 reads consensus_map_in[i] via unchecked noexcept operator[]. The size guard is OPENMS_PRECONDITION, which is empty in release (Macros.h L94). Severity is medium (not high): the loud OOB crash needs an unusual larger-out-map, while the natural empty-out misuse fails silently but recoverably (no data corruption, just no work done) and is fully avoidable by the documented/internal usage; the doc omission plus debug-only guard invites the misuse. Fix (throw Exception::InvalidParameter on size mismatch, or resize/copy out from in internally as quantify does) is source- and ABI-compatible.

### [ANAL-30] ItraqConstants::getIsotopeMatrixAsStringList / updateIsotopeMatrixFromStringList — Doc says the IsotopeMatrices vector holds 'the two matrices (4plex, 8plex)', but the code requires three and indexes by itraq_type
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h · _quant-tmt-itraq_

```cpp
static StringList getIsotopeMatrixAsStringList(const int itraq_type, const IsotopeMatrices & isotope_corrections)
```
- **Expectation:** The @param doc on both getIsotopeMatrixAsStringList and updateIsotopeMatrixFromStringList states 'isotope_corrections Vector of the two matrices (4plex, 8plex).', so a caller would size/build the vector with two entries.
- **Actual:** The implementation indexes isotope_corrections[itraq_type] where itraq_type ranges over FOURPLEX, EIGHTPLEX, TMT_SIXPLEX (SIZE_OF_ITRAQ_TYPES == 3), and updateIsotopeMatrixFromStringList even resize()s to SIZE_OF_ITRAQ_TYPES and fills all three (including ISOTOPECORRECTIONS_TMT_SIXPLEX). A caller who follows the 'two matrices' doc and passes a 2-element vector triggers out-of-bounds access for itraq_type==TMT_SIXPLEX.
- **Evidence:** Header doc: '@param[in] isotope_corrections Vector of the two matrices (4plex, 8plex).' vs .cpp: 'isotope_corrections.resize(SIZE_OF_ITRAQ_TYPES); isotope_corrections[2].setMatrix<double, 6, 4>(ItraqConstants::ISOTOPECORRECTIONS_TMT_SIXPLEX);' and access 'isotope_corrections[itraq_type](i, j)'.
- **Fix:** Fix the @param documentation to say the vector is indexed by ITRAQ_TYPES and must contain SIZE_OF_ITRAQ_TYPES (currently 3: 4plex, 8plex, TMT-6plex) entries. Doc-only change; no ABI impact.
- **Verifier correction:** Severity is medium rather than high: the failure mode is a loud crash (out-of-bounds on operator[]) for a caller who manually sizes a 2-element vector and then calls with TMT_SIXPLEX, not silent wrong results. The dominant workflow builds the vector via updateIsotopeMatrixFromStringList, which always resizes to SIZE_OF_ITRAQ_TYPES (3), so it self-corrects the size; the OOB path requires partially ignoring the itraq_type doc. The fix is doc-only (correct both @param lines in ItraqConstants.h:76,90 and the pyOpenMS docstrings in bind_analysis.cpp:720,732 to state the vector is indexed by ITRAQ_TYPES and must hold SIZE_OF_ITRAQ_TYPES==3 entries: 4plex, 8plex, TMT-6plex). No ABI impact.

### [ANAL-57] SpectralDeconvolution::setTargetDecoyType — Doc enumerates four decoy types (charge/noise/isotope dummy) that do not match the actual 3-value TargetDecoyType enum
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h · _topdown-deconv_

```cpp
void setTargetDecoyType(PeakGroup::TargetDecoyType target_decoy_type, const DeconvolvedSpectrum& target_dspec_for_decoy_calcualtion)
```
- **Expectation:** The Doxygen describing the parameter should match the actual enum the caller must pass. A caller reads the parameter doc to know which integer/value means what.
- **Actual:** The param doc claims values are 'target (0), charge dummy (1), noise dummy (2), or isotope dummy (3)' - four entries with different names and numbers. But PeakGroup::TargetDecoyType only defines target=0, noise_decoy=1, signal_decoy=2 (three values). A caller passing signal_decoy (value 2) believing it means 'noise dummy' per this doc gets the wrong decoy class.
- **Evidence:** SpectralDeconvolution.h:121 `... a target (0), charge dummy (1), noise dummy (2), or isotope dummy (3)`. PeakGroup.h:36-41 `enum TargetDecoyType { target = 0, noise_decoy, signal_decoy, };`. PeakGroup.h:34-35 class comment instead says `target (0), noise decoy (1), or signal decoy (2)`.
- **Fix:** Fix the SpectralDeconvolution::setTargetDecoyType Doxygen to match the actual 3-value enum (target/noise_decoy/signal_decoy) and reconcile the two conflicting class comments. Doc-only change, no ABI impact.
- **Verifier correction:** Severity is medium, not high: the parameter is strongly typed as the PeakGroup::TargetDecoyType enum, so realistic callers pass named enumerators (PeakGroup::noise_decoy) rather than raw integers, which guards against the worst silent-misuse. The harm is doc-induced confusion (wrong name/number mapping, inverted noise vs. signal) for a developer reading the param doc. Recommendation stands: fix the SpectralDeconvolution::setTargetDecoyType Doxygen to match the 3-value enum and reconcile the conflicting class/param comments. Also note the @param[out] tags are mislabeled (both are inputs), a secondary doc defect. Doc-only, no ABI impact.

### [ANAL-59] PeakGroup::getMzRange — getMzRange returns an inverted sentinel {-1, -10} when the charge has no peaks instead of signaling absence
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h · _topdown-deconv_

```cpp
std::tuple<double, double> getMzRange(int abs_charge) const
```
- **Expectation:** A getter returning an m/z range for a charge state should either return a valid [low, high] range or clearly indicate 'no range'. A caller naively unpacking {start,end} expects start <= end.
- **Actual:** When abs_charge is out of [min,max] range or no logMzPeak has that charge, the method returns the initializer values {mz_start=-1, mz_end=-10}, i.e. an inverted negative range, with no other signal. getRepMzRange() forwards this directly, so the same sentinel leaks out.
- **Evidence:** PeakGroup.cpp getMzRange: `double mz_start = -1; double mz_end = -10; ... return std::tuple<double,double>{mz_start, mz_end};` (the loop is skipped entirely when the charge guard fails or no peak matches).
- **Fix:** Document the {-1,-10} sentinel in the header, or return std::optional / a boolean-validity flag. Doc-only is non-breaking; changing the return type is ABI-breaking, so prefer an additional `bool hasMzRange(int)` helper plus a doc note.
- **Verifier correction:** getMzRange returns an inverted, negative sentinel {-1, -10} (start=-1 > end=-10) when abs_charge is out of [min,max] or no peak matches that charge, and getRepMzRange forwards it. The sentinel is undocumented and is reachable via the default code path (getRepAbsCharge() returns the default max_snr_abs_charge_ == -1, which fails the charge guard). It leaks unguarded into the FLASHDeconv StartMz/EndMz TSV columns (FLASHDeconvSpectrumFile.cpp:198). Severity is medium rather than high: the bad output is loud and recoverable (obviously-invalid negative/inverted m/z), not a silent numerical corruption, and the common production path usually supplies a valid rep charge. Minimal fix (document the sentinel in the header and/or add an additive bool hasMzRange(int) helper) is ABI-neutral; only changing the return type to std::optional would be ABI-breaking.

### [ANAL-66] FLASHHelperClasses::PrecalculatedAveragine::getMaxIsotopeIndex — getMaxIsotopeIndex returns an uninitialized member when setMaxIsotopeIndex was never called
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h · _topdown-misc_

```cpp
size_t getMaxIsotopeIndex() const
```
- **Expectation:** A getter on a fully-constructed object returns a well-defined value (e.g. 0 or a constructor-computed default). After building a PrecalculatedAveragine via its parameterized constructor, getMaxIsotopeIndex() should be safe to read.
- **Actual:** max_isotope_index_ is NOT initialized by either constructor (no default initializer, not in the ctor init list). Reading it before setMaxIsotopeIndex returns indeterminate garbage. MassFeatureTrace.cpp:111 uses it to size a vector: 'std::vector<float>(averagine.getMaxIsotopeIndex(), .0f)', so an unset cache yields a huge/garbage allocation.
- **Evidence:** Member decl: 'Size max_isotope_index_;' (no initializer, FLASHHelperClasses.h:75). Ctor init list (FLASHHelperClasses.cpp:15-18) sets only mass_interval_ and min_mass_. Header @warning confirms: 'The constructor does not initialise max_isotope_index_; reading this before setMaxIsotopeIndex ... returns an uninitialised value.'
- **Fix:** Initialize the member: 'Size max_isotope_index_ = 0;' (additive, no ABI break for an in-class default initializer on a non-virtual data member with unchanged layout). This turns undefined behavior into a defined (if conservative) default and removes the documented landmine.
- **Verifier correction:** The member is genuinely uninitialized by both constructors and reading it before setMaxIsotopeIndex is UB, but the claim overstates severity as "high". In all in-tree production code the value is set: SpectralDeconvolution::calculateAveragine() (SpectralDeconvolution.cpp:267-269) constructs the averagine and calls setMaxIsotopeIndex on the very next line, and the only consumer (MassFeatureTrace.cpp:111) reads an averagine obtained from getAveragine() after calculateAveragine ran. The defect is therefore reachable only through documented-against misuse (default-construct or parameterized-construct then read before set, or direct pyOpenMS construction), which together with the explicit @warning at the getter makes this medium severity, not high. Fix `Size max_isotope_index_ = 0;` is correct and ABI/layout-safe.

### [ANAL-10] MapAlignmentAlgorithmIdentification::align — align() silently clears a previously set external reference when an internal reference_index is given, but quietly keeps it otherwise
`low` · `documentation-gap` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h · _mm-alignment_

```cpp
template <typename DataType> void align(const std::vector<DataType>& data, std::vector<TransformationDescription>& transformations, Int reference_index = -1)
```
- **Expectation:** Having called setReference(externalData), a caller expects align(data, trafos, idx) to behave consistently with respect to that external reference regardless of whether idx is supplied.
- **Actual:** align() conditionally discards the previously configured reference_: if reference_index >= 0 it does `reference_.clear()` and repopulates from data[reference_index]; if reference_index < 0 it silently keeps whatever setReference() left. The reference source thus depends on the value of a defaulted integer argument, and an explicitly set external reference is silently dropped when an index is passed. This dual behavior is non-obvious from the signature.
- **Evidence:** Lines 98-118: `bool use_internal_reference = (reference_index >= 0); ... if (use_internal_reference) reference_.clear(); ... if (use_internal_reference) { ... setReference(data[reference_index]); }` with the comment 'External references set explicitly via setReference() are preserved (reference_index < 0).'
- **Fix:** Document the precedence rule on the align() declaration (passing reference_index >= 0 overrides and clears any external reference set via setReference). Doc-only; ABI-safe.

### [ANAL-2] FeatureGroupingAlgorithmUnlabeled::getResultMap — Non-const getResultMap() returns a mutable reference into the algorithm's internal pairfinder state
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h · _mm-grouping_

```cpp
ConsensusMap & getResultMap()
```
- **Expectation:** A getter named getResultMap returning a 'result' would be expected to be const (read-only access to the computed result) or to hand back ownership of a finished object.
- **Actual:** getResultMap() is non-const and returns 'pairfinder_input_[0]' by mutable reference — i.e. the live internal working buffer that addToGroup continues to swap into. The caller can mutate algorithm internals through the returned reference, and a subsequent addToGroup() will overwrite what the caller is holding. There is no const overload, so even read-only callers cannot use it on a const object.
- **Evidence:** Header: 'ConsensusMap & getResultMap() { return pairfinder_input_[0]; }'. The private member comment: 'the first element is the currently computed consensus map... after adding all maps through addToGroup it will contain the final result'.
- **Fix:** Add a const overload returning 'const ConsensusMap&' (additive, source-compatible). Document that the returned reference aliases internal state invalidated by further addToGroup calls. Removing the non-const overload would be breaking, so keep it.
- **Verifier correction:** The core verifiable surprise is the missing const overload plus an undocumented aliasing getter: getResultMap() is non-const-only and returns a mutable reference into the live internal pairfinder_input_[0] buffer, so (1) it cannot be called on a const object for read-only access, and (2) the getter's doc comment never warns the reference aliases internal state. The claim's "a subsequent addToGroup overwrites what the caller is holding" is true but only under out-of-contract use (documented order is setReference -> addToGroup -> getResultMap), and the sole in-tree caller copies the result immediately, so this is a mild const-correctness/API-hygiene surprise (low), not a high-severity silent-corruption hazard. Fix: add `const ConsensusMap& getResultMap() const`, document the aliasing/invalidation, keep the non-const overload (additive, source-compatible).

### [ANAL-16] ConsensusMapNormalizerAlgorithmThreshold::computeCorrelation — computeCorrelation does not compute any correlation; it returns mean intensity ratios (normalization factors)
`low` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h · _mm-normalizer_

```cpp
static std::vector<double> computeCorrelation(const ConsensusMap & map, const double & ratio_threshold, const std::string & acc_filter, const std::string & desc_filter)
```
- **Expectation:** A method named computeCorrelation returns correlation coefficient(s) (e.g. Pearson/Spearman in [-1,1]) between maps.
- **Actual:** It computes, per map j, the mean of per-feature intensity ratios reference_intensity/intensity_j (after outlier rejection) and returns these as per-map scaling factors. No correlation is ever computed. The class doc and the .cpp doc literally call it 'Compute one normalisation factor per map.' and the pyOpenMS docstring says 'Determines the ratio of all maps to the map with the most features'.
- **Evidence:** Header doc: '@brief Compute one normalisation factor per map.' Body: 'double ratio = feature_int[map_with_most_features_idx][k] / feature_int[j][k]; ... ratio_vector[j] = Math::mean(ratios.begin(), ratios.end());'. pyOpenMS bind_analysis.cpp: "Determines the ratio of all maps to the map with the most features".
- **Fix:** Add a correctly-named additive alias such as computeNormalizationFactors(...) delegating to the existing method, and mark computeCorrelation [[deprecated("renamed to computeNormalizationFactors; it returns mean intensity ratios, not a correlation")]]. Keeping the old symbol preserves ABI; ideal removal is breaking.
- **Verifier correction:** The misleading-name finding is real and accurately described, but its severity is low rather than high. The function is fully and correctly documented at its declaration (header @brief 'Compute one normalisation factor per map.' and a @return clause stating it returns one normalization factor per map), and it produces correct normalization output with no silently-wrong results, crashes, or data loss. The surprise is confined to the symbol name contradicting its documented/actual semantics (returns mean intensity ratios = per-map scaling factors, not a correlation coefficient). Recommended fix (additive, source-compatible): add computeNormalizationFactors(...) delegating to the existing body and mark computeCorrelation [[deprecated]]; this preserves ABI. Outright rename/removal would be breaking.

### [ANAL-19] ConsensusMapNormalizerAlgorithmMedian::passesFilters_ — Public API method carries the trailing-underscore private-helper naming convention
`low` · `inconsistent-convention` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h · _mm-normalizer_

```cpp
static bool passesFilters_(ConsensusMap::ConstIterator cf_it, const ConsensusMap& map, const std::string& acc_filter, const std::string& desc_filter)
```
- **Expectation:** In OpenMS, a trailing underscore (foo_) marks a private/internal member; a competent caller would not expect to call it, and would not expect it in the public section.
- **Actual:** passesFilters_ is declared in the public: section, is part of OPENMS_DLLAPI, and is intentionally called cross-class by ConsensusMapNormalizerAlgorithmThreshold (and could be by any external code). Its name signals 'do not call me' while its access signals the opposite.
- **Evidence:** Declared under 'public:' as 'static bool passesFilters_(...)'; called externally in ConsensusMapNormalizerAlgorithmThreshold.cpp: 'ConsensusMapNormalizerAlgorithmMedian::passesFilters_(cf_it, map, acc_filter, desc_filter)'. Class comment: 'our friend can use our passesWhitelist_() method'.
- **Fix:** Provide a public, conventionally-named alias (e.g. passesFilters) that forwards to passesFilters_, and keep passesFilters_ for ABI. Renaming outright would break the exported symbol, so treat alias as additive; full rename is breaking.

### [ANAL-20] ConsensusMapNormalizerAlgorithmMedian::passesFilters_ — Empty filter string means 'match everything' (always pass), inverting the natural 'filter rejects nothing configured' expectation only after reading code
`low` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h · _mm-normalizer_

```cpp
static bool passesFilters_(ConsensusMap::ConstIterator cf_it, const ConsensusMap& map, const std::string& acc_filter, const std::string& desc_filter)
```
- **Expectation:** Header doc says the method 'returns whether consensus feature passes accession regexp and description regexp'; a caller passing an empty regex might expect it to match the empty string only, or behave as a real regex.
- **Actual:** An empty acc_filter or desc_filter short-circuits to 'pass everything' (and also a regex that matches the empty string passes everything, even features with no identifications). The empty-string and 'regex matches ""' special cases are not visible in the signature; the header for the Median class does not state the empty-string semantics (only the Threshold header documents 'empty means match anything').
- **Evidence:** 'if ((acc_filter.empty() || boost::regex_search("", m, acc_regexp)) && (desc_filter.empty() || boost::regex_search("", m, desc_regexp))) { // feature passes (even if it has no identification!) return true; }'
- **Fix:** Document the 'empty == match anything' and 'regex matching the empty string == match anything' semantics in the Median header doc (already documented in the Threshold header). Doc-only; ABI-safe.
- **Verifier correction:** The empty-filter "match anything" behavior is a documented convention in the sibling Threshold header and a common idiom, so it is not itself the surprise. The real, entirely-undocumented surprise is that a non-empty regex which matches the empty string (e.g. `.*`, `^`, `a*`) takes the same short-circuit and passes ALL features — including features with no identifications at all — via `boost::regex_search("", m, acc_regexp)`. This silently defeats accession/description filtering for such patterns. Fix is doc-only: document both the empty-string and the regex-matches-empty-string semantics in the Median header (and ideally the regex-matches-empty case in the Threshold header too). ABI-safe.

### [ANAL-23] ConsensusMapNormalizerAlgorithmMedian (constructors/destructor doc comments) — Doc comments claim ctor/dtor/copy 'is not implemented -> private' but they are public and implemented
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h · _mm-normalizer_

```cpp
ConsensusMapNormalizerAlgorithmMedian(); virtual ~ConsensusMapNormalizerAlgorithmMedian();
```
- **Expectation:** A reader trusts the Doxygen comment; the comments state the default constructor and destructor are 'not implemented -> private'.
- **Actual:** The default constructor and destructor are declared public and ARE implemented (defaulted in the .cpp). Only the copy ctor and assignment are actually private. The comments are copy-paste errors that directly contradict the access/implementation a caller relies on.
- **Evidence:** Header (public section): '/// default constructor is not implemented -> private\n ConsensusMapNormalizerAlgorithmMedian();' and '/// destructor is not implemented -> private\n virtual ~...();'. .cpp: 'ConsensusMapNormalizerAlgorithmMedian::ConsensusMapNormalizerAlgorithmMedian() = default;' and '... ~...() = default;'.
- **Fix:** Fix the misleading Doxygen comments to describe the actual public/defaulted ctor and dtor. Comment-only change; ABI-safe.
- **Verifier correction:** The doc comments on the public default constructor (header line 34) and destructor (header line 37) falsely state they are "not implemented -> private". In fact both are declared public (header line 33 `public:`) and implemented as `= default;` in the .cpp (lines 22 and 24). The "not implemented -> private" pattern is accurate only for the truly-private copy ctor and assignment operator (header lines 27-31). Recommended fix is to correct the two stray comments; comment-only, ABI-safe. Severity is low (harmless confusion, not a correctness/safety hazard).

### [ANAL-5] FeatureDistance::FeatureDistance — Class documentation refers to a constructor parameter 'check_constraints' that does not exist (the parameter is 'force_constraints')
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h · _mm-pairfinder_

```cpp
FeatureDistance(double max_intensity = 1.0, bool force_constraints = false)
```
- **Expectation:** The constructor's bare bool argument is unreadable at call sites (FeatureDistance(maxInt, true)), so a caller relies on the class doc to know which flag they are passing. The doc names the flag they should look for.
- **Actual:** The class-level doc says the behavior 'depends on the value used for @p check_constraints in the constructor' (FeatureDistance.h:34), but the actual constructor parameter is named force_constraints (FeatureDistance.h:77) and stored as force_constraints_. There is no check_constraints anywhere. A reader grepping for check_constraints finds nothing; the bare bool at call sites already obscures intent.
- **Evidence:** FeatureDistance.h:34: 'the behavior depends on the value used for @p check_constraints in the constructor'. FeatureDistance.h:76-77: 'FeatureDistance(double max_intensity = 1.0, bool force_constraints = false);'. The @param doc at line 74 itself calls it force_constraints.
- **Fix:** Fix the class doc to say force_constraints (the @param block already uses that name). Doc-only fix is ABI-neutral and source-compatible; do not rename the parameter to avoid churn.
- **Verifier correction:** The claim is accurate but over-stated on severity. The class doc at FeatureDistance.h:34 references @p check_constraints, which does not exist; the constructor parameter is force_constraints (line 77 / @param line 74) and is stored as force_constraints_ (FeatureDistance.cpp). The only repo occurrence of check_constraints is this stale doc line. This is a low-severity, doc-only inconsistency (the correct name force_constraints already appears in the same comment's @param block), not a medium-severity misuse hazard. Fix: rename @p check_constraints to @p force_constraints on line 34. ABI: none; source-compatible.

### [ANAL-6] FeatureMapping::assignMS2IndexToFeature — Two adjacent same-typed double tolerances (each a half-width) plus a trailing bare bool that governs only the m/z unit are swap/ambiguity hazards at call sites
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h · _mm-pairfinder_

```cpp
static FeatureToMs2Indices assignMS2IndexToFeature(const MSExperiment& spectra, const FeatureMappingInfo& fm_info, const double& precursor_mz_tolerance, const double& precursor_rt_tolerance, bool ppm)
```
- **Expectation:** From the signature alone a caller cannot tell that the two doubles are half-widths (window is +/- tol), that 'ppm' switches only the m/z tolerance between ppm and Th while RT is always seconds, or that the order must be mz-then-rt. Adjacent const double& params plus a trailing unlabeled bool are easy to swap/misread.
- **Actual:** The implementation builds the m/z window via Math::getTolWindow(mz, precursor_mz_tolerance, ppm) and the RT window as [rt - precursor_rt_tolerance, rt + precursor_rt_tolerance] (FeatureMapping.cpp:43-44), confirming both are half-widths and ppm applies only to m/z. The header doc explains this, but the signature gives no protection against swapping the two doubles.
- **Evidence:** FeatureMapping.cpp:43-44: 'Math::getTolWindow(mz, precursor_mz_tolerance, ppm); fm_info.kd_tree.queryRegion(rt - precursor_rt_tolerance, rt + precursor_rt_tolerance, ...)'. Header doc FeatureMapping.h:99-101 confirms half-width and ppm-vs-Th semantics.
- **Fix:** Documentation is thorough; for safety add (additively) an overload taking a small tolerance struct (mz, rt, ppm) so call sites self-describe, keeping the existing signature for ABI stability. Pure addition is source- and ABI-compatible.
- **Verifier correction:** The surprise is real at the signature level but should be graded low, not medium/high: the half-width, RT-seconds, and ppm-only-applies-to-m/z semantics are all documented explicitly at the declaration (FeatureMapping.h:99-101), and the sole caller passes the arguments through named getters, so the swap/misread hazard is mitigated in practice. The recommended additive tolerance-struct overload is source- and ABI-compatible; abi_impact of the finding itself is none since no signature is changed.

### [ANAL-12] TransformationDescription::estimateWindow — const estimateWindow() with invert=true silently deep-copies and refits the entire transformation model on every call
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h · _mm-transform_

```cpp
double estimateWindow(double quantile = 0.99, bool invert = true, bool full_window = true, double padding_factor = 1.0) const
```
- **Expectation:** A const method named estimateWindow() that takes a quantile reads the stored residuals and computes a window cheaply, with no model re-fitting.
- **Actual:** With the default invert=true, the implementation constructs a full copy of the TransformationDescription (TransformationDescription tmp(*this);) and calls tmp.invert(), which swaps all data points AND re-runs fitModel(), re-estimating the whole (possibly b-spline/lowess) model just to measure residuals in x-units. This is an expensive hidden side effect triggered by a getter-style const call.
- **Evidence:** TransformationDescription.cpp: "TransformationDescription tmp(*this); if (invert) { ... tmp.invert(); }" and invert() calls "fitModel(model_type_, params);". Header doc only says it returns a window; nothing warns the model is refit.
- **Fix:** Document the refit cost in the header. ABI-safe: keep signature, add a @note that invert=true triggers a full copy+refit; consider an overload that accepts a pre-inverted description or a cached inverse to avoid per-call refits in hot paths (e.g. SwathMapMassCorrection).
- **Verifier correction:** const estimateWindow() always deep-copies the TransformationDescription and refits the model via the copy constructor (unconditionally, even with invert=false); with the default invert=true it refits a SECOND time inside invert(). Results are correct, so this is a hidden-cost/performance surprise, not a wrong-result bug. Severity is low (mild surprise, recoverable, callers are bounded calibration/setup paths, not inner loops), not high. Fix is documentation-only plus an optional overload accepting a pre-inverted/cached description; signature unchanged so ABI impact is none.

### [ANAL-13] TransformationDescription::estimateWindow — estimateWindow returns 0.0 (a valid window) when there are no data points, indistinguishable from a genuinely zero window
`low` · `api-contract` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h · _mm-transform_

```cpp
double estimateWindow(...) const
```
- **Expectation:** When asked to estimate an extraction window from data that does not exist, a caller would expect an error/throw or an obviously-invalid sentinel, not a plausible numeric window.
- **Actual:** With no/empty residuals the method returns 0.0, which a downstream consumer will treat as a real (degenerate) extraction half/full-width. A zero extraction window silently disables matching rather than flagging missing calibration data.
- **Evidence:** TransformationDescription.cpp: "if (diffs.empty()) { OPENMS_LOG_DEBUG ... return 0.0; }" Header @return: "If no data points are available, returns 0.0."
- **Fix:** ABI-safe: at minimum log at WARNING not DEBUG. Better: document that callers must guard against 0.0, or add a bool* success out-param overload / throw on empty so the 'no data' case is not confused with a real zero window.

### [ANAL-15] TransformationDescription::invert — invert() refits the model and (for empty-data linear) mutates the model in place via a downcast 'ugly hack'
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h · _mm-transform_

```cpp
void invert()
```
- **Expectation:** invert() inverts the transformation, leaving the stored data points conceptually consistent and the model coherent.
- **Actual:** invert() swaps first/second of every stored DataPoint AND then re-fits the model from scratch via fitModel(); for the linear-with-explicit-params/empty-data case it instead dynamic_casts the owned model to TransformationModelLinear and calls lm->invert(). So the same call has two structurally different code paths, mutates the data points, and re-runs a potentially expensive non-linear fit. The header only says 'Computes an (approximate) inverse'.
- **Evidence:** TransformationDescription.cpp invert(): swaps every DataPoint, then "// ugly hack for linear model with explicit slope/intercept parameters: if ((model_type_ == \"linear\") && data_.empty()) { TransformationModelLinear* lm = dynamic_cast<...>(model_.get()); lm->invert(); } else { ... fitModel(model_type_, params); }".
- **Fix:** ABI-safe: document that invert() rewrites stored data points and triggers a refit (cost + that the inverse of non-linear models is only approximate, as the doc hints). No signature change required.
- **Verifier correction:** Claim is factually accurate in full. Severity adjusted from the implied higher tier to low: invert() does mutate stored data points in place and trigger a refit (potentially an expensive b_spline/lowess/interpolated fit), and uses a self-described "ugly hack" dynamic_cast path for the linear+empty-data case — but the most important consequence, that non-linear inverses are only approximate, is already signposted by the word "(approximate)" in the header doc. The behavior produces correct (if approximate for non-linear) results, no data loss/crash for reasonable use; the surprises are an undocumented in-place data rewrite plus refit cost. Fix is documentation-only (note the data-point rewrite and refit cost); no signature change, ABI impact none.

### [ANAL-44] NuXLFDR::calculatePeptideAndXLQValueAndFilterAtPSMLevel — A 0.0 q-value threshold means '100% FDR / filter disabled', the opposite of the most stringent cut a caller would assume
`low` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFDR.h · _nuxl-core_

```cpp
..., std::vector<double> xl_PSM_qvalue_thresholds, std::vector<double> xl_peptidelevel_qvalue_thresholds, const std::string& out_idxml, int decoy_factor
```
- **Expectation:** For a q-value/FDR threshold parameter, 0.0 is the most stringent possible cut (accept nothing above 0% FDR); a caller passing 0.0 expects maximal stringency.
- **Actual:** `std::replace(xl_PSM_qvalue_thresholds.begin(), xl_PSM_qvalue_thresholds.end(), 0.0, 1.0);` rewrites every 0.0 entry to 1.0 (NuXLFDR.cpp line 228), i.e. treats 0.0 as 'filter disabled = 100% FDR'. The exact value that looks most stringent is silently converted to the least stringent. The passed vector is also taken by value and sorted descending, so the caller's ordering is not honored.
- **Evidence:** NuXLFDR.cpp: `// treat disabled filtering as 100% FDR` `std::replace(xl_PSM_qvalue_thresholds.begin(), xl_PSM_qvalue_thresholds.end(), 0.0, 1.0);`
- **Fix:** Use a distinct sentinel (e.g. negative or -1) for 'disabled' instead of overloading 0.0, or document this inversion at the parameter (the header does note it, but the convention itself is the surprise). ABI: behavior/contract clarification, source-compatible.
- **Verifier correction:** The 0.0->1.0 rewrite and descending in-place sort exist and are accurately quoted, but they are a documented, internally-consistent convention (every threshold in the method treats values outside the open interval (0,1) as "filter disabled"), matching the user-facing report:xlFDR semantics where 1.0 = no filtering (OpenNuXL.cpp:5041 injects 1.0 by default). The header documents the rewrite and the sort at the exact @param (lines 144-145, 166). The remaining issue is only the mild API-taste choice of overloading 0.0 rather than a distinct sentinel (e.g. -1); this is a low-severity documentation/clarity nit, not a silently-wrong-results bug, and the header already notes the behavior.

### [ANAL-46] NuXLDeisotoper::deisotopeAndSingleCharge — preserve_low_mz_peaks_threshold acts independently of preserve_high_intensity_peaks, contradicting its documented coupling
`low` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLDeisotoper.h · _nuxl-core_

```cpp
..., bool preserve_high_intensity_peaks = false, double preserve_low_mz_peaks_threshold = -1e10
```
- **Expectation:** The header doc says 'If preserve_high_intensity_peaks is set, all peaks with smaller m/z will never be filtered', so a caller expects preserve_low_mz_peaks_threshold to take effect only when preserve_high_intensity_peaks is true.
- **Actual:** The low-m/z preservation block runs unconditionally whenever `preserve_low_mz_peaks_threshold > 0.0`, with no check of preserve_high_intensity_peaks (NuXLDeisotoper.cpp lines 385-401). The two flags are independent in the code; the documented dependency is false, so a caller who sets only preserve_low_mz_peaks_threshold gets behavior the doc says requires the other flag.
- **Evidence:** NuXLDeisotoper.cpp: `if (preserve_low_mz_peaks_threshold > 0.0) { for (...) { if (spec[i].getMZ() < preserve_low_mz_peaks_threshold) {...} else { break; } } }` -- no reference to preserve_high_intensity_peaks.
- **Fix:** Fix the header doc to describe preserve_low_mz_peaks_threshold as an independent threshold; the default -1e10 (rather than 0) also obscures that the feature is off by default. ABI: doc-only, non-breaking.
- **Verifier correction:** The parameter name is fine; the header doc is what's misleading. NuXLDeisotoper.cpp:385 gates the low-m/z preservation solely on `preserve_low_mz_peaks_threshold > 0.0`, independent of preserve_high_intensity_peaks, whereas the header doc (line 71) claims the threshold only acts when preserve_high_intensity_peaks is set. Fix is doc-only (non-breaking): document preserve_low_mz_peaks_threshold as an independent threshold (peaks below it are always retained when it is > 0.0; disabled by the -1e10 default).

### [ANAL-47] NuXLReport::annotate — annotate() returns report rows but also mutates the passed-in peptide_ids by writing many meta values
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLReport.h · _nuxl-core_

```cpp
static std::vector<NuXLReportRow> annotate(const PeakMap& spectra, PeptideIdentificationList& peptide_ids, const StringList& meta_values_to_export, double marker_ions_tolerance)
```
- **Expectation:** A function with a return type std::vector<NuXLReportRow> reads as a pure builder that produces rows from inputs; the non-const PeptideIdentificationList& looks like it could be read-only-by-omission-of-const, and the header gives no documentation at all.
- **Actual:** Besides returning rows, annotate writes numerous meta values back into each PeptideHit in peptide_ids: NuXL:peptide_mass_z0, NuXL:xl_mass_z0, marker-ion intensities, NuXL:Da difference, PRECURSOR_ERROR_PPM_USERPARAM, NuXL:z1..z4 mass (NuXLReport.cpp lines 232-265). The mutation of the input is a substantial, undocumented side effect for a method whose name and return type suggest report generation.
- **Evidence:** NuXLReport.cpp: `ph.setMetaValue("NuXL:peptide_mass_z0", DataValue(peptide_weight)); ph.setMetaValue("NuXL:xl_mass_z0", xl_weight); ... ph.setMetaValue("NuXL:Da difference", (double)absolute_difference); ph.setMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM, ...); ph.setMetaValue("NuXL:z1 mass", ...);`
- **Fix:** Document in the header that peptide_ids is annotated in place (the by-value/by-ref distinction and the meta keys written), or split the meta-annotation into an explicitly named step. ABI: doc-only, non-breaking.
- **Verifier correction:** annotate() does mutate the passed-in peptide_ids in place (adding ~10 NuXL: meta values plus PRECURSOR_ERROR_PPM_USERPARAM per hit), as claimed. But the verb `annotate` is the established convention for in-place enrichment via a non-const reference (cf. the sibling annotateProteinModificationForTopHits in the same header), and the side effect is intentional/load-bearing: the only caller (OpenNuXL.cpp) reads the written ppm-error back at line 6085 and persists the enriched peptide_ids to idXML at 6338. The real, minor issue is the undocumented dual responsibility (returning report rows AND annotating hits) — a one-line header comment noting peptide_ids is annotated in place and listing the meta keys would resolve it. Severity is low (no silent wrong results, data loss, or crash; mutations are additive and surface in output), not medium.

### [ANAL-50] NuXLFragmentAdductDefinition::operator< / operator== — operator< compares by (mass, formula, name) but operator== ignores mass — set/unordered-set dedup disagree
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFragmentAdductDefinition.h · _nuxl-frag_

```cpp
bool operator<(const NuXLFragmentAdductDefinition& other) const; bool operator==(const NuXLFragmentAdductDefinition& other) const;
```
- **Expectation:** For a value type, operator< and operator== should induce a consistent equivalence: !(a<b)&&!(b<a) should mean a==b. A developer putting these in a std::set (the public alias NucleotideToFragmentAdductMap = map<char, set<NuXLFragmentAdductDefinition>>) and also comparing them with == or std::hash expects both to define 'duplicate' the same way.
- **Actual:** operator< orders by std::tie(mass, fa, name) (formula stringified), so it requires all three of mass/formula/name to be equal for equivalence. operator== compares only std::tie(formula, name) — mass is completely ignored. The std::hash specialization in the same header also hashes only formula+name and documents itself as 'consistent with operator==', but operator< (used by std::set) is NOT consistent with that equality. Two definitions with identical formula+name but slightly different mass are == and hash-equal, yet are kept as two distinct elements in a std::set.
- **Evidence:** operator<: 'return std::tie(mass, fa, name) < std::tie(other.mass, fb, other.name);' (NuXLFragmentAdductDefinition.cpp:19) vs operator==: 'return std::tie(formula, name) == std::tie(other.formula, other.name);' (NuXLFragmentAdductDefinition.cpp:24). Header doc on hash: '@note Hash is consistent with operator==.' Public alias 'using NucleotideToFragmentAdductMap = std::map<char, std::set<NuXLFragmentAdductDefinition> >;' (NuXLParameterParsing.h:43).
- **Fix:** Make the orderings consistent: have operator< order by exactly the fields operator== compares (formula, name) — i.e. drop mass from operator< (mass is derived from formula via getMonoWeight, so it is redundant for identity anyway). This is source-compatible and fixes the std::set-vs-hash dedup divergence. If preserving mass-first ordering is required for an existing on-disk/scoring order, instead add mass to operator== (and the hash). Either way pick one field set for both; do not leave them divergent.
- **Verifier correction:** The orderings are genuinely inconsistent (operator< keys on mass,formula,name; operator==/hash key on formula,name), but the consequence ('kept as two distinct std::set elements'/'set dedup disagrees') is only achievable when mass is hand-forged inconsistent with formula. In all real code paths mass = formula.getMonoWeight(), so mass is a deterministic function of formula and the two orderings coincide on every actually-constructed value — there is no current silent data loss or wrong scoring. The fix is still correct (drop mass from operator< so it keys on exactly formula+name, matching == and the documented hash), and it is source-compatible with no ABI break since signatures are unchanged. Severity is low (latent landmine / mild surprise), not high.

### [ANAL-54] NuXLFragmentAnnotationHelper::FragmentAnnotationDetail_::operator< / operator== — FragmentAnnotationDetail_ operator< uses exact mz/intensity while operator== uses a 1e-6 epsilon — ordering equivalence disagrees with equality
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFragmentAnnotationHelper.h · _nuxl-frag_

```cpp
bool operator<(const FragmentAnnotationDetail_&) const; bool operator==(const FragmentAnnotationDetail_&) const;
```
- **Expectation:** When a struct provides both < and ==, the strict-weak-ordering equivalence (!(a<b)&&!(b<a)) and == are expected to agree, so containers and std::unique/std::sort behave consistently.
- **Actual:** operator< compares std::tie(charge, shift, mz, intensity) using exact double comparison, whereas operator== treats mz and intensity as equal when their absolute difference is < 1e-6. Two details differing by, say, 1e-7 in mz are == but are strictly ordered by <, so the two relations contradict each other. The author even notes the mz/intensity terms in == are 'actually not needed' — a sign the two were not designed to match.
- **Evidence:** operator<: 'return std::tie(charge, shift, mz, intensity) < std::tie(other.charge, other.shift, other.mz, other.intensity);' (NuXLFragmentAnnotationHelper.h:49-50). operator==: 'double mz_diff = fabs(mz - other.mz); ... return (charge==other.charge && shift==other.shift && mz_diff < 1e-6 && intensity_diff < 1e-6);' with trailing comment '// mz and intensity difference comparison actually not needed but kept for completeness' (NuXLFragmentAnnotationHelper.h:55-57).
- **Fix:** Pick one comparison basis. Since the comment says mz/intensity are not needed for identity, simplify operator== to (charge, shift) and order operator< by the same key, or make both use the epsilon-tolerant form. Source-compatible (inline struct in header; only behavior of comparisons changes — audit any std::set<FragmentAnnotationDetail_> usage first).
- **Verifier correction:** The inconsistency between operator< (exact) and operator== (1e-6 epsilon) is real, but its practical impact is nil: operator== is never called anywhere in the codebase, and the struct is never stored in a std::set or run through std::unique/std::sort together with ==. The only sort in the helper (NuXLFragmentAnnotationHelper.cpp:81) is on PeptideHit::PeakAnnotation, not FragmentAnnotationDetail_. The defect is a latent style/convention smell (already noted in the source comment), not a source of silently-wrong results. abi_impact is none (inline member functions of an inline header struct; no exported symbol or layout change), and any fix would be source-compatible and behavior-neutral.

### [ANAL-55] NuXLFragmentIonGenerator::addShiftedImmoniumIons — The 'methionine-related CH5S fragment' is emitted with the same iM immonium annotation as the real M immonium, making the two peaks indistinguishable
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFragmentIonGenerator.h · _nuxl-frag_

```cpp
static void addShiftedImmoniumIons(..., PeakSpectrum::StringDataArray& partial_loss_spectrum_annotation);
```
- **Expectation:** Each distinct ion the generator appends carries a distinct, identifying annotation, so a downstream consumer can tell which chemical species produced a given peak (the class summary stresses annotation arrays parallel to the peaks).
- **Actual:** For residue M the method appends two peaks at different m/z (the C8-style immonium at 104.05285 and a separate CH5S methionine-related fragment) but annotates BOTH via getAnnotatedImmoniumIon('M', fragment_shift_name), producing the identical string 'iM+{shift}+' for both. The CH5S fragment is not labelled as CH5S, so the annotation cannot distinguish the two M-derived peaks.
- **Evidence:** Both branches call the same annotator: line 159 'getAnnotatedImmoniumIon('M', fragment_shift_name)' for the 104.05285 immonium and line 166 'getAnnotatedImmoniumIon('M', fragment_shift_name)' for 'EmpiricalFormula("CH5S").getMonoWeight() ... // methionine related fragment' (NuXLFragmentIonGenerator.cpp:153-167). Header doc: 'M additionally emits a methionine-related CH5S fragment.'
- **Fix:** Give the CH5S fragment its own annotation (e.g. 'iM(CH5S)+{shift}+', mirroring the iK(C5H10N1) variant convention used a few lines above) so peaks remain distinguishable. Behavioral/string-only change; ABI-neutral.

### [ANAL-31] PeptideAndProteinQuant::getStatistics / getPeptideResults / getProteinResults — Read-only result getters are non-const
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h · _quant-absolute-mdv_

```cpp
const Statistics& getStatistics(); const PeptideQuant& getPeptideResults(); const ProteinQuant& getProteinResults();
```
- **Expectation:** A 'get...Results()'/'getStatistics()' accessor that returns a const reference to internal data and does not modify the object should itself be const, so it can be called on a const PeptideAndProteinQuant.
- **Actual:** All three methods are declared (and defined) non-const. Implementation in PeptideAndProteinQuant.cpp:798-816 simply 'return stats_;' / 'return pep_quant_;' / 'return prot_quant_;' with no mutation, yet the methods lack the const qualifier.
- **Evidence:** Header: 'const Statistics& getStatistics();' / 'const PeptideQuant& getPeptideResults();' / 'const ProteinQuant& getProteinResults();'. Source: 'const PeptideAndProteinQuant::Statistics& PeptideAndProteinQuant::getStatistics() { return stats_; }' (no const).
- **Fix:** Add the const qualifier to all three methods. Source-compatible for callers; it does change the mangled symbol so it is technically an ABI break for the affected member functions. If strict ABI must be preserved, add new const overloads alongside the existing ones.
- **Verifier correction:** The code facts are exactly as claimed (header 182/185/188 and source 798-816 confirmed). The only adjustment is severity: this should be low, not implied-higher. The defect produces no wrong/silent results, data loss, or crashes — it merely prevents read-only getter access through a const PeptideAndProteinQuant and would fail loudly at compile time. All existing callers operate on non-const objects, so there is no current functional impact. Fix is correct (add const to all three); note it changes the mangled symbols and is therefore an ABI break for these members, though source-compatible for callers.

### [ANAL-32] AbsoluteQuantitation::getQuantMethods / getQuantMethodsAsMap — Quant-method getters are non-const
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h · _quant-absolute-mdv_

```cpp
std::vector<AbsoluteQuantitationMethod> getQuantMethods(); std::map<std::string, AbsoluteQuantitationMethod> getQuantMethodsAsMap();
```
- **Expectation:** Pure getters that only read internal state and return copies should be const and callable on a const AbsoluteQuantitation.
- **Actual:** Both are non-const. Source (AbsoluteQuantitation.cpp:95-108) only reads quant_methods_ and returns copies, but the methods are not const, so they cannot be called on a const instance.
- **Evidence:** 'std::vector<AbsoluteQuantitationMethod> AbsoluteQuantitation::getQuantMethods() { ... return quant_methods; }' and 'std::map<...> AbsoluteQuantitation::getQuantMethodsAsMap() { return quant_methods_; }' — neither marked const.
- **Fix:** Mark both methods const. Source-compatible; changes the member symbol (technically ABI). If ABI must hold, add const overloads.
- **Verifier correction:** Both getters are indeed non-const and only read state, so the defect is real. However, severity should be low rather than high/medium: the omission causes a loud compile-time error when invoked on a const instance, never a silent runtime fault. Adding const is source-compatible for all callers but changes the mangled member-function symbol, so it is technically ABI-breaking for the shared library; if strict ABI must be preserved, add const overloads alongside the existing non-const ones.

### [ANAL-37] AbsoluteQuantitation::setQuantMethods — Setter takes a non-const reference but does not modify it
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h · _quant-absolute-mdv_

```cpp
void setQuantMethods(std::vector<AbsoluteQuantitationMethod>& quant_methods)
```
- **Expectation:** A setter named setX(input) takes its input by const reference (or value). A non-const reference parameter signals to the caller that the argument may be modified.
- **Actual:** quant_methods is taken by non-const reference but the implementation only reads it (copies entries into quant_methods_ keyed by component name). This rejects rvalues/const vectors at call sites and falsely implies mutation.
- **Evidence:** Header: 'void setQuantMethods(std::vector<AbsoluteQuantitationMethod>& quant_methods);'. Source: loops 'quant_methods_[component_name] = quant_methods[i];' — read-only access.
- **Fix:** Change the parameter to 'const std::vector<AbsoluteQuantitationMethod>&'. This is source-compatible for existing callers and widens accepted arguments, but changes the mangled symbol (member ABI). If ABI must be preserved, add a const-ref overload.
- **Verifier correction:** The defect is real: setQuantMethods takes a non-const lvalue reference yet the implementation only reads it (and the header documents it @param[in]). The recommendation to switch to `const std::vector<AbsoluteQuantitationMethod>&` is correct and source-compatible for existing in-tree callers (all pass lvalues; the pyOpenMS binding already copies a const input to call it). However, it is an ABI break (the member function's mangled symbol changes), not 'source-compatible' at the ABI level — so this is best categorized as a const-correctness fix, severity low (loud compile-time rejection of rvalues/const args, no silent wrong results).

### [ANAL-40] DDAWorkflowCommons::recalibrateMS1 — peptide_ids documented '@param[in]' but passed by non-const reference
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/DDAWorkflowCommons.h · _quant-absolute-mdv_

```cpp
static void recalibrateMS1(MSExperiment & ms_centroided, PeptideIdentificationList& peptide_ids, const std::string & id_file_abs_path = "")
```
- **Expectation:** An '@param[in]' parameter should be const; a non-const reference tells the caller the argument may be mutated.
- **Actual:** peptide_ids is a non-const 'PeptideIdentificationList&' yet documented '@param[in]'. The function only forwards it to InternalCalibration::fillCalibrants (which takes it non-const for reading). The non-const reference rejects const arguments and contradicts the @param[in] tag.
- **Evidence:** Header: '@param[in] peptide_ids ... PeptideIdentificationList& peptide_ids'. Source: 'ic.fillCalibrants(peptide_ids, 25.0);' as the only use.
- **Fix:** If fillCalibrants can accept const (verify), make peptide_ids const-ref to match the @param[in]. Otherwise document it as in/out. Const-ref change is source-compatible but alters the function symbol (ABI).

### [ANAL-28] IsobaricChannelExtractor::getStats — getStats() is non-const and hands out a mutable reference to internal QC state
`low` · `const-correctness` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h · _quant-isobaric_

```cpp
ChannelQCSet& getStats()
```
- **Expectation:** A getter named getStats that 'Returns a reference to the channel statistics' for inspection would be expected to be const and to return a const reference, so callers cannot accidentally corrupt accumulated QC.
- **Actual:** getStats() is non-const and returns a non-const ChannelQCSet& aliasing the private member channel_mz_delta, so any caller can mutate (or clear) the extractor's internal statistics through the returned reference. There is no const overload, so even read-only callers must take a mutable handle.
- **Evidence:** IsobaricChannelExtractor.h:125 `ChannelQCSet& getStats();`; IsobaricChannelExtractor.cpp:935-938 `ChannelQCSet& IsobaricChannelExtractor::getStats() { return channel_mz_delta; }`. (The .cpp doc comment at :932-934 is also a wrong copy-paste of clearStats' 'Clears channel statistics'.)
- **Fix:** Add a `const ChannelQCSet& getStats() const;` overload (additive, ABI-safe) for read-only access, and fix the stale .cpp doc comment. Changing the existing non-const signature would be ABI-breaking, so keep it but prefer the const accessor in docs.
- **Verifier correction:** getStats() is a non-const-only getter returning a non-const ChannelQCSet& aliasing the private member channel_mz_delta; it lacks a const overload and is documented for inspection, so it is const-incorrect (cannot be called on a const object, and hands out a mutable handle to internal QC). However it has no callers in the codebase, so it invites misuse rather than causing silent corruption — low severity. The fix is additive: add `const ChannelQCSet& getStats() const;` (ABI-safe) and correct the copy-pasted .cpp doc comment at lines 932-934 (currently duplicates clearStats' "Clears channel statistics" text).

### [ANAL-56] PeakGroup::operator== — PeakGroup equality compares only mono mass + intensity, ignoring charge/score/decoy state
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h · _topdown-deconv_

```cpp
bool operator==(const PeakGroup& a) const
```
- **Expectation:** Two PeakGroup objects compare equal only if they represent the same deconvolution result (same mass, charge range, isotope cosine, Qscore, target/decoy type, scan, etc.). A competent caller using PeakGroup in std::set/std::unordered_set, std::find, or std::unique expects value identity.
- **Actual:** operator== returns true whenever monoisotopic_mass_ and intensity_ are bitwise-equal, regardless of charge range, scores, scan number, or TargetDecoyType. std::hash<PeakGroup> is deliberately built on the same two fields, so a target and a decoy with the same mass+intensity collide and hash-equal.
- **Evidence:** PeakGroup.cpp: `return this->monoisotopic_mass_ == a.monoisotopic_mass_ && this->intensity_ == a.intensity_;` and PeakGroup.h std::hash comment: `/// Hashes all fields used in operator== (monoisotopic mass and intensity)`.
- **Fix:** Document explicitly on operator== that equality is mass+intensity only (mass-identity semantics), or rename the concept (e.g. provide `sameMass()` and make operator== full value equality). Additive/safe fix: add a doc comment and a named predicate; changing operator== semantics is source-compatible but behaviorally breaking for existing containers, so guard it.
- **Verifier correction:** operator== compares only monoisotopic_mass_ and intensity_ (mass+intensity identity), consistent with operator</>/operator< and std::hash<PeakGroup>, rather than full value equality (charge range, scan number, scores, and target_decoy_type_ are ignored). The restricted semantics are documented at the std::hash site but not at the operator== declaration. The "target vs decoy hash-collision" is theoretical: in the actual std::set/std::map/unordered_set usages targets and decoys are not co-inserted, so no current caller gets wrong results. Recommended additive fix (add a doc comment on operator==, or a named sameMass() predicate) is reasonable; this is a documentation/naming issue, not a behavioral bug.

### [ANAL-60] DeconvolvedSpectrum::sortByQscore — sortByQscore sorts descending and keys on the 2D Q-score, not the 1D Qscore implied by the name
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h · _topdown-deconv_

```cpp
void sortByQscore()
```
- **Expectation:** By analogy with sort() ('sort by deconvolved monoisotopic masses', ascending) a caller expects sortByQscore() to order peak groups ascending by getQscore() (the 1D score the name references).
- **Actual:** It sorts in DESCENDING order (best first) using getQscore2D(), which returns max(qscore_, qscore2D_). So the ordering key is the 2D feature-level score and the direction is opposite to sort(). A caller relying on ascending order, or on the 1D score, gets surprising results.
- **Evidence:** DeconvolvedSpectrum.cpp sortByQscore: `std::sort(..., [](const PeakGroup& p1, const PeakGroup& p2){ return p1.getQscore2D() > p2.getQscore2D(); });`.
- **Fix:** Rename or document as sortByQscore2DDescending, and/or note the direction and that the key is getQscore2D. Doc-only fix is non-breaking; renaming is source-breaking.
- **Verifier correction:** sortByQscore() does sort descending (best-first) keyed on getQscore2D() rather than getQscore(), and the header doc only says "sort by Qscore", so the name/direction/key are imprecise. But getQscore2D() == max(1D,2D) degrades to the 1D score when no 2D score is set, descending-by-quality is a normal convention, and the only caller (FLASHDeconvAlgorithm.cpp:278) immediately re-sorts by mass (line 293), so no surprising results actually reach consumers. Recommend a doc-comment clarification (key = getQscore2D, order = descending); a rename is optional and would be source-breaking.

### [ANAL-61] SpectralDeconvolution::getIsotopeCosineAndIsoOffset — 'offset' is documented @param[in] but is actually an output (int&); other params are mislabeled too
`low` · `doc-direction-tag-mismatch` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h · _topdown-deconv_

```cpp
static float getIsotopeCosineAndIsoOffset(double mono_mass, const std::vector<float>& per_isotope_intensities, int& offset, const PrecalculatedAveragine& avg, int iso_int_shift, int window_width, const std::vector<double>& excluded_masses)
```
- **Expectation:** Parameter direction tags should tell the caller which arguments are inputs vs results. offset is passed by non-const reference and the doc/name say the function 'determines' the offset.
- **Actual:** offset is `int&` and is written by the function (it is the determined isotope offset), yet the Doxygen tags it `@param[in] offset output offset ...` (contradictory: tagged input, described as output). Similarly the sibling getCosine documents `@param[out] offset` for an int passed by value (an input). The direction tags are unreliable across this scoring API.
- **Evidence:** SpectralDeconvolution.h:108 `@param[in] offset output offset between input monoisotopic mono_mass and determined monoisotopic mono_mass`, signature `int& offset`. getCosine doc:99 `@param[out] offset element index offset` for a by-value `int offset`.
- **Fix:** Fix the direction tags: offset in getIsotopeCosineAndIsoOffset is [out] (or [in,out]); offset in getCosine is [in]. Doc-only, no ABI impact.
- **Verifier correction:** The Doxygen parameter-direction tags are inverted: getIsotopeCosineAndIsoOffset's `offset` is `int&` and is written (output) but tagged `@param[in]` (should be [out] or [in,out]); getCosine's `offset` is a by-value `int` read-only input but tagged `@param[out]` (should be [in]). Doc-only fix, no ABI impact. Severity is low, not a functional or API-reliability defect, because the correct direction is unambiguous from the C++ signatures themselves.

### [ANAL-62] DeconvolvedSpectrum::setPeakGroups — setPeakGroups takes a non-const lvalue reference but only copies it, implying a mutation/move that never happens
`low` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h · _topdown-deconv_

```cpp
void setPeakGroups(std::vector<PeakGroup>& x)
```
- **Expectation:** A setter taking `std::vector<PeakGroup>&` (non-const) signals that it will either modify the caller's vector or move from it; callers cannot pass a const vector or a temporary.
- **Actual:** The implementation does `std::vector<PeakGroup>().swap(peak_groups_); peak_groups_ = x;` - it makes a full copy and leaves x untouched. The non-const reference is gratuitous: it needlessly forbids passing const vectors/temporaries and misleads about ownership/mutation.
- **Evidence:** DeconvolvedSpectrum.cpp:325-329 `void DeconvolvedSpectrum::setPeakGroups(std::vector<PeakGroup>& x){ std::vector<PeakGroup>().swap(peak_groups_); peak_groups_ = x; }`.
- **Fix:** Change the parameter to `const std::vector<PeakGroup>&` (or add a `std::vector<PeakGroup>&&` move overload). Changing to const-ref is source-compatible for existing callers passing lvalues and removes the surprise; it is a signature change so flag ABI as source-compatible.
- **Verifier correction:** The claim is accurate about the code but overstates severity. setPeakGroups(std::vector<PeakGroup>& x) takes a non-const lvalue reference yet only copies from x (peak_groups_ = x) after clearing peak_groups_ via swap-with-empty; x is neither mutated nor moved-from. The non-const ref needlessly forbids const/temporary arguments and falsely implies mutation, making it the only setter in the class not using const& (cf. setOriginalSpectrum(const MSSpectrum&)). However this is a low-severity, compile-time-loud nuisance, not a source of silently wrong results or data loss. Recommended fix: change to const std::vector<PeakGroup>& (source-compatible for all current lvalue callers; ABI-breaking but source-compatible).

### [ANAL-63] FLASHDeconvAlgorithm::getDecoyAveragine — getDecoyAveragine silently returns the TARGET averagine when decoy reporting is disabled
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h · _topdown-deconv_

```cpp
const FLASHHelperClasses::PrecalculatedAveragine& getDecoyAveragine()
```
- **Expectation:** A getter named getDecoyAveragine returns the decoy averagine model (or signals that none exists). A caller using it to score/compare decoys expects a decoy-specific model.
- **Actual:** When report_decoy_ is false, the method returns getAveragine() (the target averagine) instead of a decoy model or an error. The substitution is silent and not stated in the header, so a caller can unknowingly score with the wrong (target) averagine.
- **Evidence:** FLASHDeconvAlgorithm.cpp:361-364 `... getDecoyAveragine() { if (!report_decoy_) return getAveragine(); return sd_noise_decoy_.getAveragine(); }`. Header doc only says `/// get calculated decoy averagine. Call after calculateAveragine is called.`
- **Fix:** Document the fallback (returns the target averagine when report_decoy_ is off) so callers know the returned model is not decoy-specific in that mode. Doc-only; non-breaking.
- **Verifier correction:** getDecoyAveragine() does silently return the target averagine (via getAveragine()) when report_decoy_ is false, because the noise-decoy SpectralDeconvolution (sd_noise_decoy_) is never initialized in that mode and would otherwise return an empty averagine. The header omits this. However, the only in-tree caller (FLASHDeconv.cpp) is unaffected since no decoy masses are produced when FDR reporting is off, so realistic blast radius is limited to external/pyOpenMS callers — severity low rather than medium. Fix is doc-only (note the fallback in the header), non-breaking.

### [ANAL-64] FLASHDeconvAlgorithm::getAveragine / SpectralDeconvolution::getAveragine — Read-only averagine getter is non-const, blocking use on const objects
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h · _topdown-deconv_

```cpp
const FLASHHelperClasses::PrecalculatedAveragine& getAveragine()
```
- **Expectation:** A pure getter that just returns an internal member (and per its own doc requires the caller to have already run calculateAveragine) should be const so it can be called on a const FLASHDeconvAlgorithm/SpectralDeconvolution.
- **Actual:** Both getAveragine() overloads are non-const even though SpectralDeconvolution::getAveragine simply `return avg_;` and FLASHDeconvAlgorithm::getAveragine simply forwards it. No lazy computation occurs (calculateAveragine is a separate call). The missing const needlessly prevents const usage and is inconsistent with the const getDecoyAveragine-style read-only intent.
- **Evidence:** SpectralDeconvolution.cpp:253-256 `const PrecalculatedAveragine& SpectralDeconvolution::getAveragine() { return avg_; }`; FLASHDeconvAlgorithm.cpp:356-359 forwards `sd_.getAveragine();`. Both declared without trailing const.
- **Fix:** Add const to both getAveragine() (and getDecoyAveragine()) since they perform no mutation. Adding const is source-compatible for callers and improves usability; mark ABI as source-compatible (mangled name changes).
- **Verifier correction:** All three overloads — SpectralDeconvolution::getAveragine(), FLASHDeconvAlgorithm::getAveragine(), and FLASHDeconvAlgorithm::getDecoyAveragine() — are non-const pure getters that perform no mutation (averagine is populated by the separate calculateAveragine() call). They should be marked const for consistency with the class's existing const getters (getTolerances() const, getNoiseDecoyWeight() const) and to permit use on const instances. The claim's statement that getDecoyAveragine reflects "const ... read-only intent" is inaccurate: getDecoyAveragine is also non-const. Practical impact is low — no const reference to these classes exists in the codebase, so nothing is currently blocked; this is a mild idiom inconsistency, not a correctness/safety issue. Adding const is source-compatible (callers unaffected); only mangled names change.

### [ANAL-65] MassFeatureTrace::findFeaturesAndUpdateQscore2D — is_decoy is an INPUT mode flag but is documented as @param[out]
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h · _topdown-misc_

```cpp
std::vector<FLASHHelperClasses::MassFeature> findFeaturesAndUpdateQscore2D(const PrecalculatedAveragine& averagine, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, int ms_level = 1, bool is_decoy = false)
```
- **Expectation:** A parameter tagged '@param[out] is_decoy' is an output the function writes to; a caller would expect to read it back after the call (e.g. whether decoys were found).
- **Actual:** is_decoy is a plain input bool that SELECTS which peak groups to process. The body uses it as input: 'if (is_decoy && pg.getTargetDecoyType() == ...::target) continue; if (!is_decoy && pg.getTargetDecoyType() != ...::target) continue;' and stores 'mass_feature.is_decoy = is_decoy;'. It is passed by value, so it cannot be an out-parameter at all.
- **Evidence:** Header doc: '@param[out] is_decoy if set, only process decoy spectra. otherwise only target spectra'. Impl MassFeatureTrace.cpp:71-72 and 133-134 read is_decoy to filter; line 222 copies it into the feature.
- **Fix:** Fix the Doxygen tag to '@param[in] is_decoy'. The wording 'if set, only process decoy spectra' already describes input behavior; only the [out] tag is wrong. Pure doc fix, no ABI impact. Optionally rename to process_decoys for clarity (source-compatible only if no callers use the name, but it is positional so safe).
- **Verifier correction:** The Doxygen direction tag is incorrect: is_decoy is an INPUT mode flag (pass-by-value bool, default false) that selects whether decoy or target peak groups are processed, not an output. Fix '@param[out] is_decoy' to '@param[in] is_decoy'. Pure documentation fix, no ABI impact. This is a low-severity doc inconsistency rather than a behavioral bug.

### [ANAL-67] FLASHHelperClasses::PrecalculatedAveragine::get — get(mass) returns a full IsotopeDistribution by value (copy) where a const ref is expected for a cache lookup
`low` · `return-value` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h · _topdown-misc_

```cpp
IsotopeDistribution get(double mass) const
```
- **Expectation:** An O(1) lookup into a precomputed cache named 'get' would return a const reference to the stored object; the class brief calls the accessors 'O(1) table lookups'.
- **Actual:** It returns a by-value copy of the cached IsotopeDistribution ('return isotopes_[massToIndex_(mass)];'). For a class whose entire purpose is fast averagine scoring in hot deconvolution loops, every call copies the whole isotope vector, which is a surprising and avoidable cost.
- **Evidence:** Decl 'IsotopeDistribution get(double mass) const;' (FLASHHelperClasses.h:143). Impl 'return isotopes_[massToIndex_(mass)];' (FLASHHelperClasses.cpp:143-146) returns a copy.
- **Fix:** Change return type to 'const IsotopeDistribution&' (the storage outlives the call and is never invalidated by const accessors). This is an ABI break (mangled return changes) and a potential source break for callers binding to a non-const value, so do it at a major version; alternatively add a getRef(mass) const& overload (purely additive).
- **Verifier correction:** Claim is factually correct: get(mass) returns the cached IsotopeDistribution by value (copying its peak vector) rather than as a const reference, and is called in deconvolution scoring hot loops. The correction is to severity, not substance: this is a low-severity avoidable-performance surprise, not a correctness issue — results are identical, nothing crashes or is silently wrong. Returning const IsotopeDistribution& (storage outlives the call and is never invalidated by const accessors) is the clean fix but an ABI break of the mangled return type; an additive const-ref getRef(mass) overload is the source/ABI-compatible alternative.

### [ANAL-68] FLASHHelperClasses::MassFeature::operator== (and std::hash specialization) — MassFeature equality and hashing compare only avg_mass by exact double-equality, ignoring all other fields
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h · _topdown-misc_

```cpp
bool operator==(const MassFeature& other) const { return avg_mass == other.avg_mass; }
```
- **Expectation:** operator== on a rich feature struct (charge, scans, intensities, index, ms_level, decoy flag, ...) compares feature identity; two distinct features at coincidentally equal average mass are not 'equal'.
- **Actual:** Two completely different mass features (different index, scans, charges, decoy status) compare equal whenever their avg_mass doubles are bit-for-bit equal, and they collide in std::hash. Exact == on floating-point avg_mass is also fragile. This makes MassFeature unsafe as a key in hash/ordered containers expecting identity semantics.
- **Evidence:** Header: 'bool operator==(...) const { return avg_mass == other.avg_mass; }' (lines 260-263); operator< / operator> also key on avg_mass only. std::hash 'Hash based on avg_mass (the field used in operator==)' returns OpenMS::hash_float(mf.avg_mass) (lines 348-351).
- **Fix:** Document that ordering/equality are mass-keyed (intentional for mass-trace merging) so callers do not assume identity semantics; if identity is ever needed, add an explicit sameFeature()/index-based comparator. Changing operator== semantics would be a behavioral/ABI-affecting break, so prefer the additive documented-comparator route.
- **Verifier correction:** MassFeature::operator== (and operator</>, std::hash) compare/key only on avg_mass by exact double equality, ignoring all other fields — a real but low-impact naming surprise. Crucially: (1) the std::hash is deliberately consistent with operator== (equal objects hash equal, the required invariant), so it is NOT a broken hash key; (2) these operators/hash are unused as identity keys anywhere in production — MassFeature lives in a plain std::vector and is never sorted/dedup'd/set-inserted/looked-up by them; they are exercised only in unit tests; (3) the hash was added explicitly for pyOpenMS dict compatibility, and mass-keyed ordering is a defensible mass-trace-merging convention. The recommended fix (document that ==/ordering are mass-keyed, and add an additive sameFeature()/index comparator if identity is ever needed) is additive and ABI-neutral. The claim's "unsafe hash key / two different features compare equal causing silently wrong results" framing overstates an impact that does not occur in this codebase.

### [ANAL-69] FLASHHelperClasses::getChargeMass — getChargeMass returns the proton mass narrowed to float in a mass-critical code path
`low` · `unit-or-index` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h · _topdown-misc_

```cpp
static float getChargeMass(bool positive_ioniziation_mode)
```
- **Expectation:** A charge-carrier (proton) mass used to convert m/z to neutral mass should be returned at full double precision, consistent with Constants::PROTON_MASS_U and the double-valued mass fields it feeds.
- **Actual:** It returns float, narrowing Constants::PROTON_MASS_U (a double, ~1.007276466 Da). getLogMz and LogMzPeak::getUnchargedMass then use this float in mass computations on doubles, introducing avoidable single-precision rounding (~1e-7 Da relative) into top-down masses where sub-ppm accuracy matters.
- **Evidence:** Decl 'static float getChargeMass(bool ...)' (FLASHHelperClasses.h:338); impl 'return (float)(positive_ioniziation_mode ? Constants::PROTON_MASS_U : -Constants::PROTON_MASS_U);' (FLASHHelperClasses.cpp:236-239). Used in getLogMz (cpp:241-244) and getUnchargedMass (cpp:207).
- **Fix:** Change the return type to double to match the constant and the double mass arithmetic. This is an ABI break (return type mangling) so schedule for a major version; the value change is tiny but removes a surprising precision loss in a precision-sensitive domain.
- **Verifier correction:** The code reading is correct, but the impact framing ("sub-ppm accuracy matters", precision-critical) overstates it. The narrowing yields a fixed ~5e-8 Da (~7e-5 ppm) error on the proton mass, scaling linearly with charge in Da but constant at ~7e-5 ppm relative. With the pipeline running at 5-10 ppm tolerances (SpectralDeconvolution.cpp:27) and Da-scale isotope spacing, this is ~5 orders of magnitude below any decision threshold and never changes a reported mass meaningfully or flips a result. Real but cosmetic: a tidy fix (return double to match PROTON_MASS_U) worth scheduling for a major version due to the ABI/mangling change, but low severity rather than mass-critical.
