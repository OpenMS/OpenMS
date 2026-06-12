# OpenMS POLS Audit — Batch 6: MATH + ML + COMPARISON

**Confirmed findings:** 64 (MATH 24, ML 25, COMPARISON 15) · 8 high · 33 medium · 23 low.

### [COMP-1] BinnedSharedPeakCount::operator()(const BinnedSpectrum&, const BinnedSpectrum&) — Documented @throw IncompatibleBinning never happens; mismatched binning is silently unchecked in release builds
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/COMPARISON/BinnedSharedPeakCount.h · _comp-binned_

```cpp
double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const override; // @throw IncompatibleBinning is thrown if the binning of the two input spectra are not the same
```
- **Expectation:** Per the header's @throw documentation, passing two spectra with different binning throws IncompatibleBinning, so a caller can rely on the call to validate compatibility and fail loudly on incompatible input.
- **Actual:** The implementation only checks binning via OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(...), ...). OPENMS_PRECONDITION is compiled out unless OPENMS_ASSERTIONS is defined (debug-only), so in a normal release build there is NO check at all. Even when active, it throws Exception::Precondition, not IncompatibleBinning. With incompatible bins the code proceeds to cwiseProduct on mismatched sparse vectors and returns a meaningless score.
- **Evidence:** Header: "@throw IncompatibleBinning is thrown if the binning of the two input spectra are not the same". Source (BinnedSharedPeakCount.cpp:52): OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(spec1, spec2), "Binned spectra have different bin size or spread"); Macros.h:70-94: OPENMS_PRECONDITION is defined to throw Exception::Precondition only under #ifdef OPENMS_ASSERTIONS, else expands to nothing. grep across src finds no use of IncompatibleBinning anywhere except these doc comments.
- **Fix:** Either (a) restore the documented contract by throwing Exception::IncompatibleBinning (or a real exception) unconditionally when !isCompatible(...), independent of OPENMS_ASSERTIONS, or (b) correct the header to drop the @throw IncompatibleBinning line and state that compatibility is only checked in assertion-enabled builds. Option (a) is source-compatible (adds a throw on a currently-UB path); (b) is doc-only.
- **Verifier correction:** Refinement only (not a contradiction): the claim says OPENMS_PRECONDITION is "debug-only." More precisely, on Linux/non-MSVC it is disabled even in Debug builds because config.h gates OPENMS_ASSERTIONS behind `#if (0)`; it is active only in MSVC Debug. Net effect for the documented contract is the same or worse: the @throw is never honored in any standard Linux build, release or debug. Recommended fix (a) — throw a real exception unconditionally when !isCompatible(...) — is source-compatible (adds a throw on a currently-UB path).

### [COMP-3] BinnedSumAgreeingIntensities::operator()(const BinnedSpectrum&, const BinnedSpectrum&) — Documented @throw IncompatibleBinning never thrown; incompatible bins silently produce a garbage score in release builds
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/COMPARISON/BinnedSumAgreeingIntensities.h · _comp-binned_

```cpp
double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const override; // @throw IncompatibleBinning ...
```
- **Expectation:** Header says passing spectra with different binning throws IncompatibleBinning, so callers can pass arbitrary pairs and trust the call to validate.
- **Actual:** Only OPENMS_PRECONDITION guards compatibility, which is a no-op in release builds (no OPENMS_ASSERTIONS) and throws Exception::Precondition (not IncompatibleBinning) in debug. With mismatched bins the elementwise add/subtract on differently-sized sparse vectors yields an undefined/meaningless score with no error.
- **Evidence:** Header BinnedSumAgreeingIntensities.h:59: "@throw IncompatibleBinning ...". Source BinnedSumAgreeingIntensities.cpp:52: OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(spec1, spec2), ...). Macros.h:87-94: OPENMS_PRECONDITION expands to empty when OPENMS_ASSERTIONS undefined.
- **Fix:** Throw a real exception (Exception::IncompatibleBinning / IllegalArgument) unconditionally on !isCompatible, or fix the header to remove the false @throw guarantee. Adding an unconditional throw on a currently-UB path is source-compatible.

### [COMP-4] BinnedSpectralContrastAngle::operator()(const BinnedSpectrum&, const BinnedSpectrum&) — Documented @throw IncompatibleBinning never thrown; division by sqrt of empty-spectrum norms can return NaN with no error
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/COMPARISON/BinnedSpectralContrastAngle.h · _comp-binned_

```cpp
double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const override; // @throw IncompatibleBinning ...
```
- **Expectation:** Header promises IncompatibleBinning on mismatched binning, and a similarity functor would be expected to either validate input or return a defined score for degenerate input.
- **Actual:** Compatibility is only checked via the debug-only OPENMS_PRECONDITION (no-op in release; throws Exception::Precondition not IncompatibleBinning in debug). Additionally, for an empty spectrum (all-zero bins) sum1 or sum2 is 0, so score = numerator/sqrt(0) = 0/0 = NaN is returned silently instead of a defined value or an error.
- **Evidence:** Header BinnedSpectralContrastAngle.h:53. Source BinnedSpectralContrastAngle.cpp:55-58: sum1=dot(a,a); sum2=dot(b,b); score = numerator/(sqrt(sum1*sum2)); — no guard against sum1==0 or sum2==0. Macros.h:94 shows OPENMS_PRECONDITION is empty in release.
- **Fix:** Throw unconditionally on incompatible binning, and guard the normalization (return 0.0 when either norm is 0) so the score never silently becomes NaN. Both are additive/source-compatible. Also correct the header @throw line.
- **Verifier correction:** Minor refinement: the documented IncompatibleBinning is not merely "thrown as a different exception type" — that named exception class does not exist anywhere in OpenMS (it appears only in three @throw doc comments). So the contract is unfulfillable as written. The compatibility check is debug-only (OPENMS_PRECONDITION, no-op in release; throws Exception::Precondition, not IncompatibleBinning, in debug). The unguarded normalization at BinnedSpectralContrastAngle.cpp:58 returns NaN (0/0) for an empty/all-zero spectrum. Recommended fix is source-compatible: throw a real exception unconditionally on incompatible binning, return 0.0 when either norm is 0, and correct/remove the bogus @throw line in the header (same issue exists in sibling headers BinnedSumAgreeingIntensities.h:59 and BinnedSharedPeakCount.h:55).

### [COMP-6] SpectraSTSimilarityScore::dot_bias — dot_bias default -1 does NOT trigger recomputation; only 0 does
`high` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/COMPARISON/SpectraSTSimilarityScore.h · _comp-spectra_

```cpp
double dot_bias(const BinnedSpectrum & bin1, const BinnedSpectrum & bin2, double dot_product = -1) const
```
- **Expectation:** The Doxygen says '@param[in] dot_product if -1 this value will be calculated as well.' and the default argument is -1, so a caller passing the default (or explicitly -1) expects the dot product to be computed internally.
- **Actual:** The implementation branches on `if (dot_product != 0)` and only recomputes via `(*this)(bin1, bin2)` when dot_product == 0. With the documented/default sentinel -1, it takes the `dot_product != 0` branch and divides by -1, returning a negative, meaningless dot_bias. The recompute path is only reached for dot_product == 0, the value the doc does not mention.
- **Evidence:** Header: `double dot_bias(... double dot_product = -1) const;` doc '@param[in] dot_product if -1 this value will be calculated as well.' vs cpp: `if (dot_product != 0) { return (double)numerator / dot_product; } else { return (double)numerator / (*this)(bin1, bin2); }`
- **Fix:** Fix the implementation to honor the documented sentinel: branch on `if (dot_product >= 0)` (or `!= -1`) and recompute otherwise; or change the default and doc to 0. ABI-neutral (implementation/doc fix). If changing the default value, that is source-compatible only.

### [MATH-7] Math::absdev — absdev ("absolute deviation") sums SIGNED deviations, not absolute ones - always ~0 for the default mean
`high` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/StatisticFunctions.h · _math-stats_

```cpp
template <typename IteratorType> static double absdev(IteratorType begin, IteratorType end, double mean = std::numeric_limits<double>::max())
```
- **Expectation:** A function named 'absdev' / documented as 'absolute deviation' returns the mean absolute deviation, i.e. mean(|x_i - mean|), a positive spread measure (like MeanAbsoluteDeviation a few lines above which DOES use fabs).
- **Actual:** The loop does `sum_value += *iter - mean;` with NO fabs. When `mean` is the actual mean (the default path computes it internally), the signed deviations cancel and the result is ~0 (up to floating point error). It is mathematically a near-zero number, never the absolute deviation.
- **Evidence:** Line 581-584: `for (IteratorType iter=begin; iter!=end; ++iter) { sum_value += *iter - mean; } return sum_value / std::distance(begin, end);` -- contrast with MeanAbsoluteDeviation at line 216 `mean_value += fabs(*it - mean_of_numbers);`. Caller src/topp/MRMPairFinder.cpp:286 uses it as a +/- spread estimate.
- **Fix:** This is a latent bug, not just a naming issue. Add fabs: `sum_value += std::fabs(*iter - mean);` so it actually computes mean absolute deviation. Since the current result is essentially meaningless (~0), fixing it is source/ABI-compatible at the signature level; only the (broken) numeric output changes. Alternatively deprecate in favor of the already-correct MeanAbsoluteDeviation.

### [MATH-14] Math::MultipleTesting::computeModelFDR — computeModelFDR returns an all-NaN vector (no throw) when any input PEP is NaN
`high` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/MATH/STATISTICS/MultipleTesting.h · _math-stats_

```cpp
template <class T> static std::vector<double> computeModelFDR(const std::vector<T>& data_in)
```
- **Expectation:** A method named computeModelFDR taking posterior error probabilities would, on a malformed input (a single NaN PEP), either skip it or signal an error.
- **Actual:** If ANY element is NaN, the function returns immediately with a vector of all-NaN of the same length and no diagnostic - one bad value silently wipes every FDR estimate. A caller that does not check for NaN gets silently corrupted results for the whole set.
- **Evidence:** Lines 200-208: `bool any_nan = false; ... if (any_nan) { return fdr; // propagate }` where `fdr` was initialized to all quiet_NaN at line 191.
- **Fix:** Document the all-or-nothing NaN propagation in the brief, or filter NaNs like qValue() does instead of poisoning the whole output. Doc-only fix is ABI-safe; filtering changes output values (source-compatible).

### [ML-10] ClusteringGrid::getClusters — getClusters dereferences map end() (UB / crash) for an empty cell, undocumented
`high` · `surprising-throw` · ABI: `source-compatible` · src/openms/include/OpenMS/ML/CLUSTERING/ClusteringGrid.h · _ml-clustering-b_

```cpp
std::list<int> getClusters(const CellIndex &cell_index) const
```
- **Expectation:** A const getter named getClusters returning a std::list<int> for a cell should, for a cell with no clusters, return an empty list (the natural reading of '@return list of cluster indices ... centred in this cell').
- **Actual:** Implementation is `return cells_.find(cell_index)->second;` with no existence check. For a cell that is not present in `cells_`, `find` returns `end()` and the `->second` dereference is undefined behavior (typically a crash). The only reason internal callers are safe is that every call site first guards with isNonEmptyCell(); a public caller has no way to know this from the signature.
- **Evidence:** ClusteringGrid.cpp:66-69 `std::list<int> ClusteringGrid::getClusters(const CellIndex &cell_index) const { return cells_.find(cell_index)->second; }`. Internal call sites in GridBasedClustering.h:296-298 and 592-594 always pre-check isNonEmptyCell.
- **Fix:** Guard the lookup: return an empty list when the cell is absent (additive, ABI-safe, no signature change). At minimum document the precondition that the cell must be non-empty.
- **Verifier correction:** The behavior is UB/crash (dereferencing map end()), not strictly a thrown exception; the category label "surprising-throw" is the closest available bucket but the actual failure mode is undefined behavior rather than a documented throw. Severity high because a public const getter crashes for a reasonable input (an empty/absent cell) with no loud signal. The recommended fix (guard and return an empty list) changes only the function body and is source-compatible — no signature/ABI change.

### [ML-21] ROCCurve::cutoffNeg — cutoffNeg() iterates only over positive samples but divides by the negative count
`high` · `misleading-name` · ABI: `source-compatible` · src/openms/source/ML/ROCCURVE/ROCCurve.cpp · _ml-interp-solvers_

```cpp
double cutoffNeg(double fraction = 0.95)
```
- **Expectation:** cutoffNeg(fraction) -- the symmetric partner of cutoffPos -- should compute a cutoff based on the negative class (true negatives among negatives), as the name and 'trueNeg' variable suggest.
- **Actual:** The loop only enters its body when 'cit->second' is true (i.e. the sample is a POSITIVE), yet increments a counter named trueNeg and divides it by neg_ (the negative count). It therefore mixes positive samples with the negative population size; the returned score is computed from the wrong class. cutoffPos uses the identical 'if (cit->second)' guard with truePos/pos_, which is correct there, making cutoffNeg's copy-paste error clear.
- **Evidence:** double ROCCurve::cutoffNeg(double fraction) { ... UInt trueNeg = 0; for (...) { if (cit->second) { if ((double)trueNeg++ / neg_ > 1 - fraction) { return cit->first; } } } return -1; }
- **Fix:** The negative-class loop should test '!cit->second'. Fix the predicate to '!cit->second' (and verify against the rocN convention). Source-compatible behavioral fix; document the threshold semantics. The class is already marked '[buggy and usage is discouraged]', which corroborates this.

### [COMP-7] SpectraSTSimilarityScore::operator()(const BinnedSpectrum&, const BinnedSpectrum&) — BinnedSpectrum operator() returns un-normalized dot while PeakSpectrum overload returns normalized dot
`medium` · `asymmetric-api` · ABI: `none` · src/openms/source/COMPARISON/SpectraSTSimilarityScore.cpp · _comp-spectra_

```cpp
double operator()(const BinnedSpectrum & bin1, const BinnedSpectrum & bin2) const
```
- **Expectation:** Both `operator()` overloads are documented identically ('calculates the dot product of the two spectra') and a caller would expect the same score space (a cosine/normalized dot in [0,1]) regardless of whether they pass PeakSpectrum or BinnedSpectrum.
- **Actual:** The PeakSpectrum overload normalizes each spectrum's bins to unit norm before the dot product, yielding a cosine in [0,1]. The BinnedSpectrum overload returns the raw dot of the bins with no normalization, so its magnitude depends on the (un)normalized input vectors. Callers mixing the two overloads get inconsistent score ranges.
- **Evidence:** PeakSpectrum overload: `*bin1.getBins() /= bin1.getBins()->norm(); *bin2.getBins() /= bin2.getBins()->norm(); return bin1.getBins()->dot(*bin2.getBins());` vs BinnedSpectrum overload: `return bin1.getBins()->dot(*bin2.getBins());`
- **Fix:** Document explicitly that the BinnedSpectrum overload assumes pre-normalized bins (as produced by transform()), and rename or annotate to make the precondition obvious; or normalize inside the BinnedSpectrum overload for symmetry. Doc/impl change; normalizing would change returned values (source-compatible behavior change).
- **Verifier correction:** The two operator() overloads have identical, precondition-silent docstrings but different score spaces: the PeakSpectrum overload normalizes internally (cosine in [0,1]); the BinnedSpectrum overload performs a raw, un-normalized dot. It implicitly assumes its input was produced by transform() (which normalizes) — the documented workflow and the unit test both feed it transform() output, so the score spaces coincide along the intended path. A caller who builds a BinnedSpectrum directly and passes it gets a silently un-normalized score, inconsistent with the PeakSpectrum overload and with the class doc's stated [0,1] range. Fix: document the pre-normalization precondition on the BinnedSpectrum overload (or normalize inside it for symmetry).

### [COMP-8] SpectrumPrecursorComparator::operator()(const PeakSpectrum&, const PeakSpectrum&) — Precursor comparator returns a Da-valued 'window - |delta|', not a normalized similarity
`medium` · `return-value` · ABI: `none` · src/openms/source/COMPARISON/SpectrumPrecursorComparator.cpp · _comp-spectra_

```cpp
double operator()(const PeakSpectrum & a, const PeakSpectrum & b) const
```
- **Expectation:** A PeakSpectrumCompareFunctor 'similarity' (base doc says value should be >= 0, and siblings like Zhang/SpectraST return scores in [0,1]) — a caller expects a comparable similarity, ideally normalized, where identical inputs give a fixed maximal score.
- **Actual:** It returns `window - fabs(mz1 - mz2)`, i.e. a value in Da that ranges in (0, window], where the maximum (perfect match) equals the arbitrary 'window' parameter (default 2). The self-similarity operator()(spec) returns `window`, not 1. The score scale silently depends on a tolerance parameter and is in mass units, not a normalized similarity.
- **Evidence:** `if (fabs(mz1 - mz2) > window) { return 0; } return window - fabs(mz1 - mz2);` and constructor `defaults_.setValue("window", 2, ...)`
- **Fix:** Document the score space explicitly (returns window-minus-distance in Da, max == window) or normalize by window to give [0,1]. Doc fix is ABI-neutral; normalization changes returned values (source-compatible behavior change).

### [COMP-9] SpectrumPrecursorComparator::operator()(const PeakSpectrum&, const PeakSpectrum&) — Missing precursor silently treated as m/z 0, producing a meaningless comparison
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/source/COMPARISON/SpectrumPrecursorComparator.cpp · _comp-spectra_

```cpp
double operator()(const PeakSpectrum & a, const PeakSpectrum & b) const
```
- **Expectation:** Comparing precursors of two spectra where one or both lack a precursor should signal that the comparison is undefined (throw or documented sentinel), since a precursor m/z of 0 is not a real value.
- **Actual:** If a spectrum has no precursors, its m/z defaults to 0.0 and the comparison proceeds as if the precursor were at 0 Da. Two spectra both lacking precursors compare as a perfect match (|0-0|=0 -> returns window); a real spectrum vs. a precursor-less one returns 0 only because the distance exceeds the window. The caller is never told a precursor was missing.
- **Evidence:** `double mz1 = 0.0; if (!x.getPrecursors().empty()) { mz1 = x.getPrecursors()[0].getMZ(); } double mz2 = 0.0; if (!y.getPrecursors().empty()) { mz2 = y.getPrecursors()[0].getMZ(); }`
- **Fix:** Document the no-precursor behavior, or throw/return a clear sentinel when a precursor is absent. ABI-neutral; throwing would be a source-compatible behavior change.
- **Verifier correction:** Severity adjusted from high to medium. The claim is otherwise accurate: a missing precursor is silently treated as m/z 0, and because this is a similarity functor (higher = more similar), two precursor-less spectra are scored as MAXIMALLY similar (the full window value), while a real-vs-missing pair returns 0 only incidentally. The behavior is undocumented and untested. Recommend documenting the no-precursor behavior or throwing/returning a clear sentinel when a precursor is absent; this would be a source-compatible behavior change (signature unchanged).

### [COMP-11] SpectrumCheapDPCorr::operator()(const PeakSpectrum&, const PeakSpectrum&) — const comparison operator mutates internal consensus/peak-map state and resets the weighting factor
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/COMPARISON/SpectrumCheapDPCorr.h · _comp-spectra_

```cpp
double operator()(const PeakSpectrum & a, const PeakSpectrum & b) const override
```
- **Expectation:** A const `operator()` that computes a similarity score is expected to be a pure read-only computation with no observable side effects on the functor.
- **Actual:** The const call rebuilds `lastconsensus_`, fills `peak_map_`, and at the end resets `factor_` back to 0.5 — all via `mutable` members. The accessors `lastconsensus()` / `getPeakMap()` depend on the last call, and a `setFactor()` value is silently consumed and discarded by the very next `operator()`. This is surprising hidden state for a const, score-returning call.
- **Evidence:** Members `mutable PeakSpectrum lastconsensus_; mutable double factor_; mutable std::map<UInt, UInt> peak_map_;`; cpp: `lastconsensus_ = PeakSpectrum(); ... peak_map_.clear(); ... factor_ = 0.5; return score;`
- **Fix:** Document the side effects on lastconsensus()/getPeakMap() and the one-shot nature of setFactor() prominently in the header; ideally provide an overload that returns the consensus explicitly rather than via mutable state. Doc fix ABI-neutral; API redesign would be breaking.

### [COMP-12] SpectrumCheapDPCorr::setFactor — setFactor only affects the next single comparison, then auto-resets to 0.5
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/COMPARISON/SpectrumCheapDPCorr.h · _comp-spectra_

```cpp
void setFactor(double f)
```
- **Expectation:** A setter named `setFactor` is expected to configure a persistent weighting used by all subsequent comparisons.
- **Actual:** The value set is used by exactly one following `operator()` call and is then reset to 0.5 at the end of that call, so a second comparison silently reverts to the default. The header comment '/// set weighting of the second spectrum for consensus from next function call operator' hints at this but the method name does not.
- **Evidence:** cpp `operator()` ends with `factor_ = 0.5;`; header member doc `/// weighting factor for the next consensus spectrum`
- **Fix:** Rename or clearly document as a one-shot setting (e.g. setFactorForNextCall), or make the factor persistent. Doc fix ABI-neutral; rename is source-compatible if old name retained.
- **Verifier correction:** setFactor sets a one-shot weighting: the value is consumed by the next operator() call and then reset to 0.5 at the end of that call (cpp:199), so subsequent comparisons silently revert to 0.5. However, factor_ only influences the intensities of the produced consensus spectrum (lastconsensus()), not the returned similarity score, which is computed independently in comparepeaks_/dynprog_. The reset-to-0.5 behavior is not documented (only "next" hints at one-shot scope). Recommend documenting it as one-shot (and/or renaming to setFactorForNextConsensus with the old name retained); a doc clarification is ABI-neutral.

### [COMP-14] SpectrumAlignment::getSpectrumAlignment — is_relative_tolerance silently switches between symmetric DP alignment and asymmetric nearest-neighbor matching
`medium` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/COMPARISON/SpectrumAlignment.h · _comp-spectra_

```cpp
template <...> void getSpectrumAlignment(...) const
```
- **Expectation:** A single 'tolerance' parameter would be expected to only change the width of one consistent alignment algorithm; the same method should produce comparable alignment semantics regardless of unit.
- **Actual:** When is_relative_tolerance is false (absolute Da), a banded Needleman-Wunsch-style DP minimizing m/z distance is used (each s1 peak aligns to at most one s2 peak and vice versa). When true (ppm), the code instead does a MatchedIterator nearest-neighbor match where a peak in s2 may be matched to multiple peaks in s1 and the relationship is asymmetric. The output pairing semantics change fundamentally based on a tolerance-unit flag.
- **Evidence:** Header class doc 'Method 1: ... banded alignment ... Method 2: ... a peak in s2 can be matched to none, one or multiple peaks in s1'; cpp branch `if (!param_.getValue("is_relative_tolerance").toBool()) { ... DP ... } else { MatchedIterator<...> it(s1, s2, tolerance); for (...) alignment.emplace_back(it.refIdx(), it.tgtIdx()); }`
- **Fix:** The class doc does describe both methods; surface this on the method itself and consider distinct method names for the two semantics to avoid silent algorithm swaps. Doc-only; ABI-neutral.

### [COMP-15] SpectrumCheapDPCorr::comparepeaks_ / int_cnt parameter — Invalid int_cnt value yields a -1 peak score silently folded into the returned similarity
`medium` · `silent-failure` · ABI: `none` · src/openms/source/COMPARISON/SpectrumCheapDPCorr.cpp · _comp-spectra_

```cpp
double operator()(const PeakSpectrum&, const PeakSpectrum&) const  (param int_cnt)
```
- **Expectation:** An out-of-range mode parameter should be rejected (validated/throw) rather than silently producing a negative per-peak contribution to a similarity that is otherwise non-negative.
- **Actual:** comparepeaks_ returns -1 (a sentinel, with a TODO acknowledging a missing exception) for any int_cnt outside [0,3]; this -1 is summed into the score returned by operator(), corrupting the result with no error. The int_cnt parameter is a plain integer with no valid-strings restriction.
- **Evidence:** cpp `else { // TODO exception ... return -1; }` and constructor `defaults_.setValue("int_cnt", 0, "...0 = product...3 = agreeing intensity");` (no setValidStrings / range check)
- **Fix:** Validate int_cnt at parameter-set time (restrict to {0,1,2,3}) or throw on the invalid branch instead of returning -1. ABI-neutral implementation fix.
- **Verifier correction:** int_cnt is not a parameter of operator(); it is a Param-stored value read inside the private comparepeaks_ helper. Otherwise the claim is accurate: the out-of-range else branch returns the sentinel -1 (with '// TODO exception'), which is summed into the score via operator() line 191 and dynprog_ line 223, and the int_cnt Param has no setValidStrings/min-max restriction (line 31). Triggering requires deliberately setting int_cnt outside the documented {0,1,2,3}, hence medium rather than a higher tier.

### [MATH-16] GaussFitter::GaussFitResult::log_eval_no_normalize — log_eval_no_normalize doc claims it returns scaled 'intensities' (PDF*A) but it returns an UN-scaled log-density that ignores A
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/MATH/STATISTICS/GaussFitter.h · _math-fitters_

```cpp
double log_eval_no_normalize(double x) const
```
- **Expectation:** A method named log_eval_no_normalize, documented as 'Returns the intensities (i.e. probabilities scaled by the factor A) of the PDF', should return log of the A-scaled intensity that eval() produces — i.e. log(eval(x)).
- **Actual:** The implementation returns -log(sigma) - 0.5*log(2*pi) - 0.5*((x-x0)/sigma)^2, the log of a standard-normalized Gaussian PDF. It does NOT multiply by A (the amplitude member is never read), so it is not the log of eval(). The doc block is a verbatim copy-paste of eval()'s and is wrong for this function.
- **Evidence:** Header doc: 'Returns the intensities (i.e. probabilities scaled by the factor 'A') of the PDF'. Impl (GaussFitter.cpp:139-145): `return -log(sigma) - halflogtwopi - 0.5 * pow((x - x0) / sigma, 2.0);` — A is absent. By contrast eval() (line 132-137) computes `A / pdf(ndf,x0)` scaling.
- **Fix:** Fix the Doxygen block to state it returns the natural-log density of the *normalized* Gaussian (A is ignored), and that it is NOT log(eval()). Rename ideally to log_pdf or logDensityNormalized; that rename is ABI-breaking (additive deprecated alias is safe). The doc-only fix is fully compatible.
- **Verifier correction:** The recommended remedy is a Doxygen-only fix (ABI-impact none): state that the method returns the natural log of the NORMALIZED Gaussian density (A is ignored) and that it is NOT log(eval()). A rename to log_pdf / logDensityNormalized would be ABI-breaking but is optional; an additive deprecated alias is the safe path. Severity downgraded from any "silently-wrong-results" framing to medium: the sole callers use it correctly, the function name and inline TODOs partly signal the A-independence, and no existing behavior is broken — but the copy-pasted doc actively invites silent misuse by a developer who trusts it to mean log(A-scaled intensity).

### [MATH-19] PosteriorErrorProbabilityModel::computeProbability — computeProbability silently applies a hidden score offset and returns garbage if fit() was not run, instead of signaling
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h · _math-fitters_

```cpp
double computeProbability(double score) const
```
- **Expectation:** computeProbability(score) takes the caller's raw score and returns a posterior error probability in [0,1]; calling it on an unfitted model would error or be obviously invalid.
- **Actual:** It first does `score = score + fabs(smallest_score_) + 0.001;` (a hidden re-application of the fit-time shift), then divides densities. If fit() was never called, smallest_score_=0 and the fit params are default/garbage (sigma=-1), so it silently returns a meaningless number with no error. The required pre-call ordering is documented only in a @note.
- **Evidence:** PEP cpp:576-602: `score = score + fabs(smallest_score_) + 0.001;` then `return (negative_prior_ * x_neg) / (...)`. Header line 208-209 @note: 'fit has to be used before using this function. Otherwise this function will compute nonsense.' Default ctor (cpp:41) sets smallest_score_=0; GaussFitResult default sigma=-1.
- **Fix:** Track a 'fitted_' flag and throw Exception::Precondition (or return NaN documented as such) when computeProbability is called before a successful fit. The implicit +fabs(smallest_score_)+0.001 transform should be documented at the method (not just the @note). A flag member is ABI-breaking; a thrown precondition using existing members is source-compatible.
- **Verifier correction:** computeProbability(score) is documented (header @note) to require a prior successful fit(); the `score + fabs(smallest_score_) + 0.001` line is the intended inverse of fit()'s domain shift, not a hidden corruption -- for a fitted model it is correct. The real, narrower surprise: on a default-constructed/unfitted model (smallest_score_=0, sigma=-1) there is no fitted_ flag or precondition check, so the method silently returns a meaningless value (typically ~1.0 or NaN) in [0,1] with no error, the only safeguard being an @note that says it 'will compute nonsense.' Recommend a precondition throw guarding the unfitted state (e.g. sigma<0) in the .cpp body using existing members (source-compatible); a new fitted_ member would be ABI-breaking.

### [MATH-21] PosteriorErrorProbabilityModel::fillDensities — fillDensities fills the 'incorrect' density with the Gaussian model, not the Gumbel, regardless of which fit() was run
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h · _math-fitters_

```cpp
void fillDensities(const std::vector<double>& x_scores, std::vector<double>& incorrect_density, std::vector<double>& correct_density)
```
- **Expectation:** After a Gumbel+Gauss fit (fitGumbelGauss), fillDensities should evaluate the incorrect-component density using the fitted Gumbel parameters.
- **Actual:** fillDensities always evaluates incorrect_density via incorrectly_assigned_fit_param_.eval() (a Gaussian), even when the incorrect component was fitted as a Gumbel; a code comment admits 'incorrect is currently filled with gauss as fitting gumble is not supported'. The name promises 'the densities' but silently substitutes a Gaussian for the negative component.
- **Evidence:** PEP cpp:454-461: `*incorrect = incorrectly_assigned_fit_param_.eval(score);` with comment 'TODO: incorrect is currently filled with gauss as fitting gumble is not supported'. Only the separate fillLogDensitiesGumbel uses incorrectly_assigned_fit_gumbel_param_.
- **Fix:** Document at the declaration that fillDensities/fillLogDensities always use the Gaussian negative model and that the Gumbel densities require fillLogDensitiesGumbel; or branch on the actual fitted model. Doc-only fix is fully compatible.
- **Verifier correction:** fillDensities and fillLogDensities always evaluate the incorrect/negative component via the Gaussian member incorrectly_assigned_fit_param_; they ignore the Gumbel member incorrectly_assigned_fit_gumbel_param_. This is correct/consistent for the standard Gaussian fit() (including its 'Gumbel' mode, which still fits a Gaussian and only plots a Gumbel formula). It is broken specifically after fitGumbelGauss: the fitted Gumbel is stored only in incorrectly_assigned_fit_gumbel_param_, while incorrectly_assigned_fit_param_ remains at its constructor default, so fillDensities returns densities from an uninitialized Gaussian instead of the fitted Gumbel. For Gumbel negative densities use fillLogDensitiesGumbel. Recommended fix is doc-only (document the Gaussian-only behavior at the header declaration and in the pyOpenMS docstring) or branch on the fitted model; either is source- and ABI-compatible.

### [MATH-2] EmgGradientDescent::estimateEmgParameters — Returns iteration count with no failure signal; on early numerical break the out-params can be left as DBL_MAX
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/MISC/EmgGradientDescent.h · _math-spline_

```cpp
UInt estimateEmgParameters(const std::vector<double>& xs, const std::vector<double>& ys, double& best_h, double& best_mu, double& best_sigma, double& best_tau) const
```
- **Expectation:** An estimator with out-params best_h/best_mu/best_sigma/best_tau and a UInt return should either always fill plausible parameters or signal failure. A caller reads best_* unconditionally.
- **Actual:** best_h/best_mu/best_sigma/best_tau are initialized to std::numeric_limits<double>::max(). They are only overwritten when an iteration achieves `current_E < best_E`. If the very first Loss_function evaluation is NaN/Inf (or h/mu/sigma/tau start invalid), the loop `break`s before any improvement, leaving all four out-params at DBL_MAX. The return value (iteration index) gives no indication of this; the caller cannot distinguish a 1-iteration success from a 1-iteration failure.
- **Evidence:** EmgGradientDescent.cpp line 584: `previous_E = best_h = best_mu = best_sigma = best_tau = best_E = std::numeric_limits<double>::max();` and lines 635-640 break on invalid current_E before the `if (current_E < best_E)` update at 643-651. Return is `return iter_idx;` (line 725), documented only as "number of iterations".
- **Fix:** Return a bool/enum success flag (or set out-params to NaN on failure and document it) so callers can detect divergence. Adding a separate `bool` overload or out-param is ABI-additive; changing the return type is breaking.
- **Verifier correction:** The defect is confirmed: on an early numerical break (invalid initial h/mu/sigma/tau, or a NaN/Inf first Loss_function evaluation) the out-params best_h/best_mu/best_sigma/best_tau are left at std::numeric_limits<double>::max(), and the UInt return (iter_idx, the break index) gives no way to detect this. The lone caller fitEMGPeakModel reads them unconditionally and stores them as peak metadata, silently propagating DBL_MAX. Two refinements to the claim: (1) the caller does not even inspect the return value, so 'the caller cannot distinguish 1-iteration success from failure' is moot for the only existing caller — it ignores the count entirely; (2) the high-severity DBL_MAX leak specifically requires failure on the FIRST iteration; divergence after >=1 improving iteration instead yields plausible-but-suboptimal params, so the contract-violation/silent-corruption risk is real but narrower than 'on early numerical break' suggests, warranting medium rather than high. Recommended fix (set out-params to NaN on failure, or add a bool success out-param/overload) is source-compatible/ABI-additive; only changing the UInt return type would be breaking.

### [MATH-3] EmgGradientDescent::estimateEmgParameters — const estimator writes diagnostic text to std::cout even when print_debug == 0
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/MATH/MISC/EmgGradientDescent.h · _math-spline_

```cpp
UInt estimateEmgParameters(...) const
```
- **Expectation:** With the `print_debug` parameter at its default 0 ("The level of debug information to print in the terminal"), a const compute method should be silent and not touch stdout.
- **Actual:** On any NaN/Inf in parameters or in the loss value, estimateEmgParameters writes multi-line messages to std::cout unconditionally — these prints are NOT guarded by `print_debug_`. So a quiet library call can suddenly spam stdout on numerically hard peaks regardless of the configured debug level.
- **Evidence:** EmgGradientDescent.cpp lines 625-628 (`std::cout << ... "One or more parameters are invalid." ...`) and lines 637-639 (`std::cout << ... "Bad: E value is invalid. current_E=" ...`) are outside any `if (print_debug_ ...)` guard, unlike the surrounding debug blocks.
- **Fix:** Gate these prints behind `print_debug_` (or route to OPENMS_LOG_*). ABI-safe; pure implementation change.

### [MATH-4] EmgGradientDescent::fitEMGPeakModel — left_pos/right_pos use 0.0 as a magic "unset" sentinel, so a legitimate position of 0.0 silently means "whole peak"
`medium` · `surprising-default` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/MISC/EmgGradientDescent.h · _math-spline_

```cpp
template <typename PeakContainerT> void fitEMGPeakModel(const PeakContainerT& input_peak, PeakContainerT& output_peak, const double left_pos = 0.0, const double right_pos = 0.0) const
```
- **Expectation:** The header says "If `left_pos` and `right_pos` are passed, then only that part of the peak is taken into consideration" and documents them as "RT or MZ value of the first/last point of interest". A caller passing left_pos=0.0 would expect the window to start at position 0.0.
- **Actual:** The implementation treats 0.0 as "not passed": `start_it = left_pos ? input_peak.PosBegin(left_pos) : input_peak.begin();` and `end_it = right_pos ? input_peak.PosEnd(right_pos) : input_peak.end();`. A real RT/position of exactly 0.0 is indistinguishable from the default and silently expands the window to the full peak.
- **Evidence:** EmgGradientDescent.cpp lines 737-738: `... start_it = left_pos ? input_peak.PosBegin(left_pos) : input_peak.begin(); ... end_it = right_pos ? input_peak.PosEnd(right_pos) : input_peak.end();`
- **Fix:** Document that 0.0 means "unbounded on that side", or switch the sentinel to NaN / std::optional<double> so position 0.0 is a valid bound. Adding an overload taking std::optional is ABI-additive; changing the signature is source-breaking.
- **Verifier correction:** The 0.0 magic-sentinel is real and undocumented, but the practical impact is narrower than "a legitimate position of 0.0 silently means whole peak" implies. For MSSpectrum (m/z) 0.0 is physically impossible, so no collision occurs. For MSChromatogram (RT), 0.0 is the lower bound of valid RT, so left_pos=0.0 yields essentially the same window as PosBegin(0.0); only a right_pos=0.0 (ending the window at the first point) triggers a real but pathological silent expansion to the full peak. The only in-tree caller passes genuine peak boundaries. Recommend documenting that 0.0 means "unbounded on that side", or adding an ABI-additive std::optional<double> overload; switching the existing parameter to NaN/optional in place would be source-breaking.

### [MATH-5] BSplineSmoothingSpline (constructor) / num_interior_knots() — Invalid/unsorted input does not throw and is not reported by k; failure is only observable via ok(), and num_interior_knots()/rss() return defaults that look like a valid degenerate fit
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/MISC/BSplineSmoothingSpline.h · _math-spline_

```cpp
BSplineSmoothingSpline(const std::vector<double>& x, const std::vector<double>& y, double s = -1.0, int k = 3); int num_interior_knots() const;
```
- **Expectation:** Sibling spline class CubicSpline2d throws Exception::IllegalArgument on unsorted x, mismatched sizes, or <2 points. A developer used to that convention would expect the same here, or at least that the accessors signal the failed state.
- **Actual:** The constructor logs to OPENMS_LOG_ERROR and silently returns with ok_=false on bad input (sizes, unsorted x). After a failed construct, num_interior_knots() returns 0 and rss() returns 0.0 — values that are indistinguishable from a perfect/degenerate successful fit. Only ok() reveals failure, and it is easy to skip given the other getters appear meaningful. The `k` (degree) argument is also stored but the BSpline path is always cubic, so non-3 values of k are silently ignored.
- **Evidence:** BSplineSmoothingSpline.cpp lines 28-45 (`ok_ = false; ... OPENMS_LOG_ERROR << ...; return;`), members default `num_interior_knots_ = 0; rss_ = 0.0;` (header lines 134-135). k_ is stored (header line 137) but fit_smoothing_spline/try_polynomial_fit hardcode cubic (cpp lines 265, 301-306).
- **Fix:** Either throw IllegalArgument like CubicSpline2d for consistency, or have the numeric getters throw/return NaN when !ok(). At minimum document that k!=3 is ignored. Returning NaN from getters when !ok_ is ABI-safe.
- **Verifier correction:** Two real defects, mis-stated in emphasis: (1) After a hard-failure construct, rss() returns 0.0 and num_interior_knots() returns 0 — values identical to those a SUCCESSFUL degenerate/polynomial fit produces (cpp:269, near-zero RSS), so they look like a valid fit; only ok() reveals failure. eval() is the exception: it correctly returns NaN when !ok_. So the trap is specifically in the two non-eval numeric getters, not all accessors. (2) The k (degree) argument is stored in k_ but never read; both fit paths are hardcoded cubic, so k!=3 is silently ignored and undocumented as such — concretely, MultipleTesting.cpp passes smooth_df into k expecting an effect. Recommend returning NaN (rss) / a sentinel like -1 (num_interior_knots) when !ok_, and documenting that k is currently ignored. The NaN-from-getters fix is ABI-safe (no signature change); throwing like CubicSpline2d would be source-breaking. Severity medium, not high: the primary prediction path eval() is guarded and the production caller checks ok().

### [MATH-8] Math::Histogram::maxValue / minValue — maxValue()/minValue() return the largest/smallest BIN COUNT, not the max/min data value or bin position
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/STATISTICS/Histogram.h · _math-stats_

```cpp
ValueType maxValue() const / ValueType minValue() const
```
- **Expectation:** On a class whose other accessors are minBound()/maxBound()/binSize()/centerOfBin(), a caller reading `histogram.maxValue()` would naturally expect the maximum value along the value axis (the largest measured/binned value), or at least the position of the tallest bin.
- **Actual:** It returns `*std::max_element(bins_.begin(), bins_.end())`, i.e. the highest *count* stored in any bin (a ValueType count), and minValue() the lowest count. It says nothing about which bin or what data value.
- **Evidence:** Lines 116-125: `ValueType maxValue() const { return *(std::max_element(bins_.begin(), bins_.end())); }` and `ValueType minValue() const { return *(std::min_element(bins_.begin(), bins_.end())); }`. The doc comments themselves say 'returns the bin with the highest count' / 'lowest count' - the NAME contradicts the doc.
- **Fix:** Add clearer additive accessors `maxBinCount()`/`minBinCount()` and deprecate maxValue()/minValue(). Renaming directly would break ABI; the additive aliases are source-compatible.
- **Verifier correction:** maxValue()/minValue() return the highest/lowest BIN COUNT (ValueType) across all bins via std::max_element/std::min_element — not a value-axis (BinSizeType) quantity and not the position/index of the tallest bin. The doc comments ("returns the bin with the highest count" / "lowest count") are themselves imprecise: the function returns the count itself, not the bin or its index. The misleading-name finding is genuine but should be graded medium, not high: both existing callers use it correctly as the peak/normalizing bin count, the inline doc partially signals intent, and misuse is recoverable. Recommended fix (additive maxBinCount()/minBinCount() aliases + [[deprecated]] on maxValue/minValue) is source-compatible; an outright rename would be ABI-breaking.

### [MATH-9] Math::quantile1st / quantile3rd / median — Read-only-looking quantile/median helpers SORT the caller's range in place by default (default sorted=false)
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/StatisticFunctions.h · _math-stats_

```cpp
static double quantile1st(IteratorType begin, IteratorType end, bool sorted = false)
```
- **Expectation:** A function called `median(begin,end)` / `quantile1st(begin,end)` that returns a statistic looks like a pure read of the range; a caller would not expect the underlying container to be permanently reordered.
- **Actual:** With the default `sorted = false`, each of these calls `std::sort(begin, end)` directly on the caller's iterators, mutating the source container as an undocumented-at-the-callsite side effect. Passing a const range would not compile, which masks this only partly.
- **Evidence:** median: line 140 `std::sort(begin, end);`; quantile1st: line 242 `std::sort(begin, end);`; quantile3rd: line 273 `std::sort(begin, end);`, each guarded only by `if (!sorted)` where sorted defaults to false.
- **Fix:** Document the in-place sort prominently in each brief (currently only the @param mentions it), and consider an overload taking a const range that copies internally. Behavior is long-standing; keep the default but make the side effect obvious. Additive overload is ABI-safe.

### [MATH-10] Math::quantile1st / quantile3rd — quantile1st/quantile3rd use a different (sub-median) definition than the sibling quantile(begin,end,q)
`medium` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/StatisticFunctions.h · _math-stats_

```cpp
static double quantile1st(...) / static double quantile3rd(...)
```
- **Expectation:** Within the same header, quantile1st(x) and quantile3rd(x) should agree with quantile(x, 0.25) and quantile(x, 0.75). A caller mixing the two APIs expects consistent Q1/Q3.
- **Actual:** quantile1st/quantile3rd compute the median of the lower/upper half of the data (Tukey 'hinges'-style), while quantile(begin,end,q) uses Type-7 linear interpolation. The two will return different Q1/Q3 on the same data, and SummaryStatistics uses the half-median variant while tukeyUpperFence uses the Type-7 variant - so 'lowerq' in one place and the IQR fence in another are not comparable.
- **Evidence:** quantile1st returns `median(begin, begin + (size/2)-1, true)` (line 248) vs quantile() line 325-333 doing `pos = q*(n-1); (1-frac)*x[i] + frac*x[i+1]`. SummaryStatistics line 945 uses quantile1st; tukeyUpperFence line 362-363 uses quantile(...,0.25/0.75).
- **Fix:** Document that quantile1st/3rd are hinge-style and differ from quantile(...,q); or reimplement them in terms of quantile() for one consistent definition. Reimplementation changes numeric output (source-compatible signature).

### [MATH-11] Math::BasicStatistics::update (single-range overload) — update(probabilities) computes mean/variance over implicit coordinates 0,1,2,... not over the probability values themselves
`medium` · `unit-or-index` · ABI: `none` · src/openms/include/OpenMS/MATH/STATISTICS/BasicStatistics.h · _math-stats_

```cpp
template <typename ProbabilityIterator> void update(ProbabilityIterator probability_begin, ProbabilityIterator probability_end)
```
- **Expectation:** A statistics class with mean()/variance()/sum() given a single range of numbers would, by least surprise, treat that range as the sample and report its mean/variance.
- **Actual:** The single-range update interprets the input as a probability mass function over the index positions: it computes sum=Σp_i, mean=Σ(p_i*i)/sum, variance=Σ(p_i*(i-mean)^2)/sum. So mean() is the weighted-average INDEX, not the average of the supplied numbers. This is a niche distribution convention not signalled by the generic name 'update'.
- **Evidence:** Lines 95-108: `sum_ += *iter; mean_ += *iter * pos; ... mean_ /= sum_;` with `pos` the running index; variance uses `RealType(pos) - mean_`. The two-range overload makes the coordinate explicit, confirming the first treats position as the coordinate.
- **Fix:** Document loudly that the single-range overload treats the input as weights over positions 0..n-1 (a histogram/pmf), not as raw samples; or rename to updateFromPMF(). Doc/additive rename is ABI-safe.
- **Verifier correction:** The single-range update(begin,end) interprets the input range as a probability mass / weight function over the implicit coordinates 0,1,2,...,n-1: it returns sum()=Σp_i, mean()=Σ(p_i·i)/Σp_i (a weighted-average INDEX, not the average of the input values), and variance() about that index-mean. This is intended (it is fed histograms in FeatureFinder/PoseClustering, and the explicit two-range overload supplies coordinates), and it is partially signposted by the ProbabilityIterator parameter name and the 'distribution' framing — so it is not fully silent. But the single-range overload's own doc comment does not state the position-as-coordinate semantics and mean()/variance() are undocumented, so the generic API names invite misuse by anyone treating it as a general-purpose sample-statistics class. Recommend documenting the single-range overload's PMF-over-positions semantics (and/or adding an additive updateFromPMF() alias). Doc/additive-only, ABI-safe.

### [MATH-15] Math::quantile (sorted-vector overload) — Free-function quantile(x, q) silently clamps out-of-range q instead of signaling, unlike the iterator overload which throws
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/MathFunctions.h · _math-stats_

```cpp
template<typename T1> typename T1::value_type quantile(const T1& x, double q)
```
- **Expectation:** Two functions named Math::quantile in the same library should agree on error handling for q outside [0,1].
- **Actual:** This overload silently clamps q to [0,1] (`if (q < 0.0) q = 0.; if (q > 1.0) q = 1.;`), whereas the StatisticFunctions.h Math::quantile(begin,end,q) throws Exception::InvalidValue for the same out-of-range q. A caller passing an erroneous q gets a silent endpoint here but an exception there.
- **Evidence:** MathFunctions.h lines 463-464 `if (q < 0.0) q = 0.; if (q > 1.0) q = 1.;` vs StatisticFunctions.h lines 318-322 `if (q < 0.0 || q > 1.0) { throw Exception::InvalidValue(...); }`.
- **Fix:** Make the two overloads consistent (both throw, or both clamp-and-document). Aligning error behavior is source-compatible at the signature level.
- **Verifier correction:** The two OpenMS::Math::quantile overloads disagree on out-of-range q handling: the MathFunctions.h container overload silently (and undocumentedly) clamps q to [0,1], while the StatisticFunctions.h iterator overload throws Exception::InvalidValue. Real but severity is medium, not high: the clamped result is a valid endpoint, so the failure is silent but not data-corrupting or crashing. Separately (beyond this claim), the two overloads also use different interpolation index formulas and can return different values even for valid q.

### [ML-2] SingleLinkage::operator() — threshold documented as silently ignored actually throws NotImplemented when < 1
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/ML/CLUSTERING/SingleLinkage.h · _ml-clustering-a_

```cpp
void operator()(DistanceMatrix<float> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const float threshold = 1) const override
```
- **Expectation:** Header doc says '@param threshold float value to meet Base class interface, will not be used'. The base ClusterFunctor and sibling linkages (Average/Complete) accept any threshold and use it to stop merging. A caller porting from CompleteLinkage to SingleLinkage with e.g. threshold=0.3 would expect it to be ignored (per the doc), getting a full clustering.
- **Actual:** The implementation throws Exception::NotImplemented for any threshold < 1, so passing a real threshold (which is the whole point of the base-class signature) aborts with an exception rather than being ignored as documented.
- **Evidence:** SingleLinkage.cpp:47-51 `if (threshold < 1) { OPENMS_LOG_ERROR << ...; throw Exception::NotImplemented(...); }`. Header SingleLinkage.h:51 `@param threshold ... will not be used`.
- **Fix:** Reconcile doc and behavior: either implement threshold support, or document that threshold<1 throws NotImplemented (not 'will not be used'). Additive/doc-only fix is ABI-safe; the throw itself is already shipped behavior.

### [ML-4] ClusterAnalyzer::cohesion — 'cohesion' returns average intra-cluster DISTANCE (lower = tighter), opposite of expected cohesion score
`medium` · `unit-or-index` · ABI: `none` · src/openms/include/OpenMS/ML/CLUSTERING/ClusterAnalyzer.h · _ml-clustering-a_

```cpp
std::vector<float> cohesion(const std::vector<std::vector<Size> > & clusters, const DistanceMatrix<float> & original)
```
- **Expectation:** A method named 'cohesion' is expected to return a tightness score where a HIGHER value means a more cohesive (tighter) cluster, and a perfectly tight singleton would score 0 distance / maximal cohesion.
- **Actual:** It returns the average pairwise DISTANCE within each cluster (lower value = tighter cluster, i.e. inversely related to cohesion). Moreover, for a singleton cluster it returns the GLOBAL average distance over the whole matrix, not 0, which is surprising for a per-cluster cohesion measure.
- **Evidence:** ClusterAnalyzer.cpp:667-682: accumulates `original.getValue(...)` over intra-cluster pairs into `av_c_dist`, divides by pair count, and `if (clusters[i].size() == 1) { av_c_dist = av_dist; }` (global average). Header ClusterAnalyzer.h:59 only says 'calculate the cohesions of a certain partition'.
- **Fix:** Document explicitly that 'cohesion' here is the mean intra-cluster distance (smaller = more cohesive) and that singletons are assigned the global mean distance. Renaming (e.g. meanIntraClusterDistance) would be the ideal fix but is source-breaking; prefer the doc clarification for ABI stability.
- **Verifier correction:** The claim is accurate. cohesion() returns the mean intra-cluster distance (lower = tighter), opposite to the conventional cohesion direction, and for singleton clusters it returns the global mean pairwise distance (not 0) to avoid a 0/0 NaN. Severity adjusted to medium rather than high: results are inspectable/recoverable and in a sane range (no crash or data loss), but the name-vs-value-direction inversion and undocumented singleton override can lead to silently inverted cluster rankings if a caller assumes higher=better. Recommendation stands: document that this is mean intra-cluster distance (smaller = more cohesive) and that singletons receive the global mean distance; a rename would be the ideal fix but is source-breaking, so prefer doc clarification (ABI-safe).

### [ML-5] ClusterAnalyzer::cut(const Size, const std::vector<BinaryTreeNode>&, std::vector<std::vector<Size>>&) — cut(...Size clusters) appends to the output without clearing, unlike the sibling cut(...subtrees) overload
`medium` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/ML/CLUSTERING/ClusterAnalyzer.h · _ml-clustering-a_

```cpp
void cut(const Size cluster_quantity, const std::vector<BinaryTreeNode> & tree, std::vector<std::vector<Size> > & clusters)
```
- **Expectation:** An out-parameter that is 'filled corresponding to the given cluster_quantity' (per the doc) is expected to be reset first, consistent with the sibling overload. A caller reusing a clusters vector across calls expects it replaced, not appended.
- **Actual:** The Size-vector overload never clears `clusters`; it push_back()s new clusters onto whatever was already there, then calls clusters.resize(cluster_quantity). If the caller passes a non-empty vector, pre-existing entries are kept, the result is corrupted, and resize() may even truncate the freshly-added clusters. The subtree overload, by contrast, explicitly does subtrees.clear().
- **Evidence:** ClusterAnalyzer.cpp:501-557 cut(...Size) has no clear(); line 547 `clusters.push_back(actCluster);`, line 556 `clusters.resize(cluster_quantity);`. Versus ClusterAnalyzer.cpp:570 cut(...subtrees) `subtrees.clear();`. Header doc ClusterAnalyzer.h:91 'the argument clusters is filled'.
- **Fix:** Add `clusters.clear();` at the start of the Size-vector cut() to match the subtree overload, or document that clusters must be passed empty. Adding the clear() is a behavioral fix but ABI-stable (signature unchanged).
- **Verifier correction:** Claim is mechanically accurate. Two corrections: (1) severity is medium, not high — the only in-repo caller (SpectraMerger.h:261) passes a freshly-declared empty vector and is unaffected, and the misuse is silent-but-recoverable rather than a guaranteed real-world failure. (2) abi_impact is none, not just "ABI-stable" — adding clusters.clear() changes behavior only; the signature and exported symbol are unchanged, making it both source- and binary-compatible. Note also the corroborating evidence the claim omitted: ClusterAnalyzer_test.cpp:256 itself calls clusters.clear() before reusing the vector, confirming the author was aware the function does not clear internally.

### [ML-8] ClusterAnalyzer::averageSilhouetteWidth — averageSilhouetteWidth divides each step's sum by the fixed total element count, not by the number of contributing points
`medium` · `return-value` · ABI: `none` · src/openms/include/OpenMS/ML/CLUSTERING/ClusterAnalyzer.h · _ml-clustering-a_

```cpp
std::vector<float> averageSilhouetteWidth(const std::vector<BinaryTreeNode> & tree, const DistanceMatrix<float> & original)
```
- **Expectation:** The standard average silhouette width for a clustering step is the mean of per-point silhouette values s(i) in [-1,1] over the points that participate (points in non-singleton clusters), yielding a value in [-1,1].
- **Actual:** For every step the accumulated silhouette sum (which skips singleton clusters and points with interdist 0) is divided by the constant (tree.size()+1) = total number of leaf elements, regardless of how many points actually contributed. Early steps where most points sit in singletons are therefore systematically deflated, so the returned 'average silhouette width' is not the conventional normalized average a caller expects.
- **Evidence:** ClusterAnalyzer.cpp:282-298 only adds silhouettes for clusters with size>1 and interdist!=0, but line 306 `average_silhouette_widths.push_back(average_overall_silhouette / (float)(tree.size() + 1));` divides by the fixed total element count.
- **Fix:** Either document that the per-step value is normalized by total element count (not by contributing points) or divide by the count of contributing points. Behavior change would alter outputs/tests, so confirm intended semantics before changing; at minimum clarify in the header.

### [ML-9] ClusteringGrid::isNonEmptyCell — isNonEmptyCell documents two exceptions it never throws
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ML/CLUSTERING/ClusteringGrid.h · _ml-clustering-b_

```cpp
bool isNonEmptyCell(const CellIndex &cell_index) const
```
- **Expectation:** Per the header Doxygen, isNonEmptyCell '@throw Exception::IllegalArgument if the coordinates (x,y) lie outside the grid' and '@throw Exception::InvalidValue if one of the two indices is negative.' A caller would write defensive code expecting an exception for out-of-grid / negative indices.
- **Actual:** The implementation is simply `return cells_.contains(cell_index);` — it never validates the index and never throws. A negative or out-of-grid cell index just returns false (treated as 'empty'), silently. The documented contract is entirely fictional.
- **Evidence:** ClusteringGrid.cpp:86-89 `bool ClusteringGrid::isNonEmptyCell(const CellIndex &cell_index) const { return cells_.contains(cell_index); }` vs header lines 104-106 `@throw Exception::IllegalArgument ... @throw Exception::InvalidValue ...`
- **Fix:** Fix the documentation to match behavior (remove the @throw clauses) — this is a doc-only, ABI-safe change. If validation is actually desired, add it in a new method or a clearly versioned change; do not change the existing return semantics silently.
- **Verifier correction:** isNonEmptyCell never throws either documented exception; it returns cells_.contains(cell_index). The @throw clauses are fictional and the @throw IllegalArgument clause is additionally mismatched — it describes a Point/coordinate check (copied from getIndex()) though the parameter is a CellIndex. Severity is medium (not high): existing callers rely on and require the non-throwing false-return, so there is no data loss/crash in real pipelines; the surprise only bites a defensive external caller, and it returns a benign 'false'. Recommended fix is doc-only (remove the two @throw lines), ABI-safe.

### [ML-12] GridBasedCluster::operator< / operator> / operator== — Cluster comparison/equality use only the y-coordinate of the centre, not the cluster
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/ML/CLUSTERING/GridBasedCluster.h · _ml-clustering-b_

```cpp
bool operator<(const GridBasedCluster& other) const; bool operator>(...) const; bool operator==(...) const;
```
- **Expectation:** operator== on a value type named GridBasedCluster reads as 'are these the same cluster' (same points/centre/box). operator< reads as a total ordering of clusters. A developer using clusters as keys/in a set would assume identity is by content.
- **Actual:** All three operators compare ONLY `centre_.getY()`. Two clusters with completely different points, x-coordinates and bounding boxes compare equal whenever their centre y is identical, and the ordering ignores x entirely. This actively breaks `std::map<GridBasedCluster,int> index_list` in extendClustersY: clusters sharing a centre-y collide as the same key, so `index_list.find(*c1)->second` can resolve to the WRONG original cluster index.
- **Evidence:** GridBasedCluster.cpp:50-63 `operator<` returns `centre_.getY() < other.centre_.getY();`, `operator==` returns `centre_.getY() == other.centre_.getY();`. Misused as map key at GridBasedClustering.h:303-307 `std::map<GridBasedCluster, int> index_list; ... index_list.insert(std::make_pair(..., *it));`.
- **Fix:** Document explicitly that these operators order by centre-y only (they exist to sort clusters vertically), and rename or supply a content-equality predicate; replace the GridBasedCluster-keyed map in extendClustersY with an index-keyed structure. The operator-as-y-sort is by design but the unqualified operator== is the dangerous part — at minimum compare the full centre. Changing the comparison is source-compatible but behavior-changing, so gate it carefully.
- **Verifier correction:** All three operators (GridBasedCluster.cpp:50-63) compare only centre_.getY(); this is undocumented anywhere in the class. operator== therefore is NOT content/identity equality despite the value-type name, and operator< is only a vertical (y) ordering, not a total order over clusters. The y-only operator< is plausibly by design for the "sort clusters in ascending y" step in extendClustersY, but the unqualified operator== is the dangerous part. It is used implicitly as the key comparator of std::map<GridBasedCluster,int> index_list (GridBasedClustering.h:303,307); when two distinct clusters in one x-cell share an identical centre-y the map collapses them to a single key and index_list.find(*c2)->second can return the wrong original index (h:373-384). Severity is medium (not high): the collision requires exact double equality of centre-y, so wrong results are possible but not common. Recommended fix is source-compatible (no signature change): document the y-only semantics, supply a content-equality predicate / make operator== compare the full centre (or full content), and replace the GridBasedCluster-keyed map with an index-keyed structure. Changing operator< / operator== bodies is behavior-changing but not ABI-breaking.

### [ML-15] ClusteringGrid::getGridSpacingX / getGridSpacingY — 'GridSpacing' is actually a vector of cell-boundary coordinates, not a spacing/step
`medium` · `misleading-name` · ABI: `breaking` · src/openms/include/OpenMS/ML/CLUSTERING/ClusteringGrid.h · _ml-clustering-b_

```cpp
std::vector<double> getGridSpacingX() const
```
- **Expectation:** A name like getGridSpacingX / grid_spacing_x strongly implies a spacing (delta/step between grid lines), typically a scalar or a vector of widths.
- **Actual:** The vector holds absolute boundary coordinates: the constructor sets range_x_ = (front, back) and getIndex does `lower_bound(grid_spacing_x_, position.getX(), less_equal)` — i.e. the entries are the cell-edge positions, not deltas. The min/max of the grid are the first/last entries, confirming these are coordinates. The 'spacing' name misleads about units and meaning.
- **Evidence:** ClusteringGrid.cpp:19 `range_x_(grid_spacing_x.front(),grid_spacing_x.back())`; ClusteringGrid.cpp:80 `std::lower_bound(grid_spacing_x_.begin(), grid_spacing_x_.end(), position.getX(), std::less_equal<double>())`.
- **Fix:** Document (and ideally rename in future) these as grid cell boundary coordinates, not spacings. Adding clarifying Doxygen is ABI-safe; renaming the public methods/params would be breaking.

### [ML-20] ROCCurve::AUC — AUC() returns 0.5 and prints to cerr instead of signaling an unusable dataset
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/source/ML/ROCCURVE/ROCCurve.cpp · _ml-interp-solvers_

```cpp
double AUC()
```
- **Expectation:** A method named AUC() either returns a valid area in [0,1] or signals failure (throw / sentinel) when the data cannot produce a meaningful curve.
- **Actual:** On an empty pair set it writes a message to std::cerr and returns 0.5 -- a value that is itself a legal AUC (random classifier). A caller thresholding on AUC >= 0.5 silently treats a degenerate/empty input as 'random but valid'. Note the cerr text even says 'no positives or no negatives' but the guard only checks score_clas_pairs_.empty(), so a set with positives-only or negatives-only is NOT caught here; it instead falls through to 'area /= truePos * falsePos' which divides by zero.
- **Evidence:** if (score_clas_pairs_.empty()) { std::cerr << "ROCCurve::AUC() : unsuitable dataset (no positives or no negatives)\n"; return 0.5; } ... area /= truePos * falsePos;
- **Fix:** Return a NaN or throw on an unusable set, or at minimum guard the truePos==0 || falsePos==0 cases (the only ones the cerr text claims to detect) instead of dividing by zero. Additive/source-compatible: keep the signature, change the degenerate return to std::numeric_limits<double>::quiet_NaN() and document it; or throw InvalidValue.
- **Verifier correction:** AUC() does have a genuine silent-failure smell: empty input returns the legal value 0.5 (after a cerr message), and single-class (positives-only / negatives-only) input is NOT caught by the .empty() guard despite the cerr text saying so — it instead reaches `area /= truePos * falsePos` (UInt*UInt == 0) and divides a double by zero, producing NaN/inf. However, severity is medium, not high: the class is explicitly documented as "[This class is buggy and usage is discouraged!]" (ROCCurve.h:24) and has no callers anywhere in the codebase, and the worst path (divide-by-zero) is partially loud via cerr and a propagating NaN. Fix is source-compatible: keep `double AUC()`, guard `truePos == 0 || falsePos == 0` and return std::numeric_limits<double>::quiet_NaN() (or throw Exception::InvalidValue) for degenerate sets, and make the empty case consistent with the message rather than returning the legal 0.5.

### [ML-23] NonNegativeLeastSquaresSolver::solve — Non-const Matrix/vector solve overload clobbers A and b; the 'solve' name and overload set do not warn at the call site
`medium` · `surprising-default` · ABI: `breaking` · src/openms/include/OpenMS/ML/NNLS/NonNegativeLeastSquaresSolver.h · _ml-interp-solvers_

```cpp
static Int solve(Matrix<double>& A, std::vector<double>& b, std::vector<double>& x)
```
- **Expectation:** solve(A, b, x) reads A and b to produce x; a competent caller reusing A/b afterwards expects them intact, especially because a sibling overload 'solve(const Matrix&, const Matrix&, Matrix&)' is explicitly non-destructive.
- **Actual:** The in-place overloads overwrite A and b with NNLS workspace (Q*A, Q*b). This is documented in the doxygen ('contents are unspecified on return'), but the only signature-level signal is the non-const reference; the overload that copies and the one that destroys share the same name 'solve', so an accidental overload resolution (e.g. passing a non-const Matrix) silently destroys the caller's matrix.
- **Evidence:** Doc: 'After the call the contents of @p A and @p b are clobbered by the NNLS workspace and should not be read again.' vs the const overload 'Solve ... without modifying the inputs.'
- **Fix:** Consider renaming the destructive overloads (e.g. solveInPlace) so the call site reflects the side effect; at minimum keep the strong doc. Renaming is source-compatible-ish but ABI-breaking; the conservative fix is documentation already present plus a note in the class brief.
- **Verifier correction:** The destructive overload is selected by passing std::vector-typed b and x (not merely a non-const Matrix A); the non-destructive copy overload takes Matrix-typed b and x. The shared-name hazard and silent clobbering of A and b are real and confirmed by a defensive-copy comment at an actual call site (IsobaricWorkflow.cpp:649). Severity is medium, not high: the behavior is clearly documented at the declaration, all current call sites are correct, and misuse yields silently wrong numeric results but is recoverable. The recommended remediation (rename to solveInPlace) is ABI-breaking; the documentation-only mitigation already present is ABI-none.

### [ML-25] GridSearch::evaluate — evaluate() maximizes the evaluator and uses startValue as a hard lower bound, opposite to the loss-minimizing CrossValidation sibling
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/ML/GRIDSEARCH/GridSearch.h · _ml-interp-solvers_

```cpp
auto evaluate(Functor evaluator, ... startValue, std::array<size_t,N>& resultIndices) const
```
- **Expectation:** Within the ML module a grid 'evaluate' over candidate parameters would plausibly minimize an error (as CrossValidation::gridSearch1D does: 'lower is better'); and a 'startValue' might read as a starting point.
- **Actual:** GridSearch keeps the combination with the largest evaluator result ('if (currVal > bestValue)'), i.e. it MAXIMIZES, and startValue acts as a floor: any combination scoring <= startValue is never selected and the returned bestIndices stay default-initialized (all zero) -- a silent 'no result' that looks like the first grid point. This is the opposite optimization direction from the sibling CrossValidation::gridSearch1D in the same ML cluster.
- **Evidence:** Looper: 'if (currVal > bestValue) { bestValue = currVal; bestIndices[param_index] = index; }'. CrossValidation doc: 'score(abs_errs) (lower is better)'.
- **Fix:** Document that evaluate() returns the argmax and that startValue is an exclusive lower bound (use -inf to always pick a maximum). Header-only template; doc-only fix, source/ABI compatible.
- **Verifier correction:** Severity tightened from the implied-high framing to medium: the maximization direction is genuinely undocumented and opposite to the CrossValidation sibling, and `startValue` acts as an EXCLUSIVE lower bound whose failure mode silently returns the caller's default indices (looking like grid point 0). However, the only in-tree caller (BayesianProteinInferenceAlgorithm, non-negative AUC evaluator, startValue=-1.0) is correct, no current code is broken, and the fix is doc-only — hence medium (invites silently-wrong reuse, but recoverable) rather than high. Recommended fix: document that evaluate() returns the argmax (maximizes), that `startValue` is an exclusive lower bound (pass -infinity to always select the maximum), and that resultIndices are left untouched if no combination strictly exceeds startValue. Header-only template; doc-only, source- and ABI-compatible.

### [ML-16] RansacModelQuadratic::rm_rsq_impl — rm_rsq() on the quadratic model returns chi-squared (RSS), not R-squared
`medium` · `misleading-name` · ABI: `none` · src/openms/source/ML/RANSAC/RANSACModelQuadratic.cpp · _ml-regression-ransac_

```cpp
static double RansacModelQuadratic::rm_rsq_impl(const DVecIt& begin, const DVecIt& end)
```
- **Expectation:** RansacModel::rm_rsq() is documented in the base class as 'Returns the R-squared of the data applied to the model'. A caller (or a future RANSAC loop comparing models) expects a goodness-of-fit value in [0,1] where HIGHER is better, consistent with RansacModelLinear::rm_rsq_impl which returns lin_reg.getRSquared().
- **Actual:** The quadratic override returns quad_reg.getChiSquared() (line 50: 'return quad_reg.getChiSquared();'), i.e. the weighted residual sum of squares. This is an unbounded quantity where LOWER is better -- the opposite direction and a different scale from R-squared. The two sibling models with the identical public method name thus return contradictory quantities.
- **Evidence:** RANSACModelQuadratic.cpp:37-51 'double RansacModelQuadratic::rm_rsq_impl(...) { ... QuadraticRegression quad_reg; quad_reg.computeRegression(...); return quad_reg.getChiSquared(); }' vs RANSACModelLinear.cpp:35-49 returning 'lin_reg.getRSquared();'. Base-class contract RANSACModel.h:47-59 '@brief Returns the R-squared of the data applied to the model'.
- **Fix:** Make the quadratic model actually return an R-squared (compute 1 - RSS/TSS from the points) to honor the documented contract, OR rename the base-class method and doc to a neutral 'goodness'/score and document the per-model semantics. Given OpenMS ABI policy: the additive/safe fix is to correct rm_rsq_impl to return a true R-squared (behavior fix, header signature unchanged). If a behavior change is too risky, at minimum fix the base-class Doxygen to stop promising R-squared. The current implementation is only referenced from dead DEBUG_RANSAC code, so a behavior fix has minimal blast radius.

### [ML-17] RANSAC::RANSAC — RANSAC default constructor seeds the RNG from wall-clock time, making results non-deterministic by default
`medium` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/ML/RANSAC/RANSAC.h · _ml-regression-ransac_

```cpp
explicit RANSAC(uint64_t seed = time(nullptr))
```
- **Expectation:** RANSAC is a randomized algorithm, but a default-constructed instance (the common call 'RANSAC<...> r;') is most naturally expected to be reproducible, or at least the non-determinism to be opt-in. A developer reading 'RANSAC r;' would not expect the outlier set to silently change between two runs of the same program on the same data.
- **Actual:** The default argument is 'seed = time(nullptr)', so every default-constructed RANSAC seeds from the current second. Two runs (or two instances constructed in the same second vs. different seconds) can return different inlier sets for identical input, with no indication at the call site that a time-based seed was used.
- **Evidence:** RANSAC.h:73 'explicit RANSAC(uint64_t seed = time(nullptr)): shuffler_(seed) {}'.
- **Fix:** Prefer a fixed default seed (e.g. 0 or a documented constant) for reproducibility, requiring callers to opt into time-seeding via setSeed(time(nullptr)); this is source/ABI-compatible (only the default value changes). At minimum, document prominently in the class/constructor Doxygen that the default seed is wall-clock-based and results are non-deterministic unless setSeed() is called.
- **Verifier correction:** RANSAC's default constructor seeds its RNG from wall-clock time (`seed = time(nullptr)`), so the common `RANSAC<...> r;` is non-deterministic by default — confirmed at RANSAC.h:73 and exercised in production by MRMRTNormalizer::removeOutliersRANSAC (no setSeed call). Severity is medium (not high): output remains a valid RANSAC fit and is recoverable via setSeed(), but reproducibility is silently lost, contrary to both std::mt19937's fixed-seed convention and OpenMS's own practice of fixing seeds for determinism. Fix (changing the default to a fixed constant) is ABI-safe and source-compatible since default-argument values are not part of the mangled signature.

### [ML-18] LinearRegressionWithoutIntercept::getSlope — getSlope() silently returns NaN when fewer than 2 observations were added
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ML/REGRESSION/LinearRegressionWithoutIntercept.h · _ml-regression-ransac_

```cpp
double getSlope() const
```
- **Expectation:** getSlope() is a plain const getter named like LinearRegression::getSlope(); a caller expects either a usable slope or a clear error (exception) when the fit is undefined. Nothing in the name or signature warns that the return may be a sentinel that must be checked.
- **Actual:** When fewer than two data points have been added, getSlope() returns std::numeric_limits<double>::quiet_NaN() instead of signaling failure. A caller who does not explicitly std::isnan-check the result will propagate NaN into downstream calculations (e.g. RT alignment) silently. Note also the divisor sum_xx_ can be 0 even with n_>=2 (all x==0), yielding inf/NaN with no guard.
- **Evidence:** LinearRegressionWithoutIntercept.cpp:39-47 'if (n_ < 2) { return std::numeric_limits<double>::quiet_NaN(); } return sum_xy_ / sum_xx_;'. The header doc (line 48-50) only says 'returns the slope of the estimated regression line' with no mention of NaN.
- **Fix:** Document the NaN sentinel in the getSlope() Doxygen and/or add a 'bool isValid() const' (or 'size_t getN() const') accessor so callers can test before reading; an additive accessor is fully ABI-compatible. Ideally throw Exception::UnableToFit when n_<2 or sum_xx_==0, but that is a behavior change -- flag it as the ideal fix.
- **Verifier correction:** getSlope() does return quiet_NaN() for n_<2 (and is unguarded against sum_xx_==0), and this NaN sentinel is undocumented in the public OPENMS_DLLAPI header -- a real, surprising API contract. But the asserted impact ("propagates NaN into downstream RT alignment silently") is not realized anywhere: the sole caller (FeatureFinderMultiplexAlgorithm.cpp) already guards with intensities1.size() > 5, so the sentinel path is unreachable there and no RT-alignment caller exists. Severity is medium (latent hazard inviting future misuse, recoverable), not high. The additive fix (document the sentinel; add bool isValid()/size_t getN()) is ABI-compatible (abi_impact none); throwing UnableToFit would be a behavior change but still not an ABI break.

### [ML-19] LinearRegression::computeRegression / computeRegressionWeighted — On a failed/singular fit, computeRegression mutates slope_/intercept_/chi_squared_ before (or without) throwing, leaving the object in a misleading state
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/source/ML/REGRESSION/LinearRegression.cpp · _ml-regression-ransac_

```cpp
void computeRegression(double, it x_begin, it x_end, it y_begin, bool compute_goodness=true)
```
- **Expectation:** When fitting throws Exception::UnableToFit, a caller expects either the object to be left unchanged or all members to be clearly invalid. After catching the exception they would not expect getSlope()/getIntercept()/getChiSquared() to return values computed from a partial/failed solve.
- **Actual:** computeRegression() assigns slope_, intercept_ and chi_squared_ from the gte fit BEFORE checking 'pass', then throws if !pass -- so a caught exception leaves stale fitted members readable. In computeRegressionWeighted(), when the linear system is singular ('nonsingular' false) slope_/intercept_ are NOT updated (retain prior/default 0) yet chi_squared_ IS overwritten via computeWeightedChiSquare() using those stale coefficients, and only then is UnableToFit thrown -- an inconsistent partial mutation.
- **Evidence:** LinearRegression.cpp:220-228 sets slope_/intercept_/chi_squared_ then 'if (!pass){ throw Exception::UnableToFit(...); }'. LinearRegression.cpp:285-306: 'if (nonsingular){ slope_=X[0]; intercept_=X[1]; } chi_squared_ = computeWeightedChiSquare(..., slope_, intercept_); if (nonsingular){...} else { throw ... }'.
- **Fix:** Compute into locals and only commit to members after the success check; throw before mutating any member. This is an internal-implementation fix with no header/ABI change. At minimum, reset members to a sentinel on failure so a caller that swallows the exception cannot misread stale fit values.

### [COMP-2] BinnedSharedPeakCount::precursor_mass_tolerance_ — Dead, uninitialized 'precursor_mass_tolerance_' member suggests a tolerance feature that does not exist
`low` · `misleading-name` · ABI: `breaking` · src/openms/include/OpenMS/COMPARISON/BinnedSharedPeakCount.h · _comp-binned_

```cpp
double precursor_mass_tolerance_;
```
- **Expectation:** A protected member named precursor_mass_tolerance_ implies the functor uses a precursor mass tolerance (e.g. only compares spectra whose precursors agree within a tolerance, as sibling COMPARISON class PeakAlignment does via its 'precursor_mass_tolerance' param), and that it is settable/initialized.
- **Actual:** The member is never read, never written, and never initialized in any of the three Binned* subclasses. updateMembers_() is empty, no 'precursor_mass_tolerance' parameter is registered (defaultsToParam_ has no such default), and operator() ignores precursor masses entirely. The same dead member is copy-pasted into all three subclasses. Its value is indeterminate after construction.
- **Evidence:** Declared in BinnedSharedPeakCount.h:64, BinnedSpectralContrastAngle.h:62, BinnedSumAgreeingIntensities.h:68. grep for precursor_mass_tolerance_ in the three .cpp files returns nothing (EXIT 1). updateMembers_() bodies are empty (e.g. BinnedSharedPeakCount.cpp:46-48). By contrast PeakAlignment.cpp:27/65 actually registers and uses a 'precursor_mass_tolerance' parameter.
- **Fix:** Remove the unused member (binary-incompatible: changes class layout). If ABI must be preserved short-term, at minimum initialize it (= 0.0) and document that it is unused/reserved; but the clean fix is deletion in the next ABI-breaking window since it is protected state that misleads subclassers and readers.
- **Verifier correction:** The facts are accurate, but severity should be low, not high. The member precursor_mass_tolerance_ is dead (never read/written/initialized) and misleadingly named, copy-pasted into all three Binned* functors, with an indeterminate value after construction. Because it is never used by operator() or any other method, it causes NO wrong results, data loss, or crashes — only reader/subclasser confusion (a tolerance feature is implied but absent). Clean fix is deletion (ABI-breaking layout change); short-term mitigation is value-initializing it (= 0.0) and documenting it as unused/reserved.

### [COMP-5] BinnedSpectrumCompareFunctor (base class doc) — Base-class doc promises a 'normalized' parameter that none of the concrete functors actually register
`low` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/COMPARISON/BinnedSpectrumCompareFunctor.h · _comp-binned_

```cpp
Functors normalized in the range [0,1] are identifiable at the set "normalized" parameter of the ParameterHandler
```
- **Expectation:** Per the base class documentation, a caller could inspect param_/getParameters() for a boolean 'normalized' parameter to discover whether a given functor's output is in [0,1].
- **Actual:** None of the three concrete subclasses register a 'normalized' parameter. defaultsToParam_() is called but no setValue("normalized", ...) exists anywhere in the module, and updateMembers_() is empty in every subclass. A caller querying the 'normalized' parameter will not find it.
- **Evidence:** Header BinnedSpectrumCompareFunctor.h:26: 'Functors normalized in the range [0,1] are identifiable at the set "normalized" parameter of the ParameterHandler'. grep for 'normalized'/'setValue' in the four COMPARISON Binned* .cpp files returns no parameter registration (only code comments 'score normalized to interval [0,1]').
- **Fix:** Either register the documented 'normalized' boolean parameter in each subclass's defaults (additive, source-compatible), or remove the misleading sentence from the base-class documentation (doc-only). Prefer aligning code with doc since the claim is part of the public contract.
- **Verifier correction:** The claim is accurate but its severity framing ("silently wrong results") overstates impact. The missing 'normalized' parameter is an introspection/metadata flag; its absence misleads a caller about whether output is in [0,1] but does NOT corrupt the actual similarity scores (which are correctly normalized to [0,1] in all three Binned functors). The documented discovery mechanism is real and used by PeakAlignment (PeakSpectrumCompareFunctor family) via defaults_.setValue("normalized", 1, ...), so the contract is established, not invented — the Binned* family simply fails to register it. Recommended fix (add the boolean default in each subclass, or drop the sentence) is additive/source-compatible.

### [COMP-10] SteinScottImproveScore::operator()(const PeakSpectrum&, const PeakSpectrum&) — Similarity can be negative and is silently floored to 0 below a threshold
`low` · `return-value` · ABI: `none` · src/openms/source/COMPARISON/SteinScottImproveScore.cpp · _comp-spectra_

```cpp
double operator()(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const
```
- **Expectation:** Base class documents 'The value should be greater equal 0' for a similarity functor, so a caller expects a non-negative similarity; and a 'similarity' is expected to reflect overlap, not be clamped to zero by an internal threshold without notice.
- **Actual:** The score is `(sum - z) / sqrt(sum1*sum2)` where `z = (tolerance/10000)*(sum3*sum4)` is subtracted, so the result can be negative for dissimilar spectra. Additionally any score below the 'threshold' parameter (default 0.2) is silently set to 0, so a caller cannot distinguish 'genuinely zero overlap' from 'small-but-nonzero, thresholded away'.
- **Evidence:** `score = (sum - z) / (std::sqrt((sum1 * sum2))); if (score < (float)param_.getValue("threshold")) { score = 0; }` with `constant = epsilon / 10000;` and `z = constant * (sum3 * sum4);`
- **Fix:** Document that the returned value may be negative before thresholding and that sub-threshold scores are floored to 0; consider exposing whether thresholding occurred. Doc-only fix; ABI-neutral.
- **Verifier correction:** The raw similarity score (sum - z)/sqrt(sum1*sum2) can be negative for dissimilar spectra, but with the DEFAULT threshold (0.2 > 0) all negative scores are floored to 0, so the base-class ">= 0" contract holds by default. A negative value is returned only if a caller deliberately sets threshold to a negative value. The thresholding itself is already documented in the "threshold" parameter description. The remaining issue is purely a documentation gap: the operator()/class Doxygen should state that the pre-threshold score may be negative and that the threshold parameter also serves as a configurable non-negativity floor whose lowering below 0 can produce negative returns, breaking the base-class >=0 promise. Doc-only, ABI-neutral.

### [COMP-13] SpectrumCheapDPCorr::getPeakMap / lastconsensus — getPeakMap returns a copy by value while lastconsensus returns a const reference
`low` · `api-consistency` · ABI: `none` · src/openms/include/OpenMS/COMPARISON/SpectrumCheapDPCorr.h · _comp-spectra_

```cpp
std::map<UInt, UInt> getPeakMap() const; const PeakSpectrum & lastconsensus() const
```
- **Expectation:** Sibling accessors that both expose results of the last comparison would be expected to follow a consistent convention (both by const-ref, or both by value).
- **Actual:** `lastconsensus()` returns `const PeakSpectrum&` (reference into the mutable member) while `getPeakMap()` returns `std::map<UInt,UInt>` by value (a copy). The inconsistency is easy to miss and the returned reference from lastconsensus() is invalidated by the next operator() call.
- **Evidence:** `const PeakSpectrum & lastconsensus() const;` vs `std::map<UInt, UInt> getPeakMap() const;`
- **Fix:** Make the convention consistent and document the lifetime of lastconsensus()'s reference (invalidated on next call). Changing return types is source-compatible/breaking depending on direction; documenting is ABI-neutral.
- **Verifier correction:** The claim is factually correct that lastconsensus() returns const PeakSpectrum& and getPeakMap() returns std::map<UInt,UInt> by value, and that operator() resets both members (so the reference is invalidated on the next call). But this is a minor API-consistency/style issue, not a correctness hazard: const-ref-into-member is a standard accessor idiom and the reference-lifetime caveat is the ordinary C++ rule, not a hidden trap. Severity should be low, category api-consistency. The only ABI-neutral fix (documenting lastconsensus()'s lifetime) has no ABI impact; harmonizing return types would be a source/ABI break but is not required to address the real (mild) surprise.

### [MATH-17] GaussFitter::fit — fit() takes a non-const vector& documented as [in,out] but never modifies the input points
`low` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/STATISTICS/GaussFitter.h · _math-fitters_

```cpp
GaussFitResult fit(std::vector<DPosition<2> >& points) const
```
- **Expectation:** A `std::vector&` parameter annotated `@param[in,out] points` should be mutated by the call (sorted, rescaled, filtered, etc.); a caller would defensively copy it.
- **Actual:** The Levenberg-Marquardt functor stores only a const pointer to the data and reads it via const_iterator; nothing writes back to `input`. The argument is effectively input-only, yet the signature forces a mutable lvalue and the doc says [in,out].
- **Evidence:** GaussFitter.cpp:87-116 `fit(vector<DPosition<2>>& input)`: GaussFunctor holds `const std::vector<DPosition<2>>* m_data;` and only reads `it->getX()/getY()`. No write to input. Header line 91 doc: `@param[in,out] points`. Sibling GammaDistributionFitter::fit takes `const std::vector<DPosition<2>>&` (GammaDistributionFitter.h:71).
- **Fix:** Change the parameter to `const std::vector<DPosition<2>>&` to match GammaDistributionFitter and the actual behavior (source-compatible for most callers but technically an overload/signature change; can also just add a const-ref overload). At minimum correct the doc to [in].
- **Verifier correction:** fit() reads its `points` argument as input-only; the `@param[in,out]` annotation on GaussFitter.h:91 is incorrect (should be `[in]`) and the non-const `std::vector&` is unnecessary. This is a low-severity doc/const-correctness inconsistency (the functor stores a const pointer and only reads via const_iterator), not a behavior surprise. Fix: correct the doc to `[in]` and, ideally, change the parameter to `const std::vector<DPosition<2> >&` to match GammaDistributionFitter; that recompiles cleanly for all current callers.

### [MATH-18] GumbelDistributionFitter::fit — fit() documents getGnuplotFormula() but no such method exists, and non-const& 'points' is not modified
`low` · `documentation-api-contract-mismatch` · ABI: `source-compatible` · src/openms/include/OpenMS/MATH/STATISTICS/GumbelDistributionFitter.h · _math-fitters_

```cpp
GumbelDistributionFitResult fit(std::vector<DPosition<2> >& points) const
```
- **Expectation:** Class doc says 'formula ... can be transformed into a gnuplot formula using getGnuplotFormula() after fitting' and the [in,out] points param implies mutation.
- **Actual:** There is no getGnuplotFormula() member anywhere in the class; the gnuplot formula is only emitted to stdout behind a disabled GUMBEL_DISTRIBUTION_FITTER_VERBOSE #ifdef. The non-const points vector is only read by the functor, never written.
- **Evidence:** Header lines 28-29 promise getGnuplotFormula(); class has only fit/fitWeighted/setInitialParameters. GumbelDistributionFitter.cpp:96-122 functor holds `const std::vector<DPosition<2>>* m_data;` and the formula lives inside `#ifdef GUMBEL_DISTRIBUTION_FITTER_VERBOSE`. Header line 66 doc: `@param[in,out] points`.
- **Fix:** Remove the getGnuplotFormula() promise from the class doc (it was likely copied from PosteriorErrorProbabilityModel). Make the points parameter const-ref or relabel [in]. Doc/const-ref change is source-compatible; the const-ref overload signature change is ABI-breaking.
- **Verifier correction:** The finding's underlying facts are correct but its framing is wrong. This is not a "surprising-throw" (the Exception::UnableToFit throw is documented at header line 68). It is documentation/API-contract drift: (a) the class doc (header lines 28-29) advertises a getGnuplotFormula() method that does not exist on GumbelDistributionFitter or any base — the formula is only emitted to stdout inside the disabled #ifdef GUMBEL_DISTRIBUTION_FITTER_VERBOSE block; and (b) @param[in,out] points (header line 66) is wrong because fit() only reads the vector (functor stores it as const*). The same phantom getGnuplotFormula() doc is duplicated in GaussFitter.h and GumbelMaxLikelihoodFitter.h, confirming a copy-paste artifact. Severity is low (compile-error-loud for the phantom method; merely over-restrictive for the [in,out] tag — no silent wrong results, data loss, or crash). Recommended fix is doc-only: drop the getGnuplotFormula() sentence and change @param[in,out] to @param[in]; this is source- and ABI-compatible. An optional change of points to const& would improve correctness but is ABI-breaking and not necessary.

### [MATH-20] PosteriorErrorProbabilityModel::getGumbel_ — Public static getGumbel_ computes a Gumbel density but takes a GaussFitResult, reinterpreting x0/sigma as Gumbel location/scale
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h · _math-fitters_

```cpp
static double getGumbel_(double x, const GaussFitter::GaussFitResult& params)
```
- **Expectation:** A function taking a GaussFitter::GaussFitResult evaluates a Gaussian; getGumbel_ with a trailing underscore reads as a private helper, not public API.
- **Actual:** It is declared public and computes the Gumbel density `(z*exp(-z))/sigma` where z=exp((x0-x)/sigma) — treating the Gauss result's x0 and sigma as the Gumbel location a and scale b. A caller passing an actual fitted Gaussian result gets a Gumbel value, not a Gaussian. The trailing-underscore name signals 'internal' yet it sits in the public section.
- **Evidence:** Header lines 199-204 (public): `static double getGumbel_(double x, const GaussFitter::GaussFitResult & params){ double z = exp((params.x0 - x) / params.sigma); return (z * exp(-1 * z)) / params.sigma; }`. Class-level @note line 38 explains the alpha/beta overloading but the function signature/type does not.
- **Fix:** Move getGumbel_ to the private section (it is an internal helper) or rename without the underscore and accept a GumbelDistributionFitResult to make the parameter semantics explicit. Moving to private is ABI-breaking; renaming is too. A doc note on the public declaration is source-compatible.
- **Verifier correction:** The genuine surprise is only that getGumbel_ is declared public despite a trailing underscore that conventionally marks an internal helper; it is in fact an internal helper with no external callers. The secondary claim — that the function 'reinterprets' a Gauss fit's x0/sigma as Gumbel params and thus silently mis-evaluates — is not a hidden surprise: it is explicitly documented at the class level (line 38), on the backing member variable (lines 256-257), and mirrored by getGumbelGnuplotFormula, where reusing GaussFitResult as a parameter container is an established in-class convention. The function is named and documented as computing a Gumbel density, so no Gaussian is ever promised. Severity is low (mild encapsulation/naming inconsistency, no incorrect results). The current code state has no ABI impact; the recommended fixes (move to private / rename / change parameter type) would be ABI-breaking, while a clarifying doc note on the public declaration would be source-compatible.

### [MATH-22] PosteriorErrorProbabilityModel::getIncorrectlyAssignedFitResult / getIncorrectlyAssignedGumbelFitResult — Doc comments on the incorrectly-assigned getters say they return parameters for 'correctly assigned' sequences (copy-paste)
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h · _math-fitters_

```cpp
GaussFitter::GaussFitResult getIncorrectlyAssignedFitResult() const
```
- **Expectation:** getIncorrectlyAssignedFitResult()'s documentation should describe the incorrect/negative component.
- **Actual:** All three getter doc lines read 'returns estimated parameters for correctly assigned sequences.' — the two incorrectly-assigned getters carry the wrong description, inviting a caller to mix up the positive and negative components (which is exactly the kind of sign/label error that silently corrupts a posterior).
- **Evidence:** Header lines 175-191: getCorrectlyAssignedFitResult, getIncorrectlyAssignedFitResult, and getIncorrectlyAssignedGumbelFitResult each prefixed with `///returns estimated parameters for correctly assigned sequences.`
- **Fix:** Correct the two incorrectly-assigned getters' doc comments to say 'incorrectly assigned'. Doc-only, fully ABI/source compatible.
- **Verifier correction:** The evidence is accurate, but the impact is mild, not corrupting. The two incorrectly-assigned getters' doc comments should read 'returns estimated parameters for incorrectly assigned sequences.' The fix is doc-only and fully ABI/source compatible. Because the function names themselves correctly say 'Incorrectly', the contradictory doc invites at most mild confusion rather than silent posterior corruption, so this is a low-severity documentation defect.

### [MATH-23] GumbelMaxLikelihoodFitter::GumbelMaxLikelihoodFitter(GumbelDistributionFitResult) — Single-argument constructor is non-explicit (and mislabeled 'Default constructor'), allowing implicit GumbelDistributionFitResult -> GumbelMaxLikelihoodFitter conversion
`low` · `implicit-conversion` · ABI: `none` · src/openms/include/OpenMS/MATH/STATISTICS/GumbelMaxLikelihoodFitter.h · _math-fitters_

```cpp
GumbelMaxLikelihoodFitter(GumbelDistributionFitResult init)
```
- **Expectation:** A single-arg converting constructor for a stateful fitter should be explicit; its doc should not call it the 'Default constructor'.
- **Actual:** `GumbelMaxLikelihoodFitter(GumbelDistributionFitResult init)` is not marked explicit, so a GumbelDistributionFitResult will implicitly convert to a fitter wherever one is expected. The Doxygen comment also labels it '/// Default constructor', duplicating the label of the real no-arg constructor.
- **Evidence:** Header lines 53-56: `/// Default constructor` above `GumbelMaxLikelihoodFitter();` and again `/// Default constructor` above `GumbelMaxLikelihoodFitter(GumbelDistributionFitResult init);` — the latter is single-arg and non-explicit.
- **Fix:** Mark the single-arg constructor `explicit` and fix the comment to '/// Constructor with initial parameters'. Adding explicit is source-compatible for intended (direct-init) callers and does not change ABI.
- **Verifier correction:** The header has a real defect: the single-arg converting constructor `GumbelMaxLikelihoodFitter(GumbelDistributionFitResult init)` (line 56) is non-explicit and is mislabeled '/// Default constructor' (line 55), duplicating the label on the genuine no-arg constructor (lines 53-54). It should be `explicit` with a comment like '/// Constructor with initial parameters'. But the implicit-conversion risk is low, not high: GumbelDistributionFitResult itself has only a two-arg (double,double) constructor (no single-arg converting ctor), so no scalar can chain into a fitter, and the fitter is non-copyable (private unimplemented copy ctor/assignment), which blocks pass-by-value implicit conversions. The fix is source-compatible (the sole caller already uses direct brace-init) and ABI-neutral.

### [MATH-24] PosteriorErrorProbabilityModel::fit — fit() sorts and shifts the caller's score vector in place (documented only in a trailing @note)
`low` · `documentation-clarity` · ABI: `none` · src/openms/include/OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h · _math-fitters_

```cpp
bool fit(std::vector<double>& search_engine_scores, const std::string& outlier_handling)
```
- **Expectation:** From the signature `fit(std::vector<double>& scores, ...)`, a caller might expect the vector to be read (and at most rescaled), but the sort-in-place reordering is easy to miss.
- **Actual:** fit() sorts search_engine_scores ascending in place and the smallest value drives a hidden global +fabs(smallest_score_)+0.001 shift used later by computeProbability. The reordering is only mentioned in a final @note ('the vector is sorted from smallest to biggest value!'), so a caller relying on the original order silently gets a permuted vector.
- **Evidence:** PEP cpp:70 `sort(search_engine_scores.begin(), search_engine_scores.end());`, cpp:72 `smallest_score_ = search_engine_scores[0];`, cpp:77 shift. Header lines 110/121 doc the sort only via @note.
- **Fix:** Elevate the in-place sort/shift from a @note to the @param description (it is the primary surprise of the [in,out] contract), or operate on an internal copy. Documentation change is fully compatible; copy-internally is a source-compatible behavior change.
- **Verifier correction:** fit(std::vector<double>& scores, ...) sorts the caller's vector ascending in place (cpp:248); the caller's element ORDER is changed but values are preserved. It does NOT shift the caller's vector — the +fabs(smallest_score_)+0.001 transform is applied only to an internal copy (cpp:252-258) and to the persisted member smallest_score_ later reused in computeProbability (cpp:579). The mutation is already documented both by the @param[in,out] tag and a @note at the declaration (h:107/110), so it is not a hidden side effect. Recommend lifting the 'sorted' fact into the @param line for prominence; severity is low.

### [MATH-1] Math::spline_bisection — "find the maximum point" actually brackets an arbitrary stationary point assuming a hardcoded left-edge derivative sign; can return a minimum
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/MATH/MISC/SplineBisection.h · _math-spline_

```cpp
template <class T> void spline_bisection(const T& peak_spline, double left_neighbor_mz, double right_neighbor_mz, double& max_peak_mz, double& max_peak_int, double threshold = 1e-6)
```
- **Expectation:** Name and doc ("Uses bisection to find the maximum point of a spline", out-params max_peak_mz/max_peak_int) imply it returns the location/value of the maximum of the spline on [left,right].
- **Actual:** It hardcodes `lefthand_sign = true` (assumes derivative at the left edge is positive) and then bisects on the sign change of the first derivative. It finds ANY point where the derivative crosses zero under that assumption — which is a minimum if the spline is concave-up there, and is meaningless if the derivative does not actually change sign in the interval. Nothing verifies the result is a maximum, and the assumed sign is never checked against peak_spline.derivative(left).
- **Evidence:** SplineBisection.h lines 40-69: `bool lefthand_sign = true;` ... `bool midpoint_sign = (midpoint_deriv_val < 0.0) ? false : true; if (lefthand_sign ^ midpoint_sign) { righthand = mid; } else { lefthand = mid; }` ... `max_peak_mz = (lefthand + righthand) / 2; max_peak_int = peak_spline.eval(max_peak_mz);`
- **Fix:** Document the precondition (caller must supply a bracket where the spline has a single interior maximum and positive slope at the left edge), or compute the actual sign via `peak_spline.derivative(left_neighbor_mz)` instead of hardcoding true, and/or rename to spline_find_stationary_point. ABI-safe to add a sign check; renaming is source-breaking.
- **Verifier correction:** The function brackets an interior stationary point of the spline's first derivative under the UNVERIFIED, hardcoded assumption that the slope at the left edge is positive (lefthand_sign=true, never compared to peak_spline.derivative(left_neighbor_mz)). When that precondition holds (the spline rises from the left edge to a single interior maximum), it correctly returns the maximum — which is exactly the regime both in-tree callers use (peak-apex brackets). If the precondition is violated it does NOT reliably return "a minimum"; rather the bracket collapses toward the left edge and yields a non-maximum value, with no check or warning. The genuine defect is the undocumented precondition / unchecked hardcoded sign, not a wrong result in real use. Minimal ABI-safe fix: replace `bool lefthand_sign = true;` with `bool lefthand_sign = peak_spline.derivative(left_neighbor_mz) >= 0.0;` and/or document the precondition; renaming to spline_find_stationary_point is the only source-breaking option.

### [MATH-6] BSplineSmoothingSpline::smoothing_param / num_interior_knots / rss — smoothing_param() returns the resolved/auto-computed s, not the s the caller passed; the accepted RSS may exceed it
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/MATH/MISC/BSplineSmoothingSpline.h · _math-spline_

```cpp
double smoothing_param() const; int num_interior_knots() const; double rss() const;
```
- **Expectation:** A getter named smoothing_param() paired with a constructor parameter `s` reads, to many callers, as "the s I requested". And rss() <= smoothing_param() is the implied scipy-style contract (RSS(f) <= s).
- **Actual:** When s<0 the constructor replaces it with the scipy default m - sqrt(2m) and smoothing_param() returns that resolved value (fine, but undocumented as such). More surprising: the polynomial branch accepts a fit when `rss <= s_target * 1.1` (a 10% over-budget margin), so rss() can be reported larger than smoothing_param() while ok()==true — violating the documented constraint RSS(f) <= s.
- **Evidence:** BSplineSmoothingSpline.cpp lines 48-57 (auto s_), and line 210 `if (rss <= s_target * 1.1) // Allow 10% margin` followed by accept + `ok_ = true`; the header (lines 33-35) advertises the constraint `RSS(f) <= s`.
- **Fix:** Document that smoothing_param() returns the effective (possibly auto) s and that the achieved rss() may exceed s by the internal tolerance, or tighten the acceptance to rss <= s_target. Doc/impl change only.
- **Verifier correction:** The documented contract RSS(f) <= s (header lines 33-35) is violated by the implementation: the polynomial branch accepts a fit when rss <= s_target * 1.1 (BSplineSmoothingSpline.cpp:210), and on acceptance stores the true (possibly over-budget) RSS into rss_ with ok_=true (lines 268-270). Thus rss() can exceed smoothing_param() by up to 10% while ok()==true. The auto-s resolution (s<0 -> m - sqrt(2m), returned by smoothing_param()) is real but IS documented in the header (lines 63,69) and is a normal sentinel-resolution idiom, so it is not itself a surprise. Severity is low: rss() is truthful, ok() is loud, the slack is bounded at 10%, and no current caller relies on the invariant. Fix is doc/impl only: either document the 10% tolerance at the rss()/smoothing_param() getters and in the class contract, or tighten line 210 to rss <= s_target.

### [MATH-12] Math::Histogram::getCumulativeHistogram — Two bare bool flags (complement, inclusive) with inverted/non-obvious meaning at the call site
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/MATH/STATISTICS/Histogram.h · _math-stats_

```cpp
static void getCumulativeHistogram(DataIterator begin, DataIterator end, bool complement, bool inclusive, Histogram& histogram)
```
- **Expectation:** A caller reading `getCumulativeHistogram(b, e, true, false, h)` cannot tell what the booleans mean; 'complement' would normally mean 1 - CDF (survival), and the direction of accumulation is ambiguous.
- **Actual:** complement==true accumulates from the lowest bin UP TO each value (incUntil), complement==false accumulates from each value to the TOP bin (incFrom). The mapping of 'complement' to direction is the opposite of the usual survival-function intuition and is undiscoverable at the call site, especially combined with the second bool.
- **Evidence:** Lines 267-284: `if (complement) { histogram.incUntil(*it, inclusive); } else { histogram.incFrom(*it, inclusive); }`. incUntil increments bins [0..idx] (line 232-244); incFrom increments [idx..end] (line 253-265).
- **Fix:** Replace the bool pair with named enums (e.g. enum class Direction {Ascending, Descending}; enum class Bound {Inclusive, Exclusive}) via an additive overload, and document the survival-vs-cumulative meaning of 'complement'. Additive overload is ABI-safe.
- **Verifier correction:** Two undocumented, adjacent bare bools (complement, inclusive) on getCumulativeHistogram are a real but minor readability/ergonomics concern. However, the mapping is NOT inverted: complement==true → incUntil produces the survival function (count above threshold = 1−CDF), which correctly matches the name "complement of the CDF"; complement==false → incFrom produces the CDF. The sole caller (XFDRAlgorithm.cpp:166) uses complement=true with a comment ("count the number of scores above consecutive thresholds") that confirms the intended, correct semantics. No silently-wrong results. Severity low. A named-enum additive overload plus doxygen would improve clarity and is source-compatible.

### [MATH-13] Math::MultipleTesting::pEmp — pEmp clamps every empirical p-value up to a floor of 1/m0, so it never returns 0 even for a statistic beyond all nulls
`low` · `documentation-gap` · ABI: `none` · src/openms/include/OpenMS/MATH/STATISTICS/MultipleTesting.h · _math-stats_

```cpp
template <class T> static std::vector<double> pEmp(const std::vector<T>& stat, const std::vector<T>& stat0)
```
- **Expectation:** A function returning 'empirical p-values' would be expected to return the raw fraction of null statistics >= the observed value, which can legitimately be 0 when the observed statistic exceeds every null.
- **Actual:** After computing the empirical p-values, every value <= 1/m0 is silently raised to 1/m0. A caller who expects a 0 (or who divides/log-transforms assuming the documented 'empirical p-value' definition) gets a floored value with no mention in the signature or brief doc.
- **Evidence:** Lines 309-311: `const double minp = 1.0 / static_cast<double>(m0); for (auto& vv : out) if (vv <= minp) vv = minp;`. The header doc (lines 174-180) only says 'Compute empirical p-values' with no mention of the 1/m0 floor.
- **Fix:** Document the 1/m0 lower bound in the brief (it mirrors qvalue::empPvals in R, but the floor must be stated). Doc-only; ABI-safe.
- **Verifier correction:** The 1/m0 floor is real and undocumented in the header brief (lines 174-180), but it is the standard, intended behavior of an "empirical p-value" as defined by pyprophet's pemp and R's qvalue::empPvals, which this function explicitly ports (the unit test is named "pyprophet reference vector" and asserts the floored values). It is therefore a documentation-gap, not a silent-failure: the brief should state the 1/m0 lower bound. Severity low; fix doc-only and ABI-safe.

### [ML-1] EuclideanSimilarity::operator()(const std::pair<float,float>&) — Self-similarity documented to 'yield 0' actually returns 1.0
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ML/CLUSTERING/EuclideanSimilarity.h · _ml-clustering-a_

```cpp
float operator()(const std::pair<float, float> & c) const
```
- **Expectation:** Header doc states '@brief calculates self similarity, will yield 0'. A caller reading this would treat the single-argument operator as returning 0 for identical points (i.e. distance-like semantics where 0 == identical).
- **Actual:** The implementation returns operator()(c, c) = 1 - sqrtf(0)/scale_ = 1.0, not 0. The two-argument operator is a SIMILARITY (1 - normalized_distance), so the self-similarity is the maximum value 1.0, directly contradicting the documented 0.
- **Evidence:** EuclideanSimilarity.cpp:32-35  `float EuclideanSimilarity::operator()(const std::pair<float, float> & c) const { return operator()(c, c); }` and EuclideanSimilarity.cpp:45 `return 1 - (sqrtf(...) / scale_);` -> for c==c yields 1. Header EuclideanSimilarity.h:56 `@brief calculates self similarity, will yield 0`.
- **Fix:** Fix the documentation to state the self-similarity yields 1 (maximum similarity). Do not change the return value; the class is a similarity functor and ClusterHierarchical relies on values in [0,1] with 1==identical. ABI-safe doc-only fix.
- **Verifier correction:** The claim is accurate: the single-arg operator()(const std::pair<float,float>&) returns 1.0 for self-similarity (max), not the documented 0, because it forwards to the two-arg similarity operator which returns 1 - normalized_distance and self-distance is 0. Adjustment is only to severity (high -> low): no callers exist in the codebase, the return value is the correct similarity-functor result (not silently corrupt data), and the misleading "@brief ... will yield 0" is contradicted right next to it (line 51 "[0,1]") and by the test (EuclideanSimilarity_test.cpp:67-68 assert 1-0). Recommended fix is doc-only: change "will yield 0" to "will yield the maximum similarity 1 (self-distance is 0)". ABI-safe.

### [ML-3] ClusterFunctor::operator() — cluster_tree output is annotated [in] and is cleared/overwritten by the functor
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ML/CLUSTERING/ClusterFunctor.h · _ml-clustering-a_

```cpp
virtual void operator()(DistanceMatrix<float> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const float threshold = 1) const = 0
```
- **Expectation:** ClusterFunctor.h, CompleteLinkage.h and SingleLinkage.h annotate cluster_tree as '@param[in]'. A caller reading '[in]' would believe the passed vector is read, not modified, and could pass a non-empty vector expecting it preserved or appended.
- **Actual:** Every implementation calls cluster_tree.clear() and fills it from scratch; it is purely an output. Only AverageLinkage.h documents it correctly as [out]. The [in] annotation on the base class and two of three siblings is wrong and hides the clear() side effect.
- **Evidence:** ClusterFunctor.h:63 `@param[in] cluster_tree`; CompleteLinkage.h:55 `@param[in] cluster_tree`; SingleLinkage.h:50 `@param[in] cluster_tree`. SingleLinkage.cpp:46 `cluster_tree.clear();`; AverageLinkage.cpp:50 `cluster_tree.clear();`. AverageLinkage.h:52 correctly uses `@param[out] cluster_tree`.
- **Fix:** Change the Doxygen annotations on the base class and the two siblings from [in] to [out] for cluster_tree. Doc-only, ABI-safe.
- **Verifier correction:** The factual claim is fully accurate. Only the implied severity is reduced: this is a mild, recoverable documentation surprise (an output parameter mislabeled [in]), not a high-impact silent-data-loss hazard, because the void return signature, the parameter name/description, and the correctly-labeled AverageLinkage sibling all make the output role obvious and the empty-vector usage the natural one. Recommended fix unchanged: change [in] to [out] on ClusterFunctor.h:63, SingleLinkage.h:50, and CompleteLinkage.h:55.

### [ML-6] ClusterAnalyzer::averagePopulationAberration / cut — Read-only analyzer methods are non-const and take tree by non-const reference despite not modifying it
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ML/CLUSTERING/ClusterAnalyzer.h · _ml-clustering-a_

```cpp
float averagePopulationAberration(Size cluster_quantity, std::vector<BinaryTreeNode> & tree)
```
- **Expectation:** Analysis methods that only read a clustering tree and distance matrix should be const member functions and should take the tree by const reference, so they can be called on a const ClusterAnalyzer and with const trees. averagePopulationAberration takes `std::vector<BinaryTreeNode>& tree` (non-const), implying it may modify the tree.
- **Actual:** averagePopulationAberration never modifies tree (it builds a local copy of cluster memberships), yet the parameter is a non-const reference and the method is non-const. Likewise averageSilhouetteWidth, dunnIndices, cohesion, cut, and newickTree are all non-const member functions even though they only read. This forces callers to hold a non-const ClusterAnalyzer and a non-const tree unnecessarily.
- **Evidence:** ClusterAnalyzer.h:76 `float averagePopulationAberration(Size cluster_quantity, std::vector<BinaryTreeNode> & tree);` (non-const ref). ClusterAnalyzer.cpp:599-644 only reads `tree`. None of the analyzer methods (lines 46,56,65,93,106,116) are declared const.
- **Fix:** Mark the analyzer methods const and take tree by const reference. This is source-compatible for most callers but technically an ABI/signature change (mangled names change), so schedule it for a minor-version API break or provide const overloads.
- **Verifier correction:** Only averagePopulationAberration takes the tree by non-const reference; the other tree-consuming analyzer methods (averageSilhouetteWidth, dunnIndices, both cut overloads, newickTree) already use const&. The shared, real issue is that all of these read-only methods are declared non-const on a stateless class and therefore cannot be invoked on a const ClusterAnalyzer. Fix = mark the analyzer methods const AND change averagePopulationAberration's tree parameter to const std::vector<BinaryTreeNode>&. This is source-compatible for essentially all callers but changes mangled symbol names, so it is an ABI break to schedule for a minor-version API change. Severity is low: it is a clarity/usability surprise only — no path produces wrong results, data loss, or crashes.

### [ML-7] EuclideanSimilarity::setScale — setScale Doxygen brief is a copy-paste 'clusters the indices...' unrelated to setting a scale
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ML/CLUSTERING/EuclideanSimilarity.h · _ml-clustering-a_

```cpp
void setScale(float x)
```
- **Expectation:** Documentation of setScale should describe configuring the normalization scale used to convert distance to similarity.
- **Actual:** The @brief reads 'clusters the indices according to their respective element distances' (copied from ClusterFunctor::operator()), which describes clustering, not scale-setting. The detail text below it is correct, but the brief actively misdescribes the method.
- **Evidence:** EuclideanSimilarity.h:63-71 `@brief clusters the indices according to their respective element distances ... sets the scale so that similarities can be correctly calculated from distances`.
- **Fix:** Replace the @brief with an accurate description (e.g. 'sets the normalization scale used to map distance to similarity'). Doc-only, ABI-safe.
- **Verifier correction:** The @brief for EuclideanSimilarity::setScale (EuclideanSimilarity.h:64) is a copy-paste from the linkage clustering functors and wrongly reads "clusters the indices according to their respective element distances"; setScale only does scale_ = x. Real but low severity (trivial setter; correct detail text follows; no silent wrong results). Doc-only, ABI-safe.

### [ML-11] ClusteringGrid::getIndex — getIndex throws IllegalArgument for out-of-grid positions but documents no throw
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/ML/CLUSTERING/ClusteringGrid.h · _ml-clustering-b_

```cpp
CellIndex getIndex(const Point &position) const
```
- **Expectation:** A pure coordinate-to-index mapping getter ('@return cell index (i,j) of the cell in which (x,y) lies') reads as total and exception-free; a caller would expect to be able to query the index for any position.
- **Actual:** If the position lies outside [range_x_, range_y_] the method throws Exception::IllegalArgument. The header's Doxygen for getIndex (lines 91-97) lists no @throw at all, so the throw is undocumented and surprising — a stark mirror image of isNonEmptyCell, whose documented throws don't happen here but happen there.
- **Evidence:** ClusteringGrid.cpp:73-78 throws `Exception::IllegalArgument` with message 'This position (x,y)=... is outside the range of the grid.'; header lines 91-97 contain no @throw.
- **Fix:** Document the @throw Exception::IllegalArgument on getIndex (doc-only, ABI-safe). Align the throw/no-throw contracts between getIndex and isNonEmptyCell so they are not inverted.
- **Verifier correction:** getIndex (ClusteringGrid.cpp:73-78) throws Exception::IllegalArgument for out-of-grid positions but its header Doxygen (lines 91-97) lists no @throw. The inversion with isNonEmptyCell is real but in the opposite direction from a hazard standpoint: isNonEmptyCell documents IllegalArgument + InvalidValue throws (lines 105-106) yet its body (cpp:86-89) never throws. The throw in getIndex is loud and clearly messaged, so this is a documentation-accuracy issue (undocumented throw + over-documented non-throwing sibling), not a silent-data hazard. Fix is doc-only.

### [ML-13] HashGrid::grid_dimension — grid_dimension is a public const reference that silently mutates and never shrinks (high-water mark)
`low` · `missing-documentation` · ABI: `none` · src/openms/include/OpenMS/ML/CLUSTERING/HashGrid.h · _ml-clustering-b_

```cpp
const CellIndex& grid_dimension;
```
- **Expectation:** A public member named grid_dimension documented as 'Upper-right corner of key space for cells' reads as the current bounding extent of occupied cells. A caller caching it expects it to track insert/erase, and a const reference reads as stable.
- **Actual:** grid_dimension is a const reference aliasing the private grid_dimension_, which updateGridDimension_ only ever raises (monotonic max) on insert and NEVER lowers on erase/clear. After erasing the extreme element, grid_dimension still reports the historical maximum, not the current extent. Also, being a public reference to a member that mutates under the caller is surprising for a `const` reference.
- **Evidence:** HashGrid.h:286 `const CellIndex& grid_dimension;`; updateGridDimension_ HashGrid.h:488-497 only does `if (*it_cur < *it_new) *it_cur = *it_new;`; erase (309-327) and clear (332-335) never touch grid_dimension_.
- **Fix:** Document grid_dimension as a monotonic high-water mark that is not reduced by erase/clear, and that it is the largest cell index ever inserted (not current extent). Prefer exposing it via an accessor rather than a public const reference. Doc-only fix is ABI-safe; replacing the public reference with a getter is breaking.
- **Verifier correction:** grid_dimension is a monotonic high-water mark: it records the largest (per-axis max) cell index ever inserted and is never reduced by erase or clear, so after erasing the extreme element it still reports the historical maximum rather than the current occupied extent — contradicting its doc "Upper-right corner of key space for cells." This is a documentation gap, not a hidden side effect of a misused const reference (binding a public const reference to a live private member is a standard idiom and the member is intended to grow on insert). It is low severity: no production caller reads grid_dimension (QTClusterFinder, the sole consumer, never uses it), so it cannot produce wrong results, data loss, or crashes — it only invites a mild misread for a hypothetical current-extent query. Fix: document the monotonic/never-shrinks-on-erase semantics (doc-only, ABI-safe).

### [ML-14] HashGrid::erase(iterator) — erase(iterator) leaves empty cells in the grid while insert/size assume cleanup; asymmetric with map semantics
`low` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/ML/CLUSTERING/HashGrid.h · _ml-clustering-b_

```cpp
void erase(iterator pos)
```
- **Expectation:** Modeling 'most parts of the C++ standard map interface' (header line 29), a developer expects erase to leave the container in a clean state and grid-cell bookkeeping consistent (empty cells removed, as ClusteringGrid::removeCluster does).
- **Actual:** erase(iterator) only erases the element from the cell's multimap; it never removes a now-empty CellContent from cells_. Empty cells accumulate; grid_begin()/grid_end() will still visit them, and size() iterates over them. The sibling ClusteringGrid::removeCluster explicitly removes the cell when empty, so the convention is inconsistent across the module's two grid classes.
- **Evidence:** HashGrid.h:309-313 `void erase(iterator pos){ CellContent& cell = pos.grid_it_->second; cell.erase(pos.cell_it_); }` (no empty-cell removal). Contrast ClusteringGrid.cpp:49-59 which erases the cell when its list becomes empty.
- **Fix:** Document that empty cells persist after erase (and that grid iteration may visit empty cells), or remove the cell when it becomes empty to match ClusteringGrid. Documentation fix is ABI-safe; the behavior change could affect iterator stability so gate carefully.
- **Verifier correction:** erase(iterator) (and also erase(const key_type&)) leaves now-empty cells in cells_, unlike the sibling ClusteringGrid::removeCluster which removes empty cells. The real, provable consequence is limited to the grid-cell-level introspection API (grid_begin/grid_end/grid_find/grid_at), which will visit/return stale empty cells. The element-level map interface is NOT affected: size(), empty(), and begin()/end() correctly account for and skip empty cells (verified by HashGrid_test.cpp:68-76 where size() returns 0 after erase). The claim's assertions that "size() iterates over them" causes wrong counts and that "insert assumes cleanup" are incorrect/unsubstantiated. Severity is low (mild surprise to a developer doing raw grid-cell iteration and expecting ClusteringGrid-style non-empty cells); no data loss/crash. Recommended fix: document that empty cells persist after erase and grid-level iteration may visit them, or remove the cell when empty to match ClusteringGrid (documentation fix is ABI-safe).

### [ML-22] ROCCurve::AUC / rocN / curve / cutoffPos / cutoffNeg — Query/getter methods are non-const and mutate internal state (sort + overwrite pos_/neg_)
`low` · `hidden-side-effect` · ABI: `breaking` · src/openms/include/OpenMS/ML/ROCCURVE/ROCCurve.h · _ml-interp-solvers_

```cpp
double AUC(); double rocN(Size N); std::vector<std::pair<double,double>> curve(UInt resolution = 10);
```
- **Expectation:** Methods that read out a score (AUC, rocN, curve, cutoff*) read-only summaries of the inserted data; a developer would expect them const and idempotent.
- **Actual:** All of them call the private sort() which reorders score_clas_pairs_ in place, and AUC() additionally overwrites the member counters via 'pos_ = truePos; neg_ = falsePos;'. So computing AUC() mutates pos_/neg_ that other accessors then rely on, and none can be called on a const ROCCurve. The reordering also makes the 'first encountered' semantics of cutoff* order-dependent.
- **Evidence:** void ROCCurve::sort() { if (!sorted_) { std::sort(...); sorted_ = true; } } ... double ROCCurve::AUC() { ... pos_ = truePos; neg_ = falsePos; return area; }
- **Fix:** Make these methods const with the sort cache stored in mutable members (sorted_ is already a member), and stop overwriting pos_/neg_ inside AUC() (compute locals). ABI-breaking if signatures change to const; the additive minimum is to remove the pos_/neg_ reassignment in AUC().
- **Verifier correction:** Confirmed: AUC/rocN/curve/cutoffPos/cutoffNeg are non-const and mutate score_clas_pairs_ via the private sort() (in-place std::sort), so none can be called on a const ROCCurve and they reorder inserted data. Corrected: the AUC() reassignment pos_ = truePos; neg_ = falsePos writes back values identical to the already-maintained counters (truePos/falsePos are full totals after the loop), so it is redundant, NOT a source of silently-wrong results for other accessors; the cutoff* ordering is deterministic because sort() runs first. Additionally, the class is explicitly documented at ROCCurve.h:24 as "[This class is buggy and usage is discouraged!]", which tempers the developer-expectation framing. The recommendation to mark these const with mutable cache is reasonable but signature-changing (ABI-breaking).

### [ML-24] CrossValidation::gridSearch1D — Doxygen for gridSearch1D describes a bool 'prefer_larger' parameter that does not exist
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/ML/CROSSVALIDATION/CrossValidation.h · _ml-interp-solvers_

```cpp
gridSearch1D(..., CandidateTieBreak tie_break = ...)
```
- **Expectation:** Parameter documentation should match the actual parameter list so a caller knows how to request a tie-break policy.
- **Actual:** The @brief tie-breaking note says 'choose by @p prefer_larger (true -> larger wins)', but there is no prefer_larger parameter -- the actual parameter is the enum CandidateTieBreak tie_break. A caller reading the doc would look for a missing bool and could pass arguments to the wrong slot.
- **Evidence:** Doc: 'If |score - best_score| <= @p tie_tol, choose by @p prefer_larger (true -> larger wins).' Signature: 'CandidateTieBreak tie_break = CandidateTieBreak::PreferLarger'.
- **Fix:** Fix the doc to reference tie_break / CandidateTieBreak. Pure documentation change, source- and ABI-compatible.
- **Verifier correction:** The Doxygen "Tie-breaking" note on line 103 references a non-existent @p prefer_larger with bool semantics ("true -> larger wins"), contradicting the actual parameter `CandidateTieBreak tie_break` (a 3-valued scoped enum) which is itself correctly documented one line below (line 115). The fix is to update line 103 to reference @p tie_break / CandidateTieBreak (e.g. "choose by @p tie_break: PreferLarger keeps the larger candidate, PreferSmaller the smaller, PreferAny keeps the first"). Severity is low rather than the implied higher harm: the param list is correct, the contradiction is confined to a stale prose note, and scoped-enum typing prevents silently passing a bool to the wrong slot. Pure documentation change; source- and ABI-compatible.
