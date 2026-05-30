// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Svetlana Kutuzova, Douglas McCloskey $
// $Authors: Svetlana Kutuzova, Douglas McCloskey $
// --------------------------------------------------------------------------
 
#pragma once
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabM.h>
#include <OpenMS/FORMAT/MzTabFile.h>

namespace OpenMS
{
  /**
    @brief Data processing pipeline for one Flow Injection Analysis MS (FIA-MS) sample.

    FIA-MS omits chromatographic separation; the entire sample is delivered in a short
    time window. This class consumes one @ref MSExperiment together with a target time
    cut-off and performs:
      -# truncate the experiment at @p n_seconds (drops spectra at @c RT @c >= @c n_seconds);
      -# sum the surviving spectra along the time axis with a per-band sliding bin size
         @c bin_size[i] @c = @c mzs_[i] @c / @c (resolution @c * @c 4);
      -# smooth the summed spectrum with @ref SavitzkyGolayFilter;
      -# pick peaks with @ref PeakPickerHiRes;
      -# run @ref AccurateMassSearchEngine against the configured database / adduct list
         and write the result into the caller's @c MzTab object.

    A cached picked spectrum is reused on subsequent runs of the same sample
    when @c load_cached_spectrum is @c true, skipping the cut/merge/pick steps.

    Batches of FIA-MS samples are driven by @ref FIAMSScheduler which loops over CSV
    rows and calls @ref run once per @c ;-separated time window per row.

    The workflow is inspired by Fuhrer et al. (https://pubs.acs.org/doi/10.1021/ac201267k);
    this is not an exact reimplementation.

    Configuration is via the @ref DefaultParamHandler parameter sections @c "" (root),
    @c "db:", @c "sgf:", and @c "sne:" — see the constructor for the full default map.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI FIAMSDataProcessor  :
    public DefaultParamHandler
  {
public:
    /**
      @brief Construct with the FIA-MS default parameters.

      Registers the following parameters with their defaults:
      @c filename @c = @c "fiams", @c dir_output @c = @c "", @c resolution @c = @c 120000,
      @c polarity @c in @c {positive, negative} (default @c positive),
      @c max_mz @c = @c 1500, @c bin_step @c = @c 20,
      @c db:mapping @c = @c ["CHEMISTRY/HMDBMappingFile.tsv"],
      @c db:struct @c = @c ["CHEMISTRY/HMDB2StructMapping.tsv"],
      @c positive_adducts @c = @c "CHEMISTRY/PositiveAdducts.tsv" (advanced),
      @c negative_adducts @c = @c "CHEMISTRY/NegativeAdducts.tsv" (advanced),
      @c store_progress @c in @c {true, false} (default @c true),
      @c sgf:frame_length @c = @c 11, @c sgf:polynomial_order @c = @c 4,
      @c sne:window @c = @c 10.
    */
    FIAMSDataProcessor();

    /// Default destructor.
    ~FIAMSDataProcessor() override = default;

    /// Copy constructor.
    FIAMSDataProcessor(const FIAMSDataProcessor& cp) = default;

    /// Copy assignment.
    FIAMSDataProcessor& operator=(const FIAMSDataProcessor& fdp) = default;

    /**
      @brief Run the full FIA-MS pipeline for one time window and (re)compute the accurate-mass-search result.

      Pipeline:
        -# Compute the cached picked-spectrum path as @c {dir_output}/{filename}_picked_{int(n_seconds)}.mzML.
        -# If @p load_cached_spectrum is @c true and the path exists, load it via
           @c FileHandler with @c MZML pinned and take its first spectrum as the picked
           spectrum (skip steps 3-5).
        -# Otherwise: call @ref cutForTime → @ref mergeAlongTime → @ref extractPeaks.
        -# If @c store_progress is @c "true" (and we recomputed in step 3), store the
           merged spectrum to @c {dir_output}/{filename}_merged_{int(n_seconds)}.mzML and
           the picked spectrum to the cached path.
        -# Estimate signal-to-noise via @ref trackNoise and convert the picked spectrum to
           a @c FeatureMap via @ref convertToFeatureMap.
        -# **Always** store the signal-to-noise spectrum to
           @c {dir_output}/{filename}_signal_to_noise_{int(n_seconds)}.mzML — note that
           this is independent of the @c store_progress setting.
        -# Call @ref runAccurateMassSearch and write @c {dir_output}/{filename}_{int(n_seconds)}.mzTab
           via @c MzTabFile::store.

      Progress is logged via @c OPENMS_LOG_INFO (cache load start/finish or
      calculation start/finish).

      @param[in]  experiment           Input MS experiment (centroided MS1 is assumed).
      @param[in]  n_seconds            Time cut-off (seconds). Cast to @c int when forming the per-run filename suffix, so sub-second offsets collapse.
      @param[out] output               Receives the accurate-mass-search results from @ref runAccurateMassSearch.
      @param[in]  load_cached_spectrum If @c true and the cached picked-spectrum file exists, reuse it instead of re-running the cut / merge / pick steps. Default @c true.

      @return @c true when the picked spectrum was loaded from the cached file, @c false when it was recomputed.
    */
    bool run(const MSExperiment& experiment, const float n_seconds, OpenMS::MzTab& output, const bool load_cached_spectrum = true);

    /**
      @brief Append spectra of @p experiment whose retention time is strictly less than @p n_seconds to @p output.

      @p output is @em appended to (existing entries are preserved). The RT comparison
      is strict (@c < ), so spectra at exactly @c n_seconds are dropped.

      @param[in]     experiment Source experiment.
      @param[in]     n_seconds  Upper retention-time bound (exclusive).
      @param[in,out] output     Receives the retained spectra.
    */
    void cutForTime(const MSExperiment& experiment, const float n_seconds, std::vector<MSSpectrum>& output);

    /**
      @brief Sum @p input across the time axis into a single spectrum using a per-band sliding bin size.

      For each consecutive pair @c (mzs_[i], mzs_[i+1]), calls
      @c SpectrumAddition::addUpSpectra(input, bin_sizes_[i], @c false) and keeps only the
      peaks whose m/z falls in @c [mzs_[i], mzs_[i+1]). Iteration stops one short of the
      end (@c i @c < @c mzs_.size() @c - @c 1), so the upper-most band is not emitted —
      configure @c max_mz to extend past your highest expected peak. The returned
      spectrum is sorted by m/z position.

      @param[in] input Spectra to sum (typically the output of @ref cutForTime).
      @return Single summed spectrum, sorted by m/z.
    */
    MSSpectrum mergeAlongTime(const std::vector<OpenMS::MSSpectrum>& input);

    /**
      @brief Smooth @p input with the configured @ref SavitzkyGolayFilter and pick peaks via @ref PeakPickerHiRes.

      The filter and picker are class members configured through @c sgf:* and the
      picker's own parameter section; their state is set up by the constructor and
      refreshed by @ref updateMembers_ when @c sgf:* parameters change. @p input is not
      modified — a local copy is filtered in place before picking.

      @param[in] input Smoothed/picked input spectrum.
      @return Spectrum of picked peaks.
    */
    MSSpectrum extractPeaks(const MSSpectrum& input);

    /**
      @brief Wrap @p input's peaks as @c Feature objects in a new @c FeatureMap, tagging the @c "scan_polarity" meta value.

      For every peak in @p input, emits one @c Feature with the peak's m/z and intensity
      and a @c "scan_polarity" meta value taken from the @c polarity parameter
      (@c "positive" or @c "negative"). The returned map is unsorted.

      @param[in] input Picked spectrum.
      @return FeatureMap with one feature per input peak; @c scan_polarity meta value set on every feature.
    */
    FeatureMap convertToFeatureMap(const MSSpectrum& input);

    /**
      @brief Build a parallel "noise level" spectrum: per-peak m/z carried over, intensity replaced by the local noise estimate.

      Wraps @ref SignalToNoiseEstimatorMedianRapid with window @c sne:window. Returns an
      empty spectrum when @p input is empty. The output's @c intensity field holds the
      noise estimate at the corresponding m/z (not a S/N ratio).

      @param[in] input Picked spectrum.
      @return Spectrum with the same m/z values as @p input, intensities replaced by per-peak noise estimates.
    */
    MSSpectrum trackNoise(const MSSpectrum& input);

    /**
      @brief Run @ref AccurateMassSearchEngine on @p input and write the matches into @p output.

      Configures @c AccurateMassSearchEngine with:
        - @c ionization_mode @c = @c "auto" (polarity is detected from the feature map, not taken from the class's @c polarity parameter);
        - @c mass_error_value @c = @c 1e6 @c / @c (resolution @c * @c 2) ppm;
        - @c db:mapping / @c db:struct / @c positive_adducts / @c negative_adducts forwarded from the class parameters;
        - @c keep_unidentified_masses @c = @c "false" — only identified masses are reported.
      Then calls @c init followed by @c run.

      @param[in,out] input  Feature map to search; consumed by @c AccurateMassSearchEngine::run.
      @param[out]    output Accurate-mass-search results in @c MzTab form.
    */
    void runAccurateMassSearch(FeatureMap& input, OpenMS::MzTab& output);

    /**
      @brief Return the band centres used by @ref mergeAlongTime.

      The list contains @c {bin_step, 2*bin_step, ..., (n-1)*bin_step} where
      @c n @c = @c max_mz @c / @c bin_step.

      @return Const reference to the vector of band centres. Refreshed by @ref updateMembers_ when parameters change.
    */
    const std::vector<float>& getMZs();

    /**
      @brief Return the per-band sliding bin sizes used by @ref mergeAlongTime.

      Each bin size is computed as @c bin_size[i] @c = @c mzs_[i] @c / @c (resolution @c * @c 4).

      @return Const reference to the vector of sliding bin sizes, parallel to @ref getMZs.
    */
    const std::vector<float>& getBinSizes();

protected:
    /**
      @brief Recompute the per-band @c mzs_ / @c bin_sizes_ caches and push the @c sgf:* parameters into the @ref SavitzkyGolayFilter member.

      Invoked by @ref DefaultParamHandler when any parameter changes (notably @c max_mz,
      @c bin_step, @c resolution, @c sgf:frame_length, @c sgf:polynomial_order).
    */
    void updateMembers_() override;

private:
    /**
      @brief Write @p input to @p filename as a single-spectrum mzML file via @c FileHandler.

      Used by @ref run to persist the merged spectrum, picked spectrum, and signal-to-noise
      spectrum.

      @param[in] input    Spectrum to store.
      @param[in] filename Output path for the mzML file.

      @throws Exception Propagates exceptions from @c FileHandler::storeExperiment on I/O failure (no in-class guard).
    */
    void storeSpectrum_(const MSSpectrum& input, const String& filename);

    std::vector<float> mzs_;        ///< Per-band m/z centres consumed by @ref mergeAlongTime; populated by @ref updateMembers_.
    std::vector<float> bin_sizes_;  ///< Per-band sliding bin sizes parallel to @ref mzs_; populated by @ref updateMembers_.
    SavitzkyGolayFilter sgfilter_;  ///< Smoothing filter used by @ref extractPeaks; configured from @c sgf:* parameters.
    PeakPickerHiRes picker_;        ///< Peak picker used by @ref extractPeaks; configured from the picker's own parameter section.
  };
} // namespace OpenMS
