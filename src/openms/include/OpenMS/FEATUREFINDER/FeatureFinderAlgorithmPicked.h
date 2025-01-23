// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/FEATUREFINDER/TraceFitter.h>
#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <fstream>

namespace OpenMS
{
  /**@brief The purpose of this struct is to provide definitions of classes and typedefs which are used throughout all FeatureFinder classes.  */
  struct OPENMS_DLLAPI FeatureFinderDefs
  {
    /// Index to peak consisting of two UInts (scan index / peak index)
    typedef IsotopeCluster::IndexPair IndexPair;

    /// Index to peak consisting of two UInts (scan index / peak index) with charge information
    typedef IsotopeCluster::ChargedIndexSet ChargedIndexSet;

    /// A set of peak indices
    typedef IsotopeCluster::IndexSet IndexSet;

    /// Flags that indicate if a peak is already used in a feature
    enum Flag {UNUSED, USED};

    /// Exception that is thrown if a method an invalid IndexPair is given
    class OPENMS_DLLAPI NoSuccessor :
      public Exception::BaseException
    {
public:
      NoSuccessor(const char * file, int line, const char * function, const IndexPair & index) :
        BaseException(file, line, function, "NoSuccessor", String("there is no successor/predecessor for the given Index: ") + String(index.first) + "/" + String(index.second)),
        index_(index)
      {
        Exception::GlobalExceptionHandler::setMessage(what());
      }

      ~NoSuccessor() noexcept override = default;

protected:
      IndexPair index_; // index without successor/predecessor
    };
  };


  /**
    @brief FeatureFinderAlgorithm for picked peaks.

    @htmlinclude OpenMS_FeatureFinderAlgorithmPicked.parameters

    @improvement RT model with tailing/fronting (Marc)
    @improvement More general MZ model - e.g. based on co-elution or with sulfur-averagines (Marc)

    @todo Fix output in parallel mode, change assignment of charges to threads, add parallel TOPP test (Marc)
    @todo Implement user-specified seed lists support (Marc)

    @ingroup FeatureFinder
  */
  class OPENMS_DLLAPI FeatureFinderAlgorithmPicked :
    public DefaultParamHandler, public ProgressLogger
  {
public:
    /// @name Type definitions
    //@{
    typedef MSExperiment MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef SpectrumType::FloatDataArrays FloatDataArrays;
    //@}

protected:
    typedef Peak1D PeakType;
    typedef FeatureFinderAlgorithmPickedHelperStructs::Seed Seed;
    typedef FeatureFinderAlgorithmPickedHelperStructs::MassTrace MassTrace;
    typedef FeatureFinderAlgorithmPickedHelperStructs::MassTraces MassTraces;
    typedef FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern TheoreticalIsotopePattern;
    typedef FeatureFinderAlgorithmPickedHelperStructs::IsotopePattern IsotopePattern;

public:
    /// default constructor
    FeatureFinderAlgorithmPicked();

    void setSeeds(const FeatureMap& seeds);

    void setData(const MSExperiment& map, FeatureMap& features);

    void run(PeakMap& input_map, 
      FeatureMap& features, 
      const Param& param, 
      const FeatureMap& seeds);

    virtual Param getDefaultParameters() const
    {
      return defaults_;
    }
protected:
    void run();

    /// editable copy of the map
    MapType map_;

    FeatureMap* features_;

    /// Output stream for log/debug info
    mutable std::ofstream log_;
    /// debug flag
    bool debug_;
    /// Array of abort reasons
    std::map<String, UInt> aborts_;
    /// Array of abort reasons
    std::map<Seed, String> abort_reasons_;
    /// User-specified seed list
    FeatureMap seeds_;

    /// @name Members for parameters often needed in methods
    //@{
    double pattern_tolerance_; ///< Stores mass_trace:mz_tolerance
    double trace_tolerance_; ///< Stores isotopic_pattern:mz_tolerance
    UInt min_spectra_; ///< Number of spectra that have to show the same mass (for finding a mass trace)
    UInt max_missing_trace_peaks_; ///< Stores mass_trace:max_missing
    double slope_bound_; ///< Max slope of mass trace intensities
    double intensity_percentage_; ///< Isotope pattern intensity contribution of required peaks
    double intensity_percentage_optional_; ///< Isotope pattern intensity contribution of optional peaks
    double optional_fit_improvement_; ///< Minimal improvement for leaving out optional isotope
    double mass_window_width_; ///< Width of the isotope pattern mass bins
    UInt intensity_bins_; ///< Number of bins (in RT and MZ) for intensity significance estimation
    double min_isotope_fit_; ///< Minimum isotope pattern fit for a feature
    double min_trace_score_; ///< Minimum quality of a traces
    double min_rt_span_; ///< Minimum RT range that has to be left after the fit
    double max_rt_span_; ///< Maximum RT range the model is allowed to span
    double max_feature_intersection_; ///< Maximum allowed feature intersection (if larger, that one of the feature is removed)
    String reported_mz_; ///< The mass type that is reported for features. 'maximum' returns the m/z value of the highest mass trace. 'average' returns the intensity-weighted average m/z value of all contained peaks. 'monoisotopic' returns the monoisotopic m/z value derived from the fitted isotope model.
    //@}

    /// @name Members for intensity significance estimation
    //@{
    /// RT bin width
    double intensity_rt_step_;
    /// m/z bin width
    double intensity_mz_step_;
    /// Precalculated intensity 20-quantiles (binned)
    std::vector<std::vector<std::vector<double> > > intensity_thresholds_;
    //@}

    ///Vector of precalculated isotope distributions for several mass windows
    std::vector<TheoreticalIsotopePattern> isotope_distributions_;

    // Docu in base class
    void updateMembers_() override;

    /// Writes the abort reason to the log file and counts occurrences for each reason
    void abort_(const Seed& seed, const String& reason);

    /**
     * Calculates the intersection between features.
     * The value is normalized by the size of the smaller feature, so it ranges from 0 to 1.
     */
    double intersection_(const Feature& f1, const Feature& f2) const;

    /// Returns the isotope distribution for a certain mass window
    const TheoreticalIsotopePattern& getIsotopeDistribution_(double mass) const;

    /**
      @brief Finds the best fitting position of the isotopic pattern estimate defined by @p center

      @param center the maximum peak of the isotope distribution (contains charge as well)
      @param charge The charge of the pattern
      @param best_pattern Returns the indices of the isotopic peaks. If a isotopic peak is missing -1 is returned.
    */
    double findBestIsotopeFit_(const Seed& center, UInt charge, IsotopePattern& best_pattern) const;

    /**
      Extends all mass traces of an isotope pattern in one step

      @param pattern The IsotopePattern that should be extended.
      @param traces The MassTraces datastructure where the extended mass traces will be stored in.
      @param meta_index_overall The index of the data array where the quality scores for the given charge are stored.
    */
    void extendMassTraces_(const IsotopePattern& pattern, MassTraces& traces, Size meta_index_overall) const;

    /**
      @brief Extends a single mass trace in one RT direction

      How to use this method:
        - Add the starting peak to the @p trace
        - Indicate using @c increase_rt whether to extend in downstream or upstream direction

      @param trace The trace that should be extended
      @param spectrum_index The index of the spectrum from which on the mass trace should be extended
      @param mz The mz location (center) of the trace
      @param increase_rt Indicator whether the extension is done in forward or backward direction (with respect to the current spectrum)
      @param meta_index_overall The index of the overall score
      @param min_rt The rt minimum up to which the trace will be extended.
      @param max_rt The rt maximum up to which the trace will be extended.

      @note This method assumes that it extends from a local maximum.
      @note If @c min_rt or @c max_rt are set to 0.0 no boundary is assumed in the respective direction.
    */
    void extendMassTrace_(MassTrace& trace, SignedSize spectrum_index, double mz, bool increase_rt, Size meta_index_overall, double min_rt = 0.0, double max_rt = 0.0) const;

    /// Returns the index of the peak nearest to m/z @p pos in spectrum @p spec (linear search starting from index @p start)
    Size nearest_(double pos, const MSSpectrum& spec, Size start) const;

    /**
      @brief Searches for an isotopic peak in the current spectrum and the adjacent spectra

      @param pos m/z position of the searched for peak
      @param spectrum_index index of the central spectrum
      @param pattern IsotopePattern to store found peaks
      @param pattern_index index of the isotope in the pattern
      @param peak_index starting index of the search (to avoid multiple binary searches)
    */
    void findIsotope_(double pos, Size spectrum_index, IsotopePattern& pattern, Size pattern_index, Size& peak_index) const;

    /// Calculates a score between 0 and 1 for the m/z deviation of two peaks.
    double positionScore_(double pos1, double pos2, double allowed_deviation) const;

    /// Calculates a score between 0 and 1 for the correlation between theoretical and found isotope pattern
    double isotopeScore_(const TheoreticalIsotopePattern& isotopes, IsotopePattern& pattern, bool consider_mz_distances) const;

    /**
      @brief Compute the intensity score for the peak @p peak in spectrum @p spectrum.

      The intensity score is computed by interpolating the score between the 4 nearest intensity
      bins. The scores from the different bins are weighted by the distance of the bin center to
      the peak.

      @param spectrum Index of the spectrum we are currently looking at
      @param peak Index of the peak that should be scored inside the spectrum @p spectrum
    */
    double intensityScore_(Size spectrum, Size peak) const;

    /**
      @brief Choose a the best trace fitter for the current mass traces based on the user parameter
             (symmetric, asymmetric) or based on an inspection of the mass trace (auto)

      @return A pointer to the trace fitter that should be used.
     */
    std::unique_ptr<TraceFitter> chooseTraceFitter_(double& tau);

    double intensityScore_(Size rt_bin, Size mz_bin, double intensity) const;

    /**
      @name Handling of fitted mass traces

      Methods to handle the results of the mass trace fitting process.
    */
    //@{

    /**
      @brief Creates new mass traces @p new_traces based on the fitting result and the
      original traces @p traces.

      @param fitter The TraceFitter containing the results from the rt profile fitting step.
      @param traces Original mass traces found in the experiment.
      @param new_traces Mass traces created by cropping the original mass traces.
     */
    void cropFeature_(const std::shared_ptr<TraceFitter>& fitter,
                      const MassTraces& traces,
                      MassTraces& new_traces);

    /**
      @brief Checks the feature based on different score thresholds and model constraints

      Feature can get invalid for following reasons:
      <ul>
        <li>Invalid fit: Fitted model is bigger than 'max_rt_span'</li>
        <li>Invalid feature after fit - too few traces or peaks left</li>
        <li>Invalid fit: Center outside of feature bounds</li>
        <li>Invalid fit: Less than 'min_rt_span' left after fit</li>
        <li>Feature quality too low after fit</li>
      </ul>

      @param fitter The TraceFitter containing the results from the rt profile fitting step.
      @param feature_traces Cropped feature mass traces.
      @param seed_mz Mz of the seed
      @param min_feature_score Minimal required feature score
      @param error_msg Will be filled with the error message, if the feature is invalid
      @param fit_score Will be filled with the fit score
      @param correlation Will be filled with correlation between feature and model
      @param final_score Will be filled with the final score

      @return true if the feature is valid
     */
    bool checkFeatureQuality_(const std::shared_ptr<TraceFitter>& fitter,
                              MassTraces& feature_traces,
                              const double& seed_mz, const double& min_feature_score,
                              String& error_msg, double& fit_score, double& correlation, double& final_score);

    /**
      @brief Creates several files containing plots and viewable data of the fitted mass trace

      @param fitter The TraceFitter containing the results from the rt profile fitting step.
      @param traces Original mass traces found in the spectra
      @param new_traces Cropped feature mass traces
      @param feature_ok Status of the feature
      @param error_msg If the feature is invalid, @p error_msg contains the reason
      @param final_score Final score of the feature
      @param plot_nr Index of the feature
      @param peak The Seed Peak
      @param path The path where to put the debug files (default is debug/features)
    */
    void writeFeatureDebugInfo_(const std::shared_ptr<TraceFitter>& fitter,
                                const MassTraces& traces,
                                const MassTraces& new_traces,
                                bool feature_ok, const String& error_msg, const double final_score, const Int plot_nr, const PeakType& peak,
                                const String& path  = "debug/features/");

    //@}
private:

    /// Not implemented
    FeatureFinderAlgorithmPicked& operator=(const FeatureFinderAlgorithmPicked&);
    /// Not implemented
    FeatureFinderAlgorithmPicked(const FeatureFinderAlgorithmPicked&);
  };

} // namespace OpenMS

